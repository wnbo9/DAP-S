#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <set>

using namespace Rcpp;
using namespace std;

//' Find all model combinations that have a proposal density greater than or equal to a threshold
//' @param mat A processed p x l matrix of proposal densities from SuSiE
//' @param threshold The threshold value of the proposal density
//' @return A NumericMatrix containing the unique combinations
//' @export
// [[Rcpp::export]]
List pir(NumericMatrix mat, double threshold) {
  int p = mat.nrow(); // 5001
  int L = mat.ncol(); // 3
  double logThreshold = log10(threshold);

  // open output file
  // Take log10 of the matrix, replace 0 with 10^-999 and add a row of -9999 to avoid out of bound error
  NumericMatrix logMat(p+1, L);
  for (int i = 0; i < p; i++) {
    for (int j = 0; j < L; j++) {
      if (mat(i, j) == 0) {
        logMat(i, j) = -999;
      } else {
        logMat(i, j) = log10(mat(i, j));
      }
    }
  }
  for (int j = 0; j < L; j++) {
    logMat(p, j) = -9999;
  }

  // Sort each column in descending order
  vector<vector<pair<double, int>>> sortedMat(L);
  for (int i = 0; i < L; i++) {
    for (int j = 0; j < p + 1; j++) {
      sortedMat[i].emplace_back(logMat(j, i), j);
    }
    sort(sortedMat[i].begin(), sortedMat[i].end(), greater<pair<double, int>>());
  }

  vector<int> indices(L, 0);
  set<vector<int>> Combinations;
  bool done = false;

  while (!done) {
    double currentSum = 0.0;
    vector<int> combination(L, p-1); // keep track of model name
    int pos = -1;

    // Generate the current combination
    for (int i = 0; i < L; ++i) {
        currentSum += sortedMat[i][indices[i]].first;
        // Only set combination[i] if the index is not p (the null SNP)
        if (sortedMat[i][indices[i]].second < p-1){
            combination[i] = sortedMat[i][indices[i]].second;
        }
        if (currentSum < logThreshold) {
            pos = i;
            break;
        }
    }
    
    // If the sum is greater than the threshold, save the combination
    if (currentSum >= logThreshold) {
      Combinations.insert(combination);
    }

    // Handle stopping condition
    bool increment = false;
    if (pos == 0){ // the first column already < threshold, we finish exploring
        done = true;
        break;
    } else if (pos == 1) { // the second column already < threshold
        if (indices[pos] > 0) {  // the second column does not select the first row, then we go to the previous column and indices ++
            indices[pos - 1]++;
            fill(indices.begin() + pos, indices.end(), 0);
            increment = true;
        } else { // the second column selects the first row, then we finish exploring
            done = true;
            break;
        }
    } else if (pos > 1) { // the third column or later column already < threshold
        if (indices[pos] > 0) { // the pos column does not select the first row, then we go to the previous column and indices ++
            indices[pos - 1]++;
            fill(indices.begin() + pos, indices.end(), 0);
            increment = true;
        } else { // the pos column selects the first row, then we go to the previous previous column and indices ++
            indices[pos - 2]++;
            fill(indices.begin() + pos - 1, indices.end(), 0);
            increment = true;
        }
    } else if (pos == -1) { // current model > threshold, then we move to the next model
        for (int i = L - 1; i >= 0; --i) {
            indices[i]++;
            if (indices[i] < p + 1) {
                increment = true;
                break;
            } else {
                indices[i] = 0;
            }
        }
    }

    if (!increment) {
        done = true;
    }
  }

  // Check for duplicates and track single-SNP/NULL models simultaneously
  vector<vector<int>> updatedCombinations;
  std::set<int> existingSingleSNPs;
  bool hasNullModel = false;

  for (const auto& row : Combinations) {
    // Check for duplicates
    bool isValid = true;
    unordered_set<int> seen;
    int nonNullCount = 0;
    int singleSNP = -1;

    // Process each value in the row
    for (int val : row) {
      // Check for duplicates excluding p-1=5000
      if (val != p - 1 && seen.count(val)) {
        isValid = false;
        break;
      }
      seen.insert(val);

      // Count non-NULL SNPs and track single SNP
      if (val != p - 1) {
        nonNullCount++;
        singleSNP = val;
      }
    }

    // If row is valid (no duplicates), add it and track single-SNP models
    if (isValid) {
      updatedCombinations.push_back(row);

      // Track single-SNP models and NULL model
      if (nonNullCount == 1) {
        existingSingleSNPs.insert(singleSNP);
      } else if (nonNullCount == 0) {
        hasNullModel = true;
      }
    }
  }
  
  // Convert Combinations to NumericMatrix
  NumericMatrix combo(updatedCombinations.size(), L);
  for (int i = 0; i < updatedCombinations.size(); i++) {
    for (int col = 0; col < L; col++) {
      combo(i, col) = updatedCombinations[i][col];
    }
  }

  // Create missing 1-SNP models and NULL model
  vector<int> missingModels;
  if (!hasNullModel) {
    missingModels.push_back(p - 1);
  }
  for (int snp = 0; snp < p - 1; snp++) {
    if (!existingSingleSNPs.count(snp)) {
      missingModels.push_back(snp);
    }
  }
  NumericMatrix missingModelMatrix(0, 1);
  if (!missingModels.empty()) {
    missingModelMatrix = NumericMatrix(missingModels.size(), 1);
    for (int i = 0; i < missingModels.size(); i++) {
      missingModelMatrix(i, 0) = missingModels[i];
    }
  }

  //return result; // m*p binary matrix of model configurations
  return List::create(
    Named("combo") = combo,
    Named("single") = missingModelMatrix
  );
}
