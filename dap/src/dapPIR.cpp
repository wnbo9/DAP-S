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
  int p = mat.nrow();
  int L = mat.ncol();
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
  bool hasNull = false;

  while (!done) {
    double currentSum = 0.0;
    vector<int> combination(L, -1); // keep track of model name
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
      bool isNull = true;
      for (int i = 0; i < L; i++) {
        if (combination[i] != -1){
          isNull = false;
          break;
        }
      }
      if (isNull){
        hasNull = true;
      }
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

std::vector<std::set<int>> positionElements(L);
  for (const auto &combination : Combinations) {
    for (int col = 0; col < L; col++) {
      if (combination[col] != -1) {
        positionElements[col].insert(combination[col]);
      }
    }
  }
   List positionList(L);
  for (int pos = 0; pos < L; pos++) {
    IntegerVector colElements;
    for (auto val : positionElements[pos]) {
      colElements.push_back(val);
    }
    positionList[pos] = colElements;
  }


  // Add null row and single entry rows
  for (int i = 0; i < (hasNull ? p-1 : p); i++) {
    vector<int> singleEntryRow(L, i);
    Combinations.insert(singleEntryRow);
  }


  


  // Convert Combinations to resultMatrix with unique rows
  set<vector<int>> binaryCombinations;
  for (const auto& combination : Combinations) {
    vector<int> binaryVector(p, 0);
    for (int name : combination) {
      if (name >= 0){
        binaryVector[name] = 1; // Convert the name to a binary entry
      }
    }
    binaryCombinations.insert(binaryVector);
  }
  
  // Convert binaryCombinations to NumericMatrix
  NumericMatrix result(binaryCombinations.size(), p-1); // Exclude the null row
  size_t i = 0; // Index for rows of the result matrix
  for (const auto& binaryVector : binaryCombinations) {
    for (size_t j = 0; j < p-1; j++) {
        result(i, j) = binaryVector[j];
    }
    i++;
  }
  
  //return result; // m*p binary matrix of model configurations
  return List::create(
    Named("combinations") = result,
    Named("position_elements") = positionList
  );
}
