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
vector<vector<int>> pir(const vector<vector<double>>& mat, double threshold) {

  int p = mat.size(); // 5001 including NULL SNP
  int L = mat[0].size(); // 3

  vector<vector<int>> updatedCombinations; // Store the updated combinations
  set<vector<int>> uniqueModels; // Sort the unique configurations

  // NULL model
  vector<int> nullModel(1, p - 1);
  uniqueModels.insert(vector<int>());
  updatedCombinations.push_back(nullModel); 
  
  // include all 1-SNP models and return if L == 1
  if (L == 1) {
    for (int snp = 0; snp < p - 1; snp++) {
      vector<int> model = {snp};
      uniqueModels.insert(model);
      updatedCombinations.push_back(model);
    }
    return updatedCombinations;
  }

  // L > 1, we do PIR
  double logThreshold = log10(threshold);

  // Take log10 of the matrix, replace 0 with 10^-999 and add a row of -9999 to avoid out of bound error
  vector<vector<double>> logMat(p + 1, vector<double>(L));
  for (int i = 0; i < p; i++) {
    for (int j = 0; j < L; j++) {
      logMat[i][j] = (mat[i][j] == 0) ? -999 : log10(mat[i][j]);
    }
  }
  fill(logMat[p].begin(), logMat[p].end(), -9999);


  // Sort each column in descending order
  vector<vector<pair<double, int>>> sortedMat(L);
  for (int i = 0; i < L; i++) {
    for (int j = 0; j < p + 1; j++) {
      sortedMat[i].emplace_back(logMat[j][i], j);
    }
    sort(sortedMat[i].begin(), sortedMat[i].end(), greater<pair<double, int>>());
  }

  vector<int> indices(L, 0);
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
        // Extract and sort active SNPs
        vector<int> activeSNPs;
        bool has_duplicate = false;
        set<int> seen;

        for (int val : combination) {
            if (val != p - 1) {
              if (!seen.insert(val).second) { // If the model has same SNPs, discard it
                has_duplicate = true;
                break;
              }
              activeSNPs.push_back(val);
            }
        }

        if (!has_duplicate) {
            // Only add if this set of SNPs is new
            sort(activeSNPs.begin(), activeSNPs.end());
            if (uniqueModels.insert(activeSNPs).second) {
                updatedCombinations.push_back(combination);
            }
        }
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

  for (int snp = 0; snp < p - 1; snp++) {
    vector<int> singleton = {snp};
    // If the singleton model is not already in the set, add it
    if (uniqueModels.insert(singleton).second) {
      uniqueModels.insert(singleton);
      updatedCombinations.push_back(singleton);
    }
  }


  // return results; m*p matrix of combinations and m*1 matrix of missing models; index starts from 0.
  return updatedCombinations;
}