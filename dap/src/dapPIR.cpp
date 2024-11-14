#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <set>

using namespace Rcpp;
using namespace std;

//' Find all model combinations that have a proposal density greater than or equal to a threshold
//' @param mat A p x L matrix of proposal densities from SuSiE
//' @param threshold The threshold value of the proposal density
//' @return A NumericMatrix containing the unique combinations
//' @export
// [[Rcpp::export]]
NumericMatrix pir(NumericMatrix mat, double threshold) {
  int p = mat.nrow();
  int L = mat.ncol();
  vector<int> indices(L, 0);

  // Take log10 of the matrix
  NumericMatrix logMat(p, L);
  for (int i = 0; i < p; i++) {
    for (int j = 0; j < L; j++) {
      logMat(i, j) = log10(mat(i, j));
    }
  }

  // Remove elements smaller than the threshold from each vector and keep track of original indices
  vector<vector<pair<double, int>>> sortedArr(L);
  for (int i = 0; i < L; i++) {
    for (int j = 0; j < p; j++) {
      if (logMat(j, i) >= log10(threshold)) {
        sortedArr[i].emplace_back(logMat(j, i), j);
      }
    }
    sort(sortedArr[i].begin(), sortedArr[i].end(), greater<pair<double, int>>());
  }

  vector<vector<int>> resultMatrix;
  set<vector<int>> uniqueCombinations;

  while (true) {
    double sum = 0.0;
    bool valid = true;

    // Compute the summation of the current combination
    for (int i = 0; i < L; i++) {
      sum += sortedArr[i][indices[i]].first;
      if (sum < log10(threshold)) {
        valid = false;
        break;
      }
    }

    // If the summation is greater than the threshold, save the combination as a p-vector
    if (valid && sum >= log10(threshold)) {
      vector<int> combinationVector(p, 0);
      for (int i = 0; i < L; i++) {
        combinationVector[sortedArr[i][indices[i]].second] = 1;
      }
      uniqueCombinations.insert(combinationVector);
    }

    // Find the rightmost vector that has more elements left after the current element in that vector
    int next = L - 1;
    while (next >= 0 && (indices[next] + 1 >= sortedArr[next].size())) {
      next--;
    }

    // No such vector is found, so no more combinations left
    if (next < 0) {
      break;
    }

    // Move to the next element in that vector
    indices[next]++;

    // For all vectors to the right of this vector, reset the index to the first element
    for (int i = next + 1; i < L; i++) {
      indices[i] = 0;
    }
  }

  // Add null row and single entry rows
  for (int i = 0; i < p; i++) {
    vector<int> singleEntryRow(p, 0);
    singleEntryRow[i] = 1;
    uniqueCombinations.insert(singleEntryRow);
  }

  // Convert uniqueCombinations to resultMatrix
  vector<vector<int>> uniqueResultMatrix(uniqueCombinations.begin(), uniqueCombinations.end());

  // Convert uniqueResultMatrix to NumericMatrix
  NumericMatrix result(uniqueResultMatrix.size(), p);
  for (size_t i = 0; i < uniqueResultMatrix.size(); i++) {
    for (int j = 0; j < p; j++) {
      result(i, j) = uniqueResultMatrix[i][j];
    }
  }

  return result;
}