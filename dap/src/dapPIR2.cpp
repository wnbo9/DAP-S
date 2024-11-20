#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <set>
#include <fstream>

using namespace Rcpp;
using namespace std;

//' Find all model combinations that have a proposal density greater than or equal to a threshold
//' @param mat A p x L matrix of proposal densities from SuSiE
//' @param threshold The threshold value of the proposal density
//' @return A NumericMatrix containing the unique combinations
//' @export
// [[Rcpp::export]]
NumericMatrix pir2(NumericMatrix mat, double threshold) {
  int p = mat.nrow();
  int L = mat.ncol();
  double logThreshold = log10(threshold);

  // open output file
  ofstream outFile("output.txt");
  // Take log10 of the matrix, and add a row of -9999 to avoid out of bound error
  NumericMatrix logMat(p+1, L);
  for (int i = 0; i < p; i++) {
    for (int j = 0; j < L; j++) {
      logMat(i, j) = log10(mat(i, j));
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
    vector<int> combination(L, 0); // keep track of model name
    int pos = -1;

    outFile << "---------------" << endl;
    // Generate the current combination
    for (int i = 0; i < L; ++i) {
        currentSum += sortedMat[i][indices[i]].first;
        combination[i] = sortedMat[i][indices[i]].second;
        if (currentSum < logThreshold) {
            pos = i;
            outFile << "The pos is " << pos << endl;
            break;
        }
    }

    // Print the current combination
    outFile << "Sum: " << currentSum << " -> Model Name: ";
    for (int i = 0; i < L; ++i) {
        outFile << combination[i] + 1 << " ";
    }
    outFile << endl;
    
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
            outFile << "Moving left 1 column" << endl;
            indices[pos - 1]++;
            outFile << "indices[pos - 1] is: indice[" << pos - 1 << "] = " << indices[pos - 1] << " <-> " << sortedMat[pos - 1][indices[pos - 1]].second + 1 << endl;
            fill(indices.begin() + pos, indices.end(), 0);
            increment = true;
        } else { // the second column selects the first row, then we finish exploring
            done = true;
            break;
        }
    } else if (pos > 1) { // the third column or later column already < threshold
        outFile << "indices[pos] is: " << indices[pos] << endl;
        if (indices[pos] > 0) { // the pos column does not select the first row, then we go to the previous column and indices ++
            outFile << "Moving left 1 column" << endl;
            indices[pos - 1]++;
            outFile << "indices[pos - 1] is: indice[" << pos - 1 << "] = " << indices[pos - 1] << " <-> " << sortedMat[pos - 1][indices[pos - 1]].second + 1 << endl;
            fill(indices.begin() + pos, indices.end(), 0);
            increment = true;
        } else { // the pos column selects the first row, then we go to the previous previous column and indices ++
            outFile << "Moving left 2 column" << endl;
            indices[pos - 2]++;
            outFile << "indices[pos-2] is: indice[" << pos - 2 << "] = " << indices[pos - 2] << " <-> " << sortedMat[pos - 2][indices[pos - 2]].second + 1 << endl;
            fill(indices.begin() + pos - 1, indices.end(), 0);
            increment = true;
        }
    } else if (pos == -1) { // current model > threshold, then we move to the next model
        for (int i = L - 1; i >= 0; --i) {
            indices[i]++;
            outFile << "indices[i] is: indice[" << i << "] = " << indices[i] << " <-> " << sortedMat[i][indices[i]].second + 1 << endl;
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
    outFile << "" << endl;
  }


  // Add null row and single entry rows
  for (int i = 0; i < p; i++) {
    vector<int> singleEntryRow(L, i);
    Combinations.insert(singleEntryRow);
  }

  // Convert Combinations to resultMatrix with unique rows
  set<vector<int>> binaryCombinations;
  for (const auto& combination : Combinations) {
    vector<int> binaryVector(p, 0);
    for (int name : combination) {
      binaryVector[name] = 1; // Convert the name to a binary entry
    }
    binaryCombinations.insert(binaryVector);
  }

  // Convert binaryCombinations to NumericMatrix
  NumericMatrix result(binaryCombinations.size(), p);
  size_t i = 0; // Index for rows of the result matrix
  for (const auto& binaryVector : binaryCombinations) {
    for (size_t j = 0; j < p; j++) {
        result(i, j) = binaryVector[j];
    }
    i++;
  }
  outFile.close();
  return result;
}
