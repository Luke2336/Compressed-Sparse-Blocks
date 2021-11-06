#include <iostream>
#include <vector>

#include "csb.hpp"

int main() {
  size_t N, M; // the size of matrix is N x M
  std::cin >> N >> M;
  std::vector<std::vector<double>> Matrix(N, std::vector<double>(M));
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < M; ++j) {
      std::cin >> Matrix[i][j];
    }
  }

  CSB CSB1(Matrix);
  std::vector<double> X(M), Y(N);
  for (size_t i = 0; i < M; ++i) {
    std::cin >> X[i];
  }
  CSB1.SpMV(X, Y);
  for (size_t i = 0; i < N; ++i) {
    std::cout << Y[i] << " \n"[i == N - 1];
  }
  return 0;
}