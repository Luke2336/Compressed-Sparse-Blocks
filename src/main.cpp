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
  return 0;
}