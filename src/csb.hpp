#pragma once

#include <math.h>

#include <vector>

class CSB {
private:
  const std::vector<std::vector<double>> Original_Matrix;
  const size_t N;
  const size_t M;
  const size_t Beta;
  std::vector<double> Val;
  std::vector<size_t> Row_Idx;
  std::vector<size_t> Col_Idx;
  std::vector<size_t> Blk_Ptr;

  size_t genBeta(size_t N) {
    // Beta is power of 2
    size_t sq = sqrt(N);
    size_t ret = 1;
    while (ret <= sq) {
      ret <<= 1;
    }
    return ret >> 1;
  }

  bool notZero(double Value) { return Value < -1e-6 || Value > 1e-6; }

  void genZMorton(size_t Row_Begin, size_t Row_End, size_t Col_Begin,
                  size_t Col_End, size_t Dim) {
    if (Row_Begin >= Row_End || Col_Begin >= Col_End) {
      return;
    }
    if (Dim == 1) {
      if (notZero(Original_Matrix[Row_Begin][Col_Begin])) {
        Val.emplace_back(Original_Matrix[Row_Begin][Col_Begin]);
        Row_Idx.emplace_back(Row_Begin);
        Col_Idx.emplace_back(Col_Begin);
      }
      return;
    }
    size_t Half_Dim = Dim >> 1;
    // Top-left
    genZMorton(Row_Begin, Row_Begin + Half_Dim, Col_Begin, Col_Begin + Half_Dim,
               Half_Dim);
    // Top-right
    genZMorton(Row_Begin, Row_Begin + Half_Dim, Col_Begin + Half_Dim, Col_End,
               Half_Dim);
    // Bottom-left
    genZMorton(Row_Begin + Half_Dim, Row_End, Col_Begin, Col_Begin + Half_Dim,
               Half_Dim);
    // Bottom-right
    genZMorton(Row_Begin + Half_Dim, Row_End, Col_Begin + Half_Dim, Col_End,
               Half_Dim);
  }

public:
  CSB(const std::vector<std::vector<double>> &Original_Matrix)
      : Original_Matrix(Original_Matrix), N(Original_Matrix.size()),
        M(Original_Matrix[0].size()), Beta(genBeta(Original_Matrix.size())) {
    size_t NumRowBlock = N / Beta;
    size_t NumColBlock = M / Beta;
    Blk_Ptr.emplace_back(0);
    for (size_t i = 0; i < NumRowBlock; ++i) {
      size_t Row_Begin = i * Beta;
      size_t Row_End = std::min(N, Row_Begin + Beta);
      for (size_t j = 0; j < NumColBlock; ++j) {
        size_t Col_Begin = j * Beta;
        size_t Col_End = std::min(M, Col_Begin + Beta);
        genZMorton(Row_Begin, Row_End, Col_Begin, Col_End, Beta);
        Blk_Ptr.emplace_back(Val.size());
      }
    }
  }
};