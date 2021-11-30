#pragma once

#include <math.h>

#include <iostream>
#include <vector>

class CSB {
private:
  const std::vector<std::vector<double>> Original_Matrix;
  const size_t N;
  const size_t M;
  const size_t Beta;
  const size_t NumRowBlock;
  const size_t NumColBlock;
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
        Row_Idx.emplace_back(Row_Begin % Beta);
        Col_Idx.emplace_back(Col_Begin % Beta);
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

  size_t blockIdx(int i, int j) { return i * NumColBlock + j; }

  void blockRowV(const size_t i, std::vector<int>::iterator R_Begin,
                 std::vector<int>::iterator R_End,
                 std::vector<double>::iterator X_Begin,
                 std::vector<double>::iterator X_End,
                 std::vector<double>::iterator Y_Begin,
                 std::vector<double>::iterator Y_End, std::vector<double> &X) {
    int R_Len = R_End - R_Begin;
    if (R_Len == 2) {
      size_t l = *(R_Begin) + 1;
      size_t r = *(R_Begin + 1);
      if (l == r) {
        size_t Start = Blk_Ptr.at(blockIdx(i, l));
        size_t End = Blk_Ptr.at(blockIdx(i, r) + 1) - 1;
        blockV(Start, End, Beta, X_Begin, X_End, Y_Begin, Y_End);
      } else {
        if (Blk_Ptr.at(blockIdx(i, l)) >= Blk_Ptr.at(blockIdx(i, r) + 1)) {
          return;
        }
        size_t Start = Blk_Ptr.at(blockIdx(i, l));
        size_t End = Blk_Ptr.at(blockIdx(i, r) + 1) - 1;
        for (size_t k = Start; k <= End; ++k) {
          int X_Idx = (upper_bound(Blk_Ptr.begin(), Blk_Ptr.end(), k) -
                       Blk_Ptr.begin() - 1) %
                      NumColBlock;
          *(Y_Begin + Row_Idx.at(k)) +=
              Val.at(k) * X.at(X_Idx * Beta + Col_Idx.at(k));
        }
      }
      return;
    }
    size_t Mid = (R_Len & (size_t)1) ? (R_Len >> 1) : ((R_Len >> 1) - 1);
    size_t XMid = Beta * (*(R_Begin + Mid) - *R_Begin);
    std::vector<double> Z(Y_End - Y_Begin, 0);
    blockRowV(i, R_Begin, R_Begin + Mid + 1, X_Begin, X_Begin + XMid, Y_Begin,
              Y_End, X);
    blockRowV(i, R_Begin + Mid, R_End, X_Begin + XMid, X_End, Z.begin(),
              Z.end(), X);
    for (size_t k = 0; k < Z.size(); ++k) {
      *(Y_Begin + k) += Z.at(k);
    }
  }

  void blockV(const size_t Start, const size_t End, const size_t Dim,
              std::vector<double>::iterator X_Begin,
              std::vector<double>::iterator X_End,
              std::vector<double>::iterator Y_Begin,
              std::vector<double>::iterator Y_End) {
    if (Start > End) {
      return;
    }
    if (End - Start <= Dim) {
      for (size_t k = Start; k <= End; ++k) {
        *(Y_Begin + Row_Idx.at(k)) += Val.at(k) * *(X_Begin + Col_Idx.at(k));
      }
      return;
    }
    size_t Half_Dim = Dim >> 1;
    // Binary search s2
    size_t s2 = End + 1;
    size_t Low = Start, High = End;
    while (Low <= High) {
      size_t Mid = (Low + High) >> 1;
      if (Row_Idx.at(Mid) & Half_Dim) {
        High = Mid - 1;
        s2 = Mid;
      } else {
        Low = Mid + 1;
      }
    }
    // Binary search s1
    size_t s1 = s2;
    Low = Start, High = s2 - 1;
    while (Low <= High) {
      size_t Mid = (Low + High) >> 1;
      if (Col_Idx.at(Mid) & Half_Dim) {
        High = Mid - 1;
        s1 = Mid;
      } else {
        Low = Mid + 1;
      }
    }
    // Binary search s3
    size_t s3 = End + 1;
    Low = s2, High = End;
    while (Low <= High) {
      size_t Mid = (Low + High) >> 1;
      if (Col_Idx.at(Mid) & Half_Dim) {
        High = Mid - 1;
        s3 = Mid;
      } else {
        Low = Mid + 1;
      }
    }
    blockV(Start, s1 - 1, Half_Dim, X_Begin, X_End, Y_Begin, Y_End);
    blockV(s3, End, Half_Dim, X_Begin, X_End, Y_Begin, Y_End);
    blockV(s1, s2 - 1, Half_Dim, X_Begin, X_End, Y_Begin, Y_End);
    blockV(s2, s3 - 1, Half_Dim, X_Begin, X_End, Y_Begin, Y_End);
  }

public:
  CSB(const std::vector<std::vector<double>> &Original_Matrix)
      : Original_Matrix(Original_Matrix), N(Original_Matrix.size()),
        M(Original_Matrix.at(0).size()), Beta(genBeta(N)),
        NumRowBlock(N / Beta + (N % Beta != 0)),
        NumColBlock(M / Beta + (M % Beta != 0)) {
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

  void SpMV(std::vector<double> &X, std::vector<double> &Y) {
    fill(Y.begin(), Y.end(), 0);
    for (size_t i = 0; i < NumRowBlock; ++i) {
      std::vector<int> R;
      R.emplace_back(-1);
      size_t count = 0;
      for (size_t j = 0; j < NumColBlock - 1; ++j) {
        size_t Block_Idx = blockIdx(i, j);
        count += Blk_Ptr.at(Block_Idx + 1) - Blk_Ptr.at(Block_Idx);
        if (count + Blk_Ptr.at(Block_Idx + 2) - Blk_Ptr.at(Block_Idx + 1) >
            Beta) {
          R.emplace_back(j);
          count = 0;
        }
      }
      R.emplace_back(NumColBlock - 1);
      blockRowV(i, R.begin(), R.end(), X.begin(), X.end(), Y.begin() + i * Beta,
                std::min(Y.end(), Y.begin() + (i + 1) * Beta), X);
    }
  }
};