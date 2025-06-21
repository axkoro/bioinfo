#pragma once

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

typedef struct Sequence {
  std::string header;
  std::string data;
} Sequence;

std::vector<Sequence> read_fasta(const std::string& path);

template <class T>
class Matrix {
 private:
  size_t num_rows;
  size_t num_cols;
  std::vector<T> data;

 public:
  Matrix(size_t rows, size_t cols) : num_rows(rows), num_cols(cols), data(rows * cols) {}

  T& operator()(size_t row, size_t col) { return data[row * num_cols + col]; }

  const T& operator()(size_t row, size_t col) const { return data[row * num_cols + col]; }

  size_t rows() const { return num_rows; }
  size_t cols() const { return num_cols; }

  // for debugging
  void print() const {
    if (num_rows == 0 || num_cols == 0) {
      std::cout << "[]" << std::endl;
      return;
    }

    std::vector<int> col_widths(num_cols, 0);

    for (size_t j = 0; j < num_cols; ++j) {
      for (size_t i = 0; i < num_rows; ++i) {
        std::ostringstream oss;
        oss << (*this)(i, j);
        std::string s = oss.str();
        if (static_cast<int>(s.length()) > col_widths[j]) {
          col_widths[j] = s.length();
        }
      }
    }

    for (size_t i = 0; i < num_rows; ++i) {
      for (size_t j = 0; j < num_cols; ++j) {
        std::cout << std::right << std::setw(col_widths[j] + 2) << (*this)(i, j);
      }
      std::cout << std::endl;
    }
  }
};