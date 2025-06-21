#pragma once

#include <string>
#include <unordered_map>

#include "utils.hpp"

class ScoringMatrix {
 public:
  ScoringMatrix(const std::string& path);

  int score(char row_char, char col_char) const;
  void print() const;

 private:
  Matrix<int> matrix;
  std::vector<char> alphabet;
  std::unordered_map<char, size_t> char_to_index;

  size_t get_index_for_char(char c) const;
  void load_from_stream(std::istream& stream);
};

int needleman_wunsch(const std::string& seq1, const std::string& seq2, const ScoringMatrix& mat);

void compute_guidetree(const std::vector<Sequence>& sequences, const ScoringMatrix& mat);
