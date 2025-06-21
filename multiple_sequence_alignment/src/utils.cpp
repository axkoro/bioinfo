#include "utils.hpp"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

std::vector<Sequence> read_fasta(const std::string& path) {
  std::ifstream file(path);
  if (!file.is_open()) {
    throw std::runtime_error("Could not open file at " + path);
  }

  std::vector<Sequence> sequences;
  Sequence current_sequence;
  std::string line;

  while (std::getline(file, line)) {
    if (line.empty()) {
      continue;
    }

    // Remove potential carriage return from Windows-style line endings
    if (line.back() == '\r') {
      line.pop_back();
    }

    if (line[0] == '>') {
      if (!current_sequence.data.empty()) {
        sequences.push_back(current_sequence);
      }

      // Start a new sequence
      current_sequence = Sequence();
      size_t first_char_pos = line.find_first_not_of(" \t", 1);  // skips whitespace after '>'
      current_sequence.header = line.substr(first_char_pos);
    } else {
      current_sequence.data += line;
    }
  }

  // After the loop, add the very last sequence in the file
  if (!current_sequence.data.empty()) {
    sequences.push_back(current_sequence);
  }

  return sequences;
}