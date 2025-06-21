#include "align.hpp"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "utils.hpp"

// Scoring Matrix

ScoringMatrix::ScoringMatrix(const std::string& path) : matrix(0, 0), alphabet(), char_to_index() {
  std::ifstream file(path);
  if (!file.is_open()) {
    throw std::runtime_error("Could not open file at " + path);
  }
  load_from_stream(file);
}

int ScoringMatrix::score(char row_char, char col_char) const {
  return matrix(get_index_for_char(row_char), get_index_for_char(col_char));
}

void ScoringMatrix::print() const { matrix.print(); };

size_t ScoringMatrix::get_index_for_char(char c) const {
  auto it = char_to_index.find(c);
  if (it == char_to_index.end()) {
    std::stringstream ss;
    ss << "Character '" << c << "' does not have an entry in the scoring matrix.";
    throw std::out_of_range(ss.str());
  }
  return it->second;
}

void ScoringMatrix::load_from_stream(std::istream& stream) {
  std::string line;
  size_t line_number = 0;

  // Find and parse the header line
  while (std::getline(stream, line)) {
    line_number++;
    if (line.empty() || line.front() == '#') {
      continue;
    }

    std::istringstream header_stream(line);
    std::string token;
    while (header_stream >> token) {
      if (token.length() > 1) {
        throw std::runtime_error("Header must contain single characters. Found '" + token +
                                 "' on line " + std::to_string(line_number) + ".");
      }
      char c = token.front();
      auto [_, inserted] = char_to_index.emplace(c, alphabet.size());
      if (!inserted) {
        throw std::runtime_error("Duplicate character '" + std::string(1, c) +
                                 "' in header on line " + std::to_string(line_number) + ".");
      }
      alphabet.push_back(c);
    }
    break;
  }

  if (alphabet.empty()) {
    throw std::runtime_error("A valid, non-empty header line was not found.");
  }

  // Initialize matrix and load score rows
  matrix = Matrix<int>(alphabet.size(), alphabet.size());

  for (size_t i = 0; i < alphabet.size(); ++i) {
    if (!std::getline(stream, line)) {
      throw std::runtime_error("Unexpected end of file. Expected " +
                               std::to_string(alphabet.size()) + " data rows, but found only " +
                               std::to_string(i) + ".");
    }
    line_number++;

    std::istringstream row_stream(line);
    char row_char;

    row_stream >> row_char;
    if (row_stream.fail()) {
      throw std::runtime_error("Failed to read row character on line " +
                               std::to_string(line_number) + ".");
    }

    const size_t row_idx = get_index_for_char(row_char);

    for (size_t col_idx = 0; col_idx < alphabet.size(); ++col_idx) {
      int score_value;
      row_stream >> score_value;
      if (row_stream.fail()) {
        throw std::runtime_error("Invalid or missing score on line " + std::to_string(line_number) +
                                 " for row '" + std::string(1, row_char) + "'.");
      }
      matrix(row_idx, col_idx) = score_value;
    }
  }
}

// Algorithms

int needleman_wunsch(const std::string& seq1, const std::string& seq2, const ScoringMatrix& mat) {
  const size_t num_cols = seq1.size();
  const size_t num_rows = seq2.size();

  Matrix<int> similarity_matrix(num_rows + 1, num_cols + 1);

  // Initialize first row and column
  for (size_t i = 0; i <= num_rows; i++) {
    similarity_matrix(i, 0) = i;
  }
  for (size_t i = 0; i <= num_cols; i++) {
    similarity_matrix(0, i) = i;
  }

  // Compute similarity matrix
  for (size_t row = 1; row <= num_rows; ++row) {
    for (size_t col = 1; col <= num_cols; ++col) {
      char a = seq1.at(col - 1);
      char b = seq2.at(row - 1);

      const int diagonal_score = similarity_matrix(row - 1, col - 1) + mat.score(a, b);
      const int down_score = similarity_matrix(row - 1, col) + mat.score(b, '*');
      const int right_score = similarity_matrix(row, col - 1) + mat.score(a, '*');

      int score = std::max({diagonal_score, down_score, right_score});
      similarity_matrix(row, col) = score;
    }
  }

  return similarity_matrix(num_rows, num_cols);
}

void compute_guidetree(const std::vector<Sequence>& sequences, const ScoringMatrix& mat) {
  size_t num_sequences = sequences.size();
  if (num_sequences <= 1) return;

  // Compute and print the initial pairwise similarity matrix
  Matrix<int> sim_matrix(num_sequences, num_sequences);
  std::cout << std::left << std::setw(8) << "";
  for (const auto& seq : sequences) {
    std::cout << std::setw(8) << seq.header;
  }
  std::cout << std::endl;

  for (size_t i = 0; i < num_sequences; i++) {
    std::cout << std::setw(8) << sequences[i].header;
    for (size_t j = 0; j < num_sequences; j++) {
      if (j < i) {
        std::cout << std::setw(8) << "";
      } else if (i == j) {
        sim_matrix(i, j) = std::numeric_limits<double>::min();  // to prevent self-merging
        std::cout << std::setw(8) << "";
      } else {  // j > i (right of matrix diagonal)
        int sim_score = needleman_wunsch(sequences[i].data, sequences[j].data, mat);
        sim_matrix(i, j) = sim_score;
        sim_matrix(j, i) = sim_score;
        std::cout << std::setw(8) << sim_score;
      }
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  // Hierarchical clustering
  std::vector<std::string> node_labels(num_sequences);
  std::vector<bool> active(num_sequences, true);

  // Initialize node labels with sequence numbers
  for (size_t i = 0; i < num_sequences; i++) {
    node_labels[i] = std::to_string(i + 1);
  }

  size_t active_count = num_sequences;
  while (active_count > 1) {  // Continue until only one cluster remains
    // Find the maximum similarity (smallest distance) among active nodes
    double max_sim = std::numeric_limits<double>::min();
    size_t best_i = 0, best_j = 0;

    for (size_t i = 0; i < num_sequences; i++) {
      if (!active[i]) continue;
      for (size_t j = i + 1; j < num_sequences; j++) {
        if (!active[j]) continue;
        if (sim_matrix(i, j) > max_sim) {
          max_sim = sim_matrix(i, j);
          best_i = i;
          best_j = j;
        }
      }
    }

    // if (max_sim == std::numeric_limits<double>::min()) break;

    // Output the merge in alphabetical order of labels
    std::string label_i = node_labels[best_i];
    std::string label_j = node_labels[best_j];

    if (label_i > label_j) {
      std::swap(label_i, label_j);
      std::swap(best_i, best_j);
    }

    std::cout << "(" << label_i << "," << label_j << ")" << std::endl;

    // Merge i+j -> i (and deactivate j)
    std::string new_label = label_i + "+" + label_j;

    // Update distances
    for (size_t k = 0; k < num_sequences; k++) {
      if (!active[k] || k == best_i || k == best_j) continue;

      double new_distance = (sim_matrix(best_i, k) + sim_matrix(best_j, k)) / 2.0;
      sim_matrix(best_i, k) = new_distance;
      sim_matrix(k, best_i) = new_distance;
    }

    node_labels[best_i] = new_label;
    active[best_j] = false;

    active_count = std::count(active.begin(), active.end(), true);
  }
}