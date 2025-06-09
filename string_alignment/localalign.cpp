#include <algorithm>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <vector>

// Data Structures

struct AlignmentData {
  std::string aligned_a;  // The gapped alignment core for sequence A
  std::string aligned_b;  // The gapped alignment core for sequence B
  size_t start_pos_a;     // Start index of the core in the original sequence A
  size_t end_pos_a;       // End index (exclusive) of the core in sequence A
  size_t start_pos_b;     // Start index of the core in the original sequence B
  size_t end_pos_b;       // End index (exclusive) of the core in sequence B
};

struct SmithWatermanResult {
  int alignment_score;
  std::vector<AlignmentData> alignments;
};

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
};

// FASTA Reader

std::vector<std::string> read_fasta(const std::string& path) {
  std::ifstream infile(path);
  if (!infile.is_open()) {
    throw std::runtime_error("Could not open file at " + path);
  }

  std::vector<std::string> sequences;
  std::string current_sequence;
  std::string line;

  while (std::getline(infile, line)) {
    if (line.empty()) continue;

    if (line[0] == '>') {
      if (!current_sequence.empty()) {
        sequences.push_back(current_sequence);
      }
      current_sequence.clear();
    } else {
      // Remove potential carriage return from Windows-style line endings
      if (!line.empty() && line.back() == '\r') line.pop_back();
      current_sequence += line;
    }
  }

  sequences.push_back(current_sequence);  // Add the last sequence in the file

  return sequences;
}

// Smith-Waterman Alignment

enum class Direction { DIAGONAL, UP, LEFT };

SmithWatermanResult smith_waterman(const std::string& a, const std::string& b) {
  struct Costs {
    const int match = 1;
    const int insert = -1;
    const int del = -1;
    const int replace = -1;
  } costs;

  const size_t num_rows = b.size();
  const size_t num_cols = a.size();

  // The matrices are size N+1 to accommodate the initial row/column of zeros
  Matrix<int> similarity_matrix(num_rows + 1, num_cols + 1);
  Matrix<std::vector<Direction>> backtrace_pointers(num_rows + 1, num_cols + 1);

  int max_score = 0;
  std::vector<std::pair<size_t, size_t>> max_score_entries;

  // Compute similarity matrix
  for (size_t row = 1; row <= num_rows; ++row) {
    for (size_t col = 1; col <= num_cols; ++col) {
      const int diagonal_score = similarity_matrix(row - 1, col - 1) +
                                 (a[col - 1] == b[row - 1] ? costs.match : costs.replace);
      const int up_score = similarity_matrix(row - 1, col) + costs.del;
      const int left_score = similarity_matrix(row, col - 1) + costs.insert;

      int current_max = std::max({0, diagonal_score, up_score, left_score});
      similarity_matrix(row, col) = current_max;

      // Store pointers for traceback
      if (current_max > 0) {
        if (current_max == diagonal_score)
          backtrace_pointers(row, col).push_back(Direction::DIAGONAL);
        if (current_max == up_score) backtrace_pointers(row, col).push_back(Direction::UP);
        if (current_max == left_score) backtrace_pointers(row, col).push_back(Direction::LEFT);
      }

      if (current_max > max_score) {
        max_score = current_max;
        max_score_entries.clear();
        max_score_entries.emplace_back(row, col);
      } else if (current_max == max_score) {
        max_score_entries.emplace_back(row, col);
      }
    }
  }

  std::vector<AlignmentData> alignments;
  // Use a set to store unique alignments, preventing duplicates if different
  // max-score cells trace back to the same alignment.
  std::set<std::pair<std::string, std::string>> unique_alignments;

  if (max_score > 0) {
    for (const auto& entry : max_score_entries) {
      std::string path_a;
      std::string path_b;
      size_t row = entry.first;
      size_t col = entry.second;

      const size_t end_c = col;
      const size_t end_r = row;

      while (similarity_matrix(row, col) != 0) {
        Direction dir = backtrace_pointers(row, col)[0];

        switch (dir) {
          case Direction::DIAGONAL:
            path_a += a[col - 1];
            path_b += b[row - 1];
            row--;
            col--;
            break;
          case Direction::UP:
            path_a += '_';
            path_b += b[row - 1];
            row--;
            break;
          case Direction::LEFT:
            path_a += a[col - 1];
            path_b += '_';
            col--;
            break;
        }
      }

      std::reverse(path_a.begin(), path_a.end());
      std::reverse(path_b.begin(), path_b.end());

      bool duplicate = unique_alignments.find({path_a, path_b}) != unique_alignments.end();
      if (!duplicate) {
        alignments.emplace_back(AlignmentData{path_a, path_b, col, end_c, row, end_r});
        unique_alignments.insert({path_a, path_b});
      }
    }
  }

  return {max_score, alignments};
}

void print_alignment(const std::string& seq_a, const std::string& seq_b,
                     const AlignmentData& alignment) {
  size_t prefix_len_a = alignment.start_pos_a;
  size_t prefix_len_b = alignment.start_pos_b;
  size_t suffix_len_a = seq_a.length() - alignment.end_pos_a;
  size_t suffix_len_b = seq_b.length() - alignment.end_pos_b;

  size_t max_prefix = std::max(prefix_len_a, prefix_len_b);
  size_t max_suffix = std::max(suffix_len_a, suffix_len_b);

  size_t core_len = alignment.aligned_a.length();

  size_t final_len = max_prefix + core_len + max_suffix;

  // Create two strings of the final length, filled entirely with '*' padding
  std::string display_a(final_len, '*');
  std::string display_b(final_len, '*');

  // Place the alignment cores into the display strings at an offset
  display_a.replace(max_prefix, core_len, alignment.aligned_a);
  display_b.replace(max_prefix, core_len, alignment.aligned_b);

  std::cout << display_a << "\n";
  std::cout << display_b << "\n";
}

int main(int argc, char* argv[]) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <fasta_file_path>" << std::endl;
    return 1;
  }

  try {
    std::vector<std::string> sequences = read_fasta(argv[1]);

    if (sequences.size() < 2) {
      std::cerr << "Error: FASTA file must contain at least two sequences." << std::endl;
      return 1;
    }

    if (sequences.size() > 2) {
      std::cout << "Warning: File contains more than two sequences. "
                << "Aligning the first two only." << std::endl;
    }

    SmithWatermanResult result = smith_waterman(sequences[0], sequences[1]);

    std::cout << result.alignment_score << "\n\n";

    for (const auto& alignment : result.alignments) {
      print_alignment(sequences[0], sequences[1], alignment);
      std::cout << "\n";
    }

  } catch (const std::runtime_error& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}