#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <vector>

std::vector<std::string> read_fasta(const std::string& path) {
  std::vector<std::string> sequences;
  std::ifstream infile(path);

  if (!infile.is_open()) {
    std::string error_message = "Could not open file at " + path + "\n";
    throw std::runtime_error(error_message);
  }

  std::string line;
  std::string current_sequence;

  while (std::getline(infile, line)) {
    if (line.empty()) continue;

    if (line[0] == '>') {  // A new sequence begins with '>'.
      if (!current_sequence.empty()) {
        sequences.push_back(current_sequence);
      }

      current_sequence.clear();
    } else {                      // This line is part of the current sequence, so append it.
      if (line.back() == '\r') {  // Remove the carriage return if it exists (for CRLF files).
        line.pop_back();
      }

      current_sequence += line;
    }
  }

  sequences.push_back(current_sequence);

  infile.close();
  return sequences;
}

template <class T>
class Matrix {
 private:
  size_t num_rows;
  size_t num_cols;
  std::vector<T> data;

 public:
  Matrix(size_t rows, size_t cols) : num_rows(rows), num_cols(cols), data(rows * cols) {}

  // Accessor for modifying elements
  T& operator()(size_t row, size_t col) {
    if (row >= num_rows || col >= num_cols) {
      throw std::out_of_range("Matrix access out of bounds");
    }
    return data[row * num_cols + col];
  }

  // Accessor for reading elements
  const T& operator()(size_t row, size_t col) const {
    if (row >= num_rows || col >= num_cols) {
      throw std::out_of_range("Matrix access out of bounds");
    }
    return data[row * num_cols + col];
  }

  size_t rows() const { return num_rows; }
  size_t cols() const { return num_cols; }
};

typedef struct {
  int offset;
  int length;
} Alignment;

typedef struct {
  int alignment_score;
  std::vector<Alignment> alignments;
} AlignmentResult;

std::vector<size_t> all_max_elements(const std::vector<int>& elements) {
  int max = elements[0];
  std::vector<size_t> max_elements;
  for (size_t i = 0; i < elements.size(); i++) {
    if (elements[i] > max) {
      max_elements.clear();
      max_elements.push_back(i);
      max = elements[i];
    } else if (elements[i] == max) {
      max_elements.push_back(i);
    }
  }
  return max_elements;
}

AlignmentResult smith_waterman(const std::string& a, const std::string& b) {
  struct cost_table {
    int match = 1;
    int insert = -1;
    int del = -1;
    int replace = -1;
  } cost_table;

  size_t num_cols = a.size();
  size_t num_rows = b.size();
  Matrix<int> similarity_matrix(num_rows + 1, num_cols + 1);

  Matrix<std::vector<size_t>> backtrace_pointers(
      num_rows + 1, num_cols + 1);  // contains 0 for row, 1 for col, 2 for diag

  // Compute Similarity Matrix
  // col 0 and row 0 are initialized to 0 (as are all other entries)
  for (size_t row = 1; row < num_rows; row++) {
    for (size_t col = 1; col < num_cols; col++) {
      std::vector<int> costs(3);
      costs[0] = similarity_matrix(row - 1, col) + cost_table.insert;
      costs[1] = similarity_matrix(row, col - 1) + cost_table.del;
      costs[2] =
          similarity_matrix(row - 1, col - 1) + (a == b ? cost_table.match : cost_table.replace);

      std::vector<size_t> max_elements = all_max_elements(costs);
      int max = std::max(0, costs[max_elements[0]]);
      similarity_matrix(row, col) = max;

      if (max > 0) {
        backtrace_pointers(row, col) = max_elements;
      }
    }
  }

  // Find maximum entries
  int max = 0;
  std::vector<std::pair<size_t, size_t>> max_entries;
  for (size_t row = 1; row < num_rows; row++) {
    for (size_t col = 1; col < num_cols; col++) {
      int entry = similarity_matrix(row, col);
      if (entry > max) {
        max = entry;
        max_entries.clear();
        max_entries.emplace_back(row, col);
      } else if (entry == max) {
        max_entries.emplace_back(row, col);
      }
    }
  }

  // TODO: Backtrace
  std::vector<Alignment> alignments;
  for (auto&& entry : max_entries) {
  }

  return AlignmentResult{max, alignments};
}

void print_alignment(const std::vector<std::string>& sequences,
                     const Alignment& alignment) {  // TODO:
}

int main(int argc, char* argv[]) {
  if (argc != 2) throw std::invalid_argument("Please pass exactly one file path as argument");
  std::string path = argv[1];
  std::vector<std::string> sequences = read_fasta(path);

  int num_sequences = sequences.size();
  if (num_sequences < 2) {
    throw std::invalid_argument(
        "The provided file contains less than two sequences. To perform alignment, exactly two "
        "sequences are required.");
  } else if (sequences.size() > 2) {
    std::cout << "WARNING: The provided file contains more than two sequences. Performing "
                 "alignment on the first two and ignoring the rest ..."
              << std::endl;
  }

  AlignmentResult result = smith_waterman(sequences[0], sequences[1]);

  std::cout << result.alignment_score << "\n";
  for (auto alignment : result.alignments) {
    print_alignment(sequences, alignment);
  }
}
