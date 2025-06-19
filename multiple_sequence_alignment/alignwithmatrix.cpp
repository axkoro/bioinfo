#include <fstream>
#include <iostream>
#include <set>
#include <sstream>

std::vector<std::string> read_fasta(const std::string& path) {
  std::ifstream file(path);
  if (!file.is_open()) {
    throw std::runtime_error("Could not open file at " + path);
  }

  std::vector<std::string> sequences;
  std::string current_sequence;
  std::string line;

  while (std::getline(file, line)) {
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

class ScoringMatrix {
 public:
  ScoringMatrix(const std::string& path) : matrix(0, 0) {
    std::ifstream file(path);
    if (!file.is_open()) {
      throw std::runtime_error("Could not open file at " + path);
    }
    load_from_stream(file);
  }

  int score(char row_char, char col_char) const {
    return matrix(get_index_for_char(row_char), get_index_for_char(col_char));
  }

  const std::vector<char>& get_alphabet() const { return alphabet; }

  void print() const { matrix.print(); };

 private:
  Matrix<int> matrix;
  std::vector<char> alphabet;
  std::unordered_map<char, size_t> char_to_index;

  size_t get_index_for_char(char c) const {
    auto it = char_to_index.find(c);
    if (it == char_to_index.end()) {
      std::stringstream ss;
      ss << "Character '" << c << "' does not have an entry in the scoring matrix.";
      throw std::out_of_range(ss.str());
    }
    return it->second;
  }

  void load_from_stream(std::istream& stream) {
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
          throw std::runtime_error("Invalid or missing score on line " +
                                   std::to_string(line_number) + " for row '" +
                                   std::string(1, row_char) + "'.");
        }
        matrix(row_idx, col_idx) = score_value;
      }
    }
  }
};

int needleman_wunsch(const std::string& seq1, const std::string& seq2, const ScoringMatrix& mat) {
  const size_t num_cols = seq1.size();
  const size_t num_rows = seq2.size();

  Matrix<int> similarity_matrix(num_rows + 1, num_cols + 1);

  // Initialize first row and column
  for (size_t i = 0; i < num_rows; i++) {
    similarity_matrix(i, 0) = i;
  }
  for (size_t i = 0; i < num_cols; i++) {
    similarity_matrix(0, i) = i;
  }

  // Compute similarity matrix
  for (size_t row = 1; row <= num_rows; ++row) {
    for (size_t col = 1; col <= num_cols; ++col) {
      char a = seq1.at(col - 1);
      char b = seq2.at(row - 1);

      const int diagonal_score = similarity_matrix(row - 1, col - 1) + mat.score(a, b);
      const int up_score = similarity_matrix(row - 1, col) + mat.score(a, '*');
      const int left_score = similarity_matrix(row, col - 1) + mat.score(b, '*');

      int score = std::max({diagonal_score, up_score, left_score});
      similarity_matrix(row, col) = score;
    }
  }

  return similarity_matrix(num_rows, num_cols);
}

int main(int argc, char const* argv[]) {
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0] << " <fasta_file_path> <scoring_matrix_path>" << std::endl;
    return 1;
  }

  auto sequences = read_fasta(argv[1]);
  if (sequences.size() < 2) {
    std::cerr << "Error: The FASTA file needs to contain at least 2 sequences." << std::endl;
    return 1;
  }

  ScoringMatrix mat(argv[2]);

  int score = needleman_wunsch(sequences[0], sequences[1], mat);

  std::cout << score << std::endl;

  return 0;
}
