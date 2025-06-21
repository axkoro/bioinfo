#include "align.hpp"

int main(int argc, char const* argv[]) {
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0] << " <fasta_file_path> <scoring_matrix_path>" << std::endl;
    return 1;
  }

  std::string fasta_path = argv[1];
  auto sequences = read_fasta(fasta_path);
  if (sequences.size() < 2) {
    std::cerr << "Error: The FASTA file needs to contain at least 2 sequences." << std::endl;
    return 1;
  }

  std::string matrix_path = argv[2];
  ScoringMatrix mat(matrix_path);

  int score = needleman_wunsch(sequences[0].data, sequences[1].data, mat);

  std::cout << score << std::endl;

  return 0;
}
