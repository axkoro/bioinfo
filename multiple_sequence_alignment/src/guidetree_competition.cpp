#include <iostream>

#include "align.hpp"

int main(int argc, char const* argv[]) {
  std::string fasta_path = argv[1];
  auto sequences = read_fasta(fasta_path);

  std::string matrix_path = argv[2];
  ScoringMatrix mat(matrix_path);

  auto result = compute_guidetree(sequences, mat);
  std::cout << result.guidetree_steps;

  return 0;
}
