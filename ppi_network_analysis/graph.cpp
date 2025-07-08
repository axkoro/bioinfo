#include "graph.hpp"
#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>

bool Graph::is_number(const std::string &str) {
  try {
    std::stod(str);
    return true;
  } catch (...) {
    return false;
  }
}

Graph::Graph(const std::string &filename) {
  std::ifstream file(filename);

  if (!file.is_open()) {
    throw std::runtime_error("Could not open file: " + filename);
  }

  int next_id = 0;
  std::string line;
  bool first_line = true;

  while (std::getline(file, line)) {
    if (first_line) {
      first_line = false;
      // Check if it's a header line
      std::istringstream check_stream(line);
      std::string p1, p2, val;
      check_stream >> p1 >> p2 >> val;

      if (!is_number(val)) {
        continue; // Skip header
      }
    }

    std::istringstream iss(line);
    std::string protein1, protein2, value;

    if (!(iss >> protein1 >> protein2 >> value)) {
      continue; // Skip bad lines
    }

    // Get or create IDs for proteins
    int id1, id2;

    auto it1 = protein_to_id.find(protein1);
    if (it1 == protein_to_id.end()) {
      id1 = next_id++;
      protein_to_id[protein1] = id1;
      id_to_protein.push_back(protein1);
      adjacency_list.resize(next_id);
    } else {
      id1 = it1->second;
    }

    auto it2 = protein_to_id.find(protein2);
    if (it2 == protein_to_id.end()) {
      id2 = next_id++;
      protein_to_id[protein2] = id2;
      id_to_protein.push_back(protein2);
      adjacency_list.resize(next_id);
    } else {
      id2 = it2->second;
    }

    // Add undirected edge (add to both adjacency lists)
    adjacency_list[id1].push_back(id2);
    adjacency_list[id2].push_back(id1);
  }

  file.close();
}

void Graph::clean_adjacency_lists() {
  for (size_t i = 0; i < adjacency_list.size(); i++) {
    std::sort(adjacency_list[i].begin(), adjacency_list[i].end());
    adjacency_list[i].erase(
        std::unique(adjacency_list[i].begin(), adjacency_list[i].end()),
        adjacency_list[i].end());
  }
}