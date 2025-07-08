#pragma once

#include <string>
#include <unordered_map>
#include <vector>

class Graph {
private:
  std::unordered_map<std::string, int> protein_to_id;
  std::vector<std::string> id_to_protein;

  // Check if a string is a valid number
  static bool is_number(const std::string &str);

public:
  std::vector<std::vector<int>> adjacency_list;

  explicit Graph(const std::string &filename);

  // Remove duplicates and sort for binary search
  void clean_adjacency_lists();

  // Get the number of nodes in the graph
  size_t size() const { return adjacency_list.size(); }

  // Get the degree of a node
  int degree(int node_id) const { return adjacency_list[node_id].size(); }

  // Get protein name from node ID
  const std::string &get_protein_name(int node_id) const {
    return id_to_protein[node_id];
  }

  // Get node ID from protein name (returns -1 if not found)
  int get_node_id(const std::string &protein_name) const {
    auto it = protein_to_id.find(protein_name);
    return (it != protein_to_id.end()) ? it->second : -1;
  }
};