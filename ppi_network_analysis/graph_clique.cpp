#include "graph.hpp"
#include <algorithm>
#include <iostream>
#include <set>
#include <vector>

struct Clique {
  std::set<int> nodes;

  bool operator<(const Clique &other) const { return nodes < other.nodes; }
};

// Check if a set of nodes forms a clique
bool is_clique(const std::vector<std::vector<int>> &adjacency_list,
               const std::set<int> &nodes) {
  std::vector<int> node_vec(nodes.begin(), nodes.end());

  for (size_t i = 0; i < node_vec.size(); i++) {
    for (size_t j = i + 1; j < node_vec.size(); j++) {
      int node1 = node_vec[i];
      int node2 = node_vec[j];

      // Check if node2 is in adjacency list of node1
      bool found = std::binary_search(adjacency_list[node1].begin(),
                                      adjacency_list[node1].end(), node2);
      if (!found) {
        return false;
      }
    }
  }
  return true;
}

// Compute intersection of two sorted vectors
std::vector<int> intersection(const std::vector<int> &a,
                              const std::vector<int> &b) {
  std::vector<int> result;
  std::set_intersection(a.begin(), a.end(), b.begin(), b.end(),
                        std::back_inserter(result));
  return result;
}

std::set<int> set_union(const std::set<int> &a, const std::set<int> &b) {
  std::set<int> result = a;
  result.insert(b.begin(), b.end());
  return result;
}

// Build initial set of 2-cliques (edges)
std::set<Clique> build_initial_cliques(const Graph &graph) {
  std::set<Clique> cliques;

  for (size_t i = 0; i < graph.adjacency_list.size(); i++) {
    for (int neighbor : graph.adjacency_list[i]) {
      if (static_cast<int>(i) < neighbor) { // Avoid duplicates
        Clique c;
        c.nodes.insert(i);
        c.nodes.insert(neighbor);
        cliques.insert(c);
      }
    }
  }

  return cliques;
}

// Find all maximal cliques (algo from the slides)
std::set<Clique> find_maximal_cliques(const Graph &graph) {
  // Build set S2 of all cliques of size 2 (edges)
  std::set<Clique> current_cliques = build_initial_cliques(graph);
  // std::cout << "Initial edges (2-cliques): " << current_cliques.size()
  //           << std::endl;

  int clique_size = 2;
  std::set<Clique> maximal_cliques;

  // Iteratively build larger cliques
  while (!current_cliques.empty()) {
    // std::cout << "Processing " << clique_size
    //           << "-cliques: " << current_cliques.size() << " cliques"
    //           << std::endl;

    std::set<Clique> next_cliques;
    std::vector<Clique> clique_vec(current_cliques.begin(),
                                   current_cliques.end());

    for (size_t j = 0; j < clique_vec.size(); j++) {
      for (size_t k = j + 1; k < clique_vec.size(); k++) {
        // Compute intersection of two cliques
        std::vector<int> vec1(clique_vec[j].nodes.begin(),
                              clique_vec[j].nodes.end());
        std::vector<int> vec2(clique_vec[k].nodes.begin(),
                              clique_vec[k].nodes.end());

        std::vector<int> intersect = intersection(vec1, vec2);

        // Check if they share exactly (clique_size - 2) nodes
        if (intersect.size() == static_cast<size_t>(clique_size - 2)) {
          // Compute union
          std::set<int> union_set =
              set_union(clique_vec[j].nodes, clique_vec[k].nodes);

          // Check if it forms a clique
          if (is_clique(graph.adjacency_list, union_set)) {
            Clique new_clique;
            new_clique.nodes = union_set;
            next_cliques.insert(new_clique);
          }
        }
      }
    }

    if (next_cliques.empty()) {
      // Current cliques are maximal
      maximal_cliques = current_cliques;
    }

    current_cliques = next_cliques;
    clique_size++;
  }

  return maximal_cliques;
}

void print_maximal_cliques(const std::set<Clique> &maximal_cliques,
                           const Graph &graph, int maximal_size) {
  std::cout << "\nMaximal clique size: " << maximal_size << std::endl;
  // std::cout << "Number of maximal cliques: " << maximal_cliques.size()
  //           << std::endl;

  int count = 0;
  for (const auto &clique : maximal_cliques) {
    // std::cout << "\nMaximal clique " << ++count << ":" << std::endl;
    // std::cout << "Proteins: ";
    bool first = true;
    for (int node_id : clique.nodes) {
      if (!first)
        std::cout << ", ";
      std::cout << graph.get_protein_name(node_id);
      first = false;
    }
    std::cout << std::endl;
  }
}

int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <edge_file_path>" << std::endl;
    return 1;
  }

  try {
    std::string filename = argv[1];

    // std::cout << "Reading graph from " << filename << "..." << std::endl;
    Graph graph(filename);

    graph.clean_adjacency_lists();

    // std::cout << "Graph loaded: " << graph.size() << " nodes" << std::endl;

    std::set<Clique> maximal_cliques = find_maximal_cliques(graph);

    int maximal_size = 0;
    if (!maximal_cliques.empty()) {
      maximal_size = maximal_cliques.begin()->nodes.size();
    }

    print_maximal_cliques(maximal_cliques, graph, maximal_size);

  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}