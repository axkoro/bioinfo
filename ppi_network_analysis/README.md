This directory contains the implementation for analyzing protein-protein interaction (PPI) networks.

## Files

- **Graph Class:**
  - `graph.hpp` - Header file defining the Graph class
  - `graph.cpp` - Implementation of the Graph class with file I/O and graph operations
- **Main Programs:**
  - `graph_degree.cpp` - Computes degree distribution of a PPI network
  - `graph_clique.cpp` - Finds maximal cliques in a PPI network
- **Utilities:**
  - `plot_degree_distribution.py` - Plots and analyzes the degree distribution
  - `Makefile` - Build script for compiling the C++ programs


## Usage

Compile both C++ programs (and the Graph library):
```bash
make all
```

### Degree Distribution

```bash
./graph_degree <path_to_edge_file>
```

Example:
```bash
./graph_degree data/9606.protein.physical.links.v11.5.txt
```

This will:
- Read the PPI network file
- Compute the [degree distribution](https://en.wikipedia.org/wiki/Degree_distribution)
- Output the distribution to console (format: `degree: count`)
- Save the distribution to `degree_distribution.txt` for plotting

#### Plot and Analyze Distribution

After running `graph_degree`, use the Python script to analyze:

```bash
python3 plot_degree_distribution.py
```

This will:
- Load the degree distribution from `degree_distribution.txt`
- Create multiple plots in `degree_distribution.pdf`

### Find Maximal Cliques

```bash
./graph_clique <path_to_edge_file>
```

Example:
```bash
# Extract first 1000 lines for testing
head -n 1001 9606.protein.physical.links.v11.5.txt > test_1000.txt
./graph_clique test_1000.txt
```

This will:
- Find all maximal cliques in the graph
- Output the size of maximal clique(s)
- List all proteins in each maximal clique

## Important Notes

1. **Memory Usage**: The full PPI network file is large (~2 million edges). The clique-finding algorithm may take significant time and memory for the full dataset. Start with smaller subsets for testing.

2. **File Format**: The input file should have the format:
   ```
   protein1_id protein2_id confidence_score
   ```
   The third column (confidence score) is ignored.

3. **Header Line**: The programs automatically detect and skip header lines.

4. **Graph Representation**: 
   - Edges are treated as undirected
   - Adjacency lists are kept sorted for efficiency
   - Internal integer IDs are used for performance

## Dependencies

- g++
- Python 3 with:
  - matplotlib
  - numpy
  - scipy
  - pandas

Install Python dependencies:
```bash
pip install matplotlib numpy scipy pandas
```