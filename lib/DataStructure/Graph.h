#pragma once
#include <vector>
/*
 * Basic graph class. All graph class should inherit from this.
 * Assume unlabeled graph with only adjacency list and edge list.
 */
namespace GraphLib {
class Graph {
 protected:
  std::vector<std::vector<int>> adj_list;
  std::vector<std::pair<int, int>> edge_list;
  int num_vertex, num_edge, max_degree, id;

 public:
  ~Graph() = default;
  Graph() {};
  std::vector<int> &GetNeighbors(int v) { return adj_list[v]; }

  std::vector<std::pair<int, int>> &GetEdges() { return edge_list; }

  inline std::pair<int, int> GetEdge(int i) { return edge_list[i]; }

  inline int GetNumVertices() const { return num_vertex; }

  inline int GetNumEdges() const { return num_edge; }

  inline int GetMaxDegree() const { return max_degree; }

  inline int GetDegree(int v) const { return adj_list[v].size(); }

  inline int GetId() const { return id; }
};
}  // namespace GraphLib