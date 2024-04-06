#pragma once
#include <fstream>
#include <iostream>
#include <unordered_map>

#include "Base/Base.h"
#include "Base/BasicAlgorithms.h"
#include "Graph.h"

namespace GraphLib {
    static std::unordered_map<std::string, int> global_vertex_labels, global_edge_labels;
    struct Branch {
        int vertex_label;
        std::vector<int> edge_labels;

        int GetVertexLabel() const { return vertex_label; }
        const std::vector<int> &GetEdgeLabels() const { return edge_labels; }

        bool operator<(const Branch &other) const {
            if (vertex_label != other.vertex_label)
                return vertex_label < other.vertex_label;
            return edge_labels < other.edge_labels;
        }
    };

    int BranchEditDistance(const Branch &a, const Branch &b) {
        int ed = 0;
        if (a.vertex_label != b.vertex_label) ed += 2;
        ed += std::max(a.edge_labels.size(), b.edge_labels.size());
        ed -= VectorIntersectionSize(a.edge_labels, b.edge_labels);
        return ed;
    }

    int BranchEditDistanceFromNull(const Branch &a) {
        return 2 + a.edge_labels.size();
    }

    class LabeledGraph : public Graph {
    public:
        std::vector<int> vertex_label, edge_label;
        std::unordered_map<int, int> vertex_label_frequency, edge_label_frequency;
        std::vector<std::vector<int>> adj_matrix;
        std::vector<int> degree_sequence;
        std::vector<int> branch_ids;
        std::vector<Branch> branches;

        std::vector<std::vector<int>> apsp;
        std::vector<int> eccentricity;

        LabeledGraph() {}
        ~LabeledGraph() {}

        LabeledGraph &operator=(const LabeledGraph &) = delete;

        int GetEccentricity(int v) { return eccentricity[v]; }
        int GetVertexLabel(int v) const { return vertex_label[v]; }
        int GetEdgeLabel(int edge_id) const {
            return edge_label[edge_id];
        }
        int GetEdgeLabel(int u, int v) const {
            return adj_matrix[u][v];
        }

        std::unordered_map<int, int> &GetVertexLabelFrequency() {
            return vertex_label_frequency;
        }
        std::unordered_map<int, int> &GetEdgeLabelFrequency() {
            return edge_label_frequency;
        }

        int GetVertexLabelFrequency(int l) { return vertex_label_frequency[l]; }
        int GetEdgeLabelFrequency(int l) { return edge_label_frequency[l]; }

        // Read/Write graph from/to file.
        void LoadFromFile(const std::string &filename);
        void WriteToFile(const std::string &filename);

        void PrintGraph();

        void ComputeAPSP();

        bool LoadFromStream(std::istream &in, std::string &buffer);
    };

    void LabeledGraph::LoadFromFile(const std::string &filename) {
        num_vertex = num_edge = 0;
        std::ifstream fin(filename);
        if (!fin.is_open()) {
            std::cerr << "Error opening file: " << filename << std::endl;
            exit(1);
        }
        std::string ignore, type, line;
        std::unordered_map<int, int> tmp_vertex_labels;
        while (getline(fin, line)) {
            auto tok = parse(line, " ");
            type = tok[0];
            tok.pop_front();
            if (type[0] == 'v') {
                int id = std::stoi(tok.front());
                tok.pop_front();
                int l;
                if (tok.empty())
                    l = 0;
                else {
                    auto it = global_vertex_labels.find(tok.front());
                    if (it == global_vertex_labels.end()) {
                        global_vertex_labels[tok.front()] = global_vertex_labels.size();
                    }
                    l = global_vertex_labels[tok.front()];
                    tok.pop_front();
                }
                num_vertex++;
                tmp_vertex_labels[id] = l;
                vertex_label_frequency[l]++;
            } else if (type[0] == 'e') {
                if (adj_matrix.size() != num_vertex) {
                    vertex_label.resize(num_vertex);
                    for (int i = 0; i < num_vertex; i++)
                        vertex_label[i] = tmp_vertex_labels[i];
                    adj_matrix.resize(num_vertex, std::vector<int>(num_vertex, -1));
                    adj_list.resize(num_vertex);
                    vertex_label.resize(num_vertex);
                }
                int v1, v2;
                v1 = std::stoi(tok.front());
                tok.pop_front();
                v2 = std::stoi(tok.front());
                tok.pop_front();
                adj_list[v1].push_back(v2);
                adj_list[v2].push_back(v1);
                edge_list.emplace_back(v1, v2);
                int el = 0;
                if (!tok.empty()) {
                    auto it = global_edge_labels.find(tok.front());
                    if (it == global_edge_labels.end()) {
                        global_edge_labels[tok.front()] = global_edge_labels.size() + 1;
                    }
                    el = global_edge_labels[tok.front()];
                    tok.pop_front();
                }
                edge_label.emplace_back(el);
                edge_label_frequency[el]++;
                adj_matrix[v1][v2] = el;
                adj_matrix[v2][v1] = el;
                max_degree = std::max(max_degree, (int)std::max(adj_list[v1].size(), adj_list[v2].size()));
            } else if (type[0] == 't') {
                this->id = std::stoi(tok.back());
            }
        }
        num_edge = edge_list.size();
    }

    // Read one graph from the stream. The "buffer" parameter stores a header
    // line (starting with 't') that has been read ahead for the next graph.
    // Returns true if a graph was successfully loaded.
    bool LabeledGraph::LoadFromStream(std::istream &in, std::string &buffer) {
        // Reset the graph data
        num_vertex = num_edge = max_degree = 0;
        vertex_label.clear();
        adj_matrix.clear();
        adj_list.clear();
        edge_list.clear();
        edge_label.clear();
        vertex_label_frequency.clear();
        edge_label_frequency.clear();

        // Temporary container for vertex labels (mapping local id to label)
        std::unordered_map<int, int> tmp_vertex_labels;

        std::string line;
        // If a header was already buffered from a previous call, use it.
        if (!buffer.empty()) {
            line = buffer;
            buffer.clear();
        } else {
            if (!std::getline(in, line))
                return false;  // No more data.
        }

        // The first expected line for a graph should be a header (starts with 't')
        auto tokens = parse(line, " ");
        if (tokens.empty() || tokens.front()[0] != 't') {
            std::cerr << "Expected a graph header (line starting with 't'), got: " << line << std::endl;
            exit(1);
        }
        // Set the graph id (assumed to be the last token in the header line)
        id = std::stoi(tokens.back());

        // Read lines until the next header is encountered (or end-of-file)
        while (std::getline(in, line)) {
            if (line.empty())
                continue;
            // If a new graph header is encountered, buffer it and stop reading further.
            if (line[0] == 't') {
                buffer = line;
                break;
            }
            auto tok = parse(line, " ");
            if (tok.empty())
                continue;
            std::string type = tok.front();
            tok.pop_front();
            if (type[0] == 'v') {
                // Vertex line: expected format: "v <id> [label]"
                if (tok.empty()) {
                    std::cerr << "Invalid vertex line: " << line << std::endl;
                    continue;
                }
                int vid = std::stoi(tok.front());
                tok.pop_front();
                int lab = 0;
                if (!tok.empty()) {
                    std::string lab_str = tok.front();
                    tok.pop_front();
                    // Create new label if necessary
                    if (global_vertex_labels.find(lab_str) == global_vertex_labels.end())
                        global_vertex_labels[lab_str] = global_vertex_labels.size();
                    lab = global_vertex_labels[lab_str];
                }
                tmp_vertex_labels[vid] = lab;
                vertex_label_frequency[lab]++;
                num_vertex++;
            } else if (type[0] == 'e') {
                // When processing the first edge, ensure that our data structures
                // are initialized based on the number of vertices.
                if (adj_matrix.size() != static_cast<size_t>(num_vertex)) {
                    vertex_label.resize(num_vertex);
                    for (int i = 0; i < num_vertex; i++)
                        vertex_label[i] = tmp_vertex_labels[i];
                    adj_matrix.resize(num_vertex, std::vector<int>(num_vertex, -1));
                    adj_list.resize(num_vertex);
                }
                // Edge line: expected format: "e <v1> <v2> [label]"
                if (tok.size() < 2) {
                    std::cerr << "Invalid edge line: " << line << std::endl;
                    continue;
                }
                int v1 = std::stoi(tok.front());
                tok.pop_front();
                int v2 = std::stoi(tok.front());
                tok.pop_front();
                // Update adjacency list and matrix.
                if (v1 < static_cast<int>(adj_list.size()) && v2 < static_cast<int>(adj_list.size())) {
                    adj_list[v1].push_back(v2);
                    adj_list[v2].push_back(v1);
                }
                edge_list.emplace_back(v1, v2);
                int elab = 0;
                if (!tok.empty()) {
                    std::string lab_str = tok.front();
                    tok.pop_front();
                    if (global_edge_labels.find(lab_str) == global_edge_labels.end())
                        global_edge_labels[lab_str] = global_edge_labels.size() + 1;
                    elab = global_edge_labels[lab_str];
                }
                edge_label.push_back(elab);
                edge_label_frequency[elab]++;
                if (v1 < static_cast<int>(adj_matrix.size()) && v2 < static_cast<int>(adj_matrix.size())) {
                    adj_matrix[v1][v2] = elab;
                    adj_matrix[v2][v1] = elab;
                }
                // Update maximum degree for the graph.
                max_degree = std::max(max_degree, static_cast<int>(std::max(adj_list[v1].size(), adj_list[v2].size())));
            }
            // You can add handling for other line types here if needed.
        }
        num_edge = edge_list.size();
        return true;
    }

    void LabeledGraph::WriteToFile(const std::string &filename) {
        fs::path filepath = filename;
        create_directories(filepath.parent_path());
        std::ofstream out(filename);
        out << "t " << GetNumVertices() << ' ' << GetNumEdges() / 2 << '\n';
        for (int i = 0; i < GetNumVertices(); i++) {
            out << "v " << i << ' ' << GetVertexLabel(i) << ' ' << GetDegree(i) << '\n';
        }
        int idx = 0;
        for (auto &e : edge_list) {
            if (e.first < e.second) {
                out << "e " << e.first << ' ' << e.second << ' ' << GetEdgeLabel(idx)
                    << '\n';
            }
            idx++;
        }
    }

    void LabeledGraph::PrintGraph() {
        fprintf(stderr, "Graph %d with %d vertices and %d edges\n", GetId(), GetNumVertices(), GetNumEdges());
        for (int i = 0; i < GetNumVertices(); i++) {
            fprintf(stderr, "Vertex %d: Label %d, Degree %d\n  Nbrs: ", i, GetVertexLabel(i), GetDegree(i));
            for (int j : adj_list[i]) {
                fprintf(stderr, "%d(el=%d) ", j, GetEdgeLabel(i, j));
            }
            fprintf(stderr, "\n");
        }
    }

    void LabeledGraph::ComputeAPSP() {
        printf("APSP Called! num_vertices = %d\n", GetNumVertices());
        apsp.resize(GetNumVertices(), std::vector<int>(GetNumVertices(), INF));
        for (int i = 0; i < GetNumVertices(); i++) {
            apsp[i][i] = 0;
        }
        for (auto &e : edge_list) {
            apsp[e.first][e.second] = apsp[e.second][e.first] = 1;
        }
        for (int k = 0; k < GetNumVertices(); k++) {
            for (int i = 0; i < GetNumVertices(); i++) {
                for (int j = 0; j < GetNumVertices(); j++) {
                    apsp[i][j] = std::min(apsp[i][j], apsp[i][k] + apsp[k][j]);
                }
            }
        }
        eccentricity.resize(GetNumVertices(), 0);
        for (int i = 0; i < GetNumVertices(); i++) {
            for (int j = 0; j < GetNumVertices(); j++) {
                if (apsp[i][j] != INF)
                    eccentricity[i] = std::max(eccentricity[i], apsp[i][j]);
            }
        }
    }
}  // namespace GraphLib
// note: research log marker 12
