#pragma once
#include "Base/Timer.h"
#include "GED/DynamicHungarianGED.h"

namespace GraphLib::GraphSimilarity {
    class DFSDH : public DynamicHungarianGED {
    public:
        std::vector<int> matching_order, inv_matching_order;

#ifdef AUTO
        std::vector<int> u_list = {};
        std::vector<int> v_list = {};
        virtual std::pair<int, int> GetNextPairFirst(std::vector<int>& u_list, std::vector<int>& v_list) { return {-1, -1}; };

        std::vector<int> u_vertex_to_group = {};
        std::vector<int> v_vertex_to_group = {};

        std::vector<std::vector<int>> u_groups = {};
        std::vector<std::vector<int>> v_groups = {};
        std::vector<int> u_group_sizes = {};
        std::vector<int> v_group_sizes = {};

        // std::vector<std::vector<int>> v_partitions = {};
        std::vector<std::vector<int>> u_auto_stacks = {};

        // bool** uv_auto_visited = nullptr;
        // bool* v_auto_visited = nullptr;

        void read_auto(std::string graph_id, bool is_g1) {
            std::string graph_path = "../data/auto/pubchem_large100/" + graph_id + ".txt";
            std::ifstream infile(graph_path);
            std::string line;
            if (std::getline(infile, line)) {
                std::istringstream iss(line);
                int num;
                while (iss >> num) {
                    if (is_g1)
                        u_list.push_back(num);
                    else
                        v_list.push_back(num);
                }
            }
            infile.close();

            if (is_g1) {
                u_vertex_to_group = std::vector<int>(G1->GetNumVertices(), -1);
            } else {
                v_vertex_to_group = std::vector<int>(G2->GetNumVertices(), -1);
            }

            graph_path = "../data/auto/pubchem_large100-groups/" + graph_id + ".txt";
            infile.open(graph_path);
            if (!infile.is_open()) {
                std::cerr << "Error opening file: " << graph_path << std::endl;
                return;
            }
            while (std::getline(infile, line)) {
                std::istringstream iss(line);
                int num;
                std::vector<int> group;
                while (iss >> num) {
                    group.push_back(num);
                }
                if (is_g1) {
                    u_groups.push_back(group);
                    u_group_sizes.push_back(group.size());
                    for (int v : group) {
                        u_vertex_to_group[v] = u_groups.size() - 1;
                    }
                    std::vector<int> auto_stack;
                    u_auto_stacks.push_back(auto_stack);
                } else {
                    v_groups.push_back(group);
                    v_group_sizes.push_back(group.size());
                    for (int v : group) {
                        v_vertex_to_group[v] = v_groups.size() - 1;
                    }
                    // v_auto_stacks.push_back(partition);
                }
            }
            infile.close();

            // if (is_g1) {
            //     // uv_auto_visited = new bool*[u_groups.size()];
            // } else {
            //     // for (int i = 0; i < v_groups.size(); i++) {
            //     //     uv_auto_visited[i] = new bool[v_groups.size()];
            //     //     std::fill(uv_auto_visited[i], uv_auto_visited[i] + v_groups.size(), false);
            //     // }
            //     v_auto_visited = new bool[v_groups.size()];
            //     std::fill(v_auto_visited, v_auto_visited + v_groups.size(), false);
            // }
        }
#endif

        void DFS(int u, int v) {
            cnt++;
            depth++;
            int chcost = 0;
            if (depth == G1->GetNumVertices()) {
                depth--;
                return;
            }
            hungarian_states[depth + 1] = hungarian_states[depth];
            if (depth != 0) {
                mapping[u] = v;
                inverse_mapping[v] = u;
                Match(u, v);
                chcost = ChildEditCost(u, v);
                cost += chcost;
                ComputeBranchDistanceMatrixDynamic(u, v, 1);
                UpdateParikhVector(u, v);
            }
            int lb, ub;
            int uprime = matching_order[depth];
            for (int i = 0; i < G2->GetNumVertices() - depth; i++) {
                std::tie(lb, ub) = LowerBound();
                int vprime = hungarian_states[depth + 1].assignment[uprime];
                if (lb > tau) {
                    break;
                }
                if (ub < current_best) {
                    current_best = ub;
                    if (current_best <= tau) {
                        flag = true;
                        return;
                    }
                }

                // if (!u_auto_stacks[u_vertex_to_group[uprime]].empty()) {
                //     if (u_auto_stacks[u_vertex_to_group[uprime]].back() < vprime) {
                //         continue;
                //     } else {
                //         u_auto_stacks[u_vertex_to_group[uprime]].push_back(vprime);
                //     }
                // } else {
                //     u_auto_stacks[u_vertex_to_group[uprime]].push_back(vprime);
                // }

#ifdef AUTO
                // if (depth == 0) {
                //     // for (int j = 0; j < v_groups.size(); j++) {
                //     //     std::cout << "v_groups[" << j << "]: " << v_auto_visited[j] << std::endl;
                //     // }
                //     // std::cout << "\n";
                //     if (v_auto_visited[v_vertex_to_group[vprime]]) {
                //         matrix[uprime][vprime] += INF2;
                //         hungarian_states[depth + 1].assignment[uprime] = -1;
                //         hungarian_states[depth + 1].inverse_assignment[vprime] = -1;
                //         continue;
                //     }
                // }
#endif

                // const int prev_v_group_size = v_groups.size();
                // for (int nv : G2->GetNeighbors(vprime)) {
                //     v_partitions[v_vertex_to_group[nv]].push_back(nv);
                // }
                // for (int j = 0; j < v_partitions.size(); j++) {
                //     const int partition_size = v_partitions[j].size();
                //     if (partition_size == 0) continue;
                //     if (partition_size < v_group_sizes[j]) {
                //         std::vector<int> partition = v_partitions[j];
                //         for (int k = 0; k < partition_size; k++) {
                //             v_groups[j].erase(std::find(v_groups[j].begin(), v_groups[j].end(), v_partitions[j][k]));
                //             v_vertex_to_group[v_partitions[j][k]] = v_groups.size();
                //         }
                //         v_group_sizes[j] -= partition_size;
                //         v_groups.push_back(partition);
                //         v_group_sizes.push_back(partition_size);
                //     } else {
                //         v_partitions[j].clear();
                //     }
                // }

                DFS(uprime, vprime);

                // u_auto_stacks[u_vertex_to_group[uprime]].pop_back();

                // for (int j = 0; j < v_partitions.size(); j++) {
                //     const int partition_size = v_partitions[j].size();
                //     if (partition_size == 0) continue;
                //     for (int k = 0; k < partition_size; k++) {
                //         v_groups[j].push_back(v_partitions[j][k]);
                //         v_vertex_to_group[v_partitions[j][k]] = j;
                //     }
                //     v_group_sizes[j] += partition_size;
                //     v_partitions[j].clear();
                // }
                // while (v_groups.size() > prev_v_group_size) {
                //     v_groups.pop_back();
                //     v_group_sizes.pop_back();
                // }

#ifdef AUTO
                // if (depth == 0) {
                //     v_auto_visited[v_vertex_to_group[vprime]] = true;
                // }
#endif

                if (flag) {
                    // depth--;
                    // hungarian_states[depth + 1].assignment[uprime] = -1;
                    // hungarian_states[depth + 1].inverse_assignment[vprime] = -1;
                    // for (int i = 0; i < G2->GetNumVertices(); i++) {
                    //     if (matrix[uprime][i] >= INF2) {
                    //         matrix[uprime][i] -= INF2;
                    //     }
                    // }
                    // if (depth != 0) {
                    //     RestoreParikhVector(u, v);
                    //     ComputeBranchDistanceMatrixDynamic(u, v, 0);
                    //     mapping[u] = -1;
                    //     inverse_mapping[v] = -1;
                    //     RemoveMatch(u, v);
                    //     cost -= chcost;
                    // }
                    return;
                }
#ifdef AUTO
                if (depth == 0) {
                    for (int v_auto : v_groups[v_vertex_to_group[vprime]]) {
                        matrix[uprime][v_auto] += INF2;
                        i++;
                    }
                    i--;
                } else
#endif
                {
                    matrix[uprime][vprime] += INF2;
                }
                // for (int u_auto : u_groups[u_vertex_to_group[uprime]]) {
                //     for (int v_auto : v_groups[v_vertex_to_group[vprime]]) {
                //         matrix[u_auto][v_auto] += INF2;
                //     }
                // }

                hungarian_states[depth + 1].assignment[uprime] = -1;
                hungarian_states[depth + 1].inverse_assignment[vprime] = -1;
            }

            for (int i = 0; i < G2->GetNumVertices(); i++) {
                if (matrix[uprime][i] >= INF2) {
                    matrix[uprime][i] -= INF2;
                }
            }

            // for (int i = 0; i < G1->GetNumVertices(); i++) {
            //     for (int j = 0; j < G2->GetNumVertices(); j++) {
            //         if (matrix[i][j] >= INF2) {
            //             matrix[i][j] -= INF2;
            //         }
            //     }
            // }

            if (depth != 0) {
                RestoreParikhVector(u, v);
                ComputeBranchDistanceMatrixDynamic(u, v, 0);
                mapping[u] = -1;
                inverse_mapping[v] = -1;
                RemoveMatch(u, v);
                cost -= chcost;
            }
            depth--;
        }

        void
        ComputeMatchingOrder() {
            OldMatchingOrder();
        }

        void NewMatchingOrder();
        void OldMatchingOrder();

        int ComputeGED() {
            Timer timer;
            timer.Start();

            matching_order.clear();
            ComputeMatchingOrder();
            num_nodes = 0;
            current_best = 1e9;

            Initialize();

#ifdef AUTO
            read_auto(std::to_string(G1->GetId()), true);
            read_auto(std::to_string(G2->GetId()), false);

            // print u_list, v_list, u_groups, v_groups, u_group_sizes, v_group_sizes, u_vertex_to_group, v_vertex_to_group;
            // for (int i = 0; i < u_list.size(); i++) {
            //     std::cout << "u_list[" << i << "] = " << u_list[i] << "\n";
            // }
            // for (int i = 0; i < v_list.size(); i++) {
            //     std::cout << "v_list[" << i << "] = " << v_list[i] << "\n";
            // }
            // for (int i = 0; i < u_groups.size(); i++) {
            //     std::cout << "u_groups[" << i << "] = ";
            //     for (int v : u_groups[i]) {
            //         std::cout << v << " ";
            //     }
            //     std::cout << "\n";
            // }
            // for (int i = 0; i < v_groups.size(); i++) {
            //     std::cout << "v_groups[" << i << "] = ";
            //     for (int v : v_groups[i]) {
            //         std::cout << v << " ";
            //     }
            //     std::cout << "\n";
            // }
            // for (int i = 0; i < u_vertex_to_group.size(); i++) {
            //     std::cout << "u_vertex_to_group[" << i << "] = " << u_vertex_to_group[i] << "\n";
            // }
            // for (int i = 0; i < v_vertex_to_group.size(); i++) {
            //     std::cout << "v_vertex_to_group[" << i << "] = " << v_vertex_to_group[i] << "\n";
            // }
#endif

            ComputeBranchDistanceMatrixInitial(matrix);
            ComputeParikhVector();
            tau = 0;
            depth = -1;
            auto [lb, ub] = LowerBound();
            current_best = ub;
            cnt++;
            tau = lb;
            while (true) {
                if (current_best <= tau) break;
                cnt--;
                depth = -1;
#ifdef AUTO
                // std::fill(v_auto_visited, v_auto_visited + v_groups.size(), false);
#endif
                DFS(-1, -1);
                tau++;
            }

            // tau = 501;
            // --cnt;
            // depth = -1;
            // DFS(-1, -1);
            // fprintf(stderr, "tau: %d, current_best: %d\n", tau, current_best);

            // flag = false;

            // tau = 26;
            // current_best = 29;
            // --cnt;
            // depth = -1;
            // DFS(-1, -1);
            // fprintf(stderr, "tau: %d, current_best: %d\n", tau, current_best);
            // exit(-1);

            // int l = lb;
            // int r = ub;
            // while (l <= r) {
            //     flag = false;
            //     tau = (l + r) / 2;
            //     fprintf(stderr, "[%d, %d] tau: %d\n", l, r, tau);
            //     cnt--;
            //     depth = -1;
            //     DFS(-1, -1);

            //     const bool found = (current_best <= tau);
            //     fprintf(stderr, "(%d) tau: %d, current_best: %d\n", found, tau, current_best);
            //     if (found) {
            //         r = current_best - 1;
            //     } else {
            //         l = tau + 1;
            //     }
            // }
            timer.Stop();

            log.AddResult("GED", current_best, RESULT_INT);
            log.AddResult("ElapsedTime", timer.GetTime(), RESULT_DOUBLE_FIXED);
            log.AddResult("SearchTreeNodes", cnt, RESULT_INT64);
            return current_best;
        }
    };  // namespace GraphLib::GraphSimilarity

    void DFSDH::NewMatchingOrder() {
        int N = G2->GetNumVertices();
        auto vlabel_freq = G2->GetVertexLabelFrequency();
        auto elabel_freq = G2->GetEdgeLabelFrequency();
        int root = 0;
        int cnt = 0;
        double root_weight = 0;
        for (int i = 0; i < G1->GetNumVertices(); i++) {
            double weight = 1 - vlabel_freq[G1->GetVertexLabel(i)] / double(G2->GetNumVertices());
            for (auto j : G1->GetNeighbors(i)) {
                weight += 1 - elabel_freq[G1->GetEdgeLabel(i, j)] / double(G2->GetNumEdges() * 2);
            }
            if (weight > root_weight) {
                root = i;
                root_weight = weight;
            }
        }

        std::vector<double> w(G1->GetNumVertices(), 0);
        std::vector<int> T(G1->GetNumVertices(), 0);
        w[root] = root_weight;

        std::priority_queue<std::pair<double, int>> weighted_queue;
        weighted_queue.push({w[root], root});

        while (!weighted_queue.empty()) {
            int u = weighted_queue.top().second;
            weighted_queue.pop();
            if (T[u] == 1) {
                continue;
            }

            matching_order.push_back(u);
            T[u] = 1;
            for (auto v : G1->GetNeighbors(u)) {
                // std::cout << T[v] << "\n";
                if (T[v] == 1) {
                    continue;
                }
                if (w[v] < 10e-6) {
                    w[v] += 1 - vlabel_freq[G1->GetVertexLabel(v)] / double(G2->GetNumVertices());
                }
                w[v] += 1 - elabel_freq[G1->GetEdgeLabel(u, v)] / double(G2->GetNumEdges() * 2);
                weighted_queue.push({w[v], v});
                // std::cout << "push  : " << v << " " << w[v] << "\n";
            }
        }
        for (int i = 0; i < G1->GetNumVertices(); i++) {
            if (T[i] == 0) {
                matching_order.push_back(i);
            }
        }
    }

    void DFSDH::OldMatchingOrder() {
        int N = G1->GetNumVertices();
        std::vector<int> T(N, 0), w(N, 0);
        auto vlabel_freq = G1->GetVertexLabelFrequency();
        auto elabel_freq = G1->GetEdgeLabelFrequency();
        for (int i = 0; i < N; i++) {
            w[i] -= 2 * vlabel_freq[G1->GetVertexLabel(i)];
            for (int x : G1->GetNeighbors(i)) {
                w[i] -= elabel_freq[G1->GetEdgeLabel(i, x)];
            }
        }
        std::fill(T.begin(), T.end(), 0);
        int max_idx = std::max_element(w.begin(), w.end()) - w.begin();
        std::priority_queue<std::pair<int, int>> weighted_queue;
        weighted_queue.push({w[max_idx], max_idx});
        T[max_idx] = 1;
        while (!weighted_queue.empty()) {
            int u = weighted_queue.top().second;
            weighted_queue.pop();
            matching_order.push_back(u);
            for (int x : G1->GetNeighbors(u)) {
                if (T[x] == 1)
                    continue;
                T[x] = 1;
                weighted_queue.push({w[x], x});
            }
        }
        for (int i = 0; i < N; i++) {
            if (T[i] != 1) {
                matching_order.push_back(i);
            }
        }

        inv_matching_order.resize(G1->GetNumVertices(), -1);
        for (int i = 0; i < G1->GetNumVertices(); i++) {
            inv_matching_order[matching_order[i]] = i;
        }
    }
}  // namespace GraphLib::GraphSimilarity
// note: research log marker 17
