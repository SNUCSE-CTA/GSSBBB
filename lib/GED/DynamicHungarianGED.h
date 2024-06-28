#pragma once
#include <stack>

#include "Base/Timer.h"
#include "DynamicHungarian/DynamicHungarian.h"
#include "GED/EditDistance.h"

namespace GraphLib::GraphSimilarity {
    class DynamicHungarianGED : public GraphEditDistanceSolver {
    public:
        DynamicHungarian hungarian_solver;
        std::vector<HungarianState> hungarian_states;
        std::vector<std::vector<int>> matrix;
        std::vector<int> matching_order, inv_matching_order;
        std::vector<char> row;
        std::vector<char> col;
        bool flag = false;
        int current_best = 1e3;

        const int INF = 1e3;
        const int INF2 = 1e5;
        int N = 0;
        int total_cost = 0;

        int depth = -1;
        int tau = -1;

        int cost = 0;  // mapping cost

        bool init = false;

        int64_t cnt = 0;
        int filtering_lb = 0;

        void Initialize() {
            matrix = std::vector(G2->GetNumVertices(), std::vector(G2->GetNumVertices(), 0));
            mapping = std::vector<int>(G2->GetNumVertices(), -1);
            inverse_mapping = std::vector<int>(G2->GetNumVertices(), -1);
            row = std::vector<char>(G1->GetNumVertices(), false);
            col = std::vector<char>(G2->GetNumVertices(), false);
            parikh1.assign(G1->GetNumVertices() + 1,
                           std::vector<int>(global_edge_labels.size() + 1, 0));
            parikh2.assign(G2->GetNumVertices() + 1,
                           std::vector<int>(global_edge_labels.size() + 1, 0));
            N = G2->GetNumVertices();
            hungarian_states = std::vector(G1->GetNumVertices() + 2, HungarianState(N));
            hungarian_solver.Initialize(N);
            depth = -1;
            total_cost = 0;
            cost = 0;  // mapping cost
            flag = false;
            cnt = 0;
            init = false;
            current_best = 1e3;
            tau = -1;
        }

        void ChangeAlphaBeta(std::vector<char> &row, std::vector<char> &col) {
            for (int i = 0; i < row.size(); i++) {
                if (row[i] == true) {
                    hungarian_states[depth + 1].alpha[i] = INF;
                    for (int j = 0; j < G2->GetNumVertices(); j++) {
                        hungarian_states[depth + 1].alpha[i] = std::min(hungarian_states[depth + 1].alpha[i], matrix[i][j] - hungarian_states[depth + 1].beta[j]);
                    }
                }
            }
            for (int j = 0; j < col.size(); j++) {
                if (col[j] == true) {
                    hungarian_states[depth + 1].beta[j] = INF;
                    for (int i = 0; i < G2->GetNumVertices(); i++) {
                        hungarian_states[depth + 1].beta[j] = std::min(hungarian_states[depth + 1].beta[j], matrix[i][j] - hungarian_states[depth + 1].alpha[i]);
                    }
                }
            }
            for (int i = 0; i < G2->GetNumVertices(); i++) {
                int u = i;
                int v = hungarian_states[depth + 1].assignment[u];
                if (v != -1 && hungarian_states[depth + 1].alpha[u] + hungarian_states[depth + 1].beta[v] != matrix[u][v]) {
                    hungarian_states[depth + 1].assignment[u] = -1;
                    hungarian_states[depth + 1].inverse_assignment[v] = -1;
                }
            }
        }

        void ComputeBranchDistanceMatrixDynamic(const int u, const int v, bool update) {
            // fprintf(stderr, "ComputeBranchDistanceMatrixDynamic: %d %d\n", u, v);
            auto &u_nbrs = G1->GetNeighbors(u);
            auto &v_nbrs = G2->GetNeighbors(v);
            fill(row.begin(), row.end(), 0);
            fill(col.begin(), col.end(), 0);
            /* Update cost matrix */
            // CASE1 and CASE2: _u in Nbr(u)
            for (int i = 0; i < (int)u_nbrs.size(); ++i) {
                const int _u = u_nbrs[i];
                if (mapping[_u] != -1) {
                    continue;
                }
                for (int _v = 0; _v < (int)G2->GetNumVertices(); ++_v) {
                    if (inverse_mapping[_v] != -1) {
                        // || G2->GetEdgeLabel(v, _v) != -1) {
                        continue;
                    }

                    const int l1 = G1->GetEdgeLabel(u, _u);
                    const int l2 = G2->GetEdgeLabel(v, _v);

                    if (l1 != l2) {
                        int delta = 0;
                        delta += 2;

                        delta -= ~(l1 | l2) || (parikh1[_u][0] > parikh2[_v][0]);
                        if (l1 != -1 && parikh1[_u][l1] <= parikh2[_v][l1]) {
                            delta += 1;
                        }
                        if (l2 != -1 && parikh1[_u][l2] >= parikh2[_v][l2]) {
                            delta += 1;
                        }
                        if (update) {
                            matrix[_u][_v] += delta;
                        } else {
                            matrix[_u][_v] -= delta;
                        }
                    }
                }
                row[_u] = true;
            }
            // CASE3: _u not in Nbr(u) and _v in Nbr(v)
            for (int _u = 0; _u < (int)G1->GetNumVertices(); ++_u) {
                if (mapping[_u] != -1 || G1->GetEdgeLabel(u, _u) != -1) {
                    continue;
                }
                for (int j = 0; j < (int)v_nbrs.size(); ++j) {
                    const int _v = v_nbrs[j];
                    if (inverse_mapping[_v] != -1) {
                        continue;
                    }
                    const int l1 = G1->GetEdgeLabel(u, _u);
                    const int l2 = G2->GetEdgeLabel(v, _v);
                    if (l1 != l2) {
                        int delta = 0;
                        delta += 2;

                        delta -= (parikh1[_u][0] < parikh2[_v][0]);
                        if (parikh1[_u][l2] >= parikh2[_v][l2]) {
                            delta += 1;
                        }
                        if (update) {
                            matrix[_u][_v] += delta;
                        } else {
                            matrix[_u][_v] -= delta;
                        }
                    }
                }
            }

            for (int j = 0; j < (int)v_nbrs.size(); ++j) {
                const int _v = v_nbrs[j];
                if (inverse_mapping[_v] != -1) {
                    continue;
                }
                for (int _u = G1->GetNumVertices(); _u < G2->GetNumVertices(); ++_u) {
                    if (update) {
                        matrix[_u][_v] += 1;
                    } else {
                        matrix[_u][_v] -= 1;
                    }
                }
                col[_v] = true;
            }
            if (update) {
                ChangeAlphaBeta(row, col);
            }
        }

        std::pair<int, int> LowerBound() {
            int ub = current_best + 1, lb = 0;
            if (!init) {
                total_cost = hungarian_solver.Hungarian(true, N, hungarian_states[depth + 1], matrix,
                                                        mapping, inverse_mapping, cost, tau);
                init = true;
            } else {
                total_cost = hungarian_solver.Hungarian(false, N, hungarian_states[depth + 1], matrix,
                                                        mapping, inverse_mapping, cost, tau);
            }
            lb = cost + ((total_cost + 1) / 2);
            if (lb <= tau)
                ub = ComputeDistance(hungarian_states[depth + 1].assignment, hungarian_states[depth + 1].inverse_assignment);
            return {lb, ub};
        }

        void ChangeCostINF(int u, int v) {
            matrix[u][v] += INF;
            hungarian_states[depth + 1].Disconnect(u, v, true);
        }

        void Match(int u, int v) {
            for (int j = 0; j < N; j++) {
                if (j != v) {
                    ChangeCostINF(u, j);
                }
            }
            for (int i = 0; i < N; i++) {
                if (i != u) {
                    ChangeCostINF(i, v);
                }
            }
        }

        void RemoveMatch(int u, int v) {
            // have to fix
            for (int j = 0; j < N; j++) {
                if (j != v) {
                    matrix[u][j] -= INF;
                }
            }
            for (int i = 0; i < N; i++) {
                if (i != u) {
                    matrix[i][v] -= INF;
                }
            }
        }

        void PrintMatrix() {
            for (int i = 0; i < G2->GetNumVertices(); i++) {
                std::cerr << "[";
                for (int j = 0; j < G2->GetNumVertices(); j++) {
                    if (matrix[i][j] >= 1000)
                        std::cerr << "X" << (j < G2->GetNumVertices() - 1 ? ", " : "");
                    else
                        std::cerr << matrix[i][j] << (j < G2->GetNumVertices() - 1 ? ", " : "");
                }
                std::cerr << "]\n";
            }
        };
    };
}  // namespace GraphLib::GraphSimilarity
// note: research log marker 21
