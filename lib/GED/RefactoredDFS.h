#pragma once
#include <stack>
#include "Base/Timer.h"
#include "GED/EditDistance.h"
#include "DynamicHungarian/DynamicHungarian.h"

namespace GraphLib::GraphSimilarity {
    class DFSRefactored : public GraphEditDistanceSolver {
    public:
        DynamicHungarianSolver *hungarian_solver;

        bool flag = false;
        int current_best = 1e3;
        int N = 0;
        int acc = 0;
        int total_cost = 0;
        int depth = -1;
        int tau = -1;

        int cost = 0; // mapping cost
        int64_t cnt = 0;


        int filtering_lb = 0;

        void Initialize() {
            mapping = std::vector<int>(G2->GetNumVertices(), -1);
            inverse_mapping = std::vector<int>(G2->GetNumVertices(), -1);
            parikh1.assign(G1->GetNumVertices() + 1,
                           std::vector<int>(global_edge_labels.size(), 0));
            parikh2.assign(G2->GetNumVertices() + 1,
                           std::vector<int>(global_edge_labels.size(), 0));
            N = G2->GetNumVertices();

            depth = -1;
            acc = 0;
            total_cost = 0;
            cost = 0; // mapping cost
            flag = false;
            cnt = 0;
            current_best = 1e3;
            tau = -1;
        }

        void ComputeBranchDistanceMatrixDynamic(const int u, const int v) {
            auto &u_nbrs = G1->GetNeighbors(u);
            auto &v_nbrs = G2->GetNeighbors(v);
            /* Update cost matrix */
            // CASE1 and CASE2: _u in Nbr(u)
            for (int i = 0; i < (int) u_nbrs.size(); ++i) {
                const int _u = u_nbrs[i];
                if (mapping[_u] != -1) {
                    continue;
                }
                for (int _v = 0; _v < (int) G2->GetNumVertices(); ++_v) {
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
                        hungarian_solver->UpdateCost(_u, _v, delta, depth);
                    }
                }
                hungarian_solver->MarkChange(_u, true);
            }
            // CASE3: _u not in Nbr(u) and _v in Nbr(v)
            for (int _u = 0; _u < (int) G1->GetNumVertices(); ++_u) {
                if (mapping[_u] != -1 || G1->GetEdgeLabel(u, _u) != -1) {
                    continue;
                }
                for (int j = 0; j < (int) v_nbrs.size(); ++j) {
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
                        hungarian_solver->UpdateCost(_u, _v, delta, depth);
                    }
                }
            }

            for (int j = 0; j < (int) v_nbrs.size(); ++j) {
                const int _v = v_nbrs[j];
                if (inverse_mapping[_v] != -1) {
                    continue;
                }
                for (int _u = G1->GetNumVertices(); _u < G2->GetNumVertices(); ++_u) {
                    hungarian_solver->UpdateCost(_u, _v, 1, depth);
                }
                hungarian_solver->MarkChange(_v, false);
            }
            hungarian_solver->ChangeAlphaBeta();
        }

        std::pair<int, int> LowerBound(bool init, int hungarian_threshold=-1) {
            int ub = current_best + 1, lb = 0;
            if (hungarian_threshold == -1) hungarian_threshold = 2 * (tau - cost);
            total_cost = hungarian_solver->Solve(
                init, mapping, inverse_mapping,
                2 * (tau - cost));
            // printf("LB by Hungarian: %d (init = %d)\n", total_cost, init);
            lb = cost + ((total_cost + 1) / 2);
            if (lb <= tau) {
                ub = ComputeDistance(hungarian_solver->GetAssignment(),
                                     hungarian_solver->GetInverseAssignment());
            }
            return {lb, ub};
        }

        // cost + (t + 1) / 2 <= tau
        // (t + 1) / 2 <= tau - cost
        // (t + 1) <= 2 * (tau - cost)
        // IDS dfs
        void DFS(int u, int v) {
            cnt++;
            depth++;
            int chcost = 0;
            if (depth == G1->GetNumVertices()) {
                depth--;
                return;
            }
            hungarian_solver->MemorizeCurrentState();
            if (depth != 0) {
                mapping[u] = v;
                inverse_mapping[v] = u;
                hungarian_solver->Match(u, v);
                chcost = ChildEditCost(u, v);
                cost += chcost;
                ComputeBranchDistanceMatrixDynamic(u, v);
                UpdateParikhVector(u, v);
            }
            int lb, ub;
            int uprime = matching_order[depth];
            for (int i = 0; i < G2->GetNumVertices() - depth; i++) {
                std::tie(lb, ub) = LowerBound(false);
                int vprime = hungarian_solver->GetAssigned(uprime);
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

                DFS(uprime, vprime);

                if (flag)
                    return;
                hungarian_solver->Block(uprime, vprime);
            }

            for (int i = 0; i < G2->GetNumVertices(); i++)
                hungarian_solver->UnBlock(uprime, i);
            if (depth != 0) {
                RestoreParikhVector(u, v);
                hungarian_solver->RevertUpdateCost(depth);
                mapping[u] = -1;
                inverse_mapping[v] = -1;
                hungarian_solver->RemoveMatch(u, v);
                cost -= chcost;
            }
            depth--;
            hungarian_solver->RestoreFromMemorizedState();
        }

        void PrepareGED() {
            matching_order.clear();
            ComputeMatchingOrder();
            current_best_mapping.clear();
            current_best_mapping.resize(G1->GetNumVertices(), -1);
            num_nodes = 0;
            current_best = 1e9;
        }

        int ComputeGED() {
            Timer timer;
            timer.Start();
            PrepareGED();
            Initialize();

            hungarian_solver = new DynamicHungarianSolver(G2->GetNumVertices());
            ComputeBranchDistanceMatrixInitial(hungarian_solver->cost_matrix);
            ComputeParikhVector();
            tau = 0;
            depth = -1;
            auto [lb, ub] = LowerBound(true, INF);
            current_best = ub;
            cnt++;
            tau = lb;
            // tau = threshold;
            // current_best = std::min(current_best, tau+1);
            // IDS (Main loop)
            while (true) {
                if (current_best <= tau) break;
                cnt--;
                depth = -1;
                DFS(-1, -1);
                tau++;
            }
            timer.Stop();
            log.AddResult("GED", current_best, RESULT_INT);
            log.AddResult("ElapsedTime", timer.GetTime(), RESULT_DOUBLE_FIXED);
            log.AddResult("SearchTreeNodes", cnt, RESULT_INT64);
            return tau;
        }

        void PrintMatrix() {
            for (int i = 0; i < G2->GetNumVertices(); i++) {
                for (int j = 0; j < G2->GetNumVertices(); j++) {
                    std::cout << hungarian_solver->cost_matrix[i][j] << "   ";
                }
                std::cout << "\n";
            }
            std::cout << "-------------------------------------------------\n";
        };
    };
} // namespace GraphLib::GED
// note: research log marker 23
