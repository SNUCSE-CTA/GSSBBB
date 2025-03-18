#include "GED/EditDistance.h"
#include "DynamicHungarian/DynamicHungarian.h"

namespace GraphLib::GraphSimilarity {
    using Pair = std::pair<int, int>;
    static const Pair DEFAULT_PAIR = std::make_pair(-1, -1);
    enum Direction {
        PROHIBIT,
        MATCH
    };
    class AStarState {
    public:
        int lower_bound = 0, cost = 0;
        int parent_id = -1, id = -1;
        Direction direction = PROHIBIT;
        Pair p = DEFAULT_PAIR;
        int num_matched = 0, num_prohibited = 0;

        AStarState(int id, int parent_id,
            Pair p = DEFAULT_PAIR, Direction direction = PROHIBIT,
            int lower_bound = 0, int cost = 0)
           : id(id), parent_id(parent_id), p(p), direction(direction),
            lower_bound(lower_bound), cost(cost) {
        }
    };


    class AStarBinaryBranching: public GraphEditDistanceSolver {
        int N = 0, tau = 0;
        int64_t cnt = 0;
        bool flag = false;
        DynamicHungarian hungarian_solver;
        std::vector<std::vector<int>> matrix;
        std::vector<Pair> prohibited_pairs;
        std::vector<AStarState> states;
        HungarianState hungarian_state;

        struct Comparator {
            std::vector<AStarState>& states;
            Comparator(std::vector<AStarState>& s) : states(s) {}
            bool operator()(int a, int b) const {
                if (states[a].lower_bound != states[b].lower_bound)
                    return states[a].lower_bound > states[b].lower_bound;
                if (states[a].num_matched != states[b].num_matched) {
                    return states[a].num_matched < states[b].num_matched;
                }
                return states[a].num_prohibited < states[b].num_prohibited;
            }
        };

        std::priority_queue<int, std::vector<int>, Comparator> priority_queue;

    public:
        AStarBinaryBranching() : priority_queue(Comparator(states)) {}
        void Initialize() {
            N = G2->GetNumVertices();
            hungarian_state.Initialize(N);
            mapping = std::vector<int>(G2->GetNumVertices(), -1);
            inverse_mapping = std::vector<int>(G2->GetNumVertices(), -1);
            matrix = std::vector<std::vector<int>>(N, std::vector<int>(N, -1));
            hungarian_solver.Initialize(N);
            current_best = INF;
            tau = INF;
            cnt = 0;
        }

        Pair GetNextPair(int state_id) {
            int uprime = -1, vprime = -1;
            Pair score = {-1000, -1000};
            for (int i = 0; i < G1->GetNumVertices(); i++) {
                if (mapping[i] != -1) continue;
                int assigned_j = hungarian_state.assignment[i];
                int minimum_weight = INF;
                for (int j = 0; j < G2->GetNumVertices(); j++) {
                    if (j == assigned_j) continue;
                    if (inverse_mapping[j] != -1) continue;
                    int ii = hungarian_state.inverse_assignment[j];
                    minimum_weight = std::min(minimum_weight,
                        -matrix[i][assigned_j]
                        -matrix[ii][j]
                        +matrix[i][j]
                        +matrix[ii][assigned_j]
                    );
                }
                Pair new_score = {minimum_weight,
                    std::min(G1->GetDegree(i), G2->GetDegree(assigned_j))};
                if (new_score > score) {
                    score = new_score;
                    uprime = i;
                    vprime = assigned_j;
                }
            }
            return {uprime, vprime};
        }

        void BuildBranchDistanceMatrix(int state_id) {
            DifferenceVector diff;
            diff.init(20);
            prohibited_pairs.clear();
            hungarian_state.Initialize(N);
            std::fill(mapping.begin(), mapping.end(), -1);
            std::fill(inverse_mapping.begin(), inverse_mapping.end(), -1);
            int current_state_id = states[state_id].parent_id;
            while (current_state_id > 0) {
                // if (state_id == 30) {
                //     fprintf(stderr, "Currently Looking at state %d: ", current_state_id);
                //     fprintf(stderr, "%s, %d, %d\n", states[current_state_id].direction == MATCH? "MATCH" : "PROHIBIT",
                //         states[current_state_id].p.first, states[current_state_id].p.second);
                // }
                if (states[current_state_id].direction == MATCH) {
                    int u = states[current_state_id].p.first;
                    int v = states[current_state_id].p.second;
                    mapping[u] = v;
                    inverse_mapping[v] = u;
                }
                else {
                    prohibited_pairs.push_back(states[current_state_id].p);
                }
                current_state_id = states[current_state_id].parent_id;
            }
            current_state_id = state_id;
            if (states[current_state_id].direction == MATCH) {
                int u = states[current_state_id].p.first;
                int v = states[current_state_id].p.second;
                states[current_state_id].cost += ChildEditCost(u, v);
                mapping[u] = v;
                inverse_mapping[v] = u;
            }
            else {
                if (current_state_id > 0)
                    prohibited_pairs.push_back(states[current_state_id].p);
            }
            for (int u = 0; u < G1->GetNumVertices(); u++) {
                if (mapping[u] != -1) {
                    std::fill(matrix[u].begin(), matrix[u].end(), INF);
                    matrix[u][mapping[u]] = 0;
                    continue;
                }
                for (int v = 0; v < G2->GetNumVertices(); v++) {
                    if (inverse_mapping[v] != -1)
                        matrix[u][v] = INF;
                    else {
                        matrix[u][v] = 0;
                        if (G1->GetVertexLabel(u) != G2->GetVertexLabel(v)) {
                            // if (state_id == 1236 and u == 16 and v == 16) {
                            //     fprintf(stderr, "Label Mismatch cost add +2\n");
                            // }
                            matrix[u][v] += 2;
                        }
                        diff.reset();
                        for (int u_nbr : G1->GetNeighbors(u)) {
                            int u_el = G1->GetEdgeLabel(u, u_nbr);
                            if (mapping[u_nbr] == -1) {
                                // if  (state_id == 1236 and u == 16 and v == 16) {
                                //     fprintf(stderr, "Unbounded Nbr u_%d with label %d found\n", u_nbr, u_el);
                                // }
                                diff.update(u_el, 1);
                            }
                            else {
                                int v_mapping_el = G2->GetEdgeLabel(v, mapping[u_nbr]);
                                if (v_mapping_el != u_el) {
                                    // if (state_id == 1236 and u == 16 and v == 16) {
                                    //     fprintf(stderr, "Bounded Nbr u_%d (mapped: v_%d) with label %d is mismatched, +2\n", u_nbr, mapping[u_nbr], u_el);
                                    // }
                                    matrix[u][v] += 2;
                                }
                            }
                        }
                        for (int v_nbr : G2->GetNeighbors(v)) {
                            int v_el = G2->GetEdgeLabel(v, v_nbr);
                            if (inverse_mapping[v_nbr] == -1) {
                                // if (state_id == 1236 and u == 16 and v == 16) {
                                //     fprintf(stderr, "Unbounded Nbr v_%d with label %d found\n", v_nbr, v_el);
                                // }
                                diff.update(v_el, -1);
                            }
                            else {
                                int u_mapping_el = G1->GetEdgeLabel(u, inverse_mapping[v_nbr]);
                                if (u_mapping_el == -1) {
                                    // if(state_id == 1236 and u == 16 and v == 16) {
                                    //     fprintf(stderr, "Bounded Nbr v_%d (mapped: u_%d) with label %d is mismatched, +2\n", v_nbr, inverse_mapping[v_nbr], v_el);
                                    // }
                                    matrix[u][v] += 2;
                                }
                            }
                        }
                        int inner_distance = diff.GetDifference();
                        matrix[u][v] += inner_distance;
                    }
                }
            }
            if (G1->GetNumVertices() < G2->GetNumVertices()) {

            }
            for (int v = 0; v < G2->GetNumVertices(); v++) {
                int num_mapped_nbrs = 0;
                for (int v_nbr : G2->GetNeighbors(v)) {
                    if (inverse_mapping[v_nbr] != -1) {
                        num_mapped_nbrs++;
                    }
                }
                for (int u = G1->GetNumVertices(); u < G2->GetNumVertices(); u++) {
                    matrix[u][v] = 2 + G2->GetDegree(v) + num_mapped_nbrs;
                }
            }
            for (auto [u, v] : prohibited_pairs) {
                matrix[u][v] = INF;
            }
        }

        void ProcessState(int state_id) {
            // build matrix
            BuildBranchDistanceMatrix(state_id);

            AStarState state = states[state_id];
            // if (state_id == 1236) {
            // PrintMatrix();
            //     exit(1);
            // }

            int hungarian_cost = hungarian_solver.Hungarian(
                true, N, hungarian_state, matrix, mapping, inverse_mapping, state.cost, current_best
            );
            // fprintf(stderr, "ProcessState %d (%s, %d, %d) of LB = %d: "
            //                 "Hungarian cost = %d, Current cost = %d\n",
            //                 state_id, state.direction == MATCH? "MATCH" : "PROHIBIT",
            //                 state.p.first, state.p.second, state.lower_bound,
            //                 hungarian_cost, state.cost);
            // std::cerr << "MAPPING: " <<  mapping << '\n';
            // std::cerr << "INVERSE MAPPING: " <<  inverse_mapping << '\n';
            // std::cerr << "PROHIBITED PAIRS: " <<  prohibited_pairs << '\n';
            // std::cerr << "Hung Assignment: " << hungarian_state.assignment << '\n';
            // std::cerr << "Inverse Assignment: " << hungarian_state.inverse_assignment << '\n';

            int lb = state.cost + ((hungarian_cost + 1) / 2);
            if (lb >= current_best) {
                // fprintf(stderr, "LB %d >= curbest %d, return\n\n", lb, current_best);
                return;
            }

            int ub = ComputeDistance(hungarian_state.assignment, hungarian_state.inverse_assignment);
            current_best = std::min(current_best, ub);
            // fprintf(stderr, "Update current best with %d (best = %d)\n", ub, current_best);
            Pair next_match = GetNextPair(state_id);
            // fprintf(stderr, "NextPair Chosen (%d, %d)\n", next_match.first, next_match.second);

            if (next_match.first == -1) {
                flag = true;
                return;
            }


            // match-branch
            AStarState match_state = AStarState(
                states.size(), state.id, next_match, MATCH,
                lb, state.cost);
            match_state.num_matched = state.num_matched + 1;
            match_state.num_prohibited = state.num_prohibited;
            // if (match_state.num_matched == G1->GetNumVertices()) {
            //     match_state.lower_bound = state.cost + ChildEditCost(next_match.first, next_match.second);
            // }
            states.push_back(match_state);
            priority_queue.push(match_state.id);
            // fprintf(stderr, "[%d] Matchstate %d is with lb = %d, cost %d\n",
            //     state_id, match_state.id, match_state.lower_bound, match_state.cost);

            // prohibit-branch
            matrix[next_match.first][next_match.second] += INF;
            hungarian_state.Initialize(N);
            // PrintMatrix();
            hungarian_cost = hungarian_solver.Hungarian(
                true, N, hungarian_state, matrix, mapping, inverse_mapping, state.cost, current_best
            );
            // std::cerr << "Assignment: " << hungarian_state.assignment << '\n';
            int prohibit_lower_bound = state.cost + ((hungarian_cost + 1) / 2);
            if (prohibit_lower_bound >= current_best) {
                // fprintf(stderr, "Prohibit LB %d >= curbest %d, return\n\n", prohibit_lower_bound, current_best);
                return;
            }
            AStarState prohibit_state = AStarState(
                states.size(), state.id, next_match, PROHIBIT,
                prohibit_lower_bound, state.cost);
            prohibit_state.num_matched = state.num_matched;
            prohibit_state.num_prohibited = state.num_prohibited + 1;
            states.push_back(prohibit_state);
            priority_queue.push(prohibit_state.id);


            // fprintf(stderr, "[%d] ProhibitState %d is with lb = %d, cost %d\n",
            //     state_id,
            //     prohibit_state.id, prohibit_state.lower_bound, prohibit_state.cost);
            // fprintf(stderr, "\n");
            // if (state_id >= 0) exit(0);
        }

        int ComputeGED() {
            Timer timer;
            timer.Start();
            Initialize();
            AStarState root = AStarState(0, -1);
            states.push_back(root);
            priority_queue.push(0);
            while (!priority_queue.empty()) {
                cnt++;
                int state_id  = priority_queue.top();
                priority_queue.pop();
                // fprintf(stderr, "Pop State [%d]: LB = %d, Cost = %d, Curbest = %d\n", state_id, lower_bound,
                //     states[state_id].cost,
                //     current_best);
                if (states[state_id].lower_bound >= current_best) {
                    break;
                }
                ProcessState(state_id);
                if (flag) break;
            }
            // while (!priority_queue.empty()) {
            //     priority_queue.pop();
            // }
            timer.Stop();

            log.AddResult("GED", current_best, RESULT_INT);
            log.AddResult("ElapsedTime", timer.GetTime(), RESULT_DOUBLE_FIXED);
            log.AddResult("SearchTreeNodes", cnt, RESULT_INT64);
            return current_best;
        }


        void PrintMatrix() {
            for (int i = 0; i < G2->GetNumVertices(); i++) {
                std::cerr << "[";
                for (int j = 0; j < G2->GetNumVertices(); j++) {
                    if (matrix[i][j] >= 1000)
                        std::cerr << "X" << (j < G2->GetNumVertices()-1?", ":"");
                    else
                        std::cerr << matrix[i][j] << (j < G2->GetNumVertices()-1?", ":"");
                }
                std::cerr << "]\n";
            }
        };

    };
}
// note: research log marker 61
