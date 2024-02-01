#pragma once
#include <fcntl.h>
#include <sys/wait.h>
#include <unistd.h>

#include "Base/PrettyPrint.h"
#include "Base/Timer.h"
#include "GED/DynamicHungarianGED.h"

namespace GraphLib::GraphSimilarity {
    class BinaryBranching : public DynamicHungarianGED {
    public:
        std::vector<int> temp_visited;
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

        int run_with_stdin(const char* file, const std::vector<const char*>& args,
                           const char* input_path) {
            // argv 준비
            std::vector<char*> argv;
            argv.push_back(const_cast<char*>(file));
            for (auto* a : args) argv.push_back(const_cast<char*>(a));
            argv.push_back(nullptr);

            pid_t pid = fork();
            if (pid < 0) {
                std::cerr << "fork: " << std::strerror(errno) << "\n";
                return -1;
            }
            if (pid == 0) {
                // 자식: stdin을 input_path로 리다이렉트
                if (input_path) {
                    int fd = open(input_path, O_RDONLY);
                    if (fd < 0) {
                        std::perror("open");
                        _exit(127);
                    }
                    if (dup2(fd, STDIN_FILENO) < 0) {
                        std::perror("dup2");
                        _exit(127);
                    }
                    close(fd);  // fd는 더 이상 필요 없음
                }
                execvp(file, argv.data());
                std::perror("execvp");
                _exit(127);
            }

            int status = 0;
            if (waitpid(pid, &status, 0) < 0) {
                std::cerr << "waitpid: " << std::strerror(errno) << "\n";
                return -1;
            }
            return WIFEXITED(status)     ? WEXITSTATUS(status)
                   : WIFSIGNALED(status) ? 128 + WTERMSIG(status)
                                         : -1;
        }

        void read_auto(std::string graph_id, bool is_g1) {
            std::string path = "../data/aids100/" + graph_id + ".txt";
            const char* cpath = path.c_str();
            Timer timer;
            timer.Start();
            int code = run_with_stdin("../driver/auto", {}, cpath);
            timer.Stop();
            std::cout << "exit=" << code << ", time=" << timer.GetTime() << std::endl;

            std::string graph_path = "../data/auto/aids100/" + graph_id + ".txt";
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

            graph_path = "../data/auto/aids100-groups/" + graph_id + ".txt";
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
                } else {
                    v_groups.push_back(group);
                    v_group_sizes.push_back(group.size());
                    for (int v : group) {
                        v_vertex_to_group[v] = v_groups.size() - 1;
                    }
                }
            }
            infile.close();
        }
#endif
        virtual std::pair<int, int> GetNextPair() { return {-1, -1}; };

        void AddMapping(int u, int v) {
            mapping[u] = v;
            inverse_mapping[v] = u;
            Match(u, v);
            ComputeBranchDistanceMatrixDynamic(u, v, 1);
            ComputeParikhVector();
        }

        virtual int Verify() {
            Timer timer;
            timer.Start();

            num_nodes = 0;
            current_best = 1e9;
            temp_visited = std::vector<int>(N, 0);

            Initialize();

            ComputeBranchDistanceMatrixInitial(matrix);
            ComputeParikhVector();
            depth = -1;
            tau = this->threshold;
            auto [lb, ub] = LowerBound();
            current_best = ub;
            if (lb > tau or ub <= tau) {
                goto finished;
            }
            cnt++;
            DFS(-1, -1);
        finished:
            timer.Stop();
            log.AddResult("GED", current_best, RESULT_INT);
            log.AddResult("ElapsedTime", timer.GetTime(), RESULT_DOUBLE_FIXED);
            log.AddResult("SearchTreeNodes", cnt, RESULT_INT64);
            return current_best;
        }

        virtual int ComputeGED() {
            Timer timer;
            timer.Start();

            num_nodes = 0;
            current_best = 1e9;
            temp_visited = std::vector<int>(N, 0);

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
                DFS(-1, -1);
                tau++;
            }
            timer.Stop();
            log.AddResult("GED", current_best, RESULT_INT);
            log.AddResult("ElapsedTime", timer.GetTime(), RESULT_DOUBLE_FIXED);
            log.AddResult("SearchTreeNodes", cnt, RESULT_INT64);
            return current_best;
        }

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
                // std::cout << u << " " << v << " "<<chcost << "\n";
                cost += chcost;
                ComputeBranchDistanceMatrixDynamic(u, v, 1);
                UpdateParikhVector(u, v);
            }
            int lb, ub;
            std::stack<std::pair<int, int>> uv;
            // std::vector<std::vector<int>> u_adj_vertices(u_list.size());
            // std::vector<std::vector<int>> v_adj_vertices(v_list.size());
            // int prev_u_group_size = 0, prev_v_group_size = 0;
            while (true) {
                std::tie(lb, ub) = LowerBound();
                bool valid = true;
                if (lb > tau) {
                    break;
                }
                if (ub < current_best) {
                    current_best = ub;
                    // std::cerr << ub << '\n';
                    // std::cerr << hungarian_states[depth+1].assignment << '\n';
                    // std::cerr << hungarian_states[depth+1].inverse_assignment << '\n';
                    if (current_best <= tau) {
                        flag = true;
                        return;
                    }
                }
                fill(temp_visited.begin(), temp_visited.end(), 0);

                int uprime, vprime;
#ifdef AUTO
                // if (depth == 0) {
                //     std::tie(uprime, vprime) = GetNextPairFirst(u_list, v_list);
                //     std::cout << "uprime: " << uprime << ", vprime: " << vprime << "\n";
                //     std::cout << "u_vertex_to_group[uprime]: " << u_vertex_to_group[uprime] << ", v_vertex_to_group[vprime]: " << v_vertex_to_group[vprime] << "\n";
                // } else
#endif
                {
                    std::tie(uprime, vprime) = GetNextPair();
                }
                // if (depth == 0) {
                //     std::cout << "uprime: " << uprime << ", vprime: " << vprime << "\n";
                // }

                // prev_u_group_size = u_groups.size();
                // for (int nu : G1->GetNeighbors(uprime)) {
                //     u_adj_vertices[u_vertex_to_group[nu]].push_back(nu);
                // }
                // for (int i = 0; i < u_adj_vertices.size(); i++) {
                //     if (u_adj_vertices[i].size() != u_group_sizes[i]) {
                //         std::vector<int> new_group;
                //         for (int nu : u_adj_vertices[i]) {
                //             new_group.push_back(nu);
                //             u_vertex_to_group[nu] = u_groups.size();
                //             auto it = std::find(u_groups[i].begin(), u_groups[i].end(),
                //                                 nu);
                //             if (it != u_groups[i].end()) {
                //                 u_groups[i].erase(it);
                //             }
                //         }
                //         u_group_sizes[i] -= new_group.size();
                //         u_groups.push_back(new_group);
                //         u_group_sizes.push_back(new_group.size());
                //     }
                // }

                // prev_v_group_size = v_groups.size();
                // for (int nu : G2->GetNeighbors(vprime)) {
                //     v_adj_vertices[v_vertex_to_group[nu]].push_back(nu);
                // }
                // for (int i = 0; i < v_adj_vertices.size(); i++) {
                //     if (v_adj_vertices[i].size() != v_group_sizes[i]) {
                //         std::vector<int> new_group;
                //         for (int nu : v_adj_vertices[i]) {
                //             new_group.push_back(nu);
                //             v_vertex_to_group[nu] = v_groups.size();
                //             auto it = std::find(v_groups[i].begin(), v_groups[i].end(),
                //                                 nu);
                //             if (it != v_groups[i].end()) {
                //                 v_groups[i].erase(it);
                //             }
                //         }
                //         v_group_sizes[i] -= new_group.size();
                //         v_groups.push_back(new_group);
                //         v_group_sizes.push_back(new_group.size());
                //     }
                // }

                DFS(uprime, vprime);

                // for (int i = 0; i < u_adj_vertices.size(); i++) {
                //     if (u_adj_vertices[i].size() != 0 && u_vertex_to_group[u_adj_vertices[i][0]] >= prev_u_group_size) {
                //         for (int nu : u_adj_vertices[i]) {
                //             u_groups[i].push_back(nu);
                //             u_vertex_to_group[nu] = i;
                //         }
                //         u_group_sizes[i] += u_adj_vertices[i].size();
                //     }
                //     u_adj_vertices[i].clear();
                // }
                // u_groups.erase(u_groups.begin() + prev_u_group_size, u_groups.end());
                // u_group_sizes.erase(u_group_sizes.begin() + prev_u_group_size, u_group_sizes.end());

                // for (int i = 0; i < v_adj_vertices.size(); i++) {
                //     if (v_adj_vertices[i].size() != 0 && v_vertex_to_group[v_adj_vertices[i][0]] >= prev_v_group_size) {
                //         for (int nu : v_adj_vertices[i]) {
                //             v_groups[i].push_back(nu);
                //             v_vertex_to_group[nu] = i;
                //         }
                //         v_group_sizes[i] += v_adj_vertices[i].size();
                //     }
                //     v_adj_vertices[i].clear();
                // }
                // v_groups.erase(v_groups.begin() + prev_v_group_size, v_groups.end());
                // v_group_sizes.erase(v_group_sizes.begin() + prev_v_group_size, v_group_sizes.end());

                if (flag)
                    return;
                uv.emplace(uprime, vprime);
                matrix[uprime][vprime] += INF2;
#ifdef AUTO
                if (depth == 0) {
                    for (int& u_auto : u_groups[u_vertex_to_group[uprime]]) {
                        for (int& v_auto : v_groups[v_vertex_to_group[vprime]]) {
                            matrix[u_auto][v_auto] += INF2;
                        }
                    }
                }
#endif
                hungarian_states[depth + 1].assignment[uprime] = -1;
                hungarian_states[depth + 1].inverse_assignment[vprime] = -1;
            }

            while (!uv.empty()) {
                auto [uprime, vprime] = uv.top();
                uv.pop();
                matrix[uprime][vprime] -= INF2;
#ifdef AUTO
                if (depth == 0) {
                    for (int& u_auto : u_groups[u_vertex_to_group[uprime]]) {
                        for (int& v_auto : v_groups[v_vertex_to_group[vprime]]) {
                            matrix[u_auto][v_auto] -= INF2;
                        }
                    }
                }
#endif
            }

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
    };

    class BinaryBranchingHighestContribution : public BinaryBranching {
        std::pair<int, int> GetNextPair() {
            int uprime = -1, vprime = -1, lb_prime = -1;

            for (int i = 0; i < G1->GetNumVertices(); i++) {
                if (mapping[i] != -1 || temp_visited[i]) continue;
                int j = hungarian_states[depth + 1].assignment[i];
                matrix[i][j] += INF2;
                depth++;
                hungarian_states[depth + 1] = hungarian_states[depth];
                hungarian_states[depth + 1].Disconnect(i, j);
                int hung = hungarian_solver.Hungarian(0, N, hungarian_states[depth + 1],
                                                      matrix, mapping, inverse_mapping, cost, tau);
                int new_lb = cost + ((hung + 1) / 2);
                if (new_lb > lb_prime) {
                    lb_prime = new_lb;
                    uprime = i;
                    vprime = j;
                }
                matrix[i][j] -= INF2;
                depth--;
                if (lb_prime > tau) break;
                for (int k = i + 1; k < G1->GetNumVertices(); k++) {
                    if (hungarian_states[depth + 2].assignment[k] != hungarian_states[depth + 1].assignment[k]) {
                        temp_visited[k] = true;
                    }
                }
            }
            return {uprime, vprime};
        }
    };

    class BinaryBranchingHighestWeightGap : public BinaryBranching {
        std::pair<int, int> GetNextPair() {
            int uprime = -1, vprime = -1;
            int delta = -1000;
            std::pair<int, int> score = {-1000, -1000};
            for (int i = 0; i < G1->GetNumVertices(); i++) {
                if (mapping[i] != -1) continue;
                int assigned_j = hungarian_states[depth + 1].assignment[i];
                int minimum_weight = INF;
                for (int j = 0; j < G2->GetNumVertices(); j++) {
                    if (j == assigned_j) continue;
                    if (inverse_mapping[j] != -1) continue;
                    minimum_weight = std::min(minimum_weight, matrix[i][j]);
                }
                std::pair<int, int> new_score = {minimum_weight,
                                                 std::min(G1->GetDegree(i), G2->GetDegree(assigned_j))};
                if (new_score > score) {
                    score = new_score;
                    uprime = i;
                    vprime = assigned_j;
                }
            }
            return {uprime, vprime};
        }
    };

    class BinaryBranchingBestSwap : public BinaryBranching {
        std::pair<int, int> GetNextPair() {
            int uprime = -1, vprime = -1;
            int best_weight = -1000;
            int best_degree = -1000;
            const int numVerticesG1 = G1->GetNumVertices();
            const int numVerticesG2 = G2->GetNumVertices();
            auto& state = hungarian_states[depth + 1];

            for (int i = 0; i < numVerticesG1; i++) {
                if (mapping[i] != -1) continue;
                const int assigned_j = state.assignment[i];
                int candidate_weight = INF;
                for (int j = 0; j < numVerticesG2; j++) {
                    if (j == assigned_j) continue;
                    if (inverse_mapping[j] != -1) continue;

                    int ii = state.inverse_assignment[j];
                    int candidate = -matrix[ii][j] + matrix[i][j] + matrix[ii][assigned_j];
                    if (candidate < candidate_weight) {
                        candidate_weight = candidate;
                    }
                }
                candidate_weight -= matrix[i][assigned_j];

                int deg1 = G1->GetDegree(i);
                int deg2 = G2->GetDegree(assigned_j);
                int current_min_degree = (deg1 < deg2) ? deg1 : deg2;
                if (candidate_weight > best_weight) {
                    best_weight = candidate_weight;
                    best_degree = current_min_degree;
                    uprime = i;
                    vprime = assigned_j;
                } else if (candidate_weight == best_weight && current_min_degree > best_degree) {
                    best_degree = current_min_degree;
                    uprime = i;
                    vprime = assigned_j;
                }
            }

            return {uprime, vprime};
        }

#ifdef AUTO
        std::pair<int, int> GetNextPairFirst(std::vector<int>& u_list, std::vector<int>& v_list) {
            int uprime = -1, vprime = -1;
            int best_weight = -1000;
            int best_degree = -1000;
            auto& state = hungarian_states[depth + 1];

            for (int i : u_list) std::cout << i << " ";
            std::cout << "\n";
            for (int j : v_list) std::cout << j << " ";
            std::cout << "\n";

            for (int i = 0; i < G1->GetNumVertices(); i++) {
                for (int j = 0; j < G2->GetNumVertices(); j++) {
                    if (matrix[i][j] >= INF2)
                        std::cout << "INF ";
                    else
                        std::cout << matrix[i][j] << " ";
                }
                std::cout << "\n";
            }

            for (int i : u_list) {
                if (mapping[i] != -1) continue;
                const int assigned_j = state.assignment[i];
                int candidate_weight = INF;
                for (int j : v_list) {
                    if (j == assigned_j) continue;
                    // if (inverse_mapping[j] != -1) continue;

                    int ii = state.inverse_assignment[j];
                    int candidate = -matrix[ii][j] + matrix[i][j] + matrix[ii][assigned_j];
                    if (candidate < candidate_weight) {
                        candidate_weight = candidate;
                    }
                }
                candidate_weight -= matrix[i][assigned_j];

                int deg1 = G1->GetDegree(i);
                int deg2 = G2->GetDegree(assigned_j);
                int current_min_degree = (deg1 < deg2) ? deg1 : deg2;
                if (candidate_weight > best_weight) {
                    best_weight = candidate_weight;
                    best_degree = current_min_degree;
                    uprime = i;
                    vprime = assigned_j;
                } else if (candidate_weight == best_weight && current_min_degree > best_degree) {
                    best_degree = current_min_degree;
                    uprime = i;
                    vprime = assigned_j;
                }
            }

            return {uprime, vprime};
        }
#endif
    };
}  // namespace GraphLib::GraphSimilarity