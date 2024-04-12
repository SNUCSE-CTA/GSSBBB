#pragma once
#include "Base/Base.h"

#include "Base/PrettyPrint.h"
struct HungarianState {
    int n = 0; 
    std::vector<int> assignment, inverse_assignment;
    std::vector<int> alpha, beta;

    HungarianState(int N = 0) {
        Initialize(N);
    }

    void Initialize(int N) {
        if (n != N) {
            n = N;
            assignment.resize(N, -1);
            inverse_assignment.resize(N, -1);
            alpha.resize(N, 0);
            beta.resize(N, 0);
        }
        else {
            std::fill(assignment.begin(), assignment.end(), -1);
            std::fill(inverse_assignment.begin(), inverse_assignment.end(), -1);
            std::fill(alpha.begin(), alpha.end(), 0);
            std::fill(beta.begin(), beta.end(), 0);
        }
    }

    void Disconnect(int u, int v, bool check = false) {
        if (!check or assignment[u] == v) {
            assignment[u] = -1;
            inverse_assignment[v] = -1;
        }
    }
};

struct DynamicHungarian {
    char *visX = nullptr;
    char *visY = nullptr;
    int *slack = nullptr;
    int *slackmy = nullptr;
    unsigned int *prev = nullptr;
    unsigned int *queue = nullptr;

    void Initialize(int N) {
        visX = new char[N]();
        std::memset(visX, 0, sizeof(char) * N);
        visY = new char[N]();
        std::memset(visY, 0, sizeof(char) * N);
        slack = new int[N]();
        std::memset(slack, 0, sizeof(int) * N);
        slackmy = new int[N]();
        std::memset(slackmy, 0, sizeof(int) * N);
        prev = new unsigned int[N]();
        std::memset(prev, 0, sizeof(unsigned int) * N);
        queue = new unsigned int[N]();
        std::memset(queue, 0, sizeof(unsigned int) * N);
    }

    int Hungarian(bool initialization, int n,
        HungarianState &hungarian_state,
        std::vector<std::vector<int>> &matrix,
        std::vector<int> &mapping, std::vector<int> &inverse_mapping,
        int cost, int tau
    ) {
        int *mx = hungarian_state.assignment.data();
        int *my = hungarian_state.inverse_assignment.data();
        int *lx = hungarian_state.alpha.data();
        int *ly = hungarian_state.beta.data();
        if (initialization) {
            // Initialization
            memset(mx, -1, sizeof(int) * n);
            memset(my, -1, sizeof(int) * n);
            memset(ly, 0, sizeof(int) * n);
            for (int i = 0; i < n; i++) {
                lx[i] = INF;
                for (int j = 0; j < n; j++)
                    if (matrix[i][j] < lx[i])
                        lx[i] = matrix[i][j];
                for (int j = 0; j < n; j++)
                    if (my[j] == -1 && matrix[i][j] == lx[i]) {
                        mx[i] = j;
                        my[j] = i;
                        // fprintf(stderr, "OK found init %d->%d (%d)\n",i, j, lx[i]);
                        break;
                    }
            }
        }
        int hungarian_cost = 0;
        for (int i = 0; i < n; i++) {
            // fprintf(stderr, "Mapping %d = %d\n", i, mapping[i]);
            if (mapping[i] == -1) {
                hungarian_cost += lx[i];
                // fprintf(stderr, "lx[%d] = %d\n", i, lx[i]);
            }
        }
        for (int j = 0; j < n; j++) {
            if (inverse_mapping[j] == -1) {
                hungarian_cost += ly[j];
            }
        }
        // fprintf(stderr, "after init hung cost = %d\n", hungarian_cost);
        int lb = cost + ((hungarian_cost + 1) / 2);
        if (lb > tau)
        // if (lb >= current_best)
        {
            return hungarian_cost;
        }

        for (int u = n - 1; u >= 0; u--)
            if (mx[u] == -1) {
                // Augmentation
                memset(visX, 0, sizeof(char) * n);
                memset(visY, 0, sizeof(char) * n);
                int q_n = 1;
                queue[0] = u;
                visX[u] = 1;
                for (int i = 0; i < n; i++) {
                    slack[i] = matrix[u][i] - lx[u] - ly[i];
                    slackmy[i] = u;
                }
                int target = n, X;
                while (true) {
                    for (int i = 0; i < q_n && target == n; i++) {
                        int v = queue[i];
                        for (int j = 0; j < n; j++)
                            if (!visY[j] && matrix[v][j] == lx[v] + ly[j]) {
                                if (my[j] == -1) {
                                    X = v;
                                    target = j;
                                    break;
                                }
                                visY[j] = 1;
                                X = my[j];
                                visX[X] = 1;
                                prev[X] = v;
                                queue[q_n++] = X;

                                for (int k = 0; k < n; k++)
                                    if (!visY[k] && matrix[X][k] - lx[X] - ly[k] < slack[k]) {
                                        slack[k] = matrix[X][k] - lx[X] - ly[k];
                                        slackmy[k] = X;
                                    }
                            }
                    }
                    if (target != n)
                        break;

                    q_n = 0;
                    int delta = INF;
                    for (int i = 0; i < n; i++)
                        if (!visY[i] && slack[i] < delta)
                            delta = slack[i];
                    for (int i = 0; i < n; i++) {
                        if (visX[i]) {
                            lx[i] += delta;
                            hungarian_cost += delta;
                        }
                        if (visY[i]) {
                            ly[i] -= delta;
                            hungarian_cost -= delta;
                        } else
                            slack[i] -= delta;
                    }
                    // fprintf(stderr, "Hung cost = %d after finding aug-path\n", hungarian_cost);
                    lb = cost + ((hungarian_cost + 1) / 2);
                    if (lb > tau)
                    // if (lb >= current_best)
                    {
                        return hungarian_cost;
                    }

                    for (int i = 0; i < n; i++)
                        if (!visY[i] && slack[i] == 0) {
                            if (my[i] == -1) {
                                X = slackmy[i];
                                target = i;
                                break;
                            }
                            visY[i] = 1;
                            if (!visX[my[i]]) {
                                X = my[i];
                                visX[X] = 1;
                                prev[X] = slackmy[i];
                                queue[q_n++] = X;

                                // int *tt_array = cost + X*n;
                                for (int k = 0; k < n; k++)
                                    if (!visY[k] && matrix[X][k] - lx[X] - ly[k] < slack[k]) {
                                        slack[k] = matrix[X][k] - lx[X] - ly[k];
                                        slackmy[k] = X;
                                    }
                            }
                        }
                }

                while (true) {
                    int ty = mx[X];
                    mx[X] = target;
                    my[target] = X;
                    if (X == u)
                        break;

                    X = prev[X];
                    target = ty;
                }
            }
        return hungarian_cost;
    }
};
// note: research log marker 14
