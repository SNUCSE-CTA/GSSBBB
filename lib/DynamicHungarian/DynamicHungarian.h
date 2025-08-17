#pragma once
#include "Base/Base.h"
#include "Base/PrettyPrint.h"
struct HungarianState {
    std::vector<int> assignment, inverse_assignment;
    std::vector<int> alpha, beta;

    HungarianState(int N = 0) {
        Initialize(N);
    }

    void Initialize(int N) {
        assignment.assign(N, -1);
        inverse_assignment.assign(N, -1);
        alpha.assign(N, 0);
        beta.assign(N, 0);
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
        visX = new char[N];
        visY = new char[N];
        slack = new int[N];
        slackmy = new int[N];
        prev = new unsigned int[N];
        queue = new unsigned int[N];
    }
    ~DynamicHungarian() {
        delete[] visX;
        delete[] visY;
        delete[] slack;
        delete[] slackmy;
        delete[] prev;
        delete[] queue;
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
                while (1) {
                    bool found = false;
                    for (int i = 0; i < q_n && target == n; i++) {
                        int curRow = queue[i];
                        int lx_cur = lx[curRow];
                        for (int j = 0; j < n; j++) {
                            if (!visY[j]) {
                                int newSlack = matrix[curRow][j] - lx_cur - ly[j];
                                if (newSlack < slack[j]) { slack[j] = newSlack; slackmy[j] = curRow; }
                                if (slack[j] == 0) {
                                    visY[j] = 1;
                                    if (my[j] == -1) { target = j; X = slackmy[j]; found = true; break; }
                                    else {
                                        int nextRow = my[j];
                                        if (!visX[nextRow]) { visX[nextRow] = 1; prev[nextRow] = slackmy[j]; queue[q_n++] = nextRow; }
                                    }
                                }
                            }
                        }
                        if (found)
                            break;
                    }
                    if (target != n)
                        break;
                    int delta = INF;
                    for (int j = 0; j < n; j++)
                        if (!visY[j] && slack[j] < delta)
                            delta = slack[j];
                    for (int i = 0; i < n; i++)
                        if (visX[i]) { lx[i] += delta; hungarian_cost += delta; }
                    for (int j = 0; j < n; j++) {
                        if (visY[j]) { ly[j] -= delta; hungarian_cost -= delta; }
                        else { slack[j] -= delta; }
                    }
                    lb = cost + ((hungarian_cost + 1) / 2);
                    if (lb > tau)
                        return hungarian_cost;
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
