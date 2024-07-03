#pragma once
#include <list>

#include "Base/Logger.h"
#include "Base/Timer.h"
#include "DataStructure/LabeledGraph.h"
#include "DifferenceVector.h"

namespace GraphLib::GraphSimilarity {
    bool verbosity = true;
    static int32_t LOG_EVERY = 50000;

    class GraphEditDistanceSolver {
    protected:
        LabeledGraph *G1, *G2;
        int threshold = -1, current_best = INF;
        int64_t num_nodes = 0;
        std::vector<int> mapping, inverse_mapping;
        ResultLogger log;

    public:
        virtual ~GraphEditDistanceSolver() = default;

        int GEDVerification(LabeledGraph *G1_, LabeledGraph *G2_, int threshold_) {
            this->threshold = threshold_;
            // ensure V(G1) <= V(G2)
            if (G1_->GetNumVertices() > G2_->GetNumVertices()) {
                this->G1 = G2_;
                this->G2 = G1_;
            } else {
                this->G1 = G1_;
                this->G2 = G2_;
            }
            return Verify();
        }

        virtual int Verify() {
            return 0;
        }

        int GED(LabeledGraph *G1_, LabeledGraph *G2_) {
            log.clear();
            // ensure V(G1) <= V(G2)
            if (G1_->GetNumVertices() > G2_->GetNumVertices()) {
                this->G1 = G2_;
                this->G2 = G1_;
            } else {
                this->G1 = G1_;
                this->G2 = G2_;
            }
            return ComputeGED();
        }

        virtual int ComputeGED() {
            return 0;
        }

        ResultLogger GetLog() {
            return log;
        }

        void ComputeBranchDistanceMatrixInitial(std::vector<std::vector<int> > &matrix) {
            DifferenceVector diff(global_edge_labels.size() + 1);
            for (int v = 0; v < G2->GetNumVertices(); v++) {
                auto &v_nbrs = G2->GetNeighbors(v);
                int u = 0;
                for (u = 0; u < G1->GetNumVertices(); u++) {
                    auto &u_nbrs = G1->GetNeighbors(u);
                    diff.reset();
                    if (G1->GetVertexLabel(u) != G2->GetVertexLabel(v)) {
                        matrix[u][v] += 2;
                    }
                    for (int l = 0; l < u_nbrs.size(); l++) {
                        int u_nbr = u_nbrs[l];
                        int u_el = G1->GetEdgeLabel(u, u_nbr);
                        diff.update(u_el, 1);
                    }
                    for (int r = 0; r < v_nbrs.size(); r++) {
                        int v_nbr = v_nbrs[r];
                        int v_el = G2->GetEdgeLabel(v, v_nbr);
                        diff.update(v_el, -1);
                    }
                    int inner_distance = diff.GetDifference();
                    matrix[u][v] += inner_distance;
                }
                int from_null = 2 + G2->GetDegree(v);
                for (; u < G2->GetNumVertices(); u++) {
                    matrix[u][v] = from_null;
                }
            }
        }

        int ChildEditCost(int u, int v) {
            // int child_cost = cost;
            int child_cost = 0;
            int u_label = G1->GetVertexLabel(u), v_label = G2->GetVertexLabel(v);
            if (u_label != v_label) {
                child_cost++;
            }
            int num_u_edges = 0;
            for (int u_nbr : G1->GetNeighbors(u)) {
                if (mapping[u_nbr] != -1) {
                    num_u_edges++;
                }
            }
            int ec = num_u_edges;
            for (int vprime : G2->GetNeighbors(v)) {
                int uprime = inverse_mapping[vprime];
                if (uprime == -1)
                    continue;
                ec++;
                int l1 = G1->GetEdgeLabel(u, uprime);
                int l2 = G2->GetEdgeLabel(v, vprime);
                if (l1 == -1)
                    continue;
                if (l1 == l2)
                    ec -= 2;
                else
                    ec--;
            }
            child_cost += ec;
            return child_cost;
        }

        int ComputeDistance(std::vector<int> &mapping, bool verbose = false) {
            std::vector<int> inverse_mapping(G2->GetNumVertices(), -1);
            for (int i = 0; i < G1->GetNumVertices(); i++) {
                inverse_mapping[mapping[i]] = i;
            }
            return ComputeDistance(mapping, inverse_mapping, verbose);
        }

        int ComputeDistance(std::vector<int> &mapping, std::vector<int> &inverse_mapping, bool verbose = false) {
            int cost = 0;
            // vertex re-labeling cost
            for (int i = 0; i < G1->GetNumVertices(); i++) {
                int l1 = G1->GetVertexLabel(i);
                int l2 = G2->GetVertexLabel(mapping[i]);
                if (l1 != l2) {
                    if (verbose)
                        printf("Vertex %d(%d)-%d(%d) re-labeling cost!\n", i, l1, mapping[i], l2);
                    cost++;
                }
            }
            if (verbose)
                printf("#Missing Vertices: %d\n", (G2->GetNumVertices() - G1->GetNumVertices()));
            cost += (G2->GetNumVertices() - G1->GetNumVertices());
            for (auto &[u, v] : G1->GetEdges()) {
                int fu = mapping[u], fv = mapping[v];
                int l1 = G1->GetEdgeLabel(u, v);
                int l2 = G2->GetEdgeLabel(fu, fv);
                if (l1 != l2) {
                    if (verbose)
                        printf("Edge (%d, %d)(%d)-(%d, %d)(%d) re-labeling cost!\n", u, v, l1, fu, fv, l2);
                    cost++;
                }
            }
            for (auto &[u, v] : G2->GetEdges()) {
                int inv_u = inverse_mapping[u], inv_v = inverse_mapping[v];
                // if (inv_u == -1 || inv_v == -1) {
                if (inv_u == -1 || inv_v == -1 || inv_u >= G1->GetNumVertices() || inv_v >= G1->GetNumVertices()) {
                    if (verbose)
                        printf("Edge (%d, %d) in G2 is nonexistent as G1 is (%d, %d)\n", u, v, inv_u, inv_v);
                    cost++;
                } else {
                    // printf("%d %d\n", inv_u, inv_v);
                    // printf("%d\n", G1->GetNumVertices());
                    int l = G1->GetEdgeLabel(inv_u, inv_v);
                    if (l == -1) {
                        if (verbose)
                            printf("Edge (%d, %d) in G2 is nonexistent as G1 is (%d, %d)\n", u, v, inv_u, inv_v);
                        cost++;
                    }
                }
            }
            if (verbose)
                printf("Total ED Cost: %d\n", cost);
            current_best = std::min(cost, current_best);
            return cost;
        }

        std::vector<std::vector<int> > parikh1, parikh2;

        void ComputeParikhVector() {
            const int n1 = G1->GetNumVertices();
            const int n2 = G2->GetNumVertices();
            const int np = global_edge_labels.size();
            for (int u = 0; u < n1; ++u) {
                for (const int _u : G1->GetNeighbors(u))
                    if (mapping[_u] == -1) {
                        const int u_u = G1->GetEdgeLabel(u, _u);
                        ++parikh1[u][u_u];
                        ++parikh1[u][0];
                    }
            }
            for (int v = 0; v < n2; ++v) {
                for (const int _v : G2->GetNeighbors(v))
                    if (inverse_mapping[_v] == -1) {
                        const int v_v = G2->GetEdgeLabel(v, _v);
                        ++parikh2[v][v_v];
                        ++parikh2[v][0];
                    }
            }
        }

        void UpdateParikhVector(const int u, const int v) {
            auto &u_nbrs = G1->GetNeighbors(u);
            auto &v_nbrs = G2->GetNeighbors(v);
            for (auto _u : u_nbrs) {
                if (mapping[_u] == -1) {
                    const int u_u = G1->GetEdgeLabel(u, _u);
                    parikh1[_u][u_u]--;
                    parikh1[_u][0]--;
                }
            }
            for (auto _v : v_nbrs) {
                if (inverse_mapping[_v] == -1) {
                    const int v_v = G2->GetEdgeLabel(v, _v);
                    parikh2[_v][v_v]--;
                    parikh2[_v][0]--;
                }
            }
        }

        void RestoreParikhVector(const int u, const int v) {
            auto &u_nbrs = G1->GetNeighbors(u);
            auto &v_nbrs = G2->GetNeighbors(v);
            for (auto _u : u_nbrs) {
                if (mapping[_u] == -1) {
                    const int u_u = G1->GetEdgeLabel(u, _u);
                    parikh1[_u][u_u]++;
                    parikh1[_u][0]++;
                }
            }
            for (auto _v : v_nbrs) {
                if (inverse_mapping[_v] == -1) {
                    const int v_v = G2->GetEdgeLabel(v, _v);
                    parikh2[_v][v_v]++;
                    parikh2[_v][0]++;
                }
            }
        }
    };
}  // namespace GraphLib::GraphSimilarity
// note: research log marker 22
