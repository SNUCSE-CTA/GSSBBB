#include <iostream>

#include "Base/Timer.h"
#include "DataStructure/Graph.h"
#include "DataStructure/LabeledGraph.h"
#include "GED/AStarBinaryBranching.h"
#include "GED/BinaryBranching.h"
#include "GED/DFSDH.h"
#include "GED/EditDistance.h"
using namespace std;
using namespace GraphLib;
using namespace GraphLib::GraphSimilarity;
Timer timer;

std::vector<std::string> log_entries = {
    "EditDistance", "AStarNodes", "MaxQueueSize", "InitializeTime", "AStarTime"};
int main(int argc, char* argv[]) {
    std::string g1_path, g2_path;
    g1_path = argv[1];
    g2_path = argv[2];
    LabeledGraph* G1 = new LabeledGraph();
    G1->LoadFromFile(g1_path);
    LabeledGraph* G2 = new LabeledGraph();
    G2->LoadFromFile(g2_path);
    // AStarBinaryBranching solver;
    // DFSDH solver;
    BinaryBranchingBestSwap solver;
    // BinaryBranchingHighestWeightGap solver;
    solver.GED(G1, G2);
    auto log = solver.GetLog();
    log.PrintResults();
    return 0;
}
// note: research log marker 49
