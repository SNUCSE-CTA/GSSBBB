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
    std::string query_path, database_path;
    query_path = argv[1];
    LabeledGraph* query_graph = new LabeledGraph();
    query_graph->LoadFromFile(query_path);
    database_path = argv[2];
    int tau = std::stoi(argv[3]);
    std::vector<LabeledGraph*> database_graphs;
    std::ifstream fin(database_path);
    if (!fin.is_open()) {
        std::cerr << "Error opening file: " << database_path << std::endl;
        exit(1);
    }
    std::string bufferedLine;

    int num_answer = 0, num_total = 0;
    long long search_space = 0;
    double total_time = 0.0;
    while (true) {
        LabeledGraph* data_graph = new LabeledGraph();
        if (!data_graph->LoadFromStream(fin, bufferedLine))
            break;
        BinaryBranchingBestSwap solver;
        int found_ged = solver.GEDVerification(query_graph, data_graph, tau);
        if (found_ged <= tau) {
            num_answer++;
        }
        num_total++;
        auto log = solver.GetLog();
        total_time += log.GetNumericResult("ElapsedTime");
        search_space += (long long)log.GetNumericResult("SearchTreeNodes");
        delete data_graph;
    }
    fprintf(stderr, "[Result] ElapsedTime : %.04lf\n", total_time);
    fprintf(stderr, "[Result] Answer : %d\n", num_answer);
    fprintf(stderr, "[Result] SearchSpace : %d\n", search_space);
    delete query_graph;
    return 0;
}