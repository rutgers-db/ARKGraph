/**
 * @file exp.cc
 * @author Chaoji Zuo (chaoji.zuo@rutgers.edu)
 * @brief Example running code to run ARKGraph
 * @date 2023-04-27
 *
 * @copyright Copyright (c) 2023
 */

#include <algorithm>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

#include "../src/baselines/knn_first_baseline.h"
#include "CompactKNNG/arkgraph.h"
#include "CompactKNNG/arkgraph_delta.h"
#include "CompactKNNG/arkgraph_delta_heap_version.h"
#include "CompactKNNG/arkgraph_delta_raw.h"
#include "CompactKNNG/arkgraph_delta_st.h"
#include "CompactKNNG/compact_graph_index.h"
#include "myHNSW/hnswlib-incrementally/hnswlib.h"

void SaveToCSVRow(const string &path, const int idx, const int l_bound,
                  const int r_bound, const int range, const int K_neighbor,
                  const int initial_graph_size, const int index_graph_size,
                  const string &method, const int search_ef,
                  const double &precision, const double &appr_ratio,
                  const double &search_time) {
  ofstream file;
  file.open(path, std::ios_base::app);
  if (file) {
    file << idx << "," << l_bound << "," << r_bound << "," << range << ","
         << K_neighbor << "," << initial_graph_size << "," << index_graph_size
         << "," << method << "," << search_ef << "," << precision << ","
         << appr_ratio << "," << search_time;
    file << "\n";
  }
  file.close();
}

void keepTopKNN(const vector<vector<float>> &nodes, vector<vector<int>> &knng,
                const int K) {
#pragma omp parallel for
  for (unsigned n = 0; n < knng.size(); n++) {
    auto &knn = knng[n];
    vector<int> new_knn;
    priority_queue<pair<float, int>> candidates;
    for (auto ele : knn) {
      candidates.emplace(-EuclideanDistance(nodes[n], nodes[ele]), ele);
    }
    while (new_knn.size() < K && !candidates.empty()) {
      new_knn.emplace_back(candidates.top().second);
      candidates.pop();
    }
    knn.swap(new_knn);
  }
}

int main(int argc, char **argv) {
  string dataset = "local";
  int data_size = 100000;
  string path = "";
  int index_k = 16;
  timeval t1, t2;

  for (int i = 0; i < argc; i++) {
    string arg = argv[i];
    if (arg == "-dataset") dataset = string(argv[i + 1]);
    if (arg == "-N") data_size = atoi(argv[i + 1]);
    if (arg == "-data_path") path = string(argv[i + 1]);
    if (arg == "-K") index_k = atoi(argv[i + 1]);
  }

  int query_num = 10;

  vector<vector<float>> nodes;
  vector<int> search_keys;
  ReadDataWrapper(nodes, search_keys, dataset, path, data_size);

  unsigned range = nodes.size() / 2;

  cout << "Load vecs from: " << path << endl;
  cout << "# of vecs: " << nodes.size() << endl;
  cout << "Range of KNNG: " << range << endl;

  // groundtruth
  vector<vector<float>> sub_nodes(nodes);
  sub_nodes.resize(range);
  vector<vector<int>> groundtruth;
  groundtruth.resize(range);
  for (unsigned n = 0; n < range; n++) {
    groundtruth[n] = greedyNearest(sub_nodes, sub_nodes[n], 100);
  }

  // baseline: nndescent (kgraph)
  {
    cout << endl << "Baseline: KGraph" << endl;
    auto distance = [](vector<float> a, vector<float> b) -> double {
      float ans = 0.0;
      for (int i = 0; i < a.size(); ++i) {
        ans += (a[i] - b[i]) * (a[i] - b[i]);
      }
      return ans;
    };
    typedef kgraph::VectorOracle<vector<vector<float>>, vector<float>> MyOracle;
    MyOracle oracle(sub_nodes, distance);
    kgraph::verbosity = 0;
    gettimeofday(&t1, NULL);
    kgraph::KGraph::IndexInfo info;
    kgraph::KGraph *kgraph_index = kgraph::KGraph::create();
    kgraph::KGraph::IndexParams params;
    params.K = index_k;
    kgraph_index->build(oracle, params, &info);
    vector<vector<int>> kgraph_vector;
    kgraph_vector.resize(sub_nodes.size());
    for (unsigned i = 0; i < sub_nodes.size(); ++i) {
      vector<unsigned> one_nns(130);
      unsigned pM, pL;
      kgraph_index->get_nn(i, &one_nns[0], &pM, &pL);
      for (unsigned j = 0; j < index_k; j++) {
        int point = (int)one_nns.at(j);
        kgraph_vector.at(i).emplace_back((int)one_nns.at(j));
      }
    }
    logTime(t1, t2, "Build KGraph Time");
    double recall, precision;
    evaluateKNNG(groundtruth, kgraph_vector, index_k, recall, precision);
    cout << "KGraph Recall: " << recall << endl;
    cout << "KGraph Precision: " << precision << endl;
  }

  // ARKGraph method
  {
    cout << endl << "ARKGraph:" << endl;
    gettimeofday(&t1, NULL);
    PartialIndex partial_index;
    partial_index.buildIndex(nodes, index_k, index_k);
    logTime(t1, t2, "Build ARKGraph Index Time");
    gettimeofday(&t1, NULL);
    vector<vector<int>> sub_knng;
    double recall, precision;
    partial_index.rebuildRangeKNNG(0, range - 1, sub_knng);
    keepTopKNN(nodes, sub_knng, index_k);
    logTime(t1, t2, "Rebuild Graph Time");
    evaluateKNNG(groundtruth, sub_knng, index_k, recall, precision);
    cout << "ARKGraph Recall: " << recall << endl;
    cout << "ARKGraph Precision: " << precision << endl;
  }

  // delta ARKGraph method
  {
    cout << endl << "Delta ARKGraph:" << endl;
    gettimeofday(&t1, NULL);
    deltaindex::DeltaIndex delta_index;
    delta_index.buildIndex(nodes, index_k, index_k);
    logTime(t1, t2, "Build Delta ARKGraph Index Time");
    gettimeofday(&t1, NULL);
    vector<vector<int>> sub_knng;
    double recall, precision;
    delta_index.rebuildRangeKNNG(0, range - 1, sub_knng);
    keepTopKNN(nodes, sub_knng, index_k);
    logTime(t1, t2, "Rebuild Graph Time");
    evaluateKNNG(groundtruth, sub_knng, index_k, recall, precision);
    cout << "Delta ARKGraph Recall: " << recall << endl;
    cout << "Delta ARKGraph Precision: " << precision << endl;
  }

  // delta ARKGraph Raw method
  {
    cout << endl << "Delta ARKGraph Raw:" << endl;
    gettimeofday(&t1, NULL);
    deltaindex::DeltaRawIndex delta_index;
    delta_index.buildIndex(nodes, index_k, index_k);
    logTime(t1, t2, "Build Delta ARKGraph Raw Index Time");
    gettimeofday(&t1, NULL);
    vector<vector<int>> sub_knng;
    double recall, precision;
    delta_index.rebuildRangeKNNG(0, range - 1, sub_knng);
    keepTopKNN(nodes, sub_knng, index_k);
    logTime(t1, t2, "Rebuild Graph Time");
    evaluateKNNG(groundtruth, sub_knng, index_k, recall, precision);
    cout << "Delta ARKGraph Raw Recall: " << recall << endl;
    cout << "Delta ARKGraph Raw Precision: " << precision << endl;
  }
}