
/**
 * Groundtruth: Greedy search among the range.
 *
 * Author: Chaoji Zuo
 * Date:   Nov 13, 2021
 * Email:  chaoji.zuo@rutgers.edu
 */

#include <algorithm>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

#include "baselines/greedy.h"
#include "baselines/knn_first_baseline.h"
#include "kgraph/kgraph-data.h"
#include "kgraph/kgraph.h"
#include "myHNSW/hnswlib-incrementally/hnswlib.h"
#include "range_index_base.h"

#ifdef __linux__
#include "sys/sysinfo.h"
#include "sys/types.h"
#endif

using std::cout;
using std::endl;
using std::string;
using std::vector;

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

void SaveToIndexCSVRow(const string &path, const string &method,
                       const int data_size, const int initial_graph_size,
                       const int index_graph_size, const double nn_build_time,
                       const double sort_time, const double build_time,
                       const double memory, const int node_amount,
                       const int window_count, const double index_size) {
  ofstream file;
  file.open(path, std::ios_base::app);
  if (file) {
    file << method << "," << data_size << "," << initial_graph_size << ","
         << index_graph_size << "," << nn_build_time << "," << sort_time << ","
         << build_time << "," << memory << "," << node_amount << ","
         << window_count << "," << index_size;
    file << "\n";
  }
  file.close();
}

int window_count = 0;
timeval tt1, tt2;
double tree_search_time, neighbor_search_time;
long long before_memory, after_memory;
int query_num = 1000;

string target_csv_path;
string target_index_csv_path;

void execute_one_query(BaseIndex &method, const vector<vector<float>> &data,
                       const vector<vector<float>> &querys,
                       const vector<vector<int>> &groundtruths,
                       const vector<pair<int, int>> &ranges,
                       const string method_name, const int K_query,
                       const int search_ef) {
  for (int idx = 0; idx < querys.size(); idx++) {
    timeval temp1, temp2;
    gettimeofday(&temp1, NULL);
    auto res = method.indexNeighborSearch(
        data, querys.at(idx), ranges.at(idx).first, ranges.at(idx).second,
        K_query, search_ef, vector<int>());
    gettimeofday(&temp2, NULL);
    double query_time = 0;
    CountTime(temp1, temp2, query_time);
    auto precision = countPrecision(groundtruths.at(idx), res);
    auto appr_ratio = countApproximationRatio(data, groundtruths.at(idx), res,
                                              querys.at(idx));
    SaveToCSVRow(target_csv_path, idx % query_num, ranges.at(idx).first,
                 ranges.at(idx).second,
                 ranges.at(idx).second - ranges.at(idx).first, K_query,
                 method.K1_, method.K2_, method_name, search_ef, precision,
                 appr_ratio, query_time);
  }
}

// For GroundTruth
void SaveGroundTruth(const string &path, vector<int> &gt) {
  ofstream file;
  file.open(path, std::ios_base::app);
  for (auto ele : gt) {
    file << ele << " ";
  }
  file << "\n";

  file.close();
}

void LoadGroundTruth(vector<vector<int>> &gt, string gt_path) {
  ifstream infile;
  string bline;
  string delim = " ";

  int numCols = 0;
  infile.open(gt_path, ios::in);
  int counter = 0;
  while (getline(infile, bline, '\n')) {
    counter++;
    vector<int> one_gt;
    std::pair<int, int> one_range;
    int one_id;
    vector<string> ret;
    Split(bline, delim, &ret);
    for (auto ele : ret) {
      one_gt.emplace_back(std::stoi(ele));
    }
    gt.emplace_back(one_gt);
  }
  cout << "# of GroundTruth: " << gt.size() << endl;
}

int exp1(string dataset, int data_size) {
  cout << "Start exp 1" << endl;
  print_memory();

  target_csv_path =
      "../exp_result/exp_index_storage_methods-query-" + dataset + ".csv";

  string path = "";

  vector<vector<float>> nodes;
  vector<int> search_keys;
  ReadDataWrapper(nodes, search_keys, dataset, path, data_size);
  cout << nodes.size() << endl;
  cout << "Load vecs from: " << path << endl;
  cout << "# of vecs: " << nodes.size() << endl;

  default_random_engine e;
  uniform_int_distribution<int> u(0, nodes.size());

  // calculate groundtruth.

  vector<unsigned> ranges = {2500, 5000, 7500, 10000};
  vector<int> index_k_paras = {8, 16, 32, 64};
  for (auto one_range : ranges) {
    vector<vector<float>> sub_nodes(nodes);
    sub_nodes.resize(one_range);
    vector<vector<int>> groundtruth;
    groundtruth.resize(one_range);
#pragma omp parallel for
    for (unsigned n = 0; n < one_range; n++) {
      groundtruth[n] = greedyNearest(sub_nodes, sub_nodes[n], 100);
    }

    string gt_path = "../groundtruth/knng100-10k-" + to_string(one_range) +
                     "-" + dataset + ".txt";
    for (unsigned n = 0; n < one_range; n++) {
      SaveGroundTruth(gt_path, groundtruth[n]);
    }

    // evalualte baseline recall
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

    for (auto para : index_k_paras) {
      cout << endl << "Range: " << one_range << endl;
      cout << "index_k: " << para << endl;
      gettimeofday(&tt1, NULL);
      kgraph::KGraph::IndexInfo info;
      kgraph::KGraph *kgraph_index = kgraph::KGraph::create();
      kgraph::KGraph::IndexParams params;
      params.K = para;
      kgraph_index->build(oracle, params, &info);

      vector<vector<int>> kgraph_vector;
      kgraph_vector.resize(sub_nodes.size());
      for (unsigned i = 0; i < sub_nodes.size(); ++i) {
        vector<unsigned> one_nns(130);
        unsigned pM, pL;
        kgraph_index->get_nn(i, &one_nns[0], &pM, &pL);
        for (unsigned j = 0; j < para; j++) {
          int point = (int)one_nns.at(j);
          kgraph_vector.at(i).emplace_back((int)one_nns.at(j));
        }
      }

      double build_time;
      gettimeofday(&tt2, NULL);
      CountTime(tt1, tt2, build_time);
      logTime(tt1, tt2, "Build KGraph Time");

      double recall, precision;
      evaluateKNNG(groundtruth, kgraph_vector, para, recall, precision);
      cout << "Recall: " << recall << endl;
      cout << "Precision: " << precision << endl;
    }
  }

  return 0;
}

// run exp wrapper
int main(int argc, char **argv) {
  string dataset = "biggraph";
  int data_size = 100000;

  for (int i = 0; i < argc; i++) {
    string arg = argv[i];
    if (arg == "-dataset") dataset = string(argv[i + 1]);
    if (arg == "-N") data_size = atoi(argv[i + 1]);
  }

  exp1(dataset, data_size);
  return 0;
}