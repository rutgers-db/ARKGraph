/**
 * Index Sorted Range KNN
 *
 * Author: Chaoji Zuo
 * Email:  chaoji.zuo@rutgers.edu
 * Date:   Oct 1, 2021
 */
#pragma once

#include <algorithm>
#include <iostream>
#include <numeric>
#include <queue>
#include <vector>

#include "baselines/greedy.h"
#include "hnswlib/hnswlib.h"
#include "utils/cw.h"
#include "utils/segmentTree.hpp"
#include "utils/utils.h"

using std::cout;
using std::endl;
using std::vector;

class BaseIndex {
 public:
  int window_count;
  int nodes_amount;
  double sort_time;
  double nn_build_time;
  int K1_;
  int K2_;
  bool isLog = true;
  virtual vector<int> indexNeighborSearch(
      const vector<vector<float>> &data_nodes, const vector<float> &query,
      const int l_bound, const int r_bound, const int K_neigbhors,
      const int ef) = 0;
  virtual vector<int> indexNeighborSearch(
      const vector<vector<float>> &data_nodes, const vector<float> &query,
      const int l_bound, const int r_bound, const int K_neigbhors, const int ef,
      const vector<int> &enter_points) = 0;
  virtual pair<vector<int>, vector<int>> checkIndexQuality(
      const int idx, const int lbound, const int rbound) = 0;
  virtual void buildIndex(const vector<vector<float>> &nodes,
                          const int build_knn_k, const int index_K) = 0;
  virtual ~BaseIndex() {}
};

#define SUB_NONE -1;
#define ADD_NONE -2;

struct DeltaNode {
  DeltaNode(int s, int a) : sub(s), add(a){};
  DeltaNode() {
    sub = SUB_NONE;
    add = ADD_NONE;
  };
  int sub;
  int add;
};