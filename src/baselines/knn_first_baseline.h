/**
 * baseline #1, calculate nearest range query first, then search among the
 * range.
 *
 * Author: Chaoji Zuo
 * Date:   Nov 13, 2021
 * Email:  chaoji.zuo@rutgers.edu
 */
#pragma once

#include <algorithm>
#include <iostream>
#include <numeric>
#include <vector>

#include "myHNSW/hnswlib-incrementally/hnswlib.h"
#include "../utils/utils.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;

void buildKNNFirstGraph(const vector<vector<float>> &nodes,
                        hnswlib_incre::HierarchicalNSW<float> &alg_hnsw) {
#pragma omp parallel for
  for (size_t i = 0; i < nodes.size(); ++i) {
    alg_hnsw.addPoint(nodes[i].data(), i);
  }
}

void addHNSWPointsSubgraph(const vector<vector<float>> &nodes,
                           hnswlib_incre::HierarchicalNSW<float> *alg_hnsw,
                           const int start, const int end) {
#pragma omp parallel for
  for (size_t i = start; i <= end; ++i) {
    alg_hnsw->addPoint(nodes[i].data(), i);
  }
}

void buildKNNFirstGraphSingleThread(
    const vector<vector<float>> &nodes,
    hnswlib_incre::HierarchicalNSW<float> &alg_hnsw) {
  for (size_t i = 0; i < nodes.size(); ++i) {
    alg_hnsw.addPoint(nodes[i].data(), i);
  }
}

vector<int> KNNFirstRangeSearch(
    const hnswlib_incre::HierarchicalNSW<float> &alg_hnsw,
    const vector<float> &query, const int l_bound, const int r_bound,
    const int query_k) {
  // K to run nndescent

  vector<int> result_in_range;
  auto res =
      alg_hnsw.searchKnnCloserFirst(query.data(), query_k, l_bound, r_bound);
  for (size_t j = 0; j < res.size(); j++) {
    int val = res[j].second;
    result_in_range.emplace_back(val);
  }
  return result_in_range;
}

vector<int> KNNFirstRangeSearchFixedEF(
    const hnswlib_incre::HierarchicalNSW<float> &alg_hnsw,
    const vector<float> &query, const int l_bound, const int r_bound,
    const int query_k) {
  // K to run nndescent

  vector<int> result_in_range;
  auto res = alg_hnsw.searchKnnCloserFirst(query.data(), query_k, l_bound,
                                           r_bound, true);
  for (size_t j = 0; j < res.size(); j++) {
    int val = res[j].second;
    result_in_range.emplace_back(val);
  }
  return result_in_range;
}

vector<int> KNNFirstRangeSearchFixedEF(
    hnswlib_incre::HierarchicalNSW<float> *alg_hnsw, const vector<float> &query,
    const int l_bound, const int r_bound, const int query_k,
    const int search_ef) {
  // K to run nndescent

  alg_hnsw->setEf(search_ef);

  vector<int> result_in_range;
  auto res = alg_hnsw->searchKnnCloserFirst(query.data(), query_k, l_bound,
                                           r_bound, true);
  for (size_t j = 0; j < res.size(); j++) {
    int val = res[j].second;
    result_in_range.emplace_back(val);
  }
  return result_in_range;
}