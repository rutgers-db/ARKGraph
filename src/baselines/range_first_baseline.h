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

#include "../utils/utils.h"

using std::cout;
using std::endl;
using std::vector;

// vector<int> rangeFirstRangeSearch(const vector<vector<float>> &nodes,
//                                   const vector<float> &query, const int
//                                   l_bound, const int r_bound, const int K3) {
//   auto nodes_in_range =
//       vector<vector<float>>(nodes.begin() + l_bound, nodes.begin() +
//       r_bound);
//   // The number of subspace
//   int m = 16;
//   // The number of codewords for each subspace
//   int ks = 256;

//   // Initialize PQ
//   pqdescent::PQ pq(pqdescent::PQ::Learn(nodes_in_range, m, ks));
//   pqdescent::UcharVecs encodedNodes = pq.Encode(nodes_in_range);
//   auto codewords = pq.GetCodewords();
//   auto res = pq.AD(pq.DTable(query), encodedNodes);

//   auto idx_sorted = sort_indexes(res);
//   // skip element itself
//   idx_sorted.resize(K3);
//   for (auto &idx : idx_sorted) {
//     idx += l_bound;
//   }
//   return idx_sorted;
// }

// class PQMethod {
// private:
// public:
//   void init(const vector<vector<float>> nodes);
// };

// void PQMethod::init(const vector<vector<float>> nodes){
//   // The number of subspace
//   int m = 16;
//   // The number of codewords for each subspace
//   int ks = 256;

//   // Initialize PQ
//   pqdescent::PQ pq(pqdescent::PQ::Learn(nodes, m, ks));
//   pqdescent::UcharVecs encodedNodes = pq.Encode(nodes);
// }

vector<int> rangeFirstRangeSearch(pqdescent::PQ &pq,
                                  const pqdescent::UcharVecs &encodeNodes,
                                  const vector<float> &query, const int l_bound,
                                  const int r_bound, const int K3) {
  // naive
  // vector<vector<uchar>> some_nodes;
  // for (size_t i = l_bound; i <= r_bound; i++) {
  //   some_nodes.emplace_back(encodeNodes.GetVec(i));
  // }

  // update version
  // vector<uchar> some_nodes = encodeNodes.GetMVec1D(l_bound, r_bound -
  // l_bound);
  // pqdescent::UcharVecs code_in_range(r_bound - l_bound, encodeNodes.Dim(),
  //                                    some_nodes);

  auto res = pq.AD(pq.DTable(query), encodeNodes, l_bound, r_bound);
  auto idx_sorted = sort_indexes(res);
  // skip element itself
  idx_sorted.resize(K3);
  for (auto &idx : idx_sorted) {
    idx += l_bound;
  }
  return idx_sorted;
}