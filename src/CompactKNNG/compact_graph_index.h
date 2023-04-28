/**
 * CAL
 * Compact Graph Index
 * Compact Adjacent List
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

#include "../NearestSearch_kgraph.h"
#include "../baselines/greedy.h"
#include "../range_index_base.h"
#include "../utils/cw.h"
#include "../utils/segmentTree.hpp"
#include "../utils/utils.h"
#include "kgraph/kgraph-data.h"
#include "kgraph/kgraph.h"

using std::cout;
using std::endl;
using std::vector;

namespace compact_index {

class CompactIndex : public BaseIndex {
 private:
  vector<segTree::TreeNode *> indexed_tree;

 public:
  CompactIndex() {
    window_count = 0;
    nodes_amount = 0;
    sort_time = 0;
  }

  void collectWindows(vector<CompactWindow *> &windows, vector<int> &l_arr,
                      vector<int> &r_arr, vector<int> &visited_left,
                      vector<int> &visited_right, const vector<int> &idx_sorted,
                      const int k_smallest, const int query_pos);
  void indexNodes(const int query_pos, const vector<vector<float>> &data_nodes,
                  const vector<int> &idx_sorted, const int k_smallest,
                  const int l_bound, const int r_bound, segTree::TreeNode **);
  void buildIndex(const vector<vector<float>> &nodes, const int, const int);
  vector<int> indexNeighborSearch(const vector<vector<float>> &data_nodes,
                                  const vector<float> &query, const int l_bound,
                                  const int r_bound, const int K_neigbhors,
                                  const int ef);
  vector<int> indexNeighborSearch(const vector<vector<float>> &data_nodes,
                                  const vector<float> &query, const int l_bound,
                                  const int r_bound, const int K_neigbhors,
                                  const int ef,
                                  const vector<int> &enter_points);
  pair<vector<int>, vector<int>> checkIndexQuality(const int idx,
                                                   const int lbound,
                                                   const int rbound);
  double calSize();

  void clear() {
    vector<segTree::TreeNode *>().swap(indexed_tree);
    window_count = 0;
    indexed_tree.clear();
  }
};

void CompactIndex::collectWindows(vector<CompactWindow *> &windows,
                                  vector<int> &l_arr, vector<int> &r_arr,
                                  vector<int> &visited_left,
                                  vector<int> &visited_right,
                                  const vector<int> &idx_sorted,
                                  const int k_smallest, const int query_pos) {
  for (auto idx : idx_sorted) {
    // cout << l_arr.size() << " and " << r_arr.size() << endl;
    // cout << idx << " : " << endl;
    if (l_arr.size() > (k_smallest) && idx < l_arr[1]) {
      continue;
    }
    if (r_arr.size() > (k_smallest) && idx > r_arr[r_arr.size() - 2]) {
      continue;
    }
    int insert_pos = 0;
    if (idx < query_pos) {
      bool sign = 0;
      for (auto i = l_arr.begin(); i != l_arr.end(); i++) {
        if (*i > idx) {
          l_arr.insert(i, idx);
          sign = 1;
          break;
        }
        insert_pos++;
      }
      if (!sign) {
        l_arr.emplace_back(idx);
      }
      if (l_arr.size() > k_smallest + 1) {
        l_arr.erase(l_arr.begin());
      }
      if (l_arr.size() + r_arr.size() >= k_smallest + 2) {
        // cout << "insert at " << insert_pos << endl;
        for (auto pos = 0; pos < insert_pos; pos++) {
          if (pos >= k_smallest) {
            break;
          }
          if (l_arr.size() + r_arr.size() - pos < k_smallest + 2) {
            break;
          }

#ifndef ONLY_COUNT_WINDOWS
          auto it = l_arr.begin() + pos + 1;
          CompactWindow *new_cw = new CompactWindow(k_smallest);
          new_cw->lbound = *(it - 1);
          auto r_it = r_arr.begin();
          int counter = 0;
          while (it + counter != l_arr.end()) {
            new_cw->pts[counter] = *(it + counter);
            counter++;
          }
          while (counter < k_smallest) {
            new_cw->pts[counter++] = *(r_it++);
          }
          new_cw->rbound = *(r_it);
          // new_cw->print();
          // if (query_pos == 42)
          // {
          //   cout << "insert among: " << new_cw->lbound << "," <<
          //   new_cw->pts[0] << "  ->  " << new_cw->pts.back() << "," <<
          //   new_cw->rbound << endl;
          //   // new_cw->print();
          // }
          // if (new_cw->pts[0] > query_pos)
          // {
          //   new_cw->pts[0] = query_pos;
          // }
          // if (new_cw->pts.back() < query_pos)
          // {
          //   new_cw->pts.back() = query_pos;
          // }
          windows.emplace_back(new_cw);
          visited_left.emplace_back(new_cw->lbound);
          visited_left.emplace_back(new_cw->pts[0]);
          visited_right.emplace_back(new_cw->pts.back());
          visited_right.emplace_back(new_cw->rbound);
          // segTree::insertSegment(l_bound, r_bound, new_cw->lbound,
          // new_cw->pts[0], root, new_cw, true, new_cw->pts.back(),
          // new_cw->rbound);

#endif
        }
      }
    } else if (idx > query_pos) {
      bool sign = 0;
      // cout << idx << endl;
      for (auto i = r_arr.begin(); i != r_arr.end(); i++) {
        // cout << insert_pos << endl;
        if (*i > idx) {
          r_arr.insert(i, idx);
          sign = 1;
          break;
        }
        insert_pos++;
      }
      if (!sign) {
        r_arr.insert(r_arr.end() - 1, idx);
      }
      if (r_arr.size() > k_smallest + 1) {
        r_arr.erase(r_arr.end() - 1);
      }
      if (l_arr.size() + r_arr.size() >= k_smallest + 2) {
        // cout << "insert at " << insert_pos << endl;
        // pos means at least #pos elements in right side
        for (int pos = insert_pos + 1; pos < k_smallest + 1; pos++) {
          int l_bias = k_smallest - pos;
          // cout << "insert at " << pos << endl;
          // cout << "pos: " << pos << " l_bias: " << l_bias << endl;

          if (l_bias >= l_arr.size()) {
            continue;
          }
          if (r_arr.size() + l_bias <= k_smallest) {
            break;
          }
          auto it = l_arr.end() - l_bias - 1;
          if (it == l_arr.end()) {
            break;
          }

#ifndef ONLY_COUNT_WINDOWS
          CompactWindow *new_cw = new CompactWindow(k_smallest);
          new_cw->lbound = *(it++);
          auto r_it = r_arr.begin();
          int counter = 0;
          while (it + counter != l_arr.end()) {
            new_cw->pts[counter] = *(it + counter);
            counter++;
          }
          while (counter < k_smallest) {
            new_cw->pts[counter++] = *(r_it++);
          }
          new_cw->rbound = *(r_it);
          // cout << idx << " : " << endl;
          // new_cw->print();
          // if (query_pos == 42)
          // cout << "insert among: " << new_cw->lbound << "," << new_cw->pts[0]
          // << "  ->  " << new_cw->pts.back() << "," << new_cw->rbound << endl;
          // if (new_cw->pts[0] > query_pos)
          // {
          //   new_cw->pts[0] = query_pos;
          // }
          // if (new_cw->pts.back() < query_pos)
          // {
          //   new_cw->pts.back() = query_pos;
          // }

          windows.emplace_back(new_cw);
          visited_left.emplace_back(new_cw->lbound);
          visited_left.emplace_back(new_cw->pts[0]);
          visited_right.emplace_back(new_cw->pts.back());
          visited_right.emplace_back(new_cw->rbound);
          // segTree::insertSegment(l_bound, r_bound, new_cw->lbound,
          // new_cw->pts[0], root, new_cw, true, new_cw->pts.back(),
          // new_cw->rbound);

#endif
        }
      }
    }
  }
}

// index sorted vector
// idx_sorted: index vector sorted by distance to query point
void CompactIndex::indexNodes(const int query_pos,
                              const vector<vector<float>> &data_nodes,
                              const vector<int> &idx_sorted,
                              const int k_smallest, const int l_bound,
                              const int r_bound, segTree::TreeNode **root) {
  vector<int> l_arr;
  vector<int> r_arr;

  // store visited indexes to generate tree.
  vector<int> visited_left;
  vector<int> visited_right;

  vector<CompactWindow *> windows;

  l_arr.emplace_back(-1);
  r_arr.emplace_back(r_bound + 1);

  collectWindows(windows, l_arr, r_arr, visited_left, visited_right, idx_sorted,
                 k_smallest, query_pos);

  // doesn't care about nearest negibhors
  // int l_position = 0;
  // int r_position = r_bound;
  // if (l_arr.size() >= 1)
  //   l_position = l_arr.front();
  // if (r_arr.size() >= 1)
  //   r_position = r_arr.back();
  // vector<float> dist_arr;
  // for (size_t i = l_position + 1; i < r_position; i++) {

  //   dist_arr.emplace_back(
  //       EuclideanDistance(data_nodes[query_pos], data_nodes[i], 0, 0));
  // }
  // vector<int> idx_re_sorted = sort_indexes(dist_arr);

  // vector<int> early_end_idx_re_sorted;
  // for (auto &idx : idx_re_sorted) {
  //   idx = idx + l_position + 1;
  //   early_end_idx_re_sorted.emplace_back(idx);
  // }

  // for (auto pos : l_arr) {
  //   early_end_idx_re_sorted.erase(std::remove(early_end_idx_re_sorted.begin(),
  //                                             early_end_idx_re_sorted.end(),
  //                                             pos),
  //                                 early_end_idx_re_sorted.end());
  // }
  // for (auto pos : r_arr) {
  //   early_end_idx_re_sorted.erase(std::remove(early_end_idx_re_sorted.begin(),
  //                                             early_end_idx_re_sorted.end(),
  //                                             pos),
  //                                 early_end_idx_re_sorted.end());
  // }

  // collectWindows(windows, l_arr, r_arr, visited_left, visited_right,
  //                early_end_idx_re_sorted, k_smallest, query_pos);

#ifndef ONLY_COUNT_WINDOWS
  // sort visited vector
  visited_left.emplace_back(-1);
  visited_right.emplace_back(-1);

  visited_left.emplace_back(r_bound + 1);
  visited_right.emplace_back(r_bound + 1);

  std::sort(visited_right.begin(), visited_right.end());
  visited_right.erase(std::unique(visited_right.begin(), visited_right.end()),
                      visited_right.end());
  std::sort(visited_left.begin(), visited_left.end());
  visited_left.erase(std::unique(visited_left.begin(), visited_left.end()),
                     visited_left.end());

  *root = segTree::createSegmentTree(visited_left);
  createSegmentTreeSecondDim(*root, visited_right);

  for (auto &cw : windows) {
    segTree::insertSegment(l_bound, r_bound, cw->lbound, cw->pts[0], *root, cw,
                           true, cw->pts.back(), cw->rbound);
  }
#endif
  window_count += windows.size();
}

void CompactIndex::buildIndex(const vector<vector<float>> &nodes,
                              const int initial_build_k, const int index_k) {
  auto distance = [](vector<float> a, vector<float> b) -> double {
    float ans = 0.0;
    for (int i = 0; i < a.size(); ++i) {
      ans += (a[i] - b[i]) * (a[i] - b[i]);
    }
    return ans;
  };
  indexed_tree.resize(nodes.size());
  K1_ = initial_build_k;
  K2_ = index_k;
  timeval tt1, tt2;
  gettimeofday(&tt1, NULL);
  // use kgraph to build knn graph
  kgraph::Matrix<float> data_mat;
  data_mat.resize(nodes.size(), nodes.front().size());
  for (unsigned i = 0; i < nodes.size(); ++i) {
    float *row = data_mat[i];
    for (unsigned j = 0; j < nodes.front().size(); ++j) {
      row[j] = nodes.at(i).at(j);
    }
  }
  float avg_neighbor;
  kgraph::MatrixOracle<float, kgraph::metric::l2sqr> oracle(data_mat);
  kgraph::verbosity = 0;
  kgraph::KGraph::IndexInfo info;
  kgraph::KGraph *kgraph_index = kgraph::KGraph::create();
  {
    kgraph::KGraph::IndexParams params;
    params.K = initial_build_k;
    kgraph_index->build(oracle, params, &info);
    avg_neighbor = info.M;
  }
  gettimeofday(&tt2, NULL);
  if (isLog) {
    logTime(tt1, tt2, "Compute Graph Time");
    cout << "# of average neighbors: " << avg_neighbor << endl;
  }
  int counter = 0;
  gettimeofday(&tt1, NULL);
  int hop_amount = 0;

#pragma omp parallel for
  for (int idx = 0; idx < nodes.size(); idx++) {
    timeval tt3, tt4;
    gettimeofday(&tt3, NULL);
    IncrementalNeighbors_Kgraph IN(idx, &nodes, kgraph_index, avg_neighbor);
    IN.findAllNeighbors(1, index_k);
    auto idx_sorted = IN.sorted_idxes;
    hop_amount += IN.getHop();
    gettimeofday(&tt4, NULL);
    AccumulateTime(tt3, tt4, sort_time);
    segTree::TreeNode *tree;
    indexNodes(idx, nodes, idx_sorted, index_k, 0, nodes.size() - 1, &tree);
    indexed_tree[idx] = tree;
    counter++;
  }
  nodes_amount = window_count;
  delete kgraph_index;

  if (isLog) {
    cout << "Sort Time: " << sort_time << endl;
    logTime(tt1, tt2, "Index Time");
    cout << "average hop: " << (float)hop_amount / nodes.size() << endl;
    cout << "index amount: " << window_count << endl;
  }
}

vector<int> CompactIndex::indexNeighborSearch(
    const vector<vector<float>> &data_nodes, const vector<float> &query,
    const int l_bound, const int r_bound, const int K_neigbhors, const int ef) {
  vector<int> temp;
  return indexNeighborSearch(data_nodes, query, l_bound, r_bound, K_neigbhors,
                             ef, temp);
}

// on the fly search cw
vector<int> CompactIndex::indexNeighborSearch(
    const vector<vector<float>> &data_nodes, const vector<float> &query,
    const int l_bound, const int r_bound, const int K_neigbhors, const int ef,
    const vector<int> &enter_points) {
  // assert(query.size() == data_nodes.front().size());
  unordered_map<int, bool> visited_list;
  float lower_bound = INT_MAX;
  priority_queue<pair<float, int>> top_candidates;
  priority_queue<pair<float, int>> candidate_set;
  timeval tt1, tt2;

  // multiple enter points: 10 enter points
  vector<int> enter_list;
  if (enter_points.size() == 0) {
    int interval = (r_bound - l_bound) / 11;
    for (size_t i = 1; i <= 10; i++) {
      int point = l_bound + interval * i;
      float dist = EuclideanDistance(data_nodes[point], query);
      candidate_set.push(make_pair(-dist, point));
      enter_list.emplace_back(point);
    }
  } else {
    for (auto point : enter_points) {
      float dist = EuclideanDistance(data_nodes[point], query);
      candidate_set.push(make_pair(-dist, point));
      enter_list.emplace_back(point);
    }
  }

  while (!candidate_set.empty()) {
    std::pair<float, int> current_node_pair = candidate_set.top();
    int current_node_id = current_node_pair.second;

    if (-current_node_pair.first > lower_bound) {
      break;
    }

    // cout << "current node: " << current_node_pair.second << "   ";
    candidate_set.pop();

    // only search when candidate point is inside the range
    if (current_node_id < l_bound || current_node_id > r_bound) {
      // cout << "no satisfied range point" << endl;
      continue;
    }

    // search cw on the fly
    // gettimeofday(&tt1, NULL);
    CompactWindow *cw =
        segTree::searchCW(l_bound, indexed_tree[current_node_id], r_bound);
    vector<int> current_neighbors;
    // gettimeofday(&tt2, NULL);
    // AccumulateTime(tt1, tt2, tree_search_time);

    // if cannot find in index, greedy search neigbhors
    if (cw == nullptr) {
      current_neighbors = greedyNearest(data_nodes, data_nodes[current_node_id],
                                        l_bound, r_bound, K_neigbhors);
    } else {
      current_neighbors = cw->pts;
    }

    if (current_neighbors.size() == 0) {
      cout << "current node: " << current_node_pair.second << "   ";
      cout << "size 0 neighbors" << endl;
      continue;
    }
    // print_set(current_neighbors);
    // gettimeofday(&tt1, NULL);
    for (size_t i = 0; i < current_neighbors.size(); i++) {
      int candidate_id = current_neighbors[i];
      if (!visited_list[candidate_id]) {
        visited_list[candidate_id] = true;
        float dist = EuclideanDistance(query, data_nodes[candidate_id]);
        // cout << "candidate: " << candidate_id << "  dist: " << dist << endl;
        if (top_candidates.size() < ef || lower_bound > dist) {
          candidate_set.push(make_pair(-dist, candidate_id));
          top_candidates.push(make_pair(dist, candidate_id));
          if (top_candidates.size() > ef) {
            top_candidates.pop();
          }
          if (!top_candidates.empty()) {
            lower_bound = top_candidates.top().first;
          }
        }
      }
    }
    // gettimeofday(&tt2, NULL);
    // AccumulateTime(tt1, tt2, neighbor_search_time);
  }
  vector<int> res;
  while (top_candidates.size() > K_neigbhors) {
    top_candidates.pop();
  }
  while (!top_candidates.empty()) {
    res.emplace_back(top_candidates.top().second);
    top_candidates.pop();
  }
  // cout << "tree search time: " << tree_search_time << endl;
  // cout << "neighbor search time: " << neighbor_search_time << endl;
  return res;
}

pair<vector<int>, vector<int>> CompactIndex::checkIndexQuality(
    const int idx, const int lbound, const int rbound) {
  vector<int> res = segTree::searchCW(lbound, indexed_tree[idx], rbound)->pts;
  print_set(res);
  return make_pair(res, res);
}

double CompactIndex::calSize() {
  double unit_size = sizeof(this->indexed_tree.front());
  return unit_size * window_count;
}

}  // namespace compact_index