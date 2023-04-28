/**
 * ARKGraph
 * Compact Adjacent List + Partial Ranges
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

#include "NearestSearch.h"
#include "NearestSearch_kgraph.h"
#include "baselines/greedy.h"
#include "baselines/knn_first_baseline.h"
#include "kgraph/kgraph-data.h"
#include "kgraph/kgraph.h"
#include "myHNSW/hnswlib-incrementally/hnswlib.h"
#include "myHNSW/hnswlib-nsw/hnswlib.h"
#include "range_index_base.h"
#include "utils/compactArr.hpp"
#include "utils/cw.h"
#include "utils/segmentTree.hpp"
#include "utils/utils.h"

using std::cout;
using std::endl;
using std::vector;

class PartialIndex : public BaseIndex {
 private:
  vector<pair<segTree::TreeNode *, segTree::TreeNode *>> tree_arr;

  vector<pair<compactArr::CompactArray, compactArr::CompactArray>> compact_arr;

 public:
  PartialIndex() {
    window_count = 0;
    nodes_amount = 0;
    sort_time = 0;
    hop_limit = 10;
  }
  int hop_limit;
  float avg_hop;

  void collectWindows(vector<CompactWindow *> &left_windows,
                      vector<CompactWindow *> &right_windows,
                      vector<int> &l_arr, vector<int> &r_arr,
                      vector<int> &visited_left, vector<int> &visited_right,
                      const vector<int> &idx_sorted, const int k_smallest,
                      const int query_pos);
  void indexNodes(const int query_pos, const vector<vector<float>> &data_nodes,
                  const vector<int> &idx_sorted, const int k_smallest,
                  const int l_bound, const int r_bound,
                  compactArr::CompactArray &left_arr,
                  compactArr::CompactArray &right_arr);
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
  double calSize();
  pair<vector<int>, vector<int>> checkIndexQuality(const int idx,
                                                   const int lbound,
                                                   const int rbound);

  void rebuildRangeKNNG(const int left, const int right,
                        vector<vector<int>> &knng);
  void clear() {
    for (auto &ele : compact_arr) {
      ele.first.clear();
      ele.second.clear();
    }
    vector<pair<compactArr::CompactArray, compactArr::CompactArray>>().swap(
        compact_arr);

    compact_arr.clear();
    window_count = 0;
    tree_arr.clear();
  }

  ~PartialIndex() {
    for (auto &ele : compact_arr) {
      ele.first.clear();
      ele.second.clear();
    }
    vector<pair<compactArr::CompactArray, compactArr::CompactArray>>().swap(
        compact_arr);
  }
};

void PartialIndex::collectWindows(vector<CompactWindow *> &left_windows,
                                vector<CompactWindow *> &right_windows,
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
      if (l_arr.size() >= k_smallest + 1) {
        // visited_left.emplace_back(idx);
        window_count += 1;
#ifndef ONLY_COUNT_WINDOWS
        CompactWindow *new_cw = new CompactWindow(k_smallest);
        new_cw->lbound = l_arr.front();
        for (size_t pos = 0; pos < k_smallest; pos++) {
          new_cw->pts[pos] = l_arr[pos + 1];
        }
        new_cw->rbound = query_pos;
        left_windows.emplace_back(new_cw);
        // visited_left.emplace_back(new_cw->pts.front());
        // visited_left.emplace_back(new_cw->lbound);
#endif
        // print_set(l_arr);
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
      if (r_arr.size() >= k_smallest + 1) {
        // visited_right.emplace_back(idx);
        window_count += 1;
#ifndef ONLY_COUNT_WINDOWS
        CompactWindow *new_cw = new CompactWindow(k_smallest);
        new_cw->lbound = query_pos;
        for (size_t pos = 0; pos < k_smallest; pos++) {
          new_cw->pts[pos] = r_arr[pos];
        }
        new_cw->rbound = r_arr.back();
        right_windows.emplace_back(new_cw);
        // visited_right.emplace_back(new_cw->pts.back());
        // visited_right.emplace_back(new_cw->rbound);
#endif
        // print_set(r_arr);
      }
    }
  }
}

// index sorted vector
// idx_sorted: index vector sorted by distance to query point
void PartialIndex::indexNodes(const int query_pos,
                            const vector<vector<float>> &data_nodes,
                            const vector<int> &idx_sorted, const int k_smallest,
                            const int l_bound, const int r_bound,
                            compactArr::CompactArray &left_arr,
                            compactArr::CompactArray &right_arr) {
  vector<int> l_arr;
  vector<int> r_arr;

  // store visited indexes to generate tree.
  vector<int> visited_left;
  vector<int> visited_right;

  vector<CompactWindow *> left_windows;
  vector<CompactWindow *> right_windows;

  l_arr.emplace_back(-1);
  r_arr.emplace_back(r_bound + 1);

  collectWindows(left_windows, right_windows, l_arr, r_arr, visited_left,
                 visited_right, idx_sorted, k_smallest, query_pos);

  // add dummy for not enough case
  if (left_windows.size() == 0) {
    CompactWindow *new_cw = new CompactWindow(l_arr.size() - 1);
    new_cw->lbound = l_arr.front();
    size_t pos = 0;
    for (; pos < l_arr.size() - 1; pos++) {
      new_cw->pts[pos] = l_arr[pos + 1];
    }
    // for (; pos < k_smallest; pos++) {
    //   new_cw->pts[pos] = query_pos;
    // }
    new_cw->rbound = query_pos;
    left_windows.emplace_back(new_cw);
  }

  if (right_windows.size() == 0) {
    CompactWindow *new_cw = new CompactWindow(r_arr.size() - 1);
    new_cw->lbound = query_pos;
    size_t pos = 0;
    for (; pos < r_arr.size() - 1; pos++) {
      new_cw->pts[pos] = r_arr[pos];
    }
    // new_cw->pts.resize(r_arr.size());
    new_cw->rbound = r_arr.back();
    right_windows.emplace_back(new_cw);
  }

  for (auto &cw : left_windows) {
    left_arr.insertCW(cw->lbound, cw);
  }

  reverse(right_windows.begin(), right_windows.end());

  for (auto &cw : right_windows) {
    right_arr.insertCW(cw->rbound, cw);
  }
  nodes_amount += left_windows.size();
  nodes_amount += right_windows.size();

  return;
}

void PartialIndex::buildIndex(const vector<vector<float>> &nodes,
                            const int initial_build_k, const int index_k) {
  auto distance = [](vector<float> a, vector<float> b) -> double {
    float ans = 0.0;
    for (int i = 0; i < a.size(); ++i) {
      ans += (a[i] - b[i]) * (a[i] - b[i]);
    }
    return ans;
  };
  compact_arr.clear();
  compact_arr.resize(nodes.size());
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
  double sort_time_ = 0;

#pragma omp parallel for reduction(+ : sort_time_)
  for (int idx = 0; idx < nodes.size(); idx++) {
    timeval tt3, tt4;
    gettimeofday(&tt3, NULL);
    IncrementalNeighbors_Kgraph IN(idx, &nodes, kgraph_index, avg_neighbor);
    IN.findAllNeighbors(hop_limit, index_k);
    auto idx_sorted = IN.sorted_idxes;
    hop_amount += IN.getHop();

    gettimeofday(&tt4, NULL);
    double one_sort_time;
    CountTime(tt3, tt4, one_sort_time);
    sort_time_ += one_sort_time;

    compactArr::CompactArray left_arr;
    compactArr::CompactArray right_arr;

    indexNodes(idx, nodes, idx_sorted, index_k, 0, nodes.size() - 1, left_arr,
               right_arr);
    compact_arr[idx].first = left_arr;
    compact_arr[idx].second = right_arr;
    counter++;
  }
  delete kgraph_index;
  sort_time = sort_time_;
  cout << "Sort Time: " << sort_time << endl;
  logTime(tt1, tt2, "Index Time");
  cout << "average hop: " << (float)hop_amount / nodes.size() << endl;
  cout << "index amount: " << window_count << endl;

  avg_hop = (float)hop_amount / nodes.size();
}

vector<int> PartialIndex::indexNeighborSearch(
    const vector<vector<float>> &data_nodes, const vector<float> &query,
    const int l_bound, const int r_bound, const int K_neigbhors, const int ef) {
  vector<int> temp;
  return indexNeighborSearch(data_nodes, query, l_bound, r_bound, K_neigbhors,
                             ef, temp);
}

// TODO: add a parameter like ef to control the neighbor num in searching
// on the fly search cw
vector<int> PartialIndex::indexNeighborSearch(
    const vector<vector<float>> &data_nodes, const vector<float> &query,
    const int l_bound, const int r_bound, const int K_neigbhors, const int ef,
    const vector<int> &enter_points) {
  if (compact_arr.empty()) {
    cout << "ERROR: haven't build index yet!" << endl;
    return vector<int>();
  }
  // assert(query.size() == data_nodes.front().size());
  unordered_map<int, bool> visited_list;
  float lower_bound = INT_MAX;
  priority_queue<pair<float, int>> top_candidates;
  priority_queue<pair<float, int>> candidate_set;

  // single enter point
  // lower_bound = EuclideanDistance(data_nodes[enter_point], query);
  // candidate_set.push(make_pair(-lower_bound, enter_point));

  // multiple enter points: 5 enter points
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
      // if (find(enter_list.begin(), enter_list.end(), current_node_id) ==
      // enter_list.end())
      // {
      // cout << "break in advance" << endl;
      break;
      // }
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
    CompactWindow *left_cw =
        compact_arr[current_node_id].first.searchCW(l_bound, true);
    CompactWindow *right_cw =
        compact_arr[current_node_id].second.searchCW(r_bound);

    vector<int> current_neighbors;
    // gettimeofday(&tt2, NULL);
    // AccumulateTime(tt1, tt2, tree_search_time);
    // if (cw == nullptr) {
    //   cout << current_node_id << " no cw found in (" << l_bound << ","
    //        << r_bound << ")" << endl;
    //   // segTree::print(tree_arr[current_node_id]);
    //   // assert(false);
    //   continue;
    // }
    // if cannot find in index, greedy search neigbhors
    if (left_cw == nullptr && right_cw == nullptr) {
      // current_neighbors = greedyNearest(data_nodes,
      // data_nodes[current_node_id],
      //                                   l_bound, r_bound, 10);

      // cout << current_node_id << " no cw found in (" << l_bound << ","
      //      << r_bound << ")" << endl;
      // segTree::print(tree_arr[current_node_id].first);
      // segTree::print(tree_arr[current_node_id].second);
      // assert(false);
      continue;

    } else {
      if (left_cw != nullptr)
        current_neighbors.insert(current_neighbors.end(), left_cw->pts.begin(),
                                 left_cw->pts.end());
      if (right_cw != nullptr)
        current_neighbors.insert(current_neighbors.end(), right_cw->pts.begin(),
                                 right_cw->pts.end());
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
      if (candidate_id < l_bound || candidate_id > r_bound) {
        continue;
      }
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

double PartialIndex::calSize() {
  double unit_size = sizeof(this->compact_arr.front().first);
  // cout << "unit size: " << unit_size << endl;
  return unit_size * nodes_amount;
}

pair<vector<int>, vector<int>> PartialIndex::checkIndexQuality(const int idx,
                                                             const int lbound,
                                                             const int rbound) {
  vector<int> lres, rres;
  CompactWindow *left_cw = compact_arr[idx].first.searchCW(lbound, true);
  CompactWindow *right_cw = compact_arr[idx].second.searchCW(rbound);

  // if (left_cw == nullptr) {
  //   auto temp = compact_arr[idx].first.vecs.front();
  //   print_set(temp.second->pts);
  // }
  // cout << compact_arr[idx].first.vecs.size() << endl;

  // for (auto cw : compact_arr[idx].first.vecs) {
  //   cout << cw.first << ":";
  //   cw.second->print();
  // }

  if (left_cw != nullptr)
    lres.insert(lres.end(), left_cw->pts.begin(), left_cw->pts.end());
  if (right_cw != nullptr)
    rres.insert(rres.end(), right_cw->pts.begin(), right_cw->pts.end());
  return make_pair(lres, rres);
}

// rebuild KNNG from [left,right]
void PartialIndex::rebuildRangeKNNG(const int left, const int right,
                                  vector<vector<int>> &knng) {
  knng.clear();
  knng.resize(right - left + 1);
#pragma omp parallel for
  for (unsigned n = left; n <= right; n++) {
    CompactWindow *left_cw = compact_arr[n].first.searchCW(left, true);
    CompactWindow *right_cw = compact_arr[n].second.searchCW(right);
    auto &knn = knng[n - left];
    if (left_cw != nullptr) {
      knn.insert(knn.end(), left_cw->pts.begin(), left_cw->pts.end());
    }
    if (right_cw != nullptr) {
      knn.insert(knn.end(), right_cw->pts.begin(), right_cw->pts.end());
    }
  }
}
