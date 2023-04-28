

/**
 * ARKGraph Delta
 * Compact Adjacent List + Partial Ranges + Delta Compression (Delta Node)
 * Heap version. (Heap in collectwindows function)
 *
 * Author: Chaoji Zuo
 * Email:  chaoji.zuo@rutgers.edu
 * Date:   June 20, 2022
 */
#pragma once
#pragma once

#include <algorithm>
#include <iostream>
#include <numeric>
#include <queue>
#include <vector>

#include "../NearestSearch.h"
#include "../NearestSearch_kgraph.h"
#include "../baselines/greedy.h"
#include "../baselines/knn_first_baseline.h"
#include "../range_index_base.h"
#include "../utils/cw.h"
#include "../utils/segmentTree.hpp"
#include "../utils/utils.h"
#include "hnswlib/hnswlib.h"
#include "kgraph/kgraph-data.h"
#include "kgraph/kgraph.h"
#include "myHNSW/hnswlib-incrementally/hnswlib.h"
#include "myHNSW/hnswlib-nsw/hnswlib.h"

using std::cout;
using std::endl;
using std::vector;
namespace deltaindex_heap {

class DeltaIndex : public BaseIndex {
 private:
  vector<pair<vector<DeltaNode>, vector<DeltaNode>>> indexed_arr;
  int ltemp;
  int rtemp;

 public:
  DeltaIndex() {
    window_count = 0;
    nodes_amount = 0;
    sort_time = 0;
  };
  void collectWindows(vector<DeltaNode> &left_pairs,
                      vector<DeltaNode> &right_pairs,
                      priority_queue<int, vector<int>, greater<int>> &l_arr,
                      priority_queue<int, vector<int>, less<int>> &r_arr,
                      vector<int> &visited_left, vector<int> &visited_right,
                      const vector<int> &idx_sorted, const int k_smallest,
                      const int query_pos);
  void indexNodes(const int query_pos, const vector<vector<float>> &data_nodes,
                  const vector<int> &idx_sorted, const int k_smallest,
                  const int l_bound, const int r_bound,
                  vector<DeltaNode> &left_pairs,
                  vector<DeltaNode> &right_pairs);

  void buildIndex(const vector<vector<float>> &nodes, const int build_knn_k,
                  const int hop_limit);

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

  void clear() {
    indexed_arr.clear();
    window_count = 0;
  }
};

void DeltaIndex::collectWindows(
    vector<DeltaNode> &left_pairs, vector<DeltaNode> &right_pairs,
    priority_queue<int, vector<int>, greater<int>> &l_arr,
    priority_queue<int, vector<int>, less<int>> &r_arr,
    vector<int> &visited_left, vector<int> &visited_right,
    const vector<int> &idx_sorted, const int k_smallest, const int query_pos) {
  for (auto idx : idx_sorted) {
    DeltaNode a_node;
    if (l_arr.size() >= (k_smallest) && idx < l_arr.top()) {
      continue;
    }
    if (r_arr.size() >= (k_smallest) && idx > r_arr.top()) {
      continue;
    }
    int insert_pos = 0;
    if (idx < query_pos) {
      l_arr.push(idx);

      a_node.add = idx;
      if (l_arr.size() > k_smallest) {
        a_node.sub = l_arr.top();
        l_arr.pop();
      }
      if (idx > ltemp) {
        ltemp = idx;
      }
      left_pairs.emplace_back(a_node);

    } else if (idx > query_pos) {
      bool sign = 0;
      a_node.sub = INT_MAX;
      r_arr.push(idx);

      a_node.add = idx;
      if (r_arr.size() > k_smallest) {
        a_node.sub = r_arr.top();
        r_arr.pop();
      }
      if (idx < rtemp) {
        rtemp = idx;
      }

      right_pairs.emplace_back(a_node);
    }
  }
}

// index sorted vector
// idx_sorted: index vector sorted by distance to query point
void DeltaIndex::indexNodes(const int query_pos,
                            const vector<vector<float>> &data_nodes,
                            const vector<int> &idx_sorted, const int k_smallest,
                            const int l_bound, const int r_bound,
                            vector<DeltaNode> &left_pairs,
                            vector<DeltaNode> &right_pairs) {
  // vector<int> l_arr;
  // vector<int> r_arr;

  priority_queue<int, vector<int>, greater<int>> l_arr;
  priority_queue<int, vector<int>, less<int>> r_arr;

  // store visited indexes to generate tree.
  vector<int> visited_left;
  vector<int> visited_right;

  l_arr.push(-1);
  r_arr.push(r_bound + 1);
  ltemp = -1;
  rtemp = r_bound + 1;

  collectWindows(left_pairs, right_pairs, l_arr, r_arr, visited_left,
                 visited_right, idx_sorted, k_smallest, query_pos);
  // cout << "first collect" << endl;

  timeval temp1, temp2;
  gettimeofday(&temp1, NULL);
  // doesn't care about nearest negibhors
  // second select.
  int l_position = 0;
  int r_position = r_bound;
  l_position = ltemp;
  r_position = rtemp;
  // if (l_arr.size() > 1)
  //   // l_position = l_arr.back();
  //   l_position = *max_element(l_arr.begin(), l_arr.end());
  // if (r_arr.size() > 1)
  //   // r_position = r_arr.front();
  //   r_position = *min_element(r_arr.begin(), r_arr.end());
  vector<float> dist_arr;
  for (size_t i = l_position + 1; i < r_position; i++) {
    dist_arr.emplace_back(
        EuclideanDistance(data_nodes[query_pos], data_nodes[i], 0, 0));
  }
  vector<int> idx_re_sorted = sort_indexes(dist_arr);
  vector<int> early_end_idx_re_sorted;
  for (auto &idx : idx_re_sorted) {
    idx = idx + l_position + 1;
    early_end_idx_re_sorted.emplace_back(idx);
  }
  gettimeofday(&temp2, NULL);
  AccumulateTime(temp1, temp2, sort_time);

  collectWindows(left_pairs, right_pairs, l_arr, r_arr, visited_left,
                 visited_right, early_end_idx_re_sorted, k_smallest, query_pos);
  // cout << "second collect" << endl;
  // if (query_pos == 2600) {
  //   for (auto ele : l_arr) {
  //     cout << ele << ",";
  //   }
  // }
  // cout << l_arr.size() << "," << r_arr.size() << endl;
  window_count = window_count + left_pairs.size() + right_pairs.size();

  return;
}

void DeltaIndex::buildIndex(const vector<vector<float>> &nodes,
                            const int initial_build_k = 20,
                            const int index_k = 10) {
  auto distance = [](vector<float> a, vector<float> b) -> double {
    float ans = 0.0;
    for (int i = 0; i < a.size(); ++i) {
      ans += (a[i] - b[i]) * (a[i] - b[i]);
    }
    return ans;
  };
  indexed_arr.clear();
  indexed_arr.resize(nodes.size());
  K1_ = initial_build_k;
  K2_ = index_k;
  timeval tt1, tt2;
  gettimeofday(&tt1, NULL);
  kgraph::Matrix<float> data_mat;
  data_mat.resize(nodes.size(), nodes.front().size());
  for (unsigned i = 0; i < nodes.size(); ++i) {
    float *row = data_mat[i];
    for (unsigned j = 0; j < nodes.front().size(); ++j) {
      row[j] = nodes.at(i).at(j);
    }
  }
  kgraph::MatrixOracle<float, kgraph::metric::l2sqr> oracle(data_mat);

  float avg_neighbor = initial_build_k;
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
    vector<DeltaNode> left_pair_vec, right_pair_vec;
    indexNodes(idx, nodes, idx_sorted, index_k, 0, nodes.size() - 1,
               left_pair_vec, right_pair_vec);
    indexed_arr[idx].first = left_pair_vec;
    indexed_arr[idx].second = right_pair_vec;
    counter++;
  }
  delete kgraph_index;
  cout << "Sort Time: " << sort_time << endl;
  logTime(tt1, tt2, "Index Time");
  cout << "average hop: " << (float)hop_amount / nodes.size() << endl;
  cout << "index amount: " << window_count << endl;
}

// direction: true means from left to right, false means from right to left
vector<int> decompressDeltaPath(const vector<DeltaNode> &path,
                                const int target_position,
                                const int target_size, bool direction = true) {
  vector<int> neighbors;

  int length = path.size();
  int current_pos = length / 2;
  int current_left = 0;
  int current_right = length;
  // cout << "current visiting: " << target_position << endl;
  // for (auto ele : path) {
  //   cout << "-" << ele.sub << "\t +" << ele.add << endl;
  // }

  if (direction == true) {
    while (current_pos > 1 && current_pos < (length - 1)) {
      auto cur = path.at(current_pos);
      if (cur.sub == target_position) {
        current_pos -= 1;
        break;
      }
      if ((cur.sub < target_position) &&
          target_position <= path.at(current_pos + 1).sub) {
        break;
      }
      if (current_right - current_left < 2) {
        break;
      }
      if (cur.sub < target_position) {
        current_left = current_pos;
        current_pos += (current_right - current_left) / 2;
        continue;
      }
      if (cur.sub > target_position) {
        current_right = current_pos;
        current_pos -= (current_right - current_left) / 2;
        continue;
      }
    }
    while (current_pos >= 0 && current_pos < length) {
      auto cur = path.at(current_pos);
      if (cur.add >= target_position) neighbors.emplace_back(cur.add);
      current_pos--;
      if (neighbors.size() >= target_size) {
        break;
      }
    }
  } else {
    while (current_pos > 1 && current_pos < (length - 1)) {
      auto cur = path.at(current_pos);
      if (cur.sub == target_position) {
        current_pos -= 1;
        break;
      }
      if ((cur.sub > target_position) &&
          (target_position >= path.at(current_pos + 1).sub)) {
        break;
      }
      if (current_right - current_left < 2) {
        break;
      }
      if (cur.sub < target_position) {
        current_right = current_pos;
        current_pos -= (current_right - current_left) / 2;
        continue;
      }
      if (cur.sub > target_position) {
        current_left = current_pos;
        current_pos += (current_right - current_left) / 2;
        continue;
      }
    }
    while (current_pos >= 0 && current_pos < length) {
      auto cur = path.at(current_pos);
      if (cur.add <= target_position && cur.add > 0)
        neighbors.emplace_back(cur.add);
      current_pos--;
      if (neighbors.size() >= target_size) {
        break;
      }
    }
  }
  return neighbors;
}

vector<int> DeltaIndex::indexNeighborSearch(
    const vector<vector<float>> &data_nodes, const vector<float> &query,
    const int l_bound, const int r_bound, const int K_neigbhors, const int ef) {
  vector<int> temp;
  return indexNeighborSearch(data_nodes, query, l_bound, r_bound, K_neigbhors,
                             ef, temp);
}

// on the fly search cw
// TODO: Verify search path result
vector<int> DeltaIndex::indexNeighborSearch(
    const vector<vector<float>> &data_nodes, const vector<float> &query,
    const int l_bound, const int r_bound, const int K_neigbhors, const int ef,
    const vector<int> &enter_points) {
  unordered_map<int, bool> visited_list;
  float lower_bound = INT_MAX;
  priority_queue<pair<float, int>> top_candidates;
  priority_queue<pair<float, int>> candidate_set;

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
    enter_list.assign(enter_points.begin(), enter_points.end());
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
    auto left_neighbors = decompressDeltaPath(
        indexed_arr[current_node_id].first, l_bound, K2_, true);
    auto right_neighbors = decompressDeltaPath(
        indexed_arr[current_node_id].second, r_bound, K2_, false);

    // print_set(left_neighbors);
    // print_set(right_neighbors);

    vector<int> current_neighbors(left_neighbors);
    current_neighbors.insert(current_neighbors.end(), right_neighbors.begin(),
                             right_neighbors.end());
    // gettimeofday(&tt2, NULL);
    // AccumulateTime(tt1, tt2, tree_search_time);

    if (current_neighbors.size() == 0) {
      // cout << "current node: " << current_node_pair.second << "   ";
      // cout << "size 0 neighbors" << endl;
      continue;
    }
    // print_set(current_neighbors);
    // gettimeofday(&tt1, NULL);
    for (size_t i = 0; i < current_neighbors.size(); i++) {
      int candidate_id = current_neighbors[i];
      if (!visited_list[candidate_id]) {
        visited_list[candidate_id] = true;
        float dist = EuclideanDistance(query, data_nodes[candidate_id]);
        // cout << "candidate: " << candidate_id << "  dist: " << dist <<
        // endl;
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

double DeltaIndex::calSize() {
  double unit_size = sizeof(this->indexed_arr.front().first.front());
  // cout << "unit size: " << unit_size << endl;
  return unit_size * window_count;
}

pair<vector<int>, vector<int>> DeltaIndex::checkIndexQuality(const int idx,
                                                             const int lbound,
                                                             const int rbound) {
  vector<int> lres, rres;
  lres = decompressDeltaPath(indexed_arr[idx].first, lbound, K2_, true);
  rres = decompressDeltaPath(indexed_arr[idx].second, rbound, K2_, false);
  // for (auto ele : indexed_arr[idx].first) {
  //   cout << ele.add << "," << ele.sub << endl;
  // }

  return make_pair(lres, rres);
}

}  // namespace deltaindex_heap
