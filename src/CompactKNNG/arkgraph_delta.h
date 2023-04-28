
/**
 * ARKGraph Delta
 * Compact Adjacent List + Partial Ranges + Delta Compression (Delta Node)
 *
 * Author: Chaoji Zuo
 * Email:  chaoji.zuo@rutgers.edu
 * Date:   June 20, 2022
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
// extern double tree_search_time, neighbor_search_time;
// #define ONLY_COUNT_WINDOWS

namespace deltaindex {

class DeltaIndex : public BaseIndex {
 public:
  vector<pair<vector<DeltaNode>, vector<DeltaNode>>> indexed_arr;

 public:
  int isSparse = false;
  int hop_limit;
  float avg_hop;

 public:
  DeltaIndex() {
    window_count = 0;
    nodes_amount = 0;
    sort_time = 0;
    hop_limit = 10;
  };
  void collectWindows(vector<DeltaNode> &left_pairs,
                      vector<DeltaNode> &right_pairs, vector<int> &l_arr,
                      vector<int> &r_arr, vector<int> &visited_left,
                      vector<int> &visited_right, const vector<int> &idx_sorted,
                      const int k_smallest, const int query_pos);
  void indexNodes(const int query_pos, const vector<vector<float>> &data_nodes,
                  const vector<int> &idx_sorted, const int k_smallest,
                  const int l_bound, const int r_bound,
                  vector<DeltaNode> &left_pairs,
                  vector<DeltaNode> &right_pairs);

  void buildIndex(const vector<vector<float>> &nodes, const int build_knn_k,
                  const int index_K);


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
    for (auto &ele : indexed_arr) {
      vector<DeltaNode>().swap(ele.first);
      vector<DeltaNode>().swap(ele.second);
    }
    vector<pair<vector<DeltaNode>, vector<DeltaNode>>>().swap(indexed_arr);
    indexed_arr.clear();
    window_count = 0;
  }

  void saveIndex(const string &filename);

  ~DeltaIndex() {
    for (auto &ele : indexed_arr) {
      vector<DeltaNode>().swap(ele.first);
      vector<DeltaNode>().swap(ele.second);
    }
    vector<pair<vector<DeltaNode>, vector<DeltaNode>>>().swap(indexed_arr);
  }
};

void DeltaIndex::saveIndex(const string &filepath) {
  ofstream file;
  cout << "Saving index to " << filepath << "  ..." << endl;
  file.open(filepath, std::ios_base::app);
  if (file) {
    for (auto vec : indexed_arr) {
      for (auto lls : vec.first) {
        file << lls.add << " " << lls.sub << ",";
      }
      file << ";";
      for (auto rrs : vec.second) {
        file << rrs.add << " " << rrs.sub << ",";
      }
      file << "\n";
    }
  }
  file.close();
  cout << "Done!" << endl;
}

void DeltaIndex::collectWindows(vector<DeltaNode> &left_pairs,
                                vector<DeltaNode> &right_pairs,
                                vector<int> &l_arr, vector<int> &r_arr,
                                vector<int> &visited_left,
                                vector<int> &visited_right,
                                const vector<int> &idx_sorted,
                                const int k_smallest, const int query_pos) {
  for (auto idx : idx_sorted) {
    DeltaNode a_node;
    if (l_arr.size() >= (k_smallest) && idx < l_arr[0]) {
      continue;
    }
    if (r_arr.size() >= (k_smallest) && idx > r_arr[r_arr.size() - 1]) {
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
      a_node.add = idx;
      if (l_arr.size() > k_smallest) {
        a_node.sub = l_arr.front();
        l_arr.erase(l_arr.begin());
      }
      left_pairs.emplace_back(a_node);

    } else if (idx > query_pos) {
      bool sign = 0;
      a_node.sub = INT_MAX;
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
      a_node.add = idx;
      if (r_arr.size() > k_smallest) {
        a_node.sub = r_arr.back();
        r_arr.erase(r_arr.end() - 1);
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
  vector<int> l_arr;
  vector<int> r_arr;

  // store visited indexes to generate tree.
  vector<int> visited_left;
  vector<int> visited_right;

  l_arr.emplace_back(-1);
  r_arr.emplace_back(r_bound + 1);

  collectWindows(left_pairs, right_pairs, l_arr, r_arr, visited_left,
                 visited_right, idx_sorted, k_smallest, query_pos);

  // timeval temp1, temp2;
  // gettimeofday(&temp1, NULL);
  // // doesn't care about nearest negibhors
  // // second select.
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
  // gettimeofday(&temp2, NULL);
  // AccumulateTime(temp1, temp2, sort_time);

  // collectWindows(left_pairs, right_pairs, l_arr, r_arr, visited_left,
  //                visited_right, early_end_idx_re_sorted, k_smallest,
  //                query_pos);
  // // cout << "second collect" << endl;

  window_count = window_count + left_pairs.size() + right_pairs.size();

  return;
}

void DeltaIndex::buildIndex(const vector<vector<float>> &nodes,
                            const int initial_build_k, const int index_k = 10) {
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
  double avg_neighbor = 0;
  timeval tt1, tt2;
  gettimeofday(&tt1, NULL);

  // use hierachical nn descent method
  // HVicDescent<vector<float>> hd;
  // hd.setK(initial_build_k);
  // auto idxes = hd.HierachicalcomputeGraph(nodes, distance);
  // avg_neighbor = hd.avg_neighbor;

  // use kgraph to build knn graph
  kgraph::Matrix<float> data_mat;
  data_mat.resize(nodes.size(), nodes.front().size());
  for (unsigned i = 0; i < nodes.size(); ++i) {
    float *row = data_mat[i];
    for (unsigned j = 0; j < nodes.front().size(); ++j) {
      row[j] = nodes.at(i).at(j);
    }
  }

  kgraph::MatrixOracle<float, kgraph::metric::l2sqr> oracle(data_mat);

  // use vector oracle rather than matrix, no simd optimization
  // typedef kgraph::VectorOracle<vector<vector<float>>, vector<float>>
  // MyOracle; MyOracle oracle(nodes, distance);

  avg_neighbor = initial_build_k;
  kgraph::verbosity = 0;
  kgraph::KGraph::IndexInfo info;
  kgraph::KGraph *kgraph_index = kgraph::KGraph::create();

  {
    kgraph::KGraph::IndexParams params;
    params.K = initial_build_k;
    // params.K = 2;
    kgraph_index->build(oracle, params, &info);
    avg_neighbor = info.M;
  }
  // vector<vector<int>> idxes;
  // idxes.resize(nodes.size());
  // for (unsigned i = 0; i < nodes.size(); ++i) {
  //   vector<unsigned> one_nns(100);
  //   unsigned pM, pL;
  //   kgraph_index->get_nn(i, &one_nns[0], &pM, &pL);
  //   for (unsigned j = 0; j < pM; j++) {
  //     idxes.at(i).emplace_back((int)one_nns.at(j));
  //   }
  // }
  // avg_neighbor = initial_build_k;

  gettimeofday(&tt2, NULL);
  if (isLog) {
    logTime(tt1, tt2, "Compute Graph Time");
    cout << "# of average neighbors: " << avg_neighbor << endl;
  }
  gettimeofday(&tt1, NULL);
  int hop_amount = 0;

  double sort_time_ = 0;

#pragma omp parallel for schedule(dynamic) reduction(+ : sort_time_)
  for (int idx = 0; idx < nodes.size(); idx++) {
    timeval tt3, tt4;
    gettimeofday(&tt3, NULL);
    // IncrementalNeighbors IN(idx, &nodes, &idxes, avg_neighbor);
    IncrementalNeighbors_Kgraph IN(idx, &nodes, kgraph_index, avg_neighbor);
    IN.findAllNeighbors(hop_limit, index_k);

    auto idx_sorted = IN.sorted_idxes;
    hop_amount += IN.getHop();
    gettimeofday(&tt4, NULL);
    double one_sort_time;
    CountTime(tt3, tt4, one_sort_time);
    sort_time_ += one_sort_time;

    vector<DeltaNode> left_pair_vec, right_pair_vec;
    indexNodes(idx, nodes, idx_sorted, index_k, 0, nodes.size() - 1,
               left_pair_vec, right_pair_vec);
    indexed_arr[idx].first = left_pair_vec;
    indexed_arr[idx].second = right_pair_vec;
  }

  delete kgraph_index;
  sort_time = sort_time_;
  if (isLog) {
    cout << "Sort Time: " << sort_time << endl;
    logTime(tt1, tt2, "Index Time");
    cout << "average hop: " << (float)hop_amount / nodes.size() << endl;
    cout << "index amount: " << window_count << endl;
  }
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
    if (length <= 3) {
      current_pos = length - 1;
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
    if (length <= 3) {
      current_pos = length - 1;
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
  timeval tt1, tt2;
  // int query_range = r_bound - l_bound;
  // float range_scale = (float)query_range / data_nodes.size();

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
    vector<int> left_neighbors, right_neighbors;

    left_neighbors = decompressDeltaPath(indexed_arr[current_node_id].first,
                                         l_bound, K2_, true);
    right_neighbors = decompressDeltaPath(indexed_arr[current_node_id].second,
                                          r_bound, K2_, false);
    vector<int> current_neighbors(left_neighbors);
    current_neighbors.insert(current_neighbors.end(), right_neighbors.begin(),
                             right_neighbors.end());

    // improve the connectivity
    // if ((0.2 * data_nodes.size() < (current_node_id - l_bound))) {
    //   left_neighbors =
    //   decompressDeltaPath(indexed_arr[current_node_id].first,
    //                                        l_bound, 5 * K2_, true);
    //   current_neighbors.insert(current_neighbors.end(),
    //   left_neighbors.begin(),
    //                            left_neighbors.end());
    // }
    // if ((0.2 * data_nodes.size() > (r_bound - current_node_id))) {
    //   right_neighbors =
    //   decompressDeltaPath(indexed_arr[current_node_id].second,
    //                                         r_bound, 5 * K2_, false);
    //   current_neighbors.insert(current_neighbors.end(),
    //   right_neighbors.begin(),
    //                            right_neighbors.end());
    // }

    // if ((0.4 * data_nodes.size() < (current_node_id - l_bound))) {
    //   left_neighbors =
    //   decompressDeltaPath(indexed_arr[current_node_id].first,
    //                                        l_bound, 100 * K2_, true);
    //   current_neighbors.insert(current_neighbors.end(),
    //   left_neighbors.begin(),
    //                            left_neighbors.end());
    // }
    // if ((0.4 * data_nodes.size() > (r_bound - current_node_id))) {
    //   right_neighbors =
    //   decompressDeltaPath(indexed_arr[current_node_id].second,
    //                                         r_bound, 100 * K2_, false);
    //   current_neighbors.insert(current_neighbors.end(),
    //   right_neighbors.begin(),
    //                            right_neighbors.end());
    // }

    // print_set(left_neighbors);
    // print_set(right_neighbors);

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
  sort(lres.begin(), lres.end());
  sort(rres.begin(), rres.end());
  print_set(lres);
  print_set(rres);

  return make_pair(lres, rres);
}

// rebuild KNNG from [left,right]
void DeltaIndex::rebuildRangeKNNG(const int left, const int right,
                                  vector<vector<int>> &knng) {
  knng.clear();
  knng.resize(right - left + 1);
#pragma omp parallel for
  for (unsigned n = left; n <= right; n++) {
    vector<int> lres, rres;
    lres = decompressDeltaPath(indexed_arr[n].first, left, K2_, true);
    rres = decompressDeltaPath(indexed_arr[n].second, right, K2_, false);
    auto &knn = knng[n - left];
    knn.insert(knn.end(), lres.begin(), lres.end());
    knn.insert(knn.end(), rres.begin(), rres.end());
  }
}

// class SparseLayer {
// public:
//   SparseLayer(){};
//   DeltaIndex index;
//   void buildSparseIndex(const vector<vector<float>> &nodes,
//                         const int initial_build_k, const int index_k);
//   vector<int> findEnters(const vector<float> &query, const int l_bound,
//                          const int r_bound);

// private:
//   int step_size = 10;
//   vector<vector<float>> sub_nodes;
// };

// void SparseLayer::buildSparseIndex(const vector<vector<float>> &nodes,
//                                    const int initial_build_k = 20,
//                                    const int index_k = 10) {
//   sub_nodes.clear();
//   auto start = nodes.begin();
//   size_t counter = 0;
//   while (counter < nodes.size()) {
//     sub_nodes.emplace_back(*start);
//     advance(start, step_size);
//     counter += step_size;
//   }
//   index.isLog = false;
//   index.buildIndex(sub_nodes, 20, 12);
// }

// vector<int> SparseLayer::findEnters(const vector<float> &query,
//                                     const int l_bound, const int r_bound) {
//   auto res = index.indexNeighborSearch(sub_nodes, query, l_bound /
//   step_size,
//                                        r_bound / step_size, 12, 80);
//   for (auto &ele : res) {
//     ele = ele * step_size;
//   }
//   return res;
// }

}  // namespace deltaindex