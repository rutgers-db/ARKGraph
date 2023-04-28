/**
 * Use kgraph to find more neighbors
 *
 * Author: Chaoji Zuo
 * Email:  chaoji.zuo@rutgers.edu
 * Date:   Sep 3, 2022
 */

#pragma once

#include <algorithm>
#include <iostream>
#include <unordered_set>
#include <vector>

#include "kgraph/kgraph.h"
#include "utils/utils.h"

using std::cout;
using std::endl;
using std::pair;
using std::string;
using std::unordered_set;
using std::vector;

#define _INT_MAX 2147483640

class IncrementalNeighbors_Kgraph {
 public:
  vector<int> sorted_idxes;
  // vector<int> old_candidate_idxes;
  double nn_time;
  double scan_time;
  int nn_amount;
  int scan_amount;
  vector<float> old_candidate_dists;

 private:
  int query_point;
  int cur_neighbors;
  int hop = 0;
  unordered_set<int> candidates;
  // queue<int> neigbors_queue;
  vector<int> neighbors_list;
  const vector<vector<float>> *p_nodes;
  // vector<vector<int>> *p_idxes;
  kgraph::KGraph *p_index;

  vector<float> left_candidate_dists;
  vector<int> left_candidate_idxes;
  vector<float> right_candidate_dists;
  vector<int> right_candidate_idxes;
  int nodes_size;
  int left_min;
  int right_max;
  double avg_neighbor;
  int dim;

 public:
  IncrementalNeighbors_Kgraph(const int query_point,
                              const vector<vector<float>> *nodes,
                              kgraph::KGraph *index,
                              const double average_neighbor)
      : query_point(query_point) {
    // neigbors_queue.push(query_point);
    neighbors_list.emplace_back(query_point);
    p_nodes = nodes;
    cur_neighbors = 0;
    nodes_size = nodes->size();
    left_min = 0;
    p_index = index;
    right_max = nodes_size;

    avg_neighbor = average_neighbor;
    dim = nodes->front().size();
  };
  void findMoreNeighbors(const int next);
  int size() { return cur_neighbors; }
  int nextNeighbor(const int next);
  pair<int, int> nndecentNeighbors(const int, const int);
  bool DetermineEarlyEnd(const int);
  int getHop() { return hop; }
  vector<int> findAllNeighbors(const int, const int K_index);
  vector<int> findHNSWNeighbors(const int K_index);
  // for sparse hierarchical layers
  pair<int, int> nndecentNeighbors(const vector<float> &query,
                                   const int hop_limit, const int k_index);
  vector<int> findAllNeighbors(const vector<float> &query, const int hop_limit,
                               const int K_index);
  void sort_candidates(vector<float> &v, vector<int> &original_idx) {
    vector<int> idx(v.size());
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    // using std::stable_sort instead of std::sort
    // to avoid unnecessary index re-orderings
    // when v contains elements of equal values
    stable_sort(idx.begin(), idx.end(),
                [&v](size_t i1, size_t i2) { return v[i1] < v[i2]; });

    vector<int> sorted_idx;
    vector<float> sorted_dists;
    for (auto &ele : idx) {
      sorted_idx.emplace_back(original_idx[ele]);
      sorted_dists.emplace_back(v[ele]);
    }
    v = sorted_dists;
    original_idx = sorted_idx;
  }
  ~IncrementalNeighbors_Kgraph() {
    candidates.clear();
    neighbors_list.clear();
    // free(p_idxes);
    // old_candidate_dists.clear();
  }
};

bool IncrementalNeighbors_Kgraph::DetermineEarlyEnd(const int k_index) {
  int left_max = 0;
  int right_min = nodes_size;
  left_min = query_point;
  right_max = query_point;
  for (auto ele : candidates) {
    if (ele < query_point && ele > left_max) {
      left_max = ele;
    } else if (ele > query_point && ele < right_min) {
      right_min = ele;
    }
    // if (ele < left_min) {
    //   left_min = ele;
    // } else if (ele > right_max) {
    //   right_max = ele;
    // }
  }
  left_min = *min_element(candidates.begin(), candidates.end());
  right_max = *max_element(candidates.begin(), candidates.end());

  // // update update cost model
  // // unordered_set + next_max(even distribution) + total_visit_number
  // float next_hop_complexity = neighbors_list.size() * avg_neighbor *
  //                             (right_max - left_min) / nodes_size * dim;
  // int unvisit_range = (right_min - left_max) * dim;

  // if (neighbors_list.size() > 0 && p_idxes->front().size() > 0) {
  //   int step_size =
  //       nodes_size / (neighbors_list.size() * p_idxes->front().size());
  //   if (step_size > 0) {
  //     float greedy_after_next_hop =
  //         query_point % step_size +
  //         ((query_point / step_size) + 1) * step_size - query_point;
  //     next_hop_complexity += greedy_after_next_hop;
  //   }

  //   // cout << "unvisited range: " << unvisit_range << endl;
  //   // cout << "next hop amount: " << next_hop_complexity << endl;

  //   if (unvisit_range > next_hop_complexity) {
  //     return false;
  //   } else
  //     return true;
  // }
  // return true;

  // // update cost model
  // // unordered_set + next_max (even distribution) + total_visit_number
  // //   float next_hop_complexity = neighbors_list.size() * avg_degree *
  // //                             (right_max - left_min) / nodes_size *d +
  // //                             outside_size;
  // // int unvisit_range = right_min - left_max *d;
  // float next_hop_complexity = neighbors_list.size() *
  // p_idxes->front().size()
  // *
  //                             (right_max - left_min) / nodes_size;
  // int unvisit_range = right_min - left_max;
  // if (neighbors_list.size() > 0 && p_idxes->front().size() > 0) {
  //   // int step_size =
  //   //     nodes_size / (neighbors_list.size() *
  //   p_idxes->front().size());
  //   // if (step_size > 0) {
  //   //   float greedy_after_next_hop =
  //   //       query_point % step_size +
  //   //       ((query_point / step_size) + 1) * step_size - query_point;
  //   //   next_hop_complexity += greedy_after_next_hop;
  //   // }

  int unvisit_range = right_min - left_max;
  int next_hop_complexity = neighbors_list.size() * avg_neighbor;
  if (unvisit_range > next_hop_complexity) {
    return false;
  } else
    return true;
}

// get the neighbors around the query point
pair<int, int> IncrementalNeighbors_Kgraph::nndecentNeighbors(
    const int hop_limit, const int k_index) {
  // int left_bound = query_point - 1;
  // int right_bound = query_point + 1;
  // if (left_bound < 0) {
  //   left_bound = right_bound;
  // }
  // if (right_bound == p_nodes->size()) {
  //   right_bound = left_bound;
  // }
  // auto checkBound = [this]() {
  //   int left_bound = query_point - 1;
  //   int right_bound = query_point + 1;
  //   if (candidates.find(left_bound) != candidates.end() &&
  //       candidates.find(right_bound) != candidates.end()) {
  //     return true;
  //   }
  //   return false;
  // };
  left_candidate_idxes.clear();
  left_candidate_dists.clear();
  right_candidate_idxes.clear();
  right_candidate_dists.clear();

  vector<int> left_neighbors_list;
  vector<int> right_neighbors_list;

  int left_pos = -1;
  int right_pos = nodes_size;

  if (hop_limit == 0) {
    return std::make_pair(left_pos, right_pos);
  }
  assert(neighbors_list.size() > 0);
  vector<float> query(p_nodes->at(query_point));
  while (hop < hop_limit) {
    // while (!checkBound()) {
    left_neighbors_list.clear();
    right_neighbors_list.clear();
    if (neighbors_list.empty()) {
      break;
    }

    for (auto cur_node : neighbors_list) {
      vector<unsigned> cur_node_neighbors(501);
      vector<float> cur_node_dists(501);
      unsigned pM, pL;
      p_index->get_nn((unsigned)cur_node, &cur_node_neighbors[0],
                      &cur_node_dists[0], &pM, &pL);
      // assert(pL > 500);
      // if (pL > 200) {
      //   pL = 200;
      // }
      for (unsigned i = 0; i < pM; i++) {
        int point = (int)cur_node_neighbors.at(i);
        assert(point < p_nodes->size());  // TODO: check why this happen
        if (candidates.find(point) == candidates.end()) {
          if (point == query_point) {
            continue;
          }
          candidates.insert(point);
          float dist = EuclideanDistance(query, p_nodes->at(point));
          if (point < query_point) {
            left_neighbors_list.emplace_back(point);
            left_candidate_idxes.emplace_back(point);
            left_candidate_dists.emplace_back(dist);
          } else {
            right_neighbors_list.emplace_back(point);
            right_candidate_idxes.emplace_back(point);
            right_candidate_dists.emplace_back(dist);
          }
        }
      }
    }
    hop += 1;
    neighbors_list.assign(left_neighbors_list.begin(),
                          left_neighbors_list.end());
    neighbors_list.insert(neighbors_list.end(), right_neighbors_list.begin(),
                          right_neighbors_list.end());
    cur_neighbors = candidates.size();

    // determineEarlyEnd
    if (true) {
      vector<int> left_sort_idx(left_candidate_idxes);
      vector<int> right_sort_idx(right_candidate_idxes);
      if (left_sort_idx.size() >= k_index) {
        nth_element(left_sort_idx.begin(), left_sort_idx.begin() + k_index - 1,
                    left_sort_idx.end(), std::greater<int>{});
        left_pos = left_sort_idx[k_index - 1];
      }
      if (right_sort_idx.size() >= k_index) {
        nth_element(right_sort_idx.begin(),
                    right_sort_idx.begin() + k_index - 1, right_sort_idx.end());
        right_pos = right_sort_idx.at(k_index - 1);
      }
      int unvisit_range = right_pos - left_pos;
      int next_hop_complexity = neighbors_list.size() * avg_neighbor;
      if (hop_limit == 10) {
        if (unvisit_range < next_hop_complexity) {
          break;
        } else
          continue;
      }
    }
  }

  return std::make_pair(left_pos, right_pos);

  // sort_candidates(left_candidate_dists, left_candidate_idxes);
  // old_candidate_dists.assign(left_candidate_dists.begin(),
  //                            left_candidate_dists.end());
  // sorted_idxes.assign(left_candidate_idxes.begin(),
  // left_candidate_idxes.end());

  // sort_candidates(right_candidate_dists, right_candidate_idxes);
  // old_candidate_dists.insert(old_candidate_dists.end(),
  //                            right_candidate_dists.begin(),
  //                            right_candidate_dists.end());
  // sorted_idxes.insert(sorted_idxes.end(), right_candidate_idxes.begin(),
  //                     right_candidate_idxes.end());
}

vector<int> IncrementalNeighbors_Kgraph::findAllNeighbors(const int hop_limit,
                                                          const int K_index) {
  timeval t1, t2;
  gettimeofday(&t1, NULL);
  auto nn_result = this->nndecentNeighbors(hop_limit, K_index);
  gettimeofday(&t2, NULL);
  CountTime(t1, t2, nn_time);

  nn_amount = left_candidate_idxes.size() + right_candidate_idxes.size();

  gettimeofday(&t1, NULL);
  // vector<int> left_sort_idx(left_candidate_idxes);
  // vector<int> right_sort_idx(right_candidate_idxes);
  // int left_pos = -1;
  // int right_pos = nodes_size;
  // if (left_sort_idx.size() >= K_index) {
  //   nth_element(left_sort_idx.begin(), left_sort_idx.begin() + K_index - 1,
  //               left_sort_idx.end(), greater<int>{});
  //   left_pos = left_sort_idx[K_index - 1];
  // }
  // if (right_sort_idx.size() >= K_index) {
  //   nth_element(right_sort_idx.begin(), right_sort_idx.begin() + K_index - 1,
  //               right_sort_idx.end());
  //   right_pos = right_sort_idx.at(K_index - 1);
  // }
  int left_pos = nn_result.first;
  int right_pos = nn_result.second;
  // print_set(right_sort_idx);
  // cout << left_pos << "->" << right_pos << endl;
  vector<float> dist_arr;
  vector<int> second_sorted_idxes;
  for (auto i = left_pos + 1; i < right_pos; i++) {
    if (find(left_candidate_idxes.begin(), left_candidate_idxes.end(), i) ==
            left_candidate_idxes.end() &&
        find(right_candidate_idxes.begin(), right_candidate_idxes.end(), i) ==
            right_candidate_idxes.end()) {
      dist_arr.emplace_back(
          EuclideanDistance(p_nodes->at(query_point), p_nodes->at(i), 0, 0));
      second_sorted_idxes.emplace_back(i);
    }
  }
  gettimeofday(&t2, NULL);
  CountTime(t1, t2, scan_time);
  scan_amount = right_pos - left_pos + 1;

  sorted_idxes.assign(left_candidate_idxes.begin(), left_candidate_idxes.end());
  sorted_idxes.insert(sorted_idxes.end(), right_candidate_idxes.begin(),
                      right_candidate_idxes.end());
  sorted_idxes.insert(sorted_idxes.end(), second_sorted_idxes.begin(),
                      second_sorted_idxes.end());

  old_candidate_dists.assign(left_candidate_dists.begin(),
                             left_candidate_dists.end());
  old_candidate_dists.insert(old_candidate_dists.end(),
                             right_candidate_dists.begin(),
                             right_candidate_dists.end());
  old_candidate_dists.insert(old_candidate_dists.end(), dist_arr.begin(),
                             dist_arr.end());

  IncrementalNeighbors_Kgraph::sort_candidates(old_candidate_dists,
                                               sorted_idxes);

  return sorted_idxes;
}

// use HNSW to build base graph
vector<int> IncrementalNeighbors_Kgraph::findHNSWNeighbors(const int K_index) {
  timeval t1, t2;
  gettimeofday(&t1, NULL);
  // auto nn_result = this->nndecentNeighbors(1, K_index);
  gettimeofday(&t2, NULL);
  CountTime(t1, t2, nn_time);
  int left_pos = -1;
  int right_pos = nodes_size;

  nn_amount = left_candidate_idxes.size() + right_candidate_idxes.size();
  sort(left_candidate_idxes.begin(), left_candidate_idxes.end());
  if (left_candidate_idxes.size() >= K_index) {
    left_candidate_idxes.resize(left_candidate_idxes.size() - K_index + 1);
    left_pos = left_candidate_idxes.back();
  }
  sort(right_candidate_idxes.begin(), right_candidate_idxes.end());
  if (right_candidate_idxes.size() >= K_index) {
    right_candidate_idxes.erase(right_candidate_idxes.begin(),
                                right_candidate_idxes.begin() + K_index - 1);
    right_pos = right_candidate_idxes.front();
  }
  // cout << left_pos << " , " << right_pos << endl;
  scan_amount = right_pos - left_pos + 1;

  gettimeofday(&t1, NULL);
  vector<float> dist_arr;
  vector<int> second_sorted_idxes;
  for (auto i = left_pos + 1; i < right_pos; i++) {
    dist_arr.emplace_back(
        EuclideanDistance(p_nodes->at(query_point), p_nodes->at(i), 0, 0));
    second_sorted_idxes.emplace_back(i);
  }
  IncrementalNeighbors_Kgraph::sort_candidates(dist_arr, second_sorted_idxes);
  sorted_idxes.assign(left_candidate_idxes.begin(), left_candidate_idxes.end());
  sorted_idxes.insert(sorted_idxes.end(), second_sorted_idxes.begin(),
                      second_sorted_idxes.end());
  sorted_idxes.insert(sorted_idxes.end(), right_candidate_idxes.begin(),
                      right_candidate_idxes.end());
  return sorted_idxes;
}

// get the neighbors around the query point, for sparse layer
pair<int, int> IncrementalNeighbors_Kgraph::nndecentNeighbors(
    const vector<float> &query, const int hop_limit, const int k_index) {
  left_candidate_idxes.clear();
  left_candidate_dists.clear();
  right_candidate_idxes.clear();
  right_candidate_dists.clear();

  vector<int> left_neighbors_list;
  vector<int> right_neighbors_list;

  int left_pos = -1;
  int right_pos = nodes_size;

  if (hop_limit == 0) {
    return std::make_pair(left_pos, right_pos);
  }
  assert(neighbors_list.size() > 0);
  while (hop < hop_limit) {
    // while (!checkBound()) {
    left_neighbors_list.clear();
    right_neighbors_list.clear();
    if (neighbors_list.empty()) {
      break;
    }
    for (auto cur_node : neighbors_list) {
      vector<unsigned> cur_node_neighbors(200);
      vector<float> cur_node_dists(200);
      unsigned pM, pL;

      p_index->get_nn((unsigned)cur_node, &cur_node_neighbors[0],
                      &cur_node_dists[0], &pM, &pL);

      for (unsigned i = 0; i < pM; i++) {
        auto point = cur_node_neighbors.at(i);

        if (candidates.find(point) == candidates.end()) {
          candidates.insert(point);
          float dist = EuclideanDistance(query, p_nodes->at(point));
          right_neighbors_list.emplace_back(point);
          right_candidate_idxes.emplace_back(point);
          right_candidate_dists.emplace_back(dist);
        }
      }
    }
    hop += 1;
    neighbors_list.assign(left_neighbors_list.begin(),
                          left_neighbors_list.end());
    neighbors_list.insert(neighbors_list.end(), right_neighbors_list.begin(),
                          right_neighbors_list.end());
    cur_neighbors = candidates.size();

    // determineEarlyEnd
    if (true) {
      vector<int> left_sort_idx(left_candidate_idxes);
      vector<int> right_sort_idx(right_candidate_idxes);
      if (left_sort_idx.size() >= k_index) {
        nth_element(left_sort_idx.begin(), left_sort_idx.begin() + k_index - 1,
                    left_sort_idx.end(), std::greater<int>{});
        left_pos = left_sort_idx[k_index - 1];
      }
      if (right_sort_idx.size() >= k_index) {
        nth_element(right_sort_idx.begin(),
                    right_sort_idx.begin() + k_index - 1, right_sort_idx.end());
        right_pos = right_sort_idx.at(k_index - 1);
      }
      int unvisit_range = right_pos - left_pos;
      int next_hop_complexity = neighbors_list.size() * avg_neighbor;
      if (hop_limit == 10) {
        if (unvisit_range < next_hop_complexity) {
          // cout << "early end" << endl;
          break;
        } else
          continue;
      }
    }
  }

  return std::make_pair(left_pos, right_pos);
}

vector<int> IncrementalNeighbors_Kgraph::findAllNeighbors(
    const vector<float> &query, const int hop_limit, const int K_index) {
  neighbors_list.clear();
  // add initial points
  int interval = (p_nodes->size()) / 3;
  for (size_t i = 1; i <= 2; i++) {
    int point = 0 + interval * i;
    neighbors_list.emplace_back(point);
  }

  timeval t1, t2;
  gettimeofday(&t1, NULL);
  auto nn_result = this->nndecentNeighbors(query, hop_limit, K_index);
  gettimeofday(&t2, NULL);
  CountTime(t1, t2, nn_time);
  // cout << "nn time:" << nn_time << endl;
  // cout << nn_result.first << "," << nn_result.second << endl;
  // assert(false);

  nn_amount = left_candidate_idxes.size() + right_candidate_idxes.size();

  gettimeofday(&t1, NULL);
  int left_pos = nn_result.first;
  int right_pos = nn_result.second;
  vector<float> dist_arr;
  vector<int> second_sorted_idxes;
  for (auto i = left_pos + 1; i < right_pos; i++) {
    if (find(left_candidate_idxes.begin(), left_candidate_idxes.end(), i) ==
            left_candidate_idxes.end() &&
        find(right_candidate_idxes.begin(), right_candidate_idxes.end(), i) ==
            right_candidate_idxes.end()) {
      dist_arr.emplace_back(EuclideanDistance(query, p_nodes->at(i), 0, 0));
      second_sorted_idxes.emplace_back(i);
    }
  }
  gettimeofday(&t2, NULL);
  CountTime(t1, t2, scan_time);
  scan_amount = right_pos - left_pos + 1;

  sorted_idxes.assign(left_candidate_idxes.begin(), left_candidate_idxes.end());
  sorted_idxes.insert(sorted_idxes.end(), right_candidate_idxes.begin(),
                      right_candidate_idxes.end());
  sorted_idxes.insert(sorted_idxes.end(), second_sorted_idxes.begin(),
                      second_sorted_idxes.end());

  old_candidate_dists.assign(left_candidate_dists.begin(),
                             left_candidate_dists.end());
  old_candidate_dists.insert(old_candidate_dists.end(),
                             right_candidate_dists.begin(),
                             right_candidate_dists.end());
  old_candidate_dists.insert(old_candidate_dists.end(), dist_arr.begin(),
                             dist_arr.end());

  IncrementalNeighbors_Kgraph::sort_candidates(old_candidate_dists,
                                               sorted_idxes);

  return sorted_idxes;
}