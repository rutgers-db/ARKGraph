
/**
 * ARKGraph Delta ST
 * Compact Adjacent List + Partial Ranges + Delta Compression (Segment Tree)
 *
 * Author: Chaoji Zuo
 * Email:  chaoji.zuo@rutgers.edu
 * Date:   Mar 8, 2022
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

struct TreeNode {
  int lid;
  int rid;
  TreeNode *leftChild;
  TreeNode *rightChild;
  vector<int> data;
  TreeNode() : lid(0), rid(0), leftChild(nullptr), rightChild(nullptr) {}
};

TreeNode *createTreeTopDown(const vector<int> &pts, int min, int max) {
  TreeNode *rootNode = new TreeNode();
  rootNode->lid = pts[min];
  rootNode->rid = pts[max];
  rootNode->leftChild = nullptr;
  rootNode->rightChild = nullptr;
  if (max - min == 1) {
    return rootNode;
  } else if (max - min > 1) {
    int mid = (max + min) >> 1;
    rootNode->leftChild = createTreeTopDown(pts, min, mid);
    rootNode->rightChild = createTreeTopDown(pts, mid, max);
  }
  return rootNode;
}

// create segment tree only keep endpoints
TreeNode *createSegmentTree(const vector<int> &pts) {
  const int leaf_node_size = pts.size() - 1;
  TreeNode *tree = createTreeTopDown(pts, 0, leaf_node_size);
  return tree;
}

void insertSegment(const int l, const int r, TreeNode *tree, const int idx) {
  if (tree == nullptr) {
    return;
  }
  if ((tree->rid - tree->lid == 1) || (tree->lid == l && tree->rid == r)) {
    tree->data.emplace_back(idx);
    return;
  }

  int mid = tree->rightChild->lid;

  if (r <= mid) {
    insertSegment(l, r, tree->leftChild, idx);
  } else if (l >= mid) {
    insertSegment(l, r, tree->rightChild, idx);
  } else if (l != mid) {
    insertSegment(l, mid, tree->leftChild, idx);
    insertSegment(mid, r, tree->rightChild, idx);
  }
}

vector<int> searchCW(const int pos, TreeNode *tree) {
  vector<int> res;
  // if (tree->leftChild == nullptr || tree->rightChild == nullptr) {
  //   res.insert(res.end(), tree->data.begin(), tree->data.end());
  //   return res;
  // }
  if (tree->lid <= pos && pos <= tree->rid) {
    res.insert(res.end(), tree->data.begin(), tree->data.end());
  }
  int mid = 0;
  if (tree->rightChild != nullptr)
    mid = tree->rightChild->lid;
  else if (tree->leftChild != nullptr)
    mid = tree->leftChild->rid;
  if (pos <= mid && tree->leftChild != nullptr) {
    auto tres = searchCW(pos, tree->leftChild);
    res.insert(res.end(), tres.begin(), tres.end());
  }
  if (pos >= mid && tree->rightChild != nullptr) {
    auto tres = searchCW(pos, tree->rightChild);
    res.insert(res.end(), tres.begin(), tres.end());
  }

  return res;
}

class DeltaSTIndex : public BaseIndex {
 private:
  // vector<pair<vector<DeltaNode>, vector<DeltaNode>>> indexed_arr;
  vector<pair<TreeNode *, TreeNode *>> indexed_tree;

 public:
  DeltaSTIndex() {
    window_count = 0;
    nodes_amount = 0;
    sort_time = 0;
  }
  void collectWindows(vector<DeltaNode> &left_pairs,
                      vector<DeltaNode> &right_pairs, vector<int> &l_arr,
                      vector<int> &r_arr, vector<int> &visited_left,
                      vector<int> &visited_right, const vector<int> &idx_sorted,
                      const int k_smallest, const int query_pos);
  void indexNodes(const int query_pos, const vector<vector<float>> &data_nodes,
                  const vector<int> &idx_sorted, const int k_smallest,
                  const int l_bound, const int r_bound, TreeNode **,
                  TreeNode **);

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

  void rebuildRangeKNNG(const int left, const int right,
                        vector<vector<int>> &knng);

  void clear() {
    indexed_tree.clear();
    window_count = 0;
    vector<pair<TreeNode *, TreeNode *>>().swap(indexed_tree);
  }
  ~DeltaSTIndex() { vector<pair<TreeNode *, TreeNode *>>().swap(indexed_tree); }
};

void DeltaSTIndex::collectWindows(vector<DeltaNode> &left_pairs,
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
      // else if (l_arr.size() == k_smallest + 1) {
      // }
      left_pairs.emplace_back(a_node);
      visited_left.emplace_back(idx);

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
      // else if (l_arr.size() == k_smallest + 1) {
      // }
      right_pairs.emplace_back(a_node);
      visited_right.emplace_back(idx);
    }
  }
}

// index sorted vector
// idx_sorted: index vector sorted by distance to query point
void DeltaSTIndex::indexNodes(const int query_pos,
                              const vector<vector<float>> &data_nodes,
                              const vector<int> &idx_sorted,
                              const int k_smallest, const int l_bound,
                              const int r_bound, TreeNode **left_root,
                              TreeNode **right_root) {
  vector<int> l_arr;
  vector<int> r_arr;
  vector<DeltaNode> left_pairs, right_pairs;

  // store visited indexes to generate tree.
  vector<int> visited_left;
  vector<int> visited_right;

  l_arr.emplace_back(-1);
  r_arr.emplace_back(r_bound + 1);

  collectWindows(left_pairs, right_pairs, l_arr, r_arr, visited_left,
                 visited_right, idx_sorted, k_smallest, query_pos);
  // cout << "first collect" << endl;

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

  visited_left.emplace_back(-1);
  visited_right.emplace_back(query_pos);
  visited_left.emplace_back(query_pos);
  visited_right.emplace_back(INT_MAX);
  visited_right.emplace_back(r_bound + 1);

  std::sort(visited_right.begin(), visited_right.end());
  visited_right.erase(std::unique(visited_right.begin(), visited_right.end()),
                      visited_right.end());
  std::sort(visited_left.begin(), visited_left.end());
  visited_left.erase(std::unique(visited_left.begin(), visited_left.end()),
                     visited_left.end());
  // build segment tree
  *left_root = createSegmentTree(visited_left);
  *right_root = createSegmentTree(visited_right);

  for (auto &ele : left_pairs) {
    if (ele.add == INT_MAX)
      insertSegment(ele.sub, query_pos, *left_root, query_pos);
    else
      insertSegment(ele.sub, ele.add, *left_root, ele.add);

    // if (query_pos == 290) {
    //   cout << ele.sub << "," << ele.add << endl;
    // }
  }

  for (auto &ele : right_pairs) {
    insertSegment(ele.add, ele.sub, *right_root, ele.add);
  }

  nodes_amount += (2 * pow(2, floor(log2(visited_left.size())) + 1) - 1);
  nodes_amount += (2 * pow(2, floor(log2(visited_right.size())) + 1) - 1);

  return;
}

void DeltaSTIndex::buildIndex(const vector<vector<float>> &nodes,
                              const int initial_build_k = 20,
                              const int index_k = 10) {
  auto distance = [](vector<float> a, vector<float> b) -> double {
    float ans = 0.0;
    for (int i = 0; i < a.size(); ++i) {
      ans += (a[i] - b[i]) * (a[i] - b[i]);
    }
    return ans;
  };
  indexed_tree.clear();
  indexed_tree.resize(nodes.size());
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

    IN.findAllNeighbors(10, index_k);

    auto idx_sorted = IN.sorted_idxes;
    hop_amount += IN.getHop();
    gettimeofday(&tt4, NULL);
    AccumulateTime(tt3, tt4, sort_time);

    TreeNode *left_root;
    TreeNode *right_root;
    indexNodes(idx, nodes, idx_sorted, K2_, 0, nodes.size() - 1, &left_root,
               &right_root);
    indexed_tree[idx].first = left_root;
    indexed_tree[idx].second = right_root;

    counter++;
  }
  delete kgraph_index;
  cout << "Sort Time: " << sort_time << endl;
  logTime(tt1, tt2, "Index Time");
  cout << "average hop: " << (float)hop_amount / nodes.size() << endl;
  cout << "index amount: " << window_count << endl;
}

vector<int> DeltaSTIndex::indexNeighborSearch(
    const vector<vector<float>> &data_nodes, const vector<float> &query,
    const int l_bound, const int r_bound, const int K_neigbhors, const int ef) {
  vector<int> temp;
  return indexNeighborSearch(data_nodes, query, l_bound, r_bound, K_neigbhors,
                             ef, temp);
}

// on the fly search cw
// TODO: Verify search path result
vector<int> DeltaSTIndex::indexNeighborSearch(
    const vector<vector<float>> &data_nodes, const vector<float> &query,
    const int l_bound, const int r_bound, const int K_neigbhors, const int ef,
    const vector<int> &enter_points) {
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
    gettimeofday(&tt1, NULL);
    auto left_neighbors =
        searchCW(l_bound, indexed_tree[current_node_id].first);
    auto right_neighbors =
        searchCW(r_bound, indexed_tree[current_node_id].second);

    // print_set(left_neighbors);
    // print_set(right_neighbors);

    vector<int> current_neighbors(left_neighbors);
    current_neighbors.insert(current_neighbors.end(), right_neighbors.begin(),
                             right_neighbors.end());
    // gettimeofday(&tt2, NULL);
    // AccumulateTime(tt1, tt2, tree_search_time);

    if (current_neighbors.size() == 0) {
      cout << "current node: " << current_node_pair.second << "   ";
      cout << "size 0 neighbors" << endl;
      continue;
    }
    // print_set(current_neighbors);
    gettimeofday(&tt1, NULL);
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

double DeltaSTIndex::calSize() {
  // double unit_size = sizeof(*this->indexed_tree.front().first);
  // cout << "unit size: " << unit_size << endl;
  return 24 * nodes_amount;
}

pair<vector<int>, vector<int>> DeltaSTIndex::checkIndexQuality(
    const int idx, const int lbound, const int rbound) {
  vector<int> lres, rres;
  lres = searchCW(lbound, indexed_tree[idx].first);
  rres = searchCW(rbound, indexed_tree[idx].second);
  // vector<TreeNode *> q;
  // q.push_back(indexed_tree[idx].second);
  // while (!q.empty()) {
  //   auto ptr = q.front();
  //   if (ptr->leftChild != nullptr) {
  //     q.push_back(ptr->leftChild);
  //   }
  //   if (ptr->rightChild != nullptr) {
  //     q.push_back(ptr->rightChild);
  //   }
  //   cout << "[" << ptr->lid << "," << ptr->rid << "] : ";
  //   print_set(ptr->data);
  //   q.erase(q.begin());
  // }
  sort(lres.begin(), lres.end());
  sort(rres.begin(), rres.end());
  print_set(lres);
  print_set(rres);
  return make_pair(lres, rres);
}

void DeltaSTIndex::rebuildRangeKNNG(const int left, const int right,
                                    vector<vector<int>> &knng) {
  knng.clear();
  knng.resize(right - left + 1);
#pragma omp parallel for
  for (unsigned n = left; n <= right; n++) {
    vector<int> lres, rres;
    lres = searchCW(left, indexed_tree[n].first);
    rres = searchCW(right, indexed_tree[n].second);
    auto &knn = knng[n - left];
    knn.insert(knn.end(), lres.begin(), lres.end());
    knn.insert(knn.end(), rres.begin(), rres.end());
  }
}

}  // namespace deltaindex