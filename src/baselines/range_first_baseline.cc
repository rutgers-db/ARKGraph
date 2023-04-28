/**
 * baseline #1, calculate nearest range query first, then search among the range.
 *
 * Author: Chaoji Zuo
 * Date:   Nov 13, 2021
 * Email:  chaoji.zuo@rutgers.edu
 */

#include <algorithm>
#include <iostream>
#include <numeric>
#include <vector>

#include "lib_pq/pq.hpp"
#include "../utils/utils.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;


int main(int argc, char **argv)
{

  string path = "../../data/siftsmall_base.fvecs";
  vector<vector<float>> nodes = pqdescent::ReadTopN(path, "fvecs", 1000);

  int l_bound = 190;
  int r_bound = 770;
  int query_point = 500;
  int acutal_K = 10;
  auto nodes_in_range = vector<vector<float>>(nodes.begin() + l_bound, nodes.begin() + r_bound);

  timeval t1, t2;
  gettimeofday(&t1, NULL);
  // vd.setK(10);

  // The number of subspace
  int m = 16;
  // The number of codewords for each subspace
  int ks = 256;

  // Initialize PQ
  pqdescent::PQ pq(pqdescent::PQ::Learn(nodes_in_range, m, ks));
  pqdescent::UcharVecs encodedNodes = pq.Encode(nodes_in_range);
  auto codewords = pq.GetCodewords();
  logTime(t1, t2, "index construction time");
  gettimeofday(&t1, NULL);

  auto res = pq.AD(pq.DTable(nodes_in_range[query_point - l_bound]), encodedNodes);
  auto idx_sorted = sort_indexes(res);
  // skip element itself
  idx_sorted.erase(idx_sorted.begin());
  idx_sorted.resize(acutal_K);

  print_set(idx_sorted);
  for (auto &ele : idx_sorted)
  {
    ele += l_bound;
  }
  sort(idx_sorted.begin(), idx_sorted.end());
  logTime(t1, t2, "search time");
  print_set(idx_sorted);

  return 0;
}