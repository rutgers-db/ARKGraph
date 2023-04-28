/**
**
*/

#include <algorithm>
#include <iostream>
#include <numeric>
#include <vector>
#include <sys/time.h>

#include "lib_util/utils.hpp"
#include "utils.h"

using namespace std;

void greedyNearest(const int query_pos, const vector<vector<float>> &dpts, const int k_smallest, const int l_bound, const int r_bound)
{

  vector<float> dist_arr;
  for (size_t i = l_bound; i < r_bound; i++)
  {
    dist_arr.emplace_back(EuclideanDistance(dpts[query_pos], dpts[i], 0, 0));
  }
  vector<int> sorted_idxes = sort_indexes(dist_arr);

  // skip the point itself
  if (sorted_idxes[0] == query_pos)
  {
    sorted_idxes.erase(sorted_idxes.begin());
  }
  sorted_idxes.resize(k_smallest);
  // print_set(sorted_idxes);
}

void rangeGreedy(const vector<vector<float>> &nodes, const int k_smallest, const int l_bound, const int r_bound)
{

  for (size_t i = l_bound; i < r_bound; i++)
  {
    greedyNearest(i, nodes, k_smallest, l_bound, r_bound);
  }
}


int main()
{
  timeval t1, t2;
  string path = "../../data/siftsmall_base.fvecs";
  vector<vector<float>> nodes = pqdescent::ReadTopN(path, "fvecs", 600);

  cout << nodes.size() << " " << nodes[0].size() << endl;

  gettimeofday(&t1, NULL);
  // greedyNearest(0, nodes, 10, 0, 599);
  rangeGreedy(nodes, 10, 0, 599);
  logTime(t1, t2, "cost time");
  return 0;
}