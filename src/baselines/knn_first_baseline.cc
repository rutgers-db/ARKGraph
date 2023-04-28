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


#include "../utils/utils.h"
#include "hnswlib/hnswlib.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;


int main(int argc, char **argv)
{

  string path = "../data/siftsmall_base.fvecs";
  vector<vector<float>> nodes = pqdescent::ReadTopN(path, "fvecs", 1000);

  int l_bound = 190;
  int r_bound = 770;
  int query_point = 699;
  timeval t1, t2;
  gettimeofday(&t1, NULL);
  int efConstruction = 40;
  int M = 16;

  hnswlib::L2Space space(nodes.front().size());
  hnswlib::AlgorithmInterface<float> *alg_hnsw = new hnswlib::HierarchicalNSW<float>(&space, 2 * nodes.size(), M, efConstruction);
  // hnswlib::AlgorithmInterface<float> *alg_hnsw = new hnswlib::BruteforceSearch<float>(&space, 2 * nodes.size());
  for (size_t i = 0; i < nodes.size(); ++i)
  {
    alg_hnsw->addPoint(nodes[i].data(), i);
  }
  logTime(t1, t2, "graph construction time");

  gettimeofday(&t1, NULL);
  // K to run nndescent
  int init_K = 100;
  int actual_K = 10;
  vector<vector<int>> result_in_range(r_bound - l_bound, vector<int>());
  for (size_t i = l_bound; i < r_bound; i++)
  {
    const void *p = nodes[i].data();

    auto res = alg_hnsw->searchKnnCloserFirst(p, init_K + 1); // + 1 to avoid itself

    for (size_t j = 0; j < init_K; j++)
    {
      if (result_in_range[i - l_bound].size() >= actual_K)
      {
        // sort(result_in_range[ i - l_bound].begin(), result_in_range[i - l_bound].end());
        // print_set(result_in_range[i - l_bound]);
        break;
      }
      int val = res[j].second;

      if (val < r_bound && val > l_bound)
      {
        result_in_range[i - l_bound].emplace_back(val);
      }
    }
  }

  cout << result_in_range.size() << endl;
  print_set(result_in_range[query_point - l_bound]);
  logTime(t1, t2, "search time");
  sort(result_in_range[query_point - l_bound].begin(), result_in_range[query_point - l_bound].end());
  print_set(result_in_range[query_point - l_bound]);

  // for test just one element
  // const void *p = nodes[query_point].data();
  // auto res = alg_hnsw->searchKnnCloserFirst(p, init_K + 1); // + 1 to avoid itself
  // vector<int> result_in_range;
  // for (size_t j = 1; j < init_K + 1; j++)
  // {
  //   if (result_in_range.size() >= actual_K)
  //   {
  //     break;
  //   }
  //   int val = res[j].second;

  //   if (val < r_bound && val > l_bound)
  //   {
  //     result_in_range.emplace_back(val);
  //   }
  // }
  // print_set(result_in_range);
  // logTime(t1, t2, "cost time");
  // sort(result_in_range.begin(), result_in_range.end());
  // print_set(result_in_range);

  return 0;
}