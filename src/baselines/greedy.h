/**
**
*/
#pragma once

#include <algorithm>
#include <iostream>
#include <numeric>
#include <sys/time.h>
#include <vector>

#include "../utils/utils.h"

#define _INT_MAX 2147483640
using namespace std;

vector<int> greedyNearest(const vector<vector<float>> &dpts,
                          const vector<float> query, const int l_bound,
                          const int r_bound, const int k_smallest) {

  priority_queue<pair<float, int>> top_candidates;
  float lower_bound = _INT_MAX;
  for (size_t i = l_bound; i <= r_bound; i++) {
    float dist = EuclideanDistance(query, dpts[i]);
    if (top_candidates.size() < k_smallest || dist < lower_bound) {
      top_candidates.push(make_pair(dist, i));
      if (top_candidates.size() > k_smallest) {
        top_candidates.pop();
      }

      lower_bound = top_candidates.top().first;
    }
  }
  vector<int> res;
  while (!top_candidates.empty()) {
    res.emplace_back(top_candidates.top().second);
    top_candidates.pop();
  }
  return res;
}
