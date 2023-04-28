/**
 * Compact Array for split method
 *
 * Author: Chaoji Zuo
 * Email:  chaoji.zuo@rutgers.edu
 * Date:   Feb 24, 2022
 */

#pragma once

#include "cw.h"
#include <iostream>
#include <unordered_map>
using std::cout;
using std::string;
using std::vector;

namespace compactArr {

class CompactArray {
public:
  vector<pair<int, CompactWindow *>> vecs;
  CompactArray(){};
  bool insertCW(const int, CompactWindow *);
  CompactWindow *searchCW(const int pos, bool left);
  void clear() { vector<pair<int, CompactWindow *>>().swap(vecs); }
  ~CompactArray() { vector<pair<int, CompactWindow *>>().swap(vecs); }
};

bool CompactArray::insertCW(const int pos, CompactWindow *cw) {
  if (vecs.size() == 0) {
    vecs.emplace_back(pos, cw);
    return true;
  }
  if (vecs.back().first < pos) {
    vecs.emplace_back(pos, cw);
    return true;
  } else
    return false;
}

struct compare {
  bool operator()(const pair<int, CompactWindow *> &value, const int &key) {
    return (value.first < key);
  }
  bool operator()(const int &key, const pair<int, CompactWindow *> &value) {
    return (key < value.first);
  }
};

CompactWindow *CompactArray::searchCW(const int pos, bool left = false) {

  auto res = lower_bound(vecs.begin(), vecs.end(), pos, compare());
  if (left && res != vecs.begin())
    res--;
  if (res != vecs.end()) {
    return res->second;
  } else {
    return nullptr;
  }
}

} // namespace compactArr
