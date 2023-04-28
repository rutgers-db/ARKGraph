#pragma once

#include <ctime>
#include <fstream>
#include <functional>
#include <iostream>
#include <vector>
// #include <omp.h>
#include <assert.h>
#include <sys/time.h>

#include <algorithm>
#include <iostream>
#include <numeric>
#include <queue>
#include <random>
#include <unordered_map>
#include <unordered_set>

#ifdef __linux__
#include "sys/sysinfo.h"
#include "sys/types.h"
#elif __APPLE__
#include <mach/mach_host.h>
#include <mach/mach_init.h>
#include <mach/mach_types.h>
#include <mach/vm_statistics.h>
#endif

using std::cout;
using std::endl;
using std::getline;
using std::ifstream;
using std::ios;
using std::string;
using std::vector;

void testUTIL2();

float EuclideanDistance(const vector<float> &lhs, const vector<float> &rhs,
                        const int &startDim, int lensDim);

float EuclideanDistance(const vector<float> &lhs, const vector<float> &rhs);

void AccumulateTime(timeval &t2, timeval &t1, double &val_time);
void CountTime(timeval &t1, timeval &t2, double &val_time);

// the same to sort_indexes
template <typename T>
std::vector<std::size_t> sort_permutation(const std::vector<T> &vec) {
  std::vector<std::size_t> p(vec.size());
  std::iota(p.begin(), p.end(), 0);
  std::sort(p.begin(), p.end(),
            [&](std::size_t i, std::size_t j) { return vec[i] < vec[j]; });
  return p;
}

// apply permutation
template <typename T>
void apply_permutation_in_place(std::vector<T> &vec,
                                const std::vector<std::size_t> &p) {
  std::vector<bool> done(vec.size());
  for (std::size_t i = 0; i < vec.size(); ++i) {
    if (done[i]) {
      continue;
    }
    done[i] = true;
    std::size_t prev_j = i;
    std::size_t j = p[i];
    while (i != j) {
      std::swap(vec[prev_j], vec[j]);
      done[j] = true;
      prev_j = j;
      j = p[j];
    }
  }
}

template <typename T>
vector<int> sort_indexes(const vector<T> &v) {
  // initialize original index locations
  vector<int> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values
  stable_sort(idx.begin(), idx.end(),
              [&v](size_t i1, size_t i2) { return v[i1] < v[i2]; });

  return idx;
}

template <typename T>
vector<int> sort_indexes(const vector<T> &v, const int begin_bias,
                         const int end_bias) {
  // initialize original index locations
  vector<int> idx(end_bias - begin_bias);
  iota(idx.begin() + begin_bias, idx.begin() + end_bias, 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values
  stable_sort(idx.begin() + begin_bias, idx.begin() + end_bias,
              [&v](size_t i1, size_t i2) { return v[i1] < v[i2]; });

  return idx;
}

template <typename T>
void print_set(const vector<T> &v) {
  if (v.size() == 0) {
    cout << "ERROR: EMPTY VECTOR!" << endl;
    return;
  }
  cout << "vertex in set: {";
  for (size_t i = 0; i < v.size() - 1; i++) {
    cout << v[i] << ", ";
  }
  cout << v.back() << "}" << endl;
}

void Split(std::string &s, std::string &delim, std::vector<std::string> *ret);
void ReadMatFromTxt(const std::string &path,
                    std::vector<std::vector<float>> &data,
                    const int length_limit);

void ReadMatFromTxtTwitter(const std::string &path,
                           std::vector<std::vector<float>> &data,
                           const int length_limit);
void ReadMatFromTsv(const std::string &path,
                    std::vector<std::vector<float>> &data,
                    const int length_limit);

void ReadDataWrapper(vector<vector<float>> &raw_data, vector<int> &search_keys,
                     const string &dataset, string &dataset_path,
                     const int item_num);
void ReadDataWrapper(const string &dataset, string &dataset_path,
                     vector<vector<float>> &raw_data, const int data_size,
                     string &query_path, vector<vector<float>> &querys,
                     const int query_size, vector<int> &search_keys);
void SynthesizeQuerys(const vector<vector<float>> &nodes,
                      vector<vector<float>> &querys, const int query_num);
void logTime(timeval &begin, timeval &end, const string &log);

double countPrecision(const vector<int> &truth, const vector<int> &pred);
double countApproximationRatio(const vector<vector<float>> &raw_data,
                               const vector<int> &truth,
                               const vector<int> &pred,
                               const vector<float> &query);

int YT8M2Int(const string id);
void ReadMatFromTsvYT8M(const string &path, vector<vector<float>> &data,
                        vector<int> &search_keys, const int length_limit);

void print_memory();
void record_memory(long long &);
#define _INT_MAX 2147483640

vector<int> greedyNearest(const vector<vector<float>> &dpts,
                          const vector<float> query, const int k_smallest);

void ReadGroundtruthQuery(std::vector<std::vector<int>> &gt,
                          std::vector<std::pair<int, int>> &query_ranges,
                          std::vector<int> &query_ids, std::string gt_path);

void fvecs2csv(const string &output_path, const vector<vector<float>> &nodes);
void evaluateKNNG(const vector<vector<int>> &gt,
                  const vector<vector<int>> &knng, const int K, double &recall,
                  double &precision);

namespace pqdescent {
using std::vector;
// Iterative reader class for reading .bvecs or .fvecs files.
// The next vector (std::vector<float>) is read by Next() function.
//
// Usage:
//   ItrReader reader("data.fvecs", "fvecs");
//   while(!reader.IsEnd()){
//     std::vector<float> vec = reader.Next();
//     /* some stuff for vec */
//   }
//
// Optional wrapper interface:
//   int top_n = 100;
//   std::vector<std::vector<float> > vecs = ReadTopN("data.fvecs", "fvecs",
//   top_n);

// Interface (abstract basic class) of iterative reader
class I_ItrReader {
 public:
  virtual ~I_ItrReader() {}
  virtual bool IsEnd() = 0;
  virtual std::vector<float> Next() = 0;
};

// Iterative reader for fvec file
class FvecsItrReader : I_ItrReader {
 public:
  FvecsItrReader(std::string filename);
  bool IsEnd();
  std::vector<float> Next();

 private:
  FvecsItrReader();  // prohibit default construct
  std::ifstream ifs;
  std::vector<float> vec;  // store the next vec
  bool eof_flag;
};

// Iterative reader for bvec file
class BvecsItrReader : I_ItrReader {
 public:
  BvecsItrReader(std::string filename);
  bool IsEnd();
  std::vector<float> Next();  // Read bvec, but return vec<float>
 private:
  BvecsItrReader();  // prohibit default construct
  std::ifstream ifs;
  std::vector<float> vec;  // store the next vec
  bool eof_flag;
};

// Proxy class
class ItrReader {
 public:
  // ext must be "fvecs" or "bvecs"
  ItrReader(std::string filename, std::string ext);
  ~ItrReader();

  bool IsEnd();
  std::vector<float> Next();

 private:
  ItrReader();
  I_ItrReader *m_reader;
};

// Wrapper. Read top-N vectors
// If top_n = -1, then read all vectors
std::vector<std::vector<float>> ReadTopN(std::string filename, std::string ext,
                                         int top_n = -1);

}  // namespace pqdescent