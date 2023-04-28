#include "utils.h"

float EuclideanDistance(const vector<float> &lhs, const vector<float> &rhs,
                        const int &startDim, int lensDim) {
  float ans = 0.0;
  if (lensDim == 0) {
    lensDim = lhs.size();
  }

  for (int i = startDim; i < startDim + lensDim; ++i) {
    ans += (lhs[i] - rhs[i]) * (lhs[i] - rhs[i]);
  }
  return ans;
}

void testUTIL2() { cout << "hello" << endl; }

float EuclideanDistance(const vector<float> &lhs, const vector<float> &rhs) {
  return EuclideanDistance(lhs, rhs, 0, 0);
}

// t1:begin, t2:end
void AccumulateTime(timeval &t1, timeval &t2, double &val_time) {
  val_time += (t2.tv_sec - t1.tv_sec +
               (t2.tv_usec - t1.tv_usec) * 1.0 / CLOCKS_PER_SEC);
}

void CountTime(timeval &t1, timeval &t2, double &val_time) {
  val_time = 0;
  val_time += (t2.tv_sec - t1.tv_sec +
               (t2.tv_usec - t1.tv_usec) * 1.0 / CLOCKS_PER_SEC);
}

void Split(std::string &s, std::string &delim, std::vector<std::string> *ret) {
  size_t last = 0;
  size_t index = s.find_first_of(delim, last);
  while (index != std::string::npos) {
    ret->push_back(s.substr(last, index - last));
    last = index + 1;
    index = s.find_first_of(delim, last);
  }
  if (index - last > 0) {
    ret->push_back(s.substr(last, index - last));
  }
}

void ReadDataWrapper(vector<vector<float>> &raw_data, vector<int> &search_keys,
                     const string &dataset, string &dataset_path,
                     const int item_num) {
  raw_data.clear();
  if (dataset == "glove") {
    ReadMatFromTxtTwitter(dataset_path, raw_data, item_num);
  } else if (dataset == "ml25m") {
    ReadMatFromTxt(dataset_path, raw_data, item_num);
  } else if (dataset == "sift") {
    raw_data = pqdescent::ReadTopN(dataset_path, "bvecs", item_num);
  } else if (dataset == "biggraph") {
    ReadMatFromTsv(dataset_path, raw_data, item_num);
  } else if (dataset == "local") {
    if (dataset_path == "")
      dataset_path =
          "../data/siftsmall_base.fvecs";
    cout << dataset_path << endl;
    raw_data = pqdescent::ReadTopN(dataset_path, "fvecs", item_num);
  } else if (dataset == "deep1b") {
    raw_data = pqdescent::ReadTopN(dataset_path, "fvecs", item_num);
  } else if (dataset == "deep10m") {
    raw_data = pqdescent::ReadTopN(dataset_path, "fvecs", item_num);
  } else if (dataset == "yt8m") {
    ReadMatFromTsvYT8M(dataset_path, raw_data, search_keys, item_num);
  } else {
    std::cerr << "Wrong Datset!" << endl;
    assert(false);
  }
}

void SynthesizeQuerys(const vector<vector<float>> &nodes,
                      vector<vector<float>> &querys, const int query_num) {
  int dim = nodes.front().size();
  std::default_random_engine e;
  std::uniform_int_distribution<int> u(0, nodes.size() - 1);
  querys.clear();
  querys.resize(query_num);

  for (unsigned n = 0; n < query_num; n++) {
    for (unsigned i = 0; i < dim; i++) {
      int select_idx = u(e);
      querys[n].emplace_back(nodes[select_idx][i]);
    }
  }
}

// load txt matrix data
void ReadMatFromTxt(const string &path, vector<vector<float>> &data,
                    const int length_limit = -1) {
  ifstream infile;
  string bline;
  string delim = " ";
  int numCols = 0;
  infile.open(path, ios::in);
  if (getline(infile, bline, '\n')) {
    vector<string> ret;
    Split(bline, delim, &ret);
    numCols = ret.size();
  }
  infile.close();
  // cout << "Reading " << path << " ..." << endl;
  // cout << "# of columns: " << numCols << endl;

  int counter = 0;
  if (length_limit == -1) counter = -9999999;
  // TODO: read sparse matrix
  infile.open(path, ios::in);
  while (getline(infile, bline, '\n')) {
    if (counter >= length_limit) break;
    counter++;

    vector<string> ret;
    Split(bline, delim, &ret);
    vector<float> arow(numCols);
    assert(ret.size() == numCols);
    for (int i = 0; i < ret.size(); i++) {
      arow[i] = static_cast<float>(stod(ret[i]));
    }
    data.emplace_back(arow);
  }
  infile.close();
  // cout << "# of rows: " << data.size() << endl;
}

void ReadMatFromTxtTwitter(const string &path, vector<vector<float>> &data,
                           const int length_limit = -1) {
  ifstream infile;
  string bline;
  string delim = " ";
  int numCols = 0;
  infile.open(path, ios::in);
  if (getline(infile, bline, '\n')) {
    vector<string> ret;
    Split(bline, delim, &ret);
    numCols = ret.size() - 1;
  }
  infile.close();
  cout << "Reading " << path << " ..." << endl;

  cout << "# of columns: " << numCols << endl;

  int counter = 0;
  if (length_limit == -1) counter = -9999999;
  // TODO: read sparse matrix
  infile.open(path, ios::in);
  while (getline(infile, bline, '\n')) {
    if (counter >= length_limit) break;
    counter++;

    vector<string> ret;
    Split(bline, delim, &ret);
    vector<float> arow(numCols);
    assert(ret.size() == numCols + 1);
    for (int i = 1; i < ret.size(); i++) {
      arow[i - 1] = static_cast<float>(stod(ret[i]));
    }
    data.emplace_back(arow);
  }
  infile.close();
  cout << "# of rows: " << data.size() << endl;
}

void ReadMatFromTsv(const string &path, vector<vector<float>> &data,
                    const int length_limit = -1) {
  ifstream infile;
  string bline;
  string delim = "\t";
  int numCols = 0;
  infile.open(path, ios::in);
  getline(infile, bline, '\n');
  if (getline(infile, bline, '\n')) {
    vector<string> ret;
    Split(bline, delim, &ret);
    numCols = ret.size();
  }
  infile.close();
  cout << "Reading " << path << " ..." << endl;
  cout << "# of columns: " << numCols << endl;

  int counter = 0;
  if (length_limit == -1) counter = -9999999;
  infile.open(path, ios::in);
  // skip the first line
  getline(infile, bline, '\n');
  while (getline(infile, bline, '\n')) {
    if (counter >= length_limit) break;
    counter++;

    vector<string> ret;
    Split(bline, delim, &ret);
    vector<float> arow(numCols - 1);
    assert(ret.size() == numCols);
    for (int i = 0; i < ret.size() - 1; i++) {
      arow[i] = static_cast<float>(stod(ret[i + 1]));
    }
    data.emplace_back(arow);
  }
  infile.close();
  cout << "# of rows: " << data.size() << endl;
}

int YT8M2Int(const string id) {
  int res = 0;
  for (size_t i = 0; i < 4; i++) {
    res *= 100;
    res += (int)id[i] - 38;
  }
  return res;
}

void ReadMatFromTsvYT8M(const string &path, vector<vector<float>> &data,
                        vector<int> &search_keys, const int length_limit) {
  ifstream infile;
  string bline;
  string delim = ",";
  int numCols = 0;
  infile.open(path, ios::in);
  getline(infile, bline, '\n');
  if (getline(infile, bline, '\n')) {
    vector<string> ret;
    Split(bline, delim, &ret);
    numCols = ret.size();
  }
  infile.close();
  cout << "Reading " << path << " ..." << endl;
  cout << "# of columns: " << numCols << endl;

  int counter = 0;
  if (length_limit == -1) counter = -9999999;
  infile.open(path, ios::in);
  string delim_embed = " ";

  while (getline(infile, bline, '\n')) {
    if (counter >= length_limit) break;
    counter++;

    vector<string> ret;
    Split(bline, delim, &ret);
    assert(ret.size() == numCols);

    // str 'id' to int 'id'
    // int one_search_key = YT8M2Int(ret[0]);
    int one_search_key = (int)stod(ret[1]);

    // add embedding
    string embedding_str = ret[2];
    vector<string> embedding_vec;
    vector<float> arow(1024);
    Split(embedding_str, delim_embed, &embedding_vec);
    assert(embedding_vec.size() == 1024);
    for (int i = 0; i < embedding_vec.size() - 1; i++) {
      arow[i] = static_cast<float>(stod(embedding_vec[i + 1]));
    }
    search_keys.emplace_back(one_search_key);
    data.emplace_back(arow);
  }
  infile.close();
  cout << "# of rows: " << data.size() << endl;
}

void logTime(timeval &begin, timeval &end, const string &log) {
  gettimeofday(&end, NULL);
  fprintf(stdout, ("# " + log + ": %.7fs\n").c_str(),
          end.tv_sec - begin.tv_sec +
              (end.tv_usec - begin.tv_usec) * 1.0 / CLOCKS_PER_SEC);
};

double countPrecision(const vector<int> &truth, const vector<int> &pred) {
  double num_right = 0;
  for (auto one : truth) {
    if (find(pred.begin(), pred.end(), one) != pred.end()) {
      num_right += 1;
    }
  }
  return num_right / truth.size();
}

double countApproximationRatio(const vector<vector<float>> &raw_data,
                               const vector<int> &truth,
                               const vector<int> &pred,
                               const vector<float> &query) {
  if (pred.size() == 0) {
    return 0;
  }
  vector<float> truth_dist;
  vector<float> pred_dist;
  for (auto vec : truth) {
    truth_dist.emplace_back(EuclideanDistance(query, raw_data[vec]));
  }
  for (auto vec : pred) {
    pred_dist.emplace_back(EuclideanDistance(query, raw_data[vec]));
  }
  auto max_truth = *max_element(truth_dist.begin(), truth_dist.end());
  auto max_pred = *max_element(pred_dist.begin(), pred_dist.end());
  if (pred.size() < truth.size()) {
    nth_element(truth_dist.begin(), truth_dist.begin() + pred.size() - 1,
                truth_dist.end());
    max_truth = truth_dist[pred.size() - 1];
  }
  if (max_truth != 0) return max_pred / max_truth;
  cout << "ERROR: empty pred!" << endl;
  return -1;
}

void print_memory() {
#ifdef __linux__
  struct sysinfo memInfo;

  sysinfo(&memInfo);
  // long long totalVirtualMem = memInfo.totalram;
  // // Add other values in next statement to avoid int overflow on right hand
  // // side...
  // totalVirtualMem += memInfo.totalswap;
  // totalVirtualMem *= memInfo.mem_unit;

  // long long virtualMemUsed = memInfo.totalram - memInfo.freeram;
  // // Add other values in next statement to avoid int overflow on right hand
  // // side...
  // virtualMemUsed += memInfo.totalswap - memInfo.freeswap;
  // virtualMemUsed *= memInfo.mem_unit;
  // cout << "Total Virtual Memory: " << totalVirtualMem << endl;
  // cout << "Used Virtual Memory: " << virtualMemUsed << endl;

  long long totalPhysMem = memInfo.totalram;
  // Multiply in next statement to avoid int overflow on right hand side...
  totalPhysMem *= memInfo.mem_unit;

  long long physMemUsed = memInfo.totalram - memInfo.freeram;
  // Multiply in next statement to avoid int overflow on right hand side...
  physMemUsed *= memInfo.mem_unit;

  // cout << "Total Physical Memory: " << totalPhysMem << endl;
  cout << "Used Physical Memory: " << physMemUsed << endl;
#elif __APPLE__
  vm_size_t page_size;
  mach_port_t mach_port;
  mach_msg_type_number_t count;
  vm_statistics64_data_t vm_stats;

  mach_port = mach_host_self();
  count = sizeof(vm_stats) / sizeof(natural_t);
  if (KERN_SUCCESS == host_page_size(mach_port, &page_size) &&
      KERN_SUCCESS == host_statistics64(mach_port, HOST_VM_INFO,
                                        (host_info64_t)&vm_stats, &count)) {
    long long free_memory = (int64_t)vm_stats.free_count * (int64_t)page_size;

    long long used_memory =
        ((int64_t)vm_stats.active_count + (int64_t)vm_stats.inactive_count +
         (int64_t)vm_stats.wire_count) *
        (int64_t)page_size;
    printf("free memory: %lld\nused memory: %lld\n", free_memory, used_memory);
  }
#endif
}

void record_memory(long long &memory) {
#ifdef __linux__
  struct sysinfo memInfo;
  sysinfo(&memInfo);
  long long physMemUsed = memInfo.totalram - memInfo.freeram;
  physMemUsed *= memInfo.mem_unit;
  memory = physMemUsed;
#elif __APPLE__
  vm_size_t page_size;
  mach_port_t mach_port;
  mach_msg_type_number_t count;
  vm_statistics64_data_t vm_stats;

  mach_port = mach_host_self();
  count = sizeof(vm_stats) / sizeof(natural_t);
  if (KERN_SUCCESS == host_page_size(mach_port, &page_size) &&
      KERN_SUCCESS == host_statistics64(mach_port, HOST_VM_INFO,
                                        (host_info64_t)&vm_stats, &count)) {
    memory = ((int64_t)vm_stats.active_count +
              (int64_t)vm_stats.inactive_count + (int64_t)vm_stats.wire_count) *
             (int64_t)page_size;
  }
#endif
}

vector<int> greedyNearest(const vector<vector<float>> &dpts,
                          const vector<float> query, const int k_smallest) {
  std::priority_queue<std::pair<float, int>> top_candidates;
  float lower_bound = _INT_MAX;
  for (size_t i = 0; i < dpts.size(); i++) {
    float dist = EuclideanDistance(query, dpts[i]);
    if (top_candidates.size() < k_smallest || dist < lower_bound) {
      top_candidates.push(std::make_pair(dist, i));
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
  std::reverse(res.begin(), res.end());
  return res;
}

void ReadGroundtruthQuery(vector<vector<int>> &gt,
                          vector<std::pair<int, int>> &query_ranges,
                          vector<int> &query_ids, string gt_path) {
  ifstream infile;
  string bline;
  string delim = ",";
  string space_delim = " ";

  // int numCols = 0;
  infile.open(gt_path, ios::in);
  int counter = 0;
  while (getline(infile, bline, '\n')) {
    counter++;
    vector<int> one_gt;
    std::pair<int, int> one_range;
    int one_id;
    vector<string> ret;
    Split(bline, delim, &ret);
    one_id = std::stoi(ret[0]);
    one_range.first = std::stoi(ret[1]);
    one_range.second = std::stoi(ret[2]);
    vector<string> str_gt;
    Split(ret[7], space_delim, &str_gt);
    str_gt.pop_back();
    for (auto ele : str_gt) {
      one_gt.emplace_back(std::stoi(ele));
    }
    gt.emplace_back(one_gt);
    query_ranges.emplace_back(one_range);
    query_ids.emplace_back(one_id);
  }
}

void fvecs2csv(const string &output_path, const vector<vector<float>> &nodes) {
  std::ofstream file;
  file.open(output_path, std::ios_base::app);
  for (auto row : nodes) {
    if (file) {
      for (auto ele : row) {
        file << ele << " ";
      }
      file << "\n";
    }
  }
  file.close();
}

void evaluateKNNG(const vector<vector<int>> &gt,
                  const vector<vector<int>> &knng, const int K, double &recall,
                  double &precision) {
  assert(gt.size() == knng.size());

  double all_right = 0;
  int knng_amount = 0;

#pragma omp parallel for reduction(+ : all_right) reduction(+ : knng_amount)
  for (unsigned n = 0; n < gt.size(); n++) {
    double num_right = 0;
    // skip first, itself
    for (unsigned i = 1; i < K + 1; i++) {
      int one = gt[n][i];
      if (find(knng[n].begin(), knng[n].end(), one) != knng[n].end()) {
        num_right += 1;
      }
    }
    all_right += num_right;
    knng_amount += knng[n].size();
  }
  recall = (double)all_right / (K * gt.size());
  precision = (double)all_right / (float)knng_amount;
}

namespace pqdescent {
using std::vector;
FvecsItrReader::FvecsItrReader(std::string filename) {
  ifs.open(filename, std::ios::binary);
  assert(ifs.is_open());
  Next();
}

bool FvecsItrReader::IsEnd() { return eof_flag; }

std::vector<float> FvecsItrReader::Next() {
  std::vector<float> prev_vec = vec;  // return the currently stored vec
  int D;
  if (ifs.read((char *)&D, sizeof(int))) {  // read "D"
    // Then, read a D-dim vec
    vec.resize(D);                                            // allocate D-dim
    assert(ifs.read((char *)vec.data(), sizeof(float) * D));  // Read D * float.
    eof_flag = false;
  } else {
    vec.clear();
    eof_flag = true;
  }
  return prev_vec;
}

BvecsItrReader::BvecsItrReader(std::string filename) {
  ifs.open(filename, std::ios::binary);
  assert(ifs.is_open());
  Next();
}

bool BvecsItrReader::IsEnd() { return eof_flag; }

std::vector<float> BvecsItrReader::Next() {
  std::vector<float> prev_vec = vec;  // return the currently stored vec
  int D;
  if (ifs.read((char *)&D, sizeof(int))) {  // read "D"
    // Then, read a D-dim vec
    vec.resize(D);  // allocate D-dim
    std::vector<unsigned char> buff(D);

    assert(ifs.read((char *)buff.data(),
                    sizeof(unsigned char) * D));  // Read D * uchar.

    // Convert uchar to float
    for (int d = 0; d < D; ++d) {
      vec[d] = static_cast<float>(buff[d]);
    }

    eof_flag = false;
  } else {
    vec.clear();
    eof_flag = true;
  }
  return prev_vec;
}

ItrReader::ItrReader(std::string filename, std::string ext) {
  if (ext == "fvecs") {
    m_reader = (I_ItrReader *)new FvecsItrReader(filename);
  } else if (ext == "bvecs") {
    m_reader = (I_ItrReader *)new BvecsItrReader(filename);
  } else {
    std::cerr << "Error: strange ext type: " << ext << "in ItrReader"
              << std::endl;
    exit(1);
  }
}

ItrReader::~ItrReader() { delete m_reader; }

bool ItrReader::IsEnd() { return m_reader->IsEnd(); }

std::vector<float> ItrReader::Next() { return m_reader->Next(); }

std::vector<std::vector<float>> ReadTopN(std::string filename, std::string ext,
                                         int top_n) {
  std::vector<std::vector<float>> vecs;
  if (top_n != -1) {
    vecs.reserve(top_n);
  }
  ItrReader reader(filename, ext);
  while (!reader.IsEnd()) {
    if (top_n != -1 && top_n <= (int)vecs.size()) {
      break;
    }
    vecs.push_back(reader.Next());
  }
  return vecs;
}
}  // namespace pqdescent