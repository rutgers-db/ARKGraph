/**
 * Compact Window
 *
 * Author: Chaoji Zuo & ZhiZhi Wang
 * Email:  chaoji.zuo@rutgers.edu
 * Date:   Feb 13, 2020
 */
#ifndef COMPACTWINDOW_H
#define COMPACTWINDOW_H

#include <math.h>

#include <algorithm>
#include <iostream>
#include <tuple>
#include <unordered_map>
#include <vector>
using namespace std;

// class CompactWindow {
// public:
//   float bottomk[K];
//   int c1;
//   int ck;
//   int beg;
//   int end;
//   CompactWindow(int b, int e) : beg(b), end(e) {}
// };

class CompactWindow {
 public:
  int K_val;   // k smallest
  int lbound;  // left bound
  int rbound;  // right bound
  vector<int> pts;
  CompactWindow() {}
  CompactWindow(const int k) {
    K_val = k;
    pts.resize(K_val);
  }
  CompactWindow(const int k, const int l, const int r) {
    K_val = k;
    pts.resize(K_val);
    lbound = l;
    rbound = r;
  }
  CompactWindow(const int k, const vector<int> &vec) {
    K_val = k;
    for (int i = 1; i <= k; i++) {
      pts.emplace_back(vec[i]);
    }
    lbound = vec[0];
    rbound = vec[k + 1];
  }
  void print() {
    cout << lbound << " ";
    for (int i = 0; i < K_val; i++) {
      cout << pts[i] << " ";
    }
    cout << rbound << endl;
  }
};

// class TreeNode {
// public:
//   int next; // for leaf node, the are pointers
//   int prev; // for inner node, they are segments
//   int word_pos = -1;
//   bool has_visited = false;

//   TreeNode() {}
//   TreeNode(int l, int r) {
//     next = l;
//     prev = r;
//   }
// };

// // This is a segement tree, root is 1, left child is *2 and right *2 + 1
// // for each inner node, word_pos is -1 and [next, prev] is the segment
// // for each leaf node, word_pos is the word position in the original text,
// and [next, prev] are pointers
// // at begining the chain contains only header and tail.
// class Tree {
//     public:
//     int root = 1;
//     int leafNum; // the length of the original text
//     int HEAD;
//     int TAIL;

//     vector<TreeNode> nodes; //use one vector to represent a tree and only the
//     leaf nodes store the index value

//     Tree(int size_data) {

//         leafNum = size_data;
//         nodes.clear();
//         nodes.resize(4 * size_data);
//         buildTree(leafNum);

//         // initialize the chain
//         HEAD = 0;
//         TAIL = ((int)nodes.size()) - 1;
//         nodes[HEAD].next = TAIL;
//         nodes[HEAD].word_pos = -1;
//         nodes[TAIL].prev = HEAD;
//         nodes[TAIL].word_pos = leafNum;
//     }

//     void buildTree(int leafNum) {
//         buildTreeHelper(root, 0, leafNum - 1);
//     }

//     void buildTreeHelper(const int node, const int segSt, const int segEn) {
//         if (segSt == segEn) {
//             // This is the leaf node
//             nodes[node].prev = segSt;
//             nodes[node].next = segSt;
//             nodes[node].word_pos = segSt;
//         } else {
//             // Not a leaf node
//             int mid = (segSt + segEn) >> 1;
//             int leftChild = (node << 1);
//             int rightChild = (node << 1) + 1;

//             // First build left and then right
//             buildTreeHelper(leftChild, segSt, mid);
//             buildTreeHelper(rightChild, mid + 1, segEn);

//              // Assertion !
//             assert(nodes[leftChild].next + 1 == nodes[rightChild].prev);

//             // Set left and right value
//             nodes[node].prev = nodes[leftChild].prev;
//             nodes[node].next = nodes[rightChild].next;
//         }
//     }

//     inline int GetParent(int n) const {
//         return n >> 1;
//     }

//     inline int GetLftChild(int n) const {
//         return n << 1;
//     }

//     inline int GetRgtChild(int n) const {
//         return (n << 1) + 1;
//     }

//     /**
//     * given a word, update the chain using the tree by finding either its
//     next or prev visited nodes
//     * @param {Number} params
//     */

//     void printChain() {
//         int next = nodes[HEAD].next;
//         cout << "chain: ";
//         while (next != TAIL)
//         {
//             cout << nodes[next].word_pos << " ";
//             next = nodes[next].next;
//         }
//         cout << endl;
//     }

//     int Word2Node(int node, int word_pos, int &subtree) {
//         // lowest ancestor visited
//         if (nodes[node].has_visited)
//             subtree = node;

//         if (nodes[node].word_pos != word_pos)
//         {
//             // inner node
//             int leftChild = (node << 1);
//             int rightChild = (node << 1) + 1;
//             int mid = (nodes[node].prev + nodes[node].next) >> 1;
//             //if (nodes[leftChild].prev <= word_pos && word_pos <=
//             nodes[leftChild].next) if (word_pos <= mid)
//                 return Word2Node(leftChild, word_pos, subtree);
//             else
//                 return Word2Node(rightChild, word_pos, subtree);
//         } else {
//             // leaf node
//             return node;
//         }
//     }

//     int UpdateChain(int word_pos) {
//         int subtree = 0; // the root of the first visited subtree
//         int tree_node = Word2Node(root, word_pos, subtree);
//         // get previous and next visited node
//         int prev_node;
//         int next_node;
//         if (subtree == 0) // did not find any visited subtree
//         {
//             prev_node = HEAD;
//             next_node = TAIL;
//         } else if (nodes[GetLftChild(subtree)].has_visited) {
//             // find a subtree whose left child has visited nodes
//             prev_node = FindPrevVisited(GetLftChild(subtree));
//             next_node = nodes[prev_node].next;
//         } else if (nodes[GetRgtChild(subtree)].has_visited) {
//             // find a subtree whose right child has visited nodes
//             next_node = FindNextVisited(GetRgtChild(subtree));
//             prev_node = nodes[next_node].prev;
//         }

//         // update chain
//         nodes[prev_node].next = tree_node;
//         nodes[next_node].prev = tree_node;
//         nodes[tree_node].next = next_node;
//         nodes[tree_node].prev = prev_node;

//         // update the tree
//         nodes[tree_node].has_visited = true;
//         int parent = GetParent(tree_node);
//         while (parent > 0 && !nodes[parent].has_visited) {
//             nodes[parent].has_visited = true;
//             parent = GetParent(parent);
//         }
//         return tree_node;
//     }

//     int FindNextVisited(int tree_node) {
//         while (nodes[tree_node].word_pos == -1) {
//             // inner node
//             if (nodes[GetLftChild(tree_node)].has_visited)
//                 tree_node = GetLftChild(tree_node);
//             else // use right child instead
//                 tree_node = GetRgtChild(tree_node);
//         }
//         return tree_node;
//     }

//     int FindPrevVisited(int tree_node) {
//         while (nodes[tree_node].word_pos == -1) {
//             // inner node
//             if (nodes[GetRgtChild(tree_node)].has_visited)
//                 tree_node = GetRgtChild(tree_node);
//             else // use left child instead
//                 tree_node = GetLftChild(tree_node);
//         }
//         return tree_node;
//     }

//     bool FindCompactWindows(vector<float> &data, int word_pos,
//     vector<CompactWindow> &cws) {
//         cout << "inserting " << word_pos << endl;
//         int tree_node = UpdateChain(word_pos);
//         printChain();

//         // sliding windows, add words in left and right, at most K of them
//         vector<int> prev_words;
//         vector<int> next_words;

//         int prev_node = nodes[tree_node].prev;
//         int next_node = nodes[tree_node].next;

//         for (int i = 0; i < K; i++) {
//             if (prev_node != HEAD) {
//                 prev_words.push_back(nodes[prev_node].word_pos);
//                 prev_node = nodes[prev_node].prev;
//             } else {
//                 prev_words.push_back(nodes[prev_node].word_pos);
//                 break;
//             }
//         }

//         for (int i = 0; i < K; i++) {
//             if (next_node != TAIL) {
//                 next_words.push_back(nodes[next_node].word_pos);
//                 next_node = nodes[next_node].next;
//             } else {
//                 next_words.push_back(nodes[next_node].word_pos);
//                 break;
//             }
//         }

//         // not enough
//         if (next_words.size() + prev_words.size() + 1 < K + 2)
//             return false;

//         // enought, slide the window of size K + 2 (including the starting
//         and ending positions) int beg_idx = ((int)prev_words.size()) - 1;  //
//         left bound of the winodw int end_idx = K - beg_idx - 1; // right
//         bound of the window

//         int num = cws.size(); // for printing

//         // slide until left or right bound invalid
//         while (beg_idx >= 0 && end_idx < next_words.size()) {

//             // add a compact window
//             cws.emplace_back(prev_words[beg_idx] + 1, next_words[end_idx] -
//             1);

//             // add prev nodes to the bottomk
//             int cnt = 0;
//             if (beg_idx == 0) {
//                 cws.back().c1 = word_pos;
//             } else {
//                 cws.back().c1 = prev_words[beg_idx - 1];
//                 for (int i = beg_idx - 1; i >= 0; i--)
//                     cws.back().bottomk[cnt++] = data[prev_words[i]];
//             }

//             // add current node to the bottoml
//             cws.back().bottomk[cnt++] = data[word_pos];

//             // add next nodes to the bottomk
//             if (end_idx == 0) {
//                 cws.back().ck = word_pos;
//             } else {
//                 cws.back().ck = next_words[end_idx - 1];
//                 for (int i = 0; i < end_idx; i++)
//                     cws.back().bottomk[cnt++] = data[next_words[i]];
//             }

//             assert(cnt == K);

//             beg_idx--;
//             end_idx++;
//         }

//         for (int i = num; i < cws.size(); i++)
//         {
//             cout << "l, c1, ck, r: " << cws[i].beg << " " << cws[i].c1 << " "
//             << cws[i].ck << " " << cws[i].end << endl; cout << "hash values
//             are: "; for (int j = 0; j < K; j++)
//                 cout << cws[i].bottomk[j] << " ";
//             cout << endl;
//         }

//         return true;
//     }
// };

// void optScanMethod(vector<float> &data, vector<CompactWindow> &cws){
//     vector<pair<float, int>> hash_pos; // hash value and word position pairs;
//     for (auto pos = 0; pos < data.size(); pos++)
//         hash_pos.emplace_back(data[pos], pos);

//     // first sort by hash value and second by position
//     sort(hash_pos.begin(), hash_pos.end(), [](const pair<float, int> &p1,
//     const pair<float, int> &p2){
//         if (p1.first < p2.first)
//             return true;
//         else if (p1.first > p2.first)
//             return false;
//         else
//             return p1.second < p2.second;
//     });

//     for (auto &entry : hash_pos)
//         cout << entry.first << " " << entry.second << endl;
//     cout << "sorted by hash" << endl;

//     Tree tree(data.size());
//     for (auto &entry : hash_pos)
//     {
//         int word_pos = entry.second;
//         tree.FindCompactWindows(data, word_pos, cws);
//     }
// }

#endif