/**
 * 2D segment tree to store compact windows
 *
 * Author: Chaoji Zuo
 * Date:   Nov 19 2021
 * Email:  chaoji.zuo@rutgers.edu
 */
#ifndef SEGMENTTREE_H
#define SEGMENTTREE_H
#include <iostream>
#include <unordered_map>

#include "cw.h"
using std::cout;
using std::endl;

#define MAX_IDX 10000000

namespace segTree {

class TreeNode {
 public:
  int nid;
  int lid;
  int rid;
  bool is_left;
  TreeNode *leftChild;
  TreeNode *rightChild;
  TreeNode *secondRoot;
  CompactWindow *cw;
  TreeNode()
      : nid(0),
        lid(0),
        rid(0),
        leftChild(nullptr),
        rightChild(nullptr),
        secondRoot(nullptr),
        cw(nullptr) {}
};

void print(const TreeNode *tree, bool print_cw = false) {
  queue<TreeNode *> q;
  q.push(new TreeNode(*tree));
  while (!q.empty()) {
    int count = q.size();
    while (count > 0) {
      auto tmp = q.front();
      q.pop();
      if (tmp->leftChild != nullptr) {
        q.push(new TreeNode(*tmp->leftChild));
        cout << "-";
      }
      if (tmp->rightChild != nullptr) {
        q.push(new TreeNode(*tmp->rightChild));
        cout << "-";
      }
      cout << "(" << tmp->lid << "->" << tmp->rid << ")\t";

      if (tmp->lid == 498 && tmp->rid == 2000) {
        print(tmp->secondRoot, true);
      }
      if (print_cw) {
        cout << "   ";
        if (tmp->cw != nullptr)
          tmp->cw->print();
        else {
          cout << endl;
        }
      }

      count--;
    }
    cout << endl;
  }
}

TreeNode *createTreeTopDown(int min, int max) {
  TreeNode *rootNode = new TreeNode();
  rootNode->lid = min;
  rootNode->rid = max;
  rootNode->leftChild = nullptr;
  rootNode->rightChild = nullptr;
  if (max - min == 1) {
    return rootNode;
  } else if (max - min > 1) {
    int mid = (max + min) >> 1;
    rootNode->leftChild = createTreeTopDown(min, mid);
    rootNode->rightChild = createTreeTopDown(mid, max);
  }
  return rootNode;
}

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

void createSegmentTreeSecondDim(TreeNode *tree, const vector<int> &pts) {
  if (tree->leftChild != nullptr) {
    createSegmentTreeSecondDim(tree->leftChild, pts);
  }
  if (tree->rightChild != nullptr) {
    createSegmentTreeSecondDim(tree->rightChild, pts);
  }
  // if (tree->secondRoot == nullptr)
  // {
  // cout << tree->lid << " , " << tree->rid << endl;
  tree->secondRoot = createSegmentTree(pts);
  // cout << tree->secondRoot->rightChild->lid << " , " <<
  // tree->secondRoot->rightChild->rid << endl;

  // }
}

void createTreeSecondDim(const int min, const int max, TreeNode *tree) {
  if (tree->leftChild != nullptr) {
    createTreeSecondDim(min, max, tree->leftChild);
  }
  if (tree->rightChild != nullptr) {
    createTreeSecondDim(min, max, tree->rightChild);
  }
  if (tree->leftChild == nullptr && tree->rightChild == nullptr) {
    tree->secondRoot = createTreeTopDown(min, max);
  }
}

void insertSegment(const int rmin, const int rmax, const int la, const int lb,
                   TreeNode *tree, CompactWindow *cw, bool is_left = true,
                   int ra = -1, int rb = -1) {
  // if (tree != nullptr)
  // {
  //   cout << "inserting: " << (is_left ? "left tree " : "right tree  ") << la
  //   << "  " << lb << endl; cout << "tree: " << tree->lid << "  " << tree->rid
  //   << endl;
  // }
  if (tree == nullptr || la < tree->lid || lb > tree->rid) {
    // cout << "error! wrong position!" << endl;
    return;
  }
  if ((tree->lid == la && tree->rid == lb) || (tree->rid - tree->lid) == 1) {
    // if it is ltree, start to find rtree;
    //  get the node;
    if (is_left) {
      if (tree->secondRoot == nullptr) {
        cout << "error! should not get here!" << endl;
        tree->secondRoot = createTreeTopDown(rmin, rmax);
      }
      assert(tree->secondRoot != nullptr);
      // cout << "trying" << endl;

      insertSegment(rmin, rmax, ra, rb, tree->secondRoot, cw, false);
      // cout << "after insert right " << endl;
    } else {
      tree->cw = cw;
    }

    return;
  }
  // int mid = (tree->lid + tree->rid) >> 1;
  int mid = tree->rightChild->lid;
  // cout << "mid: " << mid << endl;
  if (lb <= mid) {
    insertSegment(rmin, rmax, la, lb, tree->leftChild, cw, is_left, ra, rb);
  } else if (la >= mid) {
    insertSegment(rmin, rmax, la, lb, tree->rightChild, cw, is_left, ra, rb);
  } else if (la != mid) {
    insertSegment(rmin, rmax, la, mid, tree->leftChild, cw, is_left, ra, rb);
    insertSegment(rmin, rmax, mid, lb, tree->rightChild, cw, is_left, ra, rb);
  }
}

CompactWindow *searchCW(const int pos, TreeNode *tree, const int rpos = -2) {
  // cout << "search in : " << (rpos != -2 ? "left tree " : "right tree  ") <<
  // "for "
  //      << "  " << pos << endl;
  // cout << "tree: " << tree->lid << "  " << tree->rid << endl;
  // if (rpos == -2)
  // {
  //   print(tree);
  // }
  if (rpos == -2 && tree->cw != nullptr) {
    // cout << "find!" << endl;
    return tree->cw;
  }
  if (rpos != -2 && tree->secondRoot != nullptr) {
    // cout << "current position: " << pos << endl;
    auto searchRes = searchCW(rpos, tree->secondRoot);
    if (searchRes != nullptr) {
      return searchRes;
    } else {
      // cout << "fail to search in this position" << endl;
    }
  }

  if (tree->leftChild == nullptr || tree->rightChild == nullptr) {
    return nullptr;
  }

  // int mid = (tree->lid + tree->rid) >> 1;
  int mid = tree->rightChild->lid;

  if (pos <= mid) {
    return searchCW(pos, tree->leftChild, rpos);
  } else {
    return searchCW(pos, tree->rightChild, rpos);
  }
  return nullptr;
}

}  // namespace segTree

#endif