#pragma once
#include <vector>
#include <deque>
#include <set>
using namespace std;

#include "lazy-reversible-rbst.hpp"
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
#include <ext/pb_ds/tag_and_trait.hpp>
using namespace __gnu_pbds;

template< class T >
class OnlineSpearmanBase {
  public:
  virtual void push_back(T x_val);
  virtual void pop_front();
  virtual int spearman_d();
  virtual int size();
  double spearman_r() {
    int d = spearman_d(), n = size();
    return 100*(1.0 - 6.0*d / (n*(n*n-1)));
  }
};

struct SNode {
  int d1, d2;
  int cnt;
  SNode() {
    d1 = d2 = 0;
    cnt = 1;
  }
  SNode(int x1, int x2, int c) {
    d1 = x1;
    d2 = x2;
    cnt = c;
  }
};
inline SNode F(SNode x, SNode y) {
  SNode ret(x.d1 + y.d1, x.d2 + y.d2, x.cnt + y.cnt);
  return ret;
}
inline SNode addall(SNode x, int a) {
  if (a == 0) return x;
  int size = x.cnt;
  x.d2 += a*a*size + 2*a*x.d1;
  x.d1 += a*size;
  return x;
}
inline int none_h(int x, int y) { return x+y; }
SNode ts(SNode a) { return a; }

//< class D, class L, D (*f)(D, D), D (*g)(D, L), L (*h)(L, L), L (*p)(L, int) >

// Efficient Online Implementation of Spearman's rank correlation with Binary Search Tree
// adding/removing a value requires O(logN) timespace
template< class T >
class OnlineSpearman : public OnlineSpearmanBase<T> {
  using RBTree = LazyReversibleRBST< SNode, int, F, addall, none_h, ts >;
  RBTree tree;
  //using RBCountingTree = RBSTBase<T>;
  using CountingTree = __gnu_pbds::tree<T,null_type,less<T>,rb_tree_tag,tree_order_statistics_node_update>;
  CountingTree ctr_tree;
  int N = 0;

  public:
    deque<T> real_vals; // temporary
    RBTree::Node *root = tree.make_tree();

    // add a new element
    void push_back(T x_val) {
      real_vals.push_back(x_val);
      ctr_tree.insert(x_val);
      int z = ctr_tree.order_of_key(x_val)+1;

      // new element with (x=N+1, y=z) -> d1 = z-(N+1)
      int new_d1 = z-(N+1);
      SNode newelem = SNode(new_d1, new_d1*new_d1, 1);

      auto ptree = tree.build({newelem});
      if (root == NULL) {
        root = ptree;
      }
      else {
        if (z-1 < tree.size(root)) {
          tree.apply(root, z-1, tree.size(root), +1); // add +1 to [z-1, )
        }
        tree.insert(root, z-1, newelem); // insert new element
      }
      N += 1;
    }
    // remove an oldest element
    void pop_front() {
      T x_val = real_vals.front();
      real_vals.pop_front();

      int z = ctr_tree.order_of_key(x_val)+1;
      ctr_tree.erase(x_val);

      if (z>1) tree.apply(root, 0, z-1, +1);
      tree.erase(root, z-1);
      N -= 1;
    }
    int spearman_d() {
      if (root == NULL) return 0;
      return root->sum.d2;
    }
    int size() { return N; }
};
// Online Implementation of Spearman's rank correlation without Binary Search Tree
// adding/removing a value requires O(N) timespace
template< class T >
class OnlineSpearmanLinear : public OnlineSpearmanBase<T> {
  int N = 0;
  vector<int> D;
  deque<T> X_val;
  public:
    /*
    OnlineSpearmanLinear() {}
    OnlineSpearmanLinear(vector< pair<T, T> > arr) {
      N = arr.size();
      X_val = deque<T>(arr);
      vector<T> first, second;
      for (auto p : arr) first.push_back(arr.first);
      for (auto p : arr) second.push_back(arr.second);
      vector<int> xRank = convert_to_rank(first);
      vector<int> yRank = convert_to_rank(second);
      for (int i=0; i<N; i++) D[i] = xRank[i] - yRank[i];
    }
    */

    // add a new element
    void push_back(T x_val) {
      int z = 0;
      for (T x : X_val) if (x < x_val) z++;
      X_val.push_back(x_val);

      D.push_back(0);
      for (int i=N-1; i>=z; i--) D[i+1] = D[i] + 1; // [z, N) += 1 & insert z
      D[z] = z - N; // D_i := X_i - Y_i
      N += 1;
    }
    // remove an oldest element
    void pop_front() {
      T x_val = X_val.front();
      int z = 0;
      for (T x : X_val) if (x < x_val) z++;
      for (int i=z; i<N-1; i++) D[i] = D[i+1];
      for (int i=0; i<z; i++) D[i] += 1; // [0, k) += 1
      N -= 1;
      X_val.pop_front();
      D.pop_back();
    }
    int spearman_d() {
      int d = 0;
      for (int i=0; i<N; i++) d += D[i]*D[i];
      //cout<<"d="<<d<<"\n";
      return d;
    }
    int size() { return N; }
};


// Traditional Implementation of Spearman's rank correlation:
// adding/removing a value requires O(NlogN) timespace
// (this can be done in O(N) if implemented properly)
template< class T >
class OfflineSpearman : public OnlineSpearmanBase<T> {
  int N = 0;
  deque<T> X_val;
  public:
    // add a new element
    void push_back(T x_val) {
      N += 1;
      X_val.push_back(x_val);
    }
    // remove an oldest element
    void pop_front() {
      N -= 1;
      X_val.pop_front();
    }
    int spearman_d() {
      int n = X_val.size();
      // convert deque to vector
      vector<T> X_val2;
      for (T x : X_val) X_val2.push_back(x);
      vector<int> X = convert_to_rank(X_val2);
      vector<int> Y(n);
      for (int i=1; i<=n; i++) Y[i-1] = i;
      int d = 0;
      for (int i=0; i<n; i++) d += (X[i]-Y[i])*(X[i]-Y[i]);
      //cout<<"d="<<d<<"\n";
      return d;
    }
    int size() { return N; }

  private:
  // e.g. [3, 12123, 0] -> [2, 3, 1]
  vector<int> convert_to_rank(vector<T> X) {
    vector<pair<T, int> > X2(X.size());
    for (int i=0; i<X.size(); i++) X2[i] = pair<T, int>(X[i], i);
    sort(X2.begin(), X2.end());
    vector<int> ret(X.size());
    for (int i=0; i<X.size(); i++) ret[X2[i].second] = i+1;
    return ret;
  }
};

