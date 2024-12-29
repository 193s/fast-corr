#pragma once
#include <vector>
#include <deque>
#include <stack>
#include <map>

#include "lazy-reversible-rbst.hpp"

// G++ extensions
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
#include <ext/pb_ds/tag_and_trait.hpp>

// assuming N<INT_MAX
#define d1_type long long // sum(d_i) can be as large as O(N^2)
#define d2_type long long // sum(d_i^2) can be as large as (N^3-N)/6 <=> rho=0
// Sx = sum(t_i^3 - t_i) can be as large as N^3-N so d2_type can be used as well

const int TWO = 2; // multiply all ranks by 2 so that they can be treated as integer

template< class T >
class OnlineSpearmanBase {
  public:
  virtual void push_back(T x_val);
  virtual void push_front(T x_val);
  virtual void pop_front();
  virtual void pop_back();
  virtual d2_type spearman_d();
  virtual int size();
  double spearman_r() {
    d2_type d = spearman_d();
    // d = sum((2d_i)^2) = 4*actual_D
    int n = size();
    if (n == 1) return NAN;
    d2_type n3 = (d2_type)n*((d2_type)n*n-1);
    //cout<<"n="<<n<<", d="<<d<<"/4, Sx="<<Sx<<"\n";
    if (Sx == 0) {
      return 1.0 - (6.0/(TWO*TWO))*d / n3;
    }
    else {
      // general formula:
      // (1-(6.0/n3)*(D + Sx/12 + Sy/12)) / (sqrt(1 - Sx/n3) * sqrt(1-Sy/n3))
      // when Sy = 0,
      // = (1-(6.0/n3)*(D + Sx/12)) / (sqrt(1 - Sx/n3))
      // = (n3-6.0*(D + Sx/12)) / sqrt(n3*n3 - Sx*n3)
      // = (n3 - Sx/2.0 - 6.0*D) / sqrt(n3*n3 - Sx*n3)
      // = (2*n3 - Sx - 3*(4D)) / (2*sqrt(n3*n3 - Sx*n3))
      return (2.0*n3 - Sx - (12.0/(TWO*TWO))*d) / (2.0*n3*sqrt(1.0 - (double)Sx/(double)n3));
    }
  }
  protected:
  d2_type Sx = 0; // sum(t_i^3 - t_i)
  // Sy = 0 under the sliding constraint
};

// multiply all ranks by 2 to treat them as integer
// e.g. [3, 12123, 0] -> [2, 3, 1]*2
// e.g. [1, 2, 2, 2, 5, 5, 7] -> [1, 3, 3, 3, 5.5, 5.5, 7]*2
template< class T >
std::vector<int> convert_array_to_rank(std::vector<T> X) {
  int n = X.size();
  if (n == 0) return {};
  if (n == 1) return {1*TWO};
  // n>=2
  std::vector<std::pair<T, int> > X2(n);
  for (int i=0; i<n; i++) X2[i] = std::pair<T, int>(X[i], i);
  std::sort(X2.begin(), X2.end());
  std::vector<int> ret(n);
  //for (int i=0; i<n; i++) ret[X2[i].second] = TWO*(i+1); // works only on unique arrays
  int z = 0;
  T last_seen;
  std::stack<int> st;
  for (int i=0; i<n; i++) {
    int pos = X2[i].second;
    if (i > 0 && last_seen != X2[i].first) {
      // finalize rank
      int rank = z*TWO + 1 + st.size();// 1 + z + (st.size()-1)/2
      z += st.size();
      while (!st.empty()) {
        ret[st.top()] = rank;
        st.pop();
      }
    }
    st.push(pos);
    last_seen = X2[i].first;
  }
  if (!st.empty()) {
    // finalize rank
    int rank = z*TWO + 1 + st.size();// 1 + z + (st.size()-1)/2
    z += st.size();
    while (!st.empty()) {
      ret[st.top()] = rank;
      st.pop();
    }
  }
  //cout<<"{";for (T x:X)cout<<x<<",";cout<<"} -> ";
  //cout<<"{";for (int x:ret)cout<<x/2.0<<",";cout<<"}\n";
  return ret;
}

struct SNode {
  d1_type d1;
  d2_type d2;
  int cnt;
  SNode() {
    d1 = d2 = 0;
    cnt = 1;
  }
  SNode(d1_type x1, d2_type x2, int c) {
    d1 = x1;
    d2 = x2;
    cnt = c;
  }
};
inline SNode F(SNode x, SNode y) {
  return SNode(x.d1 + y.d1, x.d2 + y.d2, x.cnt + y.cnt);
}
inline SNode addall(SNode x, int a) {
  if (a == 0) return x;
  int size = x.cnt;
  x.d2 += (d2_type)a*(a*size + 2*x.d1); // d2 += a*a*size + 2*a*d1
  x.d1 += (d1_type)a*size;
  return x;
}
inline int none_h(int x, int y) { return x+y; }
inline SNode ts(SNode a) { return a; }

//< class D, class L, D (*f)(D, D), D (*g)(D, L), L (*h)(L, L), L (*p)(L, int) >

// Efficient Online Implementation of Spearman's rank correlation with Binary Search Tree
// adding/removing a value requires O(logN) time complexity and O(N) space complexity
template< class T >
class OnlineSpearman : public OnlineSpearmanBase<T> {
  using RBTree = LazyReversibleRBST< SNode, int, F, addall, none_h, ts >;
  RBTree tree;
  using CountingTree = __gnu_pbds::tree<
    std::pair<T, int>, __gnu_pbds::null_type, std::less<std::pair<T, int> >,
    __gnu_pbds::rb_tree_tag, __gnu_pbds::tree_order_statistics_node_update>;
  CountingTree ctr_tree;
  int N = 0;
  int id_for_tree = 0; // add unique ids to allow for duplicate values in CountingTree
  inline std::pair<int, int> _add_value(T x_val) {
    int z = ctr_tree.order_of_key(std::make_pair(x_val, -1)); // # of < x_val
    int dup = ctr_tree.order_of_key(std::make_pair(x_val, id_for_tree+111)) - z;
    OnlineSpearmanBase<T>::Sx += (d2_type)3*dup*(dup+1);
    ctr_tree.insert(std::make_pair(x_val, id_for_tree++));
    return std::make_pair(z, dup);
  }
  inline std::pair<int, int> _remove_value(T x_val) {
    int z = ctr_tree.order_of_key(std::make_pair(x_val, -1));
    int dup = ctr_tree.order_of_key(std::make_pair(x_val, id_for_tree+111)) - z - 1;
    ctr_tree.erase(ctr_tree.find_by_order(z)); // erase @z (@z+dup) will also work
    // ctr_tree.erase(std::make_pair(x_val, <id to erase>));
    assert(dup >= 0);
    OnlineSpearmanBase<T>::Sx -= (d2_type)3*dup*(dup+1);
    return std::make_pair(z, dup);
  }

  public:
    std::deque<T> real_vals;
    RBTree::Node *root = tree.make_tree();
    OnlineSpearman() {}
    OnlineSpearman(std::vector<T> x_vals) {
      for (auto &x : x_vals) push_back(x);
    }

    // add a new element
    void push_back(T x_val) {
      real_vals.push_back(x_val);
      auto p = _add_value(x_val);
      int z = p.first, dup = p.second;

      // new element with (x=z+dup/2, y=N) -> d1 = z-N+(dup/2)
      d1_type new_d1 = (z-N)*TWO + dup; // *= 2
      SNode newelem = SNode(new_d1, (d2_type)new_d1*new_d1, 1);

      if (root == NULL) root = tree.build({newelem});
      else {
        assert(tree.size(root) == N);
        // add +1(*2) to [z+dup, N) & add  +0.5(*2) to [z, z+dup)
        if (z+dup < N) tree.apply(root, z+dup, N, +1*TWO); // +1(*2) to [z+dup, N)
        if (dup > 0)   tree.apply(root, z, z+dup, +1); // +0.5(*2) to [z, z+dup)
        tree.insert(root, z+dup, newelem); // insert new element @ z+dup
      }
      N += 1;
    }
    // add a new element
    void push_front(T x_val) {
      real_vals.push_front(x_val);
      auto p = _add_value(x_val);
      int z = p.first, dup = p.second;

      // new element with (x=z+dup/2, y=0) -> d1 = z-0+(dup/2)
      d1_type new_d1 = (z-0)*TWO + dup; // *= 2
      SNode newelem = SNode(new_d1, (d2_type)new_d1*new_d1, 1);

      if (root == NULL) root = tree.build({newelem});
      else {
        assert(tree.size(root) == N);
        if (dup > 0) tree.apply(root, z, z+dup, -1); // [z, z+dup) -= 0.5(*2)
        if (z > 0)   tree.apply(root, 0, z, -1*TWO); // [0, z) -= 1(*2) (y += 1)
        tree.insert(root, z, newelem); // insert new element @ z
      }
      N += 1;
    }
    // remove an oldest element
    void pop_front() {
      T x_val = real_vals.front();
      real_vals.pop_front();
      auto p = _remove_value(x_val);
      int z = p.first, dup = p.second;

      tree.erase(root, z); // erase @z
      if (z>0)   tree.apply(root, 0, z, +1*TWO); // [0, z) += 1 (*2)
      if (dup>0) tree.apply(root, z, z+dup, +1); // [z, z+dup) += 0.5(*2) (y-=1, x-=0.5)
      N -= 1;
    }
    // remove an oldest element
    void pop_back() {
      T x_val = real_vals.back();
      real_vals.pop_back();
      auto p = _remove_value(x_val);
      int z = p.first, dup = p.second;

      tree.erase(root, z+dup); // erase @z+dup
      if (dup>0)       tree.apply(root, z, z+dup, -1);       // [z, z+dup) -= 0.5(*2) (y-=0, x-=0.5)
      if (z+dup < N-1) tree.apply(root, z+dup, N-1, -1*TWO); // [z+dup, )  -= 1(*2)
      N -= 1;
    }
    d2_type spearman_d() {
      if (root == NULL) return 0;
      return root->sum.d2;
    }
    int size() { return N; }
};
// Online Implementation of Spearman's rank correlation without Binary Search Tree
// adding/removing a value requires O(N) time complexity and O(N) space complexity
template< class T >
class OnlineSpearmanLinear : public OnlineSpearmanBase<T> {
  int N = 0;
  std::vector<int> D;
  std::deque<T> X_val;
  std::map<T, int> duplicate_counter;
  // internal function: returns the number of existing pairs with the same x-values (>=0)
  inline int _add_value(T x_val) {
    int dup = duplicate_counter[x_val]++;
    OnlineSpearmanBase<T>::Sx += (d2_type)3*dup*(dup+1);
    return dup;
  }
  // internal function: returns the number of remaining pairs with the same x-values (excluding the pair being erased, >=0)
  inline int _remove_value(T x_val) {
    int dup = --duplicate_counter[x_val]; //assert dup>=0
    if (dup == 0) duplicate_counter.erase(x_val);
    OnlineSpearmanBase<T>::Sx -= (d2_type)3*dup*(dup+1);
    return dup;
  }
  public:
    OnlineSpearmanLinear() {}
    OnlineSpearmanLinear(std::vector<T> x_vals) {
      for (auto &x : x_vals) push_back(x);
    }
    // add a new element
    void push_back(T x_val) {
      int dup = _add_value(x_val);
      int z = 0;
      for (T x : X_val) if (x < x_val) z++;
      X_val.push_back(x_val);
      D.push_back(0);
      for (int i=N-1; i>=z+dup; i--) D[i+1] = D[i] + 1*TWO; // [z+dup, N) += 1(*2)
      for (int i=z+dup-1; i>=z; i--) D[i] = D[i] + 1;     // [z, z+dup) += 0.5(*2) (& insert@z+dup)
      D[z+dup] = (z - N)*TWO + dup; // D_i := X_i - Y_i (*2)
      N += 1;
    }
    void push_front(T x_val) {
      int dup = _add_value(x_val);
      int z = 0;
      for (T x : X_val) if (x < x_val) z++;
      X_val.push_front(x_val);
      D.push_back(0);
      for (int i=N-1; i>=z+dup; i--) D[i+1] = D[i];     // [z+dup, N) += 0 (*2) (x += 1, y += 1)
      for (int i=z+dup-1; i>=z; i--) D[i+1] = D[i] - 1; // [z, z+dup) += 0.5-1(*2) (& insert@z)
      for (int i=z-1; i>=0; i--)     D[i] -= 1*TWO;     // [0, z) -= 1(*2) (y += 1)
      D[z] = (z - 0)*TWO + dup; // D_i := X_i - Y_i (*2)
      N += 1;
    }
    // remove an oldest element
    void pop_front() {
      T x_val = X_val.front();
      int dup = _remove_value(x_val);
      int z = 0;
      for (T x : X_val) if (x < x_val) z++;
      for (int i=0; i<z; i++)       D[i] += 1*TWO;     // [0, z) += 1 (*2)
      for (int i=z; i<z+dup; i++)   D[i] = D[i+1] + 1; // erase @z & [z, z+dup) += 0.5(*2)
      for (int i=z+dup; i<N-1; i++) D[i] = D[i+1];
      N -= 1;
      X_val.pop_front();
      D.pop_back();
    }
    void pop_back() {
      T x_val = X_val.back();
      int dup = _remove_value(x_val);
      int z = 0;
      for (T x : X_val) if (x < x_val) z++;
      //for (int i=0; i<z; i++) D[i] += 1*TWO; // [0, z) += 0-0 (*2)
      for (int i=z; i<z+dup; i++)   D[i] -= 1;             // [z, z+dup) -= 0.5(*2) (y-=0, x-=0.5)
      for (int i=z+dup; i<N-1; i++) D[i] = D[i+1] - 1*TWO; // erase@z+dup & [z+dup, ) -= 1(*2)
      N -= 1;
      X_val.pop_back();
      D.pop_back();
    }
    d2_type spearman_d() {
      d2_type d = 0;
      for (int i=0; i<N; i++) d += (d2_type)D[i]*D[i];
      return d;
    }
    int size() { return N; }
};


// Traditional Implementation of Spearman's rank correlation:
// adding/removing a value requires O(NlogN) time complexity and O(M) space complexity
template< class T >
class OfflineSpearman : public OnlineSpearmanBase<T> {
  int N = 0;
  std::deque<T> X_val;
  std::map<T, int> duplicate_counter;
  // internal function: returns the number of existing pairs with the same x-values (>=0)
  inline int _add_value(T x_val) {
    int dup = duplicate_counter[x_val]++;
    OnlineSpearmanBase<T>::Sx += (d2_type)3*dup*(dup+1);
    return dup;
  }
  // internal function: returns the number of remaining pairs with the same x-values (excluding the pair being erased, >=0)
  inline int _remove_value(T x_val) {
    int dup = --duplicate_counter[x_val]; //assert dup>=0
    if (dup == 0) duplicate_counter.erase(x_val);
    OnlineSpearmanBase<T>::Sx -= (d2_type)3*dup*(dup+1);
    return dup;
  }
  public:
    OfflineSpearman() {}
    OfflineSpearman(std::vector<T> x_vals) {
      for (auto &x : x_vals) push_back(x);
    }
    // add a new element
    void push_back(T x_val) {
      _add_value(x_val);
      N += 1;
      X_val.push_back(x_val);
    }
    // add a new element
    void push_front(T x_val) {
      _add_value(x_val);
      N += 1;
      X_val.push_front(x_val);
    }
    // remove an oldest element
    void pop_front() {
      T x_val = X_val.front();
      _remove_value(x_val);
      N -= 1;
      X_val.pop_front();
    }
    // remove an newest element
    void pop_back() {
      T x_val = X_val.back();
      _remove_value(x_val);
      N -= 1;
      X_val.pop_back();
    }
    d2_type spearman_d() {
      int n = X_val.size();
      // convert deque to vector
      std::vector<T> X_val2;
      for (T &x : X_val) X_val2.push_back(x);
      std::vector<int> X = convert_array_to_rank(X_val2);
      std::vector<int> Y(n);
      for (int i=1; i<=n; i++) Y[i-1] = i*TWO; // *2
      d2_type d = 0;
      for (int i=0; i<n; i++) d += (d2_type)(X[i]-Y[i])*(X[i]-Y[i]);
      return d;
    }
    int size() { return N; }
};
