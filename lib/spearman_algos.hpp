#pragma once
#include <vector>
#include <deque>
#include <map>
#include "fast_corr_base.hpp"
#include "lazy_reversible_rbst.hpp"

namespace FastCorr::MonotonicOnlineCorr {
  /**
   * Internal node class defined for Spearman algorithm
   */
  struct SNode {
    sp_d1_type d1;
    sp_d2_type d2;
    int cnt;
    SNode() { d1 = d2 = 0, cnt = 1; }
    SNode(sp_d1_type x1, sp_d2_type x2, int c) { d1 = x1, d2 = x2, cnt = c; }
  };
  inline SNode F(SNode x, SNode y) {
    return SNode(x.d1 + y.d1, x.d2 + y.d2, x.cnt + y.cnt);
  }
  inline SNode addall(SNode x, int a) {
    if (a == 0) return x;
    int size = x.cnt;
    x.d2 += (sp_d2_type)a*(a*size + 2*x.d1); // d2 += a*a*size + 2*a*d1
    x.d1 += (sp_d1_type)a*size;
    return x;
  }
  inline int none_h(int x, int y) { return x+y; }
  inline SNode ts(SNode a) { return a; }
  //  < class D, class L, D (*f)(D, D), D (*g)(D, L), L (*h)(L, L), L (*p)(L, int) >

  /**
   * Efficient Online Implementation of Spearman's rank correlation with Binary Search Tree
   * adding/removing a value requires O(logN) time complexity and O(N) space complexity
   */
  template< class T >
  class Spearman : public SpearmanBase<T> {
    using RBTree = LazyReversibleRBST< SNode, int, F, addall, none_h, ts >;
    RBTree tree;
    CountingTree< std::pair<T, int> > ctr_tree;
    int N = 0;
    int id_for_tree = 0; // add unique ids to allow for duplicate values in CountingTree
    inline std::pair<int, int> _add_value(const T &x_val) {
      int z = ctr_tree.order_of_key(std::make_pair(x_val, -1)); // # of < x_val
      int dup = ctr_tree.order_of_key(std::make_pair(x_val, id_for_tree)) - z;
      SpearmanBase<T>::Sx += (sp_d2_type)3*dup*(dup+1);
      ctr_tree.insert(std::make_pair(x_val, id_for_tree++));
      return std::make_pair(z, dup);
    }
    inline std::pair<int, int> _remove_value(const T &x_val) {
      int z = ctr_tree.order_of_key(std::make_pair(x_val, -1));
      int dup = ctr_tree.order_of_key(std::make_pair(x_val, id_for_tree)) - z - 1;
      ctr_tree.erase_kth(z); // erase @z (@z+dup) will also work
      // ctr_tree.erase(std::make_pair(x_val, <id to erase>));
      assert(dup >= 0);
      SpearmanBase<T>::Sx -= (sp_d2_type)3*dup*(dup+1);
      return std::make_pair(z, dup);
    }

    public:
      std::deque<T> real_vals;
      RBTree::Node *root = tree.make_tree();
      Spearman() {}
      Spearman(const std::vector<T> &x_vals) {
        for (auto &x : x_vals) push_back(x);
      }

      // add a new element
      void push_back(const T &x_val) {
        real_vals.push_back(x_val);
        auto p = _add_value(x_val);
        int z = p.first, dup = p.second;

        // new element with (x=z+dup/2, y=N) -> d1 = z-N+(dup/2)
        sp_d1_type new_d1 = (z-N)*2 + dup; // *= 2
        SNode newelem = SNode(new_d1, (sp_d2_type)new_d1*new_d1, 1);

        if (root == NULL) root = tree.build({newelem});
        else {
          assert(tree.size(root) == N);
          // add +1(*2) to [z+dup, N) & add  +0.5(*2) to [z, z+dup)
          if (z+dup < N) tree.apply(root, z+dup, N, +1*2); // +1(*2) to [z+dup, N)
          if (dup > 0)   tree.apply(root, z, z+dup, +1); // +0.5(*2) to [z, z+dup)
          tree.insert(root, z+dup, newelem); // insert new element @ z+dup
        }
        N += 1;
      }
      // add a new element
      void push_front(const T &x_val) {
        real_vals.push_front(x_val);
        auto p = _add_value(x_val);
        int z = p.first, dup = p.second;

        // new element with (x=z+dup/2, y=0) -> d1 = z-0+(dup/2)
        sp_d1_type new_d1 = (z-0)*2 + dup; // *= 2
        SNode newelem = SNode(new_d1, (sp_d2_type)new_d1*new_d1, 1);

        if (root == NULL) root = tree.build({newelem});
        else {
          assert(tree.size(root) == N);
          if (dup > 0) tree.apply(root, z, z+dup, -1); // [z, z+dup) -= 0.5(*2)
          if (z > 0)   tree.apply(root, 0, z, -1*2); // [0, z) -= 1(*2) (y += 1)
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
        if (z>0)   tree.apply(root, 0, z, +1*2); // [0, z) += 1 (*2)
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
        if (z+dup < N-1) tree.apply(root, z+dup, N-1, -1*2); // [z+dup, )  -= 1(*2)
        N -= 1;
      }
      sp_d2_type spearman_d() const {
        if (root == NULL) return 0;
        return root->sum.d2;
      }
      size_t size() const { return N; }
  };
  /**
   * Online Implementation of Spearman's rank correlation without Binary Search Tree
   * adding/removing a value requires O(N) time complexity and O(N) space complexity
   */
  template< class T >
  class SpearmanLinear : public SpearmanBase<T> {
    int N = 0;
    std::vector<int> D;
    std::deque<T> X_val;
    std::map<T, int> duplicate_counter;
    // internal function: returns the number of existing pairs with the same x-values (>=0)
    inline int _add_value(const T &x_val) {
      int dup = duplicate_counter[x_val]++;
      SpearmanBase<T>::Sx += (sp_d2_type)3*dup*(dup+1);
      return dup;
    }
    // internal function: returns the number of remaining pairs with the same x-values (excluding the pair being erased, >=0)
    inline int _remove_value(const T &x_val) {
      int dup = --duplicate_counter[x_val]; //assert dup>=0
      if (dup == 0) duplicate_counter.erase(x_val);
      SpearmanBase<T>::Sx -= (sp_d2_type)3*dup*(dup+1);
      return dup;
    }
    public:
      SpearmanLinear() {}
      SpearmanLinear(const std::vector<T> &x_vals) {
        for (auto &x : x_vals) push_back(x);
      }
      // add a new element
      void push_back(const T &x_val) {
        int dup = _add_value(x_val);
        int z = 0;
        for (T &x : X_val) if (x < x_val) z++;
        X_val.push_back(x_val);
        D.push_back(0);
        for (int i=N-1; i>=z+dup; i--) D[i+1] = D[i] + 2; // [z+dup, N) += 1(*2)
        for (int i=z+dup-1; i>=z; i--) D[i] = D[i] + 1;   // [z, z+dup) += 0.5(*2) (& insert@z+dup)
        D[z+dup] = (z - N)*2 + dup; // D_i := X_i - Y_i (*2)
        N += 1;
      }
      void push_front(const T &x_val) {
        int dup = _add_value(x_val);
        int z = 0;
        for (T &x : X_val) if (x < x_val) z++;
        X_val.push_front(x_val);
        D.push_back(0);
        for (int i=N-1; i>=z+dup; i--) D[i+1] = D[i];     // [z+dup, N) += 0 (*2) (x += 1, y += 1)
        for (int i=z+dup-1; i>=z; i--) D[i+1] = D[i] - 1; // [z, z+dup) += 0.5-1(*2) (& insert@z)
        for (int i=z-1; i>=0; i--)     D[i] -= 2;         // [0, z) -= 1(*2) (y += 1)
        D[z] = (z - 0)*2 + dup; // D_i := X_i - Y_i (*2)
        N += 1;
      }
      // remove an oldest element
      void pop_front() {
        T x_val = X_val.front();
        int dup = _remove_value(x_val);
        int z = 0;
        for (T &x : X_val) if (x < x_val) z++;
        for (int i=0; i<z; i++)       D[i] += 2;         // [0, z) += 1 (*2)
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
        for (T &x : X_val) if (x < x_val) z++;
        for (int i=z; i<z+dup; i++)   D[i] -= 1;         // [z, z+dup) -= 0.5(*2) (y-=0, x-=0.5)
        for (int i=z+dup; i<N-1; i++) D[i] = D[i+1] - 2; // erase@z+dup & [z+dup, ) -= 1(*2)
        N -= 1;
        X_val.pop_back();
        D.pop_back();
      }
      sp_d2_type spearman_d() const {
        sp_d2_type d = 0;
        for (int i=0; i<N; i++) d += (sp_d2_type)D[i]*D[i];
        return d;
      }
      size_t size() const { return N; }
  };


  /**
   * Traditional Implementation of Spearman's rank correlation:
   * adding/removing a value requires O(NlogN) time complexity and O(M) space complexity
   */
  template< class T >
  class OfflineSpearman : public SpearmanBase<T> {
    int N = 0;
    std::deque<T> X_val;
    std::map<T, int> duplicate_counter;
    // internal function: returns the number of existing pairs with the same x-values (>=0)
    inline int _add_value(const T &x_val) {
      int dup = duplicate_counter[x_val]++;
      SpearmanBase<T>::Sx += (sp_d2_type)3*dup*(dup+1);
      return dup;
    }
    // internal function: returns the number of remaining pairs with the same x-values (excluding the pair being erased, >=0)
    inline int _remove_value(const T &x_val) {
      int dup = --duplicate_counter[x_val]; //assert dup>=0
      if (dup == 0) duplicate_counter.erase(x_val);
      SpearmanBase<T>::Sx -= (sp_d2_type)3*dup*(dup+1);
      return dup;
    }
    public:
      OfflineSpearman() {}
      OfflineSpearman(const std::vector<T> &x_vals) {
        for (auto &x : x_vals) push_back(x);
      }
      // add a new element
      void push_back(const T &x_val) {
        _add_value(x_val);
        N += 1;
        X_val.push_back(x_val);
      }
      // add a new element
      void push_front(const T &x_val) {
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
      sp_d2_type spearman_d() const {
        int n = X_val.size();
        // convert deque to vector
        std::vector<int> X = convert_array_to_rank(std::vector<T>(X_val.begin(), X_val.end()));
        std::vector<int> Y(n);
        for (int i=0; i<n; i++) Y[i] = (i+1)*2; // *2
        sp_d2_type d = 0;
        for (int i=0; i<n; i++) d += (sp_d2_type)(X[i]-Y[i])*(X[i]-Y[i]);
        return d;
      }
      size_t size() const { return N; }
  };
}
