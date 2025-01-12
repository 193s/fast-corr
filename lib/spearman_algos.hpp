#pragma once
#include <cassert>
#include <vector>
#include <deque>
#include <map>
#include <algorithm>
#include <cmath>
#include <climits>
#include "fast_corr_base.hpp"
#include "lazy_rbst.hpp"

namespace FastCorr {
  namespace OfflineCorr {
    // O(NlogN) offline algorithm (spearman-r)
    template< class TX, class TY >
    corr_type spearman_r(const std::vector< std::pair<TX, TY> > &vals);
    template< class TX, class TY >
    corr_type spearman_r(const std::vector<TX> &x_vals, const std::vector<TY> &y_vals);
    template< class TX, class TY >
    corr_type spearman_r(const std::deque< std::pair<TX, TY> > &vals);

    template< class TX, class TY >
    sp_d2_type spearman_d(const std::vector<TX> &x_vals, const std::vector<TY> &y_vals);
    template< class TX, class TY >
    sp_d2_type spearman_d(const std::vector< std::pair<TX, TY> > &vals);
    template< class TX, class TY >
    sp_d2_type spearman_d(const std::deque< std::pair<TX, TY> > &vals);

    // Calculating sum[i] t_i*(t_i^2-1) - this is faster than using std::map
    template< class T >
    sp_d2_type offline_n3_counter(std::vector<T> xs) {
      std::sort(xs.begin(), xs.end());
      int ctr_same = 1, n = xs.size();
      sp_d2_type ret = 0;
      for (int i=1; i<n; ++i) {
        if (xs[i-1] != xs[i]) ctr_same = 0;
        ret += (sp_d2_type)3*ctr_same*(ctr_same+1);
        ++ctr_same;
      }
      return ret;
    };
  }
}
// ===== Helper functions =========================================================== //

namespace FastCorr {
  corr_type spearman_r_from_n_4d_Gx_and_Gy(int n, sp_d2_type d, sp_d2_type Gx, sp_d2_type Gy) {
    sp_d2_type n3 = (sp_d2_type)n*((sp_d2_type)n*n-1);
    if (FAST_CORR_UNLIKELY(Gx == n3 || Gy == n3)) return corr_NAN; // rank X_i can not be defined in this case
    else if (Gx == 0 && Gy == 0)
      return (corr_type)1.0 - (corr_type)1.5*(corr_type)d / (corr_type)n3;
    else {
      // general formula:
      // (Sx + Sy - D) / (2 * sqrt(Sx) * sqrt(Sy)), Sx := (n3-Gx)/12, Sy := (n3-Gy)/12
      // = (12*Sx + 12*Sy - 12*D) / (2 * sqrt(n3-Gx) * sqrt(n3-Gy))
      // = (2*n3 - Gx - Gy - 12*D) / (2 * sqrt(n3-Gx) * sqrt(n3-Gy))
      // = (2*n3 - Gx - Gy - 3*(4D)) / (2 * n3 * sqrt(1-Gx/n3) * sqrt(1-Gy/n3))
      return ((corr_type)2.0*(corr_type)n3 - (corr_type)Gx - (corr_type)Gy - (corr_type)3.0*d) / (
          (corr_type)2.0 * (corr_type)n3
          * sqrt((corr_type)1.0 - (corr_type)Gx/(corr_type)n3)
          * sqrt((corr_type)1.0 - (corr_type)Gy/(corr_type)n3));
    }
  }
  // special case when Gy = 0
  corr_type spearman_r_from_n_4d_and_Gx(int n, sp_d2_type d, sp_d2_type Gx) {
    sp_d2_type n3 = (sp_d2_type)n*((sp_d2_type)n*n-1);
    if  (FAST_CORR_UNLIKELY(Gx == n3)) return corr_NAN; // rank X_i can not be defined in this case
    else if (Gx == 0) return (corr_type)1.0 - (corr_type)1.5*(corr_type)d / (corr_type)n3;
    else {
      // when Gy = 0,
      // = (2*n3 - 3*(4D) - Gx - Gy) / (2*n3*(sqrt(1 - Gx/n3) * sqrt(1-Gy/n3)))
      // = (2*n3 - 3*(4D) - Gx) / (2*n3*(sqrt(1 - Gx/n3))
      return ((corr_type)2.0*(corr_type)n3 - (corr_type)Gx - (corr_type)3.0*d) / (
          (corr_type)2.0 * (corr_type)n3
          * sqrt((corr_type)1.0 - (corr_type)Gx/(corr_type)n3));
    }
  }

  /**
   * returns the ranks of given vector of any type
   * ranks are multiplied by 2 so that all ranks will be integers
   * e.g. [3, 12123, 0] -> [2, 3, 1]*2
   * e.g. [1, 2, 2, 2, 5, 5, 7] -> [1, 3, 3, 3, 5.5, 5.5, 7]*2
   * time complexity: O(NlogN)
   */
  template< class T >
  std::vector<int> rankdata(const std::vector<T> &X) {
    int n = X.size();
    if (FAST_CORR_UNLIKELY(n == 0)) return {};
    if (FAST_CORR_UNLIKELY(n == 1)) return {2};
    // n>=2
    std::vector<std::pair<T, int> > X2(n);
    for (int i=0; i<n; ++i) X2[i] = std::pair<T, int>(X[i], i);
    std::sort(X2.begin(), X2.end(), // O(nlogn)
        [&](std::pair<T, int> i, std::pair<T, int> j) {
          return i.first < j.first; // this is faster: order of the second values does not matter
        });
    std::vector<int> ret(n);
    //for (int i=0; i<n; ++i) ret[X2[i].second] = 2*(i+1); // works only on unique arrays
    int z = 0, head = 0;
    for (int i=1; i<n; ++i) {
      if (X2[i-1].first != X2[i].first) {
        // finalize rank
        int rank = z*2 + 1 + (i-head); // 1 + z + (number of same values - 1)/2
        z += i-head;
        for (; head<i; ++head) ret[X2[head].second] = rank;
      }
    }
    int last_rank = z*2 + 1 + (n-head); // 1 + z + (number of same values - 1)/2
    for (; head<n; ++head) ret[X2[head].second] = last_rank;
    return ret;
  }

  // ===== Online Algorithms ========================================================== //
  namespace MonotonicOnlineCorr {
    template< class T >
    class SpearmanBase : public Base<T> {
      public:
        virtual sp_d2_type spearman_d() const noexcept = 0;
        virtual size_t size() const noexcept = 0;
        corr_type spearman_r() const noexcept {
          sp_d2_type d = spearman_d(); // d = sum[i=1..n]((2d_i)^2) = 4*actual_D
          return spearman_r_from_n_4d_and_Gx(size(), d, Gx);
        }
        corr_type r() const noexcept override { return spearman_r(); } // r() is an alias for spearman_r()
      protected:
        mutable sp_d2_type Gx = 0; // sum(t_i^3 - t_i)
        // Gy = 0 under monotonic constraints
    };
    //  < class D, class L, D (*f)(D, D), D (*g)(D, L), L (*h)(L, L), L (*p)(L, int) >

    /**
     * Efficient Online Implementation of Spearman's rank correlation with Binary Search Tree
     * adding/removing a value requires O(logN) time complexity and O(N) space complexity
     */
    template< class T >
    class Spearman : public SpearmanBase<T> {
      D2LazyRBST tree;
      CountingTree< std::pair<T, unsigned int> > ctr_tree;
      int N = 0;
      unsigned int id_for_tree = 0u; // add unique ids to allow for duplicate values in CountingTree
                                     // this is assuming # of add query < UINT_MAX
      std::deque<T> real_vals;
      D2LazyRBST::Node *root = tree.make_tree();
      inline std::pair<int, int> _add_value(const T &x_val) {
        int z = ctr_tree.order_of_key(std::make_pair(x_val, 0u)); // # of < x_val
        int dup = ctr_tree.order_of_key(std::make_pair(x_val, UINT_MAX)) - z;
        SpearmanBase<T>::Gx += (sp_d2_type)3*dup*(dup+1);
        ctr_tree.insert(std::make_pair(x_val, ++id_for_tree));
        return std::make_pair(z, dup);
      }
      inline std::pair<int, int> _remove_value(const T &x_val) {
        int z = ctr_tree.order_of_key(std::make_pair(x_val, 0u));
        int dup = ctr_tree.order_of_key(std::make_pair(x_val, UINT_MAX)) - z - 1;
        ctr_tree.erase_kth(z); // erase @z (@z+dup) will also work
        // ctr_tree.erase(std::make_pair(x_val, <id to erase>));
        assert(dup >= 0);
        SpearmanBase<T>::Gx -= (sp_d2_type)3*dup*(dup+1);
        return std::make_pair(z, dup);
      }

      public:
        Spearman() {}
        Spearman(const std::vector<T> &x_vals) {
          for (auto &x : x_vals) push_back(x);
        }

        // add a new element
        void push_back(const T &x_val) override {
          real_vals.push_back(x_val);
          auto p = _add_value(x_val);
          int z = p.first, dup = p.second;

          //assert(tree.size(root) == N);
          // [0, z) | [z, z+dup)  |  insert  | [z+dup, N)
          //   +0   |     +1      | (new_d1) |     +2
          auto pairs = tree.split3(root, z, z+dup);
          auto p1 = pairs.first, p2 = pairs.second.first, p3 = pairs.second.second;
          if (p2) tree.propagate_const<+1>(p2);
          if (p3) tree.propagate_const<+2>(p3);
          // new element with (x=z+dup/2, y=N) -> d1 = z-N+(dup/2)
          sp_d1_type new_d1 = (z-N)*2 + dup; // *= 2
          auto p_new = tree.my_new(new_d1, (sp_d2_type)new_d1*new_d1);
          root = tree.merge4(p1, p2, p_new, p3);
          N += 1;
        }
        // add a new element
        void push_front(const T &x_val) override {
          real_vals.push_front(x_val);
          auto p = _add_value(x_val);
          int z = p.first, dup = p.second;

          // [0, z) |  insert  | [z, z+dup) | [z+dup, N)
          //   -2   | (new_d1) |     -1     |     +0
          auto pairs = tree.split3(root, z, z+dup);
          auto p1 = pairs.first, p2 = pairs.second.first, p3 = pairs.second.second;
          if (p1) tree.propagate_const<-2>(p1);
          if (p2) tree.propagate_const<-1>(p2);
          // new element with (x=z+dup/2, y=0) -> d1 = z-0+(dup/2)
          sp_d1_type new_d1 = (z-0)*2 + dup; // *= 2
          auto p_new = tree.my_new(new_d1, (sp_d2_type)new_d1*new_d1);
          root = tree.merge4(p1, p_new, p2, p3);
          N += 1;
        }
        // remove an oldest element
        void pop_front() override {
          T x_val = real_vals.front();
          real_vals.pop_front();
          auto p = _remove_value(x_val);
          int z = p.first, dup = p.second;

          // [0, z) |  [z]  | [z+1, z+dup+1) | [z+dup+1, N)
          //   +2   | erase |        +1      |      +0
          auto pairs = tree.split3(root, z, z+dup+1);
          auto p1 = pairs.first, p2 = pairs.second.first, p3 = pairs.second.second;
          auto pp = tree.split_by_first_element(p2); p2 = pp.second;
          tree.my_del(pp.first); // erase @z
          if (p1) tree.propagate_const<+2>(p1);
          if (p2) tree.propagate_const<+1>(p2); // y-=1, x-=0.5
          root = tree.merge3(p1, p2, p3);
          N -= 1;
        }
        // remove an oldest element
        void pop_back() override {
          T x_val = real_vals.back();
          real_vals.pop_back();
          auto p = _remove_value(x_val);
          int z = p.first, dup = p.second;

          // [0, z) | [z, z+dup) | [z+dup] | [z+dup+1, N)
          //   +0   |     -1     |  erase  |      +0
          auto pairs = tree.split3(root, z, z+dup);
          auto p1 = pairs.first, p2 = pairs.second.first, p3 = pairs.second.second;
          auto pp = tree.split_by_first_element(p3); p3 = pp.second;
          tree.my_del(pp.first); // erase @z+dup
          if (p2) tree.propagate_const<-1>(p2); // (y-=0, x-=0.5)
          if (p3) tree.propagate_const<-2>(p3);
          root = tree.merge3(p1, p2, p3);
          N -= 1;
        }
        sp_d2_type spearman_d() const noexcept override { return tree.sum_d2(root); }
        size_t size() const noexcept override { return N; }
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
        SpearmanBase<T>::Gx += (sp_d2_type)3*dup*(dup+1);
        return dup;
      }
      // internal function: returns the number of remaining pairs with the same x-values (excluding the pair being erased, >=0)
      inline int _remove_value(const T &x_val) {
        int dup = --duplicate_counter[x_val]; //assert dup>=0
        if (dup == 0) duplicate_counter.erase(x_val);
        SpearmanBase<T>::Gx -= (sp_d2_type)3*dup*(dup+1);
        return dup;
      }
      public:
        SpearmanLinear() {}
        SpearmanLinear(const std::vector<T> &x_vals) {
          D.reserve(N); // optional
          for (auto &x : x_vals) push_back(x);
        }
        // add a new element
        void push_back(const T &x_val) override {
          int dup = _add_value(x_val);
          int z = 0;
          for (T &x : X_val) if (x < x_val) ++z;
          X_val.push_back(x_val);
          D.emplace_back(0);
          for (int i=N-1; i>=z+dup; --i) D[i+1] = D[i] + 2; // [z+dup, N) += 1(*2)
          for (int i=z+dup-1; i>=z; --i) D[i] = D[i] + 1;   // [z, z+dup) += 0.5(*2) (& insert@z+dup)
          D[z+dup] = (z - N)*2 + dup; // D_i := X_i - Y_i (*2)
          N += 1;
        }
        void push_front(const T &x_val) override {
          int dup = _add_value(x_val);
          int z = 0;
          for (T &x : X_val) if (x < x_val) ++z;
          X_val.push_front(x_val);
          D.emplace_back(0);
          for (int i=N-1; i>=z+dup; --i) D[i+1] = D[i];     // [z+dup, N) += 0 (*2) (x += 1, y += 1)
          for (int i=z+dup-1; i>=z; --i) D[i+1] = D[i] - 1; // [z, z+dup) += 0.5-1(*2) (& insert@z)
          for (int i=z-1; i>=0; --i)     D[i] -= 2;         // [0, z) -= 1(*2) (y += 1)
          D[z] = (z - 0)*2 + dup; // D_i := X_i - Y_i (*2)
          N += 1;
        }
        // remove an oldest element
        void pop_front() override {
          T x_val = X_val.front();
          int dup = _remove_value(x_val);
          int z = 0;
          for (T &x : X_val) if (x < x_val) ++z;
          for (int i=0; i<z; ++i)       D[i] += 2;         // [0, z) += 1 (*2)
          for (int i=z; i<z+dup; ++i)   D[i] = D[i+1] + 1; // erase @z & [z, z+dup) += 0.5(*2)
          for (int i=z+dup; i<N-1; ++i) D[i] = D[i+1];
          N -= 1;
          X_val.pop_front();
          D.pop_back();
        }
        void pop_back() override {
          T x_val = X_val.back();
          int dup = _remove_value(x_val);
          int z = 0;
          for (T &x : X_val) if (x < x_val) ++z;
          for (int i=z; i<z+dup; ++i)   D[i] -= 1;         // [z, z+dup) -= 0.5(*2) (y-=0, x-=0.5)
          for (int i=z+dup; i<N-1; ++i) D[i] = D[i+1] - 2; // erase@z+dup & [z+dup, ) -= 1(*2)
          N -= 1;
          X_val.pop_back();
          D.pop_back();
        }
        sp_d2_type spearman_d() const noexcept override {
          sp_d2_type d = 0;
          for (int i=0; i<N; ++i) d += (sp_d2_type)D[i]*D[i];
          return d;
        }
        size_t size() const noexcept override { return N; }
    };


    /**
     * Traditional Implementation of Spearman's rank correlation:
     * adding/removing a value requires O(NlogN) time complexity and O(M) space complexity
     */
    template< class T >
    class OfflineSpearmanForBenchmark : public SpearmanBase<T> {
      int min_y_ctr = 0, max_y_ctr = 0;
      std::deque<std::pair<T, int> > vals;

      public:
        OfflineSpearmanForBenchmark() {}
        OfflineSpearmanForBenchmark(const std::vector<T> &x_vals) {
          for (auto &x : x_vals) push_back(x);
        }
        void push_back(const T &x_val) override {
          vals.push_back(std::make_pair(x_val, max_y_ctr++));
        }
        void push_front(const T &x_val) override {
          vals.push_front(std::make_pair(x_val, --min_y_ctr));
        }
        void pop_front() override {
          // y_i should all decrease by 1, but we can simply ignore that as it won't affect the result
          ++min_y_ctr;
          vals.pop_front();
        }
        void pop_back() override {
          --max_y_ctr;
          vals.pop_back();
        }
        size_t size() const noexcept override { return vals.size(); }
        sp_d2_type spearman_d() const noexcept override {
          int n = vals.size();
          std::vector<T> xs(n);
          for (int i=0; i<n; ++i) xs[i] = vals[i].first;
          SpearmanBase<T>::Gx = OfflineCorr::offline_n3_counter(xs);
          return OfflineCorr::spearman_d(vals);
        }
    };
  }

  // ===== Offline Algorithms ========================================================= //

  namespace OfflineCorr {
    // O(NlogN) spearman algorithm
    template< class TX, class TY >
    corr_type spearman_r(const std::vector<TX> &x_vals, const std::vector<TY> &y_vals) {
      int n = x_vals.size();
      sp_d2_type d = spearman_d(x_vals, y_vals);

      // calculate Gx and Gy
      sp_d2_type Gx = OfflineCorr::offline_n3_counter(x_vals);
      sp_d2_type Gy = OfflineCorr::offline_n3_counter(y_vals);
      return spearman_r_from_n_4d_Gx_and_Gy(n, d, Gx, Gy);
    }
    // wrapper for list of pairs
    template< class TX, class TY >
    corr_type spearman_r(const std::vector< std::pair<TX, TY> > &vals) {
      int n = vals.size();
      std::vector<TX> x_vals(n);
      std::vector<TY> y_vals(n);
      for (int i=0; i<n; ++i) x_vals[i] = vals[i].first;
      for (int i=0; i<n; ++i) y_vals[i] = vals[i].second;
      return spearman_r(x_vals, y_vals);
    }
    // O(NlogN) offline spearman: wrapper for deque
    template< class TX, class TY >
    corr_type spearman_r(const std::deque< std::pair<TX, TY> > &vals) {
      std::vector< std::pair<TX, TY> > vals2(vals.begin(), vals.end());
      return spearman_r(vals2);
    }

    // spearman_d (used by spearman_r)
    template< class TX, class TY >
    sp_d2_type spearman_d(const std::vector<TX> &x_vals, const std::vector<TY> &y_vals) {
      int n = x_vals.size();
      // calculate ranks and 4d
      std::vector<int> X = rankdata(x_vals);
      std::vector<int> Y = rankdata(y_vals);
      sp_d2_type d = 0;
      for (int i=0; i<n; ++i) d += (sp_d2_type)(X[i]-Y[i])*(X[i]-Y[i]);
      return d;
    }
    // wrapper for list of pairs
    template< class TX, class TY >
    sp_d2_type spearman_d(const std::vector< std::pair<TX, TY> > &vals) {
      int n = vals.size();
      std::vector<TX> x_vals(n);
      std::vector<TY> y_vals(n);
      for (int i=0; i<n; ++i) x_vals[i] = vals[i].first;
      for (int i=0; i<n; ++i) y_vals[i] = vals[i].second;
      return spearman_d(x_vals, y_vals);
    }
    // wrapper for deque
    template< class TX, class TY >
    sp_d2_type spearman_d(const std::deque< std::pair<TX, TY> > &vals) {
      std::vector< std::pair<TX, TY> > vals2(vals.begin(), vals.end());
      return spearman_d(vals2);
    }
  }
}
