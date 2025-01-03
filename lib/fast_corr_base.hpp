#pragma once
#include <cassert>
#include <vector>
#include <map>
#include <deque>
#include <algorithm>
#include <cmath>
// CountingTree<T> is currently dependent on G++ extensions
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
#include <ext/pb_ds/tag_and_trait.hpp>

namespace FastCorr {
  // -- type definitions for spearman: change accordingly --
  typedef long long sp_d1_type; // sum(d_i) can be negative for some intervals, and is at max N^2
  typedef unsigned long long sp_d2_type; // sum(d_i^2) (>=0) can be as large as (N^3-N)/6 <=> rho=0
    // d2_type is used for Sx as well: Sx = sum(t_i^3 - t_i) (>=0) can be as large as N^3-N
    // under N <= 2642245, d2_type won't exceed ULLONG_MAX=2^64-1

  // -- type definitions for kendall: change accordingly --
  typedef long long kd_n2_type; // data type used to store K,L,N*(N-1)/2 >= 0
    // under N <= 4294967296 (aprox. 4*10^9), kd_n2_type won't exceed LONG_MAX=2^63-1
    // also assuming maximum number of operations <= INT_MAX by "int id_for_tree"
  // typedef long corr_result_type; // type used for correlation coefficients ([-1, 1] or NAN)
  //                                // options are: double, float, ...

  double spearman_r_from_n_4d_and_Sx(int n, sp_d2_type d, sp_d2_type Sx);
  /**
   * @brief A module for online correlation algorithms on partial monotonicity constraints
   */
  namespace MonotonicOnlineCorr {
    template< class T >
    class KendallBase {
      public:
        virtual void push_front(const T &x_val) = 0;
        virtual void push_back(const T &x_val) = 0;
        virtual void pop_front() = 0;
        virtual void pop_back() = 0;
        virtual double kendall_tau() const = 0;
        virtual size_t size() const = 0;
        double r() const { return kendall_tau(); } // r() is an alias for kendall_tau()
    };

    template< class T >
    class SpearmanBase {
      public:
        virtual void push_back(const T &x_val) = 0;
        virtual void push_front(const T &x_val) = 0;
        virtual void pop_front() = 0;
        virtual void pop_back() = 0;
        virtual sp_d2_type spearman_d() const = 0;
        virtual size_t size() const = 0;
        double spearman_r() const {
          sp_d2_type d = spearman_d();
          // d = sum((2d_i)^2) = 4*actual_D
          int n = size();
          return spearman_r_from_n_4d_and_Sx(n, d, Sx);
        }
        double r() { return spearman_r(); } // r() is an alias for spearman_r()
      protected:
        sp_d2_type Sx = 0; // sum(t_i^3 - t_i)
        // Sy = 0 under monotonic constraints
    };
  }

  /**
   * @brief A module for online correlation algorithms (with no constraints)
   */
  namespace OnlineCorr {
    template< class TX, class TY >
    class Base {
      public:
        virtual void add(const TX &x_val, const TY &y_val) = 0;
        virtual void remove(const TX &x_val, const TY &y_val) = 0;
        virtual double r() const = 0;
        virtual size_t size() const = 0;
    };

    template< class TX, class TY >
    class Pearson : public Base<TX, TY> {
      public:
        virtual void add(const TX &x_val, const TY &y_val) = 0;
        virtual void remove(const TX &x_val, const TY &y_val) = 0;
        virtual double r() const = 0;
        virtual size_t size() const = 0;
    };

    template< class TX, class TY >
    class Kendall : public Base<TX, TY> {
      public:
        virtual void add(const TX &x_val, const TY &y_val) = 0;
        virtual void remove(const TX &x_val, const TY &y_val) = 0;
        virtual double r() const = 0;
        virtual size_t size() const = 0;
    };
  }

  /**
   * CountingTree<T>: set of type T with order_of_key(val) and erase_kth(k) implemented in O(logN)
   */
  template< class T >
  class CountingTree {
    using GPPTree = __gnu_pbds::tree< T, __gnu_pbds::null_type, std::less<T>,
    __gnu_pbds::rb_tree_tag, __gnu_pbds::tree_order_statistics_node_update>;
    GPPTree gpp_tree;

    public:
      void insert(const T &val) {
        gpp_tree.insert(val);
      }
      int order_of_key(const T &val) const {
        return gpp_tree.order_of_key(val);
      }
      void erase_kth(const int &k) {
        gpp_tree.erase(gpp_tree.find_by_order(k));
      }
      size_t size() const {
        return gpp_tree.size();
      }
  };

  /**
   * returns the ranks of given vector of any type
   * ranks are multiplied by 2 so that all ranks will be integers
   * e.g. [3, 12123, 0] -> [2, 3, 1]*2
   * e.g. [1, 2, 2, 2, 5, 5, 7] -> [1, 3, 3, 3, 5.5, 5.5, 7]*2
   * time complexity: O(NlogN)
   */
  template< class T >
  std::vector<int> convert_array_to_rank(const std::vector<T> &X) {
    int n = X.size();
    if (n == 0) return {};
    if (n == 1) return {2};
    // n>=2
    std::vector<std::pair<T, int> > X2(n);
    for (int i=0; i<n; i++) X2[i] = std::pair<T, int>(X[i], i);
    std::sort(X2.begin(), X2.end()); // O(nlogn)
    std::vector<int> ret(n);
    //for (int i=0; i<n; i++) ret[X2[i].second] = 2*(i+1); // works only on unique arrays
    int z = 0, head = 0;
    for (int i=1; i<n; i++) {
      if (X2[i-1].first != X2[i].first) {
        // finalize rank
        int rank = z*2 + 1 + (i-head); // 1 + z + (number of same values - 1)/2
        z += i-head;
        for (; head<i; head++) ret[X2[head].second] = rank;
      }
    }
    int last_rank = z*2 + 1 + (n-head); // 1 + z + (number of same values - 1)/2
    for (; head<n; head++) ret[X2[head].second] = last_rank;
    return ret;
  }

  namespace OfflineCorr {
    // O(NlogN) offline algorithm (spearman-r)
    template< class T >
    double spearman_r(const std::vector< std::pair<T, T> > &vals);
    template< class T >
    double spearman_r(const std::vector<T> &x_vals, const std::vector<T> &y_vals);

    // O(NlogN) efficient offline algorithm (tau-b)
    template< class T >
    double kendall_tau(const std::vector< std::pair<T, T> > &vals);
    // O(NlogN) efficient offline algorithm (tau-b): wrapper for deque
    template< class T >
    double kendall_tau(const std::deque< std::pair<T, T> > &vals);

    /**
     * straightforward pearson implementation: O(N) time complexity
     */
    double straightforward_pearson_r(const std::vector<double> &X, const std::vector<double> &Y) {
      assert(X.size() == Y.size());
      int n = X.size();
      double sum_X = 0, sum_Y = 0, sum_XY = 0;
      double sum_X2 = 0, sum_Y2 = 0;

      for (int i=0; i<n; i++) {
          sum_X += X[i];
          sum_Y += Y[i];
          sum_XY += X[i]*Y[i];
          sum_X2 += X[i]*X[i];
          sum_Y2 += Y[i]*Y[i];
      }
      return (double)(n*sum_XY - sum_X*sum_Y)
              / sqrt((n*sum_X2 - sum_X*sum_X) * (n*sum_Y2 - sum_Y*sum_Y));
    }
  }
}
