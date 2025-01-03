#pragma once
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
    class Base {
      public:
        virtual void push_front(const T &x_val) = 0;
        virtual void push_back(const T &x_val) = 0;
        virtual void pop_front() = 0;
        virtual void pop_back() = 0;
        virtual size_t size() const = 0;
        virtual double r() const = 0;
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
  }

  /**
   * @brief set of type T with order_of_key(val) and erase_kth(k) implemented in O(logN)
   * Currently this is used both by spearman and kendall
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
}
