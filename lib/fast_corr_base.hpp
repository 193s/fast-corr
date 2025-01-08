#pragma once

#include <type_traits>
#if __has_include(<ext/pb_ds/assoc_container.hpp>) // G++ environment
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
#include <ext/pb_ds/tag_and_trait.hpp>
#else // non-G++ environment
#include "rbst_base.hpp"
#endif

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

  // -- type definitions for correlation results --
  typedef double corr_type; // data type used to store correlation coefficients
                            // options are: double, float, long double
                            // (std::sqrt needs to be defined on corr_type)
  static_assert(std::is_floating_point<corr_type>::value, "corr_type must be a floating point type");

  corr_type spearman_r_from_n_4d_and_Sx(int n, sp_d2_type d, sp_d2_type Sx);
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
        virtual corr_type r() const = 0;
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
        virtual corr_type r() const = 0;
        virtual size_t size() const = 0;
    };
  }


  /**
   * @brief set of type T with order_of_key(val) and erase_kth(k) implemented in O(logN)
   * Currently this is used both by spearman and kendall
   */
  #if __has_include(<ext/pb_ds/assoc_container.hpp>) // G++ environment
  // Still using the G++ extension tree if it's available, as currently this seems to be faster
  template< class T >
  class CountingTree {
    using GPPTree = __gnu_pbds::tree< T, __gnu_pbds::null_type, std::less<T>,
    __gnu_pbds::rb_tree_tag, __gnu_pbds::tree_order_statistics_node_update>;
    GPPTree gpp_tree;

    public:
      int order_of_key(const T &val) const {
        return gpp_tree.order_of_key(val);
      }
      void insert(const T &val) {
        gpp_tree.insert(val);
      }
      void erase_kth(const int &k) {
        gpp_tree.erase(gpp_tree.find_by_order(k));
      }
      size_t size() const {
        return gpp_tree.size();
      }
  };

  #else // non G++ environment
  /**
   * Internal node class defined for Spearman algorithm
   */
  template <typename T>
  struct CountingNode {
    typename RBSTBase<CountingNode<T>>::Ptr l, r;
    T key;
    int cnt;
    CountingNode(const T &t = T()) : l(), r(), key(t), cnt(1) {}
  };

  template< class T >
  struct InternalCountingTree : RBSTBase<CountingNode<T>> {
    using Node = CountingNode<T>;
    using base = RBSTBase<CountingNode<T>>;
    using base::merge;
    using base::split;
    using base::count;
    using typename base::Ptr;

    InternalCountingTree() = default;

    T fold(Ptr &t, int a, int b) {
      auto x = split(t, a);
      auto y = split(x.second, b - a);
      auto ret = count(y.first);
      t = merge(x.first, merge(y.first, y.second));
      return ret;
    }
    std::pair<Ptr, Ptr> split_by_key(Ptr t, const T &key) {
      if (!t) return {nullptr, nullptr};
      push(t);
      if (key < t->key) { // assuming no duplicates
        auto s = split_by_key(t->l, key);
        t->l = s.second;
        return {s.first, update(t)};
      } else {
        auto s = split_by_key(t->r, key);
        t->r = s.first;
        return {update(t), s.second};
      }
    }
    template <typename... Args>
    void insert_by_key(Ptr &t, const T &key, const Args &... args) {
      auto x = split_by_key(t, key);
      t = merge(merge(x.first, my_new(args...)), x.second);
    }
    int order_of_key(Ptr t, const T &val) const {
      if (!t) return 0;
      if      (val >  t->key) return 1 + count(t->l) + order_of_key(t->r, val);
      else if (val == t->key) return     count(t->l);
      else return order_of_key(t->l, val);
    }

   protected:
    Ptr update(Ptr t) override {
      t->cnt = 1;
      if (t->l) t->cnt += t->l->cnt;
      if (t->r) t->cnt += t->r->cnt;
      return t;
    }
    void push(Ptr) override {};
  };
  template< class T >
  class CountingTree {
    InternalCountingTree<T> ctr_tree;
    typename InternalCountingTree<T>::Node *root = ctr_tree.make_tree();

    public:
      int order_of_key(const T &val) const {
        return ctr_tree.order_of_key(root, val);
      }
      void insert(const T &val) {
        ctr_tree.insert(root, order_of_key(val), val);
      }
      void erase_kth(const int &k) {
        ctr_tree.erase(root, k);
      }
      size_t size() const {
        return ctr_tree.size(root);
      }
  };
  #endif
}
