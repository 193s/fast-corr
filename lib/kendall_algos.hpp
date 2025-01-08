#pragma once
#include <cassert>
#include <vector>
#include <deque>
#include <map>
#include <algorithm>
#include <cmath>
#include <type_traits>
#include <climits>
#include "fast_corr_base.hpp"

namespace FastCorr {
  namespace OfflineCorr {
    // O(NlogN) efficient offline algorithm (tau-b)
    template< class TX, class TY >
    corr_type kendall_tau(const std::vector< std::pair<TX, TY> > &vals);
    // O(NlogN) efficient offline algorithm (tau-b): wrapper for deque
    template< class TX, class TY >
    corr_type kendall_tau(const std::deque< std::pair<TX, TY> > &vals);

    /* currently this is not used
    // Calculating sum[i] t_i*(t_i-1)/2 - this is faster than using std::map
    template< class T >
    kd_n2_type offline_nC2_counter(std::vector<T> xs, bool already_sorted=false) {
      if (!already_sorted) std::sort(xs.begin(), xs.end());
      int ctr_same = 1, n = xs.size();
      kd_n2_type ret = 0;
      for (int i=1; i<n; ++i) {
        if (xs[i-1] != xs[i]) ctr_same = 0;
        ret += ctr_same++;
      }
      return ret;
    };
    */
    struct BinaryIndexedTree {
      int n;
      std::vector<int> xs;
      BinaryIndexedTree(int n) : n(n), xs(n+1, 0) {}
      void add(int i, int v) {
        for (int x=i+1; x<=n; x+=x&-x) xs[x] += v;
      }
      int sum(int i) const {
        int s = 0;
        for (int x=i+1; x>0; x-=x&-x) s += xs[x];
        return s;
      }
      int sum(int a, int b) const {
        return sum(b) - sum(a-1);
      }
    };
  }

  // ===== Online Algorithms ========================================================== //

  namespace MonotonicOnlineCorr {
    template< class T >
    class KendallBase : public Base<T> {
      public:
        virtual corr_type kendall_tau() const noexcept = 0;
        corr_type r() const noexcept override { return kendall_tau(); } // r() is an alias for kendall_tau()
    };

    template< class T >
    class Kendall : public KendallBase<T> {
      CountingTree< std::pair<T, unsigned int> > ctr_tree;
      std::deque<std::pair<T, T> > vals;
      std::map<T, int> ctr_X;
      int min_y_ctr = 0, max_y_ctr = 0;
      unsigned int id_for_tree = 0u; // add unique id to allow for duplicate values in CountingTree
                                    // this is assuming # of add queries < UINT_MAX
      kd_n2_type K = 0, L = 0, n1 = 0;

      inline void _add_value(const T &x_val) {
        n1 += ctr_X[x_val]++;
      }
      inline void _remove_value(const T &x_val) {
        // equivalent to n1 -= --ctr_X[x_val]
        int tmp = ctr_X[x_val] - 1;
        n1 -= tmp;
        if (tmp == 0) ctr_X.erase(x_val);
        else ctr_X[x_val] = tmp;
      }
      public:
        Kendall() {}
        Kendall(const std::vector<T> &x_vals) {
          for (auto &x : x_vals) push_back(x);
        }
        void push_back(const T &x_val) override {
          vals.push_back(std::make_pair(x_val, max_y_ctr++));
          K += ctr_tree.order_of_key(std::make_pair(x_val, 0u));
          L += ctr_tree.size() - ctr_tree.order_of_key(std::make_pair(x_val, UINT_MAX));
          ctr_tree.insert(std::make_pair(x_val, ++id_for_tree));
          _add_value(x_val);
        }
        void pop_front() override {
          auto z = vals.front();
          T x_val = z.first;
          ++min_y_ctr;
          // y_i should all decrease by 1, but we can simply ignore that as it won't affect the result
          vals.pop_front();
          int num_lower = ctr_tree.order_of_key(std::make_pair(x_val, 0u));
          ctr_tree.erase_kth(num_lower);

          K -= ctr_tree.size() - ctr_tree.order_of_key(std::make_pair(x_val, UINT_MAX));
          L -= num_lower;
          _remove_value(x_val);
        }
        void push_front(const T &x_val) override {
          vals.push_front(std::make_pair(x_val, --min_y_ctr));
          K += ctr_tree.size() - ctr_tree.order_of_key(std::make_pair(x_val, UINT_MAX));
          L += ctr_tree.order_of_key(std::make_pair(x_val, 0u));
          ctr_tree.insert(std::make_pair(x_val, ++id_for_tree));
          _add_value(x_val);
        }
        void pop_back() override {
          auto z = vals.back();
          --max_y_ctr;
          T x_val = z.first;
          vals.pop_back();
          int num_lower = ctr_tree.order_of_key(std::make_pair(x_val, 0u));
          ctr_tree.erase_kth(num_lower);

          K -= num_lower;
          L -= ctr_tree.size() - ctr_tree.order_of_key(std::make_pair(x_val, UINT_MAX));
          _remove_value(x_val);
        }
        corr_type kendall_tau() const noexcept override {
          int N = vals.size();
          kd_n2_type n0 = (kd_n2_type)N*(N-1)/2;
          // n2 = 0 because y_i has no duplicate values
          if (FAST_CORR_UNLIKELY(n1 == n0)) return NAN; // denominator will be 0 on tau-b/c
          return (corr_type)(K-L) / (
              (corr_type)n0 * sqrt((corr_type)1.0-(corr_type)n1/(corr_type)n0)); // tau-b (n2=0)
        }
        size_t size() const noexcept override { return vals.size(); }
    };

    template< class T >
    class OfflineKendallForBenchmark : public KendallBase<T> {
      int min_y_ctr = 0, max_y_ctr = 0;
      std::deque<std::pair<T, T> > vals;

      public:
        OfflineKendallForBenchmark() {}
        OfflineKendallForBenchmark(const std::vector<T> &x_vals) {
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
        // O(NlogN) efficient offline algorithm (tau-b)
        corr_type kendall_tau() const noexcept override {
          return OfflineCorr::kendall_tau(vals);
          //return OfflineCorr::slow_kendall_tau<T>(vals);
        }
        size_t size() const noexcept override { return vals.size(); }
    };
  }
  namespace OnlineCorr {
    // TODO
    template< class TX, class TY >
    class Kendall : public Base<TX, TY> {
      public:
        virtual void add(const TX &x_val, const TY &y_val) = 0;
        virtual void remove(const TX &x_val, const TY &y_val) = 0;
        virtual corr_type r() const = 0;
        virtual size_t size() const = 0;
    };

    /** Online Kendall Algorithm when X_i of added pairs are all positive integers from 0 to MAX_X
     * time complexity is O(log N * log MAX_X)
     * space complexity is O(N log N * log MAX_X) overall
     */
    template<class TX, class TY>
    class KendallOnBoundedX : public Base<TX, TY> {
      static_assert(std::is_integral<TX>::value, "TX must be an integral type (int, long long, etc)");
      public:
        const TX MAX_X;
        int N = 0;
        KendallOnBoundedX(TX MAX_X) : MAX_X(MAX_X) {
        }
        void add(const TX &x_val, const TY &y_val) override {
          assert(0 <= x_val && x_val <= MAX_X);
          ++N;
        }
        void remove(const TX &x_val, const TY &y_val) override {
          assert(0 <= x_val && x_val <= MAX_X);
          --N;
        }
        corr_type r() const noexcept override {
          return 0;
        }
        size_t size() const noexcept override {
          return N;
        }
    };
    /** Online Kendall Algorithm when Y_i of added pairs are all positive integers from 0 to MAX_X
     * time complexity is O(log N * log MAX_X)
     * space complexity is O(N log N * log MAX_X) overall
     * Internally this is equivalent to KendallOnBoundedX, with x and y being exchanged
     */
    template<class TX, class TY>
    class KendallOnBoundedY : public Base<TX, TY> {
      static_assert(std::is_integral<TY>::value, "TY must be an integral type (int, long long, etc)");
      KendallOnBoundedX<TY, TX> algo;
      public:
        KendallOnBoundedY(TY MAX_Y) : algo(MAX_Y) {}
        inline void add   (const TX &x_val, const TY &y_val) override { algo.add   (y_val, x_val); }
        inline void remove(const TX &x_val, const TY &y_val) override { algo.remove(y_val, x_val); }
        inline corr_type r() const noexcept override { return algo.r(); }
        inline size_t size() const noexcept override { return algo.size(); }
    };
  }

  // ===== Offline Algorithms ========================================================= //

  namespace OfflineCorr {
    // internal function
    template< class TX, class TY >
    corr_type internal_kendall_tau_on_sorted_pairs(const std::vector< std::pair<TX, TY> > &sorted);

    // O(NlogN) efficient offline algorithm (tau-b)
    template< class TX, class TY >
    corr_type kendall_tau(const std::vector< std::pair<TX, TY> > &vals) {
      std::vector< std::pair<TX, TY> > sorted(vals.begin(), vals.end());
      std::sort(sorted.begin(), sorted.end()); // O(nlogn)
      return internal_kendall_tau_on_sorted_pairs<TX, TY>(sorted);
    }
    // O(NlogN) efficient offline algorithm (tau-b): wrapper for deque
    template< class TX, class TY >
    corr_type kendall_tau(const std::deque< std::pair<TX, TY> > &vals) {
      std::vector< std::pair<TX, TY> > sorted(vals.begin(), vals.end());
      std::sort(sorted.begin(), sorted.end()); // O(nlogn)
      return internal_kendall_tau_on_sorted_pairs<TX, TY>(sorted);
    }
    // wrapper for two vectors
    template< class TX, class TY >
    corr_type kendall_tau(const std::vector<TX> &x_vals, const std::vector<TY> &y_vals) {
      assert(x_vals.size() == y_vals.size());
      int n = x_vals.size();
      std::vector<std::pair<TX, TY> > vals;
      for (int i=0; i<n; ++i) vals.push_back(std::make_pair(x_vals[i], y_vals[i]));
      return kendall_tau(vals);
    }

    // internal function
    template< class TX, class TY >
    corr_type internal_kendall_tau_on_sorted_pairs(const std::vector< std::pair<TX, TY> > &sorted) {
      const int n = sorted.size();
      if (FAST_CORR_UNLIKELY(n <= 1)) return NAN;
      kd_n2_type K = 0, L = 0;
      const kd_n2_type n0 = (kd_n2_type)n*(n-1)/2;

      std::vector<std::pair<TY, int>> ys(n);
      for (int i=0; i<n; ++i) ys[i] = std::make_pair(sorted[i].second, i);
      sort(ys.begin(), ys.end()); // O(nlogn)

      std::vector<int> sorted_cmp(n); // compress TY -> int [0,H)
      int H = 0;
      for (int i=0; i<n; ++i) {
        if (i > 0 && ys[i-1].first != ys[i].first) ++H;
        sorted_cmp[ys[i].second] = H;
      }
      ++H;
      kd_n2_type n1 = 0, n2 = 0;

      BinaryIndexedTree bit(H);
      std::vector<int> ctr(H, 0);
      int head = 0;
      for (int i=0; i<n; ++i) {
        if (i > 0 && sorted[i-1].first != sorted[i].first) {
          int c = i-head;
          n1 += (kd_n2_type)c*(c-1)/2;
          for (; head<i; ++head) {
            int yi = sorted_cmp[head];
            bit.add(yi, 1); // add y
            ++ctr[yi];
          }
          // head = i
        }
        int yi = sorted_cmp[i], s = bit.sum(yi-1);
        K += s;              // K += #{yj < yi}
        L += head-s-ctr[yi]; // L += #{yj > yi}
      }
      int last_c = n-head;
      n1 += (kd_n2_type)last_c*(last_c-1)/2;
      for (int i=0; i<H; ++i) n2 += (kd_n2_type)ctr[i]*(ctr[i]-1);
      n2 /= 2;

      if (FAST_CORR_UNLIKELY(n1 == n0 || n2 == n0)) return NAN; // denominator will be 0 on tau-b and tau-c
      // return (double)(K-L) / (double)n0; // tau-a
      //for (auto p : vals) { cout<<"("<<p.first<<", "<<p.second<<"),"; } cout<<" -> tau = "<< (double)(K-L) <<"/"<< sqrt((n0-n1)*(n0-n2)) << "\n";
      return (corr_type)(K-L) / ((corr_type)n0 * sqrt(
             ((corr_type)1.0-(corr_type)n1/(corr_type)n0)
            *((corr_type)1.0-(corr_type)n2/(corr_type)n0))); // tau-b
      //int m = min(ctr_X.size(), ctr_Y.size());
      //return 2.0*(double)(K-L) / (n*n * (double)(m-1) / (double)m); // tau-c
    }

  }
}
