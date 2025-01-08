#pragma once
#include <cassert>
#include <vector>
#include <deque>
#include <map>
#include <algorithm>
#include <cmath>
#include "fast_corr_base.hpp"

namespace FastCorr {
  namespace OfflineCorr {
    // O(NlogN) efficient offline algorithm (tau-b)
    template< class TX, class TY >
    double kendall_tau(const std::vector< std::pair<TX, TY> > &vals);
    // O(NlogN) efficient offline algorithm (tau-b): wrapper for deque
    template< class TX, class TY >
    double kendall_tau(const std::deque< std::pair<TX, TY> > &vals);

    // Calculating sum[i] t_i*(t_i-1)/2 - this is faster than using std::map
    template< class T >
    kd_n2_type offline_nC2_counter(std::vector<T> xs, bool already_sorted=false) {
      if (!already_sorted) std::sort(xs.begin(), xs.end());
      int ctr_same = 1, n = xs.size();
      kd_n2_type ret = 0;
      for (int i=1; i<n; i++) {
        if (xs[i-1] != xs[i]) ctr_same = 0;
        ret += ctr_same++;
      }
      return ret;
    };
    struct BinaryIndexedTree {
      int n;
      std::vector<int> xs;
      BinaryIndexedTree(int n) : n(n) {
        xs.resize(n+1);
      }
      void add(int i, int v) {
        for (int x=i+1; x<=n; x+=x&-x) xs[x] += v;
      }
      int sum(int i) {
        int s = 0;
        for (int x=i+1; x>0; x-=x&-x) s += xs[x];
        return s;
      }
      int sum(int a, int b) {
        return sum(b) - sum(a-1);
      }
    };
  }

  // ===== Online Algorithms ========================================================== //

  namespace MonotonicOnlineCorr {
    template< class T >
    class KendallBase : public Base<T> {
      public:
        virtual double kendall_tau() const = 0;
        double r() const override { return kendall_tau(); } // r() is an alias for kendall_tau()
    };

    template< class T >
    class Kendall : public KendallBase<T> {
      CountingTree< std::pair<T, int> > ctr_tree;
      std::deque<std::pair<T, T> > vals;
      std::map<T, int> ctr_X;
      int min_y_ctr = 0, max_y_ctr = 0;
      int id_for_tree = 0; // add unique id to allow for duplicate values in CountingTree
                           // assuming maximum number of operations <= INT_MAX
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
          K += ctr_tree.order_of_key(std::make_pair(x_val, -1));
          L += ctr_tree.size() - ctr_tree.order_of_key(std::make_pair(x_val, id_for_tree));
          ctr_tree.insert(std::make_pair(x_val, id_for_tree++));
          _add_value(x_val);
        }
        void pop_front() override {
          auto z = vals.front();
          T x_val = z.first;
          min_y_ctr++;
          // y_i should all decrease by 1, but we can simply ignore that as it won't affect the result
          vals.pop_front();
          int num_lower = ctr_tree.order_of_key(std::make_pair(x_val, -1));
          ctr_tree.erase_kth(num_lower);

          K -= ctr_tree.size() - ctr_tree.order_of_key(std::make_pair(x_val, id_for_tree));
          L -= num_lower;
          _remove_value(x_val);
        }
        void push_front(const T &x_val) override {
          vals.push_front(std::make_pair(x_val, --min_y_ctr));
          K += ctr_tree.size() - ctr_tree.order_of_key(std::make_pair(x_val, id_for_tree));
          L += ctr_tree.order_of_key(std::make_pair(x_val, -1));
          ctr_tree.insert(std::make_pair(x_val, id_for_tree++));
          _add_value(x_val);
        }
        void pop_back() override {
          auto z = vals.back();
          max_y_ctr--;
          T x_val = z.first;
          vals.pop_back();
          int num_lower = ctr_tree.order_of_key(std::make_pair(x_val, -1));
          ctr_tree.erase_kth(num_lower);

          K -= num_lower;
          L -= ctr_tree.size() - ctr_tree.order_of_key(std::make_pair(x_val, id_for_tree));
          _remove_value(x_val);
        }
        double kendall_tau() const override {
          int N = vals.size();
          kd_n2_type n0 = (kd_n2_type)N*(N-1)/2;
          // n2 = 0 because y_i has no duplicate values
          if (n1 == n0) return NAN; // denominator will be 0 on tau-b/c
          return (double)(K-L) / sqrt((double)(n0-n1)*(double)n0); // tau-b (n2=0)
        }
        size_t size() const override { return vals.size(); }
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
          min_y_ctr++;
          vals.pop_front();
        }
        void pop_back() override {
          max_y_ctr--;
          vals.pop_back();
        }
        // O(NlogN) efficient offline algorithm (tau-b)
        double kendall_tau() const override {
          return OfflineCorr::kendall_tau(vals);
          //return OfflineCorr::slow_kendall_tau<T>(vals);
        }
        size_t size() const override { return vals.size(); }
    };
  }
  namespace OnlineCorr {
    // TODO
    template< class TX, class TY >
    class Kendall : public Base<TX, TY> {
      public:
        virtual void add(const TX &x_val, const TY &y_val) = 0;
        virtual void remove(const TX &x_val, const TY &y_val) = 0;
        virtual double r() const = 0;
        virtual size_t size() const = 0;
    };

    template<class TY>
    class KendallOnBoundedX : public Base<int, TY> {
      public:
        const int MAX_X;
        int N = 0;
        KendallOnBoundedX(int MAX_X) : MAX_X(MAX_X) {
        }
        void add(const int &x_val, const TY &y_val) override {
          N++;
        }
        void remove(const int &x_val, const TY &y_val) override {
          N--;
        }
        double r() const override {
          return 0;
        }
        size_t size() const override {
          return N;
        }
    };
  }

  // ===== Offline Algorithms ========================================================= //

  namespace OfflineCorr {
    // internal function
    template< class TX, class TY >
    double internal_kendall_tau_on_sorted_pairs(const std::vector< std::pair<TX, TY> > &sorted);

    // O(NlogN) efficient offline algorithm (tau-b)
    template< class TX, class TY >
    double kendall_tau(const std::vector< std::pair<TX, TY> > &vals) {
      std::vector< std::pair<TX, TY> > sorted(vals.begin(), vals.end());
      std::sort(sorted.begin(), sorted.end()); // O(nlogn)
      return internal_kendall_tau_on_sorted_pairs<TX, TY>(sorted);
    }
    // O(NlogN) efficient offline algorithm (tau-b): wrapper for deque
    template< class TX, class TY >
    double kendall_tau(const std::deque< std::pair<TX, TY> > &vals) {
      std::vector< std::pair<TX, TY> > sorted(vals.begin(), vals.end());
      std::sort(sorted.begin(), sorted.end()); // O(nlogn)
      return internal_kendall_tau_on_sorted_pairs<TX, TY>(sorted);
    }
    // wrapper for two vectors
    template< class TX, class TY >
    double kendall_tau(const std::vector<TX> &x_vals, const std::vector<TY> &y_vals) {
      assert(x_vals.size() == y_vals.size());
      int n = x_vals.size();
      std::vector<std::pair<TX, TY> > vals;
      for (int i=0; i<n; i++) vals.push_back(std::make_pair(x_vals[i], y_vals[i]));
      return kendall_tau(vals);
    }

    // internal function
    template< class TX, class TY >
    double internal_kendall_tau_on_sorted_pairs(const std::vector< std::pair<TX, TY> > &sorted) {
      const int n = sorted.size();
      if (n <= 1) return NAN;
      kd_n2_type K = 0, L = 0;
      const kd_n2_type n0 = (kd_n2_type)n*(n-1)/2;

      std::vector<std::pair<TY, int>> ys(n);
      for (int i=0; i<n; i++) ys[i] = std::make_pair(sorted[i].second, i);
      sort(ys.begin(), ys.end()); // O(nlogn)

      std::vector<int> sorted_cmp(n); // compress TY -> int [0,H)
      int H = 0;
      for (int i=0; i<n; i++) {
        if (i > 0 && ys[i-1].first != ys[i].first) H++;
        sorted_cmp[ys[i].second] = H;
      }
      H++;
      kd_n2_type n1 = 0, n2 = 0;

      BinaryIndexedTree bit(H);
      std::vector<int> ctr(H, 0);
      int head = 0;
      for (int i=0; i<n; i++) {
        if (i > 0 && sorted[i-1].first != sorted[i].first) {
          int c = i-head;
          n1 += (kd_n2_type)c*(c-1)/2;
          for (; head<i; head++) {
            int yi = sorted_cmp[head];
            bit.add(yi, 1); // add y
            ctr[yi]++;
          }
          // head = i
        }
        int yi = sorted_cmp[i], s = bit.sum(yi-1);
        K += s;              // K += #{yj < yi}
        L += head-s-ctr[yi]; // L += #{yj > yi}
      }
      int last_c = n-head;
      n1 += (kd_n2_type)last_c*(last_c-1)/2;
      for (int i=0; i<H; i++) n2 += (kd_n2_type)ctr[i]*(ctr[i]-1);
      n2 /= 2;

      if (n1 == n0 || n2 == n0) return NAN; // denominator will be 0 on tau-b and tau-c
      // return (double)(K-L) / (double)n0; // tau-a
      //for (auto p : vals) { cout<<"("<<p.first<<", "<<p.second<<"),"; } cout<<" -> tau = "<< (double)(K-L) <<"/"<< sqrt((n0-n1)*(n0-n2)) << "\n";
      return (double)(K-L) / sqrt((double)(n0-n1)*(double)(n0-n2)); // tau-b
      //int m = min(ctr_X.size(), ctr_Y.size());
      //return 2.0*(double)(K-L) / (n*n * (double)(m-1) / (double)m); // tau-c
    }

  }
}
