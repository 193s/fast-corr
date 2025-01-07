#pragma once
#include <cassert>
#include <vector>
#include <deque>
#include <map>
#include <algorithm>
#include <cmath>
#include "fast_corr_base.hpp"

namespace FastCorr::OfflineCorr {
  // O(NlogN) efficient offline algorithm (tau-b)
  template< class TX, class TY >
  double kendall_tau(const std::vector< std::pair<TX, TY> > &vals);
  // O(NlogN) efficient offline algorithm (tau-b): wrapper for deque
  template< class TX, class TY >
  double kendall_tau(const std::deque< std::pair<TX, TY> > &vals);

  // Calculating sum[i] t_i*(t_i-1)/2 - this is faster than using std::map
  template< class T >
  kd_n2_type offline_nC2_counter(std::vector<T> xs) {
    std::sort(xs.begin(), xs.end());
    int ctr_same = 1, n = xs.size();
    kd_n2_type ret = 0;
    for (int i=1; i<n; i++) {
      if (xs[i-1] != xs[i]) ctr_same = 0;
      ret += ctr_same++;
    }
    return ret;
  };
}

// ===== Online Algorithms ========================================================== //

namespace FastCorr::MonotonicOnlineCorr {
  template< class T >
  class KendallBase : public Base<T> {
    public:
      virtual double kendall_tau() const = 0;
      double r() const { return kendall_tau(); } // r() is an alias for kendall_tau()
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
    public:
      Kendall() {}
      Kendall(const std::vector<T> &x_vals) {
        for (auto &x : x_vals) push_back(x);
      }
      void push_back(const T &x_val) {
        vals.push_back(std::make_pair(x_val, max_y_ctr++));
        K += ctr_tree.order_of_key(std::make_pair(x_val, -1));
        L += ctr_tree.size() - ctr_tree.order_of_key(std::make_pair(x_val, id_for_tree));
        ctr_tree.insert(std::make_pair(x_val, id_for_tree++));
        n1 += ctr_X[x_val]++;
      }
      void pop_front() {
        auto z = vals.front();
        T x_val = z.first;
        min_y_ctr++;
        // y_i should all decrease by 1, but we can simply ignore that as it won't affect the result
        vals.pop_front();
        int num_lower = ctr_tree.order_of_key(std::make_pair(x_val, -1));
        ctr_tree.erase_kth(num_lower);

        K -= ctr_tree.size() - ctr_tree.order_of_key(std::make_pair(x_val, id_for_tree));
        L -= num_lower;
        // equivalent to n1 -= --ctr_X[x_val]
        int tmp = ctr_X[x_val] - 1;
        n1 -= tmp;
        if (tmp == 0) ctr_X.erase(x_val);
        else ctr_X[x_val] = tmp;
      }
      void push_front(const T &x_val) {
        vals.push_front(std::make_pair(x_val, --min_y_ctr));
        K += ctr_tree.size() - ctr_tree.order_of_key(std::make_pair(x_val, id_for_tree));
        L += ctr_tree.order_of_key(std::make_pair(x_val, -1));
        ctr_tree.insert(std::make_pair(x_val, id_for_tree++));
        n1 += ctr_X[x_val]++;
      }
      void pop_back() {
        auto z = vals.back();
        max_y_ctr--;
        T x_val = z.first;
        vals.pop_back();
        int num_lower = ctr_tree.order_of_key(std::make_pair(x_val, -1));
        ctr_tree.erase_kth(num_lower);

        K -= num_lower;
        L -= ctr_tree.size() - ctr_tree.order_of_key(std::make_pair(x_val, id_for_tree));
        // equivalent to n1 -= --ctr_X[x_val]
        int tmp = ctr_X[x_val] - 1;
        n1 -= tmp;
        if (tmp == 0) ctr_X.erase(x_val);
        else ctr_X[x_val] = tmp;
      }
      double kendall_tau() const {
        int N = vals.size();
        kd_n2_type n0 = (kd_n2_type)N*(N-1)/2;
        // n2 = 0 because y_i has no duplicate values
        if (n1 == n0) return NAN; // denominator will be 0 on tau-b/c
        return (double)(K-L) / sqrt((double)(n0-n1)*(double)n0); // tau-b (n2=0)
      }
      size_t size() const { return vals.size(); }
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
      void push_back(const T &x_val) {
        vals.push_back(std::make_pair(x_val, max_y_ctr++));
      }
      void push_front(const T &x_val) {
        vals.push_front(std::make_pair(x_val, --min_y_ctr));
      }
      void pop_front() {
        // y_i should all decrease by 1, but we can simply ignore that as it won't affect the result
        min_y_ctr++;
        vals.pop_front();
      }
      void pop_back() {
        max_y_ctr--;
        vals.pop_back();
      }
      // O(NlogN) efficient offline algorithm (tau-b)
      double kendall_tau() const {
        return FastCorr::OfflineCorr::kendall_tau(vals);
        //return OfflineCorr::slow_kendall_tau<T>(vals);
      }
      size_t size() const { return vals.size(); }
  };
}
namespace FastCorr::OnlineCorr {
  // TODO
  template< class TX, class TY >
  class Kendall : public Base<TX, TY> {
    public:
      virtual void add(const TX &x_val, const TY &y_val) = 0;
      virtual void remove(const TX &x_val, const TY &y_val) = 0;
      virtual double r() const = 0;
      virtual size_t size() const = 0;
  };
}

// ===== Offline Algorithms ========================================================= //

namespace FastCorr::OfflineCorr {
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
    int n = sorted.size();
    if (n <= 1) return NAN;
    kd_n2_type K = 0, L = 0;
    const kd_n2_type n0 = (kd_n2_type)n*(n-1)/2;

    std::vector<TY> ys(n);
    for (int i=0; i<n; i++) ys[i] = sorted[i].second;
    kd_n2_type n1 = 0, n2 = offline_nC2_counter(ys); // O(nlogn) but faster
    /*
    std::map<T, int> ctr_Y;
    // O(nlogn)
    for (int i=0; i<n; i++) ctr_Y[sorted[i].second]++;
    kd_n2_type n1 = 0, n2 = 0;
    for (auto &p : ctr_Y) n2 += (kd_n2_type)p.second * (p.second-1) / 2;
    */
    // the second int is an unique id to allow for duplicate values in tree
    CountingTree< std::pair<TY, int> > ctr_tree;
    int head = 0;
    for (int i=0; i<n; i++) {
      if (i > 0 && sorted[i-1].first != sorted[i].first) {
        int c = i-head;
        n1 += (kd_n2_type)c*(c-1)/2;
        for (; head<i; head++) {
          ctr_tree.insert(std::make_pair(sorted[head].second, head)); // add y
        }
        // head = i
      }
      // K += #{yj < yi}
      // L += #{yj > yi}
      TY yi = sorted[i].second;
      K += ctr_tree.order_of_key(std::make_pair(yi, -1)); // O(logn)
      L += ctr_tree.size() - ctr_tree.order_of_key(std::make_pair(yi, n)); // O(logn)
    }
    int last_c = n-head;
    n1 += (kd_n2_type)last_c*(last_c-1)/2;

    if (n1 == n0 || n2 == n0) return NAN; // denominator will be 0 on tau-b and tau-c
    // return (double)(K-L) / (double)n0; // tau-a
    //for (auto p : vals) { cout<<"("<<p.first<<", "<<p.second<<"),"; } cout<<" -> tau = "<< (double)(K-L) <<"/"<< sqrt((n0-n1)*(n0-n2)) << "\n";
    return (double)(K-L) / sqrt((double)(n0-n1)*(double)(n0-n2)); // tau-b
    //int m = min(ctr_X.size(), ctr_Y.size());
    //return 2.0*(double)(K-L) / (n*n * (double)(m-1) / (double)m); // tau-c
  }

  // O(N^2) implementation for validation (tau-b)
  template< class TX, class TY >
  double slow_kendall_tau(const std::deque<std::pair<TX, TY> > &vals) {
    int n = vals.size();
    if (n <= 1) return NAN;
    kd_n2_type K = 0, L = 0;
    const kd_n2_type n0 = (kd_n2_type)n*(n-1)/2;
    /*
    std::vector<T> xs(n);
    for (int i=0; i<n; i++) xs[i] = vals[i].first;
    kd_n2_type n1 = offline_nC2_counter(xs);
    for (int i=0; i<n; i++) xs[i] = vals[i].second;
    kd_n2_type n2 = offline_nC2_counter(xs);
    */

    std::map<TX, int> ctr_X;
    std::map<TY, int> ctr_Y;
    // O(nlogn)
    for (int i=0; i<n; i++) ctr_X[vals[i].first]++;
    for (int i=0; i<n; i++) ctr_Y[vals[i].second]++;
    kd_n2_type n1 = 0, n2 = 0;
    for (auto &p : ctr_X) n1 += (kd_n2_type)p.second * (p.second-1) / 2;
    for (auto &p : ctr_Y) n2 += (kd_n2_type)p.second * (p.second-1) / 2;
    // O(N^2)
    for (int i=0; i<n; i++) {
      for (int j=0; j<i; j++) {
        TX xi = vals[i].first,  xj = vals[j].first;
        TY yi = vals[i].second, yj = vals[j].second;
        if ((xi < xj && yi < yj) || (xi > xj && yi > yj)) K++;
        if ((xi < xj && yi > yj) || (xi > xj && yi < yj)) L++;
      }
    }
    if (n1 == n0 || n2 == n0) return NAN; // denominator will be 0 on tau-b and tau-c
    // return (double)(K-L) / (double)n0; // tau-a
    return (double)(K-L) / sqrt((double)(n0-n1)*(double)(n0-n2)); // tau-b
    //int m = min(ctr_X.size(), ctr_Y.size());
    //return 2.0*(double)(K-L) / (n*n * (double)(m-1) / (double)m); // tau-c
  }
}
