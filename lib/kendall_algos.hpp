#pragma once
#include <vector>
#include <deque>
#include <cmath>
#include <algorithm>
#include "fast_corr_base.hpp"
#include "lazy_reversible_rbst.hpp"

namespace FastCorr::MonotonicOnlineCorr {
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
        n1 -= --ctr_X[x_val];
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
        n1 -= --ctr_X[x_val];
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
  class OfflineKendall : public KendallBase<T> {
    int min_y_ctr = 0, max_y_ctr = 0;
    std::deque<std::pair<T, T> > vals;

    public:
      OfflineKendall() {}
      OfflineKendall(const std::vector<T> &x_vals) {
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
        return offline_kendall_tau(vals);
        //return offline_slow_kendall_tau<T>(vals);
      }
      size_t size() const { return vals.size(); }
  };
}
