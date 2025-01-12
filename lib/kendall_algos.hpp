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
    /** Base class for Monotonic Kendall */
    template< class T >
    class KendallBase : public Base<T> {
      public:
        virtual corr_type kendall_tau() const noexcept = 0;
        corr_type r() const noexcept override { return kendall_tau(); } // r() is an alias for kendall_tau()
    };

    /** O(logN) Monotonic Online Kendall
     * Time complexity: O(logN) each
     * Space complexity: O(N) overall
     */
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
          if (FAST_CORR_UNLIKELY(n1 == n0)) return corr_NAN; // denominator will be 0 on tau-b/c
          return (corr_type)(K-L) / (
              (corr_type)n0 * sqrt((corr_type)1.0-(corr_type)n1/(corr_type)n0)); // tau-b (n2=0)
        }
        size_t size() const noexcept override { return vals.size(); }
    };
  }

  namespace OnlineCorr {
    template< class TX, class TY >
    class Kendall : public Base<TX, TY> {
      public:
        virtual void add(const TX &x_val, const TY &y_val) = 0;
        virtual void remove(const TX &x_val, const TY &y_val) = 0;
        virtual corr_type r() const noexcept = 0;
        corr_type kendall_tau() const noexcept { return r(); };
        virtual size_t size() const noexcept = 0;
    };
    template<class TX, class TY>
    struct SegmentTreeNode {
      CountingTree<TY> ctr_tree;
      TX width;
      SegmentTreeNode<TX, TY> *l = NULL, *r = NULL;
      // width needs to be 2^n for some integer n
      SegmentTreeNode(const TX &width) : width(width) {}
      // 0 <= pos < width is assumed
      void add(const TX &pos, const TY &add_val) {
        ctr_tree.insert(add_val);
        if (FAST_CORR_UNLIKELY(width == 1)) assert(pos == 0);
        if (FAST_CORR_UNLIKELY(width == 1)) return;
        else {
          if (pos < (width>>1)) {
            if (!l) l = new SegmentTreeNode(width>>1);
            l->add(pos, add_val);
          }
          else {
            if (!r) r = new SegmentTreeNode(width>>1);
            r->add(pos - (width>>1), add_val);
          }
        }
      }
      // 0 <= pos < width is assumed
      void erase_lowerbound(const TX &pos, const TY &add_val) {
        int z = ctr_tree.order_of_key(add_val);
        ctr_tree.erase_kth(z); // erase minimum value greater than add_val
        if (FAST_CORR_UNLIKELY(width == 1)) return;
        else {
          if (pos < (width>>1)) {
            assert(FAST_CORR_LIKELY(l != NULL));
            l->erase_lowerbound(pos, add_val);
            if (l->ctr_tree.size() == 0) delete l, l = NULL;
          }
          else {
            assert(FAST_CORR_LIKELY(r != NULL));
            r->erase_lowerbound(pos - (width>>1), add_val);
            if (r->ctr_tree.size() == 0) delete r, r = NULL;
          }
        }
      }
      // calculate sum(ctr_tree.size()) of [0, pos): 0 <= pos < width is assumed
      int size_sum(const TX &pos) const {
        if (pos == 0) return 0;
        if (FAST_CORR_UNLIKELY(width == pos)) return ctr_tree.size();
        // width > 1
        if (pos < (width>>1)) {
          if (!l) return 0;
          else return l->size_sum(pos);
        }
        else {
          int lsize = l?(int)l->ctr_tree.size():0;
          if (!r) return lsize;
          else return lsize + r->size_sum(pos - (width>>1));
        }
      }
      // count sum(ctr_tree.order_of_key(y_val)) of [0, pos): 0 <= pos < width is assumed
      int count_sum(const TX &pos, const TY &y_val) const {
        if (pos == 0) return 0;
        if (FAST_CORR_UNLIKELY(width == pos)) return ctr_tree.order_of_key(y_val);
        if (pos < (width>>1)) {
          if (!l) return 0;
          else return l->count_sum(pos, y_val);
        }
        else {
          int lsum = l?(int)l->ctr_tree.order_of_key(y_val):0;
          if (!r) return lsum;
          else return lsum + r->count_sum(pos - (width>>1), y_val);
        }
      }
      // get the reference to the node of specified position: 0 <= pos < width is assumed
      SegmentTreeNode<TX, TY>* at(const TX &pos) {
        if (FAST_CORR_UNLIKELY(width == 1)) return this;
        else {
          if (pos < (width>>1)) {
            if (!l) return NULL;
            else return l->at(pos);
          }
          else {
            if (!r) return NULL;
            else return r->at(pos - (width>>1));
          }
        }
      }
    };
    // smallest 2^n such that 2^n >= x
    template<class TX>
    inline TX smallest_power2_larger_than(const TX &x) {
      TX e = 1;
      while (e<x) e<<=1;
      return e;
    }
    /** Online Kendall Algorithm when X_i of added pairs are all positive integers from 0 to MAX_X
     * Time complexity: O(log N * log MAX_X)
     * Space complexity: O(N log N * log MAX_X) overall
     */
    template<class TX, class TY>
    class KendallOnBoundedX : public Kendall<TX, TY> {
      static_assert(std::is_integral<TX>::value, "TX must be an integral type (int, long long, etc)");
      const TX MAX_X;
      int N = 0;
      SegmentTreeNode<TX, std::pair<TY, unsigned int>> *root = NULL;
      unsigned int id_for_tree = 0u; // add unique id to allow for duplicate values in CountingTree
                                    // this is assuming # of add queries < UINT_MAX
      std::map<TX, int> ctr_X;
      std::map<TY, int> ctr_Y;
      kd_n2_type K = 0, L = 0, n1 = 0, n2 = 0;
      inline void _add_value(const TX &x_val, const TY &y_val) {
        n1 += ctr_X[x_val]++;
        n2 += ctr_Y[y_val]++;
      }
      inline void _remove_value(const TX &x_val, const TY &y_val) {
        // equivalent to n1 -= --ctr_X[x_val]
        int tmp = ctr_X[x_val] - 1;
        n1 -= tmp;
        if (tmp == 0) ctr_X.erase(x_val);
        else ctr_X[x_val] = tmp;
        // equivalent to n2 -= --ctr_Y[y_val]
        tmp = ctr_Y[y_val] - 1;
        n2 -= tmp;
        if (tmp == 0) ctr_Y.erase(y_val);
        else ctr_Y[y_val] = tmp;
      }
      // add: SGN=1, remove: SGN=-1
      template<const int SGN=1>
      inline void update_KL(const TX &x_val, const TY &y_val) {
        auto lower_y = std::make_pair(y_val, 0u), upper_y = std::make_pair(y_val, UINT_MAX);
        // time cost: [ ... ] -> O(log^2N) query / ( ...) -> O(logN) query
        //
        // [a3]|[a4]|  K
        // ----+----+----
        //     |    |     -> sum: ctr_Y[y_val]
        // ----+----+----
        // [a1]|(a2)|(..) -> sum: s1
        // left side: x in [0, x_val)
        int a1 = root->count_sum(x_val, lower_y); // #{ (*<x, *<y) }
        int a3 = root->size_sum(x_val) - root->count_sum(x_val, upper_y); // #{ (*<x, *>y) }
        K += SGN*a1;
        L += SGN*a3;
        // right side: x in [x_val+1, MAX_X+1)
        //K += SGN*((root->ctr_tree.size() - root->ctr_tree.order_of_key(upper_y))
        //        - (root->size_sum(x_val+1) - root->count_sum(x_val+1, upper_y))); // #{ (*>x, *>y) }
        auto exact_x = root->at(x_val);
        int a2 = exact_x ? exact_x->ctr_tree.order_of_key(lower_y) : 0; // #{ (x, *<y) }
        int s1 = root->ctr_tree.order_of_key(lower_y);
        int a4 = exact_x ? (exact_x->ctr_tree.size() - exact_x->ctr_tree.order_of_key(upper_y)) : 0;
        auto it = ctr_Y.find(y_val);
        K += SGN*(root->ctr_tree.size() - (s1 + ((it == ctr_Y.end()) ? 0 : it->second) + a3 + a4));
        L += SGN*(s1 - a1 - a2); // #{ (*>x, *<y) }
      }
      public:
        KendallOnBoundedX(TX MAX_X) : MAX_X(MAX_X) {
          root = new SegmentTreeNode<TX, std::pair<TY, unsigned int>>
                                    (smallest_power2_larger_than<TX>(MAX_X+1));
        }
        void add(const TX &x_val, const TY &y_val) override {
          assert(0 <= x_val && x_val <= MAX_X);
          update_KL<+1>(x_val, y_val);
          _add_value(x_val, y_val);
          root->add(x_val, std::make_pair(y_val, ++id_for_tree));
          ++N;
        }
        void remove(const TX &x_val, const TY &y_val) override {
          assert(0 <= x_val && x_val <= MAX_X);
          _remove_value(x_val, y_val);
          // id to be erased doesn't necessarily match for different nodes
          root->erase_lowerbound(x_val, std::make_pair(y_val, 0u));
          update_KL<-1>(x_val, y_val);
          --N;
        }
        corr_type r() const noexcept override {
          kd_n2_type n0 = (kd_n2_type)N*(N-1)/2;
          if (FAST_CORR_UNLIKELY(n1 == n0 || n2 == n0)) return corr_NAN; // denominator will be 0 on tau-b/c
          return (corr_type)(K-L) / (
              (corr_type)n0 * sqrt(
                ((corr_type)1.0-(corr_type)n1/(corr_type)n0) *
                ((corr_type)1.0-(corr_type)n2/(corr_type)n0)
                )); // tau-b (n2=0)
          return 0;
        }
        size_t size() const noexcept override {
          return N;
        }
    };
    /** Online Kendall Algorithm when Y_i of added pairs are all positive integers from 0 to MAX_Y
     * Time complexity: O(log N * log MAX_X)
     * Space complexity: O(N log N * log MAX_X) overall
     * Internally this is equivalent to KendallOnBoundedX, with x and y being exchanged
     */
    template<class TX, class TY>
    class KendallOnBoundedY : public Kendall<TX, TY> {
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
      std::vector<std::pair<TX, TY> > vals(n);
      for (int i=0; i<n; ++i) vals[i] = std::make_pair(x_vals[i], y_vals[i]);
      return kendall_tau(vals);
    }

    // internal function
    template< class TX, class TY >
    corr_type internal_kendall_tau_on_sorted_pairs(const std::vector< std::pair<TX, TY> > &sorted) {
      const int n = sorted.size();
      if (FAST_CORR_UNLIKELY(n <= 1)) return corr_NAN;
      kd_n2_type K = 0, L = 0;
      const kd_n2_type n0 = (kd_n2_type)n*(n-1)/2;

      std::vector<std::pair<TY, int>> ys(n);
      for (int i=0; i<n; ++i) ys[i] = std::make_pair(sorted[i].second, i);
      sort(ys.begin(), ys.end(), // O(nlogn)
          [&](std::pair<TY, int> i, std::pair<TY, int> j) {
            return i.first < j.first; // this is faster: order of the second values does not matter
          });

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

      if (FAST_CORR_UNLIKELY(n1 == n0 || n2 == n0)) return corr_NAN; // denominator will be 0 on tau-b and tau-c
      // return (double)(K-L) / (double)n0; // tau-a
      //for (auto p : vals) { cout<<"("<<p.first<<", "<<p.second<<"),"; } cout<<" -> tau = "<< (double)(K-L) <<"/"<< sqrt((n0-n1)*(n0-n2)) << "\n";
      return (corr_type)(K-L) / ((corr_type)n0 * sqrt(
             ((corr_type)1.0-(corr_type)n1/(corr_type)n0)
            *((corr_type)1.0-(corr_type)n2/(corr_type)n0))); // tau-b
      //int m = min(ctr_X.size(), ctr_Y.size());
      //return 2.0*(double)(K-L) / (n*n * (double)(m-1) / (double)m); // tau-c
    }

  }

  // ===== Wrapper Classes for benchmark ============================================== //
  namespace MonotonicOnlineCorr {
    template< class T >
    class OnlineNoLimKendallForBenchmark : public KendallBase<T> {
      int MAX_Q, min_y_ctr, max_y_ctr;
      OnlineCorr::KendallOnBoundedY<T, int> kd;
      std::deque<std::pair<T, T> > vals;

      // assuming 0 <= min_y_ctr, max_y_ctr <= 2*MAX_Q
      public:
        OnlineNoLimKendallForBenchmark(int MAX_Q) :
          MAX_Q(MAX_Q), min_y_ctr(MAX_Q), max_y_ctr(MAX_Q), kd(2*MAX_Q) {}
        OnlineNoLimKendallForBenchmark(int MAX_Q, const std::vector<T> &x_vals) :
          MAX_Q(MAX_Q), min_y_ctr(MAX_Q), max_y_ctr(MAX_Q), kd(2*MAX_Q) {
          for (auto &x : x_vals) push_back(x);
        }
        void push_back(const T &x_val) override {
          kd.add(x_val, max_y_ctr);
          vals.push_back(std::make_pair(x_val, max_y_ctr++));
        }
        void push_front(const T &x_val) override {
          vals.push_front(std::make_pair(x_val, --min_y_ctr));
          kd.add(x_val, min_y_ctr);
        }
        void pop_front() override {
          // y_i should all decrease by 1, but we can simply ignore that as it won't affect the result
          ++min_y_ctr;
          kd.remove(vals.front().first, vals.front().second);
          vals.pop_front();
        }
        void pop_back() override {
          --max_y_ctr;
          kd.remove(vals.back().first, vals.back().second);
          vals.pop_back();
        }
        void reset() {
          assert(size() == 0);
          min_y_ctr = max_y_ctr = MAX_Q;
        }
        corr_type kendall_tau() const noexcept override { return kd.r(); }
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
}
