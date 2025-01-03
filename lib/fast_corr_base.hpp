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
          //if (n <= 1) return NAN;
          sp_d2_type n3 = (sp_d2_type)n*((sp_d2_type)n*n-1);
          //cout<<"n="<<n<<", d="<<d<<"/4, Sx="<<Sx<<"\n";
          if (Sx == n3) return NAN; // rank X_i can not be defined in this case
          else if (Sx == 0) return 1.0 - 1.5*d / n3;
          else {
            // general formula:
            // (1-(6.0/n3)*(D + Sx/12 + Sy/12)) / (sqrt(1 - Sx/n3) * sqrt(1-Sy/n3))
            // when Sy = 0,
            // = (1-(6.0/n3)*(D + Sx/12)) / (sqrt(1 - Sx/n3))
            // = (n3-6.0*(D + Sx/12)) / sqrt(n3*n3 - Sx*n3)
            // = (n3 - Sx/2.0 - 6.0*D) / sqrt(n3*n3 - Sx*n3)
            // = (2*n3 - Sx - 3*(4D)) / (2*sqrt(n3*n3 - Sx*n3))
            return (2.0*n3 - Sx - 3.0*d) / (2.0*n3*sqrt(1.0 - (double)Sx/(double)n3));
          }
        }
        double r() { return spearman_r(); } // r() is an alias for spearman_r()
      protected:
        sp_d2_type Sx = 0; // sum(t_i^3 - t_i)
        // Sy = 0 under the sliding constraint
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
   * CountingTree<T>: set<T> with order_of_key(T val) and erase_kth(int z) implemented
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
      void erase_kth(const int &z) {
        gpp_tree.erase(gpp_tree.find_by_order(z));
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
   */
  template< class T >
  std::vector<int> convert_array_to_rank(const std::vector<T> &X) {
    int n = X.size();
    if (n == 0) return {};
    if (n == 1) return {2};
    // n>=2
    std::vector<std::pair<T, int> > X2(n);
    for (int i=0; i<n; i++) X2[i] = std::pair<T, int>(X[i], i);
    std::sort(X2.begin(), X2.end());
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
    // finalize rank
    int rank = z*2 + 1 + (n-head); // 1 + z + (number of same values - 1)/2
    for (; head<n; head++) ret[X2[head].second] = rank;
    //cout<<"{";for (T x:X)cout<<x<<",";cout<<"} -> ";
    //cout<<"{";for (int x:ret)cout<<x/2.0<<",";cout<<"}\n";
    return ret;
  }

  // O(NlogN) efficient offline algorithm (tau-b)
  // TODO: pair<T, int>, iterator, ...
  template< class T >
  double offline_kendall_tau(const std::deque< std::pair<T, T> > &vals) {
    int N = vals.size();
    if (N <= 1) return NAN;
    kd_n2_type K = 0, L = 0, n0 = (kd_n2_type)N*(N-1)/2;
    std::vector< std::pair<T, T> > sorted;
    for (auto &p : vals) sorted.push_back(p);
    // O(NlogN)
    std::sort(sorted.begin(), sorted.end());

    //std::map<T, int> ctr_X;
    std::map<T, int> ctr_Y;
    // O(NlogN)
    //for (int i=0; i<N; i++) ctr_X[sorted[i].second]++;
    for (int i=0; i<N; i++) ctr_Y[sorted[i].second]++;
    kd_n2_type n1 = 0, n2 = 0;
    //for (auto &p : ctr_X) n1 += (kd_n2_type)p.second * (p.second-1) / 2;
    for (auto &p : ctr_Y) n2 += (kd_n2_type)p.second * (p.second-1) / 2;
    // the second int is an unique id to allow for duplicate values in tree
    CountingTree< std::pair<T, int> > ctr_tree;
    int head = 0;
    for (int i=0; i<N; i++) {
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
      T yi = sorted[i].second;
      K += ctr_tree.order_of_key(std::make_pair(yi, -1)); // O(logN)
      L += ctr_tree.size() - ctr_tree.order_of_key(std::make_pair(yi, N)); // O(logN)
    }
    int last_c = N-head;
    n1 += (kd_n2_type)last_c*(last_c-1)/2;

    if (n1 == n0 || n2 == n0) return NAN; // denominator will be 0 on tau-b and tau-c
    // return (double)(K-L) / (double)n0; // tau-a
    //for (auto p : vals) { cout<<"("<<p.first<<", "<<p.second<<"),"; } cout<<" -> tau = "<< (double)(K-L) <<"/"<< sqrt((n0-n1)*(n0-n2)) << "\n";
    return (double)(K-L) / sqrt((double)(n0-n1)*(double)(n0-n2)); // tau-b
    //int m = min(ctr_X.size(), ctr_Y.size());
    //return 2.0*(double)(K-L) / (N*N * (double)(m-1) / (double)m); // tau-c
  }

  // O(N^2) implementation for validation (tau-b)
  template< class T >
  double offline_slow_kendall_tau(const std::deque<std::pair<T, T> > &vals) {
    int N = vals.size();
    if (N <= 1) return NAN;
    kd_n2_type K = 0, L = 0;
    kd_n2_type n0 = (kd_n2_type)N*(N-1)/2;
    std::map<T, int> ctr_X, ctr_Y;
    // O(NlogN)
    for (int i=0; i<N; i++) ctr_X[vals[i].first]++;
    for (int i=0; i<N; i++) ctr_Y[vals[i].second]++;
    kd_n2_type n1 = 0, n2 = 0;
    for (auto &p : ctr_X) n1 += (kd_n2_type)p.second * (p.second-1) / 2;
    for (auto &p : ctr_Y) n2 += (kd_n2_type)p.second * (p.second-1) / 2;
    // O(N^2)
    for (int i=0; i<N; i++) {
      for (int j=0; j<i; j++) {
        T xi = vals[i].first,  xj = vals[j].first;
        T yi = vals[i].second, yj = vals[j].second;
        if ((xi < xj && yi < yj) || (xi > xj && yi > yj)) K++;
        if ((xi < xj && yi > yj) || (xi > xj && yi < yj)) L++;
      }
    }
    if (n1 == n0 || n2 == n0) return NAN; // denominator will be 0 on tau-b and tau-c
    // return (double)(K-L) / (double)n0; // tau-a
    return (double)(K-L) / sqrt((double)(n0-n1)*(double)(n0-n2)); // tau-b
    //int m = min(ctr_X.size(), ctr_Y.size());
    //return 2.0*(double)(K-L) / (N*N * (double)(m-1) / (double)m); // tau-c
  }
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
