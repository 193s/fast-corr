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

  namespace MonotonicOnlineCorr {
    template< class T >
    class KendallBase {
      public:
        virtual void push_front(T x_val) = 0;
        virtual void push_back(T x_val) = 0;
        virtual void pop_front() = 0;
        virtual void pop_back() = 0;
        virtual double kendall_tau() = 0;
        virtual int size() = 0;
    };

    template< class T >
    class SpearmanBase {
      public:
        virtual void push_back(T x_val) = 0;
        virtual void push_front(T x_val) = 0;
        virtual void pop_front() = 0;
        virtual void pop_back() = 0;
        virtual sp_d2_type spearman_d() = 0;
        virtual int size() = 0;
        double spearman_r() {
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
      protected:
        sp_d2_type Sx = 0; // sum(t_i^3 - t_i)
        // Sy = 0 under the sliding constraint
    };
  }

  namespace OnlineCorr {
    template< class TX, class TY >
    class Base {
      public:
        virtual void add(TX x_val, TY y_val) = 0;
        virtual void remove(TX x_val, TY y_val) = 0;
        virtual double r() = 0;
        virtual int size() = 0;
    };

    template< class TX, class TY >
    class Pearson : public Base<TX, TY> {
      public:
        virtual void add(TX x_val, TY y_val) = 0;
        virtual void remove(TX x_val, TY y_val) = 0;
        virtual double r() = 0;
        virtual int size() = 0;
    };

    template< class TX, class TY >
    class Kendall : public Base<TX, TY> {
      public:
        virtual void add(TX x_val, TY y_val) = 0;
        virtual void remove(TX x_val, TY y_val) = 0;
        virtual double r() = 0;
        virtual int size() = 0;
    };
  }

  // CountingTree<T>: set<T> with order_of_key(T val) and erase_kth(int z) implemented
  template< class T >
  class CountingTree {
    using GPPTree = __gnu_pbds::tree< T, __gnu_pbds::null_type, std::less<T>,
    __gnu_pbds::rb_tree_tag, __gnu_pbds::tree_order_statistics_node_update>;
    GPPTree gpp_tree;

    public:
      void insert(T val) {
        gpp_tree.insert(val);
      }
      int order_of_key(T val) {
        return gpp_tree.order_of_key(val);
      }
      void erase_kth(int z) {
        gpp_tree.erase(gpp_tree.find_by_order(z));
      }
      int size() {
        return gpp_tree.size();
      }
  };

  // multiply all ranks by 2 to treat them as integer
  // e.g. [3, 12123, 0] -> [2, 3, 1]*2
  // e.g. [1, 2, 2, 2, 5, 5, 7] -> [1, 3, 3, 3, 5.5, 5.5, 7]*2
  template< class T >
  std::vector<int> convert_array_to_rank(std::vector<T> X) {
    int n = X.size();
    if (n == 0) return {};
    if (n == 1) return {2};
    // n>=2
    std::vector<std::pair<T, int> > X2(n);
    for (int i=0; i<n; i++) X2[i] = std::pair<T, int>(X[i], i);
    std::sort(X2.begin(), X2.end());
    std::vector<int> ret(n);
    //for (int i=0; i<n; i++) ret[X2[i].second] = 2*(i+1); // works only on unique arrays
    int z = 0;
    T last_seen;
    std::stack<int> st;
    for (int i=0; i<n; i++) {
      int pos = X2[i].second;
      if (i > 0 && last_seen != X2[i].first) {
        // finalize rank
        int rank = z*2 + 1 + st.size();// 1 + z + (st.size()-1)/2
        z += st.size();
        while (!st.empty()) {
          ret[st.top()] = rank;
          st.pop();
        }
      }
      st.push(pos);
      last_seen = X2[i].first;
    }
    if (!st.empty()) {
      // finalize rank
      int rank = z*2 + 1 + st.size();// 1 + z + (st.size()-1)/2
      z += st.size();
      while (!st.empty()) {
        ret[st.top()] = rank;
        st.pop();
      }
    }
    //cout<<"{";for (T x:X)cout<<x<<",";cout<<"} -> ";
    //cout<<"{";for (int x:ret)cout<<x/2.0<<",";cout<<"}\n";
    return ret;
  }

  // O(NlogN) efficient offline algorithm (tau-b)
  // TODO: pair<T, int>, iterator, ...
  template< class T >
  double offline_kendall_tau(std::deque< std::pair<T, T> > &vals) {
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
    // O(NlogN)
    std::vector<T> cur_set;
    CountingTree< std::pair<T, int> > ctr_tree;
    int id = 0; // add unique id to allow for duplicate values in tree
    for (int i=0; i<N; i++) {
      T xi = sorted[i].first, yi = sorted[i].second;
      // K += #{yj < yi}
      // L += #{yj > yi}
      K += ctr_tree.order_of_key(std::make_pair(yi, -1));
      L += ctr_tree.size() - ctr_tree.order_of_key(std::make_pair(yi, N)); // assuming id < N
      cur_set.push_back(yi);
      if (i+1 < N && sorted[i+1].first != xi) {
        int c = cur_set.size();
        n1 += (kd_n2_type)c*(c-1)/2;
        while (cur_set.size() > 0) {
          T y = cur_set.back();
          ctr_tree.insert(std::make_pair(y, id++)); // add y
          cur_set.pop_back();
        }
      }
    }
    if (cur_set.size() > 0) {
      int c = cur_set.size();
      n1 += (kd_n2_type)c*(c-1)/2;
    }
    if (n1 == n0 || n2 == n0) return NAN; // denominator will be 0 on tau-b and tau-c
    // return (double)(K-L) / (double)n0; // tau-a
    //for (auto p : vals) { cout<<"("<<p.first<<", "<<p.second<<"),"; } cout<<" -> tau = "<< (double)(K-L) <<"/"<< sqrt((n0-n1)*(n0-n2)) << "\n";
    return (double)(K-L) / sqrt((double)(n0-n1)*(double)(n0-n2)); // tau-b
    //int m = min(ctr_X.size(), ctr_Y.size());
    //return 2.0*(double)(K-L) / (N*N * (double)(m-1) / (double)m); // tau-c
  }

  // O(N^2) implementation for validation (tau-b)
  template< class T >
  double offline_slow_kendall_tau(std::deque<std::pair<T, T> > vals) {
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
}
