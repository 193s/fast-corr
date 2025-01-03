#pragma once
// CountingTree is currently dependent on G++ extensions
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
#include <ext/pb_ds/tag_and_trait.hpp>

namespace FastCorr {
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
}
