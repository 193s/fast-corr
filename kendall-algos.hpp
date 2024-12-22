#pragma once
#include <vector>
#include <deque>
#include <set>
using namespace std;

#include "lazy-reversible-rbst.hpp"
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
#include <ext/pb_ds/tag_and_trait.hpp>
using namespace __gnu_pbds;
// O(NlogN) efficient offline algorithm (tau-b)
// TODO: pair<T, int>, iterator, ...
template< class T >
double offline_kendall_tau(deque<pair<T, T> > &vals) {
  using CountingTree = __gnu_pbds::tree<pair<T, int> ,null_type,less<pair<T, int> >,rb_tree_tag,tree_order_statistics_node_update>;
  int N = vals.size();
  int K = 0, L = 0;
  int n0 = N*(N-1)/2;
  vector< pair<T, T> > sorted;
  for (auto p : vals) sorted.push_back(p);
  // O(NlogN)
  sort(sorted.begin(), sorted.end());


  //map<T, int> ctr_X;
  map<T, int> ctr_Y;
  // O(NlogN)
  //for (int i=0; i<N; i++) ctr_X[sorted[i].second]++;
  for (int i=0; i<N; i++) ctr_Y[sorted[i].second]++;
  int n1 = 0, n2 = 0;
  //for (auto p : ctr_X) n1 += p.second * (p.second-1) / 2;
  for (auto p : ctr_Y) n2 += p.second * (p.second-1) / 2;
  // O(NlogN)
  vector<T> cur_set;
  CountingTree ctr_tree;
  int id = 0; // add unique id to allow for duplicate values in tree
  for (int i=0; i<N; i++) {
    T xi = sorted[i].first, yi = sorted[i].second;
    // K += #{yj < yi}
    // L += #{yj > yi}
    K += ctr_tree.order_of_key(make_pair(yi, -1));
    L += ctr_tree.size() - ctr_tree.order_of_key(make_pair(yi, N)); // assuming id < N
    cur_set.push_back(yi);
    if (i+1 < N && sorted[i+1].first != xi) {
      int c = cur_set.size();
      n1 += c*(c-1)/2;
      while (cur_set.size() > 0) {
        T y = cur_set.back();
        ctr_tree.insert(make_pair(y, id++)); // add y
        cur_set.pop_back();
      }
    }
  }
  if (cur_set.size() > 0) {
    int c = cur_set.size();
    n1 += c*(c-1)/2;
  }
  // return (double)(K-L) / (double)n0; // tau-a
  //for (auto p : vals) { cout<<"("<<p.first<<", "<<p.second<<"),"; } cout<<" -> tau = "<< (double)(K-L) <<"/"<< sqrt((n0-n1)*(n0-n2)) << "\n";
  return (double)(K-L) / sqrt((n0-n1)*(n0-n2)); // tau-b
  //int m = min(ctr_X.size(), ctr_Y.size());
  //return 2.0*(double)(K-L) / (N*N * (double)(m-1) / (double)m); // tau-c
}

// O(N^2) implementation for validation (tau-b)
template< class T >
double offline_slow_kendall_tau(deque<pair<T, T> > vals) {
  int N = vals.size();
  int K = 0, L = 0;
  int n0 = N*(N-1)/2;
  map<T, int> ctr_X, ctr_Y;
  // O(NlogN)
  for (int i=0; i<N; i++) ctr_X[vals[i].first]++;
  for (int i=0; i<N; i++) ctr_Y[vals[i].second]++;
  int n1 = 0, n2 = 0;
  for (auto p : ctr_X) n1 += p.second * (p.second-1) / 2;
  for (auto p : ctr_Y) n2 += p.second * (p.second-1) / 2;
  // O(N^2)
  for (int i=0; i<N; i++) {
    for (int j=0; j<i; j++) {
      T xi = vals[i].first,  xj = vals[j].first;
      T yi = vals[i].second, yj = vals[j].second;
      if ((xi < xj && yi < yj) || (xi > xj && yi > yj)) K++;
      if ((xi < xj && yi > yj) || (xi > xj && yi < yj)) L++;
    }
  }
  // return (double)(K-L) / (double)n0; // tau-a
  return (double)(K-L) / sqrt((n0-n1)*(n0-n2)); // tau-b
  //int m = min(ctr_X.size(), ctr_Y.size());
  //return 2.0*(double)(K-L) / (N*N * (double)(m-1) / (double)m); // tau-c
}


template< class T >
class OnlineKendallBase {
  public:
  virtual void push_back(T x_val);
  virtual void pop_front();
  virtual double kendall_tau();
  virtual int size();
};

template< class T >
class OnlineKendall : public OnlineKendallBase<T> {
  using CountingTree = __gnu_pbds::tree<pair<T, int> ,null_type,less<pair<T, int> >,rb_tree_tag,tree_order_statistics_node_update>;

  CountingTree ctr_tree;
  deque<pair<T, T> > vals;
  map<T, int> ctr_X, ctr_Y;
  int min_y_ctr = 0, max_y_ctr = 0;
  int id_for_tree = 0; // add unique id to allow for duplicate values in CountingTree
  int K = 0, L = 0;
  int n1 = 0;
  public:
  void push_back(T x_val) {
    vals.push_back(make_pair(x_val, max_y_ctr++));
    K += ctr_tree.order_of_key(make_pair(x_val, -1));
    L += ctr_tree.size() - ctr_tree.order_of_key(make_pair(x_val, id_for_tree+111));
    ctr_tree.insert(make_pair(x_val, id_for_tree++));
    n1 += ctr_X[x_val]++;
  }
  void pop_front() {
    auto z = vals.front();
    T x_val = z.first;
    min_y_ctr++;
    // y_i should all decrease by 1, but we can simply ignore that as it won't affect the result
    vals.pop_front();
    ctr_tree.erase(z); // FIXME
    K -= ctr_tree.size() - ctr_tree.order_of_key(make_pair(x_val, id_for_tree+111));
    L -= ctr_tree.order_of_key(make_pair(x_val, -1));
    n1 -= --ctr_X[x_val];
  }
  void push_front(T x_val) {
    vals.push_front(make_pair(x_val, --min_y_ctr));
    K += ctr_tree.size() - ctr_tree.order_of_key(make_pair(x_val, id_for_tree+111));
    L += ctr_tree.order_of_key(make_pair(x_val, -1));
    ctr_tree.insert(make_pair(x_val, id_for_tree++));
    n1 += ctr_X[x_val]++;
  }
  void pop_back() {
    auto z = vals.back();
    max_y_ctr--;
    T x_val = z.first;
    vals.pop_back();
    ctr_tree.erase(z); //FIXME
    K -= ctr_tree.order_of_key(make_pair(x_val, -1));
    L -= ctr_tree.size() - ctr_tree.order_of_key(make_pair(x_val, id_for_tree+111));
    n1 -= --ctr_X[x_val];
  }
  double kendall_tau() {
    int N = vals.size();
    int n0 = N*(N-1)/2;
    const int n2 = 0; // y should have no duplicate values
    return (double)(K-L) / sqrt((n0-n1)*(n0-n2)); // tau-b
  }
  int size() { return vals.size(); }
};

template< class T >
class OfflineKendall : public OnlineKendallBase<T> {
  int min_y_ctr = 0, max_y_ctr = 0;
  deque<pair<T, T> > vals;
  //int N = 0;

  public:
  void push_back(T x_val) {
    vals.push_back(make_pair(x_val, max_y_ctr++));
  }
  void push_front(T x_val) {
    vals.push_front(make_pair(x_val, --min_y_ctr));
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
  double kendall_tau() {
    return offline_kendall_tau(vals);
    //return offline_slow_kendall_tau<T>(vals);
  }
  int size() { return vals.size(); }
};

//using CountingTree = __gnu_pbds::tree<pair<T, int> ,null_type,less<pair<T, int> >,rb_tree_tag,tree_order_statistics_node_update>;

