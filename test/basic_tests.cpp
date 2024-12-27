#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../lib/doctest.h"
//#include <ext/pb_ds/assoc_container.hpp>
//#include <ext/pb_ds/tree_policy.hpp>
//#include <ext/pb_ds/tag_and_trait.hpp>
//using namespace __gnu_pbds;

#include <iostream>
#include <vector>
#include <iomanip>
#include <limits>
#include <chrono>
#include <random>
#include <ctime>
#include <cassert>
using namespace std;

#include "../lib/spearman-algos.hpp"
//#include "kendall-algos.hpp"
#define assertmsg(expr, msg) assert(((void)msg, expr))

const int LOOP = 3;
const double EPS = 1e-9;

bool assert_eq(vector<int> x, vector<int> y) {
  assertmsg(x.size() == y.size(), "assertion failed: size does not match");
  for (int i=0; i<x.size(); i++) {
    assertmsg(x[i] == y[i], "assertion failed");
    if (x[i] != y[i]) {
      cout<<"x={";for(int i:x)cout<<i<<",";cout<<"}, y={";for(int i:y)cout<<i<<",";cout<<"}\n";
      cout <<"x["<<i<<"]="<<x[i]<<", y["<<i<<"]="<<y[i]<<"\n";
      cout << "assertion failed\n";
      exit(1);
    }
  }
  return true;
}

// tests on helper functions
TEST_CASE("testing helper functions") {
  // vector<int> convert_array_to_rank(vector<T> arr)
  // [3, 12123, 0] -> [2, 3, 1]*2
  assert_eq(convert_array_to_rank(vector<double>({3, 123123, 0})), vector<int>({4, 6, 2}));
  // [1, 2, 2, 2, 5, 5, 7] -> [1, 3, 3, 3, 5.5, 5.5, 7]*2
  assert_eq(convert_array_to_rank(vector<double>({1, 2, 2, 2, 5, 5, 7})), vector<int>({2, 6, 6, 6, 11, 11, 14}));
  // [1, 1, 1, 1] -> [2.5, 2.5, 2.5, 2.5]
  assert_eq(convert_array_to_rank(vector<double>({1, 1, 1, 1})), vector<int>({5, 5, 5, 5}));
}

TEST_CASE("spearman basic testing") {
  OnlineSpearmanBase<double> *sp;

  for (int repeat=0; repeat<3; repeat++) {
    if (repeat == 0) sp = new OnlineSpearman<double>();
    else if (repeat == 1) sp = new OnlineSpearmanLinear<double>();
    else sp = new OfflineSpearman<double>();

    assert(isnan(sp->spearman_r())); // spearman({}) = nan
    sp->push_back(0);
    assert(isnan(sp->spearman_r())); // spearman({0}) = nan
    sp->push_back(1);
    assert(sp->spearman_r() == 1.0); // spearman({0, 1}) = 1

    // spearman({0, 1, 2, ..., n}) = 1
    for (int i=0; i<123; i++) {
      sp->push_back(2+i);
      if (sp->spearman_r() != 1.0) {
        cout<<"n="<<sp->size()<<", spearman_r="<<sp->spearman_r() <<"\n";
      }
      assert(sp->spearman_r() == 1.0);
    }
    for (int i=0; i<123; i++) {
      sp->pop_front();
      if (sp->spearman_r() != 1.0) {
        cout<<"n="<<sp->size()<<", spearman_r="<<sp->spearman_r() <<"\n";
      }
      assert(sp->spearman_r() == 1.0);
    }
  }
}
