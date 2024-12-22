#include <iostream>
#include <vector>
#include <iomanip>
#include <limits>
#include <chrono>
#include <random>
#include <ctime>
#include <cassert>
#include "spearman-algos.hpp"
#include "kendall-algos.hpp"
using namespace std;

const int LOOP = 3;
const double EPS = 1e-9;

bool assert_eq(vector<int> x, vector<int> y) {
  if (x.size() != y.size()) {
    cout << "assertion failed: size does not match\n";
    exit(1);
  }
  for (int i=0; i<x.size(); i++) if (x[i] != y[i]) {
    cout<<"x={";for(int i:x)cout<<i<<",";cout<<"}, y={";for(int i:y)cout<<i<<",";cout<<"}\n";
    cout <<"x["<<i<<"]="<<x[i]<<", y["<<i<<"]="<<y[i]<<"\n";
    cout << "assertion failed\n";
    exit(1);
  }
  return true;
}

// tests on helper functions
void tests_on_helper_functions() {
  // vector<int> convert_array_to_rank(vector<T> arr)
  // [3, 12123, 0] -> [2, 3, 1]*2
  assert_eq(convert_array_to_rank(vector<double>({3, 123123, 0})), vector<int>({4, 6, 2}));
  // [1, 2, 2, 2, 5, 5, 7] -> [1, 3, 3, 3, 5.5, 5.5, 7]*2
  assert_eq(convert_array_to_rank(vector<double>({1, 2, 2, 2, 5, 5, 7})), vector<int>({2, 6, 6, 6, 11, 11, 14}));
  // [1, 1, 1, 1] -> [2.5, 2.5, 2.5, 2.5]
  assert_eq(convert_array_to_rank(vector<double>({1, 1, 1, 1})), vector<int>({5, 5, 5, 5}));

  //
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


int main(int argc, char *argv[]) {
  tests_on_helper_functions();
  // tests
  bool random_input = false;
  bool duplicate_test = false;
  if (argc >= 2) {
    if (argv[1][0] == 'r') {
      random_input = true;
    }
    else if (argv[1][0] == 'd') {
      random_input = true;
      duplicate_test = true;
    }
  }
  int N, T;
  cin >> T >> N;
  cout << "T=" << T << ", N="<<N<<": iteration*"<<LOOP<<"\n";
  vector<double> A(T);
  if (random_input) {
    int seed = time(NULL);
    mt19937 mt(seed);
    for (int i=0; i<T; i++) A[i] = mt();
    if (duplicate_test) {
      for (int i=0; i<T/2; i++) A[i] = A[T-1-i];
      shuffle(A.begin(), A.end(), mt);
      cout << "testing with random arrays with duplicate values... seeed="<<seed<<"\n";
    }
    else {
      cout << "testing with random arrays without duplicate values... seed="<<seed<<"\n";
      for (int i=0; i<T; i++) {
        for (int j=i+1; j<T; j++) {
          if (A[i] == A[j]) {
            cout << "error: duplicate\n";
            exit(1);
          }
        }
      }
    }
  }
  else {
    for (int i=0; i<T; i++) cin >> A[i];

    for (int i=0; i<T; i++) {
      for (int j=i+1; j<T; j++) {
        if (A[i] == A[j]) {
          duplicate_test = true;
          break;
        }
      }
    }
    if (duplicate_test) cout << "input contains duplicate values\n";
  }

  // required for performance measurement
  using std::chrono::high_resolution_clock;
  using std::chrono::duration_cast;
  using std::chrono::duration;
  using std::chrono::milliseconds;


  if (true) {
    cout << "========= SPEARMAN =========\n";
    OnlineSpearmanBase<double> *sp;
    //imp_list.push_back(*(new OnlineSpearman<double>()));
    //  new OnlineSpearmanLinear<double>(), new OfflineSpearman<double>(), };

    //vector<long long> ds;
    vector<double> rs;
    for (int repeat=0; repeat<3; repeat++) {
      auto t1 = high_resolution_clock::now();
      if (repeat == 0) {
        // O(logN) efficient algorithm
        cout << "[OnlineSpearman]\n";
        sp = new OnlineSpearman<double>();
      }
      else if (repeat == 1) {
        // O(N) insert sort implementation
        cout << "[OnlineSpearmanLinear]\n";
        sp = new OnlineSpearmanLinear<double>();
      }
      else {
        // O(NlogN) straight forward implementation
        cout << "[OfflineSpearman]\n";
        sp = new OfflineSpearman<double>();
      }
      for (int _=0; _<LOOP; _++) {
        for (int i=0; i<N-1; i++) sp->push_back(A[i]);
        for (int i=N-1; i<T; i++) {
          // sliding windows: add new & remove old value
          sp->push_back(A[i]);
          sp->pop_front();
          // check results
          double r = sp->spearman_r();
          if (_ == 0 && repeat == 0) rs.push_back(r);
          else if (abs(r - rs[i-(N-1)]) > EPS) {
            cout << "verify error: (out) r="<<rs[i-(N-1)]<<" != (correct) "<<r<< endl;
            exit(1);
          }
          /*
          long long d = sp->spearman_d();
          if (_ == 0 && repeat == 0) ds.push_back(d);
          else if (d != ds[i-(N-1)]) {
            cout << "verify error: (out) d="<<ds[i-(N-1)]<<" != (correct) "<<d<< endl;
            exit(1);
          }
          */
        }
        while (sp->size()) sp->pop_front(); // clean up
      }
      auto t2 = high_resolution_clock::now();
      duration<double, std::milli> ms_double = (t2 - t1)/LOOP;
      cout << "average execution time: " << ms_double.count() << "ms\n";
    }
    cout << "tests all passed\n";
  }
  if (true) {
    cout << "========= KENDALL =========\n";
    OnlineKendallBase<double> *kd;
    vector<double> ds;
    for (int repeat=0; repeat<2; repeat++) {
      auto t1 = high_resolution_clock::now();
      if (repeat == 0) {
        // O(logN) efficient algorithm
        cout << "[OnlineKendall]\n";
        kd = new OnlineKendall<double>();
      }
      else {
        // O(NlogN) offline implementation
        cout << "[OfflineKendall]\n";
        kd = new OfflineKendall<double>();
      }
      for (int _=0; _<LOOP; _++) {
        for (int i=0; i<N-1; i++) kd->push_back(A[i]);
        for (int i=N-1; i<T; i++) {
          // sliding windows: add new & remove old value
          kd->push_back(A[i]);
          kd->pop_front();
          // check results
          double d = kd->kendall_tau();
          if (_ == 0 && repeat == 0) ds.push_back(d);
          else if (abs(d - ds[i-(N-1)]) > EPS) {
            cout << "verify error: (out) d="<<ds[i-(N-1)]<<" != (correct) "<<d<< endl;
            exit(1);
          }
        }
        while (kd->size()) kd->pop_front(); // clean up
      }
      auto t2 = high_resolution_clock::now();
      duration<double, std::milli> ms_double = (t2 - t1)/LOOP;
      cout << "average execution time: " << ms_double.count() << "ms\n";
    }

    ////////////////

    cout << "tests all passed\n";
  }

  //cout << fixed << setprecision(10) << rci_from_d_and_n(d, n) << "\n";
  return 0;
}
