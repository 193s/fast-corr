#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../lib/doctest.h"
#include "../lib/spearman-algos.hpp"
#include "../lib/kendall-algos.hpp"
#include <iostream>
#include <vector>
#include <iomanip>
#include <limits>
#include <chrono>
#include <random>
#include <ctime>
#include <cassert>
using namespace std;

const int LOOP = 3;
const double EPS = 1e-9;

void internal_test(bool duplicate_test, int T, int N);
TEST_CASE("random test without duplicates: T=20000, N=1000") {
  internal_test(false, 20000, 1000);
}
TEST_CASE("random test with duplicates: T=2000, N=1000") {
  internal_test(true, 20000, 1000);
}
void internal_test(bool duplicate_test, int T, int N) {
  cout << "T=" << T << ", N="<<N<<": iteration*"<<LOOP<<"\n";
  vector<double> A(T);
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
}
