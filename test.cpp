#include <iostream>
#include <vector>
#include <iomanip>
#include <limits>
#include <chrono>
#include "spearman-algos.hpp"
#include "kendall-algos.hpp"
using namespace std;

const int LOOP = 3;

int main() {
  auto x = OfflineKendall<double>();
  x.push_back(0);
  x.push_back(1);
  x.push_back(1);
  x.push_back(2);
  x.push_back(2);
  x.push_back(1);
  x.push_back(321);
  cout << x.kendall_tau() << "\n";

  int N, T;
  cin >> T >> N;
  vector<double> A(T);
  for (int i=0; i<T; i++) cin >> A[i];
  cout << "T=" << T << ", N="<<N<<": iteration*"<<LOOP<<"\n";

  // input should not have duplicate values
  for (int i=0; i<N; i++) {
    for (int j=i+1; j<N; j++) {
      if (A[i] == A[j]) {
        cout <<"duplicate, error\n";
        exit(1);
      }
    }
  }

  // required for performance measurement
  using std::chrono::high_resolution_clock;
  using std::chrono::duration_cast;
  using std::chrono::duration;
  using std::chrono::milliseconds;

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
          else if (abs(d - ds[i-(N-1)]) > 1e-9) {
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

  if (true) {
    cout << "========= SPEARMAN =========\n";
    OnlineSpearmanBase<double> *sp;
    //imp_list.push_back(*(new OnlineSpearman<double>()));
    //  new OnlineSpearmanLinear<double>(), new OfflineSpearman<double>(), };

    vector<int> ds;
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
          int d = sp->spearman_d();
          if (_ == 0 && repeat == 0) ds.push_back(d);
          else if (d != ds[i-(N-1)]) {
            cout << "verify error: (out) d="<<ds[i-(N-1)]<<" != (correct) "<<d<< endl;
            exit(1);
          }
        }
        while (sp->size()) sp->pop_front(); // clean up
      }
      auto t2 = high_resolution_clock::now();
      duration<double, std::milli> ms_double = (t2 - t1)/LOOP;
      cout << "average execution time: " << ms_double.count() << "ms\n";
    }
    cout << "tests all passed\n";
  }

  //cout << fixed << setprecision(10) << rci_from_d_and_n(d, n) << "\n";
  return 0;
}
