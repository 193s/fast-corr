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
#include <set>
using namespace std;

const int LOOP = 3;
const double EPS = 1e-9;

void internal_sequential_test(int T, int N, bool duplicate_test);
void internal_random_op_test(int Q, bool duplicate_test);
void internal_test(vector< pair<int, double> > operations, bool verbose); // default: verbose=true
void internal_test(vector< pair<int, double> > operations) { internal_test(operations, true); };
TEST_CASE("random test with completely random operations: T=5000") {
  internal_random_op_test(5000, true);
}
//TEST_CASE("random test without duplicates: T=20000, N=1000") {
//  internal_sequential_test(20000, 1000, false);
//}
TEST_CASE("random test with duplicates: T=10000, N=1000 sequential operations") {
  internal_sequential_test(10000, 1000, true);
}

// duplicate check in O(NlogN)
template< class T >
bool contains_duplicates(vector<T> &ret) {
  set<T> s;
  for (T x : ret) {
    if (s.find(x) != s.end()) return true;
    s.insert(x);
  }
  return false;
}
vector<double> generate_random_double_sequence(int T, int seed, bool duplicate) {
  vector<double> A(T);
  mt19937 mt(seed);
  for (int i=0; i<T; i++) A[i] = mt();
  if (duplicate) {
    for (int i=0; i<T/2; i++) A[i] = A[T-1-i];
    shuffle(A.begin(), A.end(), mt);
  }
  else {
    if (contains_duplicates(A)) {
      return generate_random_double_sequence(T, seed+1, duplicate); // re-generate with different seeds
    }
  }
  return A;
}

const int PUSH_FRONT  = 0;
const int PUSH_BACK   = 1;
const int POP_FRONT   = 2;
const int POP_BACK    = 3;
const int CALCULATE_R = 4;

void internal_random_op_test(int Q, bool duplicate_test) {
  int seed = time(NULL);
  mt19937 mt(seed);
  vector< pair<int, double> > operations;
  vector<double> pool = generate_random_double_sequence(Q, seed, duplicate_test);
  int n = 0;
  while (operations.size() < Q) {
    switch (mt()%2) {
      case 0:
        operations.push_back(make_pair(PUSH_BACK, pool[operations.size()]));
        operations.push_back(make_pair(CALCULATE_R, 0));
        n++;
        break;
      case 1:
        if (n == 0) {
          // invalid operation: try again
          continue;
        }
        operations.push_back(make_pair(POP_FRONT, 0));
        operations.push_back(make_pair(CALCULATE_R, 0));
        n--;
        break;
    }
  }
  while (operations.size() > Q) operations.pop_back();
  while (n-- > 0) operations.push_back(make_pair(POP_FRONT, 0));
  cout << "testing with completely random operations: Q="<<Q<<", seed=" << seed << ", duplicate_test="<<duplicate_test<<"\n";
  internal_test(operations, false);
}

void internal_sequential_test(int T, int N, bool duplicate_test) {
  cout << "T=" << T << ", N="<<N<<": iteration*"<<LOOP<<"\n";
  int seed = time(NULL);
  vector<double> A = generate_random_double_sequence(T, seed, duplicate_test);
  CHECK(duplicate_test == contains_duplicates(A));
  if (duplicate_test) cout << "testing with random arrays with duplicate values... seed=" << seed << "\n";
  else cout << "testing with random arrays without duplicate values... seed=" << seed << "\n";
  vector< pair<int, double> > operations;
  // sequencial operation test
  for (int i=0; i<N-1; i++) operations.push_back(make_pair(PUSH_BACK, A[i])); // add the first N-1 items
  for (int i=N-1; i<T; i++) {
    // sliding windows: add new & remove old value
    operations.push_back(make_pair(PUSH_BACK, A[i]));
    operations.push_back(make_pair(CALCULATE_R, 0)); // -> calculate r
    operations.push_back(make_pair(POP_FRONT, 0));
  }
  for (int i=0; i<N-1; i++) operations.push_back(make_pair(POP_FRONT, 0)); // clean up
  internal_test(operations);
}

void internal_test(vector< pair<int, double> > operations, bool verbose) {
  // required for performance measurement
  using std::chrono::high_resolution_clock;
  using std::chrono::duration_cast;
  using std::chrono::duration;
  using std::chrono::milliseconds;

  if (true) {
    if (verbose) cout << "========= SPEARMAN =========\n";
    OnlineSpearmanBase<double> *sp;
    //imp_list.push_back(*(new OnlineSpearman<double>()));
    //  new OnlineSpearmanLinear<double>(), new OfflineSpearman<double>(), };

    //vector<long long> ds;
    vector<double> rs;
    for (int repeat=0; repeat<3; repeat++) {
      auto t1 = high_resolution_clock::now();
      if (repeat == 0) {
        // O(logN) efficient algorithm
        if (verbose) cout << "[OnlineSpearman] O(logN)\n";
        sp = new OnlineSpearman<double>();
      }
      else if (repeat == 1) {
        // O(N) insert sort implementation
        if (verbose) cout << "[OnlineSpearmanLinear] O(N)\n";
        sp = new OnlineSpearmanLinear<double>();
      }
      else {
        // O(NlogN) straight forward implementation
        if (verbose) cout << "[OfflineSpearman] O(NlogN)\n";
        sp = new OfflineSpearman<double>();
      }
      for (int _=0; _<LOOP; _++) {
        int rs_counter = 0;
        for (auto p : operations) {
          double x_val = p.second;
          switch (p.first) {
            //case PUSH_FRONT: sp->push_front(x_val); break;
            case PUSH_BACK:  sp->push_back(x_val);  break;
            case POP_FRONT:  sp->pop_front();       break;
            //case POP_BACK:   sp->pop_back();        break;
            case CALCULATE_R:
              double r = sp->spearman_r();
              if (_ == 0 && repeat == 0) rs.push_back(r);
              else if (abs(r - rs[rs_counter]) > EPS) {
                cout << "verify error: (out) r="<<rs[rs_counter]<<" != (correct) " << r << endl;
                CHECK(abs(r - rs[rs_counter]) <= EPS);
              }
              rs_counter++;
              /*
              long long d = sp->spearman_d();
              if (_ == 0 && repeat == 0) ds.push_back(d);
              else if (d != ds[i-(N-1)]) {
                cout << "verify error: (out) d="<<ds[i-(N-1)]<<" != (correct) "<<d<< endl;
                exit(1);
              }
              */
              break;
          }
        }
        CHECK(sp->size() == 0);
      }
      auto t2 = high_resolution_clock::now();
      duration<double, std::milli> ms_double = (t2 - t1)/LOOP;
      if (verbose) cout << "average execution time: " << ms_double.count() << "ms\n";
    }
  }
  if (true) {
    if (verbose) cout << "========= KENDALL =========\n";
    OnlineKendallBase<double> *kd;
    vector<double> rs;
    for (int repeat=0; repeat<2; repeat++) {
      auto t1 = high_resolution_clock::now();
      if (repeat == 0) {
        // O(logN) efficient algorithm
        if (verbose) cout << "[OnlineKendall] O(logN)\n";
        kd = new OnlineKendall<double>();
      }
      else {
        // O(NlogN) offline implementation
        if (verbose) cout << "[OfflineKendall] O(NlogN)\n";
        kd = new OfflineKendall<double>();
      }
      for (int _=0; _<LOOP; _++) {
        int rs_counter = 0;
        for (auto p : operations) {
          double x_val = p.second;
          switch (p.first) {
            case PUSH_FRONT: kd->push_front(x_val); break;
            case PUSH_BACK:  kd->push_back(x_val);  break;
            case POP_FRONT:  kd->pop_front();       break;
            case POP_BACK:   kd->pop_back();        break;
            case CALCULATE_R:
              double r = kd->kendall_tau();
              if (_ == 0 && repeat == 0) rs.push_back(r);
              else if (abs(r - rs[rs_counter]) > EPS) {
                cout << "verify error: (out) r="<<rs[rs_counter]<<" != (correct) " << r << endl;
                CHECK(abs(r - rs[rs_counter]) <= EPS);
              }
              rs_counter++;
              break;
          }
        }
        CHECK(kd->size() == 0);
      }
      auto t2 = high_resolution_clock::now();
      duration<double, std::milli> ms_double = (t2 - t1)/LOOP;
      if (verbose) cout << "average execution time: " << ms_double.count() << "ms\n";
    }
  }
}
