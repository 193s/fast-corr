#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../lib/doctest.h"

#include <iostream>
#include <vector>
#include <iomanip>
#include <limits>
#include <chrono>
#include <random>
#include <ctime>
#include <cassert>
using namespace std;

#include "../lib/spearman_algos.hpp"
#include "../lib/kendall_algos.hpp"
using namespace FastCorr;
#define assertmsg(expr, msg) assert(((void)msg, expr))

const int LOOP = 3;
const double EPS = 1e-9;
enum class OPERATION_TYPE { PUSH_FRONT, PUSH_BACK, POP_FRONT, POP_BACK, CALCULATE_R };

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
  for (int i=0; i<T; i++) A[i] = mt()/100000.0;
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
vector<int> generate_random_int_sequence(int T, int seed, bool duplicate) {
  vector<int> A(T);
  mt19937 mt(seed);
  for (int i=0; i<T; i++) A[i] = mt();
  if (duplicate) {
    for (int i=0; i<T/2; i++) A[i] = A[T-1-i];
    shuffle(A.begin(), A.end(), mt);
  }
  else {
    if (contains_duplicates(A)) {
      return generate_random_int_sequence(T, seed+1, duplicate); // re-generate with different seeds
    }
  }
  return A;
}

int system_exec(string cmd_str) {
  char cmd[1000];
  strcpy(cmd, cmd_str.c_str());
  return system(cmd);
}

template< class T >
string internal_stringify(vector<T> x) {
  stringstream result;
  copy(x.begin(), x.end(), ostream_iterator<T>(result, ","));
  return "["+result.str()+"]";
}

TEST_CASE("check results with Python's scipy.stats") {
  //string python_cmd = "python3";
  if (!getenv("PYTHON_CMD")) {
    MESSAGE("env PYTHON_CMD not specified: skipping tests with Python's scipy.stats");
  }
  else {
    string python_cmd = getenv("PYTHON_CMD");
    int seed = time(NULL);
    MESSAGE("PYTHON_CMD=", python_cmd, ", seed=", seed);

    SUBCASE("MonotonicOnlineCorr::Spearman<T> vs scipy.stats.spearmanr (eps=1e-6)") {
      bool duplicate_test = true;
      SUBCASE("with duplicates") { duplicate_test = true; }
      SUBCASE("without duplicates") { duplicate_test = false; }

      for (int loop=0; loop<5; loop++) {
        int n = 15;
        vector<int> xs = generate_random_int_sequence(n, seed+loop, duplicate_test);
        MESSAGE("testing spearman with ", internal_stringify(xs), "...");
        double r = MonotonicOnlineCorr::Spearman<int>(xs).spearman_r();
        CHECK(0 == system_exec(python_cmd
            + " -c 'import scipy.stats; assert(1e-6 >= abs("+to_string(r)+"-"
            + "scipy.stats.spearmanr("+internal_stringify(xs)+",range("+to_string(n)+")).statistic))'"));
      }
    }

    SUBCASE("MonotonicOnlineCorr::Kendall<T> vs scipy.stats.kendalltau (eps=1e-6)") {
      bool duplicate_test = true;
      SUBCASE("with duplicates") { duplicate_test = true; }
      SUBCASE("without duplicates") { duplicate_test = false; }

      for (int loop=0; loop<5; loop++) {
        int n = 15;
        vector<int> xs = generate_random_int_sequence(n, seed+loop, duplicate_test);
        MESSAGE("testing kendall with ", internal_stringify(xs), "...");
        double r = MonotonicOnlineCorr::Kendall<int>(xs).kendall_tau();
        CHECK(0 == system_exec(python_cmd
            + " -c 'import scipy.stats; assert(1e-6 >= abs("+to_string(r)+"-"
            + "scipy.stats.kendalltau("+internal_stringify(xs)+",range("+to_string(n)+")).statistic))'"));
      }
    }

    SUBCASE("convert_array_to_rank vs scipy.stats.rankdata") {
      bool duplicate_test = true;
      SUBCASE("with duplicates") { duplicate_test = true; }
      SUBCASE("without duplicates") { duplicate_test = false; }

      for (int loop=0; loop<3; loop++) {
        int n = 10;
        vector<int> xs = generate_random_int_sequence(n, seed+loop, duplicate_test);
        vector<int> ys = convert_array_to_rank(xs);
        MESSAGE("testing rankdata with ", internal_stringify(xs), "...");
        CHECK(0 == system_exec(python_cmd
            + " -c 'import scipy.stats; assert(all("+internal_stringify(ys)
            + " == 2*scipy.stats.rankdata("+internal_stringify(xs)+")))'"));
      }
    }

  }
}
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
  MonotonicOnlineCorr::SpearmanBase<double> *sp;

  SUBCASE("<=2") {
    for (int repeat=0; repeat<3; repeat++) {
      if (repeat == 0) sp = new MonotonicOnlineCorr::Spearman<double>();
      else if (repeat == 1) sp = new MonotonicOnlineCorr::SpearmanLinear<double>();
      else sp = new MonotonicOnlineCorr::OfflineSpearman<double>();
      REQUIRE(isnan(sp->spearman_r())); // spearman({}) = nan
      sp->push_back(0);
      REQUIRE(isnan(sp->spearman_r())); // spearman({0}) = nan
      sp->push_back(0);
      REQUIRE(isnan(sp->spearman_r())); // spearman({0, 0}) = nan
      sp->push_back(0);
      REQUIRE(isnan(sp->spearman_r())); // spearman({0, 0, 0}) = nan
    }
  }
  SUBCASE("N=125") {
    for (int repeat=0; repeat<3; repeat++) {
      if (repeat == 0) sp = new MonotonicOnlineCorr::Spearman<double>();
      else if (repeat == 1) sp = new MonotonicOnlineCorr::SpearmanLinear<double>();
      else sp = new MonotonicOnlineCorr::OfflineSpearman<double>();

      REQUIRE(isnan(sp->spearman_r())); // spearman({}) = nan
      sp->push_back(0);
      REQUIRE(isnan(sp->spearman_r())); // spearman({0}) = nan
      sp->push_back(1);
      REQUIRE(sp->spearman_r() == 1.0); // spearman({0, 1}) = 1

      // spearman({0, 1, 2, ..., n}) = 1
      for (int i=0; i<123; i++) {
        sp->push_back(2+i);
        if (sp->spearman_r() != 1.0) cout<<"n="<<sp->size()<<", spearman_r="<<sp->spearman_r() <<"\n";
        REQUIRE(sp->spearman_r() == 1.0);
      }
      for (int i=0; i<123; i++) {
        sp->pop_front();
        if (sp->spearman_r() != 1.0) cout<<"n="<<sp->size()<<", spearman_r="<<sp->spearman_r() <<"\n";
        REQUIRE(sp->spearman_r() == 1.0);
      }
    }
  }
}
TEST_CASE("kendall basic testing") {
  MonotonicOnlineCorr::KendallBase<double> *kd;

  SUBCASE("<=2") {
    for (int repeat=0; repeat<3; repeat++) {
      if (repeat == 0) kd = new MonotonicOnlineCorr::Kendall<double>();
      else kd = new MonotonicOnlineCorr::OfflineKendall<double>();
      REQUIRE(isnan(kd->kendall_tau())); // spearman({}) = nan
      kd->push_back(0);
      REQUIRE(isnan(kd->kendall_tau())); // spearman({0}) = nan
      kd->push_back(0);
      REQUIRE(isnan(kd->kendall_tau())); // spearman({0, 0}) = nan
      kd->push_back(0);
      REQUIRE(isnan(kd->kendall_tau())); // spearman({0, 0, 0}) = nan
    }
  }
  SUBCASE("N=125") {
    for (int repeat=0; repeat<2; repeat++) {
      if (repeat == 0) kd = new MonotonicOnlineCorr::Kendall<double>();
      else kd = new MonotonicOnlineCorr::OfflineKendall<double>();

      REQUIRE(isnan(kd->kendall_tau())); // kendall({}) = nan
      kd->push_back(0);
      REQUIRE(isnan(kd->kendall_tau())); // kendall({0}) = nan
      kd->push_back(1);
      REQUIRE(kd->kendall_tau() == 1.0); // kendall({0, 1}) = 1

      // spearman({0, 1, 2, ..., n}) = 1
      for (int i=0; i<123; i++) {
        kd->push_back(2+i);
        if (kd->kendall_tau() != 1.0) cout<<"n="<<kd->size()<<", tau="<<kd->kendall_tau() <<"\n";
        REQUIRE(kd->kendall_tau() == 1.0);
      }
      for (int i=0; i<123; i++) {
        kd->pop_front();
        if (kd->kendall_tau() != 1.0) cout<<"n="<<kd->size()<<", tau="<<kd->kendall_tau() <<"\n";
        REQUIRE(kd->kendall_tau() == 1.0);
      }
    }
  }
}

void internal_test(vector< pair<OPERATION_TYPE, double> > operations, bool VERBOSE); // default: verbose=true
void internal_test(vector< pair<OPERATION_TYPE, double> > operations) { internal_test(operations, true); };

void internal_test_spearman(vector< pair<OPERATION_TYPE, double> > operations, bool verbose) {
  // required for performance measurement
  using std::chrono::high_resolution_clock;
  using std::chrono::duration_cast;
  using std::chrono::duration;
  using std::chrono::milliseconds;

  if (verbose) cout << "========= SPEARMAN =========\n";
  MonotonicOnlineCorr::SpearmanBase<double> *sp;
  //imp_list.push_back(*(new MonotonicOnlineCorr::Spearman<double>()));
  //  new MonotonicOnlineCorr::SpearmanLinear<double>(), new OfflineSpearman<double>(), };

  //vector<long long> ds;
  vector<double> rs;
  for (int repeat=0; repeat<3; repeat++) {
    auto t1 = high_resolution_clock::now();
    if (repeat == 0) {
      // O(logN) efficient algorithm
      if (verbose) cout << "[MonotonicOnlineCorr::Spearman] O(logN)\n";
      sp = new MonotonicOnlineCorr::Spearman<double>();
    }
    else if (repeat == 1) {
      // O(N) insert sort implementation
      if (verbose) cout << "[MonotonicOnlineCorr::SpearmanLinear] O(N)\n";
      sp = new MonotonicOnlineCorr::SpearmanLinear<double>();
    }
    else {
      // O(NlogN) straight forward implementation
      if (verbose) cout << "[MonotonicOnlineCorr::Spearman] O(NlogN)\n";
      sp = new MonotonicOnlineCorr::OfflineSpearman<double>();
    }
    for (int _=0; _<LOOP; _++) {
      int rs_counter = 0;
      for (auto p : operations) {
        double x_val = p.second;
        switch (p.first) {
          case OPERATION_TYPE::PUSH_FRONT: sp->push_front(x_val); break;
          case OPERATION_TYPE::PUSH_BACK:  sp->push_back(x_val);  break;
          case OPERATION_TYPE::POP_FRONT:  sp->pop_front();       break;
          case OPERATION_TYPE::POP_BACK:   sp->pop_back();        break;
          case OPERATION_TYPE::CALCULATE_R:
            double r = sp->spearman_r();
            if (_ == 0 && repeat == 0) rs.push_back(r);
            else if (abs(r - rs[rs_counter]) > EPS) {
              cout << "verify error: (out) r="<<rs[rs_counter]<<" != (correct) " << r << endl;
              REQUIRE(abs(r - rs[rs_counter]) <= EPS);
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
      REQUIRE(sp->size() == 0);
    }
    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> ms_double = (t2 - t1)/LOOP;
    if (verbose) cout << "average execution time: " << ms_double.count() << "ms\n";
  }
}
void internal_test_kendall(vector< pair<OPERATION_TYPE, double> > operations, bool verbose) {
  // required for performance measurement
  using std::chrono::high_resolution_clock;
  using std::chrono::duration_cast;
  using std::chrono::duration;
  using std::chrono::milliseconds;
  if (verbose) cout << "========= KENDALL =========\n";
  MonotonicOnlineCorr::KendallBase<double> *kd;
  vector<double> rs;
  for (int repeat=0; repeat<2; repeat++) {
    auto t1 = high_resolution_clock::now();
    if (repeat == 0) {
      // O(logN) efficient algorithm
      if (verbose) cout << "[MonotonicOnlineCorr::Kendall] O(logN)\n";
      kd = new MonotonicOnlineCorr::Kendall<double>();
    }
    else {
      // O(NlogN) offline implementation
      if (verbose) cout << "[MonotonicOnlineCorr::OfflineKendall] O(NlogN)\n";
      kd = new MonotonicOnlineCorr::OfflineKendall<double>();
    }
    for (int _=0; _<LOOP; _++) {
      int rs_counter = 0;
      for (auto p : operations) {
        double x_val = p.second;
        switch (p.first) {
          case OPERATION_TYPE::PUSH_FRONT: kd->push_front(x_val); break;
          case OPERATION_TYPE::PUSH_BACK:  kd->push_back(x_val);  break;
          case OPERATION_TYPE::POP_FRONT:  kd->pop_front();       break;
          case OPERATION_TYPE::POP_BACK:   kd->pop_back();        break;
          case OPERATION_TYPE::CALCULATE_R:
            double r = kd->kendall_tau();
            if (_ == 0 && repeat == 0) rs.push_back(r);
            else if (abs(r - rs[rs_counter]) > EPS) {
              cout << "verify error: (out) r="<<rs[rs_counter]<<" != (correct) " << r << endl;
              REQUIRE(abs(r - rs[rs_counter]) <= EPS);
            }
            rs_counter++;
            break;
        }
      }
      REQUIRE(kd->size() == 0);
    }
    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> ms_double = (t2 - t1)/LOOP;
    if (verbose) cout << "average execution time: " << ms_double.count() << "ms\n";
  }
}
void internal_test(vector< pair<OPERATION_TYPE, double> > operations, bool verbose) {
  internal_test_spearman(operations, verbose);
  internal_test_kendall(operations, verbose);
}
