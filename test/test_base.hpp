#include <iostream>
#include <vector>
#include <iomanip>
#include <limits>
#include <chrono>
#include <random>
#include <ctime>
#include <cassert>
#include <iterator>
using namespace std;

#include "../lib/spearman_algos.hpp"
#include "../lib/kendall_algos.hpp"
#include "../lib/pearson_algos.hpp"
using namespace FastCorr;
#define assertmsg(expr, msg) assert(((void)msg, expr))

int LOOP = 3;
std::chrono::duration<double, std::milli> BENCHMARK_MAX_ALLOWED_MS(2*1000); // 2 sec
const int BENCHMARK_MINIMUM_TIMES = 3;

const double EPS = 1e-9;
enum class OPERATION_TYPE { PUSH_FRONT, PUSH_BACK, POP_FRONT, POP_BACK, CALCULATE_R, RESTART_TIMER };

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
    //for (int i=0; i<T; i++) A[i] = mt()%4;
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

bool assert_eq(vector<int> x, vector<int> y) {
  assertmsg(x.size() == y.size(), "assertion failed: size does not match");
  for (int i=0; i<(int)x.size(); i++) {
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

void internal_test(vector< pair<OPERATION_TYPE, double> > operations, bool VERBOSE); // default: verbose=true
void internal_test(vector< pair<OPERATION_TYPE, double> > operations) { internal_test(operations, true); };

void internal_test(vector< pair<OPERATION_TYPE, double> > operations, bool verbose) {
  // required for performance measurement
  using std::chrono::high_resolution_clock;
  using std::chrono::duration_cast;
  using std::chrono::duration;
  using std::chrono::milliseconds;

  for (int corr_type=0; corr_type<2; corr_type++) {
    if (verbose) {
      if (corr_type == 0) cout << "========= SPEARMAN =========\n";
      else if (verbose) cout << "========= KENDALL =========\n";
    }
    MonotonicOnlineCorr::Base<double> *sp;

    //vector<long long> ds;
    vector<double> rs;
    for (int repeat=0; repeat<3; repeat++) {
      int loop_counter = 0;
      auto t1 = high_resolution_clock::now(); // start timer

      if (corr_type == 0) {
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
        else if (repeat == 2) {
          // O(NlogN) straight forward implementation
          if (verbose) cout << "[Offline Spearman] O(NlogN)\n";
          sp = new MonotonicOnlineCorr::OfflineSpearmanForBenchmark<double>();
        }
        else continue;
      }
      else {
        if (repeat == 0) {
          // O(logN) efficient algorithm
          if (verbose) cout << "[MonotonicOnlineCorr::Kendall] O(logN)\n";
          sp = new MonotonicOnlineCorr::Kendall<double>();
        }
        else if (repeat == 1) {
          // O(logN logU) offline implementation (U = MAX_Q ~ len(operations))
          if (verbose) cout << "[OnlineCorr::KendallOnBoundedY] O(logN logU) (U = 1e6)\n";
          sp = new MonotonicOnlineCorr::OnlineNoLimKendallForBenchmark<double>(1e6);
        }
        else if (repeat == 2) {
          // O(NlogN) offline implementation
          if (verbose) cout << "[Offline Kendall] O(NlogN)\n";
          sp = new MonotonicOnlineCorr::OfflineKendallForBenchmark<double>();
        }
        else continue;
      }
      for (int _=0; ; _++) {
        if (!verbose && _ == LOOP) break; // verbose=True <=> benchmark mode
        loop_counter++;
        int rs_counter = 0, num_q = 0;
        for (auto p : operations) {
          double x_val = p.second;
          switch (p.first) {
            case OPERATION_TYPE::PUSH_FRONT: num_q++, sp->push_front(x_val); break;
            case OPERATION_TYPE::PUSH_BACK:  num_q++, sp->push_back(x_val);  break;
            case OPERATION_TYPE::POP_FRONT:  num_q++, sp->pop_front();       break;
            case OPERATION_TYPE::POP_BACK:   num_q++, sp->pop_back();        break;
            case OPERATION_TYPE::RESTART_TIMER:
              num_q = 0;
              t1 = high_resolution_clock::now();
              break;
            case OPERATION_TYPE::CALCULATE_R:
              double r = sp->r();
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
        if (corr_type == 1 && repeat == 1) {
          ((FastCorr::MonotonicOnlineCorr::OnlineNoLimKendallForBenchmark<double> *)sp)->reset();
        }
        auto t2 = high_resolution_clock::now();
        duration<double, std::milli> ms_double = t2 - t1;
        if (verbose && loop_counter >= BENCHMARK_MINIMUM_TIMES && (t2-t1) > BENCHMARK_MAX_ALLOWED_MS) { // benchmark mode: loop end
          break;
        }
      }
      auto t2 = high_resolution_clock::now();
      duration<double, std::milli> ms_double = (t2 - t1)/loop_counter;
      if (verbose) cout << "average execution time: " << std::setprecision(2) << fixed << ms_double.count() << "ms (loop="<<loop_counter<<")\n";
    }
  }
}
// O(N^2) implementation for validation (tau-b)
template< class TX, class TY >
double slow_kendall_tau(const std::deque<std::pair<TX, TY> > &vals) {
  int n = vals.size();
  if (n <= 1) return NAN;
  kd_n2_type K = 0, L = 0;
  const kd_n2_type n0 = (kd_n2_type)n*(n-1)/2;
  /*
  std::vector<T> xs(n);
  for (int i=0; i<n; i++) xs[i] = vals[i].first;
  kd_n2_type n1 = offline_nC2_counter(xs);
  for (int i=0; i<n; i++) xs[i] = vals[i].second;
  kd_n2_type n2 = offline_nC2_counter(xs);
  */

  std::map<TX, int> ctr_X;
  std::map<TY, int> ctr_Y;
  // O(nlogn)
  for (int i=0; i<n; i++) ctr_X[vals[i].first]++;
  for (int i=0; i<n; i++) ctr_Y[vals[i].second]++;
  kd_n2_type n1 = 0, n2 = 0;
  for (auto &p : ctr_X) n1 += (kd_n2_type)p.second * (p.second-1) / 2;
  for (auto &p : ctr_Y) n2 += (kd_n2_type)p.second * (p.second-1) / 2;
  // O(N^2)
  for (int i=0; i<n; i++) {
    for (int j=0; j<i; j++) {
      TX xi = vals[i].first,  xj = vals[j].first;
      TY yi = vals[i].second, yj = vals[j].second;
      if ((xi < xj && yi < yj) || (xi > xj && yi > yj)) K++;
      if ((xi < xj && yi > yj) || (xi > xj && yi < yj)) L++;
    }
  }
  if (n1 == n0 || n2 == n0) return NAN; // denominator will be 0 on tau-b and tau-c
  // return (double)(K-L) / (double)n0; // tau-a
  return (double)(K-L) / sqrt((double)(n0-n1)*(double)(n0-n2)); // tau-b
  //int m = min(ctr_X.size(), ctr_Y.size());
  //return 2.0*(double)(K-L) / (n*n * (double)(m-1) / (double)m); // tau-c
}
