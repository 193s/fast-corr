#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "./doctest.h"
#include "./test_base.hpp"

void internal_sequential_test(int T, int N, bool duplicate_test, bool VERBOSE);
void internal_random_op_test(int Q, bool duplicate_test);
int SEED;
int BENCHMARK_T = 10000, BENCHMARK_N = 1000;
TEST_CASE("benchmark seed") {
  SEED = time(NULL);
  if (getenv("SEED")) {
    SEED = stoi(getenv("SEED"));
    //MESSAGE("env SEED specified: SEED = ", SEED, " (default=time(NULL))");
  }
  if (getenv("BENCHMARK_T")) BENCHMARK_T = stoi(getenv("BENCHMARK_T"));
  if (getenv("BENCHMARK_N")) BENCHMARK_N = stoi(getenv("BENCHMARK_N"));
}
TEST_CASE("random test with duplicates: T=5000 random operations") {
  internal_random_op_test(5000, true);
}
TEST_CASE("random test without duplicates: T=5000 random operations") {
  internal_random_op_test(5000, false);
}
TEST_CASE("random test without duplicates: T=1000, N=100 sequential operations") {
  internal_sequential_test(1000, 100, false, false);
}
TEST_CASE("[benchmark] random test with duplicates: T=10000, N=1000 sequential operations") {
  internal_sequential_test(BENCHMARK_T, BENCHMARK_N, true, true);
}

void internal_random_op_test(int Q, bool duplicate_test) {
  mt19937 mt(SEED);
  vector< pair<OPERATION_TYPE, double> > operations;
  vector<double> pool = generate_random_double_sequence(Q, SEED, duplicate_test);
  int n = 0;
  while ((int)operations.size() < Q) {
    switch (mt()%2) {
      case 0:
        // push_back or push_front
        operations.push_back(
            make_pair(mt()%2?OPERATION_TYPE::PUSH_BACK:OPERATION_TYPE::PUSH_FRONT,
              pool[operations.size()]));
        //operations.push_back(make_pair(PUSH_BACK, pool[operations.size()]));
        operations.push_back(make_pair(OPERATION_TYPE::CALCULATE_R, 0));
        n++;
        break;
      case 1:
        if (n == 0) {
          // invalid operation: try again
          continue;
        }
        else {
          // pop_front or pop_back
          operations.push_back(
              make_pair(mt()%2?OPERATION_TYPE::POP_FRONT:OPERATION_TYPE::POP_BACK, 0));
          //operations.push_back(make_pair(POP_FRONT, 0));
          operations.push_back(make_pair(OPERATION_TYPE::CALCULATE_R, 0));
          n--;
        }
        break;
    }
  }
  while ((int)operations.size() > Q) operations.pop_back();
  while (n-- > 0) operations.push_back(make_pair(OPERATION_TYPE::POP_FRONT, 0));
  cout << "testing with completely random operations: Q="<<Q<<", seed=" << SEED << ", duplicate_test="<<duplicate_test<<"\n";
  internal_test(operations, false);
}

void internal_sequential_test(int T, int N, bool duplicate_test, bool VERBOSE=false) {
  vector<double> A = generate_random_double_sequence(T, SEED, duplicate_test);
  REQUIRE(duplicate_test == contains_duplicates(A));
  if (duplicate_test) cout << "testing with random arrays with duplicate values... seed=" << SEED << "\n";
  else cout << "testing with random arrays without duplicate values... seed=" << SEED << "\n";
  if (VERBOSE) cout << "T=" << T << ", N="<<N<<"\n";
  vector< pair<OPERATION_TYPE, double> > operations;
  // sequencial operation test
  for (int i=0; i<N-1; i++) operations.push_back(make_pair(OPERATION_TYPE::PUSH_BACK, A[i])); // add the first N-1 items
  for (int i=N-1; i<T; i++) {
    // sliding windows: add new & remove old value
    operations.push_back(make_pair(OPERATION_TYPE::PUSH_BACK, A[i]));
    operations.push_back(make_pair(OPERATION_TYPE::CALCULATE_R, 0)); // -> calculate r
    operations.push_back(make_pair(OPERATION_TYPE::POP_FRONT, 0));
  }
  for (int i=0; i<N-1; i++) operations.push_back(make_pair(OPERATION_TYPE::POP_FRONT, 0)); // clean up
  internal_test(operations, VERBOSE);
}
