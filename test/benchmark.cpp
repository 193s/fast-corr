#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../lib/doctest.h"
#include "./test_base.hpp"

void internal_sequential_test(int T, int N, bool duplicate_test, bool VERBOSE);
void internal_random_op_test(int Q, bool duplicate_test);
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
  internal_sequential_test(10000, 1000, true, true);
}

void internal_random_op_test(int Q, bool duplicate_test) {
  int seed = time(NULL);
  mt19937 mt(seed);
  vector< pair<OPERATION_TYPE, double> > operations;
  vector<double> pool = generate_random_double_sequence(Q, seed, duplicate_test);
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
  cout << "testing with completely random operations: Q="<<Q<<", seed=" << seed << ", duplicate_test="<<duplicate_test<<"\n";
  internal_test(operations, false);
}

void internal_sequential_test(int T, int N, bool duplicate_test, bool VERBOSE) {
  if (VERBOSE) cout << "T=" << T << ", N="<<N<<": iteration*"<<LOOP<<"\n";
  int seed = time(NULL);
  vector<double> A = generate_random_double_sequence(T, seed, duplicate_test);
  REQUIRE(duplicate_test == contains_duplicates(A));
  if (duplicate_test) cout << "testing with random arrays with duplicate values... seed=" << seed << "\n";
  else cout << "testing with random arrays without duplicate values... seed=" << seed << "\n";
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
