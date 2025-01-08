#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "./doctest.h"
#include "./test_base.hpp"

TEST_CASE("test settings") {
  if (getenv("LOOP")) {
    LOOP = stoi(getenv("LOOP"));
    MESSAGE("env LOOP specified: LOOP = ", LOOP, " (default=3)");
  }
}

TEST_CASE("OnlineCorr::Pearson basic testing") {
  vector<double> xs, ys;
  OnlineCorr::Pearson ps;
  CHECK(isnan(OfflineCorr::pearson_r(xs, ys))); CHECK(isnan(ps.pearson_r())); // N=0
  ps.add(1, 4), xs.push_back(1), ys.push_back(4);
  CHECK(isnan(OfflineCorr::pearson_r(xs, ys))); CHECK(isnan(ps.pearson_r())); // N=1

  // add(x, y)
  ps.add(5, 3), xs.push_back(5), ys.push_back(3);
  CHECK(abs( OfflineCorr::pearson_r(xs, ys) - ps.pearson_r()) < EPS);
  ps.add(1, 9), xs.push_back(1), ys.push_back(9);
  CHECK(abs( OfflineCorr::pearson_r(xs, ys) - ps.pearson_r()) < EPS);
  ps.add(100, 2), xs.push_back(100), ys.push_back(2);
  CHECK(abs( OfflineCorr::pearson_r(xs, ys) - ps.pearson_r()) < EPS);
  ps.add(1, 9), xs.push_back(1), ys.push_back(9);
  CHECK(abs( OfflineCorr::pearson_r(xs, ys) - ps.pearson_r()) < EPS);
  ps.add(1, 9), xs.push_back(1), ys.push_back(9);
  CHECK(abs( OfflineCorr::pearson_r(xs, ys) - ps.pearson_r()) < EPS);
  ps.add(-10, -10), xs.push_back(-10), ys.push_back(-10);
  CHECK(abs( OfflineCorr::pearson_r(xs, ys) - ps.pearson_r()) < EPS);

  // remove(x, y)
  ps.remove(-10, -10), xs.pop_back(), ys.pop_back();
  CHECK(abs( OfflineCorr::pearson_r(xs, ys) - ps.pearson_r()) < EPS);
  ps.remove(1, 9), xs.pop_back(), ys.pop_back();
  CHECK(abs( OfflineCorr::pearson_r(xs, ys) - ps.pearson_r()) < EPS);
  ps.remove(1, 9), xs.pop_back(), ys.pop_back();
  CHECK(abs( OfflineCorr::pearson_r(xs, ys) - ps.pearson_r()) < EPS);
  ps.remove(100, 2), xs.pop_back(), ys.pop_back();
  CHECK(abs( OfflineCorr::pearson_r(xs, ys) - ps.pearson_r()) < EPS);
  ps.remove(1, 9), xs.pop_back(), ys.pop_back();
  CHECK(abs( OfflineCorr::pearson_r(xs, ys) - ps.pearson_r()) < EPS);

  ps.remove(5, 3), xs.pop_back(), ys.pop_back();
  CHECK(isnan(OfflineCorr::pearson_r(xs, ys))); CHECK(isnan(ps.pearson_r())); // N=1
}

TEST_CASE("OnlineCorr::Pearson with a lot of operations") {
  OnlineCorr::Pearson ps;
  vector<double> rand_xs, rand_ys;
  int Q = 1000000;
  for (int i=0; i<Q; i++) rand_xs.push_back(rand() % 100);
  for (int i=0; i<Q; i++) rand_ys.push_back(rand() % 100);

  SUBCASE("compare with offline algorithm") {
    for (int i=0; i<Q; i++) ps.add(rand_xs[i], rand_ys[i]);
    CHECK(abs( OfflineCorr::pearson_r(rand_xs, rand_ys) - ps.pearson_r()) < EPS);
  }
  SUBCASE("add and remove") {
    ps.add(1, 9);
    ps.add(3, 5);
    ps.add(-10, -10);
    //cout << ps.pearson_r() << "\n";
    double initial_r = ps.pearson_r();

    for (int i=0; i<Q; i++) ps.add(rand_xs[i], rand_ys[i]);
    for (int i=0; i<Q; i++) ps.remove(rand_xs[i], rand_ys[i]);
    //cout << ps.pearson_r() << "\n";
    CHECK(abs( ps.pearson_r() - initial_r ) < 1e-3);
  }
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

    SUBCASE("rankdata vs scipy.stats.rankdata") {
      bool duplicate_test = true;
      SUBCASE("with duplicates") { duplicate_test = true; }
      SUBCASE("without duplicates") { duplicate_test = false; }

      for (int loop=0; loop<3; loop++) {
        int n = 10;
        vector<int> xs = generate_random_int_sequence(n, seed+loop, duplicate_test);
        vector<int> ys = rankdata(xs);
        MESSAGE("testing rankdata with ", internal_stringify(xs), "...");
        CHECK(0 == system_exec(python_cmd
            + " -c 'import scipy.stats; assert(all("+internal_stringify(ys)
            + " == 2*scipy.stats.rankdata("+internal_stringify(xs)+")))'"));
      }
    }

  }
}
// tests on helper functions
TEST_CASE("testing helper functions") {
  // vector<int> rankdata(vector<T> arr)
  // [3, 12123, 0] -> [2, 3, 1]*2
  assert_eq(rankdata(vector<double>({3, 123123, 0})), vector<int>({4, 6, 2}));
  // [1, 2, 2, 2, 5, 5, 7] -> [1, 3, 3, 3, 5.5, 5.5, 7]*2
  assert_eq(rankdata(vector<double>({1, 2, 2, 2, 5, 5, 7})), vector<int>({2, 6, 6, 6, 11, 11, 14}));
  // [1, 1, 1, 1] -> [2.5, 2.5, 2.5, 2.5]
  assert_eq(rankdata(vector<double>({1, 1, 1, 1})), vector<int>({5, 5, 5, 5}));
}

TEST_CASE("spearman basic testing") {
  MonotonicOnlineCorr::SpearmanBase<double> *sp;

  SUBCASE("<=2") {
    for (int repeat=0; repeat<3; repeat++) {
      if (repeat == 0) sp = new MonotonicOnlineCorr::Spearman<double>();
      else if (repeat == 1) sp = new MonotonicOnlineCorr::SpearmanLinear<double>();
      else sp = new MonotonicOnlineCorr::OfflineSpearmanForBenchmark<double>();
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
      else sp = new MonotonicOnlineCorr::OfflineSpearmanForBenchmark<double>();

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
      else kd = new MonotonicOnlineCorr::OfflineKendallForBenchmark<double>();
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
      else kd = new MonotonicOnlineCorr::OfflineKendallForBenchmark<double>();

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
