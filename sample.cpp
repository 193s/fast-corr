#include <iostream>
#include <vector>
#include <limits>
#include "spearman-algos.hpp"
#include "kendall-algos.hpp"
using namespace std;

int main(int argc, char *argv[]) {
  auto x = OnlineKendall<double>();
  x.push_back(0);
  x.push_back(1);
  x.push_back(1);
  x.push_back(2);
  x.push_back(2);
  x.push_back(1);
  x.push_back(321);
  cout << x.kendall_tau() << "\n";
  auto y = OnlineSpearman<double>();
  y.push_back(0);
  y.push_back(1);
  y.push_back(1);
  y.push_back(2);
  y.push_back(2);
  y.push_back(1);
  y.push_back(321);
  cout << y.spearman_r() << "\n";
  return 0;
}
