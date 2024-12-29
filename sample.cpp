#include <iostream>
#include "lib/spearman-algos.hpp"
#include "lib/kendall-algos.hpp"
using namespace std;

int main(int argc, char *argv[]) {
  auto sp = OnlineSpearman<double>({0, 1, 1, 2, 2, 1});
  sp.push_back(321);
  cout << "spearman([0,1,1,2,2,1,321], [1,2,3,4,5,6,7]) = " << sp.spearman_r() << "\n";
  sp.pop_front();
  cout << "spearman([1,1,2,2,1,321], [1,2,3,4,5,6]) = " << sp.spearman_r() << "\n";

  auto kd = OnlineKendall<double>({0, 1, 1, 2, 2, 1});
  kd.push_back(321);
  cout << "kendall([0,1,1,2,2,1,321], [1,2,3,4,5,6,7]) = " << kd.kendall_tau() << "\n";
  kd.pop_front();
  cout << "kendall([1,1,2,2,1,321], [1,2,3,4,5,6]) = " << kd.kendall_tau() << "\n";
  return 0;
}
