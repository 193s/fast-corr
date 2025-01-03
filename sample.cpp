#include <iostream>
#include "lib/spearman_algos.hpp"
#include "lib/kendall_algos.hpp"
using std::cout;

int main(int argc, char *argv[]) {
  auto sp = FastCorr::MonotonicOnlineCorr::Spearman<double>({0, 1, 1, 2, 2, 1});
  sp.push_back(321);
  cout << "spearman([0,1,1,2,2,1,321], [1,2,3,4,5,6,7]) = " << sp.spearman_r() << "\n";
  sp.pop_front();
  cout << "spearman([1,1,2,2,1,321], [1,2,3,4,5,6]) = " << sp.spearman_r() << "\n";

  auto kd = FastCorr::MonotonicOnlineCorr::Kendall<double>({0, 1, 1, 2, 2, 1});
  kd.push_back(321);
  cout << "kendall([0,1,1,2,2,1,321], [1,2,3,4,5,6,7]) = " << kd.kendall_tau() << "\n";
  kd.pop_front();
  cout << "kendall([1,1,2,2,1,321], [1,2,3,4,5,6]) = " << kd.kendall_tau() << "\n";
  return 0;
}
