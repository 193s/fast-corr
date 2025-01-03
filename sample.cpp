#include <iostream>
#include <cassert>
#include <cmath>
#include "lib/spearman_algos.hpp"
#include "lib/kendall_algos.hpp"
using std::cout;

int main(int argc, char *argv[]) {
  // online algorithm examples
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

  // offline algorithm examples
  double sp2 = FastCorr::OfflineCorr::spearman_r<double>(
      {1, 1, 2, 2, 1, 321},
      {1, 2, 3, 4, 5, 6});
  cout<<"sp2="<<sp2<<"\n";
  assert(abs(sp2 - sp.spearman_r() < 1e-9));
  double kd2 = FastCorr::OfflineCorr::kendall_tau<double>(
      {1, 1, 2, 2, 1, 321},
      {1, 2, 3, 4, 5, 6});
  cout<<"kd2="<<kd2<<"\n";
  assert(abs(kd2 - kd.kendall_tau() < 1e-9));
  return 0;
}
