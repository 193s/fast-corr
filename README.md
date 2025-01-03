# fast-nonparametric-corr
[![C/C++ CI](https://github.com/193s/fast-nonparametric-corr/actions/workflows/c-cpp.yml/badge.svg)](https://github.com/193s/fast-nonparametric-corr/actions/workflows/c-cpp.yml)  
Currently works only on g++.  

## APIs
| Class | time complexity of each operation |
| ---- | ---- |
| FastCorr::MonotonicOnlineCorr::Spearman\<T\> | O(logN) |
| FastCorr::MonotonicOnlineCorr::SpearmanLinear\<T\> | O(N) |
| FastCorr::MonotonicOnlineCorr::OnlineKendall\<T\> | O(logN) |

### FastCorr::MonotonicOnlineCorr::Spearman\<T\>
- `void push_front(T x_val)`
- `void push_back(T x_val)`
- `void pop_front()`
- `void pop_back()`
- `double spearman_r()` (alias: `double r()`)

OnlineSpearmanLinear\<T\> may work faster on smaller `N`s.  
`T` can be `int`, `double`, ... or any other type with comparison operators defined.

### FastCorr::MonotonicOnlineCorr::Kendall\<T\>
- `void push_front(T x_val)`
- `void push_back(T x_val)`
- `void pop_front()`
- `void pop_back()`
- `double kendall_tau()` (alias: `double r()`)  


## Installation
`git clone https://github.com/193s/fast-nonparametric-corr/`
### testing:
`g++-14 test/basic_tests.cpp -O3`  

### benchmark:
`g++-14 test/benchmark.cpp -O3`  
<!--
- `./a.out r <<< "20000 1000"` : testing on randomized sequence without duplicate values, T=20000, N=1000
- `./a.out d <<< "20000 1000"` : testing on randomized sequence with duplicate values, T=20000, N=1000
-->

### running sample code:
`g++-14 sample.cpp -O3`  
```c++
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
```

result:
```
spearman([0,1,1,2,2,1,321], [1,2,3,4,5,6,7]) = 0.767193
spearman([1,1,2,2,1,321], [1,2,3,4,5,6]) = 0.617213

kendall([0,1,1,2,2,1,321], [1,2,3,4,5,6,7]) = 0.688033
kendall([1,1,2,2,1,321], [1,2,3,4,5,6]) = 0.544949
```

verify this with Python's `scipy.stats`:
```python
from scipy import stats
>>> '%.6f' % stats.spearmanr([0,1,1,2,2,1,321], [1,2,3,4,5,6,7]).statistic
'0.767193'
>>> '%.6f' % stats.spearmanr([1,1,2,2,1,321], [1,2,3,4,5,6]).statistic
'0.617213'

>>> '%.6f' % stats.kendalltau([0,1,1,2,2,1,321], [1,2,3,4,5,6,7]).statistic
'0.688033'
>>> '%.6f' % stats.kendalltau([1,1,2,2,1,321], [1,2,3,4,5,6]).statistic
'0.544949'
```

`PYTHON_CMD=python3 ./basic_tests` will automatically check the library's results with Python's `scipy.stats`. (adjust this to `PYTHON_CMD='pipenv run python' ./basic_tests` etc depending on your Python environment)

### Assumptions
`N<=2642245` for spearman and `N<=4294967296` for kendall is assumed during the calculation (`N` is the current number of pairs `(x,y)` in the structure at the time of each operation).  
For larger `N`s, look for the first few lines of `lib/fast_corr_base.hpp` and modify the data types (`d1_type`, `d2_type`, `kd_n2_type`) to double, \_\_int128, etc.  
f

