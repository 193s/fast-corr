# fast-corr
[![C/C++ CI](https://github.com/193s/fast-corr/actions/workflows/c-cpp.yml/badge.svg)](https://github.com/193s/fast-corr/actions/workflows/c-cpp.yml)
[![GitHub license](https://img.shields.io/badge/license-MIT-blue.svg)](https://raw.githubusercontent.com/martinus/unordered_dense/main/LICENSE)  

A header-only C++11 library providing fast implementations of various correlation coefficients.

## Installation
`git clone https://github.com/193s/fast-corr/`

## API
### MonotonicOnlineCorr - Online correlation algorithms under partial monotonic constraints
| Class name | Time complexity of each operation |
| ---- | ---- |
| FastCorr::MonotonicOnlineCorr::Spearman\<T\>       | O(logN) |
| FastCorr::MonotonicOnlineCorr::SpearmanLinear\<T\> | O(N)    |
| FastCorr::MonotonicOnlineCorr::Kendall\<T\>        | O(logN) |

Details and applications of the constraints are discussed in my paper yet to be published.
#### FastCorr::MonotonicOnlineCorr::Spearman\<T\>
- `void push_front(const T &x_val)`
- `void push_back(const T &x_val)`
- `void pop_front()`
- `void pop_back()`
- `double spearman_r()` (alias: `double r()`)

OnlineSpearmanLinear\<T\> works a bit faster on smaller `N`s and uses fewer memory space.  
`T` can be `int`, `double`, ... or any other type with comparison operators defined.  
Currently works a bit faster on g++.  

#### FastCorr::MonotonicOnlineCorr::Kendall\<T\>
- `void push_front(const T &x_val)`
- `void push_back(const T &x_val)`
- `void pop_front()`
- `void pop_back()`
- `double kendall_tau()` (alias: `double r()`)  

`kendall_tau()` function returns tau-b (tau-c option will be added soon).

----------------

### OnlineCorr - Online correlation algorithms (no constraints)
| Class name | Time complexity of each operation | Overall memory use |
| ---- | ---- | ---- |
| FastCorr::OnlineCorr::Pearson | O(1) | O(1) |
<!--| FastCorr::OnlineCorr::Spearman\<T\> | O(logN) |
| FastCorr::OnlineCorr::Kendall\<T\> | O(N) |-->

#### FastCorr::OnlineCorr::Pearson
Online algorithm works on O(1) time complexity on each query and uses O(1) memory space.  
The result may contain errors incurred by arhithmetic operations of double datatype.  
  - `void add(double x_val, double y_val)`
  - `void remove(double x_val, double y_val)`
  - `double pearson_r()` (alias: `double r()`)

----------------

### OfflineCorr - Offline implementations
| Function name | Time complexity |
| ---- | ---- |
| FastCorr::OfflineCorr::spearman\_r\<TX, TY\>  | O(NlogN) |
| FastCorr::OfflineCorr::kendall\_tau\<TX, TY\> | O(NlogN) |
| FastCorr::OfflineCorr::pearson\_r        | O(N)     |

This is nothing new, just a very straightforward set of implementations. The Kendall implementation is comparatively more efficient.  
- `double FastCorr::OfflineCorr::spearman_r<TX, TY>(const vector<TX> &x_vals, const vector<TY> &y_vals)`
- `double FastCorr::OfflineCorr::kendall_tau<TX, TY>(const vector<TX> &x_vals, const vector<TY> &y_vals)`
- `double FastCorr::OfflineCorr::pearson_r(const vector<double> &x_vals, const vector<double> &y_vals)`

----------------

## Usage
(Makefile takes C\_COMPILER environment variable: e.g. `C_COMPILER='clang++ --std=c++14' make`)

### testing:
`make basic_tests && ./basic_tests`  

### benchmark:
`make benchmark && ./benchmark`  

```
T=10000, N=1000
========= SPEARMAN =========
[MonotonicOnlineCorr::Spearman] O(logN)
average execution time: 14.62ms (loop=137)
[MonotonicOnlineCorr::SpearmanLinear] O(N)
average execution time: 13.58ms (loop=148)
[Offline Spearman] O(NlogN)
average execution time: 667.11ms (loop=4)
========= KENDALL =========
[MonotonicOnlineCorr::Kendall] O(logN)
average execution time: 5.46ms (loop=367)
[Offline Kendall] O(NlogN)
average execution time: 708.41ms (loop=3)
```
<!--
- `./a.out r <<< "20000 1000"` : testing on randomized sequence without duplicate values, T=20000, N=1000
- `./a.out d <<< "20000 1000"` : testing on randomized sequence with duplicate values, T=20000, N=1000
-->

### running a sample code:
`g++-14 sample.cpp -O3 -o sample && ./sample`  
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

### Assumptions in MonotonicOnlineCorr
`N<=2642245` for spearman and `N<=4294967296` for kendall are assumed during the calculation (`N` is the current number of pairs `(x,y)` in the structure at the time of each operation) to avoid integer overflow.  
For larger `N`s, look for the first few lines of `lib/fast_corr_base.hpp` and modify the data types (`d1_type`, `d2_type`, `kd_n2_type`) to double, \_\_int128, etc.  

