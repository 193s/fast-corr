#pragma once
#include <cassert>
#include <vector>
#include <cmath>
#include "fast_corr_base.hpp"

namespace FastCorr {
  namespace OfflineCorr {
    /**
     * O(N) straightforward offline algorithm (pearson-r)
     */
    corr_type pearson_r(const std::vector<corr_type> &X, const std::vector<corr_type> &Y) {
      assert(X.size() == Y.size());
      int n = X.size();
      corr_type sum_X = 0, sum_Y = 0, sum_XY = 0;
      corr_type sum_X2 = 0, sum_Y2 = 0;

      for (int i=0; i<n; i++) {
        sum_X += X[i];
        sum_Y += Y[i];
        sum_XY += X[i]*Y[i];
        sum_X2 += X[i]*X[i];
        sum_Y2 += Y[i]*Y[i];
      }
      return ((corr_type)n*sum_XY - sum_X*sum_Y)
              / sqrt(((corr_type)n*sum_X2 - sum_X*sum_X) * ((corr_type)n*sum_Y2 - sum_Y*sum_Y));
    }
  }
  namespace OnlineCorr {
    /**
     * straightforward pearson implementation: O(N) time complexity
     */
    class Pearson : public Base<corr_type, corr_type> {
      public:
        int N = 0;
        corr_type x_mean = 0, y_mean = 0;
        corr_type cov = 0; //  cov  = sum[i=1..n] (x_i-x_mean)*(y_i-y_mean)
        corr_type var_x = 0; // var_x = sum[i=1..n] (x_i - x_mean)^2
        corr_type var_y = 0; // var_y = sum[i=1..n] (y_i - y_mean)^2
        void add(const corr_type &x_val, const corr_type &y_val) override {
          corr_type x_new_mean = x_mean + (x_val - x_mean) / (N+1); // = (x_val + N*x_mean) / (N+1)
          corr_type y_new_mean = y_mean + (y_val - y_mean) / (N+1);
          cov += (x_val - x_mean) * (y_val - y_new_mean);
          var_x += (x_val - x_mean) * (x_val - x_new_mean);
          var_y += (y_val - y_mean) * (y_val - y_new_mean);
          x_mean = x_new_mean, y_mean = y_new_mean;
          N += 1;
        }
        /**
         * note: this does not check if the pair (x_val, y_val) exists in the set
         */
        void remove(const corr_type &x_val, const corr_type &y_val) override {
          corr_type x_new_mean = x_mean + (x_mean - x_val) / (N-1); // = (N*x_mean - x_val) / (N-1)
          corr_type y_new_mean = y_mean + (y_mean - y_val) / (N-1);
          cov -= (x_val - x_mean) * (y_val - y_new_mean);
          var_x -= (x_val - x_mean) * (x_val - x_new_mean);
          var_y -= (y_val - y_mean) * (y_val - y_new_mean);
          x_mean = x_new_mean, y_mean = y_new_mean;
          N -= 1;
        }
        corr_type pearson_r() const noexcept {
          if (N <= 1) return NAN;
          return cov / (sqrt(var_x)*sqrt(var_y));
        }
        virtual corr_type r() const noexcept override { return pearson_r(); }
        virtual size_t size() const noexcept override { return N; }
    };
  }
}
