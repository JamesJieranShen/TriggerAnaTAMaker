#pragma once
#include <cmath>
#include <cstdint>

inline float_t
safe_divide(double a, double b, double default_val = 0.0f) {
  return b == 0 ? default_val : a / b;
}

inline float_t
compute_std(double sum2, double mean, uint32_t n) {
  if (n == 0) return 0.0f;
  int ddof = n == 1 ? 0 : 1;
  float_t variance = (sum2 - mean * mean * n) / (n - ddof);
  return std::sqrt(variance);
}

inline float_t
compute_charge_balance(double sum2, double sum, uint32_t n) {
  float_t mean = n > 0 ? static_cast<float_t>(sum) / n : 0.0f;
  float_t stddev = compute_std(sum2, mean, n);
  return mean > 0.0f ? static_cast<float_t>(stddev) / mean : 0.0f;
}
