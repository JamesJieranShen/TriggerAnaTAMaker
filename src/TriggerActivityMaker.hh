#pragma once
#include "TriggerPrimitive.hh"
#include "TriggerActivity.hh"
#include <vector>

template<typename T>
class TriggerActivityMaker {
public:
  virtual ~TriggerActivityMaker() = default;
  virtual void operator()(const TriggerPrimitive& input_tp,
                          std::vector<T>& output_ta) = 0;
  virtual void flush(std::vector<T>&) {}
  virtual void configure(/* const nlohmann::json& */) {}
  std::vector<T> get_empty_vector() { return std::vector<T>{}; }
  T get_empty_activity() { return T{}; }
};
