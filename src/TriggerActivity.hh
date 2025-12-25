#pragma once
#include <TTree.h>
#include <string>

struct TriggerActivity {
public:
  virtual void LinkTree(TTree& tree, const std::string& prefix = "") = 0;
  virtual ~TriggerActivity() = default;
};
