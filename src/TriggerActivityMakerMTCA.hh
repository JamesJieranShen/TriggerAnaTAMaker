#pragma once
#include "TriggerActivityMaker.hh"

#include <cmath>
#include <queue>
#include <sys/types.h>
#include <vector>

struct TriggerActivityMTCA : public TriggerActivity {
  uint32_t time;

  uint32_t nhits_short;
  float_t sadc_short;
  float_t adc_mean_short;
  float_t adc_std_short;
  float_t charge_balance_short;
  uint32_t cathode_nhits_short;
  float_t cathode_sadc_short;
  float_t cathode_adc_mean_short;
  float_t cathode_adc_std_short;
  float_t cathode_charge_balance_short;
  uint32_t side_nhits_short;
  float_t side_sadc_short;
  float_t side_adc_mean_short;
  float_t side_adc_std_short;
  float_t side_charge_balance_short;

  uint32_t nhits_long;
  float_t sadc_long;
  float_t adc_mean_long;
  float_t adc_std_long;
  float_t charge_balance_long;
  uint32_t cathode_nhits_long;
  float_t cathode_sadc_long;
  float_t cathode_adc_mean_long;
  float_t cathode_adc_std_long;
  float_t cathode_charge_balance_long;
  uint32_t side_nhits_long;
  float_t side_sadc_long;
  float_t side_adc_mean_long;
  float_t side_adc_std_long;
  float_t side_charge_balance_long;

  void LinkTree(TTree& tree, const std::string& prefix = "") override {
    tree.Branch((prefix + "time").c_str(), &time);
    tree.Branch((prefix + "nhits_short").c_str(), &nhits_short);
    tree.Branch((prefix + "sadc_short").c_str(), &sadc_short);
    tree.Branch((prefix + "adc_mean_short").c_str(), &adc_mean_short);
    tree.Branch((prefix + "adc_std_short").c_str(), &adc_std_short);
    tree.Branch((prefix + "charge_balance_short").c_str(),
                &charge_balance_short);
    tree.Branch((prefix + "cathode_nhits_short").c_str(), &cathode_nhits_short);
    tree.Branch((prefix + "cathode_sadc_short").c_str(), &cathode_sadc_short);
    tree.Branch((prefix + "cathode_adc_mean_short").c_str(),
                &cathode_adc_mean_short);
    tree.Branch((prefix + "cathode_adc_std_short").c_str(),
                &cathode_adc_std_short);
    tree.Branch((prefix + "cathode_charge_balance_short").c_str(),
                &cathode_charge_balance_short);
    tree.Branch((prefix + "side_nhits_short").c_str(), &side_nhits_short);
    tree.Branch((prefix + "side_sadc_short").c_str(), &side_sadc_short);
    tree.Branch((prefix + "side_adc_mean_short").c_str(), &side_adc_mean_short);
    tree.Branch((prefix + "side_adc_std_short").c_str(), &side_adc_std_short);
    tree.Branch((prefix + "side_charge_balance_short").c_str(),
                &side_charge_balance_short);

    tree.Branch((prefix + "nhits_long").c_str(), &nhits_long);
    tree.Branch((prefix + "sadc_long").c_str(), &sadc_long);
    tree.Branch((prefix + "adc_mean_long").c_str(), &adc_mean_long);
    tree.Branch((prefix + "adc_std_long").c_str(), &adc_std_long);
    tree.Branch((prefix + "charge_balance_long").c_str(), &charge_balance_long);
    tree.Branch((prefix + "cathode_nhits_long").c_str(), &cathode_nhits_long);
    tree.Branch((prefix + "cathode_sadc_long").c_str(), &cathode_sadc_long);
    tree.Branch((prefix + "cathode_adc_mean_long").c_str(),
                &cathode_adc_mean_long);
    tree.Branch((prefix + "cathode_adc_std_long").c_str(),
                &cathode_adc_std_long);
    tree.Branch((prefix + "cathode_charge_balance_long").c_str(),
                &cathode_charge_balance_long);
    tree.Branch((prefix + "side_nhits_long").c_str(), &side_nhits_long);
    tree.Branch((prefix + "side_sadc_long").c_str(), &side_sadc_long);
    tree.Branch((prefix + "side_adc_mean_long").c_str(), &side_adc_mean_long);
    tree.Branch((prefix + "side_adc_std_long").c_str(), &side_adc_std_long);
    tree.Branch((prefix + "side_charge_balance_long").c_str(),
                &side_charge_balance_long);
  }
};

struct TPCompareMTCA {
  bool operator()(TriggerPrimitive const& a, TriggerPrimitive const& b) const {
    return a.time_start + a.samples_over_threshold >
           b.time_start + b.samples_over_threshold;
  }
};

typedef std::priority_queue<TriggerPrimitive, std::vector<TriggerPrimitive>,
                            TPCompareMTCA>
    TPPriorityQueueMTCA;

struct TPBufferMTCA {
  TPBufferMTCA(uint64_t _window_size_long, uint64_t _window_size_short)
      : m_windowSizeLong(_window_size_long),
        m_windowSizeShort(_window_size_short) {}

  size_t size() const { return m_buffer.size(); }
  uint64_t currentTime() const { return m_currentTime; }

  void clear() {
    TPPriorityQueueMTCA empty;
    m_buffer.swap(empty);
    m_currentTime = 0;
  }
  void add(const TriggerPrimitive& tp);
  void expire(uint64_t current_time);

  TriggerActivityMTCA computeActivity() const;

  uint64_t m_windowSizeLong, m_windowSizeShort;
  uint64_t m_currentTime;
  TPPriorityQueueMTCA m_buffer;
};

class TriggerActivityMakerMTCA
    : public TriggerActivityMaker<TriggerActivityMTCA> {
public:
  TriggerActivityMakerMTCA(uint short_window_size, uint long_window_size)
      : m_tpBuffer(long_window_size, short_window_size) {}

  void operator()(const TriggerPrimitive& input_tp,
                  std::vector<TriggerActivityMTCA>& output_ta) override {
    if (input_tp.time_start > m_tpBuffer.currentTime()) {
      output_ta.push_back(m_tpBuffer.computeActivity());
    }
    m_tpBuffer.add(input_tp);
  }
  void flush(std::vector<TriggerActivityMTCA>& output_ta) override {
    if (m_tpBuffer.size() > 0) {
      output_ta.push_back(m_tpBuffer.computeActivity());
    }
    m_tpBuffer.clear();
  }

private:
  TPBufferMTCA m_tpBuffer;
};
