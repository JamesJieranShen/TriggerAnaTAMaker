#pragma once
#include "TriggerActivityMaker.hh"

#include <cmath>
#include <queue>
#include <sys/types.h>
#include <vector>

struct TriggerActivityCluster : public TriggerActivity {
  uint32_t time;

  uint32_t nhits;
  uint32_t ntps;
  int64_t sadc;
  float_t charge_balance;
  int64_t max_dt;

  uint32_t cathode_nhits;
  uint32_t cathode_ntps;
  int64_t cathode_sadc;
  float_t cathode_charge_balance;
  int64_t cathode_max_dt;

  uint32_t side_nhits;
  uint32_t side_ntps;
  int64_t side_sadc;
  float_t side_charge_balance;
  int64_t side_max_dt;

  void LinkTree(TTree& tree, const std::string& prefix = "") override {
    tree.Branch((prefix + "time").c_str(), &time);
    tree.Branch((prefix + "nhits").c_str(), &nhits);
    tree.Branch((prefix + "ntps").c_str(), &ntps);
    tree.Branch((prefix + "sadc").c_str(), &sadc);
    tree.Branch((prefix + "charge_balance").c_str(), &charge_balance);
    tree.Branch((prefix + "max_dt").c_str(), &max_dt);

    tree.Branch((prefix + "cathode_nhits").c_str(), &cathode_nhits);
    tree.Branch((prefix + "cathode_ntps").c_str(), &cathode_ntps);
    tree.Branch((prefix + "cathode_sadc").c_str(), &cathode_sadc);
    tree.Branch((prefix + "cathode_charge_balance").c_str(),
                &cathode_charge_balance);
    tree.Branch((prefix + "cathode_max_dt").c_str(), &cathode_max_dt);

    tree.Branch((prefix + "side_nhits").c_str(), &side_nhits);
    tree.Branch((prefix + "side_ntps").c_str(), &side_ntps);
    tree.Branch((prefix + "side_sadc").c_str(), &side_sadc);
    tree.Branch((prefix + "side_charge_balance").c_str(), &side_charge_balance);
    tree.Branch((prefix + "side_max_dt").c_str(), &side_max_dt);
  }
};

// Pablo's Cluster uses peak time, I think. Let's stick with time_start for now
struct TPCompareCluster {
  bool operator()(TriggerPrimitive const& a, TriggerPrimitive const& b) const {
    return a.time_start > b.time_start;
  }
};

typedef std::priority_queue<TriggerPrimitive, std::vector<TriggerPrimitive>,
                            TPCompareCluster>
    TPPriorityQueueCluster;

struct TPBufferCluster {
  TPBufferCluster(uint64_t max_cluster_duration)
      : m_maxClusterDuration(max_cluster_duration), m_currentTime(0) {}

  size_t size() const { return m_buffer.size(); }
  uint64_t currentTime() const { return m_currentTime; }

  void clear() {
    TPPriorityQueueCluster empty;
    m_buffer.swap(empty);
    m_currentTime = 0;
  }
  void add(const TriggerPrimitive& tp);
  void expire(uint64_t current_time);

  void formCluster(std::vector<TriggerActivityCluster>&) const;

  uint64_t m_maxClusterDuration;
  uint64_t m_currentTime;
  TPPriorityQueueCluster m_buffer;
};

class TriggerActivityMakerCluster
    : public TriggerActivityMaker<TriggerActivityCluster> {
public:
  TriggerActivityMakerCluster(uint64_t max_cluster_duration)
      : m_tpBuffer(max_cluster_duration) {}

  void operator()(const TriggerPrimitive& input_tp,
                  std::vector<TriggerActivityCluster>& output_ta) override {
    if (input_tp.time_start > m_tpBuffer.currentTime()) {
      m_tpBuffer.formCluster(output_ta);
    }
    m_tpBuffer.add(input_tp);
  }
  void flush(std::vector<TriggerActivityCluster>& output_ta) override {
    if (m_tpBuffer.size() > 0) {
      m_tpBuffer.formCluster(output_ta);
    }
    m_tpBuffer.clear();
  }

private:
  TPBufferCluster m_tpBuffer;
};
