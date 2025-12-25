#pragma once
#include "TriggerActivityMaker.hh"

#include <cmath>
#include <cstdint>
#include <queue>
#include <sys/types.h>
#include <vector>

struct TriggerActivityCluster : public TriggerActivity {
  uint32_t time_start;
  uint32_t time_end;

  uint32_t n_opdets;
  uint32_t n_tps;
  uint32_t n_tps_max_opdets;
  float_t sadc_sum;
  float_t sadc_mean_opdets;
  float_t sadc_mean_tps;
  uint32_t sadc_max_tps;
  float_t sadc_max_opdets;
  float_t charge_balance_opdets;
  float_t charge_balance_tps;

  float_t tot_sum;
  float_t tot_mean;
  uint32_t tot_max;
  float_t peak_sum;
  float_t peak_mean;
  uint32_t peak_max;

  float_t t_extent;
  float_t dt_max;
  float_t dt_mean;
  float_t dt_std;

  float_t x_extent;
  float_t dx_max;
  float_t dx_mean;
  float_t dx_std;

  float_t y_extent;
  float_t dy_max;
  float_t dy_mean;
  float_t dy_std;
  float_t z_extent;
  float_t dz_max;
  float_t dz_mean;
  float_t dz_std;

  float_t dr_mean;

  float_t wall_fraction_opdets;
  float_t wall_fraction_tps;

  void LinkTree(TTree& tree, const std::string& prefix = "") override {
    tree.Branch((prefix + "time_start").c_str(), &time_start);
    tree.Branch((prefix + "time_end").c_str(), &time_end);
    tree.Branch((prefix + "n_opdets").c_str(), &n_opdets);
    tree.Branch((prefix + "n_tps").c_str(), &n_tps);
    tree.Branch((prefix + "n_tps_max_opdets").c_str(), &n_tps_max_opdets);
    tree.Branch((prefix + "sadc_sum").c_str(), &sadc_sum);
    tree.Branch((prefix + "sadc_mean_opdets").c_str(), &sadc_mean_opdets);
    tree.Branch((prefix + "sadc_mean_tps").c_str(), &sadc_mean_tps);
    tree.Branch((prefix + "sadc_max_tps").c_str(), &sadc_max_tps);
    tree.Branch((prefix + "sadc_max_opdets").c_str(), &sadc_max_opdets);
    tree.Branch((prefix + "charge_balance_opdets").c_str(), &charge_balance_opdets);
    tree.Branch((prefix + "charge_balance_tps").c_str(), &charge_balance_tps);
    tree.Branch((prefix + "tot_sum").c_str(), &tot_sum);
    tree.Branch((prefix + "tot_mean").c_str(), &tot_mean);
    tree.Branch((prefix + "tot_max").c_str(), &tot_max);
    tree.Branch((prefix + "peak_sum").c_str(), &peak_sum);
    tree.Branch((prefix + "peak_mean").c_str(), &peak_mean);
    tree.Branch((prefix + "peak_max").c_str(), &peak_max);
    tree.Branch((prefix + "t_extent").c_str(), &t_extent);
    tree.Branch((prefix + "dt_max").c_str(), &dt_max);
    tree.Branch((prefix + "dt_mean").c_str(), &dt_mean);
    tree.Branch((prefix + "dt_std").c_str(), &dt_std);
    tree.Branch((prefix + "x_extent").c_str(), &x_extent);
    tree.Branch((prefix + "dx_max").c_str(), &dx_max);
    tree.Branch((prefix + "dx_mean").c_str(), &dx_mean);
    tree.Branch((prefix + "dx_std").c_str(), &dx_std);
    tree.Branch((prefix + "y_extent").c_str(), &y_extent);
    tree.Branch((prefix + "dy_max").c_str(), &dy_max);
    tree.Branch((prefix + "dy_mean").c_str(), &dy_mean);
    tree.Branch((prefix + "dy_std").c_str(), &dy_std);
    tree.Branch((prefix + "z_extent").c_str(), &z_extent);
    tree.Branch((prefix + "dz_max").c_str(), &dz_max);
    tree.Branch((prefix + "dz_mean").c_str(), &dz_mean);
    tree.Branch((prefix + "dz_std").c_str(), &dz_std);
    tree.Branch((prefix + "dr_mean").c_str(), &dr_mean);
    tree.Branch((prefix + "wall_fraction_opdets").c_str(), &wall_fraction_opdets);
    tree.Branch((prefix + "wall_fraction_tps").c_str(), &wall_fraction_tps);
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
  TPBufferCluster(uint64_t max_cluster_time)
      : m_maxClusterTime(max_cluster_time), m_currentTime(0) {}

  size_t size() const { return m_buffer.size(); }
  uint64_t currentTime() const { return m_currentTime; }

  void clear() {
    TPPriorityQueueCluster empty;
    m_buffer.swap(empty);
    m_currentTime = 0;
  }
  void add(const TriggerPrimitive& tp);
  void expire(uint64_t current_time);

  void formClusters(std::vector<TriggerActivityCluster>&,
                    int64_t max_hit_time_diff, int64_t min_nhits,
                    float_t max_hit_distance) const;

  TriggerActivityCluster
  makeTriggerActivity(const std::vector<TriggerPrimitive>& cluster_tps) const;

  uint64_t m_maxClusterTime;
  uint64_t m_currentTime;
  TPPriorityQueueCluster m_buffer;
};

class TriggerActivityMakerCluster
    : public TriggerActivityMaker<TriggerActivityCluster> {
public:
  TriggerActivityMakerCluster(int64_t max_cluster_time,
                              int64_t max_hit_time_diff, int64_t min_nhits,
                              float_t max_hit_distance)
      : m_tpBuffer(max_cluster_time), m_maxHitTimeDiff(max_hit_time_diff),
        m_minNhits(min_nhits), m_maxHitDistance(max_hit_distance) {}

  void operator()(const TriggerPrimitive& input_tp,
                  std::vector<TriggerActivityCluster>& output_ta) override {
    if (input_tp.time_start > m_tpBuffer.currentTime()) {
      m_tpBuffer.formClusters(output_ta, m_maxHitTimeDiff, m_minNhits,
                              m_maxHitDistance);
    }
    m_tpBuffer.add(input_tp);
  }
  void flush(std::vector<TriggerActivityCluster>& output_ta) override {
    if (m_tpBuffer.size() > 0) {
      m_tpBuffer.formClusters(output_ta, m_maxHitTimeDiff, m_minNhits,
                              m_maxHitDistance);
    }
    m_tpBuffer.clear();
  }

private:
  TPBufferCluster m_tpBuffer;
  int64_t m_maxHitTimeDiff;
  int64_t m_minNhits;
  float_t m_maxHitDistance;
};
