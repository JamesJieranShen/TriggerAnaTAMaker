#include "TriggerActivityMakerCluster.hh"

#include "PMTInfo.hh"
#include "utils.hh"
#include <TVector3.h>

// assumes tps are sorted by time_start.
std::vector<std::vector<TriggerPrimitive>>
split_on_timediff(const std::vector<TriggerPrimitive>& tps,
                  int64_t max_hit_time_diff) {
  std::vector<std::vector<TriggerPrimitive>> clusters;
  if (tps.empty()) {
    return clusters;
  }
  std::vector<TriggerPrimitive> current_cluster;
  current_cluster.push_back(tps.at(0));
  for (size_t i = 1; i < tps.size(); ++i) {
    if (static_cast<int64_t>(tps.at(i).time_start - tps.at(i - 1).time_start) <=
        max_hit_time_diff) {
      current_cluster.push_back(tps.at(i));
    } else {
      clusters.push_back(current_cluster);
      current_cluster.clear();
      current_cluster.push_back(tps.at(i));
    }
  }
  if (!current_cluster.empty()) {
    clusters.push_back(current_cluster);
  }
  return clusters;
}

void
filter_by_distance(std::vector<TriggerPrimitive>& cluster_tps,
                   float_t max_hit_distance) {
  // removes tps from the cluster_tps list if there are no neighbors within
  // max_hit_distance.
  const PMTInfo& pmt_info = PMTInfo::Instance();
  int min_neighbors = 1; // NOTE: if needed make this configurable.
  std::vector<TriggerPrimitive> filtered_tps;
  for (const TriggerPrimitive& curr_tp : cluster_tps) {
    int curr_opdet_id = channel_to_opdet(curr_tp.channel);
    TVector3 curr_pos = pmt_info.GetPositionTVec(curr_opdet_id);
    int neighbor_count = 0;

    for (const TriggerPrimitive& other_tp : cluster_tps) {
      if (curr_tp == other_tp) continue;
      int other_opdet_id = channel_to_opdet(other_tp.channel);
      float_t distance2 =
          (curr_pos - pmt_info.GetPositionTVec(other_opdet_id)).Mag2();
      if (distance2 <= max_hit_distance * max_hit_distance) {
        neighbor_count++;
        if (neighbor_count >= min_neighbors) {
          break;
        }
      }
    } // end for other_tp
    if (neighbor_count >= min_neighbors) {
      filtered_tps.push_back(curr_tp);
    }
  } // end for curr_tp
  std::swap(cluster_tps, filtered_tps);
}

template <typename T>
void
compute_metrics(std::vector<T> vals, float_t& extent_val, float_t& max_val,
                float_t& mean_val, float_t& std_val) {
  if (vals.size() < 2) {
    max_val = 0.0f;
    mean_val = 0.0f;
    return;
  }
  std::sort(vals.begin(), vals.end());
  extent_val = static_cast<float_t>(vals.back()) - static_cast<float_t>(vals.front());
  float_t sum_deltas = 0;
  max_val = 0;
  for (size_t i = 1; i < vals.size(); ++i) {
    float_t delta = static_cast<float_t>(vals[i] - vals[i - 1]);
    sum_deltas += delta;
    max_val = std::max(max_val, delta);
  }
  mean_val = compute_mean(sum_deltas, vals.size() - 1);
  std_val = compute_std(static_cast<double>(sum_deltas * sum_deltas),
                        static_cast<double>(mean_val), vals.size() - 1);
}

float_t
compute_mean_distance(const std::vector<TriggerPrimitive>& tps) {
  // assumes tps are sorted by time_start.
  if (tps.size() < 2) {
    return 0.0f;
  }
  const PMTInfo& pmt_info = PMTInfo::Instance();
  size_t n_pairs = 0;
  double total_distance = 0.0;
  for (size_t i = 1; i < tps.size(); ++i) {
    const TriggerPrimitive& prev_tp = tps.at(i - 1);
    const TriggerPrimitive& curr_tp = tps.at(i);
    float_t distance =
        (pmt_info.GetPositionTVec(channel_to_opdet(curr_tp.channel)) -
         pmt_info.GetPositionTVec(channel_to_opdet(prev_tp.channel)))
            .Mag();
    total_distance += distance;
    n_pairs += 1;
  }
  return compute_mean(total_distance, n_pairs);
}

void
TPBufferCluster::expire(uint64_t current_time) {
  int64_t expire_time = current_time - m_maxClusterTime;
  while (!m_buffer.empty() &&
         static_cast<int64_t>(m_buffer.top().time_start) < expire_time) {
    m_buffer.pop();
  }
}

void
TPBufferCluster::add(const TriggerPrimitive& tp) {
  m_currentTime = tp.time_start;
  expire(tp.time_start);
  m_buffer.push(tp);
}
void
TPBufferCluster::formClusters(std::vector<TriggerActivityCluster>& output,
                              int64_t max_hit_time_diff, int64_t min_nhits,
                              float_t max_hit_distance) const {
  if (m_buffer.empty()) return;
  std::vector<TriggerPrimitive> tps;
  TPPriorityQueueCluster temp_queue = m_buffer;
  while (!temp_queue.empty()) {
    tps.push_back(temp_queue.top());
    temp_queue.pop();
  }
  std::sort(tps.begin(), tps.end(),
            [](const TriggerPrimitive& a, const TriggerPrimitive& b) {
              return a.time_start < b.time_start;
            });
  std::vector<std::vector<TriggerPrimitive>> time_clusters =
      split_on_timediff(tps, max_hit_time_diff);
  for (std::vector<TriggerPrimitive>& cluster_tps : time_clusters) {
    if (cluster_tps.size() < static_cast<size_t>(min_nhits)) {
      continue;
    }
    filter_by_distance(cluster_tps, max_hit_distance);
    if (cluster_tps.size() < static_cast<size_t>(min_nhits)) {
      continue;
    }
    std::sort(cluster_tps.begin(), cluster_tps.end());
    output.push_back(makeTriggerActivity(cluster_tps));
  }
}

TriggerActivityCluster
TPBufferCluster::makeTriggerActivity(
    const std::vector<TriggerPrimitive>& cluster_tps) const {
  TriggerActivityCluster activity;
  activity.n_tps = cluster_tps.size();
  activity.dr_mean = compute_mean_distance(cluster_tps);
  std::map<uint32_t, int64_t> opdet_sadcs;
  std::map<uint32_t, uint32_t> opdet_ntps;
  PMTInfo& pmt_info = PMTInfo::Instance();

  double total_sadc2_tps = 0.0;
  std::vector<uint64_t> times;
  times.reserve(cluster_tps.size());
  std::vector<float_t> xs;
  xs.reserve(cluster_tps.size());
  std::vector<float_t> ys;
  ys.reserve(cluster_tps.size());
  std::vector<float_t> zs;
  zs.reserve(cluster_tps.size());
  uint32_t wall_tp_count = 0;
  for (const TriggerPrimitive& tp : cluster_tps) {
    activity.time_start =
        std::min(activity.time_start, static_cast<uint32_t>(tp.time_start));
    activity.time_end = std::max(
        activity.time_end,
        static_cast<uint32_t>(tp.time_start + tp.samples_over_threshold));
    activity.sadc_sum += tp.adc_integral;
    activity.sadc_max_tps = std::max(activity.sadc_max_tps, tp.adc_integral);
    activity.tot_sum += static_cast<float_t>(tp.samples_over_threshold);
    activity.tot_max = std::max(
        activity.tot_max, static_cast<uint32_t>(tp.samples_over_threshold));
    activity.peak_sum += static_cast<double>(tp.adc_peak);
    activity.peak_max =
        std::max(activity.peak_max, static_cast<uint32_t>(tp.adc_peak));
    total_sadc2_tps += static_cast<double>(tp.adc_integral) *
                       static_cast<double>(tp.adc_integral);
    opdet_ntps[channel_to_opdet(tp.channel)] += 1;
    opdet_sadcs[channel_to_opdet(tp.channel)] += tp.adc_integral;
    times.push_back(tp.time_start);
    TVector3 pos = pmt_info.GetPositionTVec(channel_to_opdet(tp.channel));
    xs.push_back(pos.X());
    ys.push_back(pos.Y());
    zs.push_back(pos.Z());
    if (pmt_info.is_side_opdet(channel_to_opdet(tp.channel))) {
      wall_tp_count += 1;
    }
  }
  activity.n_opdets = opdet_sadcs.size();
  activity.wall_fraction_tps =
      activity.n_tps == 0
          ? 0.0
          : wall_tp_count / static_cast<float_t>(activity.n_tps);
  double total_sadc2_opdets = 0.0;
  uint32_t wall_opdet_count = 0;
  activity.n_tps_max_opdets = 0;
  for (auto const& [opdet_id, ntps] : opdet_ntps) {
    activity.n_tps_max_opdets = std::max(activity.n_tps_max_opdets, ntps);
  }
  activity.sadc_max_opdets = 0;
  for (auto const& [opdet_id, sadc] : opdet_sadcs) {
    double sadc_d = static_cast<double>(sadc);
    total_sadc2_opdets += sadc_d * sadc_d;
    if (pmt_info.is_side_opdet(opdet_id)) {
      wall_opdet_count += 1;
    }
    activity.sadc_max_opdets =
        std::max(activity.sadc_max_opdets, static_cast<float_t>(sadc));
  }
  activity.wall_fraction_opdets =
      activity.n_opdets == 0
          ? 0.0
          : wall_opdet_count / static_cast<float_t>(activity.n_opdets);
  activity.sadc_mean_opdets =
      compute_mean(activity.sadc_sum, activity.n_opdets);
  activity.sadc_mean_tps = compute_mean(activity.sadc_sum, activity.n_tps);
  activity.charge_balance_opdets = compute_charge_balance(
      total_sadc2_opdets, static_cast<double>(activity.sadc_sum),
      activity.n_opdets);
  activity.charge_balance_tps = compute_charge_balance(
      total_sadc2_tps, static_cast<double>(activity.sadc_sum), activity.n_tps);
  activity.tot_mean =
      compute_mean(activity.tot_sum, static_cast<uint32_t>(activity.n_tps));
  activity.peak_mean =
      compute_mean(activity.peak_sum, static_cast<uint32_t>(activity.n_tps));

  compute_metrics(times, activity.t_extent, activity.dt_max, activity.dt_mean,
                  activity.dt_std);
  compute_metrics(xs, activity.x_extent, activity.dx_max, activity.dx_mean,
                  activity.dx_std);
  compute_metrics(ys, activity.y_extent, activity.dy_max, activity.dy_mean,
                  activity.dy_std);
  compute_metrics(zs, activity.z_extent, activity.dz_max, activity.dz_mean,
                  activity.dz_std);
  return activity;
}
