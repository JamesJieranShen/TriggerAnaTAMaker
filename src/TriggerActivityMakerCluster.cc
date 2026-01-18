#include "TriggerActivityMakerCluster.hh"

#include "utils.hh"

#include <TVector3.h>

// assumes tps are sorted by peak time.
std::vector<std::vector<TriggerPrimitive>>
split_on_timediff(const std::vector<TriggerPrimitive>& tps,
                  int64_t max_cluster_time, int64_t max_hit_time_diff) {
  std::vector<std::vector<TriggerPrimitive>> clusters;
  if (tps.empty()) {
    return clusters;
  }
  std::vector<TriggerPrimitive> current_cluster;
  current_cluster.push_back(tps.at(0));
  for (size_t i = 1; i < tps.size(); ++i) {
    if (
        // The current TP is too late, exceeding max time diff
        (static_cast<int64_t>(tps.at(i).PeakTime() -
                              current_cluster.back().PeakTime()) >
         max_hit_time_diff) ||
        // The current cluster will exceed max cluster time with the addition of
        // the current TP
        (static_cast<int64_t>(tps.at(i).PeakTime() -
                              current_cluster.front().PeakTime()) >
         max_cluster_time)) {
      clusters.push_back(current_cluster);
      current_cluster.clear();
    }
    current_cluster.push_back(tps.at(i));
  }
  if (!current_cluster.empty()) {
    clusters.push_back(current_cluster);
  }
  return clusters;
}

// assumes tps are sorted by peak time.
std::vector<std::vector<TriggerPrimitive>>
split_on_maxtime(const std::vector<TriggerPrimitive>& tps,
                 int64_t max_cluster_time) {
  std::vector<std::vector<TriggerPrimitive>> clusters;
  if (tps.empty()) {
    return clusters;
  }
  std::vector<TriggerPrimitive> current_cluster;
  current_cluster.push_back(tps.at(0));
  for (size_t i = 1; i < tps.size(); ++i) {
    if (static_cast<int64_t>(tps.at(i).PeakTime() -
                             current_cluster.front().PeakTime()) >
        max_cluster_time) {
      clusters.push_back(current_cluster);
      current_cluster.clear();
    }
    current_cluster.push_back(tps.at(i));
  }
  if (!current_cluster.empty()) {
    clusters.push_back(current_cluster);
  }
  return clusters;
}

std::vector<TriggerPrimitive>
filter_by_distance(std::vector<TriggerPrimitive>& cluster_tps,
                   float_t max_hit_distance, int min_neighbors) {
  // removes tps from the cluster_tps list if there are no neighbors within
  // max_hit_distance.

  // do nothing if parameters are invalid
  if (max_hit_distance <= 0.0f || min_neighbors <= 0) {
    return {};
  }
  std::vector<TriggerPrimitive> filtered_tps;
  std::vector<TriggerPrimitive> removed_tps;
  for (const TriggerPrimitive& curr_tp : cluster_tps) {
    TVector3 curr_pos = curr_tp.OpDetPosition();
    int neighbor_count = 0;

    for (const TriggerPrimitive& other_tp : cluster_tps) {
      if (curr_tp == other_tp) continue;
      float_t distance2 = (curr_pos - other_tp.OpDetPosition()).Mag2();
      if (distance2 <= max_hit_distance * max_hit_distance) {
        neighbor_count++;
        if (neighbor_count >= min_neighbors) {
          break;
        }
      }
    } // end for other_tp
    if (neighbor_count >= min_neighbors) {
      filtered_tps.push_back(curr_tp);
    } else {
      removed_tps.push_back(curr_tp);
    }
  } // end for curr_tp
  std::swap(cluster_tps, filtered_tps);
  return removed_tps;
}

template <typename T>
void
compute_metrics(std::vector<T> vals, float_t& extent_val, float_t& delta_max,
                float_t& delta_mean, float_t& delta_std) {
  if (vals.size() < 2) {
    delta_max = 0.0f;
    delta_mean = 0.0f;
    return;
  }
  std::sort(vals.begin(), vals.end());
  extent_val =
      static_cast<float_t>(vals.back()) - static_cast<float_t>(vals.front());
  float_t sum_deltas = 0;
  delta_max = 0;
  for (size_t i = 1; i < vals.size(); ++i) {
    float_t delta = static_cast<float_t>(vals[i] - vals[i - 1]);
    sum_deltas += delta;
    delta_max = std::max(delta_max, delta);
  }
  delta_mean = safe_divide(sum_deltas, vals.size() - 1);
  delta_std = compute_std(static_cast<double>(sum_deltas * sum_deltas),
                          static_cast<double>(delta_mean), vals.size() - 1);
}

float_t
compute_mean_distance(const std::vector<TriggerPrimitive>& tps) {
  if (tps.size() < 2) {
    return 0.0f;
  }
  size_t n_pairs = 0;
  double total_distance = 0.0;
  for (size_t i = 1; i < tps.size(); ++i) {
    const TriggerPrimitive& prev_tp = tps.at(i - 1);
    const TriggerPrimitive& curr_tp = tps.at(i);
    float_t distance =
        (curr_tp.OpDetPosition() - prev_tp.OpDetPosition()).Mag();
    total_distance += distance;
    n_pairs += 1;
  }
  return safe_divide(total_distance, n_pairs);
}

void
TPBufferCluster::expire(uint64_t current_time) {
  int64_t expire_time = current_time - m_bufLength;
  while (!m_buffer.empty() &&
         static_cast<int64_t>(m_buffer.top().time_start) < expire_time) {
    m_buffer.pop();
  }
}

void
TPBufferCluster::add(const TriggerPrimitive& tp) {
  expire(tp.time_start);
  m_buffer.push(tp);
}
void
TPBufferCluster::formClusters(std::vector<TriggerActivityCluster>& output,
                              int64_t max_cluster_time,
                              int64_t max_hit_time_diff, int64_t min_nhits,
                              float_t max_hit_distance, int64_t min_neighbors,
                              bool last_call) {
  if (m_buffer.empty()) return;
  std::vector<TriggerPrimitive> tps;
  while (!m_buffer.empty()) {
    tps.push_back(m_buffer.top());
    m_buffer.pop();
  }
  std::sort(tps.begin(), tps.end(),
            [](const TriggerPrimitive& a, const TriggerPrimitive& b) {
              return a.PeakTime() < b.PeakTime();
            });
  std::vector<std::vector<TriggerPrimitive>> time_clusters =
      // Pablo Style. No Cluster merging between windows.
      // split_on_maxtime(tps, max_cluster_time);
      // --------
      split_on_timediff(tps, max_cluster_time, max_hit_time_diff);
      // --------
  if (!last_call) {
    // keep the last cluster in the buffer to be processed with future tps.
    std::vector<TriggerPrimitive> last_cluster = time_clusters.back();
    time_clusters.pop_back();
    for (const TriggerPrimitive& tp : last_cluster) {
      m_buffer.push(tp);
    }
  }
  for (std::vector<TriggerPrimitive>& cluster_tps : time_clusters) {
    // TODO: gather all removed TPs and add back to buffer?
    if (cluster_tps.size() < static_cast<size_t>(min_nhits)) {
      continue;
    }
    // Pablo Style:
    // std::vector<std::vector<TriggerPrimitive>> subclusters =
    //     split_on_timediff(cluster_tps, max_cluster_time, max_hit_time_diff);
    // for (std::vector<TriggerPrimitive>& subcluster_tps : subclusters) {
    //   if (subcluster_tps.size() < static_cast<size_t>(min_nhits)) {
    //     continue;
    //   }
    //   filter_by_distance(subcluster_tps, max_hit_distance, min_neighbors);
    //   if (subcluster_tps.size() < static_cast<size_t>(min_nhits)) {
    //     continue;
    //   }
    //   std::sort(subcluster_tps.begin(), subcluster_tps.end());
    //   output.push_back(makeTriggerActivity(subcluster_tps));
    //   if (m_verbosity >= Verbosity::kDebug)
    //     std::cout << "  Adding TA: [ " << output.back().time_start << ", "
    //               << output.back().time_end << " ]" << std::endl;
    //   output.back().BackTrack(subcluster_tps);
    // }
    // --------------------
    filter_by_distance(cluster_tps, max_hit_distance, min_neighbors);
    if (cluster_tps.size() < static_cast<size_t>(min_nhits)) {
      continue;
    }
    std::sort(cluster_tps.begin(), cluster_tps.end());
    output.push_back(makeTriggerActivity(cluster_tps));
    if (m_verbosity >= Verbosity::kDebug)
      std::cout << "  Adding TA: [ " << output.back().time_start << ", "
                << output.back().time_end << " ]" << std::endl;
    output.back().BackTrack(cluster_tps);
    // --------------------
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
  std::map<uint32_t, bool> opdet_onCathode;
  // PMTInfo& pmt_info = PMTInfo::Instance();

  double total_sadc2_tps = 0.0;
  std::vector<uint64_t> times;
  times.reserve(cluster_tps.size());
  std::vector<float_t> xs;
  xs.reserve(cluster_tps.size());
  std::vector<float_t> ys;
  ys.reserve(cluster_tps.size());
  std::vector<float_t> zs;
  zs.reserve(cluster_tps.size());
  std::vector<double_t> qs;
  uint32_t wall_tp_count = 0;
  for (const TriggerPrimitive& tp : cluster_tps) {
    activity.time_start =
        std::min(activity.time_start, static_cast<uint32_t>(tp.PeakTime()));
    activity.time_end =
        std::max(activity.time_end, static_cast<uint32_t>(tp.PeakTime()));
    activity.sadc_sum += tp.adc_integral;
    qs.push_back(tp.adc_integral);
    activity.sadc_max_tps = std::max(activity.sadc_max_tps, tp.adc_integral);
    activity.tot_sum += static_cast<float_t>(tp.samples_over_threshold);
    activity.tot_max = std::max(
        activity.tot_max, static_cast<uint32_t>(tp.samples_over_threshold));
    activity.peak_sum += static_cast<double>(tp.adc_peak);
    activity.peak_max =
        std::max(activity.peak_max, static_cast<uint32_t>(tp.adc_peak));
    total_sadc2_tps += static_cast<double>(tp.adc_integral) *
                       static_cast<double>(tp.adc_integral);
    uint32_t opdet_id = tp.OpDetID();
    opdet_ntps[opdet_id] += 1;
    opdet_sadcs[opdet_id] += tp.adc_integral;
    opdet_onCathode[opdet_id] = tp.OnCathode();
    times.push_back(tp.time_start);
    TVector3 pos = tp.OpDetPosition();
    xs.push_back(pos.X());
    ys.push_back(pos.Y());
    zs.push_back(pos.Z());
    if (!tp.OnCathode()) {
      wall_tp_count += 1;
    }
  }
  activity.n_opdets = opdet_sadcs.size();
  activity.wall_fraction_tps = safe_divide(wall_tp_count, activity.n_tps);
  double total_sadc2_opdets = 0.0;
  uint32_t wall_opdet_count = 0;
  float_t wall_sadc_sum = 0.0f;
  activity.n_tps_max_opdets = 0;
  for (auto const& [opdet_id, ntps] : opdet_ntps) {
    activity.n_tps_max_opdets = std::max(activity.n_tps_max_opdets, ntps);
  }
  activity.sadc_max_opdets = 0;
  for (auto const& [opdet_id, sadc] : opdet_sadcs) {
    double sadc_d = static_cast<double>(sadc);
    total_sadc2_opdets += sadc_d * sadc_d;
    if (!opdet_onCathode.at(opdet_id)) {
      wall_opdet_count += 1;
      wall_sadc_sum += static_cast<float_t>(sadc);
    }
    activity.sadc_max_opdets =
        std::max(activity.sadc_max_opdets, static_cast<float_t>(sadc));
  }
  activity.wall_fraction_opdets =
      safe_divide(wall_opdet_count, activity.n_opdets);
  activity.wall_fraction_sadc = safe_divide(wall_sadc_sum, activity.sadc_sum);
  activity.sadc_mean_opdets = safe_divide(activity.sadc_sum, activity.n_opdets);
  activity.sadc_mean_tps = safe_divide(activity.sadc_sum, activity.n_tps);
  activity.charge_balance_opdets = compute_charge_balance(
      total_sadc2_opdets, static_cast<double>(activity.sadc_sum),
      activity.n_opdets);
  activity.charge_balance_tps = compute_charge_balance(
      total_sadc2_tps, static_cast<double>(activity.sadc_sum), activity.n_tps);
  activity.tot_mean =
      safe_divide(activity.tot_sum, static_cast<uint32_t>(activity.n_tps));
  activity.peak_mean =
      safe_divide(activity.peak_sum, static_cast<uint32_t>(activity.n_tps));

  compute_metrics(times, activity.t_extent, activity.dt_max, activity.dt_mean,
                  activity.dt_std);
  compute_metrics(xs, activity.x_extent, activity.dx_max, activity.dx_mean,
                  activity.dx_std);
  compute_metrics(ys, activity.y_extent, activity.dy_max, activity.dy_mean,
                  activity.dy_std);
  compute_metrics(zs, activity.z_extent, activity.dz_max, activity.dz_mean,
                  activity.dz_std);
  double weights = 0;
  double x_sum = 0;
  double y_sum = 0;
  double z_sum = 0;
  for (size_t i = 0; i < xs.size(); i++) {
    // if (xs.at(i) > 0) continue; // skip wall PDs
    // if (xs.at(i) < 0) continue; // skip cathode PDs
    weights += qs.at(i);
    x_sum += qs.at(i) * xs.at(i);
    y_sum += qs.at(i) * ys.at(i);
    z_sum += qs.at(i) * zs.at(i);
  }
  activity.reco_x = safe_divide(x_sum, weights, -99999.0f);
  activity.reco_y = safe_divide(y_sum, weights, -99999.0f);
  activity.reco_z = safe_divide(z_sum, weights, -99999.0f);
  return activity;
}

TriggerActivityMakerCluster::TriggerActivityMakerCluster(
    const nlohmann::json& cfg)
    : m_tpBuffer(6250), m_maxClusterTime(cfg["max_cluster_time"]),
      m_maxHitTimeDiff(cfg["max_hit_time_diff"]), m_minNhits(cfg["min_nhits"]),
      m_maxHitDistance(cfg["max_hit_distance"]),
      m_minNeighbors(cfg.value("min_neighbors", 1)) {
  if (cfg.contains("filter_rules")) {
    const nlohmann::json& filter_cfg = cfg["filter_rules"];
    if (filter_cfg.contains("min_width")) {
      m_tpMinWidth = filter_cfg["min_width"];
      std::cout << "TP min width set to " << m_tpMinWidth << " from config."
                << std::endl;
    }
    if (filter_cfg.contains("min_sadc")) {
      m_tpMinSADC = filter_cfg["min_sadc"];
      std::cout << "TP min SADC set to " << m_tpMinSADC << " from config."
                << std::endl;
    }
    if (filter_cfg.contains("include_wall_PDs")) {

      m_tpIncludeWallPDs = filter_cfg["include_wall_PDs"];
      std::cout << "TP include wall PDs set to " << m_tpIncludeWallPDs
                << " from config." << std::endl;
    }
    if (filter_cfg.contains("include_cathode_PDs")) {
      m_tpIncludeCathodePDs = !filter_cfg["include_cathode_PDs"];
      std::cout << "TP include cathode PDs set to " << m_tpIncludeCathodePDs
                << " from config." << std::endl;
    }
  }
}

void
TriggerActivityMakerCluster::operator()(
    const TriggerPrimitive& input_tp,
    std::vector<TriggerActivityCluster>& output_ta) {
  if (input_tp.samples_over_threshold < m_tpMinWidth) {
    if (m_verbosity >= Verbosity::kDebug) {
      std::cout << "  Skipping TP due to width: "
                << input_tp.samples_over_threshold << " < " << m_tpMinWidth
                << std::endl;
    }
    return;
  }
  if (input_tp.adc_integral < m_tpMinSADC) {
    if (m_verbosity >= Verbosity::kDebug) {
      std::cout << "  Skipping TP due to SADC: " << input_tp.adc_integral
                << " < " << m_tpMinSADC << std::endl;
    }
    return;
  }
  if (!m_tpIncludeWallPDs && input_tp.OnWall()) {
    if (m_verbosity >= Verbosity::kDebug) {
      std::cout << "  Skipping TP on wall PD: " << input_tp.OpDetID()
                << std::endl;
    }
    return;
  }
  if (!m_tpIncludeCathodePDs && input_tp.OnCathode()) {
    if (m_verbosity >= Verbosity::kDebug) {
      std::cout << "  Skipping TP on cathode PD: " << input_tp.OpDetID()
                << std::endl;
    }
    return;
  }
  if (m_verbosity >= Verbosity::kDebug) {
    std::cout << "Adding TP: " << input_tp.time_start << " "
              << input_tp.PeakTime() << " " << input_tp.EndTime() << std::endl;
  }
  m_tpBuffer.formClusters(output_ta, m_maxClusterTime, m_maxHitTimeDiff,
                          m_minNhits, m_maxHitDistance, m_minNeighbors);
  m_tpBuffer.add(input_tp);
}
