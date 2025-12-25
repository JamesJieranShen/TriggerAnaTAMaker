#include "TriggerActivityMakerCluster.hh"
#include "utils.hh"

int compute_max_dt(std::vector<uint64_t>& hit_times) {
  if (hit_times.size() < 2) {
    return 0;
  }
  std::sort(hit_times.begin(), hit_times.end());
  int max_diff = 0;
  for (size_t i = 1; i < hit_times.size(); ++i) {
    int diff = static_cast<int>(hit_times[i] - hit_times[i - 1]);
    if (diff > max_diff) {
      max_diff = diff;
    }
  }
  return max_diff;
}

int compute_max_dt(const std::vector<uint64_t>& hit_times){
  std::vector<uint64_t> hit_times_copy = hit_times;
  return compute_max_dt(hit_times_copy);
}

void
TPBufferCluster::expire(uint64_t current_time) {
  int64_t expire_time = current_time - m_maxClusterDuration;
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
TPBufferCluster::formCluster(
    std::vector<TriggerActivityCluster>& output) const {
  if (m_buffer.empty()) {
    return;
  }
  TriggerActivityCluster activity;
  activity.time = m_currentTime;
  std::vector<TriggerPrimitive> tps;
  TPPriorityQueueCluster temp_queue = m_buffer;
  while (!temp_queue.empty()) {
    tps.push_back(temp_queue.top());
    temp_queue.pop();
  }
  std::map<uint32_t, int64_t> opdet_sadcs;
  activity.ntps = tps.size();
  activity.cathode_ntps = 0;
  activity.side_ntps = 0;
  std::vector<uint64_t> hit_times;
  hit_times.reserve(tps.size());
  std::vector<uint64_t> cathode_hit_times;
  cathode_hit_times.reserve(tps.size());
  std::vector<uint64_t> side_hit_times;
  side_hit_times.reserve(tps.size());
  for (const TriggerPrimitive& tp : tps) {
    int opdet_id = channel_to_opdet(tp.channel);
    hit_times.push_back(tp.time_start);
    if (is_cathode_opdet(opdet_id)) {
      activity.cathode_ntps += 1;
      cathode_hit_times.push_back(tp.time_start);
    } else {
      activity.side_ntps += 1;
      side_hit_times.push_back(tp.time_start);
    }
    opdet_sadcs[opdet_id] += tp.adc_integral;
  }
  int64_t total_sadc = 0;
  int64_t total_sadc_cathode = 0;
  double total_sadc2 = 0;
  double total_sadc2_cathode = 0;
  int n_cathode_hits = 0;
  for (const auto& [opdet_id, sadc] : opdet_sadcs) {
    bool is_cathode = is_cathode_opdet(opdet_id);
    total_sadc += sadc;
    total_sadc2 += sadc * sadc;
    if (is_cathode) {
      total_sadc_cathode += sadc;
      total_sadc2_cathode += sadc * sadc;
      n_cathode_hits += 1;
    }
  }
  activity.nhits = opdet_sadcs.size();
  activity.sadc = total_sadc;
  activity.charge_balance =
      compute_charge_balance(total_sadc2, total_sadc, activity.nhits);
  activity.max_dt = compute_max_dt(hit_times);

  activity.cathode_nhits = n_cathode_hits;
  activity.cathode_sadc = total_sadc_cathode;
  activity.cathode_charge_balance =
      compute_charge_balance(total_sadc2_cathode, total_sadc_cathode,
                             activity.cathode_nhits);
  activity.cathode_max_dt = compute_max_dt(cathode_hit_times);

  activity.side_nhits = activity.nhits - n_cathode_hits;
  activity.side_sadc = total_sadc - total_sadc_cathode;
  activity.side_charge_balance = compute_charge_balance(
      total_sadc2 - total_sadc2_cathode, total_sadc - total_sadc_cathode,
      activity.side_nhits);
  activity.side_max_dt = compute_max_dt(side_hit_times);

  output.push_back(activity);
}
