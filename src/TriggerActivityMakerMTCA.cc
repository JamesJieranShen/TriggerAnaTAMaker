#include <TriggerActivityMakerMTCA.hh>
#include <utils.hh>


void
compute_metrics(double sum2, double sum, int n, uint32_t& nhits, float_t& sadc,
                float_t& adc_mean, float_t& adc_std, float_t& charge_balance) {
  nhits = n;
  sadc = sum;
  adc_mean = n > 0 ? static_cast<float_t>(sum) / n : 0.0f;
  adc_std = compute_std(sum2, adc_mean, n);
  charge_balance =
      adc_mean > 0.0f ? static_cast<float_t>(adc_std) / adc_mean : 0.0f;
}


void
TPBufferMTCA::add(const TriggerPrimitive& tp) {
  m_currentTime = tp.time_start;
  expire(tp.time_start);
  m_buffer.push(tp);
}

void
TPBufferMTCA::expire(uint64_t current_time) {
  int64_t expire_time = current_time - m_windowSizeLong + 1;
  while (!m_buffer.empty() &&
         static_cast<int64_t>(m_buffer.top().time_start +
                              m_buffer.top().samples_over_threshold) <=
             expire_time) {
    m_buffer.pop();
  }
}

TriggerActivityMTCA
TPBufferMTCA::computeActivity() const {
  TriggerActivityMTCA activity;
  activity.time = m_currentTime;
  std::vector<TriggerPrimitive> tps;
  std::map<uint32_t, float_t> opdet_sadcs_short;
  std::map<uint32_t, float_t> opdet_sadcs_long;
  TPPriorityQueueMTCA temp_queue = m_buffer;
  while (!temp_queue.empty()) {
    tps.push_back(temp_queue.top());
    temp_queue.pop();
  }

  int64_t short_time_cutoff = m_currentTime - m_windowSizeShort + 1;
  int64_t long_time_cutoff = m_currentTime - m_windowSizeLong + 1;
  // std::cout << "==========" << std::endl;
  // std::cout << "Current time: " << m_currentTime << std::endl;
  for (const TriggerPrimitive& tp : tps) {
    float mean_charge =
        static_cast<float>(tp.adc_integral) / tp.samples_over_threshold;
    bool in_short =
        static_cast<int64_t>(tp.time_start + tp.samples_over_threshold) >
        short_time_cutoff;
    uint opdet_id =
        tp.channel / 10; // each PMT has two channels, they will be XX0 and XX1

    // if (!in_long) {
    // std::cout << "Current Time: " << m_currentTime << " TP on opdet "
    //           << opdet_id << " at time " << tp.time_start
    //           << " samples over threshold " << tp.samples_over_threshold
    //           << " with mean charge " << mean_charge
    //           << " -- in short window: " << in_short << std::endl;
    // }

    {
      int64_t start =
          std::max(static_cast<int64_t>(tp.time_start), long_time_cutoff);
      int64_t end = std::min(tp.time_start + tp.samples_over_threshold,
                             m_currentTime + 1);
      opdet_sadcs_long[opdet_id] += mean_charge * (end - start);
    }
    if (in_short) {
      int64_t start =
          std::max(static_cast<int64_t>(tp.time_start), short_time_cutoff);
      int64_t end = std::min(tp.time_start + tp.samples_over_threshold,
                             m_currentTime + 1);
      opdet_sadcs_short[opdet_id] += mean_charge * (end - start);
    }
  }
  double total_sadc_long = 0;
  double total_sadc2_long = 0;
  double total_sadc_cathode_long = 0;
  double total_sadc2_cathode_long = 0;
  int n_cathode_hits_long = 0;
  for (const auto& [opdet_id, sadc] : opdet_sadcs_long) {
    bool is_cathode = opdet_id >= 40;
    total_sadc_long += sadc;
    // std::cout << "Opdet " << opdet_id << " SADC Long: " << sadc << " Sum so
    // far: " << total_sadc_long << std::endl;
    total_sadc2_long += sadc * sadc;
    if (is_cathode) {
      total_sadc_cathode_long += sadc;
      total_sadc2_cathode_long += sadc * sadc;
      n_cathode_hits_long += 1;
    }
  }

  compute_metrics(total_sadc2_long, total_sadc_long, opdet_sadcs_long.size(),
                  activity.nhits_long, activity.sadc_long,
                  activity.adc_mean_long, activity.adc_std_long,
                  activity.charge_balance_long);

  compute_metrics(total_sadc2_cathode_long, total_sadc_cathode_long,
                  n_cathode_hits_long, activity.cathode_nhits_long,
                  activity.cathode_sadc_long, activity.cathode_adc_mean_long,
                  activity.cathode_adc_std_long,
                  activity.cathode_charge_balance_long);

  compute_metrics(total_sadc2_long - total_sadc_cathode_long,
                  total_sadc_long - total_sadc_cathode_long,
                  activity.nhits_long - activity.cathode_nhits_long,
                  activity.side_nhits_long, activity.side_sadc_long,
                  activity.side_adc_mean_long, activity.side_adc_std_long,
                  activity.side_charge_balance_long);

  double total_sadc_short = 0;
  double total_sadc2_short = 0;
  double total_sadc_cathode_short = 0;
  double total_sadc2_cathode_short = 0;
  int n_cathode_hits_short = 0;
  for (const auto& [opdet_id, sadc] : opdet_sadcs_short) {
    bool is_cathode = opdet_id >= 40;
    total_sadc_short += sadc;
    total_sadc2_short += sadc * sadc;
    if (is_cathode) {
      total_sadc_cathode_short += sadc;
      total_sadc2_cathode_short += sadc * sadc;
      n_cathode_hits_short += 1;
    }
  }

  compute_metrics(total_sadc2_short, total_sadc_short, opdet_sadcs_short.size(),
                  activity.nhits_short, activity.sadc_short,
                  activity.adc_mean_short, activity.adc_std_short,
                  activity.charge_balance_short);

  compute_metrics(total_sadc2_cathode_short, total_sadc_cathode_short,
                  n_cathode_hits_short, activity.cathode_nhits_short,
                  activity.cathode_sadc_short, activity.cathode_adc_mean_short,
                  activity.cathode_adc_std_short,
                  activity.cathode_charge_balance_short);

  compute_metrics(total_sadc2_short - total_sadc_cathode_short,
                  total_sadc_short - total_sadc_cathode_short,
                  activity.nhits_short - activity.cathode_nhits_short,
                  activity.side_nhits_short, activity.side_sadc_short,
                  activity.side_adc_mean_short, activity.side_adc_std_short,
                  activity.side_charge_balance_short);
  // std::cout << "SADC Long: " << activity.sadc_long << std::endl;
  // std::cout << "==========" << std::endl;
  return activity;
}
