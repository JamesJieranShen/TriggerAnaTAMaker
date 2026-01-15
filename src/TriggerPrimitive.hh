#include <TTreeReader.h>
#include <cstdint>

#pragma once
struct TriggerPrimitive {
  uint8_t version;
  uint8_t flag;
  uint8_t detid;

  uint32_t channel;
  uint16_t samples_over_threshold;
  uint64_t time_start;
  uint16_t samples_to_peak;
  uint32_t adc_integral;
  uint16_t adc_peak;

  int bt_primary_track_id;
  double bt_primary_track_numelectron_frac;
  double bt_primary_track_energy_frac;
  double bt_primary_track_ke;
  uint64_t bt_primary_pdg;
  double bt_edep;
  double bt_numelectrons;
  double bt_x, bt_y, bt_z;
  double bt_primary_x, bt_primary_y, bt_primary_z;
  int bt_mctruth_block_id;
  std::string bt_mctruth_gen_name;

  bool operator==(const TriggerPrimitive& other) const noexcept {
    return std::tie(version, flag, detid, channel, samples_over_threshold,
                    time_start, samples_to_peak, adc_integral, adc_peak) ==
           std::tie(other.version, other.flag, other.detid, other.channel,
                    other.samples_over_threshold, other.time_start,
                    other.samples_to_peak, other.adc_integral, other.adc_peak);
  }

  bool operator!=(const TriggerPrimitive& other) const noexcept {
    return !(*this == other);
  }

  bool operator<(const TriggerPrimitive& other) const noexcept {
    return std::make_tuple(time_start + samples_to_peak, channel) <
           std::make_tuple(other.time_start + other.samples_to_peak, other.channel);
  }

  bool operator>(const TriggerPrimitive& other) const noexcept {
    return other < *this;
  }
  bool operator<=(const TriggerPrimitive& other) const noexcept {
    return !(other < *this);
  }
  bool operator>=(const TriggerPrimitive& other) const noexcept {
    return !(*this < other);
  }
};

struct TriggerPrimitiveReader {
  TTreeReaderValue<uint8_t> version;
  TTreeReaderValue<uint8_t> flag;
  TTreeReaderValue<uint8_t> detid;

  TTreeReaderValue<uint32_t> channel;
  TTreeReaderValue<uint16_t> samples_over_threshold;
  TTreeReaderValue<uint64_t> time_start;
  TTreeReaderValue<uint16_t> samples_to_peak;
  TTreeReaderValue<uint32_t> adc_integral;
  TTreeReaderValue<uint16_t> adc_peak;

  TTreeReaderValue<int> bt_primary_track_id;
  TTreeReaderValue<double> bt_primary_track_numelectron_frac;
  TTreeReaderValue<double> bt_primary_track_energy_frac;
  TTreeReaderValue<double> bt_primary_track_ke;
  TTreeReaderValue<uint64_t> bt_primary_pdg;
  TTreeReaderValue<double> bt_edep;
  TTreeReaderValue<double> bt_numelectrons;
  TTreeReaderValue<double> bt_x;
  TTreeReaderValue<double> bt_y;
  TTreeReaderValue<double> bt_z;
  TTreeReaderValue<double> bt_primary_x;
  TTreeReaderValue<double> bt_primary_y;
  TTreeReaderValue<double> bt_primary_z;
  TTreeReaderValue<int> bt_mctruth_block_id;
  TTreeReaderValue<std::string> bt_mctruth_gen_name;

  explicit TriggerPrimitiveReader(TTreeReader& reader)
      : version(reader, "version"), flag(reader, "flag"),
        detid(reader, "detid"), channel(reader, "channel"),
        samples_over_threshold(reader, "samples_over_threshold"),
        time_start(reader, "time_start"),
        samples_to_peak(reader, "samples_to_peak"),
        adc_integral(reader, "adc_integral"), adc_peak(reader, "adc_peak"),
        bt_primary_track_id(reader, "bt_primary_track_id"),
        bt_primary_track_numelectron_frac(reader,
                                          "bt_primary_track_numelectron_frac"),
        bt_primary_track_energy_frac(reader, "bt_primary_track_energy_frac"),
        bt_primary_track_ke(reader, "bt_primary_track_ke"),
        bt_primary_pdg(reader, "bt_primary_pdg"), bt_edep(reader, "bt_edep"),
        bt_numelectrons(reader, "bt_numelectrons"), bt_x(reader, "bt_x"),
        bt_y(reader, "bt_y"), bt_z(reader, "bt_z"),
        bt_primary_x(reader, "bt_primary_x"),
        bt_primary_y(reader, "bt_primary_y"),
        bt_primary_z(reader, "bt_primary_z"),
        bt_mctruth_block_id(reader, "bt_truth_block_id"),
        bt_mctruth_gen_name(reader, "bt_generator_name") {}

  TriggerPrimitive get() {
    return TriggerPrimitive{
        .version = *version,
        .flag = *flag,
        .detid = *detid,
        .channel = *channel,
        .samples_over_threshold = *samples_over_threshold,
        .time_start = *time_start,
        .samples_to_peak = *samples_to_peak,
        .adc_integral = *adc_integral,
        .adc_peak = *adc_peak,
        .bt_primary_track_id = *bt_primary_track_id,
        .bt_primary_track_numelectron_frac = *bt_primary_track_numelectron_frac,
        .bt_primary_track_energy_frac = *bt_primary_track_energy_frac,
        .bt_primary_track_ke = *bt_primary_track_ke,
        .bt_primary_pdg = *bt_primary_pdg,
        .bt_edep = *bt_edep,
        .bt_numelectrons = *bt_numelectrons,
        .bt_x = *bt_x,
        .bt_y = *bt_y,
        .bt_z = *bt_z,
        .bt_primary_x = *bt_primary_x,
        .bt_primary_y = *bt_primary_y,
        .bt_primary_z = *bt_primary_z,
        .bt_mctruth_block_id = *bt_mctruth_block_id,
        .bt_mctruth_gen_name = *bt_mctruth_gen_name};
  }
};
