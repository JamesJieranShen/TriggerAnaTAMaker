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

  bool operator==(const TriggerPrimitive&) const = default;
  auto operator<=>(const TriggerPrimitive& other) const {
    return std::tie(time_start, channel) <=>
           std::tie(other.time_start, other.channel);
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

  explicit TriggerPrimitiveReader(TTreeReader& reader)
      : version(reader, "version"), flag(reader, "flag"),
        detid(reader, "detid"), channel(reader, "channel"),
        samples_over_threshold(reader, "samples_over_threshold"),
        time_start(reader, "time_start"),
        samples_to_peak(reader, "samples_to_peak"),
        adc_integral(reader, "adc_integral"), adc_peak(reader, "adc_peak") {}

  TriggerPrimitive get() {
    return TriggerPrimitive{.version = *version,
                            .flag = *flag,
                            .detid = *detid,
                            .channel = *channel,
                            .samples_over_threshold = *samples_over_threshold,
                            .time_start = *time_start,
                            .samples_to_peak = *samples_to_peak,
                            .adc_integral = *adc_integral,
                            .adc_peak = *adc_peak};
  }
};
