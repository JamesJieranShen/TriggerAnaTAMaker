#include "TriggerActivity.hh"
#include "TriggerActivityMakerCluster.hh"
#include "TriggerActivityMakerMTCA.hh"
#include "TriggerPrimitive.hh"

#include <CLI/CLI.hpp>
#include <TFile.h>
#include <TTreeReader.h>
#include <iostream>
#include <string>

int
main(int argc, char** argv) {
  CLI::App parser{
      "Toy Trigger Activity Maker, using TriggerAnaTree as an input."};
  std::string input_file, output_file, tree_name;
  parser.add_option("--input,-i", input_file, "Input TriggerAnaTree file")
      ->required();
  parser.add_option("--tree", tree_name, "Trigger Activity Tree Name")
      ->required();
  parser.add_option("--output,-o", output_file, "Output file")->required();
  CLI11_PARSE(parser, argc, argv);

  std::cout << "Input File: " << input_file << std::endl;
  std::cout << "Tree Name: " << tree_name << std::endl;
  std::cout << "Output File: " << output_file << std::endl;

  TFile fin(input_file.c_str(), "READ");
  TTreeReader tp_tree_reader(tree_name.c_str(), &fin);
  TTreeReaderValue<uint> event_reader(tp_tree_reader, "event");
  TTreeReaderValue<uint> subrun_reader(tp_tree_reader, "subrun");
  TTreeReaderValue<uint> run_reader(tp_tree_reader, "run");
  TriggerPrimitiveReader tp_reader(tp_tree_reader);
  // TriggerActivityMakerMTCA ta_maker(1, 10);
  TriggerActivityMakerCluster ta_maker(32, 23, 11, 727, 2);

  TFile fout(output_file.c_str(), "RECREATE");
  TTree out_tree("TriggerActivities", "Trigger Activities");
  int event_counter = -1;
  int prev_event = -99999;
  int prev_subrun = -99999;
  int prev_run = -99999;
  auto output_tas = ta_maker.get_empty_vector();
  auto ta_buf = ta_maker.get_empty_activity();
  out_tree.Branch("event", &event_counter);
  out_tree.Branch("orig_event", &prev_event);
  out_tree.Branch("subrun", &prev_subrun);
  out_tree.Branch("run", &prev_run);
  ta_buf.LinkTree(out_tree);

  int curr_event = -99999;
  int curr_subrun = -99999;
  int curr_run = -99999;
  std::vector<TriggerPrimitive> tps_in_event;
  while (tp_tree_reader.Next()) {
    curr_event = *event_reader;
    curr_subrun = *subrun_reader;
    curr_run = *run_reader;
    if ((curr_event != prev_event) || (curr_subrun != prev_subrun) ||
        (curr_run != prev_run)) {
      std::sort(tps_in_event.begin(), tps_in_event.end(),
                [](const TriggerPrimitive& a, const TriggerPrimitive& b) {
                  return std::tie(a.time_start, a.channel) <
                         std::tie(b.time_start, b.channel);
                });
      for (const TriggerPrimitive& curr_tp : tps_in_event) {
        ta_maker(curr_tp, output_tas);
      }
      ta_maker.flush(output_tas);
      std::cout << "Writing TAs for event " << event_counter << std::endl;
      // std::cout << "Current Event: " << curr_event << std::endl;
      std::cout << "  Number of TPs: " << tps_in_event.size() << std::endl;
      std::cout << "  Number of TAs: " << output_tas.size() << std::endl;
      for (const auto& ta : output_tas) {
        ta_buf = ta;
        out_tree.Fill();
      }
      tps_in_event.clear();
      output_tas.clear();
      event_counter++;
      prev_event = curr_event;
      prev_subrun = curr_subrun;
      prev_run = curr_run;
    }
    TriggerPrimitive curr_tp = tp_reader.get();
    if (curr_tp.time_start > 1e10) continue;
    tps_in_event.push_back(curr_tp);
  }
  std::sort(tps_in_event.begin(), tps_in_event.end(),
            [](const TriggerPrimitive& a, const TriggerPrimitive& b) {
              return std::tie(a.time_start, a.channel) <
                     std::tie(b.time_start, b.channel);
            });
  for (const TriggerPrimitive& curr_tp : tps_in_event) {
    ta_maker(curr_tp, output_tas);
  }
  ta_maker.flush(output_tas);
  std::cout << "Writing TAs for event " << event_counter << std::endl;
  std::cout << "  Number of TAs: " << output_tas.size() << std::endl;
  for (const auto& ta : output_tas) {
    ta_buf = ta;
    out_tree.Fill();
  }
  out_tree.Write();
  fout.Close();

  return 0;
}
