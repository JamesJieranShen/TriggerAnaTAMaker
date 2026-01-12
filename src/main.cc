#include "TriggerActivity.hh"
#include "TriggerActivityMakerCluster.hh"
#include "TriggerActivityMakerMTCA.hh"
#include "TriggerPrimitive.hh"

#include <CLI/CLI.hpp>
#include <TFile.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TObjString.h>
#include <iostream>
#include <string>
#include <filesystem>
#include <nlohmann/json.hpp>

int
main(int argc, char** argv) {
  CLI::App parser{
      "Toy Trigger Activity Maker, using TriggerAnaTree as an input."};
  std::string input_file, output_file, cfg, tree_name;
  parser.add_option("--config", cfg, "JSON Config string or filename.")->required();
  parser.add_option("--input,-i", input_file, "Input TriggerAnaTree file");
  parser.add_option("--tree", tree_name, "Trigger Activity Tree Name");
  parser.add_option("--output,-o", output_file, "Output file");
  CLI11_PARSE(parser, argc, argv);

  nlohmann::json parsed_cfg;
  if (std::filesystem::exists(cfg) && std::filesystem::is_regular_file(cfg)) {
    std::ifstream cfg_file(cfg);
    if (!cfg_file) throw std::runtime_error("Could not open config file: " + cfg);
    parsed_cfg = nlohmann::json::parse(cfg_file, nullptr, true, true);
  }
  else {
    parsed_cfg = nlohmann::json::parse(cfg, nullptr, true, true);
  }
  if (input_file.empty()) {
    input_file = parsed_cfg.value("input", "");
  }
  if (tree_name.empty()) {
    tree_name = parsed_cfg.value("tree", "");
  }
  if (output_file.empty()) {
    output_file = parsed_cfg.value("output", "");
  }
  if (input_file.empty()) throw std::runtime_error("Input file not specified.");
  if (tree_name.empty()) throw std::runtime_error("Tree name not specified.");
  if (output_file.empty()) throw std::runtime_error("Output file not specified.");

  
  std::cout << "Config:\n" << parsed_cfg.dump(2) << std::endl;
  std::cout << "Input File: " << input_file << std::endl;
  std::cout << "Tree Name: " << tree_name << std::endl;
  std::cout << "Output File: " << output_file << std::endl;
  TChain *theChain = new TChain(tree_name.c_str());
  theChain->Add(input_file.c_str());
  TTreeReader tp_tree_reader(theChain);
  TTreeReaderValue<uint> event_reader(tp_tree_reader, "event");
  TTreeReaderValue<uint> subrun_reader(tp_tree_reader, "subrun");
  TTreeReaderValue<uint> run_reader(tp_tree_reader, "run");
  TriggerPrimitiveReader tp_reader(tp_tree_reader);
  // TriggerActivityMakerMTCA ta_maker(1, 10);
  nlohmann::json ta_cfg = parsed_cfg["TAMaker"];
  if (ta_cfg["tool_type"] != "cluster") {
    throw std::runtime_error("Only TriggerActivityMakerCluster is supported in this example.");
  }
  // TriggerActivityMakerCluster ta_maker(32, 23, 11, 727, 2);
  TriggerActivityMakerCluster ta_maker(ta_cfg["config"]);

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

  nlohmann::json provenance = parsed_cfg;
  provenance["input_file"] = input_file;
  provenance["tree_name"] = tree_name;
  TObjString provenance_str = provenance.dump().c_str();
  provenance_str.Write("provenance",  TObject::kOverwrite);
  fout.Close();

  return 0;
}
