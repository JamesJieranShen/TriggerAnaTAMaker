#pragma once
#include <TTree.h>
#include <TriggerPrimitive.hh>
#include <map>
#include <numeric>
#include <string>

struct TriggerActivity {
public:
  virtual void LinkTree(TTree& tree, const std::string& prefix = "") = 0;
  virtual ~TriggerActivity() = default;

  virtual void LinkTreeBackTracker(TTree& tree,
                                   const std::string& prefix = "") {
    tree.Branch((prefix + "bt_primary_x").c_str(), &bt_primary_x);
    tree.Branch((prefix + "bt_primary_y").c_str(), &bt_primary_y);
    tree.Branch((prefix + "bt_primary_z").c_str(), &bt_primary_z);
    tree.Branch((prefix + "bt_primary_pdg").c_str(), &bt_primary_pdg);
    tree.Branch((prefix + "bt_primary_edep").c_str(), &bt_primary_edep);
    tree.Branch((prefix + "bt_total_edep").c_str(), &bt_total_edep);
    tree.Branch((prefix + "bt_primary_numelectrons").c_str(),
                &bt_primary_numelectrons);
    tree.Branch((prefix + "bt_total_numelectrons").c_str(),
                &bt_total_numelectrons);
    tree.Branch((prefix + "bt_primary_ke").c_str(), &bt_primary_ke);
  }

  virtual void BackTrack(const std::vector<TriggerPrimitive>& tps) {
    std::map<int, float_t> trk_edep_map;
    std::map<int, float_t> trk_numelectrons_map;
    std::map<int, uint64_t> trk_pdg_map;
    std::map<int, float_t> trk_ke_map;
    std::map<int, std::vector<float_t>> trk_posx_map;
    std::map<int, std::vector<float_t>> trk_posy_map;
    std::map<int, std::vector<float_t>> trk_posz_map;
    for (const TriggerPrimitive& tp : tps) {
      if (tp.bt_primary_track_id == -99999) continue;
      trk_edep_map[tp.bt_primary_track_id] +=
          tp.bt_edep * tp.bt_primary_track_energy_frac;
      trk_numelectrons_map[tp.bt_primary_track_id] +=
          tp.bt_numelectrons * tp.bt_primary_track_numelectron_frac;
      trk_pdg_map[tp.bt_primary_track_id] = tp.bt_primary_pdg;
      trk_ke_map[tp.bt_primary_track_id] = tp.bt_primary_track_ke;
      trk_posx_map[tp.bt_primary_track_id].push_back(tp.bt_primary_x);
      trk_posy_map[tp.bt_primary_track_id].push_back(tp.bt_primary_y);
      trk_posz_map[tp.bt_primary_track_id].push_back(tp.bt_primary_z);
    }

    int bt_primary_track_id =
        std::max_element(
            trk_numelectrons_map.begin(), trk_numelectrons_map.end(),
            [](const auto& a, const auto& b) { return a.second < b.second; })
            ->first;
    bt_primary_x =
        std::accumulate(trk_posx_map[bt_primary_track_id].begin(),
                        trk_posx_map[bt_primary_track_id].end(), 0.0) /
        trk_posx_map[bt_primary_track_id].size();
    bt_primary_y =
        std::accumulate(trk_posy_map[bt_primary_track_id].begin(),
                        trk_posy_map[bt_primary_track_id].end(), 0.0) /
        trk_posy_map[bt_primary_track_id].size();
    bt_primary_z =
        std::accumulate(trk_posz_map[bt_primary_track_id].begin(),
                        trk_posz_map[bt_primary_track_id].end(), 0.0) /
        trk_posz_map[bt_primary_track_id].size();
    bt_primary_pdg = trk_pdg_map[bt_primary_track_id];
    bt_primary_edep = trk_edep_map[bt_primary_track_id];
    bt_total_edep = std::accumulate(
        trk_edep_map.begin(), trk_edep_map.end(), 0.0,
        [](const auto& a, const auto& b) { return a + b.second; });
    bt_primary_numelectrons = trk_numelectrons_map[bt_primary_track_id];
    bt_total_numelectrons = std::accumulate(
        trk_numelectrons_map.begin(), trk_numelectrons_map.end(), 0.0,
        [](const auto& a, const auto& b) { return a + b.second; });
    bt_primary_ke = trk_ke_map[bt_primary_track_id];
  }

  float_t bt_primary_x = -99999;
  float_t bt_primary_y = -99999;
  float_t bt_primary_z = -99999;
  uint64_t bt_primary_pdg = 0;
  float_t bt_primary_edep = -99999;
  float_t bt_total_edep = -99999;
  float_t bt_primary_numelectrons = -99999;
  float_t bt_total_numelectrons = -99999;
  float_t bt_primary_ke = -99999;
};
