// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include <CCDB/BasicCCDBManager.h>
#include <cmath>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/HistogramRegistry.h"

//#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"

#include "GFWPowerArray.h"
#include "GFW.h"
#include "GFWCumulant.h"
#include "FlowContainer.h"
#include "FlowPtContainer.h"
#include "GFWConfig.h"
#include "GFWWeights.h"
#include <TProfile.h>
#include <TRandom3.h>
#include <TF1.h>
#include <iostream>



using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct GfwTutorial {

  O2_DEFINE_CONFIGURABLE(cfgCutVertex, float, 10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMin, float, 0.2f, "Minimal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMax, float, 10.0f, "Maximal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMin, float, 0.2f, "Minimal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMax, float, 3.0f, "Maximal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutEta, float, 0.8f, "Eta range for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutChi2prTPCcls, float, 2.5, "Chi2 per TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgUseNch, bool, false, "Use Nch for flow observables")
  O2_DEFINE_CONFIGURABLE(cfgNbootstrap, int, 10, "Number of subsamples")
  //O2_DEFINE_CONFIGURABLE(cfgTPCCut, float, 2, "Cut on PID Nsigma from TPC");

  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for histograms"};
  ConfigurableAxis axisPhi{"axisPhi", {60, 0.0, constants::math::TwoPI}, "phi axis for histograms"};
  ConfigurableAxis axisEta{"axisEta", {40, -1., 1.}, "eta axis for histograms"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.5, 4, 5, 6, 8, 10}, "pt axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90}, "centrality axis for histograms"};

  std::vector<double> ptbinning = {0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.5, 4, 5, 6, 8, 10};
  int ptbins = ptbinning.size() - 1;

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPtMin) && (aod::track::pt < cfgCutPtMax) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && (aod::track::tpcChi2NCl < cfgCutChi2prTPCcls);

  // Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<int64_t> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> url{"ccdb-url", "http://ccdb-test.cern.ch:8080", "url of the ccdb repository"};

  // Define output
  HistogramRegistry registry{"registry"};


  //Define flow container?
  OutputObj<FlowContainer> fFC{FlowContainer("FlowContainer")};


  // define global variables
  GFW* fGFW = new GFW(); 
  TRandom3* fRndm = new TRandom3(0);
  TAxis* fPtAxis;
  std::vector<GFW::CorrConfig> corrconfigs;

  using aodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms>>;
  using aodCollisionsRun3 = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>>;
  using aodTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra>>;

  void init(InitContext const&)
  {
    //AxisSpec ptAxis = {axisPt};
    fPtAxis = new TAxis(ptbins, &ptbinning[0]);

    ccdb->setURL(url.value);
    ccdb->setCaching(true);
    ccdb->setCreatedNotAfter(nolaterthan.value);

    // Add some output objects to the histogram registry
    registry.add("hPhi", "", {HistType::kTH1D, {axisPhi}});
    registry.add("hEta", "", {HistType::kTH1D, {axisEta}});
    registry.add("hVtxZ", "", {HistType::kTH1D, {axisVertex}});
    registry.add("hMult", "", {HistType::kTH1D, {{3000, 0.5, 3000.5}}});
    registry.add("hCent", "", {HistType::kTH1D, {{90, 0, 90}}});
     
    registry.add("c22", "", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("ptpoi", "", {HistType::kTProfile, {axisPt}});
    

//pT flow
    fGFW->AddRegion("refFull", -0.8, 0.8, 1, 1);
    fGFW->AddRegion("refN", -0.8, -0.4, 1, 1);
    fGFW->AddRegion("refP", 0.4, 0.8, 1, 1);
    fGFW->AddRegion("poiN", -0.8, -0.4, 1+fPtAxis->GetNbins(), 2);
    fGFW->AddRegion("olN", -0.8, -0.4, 1, 4);

    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN {2} refP {-2}", "ChGap22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN {3} refP {-3}", "ChGap32", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN {4} refP {-4}", "ChGap42", kFALSE));


    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN {2 2} refP {-2 -2}", "ChGap24", kFALSE));

    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiN refN | olN {2} refP {-2}", "ChGap22", kTRUE));

   
    fGFW->CreateRegions();

    TObjArray* oba = new TObjArray();
    oba->Add(new TNamed("ChGap22", "ChGap22"));
    for (Int_t i=0;i<fPtAxis->GetNbins();i++)
      oba->Add(new TNamed(Form("ChGap22_pt_%i",i+1), "ChGap22_pTDiff"));
    fFC->SetName("FlowContainer");
    fFC->SetXAxis(fPtAxis);
    fFC->Initialize(oba,axisMultiplicity,cfgNbootstrap);
    delete oba;
  }

  template <char... chars>
  void FillProfile(const GFW::CorrConfig& corrconf, const double& cent, const int& rndm) 
  {
    double dnx, val;
    dnx = fGFW->Calculate(corrconf, 0, kTRUE).real();
    if (dnx == 0)
      return;

    if (!corrconf.pTDif) {
      val = fGFW->Calculate(corrconf, 0, kFALSE).real() / dnx;
      if (TMath::Abs(val) < 1)
        fFC->FillProfile(corrconf.Head.c_str(),cent,val,dnx,rndm);
      return;
    }

    return;
  }

  template <char... chars>
  void FillPtProfile(const GFW::CorrConfig& corrconf, const double& cent, const int& rndm) //, const int& ptBin
  {
    double dnx, val;
    for(int i=1;i<=fPtAxis->GetNbins();i++){

      //LOGF(info, "pT bin: %d", i);
    
    
      dnx = fGFW->Calculate(corrconf, i-1, kTRUE).real();
      if (dnx == 0)
        continue;

      
      val = fGFW->Calculate(corrconf, i-1, kFALSE).real() / dnx;
      if (TMath::Abs(val) < 1)
        fFC->FillProfile(Form("%s_pt_%i", corrconf.Head.c_str(), i),cent,val,dnx,rndm);
    }
    return;
  }
  //change for run 2
  void process(aodCollisionsRun3::iterator const& collision, aod::BCsWithTimestamps const&, aodTracks const& tracks)
  {
    int Ntot = tracks.size();
    if (Ntot < 1){
      //LOGF(info,"INFO: too few tracks!");
      return;
    }
    if (!collision.sel8()){ //change to sel7 for run2
      //LOGF(info,"INFO: sel7 unhappy here");
      return;
    }
    //LOGF(info,"INFO: is sel7 ever happy?");
    float vtxz = collision.posZ();
    registry.fill(HIST("hVtxZ"), vtxz);
    registry.fill(HIST("hMult"), Ntot);
    //registry.fill(HIST("hCent"), collision.centRun2V0M());
    registry.fill(HIST("hCent"), collision.centFT0C());
    fGFW->Clear();
    //const auto cent = collision.centRun2V0M();
    const auto cent = collision.centFT0C();
    //LOGF(info,"INFO: centrality is %f",cent);
    float weff = 1, wacc = 1;
    for (auto& track : tracks) {
      registry.fill(HIST("hPhi"), track.phi());
      registry.fill(HIST("hEta"), track.eta());

      //fGFW->Fill(track.eta(), 1, track.phi(), wacc * weff, 1);
      // LOGF(info, "!");      
       bool WithinPtPOI = (cfgCutPtPOIMin < track.pt()) && (track.pt() < cfgCutPtPOIMax); // within POI pT range
       bool WithinPtRef = (cfgCutPtMin < track.pt()) && (track.pt() < cfgCutPtMax); // within RF pT range
       if (WithinPtRef) fGFW->Fill(track.eta(), fPtAxis->FindBin(track.pt()) - 1, track.phi(), weff * wacc, 1);
       if (WithinPtPOI) fGFW->Fill(track.eta(), fPtAxis->FindBin(track.pt()) - 1, track.phi(), weff * wacc, 2);
       if (WithinPtPOI && WithinPtRef) fGFW->Fill(track.eta(), fPtAxis->FindBin(track.pt()) - 1, track.phi(), weff * wacc, 4);
      // LOGF(info, "!!"); 
      
    }

    // Filling Flow Container
    float subsamps = fRndm->Rndm();
    FillProfile(corrconfigs.at(0), cent, subsamps);
    FillProfile(corrconfigs.at(1), cent, subsamps);
    FillProfile(corrconfigs.at(2), cent, subsamps);
    FillProfile(corrconfigs.at(3), cent, subsamps);
    FillPtProfile(corrconfigs.at(4), cent, subsamps);
    
    

  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<GfwTutorial>(cfgc)};
}
