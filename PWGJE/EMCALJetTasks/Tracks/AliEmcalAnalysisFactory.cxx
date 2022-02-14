/************************************************************************************
 * Copyright (C) 2017, Copyright Holders of the ALICE Collaboration                 *
 * All rights reserved.                                                             *
 *                                                                                  *
 * Redistribution and use in source and binary forms, with or without               *
 * modification, are permitted provided that the following conditions are met:      *
 *     * Redistributions of source code must retain the above copyright             *
 *       notice, this list of conditions and the following disclaimer.              *
 *     * Redistributions in binary form must reproduce the above copyright          *
 *       notice, this list of conditions and the following disclaimer in the        *
 *       documentation and/or other materials provided with the distribution.       *
 *     * Neither the name of the <organization> nor the                             *
 *       names of its contributors may be used to endorse or promote products       *
 *       derived from this software without specific prior written permission.      *
 *                                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY              *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     *
 ************************************************************************************/
#include <cstring>
#include <cfloat>
#include <functional>
#include <iostream>
#include <memory>
#include <vector>

#include <TObjArray.h>
#include <TObjString.h>
#include <TString.h>

#include "AliAODTrack.h"
#include "AliESDtrackCuts.h"
#include "AliEmcalAODFilterBitCuts.h"
#include "AliEmcalESDHybridTrackCuts.h"
#include "AliEmcalTrackSelection.h"
#include "AliEmcalTrackSelectionESD.h"
#include "AliEmcalTrackSelectionAOD.h"
#include "AliEmcalTriggerOfflineSelection.h"
#include "AliEMCalTriggerExtraCuts.h"
#include "AliTrackContainer.h"

#include "AliEmcalAnalysisFactory.h"

using namespace PWGJE::EMCALJetTasks;

AliEmcalTrackSelection *AliEmcalAnalysisFactory::TrackCutsFactory(TString cutstring, Bool_t aod){
  AliEmcalTrackSelection *result = NULL;
  std::unique_ptr<TObjArray> cuts(cutstring.Tokenize(","));
  std::cout << "Creating track cuts for " << (aod ? "AODs" : "ESDs") << std::endl;
  TObjArray trackcuts;
  if(!aod){
    for(auto c : *cuts){
      TString &cut = static_cast<TObjString *>(c)->String();
      if(cut == "standard"){
        auto esdcuts = GenerateDefaultCutsESD();
        esdcuts->SetName("standardRAA");
        esdcuts->SetTitle("Standard Track cuts");
        trackcuts.Add(esdcuts);
      }
      if(cut == "standardcrossedrows"){
        auto esdcuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(true, 1);
        esdcuts->DefineHistograms(kRed);
        esdcuts->SetName("standardRAA");
        esdcuts->SetTitle("Standard Track cuts");
        esdcuts->SetMinNCrossedRowsTPC(120);
        esdcuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
        trackcuts.Add(esdcuts);
      }
      if(cut == "hybrid"){
        std::cout << "Configuring standard hybrid track cuts" << std::endl;
        auto esdcuts = GenerateLooseDCACutsESD();
        esdcuts->SetTitle("hybridglobal");
        esdcuts->SetTitle("Global Hybrid tracks, loose DCA");
        trackcuts.Add(esdcuts);
      }
      if(cut == "hybrid2010_wNoRefit"){
        std::cout << "Configuring hybrid track cuts 2010, including non-ITSrefit tracks" << std::endl;
        auto hybridcuts = new PWG::EMCAL::AliEmcalESDHybridTrackCuts("2010_wNoRefit", PWG::EMCAL::AliEmcalESDHybridTrackCuts::kDef2010);
        hybridcuts->SetUseNoITSrefitTracks(kTRUE);
        trackcuts.Add(hybridcuts);
      }
      if(cut == "hybrid2010_woNoRefit"){
        std::cout << "Configuring hybrid track cuts 2010, excluding non-ITSrefit tracks" << std::endl;
        auto hybridcuts = new PWG::EMCAL::AliEmcalESDHybridTrackCuts("2010_woNoRefit", PWG::EMCAL::AliEmcalESDHybridTrackCuts::kDef2010);
        hybridcuts->SetUseNoITSrefitTracks(kFALSE);
        trackcuts.Add(hybridcuts); 
     }
      if(cut == "hybrid2011_wNoRefit"){
        std::cout << "Configuring hybrid track cuts 2011, including non-ITSrefit tracks" << std::endl;
        auto hybridcuts = new PWG::EMCAL::AliEmcalESDHybridTrackCuts("2011_wNoRefit", PWG::EMCAL::AliEmcalESDHybridTrackCuts::kDef2011);
        hybridcuts->SetUseNoITSrefitTracks(kTRUE);
        trackcuts.Add(hybridcuts);
      }
      if(cut == "hybrid2011_woNoRefit"){
        std::cout << "Configuring hybrid track cuts 2011, excluding non-ITSrefit tracks" << std::endl;
        auto hybridcuts = new PWG::EMCAL::AliEmcalESDHybridTrackCuts("2011_woNoRefit", PWG::EMCAL::AliEmcalESDHybridTrackCuts::kDef2011);
        hybridcuts->SetUseNoITSrefitTracks(kFALSE);
        trackcuts.Add(hybridcuts);
      }

      ////////////////////////////////////////////////////////////////
      /// Systematics cuts                                         ///
      ///   These cuts are based on the default cuts varying       ///
      ///   each time one cut separately                           ///
      ////////////////////////////////////////////////////////////////
      if(cut.Contains("VarITSchi2")){
        double cutvalue = ValueDecoder(cut, "VarITSchi2");
        auto esdcuts = GenerateDefaultCutsESD();
        esdcuts->SetName(TString::Format("VarITSchi2Cut%04d", static_cast<int>(cutvalue * 10.)));
        esdcuts->SetTitle(TString::Format("Default cuts - variation ITS chi2 cut at %f", cutvalue));
        esdcuts->SetMaxChi2PerClusterITS(cutvalue);
        trackcuts.Add(esdcuts);
      }
      if(cut.Contains("VarTPCchi2") && !cut.Contains("Constrained")){
        double cutvalue = ValueDecoder(cut, "VarTPCchi2");
        auto esdcuts = GenerateDefaultCutsESD();
        esdcuts->SetName(TString::Format("VarTPCchi2Cut%04d", static_cast<int>(cutvalue * 10.)));
        esdcuts->SetTitle(TString::Format("Default cuts - variation TPC chi2 cut at %f", cutvalue));
        esdcuts->SetMaxChi2PerClusterTPC(cutvalue);
        trackcuts.Add(esdcuts);
      }
      if(cut.Contains("VarTPCchi2Constrained")){
        double cutvalue = ValueDecoder(cut, "VarTPCchi2Constrained"); // in number of sigmas
        auto esdcuts = GenerateDefaultCutsESD();
        esdcuts->SetName(TString::Format("VarTPCchi2ConstrainedCut%04d", static_cast<int>(cutvalue * 10.)));
        esdcuts->SetTitle(TString::Format("Default cuts - variation TPC chi2 constrained cut at %f", cutvalue));
        esdcuts->SetMaxChi2TPCConstrainedGlobal(cutvalue);
        trackcuts.Add(esdcuts);
      }
      if(cut.Contains("VarDCAz")){
        double cutvalue = ValueDecoder(cut, "VarDCAz");
        auto esdcuts = GenerateDefaultCutsESD();
        esdcuts->SetName(TString::Format("VarDCAzCut%04d", static_cast<int>(cutvalue * 10.)));
        esdcuts->SetTitle(TString::Format("Default cuts - variation DCAz cut at %f", cutvalue));
        esdcuts->SetMaxDCAToVertexZ(cutvalue);
        trackcuts.Add(esdcuts);
      }
      if(cut.Contains("VarDCAr")){
        double cutvalue = ValueDecoder(cut, "VarDCAr");
        double p1 = cutvalue * 0.0026, p2 = cutvalue * 0.005;
        auto esdcuts = GenerateDefaultCutsESD();
        esdcuts->SetName(TString::Format("VarDCArCut%04d", static_cast<int>(cutvalue * 10.)));
        esdcuts->SetTitle(TString::Format("Default cuts - variation DCAr cut at %f sigma", cutvalue));
        esdcuts->SetMaxDCAToVertexXYPtDep(TString::Format("%f + %f/pt^1.01", p1, p2));
        trackcuts.Add(esdcuts);
      }
      if(cut.Contains("VarRatioCrossedRowsFindable")){
        double cutvalue = ValueDecoder(cut, "VarRatioCrossedRowsFindable");
        auto esdcuts = GenerateDefaultCutsESD();
        esdcuts->SetName(TString::Format("VarRatioCRFindableCut%04d", static_cast<int>(cutvalue * 10.)));
        esdcuts->SetTitle(TString::Format("Default cuts - variation ratio crossed rows over findable cut at %f", cutvalue));
        esdcuts->SetMinRatioCrossedRowsOverFindableClustersTPC(cutvalue);
        trackcuts.Add(esdcuts);
      }
      if(cut.Contains("VarFractionTPCshared")){
        double cutvalue = ValueDecoder(cut, "VarFractionTPCshared");
        auto esdcuts = GenerateDefaultCutsESD();
        esdcuts->SetName(TString::Format("VarFractionTPCShared%04d", static_cast<int>(cutvalue * 10.)));
        esdcuts->SetTitle(TString::Format("Default cuts - variation fraction TPC shared clusters %f", cutvalue));
        esdcuts->SetMaxFractionSharedTPCClusters(cutvalue);
        trackcuts.Add(esdcuts);
      }
      if(cut.Contains("VarTrackLengthDeadArea")){
        double cutvalue = ValueDecoder(cut, "VarTrackLengthDeadArea");
        AliESDtrackCuts *esdcuts = GenerateDefaultCutsESD();
        esdcuts->SetName(TString::Format("VarTLDeadAreaCut%04d", static_cast<int>(cutvalue * 10.)));
        esdcuts->SetTitle(TString::Format("Default cuts - variation track length dead area cut at %f", cutvalue));
        esdcuts->SetCutGeoNcrNcl(cutvalue, 130., 1.5, 0.0, 0.0);
        trackcuts.Add(esdcuts);
      }
      if(cut.Contains("VarTrackLengthTPCLength")){
        double cutvalue = ValueDecoder(cut, "VarTrackLengthTPCLength");
        auto esdcuts = GenerateDefaultCutsESD();
        esdcuts->SetName(TString::Format("VarTLTPCLengthCut%04d", static_cast<int>(cutvalue * 10.)));
        esdcuts->SetTitle(TString::Format("Default cuts - variation track length TPC length cut at %f", cutvalue));
        esdcuts->SetCutGeoNcrNcl(3., cutvalue, 1.5, 0.0, 0.0);
        trackcuts.Add(esdcuts);
      }
      if(cut.Contains("VarSPDhit")){
        double cutvalue = ValueDecoder(cut, "VarSPDhit");
        auto esdcuts = GenerateDefaultCutsESD();
        if(TMath::Abs(cutvalue) > DBL_EPSILON){
          esdcuts->SetName("VarSPDhitOn");
          esdcuts->SetTitle("Default cuts - variation SPD hit requirement on");
        } else {
          esdcuts->SetName("VarSPDhitOff");
          esdcuts->SetTitle("Default cuts - variation SPD hit requirement off");
          esdcuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kOff);
        }
        trackcuts.Add(esdcuts);
      }
      ////////////////////////////////////////////////////////////////
      /// Test cuts - for tracking studies                         ///
      ///   Cuts are based on loose cuts. Those which are entangled///
      ///   are loosened furthermore                               ///
      ////////////////////////////////////////////////////////////////
      if(cut.Contains("TestITSchi2")){
        double itscut = ValueDecoder(cut, "TestITSchi2");
        std::cout << "Using ITS chi2 cut variation: " << itscut << std::endl;

        auto esdcuts = GenerateLooseDCACutsESD();
        esdcuts->SetName(Form("TestITSchi2%d", int(itscut*10.)));
        esdcuts->SetTitle(Form("Loose track cuts, ITS chi2 var %.1f", itscut));

        // Do the variation
        esdcuts->SetMaxChi2PerClusterITS(itscut);
        // Set cut on the TPC global constrained chi2 to very loose
        esdcuts->SetMaxChi2TPCConstrainedGlobal(100);
        trackcuts.Add(esdcuts);
      }
      if(cut.Contains("TestTPCchi2")){
        double tpccut = ValueDecoder(cut,"TestTPCchi2");
        std::cout << "Using TPC chi2 cut variation: " << tpccut << std::endl;

        auto esdcuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(false, 1);
        esdcuts->SetName(Form("VarTPCchi2%d", int(tpccut * 10.)));
        esdcuts->SetTitle(Form("Loose track cuts, TPC chi2 var %.1f", tpccut));

        // Do the variation
        esdcuts->SetMaxChi2PerClusterTPC(tpccut);
        trackcuts.Add(esdcuts);
      }
      if(cut.Contains("TestTPCchi2Constrained")){
        double tpcconstrainedcut = ValueDecoder(cut, "TestTPCchi2Constrained");
        std::cout << "Using TPC chi2 constrained cut variation: " << tpcconstrainedcut << std::endl;

        auto esdcuts = GenerateLooseDCACutsESD();
        esdcuts->SetName(Form("VarTPCchi2constrained%d", int(tpcconstrainedcut * 10.)));
        esdcuts->SetTitle(Form("Loose track cuts, TPC constrained chi2 variation %f", tpcconstrainedcut));

        // Do the variation
        esdcuts->SetMaxChi2TPCConstrainedGlobal(tpcconstrainedcut);
        // Set ITS chi2 cut to very loose
        esdcuts->SetMaxChi2PerClusterITS(100);
        trackcuts.Add(esdcuts);
      }
      if(cut.Contains("geo")){
        auto geocuts = new AliEMCalTriggerExtraCuts();
        geocuts->SetName("geocuts");
        geocuts->SetTitle("TPC track length cut");
        geocuts->SetMinTPCTrackLengthCut();
        trackcuts.Add(geocuts);
      }
    }
    result = new AliEmcalTrackSelectionESD;
    result->AddTrackCuts(&trackcuts);
  } else {
    AliEmcalTrackSelectionAOD *aodsel = new AliEmcalTrackSelectionAOD;
    result = aodsel;
    // C++11 Lambda: Do not create multiple extra cut objects in case of AODs. If extra cut object does already exist -
    // specify new cut in the same object.
    auto FindTrackCuts = [] (const TObjArray &cuts) -> AliEMCalTriggerExtraCuts * {
      AliEMCalTriggerExtraCuts *found = nullptr;
      for(auto cutiter : cuts){
        if(cutiter->IsA() == AliEMCalTriggerExtraCuts::Class()){
          found = static_cast<AliEMCalTriggerExtraCuts *>(cutiter);
          break;
        }
      }
      return found;
    };
    for(auto c : *cuts){
      TString &cut = static_cast<TObjString *>(c)->String();
      if(cut == "standard"){
        auto filterbitcuts = new PWG::EMCAL::AliEmcalAODFilterBitCuts("globalcuts", "Global track cuts");
        filterbitcuts->SetFilterBits(AliAODTrack::kTrkGlobal);
        std::cout << "Adding standard global track cuts" << std::endl;
        trackcuts.Add(filterbitcuts);

        AliEMCalTriggerExtraCuts *crossedrowcut = FindTrackCuts(trackcuts);
        if(!crossedrowcut){
          crossedrowcut = new AliEMCalTriggerExtraCuts;
          trackcuts.Add(crossedrowcut);
        }
        crossedrowcut->SetMinTPCCrossedRows(120);
      }
      if(cut == "hybrid"){
        aodsel->GenerateTrackCuts(AliEmcalTrackSelection::kHybridTracks, AliTrackContainer::GetDefTrackCutsPeriod());
      }
      if(cut == "hybrid2010_wNoRefit"){
        aodsel->GenerateTrackCuts(AliEmcalTrackSelection::kHybridTracks2010wNoRefit, AliTrackContainer::GetDefTrackCutsPeriod());
      }
      if(cut == "hybrid2010_woNoRefit"){
        aodsel->GenerateTrackCuts(AliEmcalTrackSelection::kHybridTracks2010woNoRefit, AliTrackContainer::GetDefTrackCutsPeriod());
      }
      if(cut == "hybrid2011_wNoRefit"){
        aodsel->GenerateTrackCuts(AliEmcalTrackSelection::kHybridTracks2011wNoRefit, AliTrackContainer::GetDefTrackCutsPeriod());
      }
      if(cut == "hybrid2011_woNoRefit"){
        aodsel->GenerateTrackCuts(AliEmcalTrackSelection::kHybridTracks2011woNoRefit, AliTrackContainer::GetDefTrackCutsPeriod());
      }
      if(cut == "geo"){
        AliEMCalTriggerExtraCuts *geocuts = FindTrackCuts(trackcuts);
        if(!geocuts){
          geocuts = new AliEMCalTriggerExtraCuts;
          trackcuts.Add(geocuts);
        }
        geocuts->SetMinTPCTrackLengthCut();
      }
    }
    result->AddTrackCuts(&trackcuts);
  }

  return result;
}

AliEmcalTriggerOfflineSelection *AliEmcalAnalysisFactory::TriggerSelectionFactory(Double_t el0, Double_t eg1, Double_t eg2, Double_t ej1, Double_t ej2, AliEmcalTriggerOfflineSelection::EmcalEnergyDefinition_t endef){
  AliEmcalTriggerOfflineSelection *result = new AliEmcalTriggerOfflineSelection;
  result->SetOfflineEnergyThreshold(AliEmcalTriggerOfflineSelection::kTrgEL0, el0);
  result->SetOfflineEnergyThreshold(AliEmcalTriggerOfflineSelection::kTrgEG1, eg1);
  result->SetOfflineEnergyThreshold(AliEmcalTriggerOfflineSelection::kTrgEG2, eg2);
  result->SetOfflineEnergyThreshold(AliEmcalTriggerOfflineSelection::kTrgEJ1, ej1);
  result->SetOfflineEnergyThreshold(AliEmcalTriggerOfflineSelection::kTrgEJ2, ej2);
  result->SetEnergyDefinition(endef);
  return result;
}

TString AliEmcalAnalysisFactory::ClusterContainerNameFactory(Bool_t isAOD){
  return isAOD ? "caloClusters" : "CaloClusters";
}

TString AliEmcalAnalysisFactory::TrackContainerNameFactory(Bool_t isAOD){
  return isAOD ? "tracks" : "Tracks";
}

AliESDtrackCuts *AliEmcalAnalysisFactory::GenerateDefaultCutsESD() {
  AliESDtrackCuts *esdcuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(true, 1);
  esdcuts->DefineHistograms(kRed);
  esdcuts->SetCutGeoNcrNcl(3., 130., 1.5, 0.0, 0.0);
  esdcuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
  esdcuts->SetMaxFractionSharedTPCClusters(0.4);
  // esdcuts->SetMinNCrossedRowsTPC(0.); // Replaced by geometrical cut
  return esdcuts;

}

AliESDtrackCuts *AliEmcalAnalysisFactory::GenerateLooseDCACutsESD(){
  AliESDtrackCuts *esdcuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);
  esdcuts->SetMaxDCAToVertexXY(2.4);
  esdcuts->SetMaxDCAToVertexZ(3.2);
  esdcuts->SetDCAToVertex2D(kTRUE);
  esdcuts->SetMaxChi2TPCConstrainedGlobal(36);
  esdcuts->SetMaxFractionSharedTPCClusters(0.4);
  return esdcuts;
}

double AliEmcalAnalysisFactory::ValueDecoder(const char *cutstring, const char *tag) {
  TString cuttstring(cutstring);
  Int_t position(cuttstring.Index(tag) + strlen(tag));
  TString valuestring = cuttstring(position, 4);
  Int_t value = valuestring.Atoi();
  return static_cast<double>(value) / 10.;
}
