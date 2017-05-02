/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
#include <cstring>
#include <cfloat>
#include <functional>
#include <iostream>
#include <vector>

#include "AliAODTrack.h"
#include "AliESDtrackCuts.h"
#include "AliEmcalTrackSelection.h"
#include "AliEmcalTrackSelectionESD.h"
#include "AliEmcalTrackSelectionAOD.h"
#include "AliEmcalTriggerOfflineSelection.h"
#include "AliEMCalTriggerExtraCuts.h"

#include "AliEmcalAnalysisFactory.h"

/// \cond CLASSIMP
ClassImp(EMCalTriggerPtAnalysis::AliEmcalAnalysisFactory)
/// \endcond
///
namespace EMCalTriggerPtAnalysis {

AliEmcalTrackSelection *AliEmcalAnalysisFactory::TrackCutsFactory(TString cut, Bool_t aod){
  AliEmcalTrackSelection *result = NULL;
  if(!aod){
    std::vector<AliVCuts *> trackcuts;
    if(cut.Contains("standard")){
      AliESDtrackCuts *esdcuts = GenerateDefaultCutsESD();
      esdcuts->SetName("standardRAA");
      esdcuts->SetTitle("Standard Track cuts");
      trackcuts.push_back(esdcuts);
    }
    if(cut.Contains("standardcrossedrows")){
      AliESDtrackCuts *esdcuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(true, 1);
      esdcuts->DefineHistograms(kRed);
      esdcuts->SetName("standardRAA");
      esdcuts->SetTitle("Standard Track cuts");
      esdcuts->SetMinNCrossedRowsTPC(120);
      esdcuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
      trackcuts.push_back(esdcuts);
    }
    if(cut.Contains("hybrid")){
      AliESDtrackCuts *esdcuts = GenerateLooseDCACutsESD();
      esdcuts->SetTitle("hybridglobal");
      esdcuts->SetTitle("Global Hybrid tracks, loose DCA");
      trackcuts.push_back(esdcuts);
    }
    ////////////////////////////////////////////////////////////////
    /// Systematics cuts                                         ///
    ///   These cuts are based on the default cuts varying       ///
    ///   each time one cut separately                           ///
    ////////////////////////////////////////////////////////////////
    if(cut.Contains("VarITSchi2")){
      double cutvalue = ValueDecoder(cut, "VarITSchi2");
      AliESDtrackCuts *esdcuts = GenerateDefaultCutsESD();
      esdcuts->SetName(TString::Format("VarITSchi2Cut%04d", static_cast<int>(cutvalue * 10.)));
      esdcuts->SetTitle(TString::Format("Default cuts - variation ITS chi2 cut at %f", cutvalue));
      esdcuts->SetMaxChi2PerClusterITS(cutvalue);
      trackcuts.push_back(esdcuts);
    }
    if(cut.Contains("VarTPCchi2") && !cut.Contains("Constrained")){
      double cutvalue = ValueDecoder(cut, "VarTPCchi2");
      AliESDtrackCuts *esdcuts = GenerateDefaultCutsESD();
      esdcuts->SetName(TString::Format("VarTPCchi2Cut%04d", static_cast<int>(cutvalue * 10.)));
      esdcuts->SetTitle(TString::Format("Default cuts - variation TPC chi2 cut at %f", cutvalue));
      esdcuts->SetMaxChi2PerClusterTPC(cutvalue);
      trackcuts.push_back(esdcuts);
    }
    if(cut.Contains("VarTPCchi2Constrained")){
      double cutvalue = ValueDecoder(cut, "VarTPCchi2Constrained"); // in number of sigmas
      AliESDtrackCuts *esdcuts = GenerateDefaultCutsESD();
      esdcuts->SetName(TString::Format("VarTPCchi2ConstrainedCut%04d", static_cast<int>(cutvalue * 10.)));
      esdcuts->SetTitle(TString::Format("Default cuts - variation TPC chi2 constrained cut at %f", cutvalue));
      esdcuts->SetMaxChi2TPCConstrainedGlobal(cutvalue);
      trackcuts.push_back(esdcuts);
    }
    if(cut.Contains("VarDCAz")){
      double cutvalue = ValueDecoder(cut, "VarDCAz");
      AliESDtrackCuts *esdcuts = GenerateDefaultCutsESD();
      esdcuts->SetName(TString::Format("VarDCAzCut%04d", static_cast<int>(cutvalue * 10.)));
      esdcuts->SetTitle(TString::Format("Default cuts - variation DCAz cut at %f", cutvalue));
      esdcuts->SetMaxDCAToVertexZ(cutvalue);
      trackcuts.push_back(esdcuts);
    }
    if(cut.Contains("VarDCAr")){
      double cutvalue = ValueDecoder(cut, "VarDCAr");
      double p1 = cutvalue * 0.0026, p2 = cutvalue * 0.005;
      AliESDtrackCuts *esdcuts = GenerateDefaultCutsESD();
      esdcuts->SetName(TString::Format("VarDCArCut%04d", static_cast<int>(cutvalue * 10.)));
      esdcuts->SetTitle(TString::Format("Default cuts - variation DCAr cut at %f sigma", cutvalue));
      esdcuts->SetMaxDCAToVertexXYPtDep(TString::Format("%f + %f/pt^1.01", p1, p2));
      trackcuts.push_back(esdcuts);
    }
    if(cut.Contains("VarRatioCrossedRowsFindable")){
      double cutvalue = ValueDecoder(cut, "VarRatioCrossedRowsFindable");
      AliESDtrackCuts *esdcuts = GenerateDefaultCutsESD();
      esdcuts->SetName(TString::Format("VarRatioCRFindableCut%04d", static_cast<int>(cutvalue * 10.)));
      esdcuts->SetTitle(TString::Format("Default cuts - variation ratio crossed rows over findable cut at %f", cutvalue));
      esdcuts->SetMinRatioCrossedRowsOverFindableClustersTPC(cutvalue);
      trackcuts.push_back(esdcuts);
    }
    if(cut.Contains("VarFractionTPCshared")){
      double cutvalue = ValueDecoder(cut, "VarFractionTPCshared");
      AliESDtrackCuts *esdcuts = GenerateDefaultCutsESD();
      esdcuts->SetName(TString::Format("VarFractionTPCShared%04d", static_cast<int>(cutvalue * 10.)));
      esdcuts->SetTitle(TString::Format("Default cuts - variation fraction TPC shared clusters %f", cutvalue));
      esdcuts->SetMaxFractionSharedTPCClusters(cutvalue);
      trackcuts.push_back(esdcuts);
    }
    if(cut.Contains("VarTrackLengthDeadArea")){
      double cutvalue = ValueDecoder(cut, "VarTrackLengthDeadArea");
      AliESDtrackCuts *esdcuts = GenerateDefaultCutsESD();
      esdcuts->SetName(TString::Format("VarTLDeadAreaCut%04d", static_cast<int>(cutvalue * 10.)));
      esdcuts->SetTitle(TString::Format("Default cuts - variation track length dead area cut at %f", cutvalue));
      esdcuts->SetCutGeoNcrNcl(cutvalue, 130., 1.5, 0.0, 0.0);
      trackcuts.push_back(esdcuts);
    }
    if(cut.Contains("VarTrackLengthTPCLength")){
      double cutvalue = ValueDecoder(cut, "VarTrackLengthTPCLength");
      AliESDtrackCuts *esdcuts = GenerateDefaultCutsESD();
      esdcuts->SetName(TString::Format("VarTLTPCLengthCut%04d", static_cast<int>(cutvalue * 10.)));
      esdcuts->SetTitle(TString::Format("Default cuts - variation track length TPC length cut at %f", cutvalue));
      esdcuts->SetCutGeoNcrNcl(3., cutvalue, 1.5, 0.0, 0.0);
      trackcuts.push_back(esdcuts);
    }
    if(cut.Contains("VarSPDhit")){
      double cutvalue = ValueDecoder(cut, "VarSPDhit");
      AliESDtrackCuts *esdcuts = GenerateDefaultCutsESD();
      if(TMath::Abs(cutvalue) > DBL_EPSILON){
        esdcuts->SetName("VarSPDhitOn");
        esdcuts->SetTitle("Default cuts - variation SPD hit requirement on");
      } else {
        esdcuts->SetName("VarSPDhitOff");
        esdcuts->SetTitle("Default cuts - variation SPD hit requirement off");
        esdcuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kOff);
      }
      trackcuts.push_back(esdcuts);
    }
    ////////////////////////////////////////////////////////////////
    /// Test cuts - for tracking studies                         ///
    ///   Cuts are based on loose cuts. Those which are entangled///
    ///   are loosened furthermore                               ///
    ////////////////////////////////////////////////////////////////
    if(cut.Contains("TestITSchi2")){
      double itscut = ValueDecoder(cut, "TestITSchi2");
      std::cout << "Using ITS chi2 cut variation: " << itscut << std::endl;

      AliESDtrackCuts *esdcuts = GenerateLooseDCACutsESD();
      esdcuts->SetName(Form("TestITSchi2%d", int(itscut*10.)));
      esdcuts->SetTitle(Form("Loose track cuts, ITS chi2 var %.1f", itscut));

      // Do the variation
      esdcuts->SetMaxChi2PerClusterITS(itscut);
      // Set cut on the TPC global constrained chi2 to very loose
      esdcuts->SetMaxChi2TPCConstrainedGlobal(100);
      trackcuts.push_back(esdcuts);
    }
    if(cut.Contains("TestTPCchi2")){
      double tpccut = ValueDecoder(cut,"TestTPCchi2");
      std::cout << "Using TPC chi2 cut variation: " << tpccut << std::endl;

      AliESDtrackCuts *esdcuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(false, 1);
      esdcuts->SetName(Form("VarTPCchi2%d", int(tpccut * 10.)));
      esdcuts->SetTitle(Form("Loose track cuts, TPC chi2 var %.1f", tpccut));

      // Do the variation
      esdcuts->SetMaxChi2PerClusterTPC(tpccut);
      trackcuts.push_back(esdcuts);
    }
    if(cut.Contains("TestTPCchi2Constrained")){
      double tpcconstrainedcut = ValueDecoder(cut, "TestTPCchi2Constrained");
      std::cout << "Using TPC chi2 constrained cut variation: " << tpcconstrainedcut << std::endl;

      AliESDtrackCuts *esdcuts = GenerateLooseDCACutsESD();
      esdcuts->SetName(Form("VarTPCchi2constrained%d", int(tpcconstrainedcut * 10.)));
      esdcuts->SetTitle(Form("Loose track cuts, TPC constrained chi2 variation %f", tpcconstrainedcut));

      // Do the variation
      esdcuts->SetMaxChi2TPCConstrainedGlobal(tpcconstrainedcut);
      // Set ITS chi2 cut to very loose
      esdcuts->SetMaxChi2PerClusterITS(100);
      trackcuts.push_back(esdcuts);
    }
    if(cut.Contains("geo")){
      AliEMCalTriggerExtraCuts *geocuts = new AliEMCalTriggerExtraCuts();
      geocuts->SetName("geocuts");
      geocuts->SetTitle("TPC track length cut");
      geocuts->SetMinTPCTrackLengthCut();
      trackcuts.push_back(geocuts);
    }
    result = new AliEmcalTrackSelectionESD;
    for(std::vector<AliVCuts *>::iterator it = trackcuts.begin(); it != trackcuts.end(); ++it)
      result->AddTrackCuts(*it);
  } else {
    AliEmcalTrackSelectionAOD *aodsel = new AliEmcalTrackSelectionAOD;
    result = aodsel;
    std::vector<AliVCuts *> trackcuts;
    // C++11 Lambda: Do not create multiple extra cut objects in case of AODs. If extra cut object does already exist -
    // specify new cut in the same object.
    std::function<AliEMCalTriggerExtraCuts *(const std::vector<AliVCuts *> &)> FindTrackCuts = [] (const std::vector<AliVCuts *> &cuts) -> AliEMCalTriggerExtraCuts * {
      AliEMCalTriggerExtraCuts *found = nullptr;
      for(std::vector<AliVCuts *>::const_iterator cutiter = cuts.begin(); cutiter != cuts.end(); ++cutiter){
        if((*cutiter)->IsA() == AliEMCalTriggerExtraCuts::Class()){
          found = static_cast<AliEMCalTriggerExtraCuts *>(*cutiter);
          break;
        }
      }
      return found;
    };
    if(cut.Contains("standard")){
      aodsel->AddFilterBit(AliAODTrack::kTrkGlobal);
      AliEMCalTriggerExtraCuts *extracuts = FindTrackCuts(trackcuts);
      if(!extracuts){
        extracuts = new AliEMCalTriggerExtraCuts;
        trackcuts.push_back(extracuts);
      }
      extracuts->SetMinTPCCrossedRows(120);
    }
    if(cut.Contains("hybrid")){
      aodsel->AddFilterBit(256);
      aodsel->AddFilterBit(512);
    }
    if(cut.Contains("geo")){
      AliEMCalTriggerExtraCuts *extracuts = FindTrackCuts(trackcuts);
      if(!extracuts){
        extracuts = new AliEMCalTriggerExtraCuts;
        trackcuts.push_back(extracuts);
      }
      extracuts->SetMinTPCTrackLengthCut();
    }
    for(std::vector<AliVCuts *>::iterator it = trackcuts.begin(); it != trackcuts.end(); ++it)
      result->AddTrackCuts(*it);
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

AliESDtrackCuts *AliEmcalAnalysisFactory::GenerateDefaultCutsESD() {
  AliESDtrackCuts *esdcuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(true, 1);
  esdcuts->DefineHistograms(kRed);
  esdcuts->SetCutGeoNcrNcl(3., 130., 1.5, 0.0, 0.0);
  esdcuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
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

} /* namespace EMCalTriggerPtAnalysis */
