/*
 * AliEmcalAnalysisFactory.cxx
 *
 *  Created on: Feb 23, 2016
 *      Author: markus
 */
#include <functional>
#include <vector>

#include "AliAODTrack.h"
#include "AliESDtrackCuts.h"
#include "AliEmcalTrackSelection.h"
#include "AliEmcalTrackSelectionESD.h"
#include "AliEmcalTrackSelectionAOD.h"
#include "AliEmcalTriggerOfflineSelection.h"
#include "AliEMCalTriggerExtraCuts.h"

#include "AliEmcalAnalysisFactory.h"

ClassImp(EMCalTriggerPtAnalysis::AliEmcalAnalysisFactory)

namespace EMCalTriggerPtAnalysis {

AliEmcalTrackSelection *AliEmcalAnalysisFactory::TrackCutsFactory(TString cut, Bool_t aod){
  AliEmcalTrackSelection *result = NULL;
  if(!aod){
    std::vector<AliVCuts *> trackcuts;
    if(cut.Contains("standard")){
      AliESDtrackCuts *esdcuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(true, 1);
      esdcuts->DefineHistograms(kRed);
      esdcuts->SetName("Standard Track cuts");
      esdcuts->SetMinNCrossedRowsTPC(120);
      esdcuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
      trackcuts.push_back(esdcuts);
    }
    if(cut.Contains("hybrid")){
      AliESDtrackCuts *esdcuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);
      esdcuts->SetName("Global Hybrid tracks, loose DCA");
      esdcuts->SetMaxDCAToVertexXY(2.4);
      esdcuts->SetMaxDCAToVertexZ(3.2);
      esdcuts->SetDCAToVertex2D(kTRUE);
      esdcuts->SetMaxChi2TPCConstrainedGlobal(36);
      esdcuts->SetMaxFractionSharedTPCClusters(0.4);
      trackcuts.push_back(esdcuts);
    }
    if(cut.Contains("geo")){
      AliEMCalTriggerExtraCuts *geocuts = new AliEMCalTriggerExtraCuts();
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

AliEmcalTriggerOfflineSelection *AliEmcalAnalysisFactory::TriggerSelectionFactory(Double_t el0, Double_t eg1, Double_t eg2, Double_t ej1, Double_t ej2){
  AliEmcalTriggerOfflineSelection *result = new AliEmcalTriggerOfflineSelection;
  result->SetOfflineEnergyThreshold(AliEmcalTriggerOfflineSelection::kTrgEL0, el0);
  result->SetOfflineEnergyThreshold(AliEmcalTriggerOfflineSelection::kTrgEG1, eg1);
  result->SetOfflineEnergyThreshold(AliEmcalTriggerOfflineSelection::kTrgEG2, eg1);
  result->SetOfflineEnergyThreshold(AliEmcalTriggerOfflineSelection::kTrgEJ1, ej1);
  result->SetOfflineEnergyThreshold(AliEmcalTriggerOfflineSelection::kTrgEJ2, ej2);
  return result;
}

} /* namespace EMCalTriggerPtAnalysis */
