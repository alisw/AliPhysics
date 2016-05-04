/*
 * AliEmcalAnalysisFactory.cxx
 *
 *  Created on: Feb 23, 2016
 *      Author: markus
 */

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
    AliESDtrackCuts *esdcuts = NULL;
    if(cut == "standard"){
      esdcuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(true, 1);
      esdcuts->DefineHistograms(kRed);
      esdcuts->SetName("Standard Track cuts");
      esdcuts->SetMinNCrossedRowsTPC(120);
      esdcuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
    } else if(cut == "hybrid"){
      esdcuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);
      esdcuts->SetName("Global Hybrid tracks, loose DCA");
      esdcuts->SetMaxDCAToVertexXY(2.4);
      esdcuts->SetMaxDCAToVertexZ(3.2);
      esdcuts->SetDCAToVertex2D(kTRUE);
      esdcuts->SetMaxChi2TPCConstrainedGlobal(36);
      esdcuts->SetMaxFractionSharedTPCClusters(0.4);
    }
    result = new AliEmcalTrackSelectionESD;
    result->AddTrackCuts(esdcuts);
  } else {
    AliEmcalTrackSelectionAOD *aodsel = new AliEmcalTrackSelectionAOD;
    result = aodsel;
    if(cut == "standard"){
      aodsel->AddFilterBit(AliAODTrack::kTrkGlobal);
      AliEMCalTriggerExtraCuts *extracuts = new AliEMCalTriggerExtraCuts;
      extracuts->SetMinTPCCrossedRows(120);
      aodsel->AddTrackCuts(extracuts);
    } else if(cut == "hybrid"){
      aodsel->AddFilterBit(256);
      aodsel->AddFilterBit(512);
    }
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
