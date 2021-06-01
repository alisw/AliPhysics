/**
 * \file ReducedEventCutVariation.C
 * \brief Cut variation for the track selection in the reduced event creator task
 *
 * \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * \date Apr 16, 2015
 */
#if !defined __CINT__ || defined __MAKECINT__
#include <TString.h>
#include "AliESDtrackCuts.h"
#include "AliEmcalTrackSelection.h"
#include "AliEmcalTrackSelectionAOD.h"
#include "AliEmcalTrackSelectionESD.h"
#include "AliEMCalTriggerExtraCuts.h"
#include "AliReducedHighPtEventCreator.h"
#endif

AliEmcalTrackSelection *CreateHybridTrackCuts(bool isAOD);
AliEmcalTrackSelection *CreateDefaultTrackCuts(bool isAOD);
AliEmcalTrackSelection *CreateDefaultTrackCutsChangeCrossRows(bool isAOD, int minCrossedRows);
AliEmcalTrackSelection *CreateDefaultTrackCutsChangeFiducialAreaCut(bool isAOD );
AliEmcalTrackSelection *CreateDefaultTrackCutsChangeDCAZ(bool isAOD, Double_t dcaZ);
AliEmcalTrackSelection *CreateDefaultTrackCutsChangeDCAR(bool isAOD, TString model);
AliEmcalTrackSelection *CreateDefaultTrackCutsChangeGoldenCut(bool isAOD, Float_t chi2cut);
AliEmcalTrackSelection *CreateDefaultTrackCutsChangeChi2ITS(bool isAOD, Float_t chi2cut);
AliEmcalTrackSelection *CreateDefaultTrackCutsChangeChi2TPC(bool isAOD, Float_t chi2cut);

/**
 * Initialize cut variations in the tree creator task. Selected sample is always a subset of the hybrid track sample -
 * relaxing cuts is not possible.
 * \param task tree creator task
 * \param isAOD Flag for AOD
 */
void ReducedEventCutVariation(HighPtTracks::AliReducedHighPtEventCreator *task,  Bool_t isAOD){
  AliEmcalTrackSelection *trackcuts;
  if((trackcuts = CreateDefaultTrackCuts(isAOD))) task->AddVirtualTrackSelection(trackcuts, 0);
  if((trackcuts = CreateHybridTrackCuts(isAOD))) task->AddVirtualTrackSelection(trackcuts, 1);
  if((trackcuts = CreateDefaultTrackCutsChangeFiducialAreaCut(isAOD))) task->AddVirtualTrackSelection(trackcuts, 2);
  if((trackcuts = CreateDefaultTrackCutsChangeCrossRows(isAOD, 70))) task->AddVirtualTrackSelection(trackcuts, 3);
  if((trackcuts = CreateDefaultTrackCutsChangeCrossRows(isAOD, 100))) task->AddVirtualTrackSelection(trackcuts, 4);
  if((trackcuts = CreateDefaultTrackCutsChangeCrossRows(isAOD, 140))) task->AddVirtualTrackSelection(trackcuts, 5);
  if((trackcuts = CreateDefaultTrackCutsChangeDCAZ(isAOD, 0.5))) task->AddVirtualTrackSelection(trackcuts, 6);
  if((trackcuts = CreateDefaultTrackCutsChangeDCAZ(isAOD, 1.))) task->AddVirtualTrackSelection(trackcuts, 7);
  if((trackcuts = CreateDefaultTrackCutsChangeDCAZ(isAOD, 4.))) task->AddVirtualTrackSelection(trackcuts, 8);
  if((trackcuts = CreateDefaultTrackCutsChangeDCAR(isAOD, "0.0382+0.0350/pt^1.01"))) task->AddVirtualTrackSelection(trackcuts, 9);
  if((trackcuts = CreateDefaultTrackCutsChangeDCAR(isAOD, "0.0082+0.0350/pt^1.01"))) task->AddVirtualTrackSelection(trackcuts, 10);
  if((trackcuts = CreateDefaultTrackCutsChangeDCAR(isAOD, "0.0182+0.0650/pt^1.01"))) task->AddVirtualTrackSelection(trackcuts, 11);
  if((trackcuts = CreateDefaultTrackCutsChangeDCAR(isAOD, "0.0182+0.0150/pt^1.01"))) task->AddVirtualTrackSelection(trackcuts, 12);
  if((trackcuts = CreateDefaultTrackCutsChangeGoldenCut(isAOD, 56))) task->AddVirtualTrackSelection(trackcuts, 13);
  if((trackcuts = CreateDefaultTrackCutsChangeGoldenCut(isAOD, 18))) task->AddVirtualTrackSelection(trackcuts, 14);
  if((trackcuts = CreateDefaultTrackCutsChangeChi2ITS(isAOD, 56))) task->AddVirtualTrackSelection(trackcuts, 15);
  if((trackcuts = CreateDefaultTrackCutsChangeChi2ITS(isAOD, 18))) task->AddVirtualTrackSelection(trackcuts, 16);
  if((trackcuts = CreateDefaultTrackCutsChangeChi2TPC(isAOD, 6))) task->AddVirtualTrackSelection(trackcuts, 17);
  if((trackcuts = CreateDefaultTrackCutsChangeChi2TPC(isAOD, 2))) task->AddVirtualTrackSelection(trackcuts, 18);
}

/**
 * Hybrid track selection
 * \param isAOD Flag for AOD
 * \return virtual track selection
 */
AliEmcalTrackSelection *CreateHybridTrackCuts(bool isAOD){
  AliEmcalTrackSelection * trackSelection(NULL);
  if(isAOD){
    // Purely use filter bits
    AliEmcalTrackSelectionAOD *aodsel = new AliEmcalTrackSelectionAOD();
    aodsel->AddFilterBit(256);
    aodsel->AddFilterBit(512);
    trackSelection = aodsel;
  } else {
    AliESDtrackCuts* hybridTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);
    hybridTrackCuts->SetName("Global Hybrid tracks, loose DCA");
    hybridTrackCuts->SetMaxDCAToVertexXY(2.4);
    hybridTrackCuts->SetMaxDCAToVertexZ(3.2);
    hybridTrackCuts->SetDCAToVertex2D(kTRUE);
    hybridTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
    hybridTrackCuts->SetMaxFractionSharedTPCClusters(0.4);
    trackSelection = new AliEmcalTrackSelectionESD(hybridTrackCuts);
  }
  return trackSelection;
}

/**
 * Default track selection
 * \param isAOD Flag for AOD
 * \return virtual track selection
 */
AliEmcalTrackSelection *CreateDefaultTrackCuts(bool isAOD){
  AliEmcalTrackSelection * trackSelection(NULL);
  if(isAOD){
    AliEmcalTrackSelectionAOD *aodsel = new AliEmcalTrackSelectionAOD();
    aodsel->AddFilterBit(AliAODTrack::kTrkGlobal);
    PWGJE::EMCALJetTasks::AliEMCalTriggerExtraCuts *extraCuts = new PWGJE::EMCALJetTasks::AliEMCalTriggerExtraCuts();
    extraCuts->SetMinTPCCrossedRows(120);
    aodsel->AddTrackCuts(extraCuts);
    trackSelection = aodsel;
  } else {
    AliESDtrackCuts *standardTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(true, 1);
    standardTrackCuts->SetName("Standard Track cuts");
    standardTrackCuts->SetMinNCrossedRowsTPC(120);
    standardTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
    trackSelection = new AliEmcalTrackSelectionESD(standardTrackCuts);
  }
  return trackSelection;
}

/**
 * Crossed rows cut variation
 * \param isAOD Flag for AOD
 * \return virtual track selection
 */
AliEmcalTrackSelection *CreateDefaultTrackCutsChangeCrossRows(bool isAOD, int minCrossedRows){
  AliEmcalTrackSelection * trackSelection(NULL);
  if(isAOD){
    AliEmcalTrackSelectionAOD *aodsel = new AliEmcalTrackSelectionAOD();
    aodsel->AddFilterBit(AliAODTrack::kTrkGlobal);
    PWGJE::EMCALJetTasks::AliEMCalTriggerExtraCuts *extraCuts = new PWGJE::EMCALJetTasks::AliEMCalTriggerExtraCuts();
    extraCuts->SetMinTPCCrossedRows(minCrossedRows);
    aodsel->AddTrackCuts(extraCuts);
    trackSelection = aodsel;
  } else {
    AliESDtrackCuts *standardTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(true, 1);
    standardTrackCuts->SetName("Standard Track cuts");
    standardTrackCuts->SetMinNCrossedRowsTPC(minCrossedRows);
    standardTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
    trackSelection = new AliEmcalTrackSelectionESD(standardTrackCuts);
  }
  return trackSelection;
}

/**
 * Default cut + fiducial area cut
 * \param isAOD Flag for AOD
 * \return virtual track selection
 */
AliEmcalTrackSelection *CreateDefaultTrackCutsChangeFiducialAreaCut(bool isAOD ){
  AliEmcalTrackSelection * trackSelection(NULL);
  if(isAOD){
    AliEmcalTrackSelectionAOD *aodsel = new AliEmcalTrackSelectionAOD();
    aodsel->AddFilterBit(AliAODTrack::kTrkGlobal);
    PWGJE::EMCALJetTasks::AliEMCalTriggerExtraCuts *extraCuts = new PWGJE::EMCALJetTasks::AliEMCalTriggerExtraCuts();
    extraCuts->SetMinTPCCrossedRows(120);
    extraCuts->SetMinTPCTrackLengthCut();
    aodsel->AddTrackCuts(extraCuts);
    trackSelection = aodsel;
  } else {
    AliESDtrackCuts *standardTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(true, 1);
    standardTrackCuts->SetName("Standard Track cuts");
    standardTrackCuts->SetMinNCrossedRowsTPC(120);
    standardTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
    trackSelection = new AliEmcalTrackSelectionESD(standardTrackCuts);
    PWGJE::EMCALJetTasks::AliEMCalTriggerExtraCuts *extraCuts = new PWGJE::EMCALJetTasks::AliEMCalTriggerExtraCuts();
    extraCuts->SetMinTPCTrackLengthCut();
    trackSelection->AddTrackCuts(extraCuts);
  }
  return trackSelection;
}

/**
 * DCAZ variation
 * \param isAOD Flag for AOD
 * \return virtual track selection
 */
AliEmcalTrackSelection *CreateDefaultTrackCutsChangeDCAZ(bool isAOD, Double_t dcaZ){
  AliEmcalTrackSelection * trackSelection(NULL);
  if(isAOD){
    return NULL;
  } else {
    AliESDtrackCuts *standardTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(true, 1);
    standardTrackCuts->SetName("Standard Track cuts");
    standardTrackCuts->SetMinNCrossedRowsTPC(120);
    standardTrackCuts->SetMaxDCAToVertexZ(dcaZ);
    standardTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
    trackSelection = new AliEmcalTrackSelectionESD(standardTrackCuts);
  }
  return trackSelection;
}

/**
 * DCAR cut variation
 * \param isAOD Flag for AOD
 * \return virtual track selection
 */
AliEmcalTrackSelection *CreateDefaultTrackCutsChangeDCAR(bool isAOD, TString model){
  AliEmcalTrackSelection * trackSelection(NULL);
  if(isAOD){
    return NULL;
  } else {
    AliESDtrackCuts *standardTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(true, 1);
    standardTrackCuts->SetName("Standard Track cuts");
    standardTrackCuts->SetMinNCrossedRowsTPC(120);
    standardTrackCuts->SetMaxDCAToVertexXYPtDep(model.Data());
    trackSelection = new AliEmcalTrackSelectionESD(standardTrackCuts);
  }
  return trackSelection;
}

/**
 * Golden cut variation
 * \param isAOD Flag for AOD
 * \return virtual track selection
 */
AliEmcalTrackSelection *CreateDefaultTrackCutsChangeGoldenCut(bool isAOD, Float_t chi2cut){
  AliEmcalTrackSelection * trackSelection(NULL);
  if(isAOD){
    return NULL;
  } else {
    AliESDtrackCuts *standardTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(true, 1);
    standardTrackCuts->SetName("Standard Track cuts");
    standardTrackCuts->SetMinNCrossedRowsTPC(120);
    standardTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
    standardTrackCuts->SetMaxChi2TPCConstrainedGlobal(chi2cut);
    trackSelection = new AliEmcalTrackSelectionESD(standardTrackCuts);
  }
  return trackSelection;
}

/**
 * Chi2/ITS cluster cut variation
 * \param isAOD Flag for AOD
 * \return virtual track selection
 */
AliEmcalTrackSelection *CreateDefaultTrackCutsChangeChi2ITS(bool isAOD, Float_t chi2cut){
  AliEmcalTrackSelection * trackSelection(NULL);
  if(isAOD){
    return NULL;
  } else {
    AliESDtrackCuts *standardTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(true, 1);
    standardTrackCuts->SetName("Standard Track cuts");
    standardTrackCuts->SetMinNCrossedRowsTPC(120);
    standardTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
    standardTrackCuts->SetMaxChi2PerClusterITS(chi2cut);
    trackSelection = new AliEmcalTrackSelectionESD(standardTrackCuts);
  }
  return trackSelection;
}

/**
 * Chi2/TPC cluster cut variation
 * \param isAOD Flag for AOD
 * \return virtual track selection
 */
AliEmcalTrackSelection *CreateDefaultTrackCutsChangeChi2TPC(bool isAOD, Float_t chi2cut){
  AliEmcalTrackSelection * trackSelection(NULL);
  if(isAOD){
    return NULL;
  } else {
    AliESDtrackCuts *standardTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(true, 1);
    standardTrackCuts->SetName("Standard Track cuts");
    standardTrackCuts->SetMinNCrossedRowsTPC(120);
    standardTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
    standardTrackCuts->SetMaxChi2PerClusterTPC(chi2cut);
    trackSelection = new AliEmcalTrackSelectionESD(standardTrackCuts);
  }
  return trackSelection;
}
