/**************************************************************************
 * Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
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
#include <bitset>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

#include <TArrayD.h>
#include <TClonesArray.h>
#include <THashList.h>
#include <THistManager.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TString.h>

#include "AliAnalysisUtils.h"
#include "AliCentrality.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliEmcalTriggerOfflineSelection.h"
#include "AliESDEvent.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliVCluster.h"
#include "AliVVertex.h"
#include "AliMultSelection.h"
#include "AliMultEstimator.h"

#include "AliAnalysisTaskEmcalClustersRef.h"

ClassImp(EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalClustersRef)

namespace EMCalTriggerPtAnalysis {

/**
 * Dummy (I/O) constructor
 */
AliAnalysisTaskEmcalClustersRef::AliAnalysisTaskEmcalClustersRef() :
    AliAnalysisTaskSE(),
    fAnalysisUtil(NULL),
    fTriggerSelection(NULL),
    fHistos(NULL),
    fGeometry(NULL),
    fClusterContainer(""),
    fRequestAnalysisUtil(kTRUE),
    fTriggerStringFromPatches(kFALSE),
    fCentralityRange(-999., 999.),
    fVertexRange(-999., 999.)
{
}

/**
 * Named constructor
 * @param name Name of the task
 */
AliAnalysisTaskEmcalClustersRef::AliAnalysisTaskEmcalClustersRef(const char *name) :
    AliAnalysisTaskSE(name),
    fAnalysisUtil(),
    fTriggerSelection(NULL),
    fHistos(NULL),
    fGeometry(NULL),
    fClusterContainer(""),
    fRequestAnalysisUtil(kTRUE),
    fTriggerStringFromPatches(kFALSE),
    fCentralityRange(-999., 999.),
    fVertexRange(-999., 999.)
{
  DefineOutput(1, TList::Class());
}

/**
 * Destructor
 */
AliAnalysisTaskEmcalClustersRef::~AliAnalysisTaskEmcalClustersRef() {
  if(fTriggerSelection) delete fTriggerSelection;
}

/**
 * Creates output histograms: distribution of cluster energy for different trigger classes and number of events
 */
void AliAnalysisTaskEmcalClustersRef::UserCreateOutputObjects(){
  AliInfo(Form("Creating histograms for task %s\n", GetName()));
  fAnalysisUtil = new AliAnalysisUtils;

  TArrayD energybinning;
  CreateEnergyBinning(energybinning);
  TArrayD smbinning(14); CreateLinearBinning(smbinning, 21, -0.5, 20.5);
  TArrayD etabinning; CreateLinearBinning(etabinning, 100, -0.7, 0.7);
  fHistos = new THistManager("Ref");
  TString triggers[18] = {
      "MB", "EMC7", "DMC7",
      "EJ1", "EJ2", "EG1", "EG2", "DJ1", "DJ2", "DG1", "DG2",
      "MBexcl", "EMC7excl", "DMC7excl", "EG2excl", "EJ2excl", "DG2excl", "DJ2excl"
  };
  Double_t encuts[5] = {1., 2., 5., 10., 20.};
  Int_t sectorsWithEMCAL[10] = {4, 5, 6, 7, 8, 9, 13, 14, 15, 16};
  for(TString *trg = triggers; trg < triggers + sizeof(triggers)/sizeof(TString); trg++){
    fHistos->CreateTH1(Form("hEventCount%s", trg->Data()), Form("Event count for trigger class %s", trg->Data()), 1, 0.5, 1.5);
    fHistos->CreateTH1(Form("hEventCentrality%s", trg->Data()), Form("Event centrality for trigger class %s", trg->Data()), 103, -2., 101.);
    fHistos->CreateTH1(Form("hClusterEnergy%s", trg->Data()), Form("Cluster energy for trigger class %s", trg->Data()), energybinning);
    fHistos->CreateTH1(Form("hClusterET%s", trg->Data()), Form("Cluster transverse energy for trigger class %s", trg->Data()), energybinning);
    fHistos->CreateTH1(Form("hClusterEnergyFired%s", trg->Data()), Form("Cluster energy for trigger class %s, firing the trigger", trg->Data()), energybinning);
    fHistos->CreateTH1(Form("hClusterETFired%s", trg->Data()), Form("Cluster transverse energy for trigger class %s, firing the trigger", trg->Data()), energybinning);
    fHistos->CreateTH2(Form("hClusterEnergySM%s", trg->Data()), Form("Cluster energy versus supermodule for trigger class %s", trg->Data()), smbinning, energybinning);
    fHistos->CreateTH2(Form("hClusterETSM%s", trg->Data()), Form("Cluster transverse energy versus supermodule for trigger class %s", trg->Data()), smbinning, energybinning);
    fHistos->CreateTH2(Form("hClusterEnergyFiredSM%s", trg->Data()), Form("Cluster energy versus supermodule for trigger class %s, firing the trigger", trg->Data()), smbinning, energybinning);
    fHistos->CreateTH2(Form("hClusterETFiredSM%s", trg->Data()), Form("Cluster transverse energy versus supermodule for trigger class %s, firing the trigger", trg->Data()), smbinning, energybinning);
    fHistos->CreateTH2(Form("hEtaEnergy%s", trg->Data()), Form("Cluster energy vs. eta for trigger class %s", trg->Data()), etabinning, energybinning);
    fHistos->CreateTH2(Form("hEtaET%s", trg->Data()), Form("Cluster transverse energy vs. eta for trigger class %s", trg->Data()), etabinning, energybinning);
    fHistos->CreateTH2(Form("hEtaEnergyFired%s", trg->Data()), Form("Cluster energy vs. eta for trigger class %s, firing the trigger", trg->Data()), etabinning, energybinning);
    fHistos->CreateTH2(Form("hEtaETFired%s", trg->Data()), Form("Cluster transverse energy vs. eta for trigger class %s, firing the trigger", trg->Data()), etabinning, energybinning);
    for(int ism = 0; ism < 20; ism++){
      fHistos->CreateTH2(Form("hEtaEnergySM%d%s", ism, trg->Data()), Form("Cluster energy vs. eta in Supermodule %d for trigger %s", ism, trg->Data()), etabinning, energybinning);
      fHistos->CreateTH2(Form("hEtaETSM%d%s", ism, trg->Data()), Form("Cluster transverse energy vs. eta in Supermodule %d for trigger %s", ism, trg->Data()), etabinning, energybinning);
      fHistos->CreateTH2(Form("hEtaEnergyFiredSM%d%s", ism, trg->Data()), Form("Cluster energy vs. eta in Supermodule %d for trigger %s, firing the trigger", ism, trg->Data()), etabinning, energybinning);
      fHistos->CreateTH2(Form("hEtaETFiredSM%d%s", ism, trg->Data()), Form("Cluster transverse energy vs. eta in Supermodule %d for trigger %s, firing the trigger", ism, trg->Data()), etabinning, energybinning);
    }
    for(int isec = 0; isec < 10; isec++){
      fHistos->CreateTH2(Form("hEtaEnergySec%d%s", sectorsWithEMCAL[isec], trg->Data()), Form("Cluster energy vs.eta in tracking sector %d for trigger %s", sectorsWithEMCAL[isec], trg->Data()), etabinning, energybinning);
      fHistos->CreateTH2(Form("hEtaETSec%d%s", sectorsWithEMCAL[isec], trg->Data()), Form("Cluster transverse energy vs.eta in tracking sector %d for trigger %s", sectorsWithEMCAL[isec], trg->Data()), etabinning, energybinning);
      fHistos->CreateTH2(Form("hEtaEnergyFiredSec%d%s", sectorsWithEMCAL[isec], trg->Data()), Form("Cluster energy vs.eta in tracking sector %d for trigger %s, firing the trigger", sectorsWithEMCAL[isec], trg->Data()), etabinning, energybinning);
      fHistos->CreateTH2(Form("hEtaETFiredSec%d%s", sectorsWithEMCAL[isec], trg->Data()), Form("Cluster transverse energy vs.eta in tracking sector %d for trigger %s, firing the trigger", sectorsWithEMCAL[isec], trg->Data()), etabinning, energybinning);
    }
    for(int ien = 0; ien < 5; ien++){
      fHistos->CreateTH2(Form("hEtaPhi%dG%s", static_cast<int>(encuts[ien]), trg->Data()), Form("cluster #eta-#phi map for clusters with energy larger than %f GeV/c for trigger class %s", encuts[ien], trg->Data()), 100, -0.7, 0.7, 200, 0, 2*TMath::Pi());
      fHistos->CreateTH2(Form("hEtaPhiFired%dG%s", static_cast<int>(encuts[ien]), trg->Data()), Form("cluster #eta-#phi map for clusters fired the trigger with energy larger than %f GeV/c for trigger class %s", encuts[ien], trg->Data()), 200, -0.7, 0.7, 200, 0, 2*TMath::Pi());
    }
  }
  PostData(1, fHistos->GetListOfHistograms());
  AliDebug(1, "End creating histograms");
}


/**
 *
 * @param
 */
void AliAnalysisTaskEmcalClustersRef::UserExec(Option_t *){
  AliDebug(1, Form("%s: UserExec start\n", GetName()));
  if(!fGeometry){
    fGeometry = AliEMCALGeometry::GetInstance();
    if(!fGeometry)
      fGeometry = AliEMCALGeometry::GetInstanceFromRunNumber(InputEvent()->GetRunNumber());
  }
  TString triggerstring = "";
  TClonesArray *triggerpatches = dynamic_cast<TClonesArray *>(fInputEvent->FindListObject("EmcalTriggers"));

  if(fTriggerStringFromPatches){
    triggerstring = GetFiredTriggerClassesFromPatches(triggerpatches);
  } else {
    triggerstring = fInputEvent->GetFiredTriggerClasses();
  }

  UInt_t selectionstatus = fInputHandler->IsEventSelected();
  std::stringstream triggerdebug;
  triggerdebug << "Offline bits: " << std::bitset<sizeof(UInt_t) * 8>(selectionstatus);
  AliDebug(2, triggerdebug.str().c_str());
  Bool_t isMinBias = selectionstatus & AliVEvent::kINT7,
      isEJ1 = (selectionstatus & AliVEvent::kEMCEJE) && triggerstring.Contains("EJ1"),
      isEJ2 = (selectionstatus & AliVEvent::kEMCEJE) && triggerstring.Contains("EJ2"),
      isEG1 = (selectionstatus & AliVEvent::kEMCEGA) && triggerstring.Contains("EG1"),
      isEG2 = (selectionstatus & AliVEvent::kEMCEGA) && triggerstring.Contains("EG2"),
      isEMC7 = (selectionstatus & AliVEvent::kEMC7) && triggerstring.Contains("EMC7"),
      isDJ1 = (selectionstatus & AliVEvent::kEMCEJE) && triggerstring.Contains("DJ1"),
      isDJ2 = (selectionstatus & AliVEvent::kEMCEJE) && triggerstring.Contains("DJ2"),
      isDG1 = (selectionstatus & AliVEvent::kEMCEGA) && triggerstring.Contains("DG1"),
      isDG2 = (selectionstatus & AliVEvent::kEMCEGA) && triggerstring.Contains("DG2"),
      isDMC7 = (selectionstatus & AliVEvent::kEMC7) && triggerstring.Contains("DMC7");
  if(triggerpatches && fTriggerSelection){
      isEJ1 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgEJ1, triggerpatches);
      isEJ2 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgEJ2, triggerpatches);
      isEG1 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgEG1, triggerpatches);
      isEG2 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgEG2, triggerpatches);
      isEMC7 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgEL0, triggerpatches);
      isDJ1 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgDJ1, triggerpatches);
      isDJ2 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgDJ2, triggerpatches);
      isDG1 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgDG1, triggerpatches);
      isDG2 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgDG2, triggerpatches);
      isDMC7 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgDL0, triggerpatches);
  }
  if(!(isMinBias || isEMC7 || isEG1 || isEG2 || isEJ1 || isEJ2 || isDMC7 || isDG1 || isDG2 || isDJ1 || isDJ2)){
    AliDebug(1, Form("%s: Reject trigger\n", GetName()));
    return;
  }
  AliDebug(1, "Event selected");
  AliMultSelection *mult = dynamic_cast<AliMultSelection *>(InputEvent()->FindListObject("MultSelection"));
  if(!mult) std::cout << "Multiplicity selection not found" << std::endl;
  double centrality =  mult ? mult->GetEstimator("V0M")->GetPercentile() : -1;
  AliDebug(1, Form("%s: Centrality %f\n", GetName(), centrality));
  if(!fCentralityRange.IsInRange(centrality)){
    AliDebug(1, Form("%s: reject centrality: %f\n", GetName(), centrality));
    return;
  } else {
    AliDebug(1, Form("%s: select centrality %f\n", GetName(), centrality));
  }
  const AliVVertex *vtx = fInputEvent->GetPrimaryVertex();
  if(!vtx) vtx = fInputEvent->GetPrimaryVertexSPD();
  //if(!fInputEvent->IsPileupFromSPD(3, 0.8, 3., 2., 5.)) return;         // reject pileup event
  if(vtx->GetNContributors() < 1){
    AliDebug(1, Form("%s: Reject contributors\n", GetName()));
    return;
  }
  // Fill reference distribution for the primary vertex before any z-cut
  if(fRequestAnalysisUtil){
    AliDebug(1, Form("%s: Reject analysis util\n", GetName()));
    if(fInputEvent->IsA() == AliESDEvent::Class() && fAnalysisUtil->IsFirstEventInChunk(fInputEvent)) return;
    if(!fAnalysisUtil->IsVertexSelected2013pA(fInputEvent)) return;       // Apply new vertex cut
    if(fAnalysisUtil->IsPileUpEvent(fInputEvent)) return;       // Apply new vertex cut
  }
  // Apply vertex z cut
  if(!fVertexRange.IsInRange(vtx->GetZ())){
    AliDebug(1, Form("%s: Reject z", GetName()));
    return;
  }
  AliDebug(1, Form("%s: Event Selected\n", GetName()));

  // Fill Event counter and reference vertex distributions for the different trigger classes
  if(isMinBias){
    fHistos->FillTH1("hEventCountMB", 1);
    fHistos->FillTH1("hEventCentralityMB", centrality);
    if(!(isEMC7 || isDMC7 || isEJ1 || isEJ2 || isEG1 || isEG2 || isDJ1 || isDJ2 || isDG1 || isDG2)){
      fHistos->FillTH1("hEventCountMBexcl", 1);
      fHistos->FillTH1("hEventCentralityMBexcl", centrality);
    }
  }

  if(isEMC7){
    fHistos->FillTH1("hEventCountEMC7", 1);
    fHistos->FillTH1("hEventCentralityEMC7", centrality);
    if(!(isEJ1 || isEJ2 || isEG1 || isEG2)){
      fHistos->FillTH1("hEventCountEMC7excl", 1);
      fHistos->FillTH1("hEventCentralityEMC7excl", centrality);
    }
  }
  if(isDMC7){
    fHistos->FillTH1("hEventCountDMC7", 1);
    fHistos->FillTH1("hEventCentralityDMC7", centrality);
    if(!(isDJ1 || isDJ2 || isDG1 || isDG2)){
      fHistos->FillTH1("hEventCountDMC7excl", 1);
      fHistos->FillTH1("hEventCentralityDMC7excl", centrality);
    }
  }

  if(isEJ2){
    fHistos->FillTH1("hEventCountEJ2", 1);
    fHistos->FillTH1("hEventCentralityEJ2", centrality);
    // Check for exclusive classes
    if(!isEJ1){
      fHistos->FillTH1("hEventCountEJ2excl", 1);
      fHistos->FillTH1("hEventCentralityEJ2excl", centrality);
    }
  }
  if(isDJ2){
    fHistos->FillTH1("hEventCountDJ2", 1);
    fHistos->FillTH1("hEventCentralityDJ2", centrality);
    // Check for exclusive classes
    if(!isDJ1){
      fHistos->FillTH1("hEventCountDJ2excl", 1);
      fHistos->FillTH1("hEventCentralityDJ2excl", centrality);
    }
  }

  if(isEJ1){
    fHistos->FillTH1("hEventCountEJ1", 1);
    fHistos->FillTH1("hEventCentralityEJ1", centrality);
  }
  if(isDJ1){
    fHistos->FillTH1("hEventCountDJ1", 1);
    fHistos->FillTH1("hEventCentralityDJ1", centrality);
  }

  if(isEG2){
    fHistos->FillTH1("hEventCountEG2", 1);
    fHistos->FillTH1("hEventCentralityEG2", centrality);
    // Check for exclusive classes
    if(!(isEG1)){
      fHistos->FillTH1("hEventCountEG2excl", 1);
      fHistos->FillTH1("hEventCentralityEG2excl", centrality);
    }
  }
  if(isDG2){
    fHistos->FillTH1("hEventCountDG2", 1);
    fHistos->FillTH1("hEventCentralityDG2", centrality);
    // Check for exclusive classes
    if(!(isDG1)){
      fHistos->FillTH1("hEventCountDG2excl", 1);
      fHistos->FillTH1("hEventCentralityDG2excl", centrality);
    }
  }

  if(isEG1){
    fHistos->FillTH1("hEventCountEG1", 1);
    fHistos->FillTH1("hEventCentralityEG1", centrality);
  }
  if(isDG1){
    fHistos->FillTH1("hEventCountDG1", 1);
    fHistos->FillTH1("hEventCentralityDG1", centrality);
  }

  /*
  TList *objects = fInputEvent->GetList();
  for(TIter objiter = TIter(objects).Begin(); objiter != TIter::End(); ++ objiter){
    printf("Object %s\n", (*objiter)->GetName());
  }
  */

  TObjArray clusterEvent(1000);

  TCollection *clusterArray = dynamic_cast<TClonesArray *>(fInputEvent->FindListObject(fClusterContainer.Data()));
  if(!clusterArray){
    AliError(Form("Cluster array with name %s not found in the event", fClusterContainer.Data()));
    for(int icl = 0; icl < fInputEvent->GetNumberOfCaloClusters(); icl++){
      clusterEvent.Add(fInputEvent->GetCaloCluster(icl));
    }
    clusterArray =  &clusterEvent;
  }

  Double_t vertexpos[3];
  fInputEvent->GetPrimaryVertex()->GetXYZ(vertexpos);

  Double_t energy, eta, phi;
  for(TIter clustIter = TIter(clusterArray).Begin(); clustIter != TIter::End(); ++clustIter){
    AliVCluster *clust = static_cast<AliVCluster *>(*clustIter);
    if(!clust->IsEMCAL()) continue;
    if(clust->GetIsExotic()) continue;

    TLorentzVector posvec;
    energy = clust->GetNonLinCorrEnergy();
    clust->GetMomentum(posvec, vertexpos);
    eta = posvec.Eta();
    phi = posvec.Phi();

    // fill histograms allEta
    if(isMinBias){
      FillClusterHistograms("MB", energy, posvec.Et(), eta, phi, NULL);
      if(!(isEMC7 || isDMC7 || isEJ1 || isEJ2 || isEG1 || isEG2 || isDJ1 || isDJ2 || isDG1 || isDG2)){
        FillClusterHistograms("MBexcl", energy, posvec.Et(), eta, phi, NULL);
      }
    }
    if(isEMC7){
      FillClusterHistograms("EMC7", energy, posvec.Et(), eta, phi, NULL);
      if(!(isEJ1 || isEJ2 || isEG1 || isEG2)){
        FillClusterHistograms("EMC7excl", energy, posvec.Et(), eta, phi, NULL);
      }
    }
    if(isDMC7){
      FillClusterHistograms("DMC7", energy, posvec.Et(), eta, phi, NULL);
      if(!(isDJ1 || isDJ2 || isDG1 || isDG2)){
        FillClusterHistograms("DMC7excl", energy, posvec.Et(), eta, phi, NULL);
      }
    }
    if(isEJ2){
      TList ej2patches;
      FindPatchesForTrigger("EJ2", triggerpatches, ej2patches);
      FillClusterHistograms("EJ2", energy, posvec.Et(), eta, phi, &ej2patches);
      // check for exclusive classes
      if(!isEJ1){
        FillClusterHistograms("EJ2excl", energy, posvec.Et(), eta, phi, &ej2patches);
      }
    }
    if(isDJ2){
      TList dj2patches;
      FindPatchesForTrigger("DJ2", triggerpatches, dj2patches);
      FillClusterHistograms("DJ2", energy, posvec.Et(), eta, phi, &dj2patches);
      // check for exclusive classes
      if(!isDJ1){
        FillClusterHistograms("DJ2excl", energy, posvec.Et(), eta, phi, &dj2patches);
      }
    }
    if(isEJ1){
      TList ej1patches;
      FindPatchesForTrigger("EJ1", triggerpatches, ej1patches);
      FillClusterHistograms("EJ1", energy, posvec.Et(), eta, phi, &ej1patches);
    }
    if(isDJ1){
      TList dj1patches;
      FindPatchesForTrigger("DJ1", triggerpatches, dj1patches);
      FillClusterHistograms("DJ1", energy, posvec.Et(), eta, phi, &dj1patches);
    }
    if(isEG2){
      TList eg2patches;
      FindPatchesForTrigger("EG2", triggerpatches, eg2patches);
      FillClusterHistograms("EG2", energy, posvec.Et(), eta, phi, &eg2patches);
      // check for exclusive classes
      if(!isEG1){
        FillClusterHistograms("EG2excl", energy, posvec.Et(), eta, phi, &eg2patches);
      }
    }
    if(isDG2){
      TList dg2patches;
      FindPatchesForTrigger("DG2", triggerpatches, dg2patches);
      FillClusterHistograms("DG2", energy, posvec.Et(), eta, phi, &dg2patches);
      // check for exclusive classes
      if(!isDG1){
        FillClusterHistograms("DG2excl", energy, posvec.Et(), eta, phi, &dg2patches);
      }
    }
    if(isEG1){
      TList eg1patches;
      FindPatchesForTrigger("EG1", triggerpatches, eg1patches);
      FillClusterHistograms("EG1", energy, posvec.Et(), eta, phi, &eg1patches);
    }
    if(isDG1){
      TList dg1patches;
      FindPatchesForTrigger("DG1", triggerpatches, dg1patches);
      FillClusterHistograms("DG1", energy, posvec.Et(), eta, phi, &dg1patches);
    }
  }
  PostData(1, fHistos->GetListOfHistograms());
}

void AliAnalysisTaskEmcalClustersRef::FillClusterHistograms(TString triggerclass, double energy, double transverseenergy, double eta, double phi, TList *triggerpatches){
  Bool_t hasTriggerPatch = triggerpatches  ? CorrelateToTrigger(eta, phi, triggerpatches) : kFALSE;
  Int_t supermoduleID = -1, sector = -1;
  fGeometry->SuperModuleNumberFromEtaPhi(eta, phi, supermoduleID);
  fHistos->FillTH1(Form("hClusterEnergy%s", triggerclass.Data()), energy);
  fHistos->FillTH1(Form("hClusterET%s", triggerclass.Data()), transverseenergy);
  fHistos->FillTH2(Form("hEtaEnergy%s", triggerclass.Data()), eta, energy);
  fHistos->FillTH2(Form("hEtaET%s", triggerclass.Data()), eta, transverseenergy);
  if(supermoduleID >= 0){
    fHistos->FillTH2(Form("hClusterEnergySM%s", triggerclass.Data()), supermoduleID, energy);
    fHistos->FillTH2(Form("hClusterETSM%s", triggerclass.Data()), supermoduleID, transverseenergy);
    fHistos->FillTH2(Form("hEtaEnergySM%d%s", supermoduleID, triggerclass.Data()), eta, energy);
    fHistos->FillTH2(Form("hEtaETSM%d%s", supermoduleID, triggerclass.Data()), eta, transverseenergy);
    if(supermoduleID < 12)
      sector = 4 + int(supermoduleID/2); // EMCAL
    else
      sector = 13 + int((supermoduleID-12)/2);  // DCAL
    fHistos->FillTH2(Form("hEtaEnergySec%d%s", sector, triggerclass.Data()), eta, energy);
    fHistos->FillTH2(Form("hEtaETSec%d%s", sector, triggerclass.Data()), eta, transverseenergy);
  }
  if(hasTriggerPatch){
    fHistos->FillTH1(Form("hClusterEnergyFired%s", triggerclass.Data()), energy);
    fHistos->FillTH1(Form("hClusterETFired%s", triggerclass.Data()), energy);
    fHistos->FillTH2(Form("hEtaEnergyFired%s", triggerclass.Data()), eta, energy);
    fHistos->FillTH2(Form("hEtaETFired%s", triggerclass.Data()), eta, energy);
    if(supermoduleID >= 0){
      fHistos->FillTH2(Form("hClusterEnergyFiredSM%s", triggerclass.Data()), supermoduleID, energy);
      fHistos->FillTH2(Form("hClusterETFiredSM%s", triggerclass.Data()), supermoduleID, transverseenergy);
      fHistos->FillTH2(Form("hEtaEnergyFiredSM%d%s", supermoduleID, triggerclass.Data()), eta, energy);
      fHistos->FillTH2(Form("hEtaETFiredSM%d%s", supermoduleID, triggerclass.Data()), eta, transverseenergy);
      fHistos->FillTH2(Form("hEtaEnergyFiredSec%d%s", sector, triggerclass.Data()), eta, energy);
      fHistos->FillTH2(Form("hEtaETFiredSec%d%s", sector, triggerclass.Data()), eta, transverseenergy);
    }
  }
  Double_t encuts[5] = {1., 2., 5., 10., 20.};
  for(int ien = 0; ien < 5; ien++){
    if(energy > encuts[ien]){
      fHistos->FillTH2(Form("hEtaPhi%dG%s", static_cast<int>(encuts[ien]), triggerclass.Data()), eta, phi);
      if(hasTriggerPatch){
        fHistos->FillTH2(Form("hEtaPhiFired%dG%s", static_cast<int>(encuts[ien]), triggerclass.Data()), eta, phi);
      }
    }
  }
}

/**
 * Create new energy binning
 * @param binning
 */
void AliAnalysisTaskEmcalClustersRef::CreateEnergyBinning(TArrayD& binning) const {
  std::vector<double> mybinning;
  std::map<double,double> definitions;
  definitions.insert(std::pair<double, double>(1, 0.05));
  definitions.insert(std::pair<double, double>(2, 0.1));
  definitions.insert(std::pair<double, double>(4, 0.2));
  definitions.insert(std::pair<double, double>(7, 0.5));
  definitions.insert(std::pair<double, double>(16, 1));
  definitions.insert(std::pair<double, double>(32, 2));
  definitions.insert(std::pair<double, double>(40, 4));
  definitions.insert(std::pair<double, double>(50, 5));
  definitions.insert(std::pair<double, double>(100, 10));
  definitions.insert(std::pair<double, double>(200, 20));
  double currentval = 0.;
  mybinning.push_back(currentval);
  for(std::map<double,double>::iterator id = definitions.begin(); id != definitions.end(); ++id){
    double limit = id->first, binwidth = id->second;
    while(currentval < limit){
      currentval += binwidth;
      mybinning.push_back(currentval);
    }
  }
  binning.Set(mybinning.size());
  int ib = 0;
  for(std::vector<double>::iterator it = mybinning.begin(); it != mybinning.end(); ++it)
    binning[ib++] = *it;
}

/**
 * Create any kind of linear binning from given ranges and stores it in the binning array.
 * @param binning output array
 * @param nbins Number of bins
 * @param min lower range
 * @param max upper range
 */
void AliAnalysisTaskEmcalClustersRef::CreateLinearBinning(TArrayD& binning, int nbins, double min, double max) const {
  double binwidth = (max-min)/static_cast<double>(nbins);
  binning.Set(nbins+1);
  binning[0] = min;
  double currentlimit = min + binwidth;
  for(int ibin = 0; ibin < nbins; ibin++){
    binning[ibin+1] = currentlimit;
    currentlimit += binwidth;
  }
}

/**
 * Check whether cluster is inside a trigger patch which has fired the trigger
 * @param etaclust \f$ \eta \f$ of the cluster at center
 * @param phiclust \f$ \phi \f$ of the cluster at center
 * @param triggerpatches List of trigger patches which have fired the trigger
 * @return True if the cluster can be correlated to a triggerpatch fired the trigger, false otherwise
 */
Bool_t AliAnalysisTaskEmcalClustersRef::CorrelateToTrigger(Double_t etaclust, Double_t phiclust, TList *triggerpatches) const {
  Bool_t hasfound = kFALSE;
  for(TIter patchIter = TIter(triggerpatches).Begin(); patchIter != TIter::End(); ++patchIter){
    Double_t boundaries[4];
    GetPatchBoundaries(*patchIter, boundaries);
    Double_t etamin = TMath::Min(boundaries[0], boundaries[1]),
        etamax = TMath::Max(boundaries[0], boundaries[1]),
        phimin = TMath::Min(boundaries[2], boundaries[3]),
        phimax = TMath::Max(boundaries[2], boundaries[3]);
    if(etaclust > etamin && etaclust < etamax && phiclust > phimin && phiclust < phimax){
      hasfound = kTRUE;
      break;
    }
  }
  return hasfound;
}

/**
 * Find all patches in an event which could have fired the trigger
 * Attention: This task groups into single shower triggers (L0, EG1, EG2) and jet triggers (EJ1 and EJ2).
 * Per convention the low threshold patch is selected. No energy cut should be applied in the trigger maker
 * @param triggerclass EMCAL trigger class firing
 * @param triggerpatches Trigger patches found in the event
 * @return List of patches which could have fired the trigger
 */
void AliAnalysisTaskEmcalClustersRef::FindPatchesForTrigger(TString triggerclass, const TClonesArray * triggerpatches, TList &foundtriggers) const {
  foundtriggers.Clear();
  if(!triggerpatches) return;
  AliEmcalTriggerOfflineSelection::EmcalTriggerClass myclass = AliEmcalTriggerOfflineSelection::kTrgEL0;
  if(triggerclass == "EG1") myclass = AliEmcalTriggerOfflineSelection::kTrgEG1;
  if(triggerclass == "EG2") myclass = AliEmcalTriggerOfflineSelection::kTrgEG2;
  if(triggerclass == "EJ1") myclass = AliEmcalTriggerOfflineSelection::kTrgEJ1;
  if(triggerclass == "EJ2") myclass = AliEmcalTriggerOfflineSelection::kTrgEJ2;
  if(triggerclass == "DMC7") myclass = AliEmcalTriggerOfflineSelection::kTrgDL0;
  if(triggerclass == "DG1") myclass = AliEmcalTriggerOfflineSelection::kTrgDG1;
  if(triggerclass == "DG2") myclass = AliEmcalTriggerOfflineSelection::kTrgDG2;
  if(triggerclass == "DJ1") myclass = AliEmcalTriggerOfflineSelection::kTrgDJ1;
  if(triggerclass == "DJ2") myclass = AliEmcalTriggerOfflineSelection::kTrgDJ2;
  AliEMCALTriggerPatchInfo *mypatch = NULL;
  for(TIter patchiter = TIter(triggerpatches).Begin(); patchiter != TIter::End(); ++patchiter){
    if(!IsOfflineSimplePatch(*patchiter)) continue;
    if(AliEmcalTriggerOfflineSelection::IsDCAL(myclass)){
      if(!SelectDCALPatch(*patchiter)) continue;
    } else {
      if(SelectDCALPatch(*patchiter)) continue;
    }
    if(AliEmcalTriggerOfflineSelection::IsSingleShower(myclass)){
      if(!SelectSingleShowerPatch(*patchiter)) continue;
    } else {
      if(!SelectJetPatch(*patchiter)) continue;
    }
    double threshold = fTriggerSelection ? fTriggerSelection->GetThresholdForTrigger(myclass) : -1;
    if(GetPatchEnergy(*patchiter) > threshold) foundtriggers.Add(mypatch);
  }
}

/**
 * Apply trigger selection using offline patches and trigger thresholds based on offline ADC Amplitude
 * @param triggerpatches Trigger patches found by the trigger maker
 * @return String with EMCAL trigger decision
 */
TString AliAnalysisTaskEmcalClustersRef::GetFiredTriggerClassesFromPatches(const TClonesArray* triggerpatches) const {
  TString triggerstring = "";
  if(!triggerpatches) return triggerstring;
  Int_t nEJ1 = 0, nEJ2 = 0, nEG1 = 0, nEG2 = 0, nDJ1 = 0, nDJ2 = 0, nDG1 = 0, nDG2 = 0;
  double  minADC_J1 = 260.,
          minADC_J2 = 127.,
          minADC_G1 = 140.,
          minADC_G2 = 89.;
  for(TIter patchIter = TIter(triggerpatches).Begin(); patchIter != TIter::End(); ++patchIter){
    AliEMCALTriggerPatchInfo *patch = dynamic_cast<AliEMCALTriggerPatchInfo *>(*patchIter);
    if(!patch->IsOfflineSimple()) continue;
    if(patch->IsJetHighSimple() && patch->GetADCOfflineAmp() > minADC_J1){
      if(patch->IsDCalPHOS()) nDJ1++;
      else nEJ1++;
    }
    if(patch->IsJetLowSimple() && patch->GetADCOfflineAmp() > minADC_J2){
      if(patch->IsDCalPHOS()) nDJ2++;
      else nEJ2++;
    }
    if(patch->IsGammaHighSimple() && patch->GetADCOfflineAmp() > minADC_G1){
      if(patch->IsDCalPHOS()) nDG1++;
      else nEG1++;
    }
    if(patch->IsGammaLowSimple() && patch->GetADCOfflineAmp() > minADC_G2){
      if(patch->IsDCalPHOS()) nDG2++;
      else nEG2++;
    }
  }
  if(nEJ1) triggerstring += "EJ1";
  if(nEJ2){
    if(triggerstring.Length()) triggerstring += ",";
    triggerstring += "EJ2";
  }
  if(nEG1){
    if(triggerstring.Length()) triggerstring += ",";
    triggerstring += "EG1";
  }
  if(nEG2){
    if(triggerstring.Length()) triggerstring += ",";
    triggerstring += "EG2";
  }
  if(nDJ1){
    if(triggerstring.Length()) triggerstring += ",";
    triggerstring += "DJ1";
  }
  if(nDJ2){
    if(triggerstring.Length()) triggerstring += ",";
    triggerstring += "DJ2";
  }
  if(nDG1){
    if(triggerstring.Length()) triggerstring += ",";
    triggerstring += "DG1";
  }
  if(nDG2){
    if(triggerstring.Length()) triggerstring += ",";
    triggerstring += "DG2";
  }
  return triggerstring;
}

void AliAnalysisTaskEmcalClustersRef::GetPatchBoundaries(TObject *o, Double_t *boundaries) const {
  AliEMCALTriggerPatchInfo *patch= dynamic_cast<AliEMCALTriggerPatchInfo *>(o);
  boundaries[0] = patch->GetEtaMin();
  boundaries[1] = patch->GetEtaMax();
  boundaries[2] = patch->GetPhiMin();
  boundaries[3] = patch->GetPhiMax();
}

bool AliAnalysisTaskEmcalClustersRef::IsOfflineSimplePatch(TObject *o) const {
  AliEMCALTriggerPatchInfo *patch = dynamic_cast<AliEMCALTriggerPatchInfo *>(o);
  return patch->IsOfflineSimple();
}

bool AliAnalysisTaskEmcalClustersRef::SelectDCALPatch(TObject *o) const {
  AliEMCALTriggerPatchInfo *patch = dynamic_cast<AliEMCALTriggerPatchInfo *>(o);
  return patch->GetRowStart() >= 64;
}

bool AliAnalysisTaskEmcalClustersRef::SelectSingleShowerPatch(TObject *o) const{
  AliEMCALTriggerPatchInfo *patch = dynamic_cast<AliEMCALTriggerPatchInfo *>(o);
  return patch->IsGammaLowSimple();
}

bool AliAnalysisTaskEmcalClustersRef::SelectJetPatch(TObject *o) const{
  AliEMCALTriggerPatchInfo *patch = dynamic_cast<AliEMCALTriggerPatchInfo *>(o);
  if(!patch->IsOfflineSimple()) return false;
  return patch->IsJetLowSimple();
}

double AliAnalysisTaskEmcalClustersRef::GetPatchEnergy(TObject *o) const {
  double energy = 0.;
  AliEMCALTriggerPatchInfo *patch = dynamic_cast<AliEMCALTriggerPatchInfo *>(o);
  energy = patch->GetPatchE();
  return energy;
}

} /* namespace EMCalTriggerPtAnalysis */
