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
#include <TLorentzVector.h>
#include <TMath.h>
#include <TString.h>

#include "AliAnalysisUtils.h"
#include "AliEMCALGeometry.h"
#include "AliEmcalTriggerPatchInfoAP.h"
#include "AliEmcalTriggerPatchInfoAPV1.h"
#include "AliEMCalHistoContainer.h"
#include "AliESDEvent.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliVCluster.h"
#include "AliVVertex.h"

#include "AliAnalysisTaskEmcalClustersRef.h"

ClassImp(EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalClustersRef)

namespace EMCalTriggerPtAnalysis {

/**
 * Dummy (I/O) constructor
 */
AliAnalysisTaskEmcalClustersRef::AliAnalysisTaskEmcalClustersRef() :
    AliAnalysisTaskSE(),
    fAnalysisUtil(NULL),
    fHistos(NULL),
    fGeometry(NULL),
    fClusterContainer(""),
    fRequestAnalysisUtil(kTRUE),
    fTriggerStringFromPatches(kFALSE)
{
  for(int itrg = 0; itrg < kECRntrig; itrg++){
    fOfflineEnergyThreshold[itrg] = -1;
  }
}

/**
 * Named constructor
 * @param name Name of the task
 */
AliAnalysisTaskEmcalClustersRef::AliAnalysisTaskEmcalClustersRef(const char *name) :
    AliAnalysisTaskSE(name),
    fAnalysisUtil(),
    fHistos(NULL),
    fGeometry(NULL),
    fClusterContainer(""),
    fRequestAnalysisUtil(kTRUE),
    fTriggerStringFromPatches(kFALSE)
{
  for(int itrg = 0; itrg < kECRntrig; itrg++){
    fOfflineEnergyThreshold[itrg] = -1;
  }
  DefineOutput(1, TList::Class());
}

/**
 * Destructor
 */
AliAnalysisTaskEmcalClustersRef::~AliAnalysisTaskEmcalClustersRef() {
}

/**
 * Creates output histograms: distribution of cluster energy for different trigger classes and number of events
 */
void AliAnalysisTaskEmcalClustersRef::UserCreateOutputObjects(){
  fAnalysisUtil = new AliAnalysisUtils;

  TArrayD energybinning;
  CreateEnergyBinning(energybinning);
  TArrayD smbinning(14); CreateLinearBinning(smbinning, 21, -0.5, 20.5);
  TArrayD etabinning; CreateLinearBinning(etabinning, 100, -0.7, 0.7);
  fHistos = new AliEMCalHistoContainer("Ref");
  TString triggers[18] = {
      "MB", "EMC7", "DMC7",
      "EJ1", "EJ2", "EG1", "EG2", "DJ1", "DJ2", "DG1", "DG2",
      "MBexcl", "EMC7excl", "DMC7excl", "EG2excl", "EJ2excl", "DG2excl", "DJ2excl"
  };
  Double_t encuts[5] = {1., 2., 5., 10., 20.};
  Int_t sectorsWithEMCAL[10] = {4, 5, 6, 7, 8, 9, 13, 14, 15, 16};
  for(TString *trg = triggers; trg < triggers + sizeof(triggers)/sizeof(TString); trg++){
    fHistos->CreateTH1(Form("hEventCount%s", trg->Data()), Form("Event count for trigger class %s", trg->Data()), 1, 0.5, 1.5);
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
    for(int isec = 0; isec < 01; isec++){
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
}


/**
 *
 * @param
 */
void AliAnalysisTaskEmcalClustersRef::UserExec(Option_t *){
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
      isEJ1 = (selectionstatus & AliVEvent::kEMCEJE) && triggerstring.Contains("EJ1") && IsOfflineSelected(kECREJ1, triggerpatches),
      isEJ2 = (selectionstatus & AliVEvent::kEMCEJE) && triggerstring.Contains("EJ2") && IsOfflineSelected(kECREJ2, triggerpatches),
      isEG1 = (selectionstatus & AliVEvent::kEMCEGA) && triggerstring.Contains("EG1") && IsOfflineSelected(kECREG1, triggerpatches),
      isEG2 = (selectionstatus & AliVEvent::kEMCEGA) && triggerstring.Contains("EG2") && IsOfflineSelected(kECREG2, triggerpatches),
      isEMC7 = (selectionstatus & AliVEvent::kEMC7) && triggerstring.Contains("EMC7") && IsOfflineSelected(kECREL0, triggerpatches),
      isDJ1 = (selectionstatus & AliVEvent::kEMCEJE) && triggerstring.Contains("DJ1") && IsOfflineSelected(kECREJ1, triggerpatches),
      isDJ2 = (selectionstatus & AliVEvent::kEMCEJE) && triggerstring.Contains("DJ2") && IsOfflineSelected(kECREJ2, triggerpatches),
      isDG1 = (selectionstatus & AliVEvent::kEMCEGA) && triggerstring.Contains("DG1") && IsOfflineSelected(kECREG1, triggerpatches),
      isDG2 = (selectionstatus & AliVEvent::kEMCEGA) && triggerstring.Contains("DG2") && IsOfflineSelected(kECREG2, triggerpatches),
      isDMC7 = (selectionstatus & AliVEvent::kEMC7) && triggerstring.Contains("DMC7") && IsOfflineSelected(kECREL0, triggerpatches);
  if(!(isMinBias || isEMC7 || isEG1 || isEG2 || isEJ1 || isEJ2 || isDMC7 || isDG1 || isDG2 || isDJ1 || isDJ2)) return;
  const AliVVertex *vtx = fInputEvent->GetPrimaryVertex();
  //if(!fInputEvent->IsPileupFromSPD(3, 0.8, 3., 2., 5.)) return;         // reject pileup event
  if(vtx->GetNContributors() < 1) return;
  // Fill reference distribution for the primary vertex before any z-cut
  if(fRequestAnalysisUtil){
    if(fInputEvent->IsA() == AliESDEvent::Class() && fAnalysisUtil->IsFirstEventInChunk(fInputEvent)) return;
    if(!fAnalysisUtil->IsVertexSelected2013pA(fInputEvent)) return;       // Apply new vertex cut
    if(fAnalysisUtil->IsPileUpEvent(fInputEvent)) return;       // Apply new vertex cut
  }
  // Apply vertex z cut
  if(vtx->GetZ() < -10. || vtx->GetZ() > 10.) return;

  // Fill Event counter and reference vertex distributions for the different trigger classes
  if(isMinBias){
    fHistos->FillTH1("hEventCountMB", 1);
    if(!(isEMC7 || isDMC7 || isEJ1 || isEJ2 || isEG1 || isEG2 || isDJ1 || isDJ2 || isDG1 || isDG2)){
      fHistos->FillTH1("hEventCountMBexcl", 1);
    }
  }

  if(isEMC7){
    fHistos->FillTH1("hEventCountEMC7", 1);
    if(!(isEJ1 || isEJ2 || isEG1 || isEG2)){
      fHistos->FillTH1("hEventCountEMC7excl", 1);
    }
  }
  if(isDMC7){
    fHistos->FillTH1("hEventCountDMC7", 1);
    if(!(isDJ1 || isDJ2 || isDG1 || isDG2)){
      fHistos->FillTH1("hEventCountEMC7excl", 1);
    }
  }

  if(isEJ2){
    fHistos->FillTH1("hEventCountEJ2", 1);
    // Check for exclusive classes
    if(!isEJ1){
      fHistos->FillTH1("hEventCountEJ2excl", 1);
    }
  }
  if(isDJ2){
    fHistos->FillTH1("hEventCountDJ2", 1);
    // Check for exclusive classes
    if(!isDJ1){
      fHistos->FillTH1("hEventCountDJ2excl", 1);
    }
  }

  if(isEJ1){
    fHistos->FillTH1("hEventCountEJ1", 1);
  }
  if(isDJ1){
    fHistos->FillTH1("hEventCountDJ1", 1);
  }

  if(isEG2){
    fHistos->FillTH1("hEventCountEG2", 1);
    // Check for exclusive classes
    if(!(isEG1)){
      fHistos->FillTH1("hEventCountEG2excl", 1);
    }
  }
  if(isDG2){
    fHistos->FillTH1("hEventCountDG2", 1);
    // Check for exclusive classes
    if(!(isDG1)){
      fHistos->FillTH1("hEventCountDG2excl", 1);
    }
  }

  if(isEG1){
    fHistos->FillTH1("hEventCountEG1", 1);
  }
  if(isDG1){
    fHistos->FillTH1("hEventCountDG1", 1);
  }

  TClonesArray *clusterArray = dynamic_cast<TClonesArray *>(fInputEvent->FindListObject(fClusterContainer.Data()));
  if(!clusterArray){
    AliError("Cluster array not found");
    return;
  }

  Double_t vertexpos[3];
  fInputEvent->GetPrimaryVertex()->GetXYZ(vertexpos);

  Double_t energy, eta, phi;
  for(TIter clustIter = TIter(clusterArray).Begin(); clustIter != TIter::End(); ++clustIter){
    AliVCluster *clust = static_cast<AliVCluster *>(*clustIter);
    if(!clust->IsEMCAL()) continue;

    TLorentzVector posvec;
    energy = clust->E();
    clust->GetMomentum(posvec, vertexpos);
    eta = posvec.Eta();
    phi = posvec.Phi();

    // fill histograms allEta
    if(isMinBias){
      FillClusterHistograms("MB", energy, posvec.Et(), eta, phi, NULL);
      if(!(isEMC7 || isDMC7 || isEJ1 || isEJ2 || isEG1 || isEG2 || isDJ1 || isDJ2 || isDG1 || isDG2)){
        FillClusterHistograms("hMBexcl", energy, posvec.Et(), eta, phi, NULL);
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
  EmcalTriggerClass myclass = kECREL0;
  if(triggerclass == "EG1") myclass = kECREG1;
  if(triggerclass == "EG2") myclass = kECREG2;
  if(triggerclass == "EJ1") myclass = kECREJ1;
  if(triggerclass == "EJ2") myclass = kECREJ2;
  if(triggerclass == "DMC7") myclass = kECRDL0;
  if(triggerclass == "DG1") myclass = kECRDG1;
  if(triggerclass == "DG2") myclass = kECRDG2;
  if(triggerclass == "DJ1") myclass = kECRDJ1;
  if(triggerclass == "DJ2") myclass = kECRDJ2;
  bool isSingleShower = (
      (myclass == kECREG1) || (triggerclass == kECREG2) || (myclass == kECREL0)
      || (myclass == kECRDG1) || (triggerclass == kECRDG2) || (myclass == kECRDL0)
  ),
      isDCAL = ((myclass == kECRDL0) || (myclass == kECRDG1) || (myclass == kECRDG2) || (myclass == kECRDJ1) || (myclass == kECRDJ2));
  AliEmcalTriggerPatchInfoAPV1 *mypatch = NULL;
  for(TIter patchiter = TIter(triggerpatches).Begin(); patchiter != TIter::End(); ++patchiter){
    if(!IsOfflineSimplePatch(*patchiter)) continue;
    if(isDCAL){
      if(!SelectDCALPatch(*patchiter)) continue;
    } else {
      if(SelectDCALPatch(*patchiter)) continue;
    }
    if(isSingleShower){
      if(!SelectSingleShowerPatch(*patchiter)) continue;
    } else {
      if(!SelectJetPatch(*patchiter)) continue;
    }
    if(GetPatchEnergy(*patchiter) > fOfflineEnergyThreshold[myclass]) foundtriggers.Add(mypatch);
  }
}

/**
 * Apply additional cut requiring at least one offline patch above a given energy (not fake ADC!)
 * Attention: This task groups into single shower triggers (L0, EG1, EG2) and jet triggers (EJ1 and EJ2).
 * Per convention the low threshold patch is selected. No energy cut should be applied in the trigger maker
 * @param trgcls Trigger class for which to apply additional offline patch selection
 * @param triggerpatches Array of trigger patches
 * @return True if at least on patch above threshold is found or no cut is applied
 */
Bool_t AliAnalysisTaskEmcalClustersRef::IsOfflineSelected(EmcalTriggerClass trgcls, const TClonesArray * const triggerpatches) const {
  if(fOfflineEnergyThreshold[trgcls] < 0) return true;
  bool isSingleShower = ((trgcls == kECREL0) || (trgcls == kECREG1) || (trgcls == kECREG2) ||
      (trgcls == kECRDL0) || (trgcls == kECRDG1) || (trgcls == kECRDG2)),
      isDCAL = ((trgcls == kECRDL0) || (trgcls == kECRDG1) || (trgcls == kECRDG2) || (trgcls == kECRDJ1) || (trgcls == kECRDJ2));
  int nfound = 0;
  for(TIter patchIter = TIter(triggerpatches).Begin(); patchIter != TIter::End(); ++patchIter){
    bool isDCALpatch = SelectDCALPatch(*patchIter);
    if(isDCAL){
      if(!isDCALpatch) continue;
    } else {
      if(isDCALpatch) continue;
    }
    if(!IsOfflineSimplePatch(*patchIter)) continue;
    if(isSingleShower){
      if(!SelectSingleShowerPatch(*patchIter)) continue;
    } else{
      if(!SelectJetPatch(*patchIter)) continue;
    }
    if(GetPatchEnergy(*patchIter) > fOfflineEnergyThreshold[trgcls]) nfound++;
  }
  return nfound > 0;
}

/**
 * Apply trigger selection using offline patches and trigger thresholds based on offline ADC Amplitude
 * @param triggerpatches Trigger patches found by the trigger maker
 * @return String with EMCAL trigger decision
 */
TString AliAnalysisTaskEmcalClustersRef::GetFiredTriggerClassesFromPatches(const TClonesArray* triggerpatches) const {
  TString triggerstring = "";
  if(!triggerpatches) return triggerstring;
  Int_t nEJ1 = 0, nEJ2 = 0, nEG1 = 0, nEG2 = 0;
  double  minADC_EJ1 = 260.,
          minADC_EJ2 = 127.,
          minADC_EG1 = 140.,
          minADC_EG2 = 89.;
  for(TIter patchIter = TIter(triggerpatches).Begin(); patchIter != TIter::End(); ++patchIter){
    AliEmcalTriggerPatchInfo *patch = dynamic_cast<AliEmcalTriggerPatchInfo *>(*patchIter);
    if(!patch->IsOfflineSimple()) continue;
    if(patch->IsJetHighSimple() && patch->GetADCOfflineAmp() > minADC_EJ1) nEJ1++;
    if(patch->IsJetLowSimple() && patch->GetADCOfflineAmp() > minADC_EJ2) nEJ2++;
    if(patch->IsGammaHighSimple() && patch->GetADCOfflineAmp() > minADC_EG1) nEG1++;
    if(patch->IsGammaLowSimple() && patch->GetADCOfflineAmp() > minADC_EG2) nEG2++;
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
  return triggerstring;
}

void AliAnalysisTaskEmcalClustersRef::GetPatchBoundaries(TObject *o, Double_t *boundaries) const {
  AliEmcalTriggerPatchInfo *patch(NULL);
  if((patch = dynamic_cast<AliEmcalTriggerPatchInfo *>(o))){
    boundaries[0] = patch->GetEtaMin();
    boundaries[1] = patch->GetEtaMax();
    boundaries[2] = patch->GetPhiMin();
    boundaries[3] = patch->GetPhiMax();
  } else {
    AliEmcalTriggerPatchInfoAPV1 *newpatch = dynamic_cast<AliEmcalTriggerPatchInfoAPV1 *>(o);
    if(newpatch){
      boundaries[0] = newpatch->GetEtaMin();
      boundaries[1] = newpatch->GetEtaMax();
      boundaries[2] = newpatch->GetPhiMin();
      boundaries[3] = newpatch->GetPhiMax();
    }
  }

}

bool AliAnalysisTaskEmcalClustersRef::IsOfflineSimplePatch(TObject *o) const {
  AliEmcalTriggerPatchInfo *patch(NULL);
  if((patch = dynamic_cast<AliEmcalTriggerPatchInfo *>(o))){
    return patch->IsOfflineSimple();
  } else {
    AliEmcalTriggerPatchInfoAPV1 *newpatch = dynamic_cast<AliEmcalTriggerPatchInfoAPV1 *>(o);
    if(newpatch) return newpatch->IsOfflineSimple();
  }
  return false;
}

bool AliAnalysisTaskEmcalClustersRef::SelectDCALPatch(TObject *o) const {
  AliEmcalTriggerPatchInfo *patch(NULL);
  if((patch = dynamic_cast<AliEmcalTriggerPatchInfo *>(o))){
    return patch->GetRowStart() >= 64;
  } else {
    AliEmcalTriggerPatchInfoAPV1 *newpatch = dynamic_cast<AliEmcalTriggerPatchInfoAPV1 *>(o);
    if(newpatch) return newpatch->GetRowStart() >= 64;
  }
  return false;
}

bool AliAnalysisTaskEmcalClustersRef::SelectSingleShowerPatch(TObject *o) const{
  AliEmcalTriggerPatchInfo *patch(NULL);
  if((patch = dynamic_cast<AliEmcalTriggerPatchInfo *>(o))){
    if(!patch->IsOfflineSimple()) return false;
    return patch->IsGammaLowSimple();
  } else {
    AliEmcalTriggerPatchInfoAPV1 *newpatch = dynamic_cast<AliEmcalTriggerPatchInfoAPV1 *>(o);
    if(newpatch) return newpatch->IsGammaLowSimple();
  }
  return false;
}

bool AliAnalysisTaskEmcalClustersRef::SelectJetPatch(TObject *o) const{
  AliEmcalTriggerPatchInfo *patch(NULL);
  if((patch = dynamic_cast<AliEmcalTriggerPatchInfo *>(o))){
    if(!patch->IsOfflineSimple()) return false;
    return patch->IsJetLowSimple();
  } else {
    AliEmcalTriggerPatchInfoAPV1 *newpatch = dynamic_cast<AliEmcalTriggerPatchInfoAPV1 *>(o);
    if(newpatch) return newpatch->IsJetLowSimple();
  }
  return false;
}


double AliAnalysisTaskEmcalClustersRef::GetPatchEnergy(TObject *o) const {
  double energy = 0.;
  AliEmcalTriggerPatchInfo *patch(NULL);
  if((patch = dynamic_cast<AliEmcalTriggerPatchInfo *>(o))){
    if(!patch->IsOfflineSimple()) return false;
    energy = patch->GetPatchE();
  } else {
    AliEmcalTriggerPatchInfoAPV1 *newpatch = dynamic_cast<AliEmcalTriggerPatchInfoAPV1 *>(o);
    if(newpatch) energy = newpatch->GetPatchE();
  }
  return energy;
}

} /* namespace EMCalTriggerPtAnalysis */
