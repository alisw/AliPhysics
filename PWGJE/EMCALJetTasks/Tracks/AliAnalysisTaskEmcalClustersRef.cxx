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
#include <TString.h>

#include "AliAnalysisUtils.h"
#include "AliEMCALGeometry.h"
#include "AliEmcalTriggerPatchInfo.h"
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
    fTriggerStringFromPatches(kFALSE)
{

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
    fTriggerStringFromPatches(kFALSE)
{
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
  TArrayD smbinning(14); CreateLinearBinning(smbinning, 14, -0.5, 13.5);
  TArrayD etabinning; CreateLinearBinning(etabinning, 100, -0.7, 0.7);
  fHistos = new AliEMCalHistoContainer("Ref");
  TString triggers[17] = {
      "MB", "EMC7",
      "EJ1", "EJ2", "EG1", "EG2",
      "EMC7excl","EG1excl", "EG2excl", "EJ1excl", "EJ2excl",
      "E1combined", "E1Jonly", "E1Gonly", "E2combined", "E2Jonly", "E2Gonly"
  };
  Double_t encuts[5] = {1., 2., 5., 10., 20.};
  for(TString *trg = triggers; trg < triggers + sizeof(triggers)/sizeof(TString); trg++){
    fHistos->CreateTH1(Form("hEventCount%s", trg->Data()), Form("Event count for trigger class %s", trg->Data()), 1, 0.5, 1.5);
    fHistos->CreateTH1(Form("hClusterEnergy%s", trg->Data()), Form("Cluster energy for trigger class %s", trg->Data()), energybinning);
    fHistos->CreateTH1(Form("hClusterEnergyFired%s", trg->Data()), Form("Cluster energy for trigger class %s, firing the trigger", trg->Data()), energybinning);
    fHistos->CreateTH2(Form("hClusterEnergySM%s", trg->Data()), Form("Cluster energy versus supermodule for trigger class %s", trg->Data()), smbinning, energybinning);
    fHistos->CreateTH2(Form("hClusterEnergyFiredSM%s", trg->Data()), Form("Cluster energy versus supermodule for trigger class %s, firing the trigger", trg->Data()), smbinning, energybinning);
    fHistos->CreateTH2(Form("hEtaEnergy%s", trg->Data()), Form("Cluster energy vs. eta for trigger class %s", trg->Data()), etabinning, energybinning);
    fHistos->CreateTH2(Form("hEtaEnergyFired%s", trg->Data()), Form("Cluster energy vs. eta for trigger class %s, firing the trigger", trg->Data()), etabinning, energybinning);
    for(int ism = 0; ism < 12; ism++){
      fHistos->CreateTH2(Form("hEtaEnergySM%d%s", ism, trg->Data()), Form("Cluster energy vs. eta in Supermodule %d for trigger %s", ism, trg->Data()), etabinning, energybinning);
      fHistos->CreateTH2(Form("hEtaEnergyFiredSM%d%s", ism, trg->Data()), Form("Cluster energy vs. eta in Supermodule %d for trigger %s, firing the trigger", ism, trg->Data()), etabinning, energybinning);
    }
    for(int isec = 4; isec < 10; isec++){
      fHistos->CreateTH2(Form("hEtaEnergySec%d%s", isec, trg->Data()), Form("Cluster energy vs.eta in tracking sector %d for trigger %s", isec, trg->Data()), etabinning, energybinning);
      fHistos->CreateTH2(Form("hEtaEnergyFiredSec%d%s", isec, trg->Data()), Form("Cluster energy vs.eta in tracking sector %d for trigger %s, firing the trigger", isec, trg->Data()), etabinning, energybinning);
    }
    for(int ien = 0; ien < 5; ien++){
      fHistos->CreateTH2(Form("hEtaPhi%dG%s", static_cast<int>(encuts[ien]), trg->Data()), Form("cluster #eta-#phi map for clusters with energy larger than %f GeV/c for trigger class %s", encuts[ien], trg->Data()), 100, -0.7, 0.7, 100, 1.4, 3.2);
      fHistos->CreateTH2(Form("hEtaPhiFired%dG%s", static_cast<int>(encuts[ien]), trg->Data()), Form("cluster #eta-#phi map for clusters fired the trigger with energy larger than %f GeV/c for trigger class %s", encuts[ien], trg->Data()), 100, -0.7, 0.7, 100, 1.4, 3.2);
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
      isEMC7 = (selectionstatus & AliVEvent::kEMC7) && triggerstring.Contains("CEMC7");
  if(!(isMinBias || isEMC7 || isEG1 || isEG2 || isEJ1 || isEJ2)) return;
  const AliVVertex *vtx = fInputEvent->GetPrimaryVertex();
  //if(!fInputEvent->IsPileupFromSPD(3, 0.8, 3., 2., 5.)) return;         // reject pileup event
  if(vtx->GetNContributors() < 1) return;
  // Fill reference distribution for the primary vertex before any z-cut
  if(fInputEvent->IsA() == AliESDEvent::Class() && fAnalysisUtil->IsFirstEventInChunk(fInputEvent)) return;
  if(!fAnalysisUtil->IsVertexSelected2013pA(fInputEvent)) return;       // Apply new vertex cut
  if(fAnalysisUtil->IsPileUpEvent(fInputEvent)) return;       // Apply new vertex cut
  // Apply vertex z cut
  if(vtx->GetZ() < -10. || vtx->GetZ() > 10.) return;

  // Fill Event counter and reference vertex distributions for the different trigger classes
  if(isMinBias){
    fHistos->FillTH1("hEventCountMB", 1);
  }
  if(isEMC7){
    fHistos->FillTH1("hEventCountEMC7", 1);
    if(!isMinBias){
      fHistos->FillTH1("hEventCountEMC7excl", 1);
    }
  }
  if(isEJ2){
    fHistos->FillTH1("hEventCountEJ2", 1);
    // Check for exclusive classes
    if(!(isMinBias)){
      fHistos->FillTH1("hEventCountEJ2excl", 1);
    }
    if(isEG1 || isEG2){
      fHistos->FillTH1("hEventCountE2combined", 1);
    } else {
      fHistos->FillTH1("hEventCountE2Jonly", 1);
    }

  }
  if(isEJ1){
    fHistos->FillTH1("hEventCountEJ1", 1);
    // Check for exclusive classes
    if(!(isMinBias || isEJ2)){
      fHistos->FillTH1("hEventCountEJ1excl", 1);
    }
    if(isEG1 || isEG2){
      fHistos->FillTH1("hEventCountE1combined", 1);
    } else {
      fHistos->FillTH1("hEventCountE1Jonly", 1);
    }
  }
  if(isEG2){
    fHistos->FillTH1("hEventCountEG2", 1);
    // Check for exclusive classes
    if(!(isMinBias)){
      fHistos->FillTH1("hEventCountEG2excl", 1);
    }
    if(!(isEJ1 || isEJ2)){
      fHistos->FillTH1("hEventCountE2Gonly", 1);
    }
  }
  if(isEG1){
    fHistos->FillTH1("hEventCountEG1", 1);
    // Check for exclusive classes
    if(!(isMinBias || isEG2)){
      fHistos->FillTH1("hEventCountEG1excl", 1);
    }
    if(!(isEJ1 || isEJ2)){
      fHistos->FillTH1("hEventCountE1Gonly", 1);
    }
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
      FillClusterHistograms("MB", energy, eta, phi, NULL);
    }
    if(isEMC7){
      FillClusterHistograms("EMC7", energy, eta, phi, NULL);
      if(!isMinBias){
        FillClusterHistograms("EMC7excl", energy, eta, phi, NULL);
      }
    }
    if(isEJ2){
      TList ej2patches;
      FindPatchesForTrigger("EJ2", triggerpatches, ej2patches);
      FillClusterHistograms("EJ2", energy, eta, phi, &ej2patches);
      // check for exclusive classes
      if(!isMinBias){
        FillClusterHistograms("EJ2excl", energy, eta, phi, &ej2patches);
      }
      if(isEG1 || isEG2){
        FillClusterHistograms("E2combined", energy, eta, phi, &ej2patches);
      } else {
        FillClusterHistograms("E2Jonly", energy, eta, phi, &ej2patches);
      }
    }
    if(isEJ1){
      TList ej1patches;
      FindPatchesForTrigger("EJ1", triggerpatches, ej1patches);
      FillClusterHistograms("EJ1", energy, eta, phi, &ej1patches);
      // check for exclusive classes
      if(!(isMinBias || isEJ2)){
        FillClusterHistograms("EJ1excl", energy, eta, phi, &ej1patches);
      }
      if(isEG1 || isEG2) {
        FillClusterHistograms("E1combined", energy, eta, phi, &ej1patches);
      } else {
        FillClusterHistograms("E1Jonly", energy, eta, phi, &ej1patches);
      }
    }
    if(isEG2){
      TList eg2patches;
      FindPatchesForTrigger("EG2", triggerpatches, eg2patches);
      FillClusterHistograms("EG2", energy, eta, phi, &eg2patches);
      // check for exclusive classes
      if(!isMinBias){
        FillClusterHistograms("EG2excl", energy, eta, phi, &eg2patches);
      }
      if(!(isEJ2 || isEJ1)){
        FillClusterHistograms("E2Gonly", energy, eta, phi, &eg2patches);
      }
    }
    if(isEG1){
      TList eg1patches;
      FindPatchesForTrigger("EG1", triggerpatches, eg1patches);
      FillClusterHistograms("EG1", energy, eta, phi, &eg1patches);
      if(!(isEG2 || isMinBias))
        FillClusterHistograms("EG1excl", energy, eta, phi, &eg1patches);
      if(!(isEJ1 || isEJ2)){
        FillClusterHistograms("E1Gonly", energy, eta, phi, &eg1patches);
      }
    }
  }
  PostData(1, fHistos->GetListOfHistograms());
}

void AliAnalysisTaskEmcalClustersRef::FillClusterHistograms(TString triggerclass, double energy, double eta, double phi, TList *triggerpatches){
  Bool_t hasTriggerPatch = triggerpatches  ? CorrelateToTrigger(eta, phi, triggerpatches) : kFALSE;
  Int_t supermoduleID = -1, sector = -1;
  fGeometry->SuperModuleNumberFromEtaPhi(eta, phi, supermoduleID);
  fHistos->FillTH1(Form("hClusterEnergy%s", triggerclass.Data()), energy);
  fHistos->FillTH2(Form("hEtaEnergy%s", triggerclass.Data()), eta, energy);
  if(supermoduleID >= 0){
    fHistos->FillTH2(Form("hClusterEnergySM%s", triggerclass.Data()), supermoduleID, energy);
    fHistos->FillTH2(Form("hEtaEnergySM%d%s", supermoduleID, triggerclass.Data()), eta, energy);
    sector = 4 + int(supermoduleID/2);
    fHistos->FillTH2(Form("hEtaEnergySec%d%s", sector, triggerclass.Data()), eta, energy);
  }
  if(hasTriggerPatch){
    fHistos->FillTH1(Form("hClusterEnergyFired%s", triggerclass.Data()), energy);
    fHistos->FillTH2(Form("hEtaEnergyFired%s", triggerclass.Data()), eta, energy);
    if(supermoduleID >= 0){
      fHistos->FillTH2(Form("hClusterEnergyFiredSM%s", triggerclass.Data()), supermoduleID, energy);
      fHistos->FillTH2(Form("hEtaEnergyFiredSM%d%s", supermoduleID, triggerclass.Data()), eta, energy);
      fHistos->FillTH2(Form("hEtaEnergyFiredSec%d%s", sector, triggerclass.Data()), eta, energy);
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
    AliEmcalTriggerPatchInfo *mypatch = static_cast<AliEmcalTriggerPatchInfo *>(*patchIter);
    Double_t etamin = TMath::Min(mypatch->GetEtaMin(), mypatch->GetEtaMax()),
        etamax = TMath::Max(mypatch->GetEtaMin(), mypatch->GetEtaMax()),
        phimin = TMath::Min(mypatch->GetPhiMin(), mypatch->GetPhiMax()),
        phimax = TMath::Max(mypatch->GetPhiMin(), mypatch->GetPhiMax());
    if(etaclust > etamin && etaclust < etamax && phiclust > phimin && phiclust < phimax){
      hasfound = kTRUE;
      break;
    }
  }
  return hasfound;
}

/**
 * Find all patches in an event which could have fired the trigger
 * @param triggerclass EMCAL trigger class firing
 * @param triggerpatches Trigger patches found in the event
 * @return List of patches which could have fired the trigger
 */
void AliAnalysisTaskEmcalClustersRef::FindPatchesForTrigger(TString triggerclass, const TClonesArray * triggerpatches, TList &foundtriggers) const {
  double  minADC_EJ1 = 260.,
          minADC_EJ2 = 127.,
          minADC_EG1 = 140.,
          minADC_EG2 = 89.;
  if(!triggerpatches) return;
  Bool_t isEG1 = (triggerclass == "EG1"),
      isEG2 = (triggerclass == "EG2"),
      isEJ1 = (triggerclass == "EJ1"),
      isEJ2 = (triggerclass == "EJ2");
  AliEmcalTriggerPatchInfo *mypatch = NULL;
  for(TIter patchiter = TIter(triggerpatches).Begin(); patchiter != TIter::End(); ++patchiter){
    mypatch = dynamic_cast<AliEmcalTriggerPatchInfo *>(*patchiter);
    if(!mypatch->IsOfflineSimple()) continue;
    if(isEG1){
      if(mypatch->IsGammaHighSimple() && mypatch->GetADCAmp() > minADC_EG1)
        foundtriggers.Add(mypatch);
    }
    if(isEG2){
      if(mypatch->IsGammaLowSimple() && mypatch->GetADCAmp() > minADC_EG2)
        foundtriggers.Add(mypatch);
    }
    if(isEJ1){
      if(mypatch->IsJetHighSimple() && mypatch->GetADCAmp() > minADC_EJ1)
        foundtriggers.Add(mypatch);
    }
    if(isEJ2){
      if(mypatch->IsJetLowSimple() && mypatch->GetADCAmp() > minADC_EJ2)
        foundtriggers.Add(mypatch);
    }
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


} /* namespace EMCalTriggerPtAnalysis */
