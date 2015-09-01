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
#include <vector>
#include <map>

#include <TArrayD.h>
#include <TClonesArray.h>
#include <THashList.h>
#include <TLorentzVector.h>
#include <TString.h>

#include "AliAnalysisUtils.h"
#include "AliEmcalTriggerPatchInfo.h"
#include "AliEMCalHistoContainer.h"
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
    fHistos(),
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
  fHistos = new AliEMCalHistoContainer("Ref");
  TString triggers[14] = {"MB", "EJ1", "EJ2", "EG1", "EG2", "EG2excl", "EJ2excl", "MBexcl", "E1combined", "E1Jonly", "E1Gonly", "E2combined", "E2Jonly", "E2Gonly"};
  Double_t encuts[5] = {1., 2., 5., 10., 20.};
  for(TString *trg = triggers; trg < triggers+14; trg++){
    fHistos->CreateTH1(Form("hEventCount%s", trg->Data()), Form("Event count for trigger class %s", trg->Data()), 1, 0.5, 1.5);
    fHistos->CreateTH1(Form("hClusterEnergy%s", trg->Data()), Form("Cluster energy for trigger class %s", trg->Data()), energybinning);
    for(int ien = 0; ien < 5; ien++){
      fHistos->CreateTH2(Form("hEtaPhi%dG%s", static_cast<int>(encuts[ien]), trg->Data()), Form("cluster #eta-#phi map for clusters with energy larger than %f GeV/c for trigger class %s", encuts[ien], trg->Data()), 100, -0.7, 0.7, 100, 1.4, 3.2);
    }
  }
  PostData(1, fHistos->GetListOfHistograms());
}


/**
 *
 * @param
 */
void AliAnalysisTaskEmcalClustersRef::UserExec(Option_t *){
  TString triggerstring = "";
  if(fTriggerStringFromPatches){
    triggerstring = GetFiredTriggerClassesFromPatches(dynamic_cast<TClonesArray *>(fInputEvent->FindListObject("EmcalTriggers")));
  } else {
    triggerstring = fInputEvent->GetFiredTriggerClasses();
  }
  Bool_t isMinBias = fInputHandler->IsEventSelected() & AliVEvent::kINT7,
      isEJ1 = triggerstring.Contains("EJ1"),
      isEJ2 = triggerstring.Contains("EJ2"),
      isEG1 = triggerstring.Contains("EG1"),
      isEG2 = triggerstring.Contains("EG2");
  if(!(isMinBias || isEG1 || isEG2 || isEJ1 || isEJ2)) return;
  const AliVVertex *vtx = fInputEvent->GetPrimaryVertex();
  //if(!fInputEvent->IsPileupFromSPD(3, 0.8, 3., 2., 5.)) return;         // reject pileup event
  if(vtx->GetNContributors() < 1) return;
  // Fill reference distribution for the primary vertex before any z-cut
  if(!fAnalysisUtil->IsVertexSelected2013pA(fInputEvent)) return;       // Apply new vertex cut
  if(fAnalysisUtil->IsPileUpEvent(fInputEvent)) return;       // Apply new vertex cut
  // Apply vertex z cut
  if(vtx->GetZ() < -10. || vtx->GetZ() > 10.) return;

  // Fill Event counter and reference vertex distributions for the different trigger classes
  if(isMinBias){
    fHistos->FillTH1("hEventCountMB", 1);
    // Check for exclusive classes
    if(!(isEG1 || isEG2 || isEJ1 || isEJ2)){
      fHistos->FillTH1("hEventCountMBexcl", 1);
    }
  }
  if(isEJ1){
    fHistos->FillTH1("hEventCountEJ1", 1);
    if(isEG1 || isEG2){
      fHistos->FillTH1("hEventCountE1combined", 1);
    } else {
      fHistos->FillTH1("hEventCountE1Jonly", 1);
    }

  }
  if(isEJ2){
    fHistos->FillTH1("hEventCountEJ2", 1);
    // Check for exclusive classes
    if(!isEJ1){
      fHistos->FillTH1("hEventCountEJ2excl", 1);
    }
    if(isEG1 || isEG2){
      fHistos->FillTH1("hEventCountE2combined", 1);
    } else {
      fHistos->FillTH1("hEventCountE2Jonly", 1);
    }

  }
  if(isEG1){
    fHistos->FillTH1("hEventCountEG1", 1);
  }
  if(isEG2){
    fHistos->FillTH1("hEventCountEG2", 1);
    // Check for exclusive classes
    if(!isEG1){
      fHistos->FillTH1("hEventCountEG2excl", 1);
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
      FillClusterHistograms("MB", energy, eta, phi);
      // check for exclusive classes
      if(!(isEG1 || isEG2 || isEJ1 || isEJ2)){
        FillClusterHistograms("MBexcl", energy, eta, phi);
      }
    }
    if(isEJ1){
      FillClusterHistograms("EJ1", energy, eta, phi);
      if(isEG1 || isEG2) {
        FillClusterHistograms("E1combined", energy, eta, phi);
      } else {
        FillClusterHistograms("E1Jonly", energy, eta, phi);
      }
    }
    if(isEJ2){
      FillClusterHistograms("EJ2", energy, eta, phi);
      // check for exclusive classes
      if(!isEJ1){
        FillClusterHistograms("EJ2excl", energy, eta, phi);
      }
      if(isEG1 || isEG2){
        FillClusterHistograms("E2combined", energy, eta, phi);
      } else {
        FillClusterHistograms("E2Jonly", energy, eta, phi);
      }
    }
    if(isEG1){
      FillClusterHistograms("EG1", energy, eta, phi);
      if(!(isEJ1 || isEJ2)){
        FillClusterHistograms("E1Gonly", energy, eta, phi);
      }
    }
    if(isEG2){
      FillClusterHistograms("EG2", energy, eta, phi);
      // check for exclusive classes
      if(!isEG1){
        FillClusterHistograms("EG2excl", energy, eta, phi);
      }
      if(!(isEJ2 || isEJ1)){
        FillClusterHistograms("E2Gonly", energy, eta, phi);
      }
    }

  }
  PostData(1, fHistos->GetListOfHistograms());
}

void AliAnalysisTaskEmcalClustersRef::FillClusterHistograms(TString triggerclass, double energy, double eta, double phi){
  fHistos->FillTH1(Form("hClusterEnergy%s", triggerclass.Data()), energy);
  Double_t encuts[5] = {1., 2., 5., 10., 20.};
  for(int ien = 0; ien < 5; ien++){
    if(energy > encuts[ien]){
      fHistos->FillTH2(Form("hEtaPhi%dG%s", static_cast<int>(encuts[ien]), triggerclass.Data()), eta, phi);
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
 * Apply trigger selection using offline patches and trigger thresholds based on offline ADC Amplitude
 * @param triggerpatches Trigger patches found by the trigger maker
 * @return String with EMCAL trigger decision
 */
TString AliAnalysisTaskEmcalClustersRef::GetFiredTriggerClassesFromPatches(const TClonesArray* triggerpatches) const {
  TString triggerstring = "";
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
