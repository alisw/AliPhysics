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
#include <map>
#include <vector>

#include <TArrayD.h>
#include <TClonesArray.h>
#include <THashList.h>
#include <TString.h>

#include "AliAnalysisUtils.h"
#include "AliEMCalHistoContainer.h"
#include "AliEmcalTriggerPatchInfoAP.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"

#include "AliAnalysisTaskEmcalPatchesRef.h"

ClassImp(EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalPatchesRef)

namespace EMCalTriggerPtAnalysis {

/**
 * Dummy (I/O) onstructor
 */
AliAnalysisTaskEmcalPatchesRef::AliAnalysisTaskEmcalPatchesRef() :
        AliAnalysisTaskSE(),
        fAnalysisUtil(NULL),
        fHistos(NULL),
        fTriggerStringFromPatches(kFALSE)
{
  for(int itrg = 0; itrg < kEPRntrig; itrg++){
    fOfflineEnergyThreshold[itrg] = -1;
  }
}

/**
 * Named constructor
 * @param name Name of the task
 */
AliAnalysisTaskEmcalPatchesRef::AliAnalysisTaskEmcalPatchesRef(const char *name):
    AliAnalysisTaskSE(name),
    fAnalysisUtil(NULL),
    fHistos(NULL),
    fTriggerStringFromPatches(kFALSE)
{
  for(int itrg = 0; itrg < kEPRntrig; itrg++){
    fOfflineEnergyThreshold[itrg] = -1;
  }
  DefineOutput(1, TList::Class());
}

/**
 * Destructor
 */
AliAnalysisTaskEmcalPatchesRef::~AliAnalysisTaskEmcalPatchesRef() {
}

/**
 * Creating output histograms:
 * + Patch (calibrated) energy spectrum - separated by patch type - for different trigger classes
 * + Patch eta-phi map - separated by patch type - for different trigger classes and different min. energies
 */
void AliAnalysisTaskEmcalPatchesRef::UserCreateOutputObjects(){
  fAnalysisUtil = new AliAnalysisUtils;

  TArrayD energybinning, etabinning;
  CreateEnergyBinning(energybinning);
  CreateLinearBinning(etabinning, 100, -0.7, 0.7);
  fHistos = new AliEMCalHistoContainer("Ref");
  TString triggers[15] = {"MB", "EMC7", "EJ1", "EJ2", "EG1", "EG2", "EG2excl", "EJ2excl", "MBexcl", "E1combined", "E1Jonly", "E1Gonly", "E2combined", "E2Jonly", "E2Gonly"};
  TString patchtype[5] = {"EG1", "EG2", "EJ1", "EJ2", "EMC7"};
  Double_t encuts[5] = {1., 2., 5., 10., 20.};
  for(TString *trg = triggers; trg < triggers+15; trg++){
    fHistos->CreateTH1(Form("hEventCount%s", trg->Data()), Form("Event count for trigger class %s", trg->Data()), 1, 0.5, 1.5);
    for(int ipatch = 0; ipatch < 5; ipatch++){
      fHistos->CreateTH1(Form("h%sPatchEnergy%s", patchtype[ipatch].Data(), trg->Data()), Form("%s-patch energy for trigger class %s", patchtype[ipatch].Data(), trg->Data()), energybinning);
      fHistos->CreateTH1(Form("h%sPatchET%s", patchtype[ipatch].Data(), trg->Data()), Form("%s-patch transverse energy for trigger class %s", patchtype[ipatch].Data(), trg->Data()), energybinning);
      fHistos->CreateTH2(Form("h%sPatchEnergyEta%s", patchtype[ipatch].Data(), trg->Data()), Form("%s-patch energy for trigger class %s", patchtype[ipatch].Data(), trg->Data()), energybinning, etabinning);
      fHistos->CreateTH2(Form("h%sPatchETEta%s", patchtype[ipatch].Data(), trg->Data()), Form("%s-patch transverse energy for trigger class %s", patchtype[ipatch].Data(), trg->Data()), energybinning, etabinning);
      for(int ien = 0; ien < 5; ien++){
        fHistos->CreateTH2(Form("h%sEtaPhi%dG%s", patchtype[ipatch].Data(), static_cast<int>(encuts[ien]), trg->Data()), Form("%s-patch #eta-#phi map for clusters with energy larger than %f GeV/c for trigger class %s", patchtype[ipatch].Data(), encuts[ien], trg->Data()), 100, -0.7, 0.7, 100, 1.4, 3.2);
      }
    }
  }
  PostData(1, fHistos->GetListOfHistograms());

}

/**
 * Run event loop
 * @param Not used
 */
void AliAnalysisTaskEmcalPatchesRef::UserExec(Option_t *){
  TClonesArray *patches = dynamic_cast<TClonesArray *>(fInputEvent->FindListObject("EmcalTriggers"));
  TString triggerstring = "";
  if(fTriggerStringFromPatches){
    triggerstring = GetFiredTriggerClassesFromPatches(patches);
  } else {
    triggerstring = fInputEvent->GetFiredTriggerClasses();
  }
  UInt_t selectionstatus = fInputHandler->IsEventSelected();
  Bool_t isMinBias = selectionstatus & AliVEvent::kINT7,
      isEMC7 = (selectionstatus & AliVEvent::kEMC7) && triggerstring.Contains("EMC7") && IsOfflineSelected(kEPREL0, patches),
      isEJ1 = (selectionstatus & AliVEvent::kEMCEJE) && triggerstring.Contains("EJ1") && IsOfflineSelected(kEPREJ1, patches),
      isEJ2 = (selectionstatus & AliVEvent::kEMCEJE) && triggerstring.Contains("EJ2") && IsOfflineSelected(kEPREJ2, patches),
      isEG1 = (selectionstatus & AliVEvent::kEMCEGA) && triggerstring.Contains("EG1") && IsOfflineSelected(kEPREG1, patches),
      isEG2 = (selectionstatus & AliVEvent::kEMCEGA) && triggerstring.Contains("EG2") && IsOfflineSelected(kEPREG2, patches);
  if(!(isMinBias || isEMC7 || isEG1 || isEG2 || isEJ1 || isEJ2)) return;
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
  if(isEMC7){
    fHistos->FillTH1("hEventCountEMC7", 1);
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

  if(!patches){
    AliError("Trigger patch container not available");
    return;
  }

  Double_t vertexpos[3];
  fInputEvent->GetPrimaryVertex()->GetXYZ(vertexpos);

  Double_t energy, eta, phi, et;
  for(TIter patchIter = TIter(patches).Begin(); patchIter != TIter::End(); ++patchIter){
    AliEmcalTriggerPatchInfo *patch = static_cast<AliEmcalTriggerPatchInfo *>(*patchIter);
    if(!patch->IsOfflineSimple()) continue;

    std::vector<TString> patchnames;
    if(patch->IsJetHighSimple())
      patchnames.push_back("EJ1");
    if(patch->IsJetLowSimple())
      patchnames.push_back("EJ2");
    if(patch->IsGammaHighSimple())
      patchnames.push_back("EG1");
    if(patch->IsGammaLowSimple()){
      patchnames.push_back("EG2");
      patchnames.push_back("EMC7");   // Treat EG2 patch also as EMC7 patch (single shower);
    }
    if(!patchnames.size()){
      // Undefined patch type - ignore
      continue;
    }

    TLorentzVector posvec;
    energy = patch->GetPatchE();
    eta = patch->GetEtaGeo();
    phi = patch->GetPhiGeo();
    et = patch->GetLorentzVectorCenterGeo().Et();

    // fill histograms allEta
    for(std::vector<TString>::iterator nameit = patchnames.begin(); nameit != patchnames.end(); ++nameit){
      if(isMinBias){
        FillPatchHistograms("MB", *nameit, energy, et, eta, phi);
        // check for exclusive classes
        if(!(isEG1 || isEG2 || isEJ1 || isEJ2)){
          FillPatchHistograms("MBexcl", *nameit, energy, et, eta, phi);
        }
      }
      if(isEMC7){
        FillPatchHistograms("EMC7", *nameit, energy, et, eta, phi);
      }
      if(isEJ1){
        FillPatchHistograms("EJ1", *nameit, energy, et, eta, phi);
        if(isEG1 || isEG2) {
          FillPatchHistograms("E1combined", *nameit, energy, et, eta, phi);
        } else {
          FillPatchHistograms("E1Jonly", *nameit, energy, et, eta, phi);
        }
      }
      if(isEJ2){
        FillPatchHistograms("EJ2", *nameit, energy, et, eta, phi);
        // check for exclusive classes
        if(!isEJ1){
          FillPatchHistograms("EJ2excl", *nameit, energy, et, eta, phi);
        }
        if(isEG1 || isEG2){
          FillPatchHistograms("E2combined", *nameit, energy, et, eta, phi);
        } else {
          FillPatchHistograms("E2Jonly", *nameit, energy, et, eta, phi);
        }
      }
      if(isEG1){
        FillPatchHistograms("EG1", *nameit, energy, et, eta, phi);
        if(!(isEJ1 || isEJ2)){
          FillPatchHistograms("E1Gonly", *nameit, energy, et, eta, phi);
        }
      }
      if(isEG2){
        FillPatchHistograms("EG2", *nameit, energy, et, eta, phi);
        // check for exclusive classes
        if(!isEG1){
          FillPatchHistograms("EG2excl", *nameit, energy, et, eta, phi);
        }
        if(!(isEJ2 || isEJ1)){
          FillPatchHistograms("E2Gonly", *nameit, energy, et, eta, phi);
        }
      }
    }
  }
  PostData(1, fHistos->GetListOfHistograms());
}

/**
 * Filling patch related histogram
 * @param triggerclass Name of the trigger class firing the event
 * @param patchname Name of the patchtype
 * @param energy Calibrated energy of the patch
 * @param eta Patch eta at the geometrical center
 * @param phi Patch phi at the geometrical center
 */
void AliAnalysisTaskEmcalPatchesRef::FillPatchHistograms(TString triggerclass, TString patchname, double energy, double transverseenergy, double eta, double phi){
  fHistos->FillTH1(Form("h%sPatchEnergy%s", patchname.Data(), triggerclass.Data()), energy);
  fHistos->FillTH1(Form("h%sPatchET%s", patchname.Data(), triggerclass.Data()), transverseenergy);
  fHistos->FillTH2(Form("h%sPatchEnergyEta%s", patchname.Data(), triggerclass.Data()), eta, energy);
  fHistos->FillTH2(Form("h%sPatchETEta%s", patchname.Data(), triggerclass.Data()), eta, transverseenergy);
  Double_t encuts[5] = {1., 2., 5., 10., 20.};
  for(int ien = 0; ien < 5; ien++){
    if(energy > encuts[ien]){
      fHistos->FillTH2(Form("h%sEtaPhi%dG%s", patchname.Data(), static_cast<int>(encuts[ien]), triggerclass.Data()), eta, phi);
    }
  }
}

/**
 * Create new energy binning
 * @param binning
 */
void AliAnalysisTaskEmcalPatchesRef::CreateEnergyBinning(TArrayD& binning) const {
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
void AliAnalysisTaskEmcalPatchesRef::CreateLinearBinning(TArrayD& binning, int nbins, double min, double max) const {
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
 * Apply additional cut requiring at least one offline patch above a given energy (not fake ADC!)
 * Attention: This task groups into single shower triggers (L0, EG1, EG2) and jet triggers (EJ1 and EJ2).
 * Per convention the low threshold patch is selected. No energy cut should be applied in the trigger maker
 * @param trgcls Trigger class for which to apply additional offline patch selection
 * @param triggerpatches Array of trigger patches
 * @return True if at least on patch above threshold is found or no cut is applied
 */
Bool_t AliAnalysisTaskEmcalPatchesRef::IsOfflineSelected(EmcalTriggerClass trgcls, const TClonesArray * const triggerpatches) const {
  if(fOfflineEnergyThreshold[trgcls] < 0) return true;
  bool isSingleShower = ((trgcls == kEPREL0) || (trgcls == kEPREG1) || (trgcls == kEPREG2));
  int nfound = 0;
  AliEmcalTriggerPatchInfo *patch = NULL;
  for(TIter patchIter = TIter(triggerpatches).Begin(); patchIter != TIter::End(); ++patchIter){
    patch = static_cast<AliEmcalTriggerPatchInfo *>(*patchIter);
    if(!patch->IsOfflineSimple()) continue;
    if(isSingleShower){
     if(!patch->IsGammaLowSimple()) continue;
    } else {
      if(!patch->IsJetLowSimple()) continue;
    }
    if(patch->GetPatchE() > fOfflineEnergyThreshold[trgcls]) nfound++;
  }
  return nfound > 0;
}


/**
 * Apply trigger selection using offline patches and trigger thresholds based on offline ADC Amplitude
 * @param triggerpatches Trigger patches found by the trigger maker
 * @return String with EMCAL trigger decision
 */
TString AliAnalysisTaskEmcalPatchesRef::GetFiredTriggerClassesFromPatches(const TClonesArray* triggerpatches) const {
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
