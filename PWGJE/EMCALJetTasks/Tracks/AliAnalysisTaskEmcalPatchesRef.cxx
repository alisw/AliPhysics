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
#include <THistManager.h>
#include <THashList.h>
#include <TString.h>

#include "AliAnalysisUtils.h"
#include "AliESDEvent.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliEmcalTriggerOfflineSelection.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliMultSelection.h"
#include "AliMultEstimator.h"

#include "AliAnalysisTaskEmcalPatchesRef.h"

ClassImp(EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalPatchesRef)

namespace EMCalTriggerPtAnalysis {

/**
 * Dummy (I/O) onstructor
 */
AliAnalysisTaskEmcalPatchesRef::AliAnalysisTaskEmcalPatchesRef() :
        AliAnalysisTaskSE(),
        fAnalysisUtil(NULL),
        fTriggerSelection(NULL),
        fHistos(NULL),
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
AliAnalysisTaskEmcalPatchesRef::AliAnalysisTaskEmcalPatchesRef(const char *name):
    AliAnalysisTaskSE(name),
    fAnalysisUtil(NULL),
    fTriggerSelection(NULL),
    fHistos(NULL),
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
AliAnalysisTaskEmcalPatchesRef::~AliAnalysisTaskEmcalPatchesRef() {
}

/**
 * Creating output histograms:
 * + Patch (calibrated) energy spectrum - separated by patch type - for different trigger classes
 * + Patch eta-phi map - separated by patch type - for different trigger classes and different min. energies
 */
void AliAnalysisTaskEmcalPatchesRef::UserCreateOutputObjects(){
  AliInfo(Form("Creating histograms for task %s\n", GetName()));
  fAnalysisUtil = new AliAnalysisUtils;

  TArrayD energybinning, etabinning;
  CreateEnergyBinning(energybinning);
  CreateLinearBinning(etabinning, 100, -0.7, 0.7);
  fHistos = new THistManager("Ref");
  TString triggers[18] = {"MB", "EMC7", "DMC7",
      "EJ1", "EJ2", "EG1", "EG2", "DJ1", "DJ2", "DG1", "DG2",
      "MBexcl", "EMC7excl", "DMC7excl", "EG2excl", "EJ2excl", "DG2excl", "DJ2excl"
  };
  TString patchtype[10] = {"EG1", "EG2", "EJ1", "EJ2", "EMC7", "DG1", "DG2", "DJ1", "DJ2", "DMC7"};
  Double_t encuts[5] = {1., 2., 5., 10., 20.};
  for(TString *trg = triggers; trg < triggers+18; trg++){
    fHistos->CreateTH1(Form("hEventCount%s", trg->Data()), Form("Event count for trigger class %s", trg->Data()), 1, 0.5, 1.5);
    for(int ipatch = 0; ipatch < 10; ipatch++){
      fHistos->CreateTH1(Form("h%sPatchEnergy%s", patchtype[ipatch].Data(), trg->Data()), Form("%s-patch energy for trigger class %s", patchtype[ipatch].Data(), trg->Data()), energybinning);
      fHistos->CreateTH1(Form("h%sPatchET%s", patchtype[ipatch].Data(), trg->Data()), Form("%s-patch transverse energy for trigger class %s", patchtype[ipatch].Data(), trg->Data()), energybinning);
      fHistos->CreateTH2(Form("h%sPatchEnergyEta%s", patchtype[ipatch].Data(), trg->Data()), Form("%s-patch energy for trigger class %s", patchtype[ipatch].Data(), trg->Data()), energybinning, etabinning);
      fHistos->CreateTH2(Form("h%sPatchETEta%s", patchtype[ipatch].Data(), trg->Data()), Form("%s-patch transverse energy for trigger class %s", patchtype[ipatch].Data(), trg->Data()), energybinning, etabinning);
      for(int ien = 0; ien < 5; ien++){
        fHistos->CreateTH2(Form("h%sEtaPhi%dG%s", patchtype[ipatch].Data(), static_cast<int>(encuts[ien]), trg->Data()), Form("%s-patch #eta-#phi map for patches with energy larger than %f GeV/c for trigger class %s", patchtype[ipatch].Data(), encuts[ien], trg->Data()), 100, -0.7, 0.7, 200, 0, TMath::TwoPi());
        fHistos->CreateTH2(Form("h%sColRow%dG%s", patchtype[ipatch].Data(), static_cast<int>(encuts[ien]), trg->Data()), Form("%s-patch col-row map for patches with energy larger than %f GeV/c for trigger class %s", patchtype[ipatch].Data(), encuts[ien], trg->Data()), 48, -0.5, 47.5, 104, -0.5, 103.5);
      }
    }
  }
  PostData(1, fHistos->GetListOfHistograms());
  AliDebug(1, "Histograms done");
}

/**
 * Run event loop
 * @param Not used
 */
void AliAnalysisTaskEmcalPatchesRef::UserExec(Option_t *){
  AliDebug(1, Form("%s: Start function\n", GetName()));
  TClonesArray *patches = dynamic_cast<TClonesArray *>(fInputEvent->FindListObject("EmcalTriggers"));
  TString triggerstring = "";
  if(fTriggerStringFromPatches){
    triggerstring = GetFiredTriggerClassesFromPatches(patches);
  } else {
    triggerstring = fInputEvent->GetFiredTriggerClasses();
  }
  UInt_t selectionstatus = fInputHandler->IsEventSelected();
  Bool_t isMinBias = selectionstatus & AliVEvent::kINT7,
      isEMC7 = (selectionstatus & AliVEvent::kEMC7) && triggerstring.Contains("EMC7"),
      isEJ1 = (selectionstatus & AliVEvent::kEMCEJE) && triggerstring.Contains("EJ1"),
      isEJ2 = (selectionstatus & AliVEvent::kEMCEJE) && triggerstring.Contains("EJ2"),
      isEG1 = (selectionstatus & AliVEvent::kEMCEGA) && triggerstring.Contains("EG1"),
      isEG2 = (selectionstatus & AliVEvent::kEMCEGA) && triggerstring.Contains("EG2"),
      isDMC7 = (selectionstatus & AliVEvent::kEMC7) && triggerstring.Contains("DMC7"),
      isDJ1 = (selectionstatus & AliVEvent::kEMCEJE) && triggerstring.Contains("DJ1"),
      isDJ2 = (selectionstatus & AliVEvent::kEMCEJE) && triggerstring.Contains("DJ2"),
      isDG1 = (selectionstatus & AliVEvent::kEMCEGA) && triggerstring.Contains("DG1"),
      isDG2 = (selectionstatus & AliVEvent::kEMCEGA) && triggerstring.Contains("DG2");
  if(patches && fTriggerSelection){
      isEMC7 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgEL0, patches);
      isEJ1 &=  fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgEJ1, patches);
      isEJ2 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgEJ2, patches);
      isEG1 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgEG1, patches);
      isEG2 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgEG2, patches);
      isDMC7 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgDL0, patches);
      isDJ1 &=  fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgDJ1, patches);
      isDJ2 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgDJ2, patches);
      isDG1 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgDG1, patches);
      isDG2 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgDG2, patches);

  }
  if(!(isMinBias || isEMC7 || isEG1 || isEG2 || isEJ1 || isEJ2 || isDMC7 || isDG1 || isDG2 || isDJ1 || isDJ2)){
    AliDebug(1, Form("%s: Reject trigger\n", GetName()));
    return;
  }
  AliMultSelection *mult = dynamic_cast<AliMultSelection *>(InputEvent()->FindListObject("MultSelection"));
  double centrality =  mult ? mult->GetEstimator("V0M")->GetPercentile() : -1;
  AliDebug(1, Form("%s: Centrality %f\n", GetName(), centrality));
  if(!fCentralityRange.IsInRange(centrality)){
    AliDebug(1, Form("%s: Reject centrality\n", GetName()));
    return;
  }
  const AliVVertex *vtx = fInputEvent->GetPrimaryVertex();
  if(!vtx) vtx = fInputEvent->GetPrimaryVertexSPD();
  //if(!fInputEvent->IsPileupFromSPD(3, 0.8, 3., 2., 5.)) return;         // reject pileup event
  if(vtx->GetNContributors() < 1){
    AliDebug(1, Form("%s: Reject contributor\n", GetName()));
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
    AliDebug(1, Form("%s: Reject Z\n", GetName()));
    return;
  }
  AliDebug(1, Form("%s: Event Selected\n", GetName()));

  // Fill Event counter and reference vertex distributions for the different trigger classes
  if(isMinBias){
    fHistos->FillTH1("hEventCountMB", 1);
    // Check for exclusive classes
    if(!(isEMC7 || isDMC7 || isEJ1 || isEJ2 || isEG1 || isEG2 || isDJ1 || isDJ2 || isDG1 || isDG2)){
      fHistos->FillTH1("hEventCountMBexcl", 1);
    }
  }

  // L0 triggers
  if(isEMC7){
    fHistos->FillTH1("hEventCountEMC7", 1);
    if(!(isEJ1 || isEJ2 || isEG1 || isEG2)){
      fHistos->FillTH1("hEventCountEMC7excl", 1);
    }
  }
  if(isDMC7){
    fHistos->FillTH1("hEventCountDMC7", 1);
    if(!(isDJ1 || isDJ2 || isDG1 || isDG2)){
      fHistos->FillTH1("hEventCountDMC7excl", 1);
    }
  }

  // L1 jet triggers
  if(isEJ1){
    fHistos->FillTH1("hEventCountEJ1", 1);
  }
  if(isDJ1){
    fHistos->FillTH1("hEventCountEJ1", 1);
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

  // L1 gamma triggers
  if(isEG1){
    fHistos->FillTH1("hEventCountEG1", 1);
  }
  if(isDG1){
    fHistos->FillTH1("hEventCountDG1", 1);
  }
  if(isEG2){
    fHistos->FillTH1("hEventCountEG2", 1);
    // Check for exclusive classes
    if(!isEG1){
      fHistos->FillTH1("hEventCountEG2excl", 1);
    }
  }
  if(isDG2){
    fHistos->FillTH1("hEventCountDG2", 1);
    // Check for exclusive classes
    if(!isEG1){
      fHistos->FillTH1("hEventCountDG2excl", 1);
    }
  }

  if(!patches){
    AliError("Trigger patch container not available");
    return;
  }

  AliDebug(1, Form("%s: Number of trigger patches %d\n", GetName(), patches->GetEntries()));

  Double_t vertexpos[3];
  fInputEvent->GetPrimaryVertex()->GetXYZ(vertexpos);

  Double_t energy, eta, phi, et;
  Int_t col, row;
  for(TIter patchIter = TIter(patches).Begin(); patchIter != TIter::End(); ++patchIter){
    if(!IsOfflineSimplePatch(*patchIter)) continue;
    AliEMCALTriggerPatchInfo *patch = static_cast<AliEMCALTriggerPatchInfo *>(*patchIter);

    bool isDCAL         = SelectDCALPatch(*patchIter),
        isSingleShower  = SelectSingleShowerPatch(*patchIter),
        isJetPatch      = SelectJetPatch(*patchIter);

    std::vector<TString> patchnames;
    if(isJetPatch){
      if(isDCAL){
        patchnames.push_back("DJ1");
        patchnames.push_back("DJ2");
      } else {
        patchnames.push_back("EJ1");
        patchnames.push_back("EJ2");
      }
    }
    if(isSingleShower){
      if(isDCAL){
        patchnames.push_back("DMC7");
        patchnames.push_back("DG1");
        patchnames.push_back("DG2");
      } else {
        patchnames.push_back("EMC7");
        patchnames.push_back("EG1");
        patchnames.push_back("EG2");
      }
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
        FillPatchHistograms("MB", *nameit, energy, et, eta, phi, col, row);
        // check for exclusive classes
        if(!(isEMC7 || isDMC7 || isEJ1 || isEJ2 || isEG1 || isEG2 || isDJ1 || isDJ2 || isDG1 || isDG2)){
          FillPatchHistograms("MBexcl", *nameit, energy, et, eta, phi, col, row);
        }
      }

      if(isEMC7){
        FillPatchHistograms("EMC7", *nameit, energy, et, eta, phi, col, row);
        if(!(isEJ1 || isEJ2 || isEG1 || isEG2)){
          FillPatchHistograms("EMC7excl", *nameit, energy, et, eta, phi, col, row);
        }
      }
      if(isDMC7){
        FillPatchHistograms("DMC7", *nameit, energy, et, eta, phi, col, row);
        if(!(isDJ1 || isDJ2 || isDG1 || isDG2)){
          FillPatchHistograms("DMC7excl", *nameit, energy, et, eta, phi, col, row);
        }
      }

      if(isEJ1){
        FillPatchHistograms("EJ1", *nameit, energy, et, eta, phi, col, row);
      }
      if(isDJ1){
        FillPatchHistograms("DJ1", *nameit, energy, et, eta, phi, col, row);
      }


      if(isEJ2){
        FillPatchHistograms("EJ2", *nameit, energy, et, eta, phi, col, row);
        // check for exclusive classes
        if(!isEJ1){
          FillPatchHistograms("EJ2excl", *nameit, energy, et, eta, phi, col, row);
        }
      }
      if(isDJ2){
        FillPatchHistograms("DJ2", *nameit, energy, et, eta, phi, col, row);
        // check for exclusive classes
        if(!isDJ1){
          FillPatchHistograms("DJ2excl", *nameit, energy, et, eta, phi, col, row);
        }
      }

      if(isEG1){
        FillPatchHistograms("EG1", *nameit, energy, et, eta, phi, col, row);
      }
      if(isDG1){
        FillPatchHistograms("DG1", *nameit, energy, et, eta, phi, col, row);
      }

      if(isEG2){
        FillPatchHistograms("EG2", *nameit, energy, et, eta, phi, col, row);
        // check for exclusive classes
        if(!isEG1){
          FillPatchHistograms("EG2excl", *nameit, energy, et, eta, phi, col, row);
        }
      }
      if(isDG2){
        FillPatchHistograms("DG2", *nameit, energy, et, eta, phi, col, row);
        // check for exclusive classes
        if(!isDG1){
          FillPatchHistograms("DG2excl", *nameit, energy, et, eta, phi, col, row);
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
void AliAnalysisTaskEmcalPatchesRef::FillPatchHistograms(TString triggerclass, TString patchname, double energy, double transverseenergy, double eta, double phi, int col, int row){
  fHistos->FillTH1(Form("h%sPatchEnergy%s", patchname.Data(), triggerclass.Data()), energy);
  fHistos->FillTH1(Form("h%sPatchET%s", patchname.Data(), triggerclass.Data()), transverseenergy);
  fHistos->FillTH2(Form("h%sPatchEnergyEta%s", patchname.Data(), triggerclass.Data()), eta, energy);
  fHistos->FillTH2(Form("h%sPatchETEta%s", patchname.Data(), triggerclass.Data()), eta, transverseenergy);
  Double_t encuts[5] = {1., 2., 5., 10., 20.};
  for(int ien = 0; ien < 5; ien++){
    if(energy > encuts[ien]){
      fHistos->FillTH2(Form("h%sEtaPhi%dG%s", patchname.Data(), static_cast<int>(encuts[ien]), triggerclass.Data()), eta, phi);
      fHistos->FillTH2(Form("h%sEtaPhi%dG%s", patchname.Data(), static_cast<int>(encuts[ien]), triggerclass.Data()), col, row);
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
 * Apply trigger selection using offline patches and trigger thresholds based on offline ADC Amplitude
 * @param triggerpatches Trigger patches found by the trigger maker
 * @return String with EMCAL trigger decision
 */
TString AliAnalysisTaskEmcalPatchesRef::GetFiredTriggerClassesFromPatches(const TClonesArray* triggerpatches) const {
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
  if(nEJ2){
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

void AliAnalysisTaskEmcalPatchesRef::GetPatchBoundaries(TObject *o, Double_t *boundaries) const {
  AliEMCALTriggerPatchInfo *patch = dynamic_cast<AliEMCALTriggerPatchInfo *>(o);
  boundaries[0] = patch->GetEtaMin();
  boundaries[1] = patch->GetEtaMax();
  boundaries[2] = patch->GetPhiMin();
  boundaries[3] = patch->GetPhiMax();
}

bool AliAnalysisTaskEmcalPatchesRef::IsOfflineSimplePatch(TObject *o) const {
  AliEMCALTriggerPatchInfo *patch = dynamic_cast<AliEMCALTriggerPatchInfo *>(o);
  return patch->IsOfflineSimple();
}

bool AliAnalysisTaskEmcalPatchesRef::SelectDCALPatch(TObject *o) const {
  AliEMCALTriggerPatchInfo *patch = dynamic_cast<AliEMCALTriggerPatchInfo *>(o);
  return patch->GetRowStart() >= 64;
}

bool AliAnalysisTaskEmcalPatchesRef::SelectSingleShowerPatch(TObject *o) const{
  AliEMCALTriggerPatchInfo *patch = dynamic_cast<AliEMCALTriggerPatchInfo *>(o);
  if(!patch->IsOfflineSimple()) return false;
  return patch->IsGammaLowSimple();
}

bool AliAnalysisTaskEmcalPatchesRef::SelectJetPatch(TObject *o) const{
  AliEMCALTriggerPatchInfo *patch = dynamic_cast<AliEMCALTriggerPatchInfo *>(o);
  if(!patch->IsOfflineSimple()) return false;
  return patch->IsJetLowSimple();
}

double AliAnalysisTaskEmcalPatchesRef::GetPatchEnergy(TObject *o) const {
  double energy = 0.;
  AliEMCALTriggerPatchInfo *patch = dynamic_cast<AliEMCALTriggerPatchInfo *>(o);
  energy = patch->GetPatchE();
  return energy;
}


} /* namespace EMCalTriggerPtAnalysis */
