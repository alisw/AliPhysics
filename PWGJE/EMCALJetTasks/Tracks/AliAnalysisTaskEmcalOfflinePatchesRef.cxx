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
#include <TClonesArray.h>
#include <THashList.h>
#include <THistManager.h>
#include <TString.h>

#include "AliAnalysisUtils.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliESDEvent.h"
#include "AliInputEventHandler.h"
#include "AliVEvent.h"
#include "AliVVertex.h"

#include "AliAnalysisTaskEmcalOfflinePatchesRef.h"

#include <array>

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalOfflinePatchesRef)

using namespace PWGJE::EMCALJetTasks;

AliAnalysisTaskEmcalOfflinePatchesRef::AliAnalysisTaskEmcalOfflinePatchesRef():
  AliAnalysisTaskSE(),
  fAnalysisUtil(nullptr),
  fGeometry(nullptr),
  fHistos(nullptr)
{
}

AliAnalysisTaskEmcalOfflinePatchesRef::AliAnalysisTaskEmcalOfflinePatchesRef(const char *name):
  AliAnalysisTaskSE(name),
  fAnalysisUtil(nullptr),
  fGeometry(nullptr),
  fHistos(nullptr)
{
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskEmcalOfflinePatchesRef::~AliAnalysisTaskEmcalOfflinePatchesRef() {
}

void AliAnalysisTaskEmcalOfflinePatchesRef::UserCreateOutputObjects(){
  fAnalysisUtil = new AliAnalysisUtils;

  fHistos = new THistManager("Ref");
  TString triggername;
  // Plots at global level:
  // Energy vs. supermodule
  // Energy vs. eta (all sectors)
  // Energy vs. eta for sector
  std::array<TString, 3> patchnames = {
      "EL0", "EG1", "EG2"
  };
  for(auto mytrg : patchnames){
    triggername = mytrg;
    fHistos->CreateTH1(Form("hEventCount%s", triggername.Data()), Form("Event counter for trigger type %s", triggername.Data()), 1, 0.5, 1.5);
    fHistos->CreateTH2(Form("hPatchEnergy%s", triggername.Data()), Form("Patch energy versus supermodule for trigger %s", triggername.Data()), 12, -0.5, 11.5, 200, 0., 200.);
    fHistos->CreateTH2(Form("hPatchET%s", triggername.Data()), Form("Patch transverse energy versus supermodule for trigger %s", triggername.Data()), 12, -0.5, 11.5, 200, 0., 200.);
    fHistos->CreateTH2(Form("hPatchADC%s", triggername.Data()), Form("Patch online ADC versus supermodule for trigger %s", triggername.Data()), 12, -0.5, 11.5, 2100, 0., 2100.);
    fHistos->CreateTH2(Form("hPatchEnergyEta%s", triggername.Data()), Form("Patch energy versus eta for trigger %s", triggername.Data()), 100, -0.7, 0.7, 200., 0., 200.);
    fHistos->CreateTH2(Form("hPatchETEta%s", triggername.Data()), Form("Patch transverse energy versus eta for trigger %s", triggername.Data()), 100, -0.7, 0.7, 200., 0., 200.);
    fHistos->CreateTH2(Form("hPatchADCEta%s", triggername.Data()), Form("Patch online ADC versus eta for trigger %s", triggername.Data()), 100, -0.7, 0.7, 2100., 0., 2100.);
    fHistos->CreateTH2(Form("hEvSelPatchEnergy%s", triggername.Data()), Form("Patch energy versus supermodule for trigger %s in selected events", triggername.Data()), 12, -0.5, 11.5, 200, 0., 200.);
    fHistos->CreateTH2(Form("hEvSelPatchET%s", triggername.Data()), Form("Patch transverse energy versus supermodule for trigger %s in selected events", triggername.Data()), 12, -0.5, 11.5, 200, 0., 200.);
    fHistos->CreateTH2(Form("hEvSelPatchADC%s", triggername.Data()), Form("Patch online ADC versus supermodule for trigger %s in selected events", triggername.Data()), 12, -0.5, 11.5, 2100, 0., 2100.);
    fHistos->CreateTH2(Form("hEvSelPatchEnergyEta%s", triggername.Data()), Form("Patch energy versus eta for trigger %s in selected events", triggername.Data()), 100, -0.7, 0.7, 200., 0., 200.);
    fHistos->CreateTH2(Form("hEvSelPatchETEta%s", triggername.Data()), Form("Patch transverse energy versus eta for trigger %s in selected events", triggername.Data()), 100, -0.7, 0.7, 200., 0., 200.);
    fHistos->CreateTH2(Form("hEvSelPatchADCEta%s", triggername.Data()), Form("Patch online ADC versus eta for trigger %s  in selected events", triggername.Data()), 100, -0.7, 0.7, 2100., 0., 2100.);
    for(int ism = 0; ism <= 9; ism++){ // small sectors do not yet contribute to trigger decision, thus they are in here for the future
      fHistos->CreateTH2(Form("hPatchEnergyEta%sSM%d", triggername.Data(), ism), Form("Patch energy versus eta for trigger %s, Supermodule %d", triggername.Data(), ism), 100, -0.7, 0.7, 200., 0., 200.);
      fHistos->CreateTH2(Form("hPatchETEta%sSM%d", triggername.Data(), ism), Form("Patch transverse energy versus eta for trigger %s, Supermodule %d", triggername.Data(), ism), 100, -0.7, 0.7, 200., 0., 200.);
      fHistos->CreateTH2(Form("hPatchADCEta%sSM%d", triggername.Data(), ism), Form("Patch online ADC versus eta for trigger %s, Supermodule %d", triggername.Data(), ism), 100, -0.7, 0.7, 2100., 0., 2100.);
      fHistos->CreateTH2(Form("hEvSelPatchEnergyEta%sSM%d", triggername.Data(), ism), Form("Patch energy versus eta for trigger %s, Supermodule %d in selected events", triggername.Data(), ism), 100, -0.7, 0.7, 200., 0., 200.);
      fHistos->CreateTH2(Form("hEvSelPatchETEta%sSM%d", triggername.Data(), ism), Form("Patch transverse energy versus eta for trigger %s, Supermodule %d in selected events", triggername.Data(), ism), 100, -0.7, 0.7, 200., 0., 200.);
      fHistos->CreateTH2(Form("hEvSelPatchADCEta%sSM%d", triggername.Data(), ism), Form("Patch online ADC versus eta for trigger %s, Supermodule %d in selected events", triggername.Data(), ism), 100, -0.7, 0.7, 2100., 0., 2100.);
    }
    for(int isec = 4; isec <= 9; isec++){ // small sectors do not yet contribute to trigger decision, thus they are in here for the future
      fHistos->CreateTH2(Form("hPatchEnergyEta%sSector%d", triggername.Data(), isec), Form("Patch energy versus eta for trigger %s Sector %d", triggername.Data(), isec), 100, -0.7, 0.7, 200., 0., 200.);
      fHistos->CreateTH2(Form("hPatchETEta%sSector%d", triggername.Data(), isec), Form("Patch transverse energy versus eta for trigger %s Sector %d", triggername.Data(), isec), 100, -0.7, 0.7, 200., 0., 200.);
      fHistos->CreateTH2(Form("hPatchADCEta%sSector%d", triggername.Data(), isec), Form("Patch online ADC versus eta for trigger %s, Sector %d", triggername.Data(), isec), 100, -0.7, 0.7, 2100., 0., 2100.);
      fHistos->CreateTH2(Form("hEvSelPatchEnergyEta%sSector%d", triggername.Data(), isec), Form("Patch energy versus eta for trigger %s Sector %d in selectedEvents", triggername.Data(), isec), 100, -0.7, 0.7, 200., 0., 200.);
      fHistos->CreateTH2(Form("hEvSelPatchETEta%sSector%d", triggername.Data(), isec), Form("Patch transverse energy versus eta for trigger %s Sector %d in selectedEvents", triggername.Data(), isec), 100, -0.7, 0.7, 200., 0., 200.);
      fHistos->CreateTH2(Form("hEvSelPatchADCEta%sSector%d", triggername.Data(), isec), Form("Patch online ADC versus eta for trigger %s, Sector %d in selectedEvents", triggername.Data(), isec), 100, -0.7, 0.7, 2100., 0., 2100.);
    }
  }
  PostData(1, fHistos->GetListOfHistograms());
}

void AliAnalysisTaskEmcalOfflinePatchesRef::UserExec(Option_t *){
  if(!fGeometry){
    fGeometry = AliEMCALGeometry::GetInstance();
    if(!fGeometry)
      fGeometry = AliEMCALGeometry::GetInstanceFromRunNumber(InputEvent()->GetRunNumber());
  }
  TClonesArray *patches = dynamic_cast<TClonesArray *>(fInputEvent->FindListObject("EmcalTriggers"));
  TString triggerstring = fInputEvent->GetFiredTriggerClasses();
  UInt_t selectionstatus = fInputHandler->IsEventSelected();
  Bool_t isMinBias = selectionstatus & AliVEvent::kINT7,
      isEG1 = (selectionstatus & AliVEvent::kEMCEGA) && triggerstring.Contains("EG1"),
      isEG2 = (selectionstatus & AliVEvent::kEMCEGA) && triggerstring.Contains("EG2"),
      isEMC7 = (selectionstatus & AliVEvent::kEMC7) && triggerstring.Contains("EMC7");
  if(!(isMinBias || isEG1 || isEG2 || isEMC7)) return;
  const AliVVertex *vtx = fInputEvent->GetPrimaryVertex();
  //if(!fInputEvent->IsPileupFromSPD(3, 0.8, 3., 2., 5.)) return;         // reject pileup event
  if(vtx->GetNContributors() < 1) return;
  if(fInputEvent->IsA() == AliESDEvent::Class() && fAnalysisUtil->IsFirstEventInChunk(fInputEvent)) return;
  // Fill reference distribution for the primary vertex before any z-cut
  if(!fAnalysisUtil->IsVertexSelected2013pA(fInputEvent)) return;       // Apply new vertex cut
  if(fAnalysisUtil->IsPileUpEvent(fInputEvent)) return;       // Apply new vertex cut
  // Apply vertex z cut
  if(vtx->GetZ() < -10. || vtx->GetZ() > 10.) return;
  if(isEMC7) fHistos->FillTH1("hEventCountEL0", 1);
  if(isEG1) fHistos->FillTH1("hEventCountEG1", 1);
  if(isEG2) fHistos->FillTH1("hEventCountEG2", 1);

  AliEMCALTriggerPatchInfo *mypatch(nullptr);
  Int_t supermoduleID = -1;
  TString patchname;
  for(TIter patchiter = TIter(patches).Begin(); patchiter != TIter::End(); ++patchiter){
    mypatch = dynamic_cast<AliEMCALTriggerPatchInfo *>(*patchiter);
    if(!mypatch) continue;
    // Select only gamma and L0 online patches
    if(!mypatch->IsOfflineSimple()) continue;
    if(!mypatch->IsGammaLowSimple()) continue;
    fGeometry->SuperModuleNumberFromEtaPhi(mypatch->GetEtaCM(), mypatch->GetPhiCM(), supermoduleID);
    Int_t sector = 4 + int(supermoduleID / 2);
    if(mypatch->IsGammaLowSimple() && mypatch->GetPatchE() > 3){
      FillTriggerPatchHistos("EL0", mypatch, supermoduleID, sector, kFALSE);
      FillTriggerPatchHistos("EG1", mypatch, supermoduleID, sector, kFALSE);
      FillTriggerPatchHistos("EG2", mypatch, supermoduleID, sector, kFALSE);
    }

    // correlate patches with real event classes selected
    if(isEMC7 && mypatch->IsGammaLowSimple()) FillTriggerPatchHistos("EL0", mypatch, supermoduleID, sector, kTRUE);
    if(isEG1 && mypatch->IsGammaLowSimple()) FillTriggerPatchHistos("EG1", mypatch, supermoduleID, sector, kTRUE);
    if(isEG2 && mypatch->IsGammaLowSimple()) FillTriggerPatchHistos("EG2", mypatch, supermoduleID, sector, kTRUE);
  }

  PostData(1, fHistos->GetListOfHistograms());
}

/**
 * Fill trigger patch histograms for a given patchtype with relevant information.
 * Information is also sorted according to sector and supermodule
 * Plots at global level:
 *  - Energy vs. supermodule
 *  - Energy vs. eta (all sectors)
 *  - Energy vs. eta for sector
 * @param patchtype Type of the reconstructed trigger patch
 * @param recpatch Reconstructed trigger patch with all information
 * @param supermoduleID ID of the supermodule
 * @param sector Sector in global numbering scheme
 * @param evsel Check whether additional event selection is applied
 */
void AliAnalysisTaskEmcalOfflinePatchesRef::FillTriggerPatchHistos(const char *patchtype, const AliEMCALTriggerPatchInfo * const recpatch, Int_t supermoduleID, Int_t sector, Bool_t evsel){
  TString fbase = evsel ? "hEvSel" : "h";
  fHistos->FillTH2(Form("%sPatchEnergy%s", fbase.Data(), patchtype), supermoduleID, recpatch->GetPatchE());
  fHistos->FillTH2(Form("%sPatchET%s", fbase.Data(), patchtype), supermoduleID, recpatch->GetLorentzVectorCenterGeo().Et());
  fHistos->FillTH2(Form("%sPatchADC%s", fbase.Data(), patchtype), supermoduleID, recpatch->GetADCAmp());
  fHistos->FillTH2(Form("%sPatchEnergyEta%s", fbase.Data(), patchtype), recpatch->GetEtaCM(), recpatch->GetPatchE());
  fHistos->FillTH2(Form("%sPatchETEta%s", fbase.Data(), patchtype), recpatch->GetEtaCM(), recpatch->GetLorentzVectorCenterGeo().Et());
  fHistos->FillTH2(Form("%sPatchADCEta%s", fbase.Data(), patchtype), recpatch->GetEtaCM(), recpatch->GetADCAmp());
  if(sector >= 4 && sector < 10){
    fHistos->FillTH2(Form("%sPatchEnergyEta%sSector%d", fbase.Data(), patchtype, sector), recpatch->GetEtaCM(), recpatch->GetPatchE());
    fHistos->FillTH2(Form("%sPatchETEta%sSector%d", fbase.Data(), patchtype, sector), recpatch->GetEtaCM(), recpatch->GetLorentzVectorCenterGeo().Et());
    fHistos->FillTH2(Form("%sPatchADCEta%sSector%d", fbase.Data(), patchtype, sector), recpatch->GetEtaCM(), recpatch->GetADCAmp());
  }
  if(supermoduleID >= 0 && supermoduleID < 10){
    fHistos->FillTH2(Form("%sPatchEnergyEta%sSM%d", fbase.Data(), patchtype, supermoduleID), recpatch->GetEtaCM(), recpatch->GetPatchE());
    fHistos->FillTH2(Form("%sPatchETEta%sSM%d", fbase.Data(), patchtype, supermoduleID), recpatch->GetEtaCM(), recpatch->GetLorentzVectorCenterGeo().Et());
    fHistos->FillTH2(Form("%sPatchADCEta%sSM%d", fbase.Data(), patchtype, supermoduleID), recpatch->GetEtaCM(), recpatch->GetADCAmp());
  }
}
