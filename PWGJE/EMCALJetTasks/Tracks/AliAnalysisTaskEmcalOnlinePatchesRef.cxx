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
#include <TString.h>

#include "AliAnalysisUtils.h"
#include "AliEMCALGeometry.h"
#include "AliEmcalTriggerPatchInfo.h"
#include "AliEMCalHistoContainer.h"
#include "AliInputEventHandler.h"
#include "AliVEvent.h"
#include "AliVVertex.h"

#include "AliAnalysisTaskEmcalOnlinePatchesRef.h"

#if __cplusplus < 201103L
/*
 * Old C++
 */
#define nullptr NULL
#include <vector>
#else
#include <array>
#endif

ClassImp(EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalOnlinePatchesRef)

namespace EMCalTriggerPtAnalysis {

AliAnalysisTaskEmcalOnlinePatchesRef::AliAnalysisTaskEmcalOnlinePatchesRef():
  AliAnalysisTaskSE(),
  fAnalysisUtil(nullptr),
  fGeometry(nullptr),
  fHistos(nullptr)
{
}

AliAnalysisTaskEmcalOnlinePatchesRef::AliAnalysisTaskEmcalOnlinePatchesRef(const char *name):
  AliAnalysisTaskSE(name),
  fAnalysisUtil(nullptr),
  fGeometry(nullptr),
  fHistos(nullptr)
{
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskEmcalOnlinePatchesRef::~AliAnalysisTaskEmcalOnlinePatchesRef() {
}

void AliAnalysisTaskEmcalOnlinePatchesRef::UserCreateOutputObjects(){
  fAnalysisUtil = new AliAnalysisUtils;

  fHistos = new AliEMCalHistoContainer("Ref");
  TString triggername;
  // Plots at global level:
  // Energy vs. supermodule
  // Energy vs. eta (all sectors)
  // Energy vs. eta for sector
#if __cplusplus >= 201103L
  /*
   * Version for beautifull C++11
   */
  std::array<TString, 3> patchnames = {
      "EL0", "EG1", "EG2"
  };
  for(auto mytrg : patchnames){
    triggername = mytrg;
#else
  /*
   * Backward compatible version for the ancient technology
   */
  std::vector<TString> patchnames;
  patchnames.push_back("EL0");
  patchnames.push_back("EG1");
  patchnames.push_back("EG2");
  for(std::vector<TString>::iterator mytrg = patchnames.begin(); mytrg != patchnames.end(); ++mytrg){
    triggername = *mytrg;
#endif
    fHistos->CreateTH2(Form("hPatchEnergy%s", triggername.Data()), Form("Patch energy versus supermodule for trigger %s", triggername.Data()), 12, -0.5, 11.5, 200, 0., 200.);
    fHistos->CreateTH2(Form("hPatchADC%s", triggername.Data()), Form("Patch online ADC versus supermodule for trigger %s", triggername.Data()), 12, -0.5, 11.5, 2100, 0., 2100.);
    fHistos->CreateTH2(Form("hPatchEnergyEta%s", triggername.Data()), Form("Patch energy versus eta for trigger %s", triggername.Data()), 100, -0.7, 0.7, 200., 0., 200.);
    fHistos->CreateTH2(Form("hPatchADCEta%s", triggername.Data()), Form("Patch energy versus eta for trigger %s", triggername.Data()), 100, -0.7, 0.7, 2100., 0., 2100.);
    for(int isec = 5; isec <= 10; isec++){
      fHistos->CreateTH2(Form("hPatchEnergyEta%sSector%d", triggername.Data(), isec), Form("Patch energy versus eta for trigger %s", triggername.Data()), 100, -0.7, 0.7, 200., 0., 200.);
      fHistos->CreateTH2(Form("hPatchADCEta%sSector%d", triggername.Data(), isec), Form("Patch energy versus eta for trigger %s", triggername.Data()), 100, -0.7, 0.7, 2100., 0., 2100.);
    }
  }
  PostData(1, fHistos->GetListOfHistograms());
}

void AliAnalysisTaskEmcalOnlinePatchesRef::UserExec(Option_t *){
  if(!fGeometry){
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
  // Fill reference distribution for the primary vertex before any z-cut
  if(!fAnalysisUtil->IsVertexSelected2013pA(fInputEvent)) return;       // Apply new vertex cut
  if(fAnalysisUtil->IsPileUpEvent(fInputEvent)) return;       // Apply new vertex cut
  // Apply vertex z cut
  if(vtx->GetZ() < -10. || vtx->GetZ() > 10.) return;

  AliEmcalTriggerPatchInfo *mypatch(nullptr);
  Int_t supermoduleID = -1;
  TString patchname;
  for(TIter patchiter = TIter(patches).Begin(); patchiter != TIter::End(); ++patchiter){
    mypatch = dynamic_cast<AliEmcalTriggerPatchInfo *>(*patchiter);
    if(!mypatch) continue;
    // Select only gamma and L0 online patches
    if(mypatch->IsOfflineSimple()) continue;
    if(!(mypatch->IsGammaHigh() || mypatch->IsGammaLow() || mypatch->IsLevel0())) continue;
    fGeometry->SuperModuleNumberFromEtaPhi(mypatch->GetEtaCM(), mypatch->GetPhiCM(), supermoduleID);
    Int_t sector = 5 + supermoduleID / 2;
    if(mypatch->IsLevel0()) FillTriggerPatchHistos("EL0", mypatch, supermoduleID, sector);
    if(mypatch->IsGammaHigh()) FillTriggerPatchHistos("EG1", mypatch, supermoduleID, sector);
    if(mypatch->IsGammaLow()) FillTriggerPatchHistos("EG2", mypatch, supermoduleID, sector);
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
 */
void AliAnalysisTaskEmcalOnlinePatchesRef::FillTriggerPatchHistos(const char *patchtype, const AliEmcalTriggerPatchInfo * const recpatch, Int_t supermoduleID, Int_t sector){
  fHistos->FillTH2(Form("hPatchEnergy%s", patchtype), supermoduleID, recpatch->GetPatchE());
  fHistos->FillTH2(Form("hPatchADC%s", patchtype), supermoduleID, recpatch->GetADCAmp());
  fHistos->FillTH2(Form("hPatchEnergyEta%s", patchtype), recpatch->GetEtaCM(), recpatch->GetPatchE());
  fHistos->FillTH2(Form("hPatchADCEta%s", patchtype), recpatch->GetEtaCM(), recpatch->GetADCAmp());
  fHistos->FillTH2(Form("hPatchEnergyEta%sSector%d", patchtype, sector), recpatch->GetEtaCM(), recpatch->GetPatchE());
  fHistos->FillTH2(Form("hPatchADCEta%sSector%d", patchtype, sector), recpatch->GetEtaCM(), recpatch->GetADCAmp());
}

} /* namespace EMCalTriggerPtAnalysis */
