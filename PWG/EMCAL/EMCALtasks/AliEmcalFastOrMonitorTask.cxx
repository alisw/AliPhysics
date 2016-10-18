/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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
#include <algorithm>
#include <vector>
#include <THashList.h>
#include <THistManager.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TVector3.h>

#include "AliEmcalFastOrMonitorTask.h"
#include "AliEMCALGeometry.h"
#include "AliInputEventHandler.h"
#include "AliVCaloTrigger.h"
#include "AliVEvent.h"
#include "AliVVertex.h"

/// \cond CLASSIMP
ClassImp(AliEmcalFastOrMonitorTask)
/// \endcond

AliEmcalFastOrMonitorTask::AliEmcalFastOrMonitorTask() :
  AliAnalysisTaskSE(),
  fHistos(nullptr),
  fGeom(nullptr),
  fLocalInitialized(false),
  fOldRun(-1),
  fRequestTrigger(AliVEvent::kAny),
  fTriggerPattern("")
{

}

AliEmcalFastOrMonitorTask::AliEmcalFastOrMonitorTask(const char *name) :
  AliAnalysisTaskSE(name),
  fHistos(nullptr),
  fGeom(nullptr),
  fLocalInitialized(false),
  fOldRun(-1),
  fRequestTrigger(AliVEvent::kAny),
  fTriggerPattern("")
{
  DefineOutput(1, TList::Class());
}

AliEmcalFastOrMonitorTask::~AliEmcalFastOrMonitorTask() {

}

void AliEmcalFastOrMonitorTask::UserCreateOutputObjects() {
  fHistos = new THistManager("fastOrHistos");

  const int kMaxCol = 48, kMaxRow = 104, kMaxFastOr = kMaxRow * kMaxCol;

  fHistos->CreateTH1("hEvents", "Number of events", 1, 0.5, 1.5);
  fHistos->CreateTH1("hFastOrFrequencyL0", "FastOr frequency at Level0", kMaxFastOr, -0.5, kMaxFastOr - 0.5);
  fHistos->CreateTH1("hFastOrFrequencyL1", "FastOr frequency at Level1", kMaxFastOr, -0.5, kMaxFastOr - 0.5);
  fHistos->CreateTH2("hFastOrAmplitude", "FastOr amplitudes", kMaxFastOr, -0.5, kMaxFastOr - 0.5, 513, -0.5, 512.5);
  fHistos->CreateTH2("hFastOrTimeSum", "FastOr time sum", kMaxFastOr, -0.5, kMaxFastOr - 0.5, 2049, -0.5, 2048.5);
  fHistos->CreateTH2("hFastOrTransverseTimeSum", "FastOr transverse time sum", kMaxFastOr, -0.5, kMaxFastOr - 0.5, 2049, -0.5, 2048.5);
  fHistos->CreateTH2("hFastOrNL0Times", "FastOr Number of L0 times", kMaxFastOr, -0.5, kMaxFastOr - 0.5, 16, -0.5, 15.5);
  fHistos->CreateTH2("hFastOrColRowFrequencyL0", "FastOr Frequency (col-row) at Level1", kMaxCol, -0.5, kMaxCol - 0.5, kMaxRow, -0.5, kMaxRow - 0.5);
  fHistos->CreateTH2("hFastOrColRowFrequencyL1", "FastOr Frequency (col-row) at Level0", kMaxCol, -0.5, kMaxCol - 0.5, kMaxRow, -0.5, kMaxRow - 0.5);

  PostData(1, fHistos->GetListOfHistograms());
}

void AliEmcalFastOrMonitorTask::ExecOnce(){
  fGeom = AliEMCALGeometry::GetInstanceFromRunNumber(InputEvent()->GetRunNumber());
}

void AliEmcalFastOrMonitorTask::RunChanged(){

}

void AliEmcalFastOrMonitorTask::UserExec(Option_t *) {
  if(!fLocalInitialized){
    ExecOnce();
    fLocalInitialized = true;
  }

  // Run change
  if(InputEvent()->GetRunNumber() != fOldRun){
    RunChanged();
    fOldRun = InputEvent()->GetRunNumber();
  }

  // Check trigger
  if(!(fInputHandler->IsEventSelected() & fRequestTrigger)) return;
  if(fTriggerPattern.Length()){
    if(!TString(InputEvent()->GetFiredTriggerClasses()).Contains(fTriggerPattern)) return;
  }

  const AliVVertex *vtx = fInputEvent->GetPrimaryVertex();
  Double_t vtxpos[3];
  vtx->GetXYZ(vtxpos);

  fHistos->FillTH1("hEvents", 1);

  AliVCaloTrigger *triggerdata = InputEvent()->GetCaloTrigger("EMCAL");
  triggerdata->Reset();
  Int_t nl0times, l1timesum, fastOrID, globCol, globRow;
  Float_t amp;
  while(triggerdata->Next()){
    triggerdata->GetAmplitude(amp);
    triggerdata->GetNL0Times(nl0times);
    triggerdata->GetL1TimeSum(l1timesum);
    triggerdata->GetPosition(globCol, globRow);
    fGeom->GetTriggerMapping()->GetAbsFastORIndexFromPositionInEMCAL(globCol, globRow, fastOrID);
    if(amp > 1e-5){
      fHistos->FillTH2("hFastOrColRowFrequencyL0", globCol, globRow);
      fHistos->FillTH1("hFastOrFrequencyL0", fastOrID);
    }
    if(l1timesum){
      fHistos->FillTH2("hFastOrColRowFrequencyL1", globCol, globRow);
      fHistos->FillTH1("hFastOrFrequencyL1", fastOrID);
    }
    if(std::find(fMaskedFastors.begin(), fMaskedFastors.end(), fastOrID) == fMaskedFastors.end()){
      fHistos->FillTH2("hFastOrAmplitude", fastOrID, amp);
      fHistos->FillTH2("hFastOrTimeSum", fastOrID, l1timesum);
      fHistos->FillTH2("hFastOrNL0Times", fastOrID, nl0times);
      fHistos->FillTH2("hFastOrTransverseTimeSum", fastOrID, GetTransverseTimeSum(fastOrID, l1timesum, vtxpos));
    }
  }

  PostData(1, fHistos->GetListOfHistograms());
}

Double_t AliEmcalFastOrMonitorTask::GetTransverseTimeSum(Int_t fastorAbsID, Double_t adc, const Double_t *vertex) const{
  Int_t cellIDs[4];
  fGeom->GetTriggerMapping()->GetCellIndexFromFastORIndex(fastorAbsID, cellIDs);
  std::vector<double> eta, phi;
  for(int i = 0l; i < 4; i++){
    double etatmp, phitmp;
    fGeom->EtaPhiFromIndex(cellIDs[i], etatmp, phitmp);
    eta.push_back(etatmp);
    phi.push_back(phitmp);
  }

  // Calculate FastOR position: for the approximation take mean eta and phi
  // Radius is taken from the geometry
  TVector3 fastorPos, vertexPos(vertex[0], vertex[1], vertex[2]);
  fastorPos.SetPtEtaPhi(fGeom->GetIPDistance(), TMath::Mean(eta.begin(), eta.end()), TMath::Mean(phi.begin(), phi.end()));
  fastorPos -= vertexPos;

  TLorentzVector evec(fastorPos, adc);
  return evec.Et();
}
