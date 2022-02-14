/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

#include <TH1D.h>
#include <TList.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>

#include "AliVEventHandler.h"
#include "AliAnalysisManager.h"

#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliVParticle.h"

#include "AliMuonTrackCuts.h"
#include "AliMuonPairCuts.h"
#include "AliAnalysisMuonUtility.h"
#include "AliMultSelection.h"

#include "AliPicoDQ.h"
#include "AliPicoDQheader.h"
#include "AliAnalysisTaskSEMultDQdev.h"

ClassImp(AliAnalysisTaskSEMultDQdev)

//_____________________________________________________________________________
AliAnalysisTaskSEMultDQdev::AliAnalysisTaskSEMultDQdev() :
AliAnalysisTaskSE(),
fHeader(nullptr),
fDQClArr(nullptr),
fCutsDimu(nullptr),
fListOutputs(nullptr)
{
//
//  AliAnalysisTaskSEMultDQdev::AliAnalysisTaskSEMultDQdev() :
//
}

//_____________________________________________________________________________
AliAnalysisTaskSEMultDQdev::AliAnalysisTaskSEMultDQdev(const char *s) :
AliAnalysisTaskSE(s),
fHeader(nullptr),
fDQClArr(nullptr),
fCutsDimu(nullptr),
fListOutputs(nullptr)
{
//
//  AliAnalysisTaskSEMultDQdev::AliAnalysisTaskSEMultDQdev(const char *s, const Bool_t b) :
//

  DefineOutput(1, TList::Class());
}

//_____________________________________________________________________________
AliAnalysisTaskSEMultDQdev::~AliAnalysisTaskSEMultDQdev()
{
//
//  AliAnalysisTaskSEMultDQdev::~AliAnalysisTaskSEMultDQdev()
//

  if (fHeader)  { delete fHeader;  fHeader  = nullptr; }
  if (fDQClArr) { delete fDQClArr; fDQClArr = nullptr; }

  if (fCutsDimu) { delete fCutsDimu; fCutsDimu = nullptr; }

  if (fListOutputs) { delete fListOutputs; fListOutputs = nullptr; }
}

//_____________________________________________________________________________
void AliAnalysisTaskSEMultDQdev::UserCreateOutputObjects()
{
//
//  void AliAnalysisTaskSEMultDQdev::UserCreateOutputObjects()
//
  if (!fCutsDimu) {
    fCutsDimu = new AliMuonPairCuts("AliDQcuts", "ALICE DQ cuts");
    fCutsDimu->GetMuonTrackCuts().SetAllowDefaultParams();
    fCutsDimu->SetFilterMask(fCutsDimu->GetFilterMask() | AliMuonPairCuts::kBothMuPdca);
//    if (fUpsilon) fCutsDimu->IsApplySharpPtCutInMatching();
  }

  fCutsDimu->Print("pair");
//=============================================================================

  if (!fHeader) {
    fHeader = new AliPicoDQheader();
    fHeader->SetName("PicoDQheader");
    AddAODBranch("AliPicoDQheader", &fHeader);
  }

  if (!fDQClArr) {
    fDQClArr = new TClonesArray("AliPicoDQ", 0);
    fDQClArr->SetName("PicoDQ");
    AddAODBranch("TClonesArray", &fDQClArr);
  }
//=============================================================================

  if (fListOutputs) { delete fListOutputs; fListOutputs = nullptr; }

  fListOutputs = new TList();
  fListOutputs->SetOwner();

  const auto bStatus(TH1::AddDirectoryStatus());
  TH1::AddDirectory(kFALSE);

  auto h(new TH1D("hEvent", "", 5, -0.5, 4.5));

  fListOutputs->Add(h);

  TH1::AddDirectory(bStatus);
  PostData(1, fListOutputs);

//=============================================================================

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEMultDQdev::UserExec(Option_t *)
{
//
//  void AliAnalysisTaskSEMultDQdev::UserExec(Option_t *)
//
  fHeader->Reset();
  fDQClArr->Delete();
//=============================================================================

  if (IsNotEventSelected()) return;
//=============================================================================

  AliAnalysisMuonUtility::SetUseSmearedTracks(kFALSE,kFALSE);
  const auto ntrks(AliAnalysisMuonUtility::GetNTracks(InputEvent()));
//=============================================================================

  auto l(fDQClArr->GetEntriesFast());
  for (auto i=0, k=1; k<ntrks; ++i, ++k) {
    const auto p1(AliAnalysisMuonUtility::GetTrack(i,InputEvent()));
    
    if (!p1) continue;
 
    const auto w(p1->Charge());
    for (auto j=k; j<ntrks; ++j) {
      const auto p2(AliAnalysisMuonUtility::GetTrack(j,InputEvent()));
      
      if (!p2) continue;

      if ((w*(p2->Charge()))>0) continue;
      if (!(fCutsDimu->IsSelected(p1,p2))) continue;
  
      auto v(AliAnalysisMuonUtility::GetTrackPair(p1,p2));

      new ((*fDQClArr)[l++]) AliPicoDQ(v);
    }
  }
//=============================================================================

  AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()->SetFillAOD(kTRUE);
//=============================================================================

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEMultDQdev::Terminate(Option_t *)
{
//
//  void AliAnalysisTaskSEMultDQdev::Terminate(Option_t *)
//

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEMultDQdev::NotifyRun()
{
//
//  void AliAnalysisTaskSEMultDQdev::NotifyRun()
//

  if (fCutsDimu) fCutsDimu->SetRun(fInputHandler);
//=============================================================================

  return;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskSEMultDQdev::IsNotEventSelected()
{
//
//  Bool_t AliAnalysisTaskSEMultDQdev::IsNotEventSelected()
//

  const auto pAOD(dynamic_cast<AliAODEvent*>(InputEvent()));
  const auto pESD(dynamic_cast<AliESDEvent*>(InputEvent()));
  if (!(pAOD || pESD)) { AliFatal("Event is neither of AOD nor ESD type"); return kTRUE; }
//=============================================================================

  auto h(static_cast<TH1D*>(fListOutputs->FindObject("hEvent")));
  h->Fill(0.);
//=============================================================================

  const auto pSPDVtx(pAOD ? dynamic_cast<const AliVVertex*>(pAOD->GetPrimaryVertexSPD()) :
                            dynamic_cast<const AliVVertex*>(pESD->GetPrimaryVertexSPD()));

  if (!pSPDVtx) return kTRUE;  
//=============================================================================

  const auto dSPDVz(pSPDVtx->GetZ());

  if ((dSPDVz<-10.) || (dSPDVz>=10.))   return kTRUE;
  h->Fill(1.);

  if ((pSPDVtx->GetNContributors())<=0)  return kTRUE;
  h->Fill(2.);

  Double_t dCov[6]; pSPDVtx->GetCovarianceMatrix(dCov);
  if (TMath::Sqrt(dCov[5])>=0.25)      return kTRUE;
  h->Fill(3.);  
//=============================================================================

  auto pms(dynamic_cast<AliMultSelection*>(InputEvent()->FindListObject("MultSelection")));
  if (!pms) return kTRUE;

  const auto dmultSPDtrkls(pms->GetMultiplicityPercentile("SPDTracklets"));
  const auto dmultV0M(pms->GetMultiplicityPercentile("V0M"));
  const auto dmultV0C(pms->GetMultiplicityPercentile("V0C"));

  const auto pm(InputEvent()->GetMultiplicity());
  if (!pm) return kTRUE;

  UInt_t wmt(0); 
  for (auto i=0; i<pm->GetNumberOfTracklets(); ++i) {
    const auto d(pm->GetEta(i));
    if ((d<-1.) || (d>=1.)) continue;
    wmt += 1;
  }
  
  h->Fill(4.);
//=============================================================================
  
  fHeader->SetEventInfo(
    fInputHandler->IsEventSelected(),
    pAOD ? pAOD->GetHeader()->GetL0TriggerInputs() : pESD->GetHeader()->GetL0TriggerInputs(),
    pAOD ? pAOD->GetFiredTriggerClasses()          : pESD->GetFiredTriggerClasses(),
    dSPDVz, dmultSPDtrkls, dmultV0M, dmultV0C, wmt);
//=============================================================================

  return kFALSE;
}

