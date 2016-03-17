// -*- C++ -*-
// $Id$

/**************************************************************************
 * Author: C. Mayer                                                       *
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

#include <TFile.h>
#include <TList.h>
#include <TH2.h>
#include <TString.h>
#include <TBits.h>

#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskADChargeMonitoring.h"
#include "AliRawEventHeaderBase.h"

#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliESDAD.h"
#include "AliESDADfriend.h"

#include "AliCDBManager.h"
#include "AliCDBEntry.h"

#include "AliAnalysisUtils.h"
#include "ADESDFriendUtils.h"
#include "AliTriggerRunScalers.h"
#include "AliTriggerScalers.h"
#include "AliTriggerScalersRecord.h"

ClassImp(AliAnalysisTaskADChargeMonitoring);

AliAnalysisTaskADChargeMonitoring::AliAnalysisTaskADChargeMonitoring(const char *name)
  : AliAnalysisTaskSE(name)
  , fFillTTree(kFALSE)
  , fTL(NULL)
  , fTE(NULL)
  , fTimeStamp(0)
  , fBC(0)
  , fClassMask(0)
  , fClassMaskNext50(0)
  , fPSInt0(0)
  , fPSInt1(0)
  , fAnalysisUtils(NULL)
  , fESDADfriendUtils(NULL) {

  for (Int_t i=0; i<AliADRawStream::kNScalers; ++i)
    fScalers[i] = 0;
  
  for (Int_t i=0; i<AliADRawStream::kNChannels; ++i) {
    fTriggerCharges[i] = 0.0f;
    fBBFlags[i]        = kFALSE;
    fHChargeTime[i]    = NULL;
  }

  fPileUp[0] = fPileUp[1] = kFALSE;

  DefineOutput(1, TList::Class());
}

AliAnalysisTaskADChargeMonitoring::~AliAnalysisTaskADChargeMonitoring() {
  if (AliAnalysisManager::GetAnalysisManager() && 
      AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() == AliAnalysisManager::kProofAnalysis)
    return;

  delete fTL;
  fTL = NULL;
  
  delete fAnalysisUtils;
  fAnalysisUtils = NULL;

  delete fESDADfriendUtils;
  fESDADfriendUtils = NULL;
}

void AliAnalysisTaskADChargeMonitoring::UserCreateOutputObjects() {
  fTL = new TList;
  fTL->SetOwner(kTRUE);

  TDirectory *owd = gDirectory;
  OpenFile(1);
  fTE = new TTree;
  fTE->SetName("TE");
  fTE->Branch("timeStamp",       &fTimeStamp);
  fTE->Branch("bc",              &fBC);
  fTE->Branch("ClassMask",       &fClassMask);
  fTE->Branch("ClassMaskNext50", &fClassMaskNext50);
  fTE->Branch("PSInt0",          &fPSInt0);
  fTE->Branch("PSInt1",          &fPSInt1);  
  fTE->Branch("PileUp",          &fPileUp, "pileup/O:SPDClusterVsTracklet");
  fTE->Branch("scalers",         &fScalers,
	      Form("val[%d]/i", AliADRawStream::kNChannels));
  fTE->Branch("BBFlags",         &fBBFlags,
	      Form("val[%d]/O", AliADRawStream::kNChannels));
  fTE->Branch("TriggerCharges",  &fTriggerCharges,
	      Form("val[%d]/F", AliADRawStream::kNChannels));
  owd->cd();
  fTL->Add(fTE);

  PostData(1, fTL);

  fAnalysisUtils    = new AliAnalysisUtils;
  fESDADfriendUtils = new ADESDFriendUtils;
}

void AliAnalysisTaskADChargeMonitoring::NotifyRun() { 
  AliCDBManager *man   = fESDADfriendUtils->Init(fCurrentRunNumber);
  AliCDBEntry   *entry = man->Get("GRP/CTP/Scalers");
  if (NULL == entry) {
    AliFatal("NULL == entry"); return;
  }
  const AliTriggerRunScalers *triggerScalers = dynamic_cast<const AliTriggerRunScalers*>(entry->GetObject());
  if (NULL == triggerScalers) {
    AliFatal("NULL == triggerScalers"); return;
  }

  const TObjArray *a = triggerScalers->GetScalersRecords();
  if (NULL == a) {
    AliFatal("NULL == a"); return;
  }
  const AliTriggerScalersRecord *sFirst = dynamic_cast<const AliTriggerScalersRecord*>(a->First());
  const AliTriggerScalersRecord *sLast  = dynamic_cast<const AliTriggerScalersRecord*>(a->Last());

  const AliTimeStamp *tFirst = sFirst->GetTimeStamp();
  const AliTimeStamp *tLast  = sLast->GetTimeStamp();

  const Double_t dt = 2*60;
  const Double_t t0 = dt*Int_t(tFirst->GetSeconds()/dt);
  const Double_t t1 = dt*Int_t(tLast->GetSeconds() /dt) + dt;

  for (Int_t ch=0; ch<AliADRawStream::kNChannels; ++ch) {
    const TString histName  = TString::Format("ch%02d", ch);
    const TString histTitle = TString::Format("ch%0d;time;trigger charge ch%02d (ADC)", ch, ch);
    fHChargeTime[ch] = new TH2F(histName,
				histTitle,
				Int_t((t1-t0)/dt+0.5), t0, t1,
				256, 0, 1024);
    fTL->Add(fHChargeTime[ch]);
  }
  

}

void AliAnalysisTaskADChargeMonitoring::UserExec(Option_t* ) {
  AliESDEvent* esdEvent = dynamic_cast<AliESDEvent*>(InputEvent());
  if (NULL == esdEvent) {
    AliError("NULL == esdEvent");
    return;
  }

  AliESDHeader* esdHeader = dynamic_cast<AliESDHeader*>(esdEvent->GetHeader());
  if (NULL == esdHeader) {
    AliError("NULL == esdHeader");
    return;
  }

  AliESDfriend *esdFriend = esdEvent->FindFriend();  
  if (NULL == esdFriend) {
    AliError("NULL == esdFriend");
    return;
  }
  AliESDADfriend* esdADfriend = esdFriend->GetADfriend();
  if (NULL == esdADfriend) {
    AliError("NULL == esdADfriend");
    return;
  }

  if (esdEvent->GetEventType() != AliRawEventHeaderBase::kPhysicsEvent)
    return;

  fTimeStamp       = esdEvent->GetTimeStamp();
  fBC              = esdEvent->GetBunchCrossNumber();
  fClassMask       = esdEvent->GetTriggerMask();
  fClassMaskNext50 = esdEvent->GetTriggerMaskNext50();

  const TBits& ir1Map = esdHeader->GetIRInt1InteractionMap();
  const TBits& ir2Map = esdHeader->GetIRInt2InteractionMap();

  fPSInt0 = fPSInt1 = 0;
  for (Int_t i=0; i<32; ++i) {
    fPSInt0 |= (1U<<i)*ir1Map.TestBitNumber(90+i-16);
    fPSInt1 |= (1U<<i)*ir2Map.TestBitNumber(90+i-16);
  }

  fPileUp[0] = fAnalysisUtils->IsPileUpEvent(esdEvent);
  fPileUp[1] = fAnalysisUtils->IsSPDClusterVsTrackletBG(esdEvent);

  fESDADfriendUtils->Update(esdADfriend);
  for (Int_t ch=0; ch<AliADRawStream::kNChannels; ++ch) {
    fTriggerCharges[ch] = fESDADfriendUtils->GetADCPedSub(ch, 10);
    fBBFlags[ch]        = esdADfriend->GetBBFlag(ch, 10);

    if (fPSInt0 == (1U<<16) &&
	fPSInt1 == (1U<<16) &&
	fBBFlags[ch]        &&
	!fPileUp[0]         &&
	!fPileUp[1])
      fHChargeTime[ch]->Fill(fTimeStamp, fTriggerCharges[ch]);
  }

  for (Int_t i=0; i<AliADRawStream::kNScalers; ++i)
    fScalers[i] = esdADfriend->GetTriggerScalers(i);

  if (fFillTTree)
    fTE->Fill();

  PostData(1, fTL);
}

void AliAnalysisTaskADChargeMonitoring::Terminate(Option_t* ) {
  // NOP
}
