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

#include "ADESDFriendUtils.h"

ClassImp(AliAnalysisTaskADChargeMonitoring);

AliAnalysisTaskADChargeMonitoring::AliAnalysisTaskADChargeMonitoring(const char *name)
  : AliAnalysisTaskSE(name)
  , fTE(NULL)
  , fTimeStamp(0)
  , fBC(0)
  , fClassMask(0)
  , fClassMaskNext50(0)
  , fPSInt0(0)
  , fPSInt1(0)
  , fESDADfriendUtils(NULL) {

  for (Int_t i=0; i<AliADRawStream::kNScalers; ++i)
    fScalers[i] = 0;
  
  for (Int_t i=0; i<AliADRawStream::kNChannels; ++i) {
    fTriggerCharges[i] = 0.0f;
    fBBFlags[i]        = kFALSE;
  }

  DefineOutput(1, TTree::Class());
}

AliAnalysisTaskADChargeMonitoring::~AliAnalysisTaskADChargeMonitoring() {
  if (AliAnalysisManager::GetAnalysisManager() && 
      AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() == AliAnalysisManager::kProofAnalysis)
    return;

  delete fTE;
  fTE = NULL;
  
  delete fESDADfriendUtils;
  fESDADfriendUtils = NULL;
}

void AliAnalysisTaskADChargeMonitoring::UserCreateOutputObjects() {
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
  fTE->Branch("scalers",         &fScalers,
	      Form("val[%d]/i", AliADRawStream::kNChannels));
  fTE->Branch("BBFlags",         &fBBFlags,
	      Form("val[%d]/O", AliADRawStream::kNChannels));
  fTE->Branch("TriggerCharges",  &fTriggerCharges,
	      Form("val[%d]/F", AliADRawStream::kNChannels));
  owd->cd();

  PostData(1, fTE);

  fESDADfriendUtils = new ADESDFriendUtils;
}

void AliAnalysisTaskADChargeMonitoring::NotifyRun() {
  fESDADfriendUtils->Init(fCurrentRunNumber);
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

  fTimeStamp         = esdEvent->GetTimeStamp();
  fBC                = esdEvent->GetBunchCrossNumber();
  fClassMask       = esdEvent->GetTriggerMask();
  fClassMaskNext50 = esdEvent->GetTriggerMaskNext50();

  const TBits& ir1Map = esdHeader->GetIRInt1InteractionMap();
  const TBits& ir2Map = esdHeader->GetIRInt2InteractionMap();

  fPSInt0 = fPSInt1 = 0;
  for (Int_t i=0; i<32; ++i) {
    fPSInt0 |= (1U<<i)*ir1Map.TestBitNumber(90+i-16);
    fPSInt1 |= (1U<<i)*ir2Map.TestBitNumber(90+i-16);
  }

  fESDADfriendUtils->Update(esdADfriend);
  for (Int_t i=0; i<AliADRawStream::kNChannels; ++i) {
    fTriggerCharges[i] = fESDADfriendUtils->GetADCPedSub(i, 10);
    fBBFlags[i]        = esdADfriend->GetBBFlag(i, 10);
  }

  for (Int_t i=0; i<AliADRawStream::kNScalers; ++i)
    fScalers[i] = esdADfriend->GetTriggerScalers(i);

  fTE->Fill();

  PostData(1, fTE);
}

void AliAnalysisTaskADChargeMonitoring::Terminate(Option_t* ) {
  // NOP
}
