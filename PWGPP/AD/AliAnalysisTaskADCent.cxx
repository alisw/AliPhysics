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
#include <TTree.h>

#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliAODEvent.h"
#include "AliAODHeader.h"

#include "AliAnalysisTaskADCent.h"

ClassImp(AliAnalysisTaskADCent);

AliAnalysisTaskADCent::AliAnalysisTaskADCent(const char *name)
  : AliAnalysisTaskSE(name)
  , fAD()
  , fIR1Map()
  , fIR2Map()
  , fCent()
{
  DefineOutput(1, TTree::Class());
}

AliAnalysisTaskADCent::~AliAnalysisTaskADCent() {
  if (AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() == AliAnalysisManager::kProofAnalysis)
    return;

  SafeDelete(fTE);
}


void AliAnalysisTaskADCent::NotifyRun() {
  // NOP for now
}

void AliAnalysisTaskADCent::UserCreateOutputObjects() {
  TDirectory *owd = gDirectory;
  OpenFile(1);
  fTE = new TTree("TE", "");
  fTE->Branch("fAD",     &fAD,     32000, 1);
  fTE->Branch("fIR1Map", &fIR1Map, 32000, 1);
  fTE->Branch("fIR2Map", &fIR2Map, 32000, 1);
  fTE->Branch("fCent",   &fCent,   32000, 1);
  PostData(1, fTE);
  owd->cd();
}

void AliAnalysisTaskADCent::UserExec(Option_t* ) {
  AliAnalysisManager* pManager(AliAnalysisManager::GetAnalysisManager());
  if (NULL == pManager) return;

  AliInputEventHandler* pInputHandler(dynamic_cast<AliInputEventHandler*>(pManager->GetInputEventHandler()));
  if (NULL == pInputHandler) return;

  if (!pInputHandler->IsEventSelected())
    return;

  AliAODEvent* pAOD(dynamic_cast<AliAODEvent*>(InputEvent()));
  if (NULL == pAOD) return;

  AliAODHeader *pAODHeader = dynamic_cast<AliAODHeader*>(pAOD->GetHeader());
  if(!pAODHeader) AliFatal("Not a standard AOD");
  if (NULL == pAODHeader) return;

  AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(InputEvent());
  if (NULL == aodEvent) {
    AliError("NULL == esdEvent");
    return;
  }
  const AliAODAD* aodAD = aodEvent->GetADData();
  if (NULL == aodAD) {
    AliError("NULL == aodAD");
    return;
  }
  const AliCentrality *mult = pAODHeader->GetCentralityP();
  if (!mult) {
    AliError("NULL == mult");
    return;
  }
  fAD     = *aodAD;
  fIR1Map =  pAODHeader->GetIRInt1InteractionMap();
  fIR2Map =  pAODHeader->GetIRInt2InteractionMap();
  fCent   = *mult;

  fTE->Fill();

  PostData(1, fTE);
}

void AliAnalysisTaskADCent::Terminate(Option_t* ) {
  // NOP
}

