/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

#include <iostream>
#include <TList.h>

#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliMCEvent.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"

#include "AlidNdPtEventCuts.h"

using namespace std;

ClassImp(AlidNdPtEventCuts)

//_____________________________________________________________________________
AlidNdPtEventCuts::AlidNdPtEventCuts(const Char_t* name,const Char_t *title) : 
AliAnalysisCuts(name, title)
, fTriggerRequired(kTRUE)
, fRecVertexRequired(kTRUE)
, fEventProcessType(AliPWG0Helper::kInvalidProcess)
, fMinNContributors(0)
, fMaxNContributors(0)
, fMaxR(0)
, fMinZv(0)
, fMaxZv(0)
, fMeanXv(0)
, fMeanYv(0)
, fMeanZv(0)
, fSigmaMeanXv(0)
, fSigmaMeanYv(0)
, fSigmaMeanZv(0)
, fRedoTPCVertex(kTRUE)
, fUseBeamSpotConstraint(kTRUE)
{
  // default constructor 
  
  // init data members with defaults
  Init();
}

//_____________________________________________________________________________
AlidNdPtEventCuts::~AlidNdPtEventCuts()  
{
  // destructor
}

//_____________________________________________________________________________
void AlidNdPtEventCuts::Init()  
{
  // set default values
  SetTriggerRequired();
  SetRecVertexRequired();
  SetEventProcessType();
  SetNContributorsRange();
  SetMaxR();
  SetZvRange();
  SetMeanXYZv();
  SetSigmaMeanXYZv();
  SetRedoTPCVertex();
  SetUseBeamSpotConstraint();
}

//_____________________________________________________________________________
Bool_t AlidNdPtEventCuts::AcceptEvent(AliESDEvent *esdEvent,AliMCEvent *mcEvent, const AliESDVertex *vtx)
{
  // Check event selection cuts
  Bool_t retValue=kTRUE;

  if(!esdEvent) return kFALSE;
  if(!IsRecVertexRequired()) return kTRUE;
  if(!vtx) return kFALSE;
  if(!vtx->GetStatus()) return kFALSE;

  if(mcEvent) {
   // check MC event conditions
   AliHeader* header = mcEvent->Header();
   if(!header) return kFALSE;
  
    // select event type (ND-non diffractive, SD-single diffractive, DD-double diffractive)
    if(fEventProcessType == AliPWG0Helper::kInvalidProcess) { 
      retValue=kTRUE;
    } 
    else if(fEventProcessType == AliPWG0Helper::kSD || fEventProcessType == AliPWG0Helper::kDD) {
      AliPWG0Helper::MCProcessType processType = AliPWG0Helper::GetEventProcessType(header);
      if(processType == AliPWG0Helper::kND) retValue=kFALSE;
      else retValue=kTRUE;
    }
    else if(fEventProcessType == AliPWG0Helper::GetEventProcessType(header)) { 
      retValue=kTRUE;
    }
    else 
      retValue=kFALSE;
  }

  /*
  Float_t R = TMath::Sqrt(vtx->GetXv()*vtx->GetXv()+vtx->GetYv()*vtx->GetYv());
  if(vtx->GetNContributors() < fMinNContributors) return kFALSE; 
  if(vtx->GetNContributors() > fMaxNContributors) return kFALSE; 
  if(R > fMaxR) return kFALSE; 
  */

  if(vtx->GetZv() < fMinZv) return kFALSE; 
  if(vtx->GetZv() > fMaxZv) return kFALSE; 

return retValue;  
}

//_____________________________________________________________________________
Bool_t AlidNdPtEventCuts::AcceptMCEvent(AliMCEvent *mcEvent)
{
  // Check event selection cuts
  if(!mcEvent) return kFALSE;

  Bool_t retValue=kTRUE;

  // check MC event conditions
  AliHeader* header = mcEvent->Header();
  if(!header) return kFALSE;

  AliGenEventHeader* genHeader = header->GenEventHeader();
  if (!genHeader) {
    AliDebug(AliLog::kError, "Could not retrieve genHeader from Header");
    return kFALSE;
  }
  TArrayF vtxMC(3);
  genHeader->PrimaryVertex(vtxMC);
  
  // select event type (ND-non diffractive, SD-single diffractive, DD-double diffractive)
  if(fEventProcessType == AliPWG0Helper::kInvalidProcess) { 
     retValue=kTRUE;
  } else {
     if(fEventProcessType == AliPWG0Helper::GetEventProcessType(header)) retValue=kTRUE;
     else retValue=kFALSE;
  }

  /*
  Float_t R = TMath::Sqrt(vtxMC[0]*vtxMC[0]+vtxMC[1]*vtxMC[1]);
  if(R > fMaxR) return kFALSE; 
  */

  if(vtxMC[2] < fMinZv) return kFALSE; 
  if(vtxMC[2] > fMaxZv) return kFALSE; 

return retValue;  
}

//_____________________________________________________________________________
Long64_t AlidNdPtEventCuts::Merge(TCollection* list) 
{
  // Merge list of objects (needed by PROOF)
  if (!list)
  return 0;

  if (list->IsEmpty())
  return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj = 0;

  Int_t count=0;
  while((obj = iter->Next()) != 0) 
  {
    AlidNdPtEventCuts* entry = dynamic_cast<AlidNdPtEventCuts*>(obj);
    if (entry == 0)  
      continue; 

  count++;
  }

return count;
}
