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
//
// Cut step class
// Select all tracks surviving cuts in one special cut step
// Used in AliHFEtrackFilter
// 
// Author:
//   Markus Fasel <M.Fasel@gsi.de>
//
#include <TObjArray.h>

#include "AliAnalysisCuts.h"
#include "AliCFCutBase.h"
#include "AliHFEcutStep.h"
#include "AliLog.h"
#include "AliMCEvent.h"

ClassImp(AliHFEcutStep)

//__________________________________________________________________
AliHFEcutStep::AliHFEcutStep(const Char_t *name):
  TNamed(name, ""),
  fCuts(NULL)
{
  //
  // Default Constructor
  //
  fCuts = new TObjArray;
}

//__________________________________________________________________
AliHFEcutStep::AliHFEcutStep(const AliHFEcutStep &o):
  TNamed(o),
  fCuts(NULL)
{
  //
  // Copy constructor
  //
  o.Copy(*this);
}

//__________________________________________________________________
AliHFEcutStep &AliHFEcutStep::operator=(const AliHFEcutStep &o){
  //
  // Assignment operator
  //
  if(this != &o)
    o.Copy(*this);
  return *this;
}
 
//__________________________________________________________________
AliHFEcutStep::~AliHFEcutStep(){
  //
  // destructor
  //
  delete fCuts;
}

//__________________________________________________________________
void AliHFEcutStep::Copy(TObject &o) const{
  //
  // Copy into content into object o
  //
  TNamed::Copy(o);
  AliHFEcutStep &target = dynamic_cast<AliHFEcutStep &>(o);

  // Make copy
  target.fCuts = dynamic_cast<TObjArray *>(fCuts->Clone());
}

//__________________________________________________________________
Bool_t AliHFEcutStep::IsSelected(TObject *o){
  //
  // Filter tracks in the given cut step
  // Apply all cut objects
  //
  AliDebug(1, Form("Cut Step %s: Number of cut objects: %d", GetName(), fCuts->GetEntriesFast()));
  if(!fCuts->GetEntriesFast()) return kTRUE;
  Bool_t isSelected = kTRUE;
  AliAnalysisCuts *cuts = NULL;
  for(Int_t iCut = 0; iCut < fCuts->GetEntriesFast(); iCut++){
    cuts = dynamic_cast<AliAnalysisCuts *>(fCuts->UncheckedAt(iCut));
    if(cuts && (!cuts->IsSelected(o))) isSelected = kFALSE;
  }
  AliDebug(1, Form("Accepted: %s", isSelected ? "yes" : "no"));
  return isSelected;
} 

//__________________________________________________________________
AliAnalysisCuts *AliHFEcutStep::GetCut(const Char_t *cutName){
  //
  // return cut object
  //
  return dynamic_cast<AliAnalysisCuts *>(fCuts->FindObject(cutName));
}

//__________________________________________________________________
void AliHFEcutStep::AddCut(AliAnalysisCuts *cut){
  //
  // Add cut object to the cut step
  //
  fCuts->Add(cut);
}

//__________________________________________________________________
void AliHFEcutStep::SetMC(const AliMCEvent *mc){
  //
  // Set MC information to the cuts in the cut step
  //
  AliCFCutBase *cfc = NULL;
  for(Int_t icut = 0; icut < fCuts->GetEntriesFast(); icut++){
    if((cfc = dynamic_cast<AliCFCutBase *>(fCuts->UncheckedAt(icut)))) cfc->SetMCEventInfo(mc);
  }
}

//__________________________________________________________________
void AliHFEcutStep::SetRecEvent(const AliVEvent *rec){
  //
  // Publish rec event to the cut step
  //
  AliCFCutBase *cfc = NULL;
  for(Int_t icut = 0; icut < fCuts->GetEntriesFast(); icut++){
    if((cfc = dynamic_cast<AliCFCutBase *>(fCuts->UncheckedAt(icut)))) cfc->SetRecEventInfo(rec);
  }
}

