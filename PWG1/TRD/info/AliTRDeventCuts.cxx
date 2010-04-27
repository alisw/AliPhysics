/*************************************************************************
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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
* provided "as is" without express or implied warranty.                  * **************************************************************************/

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// Event cut class for the TRD Performance Train                          //
//                                                                        //
// Encapsulation of events cuts for usage with TRD performance train      //
// - reconstructed vertex                                                 //
// - vertex Z position                                                    //
// - vertex number of contributors                                        //
// - trigger selection list                                               //
//                                                                        //
// *) see constructor for default values                                  //
//                                                                        //
// author                                                                 //
// Markus Fasel <m.fasel@gsi.de>                                          //
//                                                                        //
////////////////////////////////////////////////////////////////////////////


#include <TIterator.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TString.h>

#include "AliESDEvent.h"
#include "AliESDVertex.h"

#include "AliTRDeventCuts.h"

ClassImp(AliTRDeventCuts)

//______________________________________________________________
AliTRDeventCuts::AliTRDeventCuts():
  TNamed("trdEventCuts", ""),
  fTriggerNames(NULL),
  fVertexN(1),
  fVertexZ(15.)
{
  //
  // Dummy Constructor
  //
}

//______________________________________________________________
AliTRDeventCuts::AliTRDeventCuts(const Char_t *name):
  TNamed(name, ""),
  fTriggerNames(NULL),
  fVertexN(1),
  fVertexZ(15.)
{
  //
  // Default Constructor
  //
}

//______________________________________________________________
AliTRDeventCuts::~AliTRDeventCuts()
{
// Destructor

  if(fTriggerNames) fTriggerNames->Delete();
  delete fTriggerNames;
}

//______________________________________________________________
Bool_t AliTRDeventCuts::IsSelected(AliESDEvent *ev, Bool_t col)
{
// Apply cuts

  Bool_t select = kTRUE;
  if(fTriggerNames){
    Bool_t passTrigger = kFALSE; 
    TString trgname = ev->GetFiredTriggerClasses();
    TObjArray *triggers = trgname.Tokenize("");
    TIterator *trgIter = triggers->MakeIterator(); 
    TObjString *trg = NULL;
    while((trg = dynamic_cast<TObjString *>(trgIter->Next())))
      passTrigger = passTrigger || CheckTrigger(trg->String().Data());
    delete triggers; delete trgIter;
    select = select && passTrigger;
  }
  // select only physical events
  select = select && (ev->GetEventType() == 7);

  if(!col) return select;

  const AliESDVertex *primVtx = ev->GetPrimaryVertex();
  if(fVertexN > 0)
    select = select && (primVtx && primVtx->GetNContributors() >= fVertexN);
  if(fVertexZ >= 0.)
    select = select && (primVtx && TMath::Abs(primVtx->GetZ()) <= fVertexZ);
  return select;
}

//______________________________________________________________
void AliTRDeventCuts::AddTrigger(const Char_t *name)
{
// Add trigger name according to the logbook

  if(!fTriggerNames) fTriggerNames = new TObjArray;
  if(CheckTrigger(name)) return;
  fTriggerNames->Add(new TObjString(name));
}

//______________________________________________________________
Bool_t AliTRDeventCuts::CheckTrigger(const Char_t *name)
{
// check if trigger with id "name" is on the accepted trigger list

  TObjString *trg = NULL;
  TIterator *it = fTriggerNames->MakeIterator();
  Bool_t triggerExists = kFALSE;
  while((trg = dynamic_cast<TObjString *>(it))){
    if(!trg->String().CompareTo(name)) triggerExists = kTRUE;
  }
  delete it;
  return triggerExists;
}
