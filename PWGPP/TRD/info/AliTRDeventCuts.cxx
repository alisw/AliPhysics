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

#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"

#include "AliTRDeventInfo.h"
#include "AliTRDeventCuts.h"

ClassImp(AliTRDeventCuts)

//______________________________________________________________
AliTRDeventCuts::AliTRDeventCuts()
  :TNamed("trdEventCuts", "")
  ,fTriggerNames(NULL)
  ,fBunches(NULL)
  ,fEventType(7)
  ,fVertexN(1)
  ,fVertexZ(15.)
{
  //
  // Dummy Constructor
  //
  
}

//______________________________________________________________
AliTRDeventCuts::AliTRDeventCuts(const Char_t *name)
  :TNamed(name, "")
  ,fTriggerNames(NULL)
  ,fBunches(NULL)
  ,fEventType(7)
  ,fVertexN(1)
  ,fVertexZ(15.)
{
  //
  // Default Constructor
  //
}

//______________________________________________________________
AliTRDeventCuts::AliTRDeventCuts(const AliTRDeventCuts &ref)
  :TNamed((TNamed&)ref)
  ,fTriggerNames(NULL)
  ,fBunches(NULL)
  ,fEventType(ref.fEventType)
  ,fVertexN(ref.fVertexN)
  ,fVertexZ(ref.fVertexZ)
{
// Copy constructor
  if(ref.fTriggerNames){
    for(Int_t it(0); it<ref.fTriggerNames->GetEntriesFast(); it++) AddTrigger(((TObjString*)(*ref.fTriggerNames)[it])->GetName());
  }
  if(ref.fBunches) SetBunchSelection(AliTRDeventInfo::kLHCbunches, ref.fBunches);
}

//______________________________________________________________
AliTRDeventCuts::~AliTRDeventCuts()
{
// Destructor

  if(fTriggerNames) fTriggerNames->Delete();
  delete fTriggerNames;
  if(fBunches) delete [] fBunches;
}

//______________________________________________________________
Bool_t AliTRDeventCuts::IsSelected(AliESDEvent *ev, Bool_t col)
{
// Apply cuts

  Bool_t select = kTRUE;
  if(fTriggerNames){
    Bool_t passTrigger = kFALSE; 
    TString trgname = ev->GetFiredTriggerClasses();
    TObjArray *triggers = trgname.Tokenize(" ");
    TIterator *trgIter = triggers->MakeIterator(); 
    TObjString *trg(NULL);
    while((trg = dynamic_cast<TObjString *>(trgIter->Next())))
      passTrigger = passTrigger || CheckTrigger(trg->String().Data());
    delete triggers; delete trgIter;
    select = select && passTrigger;
  }
  if(!select){
    AliDebug(1, Form("Reject Ev[%d] for trigger[%s]", ev->GetEventNumberInFile(), ev->GetFiredTriggerClasses().Data()));
    return select;
  }
  // select only physical events
  select = select && (ev->GetEventType() == fEventType);
  if(!select){
    AliDebug(1, Form("Reject Ev[%d] for EvType[%d]", ev->GetEventNumberInFile(), ev->GetEventType()));
    return select;
  }

  if(!col) return select;

  // vertex selection
  const AliESDVertex *primVtx = ev->GetPrimaryVertex();
  if(fVertexN > 0)
    select = select && (primVtx && primVtx->GetNContributors() >= fVertexN);
  if(fVertexZ >= 0.)
    select = select && (primVtx && TMath::Abs(primVtx->GetZ()) <= fVertexZ);
  if(!select){
    AliDebug(1, Form("Reject Ev[%d] for Vertex[%p][%d %6.2f]", ev->GetEventNumberInFile(), (void*)primVtx, primVtx?primVtx->GetNContributors():0, primVtx?TMath::Abs(primVtx->GetZ()):0));
    return select;
  }

  // bunch cross selection
  if(fBunches){
    Int_t evBC(ev->GetBunchCrossNumber()), ibc(0);
    Bool_t kFOUND(kFALSE);
    while(fBunches[ibc]>0){
      if(evBC==fBunches[ibc]){
        kFOUND = kTRUE;
        break;
      }
      ibc++;
    }
    select = select && kFOUND;
  }
  if(!select){
    AliDebug(1, Form("Reject Ev[%d] for BunchCross[%d]", ev->GetEventNumberInFile(), ev->GetBunchCrossNumber()));
    return select;
  }
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

  if(!fTriggerNames) return kFALSE;
  Bool_t kExists(kFALSE);
  for(Int_t it(0); it<fTriggerNames->GetEntriesFast(); it++){
    if(((TObjString*)(*fTriggerNames)[it])->String().CompareTo(name)==0){
      kExists = kTRUE;
      break;
    }
  }
  return kExists;
}


//______________________________________________________________
void AliTRDeventCuts::Print(Option_t */*opt*/) const
{
// Print content of event cuts
  printf("Event Type       : %2d\n", fEventType);
  printf("Vertex  selection: N[%2d] Z[cm]=%6.2f\n", fVertexN, fVertexZ);
  if(fTriggerNames){
    printf("Trigger selection: ");
    for(Int_t it(0); it<fTriggerNames->GetEntriesFast(); it++) printf("\"%s\" ", ((TObjString*)(*fTriggerNames)[it])->GetName());
    printf("\n");
  }
  if(fBunches){
    printf("Bunches selection: ");
    for(Int_t ibc(0); ibc<AliTRDeventInfo::kLHCbunches; ibc++){
      if(fBunches[ibc]<0) break;
      printf("%4d ", fBunches[ibc]);
    }
    printf("\n");
  }
}

//______________________________________________________________
void AliTRDeventCuts::SetBunchSelection(Int_t n, Int_t bunches[])
{
// Set Bunch selection for run
  if(!fBunches) fBunches = new Int_t[AliTRDeventInfo::kLHCbunches];
  for(Int_t ibc(0); ibc<AliTRDeventInfo::kLHCbunches; ibc++) fBunches[ibc] = ibc<n?bunches[ibc]:-1;
}

