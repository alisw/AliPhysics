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

///////////////////////////////////////////////////////////////////////////
//          ----   CORRECTION FRAMEWORK   ----
// Class to cut on the number of AliTrackReference's 
// for each detector. Applies on pair of tracks (AliCFPair)
///////////////////////////////////////////////////////////////////////////
// author : R. Vernet (renaud.vernet@cern.ch)
///////////////////////////////////////////////////////////////////////////

#include "AliMCParticle.h"
#include "AliCFPairAcceptanceCuts.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"

ClassImp(AliCFPairAcceptanceCuts)

//______________________________
AliCFPairAcceptanceCuts::AliCFPairAcceptanceCuts() : 
  AliCFCutBase()
{
  //
  //ctor
  //
  fCutNeg = new AliCFAcceptanceCuts();
  fCutPos = new AliCFAcceptanceCuts();
}

//______________________________
AliCFPairAcceptanceCuts::AliCFPairAcceptanceCuts(const Char_t* name, const Char_t* title) : 
  AliCFCutBase(name,title)
{
  //
  //ctor
  //
  fCutNeg = new AliCFAcceptanceCuts(name,title);
  fCutPos = new AliCFAcceptanceCuts(name,title);
}

//______________________________
AliCFPairAcceptanceCuts::AliCFPairAcceptanceCuts(const AliCFPairAcceptanceCuts& c) : 
  AliCFCutBase(c),
  fCutNeg(c.fCutNeg),
  fCutPos(c.fCutPos)
{
  //
  //copy ctor
  //
}

//______________________________
AliCFPairAcceptanceCuts& AliCFPairAcceptanceCuts::operator=(const AliCFPairAcceptanceCuts& c)
{
  //
  // Assignment operator
  //
  if (this != &c) {
    AliCFCutBase::operator=(c) ;
    fCutNeg = c.fCutNeg ;
    fCutPos = c.fCutPos ;
  }
  return *this ;
}

//______________________________
Bool_t AliCFPairAcceptanceCuts::IsSelected(TObject* obj) {
  //
  // checks the number of track references associated to 'obj'
  // 'obj' must be an AliMCParticle
  //

  if (!obj) return kFALSE ;
  TString className(obj->ClassName());
  if (className.CompareTo("AliMCParticle") != 0) {
    Error("IsSelected","obj must point to a AliMCParticle !");
    return kFALSE ;
  }

  TParticle* part = ((AliMCParticle*)obj)->Particle() ;
  if (!part || part->GetNDaughters()!=2) return kFALSE ;
  Int_t lab0 = part->GetDaughter(0);
  Int_t lab1 = part->GetDaughter(1);
  
  AliMCParticle* negDaughter = fMCInfo->MCEvent()->GetTrack(lab0) ;
  AliMCParticle* posDaughter = fMCInfo->MCEvent()->GetTrack(lab1) ;

  if (!fCutNeg->IsSelected(negDaughter)) return kFALSE;
  if (!fCutPos->IsSelected(posDaughter)) return kFALSE; 

  return kTRUE ;
}

//______________________________
void AliCFPairAcceptanceCuts::SetEvtInfo(TObject* mcInfo) {
  //
  // Sets pointer to MC event information (AliMCEventHandler)
  //

  if (!mcInfo) {
    Error("SetEvtInfo","Pointer to MC Event Handler is null !");
    return;
  }
  
  TString className(mcInfo->ClassName());
  if (className.CompareTo("AliMCEventHandler") != 0) {
    Error("SetEvtInfo","argument must point to an AliMCEventHandler !");
    return ;
  }
  
  fMCInfo = (AliMCEventHandler*) mcInfo ;
}
