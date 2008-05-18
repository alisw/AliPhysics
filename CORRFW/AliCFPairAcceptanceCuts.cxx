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
#include "AliMCEvent.h"
#include "TBits.h"
#include "AliLog.h"

ClassImp(AliCFPairAcceptanceCuts)

//______________________________
AliCFPairAcceptanceCuts::AliCFPairAcceptanceCuts() : 
  AliCFCutBase(),
  fMCInfo(0x0),
  fCutNeg(new AliCFAcceptanceCuts()),
  fCutPos(new AliCFAcceptanceCuts()),
  fBitmap(new TBits(0))
{
  //
  //Default Constructor
  //
}

//______________________________
AliCFPairAcceptanceCuts::AliCFPairAcceptanceCuts(const Char_t* name, const Char_t* title) : 
  AliCFCutBase(name,title),
  fMCInfo(0x0),
  fCutNeg(new AliCFAcceptanceCuts(name,title)),
  fCutPos(new AliCFAcceptanceCuts(name,title)),
  fBitmap(new TBits(0))
{
  //
  //Named Constructor
  //
}

//______________________________
AliCFPairAcceptanceCuts::AliCFPairAcceptanceCuts(const AliCFPairAcceptanceCuts& c) : 
  AliCFCutBase(c),
  fMCInfo(c.fMCInfo),
  fCutNeg(c.fCutNeg),
  fCutPos(c.fCutPos),
  fBitmap(c.fBitmap)
{
  //
  //Copy Constructor
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
    fMCInfo = c.fMCInfo ;
    fCutNeg = c.fCutNeg ;
    fCutPos = c.fCutPos ;
    fBitmap = c.fBitmap ;
  }
  return *this ;
}

//__________________________________________________________
Bool_t AliCFPairAcceptanceCuts::IsSelected(TObject* obj) {
  //
  // checks the number of track references associated to 'obj'
  // 'obj' must be an AliMCParticle
  //
  //
  // check if selections on 'obj' are passed
  // 'obj' must be an AliMCParticle
  //
  
  SelectionBitMap(obj);

  //   if (fIsQAOn) FillHistograms(obj,kFALSE);
  Bool_t isSelected = kTRUE;

  for (UInt_t icut=0; icut<fBitmap->GetNbits(); icut++) {
    if (!fBitmap->TestBitNumber(icut)) {
      isSelected = kFALSE;
      break;
    }
  }  

  if (!isSelected) return kFALSE ;
  //   if (fIsQAOn) FillHistograms(obj,kTRUE);
  return kTRUE;
}

//__________________________________________________________
void AliCFPairAcceptanceCuts::SelectionBitMap(TObject* obj) 
{
  //
  // test if the track passes the single cuts
  // and store the information in a bitmap
  //

  for (UInt_t i=0; i<kNCuts; i++) fBitmap->SetBitNumber(i,kFALSE);

  if (!obj) return;
  TString className(obj->ClassName());
  if (className.CompareTo("AliMCParticle") != 0) {
    AliError("obj must point to an AliMCParticle !");
    return ;
  }

  TParticle* part = (dynamic_cast<AliMCParticle*>(obj))->Particle() ;
  if (!part || part->GetNDaughters() !=2) return ;

  Int_t lab0 = part->GetDaughter(0);
  Int_t lab1 = part->GetDaughter(1);
  AliMCParticle* negDaughter = fMCInfo->GetTrack(lab0) ;
  AliMCParticle* posDaughter = fMCInfo->GetTrack(lab1) ;

  Int_t iCutBit = 0;

  if (fCutNeg->IsSelected(negDaughter)) fBitmap->SetBitNumber(iCutBit,kTRUE);
  iCutBit++;
  
  if (fCutPos->IsSelected(posDaughter)) fBitmap->SetBitNumber(iCutBit,kTRUE);
}

//______________________________
void AliCFPairAcceptanceCuts::SetEvtInfo(TObject* mcInfo) {
  //
  // Sets pointer to MC event information (AliMCEvent)
  //

  if (!mcInfo) {
    Error("SetEvtInfo","Pointer to MC Event is null !");
    return;
  }
  
  TString className(mcInfo->ClassName());
  if (className.CompareTo("AliMCEvent") != 0) {
    Error("SetEvtInfo","argument must point to an AliMCEvent !");
    return ;
  }
  
  fMCInfo = (AliMCEvent*) mcInfo ;
}
