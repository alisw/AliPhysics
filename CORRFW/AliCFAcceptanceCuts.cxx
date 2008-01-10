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
// AliCFAcceptanceCuts implementation
// Class to cut on the number of AliTrackReference's 
// for each detector
///////////////////////////////////////////////////////////////////////////
// author : R. Vernet (renaud.vernet@cern.ch)
///////////////////////////////////////////////////////////////////////////

#include "AliLog.h"
#include "AliMCParticle.h"
#include "AliCFAcceptanceCuts.h"

ClassImp(AliCFAcceptanceCuts)

//______________________________
AliCFAcceptanceCuts::AliCFAcceptanceCuts() : 
  AliCFCutBase(),
  fMCInfo(0x0),
  fMinNHitITS(0),
  fMinNHitTPC(0),
  fMinNHitTRD(0),
  fMinNHitTOF(0),
  fMinNHitMUON(0)

{
  //
  //ctor
  //
}

//______________________________
AliCFAcceptanceCuts::AliCFAcceptanceCuts(const Char_t* name, const Char_t* title) : 
  AliCFCutBase(name,title),
  fMCInfo(0x0),
  fMinNHitITS(0),
  fMinNHitTPC(0),
  fMinNHitTRD(0),
  fMinNHitTOF(0),
  fMinNHitMUON(0)
{
  //
  //ctor
  //
}

//______________________________
AliCFAcceptanceCuts::AliCFAcceptanceCuts(const AliCFAcceptanceCuts& c) : 
  AliCFCutBase(c),
  fMCInfo(c.fMCInfo),
  fMinNHitITS(c.fMinNHitITS),
  fMinNHitTPC(c.fMinNHitTPC),
  fMinNHitTRD(c.fMinNHitTRD),
  fMinNHitTOF(c.fMinNHitTOF),
  fMinNHitMUON(c.fMinNHitMUON)
{
  //
  //copy ctor
  //
}

//______________________________
AliCFAcceptanceCuts& AliCFAcceptanceCuts::operator=(const AliCFAcceptanceCuts& c)
{
  //
  // Assignment operator
  //
  if (this != &c) {
    AliCFCutBase::operator=(c) ;
    fMCInfo=c.fMCInfo;
    fMinNHitITS=c.fMinNHitITS;
    fMinNHitTPC=c.fMinNHitTPC;
    fMinNHitTRD=c.fMinNHitTRD;
    fMinNHitTOF=c.fMinNHitTOF;
    fMinNHitMUON=c.fMinNHitMUON;
  }
  return *this ;
}

//______________________________
Bool_t AliCFAcceptanceCuts::IsSelected(TObject* obj) {
  //
  // checks the number of track references associated to 'obj'
  // 'obj' must be an AliMCParticle
  //

  if (!obj) return kFALSE ;

  TString className(obj->ClassName());
  if (className.CompareTo("AliMCParticle") != 0) {
    AliError("obj must point to a AliMCParticle !");
    return kFALSE ;
  }

  AliMCParticle* part = (AliMCParticle*) obj ;
  if(!part) return kFALSE;

  Int_t nHitsITS=0, nHitsTPC=0, nHitsTRD=0, nHitsTOF=0, nHitsMUON=0 ;
  for (Int_t iTrackRef=0; iTrackRef<part->GetNumberOfTrackReferences(); iTrackRef++) {
    AliTrackReference * trackRef = part->GetTrackReference(iTrackRef);
    if(trackRef){
      Int_t detectorId = trackRef->DetectorId();
      switch(detectorId) {
      case AliTrackReference::kITS  : nHitsITS++  ; break ;
      case AliTrackReference::kTPC  : nHitsTPC++  ; break ;
      case AliTrackReference::kTRD  : nHitsTRD++  ; break ;
      case AliTrackReference::kTOF  : nHitsTOF++  ; break ;
      case AliTrackReference::kMUON : nHitsMUON++ ; break ;
      default : break ;
      }
    }
  }
  
  if (nHitsITS  < fMinNHitITS  ) return kFALSE;
  if (nHitsTPC  < fMinNHitTPC  ) return kFALSE;
  if (nHitsTRD  < fMinNHitTRD  ) return kFALSE;
  if (nHitsTOF  < fMinNHitTOF  ) return kFALSE;
  if (nHitsMUON < fMinNHitMUON ) return kFALSE;


  return kTRUE ;
}


void AliCFAcceptanceCuts::SetEvtInfo(TObject* mcInfo) {
  //
  // Sets pointer to MC event information (AliMCEventHandler)
  //

  if (!mcInfo) {
    AliError("Pointer to MC Event Handler is null !");
    return;
  }
  
  TString className(mcInfo->ClassName());
  if (className.CompareTo("AliMCEventHandler") != 0) {
    AliError("argument must point to an AliMCEventHandler !");
    return ;
  }
  
  fMCInfo = (AliMCEventHandler*) mcInfo ;
}
