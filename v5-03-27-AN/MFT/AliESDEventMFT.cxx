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

//====================================================================================================================================================
//
//      ESD Event with MUON+MFT muon tracks (AliMuonForwardTrack)
//
//      Contact author: antonio.uras@cern.ch
//
//====================================================================================================================================================

#include "AliESDEvent.h"
#include "TClonesArray.h"
#include "AliMuonForwardTrack.h"
#include "AliESDEventMFT.h"

ClassImp(AliESDEventMFT)

//====================================================================================================================================================

AliESDEventMFT::AliESDEventMFT():
  AliESDEvent(), 
  fMuonForwardTracks(0x0)
{

  // default constructor 
  fMuonForwardTracks = new TClonesArray("AliMuonForwardTrack", 0);
  AddObject(fMuonForwardTracks);

}

//====================================================================================================================================================

AliESDEventMFT::AliESDEventMFT(AliESDEvent &esdEvent):
  AliESDEvent(esdEvent), 
  fMuonForwardTracks(0x0)
{

  AliDebug(1, "building array of muon tracks");
  fMuonForwardTracks = new TClonesArray("AliMuonForwardTrack");
  AliDebug(1, "adding array of muon tracks to list");
  AddObject(fMuonForwardTracks);
  AliDebug(1, "event created!");
    
}

//====================================================================================================================================================

AliESDEventMFT::AliESDEventMFT(const AliESDEventMFT &esdEventMFT): 
  AliESDEvent(esdEventMFT),
  fMuonForwardTracks(esdEventMFT.fMuonForwardTracks)
{

  // copy constructor
  AddObject(fMuonForwardTracks);
  
}

//====================================================================================================================================================

AliESDEventMFT& AliESDEventMFT::operator=(const AliESDEventMFT &esdEventMFT) {

  // Asignment operator

  // check assignement to self
  if (this == &esdEventMFT) return *this;

  // base class assignement
  AliESDEvent::operator=(esdEventMFT);
  
  // clear memory
  Clear();
  
  fMuonForwardTracks = esdEventMFT.fMuonForwardTracks;

  return *this;

}

//====================================================================================================================================================

AliESDEventMFT::~AliESDEventMFT() {

  // destructor

  if (fMuonForwardTracks) {
    //    fMuonForwardTracks->Delete();
    delete fMuonForwardTracks;
  }

}

//====================================================================================================================================================

void AliESDEventMFT::AddMuonForwardTrack(const AliMuonForwardTrack *muonForwardTrack) {

    TClonesArray &myMuonForwardTracks = *fMuonForwardTracks;
    new (myMuonForwardTracks[fMuonForwardTracks->GetEntriesFast()]) AliMuonForwardTrack(*muonForwardTrack);

}

//====================================================================================================================================================
