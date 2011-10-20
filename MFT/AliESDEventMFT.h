#ifndef AliESDEventMFT_H
#define AliESDEventMFT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

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

//====================================================================================================================================================

class AliMuonForwardTrack;

class AliESDEventMFT: public AliESDEvent { 

public:

  AliESDEventMFT();
  AliESDEventMFT(AliESDEvent &esdEvent);

  AliESDEventMFT(const AliESDEventMFT&);
  AliESDEventMFT &operator=(const AliESDEventMFT&);

  virtual ~AliESDEventMFT();

  AliMuonForwardTrack *GetMuonForwardTrack(Int_t i) const { 
    return (AliMuonForwardTrack*)(fMuonForwardTracks?fMuonForwardTracks->UncheckedAt(i):0x0); 
  }

  void AddMuonForwardTrack(const AliMuonForwardTrack *muonForwardTrack);

  Int_t GetNMuonForwardTracks() const { return fMuonForwardTracks?fMuonForwardTracks->GetEntriesFast():0; }

private:
 
  TClonesArray *fMuonForwardTracks;       // array of AliMuonForwardTrack

  ClassDef(AliESDEventMFT, 1) 

};

//====================================================================================================================================================

#endif
