#ifndef ALIMUONTRANSIENTDIGIT_H
#define ALIMUONTRANSIENTDIGIT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliMUONDigit.h"
#include "TObjArray.h"

class AliMUONTransientDigit : public AliMUONDigit 
{
  public:
    AliMUONTransientDigit() {fTrackList=0;}
    AliMUONTransientDigit(const AliMUONTransientDigit& digit);
    AliMUONTransientDigit(Int_t rpad, Int_t *digits);
    virtual ~AliMUONTransientDigit();

    Int_t Chamber() const {return fChamber;}
    Int_t GetNTracks() const {return fTrackList->GetEntriesFast();}
    Int_t GetTrack(Int_t i) const;
    Int_t GetCharge(Int_t i) const;
    void AddToTrackList(Int_t track, Int_t charge);
    void UpdateTrackList(Int_t track, Int_t charge);
    AliMUONTransientDigit & operator =(const AliMUONTransientDigit & rhs);
    
  protected:
    Int_t          fChamber;       // chamber number of pad
    TObjArray     *fTrackList;     // List of tracks contributing

  ClassDef(AliMUONTransientDigit,1)  // Transient digit for MUON
};
#endif

