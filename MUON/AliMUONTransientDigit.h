#ifndef ALIMUONTRANSIENTDIGIT_H
#define ALIMUONTRANSIENTDIGIT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliMUONDigit.h"
class TObjArray;


class AliMUONTransientDigit : public AliMUONDigit {
 public:
    Int_t          fChamber;       // chamber number of pad
    TObjArray     *fTrackList;     // List of tracks contributing
 public:
    AliMUONTransientDigit() {fTrackList=0;}
    AliMUONTransientDigit(const AliMUONTransientDigit& digit);
    AliMUONTransientDigit(Int_t rpad, Int_t *digits);
    virtual ~AliMUONTransientDigit();
    TObjArray  *TrackList()   {return fTrackList;}
    AliMUONTransientDigit & operator =(const AliMUONTransientDigit & rhs);
    
    ClassDef(AliMUONTransientDigit,1)  // Transient digit for set MUON
};
#endif

