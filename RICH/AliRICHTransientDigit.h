#ifndef ALIRICHTRANSIENTDIGIT_H
#define ALIRICHTRANSIENTDIGIT_H


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id $ */
#include "AliRICHDigit.h"

class AliRICHTransientDigit : public AliRICHDigit {
 public:
    Int_t          fRpad;          // r_pos of pad
    Int_t          fChamber;       // chamber number of pad
    TObjArray     *fTrackList;     // list of tracks
 public:
    AliRICHTransientDigit() {fTrackList=0;}
    AliRICHTransientDigit(Int_t ich, Int_t *digits);
    virtual ~AliRICHTransientDigit() {}
    
    TObjArray  *TrackList()   {return fTrackList;}
    
    ClassDef(AliRICHTransientDigit,1)  //Digits for set:RICH
};

#endif
