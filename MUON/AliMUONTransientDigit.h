#ifndef ALIMUONTRANSIENTDIGIT_H
#define ALIMUONTRANSIENTDIGIT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

/// \ingroup base
/// \class AliMUONTransientDigit
/// \brief MUON transient digit

#include "AliMUONDigit.h"

class TObjArray;

class AliMUONTransientDigit : public AliMUONDigit 
{
  public:
    AliMUONTransientDigit();
    AliMUONTransientDigit(Int_t rpad, Int_t *digits);
    virtual ~AliMUONTransientDigit();

    Int_t Chamber() const {return fChamber;}
    Int_t GetNTracks() const {return fTrackList->GetEntriesFast();}
    Int_t GetTrack(Int_t i) const;
    Int_t GetCharge(Int_t i) const;
    void AddToTrackList(Int_t track, Int_t charge);
    void UpdateTrackList(Int_t track, Int_t charge);
    
  protected:
    AliMUONTransientDigit(const AliMUONTransientDigit& digit);
    AliMUONTransientDigit & operator =(const AliMUONTransientDigit & rhs);

    Int_t          fChamber;       ///< chamber number of pad
    TObjArray     *fTrackList;     ///< List of tracks contributing

  ClassDef(AliMUONTransientDigit,1)  // Transient digit for MUON
};
#endif

