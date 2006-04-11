#ifndef ALIMUONDIGITMAPA1_H
#define ALIMUONDIGITMAPA1_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

/// \ingroup base
/// \class AliMUONDigitMapA1
/// \brief Implements cluster Map as a 1-dim array
/// 
/// \author Christian Finck

#include "AliHitMap.h"

#include <TObject.h>

class TObjArray;

class AliMUONDigitMapA1 : public TObject
{
 public:
    AliMUONDigitMapA1();
    AliMUONDigitMapA1(Int_t idDE, Int_t npx, Int_t npy);
    virtual ~AliMUONDigitMapA1();

    // Fill hits from list of digits into hit map
    virtual  void  FillHits(TObjArray* dig);
    // Clear the hit map
    virtual  void  Clear(const char *opt = "");
    // Set a single hit
    virtual  void  SetHit(Int_t ix, Int_t iy, Int_t idigit);
    // Delete a single hit
    virtual  void  DeleteHit(Int_t ix, Int_t iy);
    // Get index of hit in the list of digits
    virtual Int_t  GetHitIndex(Int_t ix, Int_t iy) const;
    // Get pointer to digit
    virtual TObject*  GetHit(Int_t ix, Int_t iy) const;
    // Flag a hit as used
    virtual  void  FlagHit(Int_t ix, Int_t iy);
    // Test hit status
    virtual FlagType TestHit(Int_t ix, Int_t iy);

 protected:
    AliMUONDigitMapA1(const AliMUONDigitMapA1 & hitMap);
    // Assignment operator
    AliMUONDigitMapA1& operator = (const AliMUONDigitMapA1& rhs);
    
 private:
    // Check index
    Int_t CheckedIndex(Int_t ix, Int_t iy) const;
 private:
    Int_t fIdDE;                          // id DE
    Int_t fNpx;                           // Maximum number of pads in x
    Int_t fNpy;                           // Maximum number of pads in y
    TObjArray *fDigits;                   // Pointer to digits
    Int_t fMaxIndex;                      // maximum index in hit map
    Int_t *fHitMap;                       // ! [fMaxIndex]         

    ClassDef(AliMUONDigitMapA1,0) // Implements HitMap as a 1-dim array
};
#endif	
