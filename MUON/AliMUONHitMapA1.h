#ifndef ALIMUONHITMAPA1_H
#define ALIMUONHITMAPA1_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliMUONHitMap.h"
class AliMUONSegmentation;
class TObjArray;



class AliMUONHitMapA1 :
public AliMUONHitMap 
{
 public:
    AliMUONHitMapA1(AliMUONSegmentation *seg, TObjArray *dig);
    AliMUONHitMapA1(const AliMUONHitMapA1 & hitMap);
    
    virtual ~AliMUONHitMapA1();
    // Fill hits from list of digits into hit map
    virtual  void  FillHits();
    // Clear the hit map
    virtual  void  Clear();
    // Set a single hit
    virtual  void  SetHit(Int_t ix, Int_t iy, Int_t idigit);
    // Delete a single hit
    virtual  void  DeleteHit(Int_t ix, Int_t iy);
    // Get index of hit in the list of digits
    virtual Int_t  GetHitIndex(Int_t ix, Int_t iy);
    // Get pointer to digit
    virtual TObject*  GetHit(Int_t ix, Int_t iy);
    // Flag a hit as used
    virtual  void  FlagHit(Int_t ix, Int_t iy);
    // Test hit status
    virtual FlagType TestHit(Int_t ix, Int_t iy);
    // Assignment operator
    AliMUONHitMapA1& operator = (const AliMUONHitMapA1& rhs);
    
 private:
    // Check index
    Int_t CheckedIndex(Int_t ix, Int_t iy);
 private:
    AliMUONSegmentation *fSegmentation;   // Chamber segmentation
    Int_t fNpx;                           // Maximum number of pads in x
    Int_t fNpy;                           // Maximum number of pads in y
    TObjArray *fDigits;                   // Pointer to digits
    Int_t fNdigits;                       // Number of digits
    Int_t *fHitMap;                       // [fMaxIndex]     
    Int_t fMaxIndex;                      // maximum index in hit map
    

    ClassDef(AliMUONHitMapA1,1) // Implements HitMap as a 1-dim array
};
#endif	
