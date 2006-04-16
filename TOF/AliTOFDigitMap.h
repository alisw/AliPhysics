#ifndef ALITOFDIGITMAP_H
#define ALITOFDIGITMAP_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


////////////////////////////////////////////////////////////////////////
//
// AliTOFDigitMap class
//
// digitmap enables fast check if the pad was already hit
//
// Author: F. Pierella | pierella@bo.infn.it
//
////////////////////////////////////////////////////////////////////////

#include "AliHitMap.h"
#include "TObject.h"
#include "AliTOFGeometry.h"

class TClonesArray;
class AliTOFGeometry;

class AliTOFDigitMap : public TObject
{
 public:
    AliTOFDigitMap();
    AliTOFDigitMap(TClonesArray *dig, AliTOFGeometry *tofGeom);
    AliTOFDigitMap(const AliTOFDigitMap & digitMap);
    
    virtual ~AliTOFDigitMap();
    // Clear the hit map
    virtual  void  Clear(const char *opt = "");
    // Set a single hit
    virtual  void  SetHit(Int_t *vol, Int_t idigit);
    virtual  void  SetHit(Int_t *vol);
    // Get index of hit in the list of digits
    virtual Int_t  GetHitIndex(Int_t *vol) const;
    // Get pointer to digit
    virtual TObject*  GetHit(Int_t *vol) const;
    // Test hit status
    virtual FlagType TestHit(Int_t *vol) const;
    // Assignment operator
    AliTOFDigitMap& operator = (const AliTOFDigitMap& rhs);
    
 private:
    // Check index
    Int_t CheckedIndex(Int_t *vol) const;
 private:
    Int_t fNSector;                       // Number of sectors
    Int_t fNplate;                        // Number of plates
    Int_t fNstrip;                        // Maximum number of strips
    Int_t fNpx;                           // Number of pads in x
    Int_t fNpz;                           // Number of pads in z

    TClonesArray *fDigits;                // Pointer to sdigits
    Int_t fMaxIndex;                      // maximum index in hit map
    Int_t *fDigitMap;                     // ! [fMaxIndex]         

    AliTOFGeometry *fTOFGeometry;         // Pointer to the TOF geometry

    ClassDef(AliTOFDigitMap,0) // Implements DigitMap as a 1-dim array
};
#endif	
