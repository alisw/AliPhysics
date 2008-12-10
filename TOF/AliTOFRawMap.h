#ifndef ALITOFRAWMAP_H
#define ALITOFRAWMAP_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//                                            //
//   AliTOFRawMap class                       //
//                                            //
//   It enables fast check                    //
//   if the TDC channel was already engaged   //
//   for a measurement                        //
//                                            //
////////////////////////////////////////////////

#include "TObject.h"

#include "AliHitMap.h"

class TClonesArray;

class AliTOFRawMap : public TObject
{
 public:
    AliTOFRawMap();
    AliTOFRawMap(TClonesArray *sdig);
    AliTOFRawMap(const AliTOFRawMap & rawMap);
    AliTOFRawMap &operator=(const AliTOFRawMap & rawMap);
    
    virtual ~AliTOFRawMap();
    // Clear the raw map
    virtual  void  Clear(const char *opt = "");
    // Set a single raw
    virtual  void  SetHit(Int_t *slot, Int_t idigit);
    virtual  void  SetHit(Int_t *slot);
    // Get index of hit in the list of digits
    virtual Int_t  GetHitIndex(Int_t *vol) const;
    // Get pointer to digit
    virtual TObject*  GetHit(Int_t *vol) const;
    // Test hit status
    virtual FlagType TestHit(Int_t *vol) const;
    
 private:
    // Check index
    Int_t CheckedIndex(Int_t *slot) const;
 private:
    Int_t fNtrm;            // Number of TRM
    Int_t fNtrmChain;       // Number of TRM chains per TRM
    Int_t fNtdc;            // Number of TDCs per TRM
    Int_t fNtdcChannel;     // Number of TDC channels per TDC

    TClonesArray *fRawData; // Pointer to raw data
    Int_t fMaxIndex;        // maximum index in hit map
    Int_t *fRawMap;         // ! [fMaxIndex]         

    ClassDef(AliTOFRawMap,0) // Implements RawMap as a 1-dim array
};
#endif	
