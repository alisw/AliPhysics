#ifndef ALICALOFITSUBARRAY_H
#define ALICALOFITSUBARRAY_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: $ */

#include "Rtypes.h"

// Container class to hold info from bunches/samples
// selected for signal fitting.
// Variables are:
//  Int_t   fBunchIndex;  // Index for selected bunch
//  Int_t   fMaxRev;      // Max index in reversed array
//  Int_t   fFirst;   // first index in array used for fit
//  Int_t   fLast;    // last index in array used for fit

class  AliCaloFitSubarray
{
 public:
  explicit AliCaloFitSubarray( const Int_t bunchIndex, 
			       const Int_t maxrev, 
			       const Int_t first, 
			       const Int_t last ); 

  explicit AliCaloFitSubarray(const Int_t init);

  AliCaloFitSubarray(const AliCaloFitSubarray &fitSubarray);
  AliCaloFitSubarray& operator = (const AliCaloFitSubarray& source) ;

  virtual  ~AliCaloFitSubarray();

  void SetBunchIndex(Int_t i) { fBunchIndex = i;};
  void SetMaxRev(Int_t i) { fMaxRev = i;};
  void SetFirst(Int_t i) { fFirst = i; };
  void SetLast(Int_t i) { fLast = i; };
  
  Int_t GetBunchIndex() const  { return fBunchIndex;};
  Int_t GetMaxRev() const { return fMaxRev;};
  Int_t GetFirst() const { return fFirst; };
  Int_t GetLast() const { return fLast; };
  
 private:

  Int_t   fBunchIndex;  // Index for selected bunch
  Int_t   fMaxRev;      // Max index in reversed array
  Int_t   fFirst;   // first index in array used for fit
  Int_t   fLast;    // last index in array used for fit
};

#endif
