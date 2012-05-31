#ifndef ALIMERGECOMBI_H
#define ALIMERGECOMBI_H
/* Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////
//
//  Class to make combination for merging
//  Returns combination of input event numbers                
//  Author: Jiri Chudoba (CERN), 2001
//
////////////////////////////////////////////////////////////////////////

// --- ROOT system ---

#include <TNamed.h>

// --- AliRoot header files ---

class AliMergeCombi: public TNamed {

public:
  AliMergeCombi();
  AliMergeCombi(Int_t dim, Int_t sperb);
  virtual ~AliMergeCombi();
  Bool_t Combination(Int_t evNumber[], Int_t delta[]);
  
private:  
  Int_t fDim;               //! dimension of arrays evNumber and delta
  Int_t fSperb;             //! signal per background ratio
  Int_t fCounter;           //! counter for calls
  
  ClassDef(AliMergeCombi,1)
};

#endif // ALIMERGECOMBI_H

