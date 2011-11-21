// @(#) $Id$
// Original: AliHLTInterMerger.h,v 1.4 2004/02/02 15:00:34 loizides Exp $

#ifndef ALIHLTTPCINTERMERGER_H
#define ALIHLTTPCINTERMERGER_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTTPCInterMerger.h
    @author Uli Frankenfeld, maintained by Matthias Richter
    @date   
    @brief  The HLT TPC track segment merger
*/

#ifndef __CINT__ 
#include "AliHLTTPCMerger.h"
#endif

#include "AliHLTTPCRootTypes.h"

/** 
 * @class AliHLTTPCInterMerger
 * The HLTTPC track segment merger
 *
 * @ingroup alihlt_tpc
 */
class AliHLTTPCInterMerger : public AliHLTTPCMerger {

 public:
  AliHLTTPCInterMerger();
  virtual ~AliHLTTPCInterMerger();

  void Init(Int_t *row,Int_t p){fRowMin=row[0];fRowMax=row[1];fPatch=p;}
  void SlowMerge();
  Int_t Merge();
  void MMerge();  //Loop over tracks from different subsectors
  
 private:
  AliHLTTPCInterMerger(const AliHLTTPCInterMerger&);
  AliHLTTPCInterMerger& operator=(const AliHLTTPCInterMerger&);

  Int_t fPatch;  // current patch
  Int_t fRowMin; // min row
  Int_t fRowMax; // max row

  ClassDef(AliHLTTPCInterMerger,0) //Intermerging class
};

#endif
