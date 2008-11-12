// @(#) $Id$
// Original: AliHLTGlobalMerger.h,v 1.6 2004/02/02 15:00:34 loizides 

#ifndef ALIHLTTPC_GLOBALMERGER_H
#define ALIHLTTPC_GLOBALMERGER_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTTPCGlobalMerger.h
    @author Uli Frankenfeld, maintained by Matthias Richter
    @date   
    @brief  The HLT TPC slice merger
*/

#ifndef  __CINT__
#include "AliHLTTPCMerger.h"
#endif

#include "AliHLTTPCRootTypes.h"

/** 
 * @class AliHLTTPCGlobalMerger
 * The HLTTPC Slice merger
 *
 * @ingroup alihlt_tpc
 */
class AliHLTTPCGlobalMerger : public AliHLTTPCMerger{

 public:
  AliHLTTPCGlobalMerger();
  virtual ~AliHLTTPCGlobalMerger();
  
  void Setup(Int_t first,Int_t last);
  void InitSlice(Int_t slice);
  void SlowMerge(const Char_t *path="./");
  void Merge();  //Loop over tracks from different sectors

 private:
  Double_t CheckTracks(AliHLTTPCTrack *innertrack,AliHLTTPCTrack *outertrack,Int_t slice);
  
  Int_t fNSlices; // no of slices
  Int_t fFirst;   // first slice?
  Int_t fLast;    // last slice?

  ClassDef(AliHLTTPCGlobalMerger,1) //Slice merger
};

#endif
