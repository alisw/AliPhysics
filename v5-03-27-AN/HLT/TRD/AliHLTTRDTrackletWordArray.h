// $Id$
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

#ifndef _ALIHLTTRDTRACKLETWORDARRAY_H_
#define _ALIHLTTRDTRACKLETWORDARRAY_H_

#include "AliHLTDataTypes.h"
#include "AliHLTStdIncludes.h"

struct AliHLTTRDTrackletWordArray {
  AliHLTTRDTrackletWordArray(Int_t det):fDet(det),fCount(0){}
  Int_t          fDet;
  AliHLTUInt32_t fCount;
#if defined(__HP_aCC) || defined(__DECCXX) || defined(__SUNPRO_CC)
  UInt_t fTracklets[1];
#else
  UInt_t fTracklets[];
#endif
};

//typedef struct AliTRDHLTTrackletWordArray AliTRDHLTTrackletWordArray ?plain c!

#endif
