// $Id$
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

#ifndef _ALIHLTEXTERNALTRACKPARAM_H_
#define _ALIHLTEXTERNALTRACKPARAM_H_

#include "AliHLTStdIncludes.h"

/**
 * @struct AliHLTExternalTrackParam
 * This in a struct for tracks which are closer to the format used by offline.
 * The offline parametrisation we find in AliExternalTrackParam class.
 * The covariance matrix is in Double_t since this is the format in AliExternalTrackParam.
 * This saves time in translating from Float_t to Double_t. The other values has to be copied 
 * anyway, so these can be Float_t for saveing space.The charge is now removed and the Pt has a sign.
 *
 * @ingroup alihlt_component_datatypes
 */
struct AliHLTExternalTrackParam
    {
      Float_t fAlpha;
      Float_t fX;
      Float_t fY;
      Float_t fZ;
      Float_t fLastX;
      Float_t fLastY;
      Float_t fLastZ;
      Float_t fq1Pt;
      Float_t fSinPsi;
      Float_t fTgl;
      Float_t fC[15];
      UInt_t  fNPoints;
#if defined(__HP_aCC) || defined(__DECCXX) || defined(__SUNPRO_CC)
      UInt_t  fPointIDs[1];
#else
      UInt_t  fPointIDs[0];
#endif
     };

typedef struct AliHLTExternalTrackParam AliHLTExternalTrackParam;

struct AliHLTTracksData {
	AliHLTUInt32_t fCount;
#if defined(__HP_aCC) || defined(__DECCXX) || defined(__SUNPRO_CC)
	AliHLTExternalTrackParam fTracklets[1];
#else
	AliHLTExternalTrackParam fTracklets[];
#endif
};

typedef struct AliHLTTracksData AliHLTTracksData;

#endif
