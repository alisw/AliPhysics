// $Id$
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

#ifndef ALIHLTEXTERNALTRACKPARAM_H
#define ALIHLTEXTERNALTRACKPARAM_H

#include "AliHLTDataTypes.h"
#include "AliHLTStdIncludes.h"

/**
 * @struct AliHLTExternalTrackParam
 * This in a struct for tracks which are closer to the format used by offline.
 * The offline parametrisation we find in AliExternalTrackParam class.
 * The covariance matrix is in Double_t since this is the format in AliExternalTrackParam.
 * This saves time in translating from Float_t to Double_t. The other values has to be copied 
 * anyway, so these can be Float_t for saveing space.The charge is now removed and the Pt has a sign.
 *
 * The array of points is just appended to the structure. The member of array size 0
 * is not supported by all compilers. An open issue is that the code can not work if
 * the array size is 1, this was just added in order to make it compile.
 * @ingroup alihlt_component_datatypes
 */
struct AliHLTExternalTrackParam
    {
      Float_t fAlpha;  // azimuthal angle of reference frame
      Float_t fX;      // x: radial distance
      Float_t fY;      // local Y-coordinate of a track (cm)
      Float_t fZ;      // local Z-coordinate of a track (cm)
      Float_t fSinPsi; // local sine of the track momentum azimuthal angle
      Float_t fTgl;    // tangent of the track momentum dip angle
      Float_t fq1Pt;   // 1/pt (1/(GeV/c))
      Float_t fC[15];  // covariance matrix
      Float_t fLastX;  // x of last point
      Float_t fLastY;  // y of last point
      Float_t fLastZ;  // z of last point
      Int_t   fTrackID;// track id
      UInt_t  fFlags;  // flags
      UInt_t  fNPoints;// number of points
#if defined(__HP_aCC) || defined(__DECCXX) || defined(__SUNPRO_CC)
      UInt_t  fPointIDs[1]; // array of points
#else
      UInt_t  fPointIDs[0]; // array of points
#endif
     };

typedef struct AliHLTExternalTrackParam AliHLTExternalTrackParam;

struct AliHLTTracksData {
  AliHLTUInt32_t fCount; // number of tracklets
#if defined(__HP_aCC) || defined(__DECCXX) || defined(__SUNPRO_CC)
  AliHLTExternalTrackParam fTracklets[1]; // array of tracklets
#else
  AliHLTExternalTrackParam fTracklets[]; // array of tracklets
#endif
};

typedef struct AliHLTTracksData AliHLTTracksData;

#endif
