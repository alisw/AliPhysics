// @(#) $Id$

#ifndef _ALIL3TRACKSEGMENTDATA_H_
#define _ALIL3TRACKSEGMENTDATA_H_

#include "AliL3RootTypes.h"

struct AliL3TrackSegmentData
    {
	Float_t fX;
	Float_t fY;
	Float_t fZ;
        Float_t fLastX;
        Float_t fLastY;
        Float_t fLastZ;
	Double_t fPt;
	Double_t fPsi;
        Double_t fTgl;
        Int_t fCharge;
#ifdef ROWHOUGH
        UInt_t  fWeight;
        Int_t  fTrackID;
#endif
	UInt_t  fNPoints;
#if defined(__HP_aCC) || defined(__DECCXX) || defined(__SUNPRO_CC)
	UInt_t  fPointIDs[1];
#else
	UInt_t  fPointIDs[0];
#endif
    };

typedef struct AliL3TrackSegmentData AliL3TrackSegmentData;

#endif /* _ALIL3TRACKSEGMENTDATA_H_ */
