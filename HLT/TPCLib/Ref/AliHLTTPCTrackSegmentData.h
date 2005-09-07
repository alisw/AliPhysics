// @(#) $Id$

#ifndef _ALIHLTTPCTRACKSEGMENTDATA_H_
#define _ALIHLTTPCTRACKSEGMENTDATA_H_

#include "AliHLTTPCRootTypes.h"

struct AliHLTTPCTrackSegmentData
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
	UInt_t  fNPoints;
	UInt_t  fPointIDs[0];
    };

typedef struct AliHLTTPCTrackSegmentData AliHLTTPCTrackSegmentData;

#endif /* _ALIHLTTPCTRACKSEGMENTDATA_H_ */
