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
	UInt_t  fNPoints;
	UInt_t  fPointIDs[0];
    };

typedef struct AliL3TrackSegmentData AliL3TrackSegmentData;




#endif /* _ALIL3TRACKSEGMENTDATA_H_ */
