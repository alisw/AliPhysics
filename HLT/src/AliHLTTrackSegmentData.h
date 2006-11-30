// @(#) $Id$

#ifndef _ALIL3TRACKSEGMENTDATA_H_
#define _ALIL3TRACKSEGMENTDATA_H_

#include "AliHLTRootTypes.h"

struct AliHLTTrackSegmentData
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
	Double_t fPterr;
	Double_t fPsierr;
        Double_t fTglerr;
        Int_t fCharge;
#ifdef ROWHOUGHPARAMS
      /* needed for PDC */
        UInt_t  fWeight;
        Int_t  fTrackID;
        Int_t  fRowRange1;
        Int_t  fRowRange2;
        Int_t  fSector;
        Float_t  fPID;
        Float_t  fBinX;
        Float_t  fBinY;
        Float_t  fBinXSize;
        Float_t  fBinYSize;
#endif
	UInt_t  fNPoints;
#if defined(__HP_aCC) || defined(__DECCXX) || defined(__SUNPRO_CC)
	UInt_t  fPointIDs[1];
#else
	UInt_t  fPointIDs[0];
#endif
    };

typedef struct AliHLTTrackSegmentData AliHLTTrackSegmentData;

#endif /* _ALIL3TRACKSEGMENTDATA_H_ */
