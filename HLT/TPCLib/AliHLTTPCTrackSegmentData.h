// @(#) $Id$
// Original: AliHLTTrackSegmentData.h,v 1.7 2005/03/31 04:48:59 cvetan 
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
	Double_t fPterr;
	Double_t fPsierr;
        Double_t fTglerr;
        Int_t fCharge;
#ifdef INCLUDE_TPC_HOUGH
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
#endif // INCLUDE_TPC_HOUGH
	UInt_t  fNPoints;
#if defined(__HP_aCC) || defined(__DECCXX) || defined(__SUNPRO_CC)
	UInt_t  fPointIDs[1];
#else
	UInt_t  fPointIDs[0];
#endif
    };

typedef struct AliHLTTPCTrackSegmentData AliHLTTPCTrackSegmentData;

#endif /* _ALIHLTTPCTRACKSEGMENTDATA_H_ */
