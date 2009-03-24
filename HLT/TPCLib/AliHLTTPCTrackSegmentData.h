// @(#) $Id$
// Original: AliHLTTrackSegmentData.h,v 1.7 2005/03/31 04:48:59 cvetan 

#ifndef ALIHLTTPCTRACKSEGMENTDATA_H
#define ALIHLTTPCTRACKSEGMENTDATA_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

#include "AliHLTTPCRootTypes.h"

/**
 * @struct AliHLTTPCTrackSegmentData
 * Primitive data exchange structure for TPC tracks.
 *
 * @ingroup alihlt_tpc_datastructs
 */
struct AliHLTTPCTrackSegmentData
    {
      /** x coordinate of the first point assigned to track segment */
      Float_t fX;                                                          ///
      /** y coordinate of the first point assigned to track segment */
      Float_t fY;                                                          ///
      /** z coordinate of the first point assigned to track segment */
      Float_t fZ;                                                          ///
      /** x coordinate of the last point assigned to track segment */
      Float_t fLastX;                                                      ///
      /** y coordinate of the last point assigned to track segment */
      Float_t fLastY;                                                      ///
      /** z coordinate of the last point assigned to track segment */
      Float_t fLastZ;                                                      ///
      /** transvers momentum at first point */
      Double_t fPt;                                                        ///
      /** local sine of the track momentum azimuthal angle at first point */
      Double_t fPsi;                                                       ///
      /** tangent of the track momentum dip angle at first point */
      Double_t fTgl;                                                       ///
      /** error at first point */
      Double_t fY0err;                                                     ///
      /** error at first point */
      Double_t fZ0err;                                                     ///
      /** error at first point */
      Double_t fPterr;                                                     ///
      /** error at first point */
      Double_t fPsierr;                                                    ///
      /** error at first point */
      Double_t fTglerr;                                                    ///
      /** total charge */
      Int_t fCharge;                                                       ///
#ifdef INCLUDE_TPC_HOUGH
#ifdef ROWHOUGHPARAMS
      /* needed for PDC */
      UInt_t  fWeight;                // hough tracking parameters, deprecated
      Int_t  fTrackID;                // hough tracking parameters, deprecated
      Int_t  fRowRange1;              // hough tracking parameters, deprecated
      Int_t  fRowRange2;              // hough tracking parameters, deprecated
      Int_t  fSector;                 // hough tracking parameters, deprecated
      Float_t  fPID;                  // hough tracking parameters, deprecated
      Float_t  fBinX;                 // hough tracking parameters, deprecated
      Float_t  fBinY;                 // hough tracking parameters, deprecated
      Float_t  fBinXSize;             // hough tracking parameters, deprecated
      Float_t  fBinYSize;             // hough tracking parameters, deprecated
#endif
#endif // INCLUDE_TPC_HOUGH
      /** number of points attached in the following array */
      UInt_t  fNPoints;                                                    ///
#if defined(__HP_aCC) || defined(__DECCXX) || defined(__SUNPRO_CC)
      /** array of assigned points */
      UInt_t  fPointIDs[1];                                                ///
#else
      /** array of assigned points */
      UInt_t  fPointIDs[0];                                                ///
#endif
    };

typedef struct AliHLTTPCTrackSegmentData AliHLTTPCTrackSegmentData;

/**
 * @struct AliHLTTPCTrackSegmentDataV1
 * Former structure track segments, valid until July 2008
 * revision 27415
 *
 * @ingroup alihlt_tpc_datastructs
 */
struct AliHLTTPCTrackSegmentDataV1
    {
      /** x coordinate of the first point assigned to track segment */
      Float_t fX;                                                          ///
      /** y coordinate of the first point assigned to track segment */
      Float_t fY;                                                          ///
      /** z coordinate of the first point assigned to track segment */
      Float_t fZ;                                                          ///
      /** x coordinate of the last point assigned to track segment */
      Float_t fLastX;                                                      ///
      /** y coordinate of the last point assigned to track segment */
      Float_t fLastY;                                                      ///
      /** z coordinate of the last point assigned to track segment */
      Float_t fLastZ;                                                      ///
      /** transvers momentum at first point */
      Double_t fPt;                                                        ///
      /** local sine of the track momentum azimuthal angle at first point */
      Double_t fPsi;                                                       ///
      /** tangent of the track momentum dip angle at first point */
      Double_t fTgl;                                                       ///
      /** error at first point */
      Double_t fPterr;                                                     ///
      /** error at first point */
      Double_t fPsierr;                                                    ///
      /** error at first point */
      Double_t fTglerr;                                                    ///
      /** total charge */
      Int_t fCharge;                                                       ///
      /** number of points attached in the following array */
      UInt_t  fNPoints;                                                    ///
#if defined(__HP_aCC) || defined(__DECCXX) || defined(__SUNPRO_CC)
      /** array of assigned points */
      UInt_t  fPointIDs[1];                                                ///
#else
      /** array of assigned points */
      UInt_t  fPointIDs[0];                                                ///
#endif
    };

typedef struct AliHLTTPCTrackSegmentDataV1 AliHLTTPCTrackSegmentDataV1;

#endif /* _ALIHLTTPCTRACKSEGMENTDATA_H_ */
