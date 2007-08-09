#ifndef _ALIHLTTPCTRACKLETFORMAT_HPP_
#define _ALIHLTTPCTRACKLETFORMAT_HPP_

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* AliHLTTPCTrackletData
 */

#include "AliHLTDataTypes.h"
#include "AliHLTTPCTrackSegmentData.h"

struct AliHLTTPCTrackletData
    {
	AliHLTUInt32_t fTrackletCnt;
#ifndef __SUNPRO_CC
	AliHLTTPCTrackSegmentData fTracklets[];
#else
	AliHLTTPCTrackSegmentData fTracklets[1];
#endif
	//AliHLTTPCSpacePointData fSpacePoints[];
    };


#endif // _ALIHLTTPCTRACKLETFORMAT_HPP_
