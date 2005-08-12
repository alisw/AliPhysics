#ifndef _ALIHLTTPCTRACKLETFORMAT_HPP_
#define _ALIHLTTPCTRACKLETFORMAT_HPP_

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* AliHLTTPCTrackletData
 */

#include "AliHLTDataTypes.h"
#include "AliL3TrackSegmentData.h"

struct AliHLTTPCTrackletData
    {
	AliHLTUInt32_t fTrackletCnt;
	AliL3TrackSegmentData fTracklets[];
	//AliL3SpacePointData fSpacePoints[];
    };


#endif // _ALIHLTTPCTRACKLETFORMAT_HPP_
