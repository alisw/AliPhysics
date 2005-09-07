#ifndef _ALIHLTTPCSECTORTRACKFORMAT_HPP_
#define _ALIHLTTPCSECTORTRACKFORMAT_HPP_

/*
***************************************************************************
**
** $Author$ - Initial Version by Timm Morten Steinbeck
**
** $Id$ 
**
***************************************************************************
*/

#include "AliHLTDataTypes.h"
#include "AliHLTTPCTrackSegmentData.h"

struct AliHLTTPCTrackletData
    {
	AliHLTUInt32_t fTrackletCnt;
	AliHLTTPCTrackSegmentData fTracklets[];
	//AliHLTTPCSpacePointData fSpacePoints[];
    };


/*
***************************************************************************
**
** $Author$ - Initial Version by Timm Morten Steinbeck
**
** $Id$ 
**
***************************************************************************
*/

#endif // _ALIHLTTPCSECTORTRACKFORMAT_HPP_
