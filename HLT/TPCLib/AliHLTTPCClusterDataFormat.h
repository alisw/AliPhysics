#ifndef _ALIHLTTPCCLUSTERFORMAT_H_
#define _ALIHLTTPCCLUSTERFORMAT_H_

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliHLTTPCSpacePointData.h"

/* AliHLTTPCClusterData
 */
struct AliHLTTPCClusterData
    {
	AliHLTUInt32_t fSpacePointCnt;
#ifndef __SUNPRO_CC
	AliHLTTPCSpacePointData fSpacePoints[];
#else
	AliHLTTPCSpacePointData fSpacePoints[1];
#endif
    };

#endif // _ALIHLTTPCCLUSTERFORMAT_H_
