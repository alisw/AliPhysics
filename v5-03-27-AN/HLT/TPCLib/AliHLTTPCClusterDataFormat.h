// $Id$
#ifndef _ALIHLTTPCCLUSTERFORMAT_H_
#define _ALIHLTTPCCLUSTERFORMAT_H_

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *

#include "AliHLTTPCSpacePointData.h"

/**
 * @struct AliHLTTPCClusterData
 * Primitive data exchange structure for TPC clusters.
 * The data format contains one 32bit count member and the array
 * of spacepoint data structures.
 *
 * @ingroup alihlt_tpc_datastructs
 */
struct AliHLTTPCClusterData
    {
	AliHLTUInt32_t fSpacePointCnt;
#if !defined(__SUNPRO_CC) && !defined(__clang__)
	AliHLTTPCSpacePointData fSpacePoints[];
#else
	AliHLTTPCSpacePointData fSpacePoints[1];
#endif
    };

#endif // _ALIHLTTPCCLUSTERFORMAT_H_
