// $Id$
#ifndef ALIHLTITSCLUSTERFORMAT_H
#define ALIHLTITSCLUSTERFORMAT_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *

#include "AliHLTITSSpacePointData.h"

/**
 * @struct AliHLTITSClusterData
 * Primitive data exchange structure for ITS clusters.
 * The data format contains one 32bit count member and the array
 * of spacepoint data structures.
 *
 * @ingroup alihlt_its_datastructs
 */
struct AliHLTITSClusterData
    {
	AliHLTUInt32_t fSpacePointCnt;
#ifndef __SUNPRO_CC
	AliHLTITSSpacePointData fSpacePoints[];
#else
	AliHLTITSSpacePointData fSpacePoints[1];
#endif
    };

#endif // ALIHLTITSCLUSTERFORMAT_H
