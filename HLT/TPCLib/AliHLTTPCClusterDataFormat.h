#ifndef _ALIHLTTPCCLUSTERFORMAT_H_
#define _ALIHLTTPCCLUSTERFORMAT_H_

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* AliHLTTPCClusterData
 */

struct AliHLTTPCClusterData
    {
	AliHLTUInt32_t fSpacePointCnt;
	AliL3SpacePointData fSpacePoints[];
    };

#endif // _ALIHLTTPCCLUSTERFORMAT_H_
