#ifndef _ALIHLTTPCCLUSTERFORMAT_H_
#define _ALIHLTTPCCLUSTERFORMAT_H_


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
// #include "AliHLTTPCMemHandler.h"
// #include "AliHLTTPCTransform.h"
#include "AliHLTTPCSpacePointData.h"

struct AliHLTTPCClusterData
    {
	AliHLTUInt32_t fSpacePointCnt;
	AliHLTTPCSpacePointData fSpacePoints[];
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

#endif // _ALIHLTTPCCLUSTERFORMAT_H_
