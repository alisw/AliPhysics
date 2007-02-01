#ifndef ALIHLTPHOSVALIDCELLDATASTRUCT_H
#define ALIHLTPHOSVALIDCELLDATASTRUCT_H 

#include "AliHLTPHOSCommonDefs.h"
//#include "AliHLTDataTypes.h"
//#include "Rtypes.h"

struct AliHLTPHOSValiCellDataStruct
{
  AliHLTUInt16_t fRcuX;
  AliHLTUInt16_t fRcuZ;
  AliHLTUInt16_t fModuleID;
  AliHLTUInt16_t fRow;
  AliHLTUInt16_t fCol;
  AliHLTUInt16_t fGain;
};


#endif
