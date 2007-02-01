#ifndef ALIHLTPHOSRCUCELLENERGYDATASTRUCT_H
#define ALIHLTPHOSRCUCELLENERGYDATASTRUCT_H

#include "AliHLTPHOSCommonDefs.h"
#include "AliHLTPHOSValidCellDataStruct.h"

struct AliHLTPHOSRcuCellEnergyDataStruct
{
  AliHLTUInt16_t fRcuX;
  AliHLTUInt16_t fRcuZ;
  AliHLTUInt16_t fModuleID;
  AliHLTUInt16_t fCnt;
  AliHLTPHOSValiCellDataStruct fValidData;
  Double_t fCellEnergies[ N_ROWS_RCU][ N_COLUMNS_RCU][N_GAINS];
  
};


#endif
