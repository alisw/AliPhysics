#ifndef ALIHLTPHOSNODULECELLENERGYDATASTRUCT_H
#define ALIHLTPHOSNODULECELLENERGYDATASTRUCT_H

#include "AliHLTPHOSCommonDefs.h"

struct AliHLTPHOSModuleCellEnergyDataStruct
{
  AliHLTUInt8_t fModuleID;
  AliHLTUInt16_t fCnt;
  AliHLTUInt16_t fValidData[N_ROWS_MOD*N_COLUMNS_MOD*N_GAINS];
  unsigned long cellEnergies[N_ROWS_MOD][N_COLUMNS_MOD][N_GAINS];
};


#endif
