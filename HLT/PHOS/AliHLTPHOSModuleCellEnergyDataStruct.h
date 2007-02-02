#ifndef ALIHLTPHOSNODULECELLENERGYDATASTRUCT_H
#define ALIHLTPHOSNODULECELLENERGYDATASTRUCT_H

#include "AliHLTPHOSCommonDefs.h"
#include "AliHLTPHOSValiCellDataStruct.h"

struct AliHLTPHOSModuleCellEnergyDataStruct
{
  AliHLTUInt8_t fModuleID;
  AliHLTUInt16_t fCnt;
  //  AliHLTUInt16_t fValidData[N_ROWS_MOD*N_COLUMNS_MOD*N_GAINS];
  
  AliHLTPHOSValidCellDataStruct fValidData[N_RCUS]; 

  unsigned long cellEnergies[N_ROWS_MOD][N_COLUMNS_MOD][N_GAINS];
};


#endif
