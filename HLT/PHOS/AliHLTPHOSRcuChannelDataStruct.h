#ifndef ALIHLTPHOSRCUCHANNELDATASTRUCT_H
#define ALIHLTPHOSRCUCHANNELDATASTRUCT_H

//#include "AliHLTPHOSCommonDefs.h"
#include "AliHLTPHOSValidChannelDataStruct.h"
#include "AliHLTPHOSConstants.h" 

using namespace PhosHLTConst;

struct AliHLTPHOSRcuChannelDataStruct
{
  AliHLTUInt32_t  fNValidChannels;
  AliHLTUInt8_t   fModuleID;
  AliHLTUInt8_t   fRcuX;
  AliHLTUInt8_t   fRcuZ; 
  //  AliHLTPHOSValidChannelDataStruct fValidData[N_ROWS_RCU*N_COLUMNS_RCU*N_GAINS];
  //  AliHLTPHOSValidChannelDataStruct fValidData[512];
  AliHLTPHOSValidChannelDataStruct fValidData[N_ZROWS_RCU*N_XCOLUMNS_RCU*N_GAINS];
  //  AliHLTUInt16_t  fBuffer[(ALTRO_MAX_SAMPLES+2)*N_ZROWS_RCU*N_XCOLUMNS_RCU*N_GAINS];
  

};



#endif  
