#ifndef ALIHLTPHOSRCUCHANNELDATASTRUCT_H
#define ALIHLTPHOSRCUCHANNELDATASTRUCT_H
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Author:  Per Thomas Hille  <perthi@fys.uio.no>                 *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

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
