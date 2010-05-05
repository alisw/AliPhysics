//-*- Mode: C++ -*-
// $Id: AliHLTPHOSRcuCellAccumulatedEnergyDataStruct.h 31490 2009-03-15 16:27:11Z odjuvsla $

#ifndef ALIHLTPHOSRCUCELLACCUMULATEDENERGYDATASTRUCT_H
#define ALIHLTPHOSRCUCELLACCUMULATEDENERGYDATASTRUCT_H

/***************************************************************************
 * Copyright(c) 2007, ALICE Experiment at CERN, All rights reserved.       *
 *                                                                         *
 * Author: Per Thomas Hille <perthi@fys.uio.no> for the ALICE HLT Project. *
 * Contributors are mentioned in the code where appropriate.               *
 *                                                                         *
 * Permission to use, copy, modify and distribute this software and its    *
 * documentation strictly for non-commercial purposes is hereby granted    *
 * without fee, provided that the above copyright notice appears in all    *
 * copies and that both the copyright notice and this permission notice    *
 * appear in the supporting documentation. The authors make no claims      *
 * about the suitability of this software for any purpose. It is           *
 * provided "as is" without express or implied warranty.                   *
 **************************************************************************/

//#include "AliHLTPHOSCommonDefs.h"
//#include "AliHLTCaloConstant.h"
// using namespace PhosHLTConst;
//using namespace CaloHLTConst;

#include "AliHLTEMCALConstant.h"

using namespace EmcalHLTConst;

struct AliHLTCaloRcuCellAccumulatedEnergyDataStruct
{
  AliHLTUInt8_t fModuleID;
  AliHLTUInt8_t fRcuX;
  AliHLTUInt8_t fRcuZ; 
  float fAccumulatedEnergies[NXCOLUMNSMOD][NZROWSMOD][NGAINS];
  AliHLTUInt32_t fHits[NXCOLUMNSMOD][NZROWSMOD][NGAINS];
  AliHLTUInt32_t fDeadChannelMap[NXCOLUMNSMOD][NZROWSMOD][NGAINS]; 
  
  /*
  AliHLTUInt8_t fModuleID;
  AliHLTUInt8_t fRcuX;
  AliHLTUInt8_t fRcuZ; 
  float fAccumulatedEnergies[NXCOLUMNSMOD][NZROWSMOD][NGAINS];
  AliHLTUInt32_t fHits[NXCOLUMNSMOD][NZROWSMOD][NGAINS];
  AliHLTUInt32_t fDeadChannelMap[NXCOLUMNSMOD][NZROWSMOD][NGAINS];
  */  

};


#endif
