//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTPHOSRCUCELLENERGYDATASTRUCT_H
#define ALIHLTPHOSRCUCELLENERGYDATASTRUCT_H

/**************************************************************************
 * Copyright(c) 2007, ALICE Experiment at CERN, All rights reserved.      *
 *                                                                        *
 * Author: Per Thomas Hille <perthi@fys.uio.no> for the ALICE HLT Project.*
 * Contributors are mentioned in the code where appropriate.              *
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
#include "AliHLTPHOSValidCellDataStruct.h"
#include "Rtypes.h"

#include "AliHLTPHOSConstants.h"
using namespace PhosHLTConst;

struct AliHLTPHOSRcuCellEnergyDataStruct
{
  //  Int_t fShmAddress;
  Int_t fSize;         //size in charcters
  Int_t fModuleID;
  Int_t fRcuX;
  Int_t fRcuZ;
  Int_t fCnt;
  Int_t fAlgorithm;
  Int_t fInfo;
  Bool_t fHasRawData;
  AliHLTPHOSValidCellDataStruct fValidData[NXCOLUMNSRCU*NZROWSRCU*NGAINS];
};

#endif
