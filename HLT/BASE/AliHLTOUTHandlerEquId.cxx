// $Id$

//**************************************************************************
//* This file is property of and copyright by the                          * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** @file   AliHLTOUTHandlerEquId.cxx
    @author Matthias Richter
    @date   
    @brief  HLTOUT handler returning equipment id from data type and spec.
*/

#include "AliHLTOUTHandlerEquId.h"
#include "AliHLTOUT.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTOUTHandlerEquId)

AliHLTOUTHandlerEquId::AliHLTOUTHandlerEquId()
{ 
  // Base class for DDL raw data redirection handlers.
  //
  // HLTOUT handlers of this type are used for the replacement of detector
  // reconstruction input by the corresponding data from the HLT output.
  // From the data type and specification of an HLT output block the
  // corresponding equipment id of the original detector streams is determined.
}

AliHLTOUTHandlerEquId::~AliHLTOUTHandlerEquId()
{
  // destructor
}

int AliHLTOUTHandlerEquId::ProcessData(AliHLTOUT* pData)
{
  // process data
  if (!pData) return -EINVAL;
  AliHLTComponentDataType dt=kAliHLTVoidDataType;
  AliHLTUInt32_t spec=kAliHLTVoidDataSpec;
  int iResult=pData->GetDataBlockDescription(dt, spec);
  if (iResult>=0) {
    iResult=(int)spec;
  }
  return iResult;
}
