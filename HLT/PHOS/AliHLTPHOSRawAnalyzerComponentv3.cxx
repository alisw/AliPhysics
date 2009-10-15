
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Per Thomas Hille, Oystein Djuvsland                   *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliHLTPHOSRawAnalyzerComponentv3.h"
#include "AliHLTPHOSMapper.h"
#include "AliHLTPHOSDefinitions.h"
#include "AliHLTPHOSUtilities.h"

AliHLTPHOSRawAnalyzerComponentv3::AliHLTPHOSRawAnalyzerComponentv3():
  AliHLTCaloRawAnalyzerComponentv3()
{
  //comment
}


AliHLTPHOSRawAnalyzerComponentv3::~AliHLTPHOSRawAnalyzerComponentv3()
{
  //comment
}

int 
AliHLTPHOSRawAnalyzerComponentv3::Deinit()
{
  //comment
  return 0;
}

void
AliHLTPHOSRawAnalyzerComponentv3::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  //comment
  list.clear();
  list.push_back( AliHLTPHOSDefinitions::fgkDDLPackedRawDataType | kAliHLTDataOriginPHOS );
}

int
AliHLTPHOSRawAnalyzerComponentv3::InitMapping( const int spec)
{ 

  //See base class for documentation
  // fPrintInfo = kFALSE;

  if(fMapperPtr == 0)
    {
      fMapperPtr = new AliHLTPHOSMapper();
    }
 
  if(fMapperPtr->GetIsInitializedMapping() == false)
    {
      HLTError("%d:%d, ERROR, mapping not initialized ", __FILE__, __LINE__ );
      exit(-2);
    }

  return iResult;
}


