
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
#include "AliHLTPHOSDefinitions.h"
#include "AliHLTPHOSMapper.h"


AliHLTPHOSRawAnalyzerComponentv3::AliHLTPHOSRawAnalyzerComponentv3() :
   AliHLTCaloRawAnalyzerComponentv3("PHOS")
{
  // See header file for class documentation
   InitMapping(0x1); //using 0x1 to avoid error message
}


AliHLTPHOSRawAnalyzerComponentv3::~AliHLTPHOSRawAnalyzerComponentv3()
{
  //comment
}


void
AliHLTPHOSRawAnalyzerComponentv3::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  //comment
  list.clear();
  list.push_back( AliHLTPHOSDefinitions::fgkDDLPackedRawDataType);
}


AliHLTComponentDataType 
AliHLTPHOSRawAnalyzerComponentv3::GetOutputDataType()
{
  //comment
  return AliHLTPHOSDefinitions::fgkChannelDataType;
}


void AliHLTPHOSRawAnalyzerComponentv3::InitMapping ( const int specification )
{
   // See header file for class documentation
   fMapperPtr = new AliHLTPHOSMapper;
   fMapperPtr->InitDDLSpecificationMapping();
   fMapperPtr->InitAltroMapping(specification);
}

