/**************************************************************************
 * Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved.      *
 *                                                                        *
 * Author: Per Thomas Hille for the ALICE HLT Project.                    *
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

#include "AliHLTPHOSRawAnalyzerComponent.h"
#include <iostream>

//ClassImp(AliHLTPHOSRawAnalyzerComponent) 

AliHLTPHOSRawAnalyzerComponent::AliHLTPHOSRawAnalyzerComponent()
{

} 

AliHLTPHOSRawAnalyzerComponent::~AliHLTPHOSRawAnalyzerComponent()
{

}


AliHLTPHOSRawAnalyzerComponent::AliHLTPHOSRawAnalyzerComponent(const AliHLTPHOSRawAnalyzerComponent & ) : AliHLTProcessor()
{

}

int 
AliHLTPHOSRawAnalyzerComponent::Deinit()
{
  return 0;
}

int 
AliHLTPHOSRawAnalyzerComponent::DoDeinit()
{
  return 0;
}

const char* 
AliHLTPHOSRawAnalyzerComponent::GetComponentID()
{
  return 0;
}

void
AliHLTPHOSRawAnalyzerComponent::GetInputDataTypes(std::vector<AliHLTComponentDataType, std::allocator<AliHLTComponentDataType> >&)
{

}

AliHLTComponentDataType 
AliHLTPHOSRawAnalyzerComponent::GetOutputDataType()
{
  AliHLTComponentDataType tmp;
  return tmp;
}

void
AliHLTPHOSRawAnalyzerComponent::GetOutputDataSize(long unsigned int&, double&)
{

}

void 
AliHLTPHOSRawAnalyzerComponent::GetOutputDataSize(long  int&, double&)
{

}

AliHLTComponent* 
AliHLTPHOSRawAnalyzerComponent::Spawn()
{
  return 0;
}

int 
AliHLTPHOSRawAnalyzerComponent::DoEvent(const AliHLTComponentEventData&, const AliHLTComponentBlockData*, AliHLTComponentTriggerData&, AliHLTUInt8_t*, AliHLTUInt32_t&, std::vector<AliHLTComponentBlockData, std::allocator<AliHLTComponentBlockData> >&)
{
  printf("\nPHOSHLT DoEvent, not yet implemented\n");
  return 0;
}
