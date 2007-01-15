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

#include "AliHLTPHOSAnalyzerComponent.h"
#include <iostream>

//ClassImp(AliHLTPHOSAnalyzerComponent) 

AliHLTPHOSAnalyzerComponent::AliHLTPHOSAnalyzerComponent()
{

} 

AliHLTPHOSAnalyzerComponent::~AliHLTPHOSAnalyzerComponent()
{

}


AliHLTPHOSAnalyzerComponent::AliHLTPHOSAnalyzerComponent(const AliHLTPHOSAnalyzerComponent & ) : AliHLTProcessor()
{

}

int 
AliHLTPHOSAnalyzerComponent::Deinit()
{
  return 0;
}

int 
AliHLTPHOSAnalyzerComponent::DoDeinit()
{
  return 0;
}

const char* 
AliHLTPHOSAnalyzerComponent::GetComponentID()
{
  return 0;
}

void
AliHLTPHOSAnalyzerComponent::GetInputDataTypes(std::vector<AliHLTComponentDataType, std::allocator<AliHLTComponentDataType> >&)
{

}

AliHLTComponentDataType 
AliHLTPHOSAnalyzerComponent::GetOutputDataType()
{
  AliHLTComponentDataType tmp;
  return tmp;
}

void
AliHLTPHOSAnalyzerComponent::GetOutputDataSize(long unsigned int&, double&)
{

}

void 
AliHLTPHOSAnalyzerComponent::GetOutputDataSize(long  int&, double&)
{

}

AliHLTComponent* 
AliHLTPHOSAnalyzerComponent::Spawn()
{
  return 0;
}

int 
AliHLTPHOSAnalyzerComponent::DoEvent(const AliHLTComponentEventData&, const AliHLTComponentBlockData*, AliHLTComponentTriggerData&, AliHLTUInt8_t*, AliHLTUInt32_t&, std::vector<AliHLTComponentBlockData, std::allocator<AliHLTComponentBlockData> >&)
{
  printf("\nPHOSHLT DoEvent, not yet implemented\n");
  return 0;
}
