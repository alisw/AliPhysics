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

const AliHLTComponentDataType AliHLTPHOSRawAnalyzerComponent::inputDataTypes[]={kAliHLTVoidDataType,{0,"",""}}; //'zero' terminated array
const AliHLTComponentDataType AliHLTPHOSRawAnalyzerComponent::outputDataType=kAliHLTVoidDataType;



//ClassImp(AliHLTPHOSRawAnalyzerComponent) 
AliHLTPHOSRawAnalyzerComponent::AliHLTPHOSRawAnalyzerComponent():AliHLTProcessor(), eventCount(0)
{

} 

AliHLTPHOSRawAnalyzerComponent::~AliHLTPHOSRawAnalyzerComponent()
{

}


AliHLTPHOSRawAnalyzerComponent::AliHLTPHOSRawAnalyzerComponent(const AliHLTPHOSRawAnalyzerComponent & ) : AliHLTProcessor(), eventCount(0)
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
 Logging(kHLTLogInfo, "HLT", "PHOS", ",AliHLTPHOSRawAnalyzerComponen DoDeinit");
  return 0;

}

//const char* 
//AliHLTPHOSRawAnalyzerComponent::GetComponentID()
//{
//  return 0;
//}

void
AliHLTPHOSRawAnalyzerComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  const AliHLTComponentDataType* pType=inputDataTypes;
  while (pType->fID!=0) {
    list.push_back(*pType);
    pType++;
  }
}

AliHLTComponentDataType 
AliHLTPHOSRawAnalyzerComponent::GetOutputDataType()
{
  return outputDataType;
}

void
AliHLTPHOSRawAnalyzerComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  constBase = 0;inputMultiplier = 0;
}

int 
AliHLTPHOSRawAnalyzerComponent::DoEvent(const AliHLTComponentEventData&, const AliHLTComponentBlockData*, AliHLTComponentTriggerData&, AliHLTUInt8_t*, AliHLTUInt32_t&, std::vector<AliHLTComponentBlockData, std::allocator<AliHLTComponentBlockData> >&)
{
  Logging(kHLTLogInfo, "HLT", "Sample", "PhosHLTRawAnalyzerComonent, DoEvent");
  printf("\nPHOSHLT DoEvent: processing event: %d\n", eventCount);
  //  printf("\nPHOSHLT DoEvent, not yet implemented\n");
  eventCount++;
  return 0;
}


int
AliHLTPHOSRawAnalyzerComponent::DoInit( int argc, const char** argv )
{
  printf("\nInside AliHLTPHOSRawAnalyzerComponent:DoInit\n");
  if (argc==0 && argv==NULL) {
    // this is currently just to get rid of the warning "unused parameter"
  }
  return 0;
}
