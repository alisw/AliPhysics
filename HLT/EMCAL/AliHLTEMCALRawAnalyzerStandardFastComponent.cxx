/**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        *
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Rudiger Haake <ruediger.haake@cern.ch>, Yale          *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************/

// Component structure copied from AliHLTEMCALRawAnalyzerStandardComponent

#include "AliHLTEMCALRawAnalyzerStandardFastComponent.h"
#include "AliCaloRawAnalyzerKStandardFast.h"


//AliHLTEMCALRawAnalyzerStandardFastComponent::AliHLTEMCALRawAnalyzerStandardFastComponent : AliHLTEMCALRawAnalyzerComponent()
AliHLTEMCALRawAnalyzerStandardFastComponent::AliHLTEMCALRawAnalyzerStandardFastComponent() : AliHLTEMCALRawAnalyzerComponent( kStandardFast )
{
}


AliHLTEMCALRawAnalyzerStandardFastComponent::~AliHLTEMCALRawAnalyzerStandardFastComponent()
{
  // destructor
}

int 
AliHLTEMCALRawAnalyzerStandardFastComponent::DoDeinit()
{
  //comment
  if (0 != fAnalyzerPtr)
    {
      delete fAnalyzerPtr;
      fAnalyzerPtr = 0;
    }

  return AliHLTEMCALRawAnalyzerComponent::DoDeinit();
}

const char* 
AliHLTEMCALRawAnalyzerStandardFastComponent::GetComponentID()
{
  // component id
  return "EmcalRawStandardFast";
}


AliHLTComponent* 
AliHLTEMCALRawAnalyzerStandardFastComponent::Spawn()
{
  // spawn component
  return new AliHLTEMCALRawAnalyzerStandardFastComponent();
}
 
