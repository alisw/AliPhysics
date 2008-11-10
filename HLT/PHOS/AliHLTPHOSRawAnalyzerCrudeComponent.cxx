// $Id$


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


#include "AliHLTPHOSRawAnalyzerCrudeComponent.h"
#include "AliHLTPHOSRawAnalyzerCrude.h"

AliHLTPHOSRawAnalyzerCrudeComponent gAliHLTPHOSRawAnalyzerCrudeComponent;

//___________________________________________________________________________
AliHLTPHOSRawAnalyzerCrudeComponent::AliHLTPHOSRawAnalyzerCrudeComponent()
{
  fAnalyzerPtr = new AliHLTPHOSRawAnalyzerCrude();
} 

//___________________________________________________________________________
AliHLTPHOSRawAnalyzerCrudeComponent::~AliHLTPHOSRawAnalyzerCrudeComponent()
{
  if(fAnalyzerPtr)
    {
      delete fAnalyzerPtr;
      fAnalyzerPtr = 0;
    }
}

//___________________________________________________________________________
AliHLTPHOSRawAnalyzerCrudeComponent::AliHLTPHOSRawAnalyzerCrudeComponent(const AliHLTPHOSRawAnalyzerCrudeComponent & ):AliHLTPHOSRawAnalyzerComponent()
{

}

int
AliHLTPHOSRawAnalyzerCrudeComponent::Deinit()
{
  
  if(fAnalyzerPtr)
    {
      delete fAnalyzerPtr;
      fAnalyzerPtr = 0;
    }
  Logging(kHLTLogInfo, "HLT", "PHOS", ",AliHLTPHOSRawAnalyzerCrudeComponent Deinit");
  return 0;
}

//___________________________________________________________________________
const char* 
AliHLTPHOSRawAnalyzerCrudeComponent::GetComponentID()
{
  return "PhosRawCrude";
}

//___________________________________________________________________________
AliHLTComponent*
AliHLTPHOSRawAnalyzerCrudeComponent::Spawn()
{
  return new AliHLTPHOSRawAnalyzerCrudeComponent;
}
