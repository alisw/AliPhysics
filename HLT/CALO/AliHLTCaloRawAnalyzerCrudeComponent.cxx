// $Id: AliHLTCaloRawAnalyzerCrudeComponent.cxx 29824 2008-11-10 13:43:55Z richterm $


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


#include "AliHLTCaloRawAnalyzerCrudeComponent.h"
#include "AliHLTCaloRawAnalyzerCrude.h"

//AliHLTCaloRawAnalyzerCrudeComponent gAliHLTCaloRawAnalyzerCrudeComponent;

//___________________________________________________________________________
AliHLTCaloRawAnalyzerCrudeComponent::AliHLTCaloRawAnalyzerCrudeComponent()
{
  fAnalyzerPtr = new AliHLTCaloRawAnalyzerCrude();
} 

//___________________________________________________________________________
AliHLTCaloRawAnalyzerCrudeComponent::~AliHLTCaloRawAnalyzerCrudeComponent()
{
  if(fAnalyzerPtr)
    {
      delete fAnalyzerPtr;
      fAnalyzerPtr = 0;
    }
}

//___________________________________________________________________________
AliHLTCaloRawAnalyzerCrudeComponent::AliHLTCaloRawAnalyzerCrudeComponent(const AliHLTCaloRawAnalyzerCrudeComponent & ):AliHLTCaloRawAnalyzerComponentv3()
{

}

int
AliHLTCaloRawAnalyzerCrudeComponent::Deinit()
{
  
  if(fAnalyzerPtr)
    {
      delete fAnalyzerPtr;
      fAnalyzerPtr = 0;
    }
  Logging(kHLTLogInfo, "HLT", "Calo", ",AliHLTCaloRawAnalyzerCrudeComponent Deinit");
  return 0;
}

//___________________________________________________________________________
const char* 
AliHLTCaloRawAnalyzerCrudeComponent::GetComponentID()
{
  return "CaloRawCrude";
}

//___________________________________________________________________________

/*
AliHLTComponent*
AliHLTCaloRawAnalyzerCrudeComponent::Spawn()
{
  return new AliHLTCaloRawAnalyzerCrudeComponent;
}
*/
