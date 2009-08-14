
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


#include "AliHLTPHOSRawAnalyzerCrudeComponentv3.h"
#include "AliHLTPHOSRawAnalyzerCrude.h"

AliHLTPHOSRawAnalyzerCrudeComponentv3 gAliHLTPHOSRawAnalyzerCrudeComponentv3;

//___________________________________________________________________________
AliHLTPHOSRawAnalyzerCrudeComponentv3::AliHLTPHOSRawAnalyzerCrudeComponentv3() 
{
  fAnalyzerPtr = new AliHLTPHOSRawAnalyzerCrude();
} 

//___________________________________________________________________________
AliHLTPHOSRawAnalyzerCrudeComponentv3::~AliHLTPHOSRawAnalyzerCrudeComponentv3()
{
  if(fAnalyzerPtr)
    {
      delete fAnalyzerPtr;
      fAnalyzerPtr = 0;
    }
}

//___________________________________________________________________________
AliHLTPHOSRawAnalyzerCrudeComponentv3::AliHLTPHOSRawAnalyzerCrudeComponentv3(const AliHLTPHOSRawAnalyzerCrudeComponentv3 & ):AliHLTPHOSRawAnalyzerComponentv3()
{

}

int
AliHLTPHOSRawAnalyzerCrudeComponentv3::Deinit()
{
  
  if(fAnalyzerPtr)
    {
      delete fAnalyzerPtr;
      fAnalyzerPtr = 0;
    }
  return 0;
}

//___________________________________________________________________________
const char* 
AliHLTPHOSRawAnalyzerCrudeComponentv3::GetComponentID()
{
  return "PhosRawCrudev3";
}

//___________________________________________________________________________
AliHLTComponent*
AliHLTPHOSRawAnalyzerCrudeComponentv3::Spawn()
{
  return new AliHLTPHOSRawAnalyzerCrudeComponentv3;
}
