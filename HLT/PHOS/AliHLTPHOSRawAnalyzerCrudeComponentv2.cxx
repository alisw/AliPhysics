// $Id$


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


#include "AliHLTPHOSRawAnalyzerCrudeComponentv2.h"
#include "AliHLTPHOSRawAnalyzerCrude.h"

AliHLTPHOSRawAnalyzerCrudeComponentv2 gAliHLTPHOSRawAnalyzerCrudeComponentv2;

//___________________________________________________________________________
AliHLTPHOSRawAnalyzerCrudeComponentv2::AliHLTPHOSRawAnalyzerCrudeComponentv2() 
{
  fAnalyzerPtr = new AliHLTPHOSRawAnalyzerCrude();
} 

//___________________________________________________________________________
AliHLTPHOSRawAnalyzerCrudeComponentv2::~AliHLTPHOSRawAnalyzerCrudeComponentv2()
{
  if(fAnalyzerPtr)
    {
      delete fAnalyzerPtr;
      fAnalyzerPtr = 0;
    }
}

//___________________________________________________________________________
AliHLTPHOSRawAnalyzerCrudeComponentv2::AliHLTPHOSRawAnalyzerCrudeComponentv2(const AliHLTPHOSRawAnalyzerCrudeComponentv2 & ):AliHLTPHOSRawAnalyzerComponentv2()
{

}

int
AliHLTPHOSRawAnalyzerCrudeComponentv2::Deinit()
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
AliHLTPHOSRawAnalyzerCrudeComponentv2::GetComponentID()
{
  return "PhosRawCrudev2";
}

//___________________________________________________________________________
AliHLTComponent*
AliHLTPHOSRawAnalyzerCrudeComponentv2::Spawn()
{
  return new AliHLTPHOSRawAnalyzerCrudeComponentv2;
}
