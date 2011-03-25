/**************************************************************************
 * This file is property of and copyright by the Experimental Nuclear     *
 * Physics Group, Dep. of Physics                                         *
 * University of Oslo, Norway, 2007                                       *
 *                                                                        *
 * Author: Per Thomas Hille <perthi@fys.uio.no> for the ALICE HLT Project.*
 * Contributors are mentioned in the code where appropriate.              *
 * Please report bugs to perthi@fys.uio.no                                *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliHLTEMCALRawAnalyzerFastFitComponent.h"
#include "AliCaloRawAnalyzerFastFit.h"

//Evaluation of Amplitude and Peak position using Least Mean Square (FastFit) fit
//-----------
//-----------
//-----------
//-----------


AliHLTEMCALRawAnalyzerFastFitComponent  gAliHLTEMCALRawAnalyzerFastFitComponent;

AliHLTEMCALRawAnalyzerFastFitComponent::AliHLTEMCALRawAnalyzerFastFitComponent() : AliHLTEMCALRawAnalyzerComponent(kFastFit)
{
  // fAnalyzerPtr
  //  fAnalyzerPtr = new AliCaloRawAnalyzerFastFit();
}


AliHLTEMCALRawAnalyzerFastFitComponent::~AliHLTEMCALRawAnalyzerFastFitComponent()
{
  delete fAnalyzerPtr;
}



const char* 
AliHLTEMCALRawAnalyzerFastFitComponent::GetComponentID()
{
  return "EmcalRawFastFit";
}



AliHLTComponent* 
AliHLTEMCALRawAnalyzerFastFitComponent::Spawn()
{
  return new AliHLTEMCALRawAnalyzerFastFitComponent();
}


int 
AliHLTEMCALRawAnalyzerFastFitComponent::Deinit()
{
  return 0;
}
