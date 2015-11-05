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


// Evaluation of amplitude using  just the
// max  ADC walu - pedestal, and pek position
// using the time index of the max
// amplitude
// ------------
// ------------


#include "AliHLTEMCALRawAnalyzerStandardComponent.h"
//#include "AliHLTCaloRawAnalyzerCrude.h"
#include "AliCaloRawAnalyzerKStandard.h"


//AliHLTEMCALRawAnalyzerStandardComponent::AliHLTEMCALRawAnalyzerStandardComponent : AliHLTEMCALRawAnalyzerComponent()
AliHLTEMCALRawAnalyzerStandardComponent::AliHLTEMCALRawAnalyzerStandardComponent() : AliHLTEMCALRawAnalyzerComponent( kStandard )
{
  // constructor
  //  fAnalyzerPtr = new   AliHLTCaloRawAnalyzerCrude();
  // fAnalyzerPtr = new   AliCaloRawAnalyzerCrude();
}


AliHLTEMCALRawAnalyzerStandardComponent::~AliHLTEMCALRawAnalyzerStandardComponent()
{
  // destructor
  /*
  if (0 != fAnalyzerPtr)
    {
      delete fAnalyzerPtr;
      fAnalyzerPtr = 0;
    }
  */
}

int 
AliHLTEMCALRawAnalyzerStandardComponent::DoDeinit()
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
AliHLTEMCALRawAnalyzerStandardComponent::GetComponentID()
{
  // component id
  return "EmcalRawStandard";
}


AliHLTComponent* 
AliHLTEMCALRawAnalyzerStandardComponent::Spawn()
{
  // spawn component
  return new AliHLTEMCALRawAnalyzerStandardComponent();
}
 
