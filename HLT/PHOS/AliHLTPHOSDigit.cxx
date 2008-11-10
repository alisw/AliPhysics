// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Oystein Djuvsland                                     *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliHLTPHOSDigit.h"
#include "AliHLTPHOSAltroConfig.h"

/**
 * Digit class for PHOS HLT
 *
 * @file   AliHLTPHOSDigit.cxx 
 * @author Oystein Djuvsland
 * @date
 * @brief  Digit class for PHOS HLT
 */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

ClassImp(AliHLTPHOSDigit);

AliHLTPHOSDigit::AliHLTPHOSDigit() :
  TObject(),
  AliHLTPHOSBase(),
  fX(-1),
  fZ(-1),
  fAmplitude(-1),
  fTime(-1),
  fEnergy(-1),
  fGain(-1),
  fSamples(55),
  fPreSamples(15),
  fTotalSamples(70),
  fDebugVar(-1),
  fData(0),
  fCrazyness(0),
  fBaseline(0)
{
  //See header file for documentation
  //added by PT
  fSamples = fNSamples;
  fPreSamples = fNPresamples;
  fTotalSamples = fNTotalSamples;
  //   fData = new Int_t[fNTotalSamples];
  fData = new Int_t[fNTotalSamples];

}

AliHLTPHOSDigit::~AliHLTPHOSDigit()
{
  //See header file for documentation
}

void 
AliHLTPHOSDigit::SetRawData(Int_t *dataPtr)
{
  // See header file for documentation
  //modified by PT
  //  for(Int_t i = 0; i < 70; i++)
  //    {
  //     fData[i] = dataPtr[i];
  //    }
  for(Int_t i = 0; i < fNTotalSamples; i++)
    {
      fData[i] = dataPtr[i];
    }
}
 

void 
AliHLTPHOSDigit::ResetDigit()
  {
    // See header file for documentation
    fZ = -1;
    fX = -1;
    fAmplitude = -1;
    fTime = -1;
    fEnergy =-1;
    fGain = -1;
    fSamples = 55;
    fPreSamples =15;
    fTotalSamples =70;
    fDebugVar = -1;
  }
