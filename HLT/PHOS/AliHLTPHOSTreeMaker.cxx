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

/**
 * @file   AliHLTPHOSTreeMaker.cxx
 * @author Oystein Djuvsland
 * @date
 * @brief  Tree maker  for PHOS HLT
 */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTPHOSTreeMaker.h"
#include "AliHLTPHOSBase.h"
#include "AliHLTPHOSDigitContainerDataStruct.h"
#include "AliHLTPHOSDigitDataStruct.h"
#include "AliHLTPHOSDigit.h"
#include "TClonesArray.h"
#include "TTree.h"

ClassImp(AliHLTPHOSTreeMaker);



AliHLTPHOSTreeMaker::AliHLTPHOSTreeMaker() :
  AliHLTPHOSBase(),
  fDigitArrayPtr(0),
  fDigitTreePtr(0)
{
  //See header file for documentation
  fDigitArrayPtr = new TClonesArray("AliHLTPHOSDigit", 300); //!!!!!!!!!!!!!!!!
  fDigitTreePtr = new TTree("digitTree", "Digits Tree");

  fDigitTreePtr->Branch("Digit", &fDigitArrayPtr);

}


AliHLTPHOSTreeMaker::~AliHLTPHOSTreeMaker()
{
  //See header file for documentation
}

Int_t
AliHLTPHOSTreeMaker::MakeDigitArray(AliHLTPHOSDigitContainerDataStruct *digitContainer, Int_t nDigits)
{
  /*

  //See header file for documentation
  AliHLTPHOSDigit *digit = 0;
  AliHLTPHOSDigitDataStruct *digitStruct = 0;

  for(UInt_t i = 0; i < digitContainer->fNDigits; i++)
    {
      digitStruct = &(digitContainer->fDigitDataStruct[i]);
      digit = (AliHLTPHOSDigit*)fDigitArrayPtr->New(i + nDigits);
      digit->SetX(digitStruct->fX);
      digit->SetZ(digitStruct->fZ);
      digit->SetAmplitude(digitStruct->fAmplitude);
      digit->SetTime(digitStruct->fTime);
      digit->SetGain(digitStruct->fGain);
      digit->SetRawData(digitStruct->fData);
      digit->SetCrazyness(digitStruct->fCrazyness);
      digit->SetBaseline(digitStruct->fBaseline);
    }
  return digitContainer->fNDigits;
  */
  return 0;
}

void
AliHLTPHOSTreeMaker::FillDigitTree()
{
  //See header file for documentation
  fDigitTreePtr->Fill();
  fDigitArrayPtr->Clear();
}
 
void 
AliHLTPHOSTreeMaker::SetDigitTree(TTree *tree) 
{ 
  //See header file for documentation
  fDigitTreePtr = tree; 
  fDigitTreePtr->Branch("Digit", &fDigitArrayPtr);
}
