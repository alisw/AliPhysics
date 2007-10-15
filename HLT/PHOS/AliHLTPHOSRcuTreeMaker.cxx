 
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


#include "AliHLTPHOSRcuTreeMaker.h"
#include "AliHLTPHOSBase.h"
#include "AliHLTPHOSRcuDigitContainerDataStruct.h"
#include "AliHLTPHOSDigitDataStruct.h"
#include "AliHLTPHOSDigit.h"
#include "TClonesArray.h"
#include "TTree.h"

ClassImp(AliHLTPHOSRcuTreeMaker);

AliHLTPHOSRcuTreeMaker::AliHLTPHOSRcuTreeMaker() :
  AliHLTPHOSBase(),
  fDigitArrayPtr(0),
  fDigitTreePtr(0)
{

  fDigitArrayPtr = new TClonesArray("AliHLTPHOSRcuDigit", 300); //!!!!!!!!!!!!!!!!
  fDigitTreePtr = new TTree("digitTree", "Digits Tree");

  fDigitTreePtr->Branch("Digit", &fDigitArrayPtr);

}

AliHLTPHOSRcuTreeMaker::~AliHLTPHOSRcuTreeMaker()
{
}

Int_t
AliHLTPHOSRcuTreeMaker::MakeDigitArray(AliHLTPHOSRcuDigitContainerDataStruct *digitContainer, Int_t nDigits)
{
  AliHLTPHOSDigit *digit = 0;
  AliHLTPHOSDigitDataStruct *digitStruct = 0;

  for(Int_t i = 0; i < digitContainer->fNDigits; i++)
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
}

void
AliHLTPHOSRcuTreeMaker::FillDigitTree()
{
  fDigitTreePtr->Fill();
  fDigitArrayPtr->Clear();
}
 
void 
AliHLTPHOSRcuTreeMaker::SetDigitTree(TTree *tree) 
{ 
  fDigitTreePtr = tree; 
  fDigitTreePtr->Branch("Digit", &fDigitArrayPtr);
}
