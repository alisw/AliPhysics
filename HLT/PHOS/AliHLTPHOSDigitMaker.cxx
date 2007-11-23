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
 * @file   AliHLTPHOSClusterizer.cxx
 * @author Oystein Djuvsland
 * @date 
 * @brief  Clusterizer for PHOS HLT  
 */
      

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTPHOSDigitMaker.h"
#include "AliHLTPHOSDigit.h"
#include "AliHLTPHOSConstants.h"
#include "AliHLTPHOSBaseline.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TFile.h"

#include "AliHLTPHOSValidCellDataStruct.h"
#include "AliHLTPHOSRcuCellEnergyDataStruct.h"
#include "AliHLTPHOSDigitDataStruct.h"
#include "AliHLTPHOSDigitContainerDataStruct.h"


ClassImp(AliHLTPHOSDigitMaker);

using namespace PhosHLTConst;

AliHLTPHOSDigitMaker::AliHLTPHOSDigitMaker() :
  AliHLTPHOSBase(),
  fCellDataPtr(0),
  // fDigitContainerStructPtr(0),
  fDigitArrayPtr(0),
  fDigitPtr(0),
  // fDigitStructPtr(0),
  fDigitCount(0), 
  fNrPresamples(10),
  fDigitThreshold(0)
{
  // See header file for documentation
}
  
AliHLTPHOSDigitMaker::~AliHLTPHOSDigitMaker() 
{
  //See header file for documentation
}

Int_t
AliHLTPHOSDigitMaker::MakeDigits(AliHLTPHOSRcuCellEnergyDataStruct* rcuData)
{
  //See header file for documentation
  Int_t i = 0;
  Int_t j = 0;
  Int_t x = -1;
  Int_t z = -1;
  Float_t amplitude = 0;
  for ( i = 0; i < rcuData->fCnt; i++ )
  {
    fCellDataPtr = & ( rcuData->fValidData[i] );
    x = fCellDataPtr->fX + rcuData->fRcuX * N_XCOLUMNS_RCU;
    z = fCellDataPtr->fZ + rcuData->fRcuZ * N_ZROWS_RCU;
    amplitude = fCellDataPtr->fEnergy;
    if ( amplitude > fDigitThreshold )
      {
        fDigitStructPtr = & ( fDigitContainerStructPtr->fDigitDataStruct[j + fDigitCount] );
        fDigitStructPtr->fX = ( fCellDataPtr->fX + rcuData->fRcuX * N_XCOLUMNS_RCU );
        fDigitStructPtr->fZ = ( fCellDataPtr->fZ + rcuData->fRcuZ * N_ZROWS_RCU );
        fDigitStructPtr->fAmplitude = ( amplitude );
        fDigitStructPtr->fTime = ( fCellDataPtr->fTime );
        fDigitStructPtr->fGain = ( fCellDataPtr->fGain );
        fDigitStructPtr->SetRawData ( fCellDataPtr->fData );
        fDigitStructPtr->fCrazyness = ( fCellDataPtr->fCrazyness );
        fDigitStructPtr->fBaseline = -1;
        j++;
      }
  }
  fDigitCount += j;
  return fDigitCount; 
}
/*
Int_t
AliHLTPHOSDigitMaker::SetDigitsTree(TTree *tree)
{
  TBranch * digBranch = tree->Branch("digits","TClonesArray",fDebugDigitArrayPtr); 
}
*/

void
AliHLTPHOSDigitMaker::Reset()
{ 
  //fDigitArrayPtr->Clear();
  fDigitCount = 0;
}

  

 
 
