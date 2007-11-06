
#include "AliHLTPHOSRcuDigitMaker.h"
#include "AliHLTPHOSDigit.h"
#include "AliHLTPHOSConstants.h"
#include "AliHLTPHOSBaseline.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TFile.h"

#include "AliHLTPHOSValidCellDataStruct.h"
#include "AliHLTPHOSRcuCellEnergyDataStruct.h"
//#include "AliHLTPHOSDigitContainerStruct.h"
#include "AliHLTPHOSDigitDataStruct.h"
#include "AliHLTPHOSRcuDigitContainerDataStruct.h"
       

//ClassImp(AliHLTPHOSRcuDigitMaker);

using namespace PhosHLTConst;

AliHLTPHOSRcuDigitMaker::AliHLTPHOSRcuDigitMaker() :
  AliHLTPHOSBase(),
  fCellDataPtr(0),
  // fDigitContainerStructPtr(0),
  fDigitArrayPtr(0),
  fDigitPtr(0),
  // fDigitStructPtr(0),
  fDigitCount(0), 
  fDigitThreshold(0),
  fNrPresamples(10)
{
  //comment
}
  
AliHLTPHOSRcuDigitMaker::~AliHLTPHOSRcuDigitMaker() 

{
}

Int_t
AliHLTPHOSRcuDigitMaker::MakeDigits(AliHLTPHOSRcuCellEnergyDataStruct* rcuData)
{
  //comment
  Int_t i = 0;
  Int_t j = 0;
  Int_t x = -1;
  Int_t z = -1;
  Float_t amplitude = 0;
  for ( i = 0; i < rcuData->fCnt; i++ )
  {
    fCellDataPtr = & ( rcuData->fValidData[i] );
    x = fCellDataPtr->fX;
    z = fCellDataPtr->fZ;
    amplitude = fCellDataPtr->fEnergy;
    if ( amplitude > fDigitThreshold )
      {
        fDigitStructPtr = & ( fDigitContainerStructPtr->fDigitDataStruct[j + fDigitCount] );
        fDigitStructPtr->fX = fCellDataPtr->fX;
        fDigitStructPtr->fZ = fCellDataPtr->fZ;
        fDigitStructPtr->fAmplitude = ( amplitude );
        fDigitStructPtr->fTime = fCellDataPtr->fTime ;
        fDigitStructPtr->fGain = ( fCellDataPtr->fGain );
        fDigitStructPtr->SetRawData ( fCellDataPtr->fData );
        fDigitStructPtr->fCrazyness  =   fCellDataPtr->fCrazyness;
        fDigitStructPtr->fBaseline = -1;
        j++;
      }
  }
  fDigitCount += j;
  return fDigitCount; 
}
/*
Int_t
AliHLTPHOSRcuDigitMaker::SetDigitsTree(TTree *tree)
{
  TBranch * digBranch = tree->Branch("digits","TClonesArray",fDebugDigitArrayPtr); 
}
*/

void
AliHLTPHOSRcuDigitMaker::Reset()
{ 
 // fDigitArrayPtr->Clear();
  fDigitCount = 0;
}

  

 
 
