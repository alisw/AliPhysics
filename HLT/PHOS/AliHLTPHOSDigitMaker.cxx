
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
//#include "AliHLTPHOSDigitContainerStruct.h"
#include "AliHLTPHOSDigitDataStruct.h"
#include "AliHLTPHOSDigitContainerDataStruct.h"
       

//ClassImp(AliHLTPHOSDigitMaker);

using namespace PhosHLTConst;

AliHLTPHOSDigitMaker::AliHLTPHOSDigitMaker() :
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
 
}
  
AliHLTPHOSDigitMaker::~AliHLTPHOSDigitMaker() 
{

}

Int_t
AliHLTPHOSDigitMaker::MakeDigits(AliHLTPHOSRcuCellEnergyDataStruct* rcuData)
{

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
 // fDigitArrayPtr->Clear();
  fDigitCount = 0;
}

  

 
 
