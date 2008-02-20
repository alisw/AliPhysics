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
 * @brief  Digit maker for PHOS HLT  
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
#include "TH2F.h"

#include "AliHLTPHOSValidCellDataStruct.h"
#include "AliHLTPHOSRcuCellEnergyDataStruct.h"
#include "AliHLTPHOSDigitDataStruct.h"
#include "AliHLTPHOSDigitContainerDataStruct.h"
#include "AliHLTPHOSSharedMemoryInterface.h" // added by PTH

ClassImp(AliHLTPHOSDigitMaker);

using namespace PhosHLTConst;

AliHLTPHOSDigitMaker::AliHLTPHOSDigitMaker() :
  AliHLTPHOSBase(),
  fCellDataPtr(0),
  fDigitContainerStructPtr(0),
  fDigitArrayPtr(0),
  fDigitStructPtr(0),
  fDigitCount(0)
{
  // See header file for documentation
  for(UInt_t x = 0; x < N_XCOLUMNS_MOD; x++)
    {
      for(UInt_t z = 0; z < N_ZROWS_MOD; z++)
	{
	  fHighGainFactors[x][z] = 0.005;
	  fLowGainFactors[x][z] = 0.08;
	}
    }
  
}
  
AliHLTPHOSDigitMaker::~AliHLTPHOSDigitMaker() 
{
  //See header file for documentation
}

Int_t
AliHLTPHOSDigitMaker::MakeDigits(AliHLTPHOSRcuCellEnergyDataStruct* rcuData)
{
  //See header file for documentation
  
  Int_t j = 0;
  Int_t xMod = -1;
  Int_t zMod = -1;
  Float_t amplitude = 0;
 
  for(UInt_t x = 0; x < N_XCOLUMNS_RCU; x++)
    {
      for(UInt_t z = 0; z < N_ZROWS_RCU; z++) 
	{
	  fCellDataPtr = &(rcuData->fValidData[x][z][HIGH_GAIN]);
	  xMod = x + rcuData->fRcuX * N_XCOLUMNS_RCU;
	  zMod = z + rcuData->fRcuZ * N_ZROWS_RCU;
	  amplitude = fCellDataPtr->fEnergy;
	  if(amplitude > fDigitThresholds[xMod][zMod][HIGH_GAIN] && amplitude < MAX_BIN_VALUE)
	    {
	      fDigitStructPtr = &(fDigitContainerStructPtr->fDigitDataStruct[j+fDigitCount]);
	      fDigitStructPtr->fX = xMod;
	      fDigitStructPtr->fZ = zMod;
	      fDigitStructPtr->fAmplitude = amplitude;
	      fDigitStructPtr->fEnergy = amplitude * fHighGainFactors[xMod][zMod];
	      //	      cout << "Amplitude in GeVs: " << fDigitStructPtr->fEnergy << endl;
	      //TODO: fix time
	      fDigitStructPtr->fTime = fCellDataPtr->fTime * 0.0000001;
	      fDigitStructPtr->fCrazyness = fCellDataPtr->fCrazyness;
	      fDigitStructPtr->fModule = rcuData->fModuleID;
	      j++;
	    }
	  else if(amplitude >= MAX_BIN_VALUE)
	    {
	      fCellDataPtr = & (rcuData->fValidData[x][z][LOW_GAIN]);
	      amplitude = fCellDataPtr->fEnergy;
	      if(amplitude > fDigitThresholds[xMod][zMod][LOW_GAIN])
		{	      
		  fDigitStructPtr = &(fDigitContainerStructPtr->fDigitDataStruct[j+fDigitCount]);
		  fDigitStructPtr->fX = xMod;
		  fDigitStructPtr->fZ = zMod;
		  fDigitStructPtr->fAmplitude = amplitude;
		  fDigitStructPtr->fEnergy = amplitude * fLowGainFactors[xMod][zMod];
		  //TODO: fix time
		  fDigitStructPtr->fTime = fCellDataPtr->fTime  * 0.0000001;;
		  fDigitStructPtr->fCrazyness = fCellDataPtr->fCrazyness;
		  fDigitStructPtr->fModule = rcuData->fModuleID;
		  j++;
		}
	    }
	}
    }
  
  fDigitCount += j;
  fDigitContainerStructPtr->fNDigits = fDigitCount;
  return fDigitCount; 
}


void
AliHLTPHOSDigitMaker::Reset()
{ 
  //See header file for documentation
  fDigitCount = 0;
}

void 
AliHLTPHOSDigitMaker::SetGlobalHighGainFactor(Float_t factor)
{
  //See header file for documentation
  for(UInt_t x = 0; x < N_XCOLUMNS_MOD; x++)
    {
      for(UInt_t z = 0; z < N_ZROWS_MOD; z++)
	{
	  fHighGainFactors[x][z] = factor;
	}
    }
}

void
AliHLTPHOSDigitMaker::SetGlobalLowGainFactor(Float_t factor)
{
  //See header file for documentation
  for(UInt_t x = 0; x < N_XCOLUMNS_MOD; x++)
    {
      for(UInt_t z = 0; z < N_ZROWS_MOD; z++)
	{
	  fLowGainFactors[x][z] = factor;
	}
    }
}

void 
AliHLTPHOSDigitMaker::SetDigitThresholds(const char* filepath, Int_t nSigmas)
{
  //See header file for documentation

  TFile *histFile = new TFile(filepath);
  
  TH2F *lgHist = (TH2F*)histFile->Get("RMSLGMapHist");
  TH2F *hgHist = (TH2F*)histFile->Get("RMSHGMapHist");

  for(UInt_t x = 0; x < N_XCOLUMNS_MOD; x++)
    {
      for(UInt_t z = 0; z < N_ZROWS_MOD; z++)
	{
	  fDigitThresholds[x][z][LOW_GAIN] = lgHist->GetBinContent(x, z) * nSigmas;
	  fDigitThresholds[x][z][HIGH_GAIN] = hgHist->GetBinContent(x, z) * nSigmas;
	}
    }
}
