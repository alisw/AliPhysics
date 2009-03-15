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

#include "AliHLTPHOSClusterizer.h"
#include "AliHLTPHOSBase.h"
#include "AliHLTLogging.h"
#include "TMath.h"
#include "AliHLTPHOSRecPointContainerStruct.h"
#include "AliHLTPHOSRecPointDataStruct.h"
#include "AliHLTPHOSDigitDataStruct.h"
#include "AliHLTPHOSDigitContainerDataStruct.h"
#include "TClonesArray.h"
#include "AliPHOSDigit.h"
#ifndef HAVENOT__PHOSRECOPARAMEMC // set from configure if EMC functionality not available in AliPHOSRecoParam
#include "AliPHOSRecoParam.h"
#else
#include "AliPHOSRecoParamEmc.h"
#endif
#include <iostream>

using namespace std;

ClassImp(AliHLTPHOSClusterizer);

AliHLTPHOSClusterizer::AliHLTPHOSClusterizer():
  AliHLTPHOSBase(),
  fRecPointDataPtr(0),
  fDigitDataPtr(0),
  fEmcClusteringThreshold(0),
  fEmcMinEnergyThreshold(0),
  fEmcTimeGate(0),
  fDigitsInCluster(0),
  fDigitContainerPtr(0),
  fMaxDigitIndexDiff(2*NZROWSMOD)
{
  //See header file for documentation
  fEmcClusteringThreshold = 0.2;
  fEmcMinEnergyThreshold = 0.03;
  fEmcTimeGate = 1.e-6 ;
}//end


AliHLTPHOSClusterizer::~AliHLTPHOSClusterizer()  
{
  //See header file for documentation
}

void 
AliHLTPHOSClusterizer::SetRecPointDataPtr(AliHLTPHOSRecPointDataStruct* recPointDataPtr)
{
  fRecPointDataPtr = recPointDataPtr;
}

void
AliHLTPHOSClusterizer::SetRecoParameters(AliPHOSRecoParam* params)
{
  //see header file for documentation
#ifndef HAVE_NOT_PHOSRECOPARAMEMC // set from configure if EMC functionality not available in AliPHOSRecoParam
  // the new AliPHOSRecoParam functions, available from revision
  //  fEmcClusteringThreshold = params->GetEMCClusteringThreshold();
  // fEmcMinEnergyThreshold = params->GetEMCMinE();
  //  fLogWeight = params->GetEMCLogWeight();
  params++;
  params--;
#else
  fEmcClusteringThreshold = params->GetClusteringThreshold();
  fEmcMinEnergyThreshold = params->GetMinE();
  fLogWeight = params->GetLogWeight();
#endif
}  

Int_t 
AliHLTPHOSClusterizer::ClusterizeEvent(UInt_t availableSize, UInt_t& totSize)
{
  //see header file for documentation
  Int_t nRecPoints = 0;

  UInt_t maxRecPointSize = sizeof(AliHLTPHOSRecPointDataStruct) + (sizeof(AliHLTPHOSDigitDataStruct) << 7); //Reasonable estimate... 

  //Clusterization starts
  for(UInt_t i = 0; i < fDigitContainerPtr->fNDigits; i++)
    { 
      fDigitsInCluster = 0;
      
      if(fDigitContainerPtr->fDigitDataStruct[i].fEnergy < fEmcClusteringThreshold)
	{
	  continue;
	}
      if(availableSize < (totSize + maxRecPointSize)) 
	{
	  return -1; //Might get out of buffer, exiting
	}

      // First digit is placed at the fDigits member variable in the recpoint
      fDigitDataPtr = &(fRecPointDataPtr->fDigits);

      fRecPointDataPtr->fAmp = 0;
      fRecPointDataPtr->fModule = fDigitContainerPtr->fDigitDataStruct[i].fModule;
      // Assigning digit data to the digit pointer
      fRecPointDataPtr->fDigits = fDigitContainerPtr->fDigitDataStruct[i];

      // Incrementing the pointer to be ready for new entry
      fDigitDataPtr++;

      fRecPointDataPtr->fAmp += fDigitContainerPtr->fDigitDataStruct[i].fEnergy;
      fDigitContainerPtr->fDigitDataStruct[i].fEnergy = 0;
      fDigitsInCluster++;
      nRecPoints++;

      // Scanning for the neighbours
      ScanForNeighbourDigits(i, fRecPointDataPtr);

      totSize += sizeof(AliHLTPHOSRecPointDataStruct) + (fDigitsInCluster-1)*sizeof(AliHLTPHOSDigitDataStruct);   
      fRecPointDataPtr->fMultiplicity = fDigitsInCluster;     

      fRecPointDataPtr = reinterpret_cast<AliHLTPHOSRecPointDataStruct*>(fDigitDataPtr);
    }//end of clusterization

   return nRecPoints;
}

void
AliHLTPHOSClusterizer::ScanForNeighbourDigits(Int_t index, AliHLTPHOSRecPointDataStruct* recPoint)
{
  //see header file for documentation
  Int_t max = TMath::Min((Int_t)fDigitContainerPtr->fNDigits, (Int_t)fMaxDigitIndexDiff+index);
  Int_t min = TMath::Max(0, (Int_t)(index - (Int_t)fMaxDigitIndexDiff));

  max = fDigitContainerPtr->fNDigits;
  min = 0;
  for(Int_t j = min; j < max; j++)
    {
      if(fDigitContainerPtr->fDigitDataStruct[j].fEnergy > fEmcMinEnergyThreshold)
	{
	  if(j != index)
	    {
	      if(AreNeighbours(&(fDigitContainerPtr->fDigitDataStruct[index]),
			       &(fDigitContainerPtr->fDigitDataStruct[j])))
		{
		  // Assigning value to digit ptr
		  *fDigitDataPtr = fDigitContainerPtr->fDigitDataStruct[j];
		  // Incrementing digit pointer to be ready for new entry
		  fDigitDataPtr++;

		  recPoint->fAmp += fDigitContainerPtr->fDigitDataStruct[j].fEnergy;
		  fDigitContainerPtr->fDigitDataStruct[j].fEnergy = 0;	      
		  fDigitsInCluster++;
		  ScanForNeighbourDigits(j, recPoint);
		}
	    }
	}
    }
  return;
}

Int_t 
AliHLTPHOSClusterizer::AreNeighbours(AliHLTPHOSDigitDataStruct* digit1, 
					    AliHLTPHOSDigitDataStruct* digit2)
{
  //see header file for documentation
  if ( (digit1->fModule == digit2->fModule) /*&& (coord1[1]==coord2[1])*/ ) // inside the same PHOS module
    { 
      Int_t rowdiff = TMath::Abs( digit1->fZ - digit2->fZ );  
      Int_t coldiff = TMath::Abs( digit1->fX - digit2->fX ); 
      if (( coldiff <= 1   &&  rowdiff < 1 ) || ( coldiff < 1   &&  rowdiff <= 1 ))
	{
	  if(TMath::Abs(digit1->fTime - digit2->fTime ) < fEmcTimeGate)
	    {
	      return 1; 
	    }
	}
    }
  return 0;
}
