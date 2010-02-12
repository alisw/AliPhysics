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
 * @file   AliHLTCaloClusterizer.cxx
 * @author Oystein Djuvsland
 * @date 
 * @brief  Clusterizer for PHOS HLT 
 */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTCaloClusterizer.h"
#include "AliHLTLogging.h"
#include "TMath.h"
#include "AliHLTCaloRecPointDataStruct.h"
#include "AliHLTCaloDigitDataStruct.h"
#include "AliHLTCaloDigitContainerDataStruct.h"
#include "AliHLTCaloConstantsHandler.h"

ClassImp(AliHLTCaloClusterizer);

AliHLTCaloClusterizer::AliHLTCaloClusterizer(TString det):
  AliHLTCaloConstantsHandler(det),
  fRecPointDataPtr(0),
  fDigitIndexPtr(0),
  fEmcClusteringThreshold(0),
  fEmcMinEnergyThreshold(0),
  fEmcTimeGate(0),
  fDigitsInCluster(0),
  fDigitsPointerArray(0),
  fDigitContainerPtr(0),
  fMaxDigitIndexDiff(0),
  fNDigits(0)
{
  //See header file for documentation
  fEmcClusteringThreshold = 0.2;
  fEmcMinEnergyThreshold = 0.03;
  fEmcTimeGate = 1.e-6 ;
  
  fMaxDigitIndexDiff = 2*fCaloConstants->GetNZROWSMOD();

}//end


//BALLE how do you set the right detector?
// AliHLTCaloClusterizer::AliHLTCaloClusterizer(const AliHLTCaloClusterizer &) :
//   AliHLTCaloConstantsHandler("BALLE"),
//   fRecPointDataPtr(0),
//   fDigitDataPtr(0),
//   fEmcClusteringThreshold(0),
//   fEmcMinEnergyThreshold(0),
//   fEmcTimeGate(0),
//   fDigitsInCluster(0),
//   fDigitContainerPtr(0),
//   fMaxDigitIndexDiff(0)
// {
//   // dummy copy constructor
// }//end


AliHLTCaloClusterizer::~AliHLTCaloClusterizer()  
{
  //See header file for documentation
}

void 
AliHLTCaloClusterizer::SetRecPointDataPtr(AliHLTCaloRecPointDataStruct* recPointDataPtr)
{
  // See header file for documentation
  fRecPointDataPtr = recPointDataPtr;
}

Int_t 
AliHLTCaloClusterizer::ClusterizeEvent(Int_t nDigits, UInt_t availableSize, UInt_t& totSize)
{
  //see header file for documentation
  Int_t nRecPoints = 0;
  
  fAvailableSize = availableSize;
  
  fNDigits = nDigits;

  UInt_t maxRecPointSize = sizeof(AliHLTCaloRecPointDataStruct) + (sizeof(AliHLTCaloDigitDataStruct) << 7); //Reasonable estimate... 

  //Clusterization starts
  for(Int_t i = 0; i < nDigits; i++)
    { 
      fDigitsInCluster = 0;
      //      printf("ENERGY: %f\n", fDigitsPointerArray[i]->fEnergy);
      if(fDigitsPointerArray[i]->fEnergy < fEmcClusteringThreshold)
	{
	  continue;
	}
	
	if(fAvailableSize < (sizeof(AliHLTCaloRecPointDataStruct)))
	{
	  HLTError("Out of buffer, stopping clusterisation");
	  return -1; 
	}
	
      //            printf("cluster candidate!\n");
      // First digit is placed at the fDigits member variable in the recpoint
      fDigitIndexPtr = &(fRecPointDataPtr->fDigits);

      fRecPointDataPtr->fAmp = 0;
      fRecPointDataPtr->fModule = fDigitsPointerArray[i]->fModule;

      // Assigning the digit to this rec point
      fRecPointDataPtr->fDigits = i;

      // Incrementing the pointer to be ready for new entry
      fDigitIndexPtr++;

      fRecPointDataPtr->fAmp += fDigitsPointerArray[i]->fEnergy;
      fDigitsPointerArray[i]->fEnergy = 0;
      fDigitsInCluster++;
      nRecPoints++;

      // Scanning for the neighbours
      if(ScanForNeighbourDigits(i, fRecPointDataPtr) != 0)
      {
	 return -1;
      }

      totSize += sizeof(AliHLTCaloRecPointDataStruct) + (fDigitsInCluster-1)*sizeof(AliHLTCaloDigitDataStruct);   
      fRecPointDataPtr->fMultiplicity = fDigitsInCluster;     
      //      printf("Rec point energy: %f\n", fRecPointDataPtr->fAmp);
      fRecPointDataPtr = reinterpret_cast<AliHLTCaloRecPointDataStruct*>(fDigitIndexPtr);

    }//end of clusterization

   return nRecPoints;
}

Int_t
AliHLTCaloClusterizer::ScanForNeighbourDigits(Int_t index, AliHLTCaloRecPointDataStruct* recPoint)
{
  //see header file for documentation
  Int_t max = TMath::Min(fNDigits, (Int_t)fMaxDigitIndexDiff+index);
  Int_t min = TMath::Max(0, (Int_t)(index - (Int_t)fMaxDigitIndexDiff));

  max = fNDigits;
  min = 0;
  for(Int_t j = min; j < max; j++)
    {
      if(fDigitsPointerArray[j]->fEnergy > fEmcMinEnergyThreshold)
	{
	  if(j != index)
	    {
	      if(AreNeighbours(fDigitsPointerArray[index],
			       fDigitsPointerArray[j]))
		{
		  if(fAvailableSize < (sizeof(Int_t)))
		     {
			HLTError("Out of buffer, stopping clusterisation");
			return -1; 
		     }	
		  // Assigning index to digit
		  *fDigitIndexPtr = j;
		  // Incrementing digit pointer to be ready for new entry
		  fDigitIndexPtr++;

		  recPoint->fAmp += fDigitsPointerArray[j]->fEnergy;
		  fDigitsPointerArray[j]->fEnergy = 0;	      
		  fDigitsInCluster++;
		  ScanForNeighbourDigits(j, recPoint);
		}
	    }
	}
    }
  return 0;
}

Int_t 
AliHLTCaloClusterizer::AreNeighbours(AliHLTCaloDigitDataStruct* digit1, 
					    AliHLTCaloDigitDataStruct* digit2)
{
  //see header file for documentation
  if ( (digit1->fModule == digit2->fModule) /*&& (coord1[1]==coord2[1])*/ ) // inside the same PHOS module
    { 
//       Int_t rowdiff = TMath::Abs( digit1->fZ - digit2->fZ );  
//       Int_t coldiff = TMath::Abs( digit1->fX - digit2->fX ); 
//       if (( coldiff <= 1   &&  rowdiff == 0 ) || ( coldiff == 0 &&  rowdiff <= 1 ))
// 	{
// 	  cout << "Are neighbours: digit (E = "  << digit1->fEnergy << ") with x = " << digit1->fX << " and z = " << digit1->fZ << 
// 	    " is neighbour with digit (E = " << digit2->fEnergy << ") with x = " << digit2->fX << " and z = " << digit2->fZ << endl;

// 	  if(TMath::Abs(digit1->fTime - digit2->fTime ) < fEmcTimeGate)
// 	    {
// 	      return 1; 
// 	    }
// 	}

      Float_t rowdiff = TMath::Abs( digit1->fZ - digit2->fZ );  
      Float_t coldiff = TMath::Abs( digit1->fX - digit2->fX ); 
      if (( coldiff <= 2.4   &&  rowdiff < 0.4 ) || ( coldiff < 0.4 &&  rowdiff <= 2.4 ))
	{
	  //	  cout << "Are neighbours: digit (E = "  << digit1->fEnergy << ") with x = " << digit1->fX << " and z = " << digit1->fZ << 
	  //	    " is neighbour with digit (E = " << digit2->fEnergy << ") with x = " << digit2->fX << " and z = " << digit2->fZ << endl;

	  if(TMath::Abs(digit1->fTime - digit2->fTime ) < fEmcTimeGate)
	    {
	      return 1; 
	    }
	}
      else
	{
	  //  cout << "Not neighbours: digit (E = "  << digit1->fEnergy << ") with x = " << digit1->fX << " and z = " << digit1->fZ << 
	  //  " is not neighbour with digit (E = " << digit2->fEnergy << ") with x = " << digit2->fX << " and z = " << digit2->fZ << endl;
	}
    }
  return 0;
}
