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
  fRecPointArray(0),
  fRecPointDataPtr(0),
  fFirstRecPointPtr(0),
  fArraySize(0),
  fAvailableSize(0),
  fUsedSize(0),
  fNRecPoints(0),
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
  
  
  fArraySize = 10;
  fRecPointArray = new AliHLTCaloRecPointDataStruct*[fArraySize];
  
  fAvailableSize = sizeof(AliHLTCaloRecPointDataStruct) * 20;
  fRecPointDataPtr = reinterpret_cast<AliHLTCaloRecPointDataStruct*>(new UChar_t[fAvailableSize]);
  fFirstRecPointPtr = fRecPointDataPtr;  
  printf("Start of rec point data: %x, end of rec point data: %x\n", fRecPointDataPtr, reinterpret_cast<UChar_t*>(fRecPointDataPtr) + fAvailableSize);

}//end

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
AliHLTCaloClusterizer::ClusterizeEvent(Int_t nDigits)
{
  //see header file for documentation
  Int_t nRecPoints = 0;
  fNRecPoints = 0;
  fUsedSize = 0;
  fNDigits = nDigits;
  fRecPointDataPtr = fFirstRecPointPtr;
  //Clusterization starts
  for(Int_t i = 0; i < nDigits; i++)
    { 
      fDigitsInCluster = 0;
      //      printf("ENERGY: %f\n", fDigitsPointerArray[i]->fEnergy);
      if(fDigitsPointerArray[i]->fEnergy < fEmcClusteringThreshold)
	{
	  continue;
	}
      CheckArray();
      CheckBuffer();
      //            printf("cluster candidate!\n");
      // First digit is placed at the fDigits member variable in the recpoint
      fDigitIndexPtr = &(fRecPointDataPtr->fDigits);

      fRecPointDataPtr->fAmp = 0;
      fRecPointDataPtr->fModule = fDigitsPointerArray[i]->fModule;

      // Assigning the digit to this rec point
      fRecPointDataPtr->fDigits = i;
      printf("Clusterizier: adding digit:  index pointer: %x, index: %d\n", fDigitIndexPtr, *fDigitIndexPtr);
      fUsedSize += sizeof(AliHLTCaloRecPointDataStruct);
      
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

      //fUsedSize += sizeof(AliHLTCaloRecPointDataStruct) + (fDigitsInCluster-1)*sizeof(AliHLTCaloDigitDataStruct);   
      
      fRecPointDataPtr->fMultiplicity = fDigitsInCluster;     
      printf("Rec point energy: %f\n", fRecPointDataPtr->fAmp);
      printf("Multiplicity: %d\n", fDigitsInCluster);
      fRecPointArray[fNRecPoints] = fRecPointDataPtr; 
      
      fRecPointDataPtr = reinterpret_cast<AliHLTCaloRecPointDataStruct*>(fDigitIndexPtr);
      
    }//end of clusterization
   fNRecPoints = nRecPoints;
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
// 		  if((fAvailableSize - fUsedSize) < sizeof(Int_t))
// 		     {
// 			UChar_t *tmp = new UChar_t[fAvailableSize*2];
// 			memcpy(tmp, fRecPointDataPtr, fAvailableSize);
// 			for(Int_t n = 0; n < fNRecPoints; n++)
// 			{
// 			   fRecPointArray[n] = reinterpret_cast<AliHLTCaloRecPointDataStruct*>(reinterpret_cast<UChar_t*>(fRecPointArray[n]) - reinterpret_cast<UChar_t*>(fFirstRecPointPtr) + reinterpret_cast<UChar_t*>(tmp));
// 			}
// 			fRecPointDataPtr = reinterpret_cast<AliHLTCaloRecPointDataStruct*>(tmp);
// 			fFirstRecPointPtr = fRecPointDataPtr;
// 			fUsedSize = 0;
// 		     }	
		  CheckBuffer();
		  // Assigning index to digit
		  printf("Digit index pointer: %x\n", fDigitIndexPtr);
		  *fDigitIndexPtr = j;
		  fUsedSize += sizeof(Int_t);
		  
		  printf("Clusterizier: adding digit:  index pointer: %x, index: %d\n", fDigitIndexPtr, *fDigitIndexPtr); 
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
      Int_t rowdiff = TMath::Abs( digit1->fZ - digit2->fZ );  
      Int_t coldiff = TMath::Abs( digit1->fX - digit2->fX ); 
      if (( coldiff <= 1   &&  rowdiff == 0 ) || ( coldiff == 0 &&  rowdiff <= 1 ))
	{
//	  cout << "Are neighbours: digit (E = "  << digit1->fEnergy << ") with x = " << digit1->fX << " and z = " << digit1->fZ << 
//	    " is neighbour with digit (E = " << digit2->fEnergy << ") with x = " << digit2->fX << " and z = " << digit2->fZ << endl;

	  if(TMath::Abs(digit1->fTime - digit2->fTime ) < fEmcTimeGate)
	    {
	      return 1; 
	    }
	}

   /*   Float_t rowdiff = TMath::Abs( digit1->fZ - digit2->fZ );  
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
   */
      else
	{
//	    cout << "Not neighbours: digit (E = "  << digit1->fEnergy << ") with x = " << digit1->fX << " and z = " << digit1->fZ << 
//	    " is not neighbour with digit (E = " << digit2->fEnergy << ") with x = " << digit2->fX << " and z = " << digit2->fZ << endl;
	}
    }
  return 0;
}



Int_t AliHLTCaloClusterizer::CheckArray()
{
      printf("CheckArray: fArraySize: %d, fNRecPoints: %d\n", fArraySize, fNRecPoints);
      if(fArraySize == fNRecPoints)
	{
	   printf("Expanding array...");
	   fArraySize *= 2;
	   AliHLTCaloRecPointDataStruct **tmp = new AliHLTCaloRecPointDataStruct*[fArraySize];
	   memcpy(tmp, fRecPointArray, fArraySize/2 * sizeof(AliHLTCaloRecPointDataStruct*));
	   delete fRecPointArray;
	   fRecPointArray = tmp;
	}
   return 0;
}

Int_t AliHLTCaloClusterizer::CheckBuffer()
{
   // See header file for class documentation 
	 printf("CheckBuffer: Used size %d, fAvailableSize: %d\n", fUsedSize, fAvailableSize);
	if((fAvailableSize - fUsedSize) < sizeof(AliHLTCaloRecPointDataStruct) )
	{
	   printf("Expanding buffer...\n");
	    Int_t recPointOffset = reinterpret_cast<UChar_t*>(fRecPointDataPtr) - reinterpret_cast<UChar_t*>(fFirstRecPointPtr);
	    fAvailableSize *= 2;
	    UChar_t *tmp = new UChar_t[fAvailableSize];
	    memcpy(tmp, fRecPointDataPtr, fAvailableSize/2);
	    for(Int_t n = 0; n < fNRecPoints; n++)
	    {
	       fRecPointArray[n] = reinterpret_cast<AliHLTCaloRecPointDataStruct*>(reinterpret_cast<UChar_t*>(fRecPointArray[n]) - reinterpret_cast<UChar_t*>(fFirstRecPointPtr) + reinterpret_cast<UChar_t*>(tmp));
	    }
	    delete fRecPointDataPtr;
	    fFirstRecPointPtr = reinterpret_cast<AliHLTCaloRecPointDataStruct*>(tmp);
	    fRecPointDataPtr = reinterpret_cast<AliHLTCaloRecPointDataStruct*>(tmp + recPointOffset);
	    fUsedSize = 0;
	}
   return 0;
}

Int_t AliHLTCaloClusterizer::CheckDigits(AliHLTCaloRecPointDataStruct** recArray, AliHLTCaloDigitDataStruct** digitArray, Int_t nRP)
{
  AliHLTCaloRecPointDataStruct **recpoints = recArray;
  AliHLTCaloDigitDataStruct **digits = digitArray;
  Int_t nRecPoints = nRP;
  
  if(recArray == 0)
  {
     recpoints = fRecPointArray;
  }
  if(digitArray == 0)
  {
     digits = fDigitsPointerArray;
  }
  if(nRP == 0)
  {
     nRecPoints = fNRecPoints;
  }
  printf("CL: CheckDigits: Number of rec points: %d\n", nRecPoints);
  for(Int_t i = 0; i < nRecPoints; i++)
  {
          
     AliHLTCaloRecPointDataStruct *recPoint = recpoints[i];

     //AliHLTCaloRecPointDataStruct *recPoint = fRecPointArray[0];
     Int_t multiplicity = recPoint->fMultiplicity;
     Int_t *digitIndexPtr = &(recPoint->fDigits);
     printf("CL: Rec point with energy: %f, multiplicity: %d\n", recPoint->fAmp, recPoint->fMultiplicity);
     for(Int_t j = 0; j < multiplicity; j++)
     {
	//AliHLTCaloRecPointDataStruct *recPoint = fRecPointArray[j];
	AliHLTCaloDigitDataStruct *digit = digits[*digitIndexPtr];
	printf("CL: Digit ID: %d, energy: %f, index: %d, indexpointer: %x\n", digit->fID, digit->fEnergy, *digitIndexPtr, digitIndexPtr);
	digitIndexPtr++;
	//recPoint = reinterpret_cast<AliHLTCaloRecPointDataStruct*>(digitIndexPtr);
     }
  }
     
     
     
     
     
     
     
}

Int_t AliHLTCaloClusterizer::CheckDigits(AliHLTCaloRecPointDataStruct** recArray, AliHLTCaloDigitDataStruct* digitArray, Int_t nRP)
{
  AliHLTCaloRecPointDataStruct **recpoints = recArray;
  AliHLTCaloDigitDataStruct *digits = digitArray;
  Int_t nRecPoints = nRP;
  
  if(recArray == 0)
  {
     recpoints = fRecPointArray;
  }
    if(nRP == 0)
  {
     nRecPoints = fNRecPoints;
  }
  printf("CL: CheckDigits: Number of rec points: %d\n", nRecPoints);
  for(Int_t i = 0; i < nRecPoints; i++)
  {
          
     AliHLTCaloRecPointDataStruct *recPoint = recpoints[i];

     //AliHLTCaloRecPointDataStruct *recPoint = fRecPointArray[0];
     Int_t multiplicity = recPoint->fMultiplicity;
     Int_t *digitIndexPtr = &(recPoint->fDigits);
     printf("CL: Rec point with energy: %f, multiplicity: %d\n", recPoint->fAmp, recPoint->fMultiplicity);
     for(Int_t j = 0; j < multiplicity; j++)
     {
	//AliHLTCaloRecPointDataStruct *recPoint = fRecPointArray[j];
	AliHLTCaloDigitDataStruct digit = digits[*digitIndexPtr];
	printf("CL: digits: %x, recpoints: %x, digitIndexPtr: %x\n", digits, recpoints, digitIndexPtr);
	printf("CL: Digit ID: %d, energy: %f, index: %d, indexpointer: %x\n", digit.fID, digit.fEnergy, *digitIndexPtr, digitIndexPtr);
	digitIndexPtr++;
	//recPoint = reinterpret_cast<AliHLTCaloRecPointDataStruct*>(digitIndexPtr);
     }
  }
     
     
     
     
     
     
     
}