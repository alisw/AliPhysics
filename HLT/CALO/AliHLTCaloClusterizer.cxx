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
        fCompareFunction(CompareDigitsByPosition),
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
        fNDigits(0),
        fSortedByPosition(false),
        fSortedByEnergy(false),
        fSortDigits(false)
{
    //See header file for documentation
    //fEmcClusteringThreshold = 0.2;
    //fEmcMinEnergyThreshold = 0.03;

    fEmcClusteringThreshold = 0.1;
    fEmcMinEnergyThreshold = 0.01;
    fEmcTimeGate = 1.e-6 ;

    fMaxDigitIndexDiff = 2*fCaloConstants->GetNZROWSMOD();


    fArraySize = 10;
    fRecPointArray = new AliHLTCaloRecPointDataStruct*[fArraySize];

    fAvailableSize = sizeof(AliHLTCaloRecPointDataStruct) * 20;
    fFirstRecPointPtr = reinterpret_cast<AliHLTCaloRecPointDataStruct*>(new UChar_t[fAvailableSize]);
    fRecPointDataPtr = fFirstRecPointPtr;

}//end

AliHLTCaloClusterizer::~AliHLTCaloClusterizer()
{
    //See header file for documentation
  delete [] fRecPointDataPtr;
  delete [] fRecPointArray;
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

    // Sort our digits
    SortDigits();

    //Clusterization starts
    for (Int_t i = 0; i < nDigits; i++)
    {
        fDigitsInCluster = 0;

	 HLTDebug("Digit with energy: %f", fDigitsPointerArray[i]->fEnergy);
	
        if (fDigitsPointerArray[i]->fEnergy < fEmcClusteringThreshold && fSortedByEnergy)
        {
	   // Since we have sorted by energy the next digit will have even lower energy, so we return 
	   return fNRecPoints;
	}

	if(fDigitsPointerArray[i]->fAssociatedCluster != -1)
	{
	   // The digit is added to a previous cluster, continue
	   continue;
	}

	CheckArray();
        CheckBuffer();

        // First digit is placed at the fDigits member variable in the recpoint
        fDigitIndexPtr = &(fRecPointDataPtr->fDigits);

        fRecPointDataPtr->fAmp = 0;
        fRecPointDataPtr->fModule = fDigitsPointerArray[i]->fModule;

        // Assigning the digit to this rec point
        fRecPointDataPtr->fDigits = i;
        fUsedSize += sizeof(AliHLTCaloRecPointDataStruct);

        // Incrementing the pointer to be ready for new entry
        fDigitIndexPtr++;

        fRecPointDataPtr->fAmp += fDigitsPointerArray[i]->fEnergy;
    
	
	//fDigitsPointerArray[i]->fEnergy = 0;
        fDigitsPointerArray[i]->fAssociatedCluster = fNRecPoints;
	
	
	fDigitsInCluster++;
        nRecPoints++;

        // Scanning for the neighbours
        if (ScanForNeighbourDigits(i, fRecPointDataPtr) != 0)
        {
            return -1;
        }

        //fUsedSize += sizeof(AliHLTCaloRecPointDataStruct) + (fDigitsInCluster-1)*sizeof(AliHLTCaloDigitDataStruct);

        fRecPointDataPtr->fMultiplicity = fDigitsInCluster;
        fRecPointArray[fNRecPoints] = fRecPointDataPtr;

        fRecPointDataPtr = reinterpret_cast<AliHLTCaloRecPointDataStruct*>(fDigitIndexPtr);

        fNRecPoints++;

    }//end of clusterization

    return nRecPoints;
}

Int_t
AliHLTCaloClusterizer::ScanForNeighbourDigits(Int_t index, AliHLTCaloRecPointDataStruct* recPoint)
{
    //see header file for documentation

    // The following cuts can be used if we sort by posisiton. Not tested, but it should be fine...
    Int_t max = TMath::Min(fNDigits, (Int_t)fMaxDigitIndexDiff+index);
    Int_t min = TMath::Max(0, (Int_t)(index - (Int_t)fMaxDigitIndexDiff));

    // All digits for now
    max = fNDigits;
    min = 0;

    for (Int_t j = min; j < max; j++)
    {
        if (fDigitsPointerArray[j]->fAssociatedCluster == -1 &&  fDigitsPointerArray[j]->fEnergy > fEmcMinEnergyThreshold)
        {
            if (j != index)
            {
                if (AreNeighbours(fDigitsPointerArray[index],
                                  fDigitsPointerArray[j]))
                {
                    // Check that the buffer is large enough for adding a digit (can be heavily improved wrt performance)
                    CheckBuffer();

                    // Assigning index to digit
                    *fDigitIndexPtr = j;
                    fUsedSize += sizeof(Int_t);

                    // Incrementing digit pointer to be ready for new entry
                    fDigitIndexPtr++;

                    // Adding the digit energy to the rec point
                    fRecPointDataPtr->fAmp += fDigitsPointerArray[j]->fEnergy;

                    // Setting energy to 0
		    //fDigitsPointerArray[j]->fEnergy = 0;
		    
		    // Setting the associated cluster 
		    fDigitsPointerArray[j]->fAssociatedCluster = fNRecPoints;
		    
		    HLTDebug("Added digit with index: %d, energy: %f, to associated cluster: %d", fDigitsPointerArray[j]->fID, fDigitsPointerArray[j]->fEnergy, fDigitsPointerArray[j]->fAssociatedCluster);
		    
                    fDigitsInCluster++;

                    // Scan for neighbours of this digit
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

        // As in the offline code we define neighbours as cells that share an edge, a corner is not  enough
	//        if (( coldiff <= 1   &&  rowdiff == 0 ) || ( coldiff == 0 &&  rowdiff <= 1 ))
        if (( coldiff <= 1) || ( rowdiff <= 1 ))
        {
            // Check also for time
            if (TMath::Abs(digit1->fTime - digit2->fTime ) < fEmcTimeGate)
            {
                return 1;
            }
        }
    }
    return 0;
}



Int_t AliHLTCaloClusterizer::CheckArray()
{
    // See header file for class documentation
    if (fArraySize == fNRecPoints)
    {
        fArraySize *= 2;
        AliHLTCaloRecPointDataStruct **tmp = new AliHLTCaloRecPointDataStruct*[fArraySize];
        memcpy(tmp, fRecPointArray, fArraySize/2 * sizeof(AliHLTCaloRecPointDataStruct*));
        delete [] fRecPointArray;
        fRecPointArray = tmp;
    }
    return 0;
}

Int_t AliHLTCaloClusterizer::CheckBuffer()
{
    // See header file for class documentation
    if ((fAvailableSize - fUsedSize) < (Int_t)sizeof(AliHLTCaloRecPointDataStruct))
    {
        Int_t recPointOffset = reinterpret_cast<UChar_t*>(fRecPointDataPtr) - reinterpret_cast<UChar_t*>(fFirstRecPointPtr);
        Int_t digitIndexOffset = reinterpret_cast<UChar_t*>(fDigitIndexPtr) - reinterpret_cast<UChar_t*>(fRecPointDataPtr);
        UChar_t *tmp = new UChar_t[fAvailableSize*2];

        memcpy(tmp, fFirstRecPointPtr, fUsedSize);
        fAvailableSize *= 2;
        for (Int_t n = 0; n < fNRecPoints; n++)
        {
            fRecPointArray[n] = reinterpret_cast<AliHLTCaloRecPointDataStruct*>(reinterpret_cast<UChar_t*>(fRecPointArray[n]) - reinterpret_cast<UChar_t*>(fFirstRecPointPtr) + reinterpret_cast<UChar_t*>(tmp));
        }
        delete [] fFirstRecPointPtr;
        fFirstRecPointPtr = reinterpret_cast<AliHLTCaloRecPointDataStruct*>(tmp);
        fRecPointDataPtr = reinterpret_cast<AliHLTCaloRecPointDataStruct*>(tmp + recPointOffset);
        fDigitIndexPtr = reinterpret_cast<Int_t*>(reinterpret_cast<UChar_t*>(fRecPointDataPtr) + digitIndexOffset);
        //fUsedSize = 0;
    }
    return 0;
}

void AliHLTCaloClusterizer::SetSortDigitsByPosition()
{
    // Sort the digit pointers by position
    fCompareFunction = &CompareDigitsByPosition;
    fSortDigits = true;
    fSortedByPosition = true;
}

void AliHLTCaloClusterizer::SetSortDigitsByEnergy()
{
    // See header file for class documentation
    fCompareFunction = &CompareDigitsByEnergy;
    fSortDigits = true;
    fSortedByEnergy = true;
}

void AliHLTCaloClusterizer::SortDigits()
{
    // See header file for class documentation
    if (fSortDigits) qsort(fDigitsPointerArray, fNDigits, sizeof(AliHLTCaloDigitDataStruct*), fCompareFunction);
}

Int_t
AliHLTCaloClusterizer::CompareDigitsByPosition(const void *dig0, const void *dig1)
{
    // See header file for documentation
    return (*((AliHLTCaloDigitDataStruct**)(dig0)))->fID - (*((AliHLTCaloDigitDataStruct**)(dig1)))->fID;
}

Int_t
AliHLTCaloClusterizer::CompareDigitsByEnergy(const void *dig0, const void *dig1)
{
    // See header file for documentation
  return (float)((*((AliHLTCaloDigitDataStruct**)(dig1)))->fEnergy - (*((AliHLTCaloDigitDataStruct**)(dig0)))->fEnergy);
}
