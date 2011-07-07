

/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Oystein Djuvsland <oysteind@ift.uib.no>                       *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliHLTCaloClusterizerNbyN.h"

/**
 * @file   AliHLTCaloClusterizerNbyN.cxx
 * @author Oystein Djuvsland
 * @date
 * @brief  Clusterizer for PHOS HLT
 */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
 

AliHLTCaloClusterizerNbyN::AliHLTCaloClusterizerNbyN(TString det) : AliHLTCaloClusterizer(det)
,fN(3)
{
// Constructor
}
AliHLTCaloClusterizerNbyN::~AliHLTCaloClusterizerNbyN()
{
// Destructor
}

Int_t AliHLTCaloClusterizerNbyN::ClusterizeEvent(Int_t nDigits)
{

    //see header file for documentation
    Int_t nRecPoints = 0;
    fNRecPoints = 0;
    fUsedSize = 0;
    fNDigits = nDigits;
    fRecPointDataPtr = fFirstRecPointPtr;

    // Sort our digits
    fSortedByEnergy = true;
    SortDigits();

    //Clusterization starts
    for (Int_t i = 0; i < nDigits; i++)
    {
        fDigitsInCluster = 0;

        HLTDebug("Digit with energy: %f", fDigitsPointerArray[i]->fEnergy);

        if (fDigitsPointerArray[i]->fEnergy < fEmcClusteringThreshold)
        {
            // Since we have sorted by energy the next digit will have even lower energy, so we return
            return fNRecPoints;
        }

        if (fDigitsPointerArray[i]->fAssociatedCluster != -1) // Digit is neighbour with a higher energy digit
        {
            continue;
        }

        // Check if we enough space to write to, if not we expand it automatically.
        CheckArray();
        CheckBuffer();

        // Create the rec point and add the digit.
        // First digit is placed at the fDigits member variable in the recpoint
        fDigitIndexPtr = &(fRecPointDataPtr->fDigits);

        fRecPointDataPtr->fAmp = 0;
        fRecPointDataPtr->fModule = fDigitsPointerArray[i]->fModule;
	fRecPointDataPtr->fTime = fDigitsPointerArray[i]->fTime;

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


        // Then we loop over all the other digits (stupid, I know...)
        Int_t maxDiff = fN/2;
        for (Int_t j = 0; j < nDigits; j++)
        {
            if (fDigitsPointerArray[j]->fEnergy < fEmcMinEnergyThreshold) break; // Sorted by energy

            if (TMath::Abs(fDigitsPointerArray[i]->fX - fDigitsPointerArray[j]->fX) <= maxDiff
                    && TMath::Abs(fDigitsPointerArray[i]->fZ - fDigitsPointerArray[j]->fZ) <= maxDiff) // The digit is in our grid
            {
                if (TMath::Abs(fDigitsPointerArray[i]->fX - fDigitsPointerArray[j]->fX) == 1
                        || TMath::Abs(fDigitsPointerArray[i]->fZ - fDigitsPointerArray[j]->fZ) == 1) // The digit neighbour to the seed
                {
                    // This means the digit is not a local maximum
                    fDigitsPointerArray[j]->fAssociatedCluster = fNRecPoints;
                
		    // Check that the buffer is large enough for adding a digit (can be heavily improved wrt performance)
		    CheckBuffer();

		    // Assigning index to digit
		    *fDigitIndexPtr = j;
		    fUsedSize += sizeof(Int_t);

		    // Incrementing digit pointer to be ready for new entry
		    fDigitIndexPtr++;

		    // Adding the digit energy to the rec point
		    fRecPointDataPtr->fAmp += fDigitsPointerArray[j]->fEnergy;
		
		    // Count it
		    fDigitsInCluster++;
		}
            }
        }
  
        fRecPointDataPtr->fMultiplicity = fDigitsInCluster;
        fRecPointArray[fNRecPoints] = fRecPointDataPtr;

        fRecPointDataPtr = reinterpret_cast<AliHLTCaloRecPointDataStruct*>(fDigitIndexPtr);

        fNRecPoints++;

    }//end of clusterization

    return nRecPoints;
}

