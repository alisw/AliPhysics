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
 * @file   AliHLTPHOSPhysicsAnalyzerSpectrum.cxx
 * @author Oystein Djuvsland
 * @date 
 * @brief  Invariant mass spectrum from 2 gammas  */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTPHOSPhysicsAnalyzerSpectrum.h"
#include "AliHLTPHOSRecPointContainerStruct.h"
#include "AliHLTPHOSRecPointDataStruct.h"
#include <cmath>
#include "math.h"
#include "TH1F.h"


ClassImp(AliHLTPHOSPhysicsAnalyzerSpectrum);



AliHLTPHOSPhysicsAnalyzerSpectrum::AliHLTPHOSPhysicsAnalyzerSpectrum():AliHLTPHOSPhysicsAnalyzer(), fPos0Ptr(0), fPos1Ptr(0), fThresholdPtr(0), fEnergyPtr(0)
{
  //Constructor
  //See header file for documentation
  fEnergyPtr = new Float_t[2];
  fPos0Ptr = new Float_t[3];
  fPos1Ptr = new Float_t[3];
  fThresholdPtr = new Float_t[2];
  fThresholdPtr[0] = 0;
  fThresholdPtr[1] = 0;
  
}

AliHLTPHOSPhysicsAnalyzerSpectrum::AliHLTPHOSPhysicsAnalyzerSpectrum(const AliHLTPHOSPhysicsAnalyzerSpectrum &):AliHLTPHOSPhysicsAnalyzer(), fPos0Ptr(0), fPos1Ptr(0), fThresholdPtr(0), fEnergyPtr(0)
{
  //Copy constructor not implemented
  //See header file for documentation
}

AliHLTPHOSPhysicsAnalyzerSpectrum::~AliHLTPHOSPhysicsAnalyzerSpectrum()
{
  //Destructor
  //See header file for documentation
  if(fRecPointsPtr) fRecPointsPtr = 0;

  if(fRootHistPtr) fRootHistPtr = 0;

  if(fThresholdPtr)
    {
      delete [] fThresholdPtr;
      fThresholdPtr = 0;
    }
  if(fEnergyPtr)
    {
      delete [] fEnergyPtr;
      fEnergyPtr = 0;
    }
  if(fPos0Ptr)
    {
      delete [] fPos0Ptr;
      fPos0Ptr = 0;
    }
  if(fPos1Ptr)
    {
      delete [] fPos1Ptr;
      fPos1Ptr = 0;
    }

}
void
AliHLTPHOSPhysicsAnalyzerSpectrum::Analyze(AliHLTPHOSRecPointContainerStruct* recPointsArrayPtr, Int_t nRecPoints)
{
  //Analyzing a set of recPoints
  //See header file for documentation
  Float_t cosOpeningAngle = 0;

  AliHLTPHOSRecPointDataStruct* firstRecPointPtr;
  AliHLTPHOSRecPointDataStruct* secondRecPointPtr;

  if(nRecPoints > 1)
    {
      for(Int_t i = 0; i < nRecPoints-1; i++)
	{
	  firstRecPointPtr = &(recPointsArrayPtr->fRecPointArray[i]);
	  fEnergyPtr[0] = firstRecPointPtr->fAmp;
	  if(fEnergyPtr[0] > fThresholdPtr[0])
	    {
	      
	      for(Int_t j = i+1; j < nRecPoints; j++)
		{
		  secondRecPointPtr = &(recPointsArrayPtr->fRecPointArray[j]);
		  fEnergyPtr[1] = secondRecPointPtr->fAmp;
		  if(fEnergyPtr[1] > fThresholdPtr[1])
		    {
		      GlobalPosition(firstRecPointPtr, fPos0Ptr);
		      GlobalPosition(secondRecPointPtr, fPos1Ptr);

		      cosOpeningAngle = (fPos0Ptr[0]*fPos1Ptr[0] + fPos0Ptr[1]*fPos1Ptr[1] + fPos0Ptr[2]*fPos1Ptr[2])/
			(sqrt(fPos0Ptr[0]*fPos0Ptr[0] + fPos0Ptr[1]*fPos0Ptr[1] + fPos0Ptr[2]*fPos0Ptr[2])*
			 sqrt(fPos1Ptr[0]*fPos1Ptr[0] + fPos1Ptr[1]*fPos1Ptr[1] + fPos1Ptr[2]*fPos1Ptr[2]));
		      
		      fRootHistPtr->Fill(sqrt(2*fEnergyPtr[0]*fEnergyPtr[1]*(1 - cosOpeningAngle)));
		    }
		  
		}	   	  
	      
	    }
	  
	}
      
    }
}

Float_t 
AliHLTPHOSPhysicsAnalyzerSpectrum::EvalDistance()
{
  //Evaluate the distance between the two recPoints
  //See header file for documentation
  if(fPos0Ptr && fPos1Ptr)
    return sqrt(pow(fPos1Ptr[0]-fPos0Ptr[0],2) + pow(fPos1Ptr[1]-fPos0Ptr[1],2) + pow(fPos1Ptr[2]-fPos0Ptr[2],2));
  return -1;

}

Int_t
AliHLTPHOSPhysicsAnalyzerSpectrum::SetThreshold(Float_t photonEnergy0, Float_t photonEnergy1)
{
  //Setting the cut thresholds
  //See header file for documentation
  fThresholdPtr[0] = photonEnergy0;

  fThresholdPtr[1] = photonEnergy1;

  return 0;

}
