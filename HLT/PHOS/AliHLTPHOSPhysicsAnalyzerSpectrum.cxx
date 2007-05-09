
/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Øystein Djuvsland <oysteind@ift.uib.no>                       *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


#include "AliHLTPHOSPhysicsAnalyzerSpectrum.h"
#include "AliHLTPHOSClusterDataStruct.h"
#include <cmath>
#include "math.h"
#include <iostream>

ClassImp(AliHLTPHOSPhysicsAnalyzerSpectrum);

AliHLTPHOSPhysicsAnalyzerSpectrum::AliHLTPHOSPhysicsAnalyzerSpectrum():AliHLTPHOSPhysicsAnalyzer()
{

  fEnergyPtr = new Float_t[2];
  fPos0Ptr = new Float_t[3];
  fPos1Ptr = new Float_t[3];
  fThresholdPtr = new Float_t[2];
  fThresholdPtr[0] = 0;
  fThresholdPtr[1] = 0;
  
}

AliHLTPHOSPhysicsAnalyzerSpectrum::AliHLTPHOSPhysicsAnalyzerSpectrum(const AliHLTPHOSPhysicsAnalyzerSpectrum &):AliHLTPHOSPhysicsAnalyzer()
{

  cout << "AliHLTPHOSPhysicsAnalyzerSpectrum: Copy constructor not implemented yet!" << endl;

}

AliHLTPHOSPhysicsAnalyzerSpectrum::~AliHLTPHOSPhysicsAnalyzerSpectrum()
{

  if(fClustersPtr) fClustersPtr = 0;

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
AliHLTPHOSPhysicsAnalyzerSpectrum::Analyze(AliHLTPHOSClusterDataStruct* clustersPtr[10000], Int_t nClusters)
{

  Float_t cosOpeningAngle = 0;

  cout << "number of clusters: " << nClusters << endl;

  if(nClusters > 1)
    {
      for(Int_t i = 0; i < nClusters-1; i++)
	{
	  
	  fEnergyPtr[0] = clustersPtr[i]->fClusterEnergy;

	  if(fEnergyPtr[0] > fThresholdPtr[0])
	    {
	      
	      for(Int_t j = i+1; j < nClusters; j++)
		{
		  
		  fEnergyPtr[1] = clustersPtr[j]->fClusterEnergy;

		  if(fEnergyPtr[1] > fThresholdPtr[1])
		    {
		      GlobalPosition(clustersPtr[i], fPos0Ptr);

		      GlobalPosition(clustersPtr[j], fPos1Ptr);
		      
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
  
  return sqrt(pow(fPos1Ptr[0]-fPos0Ptr[0],2) + pow(fPos1Ptr[1]-fPos0Ptr[1],2) + pow(fPos1Ptr[2]-fPos0Ptr[2],2));

}

Float_t 
AliHLTPHOSPhysicsAnalyzerSpectrum::EvalCutDistance(Float_t cutMass)
{

  Float_t cosCutOpeningAngle = 0;
  
  cosCutOpeningAngle = 1 - cutMass*cutMass/(2*fEnergyPtr[0]*fEnergyPtr[1]);

  return 1.5*2*sin(acos(cosCutOpeningAngle));
  
}

Int_t
AliHLTPHOSPhysicsAnalyzerSpectrum::SetThreshold(Float_t photonEnergy0, Float_t photonEnergy1)
{
  
  if(!fThresholdPtr) fThresholdPtr = new Float_t[2];

  fThresholdPtr[0] = photonEnergy0;

  fThresholdPtr[1] = photonEnergy1;

  return 0;

}
