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
 * @file   AliHLTPHOSClusterAnalyser.cxx
 * @author Oystein Djuvsland
 * @date 
 * @brief  Cluster analyser for PHOS HLT 
 */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTPHOSClusterAnalyser.h"
#include "AliHLTPHOSRecPointContainerStruct.h"
#include "AliHLTPHOSRecPointDataStruct.h"
#include "AliHLTPHOSCaloClusterContainerStruct.h"
#include "AliHLTPHOSCaloClusterDataStruct.h"
#include "AliHLTPHOSPhysicsAnalyzerSpectrum.h"
#include "AliPHOSGeometry.h"
#include "TMath.h"

ClassImp(AliHLTPHOSClusterAnalyser);

AliHLTPHOSClusterAnalyser::AliHLTPHOSClusterAnalyser() :
  AliHLTPHOSBase(),
  fLogWeight(0),
  fRecPointsPtr(0),
  fCaloClustersPtr(0),
  fPHOSGeometry(0),
  fAnalyzerPtr(0),
  fDoClusterFit(false),
  fHaveCPVInfo(false),
  fDoPID(false),
  fHaveDistanceToBadChannel(false)
{
  //See header file for documentation
  fLogWeight = 4.5;

  fAnalyzerPtr = new AliHLTPHOSPhysicsAnalyzerSpectrum();
  fPHOSGeometry = AliPHOSGeometry::GetInstance("noCPV");
}

AliHLTPHOSClusterAnalyser::~AliHLTPHOSClusterAnalyser() 
{
}

void 
AliHLTPHOSClusterAnalyser::SetCaloClusterContainer(AliHLTPHOSCaloClusterContainerStruct *caloClusterContainerPtr)
{ 
  //see header file for documentation
  fCaloClustersPtr = caloClusterContainerPtr; 
  fCaloClustersPtr->fNCaloClusters = 0;
}

Int_t
AliHLTPHOSClusterAnalyser::CalculateCenterOfGravity()
{
  //see header file for documentation
  Float_t wtot = 0.;
  Float_t x = 0.;
  Float_t z = 0.;
  Float_t xi = 0.;
  Float_t zi = 0.;
  AliHLTPHOSRecPointDataStruct *recPoint = 0;
  AliHLTPHOSDigitDataStruct *digit = 0;
  //AliPHOSGeometry * phosgeom =  AliPHOSGeometry::GetInstance() ;

  UInt_t iDigit = 0;
  UInt_t iRecPoint = 0;

  for(iRecPoint=0; iRecPoint < fRecPointsPtr->fNRecPoints; iRecPoint++) 
    {
      recPoint = &(fRecPointsPtr->fRecPointArray[iRecPoint]);
      for(iDigit = 0; iDigit < recPoint->fMultiplicity; iDigit++)
	{
	  digit = &(recPoint->fDigitsList[iDigit]);
	  
	  xi = digit->fX;
	  zi = digit->fZ;
	  
	  if (recPoint->fAmp > 0 && digit->fEnergy > 0) 
	    {
	      Float_t w = TMath::Max( 0., fLogWeight + TMath::Log( digit->fEnergy / recPoint->fAmp ) ) ;
	      x    += xi * w ;
	      z    += zi * w ;
	      wtot += w ;
	    }
	}
      
      if (wtot>0) 
	{
	  recPoint->fX = x/wtot ;
	  recPoint->fZ = z/wtot ;
	}
      else
	{
	  recPoint->fAmp = 0;
	}
    }
  return 0;
}


Int_t 
AliHLTPHOSClusterAnalyser::CalculateRecPointMoments()
{
  //See header file for documentation
  return 0;
}

Int_t 
AliHLTPHOSClusterAnalyser::CalculateClusterMoments(AliHLTPHOSRecPointDataStruct */*recPointPtr*/, AliHLTPHOSCaloClusterDataStruct* /*clusterPtr*/)
{
  //See header file for documentation
  return 0;
}


Int_t 
AliHLTPHOSClusterAnalyser::DeconvoluteClusters()
{
  //See header file for documentation
  return 0;
}

Int_t 
AliHLTPHOSClusterAnalyser::CreateClusters()
{
  //See header file for documentation
  
  AliHLTPHOSRecPointDataStruct* recPointPtr = 0;
  AliHLTPHOSCaloClusterDataStruct* caloClusterPtr = 0;
  AliHLTPHOSDigitDataStruct* digitPtr = 0;

  Float_t localPos[2];
  Float_t globalPos[3];
  Int_t id = -1;

  //fPHOSGeometry = AliPHOSGeometry::GetInstance("noCPV");

  for(UInt_t i = 0; i < fRecPointsPtr->fNRecPoints; i++)
    {
      
      caloClusterPtr = &(fCaloClustersPtr->fCaloClusterArray[i]);
      recPointPtr = &(fRecPointsPtr->fRecPointArray[i]);
      
      localPos[0] = recPointPtr->fX;
      localPos[1] = recPointPtr->fZ;
      
      fAnalyzerPtr->GlobalPosition( localPos, globalPos, recPointPtr->fModule);

      caloClusterPtr->fGlobalPos[0] = globalPos[0];
      caloClusterPtr->fGlobalPos[1] = globalPos[1];
      caloClusterPtr->fGlobalPos[2] = globalPos[2];

      caloClusterPtr->fNCells = recPointPtr->fMultiplicity;
  
      for(UInt_t j = 0; j < caloClusterPtr->fNCells; j++)
	{
	  digitPtr = &(recPointPtr->fDigitsList[j]);
	  //fPHOSGeometry->RelPosToAbsId((Int_t)(recPointPtr->fModule + 1), (double)(digitPtr->fX), (double)(digitPtr->fZ), id);
	  caloClusterPtr->fCellsAbsId[j] = id;
	  caloClusterPtr->fCellsAmpFraction[j] = digitPtr->fAmplitude/recPointPtr->fAmp;
	}
      caloClusterPtr->fEnergy = recPointPtr->fAmp;
      
      if(fDoClusterFit)
	{
	  FitCluster(recPointPtr);
	}
      else
	{
	  caloClusterPtr->fDispersion = 0;
	  caloClusterPtr->fFitQuality = 0;
	  caloClusterPtr->fM20 = 0;
	  caloClusterPtr->fM02 = 0;
	  caloClusterPtr->fM11 = 0;
	}
      if(fHaveCPVInfo)
	{
	  caloClusterPtr->fEmcCpvDistance = GetCPVDistance(recPointPtr);
	}
      else
	{
	  caloClusterPtr->fEmcCpvDistance = -1;
	}
      if(fDoPID)
	{
	  DoParticleIdentification(caloClusterPtr);
	}
      else
	{
	  for(Int_t k = 0; k < AliPID::kSPECIESN; k++)
	    {
	      caloClusterPtr->fPID[k] = 0;
	    }
	}

      if(fHaveDistanceToBadChannel)
	{
	  caloClusterPtr->fDistanceToBadChannel = GetDistanceToBadChannel(caloClusterPtr);
	}
      else
	{
	  caloClusterPtr->fDistanceToBadChannel = -1;
	}

      caloClusterPtr->fClusterType = '\0';
    }
  fCaloClustersPtr->fNCaloClusters = fRecPointsPtr->fNRecPoints;

 
  return fCaloClustersPtr->fNCaloClusters;

}

