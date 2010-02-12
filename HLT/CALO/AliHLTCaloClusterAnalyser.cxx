// $Id: AliHLTCaloClusterAnalyser.cxx 35107 2009-09-30 01:45:06Z phille $

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
 * @file   AliHLTCaloClusterAnalyser.cxx
 * @author Oystein Djuvsland
 * @date 
 * @brief  Cluster analyser for Calo HLT 
 */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTCaloClusterAnalyser.h"
#include "AliHLTCaloRecPointHeaderStruct.h"
#include "AliHLTCaloRecPointDataStruct.h"
#include "AliHLTCaloClusterDataStruct.h"
//#include "AliHLTCaloPhysicsAnalyzer.h"
#include "AliHLTCaloGeometry.h"
#include "AliESDCaloCluster.h"
#include "TMath.h"
#include "TVector3.h"
#include "TH1F.h"
#include "TFile.h"
#include "AliHLTCaloClusterizer.h"

ClassImp(AliHLTCaloClusterAnalyser);

AliHLTCaloClusterAnalyser::AliHLTCaloClusterAnalyser() : 
  //  AliHLTCaloBase(),
  fLogWeight(4.5),
  fRecPointArray(0),
  fDigitDataArray(0),
  fNRecPoints(0),
  fCaloClusterDataPtr(0),
  fCaloClusterHeaderPtr(0),
  //fAnalyzerPtr(0),
  fDoClusterFit(false),
  fHaveCPVInfo(false),
  fDoPID(false),
  fHaveDistanceToBadChannel(false),
  fGeometry(0),
  fClusterType(AliESDCaloCluster::kPHOSCluster)
{
  //See header file for documentation
}

AliHLTCaloClusterAnalyser::~AliHLTCaloClusterAnalyser() 
{
  // See header file for class documentation
}

void 
AliHLTCaloClusterAnalyser::SetCaloClusterData(AliHLTCaloClusterDataStruct *caloClusterDataPtr)
{ 
  //see header file for documentation
  fCaloClusterDataPtr = caloClusterDataPtr; 
}

void
AliHLTCaloClusterAnalyser::SetRecPointArray(AliHLTCaloRecPointDataStruct **recPointDataPtr, Int_t nRecPoints)
{ 
  fRecPointArray = recPointDataPtr; 
  fNRecPoints = nRecPoints;
}

void 
AliHLTCaloClusterAnalyser::SetDigitDataArray(AliHLTCaloDigitDataStruct *digits) 
{ 
//   AliHLTCaloClusterizer cl("PHOS");
  // cl.CheckDigits(fRecPointArray, digits, fNRecPoints);
   fDigitDataArray = digits; 
   //cl.CheckDigits(fRecPointArray, fDigitDataArray, fNRecPoints);
}

Int_t
AliHLTCaloClusterAnalyser::CalculateCenterOfGravity()
{
  //see header file for documentation
  Float_t wtot = 0.;
  Float_t x = 0.;
  Float_t z = 0.;
  Float_t xi = 0.;
  Float_t zi = 0.;

  AliHLTCaloDigitDataStruct *digit = 0;

  UInt_t iDigit = 0;

  for(Int_t iRecPoint=0; iRecPoint < fNRecPoints; iRecPoint++) 
    {
      AliHLTCaloRecPointDataStruct *recPoint = fRecPointArray[iRecPoint];
      //      digit = &(recPoint->fDigits);

      Int_t *digitIndexPtr = &(recPoint->fDigits);

      for(iDigit = 0; iDigit < recPoint->fMultiplicity; iDigit++)
	{

	  digit = &(fDigitDataArray[*digitIndexPtr]);

	  xi = digit->fX;
	  zi = digit->fZ;

	  if (recPoint->fAmp > 0 && digit->fEnergy > 0) 
	    {
	      Float_t w = TMath::Max( 0., fLogWeight + TMath::Log( digit->fEnergy / recPoint->fAmp ) ) ;
	      x    += xi * w ;
	      z    += zi * w ;
	      wtot += w ;
	    }
	  digitIndexPtr++;
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
AliHLTCaloClusterAnalyser::CalculateRecPointMoments()
{
  //See header file for documentation
  return 0;
}

Int_t 
AliHLTCaloClusterAnalyser::CalculateClusterMoments(AliHLTCaloRecPointDataStruct */*recPointPtr*/, AliHLTCaloClusterDataStruct* /*clusterPtr*/)
{
  //See header file for documentation
  return 0;
}


Int_t 
AliHLTCaloClusterAnalyser::DeconvoluteClusters()
{
  //See header file for documentation
  return 0;
}

Int_t 
AliHLTCaloClusterAnalyser::CreateClusters(Int_t nRecPoints, UInt_t availableSize, UInt_t& totSize)
{
  //See header file for documentation

   
  totSize += sizeof(AliHLTCaloClusterDataStruct);
  fNRecPoints = nRecPoints;

  if(fGeometry == 0)
  {
     HLTError("No geometry object is initialised, creation of clusters stopped");
  }

  CalculateCenterOfGravity();

  //  AliHLTCaloDigitDataStruct* digitPtr = &(recPointPtr->fDigits);  
  AliHLTCaloDigitDataStruct* digitPtr = 0;

  AliHLTCaloClusterDataStruct* caloClusterPtr = 0;
  UShort_t* cellIDPtr = 0;
  Float_t* cellAmpFracPtr = 0;;
  
//  Int_t id = -1;
  TVector3 globalPos;

  for(Int_t i = 0; i < fNRecPoints; i++) //TODO needs fix when we start unfolding (number of clusters not necessarily same as number of recpoints gotten from the clusterizer
    {
      if((availableSize - totSize)  < sizeof(AliHLTCaloClusterDataStruct))
      {
	 HLTError("Out of buffer");
	 return -ENOBUFS;
      }
      
      caloClusterPtr = fCaloClusterDataPtr;
     
      cellIDPtr = &(caloClusterPtr->fCellsAbsId);
      cellAmpFracPtr = &(caloClusterPtr->fCellsAmpFraction);
     
      AliHLTCaloRecPointDataStruct *recPointPtr = fRecPointArray[i];
      
      AliHLTCaloGlobalCoordinate globalCoord;
      fGeometry->GetGlobalCoordinates(*recPointPtr, globalCoord);

      caloClusterPtr->fGlobalPos[0] = globalCoord.fX;
      caloClusterPtr->fGlobalPos[1] = globalCoord.fY;
      caloClusterPtr->fGlobalPos[2] = globalCoord.fZ;

      caloClusterPtr->fNCells = 0;//recPointPtr->fMultiplicity;
      
      Int_t tmpSize = totSize + (caloClusterPtr->fNCells-1)*(sizeof(Short_t) + sizeof(Float_t));
      
      if((availableSize - totSize)  < tmpSize)
      {
	 HLTError("Out of buffer");
	 return -ENOBUFS;
      }
      
      for(UInt_t j = 0; j < caloClusterPtr->fNCells; j++)
	{
/*	  fGeometry->GetCellAbsId(recPointPtr->fModule, digitPtr->fX, digitPtr->fZ, id);
	  *cellIDPtr = id;
	  *cellAmpFracPtr = digitPtr->fEnergy/recPointPtr->fAmp;
	  digitPtr++;
	  cellIDPtr = reinterpret_cast<UShort_t*>(reinterpret_cast<char*>(cellAmpFracPtr) + sizeof(Float_t)); 
	  cellAmpFracPtr = reinterpret_cast<Float_t*>(reinterpret_cast<char*>(cellIDPtr) + sizeof(Short_t)); */
	}

      totSize += tmpSize;

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

      caloClusterPtr->fClusterType = fClusterType;
      //      totSize += sizeof(AliHLTCaloClusterDataStruct) + (caloClusterPtr->fNCells)*(sizeof(Short_t) +sizeof(Float_t)-1);   
      //totSize += sizeof(AliHLTCaloClusterDataStruct) + (caloClusterPtr->fNCells-1)*(sizeof(Short_t) + sizeof(Float_t));   

      //      caloClusterPtr = reinterpret_cast<AliHLTCaloClusterDataStruct*>(cellAmpFracPtr);
      caloClusterPtr = reinterpret_cast<AliHLTCaloClusterDataStruct*>(cellIDPtr);

      recPointPtr = reinterpret_cast<AliHLTCaloRecPointDataStruct*>(digitPtr);
      //digitPtr = &(recPointPtr->fDigits);  
    }

return fNRecPoints;

}

