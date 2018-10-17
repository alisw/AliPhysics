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
#include "AliHLTCaloRecoParamHandler.h"

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
  fClusterType(AliVCluster::kEMCALClusterv1),
  fRecoParamsPtr(0),
  fCutOnSingleCellClusters(false),
  fSingleCellEnergyCut(0.5)
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
AliHLTCaloClusterAnalyser::SetDigitDataArray(AliHLTCaloDigitDataStruct **digits) 
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
       wtot = 0;
       x = 0; 
       z = 0;

/*       Float_t maxAmp = 0;
       Int_t maxX = 0;
       Int_t maxZ = 0;*/
       if (fDigitDataArray[*digitIndexPtr])

	 for(iDigit = 0; iDigit < recPoint->fMultiplicity; iDigit++)
	   {
	     
	     digit = fDigitDataArray[*digitIndexPtr];
	     
	     xi = digit->fX;
	     zi = digit->fZ;
	     
	     //xi = digit->fX+0.5;
	     //zi = digit->fZ+0.5;
	     
	     if (recPoint->fAmp > 0 && digit->fEnergy > 0) 
	       {
		 Float_t w = TMath::Max( 0., fLogWeight + TMath::Log( digit->fEnergy / recPoint->fAmp ) ) ;
		 x    += xi * w ;
		 z    += zi * w ;
		 wtot += w ;
		 /*	      if(digit->fEnergy > maxAmp)
			      {
			      maxAmp = digit->fEnergy;
			      maxX = digit->fX;// + 0.5;
			      maxZ = digit->fZ;// + 0.5;
			      }*/
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
	   recPoint->fX = -9999;
	   recPoint->fZ =-9999;
	   // no good crashes depth with FP exception
	   //recPoint->fAmp = 0;
	 }
//	printf("Max digit: E = %f, x = %d, z= %d, cluster: E = %f, x = %f, z = %f\n" , maxAmp, maxX, maxZ, recPoint->fAmp, recPoint->fX, recPoint->fZ);
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

   if(fRecoParamsPtr)
   {
      fLogWeight = fRecoParamsPtr->GetLogWeight();
   }
   
  fNRecPoints = nRecPoints;

  if(fGeometry == 0)
  {
     HLTError("No geometry object is initialised, creation of clusters stopped");
     return  -1;
  }

  AliHLTCaloClusterDataStruct* caloClusterPtr = 0;

  TVector3 globalPos;
  for(Int_t i = 0; i < fNRecPoints; i++)
  {
    if((availableSize - totSize)  < sizeof(AliHLTCaloClusterDataStruct))
    {
        HLTError("Out of buffer: available size is: %d, total size used: %d", availableSize, totSize);
        return -ENOBUFS;
    }
    
    AliHLTCaloRecPointDataStruct *recPointPtr = fRecPointArray[i];

    if(fCutOnSingleCellClusters && recPointPtr->fAmp > fSingleCellEnergyCut && recPointPtr->fMultiplicity == 1) continue;
    
    totSize += sizeof(AliHLTCaloClusterDataStruct);
    caloClusterPtr = fCaloClusterDataPtr;
    caloClusterPtr->fChi2 = 0;
    caloClusterPtr->fClusterType = kUndef;
    caloClusterPtr->fDispersion = 0;
    caloClusterPtr->fDistanceToBadChannel = 0;
    caloClusterPtr->fDistToBadChannel = 0;
    caloClusterPtr->fEmcCpvDistance = 0;
    caloClusterPtr->fEnergy = 0;
    caloClusterPtr->fFitQuality = 0;
    caloClusterPtr->fID = recPointPtr->fLeadingCellID;
    caloClusterPtr->fM02 = recPointPtr->fM2x;
    caloClusterPtr->fM20 = recPointPtr->fM2z;
    caloClusterPtr->fNCells = 0;
    caloClusterPtr->fNExMax = 0;
    caloClusterPtr->fTOF = 0;
    caloClusterPtr->fTrackDx = -999;
    caloClusterPtr->fTrackDz = -999;
    
    AliHLTCaloGlobalCoordinate globalCoord;

    // 0 = assume photon
    fGeometry->GetGlobalCoordinates(*recPointPtr, globalCoord, 0);

    caloClusterPtr->fModule = recPointPtr->fModule;
    caloClusterPtr->fGlobalPos[0] =  globalCoord.fX;
    caloClusterPtr->fGlobalPos[1] =  globalCoord.fY;
    caloClusterPtr->fGlobalPos[2] =  globalCoord.fZ;

    HLTDebug("Cluster local position: x = %f, z = %f", recPointPtr->fX, recPointPtr->fZ);
    HLTDebug("Cluster global position: x = %f, y = %f, z = %f", globalCoord.fX, globalCoord.fY, globalCoord.fZ);
    
    caloClusterPtr->fNCells = recPointPtr->fMultiplicity;
    caloClusterPtr->fClusterType = fClusterType;

    if(fRecoParamsPtr)
    {
        caloClusterPtr->fEnergy = fRecoParamsPtr->GetCorrectedEnergy(recPointPtr->fAmp);
    }
    else
    {
        caloClusterPtr->fEnergy = recPointPtr->fAmp;
    }

    // Set the time of the maximum digit as cluster time
    caloClusterPtr->fTOF = recPointPtr->fTime;
    
    HLTDebug("Cluster global position: x = %f, y = %f, z = %f, energy: %f, time: %f, number of cells: %d, cluster pointer: %x", globalCoord.fX, globalCoord.fY, globalCoord.fZ, caloClusterPtr->fEnergy, caloClusterPtr->fTOF, caloClusterPtr->fNCells,  caloClusterPtr);

    for(Int_t k = 0; k < AliPID::kSPECIESCN; k++)
      caloClusterPtr->fPID[k] = 0;

    caloClusterPtr->fClusterType = fClusterType;
    fCaloClusterDataPtr++;
  }
  return fNRecPoints;
}

