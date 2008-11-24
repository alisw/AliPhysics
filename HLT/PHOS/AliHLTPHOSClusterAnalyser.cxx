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
#include "AliHLTPHOSRecPointHeaderStruct.h"
#include "AliHLTPHOSRecPointDataStruct.h"
#include "AliHLTPHOSCaloClusterHeaderStruct.h"
#include "AliHLTPHOSCaloClusterDataStruct.h"
#include "AliHLTPHOSPhysicsAnalyzer.h"
#include "AliPHOSGeometry.h"
#include "TMath.h"

ClassImp(AliHLTPHOSClusterAnalyser);

AliHLTPHOSClusterAnalyser::AliHLTPHOSClusterAnalyser() :
  AliHLTPHOSBase(),
  fLogWeight(0),
  fRecPointDataPtr(0),
  fNRecPoints(0),
  fCaloClusterDataPtr(0),
  fCaloClusterHeaderPtr(0),
  fPHOSGeometry(0),
  fAnalyzerPtr(0),
  fDoClusterFit(false),
  fHaveCPVInfo(false),
  fDoPID(false),
  fHaveDistanceToBadChannel(false)
{
  //See header file for documentation
  fLogWeight = 4.5;

  fAnalyzerPtr = new AliHLTPHOSPhysicsAnalyzer();
  fPHOSGeometry = AliPHOSGeometry::GetInstance("noCPV");
}

AliHLTPHOSClusterAnalyser::~AliHLTPHOSClusterAnalyser() 
{
}

void 
AliHLTPHOSClusterAnalyser::SetCaloClusterDataPtr(AliHLTPHOSCaloClusterDataStruct *caloClusterDataPtr)
{ 
  //see header file for documentation
  fCaloClusterDataPtr = caloClusterDataPtr; 
}
void
AliHLTPHOSClusterAnalyser::SetRecPointDataPtr(AliHLTPHOSRecPointHeaderStruct *recPointDataPtr)
{ 
  fNRecPoints = recPointDataPtr->fNRecPoints;
  fRecPointDataPtr = reinterpret_cast<AliHLTPHOSRecPointDataStruct*>(reinterpret_cast<Char_t*>(recPointDataPtr)+sizeof(AliHLTPHOSRecPointHeaderStruct)); 
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

  AliHLTPHOSDigitDataStruct *digit = 0;
  //AliPHOSGeometry * phosgeom =  AliPHOSGeometry::GetInstance() ;

  AliHLTPHOSRecPointDataStruct *recPoint = fRecPointDataPtr;

  UInt_t iDigit = 0;

  for(Int_t iRecPoint=0; iRecPoint < fNRecPoints; iRecPoint++) 
    {
      digit = &(recPoint->fDigits);
      for(iDigit = 0; iDigit < recPoint->fMultiplicity; iDigit++)
	{
	  
	  xi = digit->fX;
	  zi = digit->fZ;
	  //cout << "COG digits (x:z:amp:time): " << xi << " : " << zi << " : " << digit->fEnergy << " : " << digit->fTime << endl;
	  if (recPoint->fAmp > 0 && digit->fEnergy > 0) 
	    {
	      Float_t w = TMath::Max( 0., fLogWeight + TMath::Log( digit->fEnergy / recPoint->fAmp ) ) ;
	      x    += xi * w ;
	      z    += zi * w ;
	      wtot += w ;
	    }
	  digit++;
	}
      //cout << endl;
      if (wtot>0) 
	{
	  recPoint->fX = x/wtot ;
	  recPoint->fZ = z/wtot ;
	}
      else
	{
	  recPoint->fAmp = 0;
	}
      recPoint = reinterpret_cast<AliHLTPHOSRecPointDataStruct*>(digit);
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
AliHLTPHOSClusterAnalyser::CreateClusters(UInt_t availableSize, UInt_t& totSize)
{
  //See header file for documentation

  UInt_t maxClusterSize = sizeof(AliHLTPHOSCaloClusterDataStruct) + (6 << 7); //Reasonable estimate... (6 = sizeof(Short_t) + sizeof(Float_t)

  AliHLTPHOSRecPointDataStruct* recPointPtr = fRecPointDataPtr;
  AliHLTPHOSDigitDataStruct* digitPtr = &(recPointPtr->fDigits);  
 
  AliHLTPHOSCaloClusterDataStruct* caloClusterPtr = fCaloClusterDataPtr;
  UShort_t* cellIDPtr = &(caloClusterPtr->fCellsAbsId);
  Float_t* cellAmpFracPtr = &(caloClusterPtr->fCellsAmpFraction);
  
  Float_t localPos[2];

  Float_t globalPos[3];
  Int_t id = -1;
  
  //fPHOSGeometry = AliPHOSGeometry::GetInstance("noCPV");

  for(Int_t i = 0; i < fNRecPoints; i++) //FIXME needs fix when we start unfolding (number of clusters not necessarily same as number of recpoints gotten from the clusterizer
    {
      
      if(availableSize < (totSize + maxClusterSize)) 
	{
	  return -1; //Might get out of buffer, exiting
	}
      localPos[0] = recPointPtr->fX;
      localPos[1] = recPointPtr->fZ;
      
      fAnalyzerPtr->GlobalPosition( localPos, globalPos, recPointPtr->fModule);

      caloClusterPtr->fGlobalPos[0] = globalPos[0];
      caloClusterPtr->fGlobalPos[1] = globalPos[1];
      caloClusterPtr->fGlobalPos[2] = globalPos[2];

      //cout << "Local Position (x:z:module): " << localPos[0] << " : "<< localPos[1] << " : " << recPointPtr->fModule << endl;

      //cout << "Global Position (x:y:z): " << globalPos[0] << " : "<< globalPos[1] << " : " << globalPos[2] << endl << endl;

      caloClusterPtr->fNCells = recPointPtr->fMultiplicity;
  
      cellIDPtr = &(caloClusterPtr->fCellsAbsId);
      cellAmpFracPtr = &(caloClusterPtr->fCellsAmpFraction);
     
      for(UInt_t j = 0; j < caloClusterPtr->fNCells; j++)
	{
	  // fPHOSGeometry->RelPosToAbsId((Int_t)(recPointPtr->fModule + 1), (double)(digitPtr->fX), (double)(digitPtr->fZ), id);
	  *cellIDPtr = id;
	  *cellAmpFracPtr = digitPtr->fEnergy/recPointPtr->fAmp;
	  digitPtr++;
	  cellIDPtr = reinterpret_cast<UShort_t*>(reinterpret_cast<char*>(cellAmpFracPtr) + sizeof(Float_t)); 
	  cellAmpFracPtr = reinterpret_cast<Float_t*>(reinterpret_cast<char*>(cellIDPtr) + sizeof(Short_t)); 
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
	  //	  caloClusterPtr->fM11 = 0;
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
      totSize += sizeof(AliHLTPHOSCaloClusterDataStruct) + (caloClusterPtr->fNCells)*(sizeof(Short_t) +sizeof(Float_t)-1);   

      caloClusterPtr = reinterpret_cast<AliHLTPHOSCaloClusterDataStruct*>(cellAmpFracPtr);
      recPointPtr = reinterpret_cast<AliHLTPHOSRecPointDataStruct*>(digitPtr);
      digitPtr = &(recPointPtr->fDigits);  
    }
 
  return fNRecPoints;

}

