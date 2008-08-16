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
 * @file   AliHLTPHOSClusterizer.cxx
 * @author Oystein Djuvsland
 * @date 
 * @brief  Clusterizer for PHOS HLT 
 */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTPHOSClusterizer.h"
#include "AliHLTPHOSBase.h"
#include "TMath.h"
#include "AliHLTPHOSRecPointContainerStruct.h"
#include "AliHLTPHOSRecPointDataStruct.h"
#include "AliHLTPHOSDigitDataStruct.h"
#include "AliHLTPHOSDigitContainerDataStruct.h"
#include "TClonesArray.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSDigit.h"
#ifndef HAVE_NOT_PHOSRECOPARAMEMC // set from configure if EMC functionality not available in AliPHOSRecoParam
#include "AliPHOSRecoParam.h"
#else
#include "AliPHOSRecoParamEmc.h"
#endif

ClassImp(AliHLTPHOSClusterizer);

AliHLTPHOSClusterizer::AliHLTPHOSClusterizer():
  AliHLTPHOSBase(),
  fEmcClusteringThreshold(0),
  fEmcMinEnergyThreshold(0),
  fEmcTimeGate(0),
  fLogWeight(0),
  fDigitsInCluster(0),
  fOnlineMode(true),
  fDigitArrayPtr(0),
  fEmcRecPointsPtr(0),
  fDigitPtr(0),
  fDigitContainerPtr(0),
  fRecPointContainerPtr(0),
  fPHOSGeometry(0),
  fGetterPtr(0)
{
  //See header file for documentation
  fPHOSGeometry = AliPHOSGeometry::GetInstance("noCPV");
  fEmcClusteringThreshold = 0.2;
  fEmcMinEnergyThreshold = 0.03;
  fEmcTimeGate = 1.e-6 ;
  fLogWeight = 4.5;
  

}//end


AliHLTPHOSClusterizer::~AliHLTPHOSClusterizer()  
{
  //See header file for documentation
}

void 
AliHLTPHOSClusterizer::SetRecPointContainer(AliHLTPHOSRecPointContainerStruct* recPointContainerPtr)
  { 
    fRecPointContainerPtr = recPointContainerPtr; 
    fRecPointContainerPtr->fNRecPoints = 0;
  }

void
//AliHLTPHOSClusterizer::SetRecoParameters(AliPHOSRecoParam* params)
//CRAP PTH
AliHLTPHOSClusterizer::SetRecoParameters(AliPHOSRecoParam* /*params*/)
{
  //<<<<<<< HEAD:AliHLTPHOSClusterizer.cxx
  //=======



    //>>>>>>> e76e4d6d226ed1b3b26e0441a108cfe0e0bf0290:AliHLTPHOSClusterizer.cxx
  /*
  //see header file for documentation
#ifndef HAVE_NOT_PHOSRECOPARAMEMC // set from configure if EMC functionality not available in AliPHOSRecoParam
  // the new AliPHOSRecoParam functions, available from revision

    CRAP PTH
    temporarily commeneted out by PTH
    
  fEmcClusteringThreshold = params->GetEMCClusteringThreshold();
  fEmcMinEnergyThreshold = params->GetEMCMinE();
  fLogWeight = params->GetEMCLogWeight();
#else
  fEmcClusteringThreshold = params->GetClusteringThreshold();
  fEmcMinEnergyThreshold = params->GetMinE();
  fLogWeight = params->GetLogWeight();
#endif
  */
}

void
AliHLTPHOSClusterizer::SetOfflineMode(AliPHOSLoader* getter)
{
  //see header file for documentation
  fRecPointContainerPtr = new AliHLTPHOSRecPointContainerStruct();
  fDigitContainerPtr = new AliHLTPHOSDigitContainerDataStruct();
  fGetterPtr = getter;
  fDigitArrayPtr = fGetterPtr->Digits();
  fEmcRecPointsPtr = fGetterPtr->EmcRecPoints();
  fOnlineMode = false;
}

Int_t
AliHLTPHOSClusterizer::GetEvent(Int_t i)
{
  //see header file for documentation
  Int_t coord[4];

  //  fGetterPtr->Event(i, "D");
  fGetterPtr->GetEvent();
  for(Int_t j = 0; j < fDigitArrayPtr->GetEntries(); j++)
    {
      fDigitPtr = (AliPHOSDigit*)fDigitArrayPtr->At(j);
      fPHOSGeometry->AbsToRelNumbering(fDigitPtr->GetId(), coord);
      fDigitContainerPtr->fDigitDataStruct[j].fX = coord[3];
      fDigitContainerPtr->fDigitDataStruct[j].fZ = coord[2];
      fDigitContainerPtr->fDigitDataStruct[j].fModule = coord[0];
      fDigitContainerPtr->fDigitDataStruct[j].fEnergy = fDigitPtr->GetEnergy();
      fDigitContainerPtr->fDigitDataStruct[j].fTime = fDigitPtr->GetTime();
    }
  fDigitContainerPtr->fNDigits = fDigitArrayPtr->GetEntriesFast();
  return 0;
}

Int_t 
AliHLTPHOSClusterizer::GetNEvents()
{
  //see header file for documentation
  if(fOnlineMode)
    {
      //     Logging(kHLTLogWarning, __FILE__ , "information not available" , "GetNEvents()  Number of events not available in online mod");
      return -1;
    }
  return fGetterPtr->MaxEvent();
}
  

Int_t 
AliHLTPHOSClusterizer::ClusterizeEvent()
{
  //see header file for documentation
  Int_t nRecPoints = 0;
  UInt_t i = 0;

  AliHLTPHOSRecPointDataStruct *recPoint = 0;

  //Clusterization starts
  for(i = 0; i < fDigitContainerPtr->fNDigits; i++)
    { 
      
      fDigitsInCluster = 0;
      if(fDigitContainerPtr->fDigitDataStruct[i].fEnergy < fEmcClusteringThreshold)
	{
	  continue;
	}
      recPoint = &(fRecPointContainerPtr->fRecPointArray[nRecPoints]);
      recPoint->fAmp = 0;
      //TODO!!!!!!!
      recPoint->fModule = fDigitContainerPtr->fDigitDataStruct[i].fModule;
      //TODO!!!!!!!
      //recPoint->
      recPoint->fDigitsList[fDigitsInCluster] =  fDigitContainerPtr->fDigitDataStruct[i];
      recPoint->fAmp += fDigitContainerPtr->fDigitDataStruct[i].fEnergy;
      fDigitContainerPtr->fDigitDataStruct[i].fEnergy = 0;
      fDigitsInCluster++;
      nRecPoints++;
      ScanForNeighbourDigits(i, recPoint);

      recPoint->fMultiplicity = fDigitsInCluster;
      
    }//end of clusterization
  fRecPointContainerPtr->fNRecPoints = nRecPoints;
  
  return nRecPoints;
}

void
AliHLTPHOSClusterizer::ScanForNeighbourDigits(Int_t index, AliHLTPHOSRecPointDataStruct* recPoint)

{
  //see header file for documentation

  for(UInt_t j = 0; j < fDigitContainerPtr->fNDigits; j++)
    {
      if(fDigitContainerPtr->fDigitDataStruct[j].fEnergy > fEmcMinEnergyThreshold)
	{
	  switch(AreNeighbours(&(fDigitContainerPtr->fDigitDataStruct[index]),
			       &(fDigitContainerPtr->fDigitDataStruct[j])))
	    {
	    case 0:
	      break; 
	    case 1:      
	      recPoint->fDigitsList[fDigitsInCluster] =  fDigitContainerPtr->fDigitDataStruct[j];
	      recPoint->fAmp += fDigitContainerPtr->fDigitDataStruct[j].fEnergy;
	      fDigitContainerPtr->fDigitDataStruct[j].fEnergy = 0;
	      fDigitsInCluster++;
	      ScanForNeighbourDigits(j, recPoint);
	      break;
	    case 2:
	      break;
	    }
	}
    } 
  return;
}


Int_t 
AliHLTPHOSClusterizer::AreNeighbours(AliHLTPHOSDigitDataStruct* digit1, 
					    AliHLTPHOSDigitDataStruct* digit2)
{
  //see header file for documentation

  if ( (digit1->fModule == digit2->fModule) /*&& (coord1[1]==coord2[1])*/ ) // inside the same PHOS module
    { 
      Int_t rowdiff = TMath::Abs( digit1->fZ - digit2->fZ );  
      Int_t coldiff = TMath::Abs( digit1->fX - digit2->fX ); 
          
      if (( coldiff <= 1 )  && ( rowdiff <= 1 ))
	{
	  if(TMath::Abs(digit1->fTime - digit2->fTime ) < fEmcTimeGate)
	    {
	      return 1; 
	    }
	}
    }
  return 0;
}

void
AliHLTPHOSClusterizer::CalculateCenterOfGravity()
{
  //see header file for documentation
  Float_t wtot = 0.;
  Float_t x = 0.;
  Float_t z = 0.;
  Float_t xi = 0.;
  Float_t zi = 0.;

  AliHLTPHOSRecPointDataStruct *recPoint = 0;
  AliHLTPHOSDigitDataStruct *digit = 0;
  UInt_t iDigit = 0;
  UInt_t iRecPoint = 0;

  for(iRecPoint=0; iRecPoint<fRecPointContainerPtr->fNRecPoints; iRecPoint++) 
    {
      recPoint = &(fRecPointContainerPtr->fRecPointArray[iRecPoint]);
      for(iDigit = 0; iDigit < recPoint->fMultiplicity; iDigit++)
	{
	  digit = &(recPoint->fDigitsList[iDigit]);
	  
	  //fPHOSGeometry->AbsToRelNumbering(digit->fID, relid) ;
	  //  fPHOSGeometry->RelPosInModule(relid, xi, zi);
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
      
}


