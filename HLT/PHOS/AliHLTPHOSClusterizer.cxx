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


/** @file   AliHLTPHOSClusterizer.cxx
    @author Øystein Djuvsland
    @date   
    @brief  A temporary clusterizer for PHOS
*/



#include "AliHLTPHOSPhysicsDefinitions.h"
#include "AliHLTPHOSClusterizer.h"
#include "AliHLTPHOSCommonDefs.h"
#include "TVector3.h"
#include "TMath.h"
#include "AliHLTPHOSRcuCellEnergyDataStruct.h"
#include "AliHLTPHOSRecPointListDataStruct.h"
#include "AliHLTPHOSValidCellDataStruct.h"
#include "AliHLTPHOSRecPointDataStruct.h"
#include "AliHLTPHOSClusterDataStruct.h"


ClassImp(AliHLTPHOSClusterizer);

/**
* Main constructor
**/
AliHLTPHOSClusterizer::AliHLTPHOSClusterizer():fPHOSModule(-1), fThreshold(0), fClusterThreshold(0), 
					       fHighGainFactor(0.005), fLowGainFactor(0.08),
					       fArraySize(3), fMultiplicity(fArraySize*fArraySize)
{
  //Main constructor
  
}//end

AliHLTPHOSClusterizer::AliHLTPHOSClusterizer(const AliHLTPHOSClusterizer &):fPHOSModule(-1), fThreshold(0), fClusterThreshold(0), 
									    fHighGainFactor(0.005), fLowGainFactor(0.08),
									    fArraySize(3), fMultiplicity(fArraySize*fArraySize)
{
  //Copy constructor, not implemented
}//end

AliHLTPHOSClusterizer:: ~AliHLTPHOSClusterizer()  
{
  //Destructor
}

/**
* Building a 2D array of cell energies of the PHOS detector
* @param cellData object containing the cell energies from one event
* @param recPointList list to be filled with coordinates of local maxima in the detector
**/

Int_t
AliHLTPHOSClusterizer::BuildCellEnergyArray(AliHLTPHOSRcuCellEnergyDataStruct* cellData, 
					    AliHLTPHOSRecPointListDataStruct* recPointList)
{
  //Build the cell energy array of the detector

  Int_t x = 0;
  Int_t z = 0;
  Int_t gain = 0;
  Int_t xMod = 0;
  Int_t zMod = 0;
  Int_t index = 0;
  Double_t energyCount = 0;

  for(Int_t cell = 0; cell < cellData->fCnt; cell++)
    {
     
      z = (cellData->fValidData[cell]).fZ;
      x = (cellData->fValidData[cell]).fX;
      gain = (cellData->fValidData[cell]).fGain;

      zMod = z + (cellData->fRcuZ)*N_ROWS_RCU;
      xMod = x + (cellData->fRcuX)*N_COLUMNS_RCU;
      
      energyCount = (cellData->fValidData[cell]).fEnergy;

      if(gain == 0 && energyCount < 1023) 
	{
	  fEnergyArray[xMod][zMod] = fHighGainFactor * energyCount;
	 
	  if(fEnergyArray[xMod][zMod] > fClusterThreshold)
	    {	      
	      recPointList[index].fZ = zMod;
	      recPointList[index].fX = xMod;
	      
	      for(Int_t j = 0; j < index; j++)
		{
		  if(recPointList[j].fZ == zMod && recPointList[j].fX == xMod)
		    recPointList[j].fZ = -1;
		}
	      index++;
	    }
	  
	}
      else if(gain == 1 && fEnergyArray[xMod][zMod] == 0) 
	{
	  fEnergyArray[xMod][zMod] = fLowGainFactor * energyCount;
	   if(fEnergyArray[xMod][zMod] > fClusterThreshold)
	    {	      
	      recPointList[index].fZ = zMod;
	      recPointList[index].fX = xMod;
	      recPointList[index].fModule = cellData->fModuleID;
	      index++;
	    }
	}

    }  
  
  fPHOSModule = cellData->fModuleID;

  return index;
}//end BuildCellEnergyArray



/**
* Creating an array of rec points
* @param recPointStructArrayPtr array to store the rec points
* @param list list of rec points
* @param nPoints number of points
**/
Int_t 
AliHLTPHOSClusterizer::CreateRecPointStructArray(AliHLTPHOSRecPointDataStruct* recPointStructArrayPtr, 
						 AliHLTPHOSRecPointListDataStruct* list, 
						 Int_t nPoints) 

{
  //Create the rec point struct array

  Int_t flag = 0;
  Int_t edgeFlagRows = 0;
  Int_t edgeFlagCols = 0;
  Int_t k = 0;
  Int_t nRecPoints = 0;
  Int_t z = 0;
  Int_t x = 0;

  Float_t* energiesList = NULL;
  
  for(Int_t point = 0; point < nPoints; point++)
    {
      z = list[point].fZ;
      x = list[point].fX;
      if(z == -1) continue;

      if((z-fArraySize/2) < 0/*= - N_ROWS_MOD/2*/ || (z+fArraySize/2) >= N_ROWS_MOD/*/2*/) 
	{
	  edgeFlagRows = 1;
	  continue;
	}
      if((x-fArraySize/2) < 0/*= - N_COLUMNS_MOD/2*/ || (x+fArraySize/2) >= N_COLUMNS_MOD) 
	{
	  edgeFlagCols = 1;
	  continue;
	}
      

      if(!flag) energiesList = new Float_t[fMultiplicity];
      flag = 0;
      k = 0;	      
	  
      for(Int_t i = -fArraySize/2; i <= fArraySize/2; i++)
	{
	  if(flag) break;
	  for(Int_t j = -fArraySize/2; j <= fArraySize/2; j++)
	    {			  
	      
	      if(fEnergyArray[x+i][z+j] > fEnergyArray[x][z] && abs(i) < 2 && abs(j) < 2)
		{
		  flag = 1;
		  break;
		}
	      energiesList[k] = fEnergyArray[x+i][z+j];
	      k++;
	    }
	}
	  
      
      if(!flag && k)
	{
	  recPointStructArrayPtr[nRecPoints].fEnergiesListPtr = energiesList;
	  recPointStructArrayPtr[nRecPoints].fCoordinatesPtr[0] = x;
	  recPointStructArrayPtr[nRecPoints].fCoordinatesPtr[1] = z;
	  recPointStructArrayPtr[nRecPoints].fX = x;
	  recPointStructArrayPtr[nRecPoints].fZ = z;
	  recPointStructArrayPtr[nRecPoints].fMultiplicity = fMultiplicity;
	  recPointStructArrayPtr[nRecPoints].fPHOSModule = list[point].fModule;
	  nRecPoints++;
	}
    }



  energiesList = NULL;

  return nRecPoints;
  
}//end CreateRecPointStructArray


/**
* Calculating the center of gravity of a rec point
* Not working well at this point!
* @param recPointPtr pointer to the rec point
**/
Int_t
AliHLTPHOSClusterizer::CalculateCenterOfGravity(AliHLTPHOSRecPointDataStruct* recPointPtr)
{
  //Calculate the center of gravity

  Float_t xt = 0;
  Float_t zt = 0;
  Float_t xi = 0;
  Float_t zj = 0;
  Float_t w = 0;
  Float_t w0 = 4.5;
  Float_t wtot = 0;

  Int_t x = recPointPtr->fCoordinatesPtr[0];
  Int_t z = recPointPtr->fCoordinatesPtr[1];
  
  for(Int_t i = -fArraySize/2; i <= fArraySize/2; i++)
    {
      for(Int_t j = -fArraySize/2; j <= fArraySize/2; j++)
	{			  
	  xi = x + i;
	  zj = z + j;
	  w = TMath::Max( (Float_t)0., (Float_t)(w0 + log( fEnergyArray[x+i][z+j] / fEnergyArray[x][z]) ) ) ;
	  xt += xi * w ;
	  zt += zj * w ;
	  wtot += w ;
	  
	}
    }
  xt /= wtot ;
  zt /= wtot ;
  
  recPointPtr->fX = xt;
  recPointPtr->fZ = zt;
  
  return 0;

}//end CalculateCenterOfGravity
  


/**
* Create a simpler data struct for a rec point, containing
* only coordinates, energy and timing
* @param recPointPtr pointer to the rec point
* @param clusterStructPtr pointer to the cluster struct
**/
Int_t
AliHLTPHOSClusterizer::ClusterizeStruct(AliHLTPHOSRecPointDataStruct* recPointPtr, AliHLTPHOSClusterDataStruct* clusterStructPtr)
{
  //Simplify the rec points

  Float_t clusterEnergy = 0;
  Float_t* energiesListPtr = recPointPtr->fEnergiesListPtr;
  
  for(Int_t i = 0; i < recPointPtr->fMultiplicity; i++)
    {
      clusterEnergy += energiesListPtr[i];
    }

  clusterStructPtr->fLocalPositionPtr[0] = recPointPtr->fX;
  clusterStructPtr->fLocalPositionPtr[1] = recPointPtr->fZ;
  clusterStructPtr->fClusterEnergy = clusterEnergy;
  clusterStructPtr->fPHOSModule = recPointPtr->fPHOSModule;

  return 0;

}//end ClusterizeStruct




/**
* Resets the cell energy array
**/ 
Int_t 
AliHLTPHOSClusterizer::ResetCellEnergyArray()

{
  //Reset the cell energy array

  for(Int_t x = 0; x < N_ROWS_MOD; x++)
    {
      for(Int_t z = 0; z < N_COLUMNS_MOD; z++)
	{
	  fEnergyArray[x][z] = 0;
	}
    }

  return 0;

}//end ResetCellEnergyArray

