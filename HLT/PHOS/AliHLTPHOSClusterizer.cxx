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
#include "AliHLTPHOSConstants.h"

using namespace PhosHLTConst;


ClassImp(AliHLTPHOSClusterizer);

/**
* Main constructor
**/
AliHLTPHOSClusterizer::AliHLTPHOSClusterizer():fPHOSModule(0), fThreshold(0), fClusterThreshold(0), 
					       fHighGainFactor(0.005), fLowGainFactor(0.08),
					       fArraySize(3), fMultiplicity(fArraySize*fArraySize)
{
  //See header file for documentation
  
}//end

AliHLTPHOSClusterizer::AliHLTPHOSClusterizer(const AliHLTPHOSClusterizer &):fPHOSModule(0), fThreshold(0), fClusterThreshold(0), 
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

int
AliHLTPHOSClusterizer::BuildCellEnergyArray(AliHLTPHOSRcuCellEnergyDataStruct* cellData, 
					    AliHLTPHOSRecPointListDataStruct* recPointList)
{
  //BuildCellEnergyArray

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
int 
AliHLTPHOSClusterizer::CreateRecPointStructArray(AliHLTPHOSRecPointDataStruct* recPointStructArrayPtr, 
						 AliHLTPHOSRecPointListDataStruct* list, 
						 int nPoints) 

{
  //CreateRecPointStructArray

  Int_t flag = 0;
  Int_t edgeFlagRows = 0;
  Int_t edgeFlagCols = 0;
  Int_t k = 0;
  Int_t nRecPoints = 0;
  Int_t z = 0;
  Int_t x = 0;
  Int_t i = 0;
  Int_t j = 0;
  Int_t a = 0;

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
	  
      for(i = -fArraySize/2; i <= fArraySize/2; i++)
	{
	  if(flag) break;
	  for(j = -fArraySize/2; j <= fArraySize/2; j++)
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
	  for(a = 0; a < k; a++) 
	    {
	      (recPointStructArrayPtr[nRecPoints].fEnergiesListPtr)[a] = energiesList[a];
	    }
	  recPointStructArrayPtr[nRecPoints].fCoordinatesPtr[0] = x;
	  recPointStructArrayPtr[nRecPoints].fCoordinatesPtr[1] = z;
	  recPointStructArrayPtr[nRecPoints].fX = x;
	  recPointStructArrayPtr[nRecPoints].fZ = z;
	  // recPointStructArrayPtr[nRecPoints].fMultiplicity = fMultiplicity;
	  recPointStructArrayPtr[nRecPoints].fPHOSModule = list[point].fModule;
	  nRecPoints++;
	}
    }

  if(energiesList)
    {
      delete [] energiesList;
      energiesList = NULL;
    }

  return nRecPoints;
  
}//end CreateRecPointStructArray


/**
* Calculating the center of gravity of a rec point
* Not working well at this point!
* @param recPointPtr pointer to the rec point
**/
int
AliHLTPHOSClusterizer::CalculateCenterOfGravity(AliHLTPHOSRecPointDataStruct* recPointPtr)
{
  //CalculateCenterOfGravity
  //Copied from offline code, and modified

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
	  //printf("log in cog: %f\n", log( fEnergyArray[x+i][z+j] / fEnergyArray[x][z]));
	  xt += xi * w ;
	  zt += zj * w ;
	  wtot += w ;	  
	}
    }
  xt /= wtot ;
  zt /= wtot ;
  
  //printf("wtot cog: %f\n", wtot);

  recPointPtr->fX = xt;
  recPointPtr->fZ = zt;
  
  return 0;

}//end CalculateCenterOfGravity
  
int
AliHLTPHOSClusterizer::CalculateMoments(AliHLTPHOSRecPointDataStruct* recPointPtr, Bool_t axisOnly)
{
  //Copied from offline code, and modified
  
  
  // Calculate the shower moments in the eigen reference system
  // M2x, M2z, M3x, M4z
  // Calculate the angle between the shower position vector and the eigen vector
  
  Double_t wtot = 0. ;
  Double_t x    = 0.;
  Double_t z    = 0.;
  Double_t dxx  = 0.;
  Double_t dzz  = 0.;
  Double_t dxz  = 0.;
  Double_t lambda0=0, lambda1=0;
  Float_t logWeight = 4.5;
  
  //Int_t iDigit;
  
  // 1) Find covariance matrix elements:
  //    || dxx dxz ||
  //    || dxz dzz ||
  
  Int_t xi;
  Int_t zj;
  Double_t w;
  
  Int_t xc = recPointPtr->fCoordinatesPtr[0];
  Int_t zc = recPointPtr->fCoordinatesPtr[1];
  
  for(Int_t i = -fArraySize/2; i <= fArraySize/2; i++)
    {
      xi = xc + i;
      for(Int_t j = -fArraySize/2; j <= fArraySize/2; j++)
	{
	  zj = zc + j;
	  if (fEnergyArray[xi][zj] > 0) 
	    {
	      w     = TMath::Max(0.,logWeight + log( fEnergyArray[xi][zj]/fEnergyArray[xc][zc] ) ) ;
	      //printf("log in mom: %f\n", TMath::Log( fEnergyArray[xi][zj] / fEnergyArray[xc][zc]));
	      x    += w * xi ;
	      z    += w * zj ; 
	      dxx  += w * xi * xi ;
	      dzz  += w * zj * zj ;
	      dxz  += w * xi * zj ; 
	      wtot += w ;
	    }
	}
    }
  //printf("wtot: %f\n", wtot);
  if (wtot>0) 
    {
      x   /= wtot ;
      z   /= wtot ;
      dxx /= wtot ;
      dzz /= wtot ;
      dxz /= wtot ;
      dxx -= x * x ;
      dzz -= z * z ;
      dxz -= x * z ;
      
  // 2) Find covariance matrix eigen values lambda0 and lambda1

      lambda0 =  0.5 * (dxx + dzz) + TMath::Sqrt( 0.25 * (dxx - dzz) * (dxx - dzz) + dxz * dxz )  ;
      lambda1 =  0.5 * (dxx + dzz) - TMath::Sqrt( 0.25 * (dxx - dzz) * (dxx - dzz) + dxz * dxz )  ;
    }
    
  recPointPtr->fM2x   = lambda0;
  recPointPtr->fM2z   = lambda1;
  //  printf("lambda0 = %f -- lambda1 = %f\n", lambda0, lambda1);
  // 3) Find covariance matrix eigen vectors e0 and e1
  if(!axisOnly)
    {
      TVector2 e0,e1;
      if (dxz != 0)
	e0.Set(1.,(lambda0-dxx)/dxz);
      else
	e0.Set(0.,1.);
      
      e0 = e0.Unit();
      e1.Set(-e0.Y(),e0.X());
      
      // 4) Rotate cluster tensor from (x,z) to (e0,e1) system
      //    and calculate moments M3x and M4z
      
      Float_t cosPhi = e0.X();
      Float_t sinPhi = e0.Y();
      
      Float_t xiPHOS ;
      Float_t zjPHOS ;
      Double_t dx3, dz3, dz4;
      wtot = 0.;
      x    = 0.;
      z    = 0.;
      dxx  = 0.;
      dzz  = 0.;
      dxz  = 0.;
      dx3  = 0.;
      dz3  = 0.;
      dz4  = 0.;
      for(Int_t i = -fArraySize/2; i <= fArraySize/2; i++)
	{
	  for(Int_t j = -fArraySize/2; j <= fArraySize/2; j++)
	    {
	      xi = xc + i;
	      zj = zc + j;
	      xiPHOS = (Float_t)xi*cosPhi + (Float_t)zj*sinPhi;
	      zjPHOS = (Float_t)zj*cosPhi - (Float_t)xi*sinPhi;
	      //	      xi = (Float_t)xi*cosPhi + (Float_t)zj*sinPhi;
	      //	      zj = (Float_t)zj*cosPhi - (Float_t)xi*sinPhi;
	      if (fEnergyArray[xi][zj] > 0) 
		{
		  w     = TMath::Max(0.,logWeight+TMath::Log(fEnergyArray[xi][zj]/fEnergyArray[xc][zc] ) ) ;
		  //  printf("log in mom: %f\n", TMath::Log( fEnergyArray[xi][zj] / fEnergyArray[xc][zc]));
			  
		  x    += w * xi ;
		  z    += w * zj ; 
		  dxx  += w * xi * xi ;
		  dzz  += w * zj * zj ;
		  dxz  += w * xi * zj ; 
		  dx3  += w * xi * xi * xi;
		  dz3  += w * zj * zj * zj ;
		  dz4  += w * zj * zj * zj * zj ;
		  wtot += w ;
		}
	    }
	}
      if (wtot>0) 
	{
	  x   /= wtot ;
	  z   /= wtot ;
	  dxx /= wtot ;
	  dzz /= wtot ;
	  dxz /= wtot ;
	  dx3 /= wtot ;
	  dz3 /= wtot ;
	  dz4 /= wtot ;
	  dx3 += -3*dxx*x + 2*x*x*x;
	  dz4 += -4*dz3*z + 6*dzz*z*z -3*z*z*z*z;
	  dxx -= x * x ;
	  dzz -= z * z ;
	  dxz -= x * z ;
	}
  
      // 5) Find an angle between cluster center vector and eigen vector e0
      
      Float_t phi = 0;
      if ( (x*x+z*z) > 0 ) 
	phi = TMath::ACos ((x*e0.X() + z*e0.Y()) / sqrt(x*x + z*z));

      recPointPtr->fM3x   = dx3;
      recPointPtr->fM4z   = dz4;
      recPointPtr->fPhixe = phi;
      //      printf("%f\n",phi);
    }
  return 0;
}


/**
* Create a simpler data struct for a rec point, containing
* only coordinates, energy and timing
* @param recPointPtr pointer to the rec point
* @clusterStructPtr pointer to the cluster struct
**/
int
AliHLTPHOSClusterizer::ClusterizeStruct(AliHLTPHOSRecPointDataStruct* recPointPtr, AliHLTPHOSClusterDataStruct* clusterStructPtr)
{
  //ClusterizeStruct

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
int 
AliHLTPHOSClusterizer::ResetCellEnergyArray()

{
  //ResetCellEnergyArray

  for(Int_t x = 0; x < N_ROWS_MOD; x++)
    {
      for(Int_t z = 0; z < N_COLUMNS_MOD; z++)
	{
	  fEnergyArray[x][z] = 0;
	}
    }

  return 0;

}//end ResetCellEnergyArray

