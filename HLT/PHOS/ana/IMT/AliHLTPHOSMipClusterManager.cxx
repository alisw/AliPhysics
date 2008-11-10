// $Id$

/**************************************************************************
 * This file is property of and copyright by the Experimental Nuclear     *
 * Physics Group, Dep. of Physics                                         *
 * University of Oslo, Norway, 2007                                       *
 *                                                                        *
 * Author: Per Thomas Hille <perthi@fys.uio.no> for the ALICE HLT Project.*
 * Contributors are mentioned in the code where appropriate.              *
 * Please report bugs to perthi@fys.uio.no                                *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
#include "AliHLTPHOSMipClusterManager.h"
#include "AliHLTPHOSMipCluster.h"
#include "AliHLTPHOSDigit.h"


AliHLTPHOSMipClusterManager::AliHLTPHOSMipClusterManager()  
{
  
  for(int z=0; z<N_ZROWS_MOD; z++)
    {
      for(int x =0; x <N_XCOLUMNS_MOD; x++)
	{
	  fMipClusterFarm[z][x] = new  AliHLTPHOSMipCluster();
  	  fMipClusterTable[z][x] = 0;
	  //	  fMipClusterFarm[z][x]->PrintInfo();
	}
    }
}



AliHLTPHOSMipClusterManager::~AliHLTPHOSMipClusterManager()
{

}


void
AliHLTPHOSMipClusterManager::AddDigit(AliHLTPHOSDigit *inDigit)
{
  int tmpZ = inDigit->GetZ();
  int tmpX = inDigit->GetX();

  //  cout << "searching for neighbour ar Z =" << tmpZ << "  tmpX ="<< tmpX <<endl;

  AliHLTPHOSMipCluster *mipCluster = GetNeighbour(Z_CENTER_RANGE, X_CENTER_RANGE, tmpZ, tmpX);
		 
  if(mipCluster != 0)
    {
      if(mipCluster->GetCenterAmplitude() < inDigit->GetAmplitude())
	{
	  int oldZ = mipCluster->GetZ();
	  int oldX = mipCluster->GetX();
	  AliHLTPHOSMipCluster *tmpCluster = MoveCluster(oldZ, oldX, tmpZ, tmpX);
	  tmpCluster->AddDigit(inDigit);
	}
      else
	{
	  //	  cout <<"adding digit to exsising cluster"<< endl;	 	
	  mipCluster->AddDigit(inDigit);
	}
    }
		 
  else
    {
      //      cout <<"creating new  cluster"<< endl;  
      NewMipTableEntry(inDigit, tmpZ, tmpX);
      

      //     cout << "fMipCnt   " << fMipCnt  <<endl;
      fMipCnt ++;  
    }
}


void
AliHLTPHOSMipClusterManager::NewMipTableEntry(AliHLTPHOSDigit *digit, int z, int x)
{
  if(digit == 0)
    {
      cout << "AliHLTPHOSMipClusterManager::NewMipTableEntry, Error, digit = NULL" << endl;
    }

  else
    {
      //      fMipClusterTable[z][x] =  fMipClusterFarm[fMipCnt];
      fMipClusterTable[z][x] =  fMipClusterFarm[z][x];
      fMipClusterTable[z][x]->SetCenterCoordinate(z,  x); 
      fMipClusterFarm[z][x]->AddDigit(digit); 
    }

}


AliHLTPHOSMipCluster* 
AliHLTPHOSMipClusterManager::MoveCluster(int zFrom, int xFrom, int zTo, int xTo )
{
  AliHLTPHOSMipCluster *swapTo   = fMipClusterFarm[zFrom][xFrom]; 
  AliHLTPHOSMipCluster *swapFrom = fMipClusterFarm[zTo][xTo];   

  //  fMipClusterFarm[zTo][xTo] = swapFrom; 
  // fMipClusterFarm[zFrom][xFrom] = swapTo; 

  fMipClusterFarm[zTo][xTo] = swapTo; 
  fMipClusterFarm[zFrom][xFrom] = swapFrom; 

  fMipClusterTable[zFrom][xFrom] = 0; 
  fMipClusterTable[zTo][xTo] = fMipClusterFarm[zTo][xTo];

  fMipClusterFarm[zFrom][xFrom]->SetCenterCoordinate(zFrom, xFrom);
  fMipClusterFarm[zTo][xTo]->SetCenterCoordinate(zTo, xTo);


  fMipClusterTable[zFrom][xFrom] = 0; 
  fMipClusterTable[zTo][xTo] = fMipClusterFarm[zTo][xTo];

  //  fMipClusterFarm[zFrom][xFrom]->Remap();
  fMipClusterFarm[zFrom][xFrom]->ResetMipCluster();

  fMipClusterFarm[zTo][xTo]->Remap();

  return fMipClusterTable[zTo][xTo];
}


AliHLTPHOSMipCluster*
AliHLTPHOSMipClusterManager::GetNeighbour(int zRange, int xRange, int zIn, int xIn)
{
  AliHLTPHOSMipCluster* tmp = 0;

  if(fMipClusterTable[zIn][xIn] != 0)
    {
      //      cout  << "ERROR; table[" <<zIn << "]["<< xIn << "]" << "already filled for this event" << endl;
    }
  else
    {
      int tmpZmin = 0;
      int tmpZmax = 0;    
      int tmpXmin = 0;
      int tmpXmax = 0;        

      
      if(zIn -zRange < 0)
	{
	  tmpZmin = 0;
	}
      else
	{
	  tmpZmin = zIn - zRange;
	}
     
      if(zIn + zRange > N_ZROWS_MOD)
       {
	 tmpZmax = N_ZROWS_MOD ;
       }
      else
	{
	  tmpZmax = zIn + zRange;
	}


       if(xIn - xRange < 0)
	{
	  tmpXmin = 0;
	}
      else
	{
	  tmpXmin = xIn - xRange;
	}
     
      if(xIn + xRange > N_XCOLUMNS_MOD)
       {
	 tmpXmax = N_XCOLUMNS_MOD ;
       } 
      else
	{
	  tmpXmax = xIn + xRange;
	}
      
      for(int z = tmpZmin; z < tmpZmax; z ++)
	{
	  for(int x = tmpXmin; x < tmpXmax; x ++) 
	    {
	      if( fMipClusterTable[z][x] != 0 )
		{
		  tmp = fMipClusterTable[z][x];
		}
	    }
	}
    }
  
  return tmp;
  
}


void 
AliHLTPHOSMipClusterManager::AddToMipTable(AliHLTPHOSDigit *digit)
{

  int z = digit->GetZ();
  int x = digit->GetX();
  
  if(fMipClusterTable[z][x] == 0)
    {
      fMipClusterTable[z][x] = fMipClusterFarm[z][x];
    }
}


void
AliHLTPHOSMipClusterManager::PrintDigitInfo(AliHLTPHOSDigit *digit)
{
  Int_t *data = digit->GetRawData();
  //  cout << endl;
  cout <<"Gain = " << digit->GetGain()  <<endl;
  cout <<"Z = " << digit->GetZ()  <<endl;
  cout <<"X = " << digit->GetX()  <<endl;
  cout <<"amplitude = " << digit->GetAmplitude()  <<endl; 
  cout <<"crazynes =  " << digit->GetCrazyness()  <<endl;
  
  for(int i=0; i< 70; i++)
    {
      cout  << data[i]<<"\t"; 
    }
   
  void CopyDigit(AliHLTPHOSDigit *inDigit, AliHLTPHOSDigit *copiedDigit);
  cout << endl;
}


void
AliHLTPHOSMipClusterManager::ResetMipManager()
{
  ResetClusterFarm();

  for(int z=0; z<N_ZROWS_MOD; z++)
    {
      for(int x =0; x <N_XCOLUMNS_MOD; x++)
	{
	  //	  if(fMipClusterTable[z][x] != 0)
	  //	    {
	  //	      ResetClusterFarm();
	  //	    }
	  
	  fMipClusterTable[z][x] = 0;
	}
    }
}


void
AliHLTPHOSMipClusterManager::ResetClusterFarm()
{

  for(int z=0;  z < N_ZROWS_MOD; z++ )
    {
      for(int x=0; x < N_XCOLUMNS_MOD; x ++)
	{
	  fMipClusterFarm[z][x]->ResetMipCluster();
	}
    }
}



void
AliHLTPHOSMipClusterManager::PrintClusters(int minEntries)
{
  for(int z=0; z<N_ZROWS_MOD; z++)
    {
      for(int x =0; x <N_XCOLUMNS_MOD; x++)
	{
	  if(fMipClusterTable[z][x] != 0)
	    {
	      if( fMipClusterTable[z][x]->GetEntries() >= minEntries  )
	      {
		cout << endl;
		cout << "AliHLTPHOSIsolatedMipTrigger::PrintClusters, printing info for"  << endl;
		cout << "mipClusterTable[" << z << "][" << x << "]";
		fMipClusterTable[z][x]->PrintInfo();
	      }
	    }
	}
    }
}


AliHLTPHOSMipCluster*
AliHLTPHOSMipClusterManager::GetCluster(int z, int x)
{
  return fMipClusterTable[z][x];
}


