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
#include "AliHLTPHOSMipCluster.h"
#include "AliHLTPHOSDigit.h"

AliHLTPHOSMipCluster::AliHLTPHOSMipCluster() : fZ(0),
					       fX(0),
					       fEntries(0)
{
  for(int z=0; z <Z_RANGE; z++)
    {
      for(int x=0;  x < X_RANGE; x++)
	{
	  fMipcluster[z][x]= new AliHLTPHOSDigit();
	}
    }
}


AliHLTPHOSMipCluster::~AliHLTPHOSMipCluster()
{

}


int             
AliHLTPHOSMipCluster::GetEntries()
{
  return fEntries;
}


void
AliHLTPHOSMipCluster::ResetMipCluster()
{
   for(int i=0; i<5; i++)
    {
      for(int j=0; j<5; j++)
	{
	  fMipcluster[i][j]->SetX(-1);
	  fMipcluster[i][j]->SetZ(-1);
	  fMipcluster[i][j]->SetAmplitude(-1); 
	  fMipcluster[i][j]->SetCrazyness(-1);
	  //	  fMipcluster[i][j]->SetZ(0);
	  fMipcluster[i][j]->SetEnergy(-1);
	}

    } 
   fEntries = 0;
   fZ = 0;
   fX = 0;
}


void
AliHLTPHOSMipCluster::AddDigit(AliHLTPHOSDigit *digit)
{
  int relZ = digit->GetZ() -fZ +2; //Z corrdinate relative to the 5x5 clustertable
  int relX = digit->GetX() -fX +2; //X corrdinate relative to the 5x5 clustertable
  CopyDigit(digit, fMipcluster[relZ][relX]);
  fEntries  ++;
}


void 
AliHLTPHOSMipCluster::CopyDigit(AliHLTPHOSDigit *inDigit, AliHLTPHOSDigit *copiedDigit)
{
  if(inDigit == 0)
    {
      cout << "AliHLTPHOSMipCluster::CopyDigit, Fatal ERROR digit == NULL" << endl;
    }

  else if(copiedDigit == 0)
    {
      cout << "AliHLTPHOSMipCluster::CopyDigit, Fatal ERROR copydigit == NULL" << endl;
    } 

  else
    {
      copiedDigit->SetZ(inDigit->GetZ());
      copiedDigit->SetX(inDigit->GetX());
      copiedDigit->SetAmplitude(inDigit->GetAmplitude());
      copiedDigit->SetCrazyness(inDigit->GetCrazyness());
      copiedDigit->SetEnergy(inDigit->GetEnergy());
      copiedDigit->SetGain(inDigit->GetGain());
      copiedDigit->SetTime(inDigit->GetTime());
      copiedDigit->SetRawData(inDigit->GetRawData());
    }
}


void
AliHLTPHOSMipCluster::Remap()
{
  AliHLTPHOSDigit *tmpmipcluster[Z_RANGE][X_RANGE]; 

  for(int z=0; z < Z_RANGE; z++ )
    {
      for(int x=0; x < X_RANGE; x++ )
	{
	  tmpmipcluster[z][x] = fMipcluster[z][x];  
	}   
    }
  
  for(int z=0; z < Z_RANGE; z++ )
    {
      for(int x=0; x < X_RANGE; x++ )
	{
	  int relZ = tmpmipcluster[z][x]->GetZ() -fZ +2;
	  int relX = tmpmipcluster[z][x]->GetX() -fX +2;

	  if(relZ < 0 || relX <0)
	    {
	      if(  (tmpmipcluster[z][x]->GetZ() >= 0 )  &&  (tmpmipcluster[z][x]->GetX() >=0)  ) 
		{

		  /*
		  cout << endl;
		  cout << "INFO:  AliHLTPHOSMipCluster::Remap() (relZ, relX)           =  " <<"(" << relZ << ","<< relX << ")"<<endl;
		  cout << "    :  AliHLTPHOSMipCluster::Remap() (clusterZ, clusterX)  =  " <<"(" <<  tmpmipcluster[z][x]->GetZ()  << ","<<  tmpmipcluster[z][x]->GetX() << ")"<<endl;
		  cout << "    :  AliHLTPHOSMipCluster::Remap() (z, x)                 =  " <<"(" << z << ","<< x << ")"<<endl;
		  cout << "    :  AliHLTPHOSMipCluster::Remap() (fZ, fX)               =  " <<"(" << fZ << ","<< fX << ")"<<endl;
		  cout << endl;
		  */
		}
	    }


	  if( (relZ >= 0)   &&  (relZ < Z_RANGE) &&  (relX >= 0)  && (relX < X_RANGE)  &&    (z >= 0)  && (z < Z_RANGE) &&  (x >= 0)  && (x < X_RANGE))
	    {
	      if(  (relZ >= 0 )  &&   (relX >= 0)  )
		{
		  if((relZ < Z_RANGE)  && (relX < X_RANGE)  &&    (z >= 0)  && (z < Z_RANGE) &&  (x >= 0)  && (x < X_RANGE))
		    {
		   
		      int relZ = tmpmipcluster[z][x]->GetZ() -fZ +2;
		      int relX = tmpmipcluster[z][x]->GetX() -fX +2;   
		    
		      if( (relZ == 0) && (relX == 0) )
			{
			  cout <<"AliHLTPHOSMipCluster::Remap(), ERROR,(z,x) ="<< "(" <<relZ << "," << relX <<")" <<endl;
			}
		      

		      AliHLTPHOSDigit *swap = fMipcluster[relZ][relX];
		      fMipcluster[relZ][relX] = tmpmipcluster[z][x];
		      fMipcluster[z][x] = swap;
		    }  
		  
	
		  else
		    {
		      /*
		      cout << endl;
		      cout <<"AliHLTPHOSMipCluster::Remap(), fatal ERROR in coordinates" << endl;
		      cout << "AliHLTPHOSMipCluster::Remap(),  fZ ="<<  fZ  << endl;
		      cout << "AliHLTPHOSMipCluster::Remap(),  fX ="<<  fX  << endl;	  
		      cout << "AliHLTPHOSMipCluster::Remap(),  tmpmipcluster[z][x]->GetZ() ="<<  tmpmipcluster[z][x]->GetZ()  << endl;
		      cout << "AliHLTPHOSMipCluster::Remap(),  tmpmipcluster[z][x]->GetX() ="<<  tmpmipcluster[z][x]->GetX()  << endl;	  
		      cout << "AliHLTPHOSMipCluster::Remap(), relZ =  "  << relZ << endl;
		      cout << "AliHLTPHOSMipCluster::Remap(), relX =  "  << relX << endl;	
		      cout << endl;
		      */
		    }

		  
		}
	    }
	  
	}   
    }
  //  cout << "AliHLTPHOSMipCluster::Remap(), clusters after remap ="  << endl;
  //  PrintInfo();
}

void
AliHLTPHOSMipCluster::ResetDigit(int z, int x)
{
  //  fMipcluster[z][x];
}


void
AliHLTPHOSMipCluster::PrintInfo()
{
  cout << endl;
  cout << "Printing amplitudes for cluster " << endl;
  cout << "Z Center = " << fZ << endl;
  cout << "X Center = " << fX << endl;

  cout << "3x3 sum =" << Get3x3Sum() << endl;
   
  for(int z = 0; z< Z_RANGE; z ++)
    {
      for(int x = 0;  x < X_RANGE;  x ++)  
	{
	  cout  << fMipcluster[z][x]  <<"\t";;
	}

      cout << endl;
    } 
  
  
  for(int z = 0; z< Z_RANGE; z ++)
    {
      for(int x = 0;  x < X_RANGE;  x ++)  
	{
	  cout  << fMipcluster[z][x]->GetAmplitude()  <<"\t";;
	}

      cout << endl;
    }

  //  cout << endl;

  for(int z = 0; z< Z_RANGE; z ++)
    {
      for(int x = 0;  x < X_RANGE;  x ++)  
	{
	  
	  cout  << "("<< fMipcluster[z][x]->GetZ()  <<","<<fMipcluster[z][x]->GetX() << ")" <<"\t";
	}
      cout << endl;
    }

    cout << endl;
}


void
AliHLTPHOSMipCluster:: SetCenterCoordinate(int z, int x)
{
  fZ = z;
  fX = x;
}


Float_t
AliHLTPHOSMipCluster::GetCenterAmplitude()
{
  return  fMipcluster[2][2]->GetAmplitude();
}


Int_t
AliHLTPHOSMipCluster::GetZ()
{
  return fZ;
}


Int_t
AliHLTPHOSMipCluster::GetX()
{
  return fX;
}


Float_t
AliHLTPHOSMipCluster::Get3x3Sum()
{
  Float_t tmpSum = 0;
  
  Float_t tmpAmp = 0;

  for(int z =1; z <= 3; z++)
    {
      for(int x =1; x <= 3; x++)
	{
	  //	  cout <<  "AliHLTPHOSMipCluster::Get3x3Sum(), z =" << z << "x=" << x << endl;

	  tmpAmp = fMipcluster[z][x]->GetAmplitude(); 

	  if(tmpAmp > 0)
	    {
	      tmpSum += tmpAmp;
	    }
	  
	} 
      
    } 

  return tmpSum;
}
