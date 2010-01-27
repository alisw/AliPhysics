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

// The Peak-Finder algorithm
// The amplitude is extracted  as a
// weighted sum of the samples using the 
// best possible weights.
// The wights is calculated only once and the
// Actual extraction of amplitude and peak position
// Is done with a simple vector multiplication, allowing for
// Extreemely fast computations. 

#include "AliCaloRawAnalyzerPeakFinder.h"
#include "AliCaloBunchInfo.h"
#include "AliCaloFitResults.h"
#include <iostream>
#include "unistd.h"
#include "TMath.h"
#include "AliLog.h"

using namespace std;

AliCaloRawAnalyzerPeakFinder::AliCaloRawAnalyzerPeakFinder() :AliCaloRawAnalyzer(), 
								fTof(0), 
								fAmp(0)
{
  //comment

  fNsampleCut = 5;

  for(int i=0; i < MAXSTART; i++)
    {
      for(int j=0; j < SAMPLERANGE; j++ )
	{
	  fPFAmpVectors[i][j] = new double[100];
	  fPFTofVectors[i][j] = new double[100];
	  
	  for(int k=0; k < 100; k++ )
	    {
	      fPFAmpVectors[i][j][k] = 0; 
	      fPFTofVectors[i][j][k] = 0;
	    }
	}
    }
  LoadVectors();
}


AliCaloRawAnalyzerPeakFinder::~AliCaloRawAnalyzerPeakFinder()
{
  //comment
  for(int i=0; i < MAXSTART; i++)
    {
      for(int j=0; j < SAMPLERANGE; j++ )
	{
	  delete[] fPFAmpVectors[i][j];
	  delete[] fPFTofVectors[i][j];
	}
    }
}


AliCaloFitResults 
AliCaloRawAnalyzerPeakFinder::Evaluate( const vector<AliCaloBunchInfo> &bunchvector, const UInt_t altrocfg1,  const UInt_t altrocfg2 )
{
  // Extracting the amplitude using the Peak-Finder algorithm
  // The amplitude is a weighted sum of the samples using 
  // optimum weights.

  short maxampindex; //index of maximum amplitude
  short maxamp; //Maximum amplitude
  fAmp = 0;
  fAmpA[0] = 0;
  fAmpA[1] = 0;
  fAmpA[2] = 0;

  int index = SelectBunch( bunchvector,  &maxampindex,  &maxamp );
 
  if( index >= 0)
    {
      Float_t ped = ReverseAndSubtractPed( &(bunchvector.at(index))  ,  altrocfg1, altrocfg2, fReversed  );
      Float_t maxf = TMath::MaxElement(   bunchvector.at(index).GetLength(),  fReversed );
      
      if(  maxf < fAmpCut  ||  ( maxamp - ped) > 900  )	 // (maxamp - ped) > 900 = Close to saturation (use low gain then)
	{
	  //	  cout << __FILE__ << __LINE__ <<":, maxamp = " << maxamp << ", ped = "<< ped  << ",. maxf = "<< maxf << ", maxampindex = "<< maxampindex  << endl;
	  return  AliCaloFitResults( maxamp, ped,  -1, maxf,   maxampindex, -1, -1 );
 	}
      
      int first;
      int last;
      
      if ( maxf > fAmpCut )
	{
	  SelectSubarray( fReversed,  bunchvector.at(index).GetLength(),  maxampindex -  bunchvector.at(index).GetStartBin(), &first, &last);
	  int nsamples =  last - first;
	  if( ( nsamples  )  >= fNsampleCut )
	    {
	      int startbin = bunchvector.at(index).GetStartBin();  
	      int n = last -first;  
	      int pfindex = n - fNsampleCut; 
	      pfindex = pfindex > SAMPLERANGE ? SAMPLERANGE : pfindex;

	      int dt =  maxampindex - startbin -2; 
	      for(int i=0; i < SAMPLERANGE; i++ )
		{
		  //	  int dt =  maxampindex - startbin -2;
		  
		  //		  double tmp[3];
 
		  //	  tmp[0]  =  fReversed[ dt  +i -1]; 
		  //		  tmp[1]  =  fReversed[ dt  +i];  
		  //		  tmp[2]  =  fReversed[ dt  +i +1]; 
		  
		  for(int j = 0; j < 3; j++)
		    {
		      //    fAmpA[j] += fPFAmpVectors[0][pfindex][i]*tmp[j]; 
		      fAmpA[j] += fPFAmpVectors[0][pfindex][i]*fReversed[ dt  +i +j -1];

		    }
		}
	      
	      double diff = 9999;

	      int tmpindex = 0;

	      for(int k=0; k < 3; k ++)
		{
		  //		  cout << __FILE__ << __LINE__ << "amp[="<< k <<"] = " << fAmpA[k] << endl;
		  if(  TMath::Abs(fAmpA[k] - ( maxamp - ped) )  < diff)
		    {
		      diff = TMath::Abs(fAmpA[k] - ( maxamp - ped));
		      tmpindex = k; 
		    }
		}
	      
	      double tof = 0;

	      for(int k=0; k < SAMPLERANGE; k++   )
		{
		  tof +=  fPFTofVectors[0][pfindex][k]*fReversed[ dt  +k + tmpindex -1 ];   
		}
	      
	      tof = tof /  fAmpA[tmpindex];

	      //   return AliCaloFitResults( maxamp, ped , -1, fAmp, -1, -1, -1 );  
	      return AliCaloFitResults( maxamp, ped , -1, fAmpA[tmpindex], tof, -2, -3 ); 
		  
	    }

	  else
	    {
	      return AliCaloFitResults( maxamp, ped , -5, maxf, -6, -7, -8 ); 
	    }
	}
    }
 
  //  cout << __FILE__ << __LINE__ <<  "WARNING, returning amp = -1 "  <<  endl;

  return  AliCaloFitResults(-1, -1);
}


void 
AliCaloRawAnalyzerPeakFinder::LoadVectors()
{
  //Read in the Peak finder vecors from file
  for(int i = 0;  i < MAXSTART ; i++)
    {
      for( int j=0; j < SAMPLERANGE; j++)
	{
	  char filename[256];
	  int n = j+fNsampleCut;
	  double start = (double)i+0.5;
	  sprintf(filename, "%s/EMCAL/vectors-emcal/start%.1fN%dtau0.235fs10dt1.0.txt", getenv("ALICE_ROOT"), start, n);
	  FILE *fp = fopen(filename, "r");
	  
	  if(fp == 0)
	    {
	      AliFatal( Form( "could not open file: %s", filename ) );
	    }
	  else
	    {
	      for(int m = 0; m < n ; m++ )
		{
		  fscanf(fp, "%lf\t", &fPFAmpVectors[i][j][m] );
		}
	      fscanf(fp, "\n");
	      
	      for(int m = 0; m < n ; m++ )
		{
		  fPFTofVectors[i][j][m] = 1;
	 //	  fscanf(fp, "%lf\t", &fPFTofVectors[i][j][m] );
		}
	      
	      fclose (fp);
	    }
	}
    }
}
