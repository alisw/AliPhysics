/**************************************************************************
 * This file is property of and copyright by the Relativistic Heavy Ion   *
 * Group (RHIG),  Department of Physics Yale University, US, 2010         *
 *                                                                        *
 * Author: Per Thomas Hille <perthomas.hille@yale.edu> for the ALICE EMCAL*
 * project. Contributors are mentioned in the code where appropriate.     *
 * Please report bugs to perthomas.hille@yale.edu                         *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                *
 **************************************************************************/

//Container class for Peak Finder vectors

#include "AliCaloPeakFinderVectors.h"
#include <iostream>


using namespace std;



ClassImp( AliCaloPeakFinderVectors)



AliCaloPeakFinderVectors::AliCaloPeakFinderVectors()
{
  ResetVectors();
}


AliCaloPeakFinderVectors::~AliCaloPeakFinderVectors()
{
  
}


void 
AliCaloPeakFinderVectors::ResetVectors()
{
  // As implied by function name
  for(int i=0; i < MAXSTART; i++ )
    {
      for(int j=0; j < SAMPLERANGE; j++)
	{
	  if(i < MAXSTART  && j < SAMPLERANGE )
	    {
	      for(int k = 0;  k < 100; k++)
		{
		  fPFAmpV[i][j][k]   =  0 ;
		  fPFTofV[i][j][k]   =  0 ;
		  fPFAmpVC[i][j][k]  =  0 ;
		  fPFTofVC[i][j][k]  =  0 ;
		}
	    }
	}	
    }
}


void 
AliCaloPeakFinderVectors::SetVector(const int i, const int j, const Double_t  *const a, const Double_t *const t,   
				    const Double_t *const ac, const Double_t *const tc )
{
  // As implied by function name
  if(i < MAXSTART  && j < SAMPLERANGE )
    {
      for(int k = 0;  k < 100; k++)
	{
	  fPFAmpV[i][j][k] =  a[k];
	  fPFTofV[i][j][k] =  t[k];
	  fPFAmpVC[i][j][k] = ac[k];
	  fPFTofVC[i][j][k] = tc[k];
	}
    }
}


void 
AliCaloPeakFinderVectors::GetVector(const int i, const int j, Double_t *const a, Double_t *const t,   
				    Double_t *const ac, Double_t *const tc ) const
{
  // As implied by function name
  if(i < MAXSTART  && j < SAMPLERANGE )
    {
      for( int k = 0;  k < 100; k++)
	{
	  a[k]  = fPFAmpV[i][j][k];
	  t[k]  = fPFTofV[i][j][k];
	  ac[k] = fPFAmpVC[i][j][k];
	  tc[k] = fPFTofVC[i][j][k];
	}
    }
}


void 
AliCaloPeakFinderVectors::PrintVectors() const
{
  // As implied by function name
  cout << __FILE__ << __LINE__ << __FUNCTION__ << endl;
  for(int i= 0; i < MAXSTART; i++ )
    {
      for(int j=0; j < SAMPLERANGE; j++ )
	{
	  for(int k=0; k < 10; k++ )
	    {
	      cout << fPFAmpV[i][j][k] << "\t";
	    }
	  cout << endl;
	}
      cout << endl; 
    }
}  

