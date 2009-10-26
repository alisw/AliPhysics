// $Id: AliHLTCALORawAnalyzerPeakFinder.cxx 29824 2008-11-10 13:43:55Z richterm $

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




#include "AliHLTCaloRawAnalyzerPeakFinder.h"
//#include <iostream>
#include <cmath>
#include "AliHLTCaloUtilities.h" 

using std::cout;
using std::endl;

//ClassImp(AliHLTCaloRawAnalyzerPeakFinder) 



/**
 * The AliHLTCaloPeakfinder class is the class for extracting the basic signal parameters
 * "timing" and "energy" from the Calorimeter's raw data. Physical data will for a given readout channel be
 * a sequense of ADC digitized 10 bit integer values, however for performance reasons all values used in
 * calculation is of type double.
 **/
AliHLTCaloRawAnalyzerPeakFinder::AliHLTCaloRawAnalyzerPeakFinder():AliHLTCaloRawAnalyzer(), 
                                                                   fTVectorPtr(0), 
								   fAVectorPtr(0), 
								   fTVectorSize(0), 
								   fAVectorSize(0)
//  fUtilitiesPtr(0)
{
  //  fUtilitiesPtr = new  AliHLTCaloUtilities(); 
  //  cout <<"PeakFinder:You cannot invoke the Fitter without arguments"<<endl;;
}


//___________________________________________________________________
AliHLTCaloRawAnalyzerPeakFinder::~AliHLTCaloRawAnalyzerPeakFinder()
{

} //end AliHLTCaloRawAnalyzerPeakFinder



void 
AliHLTCaloRawAnalyzerPeakFinder::SetTVector(Double_t *tVec, Int_t size)
{
  fTVectorSize = size;

  if(fTVectorPtr != 0)
    {
      delete fTVectorPtr;
    }
  
  fTVectorPtr = new Double_t[size];

  for(int i=0; i< size; i++)
    {
      fTVectorPtr[i] = tVec[i];
    }
}


//___________________________________________________________________
void
AliHLTCaloRawAnalyzerPeakFinder::SetAVector(Double_t *aVec, Int_t size)
{
    
  fAVectorSize = size;

  if(fAVectorPtr != 0)
    {
      delete fAVectorPtr;
    }
  
  fAVectorPtr = new Double_t[size];

  for(int i=0; i< size; i++)
    {
      fAVectorPtr[i] = aVec[i];
    }
}



//___________________________________________________________________
void 
AliHLTCaloRawAnalyzerPeakFinder::Evaluate(Int_t /*start*/, Int_t length)
{
  fDTof = 0;
  fDAmpl = 0;
  Int_t tmpLength;

  if(fTVectorPtr == 0 || fAVectorPtr == 0)
    {

    }
  else
    {
      if(length <  fTVectorSize)
	{
	  tmpLength = length;
	}
      else
	{
	  tmpLength = fTVectorSize;
	}
      
      for(int i=0; i < tmpLength; i++)
	{  
	  //fDAmpl += fAVectorPtr[i]*fIntDataPtr[i]; removed 18 april 2008    
	  fDAmpl += fAVectorPtr[i]*fShortDataPtr[i];   
	}

      for(int i=0; i < tmpLength; i++)
	{   
	  //  fDTof += fTVectorPtr[i]*fIntDataPtr[i];  removed 18 april 2008   
	  fDTof += fTVectorPtr[i]*fShortDataPtr[i]; 
	}
      
      if(fDAmpl > 900)
	{
	  //	  Double_t tmpMax = MaxValue(const_cast<unsigned int*>(fIntDataPtr), tmpLength);  removed 18 april 2008   
	  double tmpMax = fUtilitiesPtr->MaxValue(fShortDataPtr, tmpLength); 

	  if(tmpMax == 1023)
	    {
	      fDAmpl = tmpMax;
	    }
	}

      fDTof = fDTof/fDAmpl;

    }
} //end Evaluate




