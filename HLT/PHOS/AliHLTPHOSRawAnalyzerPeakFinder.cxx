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

#include "AliHLTPHOSRawAnalyzerPeakFinder.h"
#include <iostream>
#include <cmath>

using std::cout;
using std::endl;

ClassImp(AliHLTPHOSRawAnalyzerPeakFinder) 


AliHLTPHOSRawAnalyzerPeakFinder::AliHLTPHOSRawAnalyzerPeakFinder(const AliHLTPHOSRawAnalyzerPeakFinder&):AliHLTPHOSRawAnalyzer() , fTVectorPtr(0), fAVectorPtr(0), fTVectorSize(0), fAVectorSize(0)
{


}


/**
 * The AliHLTPHOSPeakfinder class is the class for extracting the basic signal parameters
 * "timing" and "energy" from the PHOS raw data. Physical data will for a given readout channel be
 * a sequense of ADC digitized 10 bit integer values, however for performance reasons all values used in
 * calculation is of type double.
 **/
AliHLTPHOSRawAnalyzerPeakFinder::AliHLTPHOSRawAnalyzerPeakFinder():AliHLTPHOSRawAnalyzer(), fTVectorPtr(0), fAVectorPtr(0), fTVectorSize(0), fAVectorSize(0)
{
  //  cout <<"PeakFinder:You cannot invoke the Fitter without arguments"<<endl;;
}



AliHLTPHOSRawAnalyzerPeakFinder::~AliHLTPHOSRawAnalyzerPeakFinder()
{

} //end AliHLTPHOSRawAnalyzerPeakFinder

void 
AliHLTPHOSRawAnalyzerPeakFinder::SetTVector(Double_t *tVec, Int_t size)
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


void
AliHLTPHOSRawAnalyzerPeakFinder::SetAVector(Double_t *aVec, Int_t size)
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



/**
* Extraction of timing and energy using the Peakfinde Algorithm.
* The. The parameters "start" and "length" defines a sub array  of the data array
* that will be used for the the fit. If start+length must not exeed the total length
* of the Data array. "start" must be chosen as close as possible to t0.
* The baseline must also be subtracted.
* The length of "tVector" and "aVector" mus be equal to length.
* "index + length" must not exeed the length of the data array set in the constructor.
* @param start the start index of the subarray of the data array. 
* @param length the number of samples to use starting from index 
* @param tVector the peakfinder vector for timing
* @param aVector the peakfinder vector for amplitude (energy)
**/
void 
AliHLTPHOSRawAnalyzerPeakFinder::Evaluate(int start, int length)
{
  //  printf("\n AliHLTPHOSRawAnalyzerPeakFinder::Evaluat from index %d to %d\n", start, start + length);
  fDTof = 0;
  fDAmpl = 0;
  Int_t tmpLength;


  if(fTVectorPtr == 0 || fAVectorPtr == 0)
    {
      printf("\nError: the peakfinder vectors are not specified, aborting !!!\n");
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
	  fDAmpl += fAVectorPtr[i]*fFloatDataPtr[i];    
	}

      for(int i=0; i < tmpLength; i++)
	{   
	  fDTof += fTVectorPtr[i]*fFloatDataPtr[i]; 
	}
      
      if(fDAmpl > 900)
	{
	  Double_t tmpMax = GetMaxValue(fFloatDataPtr, tmpLength);
	  if(tmpMax == 1023)
	    {
	      fDAmpl = tmpMax;
	    }
	}

      fDTof = fDTof/fDAmpl;

    }
} //end Evaluate




