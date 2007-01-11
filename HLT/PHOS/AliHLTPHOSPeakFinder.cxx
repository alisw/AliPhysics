
/**************************************************************************
 * Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved.      *
 *                                                                        *
 * Author: Per Thomas Hille for the ALICE HLT Project.                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/



#include "AliHLTPHOSPeakFinder.h"
#include <iostream>

using std::cout;
using std::endl;

ClassImp(AliHLTPHOSPeakFinder) 


AliHLTPHOSPeakFinder::AliHLTPHOSPeakFinder(const AliHLTPHOSPeakFinder&):TObject(), AliHLTPHOSRawAnalyzer(), fDTofGuess(0), fDAmplGuess(0), kfMCovarPtrPtr(0), fPCovarPtrPtr(0)
{

}

//AliHLTComponent::GetComponentID()}

/**
 * The AliHLTPHOSPeakfinder class is the class for extracting the basic signal parameters
 * "timing" and "energy" from the PHOS raw data. Physical data will for a given readout channel be
 * a sequense of ADC digitized 10 bit integer values, however for performance reasons all values used in
 * calculation is of type double.
 **/
AliHLTPHOSPeakFinder::AliHLTPHOSPeakFinder():fDTofGuess(0), fDAmplGuess(0), kfMCovarPtrPtr(0), fPCovarPtrPtr(0) 
{
  cout <<"You cannot invoke the Fitter without arguments"<<endl;;
}


/**
* Main constructor
* @param dataPtr Data array for wich a subarray will be taken to perform the fit
* @param fs the sampling frequency in entities of MHz. Needed in order to calculate physical time
**/
AliHLTPHOSPeakFinder::AliHLTPHOSPeakFinder(double *dtaPtr, double fs): fDTofGuess(0), fDAmplGuess(0), kfMCovarPtrPtr(0), fPCovarPtrPtr(0) 
{
  fFloatDataPtr = dtaPtr;  
  fSampleFrequency = fs;
} //end   AliHLTPHOSPeakFinder 


AliHLTPHOSPeakFinder::~AliHLTPHOSPeakFinder()
{

} //end AliHLTPHOSPeakFinder


void
AliHLTPHOSPeakFinder::Analyze() const
{

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
AliHLTPHOSPeakFinder::FitPeakFinder(int start, int length, double *tVector, double *aVector)
{
  fDTof = 0;
  fDAmpl = 0;

  printf("\nstart = %d, length = %d\n", start, length);   

  double tmpTime[1008];

  //  double tmpTime[length];

  for(int i=0; i < length; i++)
    {  
      fDAmpl += aVector[i]*fFloatDataPtr[i];    
    }
  
  for(int i=0; i < length; i++)
    {   
      tmpTime[i] = tVector[i]*fFloatDataPtr[i];
      fDTof = fDTof + tmpTime[i]; 
   }

  fDTof = fDTof/fDAmpl;
  //thats all 
} //end FitPeakFinder



/**
 * This method finds the start index of the pulse (relative to the data array) by serching for
 * three or more continious samples above trheshold.
 * @param treshold treshold to use when searchin for the start of the pulse
 **/
int 
AliHLTPHOSPeakFinder::FindStartIndex(double treshold)
{
  printf("\ntreshold = %f \n", treshold);
  cout << "Find Start index not yet implemented" << endl;
  return 0;
} //end FindStartIndex


/**
 * This function applies only to the Chi and Least mean square fit. An initial guess is made
 * based on the average of the first 5 samples and the first value exeeding this value.
 **/
void 
AliHLTPHOSPeakFinder::MakeInitialGuess()
{
  cout << "Make initial guess not yet implemeted" << endl;
}


/**
 * This function applies only to the Chi and Least mean square fit. An initial guess is made
 * based on the average of the first 5 samples and the first value exeeding threshold + this value.
 * @param treshold The index of the first value above treshold is ntaken to be the first value.
 **/
void 
AliHLTPHOSPeakFinder::MakeInitialGuess(int treshold)
{
  printf("\ntreshold = %d\n", treshold);
  cout << "Make initial guess not yet implemeted" << endl;  
}


