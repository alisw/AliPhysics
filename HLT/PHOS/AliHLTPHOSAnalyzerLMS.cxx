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


#include "AliHLTPHOSAnalyzerLMS.h"
#include <iostream>

using std::cout;
using std::endl;

ClassImp(AliHLTPHOSAnalyzerLMS) 


AliHLTPHOSAnalyzerLMS::AliHLTPHOSAnalyzerLMS(const AliHLTPHOSAnalyzerLMS&):AliHLTPHOSAnalyzer(), kfMCovarPtrPtr(0), fPCovarPtrPtr(0)
{

}


/**
 * The AliHLTPHOSPeakfinder class is the class for extracting the basic signal parameters
 * "timing" and "energy" from the PHOS raw data. Physical data will for a given readout channel be
 * a sequense of ADC digitized 10 bit integer values, however for performance reasons all values used in
 * calculation is of type double.
 **/
AliHLTPHOSAnalyzerLMS::AliHLTPHOSAnalyzerLMS():AliHLTPHOSAnalyzer(), kfMCovarPtrPtr(0), fPCovarPtrPtr(0) 
{
  cout <<"You cannot invoke the Fitter without arguments"<<endl;;
}


/**
* Main constructor
* @param dataPtr Data array for wich a subarray will be taken to perform the fit
* @param fs the sampling frequency in entities of MHz. Needed in order to calculate physical time
**/
AliHLTPHOSAnalyzerLMS::AliHLTPHOSAnalyzerLMS(double *dtaPtr, double fs):AliHLTPHOSAnalyzer(),kfMCovarPtrPtr(0), fPCovarPtrPtr(0) 
{
  fFloatDataPtr = dtaPtr;  
  fSampleFrequency = fs;
} //end   AliHLTPHOSAnalyzerLMS 


AliHLTPHOSAnalyzerLMS::~AliHLTPHOSAnalyzerLMS()
{

} //end AliHLTPHOSAnalyzerLMS


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
**/
void 
AliHLTPHOSAnalyzerLMS::Evaluate(int start, int length)
{
  /*

  */

  //thats all 
} //end FitLMS





