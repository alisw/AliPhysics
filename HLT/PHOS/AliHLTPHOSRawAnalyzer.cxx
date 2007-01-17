#include "AliHLTPHOSRawAnalyzer.h"
#include <iostream>

using std::cout;
using std::endl;



AliHLTPHOSRawAnalyzer:: AliHLTPHOSRawAnalyzer():fFloatDataPtr(0), fSampleFrequency(10), fTau(2), fDTof(99999), fDAmpl(99999), n(99999)
{
cout <<"AliHLTPHOSRawAnalyzer Base class::Creating a new AliHLTPHOSRawAnalyzer() object"<<endl;
}

AliHLTPHOSRawAnalyzer::~AliHLTPHOSRawAnalyzer()
{

}

AliHLTPHOSRawAnalyzer::AliHLTPHOSRawAnalyzer(AliHLTPHOSRawAnalyzer const&):fFloatDataPtr(0), fSampleFrequency(10), fTau(2), fDTof(99999), fDAmpl(99999), n(99999)
{

}

/**
* Main constructor
* @param dataPtr Data array for wich a subarray will be taken to perform the fit
* @param fs the sampling frequency in entities of MHz. Needed in order to calculate physical time
**/
AliHLTPHOSRawAnalyzer::AliHLTPHOSRawAnalyzer(double *dtaPtr, double fs):fFloatDataPtr(0), fSampleFrequency(10), fTau(2), fDTof(99999), fDAmpl(99999), n(99999)
{
  fFloatDataPtr = dtaPtr;  
  fSampleFrequency = fs;
} //end  


/**
* Attemps to level the basline to zero.
* The baseline will be calculated from the pretrigger samples and subtracted from the 
* data array. 
* If pretrigger samples are not present then the basline correction will be incorrect. 
* @param dataPtr array for wich to correct the basline 
* @param N the number of pretrigger samples used to calculate the baseline.
**/
void 
AliHLTPHOSRawAnalyzer::BaselineCorrection(double *dataPtr, int N)
{
  fFloatDataPtr = dataPtr;  
  n=N;
  cout << "Baseline correction not yet implemeted" << endl;
} //end BaselineCorrection


/**
* Shifts the basline with the amount given by baselineValue
* If pretrigger samples are not present then the basline correction will be incorrect. 
* @param dataPtr array for wich to correct the basline 
* @param BaslineValue the basline value to subtract..
**/
void 
AliHLTPHOSRawAnalyzer::BaselineCorrection(double *dataPtr, double baselineValue)
{
  fFloatDataPtr = dataPtr;   
  printf("\nbaselineValue = %f\n", baselineValue);
  cout << "Baseline correction not yet implemeted" << endl;
} //end BaslineCorrection


/**
 * Gives the timing in entities of sample indexes
 * Physical time is found by multiplying  with the sampling intervall (Ts).
 **/
float
AliHLTPHOSRawAnalyzer::GetTiming()
{
  return fDTof;
} //end GetTiming


/**
 * Gives the time in entities of ADC channels (quantization levels).  
 * Absolute enrgy is found by multiplying with offline calibration constants.
 **/
float
AliHLTPHOSRawAnalyzer::GetEnergy()
{
  return fDAmpl;
} //end GetEnergy


/**
 * Set data array. Overrides data data array set in the constructor.
 **/
void 
AliHLTPHOSRawAnalyzer::SetData(double *data)
{
  //  double *dta;
  //  dta = data;
  fFloatDataPtr = data;
  //  cout << "Set data not yet implemented" << endl;
}

void 
AliHLTPHOSRawAnalyzer::SetSampleFreq(double freq)
{

}

int 
AliHLTPHOSRawAnalyzer::FindStartIndex(double treshold)
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
AliHLTPHOSRawAnalyzer::MakeInitialGuess()
{
  cout << "Make initial guess not yet implemeted" << endl;
}


/**
 * This function applies only to the Chi and Least mean square fit. An initial guess is made
 * based on the average of the first 5 samples and the first value exeeding threshold + this value.
 * @param treshold The index of the first value above treshold is ntaken to be the first value.
 **/
void 
AliHLTPHOSRawAnalyzer::MakeInitialGuess(int treshold)
{
  printf("\ntreshold = %d\n", treshold);
  cout << "Make initial guess not yet implemeted" << endl;  
}
