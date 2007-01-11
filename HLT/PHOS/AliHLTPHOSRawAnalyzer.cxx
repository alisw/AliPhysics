#include "AliHLTPHOSRawAnalyzer.h"
#include <iostream>

using std::cout;
using std::endl;



AliHLTPHOSRawAnalyzer:: AliHLTPHOSRawAnalyzer():fFloatDataPtr(0), fSampleFrequency(10), fTau(2), fDTof(99999), fDAmpl(99999), n(99999)
{


}

//int AliHLTPHOSRawAnalyzer::Deinit()
//{
//  return(0);
//}
AliHLTPHOSRawAnalyzer::~AliHLTPHOSRawAnalyzer()
{

}

AliHLTPHOSRawAnalyzer::AliHLTPHOSRawAnalyzer(AliHLTPHOSRawAnalyzer const&):fFloatDataPtr(0), fSampleFrequency(10), fTau(2), fDTof(99999), fDAmpl(99999), n(99999)
{

}

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
  double *dta;
  dta = data;
  cout << "Set data not yet implemented" << endl;
}
