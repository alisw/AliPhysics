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

#include "AliHLTPHOSRawAnalyzer.h"
//#include <iostream>

//using std::cout;
//using std::endl;

AliHLTPHOSRawAnalyzer:: AliHLTPHOSRawAnalyzer(): AliHLTPHOSBase(), fIntDataPtr(0), fSampleFrequency(10), fTau(2), fDTof(99999), fDAmpl(99999)
{
 fIntDataPtr = new UInt_t[1008];
}

AliHLTPHOSRawAnalyzer::~AliHLTPHOSRawAnalyzer()
{
  delete[] fIntDataPtr;
}



/**
* Main constructor
* @param dtaPtr Data array for wich a subarray will be taken to perform the fit
* @param fs the sampling frequency in entities of MHz. Needed in order to calculate physical time
**/
AliHLTPHOSRawAnalyzer::AliHLTPHOSRawAnalyzer(double *dtaPtr, double fs): AliHLTPHOSBase(), fIntDataPtr(0), fSampleFrequency(10), fTau(2), fDTof(99999), fDAmpl(99999), fStartIndex(0)
{
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
  cout << "Baseline correction not yet implemeted" << endl;
} //end BaselineCorrection


/**
* Shifts the baseline with the amount given by baselineValue
* If pretrigger samples are not present then the basline correction will be incorrect. 
* @param dataPtr array for wich to correct the basline 
* @param baselineValue the basline value to subtract..
**/
void 
AliHLTPHOSRawAnalyzer::BaselineCorrection(double *dataPtr, double baselineValue)
{
  printf("\nbaselineValue = %f\n", baselineValue);
  cout << "Baseline correction not yet implemeted" << endl;
} //end BaslineCorrection


/**
 * Gives the timing in entities of sample indexes
 * Physical time is found by multiplying  with the sampling intervall (Ts).
 **/
float
AliHLTPHOSRawAnalyzer::GetTiming() const 
{
  return fDTof;
} //end GetTiming


/**
 * Gives the time in entities of ADC channels (quantization levels).  
 * Absolute enrgy is found by multiplying with offline calibration constants.
 **/
float
AliHLTPHOSRawAnalyzer::GetEnergy() const
{
  return fDAmpl;
} //end GetEnergy


/**
 * Set data array. Overrides data data array set in the constructor.
 **/
void 
AliHLTPHOSRawAnalyzer::SetData(UInt_t *data)
{
  fIntDataPtr = data;
}

/**
 * Set data array. Overrides data data array set in the constructor.
 **/
void 
AliHLTPHOSRawAnalyzer::SetData(double *data)
{
  fFloatDataPtr = data;
}



void 
AliHLTPHOSRawAnalyzer::SetSampleFreq(double freq)
{
  fSampleFrequency = freq;
}

int 
AliHLTPHOSRawAnalyzer::FindStartIndex(double treshold)
{
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
  cout << "Make initial guess not yet implemeted" << endl;  
}


void
AliHLTPHOSRawAnalyzer::SetStartIndex(int index)
{
  fStartIndex = index;
}



void 
AliHLTPHOSRawAnalyzer::SetTVector(Double_t *tVector, Int_t size)
{
  cout <<"ERROR: AliHLTPHOSRawAnalyzer::SetTVector:  You cannot set the peakfindervector here, must be set in derived class peakfinder"<<endl;
}



void
AliHLTPHOSRawAnalyzer::SetAVector(Double_t *aVector, Int_t size)
{
 cout <<"ERROR: AliHLTPHOSRawAnalyzer::SetAVector:  You cannot set the peakfindervector here, must be set in derived class peakfinder"<<endl;
}

/*
UInt_t
AliHLTPHOSRawAnalyzer::GetMaxValue(UInt_t *dta, Int_t size) const
{

  Double_t tmpMax = 0;

  for(int i = 0; i < size; i++)
    {
      if(dta[i] > tmpMax)
	{
	  tmpMax = dta[i];
	}
    }
  
  return tmpMax;

}
*/
