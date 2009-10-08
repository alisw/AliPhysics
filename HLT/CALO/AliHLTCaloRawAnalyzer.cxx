// $Id: AliHLTCALORawAnalyzer.cxx 34264 2009-08-14 18:29:23Z odjuvsla $

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


#include "AliHLTCaloRawAnalyzer.h"
#include "AliHLTCaloUtilities.h" 

AliHLTCaloRawAnalyzer:: AliHLTCaloRawAnalyzer(): //AliHLTCaloBase(),   
						 //	 fDoCorrectBaselineUsingFirstFiveSamples(false),
						 fDoCorrectBaselineUsingFirstFiveSamples(true),
						 fDoubleDataPtr(0), 
						 fIntDataPtr(0), 
						 fShortDataPtr(0),
						 fSampleFrequency(10),
						 fDTofGuess(-1),
						 fDAmplGuess(-1),
						 fTau(2), 
						 fDTof(99999), 
						 fDAmpl(99999),
						 fStartIndex(0),
						 fUseShortValues(false),
						 fUtilitiesPtr(0)

{
  //  fIntDataPtr = new UInt_t[1008];

  //  fDoubleDataPtr;   
  fUtilitiesPtr = new  AliHLTCaloUtilities(); 
}

AliHLTCaloRawAnalyzer::~AliHLTCaloRawAnalyzer()
{
  //  delete[] fIntDataPtr;
}





/**
* Main constructor
* param dtaPtr Data array for wich a subarray will be taken to perform the fit
* @param fs the sampling frequency in entities of MHz. Needed in order to calculate physical time
**/
AliHLTCaloRawAnalyzer::AliHLTCaloRawAnalyzer(double * /*dtaPtr*/, double fs): //AliHLTCaloBase(), 
									      fDoCorrectBaselineUsingFirstFiveSamples(false),
									      fDoubleDataPtr(0), 
									      fIntDataPtr(0), 
									      fShortDataPtr(0),
									      fSampleFrequency(10),
									      fDTofGuess(-1),
									      fDAmplGuess(-1),
									      fTau(2), 
									      fDTof(99999), 
									      fDAmpl(99999),
									      fStartIndex(0),
									      fUseShortValues(false),
									      fUtilitiesPtr(0)
									      
{
  fSampleFrequency = fs;
} //end  



void 
AliHLTCaloRawAnalyzer::SetCorrectBaselineUsingFirstFiveSamples()
{
  fDoCorrectBaselineUsingFirstFiveSamples = true;
}
 

void 
//AliHLTCaloRawAnalyzer::CorrectBaselineUsingFirstFiveSamples(double *data, int length)
//AliHLTCaloRawAnalyzer::CorrectBaselineUsingFirstFiveSamples(int *data, int length)
AliHLTCaloRawAnalyzer::CorrectBaselineUsingFirstFiveSamples(UInt_t */*data*/, const int /*length*/)
{
  //  cout << "AliHLTCaloRawAnalyzer::CorrectBaselineUsingFirstFiveSamples" << endl;

  //CRAP PTH
  /*
  
  unsigned int sumOfFirstFiveSamples = 0;

  for(int i=0; i< 5; i++)
    {
      sumOfFirstFiveSamples += data[i];
    }

  unsigned int valueToSubtract = sumOfFirstFiveSamples/5;

  for(int j = 0; j < length; j++)
    {
      if( (int)(data[j] - valueToSubtract) > 0)
	{
	  data[j] = data[j] - valueToSubtract;
	}
      else
	{
	  data[j] = 0;
	}
    }
  */
}



// /**
// * Attemps to level the basline to zero.
// * The baseline will be calculated from the pretrigger samples and subtracted from the 
// * data array. 
// * If pretrigger samples are not present then the basline correction will be incorrect. 
// * @param dataPtr array for wich to correct the basline 
// * @param N the number of pretrigger samples used to calculate the baseline.
// **/
// void 
// AliHLTCaloRawAnalyzer::BaselineCorrection(double * /*dataPtr*/, int /*N*/)
// {

// } //end BaselineCorrection


/**
* Shifts the baseline with the amount given by baselineValue
* If pretrigger samples are not present then the basline correction will be incorrect. 
* @param dataPtr array for wich to correct the basline 
* @param baselineValue the basline value to subtract..
**/
void 
AliHLTCaloRawAnalyzer::BaselineCorrection(double * /*dataPtr*/, double /*baselineValue*/)
{

} //end BaslineCorrection


/**
 * Gives the timing in entities of sample indexes
 * Physical time is found by multiplying  with the sampling intervall (Ts).
 **/
float
AliHLTCaloRawAnalyzer::GetTiming() const 
{
  return fDTof;
} //end GetTiming


/**
 * Gives the time in entities of ADC channels (quantization levels).  
 * Absolute enrgy is found by multiplying with offline calibration constants.
 **/
float
AliHLTCaloRawAnalyzer::GetEnergy() const
{
  return fDAmpl;
} //end GetEnergy



/**
 * Set data array. Overrides data data array set in the constructor.
 **/
 void 

 AliHLTCaloRawAnalyzer::SetData(const UInt_t *data, const int length) 
 // AliHLTCaloRawAnalyzer::SetData(UInt_t *data, const int length) 
// AliHLTCaloRawAnalyzer::SetData(int *data, const int length) 
 {
   fIntDataPtr = const_cast<UInt_t *>(data);

   if(fDoCorrectBaselineUsingFirstFiveSamples == true)
     {
       CorrectBaselineUsingFirstFiveSamples(fIntDataPtr, length);
     }

   //   fIntDataPtr = data;

 }

void
 AliHLTCaloRawAnalyzer::SetData(const UShort_t *data, const int length) 
 // AliHLTCaloRawAnalyzer::SetData(UInt_t *data, const int length) 
// AliHLTCaloRawAnalyzer::SetData(int *data, const int length) 
 {

   fShortDataPtr = const_cast<UShort_t *>(data);
   fUseShortValues = true;
   if(fDoCorrectBaselineUsingFirstFiveSamples == true)
     {
       CorrectBaselineUsingFirstFiveSamples(fIntDataPtr, length);
     }

   //   fIntDataPtr = data;

 }




/**
 * Set data array. Overrides data data array set in the constructor.
 **/
// void 
// //AliHLTCaloRawAnalyzer::SetData(const double *data) 
// AliHLTCaloRawAnalyzer::SetData(double *data, const int length) 
// {
//   if(fDoCorrectBaselineUsingFirstFiveSamples == true)
//     {
//       CorrectBaselineUsingFirstFiveSamples(data, length);
//     }


//   fDoubleDataPtr = data;
// }



void 
AliHLTCaloRawAnalyzer::SetSampleFreq(double freq)
{
  fSampleFrequency = freq;
}

int 
AliHLTCaloRawAnalyzer::FindStartIndex(double /*treshold*/)
{
  // cout << "Find Start index not yet implemented" << endl;
  return 0;
} //end FindStartIndex


/**
 * This function applies only to the Chi and Least mean square fit. An initial guess is made
 * based on the average of the first 5 samples and the first value exeeding this value.
 **/
void 
AliHLTCaloRawAnalyzer::MakeInitialGuess()
{
  //  cout << "Make initial guess not yet implemeted" << endl;
}


/**
 * This function applies only to the Chi and Least mean square fit. An initial guess is made
 * based on the average of the first 5 samples and the first value exeeding threshold + this value.
 * @param treshold The index of the first value above treshold is ntaken to be the first value.
 **/
void 
AliHLTCaloRawAnalyzer::MakeInitialGuess(int /*treshold*/)
{
  //  cout << "Make initial guess not yet implemeted" << endl;  
}


void
AliHLTCaloRawAnalyzer::SetStartIndex(int index)
{
  fStartIndex = index;
}



void 
AliHLTCaloRawAnalyzer::SetTVector(Double_t * /*tVector*/, Int_t /*size*/)
{
  //  cout <<"ERROR: AliHLTCaloRawAnalyzer::SetTVector:  You cannot set the peakfindervector here, must be set in derived class peakfinder"<<endl;
}



void
AliHLTCaloRawAnalyzer::SetAVector(Double_t * /*aVector*/, Int_t /*size*/)
{
  // cout <<"ERROR: AliHLTCaloRawAnalyzer::SetAVector:  You cannot set the peakfindervector here, must be set in derived class peakfinder"<<endl;
}


