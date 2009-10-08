// $Id: AliHLTCALORawAnalyzerCrude.cxx 34622 2009-09-04 13:22:01Z odjuvsla $

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

#include "AliHLTCaloRawAnalyzerCrude.h"

#include <iostream>

using namespace std;

//ClassImp(AliHLTCaloRawAnalyzerCrude) 


/**
 * The AliHLTPHOSPeakfCinder class is the class for extracting the basic signal parameters
 * "timing" and "energy" from the PHOS raw data. Physical data will for a given readout channel be
 * a sequense of ADC digitized 10 bit integer values, however for performance reasons all values used in
 * calculation is of type double.
 **/
//____________________________________________________________________________
AliHLTCaloRawAnalyzerCrude::AliHLTCaloRawAnalyzerCrude():AliHLTCaloRawAnalyzer() 
{

}


//____________________________________________________________________________
AliHLTCaloRawAnalyzerCrude::~AliHLTCaloRawAnalyzerCrude()
{

} //end AliHLTCaloRawAnalyzerCrude


/**
* Extraction of timing and energy using Crude estimate.
* The. The parameters "start" and "length" defines a sub array  of the data array
* that will be used for the the fit. If start+length must not exeed the total length
* of the Data array. "start" must be chosen as close as possible to t0.
* The baseline must also be subtracted.                                                              .
* "index + length" must not exeed the length of the data array set in the constructor.
* @param start the start index of the subarray of the data array. 
* @param length the number of samples to use starting from index 
**/
// //____________________________________________________________________________
// void 
// AliHLTCaloRawAnalyzerCrude::Evaluate(int start, int length)
// {
//   double tmpAmplitudeMax =0; 
//   double tmpTime = 0;

//   for(int i=start; i<length; i++)
//     {
//       if(fIntDataPtr[i] >  tmpAmplitudeMax && i > 5)
// 	{
// 	  tmpAmplitudeMax = fIntDataPtr[i];
// 	  tmpTime = i;		     
// 	}
//     }
	
//   fDAmpl = tmpAmplitudeMax;
//   fDTof =  tmpTime;
 
//   //thats all 
// } //end Crude


void 
AliHLTCaloRawAnalyzerCrude::Evaluate(int start, int length)
{
  //  cout << "AliHLTCaloRawAnalyzerCrude::Evaluate TP0"  << endl;

  //DumpData(T

  //  DumpData(fDoubleDataPtr,50, 25);

  if(fUseShortValues == true)
    {
      EvaluateShort(start, length);
    }
  else
    {

      double tmpAmplitudeMax =0; 
      double tmpTime = 0;

      for(int i=start; i<length; i++)
	{
	  //   if(fDoubleDataPtr[i] >  tmpAmplitudeMax && i > 5)
	  if(fIntDataPtr[i] >  tmpAmplitudeMax && i > 5)
	    {
	      //	  tmpAmplitudeMax = fDoubleDataPtr[i];
	      tmpAmplitudeMax = fIntDataPtr[i];
	      tmpTime = i;		     
	    }
	}

      fDAmpl = tmpAmplitudeMax;
      fDTof =  tmpTime;
    }
 
  //thats all 
} //end Crude

void 
AliHLTCaloRawAnalyzerCrude::EvaluateShort(int start, int length)
{
  //  cout << "AliHLTCaloRawAnalyzerCrude::Evaluate TP0"  << endl;

  //DumpData(T

  //  DumpData(fDoubleDataPtr,50, 25);

  double tmpAmplitudeMax =0; 
  double tmpTime = 0;

  for(int i=start; i<length; i++)
    {
      //      cout << "sample index: " << i << endl;
        //   if(fDoubleDataPtr[i] >  tmpAmplitudeMax && i > 5)
      if(fShortDataPtr[i] >  tmpAmplitudeMax && i > 5)
	{
	  //	  tmpAmplitudeMax = fDoubleDataPtr[i];
	  tmpAmplitudeMax = fShortDataPtr[i];
	  tmpTime = i;		     
	}
    }

 
  fDAmpl = tmpAmplitudeMax;
  fDTof =  tmpTime;
 
  //thats all 
} //end Crude
