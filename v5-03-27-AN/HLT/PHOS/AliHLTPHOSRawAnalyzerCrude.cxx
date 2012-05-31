// $Id$

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

#include "AliHLTPHOSRawAnalyzerCrude.h"

//ClassImp(AliHLTPHOSRawAnalyzerCrude) 


/**
 * The AliHLTPHOSPeakfinder class is the class for extracting the basic signal parameters
 * "timing" and "energy" from the PHOS raw data. Physical data will for a given readout channel be
 * a sequense of ADC digitized 10 bit integer values, however for performance reasons all values used in
 * calculation is of type double.
 **/
//____________________________________________________________________________
AliHLTPHOSRawAnalyzerCrude::AliHLTPHOSRawAnalyzerCrude():AliHLTPHOSRawAnalyzer() 
{

}


//____________________________________________________________________________
AliHLTPHOSRawAnalyzerCrude::~AliHLTPHOSRawAnalyzerCrude()
{

} //end AliHLTPHOSRawAnalyzerCrude


void 
AliHLTPHOSRawAnalyzerCrude::Evaluate(int start, int length)
{
  double tmpAmplitudeMax =0; 
  double tmpTime = 0;

  for(int i=start; i<length; i++)
    {
      if(fShortDataPtr[i] >  tmpAmplitudeMax  )	
	{
	  tmpAmplitudeMax = fShortDataPtr[i];
	  tmpTime = i;		     
	}
    }

  fDAmpl = tmpAmplitudeMax;
  fDTof =  tmpTime;
  //thats all 
} //end Crude

