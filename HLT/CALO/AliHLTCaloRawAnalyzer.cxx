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

AliHLTCaloRawAnalyzer:: AliHLTCaloRawAnalyzer():fDoCorrectBaselineUsingFirstFiveSamples(false),
						fShortDataPtr(0),
						fSampleFrequency(10),
						fTau(2), 
						fDTof(99999), 
						fDAmpl(99999),
						fUtilitiesPtr(0)

{
  fUtilitiesPtr = new  AliHLTCaloUtilities(); 
}


AliHLTCaloRawAnalyzer::~AliHLTCaloRawAnalyzer()
{
  delete  fUtilitiesPtr;
}


void 
AliHLTCaloRawAnalyzer::SetCorrectBaselineUsingFirstFiveSamples()
{
  fDoCorrectBaselineUsingFirstFiveSamples = true;
}
 

void 
AliHLTCaloRawAnalyzer::CorrectBaselineUsingFirstFiveSamples(UShort_t *data, const int length)
{
  UShort_t sumOfFirstFiveSamples = 0;
  for(int i=0; i< 5; i++)
    {
      sumOfFirstFiveSamples += data[i];
    }
  UShort_t valueToSubtract = sumOfFirstFiveSamples/5;
  for(int j = 0; j < length; j++)
    {
      if( ( data[j] - valueToSubtract) > 0)
	{
	  data[j] = data[j] - valueToSubtract;
	}
      else
	{
	  data[j] = 0;
	}
    }
}


void
AliHLTCaloRawAnalyzer::SetData(const UShort_t *data, const int length) 
{
  fShortDataPtr = const_cast<UShort_t *>(data);
  if(fDoCorrectBaselineUsingFirstFiveSamples == true)
     {
       CorrectBaselineUsingFirstFiveSamples( fShortDataPtr, length );
     }
 }




