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


AliHLTCaloRawAnalyzerCrude::AliHLTCaloRawAnalyzerCrude():AliHLTCaloRawAnalyzer() 
{

}


AliHLTCaloRawAnalyzerCrude::~AliHLTCaloRawAnalyzerCrude()
{

} //end AliHLTCaloRawAnalyzerCrude


void 
AliHLTCaloRawAnalyzerCrude::Evaluate(int start, int length)
{
  fDAmpl =0;
  fDTof = 0;
  int i=0;
  
  for(i=start; i<length; ++i)
      {
	if(fShortDataPtr[i] > fDAmpl  )
	  {
	    fDAmpl   = fShortDataPtr[i];
	    fDTof  = i;		     
	  }
      }
} //end Crude

