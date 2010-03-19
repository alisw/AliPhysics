/**************************************************************************
 * This file is property of and copyright by the Experimental Nuclear     *
 * Physics Group, Dep. of Physics                                         *
 * University of Oslo, Norway, 2007                                       *
 *                                                                        *
 * Author: Per Thomas Hille <perthi@fys.uio.no> for the ALICE HLT Project.*
 * Contributors are mentioned in the code where appropriate.              *
 * Please report bugs to perthi@fys.uio.no                                *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// Evaluation of amplitude
// as max sample value - pedestal
// Not veru accurate, but very robust
// --------------
// --------------

#include "AliCaloRawAnalyzerCrude.h"
#include "AliCaloFitResults.h"
#include "AliCaloBunchInfo.h"


ClassImp(AliCaloRawAnalyzerCrude)  


AliCaloRawAnalyzerCrude::AliCaloRawAnalyzerCrude() : AliCaloRawAnalyzer("Crude", "Crude")
{

}


AliCaloRawAnalyzerCrude::~AliCaloRawAnalyzerCrude()
{

}


AliCaloFitResults
AliCaloRawAnalyzerCrude::Evaluate(const vector<AliCaloBunchInfo> &bunchvector, const UInt_t /*altrocfg1*/,  const UInt_t /*altrocfg2*/)
{
  // Evaluation of signal parameters
  if( bunchvector.size()  <=  0 )
    {
      return AliCaloFitResults(AliCaloFitResults::kInvalid, AliCaloFitResults::kInvalid);
    }

  Int_t amp = 0;
  Int_t tof = -99;
  const UShort_t *sig;
  
  double ped = EvaluatePedestal( bunchvector.at(0).GetData(), bunchvector.at(0).GetLength() ) ;

  for( unsigned int i= 0; i < bunchvector.size(); ++i)
    {
      sig = bunchvector.at(i).GetData();
      int length = bunchvector.at(i).GetLength(); 
      
      for(int j = 0; j < length; j ++)
	if( sig[j] > amp  )
	  {
	    amp   = sig[j];
	    tof   = bunchvector.at(i).GetStartBin() - j;		     
	  }
    }

  //:EvaluatePedestal(const UShort_t * const data, const int length )
  //  double ped = EvaluatePedestal(sig, length) ;
  return  AliCaloFitResults(amp, ped, AliCaloFitResults::kCrude, amp - ped, tof);
  
} //end Crude


