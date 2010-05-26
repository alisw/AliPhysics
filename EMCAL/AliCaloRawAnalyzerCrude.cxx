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
#include "TMath.h"

using namespace std;

ClassImp(AliCaloRawAnalyzerCrude)  


AliCaloRawAnalyzerCrude::AliCaloRawAnalyzerCrude() : AliCaloRawAnalyzer("Crude", "Crude")
{

}


AliCaloRawAnalyzerCrude::~AliCaloRawAnalyzerCrude()
{

}


AliCaloFitResults
AliCaloRawAnalyzerCrude::Evaluate(const vector<AliCaloBunchInfo> &bunchvector, const UInt_t altrocfg1,  const UInt_t altrocfg2)
{
  // Evaluation of signal parameters
  short maxampindex; //index of maximum amplitude
  short maxamp; //Maximum amplitude
  int index = SelectBunch( bunchvector,  &maxampindex,  &maxamp );
 
  if( index >= 0)
    {
      Float_t ped = ReverseAndSubtractPed( &(bunchvector.at(index))  ,  altrocfg1, altrocfg2, fReversed  );
      Float_t maxf = TMath::MaxElement( bunchvector.at(index).GetLength(),  fReversed );
      short timebinOffset = maxampindex - (bunchvector.at(index).GetLength()-1);

      if(  maxf < fAmpCut  ||  ( maxamp - ped) > fOverflowCut  ) // (maxamp - ped) > fOverflowCut = Close to saturation (use low gain then)
	{
	  return  AliCaloFitResults( maxamp, ped, AliCaloFitResults::kCrude, maxf, timebinOffset);
	}
      else if ( maxf >= fAmpCut ) // no if statement needed really; keep for readability
	{
	  int first = 0;
	  int last = 0;
	  int maxrev =  maxampindex -  bunchvector.at(index).GetStartBin();
	  SelectSubarray( fReversed,  bunchvector.at(index).GetLength(), maxrev , &first, &last);

	  Float_t chi2 = CalculateChi2(maxf, maxrev, first, last);
	  Int_t ndf = last - first - 1; // nsamples - 2
	  return AliCaloFitResults( maxamp, ped, AliCaloFitResults::kCrude, maxf, timebinOffset,
				    timebinOffset, chi2, ndf, AliCaloFitResults::kDummy, AliCaloFitSubarray(index, maxrev, first, last) ); 
	} // ampcut
    } // bunch index    

  return AliCaloFitResults(AliCaloFitResults::kInvalid , AliCaloFitResults::kInvalid);

} //end Crude


