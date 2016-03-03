/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
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

/*
  Author: R. GUERNANE LPSC Grenoble CNRS/IN2P3
*/


#include "AliCaloRawAnalyzerFakeALTRO.h"
#include "AliCaloBunchInfo.h"
#include "AliCaloFitResults.h"
#include "AliLog.h"
#include "TMath.h"
#include <stdexcept>
#include <iostream>
#include "TF1.h"
#include "TGraph.h"
#include "AliCaloConstants.h"

using namespace std;

ClassImp( AliCaloRawAnalyzerFakeALTRO )


AliCaloRawAnalyzerFakeALTRO::AliCaloRawAnalyzerFakeALTRO() : AliCaloRawAnalyzerFitter("Chi Square Fit", "FakeAltro")
{
  // constructor
  
  fAlgo= Algo::kFakeAltro;
}

AliCaloRawAnalyzerFakeALTRO::~AliCaloRawAnalyzerFakeALTRO()
{
  // destructor
  
  //delete fTf1;
}

AliCaloFitResults
AliCaloRawAnalyzerFakeALTRO::Evaluate( const vector<AliCaloBunchInfo>  &bunchvector,
                                      UInt_t altrocfg1, UInt_t altrocfg2 )
{
  // Extracting signal parameters using fitting
  
  short maxampindex; //index of maximum amplitude
  short maxamp; //Maximum amplitude
  int index = SelectBunch( bunchvector,  &maxampindex,  &maxamp );
  
  if( index >= 0)
  {
    Float_t  ped  = ReverseAndSubtractPed( &(bunchvector.at(index))  ,  altrocfg1, altrocfg2, fReversed  );
    Float_t maxf = TMath::MaxElement( bunchvector.at(index).GetLength(),  fReversed );
    short maxrev = maxampindex  -  bunchvector.at(index).GetStartBin();
    // timebinOffset is timebin value at maximum (maxrev)
    short timebinOffset = maxampindex - (bunchvector.at(index).GetLength()-1);
    Float_t time = (timebinOffset*TIMEBINWITH)-fL1Phase;
    if(  maxf < fAmpCut  ||  ( maxamp - ped) > fOverflowCut  ) // (maxamp - ped) > fOverflowCut = Close to saturation (use low gain then)
    {
      return  AliCaloFitResults( maxamp, ped, Ret::kCrude, maxf, time);
    }
    else if ( maxf >= fAmpCut )
    {
      int first = 0;
      int last = 0;
      SelectSubarray( fReversed,  bunchvector.at(index).GetLength(),  maxrev, &first, &last, fFitArrayCut );
      int nsamples =  last - first + 1;
      
      if( ( nsamples  )  >= fNsampleCut )
	    {
	      Float_t tmax = (maxrev - first); // local tmax estimate
	      TGraph *graph =  new TGraph(  nsamples, fXaxis,  &fReversed[first] );
	      fTf1->SetParameter(0, maxf*fkEulerSquared );
	      fTf1->SetParameter(1, tmax - fTau);
	      // set rather loose parameter limits
	      fTf1->SetParLimits(0, 0.5*maxf*fkEulerSquared, 2*maxf*fkEulerSquared );
	      fTf1->SetParLimits(1, tmax - fTau - 4, tmax - fTau + 4);
        
	      if (fFixTau) {
          fTf1->FixParameter(2, fTau);
	      }
	      else {
          fTf1->ReleaseParameter(2); // allow par. to vary
          fTf1->SetParameter(2, fTau);
	      }
        
	      Short_t tmpStatus = 0;
	      try {
          tmpStatus =  graph->Fit(fTf1, "Q0RW");
	      }
	      catch (const std::exception & e) {
          AliError( Form("TGraph Fit exception %s, fit status %d", e.what(),tmpStatus) );
          return AliCaloFitResults( maxamp, ped, Ret::kNoFit, maxf, time,
                                   (int)time, Ret::kDummy, Ret::kDummy,  Ret::kDummy, AliCaloFitSubarray(index, maxrev, first, last) );
	      }
        
	      if( fVerbose == true )
        {
          AliCaloRawAnalyzer::PrintBunch( bunchvector.at(index) );
          PrintFitResult( fTf1 ) ;
        }
	      // global tmax
	      tmax = fTf1->GetParameter(1) + timebinOffset - (maxrev - first) // abs. t0
        + fTf1->GetParameter(2); // +tau, makes sum tmax
	    Float_t timemax = (tmax*TIMEBINWITH)-fL1Phase;
        delete graph;
        return AliCaloFitResults( maxamp, ped , Ret::kFitPar,
                                 fTf1->GetParameter(0)/fkEulerSquared,
                                 timemax,
                                 (int)timemax,
                                 fTf1->GetChisquare(),
                                 fTf1->GetNDF(),
                                 Ret::kDummy, AliCaloFitSubarray(index, maxrev, first, last) );
				
        //     delete graph;
        
	    }
      else
	    {
	      Float_t chi2 = CalculateChi2(maxf, maxrev, first, last);
	      Int_t ndf = last - first - 1; // nsamples - 2
	      return AliCaloFitResults( maxamp, ped, Ret::kCrude, maxf, time,
                                 (int)time, chi2, ndf, Ret::kDummy, AliCaloFitSubarray(index, maxrev, first, last) );
	    }
    } // ampcut
  }
  return AliCaloFitResults(  Ret::kInvalid,  Ret::kInvalid );
}

