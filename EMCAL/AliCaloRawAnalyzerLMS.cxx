/**************************************************************************
 * This file is property of and copyright by                              *
 * the Relativistic Heavy Ion Group (RHIG), Yale University, US, 2009     *
 *                                                                        *
 * Primary Author: Per Thomas Hille <p.t.hille@fys.uio.no>                *
 *                                                                        *
 * Contributors are mentioned in the code where appropriate.              *
 * Please report bugs to p.t.hille@fys.uio.no                             *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


// Extraction of amplitude and peak position
// FRom CALO raw data using
// least square fit for the
// Moment assuming identical and 
// independent errors (equivalent with chi square)
// 

#include "AliCaloRawAnalyzerLMS.h"
#include "AliCaloBunchInfo.h"
#include "AliCaloFitResults.h"
#include "AliLog.h"
#include "TMath.h"
#include <iostream>
#include "TF1.h"
#include "TGraph.h"

using namespace std;


#define BAD 4  //CRAP PTH
 
ClassImp( AliCaloRawAnalyzerLMS )


AliCaloRawAnalyzerLMS::AliCaloRawAnalyzerLMS() : AliCaloRawAnalyzer("Chi Square Fit", "LMS"),
						 fkEulerSquared(7.389056098930650227),
						 fSig(0),
						 fTf1(0)
{
  //comment
  for(int i=0; i < MAXSAMPLES; i++)
    {
      fXaxis[i] = i;
    }
  
  fTf1 = new TF1( "myformula", "[0]*((x - [1])/2.35)^2*exp(-2*(x -[1])/2.35)",  0, 30 ); 
  //  fTf1 = new TF1( "myformula", "[0]*((x - [1])/[2])^[3]*exp(-[3]*(x -[1])/[2])",  0, 30 ); 

}


AliCaloRawAnalyzerLMS::~AliCaloRawAnalyzerLMS()
{
  delete fTf1;
}



AliCaloFitResults
AliCaloRawAnalyzerLMS::Evaluate( const vector<AliCaloBunchInfo>  &bunchvector, const UInt_t altrocfg1,  const UInt_t altrocfg2 )
{
  // Extracting signal parameters using fitting
  short maxampindex; //index of maximum amplitude
  short maxamp; //Maximum amplitude
  int index = SelectBunch( bunchvector,  &maxampindex,  &maxamp );
  
  if( index >= 0)
    {
      Float_t  ped  = ReverseAndSubtractPed( &(bunchvector.at(index))  ,  altrocfg1, altrocfg2, fReversed  );
      int first;
      int last;
      Float_t maxf = TMath::MaxElement( bunchvector.at(index).GetLength(),  fReversed );

      if ( maxf > fAmpCut )
	{
	  short maxrev = maxampindex  -  bunchvector.at(index).GetStartBin();
	  short timebinOffset = maxampindex - (bunchvector.at(index).GetLength()-1);
	  SelectSubarray( fReversed,  bunchvector.at(index).GetLength(),  maxrev, &first, &last);
	  int nsamples =  last - first + 1;
	  
	  if( ( nsamples  )  > fNsampleCut )
	    {
	      
	      TGraph *graph =  new TGraph(  nsamples, fXaxis,  &fReversed[first] );
	      fTf1->SetParameter(0, maxf*fkEulerSquared );
	      fTf1->SetParameter(1, 0.2);
	      //	      fTf1->SetParameter(2,  2);
	      //	      fTf1->SetParameter(2,  3); 
	    
	      //     return   AliCaloFitResults( -1 , -1 , -1, -1, -1, -1 , -1); 

	      Short_t tmpStatus =  graph->Fit(fTf1, "Q0RW");
	     
	      //       return   AliCaloFitResults( -1 , -1 , -1, -1, -1, -1 , -1); 
 
	      if( fVerbose == true )
		{
		  AliCaloRawAnalyzer::PrintBunch( bunchvector.at(index) ); 
		  PrintFitResult( fTf1 ) ;
		}  
	      
	        delete graph;
		return AliCaloFitResults( maxamp, ped ,    tmpStatus,  
					 fTf1->GetParameter(0)/fkEulerSquared, 
					 fTf1->GetParameter(1) + timebinOffset,  
					 fTf1->GetChisquare(), 
					 fTf1->GetNDF());
		
		
		//     delete graph;
	
	    }
	  else
	    {
	      return AliCaloFitResults( maxamp, ped, -1, maxf, -1, -1, -1 );
	    }
	}
      else
	{
	  return AliCaloFitResults( maxamp , ped, -1, maxf, -1, -1, -1 );
	}       
    }
  else
    {
      return AliCaloFitResults( -99, -99 );
    }

  return AliCaloFitResults( -99, -99 );
  
}


void 
AliCaloRawAnalyzerLMS::PrintFitResult(const TF1 *f) const
{
  //comment
  cout << endl;
  cout << __FILE__ << __LINE__ << "Using this samplerange we get" << endl;
  cout << __FILE__ << __LINE__ << "AMPLITUDE = " << f->GetParameter(0)/fkEulerSquared << ",.. !!!!" << endl;
  cout << __FILE__ << __LINE__ << "TOF = " << f->GetParameter(1) << ",.. !!!!" << endl;
  cout << __FILE__ << __LINE__ << "NDF = " << f->GetNDF() << ",.. !!!!" << endl;
  //  cout << __FILE__ << __LINE__ << "STATUS = " << f->GetStatus()  << ",.. !!!!" << endl << endl;
  cout << endl << endl;
}
