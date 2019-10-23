/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 *                                                                        *
 * Author: Martin Poghosyan                                               *
 * Fitting ALTRO data using gamma-2 function                              *
 **************************************************************************/

#include "TMath.h"
#include "TRandom.h"

#include "AliCaloRawAnalyzerGamma2.h"
#include "AliCaloBunchInfo.h"
#include "AliCaloFitResults.h"
#include "AliLog.h"


using namespace std;

ClassImp(AliCaloRawAnalyzerGamma2)

AliCaloRawAnalyzerGamma2::AliCaloRawAnalyzerGamma2():AliCaloRawAnalyzerFitter("Chi Square ( Gamma2 )", "KGamma2"),
fNiter(0),
fNiterationsMax(10)
{

}


AliCaloFitResults AliCaloRawAnalyzerGamma2::Evaluate( const vector<AliCaloBunchInfo>  &bunchlist,
                                       UInt_t altrocfg1, UInt_t altrocfg2 )
{
  Float_t pedEstimate  = 0;
  short maxADC = 0;
  Int_t first = 0;
  Int_t last = 0;
  Int_t bunchIndex = 0;
  Float_t ampEstimate = 0;
  short timeEstimate  = 0;
  Float_t time = 0;
  Float_t amp=0;
  Float_t chi2 = 0;
  Int_t ndf = 0;
  Bool_t fitDone = kFALSE;
  
  int nsamples = PreFitEvaluateSamples( bunchlist, altrocfg1, altrocfg2, bunchIndex, ampEstimate, 
                                       maxADC, timeEstimate, pedEstimate, first, last,   (int)fAmpCut ); 
  
  
  if (ampEstimate >= fAmpCut  ) 
  { 
    time = timeEstimate; 
    Int_t timebinOffset = bunchlist.at(bunchIndex).GetStartBin() - (bunchlist.at(bunchIndex).GetLength()-1); 
    amp = ampEstimate; 



    
    if ( nsamples >2 && maxADC< OVERFLOWCUT ) 
    { 
		DoParabolaFit(timeEstimate-1, amp, time);
		fNiter=0;
		fitDone=DoFit_1peak(first, nsamples, amp,  time, chi2);

		if(!fitDone)
			{
			amp=ampEstimate;
			time = timeEstimate;
			chi2=1.e9; 
			}

      time += timebinOffset;
      timeEstimate += timebinOffset;
      ndf = nsamples - 2;



    }
  }


  if ( fitDone ) 
  { 
    Float_t ampAsymm = (amp - ampEstimate)/(amp + ampEstimate);
    Float_t timeDiff = time - timeEstimate;
    
    if ( (TMath::Abs(ampAsymm) > 0.1) || (TMath::Abs(timeDiff) > 2) ) 
    {
      amp     = ampEstimate;
      time    = timeEstimate; 
      fitDone = kFALSE;
    } 
  }  
  if (amp >= fAmpCut ) 
  { 
    if ( ! fitDone) 
    { 
      amp += (0.5 - gRandom->Rndm()); 
    }
    time = time * TIMEBINWITH; 
    time -= fL1Phase;
    
    return AliCaloFitResults( -99, -99, fAlgo , amp, time,
                             (int)time, chi2, ndf, Ret::kDummy );
  }
  return AliCaloFitResults( Ret::kInvalid, Ret::kInvalid );
}



Bool_t AliCaloRawAnalyzerGamma2::DoFit_1peak(Int_t ifirst, Int_t nS, Float_t &A0, Float_t &t0, Float_t &chi2)
{
// fit using gamma-2 function 	(ORDER =2 assumed)

	if(nS<3) 
		return kFALSE;
	if(fNiter>fNiterationsMax) 
		return kFALSE;

	Double_t expi, ti, g_i, g_1i, g_2i, g1_1, g2_i, gp_i, g2p_i, q1_i, q2_i, D, dA, dt, delta;
	Double_t c11=0;
	Double_t c12=0;
	Double_t c21=0;
	Double_t c22=0;
	Double_t d1=0;
	Double_t d2=0;
	chi2=0;

	fNiter++;

	for(Int_t i=0; i<nS; i++)
	{

	ti=(i-t0)/TAU;
	if((ti+1)<0) continue;
	
	expi=TMath::Exp(-2*ti);
	g_2i=expi;	
	g_1i=(ti+1)*g_2i;	
	g_i =(ti+1)*g_1i;	
	gp_i=2*(g_i-g_1i);
	g2p_i=2*g_i*gp_i;
	q1_i=g_2i*(2*ti+1);
	q2_i=g_1i*g_1i*(4*ti+1);
	c11+=(GetReversed(i)-A0*2*g_i)*gp_i;
	c12+=g_i*g_i;
	c21+=GetReversed(i)*q1_i-A0*q2_i;
	c22+=g_i*g_1i;
	delta=A0*g_i-GetReversed(i);
	d1+=delta*g_i;
	d2+=delta*g_1i;
	chi2+=(delta*delta);
	}

	D=c11*c22-c12*c21;

	if(TMath::Abs(D)<1.e-9) 
		return kFALSE;

 	dt=(d1*c22-d2*c12)/D*TAU;
	dA=(d1*c21-d2*c11)/D;

 	t0+=dt;
	A0+=dA;

	Bool_t res=kTRUE;

	if(TMath::Abs(dA)>1 || TMath::Abs(dt)>0.01)
	res=DoFit_1peak(ifirst, nS, A0, t0, chi2);

return res;
}


void AliCaloRawAnalyzerGamma2::DoParabolaFit(Int_t x, Float_t &amp, Float_t &time)
{
// estimate initial values 
Double_t a=(GetReversed(x+2)+GetReversed(x)-2.*GetReversed(x+1))/2.;

if(TMath::Abs(a)<1.e09)
	{
	amp=GetReversed(x+1);
	time=x+1;
	return;
	}

Double_t b=GetReversed(x+1)-GetReversed(x)-a*(2.*x+1);
Double_t c=GetReversed(x)-b*x-a*x*x;

time=-b/2./a;
amp=a*time*time+b*time+c;

return ;
}


