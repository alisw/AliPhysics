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

#include "AliCaloRawAnalyzerKStandard.h"
#include "AliCaloBunchInfo.h"
#include "AliCaloFitResults.h"
#include "AliLog.h"
#include "TMath.h"
#include <stdexcept>
#include <iostream>
#include "TF1.h"
#include "TGraph.h"
#include "TRandom.h"

using namespace std;

ClassImp( AliCaloRawAnalyzerKStandard )


AliCaloRawAnalyzerKStandard::AliCaloRawAnalyzerKStandard() : AliCaloRawAnalyzerFitter("Chi Square ( kStandard )", "KStandard")
{
  
  fAlgo = Algo::kStandard;
}


AliCaloRawAnalyzerKStandard::~AliCaloRawAnalyzerKStandard()
{
  //  delete fTf1;
}


AliCaloFitResults
AliCaloRawAnalyzerKStandard::Evaluate( const vector<AliCaloBunchInfo>  &bunchlist, const UInt_t altrocfg1,  const UInt_t altrocfg2 )
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
					maxADC, timeEstimate, pedEstimate, first, last,   fAmpCut ); 
  
  
  if (ampEstimate >= fAmpCut  ) 
    { 
      time = timeEstimate; 
      Int_t timebinOffset = bunchlist.at(bunchIndex).GetStartBin() - (bunchlist.at(bunchIndex).GetLength()-1); 
      amp = ampEstimate; 
      
      if ( nsamples > 1 && maxADC< OVERFLOWCUT ) 
	{ 
	  FitRaw(first, last, amp, time, chi2, fitDone);
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
				time, chi2, ndf, Ret::kDummy );
     }
  return AliCaloFitResults( Ret::kInvalid, Ret::kInvalid );
}

	
//____________________________________________________________________________ 
void
 AliCaloRawAnalyzerKStandard::FitRaw(const Int_t firstTimeBin, const Int_t lastTimeBin, Float_t & amp, Float_t & time, Float_t & chi2, Bool_t & fitDone) const 
{ 
  // Fits the raw signal time distribution
  int nsamples = lastTimeBin - firstTimeBin + 1;
  fitDone = kFALSE;
  if (nsamples < 3) { return; } 
  
  TGraph *gSig =  new TGraph( nsamples); 
 
  for (int i=0; i<nsamples; i++) 
    {
      Int_t timebin = firstTimeBin + i;    
      gSig->SetPoint(i, timebin, GetReversed(timebin)); 
    }
      
  TF1 * signalF = new TF1("signal", RawResponseFunction, 0, TIMEBINS , 5);
  signalF->SetParameters(10.,5., TAU  ,ORDER,0.); //set all defaults once, just to be safe
  signalF->SetParNames("amp","t0","tau","N","ped");
  signalF->FixParameter(2,TAU); 
  signalF->FixParameter(3,ORDER); 
  signalF->FixParameter(4, 0); 
  signalF->SetParameter(1, time);
  signalF->SetParameter(0, amp);
  signalF->SetParLimits(0, 0.5*amp, 2*amp );
  signalF->SetParLimits(1, time - 4, time + 4); 
      
  try {			
    gSig->Fit(signalF, "QROW"); // Note option 'W': equal errors on all points
    amp  = signalF->GetParameter(0); 
    time = signalF->GetParameter(1);
    chi2 = signalF->GetChisquare();
    fitDone = kTRUE;
  }
  catch (const std::exception & e) {
    AliError( Form("TGraph Fit exception %s", e.what()) ); 
    // stay with default amp and time in case of exception, i.e. no special action required
    fitDone = kFALSE;
  }

  delete signalF;
  delete gSig; // delete TGraph
  return;
}


//__________________________________________________________________
void 
AliCaloRawAnalyzerKStandard::FitParabola(const TGraph *gSig, Float_t & amp) const 
{
  //BEG YS alternative methods to calculate the amplitude
  Double_t * ymx = gSig->GetX() ; 
  Double_t * ymy = gSig->GetY() ; 
  const Int_t kN = 3 ; 
  Double_t ymMaxX[kN] = {0., 0., 0.} ; 
  Double_t ymMaxY[kN] = {0., 0., 0.} ; 
  Double_t ymax = 0. ; 
  // find the maximum amplitude
  Int_t ymiMax = 0 ;  
  for (Int_t ymi = 0; ymi < gSig->GetN(); ymi++) {
    if (ymy[ymi] > ymMaxY[0] ) {
      ymMaxY[0] = ymy[ymi] ; //<========== This is the maximum amplitude
      ymMaxX[0] = ymx[ymi] ;
      ymiMax = ymi ; 
    }
  }
  // find the maximum by fitting a parabola through the max and the two adjacent samples
  if ( ymiMax < gSig->GetN()-1 && ymiMax > 0) {
    ymMaxY[1] = ymy[ymiMax+1] ;
    ymMaxY[2] = ymy[ymiMax-1] ; 
    ymMaxX[1] = ymx[ymiMax+1] ;
    ymMaxX[2] = ymx[ymiMax-1] ; 
    if (ymMaxY[0]*ymMaxY[1]*ymMaxY[2] > 0) {
      //fit a parabola through the 3 points y= a+bx+x*x*x
      Double_t sy = 0 ; 
      Double_t sx = 0 ; 
      Double_t sx2 = 0 ; 
      Double_t sx3 = 0 ; 
      Double_t sx4 = 0 ; 
      Double_t sxy = 0 ; 
      Double_t sx2y = 0 ; 
      for (Int_t i = 0; i < kN ; i++) {
        sy += ymMaxY[i] ; 
        sx += ymMaxX[i] ; 		
        sx2 += ymMaxX[i]*ymMaxX[i] ; 
        sx3 += ymMaxX[i]*ymMaxX[i]*ymMaxX[i] ; 
        sx4 += ymMaxX[i]*ymMaxX[i]*ymMaxX[i]*ymMaxX[i] ; 
        sxy += ymMaxX[i]*ymMaxY[i] ; 
        sx2y += ymMaxX[i]*ymMaxX[i]*ymMaxY[i] ; 
      }
      Double_t cN = (sx2y*kN-sy*sx2)*(sx3*sx-sx2*sx2)-(sx2y*sx-sxy*sx2)*(sx3*kN-sx*sx2); 
      Double_t cD = (sx4*kN-sx2*sx2)*(sx3*sx-sx2*sx2)-(sx4*sx-sx3*sx2)*(sx3*kN-sx*sx2) ;
      Double_t c  = cN / cD ; 
      Double_t b  = ((sx2y*kN-sy*sx2)-c*(sx4*kN-sx2*sx2))/(sx3*kN-sx*sx2) ;
      Double_t a  = (sy-b*sx-c*sx2)/kN  ;
      Double_t xmax = -b/(2*c) ; 
      ymax = a + b*xmax + c*xmax*xmax ;//<========== This is the maximum amplitude
      amp = ymax;
    }
  }
  
  Double_t diff = TMath::Abs(1-ymMaxY[0]/amp) ; 
  if (diff > 0.1) 
    {
      amp = ymMaxY[0] ; 
    }

  return;
}


//__________________________________________________________________
Double_t 
AliCaloRawAnalyzerKStandard::RawResponseFunction(Double_t *x, Double_t *par)
{
  Double_t signal = 0.;
  Double_t tau    = par[2];
  Double_t n      = par[3];
  Double_t ped    = par[4];
  Double_t xx     = ( x[0] - par[1] + tau ) / tau ;

  if (xx <= 0) 
    signal = ped ;  
  else {  
    signal = ped + par[0] * TMath::Power(xx , n) * TMath::Exp(n * (1 - xx )) ; 
  }
  return signal ;  
}

