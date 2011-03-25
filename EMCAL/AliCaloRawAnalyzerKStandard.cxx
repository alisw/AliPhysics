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


#define BAD 4  //CRAP PTH
 
ClassImp( AliCaloRawAnalyzerKStandard )


AliCaloRawAnalyzerKStandard::AliCaloRawAnalyzerKStandard() : AliCaloRawAnalyzer("Chi Square ( kStandard )", "KStandard"),
						 fkEulerSquared(7.389056098930650227),
						 fTf1(0),
						 fTau(2.35),
						 fFixTau(kTRUE)
{
  
  fAlgo = Algo::kStandard;
  //comment
  for(int i=0; i < ALTROMAXSAMPLES; i++)
    {
      fXaxis[i] = i;
    }
  
  fTf1 = new TF1( "myformula", "[0]*((x - [1])/[2])^2*exp(-2*(x -[1])/[2])",  0, 30 ); 
  if (fFixTau) 
    {
    fTf1->FixParameter(2, fTau);
    }
  else 
    {
      fTf1->ReleaseParameter(2); // allow par. to vary
      fTf1->SetParameter(2, fTau);
    }
}


AliCaloRawAnalyzerKStandard::~AliCaloRawAnalyzerKStandard()
{
  delete fTf1;
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
      //Int_t id = fGeom->GetAbsCellIdFromCellIndexes(in.GetModule(), in.GetRow(), in.GetColumn()) ;
      // lowGain  = in.IsLowGain();
     
      time = time * TIMEBINWITH; 
      
      /////////////!!!!!!!!!time -= in.GetL1Phase();

      time -= fL1Phase;

      //    AliDebug(2,Form("id %d lowGain %d amp %g", id, lowGain, amp));
      //  AddDigit(digitsArr, id, lowGain, amp, time, chi2, ndf); 
      

      return AliCaloFitResults( -99, -99, fAlgo , amp, time,
				time, chi2, ndf, Ret::kDummy );
      
      
      //	AliCaloFitSubarray(index, maxrev, first, last));
    
    }
  
  
  return AliCaloFitResults( Ret::kInvalid, Ret::kInvalid );
}




/*  
   return AliCaloFitResults( maxamp, ped, Ret::kCrude, maxf, timebinOffset,
					timebinOffset, chi2, ndf, Ret::kDummy, AliCaloFitSubarray(index, maxrev, first, last) ); 
	    }
        } // ampcut
    }
  return AliCaloFitResults( Ret::kInvalid, Ret::kInvalid );
*/


/*
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
      if(  maxf < fAmpCut  ||  ( maxamp - ped) > fOverflowCut  ) // (maxamp - ped) > fOverflowCut = Close to saturation (use low gain then)
	{
	  return  AliCaloFitResults( maxamp, ped, Ret::kCrude, maxf, timebinOffset);
 	}            
      else if ( maxf >= fAmpCut )
	{
	  int first = 0;
	  int last = 0;
	  SelectSubarray( fReversed,  bunchvector.at(index).GetLength(),  maxrev, &first, &last, fFitArrayCut);
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
		AliError( Form("TGraph Fit exception %s", e.what()) ); 
		return AliCaloFitResults( maxamp, ped, Ret::kNoFit, maxf, timebinOffset,
					  timebinOffset, Ret::kDummy, Ret::kDummy, Ret::kDummy, AliCaloFitSubarray(index, maxrev, first, last) );
	      }

	      if( fVerbose == true )
		{
		  AliCaloRawAnalyzer::PrintBunch( bunchvector.at(index) ); 
		  PrintFitResult( fTf1 ) ;
		}  
	      // global tmax
	      tmax = fTf1->GetParameter(1) + timebinOffset - (maxrev - first) // abs. t0
		+ fTf1->GetParameter(2); // +tau, makes sum tmax
	      
	        delete graph;
		return AliCaloFitResults( maxamp, ped , Ret::kFitPar,
					  fTf1->GetParameter(0)/fkEulerSquared, 
					  tmax,
					  timebinOffset,  
					  fTf1->GetChisquare(), 
					  fTf1->GetNDF(),
					  Ret::kDummy, AliCaloFitSubarray(index, maxrev, first, last) );
				
		//     delete graph;
	
	    }
	  else
	    {
	      
	      Float_t chi2 = CalculateChi2(maxf, maxrev, first, last);
	      Int_t ndf = last - first - 1; // nsamples - 2
	      return AliCaloFitResults( maxamp, ped, Ret::kCrude, maxf, timebinOffset,
					timebinOffset, chi2, ndf, Ret::kDummy, AliCaloFitSubarray(index, maxrev, first, last) ); 
	    }
        } // ampcut
    }
  return AliCaloFitResults( Ret::kInvalid, Ret::kInvalid );
  
}
*/



void 
AliCaloRawAnalyzerKStandard::PrintFitResult(const TF1 *f) const
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




	
//____________________________________________________________________________ 
void
 AliCaloRawAnalyzerKStandard::FitRaw(const Int_t firstTimeBin, const Int_t lastTimeBin, Float_t & amp, Float_t & time, Float_t & chi2, Bool_t & fitDone) const 
{ // Fits the raw signal time distribution
  
  //--------------------------------------------------
  //Do the fit, different fitting algorithms available
  //--------------------------------------------------
  
  // fprintf(fp, "%s:%d:%s\n", __FILE__, __LINE__, __FUNCTION__ );

  int nsamples = lastTimeBin - firstTimeBin + 1;
  fitDone = kFALSE;
  
  // switch(fFittingAlgorithm) 
  //   {
  //  case Algo::kStandard:
  //    {
  if (nsamples < 3) { return; } // nothing much to fit
  //printf("Standard fitter \n");

  // Create Graph to hold data we will fit 
  
  TGraph *gSig =  new TGraph( nsamples); 
 
  for (int i=0; i<nsamples; i++) 
    {
      Int_t timebin = firstTimeBin + i;    
      gSig->SetPoint(i, timebin, GetReversed(timebin)); 
    }
      
  TF1 * signalF = new TF1("signal", RawResponseFunction, 0, TIMEBINS , 5);
  signalF->SetParameters(10.,5., TAU  ,ORDER,0.); //set all defaults once, just to be safe
  signalF->SetParNames("amp","t0","tau","N","ped");
  signalF->FixParameter(2,TAU); // tau in units of time bin
  signalF->FixParameter(3,ORDER); // order
  signalF->FixParameter(4, 0); // pedestal should be subtracted when we get here 
  signalF->SetParameter(1, time);
  signalF->SetParameter(0, amp);
  // set rather loose parameter limits
  signalF->SetParLimits(0, 0.5*amp, 2*amp );
  signalF->SetParLimits(1, time - 4, time + 4); 
      
  try {			
    gSig->Fit(signalF, "QROW"); // Note option 'W': equal errors on all points
    // assign fit results
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
      
  //printf("Std   : Amp %f, time %g\n",amp, time);
  delete gSig; // delete TGraph
      
  //	break;
  //    }//kStandard Fitter
  //----------------------------
 
  /*
    case Algo::kLogFit:
    {
    if (nsamples < 3) { return; } // nothing much to fit
    //printf("LogFit \n");
      
    // Create Graph to hold data we will fit 
    TGraph *gSigLog =  new TGraph( nsamples); 
    for (int i=0; i<nsamples; i++) {
    Int_t timebin = firstTimeBin + i;    
    gSigLog->SetPoint(timebin, timebin, TMath::Log(fRawAnalyzer->GetReversed(timebin) ) ); 
    }
      
    TF1 * signalFLog = new TF1("signalLog", RawResponseFunctionLog, 0,  TIMEBINS , 5);
    signalFLog->SetParameters(2.3, 5.,TAU,ORDER,0.); //set all defaults once, just to be safe
    signalFLog->SetParNames("amplog","t0","tau","N","ped");
    signalFLog->FixParameter(2,TAU); // tau in units of time bin
    signalFLog->FixParameter(3, ORDER); // order
    signalFLog->FixParameter(4, 0); // pedestal should be subtracted when we get here 
    signalFLog->SetParameter(1, time);
    if (amp>=1) {
    signalFLog->SetParameter(0, TMath::Log(amp));
    }
      
    gSigLog->Fit(signalFLog, "QROW"); // Note option 'W': equal errors on all points
      
    // assign fit results
    Double_t amplog = signalFLog->GetParameter(0); //Not Amp, but Log of Amp
    amp = TMath::Exp(amplog);
    time = signalFLog->GetParameter(1);
    fitDone = kTRUE;
      
    delete signalFLog;
    //printf("LogFit: Amp %f, time %g\n",amp, time);
    delete gSigLog; 
    break;
    } //kLogFit 
      //----------------------------	
      //----------------------------
      }//switch fitting algorithms
  */
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
    amp = ymMaxY[0] ; 
  //printf("Yves   : Amp %f, time %g\n",amp, time);
  //END YS
  return;
}



//__________________________________________________________________
Double_t 
AliCaloRawAnalyzerKStandard::RawResponseFunction(Double_t *x, Double_t *par)
{
  // Matches version used in 2007 beam test
  //
  // Shape of the electronics raw reponse:
  // It is a semi-gaussian, 2nd order Gamma function of the general form
  //
  // xx = (t - t0 + tau) / tau  [xx is just a convenient help variable]
  // F = A * (xx**N * exp( N * ( 1 - xx) )   for xx >= 0
  // F = 0                                   for xx < 0 
  //
  // parameters:
  // A:   par[0]   // Amplitude = peak value
  // t0:  par[1]
  // tau: par[2]
  // N:   par[3]
  // ped: par[4]
  //
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

