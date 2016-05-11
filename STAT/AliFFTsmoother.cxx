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
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
/// \file AliFFTsmoother.cxx
/// \class AliFFTsmmother
///  FFT smoothing algorithm. Transform input histogram removing outlier frequency or smoothing frequencies
///  


#include "TH1.h"
#include "TGraph.h"
#include "TTreeStream.h"
#include "TVirtualFFT.h"
#include "TMath.h"
#include "TVector.h"
#include "TStatToolkit.h"
#include "AliFFTsmoother.h"


TGraph  * AliFFTsmoother::ReplaceOutlierFrequenciesMedian(TH1 *hinput, Double_t outlierCut, Int_t medianRange, Float_t smoothSigma,  Int_t lowBand, TTreeSRedirector * pcstream)
{
  /// 
  ///   FFT smoothing algorithm. Transform input histogram removing outlier frequency 
  ///   Code used for the analysis of the CE electron transparency measurement
  ///   authors: Marian + Mesut
  ///   Current algorithm is sensitive to frequency "aliasing" 
  ///       see detailed discussin in the UnitTest routines 
  ///       parasitic frequency and input histogram binning should be in "sync" 
  ///       e.g having repetetive structure 4 bins nx4 bins hsould be used
  /// 
  ///   Parameters:
  ///   hinput              - input histogram
  ///   outlierCut          - nsigma value abouve which the frequency is replaced |f_{i}-(f_{i+1)+f_{i-1})/2|<outlierCut*sigma
  ///   lowBand             - low freuency band (untoched)
  ///   pcstream            - ption to dump the intermedait results FFT into tree for further analysis
  ///
  /// Algorithm:
  ///  1.) Make FFT of the input histogram
  ///  2.) Replace outlier frequencies
  ///  3.) Make a back FFT  transformation
  ///  4.) In case specified streamer dump intermediat results for the checking purposes
  TGraph grInput(hinput);
  Int_t fftLength = hinput->GetNbinsX();  
  TH1D *hm = 0;
  TVirtualFFT::SetTransform(0);
  hm = (TH1D*)hinput->FFT(hm, "MAG");
  hm->SetTitle("Magnitude of the 1st transform");

  //
  //   1.) Make FFT of the input histogram
  //
  TVirtualFFT *fft = TVirtualFFT::GetCurrentTransform();
  TVectorD reFull(fftLength);
  TVectorD imFull(fftLength);
  TVectorD magFull(fftLength);
  TVectorD magS(fftLength);
  TVectorD phaseFull(fftLength);
  fft->GetPointsComplex(reFull.GetMatrixArray(),imFull.GetMatrixArray());
  for (Int_t ipoint=1; ipoint<fftLength/2-1; ipoint++){
    magFull[ipoint]=TMath::Sqrt(reFull[ipoint]*reFull[ipoint]+imFull[ipoint]*imFull[ipoint]);
    phaseFull[ipoint]=TMath::ATan2(imFull[ipoint], reFull[ipoint]);
  }
  //
  TVectorD reMod(fftLength);
  TVectorD imMod(fftLength);
  TVectorD vecDiff(fftLength);
  //
  //   2.) identify and replace  outlier frequencies
  // 
  TVectorD magMed(fftLength);
  for (Int_t ipoint=0; ipoint<fftLength/2-1; ipoint++){
    Int_t dedge=TMath::Min(ipoint, fftLength/2-ipoint-1);
    Int_t   window= TMath::Min(dedge, medianRange);
    magMed[ipoint]=TMath::Median(1+2*window, &(magFull.GetMatrixArray()[ipoint-window]));
    vecDiff[ipoint]=0;
    if (magMed[ipoint]>0){
      vecDiff[ipoint]=(magFull[ipoint]-magMed[ipoint]);
    }
  }  
  for (Int_t ipoint=1; ipoint<fftLength/2-1; ipoint++){
    vecDiff[ipoint]= (magFull[ipoint]-(magMed[ipoint-1]+magMed[ipoint+1])*0.5);
  }
  //
  Double_t meanT, rmsT;
  TStatToolkit::EvaluateUni(fftLength/2,vecDiff.GetMatrixArray(), meanT, rmsT,  0.95*(fftLength/2));  
  if (smoothSigma<=0){
    for (Int_t ipoint=0; ipoint<fftLength/2-1; ipoint++){
      reMod[ipoint]=reFull[ipoint];
      imMod[ipoint]=imFull[ipoint];
      if (ipoint<lowBand) continue;
      if (ipoint>=fftLength/2-lowBand) continue;
      if (TMath::Abs(vecDiff[ipoint]-meanT)>rmsT*outlierCut){
	reMod[ipoint]=magMed[ipoint]*TMath::Cos(phaseFull[ipoint]);
	imMod[ipoint]=magMed[ipoint]*TMath::Sin(phaseFull[ipoint]);
      }
    }
  }else{
    for (Int_t ipoint=0; ipoint<fftLength/2-1; ipoint++){
      reMod[ipoint]=reFull[ipoint];
      imMod[ipoint]=imFull[ipoint];
      if (ipoint<lowBand) continue;
      if (ipoint>=fftLength/2-lowBand) continue;
      //
      Double_t sumM=0, sumW=0;
      for (Int_t dpoint=-4.*smoothSigma; dpoint<4.*smoothSigma; dpoint++){
	if (dpoint+ipoint<0) continue;
	if (dpoint+ipoint>fftLength/2-1) continue;	
	if (TMath::Abs(vecDiff[ipoint+ipoint]-meanT)>rmsT*outlierCut){
	  continue;
	}
	Double_t w= TMath::Gaus(dpoint);
	sumM+=magFull[ipoint+dpoint]*w;
	sumW+=w;
      }
      sumM/=sumW;
      reMod[ipoint]=sumM*TMath::Cos(phaseFull[ipoint]);
      imMod[ipoint]=sumM*TMath::Sin(phaseFull[ipoint]);
    }    
  }
  //
  //   3.) Make a back FFT  transformation
  //
  TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &fftLength, "C2R M K");
  fft_back->SetPointsComplex(reMod.GetMatrixArray(),imMod.GetMatrixArray());
  fft_back->Transform();  
  TH1D *hb = 0;
  hb = (TH1D*)(TH1::TransformHisto(fft_back,hb,"Re"));
  hb->SetTitle("The backward transform result");
  hb->Scale(1./Double_t(fftLength));    
  TGraph *grOutput = new TGraph(grInput);
  for (Int_t i=0; i<fftLength; i++){
    //hinput->SetBinContent(i, hb->GetBinContent(i-firstBin));
    grOutput->GetY()[i]=hb->GetBinContent(i+1);
  }
  //
  //  4.) in case specified streamer dump intermediat results for the checking purposes
  //
  if (pcstream!=NULL){
    (*pcstream)<<"fft"<<
      "hinput."<<hinput<<         // input histogram
      "grInput.="<<&grInput<<     // input graph 
      "grOutput.="<<grOutput<<    // output graph with removed outlier frequencies
      //
      "reFull.="<<&reFull<<       // fft real part for original graph
      "imFull.="<<&imFull<<       // fft imaginary part for original graph
      "reMod.="<<&reMod<<
      "imMod.="<<&imMod<<
      //
      "magFull.="<<&magFull<<
      "magMed.="<<&magMed<<
      "vecDiff.="<<&vecDiff<<
      //
      "meanT="<<meanT<<          // mean difference
      "rmsT="<<rmsT<<            // rms difference
      "\n";
  }
  delete hb;
  return grOutput;  
}


