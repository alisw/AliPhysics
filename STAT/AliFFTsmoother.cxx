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


TGraph  * AliFFTsmoother::ReplaceOutlierFrequencies(TH1 *hinput, Int_t firstBin, Int_t lastBin, Double_t outlierCut, Int_t skipFreq, TTreeSRedirector * pcstream)
{
  /// 
  ///   FFT smoothing algorithm. Transform input histogram removing outlier frequency 
  ///   Code used for the analysis of the CE electron transparency measurement
  ///   authors: Marian + Mesut
  /// 
  ///   Parameters:
  ///   hinput              - input histogram
  ///   firstBin, lastBin   - range of the bins which will be touched
  ///   outlierCut          - nsigma value abouve which the frequency is replaced |f_{i}-(f_{i+1)+f_{i-1})/2|<outlierCut*sigma
  ///   skip freq           - option to skip low frequency (better to sue range?)
  ///   pcstream            - ption to dump the intermedait results FFT into tree for further analysis
  ///
  /// Algorithm:
  ///  1.) Make FFT of the input histogram
  ///  2.) Replace outlier frequencies
  ///  3.) Make a back FFT  transformation
  ///  4.) In case specified streamer dump intermediat results for the checking purposes
  TGraph grInput(hinput);
  Int_t fftLength = lastBin-firstBin-1;  
  TH1D *hsin = new TH1D("hsin", "hsin", fftLength, 0., fftLength);
  for (Int_t i=firstBin; i<lastBin; i++){
    hsin->SetBinContent(i-firstBin, hinput->GetBinContent(i));
  }
  //
  //   1.) Make FFT of the input histogram
  //
  TVirtualFFT *fft = TVirtualFFT::GetCurrentTransform();
  TVectorD reFull(fftLength);
  TVectorD imFull(fftLength);
  fft->GetPointsComplex(reFull.GetMatrixArray(),imFull.GetMatrixArray());
  TVectorD reMod(fftLength);
  TVectorD imMod(fftLength);
  TVectorD vecDiff(fftLength);
  TVectorD vecDiffL(fftLength);
  TVectorD vecDiffR(fftLength);
  //
  //   2.) Replace outlier frequencies
  // 
  for (Int_t ipoint=1; ipoint<fftLength/2-1; ipoint++){
    Double_t diffRe=reFull[ipoint]-(reFull[ipoint-1]+reFull[ipoint+1])*0.5;
    Double_t diffIm=imFull[ipoint]-(imFull[ipoint-1]+imFull[ipoint+1])*0.5;
    Double_t diff=TMath::Sqrt(diffRe*diffRe+diffIm*diffIm);
    vecDiff[ipoint]=diff;
    vecDiffL[ipoint]=TMath::Sqrt((reFull[ipoint]-reFull[ipoint-1])*(reFull[ipoint]-reFull[ipoint-1])+(reFull[ipoint]-imFull[ipoint-1])*(reFull[ipoint]-imFull[ipoint-1]));
    vecDiffR[ipoint]=TMath::Sqrt((reFull[ipoint]-reFull[ipoint+1])*(reFull[ipoint]-reFull[ipoint+1])+(reFull[ipoint]-imFull[ipoint+1])*(reFull[ipoint]-imFull[ipoint+1]));
  }
  Double_t mean90, rms90;
  TStatToolkit::EvaluateUni(fftLength/2,vecDiff.GetMatrixArray(), mean90, rms90,  0.9*(fftLength/2));
  for (Int_t ipoint=0; ipoint<fftLength/2; ipoint++){
    reMod[ipoint]=reFull[ipoint];
    imMod[ipoint]=imFull[ipoint];
    if (ipoint<skipFreq) continue;
    if (ipoint>=fftLength/2-skipFreq) continue;
    if (vecDiff[ipoint]>mean90+rms90*outlierCut){
      if (vecDiff[ipoint-1]<vecDiff[ipoint+1]){
	reMod[ipoint]=reMod[ipoint-1];
	imMod[ipoint]=imMod[ipoint-1];
      }else{
	reMod[ipoint]=reMod[ipoint+1];
	imMod[ipoint]=imMod[ipoint+1];
      }
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
  for (Int_t i=firstBin+1; i<lastBin; i++){
    hinput->SetBinContent(i, hb->GetBinContent(i-firstBin));
  }
  TGraph *grOutput = new TGraph(hinput);
  //
  //  4.) in case specified streamer dump intermediat results for the checking purposes
  //
  if (pcstream!=NULL){
    (*pcstream)<<"fft"<<
      "hinput."<<hinput<<         // input histogram
      "grInput.="<<&grInput<<     // input graph 
      "grOutput.="<<grOutput<<    // output graph with removed outlier frequencies
      "reFull.="<<&reFull<<       // fft real part for original graph
      "imFull.="<<&imFull<<       // fft imaginary part for original graph
      "reMod.="<<&reMod<<
      "imMod.="<<&imMod<<
      "\n";
  }
  delete hb;
  delete hsin;
  return grOutput;  
}

TGraph  *  AliFFTsmoother::SmoothFrequencies(TH1 */*hinput*/, Int_t /*firstBin*/, Int_t /*lastBin*/, Double_t /*outlierCut*/, Int_t /*skipFreq*/, TTreeSRedirector * /*pcstream*/){
  ///
  /// Algortihm used to remove outlier frequncies and to smooth FFT freqiencies for the ion tail analysis
  ///     (MI, Mesut)
  ///     to be ported from the compiled macro after cleaning
  /// Significantly slower algorithm compared to one above because of suage of robust fitter for frequency smoothing
  /// 
  return 0;
}
