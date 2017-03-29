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

#include "AliAnalysisMuMuSpectraCapsulePbPb.h"

ClassImp(AliAnalysisMuMuSpectraCapsulePbPb)

#include "TF1.h"
#include "TProfile.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "THashList.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TPaveText.h"
#include "TLine.h"
#include "TList.h"
#include "TMath.h"
#include "TObjArray.h"
#include "AliAnalysisMuMuBinning.h"
#include "AliLog.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TColor.h"
#include "TROOT.h"
#include "AliCounterCollection.h"
#include "AliMergeableCollection.h"
#include "AliAnalysisMuMuSpectra.h"
#include "AliAnalysisMuMuResult.h"
#include "AliAnalysisMuMuJpsiResult.h"
#include "AliAnalysisMuMuSpectraCapsulePbPb.h"
#include <fstream>
#include <string>

using std::cout;
using std::endl;
using std::ifstream;


namespace
{
  const Double_t BR                = 5.96/100; // Branching ratio
  const Double_t BRerr             = 0.03/5.96; // Branching ratio
  //Normalization factor
  const Double_t Fnorm             = 11.842;    // Normalization
  const Double_t FnormStat         = 0.00095;     // Normalization
  const Double_t FnormSyst         = 0.059;     // Normalization

  //pp Cross-section integrated for pp@5TeV 0<pt<8 , -4<y<-2.5
  // const Double_t sigmaPP        = 5.55;    // for fully integrated case
  // const Double_t dsigmaPP       = 0.08;    // idem
  // const Double_t dsigmaPPCorr   = 0.28;    // for fully integrated case

  // pp Cross-section integrated for pp@5TeV 0.3<pT<8 , -4<y<-2.5
  // const Double_t sigmaPP        = 5.47;    // for fully integrated case
  // const Double_t dsigmaPP       = 0.08;    // idem
  // const Double_t dsigmaPPCorr   = 0.27;    // for fully integrated case

  //pp Cross-section integrated for pp@5TeV 0<pT<12 , -4<y<-2.5
  const Double_t sigmaPP           = 5.61;    // for fully integrated case
  const Double_t dsigmaPP          = 0.08;    // idem
  const Double_t dsigmaPPCorr      = 0.27;    // for fully integrated case

  // Global MC sys. err. for centrality integrated in pt and Y
  const Double_t MCParamError      = 2.0;//%
  // NofMUL correspondind to signal extraction from 2015
  const Double_t Mul2015           = 126778700;
  // Syst. associated to matching
  const Double_t MatchingError     = 1.;//%
  // Syst. associated to Intrinsic trigger efficiency
  const Double_t TriggerError      = 1.5;//%
  // Syst. associated to Intrinsic trigger efficiency
  const Double_t TrackingErrorCent = 1.;//%
  // Syst. associated to Intrinsic trigger efficiency
  const Double_t TriggerErrorCent  = 1.;//%


}


//_____________________________________________________________________________
 AliAnalysisMuMuSpectraCapsulePbPb::AliAnalysisMuMuSpectraCapsulePbPb(
const AliAnalysisMuMuSpectra*  spectra,
const TString                 spectraPath,
const char                  * externFile,
const char                  * externFile2)
:
  AliAnalysisMuMuSpectraCapsule(spectra,spectraPath,externFile,externFile2)
  // fSpectra(spectra),
  // fSpectraName(spectraPath),
  // fExternFile(externFile),
  // fExternFile2(externFile2),
  // fPrintFlag(kFALSE)
{
  // //Check point
  // if (!fSpectra)
  // {
  //   AliError(Form("Cannot find spectra wih name %s Please check the name",fSpectra->GetName()));
  //   return;
  // }
  // AliDebug(1, Form(" - spectra(%s) = %p ",fSpectra->GetName(),fSpectra));


  // if (fSpectraName.IsNull())
  // {
  //   AliWarning(Form("No spectra name ! "));
  //   return;
  // }

  if(!AliAnalysisMuMuSpectraCapsule::SetConstantFromExternFile(fExternFile2,&fConstArray[0],&fSpectraName))
  {
    AliWarning(Form("No extern file readed"));
  }

}

//_____________________________________________________________________________
AliAnalysisMuMuSpectraCapsulePbPb::~AliAnalysisMuMuSpectraCapsulePbPb()
{
  // dtor
}

//_____________________________________________________________________________
TGraphErrors* AliAnalysisMuMuSpectraCapsulePbPb::ComputeYield( const char* what, const TH1* histo, const char* sResName, Double_t MUL)
{
  /// @brief Compute Yield.
  /// @argument  what  the yield nominator, i.e NofJPsi, meanPT etc. (null by default)
  /// @argument  histo histogramme of Equivalent MinBias

  if(!GetSpectra() || strcmp(what,"")==1) return 0x0;

  printf("here !\n");

  // Some constants
  const TString graphTitle = Form("%s-YIELD",GetSpectraName().Data());
  TString sres(sResName);

  // Pointers to handle results and subresults and binning
  AliAnalysisMuMuResult    * result;
  AliAnalysisMuMuBinning   ::Range* r;

 // Array to store bins for the while loop
  TObjArray * bins=GetSpectra()->Binning()->CreateBinObjArray();// (intrinseque 'new')
  if (!bins)
  {
    AliError(Form("Cannot find bins"));
    return 0x0;
  }

  // Here we define some pointers
  TGraphErrors* graphAll[2];
  TGraphErrors*graph(0x0);
  TGraphErrors*graph_sysUncorr(0x0);

  Double_t    * binArray(0x0) ;// (intrinseque 'new')
  Int_t binsX = 0;

  //________Define histo according to bin type
  if (GetSpectraName().Contains("-INTEGRATED"))
  {
    graph           = new TGraphErrors(1);
    graph_sysUncorr = new TGraphErrors(1);
    graph->SetTitle(graphTitle.Data());

    // graph_sysUncorr->SetFillColorAlpha(5,0.05);
  }
  else if (GetSpectraName().Contains("-PT") || GetSpectraName().Contains("-Y"))
  {
    binArray = GetSpectra()->Binning()->CreateBinArray();
    binsX    = GetSpectra()->Binning()->GetNBinsX();

    if (!binArray)
    {
      AliError(Form("Cannot set binArray"));
      return 0x0;
    }
    if (binsX==0)
    {
      AliError(Form("Cannot set binsX"));
      return 0x0;
    }

    graph           = new TGraphErrors(binsX);
    // graph_sysUncorr = new TGraphErrors(binsX);
    graph->SetTitle(graphTitle.Data());

    // graph_sysUncorr->SetFillColorAlpha(5,0.05);
  }
  else
  {
    cout << "Unknowned Bin type !" << endl;
    return 0x0;
  }
  //________

  //________Counters and Iterator for bin
  Int_t nofResult = 0;
  TIter nextBin(bins);
  nextBin.Reset();
  //________

  // Loop on bins
  //==============================================================================
  while ((r = static_cast<AliAnalysisMuMuBinning::Range*>(nextBin())))
  {
    //________Make bin a MuMuResult
    result = GetSpectra()->GetResultForBin(*r);
    if (!result)
    {
      AliError(Form("Cannot find result "));
      return 0x0;
    }
    AliDebug(1, Form("result(%s) = %p ",result->GetName(),result));
    //________

    // Store quantities
    Double_t NofWhattTot = result->GetValue(what,sres.Data());
    Double_t NofWhattTotError = result->GetErrorStat(what,sres.Data());

    Double_t nEqMBTot      = 0.0;
    Double_t nEqMBTotError = 0.0;
    if (histo)
    {
      nEqMBTot = histo->GetBinContent(nofResult+1);
      nEqMBTotError = histo->GetBinError(nofResult+1);
      AliDebug(1,Form("histo used    : %s",histo->GetTitle()));
    } else {
      Double_t nmul = MUL!=0. ? MUL :  Mul2015;
      nEqMBTot      = Fnorm * nmul;
      nEqMBTotError = nmul*FnormStat/Fnorm ;
    }

    if( NofWhattTot==0||NofWhattTotError==0||nEqMBTot==0||nEqMBTotError==0)
    {
      AliError("Cannot set quantities properly");
      return 0x0;
    }

    AliDebug(1,Form("%s            = %f",what,NofWhattTot));
    AliDebug(1,Form("%s error      = %f",what,NofWhattTotError));
    AliDebug(1,Form("nEqMBTot      = %f",nEqMBTot));
    AliDebug(1,Form("nEqMBTotError = %f",nEqMBTotError));

    //________Compute yield in case of fully integrated spectra
    if(GetSpectraName().Contains("-INTEGRATED"))
    {

      Double_t yieldInt = NofWhattTot/(nEqMBTot*BR*fConstArray[8]);
      Double_t yieldIntError = yieldInt*AliAnalysisMuMuResult::ErrorAB(NofWhattTot,NofWhattTotError,nEqMBTot,TMath::Sqrt(nEqMBTot));

      // Add results to TGraphs
      graph->SetPoint(nofResult,fConstArray[0],yieldInt);
      graph->SetPointError(nofResult,fConstArray[1],yieldIntError);
    } else {

      // read exterfile and get the correct value
      float valueArray[13];
      //  valueArray[0], valueArray[1], valueArray[2], valueArray[3],     valueArray[4], valueArray[5], valueArray[6], valueArray[7], valueArray[8], valueArray[9], valueArray[10], valueArray[11] valueArray[12]
      //  sigmapp         dsigmapp      dsigmappUnCorr   dsigmappCorr(%)  AccEff         dAccEff        sysMC          TrajEffError    TriggerError  MatchingError,  NofJpsi,       NofJpsiStat,    NofJpsiSys

      if(ReadFromFile(r->AsString(),&valueArray[0])==kFALSE) return 0x0;
      AliDebug(1, " Values correctly read from extern file");

      Double_t yieldInt = NofWhattTot/(nEqMBTot*valueArray[4]);
      Double_t yieldIntError = yieldInt*AliAnalysisMuMuResult::ErrorAB(NofWhattTot,NofWhattTotError,nEqMBTot,TMath::Sqrt(nEqMBTot));

      printf("bin : %s -- %s = %f +/- %f --  nEqMBTot : %f -- yield : %f +/- %f\n",r->AsString().Data(),what,NofWhattTot,NofWhattTotError,nEqMBTot,yieldInt,yieldIntError );

      //Fill graph
      Double_t binCenter = (binArray[nofResult+1]-binArray[nofResult])/2 + binArray[nofResult];
      graph->SetPoint(nofResult,binCenter,yieldInt);
      graph->SetPointError(nofResult,r->WidthX()/5,yieldIntError);
    }
    //________

    nofResult++;
  }

  // Config. graphics
  if(GetSpectraName().Contains("-INTEGRATED"))graph->GetXaxis()->SetTitle(Form("INTEGRATED"));
  else if (GetSpectraName().Contains("-PT"))graph->GetXaxis()->SetTitle(Form("PT"));
  else if (GetSpectraName().Contains("-Y"))graph->GetXaxis()->SetTitle(Form("Y"));
  graph->GetYaxis()->SetTitle("Yield");
  graph->SetMarkerColor(4);
  graph->SetMarkerStyle(21);

  delete bins;
  delete binArray;

 return graph ;
}

//_____________________________________________________________________________
void AliAnalysisMuMuSpectraCapsulePbPb::DrawResults( const char* what, const char* particle,const char* subresults) const
{
  /**
   *
   * Print fit results on a single canvas
   *
   */

  //Check point
  if(!GetSpectra() ) return ;

  // Pointers to handle results and subresults
  AliAnalysisMuMuResult    * result;
  AliAnalysisMuMuJpsiResult* subresult;
  AliAnalysisMuMuResult    * sr;
  AliAnalysisMuMuBinning   ::Range* r;
  TH1 * h = 0x0;

  const TString sres(subresults);

  //Pointer for functions
  TF1* f1 = 0x0;
  TF1* f2 = 0x0;
  TF1* f3 = 0x0;
  TF1* f4 = 0x0;

  // Arrays
  TObjArray* histos = new TObjArray(0x0);// Array to store histograms
  TObjArray* bins=GetSpectra()->Binning()->CreateBinObjArray();// Array to store bins

  if (!bins){
    AliError(Form("Cannot find bins"));
    return;
  }

  // Settings for histo
  Double_t xmin(-1);
  Double_t xmax(-1);

  if ( fSpectraName.Contains(particle)){
    xmin = 2;
    xmax = 6;
  }

  //Iterator for bin
  TIter nextBin(bins);

  std::vector<double> NofWhat;
  std::vector<double> NofWhatErr;
  std::vector<double> SoverB;
  std::vector<double> BinRange;
  std::vector<double> massJpsi;
  std::vector<double> massJpsiErr;
  std::vector<double> sigmaJpsi;
  std::vector<double> sigmaJpsiErr;
  std::vector<double> chi2perndf;

  TString BinType;

  // --- Loop on bins
  while ((r = static_cast<AliAnalysisMuMuBinning::Range*>(nextBin())))
  {
    // Make bin a MuMuResult
    BinType =Form("%s",GetSpectra()->GetName());

    result = GetSpectra()->GetResultForBin(*r);
    if (!result)
    {
      AliError(Form("Cannot find result "));
      return;
    }
    AliDebug(1, Form("result(%s) = %p ",result->GetName(),result));

    TIter nextSubResult(result->SubResults());// Iterator for subresults
    nextSubResult.Reset();
    // Loop on subresults
    //==============================================================================
    while ((sr = static_cast<AliAnalysisMuMuResult*>(nextSubResult())))
    {
      // Get our final result
      subresult = static_cast<AliAnalysisMuMuJpsiResult*>(result->SubResult(Form("%s",sr->GetName())));
      if (!subresult){
        AliError(Form("Cannot find subresult "));
        return;
      }
      AliDebug(1,Form("subresult(%s) = %p",sr->GetName(),subresult));
      if(!sres.IsNull() && !sres.Contains(sr->GetName())) continue;

      // Get histo
      if ( subresult ) h = (TH1*)subresult->Histo();
      AliDebug(1, Form(" - Histo(%s) = %p ",h->GetTitle(),h));

      // Store it
      if(h) histos->Add(h);
      else {
        AliError(Form("Cannot set histo result "));
        return;
      }

      // To go on the TLegend
      if (subresult->HasValue(what))             NofWhat.push_back( subresult->GetValue(what) );
      else NofWhat.push_back( 0.0 );

      if (subresult->HasValue(what))             NofWhatErr.push_back( subresult->GetErrorStat(what) );
      else NofWhatErr.push_back( 0.0 );

      if(subresult->HasValue("SignalOverBkg3s")) SoverB.push_back ( subresult->GetValue("SignalOverBkg3s") );
      else SoverB.push_back (0.0);

      if(subresult->HasValue("mJPsi"))           massJpsi.push_back( subresult->GetValue("mJPsi") );
      else massJpsi.push_back(0.0);

      if(subresult->HasValue("mJPsi"))           massJpsiErr.push_back( subresult->GetErrorStat("mJPsi") );
      else massJpsiErr.push_back(0.0);

      if(subresult->HasValue("sJPsi"))           sigmaJpsi.push_back( subresult->GetValue("sJPsi") );
      else sigmaJpsi.push_back(0.0);

      if(subresult->HasValue("sJPsi"))           sigmaJpsiErr.push_back( subresult->GetErrorStat("sJPsi") );
      else sigmaJpsiErr.push_back(0.0);

      if(subresult->HasValue("FitChi2PerNDF"))   chi2perndf.push_back( subresult->GetErrorStat("FitChi2PerNDF") );
      else chi2perndf.push_back(0.0);

      BinRange.push_back ( r->Xmin() );
      BinRange.push_back ( r->Xmax() );
    }
  }

  // Configure canvas
  Int_t nx(1);
  Int_t ny(1);
  Int_t nofResult = histos->GetEntries(); // # of histo
  if ( nofResult == 2 ){
    nx = 2;
    ny=0;
  }
  else if ( nofResult > 2 ){
    ny = TMath::Nint(TMath::Sqrt(nofResult));
    nx = TMath::Nint((nofResult/ny) +0.6);
  }

  TCanvas *c = new TCanvas();

  // --- Configure pad ---
  SetCanvasStyle(c);
  c->Divide(nx,ny,0,0);
  c->SetTitle(Form("%s",fSpectraName.Data()));
  c->Draw();

  AliDebug(1, Form(" Canvas divided in %dx%d",nx,ny));

  TIter nextHisto(histos);
  TH1 * h2;
  Int_t n=0;

  // Loop on Pad
  while ((h2 = static_cast<TH1*>(nextHisto())))
  {
    AliDebug(1,Form(" - subcanvas = %d",n));
    if (h2){
      n++;
      c->cd(n);// got to pad

      // --- Configure histo ---
      Double_t scale = (h2->GetNbinsX())/(h2->GetXaxis()->GetXmax()-h2->GetXaxis()->GetXmin());
      h2->GetXaxis()->SetRangeUser(1.,6.);
      h2->GetXaxis()->SetTitleOffset(1.1);
      h2->SetMinimum(1);
      h2->GetXaxis()->SetLabelFont(42);
      h2->GetXaxis()->SetTitleFont(42);
      h2->GetXaxis()->SetTitleSize(0.08);
      h2->GetXaxis()->SetTitle("m_{#mu#mu} (GeV/#it{c}^{2})");
      h2->GetYaxis()->CenterTitle();
      h2->GetYaxis()->SetLabelFont(42);
      h2->GetYaxis()->SetTitleFont(42);
      h2->SetMarkerSize(0.8);
      h2->SetMarkerStyle(20);
      h2->GetYaxis()->SetTitle(Form("Counts per %.0f MeV/#it{c}^{2}",1000*h2->GetBinWidth(4)));

      h2->Draw();

      // --- Get fitting functions and draw them ---
      if(h2->GetFunction("signal+bck") )      f1 = h2->GetFunction("signal+bck");
      else if(h2->GetFunction("fitMeanpt") )  f1 = h2->GetFunction("fitMeanpt");
      else if(h2->GetFunction("signal") )     f1 = h2->GetFunction("signal");
      if(f1)f1->SetNpx(150);
      if(f1)f1->SetLineColor(kBlue);
      Double_t chi2 =f1->GetChisquare()/f1->GetNDF();

      printf("FitChi2PerNDF = %f",chi2);


      if(h2->GetFunction("signalJPsi") ) f2 = h2->GetFunction("signalJPsi");
      if(f2)f2->SetLineColor(kRed);
      if(f2)f2->SetLineStyle(7);
      if(f2)f2->SetNpx(150);

      if(h2->GetFunction("signalPsiP") )f3 = h2->GetFunction("signalPsiP");
      if(f3)f3->SetLineColor(kRed);
      if(f3)f3->SetLineStyle(7);
      if(f3)f3->SetNpx(150);

      if(h2->GetFunction("bck") )f4 = h2->GetFunction("bck");
      if(f4)f4->SetLineColor(kGray+2);
      if(f4)f4->SetLineStyle(7);
      if(f4)f4->SetNpx(150);

      if(f1) f1->DrawCopy("same");
      if(f2) f2->DrawCopy("same");
      if(f3) f3->DrawCopy("same");
      if(f4) f4->DrawCopy("same");

      // --- Config. first legend pad ---

      TLegend* leg = new TLegend(0.5209804,0.2662884,0.7326179,0.7057458);
      leg->SetTextSize(0.05);
      leg->SetBorderSize(0);

      leg->AddEntry((TObject*)0,Form("#chi^{2}/ndf = %.2f ",chi2perndf[n-1]),"");
      leg->AddEntry((TObject*)0,Form("S/B = %.1f",SoverB[n-1]),"");
      leg->AddEntry((TObject*)0,Form("%s = %.0f +/-  %.0f",what,NofWhat[n-1],NofWhatErr[n-1]),"");
      leg->AddEntry((TObject*)0,Form("m_{J/#psi} = %.0f +/-  %.1f  MeV/#it{c}",1000*massJpsi[n-1],1000*massJpsiErr[n-1]),"");
      leg->AddEntry((TObject*)0,Form("#sigma_{J/#psi} = %.0f +/-  %.1f MeV/#it{c}",1000*sigmaJpsi[n-1],1000*sigmaJpsiErr[n-1]),"");
      leg->Draw("same");

      // --- Playground for a second pad if needed   ---

      TPaveText* pt = new TPaveText(0.55,0.75,0.75,0.95,"brNDC");
      pt->SetBorderSize(0);
      pt->SetFillStyle(0);
      pt->SetLineWidth(2);
      pt->SetTextAlign(12);
      pt->SetTextFont(42);
      pt->SetTextSize(0.058);
      // pt->AddText("ALICE Performance 20/08/2016 ");
      pt->AddText("ALICE pp #sqrt{#it{s}} = 5.02 TeV, L_{int} = 109 #pm 2.1 % nb^{-1}");
      pt->AddText("0 < #it{p}_{T} < 12 GeV/#it{c}^{2}");
      pt->AddText("2.5 < #it{y} < 4");

      // if(BinType.CompareTo("PSI-YVSPT-2DBIN1")==0){
      //   if (n == 1)pt->AddText(Form(" %.0f < #it{p}_{T } < %.0f GeV/#it{c}",BinRange[n-1],BinRange[n]));
      //   else pt->AddText(Form(" %.0f < #it{p}_{T } < %.0f GeV/#it{c}^{2}",BinRange[n],BinRange[n+1]));
      //   pt->AddText(Form("3.25 < #it{y} < 4"));
      // }
      // else if(BinType.CompareTo("PSI-YVSPT-2DBIN2")==0){
      //   if (n == 1)pt->AddText(Form(" %.0f < #it{p}_{T } < %.0f GeV/#it{c}",BinRange[n-1],BinRange[n]));
      //   else pt->AddText(Form(" %.0f < #it{p}_{T } < %.0f GeV/#it{c}^{2}",BinRange[n],BinRange[n+1]));
      //   pt->AddText(Form("2.5 < #it{y} < 3.25"));
      // }
      // else if(BinType.CompareTo("PSI-PT")==0){
      //   if (n == 1)pt->AddText(Form(" %.0f < #it{p}_{T } < %.0f GeV/#it{c}",BinRange[n-1],BinRange[n]));
      //   else pt->AddText(Form(" %.0f < #it{p}_{T } < %.0f GeV/#it{c}^{2}",BinRange[n],BinRange[n+1]));
      //   pt->AddText(Form("2.5 < #it{y} < 4"));
      // }
      // else if(BinType.CompareTo("PSI-Y")==0){
      //   if (n == 1)pt->AddText(Form(" 0 < #it{p}_{T } < 12 GeV/#it{c}^{2}, %.0f < #it{y} < %.0f",BinRange[n-1],BinRange[n]));
      //   else pt->AddText(Form(" 0 < #it{p}_{T } < 12 GeV/#it{c}^{2}, %.0f < #it{y} < %.0f",BinRange[n],BinRange[n+1]));
      // }
      // else if(BinType.CompareTo("PSI-INTEGRATED")==0){
      //   // if (n == 1)pt->AddText(Form(" 0 < #it{p}_{T } < 12 GeV/#it{c}^{2}, %.0f < #it{y} < %.0f",BinRange[n-1],BinRange[n]));
      //   // else pt->AddText(Form(" 0 < #it{p}_{T } < 12 GeV/#it{c}^{2}, %.0f < #it{y} < %.0f",BinRange[n],BinRange[n+1]));
      //   pt->AddText(Form("2.5 < #it{y} < 4"));
      //   pt->AddText(Form(" 0 < #it{p}_{T } < 12 GeV/#it{c}^{2}"));
      // }
      pt->Draw();
      // gPad->Modified();
      // gPad->Update();
    }
    else
    {
      AliError(Form("Cannot find histogram stored at %d ",n));
      continue;
    }
  }
  delete bins;
  delete histos;
}


//_____________________________________________________________________________
TList * AliAnalysisMuMuSpectraCapsulePbPb::RAAasGraphic(Double_t MUL) const
{
   ///
   /// Run over each bin, calculate RAA according to fBinType throught GetValuesFromExternFiles() :
   /// Return a graph to be deleted by owner.
   ///

  // Some constants
  const TString histoName = Form("%s",fSpectraName.Data());

  //Check point
  if(!GetSpectra() || fExternFile.IsNull() ) return 0x0 ;

  //Check point
  if (MUL==0)
  {
    AliError(Form("NofMUL is null"));
    return 0x0;
  }

  // Pointers to handle results and subresults and binning
  AliAnalysisMuMuResult    * result;
  AliAnalysisMuMuBinning   ::Range* r;

 // Array to store bins for the while loop
  TObjArray * bins=GetSpectra()->Binning()->CreateBinObjArray();// (intrinseque 'new')
  if (!bins)
  {
    AliError(Form("Cannot find bins"));
    return 0x0;
  }
  // Array for listed quantities
  Double_t num[8]={0.};
  //  num[0]   ,  num[1]   ,   num[2]   ,   num[3]  ,  num[4] ,  num[5]  ,   num[6]   ,   num[7]
  //  NofJpsi     JPsiStat     JPsiSyst     NormTot    RAA       StatErr     SystCorrErr  SystUnCorrErr
  // --------------------------

  // Here we define some pointers
  TGraphErrors*graph(0x0);
  TGraphErrors*graph_sysUncorr(0x0);

  // One graph with one point for all the bins
  TGraphErrors*graph_syscorr(0x0);
  graph_syscorr = new TGraphErrors(1);


  Double_t    * binArray(0x0) ;// (intrinseque 'new')
  Int_t binsX = 0;

  //________Define histo according to bin type
  if (fSpectraName.Contains("-INTEGRATED")){
    graph           = new TGraphErrors(1);
    graph_sysUncorr = new TGraphErrors(1);
    graph->SetTitle(histoName.Data());
    graph_sysUncorr->SetFillColorAlpha(5,0.05);
    graph_syscorr->SetFillColorAlpha(6,0.05);
  }
  else if (fSpectraName.Contains("-PT")|| fSpectraName.Contains("-Y")) {
    binArray =GetSpectra()->Binning()->CreateBinArray();
    binsX    = GetSpectra()->Binning()->GetNBinsX();

    if (!binArray){// Protection
      AliError(Form("Cannot set binArray"));
      return 0x0;
    }
    if (binsX==0){// Protection
      AliError(Form("Cannot set binsX"));
      return 0x0;
    }
    graph           = new TGraphErrors(binsX);
    graph_sysUncorr = new TGraphErrors(binsX);
    graph->SetTitle(histoName.Data());
    graph->SetMinimum(0.);
    graph->SetMaximum(1.2);
    graph_sysUncorr->SetFillColorAlpha(5,0.05);
    graph_syscorr->SetFillColorAlpha(6,0.05);
  }
  else{// protection
    cout << "Unknowned Bin type !" << endl;
    return 0x0;
  }
  //________

  //________Counters and Iterator for bin
  Int_t nofResult = 0;
  TIter nextBin(bins);
  nextBin.Reset();
  //________

  // Loop on bins
  //==============================================================================
  while ((r = static_cast<AliAnalysisMuMuBinning::Range*>(nextBin())))
  {
    //________Make bin a MuMuResult
    result = GetSpectra()->GetResultForBin(*r);
    if (!result){
      AliError(Form("Cannot find result "));
      return 0x0;
    }
    AliDebug(1, Form("result(%s) = %p ",result->GetName(),result));
    //________

    // Get a string with bin name
    TString binAsString = r->AsString();

    // Store quantities
    num[0] = result->GetValue("NofJPsi");
    num[1] = result->GetErrorStat("NofJPsi");
    num[2] = result->GetRMS("NofJPsi");

    //Main methods
    if(!ComputeRAA(binAsString,&num[0],MUL,r->WidthX())) continue;

    // Set the corr. syst. point at  x=0.5,y=1. Computing for each results, maybe a cleaver way to do...
    graph_syscorr->SetPoint(0,1.,1.);
    graph_syscorr->SetPointError(0,1.,num[6]);

    //________Compute R_AA in case of fully integrated spectra
    if(fSpectraName.Contains("-INTEGRATED")){
      //Output messages
      cout << Form("") << endl;
      printf("  |    %s    || RAA : %.3f +/-  %.3f (%.2f %%) +/- %.3f (%.2f %%)  | Npart :  %.1f   %.1f  | \n"
        ,binAsString.Data(),
        num[4],
        num[5],
        100*num[5]/num[4],
        num[7],
        100*num[7]/num[4],
        fConstArray[0],fConstArray[1]);

      // Add results to TGraphs
      graph->SetPoint(nofResult,fConstArray[0],num[4]);
      graph->SetPointError(nofResult,fConstArray[1],num[5]);
      graph_sysUncorr->SetPoint(nofResult,fConstArray[0],num[4]);
      graph_sysUncorr->SetPointError(nofResult,0.2,num[7]);
    }
    else if (fSpectraName.Contains("-PT")){
      cout << Form("") << endl;
      printf("  |    %s    || RAA : %.3f +/-  %.3f (%.2f %%) +/- %.3f (%.2f %%)  | Npart :  %.1f   %.1f  | \n"
        ,binAsString.Data(),
        num[4],
        num[5],
        100*num[5]/num[4],
        num[7],
        100*num[7]/num[4],
        fConstArray[0],fConstArray[1]);
      //Fill graph
      Double_t binCenter = (binArray[nofResult+1]-binArray[nofResult])/2 + binArray[nofResult] ;
      graph->SetPoint(nofResult,binCenter,num[4]);
      graph->SetPointError(nofResult,r->WidthX()/5,num[5]);
      graph_sysUncorr->SetPoint(nofResult,binCenter,num[4]);
      graph_sysUncorr->SetPointError(nofResult,r->WidthX()/5,num[7]);
    }
    else if (fSpectraName.Contains("-Y")){
      cout << Form("") << endl;
      printf("  |    %s    || RAA : %.3f +/-  %.3f (%.2f %%) +/- %.3f (%.2f %%)  | Npart :  %.1f   %.1f  | \n"
        ,binAsString.Data(),
        num[4],
        num[5],
        100*num[5]/num[4],
        num[7],
        100*num[7]/num[4],
        fConstArray[0],fConstArray[1]);
      //Fill graph
      Double_t binCenter = -((binArray[nofResult+1]-binArray[nofResult])/2 + binArray[nofResult]) ;
      graph->SetPoint(nofResult,binCenter,num[4]);
      graph->SetPointError(nofResult,r->WidthX()/5,num[5]);
      graph_sysUncorr->SetPoint(nofResult,binCenter,num[4]);
      graph_sysUncorr->SetPointError(nofResult,r->WidthX()/5,num[7]);
    }
    else return 0x0;
    //________

    nofResult++;
  }

  // Config. graphics
  if(fSpectraName.Contains("INTEGRATED"))graph->GetXaxis()->SetTitle(Form("<NPart>"));
  else if (fSpectraName.Contains("-PT"))graph->GetXaxis()->SetTitle(Form("PT"));
  else if (fSpectraName.Contains("-Y"))graph->GetXaxis()->SetTitle(Form("Y"));
  graph->GetYaxis()->SetTitle("R_{AA}");
  graph->SetMarkerColor(4);
  graph->SetMarkerStyle(21);

  //Add and merge all Graph
  TList* l = new TList();
  l->SetOwner(kTRUE);
  l->Add(graph);
  l->Add(graph_sysUncorr);
  l->Add(graph_syscorr);

  delete bins;
  delete binArray;

  return l ;
}


//_____________________________________________________________________________
Bool_t AliAnalysisMuMuSpectraCapsulePbPb::ComputeRAA(TString sbin, Double_t numArray[], Double_t MUL, Double_t binwidth) const
{
   ///
   /// Checks bin type and read files (or not) accordingly. Then computes and stores several results in numArray.
   ///

  AliWarning("INNER NORMALIZATION FACTOR, YOU MIGHT CHECK THE CODE !!");

  Double_t CentralityNormalization =1.;
  if (fSpectraName.Contains("V0M_00.00_90.00"))CentralityNormalization       = 1.;
  else if (fSpectraName.Contains("V0M_00.00_20.00")) CentralityNormalization = (2./9.);
  else CentralityNormalization                                               = (1./9.);

  // Taking 2015 value
  if(Mul2015!=0) MUL= Mul2015;

  //________PT and Y case
  if (fSpectraName.Contains("-PT") || fSpectraName.Contains("-Y"))
  {
    // read exterfile and get the correct value
    float valueArray[13];
    //  valueArray[0], valueArray[1], valueArray[2], valueArray[3],     valueArray[4], valueArray[5], valueArray[6], valueArray[7], valueArray[8], valueArray[9], valueArray[10], valueArray[11] valueArray[12]
    //  sigmapp         dsigmapp      dsigmappUnCorr   dsigmappCorr(%)  AccEff         dAccEff        sysMC          TrajEffError    TriggerError  MatchingError,  NofJpsi,       NofJpsiStat,    NofJpsiSys

    if(ReadFromFile(sbin,&valueArray[0])==kFALSE) return kFALSE;
    AliDebug(1, " Values correctly read from extern file");

    //Select the source of NofJpsi
    if(valueArray[10]!=0.)numArray[0]=valueArray[10];
    if(valueArray[11]!=0.)numArray[1]=valueArray[11];
    if(valueArray[12]!=0.)numArray[2]=valueArray[12];

    if(fPrintFlag){
      printf("\n");
      printf("BR                           = %f +/- %f \n", BR,BRerr);
      printf("Fnorm                        = %f +/- %f (%f %% ) +/- %f (%f %% ) \n", Fnorm,FnormStat,100*FnormStat/Fnorm,FnormSyst,100*FnormSyst/Fnorm);
      printf("sigma pp                     = %f +/- %f +/- %f (global %f %%) ub \n", valueArray[0],valueArray[1],valueArray[2],valueArray[3]);
      printf("AccEff                       = %f +/- %f (%f %% )\n", valueArray[4],valueArray[5],100*valueArray[5]/valueArray[4]);
      printf("TAA                          = %f +/- %f (%f %% ) \n", fConstArray[2],fConstArray[3],100*fConstArray[3]/fConstArray[2]);
      printf("Systematic A.P               = %f %% \n", fConstArray[4]);
      printf("MUL                          = %f\n", MUL);
      printf("NofJpsi                      = %f +/- %f (%f %% ) +/- %f (%f %% ) \n", numArray[0],numArray[1],100*numArray[1]/numArray[0],numArray[2],100*numArray[2]/numArray[0]);
      printf("systematic A.P               = %f %% \n", fConstArray[4]);
      printf("sysMC                        = %f %% \n", valueArray[6]);
      printf("TrajEffError                 = %f %% \n", valueArray[7]);
      printf("TriggerError                 = %f %% \n", valueArray[8]);
      printf("TriggerError   (local board) = %f %% \n", TriggerError);
      printf("MatchingError                = %f %% \n", valueArray[9]);
      printf("CentralityNormalization      = %f\n", CentralityNormalization);
      }

     // Compute Raa
     numArray[3] = CentralityNormalization*fConstArray[2]*BR*MUL*Fnorm*(valueArray[0]/1000.)*(valueArray[4]);
     numArray[4] = numArray[0]/numArray[3]/binwidth;

    // Divide by 1.5 in case of pt bins
    if(fSpectraName.Contains("-PT") )numArray[4] =numArray[4]/1.5;

    //Stat error (%)
    numArray[5] =TMath::Sqrt(
      pow(numArray[1]/numArray[0],2)          //signal
    + pow(valueArray[5]/valueArray[4],2)      // accxeff
    + pow(FnormStat/Fnorm,2));                // Fnorm

    numArray[5] = numArray[4]*numArray[5];

    //Corr error in %
    numArray[6]     =  TMath::Sqrt(
    + fConstArray[4]*fConstArray[4]             // syst. AP
    + valueArray[3]*valueArray[3]               // corr. pp cross-section
    + TriggerErrorCent *TriggerErrorCent         // Trigger Local board (%)
    + pow(100*fConstArray[3]/fConstArray[2],2)  // Taa
    + pow(100*FnormSyst/Fnorm,2)                // FNorm
    - pow(100*BRerr,2));                        // Br cancel out

    //UnCorr error absolut
      numArray[7] = numArray[4] * TMath::Sqrt(
      pow(numArray[2]/numArray[0],2)          //signal
      + pow(valueArray[6]/100,2)              // MC input
      + pow(valueArray[7]/100,2)              // tracking
      + pow(valueArray[8]/100,2)              // trigger
      + pow(TriggerError/100,2)               // trigger (local board)
      + pow(valueArray[9]/100,2)              // matching
      + pow(valueArray[1]/valueArray[0],2)    // stat. cross-section pp
      + pow(valueArray[2]/valueArray[0],2));  // uncorr. cross-section pp
  }
  //________Integrated case
  else if(fSpectraName.Contains("-INTEGRATED")){

    if(fConstArray[10]!=0.)numArray[0]=fConstArray[10];
    if(fConstArray[11]!=0.)numArray[1]=fConstArray[11];
    if(fConstArray[12]!=0.)numArray[2]=fConstArray[12];

    if(fPrintFlag){
      printf("\n");
      printf("BR                           = %f +/- %f \n", BR,BRerr);
      printf("Fnorm                        = %f +/- %f (%f %% ) +/- %f (%f %% ) \n", Fnorm,FnormStat,100*FnormStat/Fnorm,FnormSyst,100*FnormSyst/Fnorm);
      printf("sigma pp                     = %f +/- %f +/- %f (%f %%) (global %f %%) ub \n", sigmaPP,dsigmaPP,100*dsigmaPP/sigmaPP,dsigmaPPCorr,100*dsigmaPPCorr/sigmaPP);
      printf("AccEff                       = %f +/- %f (%f %% )\n", fConstArray[8],fConstArray[9],100*fConstArray[9]/fConstArray[8]);
      printf("TAA                          = %f +/- %f (%f %% ) \n", fConstArray[2],fConstArray[3],100*fConstArray[3]/fConstArray[2]);
      printf("MUL                          = %f\n", MUL);
      printf("NofJpsi                      = %f +/- %f (%f %% ) +/- %f (%f %% ) \n", numArray[0],numArray[1],100*numArray[1]/numArray[0],numArray[2],100*numArray[2]/numArray[0]);
      printf("systematic A.P               = %f %% \n", fConstArray[4]);
      printf("TrajEffError                 = %f %% \n", fConstArray[5]);
      printf("TriggerError                 = %f %% \n", fConstArray[6]);
      printf("TriggerError   (local board) = %f %% \n", TriggerError);
      printf("MatchingError                = %f %% \n", fConstArray[7]);
      printf("Systematic MC                = %f %% \n", MCParamError);
      printf("CentralityNormalization      = %f\n", CentralityNormalization);
    }

    //Compute Raa
    numArray[3] = CentralityNormalization*BR*fConstArray[2]*Fnorm*MUL*(sigmaPP/1000)*(fConstArray[8]);
    numArray[4] = numArray[0]/numArray[3];
    //Stat error
    numArray[5] = numArray[4] * AliAnalysisMuMuResult::ErrorAB(numArray[0],numArray[1],fConstArray[8],fConstArray[9]);
    //                               signal
    //Corr error
    numArray[6]        =  TMath::Sqrt(
    MCParamError       *MCParamError             // MC input (%)
    + fConstArray[5]   *fConstArray[5]           // Tracking (%)
    + fConstArray[6]   *fConstArray[6]           // Trigger  (%)
    + TriggerError     *TriggerError             // Trigger Local board (%)
    + TrackingErrorCent*TrackingErrorCent        // Trigger Local board (%)
    + TriggerErrorCent *TriggerErrorCent         // Trigger Local board (%)
    + fConstArray[7]   *fConstArray[7]           // matching. (%)
    + pow(100          *dsigmaPP/sigmaPP,2)      // stat. pp cross section
    + pow(100          *dsigmaPPCorr/sigmaPP,2)  // corr. pp cross-section
    + pow(100          *FnormSyst/Fnorm,2)       // FNorm
    - pow(100          *BRerr,2));               // BR cancel out

    //Uncorrelated error
    numArray[7] = numArray[4] * AliAnalysisMuMuResult::ErrorABC(numArray[0],numArray[2],fConstArray[2],fConstArray[3], 100.,fConstArray[4]);
    //                                                                  signal                       TAA                       AP
  } else {
    AliError("Unowned bin type... I Told you !");
    return kFALSE;
  }
 return kTRUE;
}

//_____________________________________________________________________________
void AliAnalysisMuMuSpectraCapsulePbPb::Print(Option_t* opt) const
{
  /**
   *
   * Print spectra
   *
   */

  //Check point
  if(!GetSpectra()) return ;
  GetSpectra()->Print(opt);
}

//_____________________________________________________________________________
void AliAnalysisMuMuSpectraCapsulePbPb::PrintConst() const
{
    ///
    /// Print member constants on the terminal
    ///

  //Check point
  if(!GetSpectra()) return ;
  else{
    cout <<      " ================================================================ " << endl;
    cout << Form(" Constants for Spectra %s",fSpectraName.Data()) << endl;
    cout <<      " ================================================================ " << endl;
    cout << Form(" -- Value of <Npart>                    = %f",fConstArray[0]) << endl;
    cout << Form(" -- Value of d<Npart>                   = %f",fConstArray[1]) << endl;
    cout << Form(" -- Value of TAA                        = %f",fConstArray[2]) << endl;
    cout << Form(" -- Value of dTAA                       = %f",fConstArray[3]) << endl;
    cout << Form(" -- Value of sys.AP                     = %f",fConstArray[4]) << endl;
    cout << Form(" -- Value of Traj. err.                 = %f",fConstArray[5]) << endl;
    cout << Form(" -- Value of Trigg. err.                = %f",fConstArray[6]) << endl;
    cout << Form(" -- Value of Pair. err.                 = %f",fConstArray[7]) << endl;
    cout << Form(" -- Value of AccEff                     = %f",fConstArray[8]) << endl;
    cout << Form(" -- Value of dAccEff                    = %f",fConstArray[9]) << endl;
    cout << Form(" -- Value of NofJpsi from exterfile     = %f",fConstArray[10]) << endl;
    cout << Form(" -- Value of StatNofJpsi from exterfile = %f",fConstArray[11]) << endl;
    cout << Form(" -- Value of SystNofJpsi from exterfile = %f",fConstArray[12]) << endl;
  }
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuSpectraCapsulePbPb::ReadFromFile(TString sbin, float valueArray[]) const
{
    ///
    /// Read extern file lines and store associated values. Exemple of line :
    ///#intervalLow intervalHight sigmapp     dsigmapp   dsigmappUncorr  dsigmappcorr(%)  AccEff    dAccEff   sysMC   TrajEffError(%)  TriggerError(%)  PairError(%)  NofPsi StatJpsi SystJpsi
    ///00           01            0.8421      0.0255      0.0463        0.052             0.13434   0.00058   1.600    4.0              4.6             1.0           36836  1292     760
    ///
    /// All white space must be single whitespace, i.e " " and not "<tab>"

    Bool_t ok =kFALSE;

    //________Open file
    ifstream infile(fExternFile.Data(),std::ios::in);
    TString line;
    TObjArray* lineArray;

    if (infile)
    {
      AliDebug(1, " ==== opening file ==== ");
      // Loop until end of file is reached
      while(infile.eof()!=kTRUE){

        //read the line
        line.ReadLine(infile,kFALSE);
        if (line.BeginsWith("#"))continue;
        AliDebug(1,Form(" Read line : %s",line.Data()));

        // Put the line in a TObjArray
        lineArray = line.Tokenize(" ");

        // Select the good interval. Since interval is written in <binAsString>, just need them to match
        TString intervalLow  = TString::Format("%.2f",static_cast<TObjString*>(lineArray->At(0))->String().Atof());
        TString intervalHigh = TString::Format("%.2f",static_cast<TObjString*>(lineArray->At(1))->String().Atof());
        if(sbin.Contains(Form("%s",intervalLow.Data())) && sbin.Contains(Form("%s",intervalHigh.Data()))){
            AliDebug(1,Form(" -- line selected -- "));
            ok = kTRUE;
            break;
        }
        else continue;
      }
      infile.close();
      AliDebug(1, " ==== closing file ==== ");

      // Store the value
        for (int i =0 ; i<13 ; i++) {
            valueArray[i]= static_cast<TObjString*>(lineArray->At(i+2))->String().Atof();
        }
        return ok;
    }
    else return ok;
}

//_____________________________________________________________________________
void  AliAnalysisMuMuSpectraCapsulePbPb::SetCanvasStyle(TCanvas *can) const {
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  int font = 42;

  gROOT->SetStyle("Plain");
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(10);
  gStyle->SetCanvasColor(10);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleBorderSize(1);
  gStyle->SetStatColor(10);
  gStyle->SetStatBorderSize(1);
  gStyle->SetLegendBorderSize(1);
  gStyle->SetDrawBorder(0);
  gStyle->SetTextFont(font);
  gStyle->SetStatFont(font);
  gStyle->SetStatFontSize(0.05);
  gStyle->SetStatX(0.97);
  gStyle->SetStatY(0.98);
  gStyle->SetStatH(0.03);
  gStyle->SetStatW(0.3);
  gStyle->SetTickLength(0.02,"y");
  gStyle->SetEndErrorSize(3);
  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetLabelFont(font,"xyz");
  gStyle->SetLabelOffset(0.01,"xyz");
  gStyle->SetTitleFont(font,"xyz");
  gStyle->SetTitleOffset(1.1,"xy");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetMarkerSize(1.3);
  gStyle->SetPalette(1,0);

  gROOT->ForceStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  gStyle->SetLegendFont(42);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(10);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetEndErrorSize(0);


  // can->SetFillColor(0);

  can->SetBorderMode(0);

  can->SetBorderSize(0);

  can->SetLeftMargin(0.18);
  can->SetRightMargin(0.1);
  can->SetBottomMargin(0.1518219);
  can->SetTopMargin(0.);
  can->SetFrameBorderMode(0);

}
