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

#include "AliAnalysisMuMuSpectraCapsulePbP.h"

ClassImp(AliAnalysisMuMuSpectraCapsulePbP)

#include "TF1.h"
#include "TProfile.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "THashList.h"
#include "TLegend.h"
#include "TLine.h"
#include "TList.h"
#include "TMath.h"
#include "TObjArray.h"
#include "AliAnalysisMuMuBinning.h"
#include "AliLog.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "AliCounterCollection.h"
#include "AliMergeableCollection.h"
#include "AliAnalysisMuMuSpectra.h"
#include "AliAnalysisMuMuResult.h"
#include "AliAnalysisMuMuJpsiResult.h"
#include "AliAnalysisMuMuSpectraCapsulePbP.h"
#include <iostream>
#include <string>

using std::cout;
using std::endl;
using std::ifstream;

namespace
{

  const Double_t BR             = 5.93/100; // Branching ratio
  //Normalization factor
  //FIXME : Fnorm store in TH1, make in general
  const Double_t Fnorm          = 27.51;    // Normalization
  const Double_t FnormStat      = 0.01;     // Normalization
  const Double_t FnormSyst      = 0.97;     // Normalization
  //pp Cross-section integrated in pt,y
  const Double_t sigmaPP        = 3.343;    // for fully integrated case
  const Double_t dsigmaPP       = 0.033;    // idem
  const Double_t dsigmaPPCorr   = 0.022;    // for fully integrated case
  const Double_t dsigmaPPUncorr = 0.021;    // idem
  // Global MC sys. err. for centrality integrated in pt and Y
  const Double_t MCParamError   = 3/100;
  // Corr. error for centrality
  const Double_t TrajCENT       = 11/100;
  const Double_t TriggCENT      = 2/100;
  const Double_t PairCENT       = 1/100;
  // Corr. error for pt case
  const Double_t TrajPT         = 1/100;
  const Double_t TriggPT        = 1/100;
  const Double_t PairPt         = 1/100;
  // Corr. error for y case
  const Double_t TrajY          = 1/100;
  const Double_t TriggY         = 1/100;
  const Double_t PairY          = 1/100;
}



//_____________________________________________________________________________
 AliAnalysisMuMuSpectraCapsulePbP::AliAnalysisMuMuSpectraCapsulePbP(
const AliAnalysisMuMuSpectra*  spectra,
const TString                 spectraPath,
const char                  * externFile,
const char                  * externFile2)
:
  AliAnalysisMuMuSpectraCapsule(spectra,spectraPath,externFile,externFile2)
  // fSpectra(spectra),
  // fSpectraName(spectraPath),
  // fExternFile(externFile),
  // fExternFile2(externFile2)
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
AliAnalysisMuMuSpectraCapsulePbP::~AliAnalysisMuMuSpectraCapsulePbP()
{
  // dtor
}

//_____________________________________________________________________________
TGraphErrors* AliAnalysisMuMuSpectraCapsulePbP::ComputeYield( const char* what, const TH1* histo, const char* sResName)
{
  /// @brief Compute Yield.
  /// @argument what  the yield nominator, i.e NofJPsi, meanPT etc. (null by default)
  /// @argument histo  histogramme of Equivalent MinBias

  if(!GetSpectra() || histo==0x0|| strcmp(what,"")==1) return 0x0;

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
  TGraphErrors*graph(0x0);
  // TGraphErrors*graph_sysUncorr(0x0);

  Double_t    * binArray(0x0) ;// (intrinseque 'new')
  Int_t binsX = 0;

  //________Define histo according to bin type
  if (GetSpectraName().Contains("-INTEGRATED"))
  {
    graph           = new TGraphErrors(1);
    // graph_sysUncorr = new TGraphErrors(1);
    graph->SetTitle(graphTitle.Data());

    // graph_sysUncorr->SetFillColorAlpha(5,0.05);
  }
  else if (GetSpectraName().Contains("-PT")|| GetSpectraName().Contains("-Y"))
  {
    binArray =GetSpectra()->Binning()->CreateBinArray();
    binsX = GetSpectra()->Binning()->GetNBinsX();

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

    Double_t nEqMBTot      = histo->GetBinContent(nofResult+1);
    Double_t nEqMBTotError = histo->GetBinError(nofResult+1);

    AliDebug(1,Form("histo used    : %s",histo->GetTitle()));
    AliDebug(1,Form("%s            = %f",what,NofWhattTot));
    AliDebug(1,Form("%s error      = %f",what,NofWhattTotError));
    AliDebug(1,Form("nEqMBTot      = %f",nEqMBTot));
    AliDebug(1,Form("nEqMBTotError = %f",nEqMBTotError));

    if( NofWhattTot==0||NofWhattTotError==0||nEqMBTot==0||nEqMBTotError==0)
    {
      AliError("Cannot set quantities properly");
      return 0x0;
    }

    //________Compute R_AA in case of fully integrated spectra
    if(GetSpectraName().Contains("-INTEGRATED"))
    {
      Double_t yieldInt = NofWhattTot/(nEqMBTot*BR);
      Double_t yieldIntError = yieldInt*AliAnalysisMuMuResult::ErrorAB(NofWhattTot,NofWhattTotError,nEqMBTot,TMath::Sqrt(nEqMBTot));

      // Add results to TGraphs
      graph->SetPoint(nofResult,fConstArray[0],yieldInt);
      graph->SetPointError(nofResult,fConstArray[1],yieldIntError);
      // graph_sysUncorr->SetPoint(nofResult,0,num[4]);
      // graph_sysUncorr->SetPointError(nofResult,0.2,num[7]);
    }
    else
    {
      Double_t yieldInt = NofWhattTot/(nEqMBTot*BR);
      Double_t yieldIntError = yieldInt*AliAnalysisMuMuResult::ErrorAB(NofWhattTot,NofWhattTotError,nEqMBTot,TMath::Sqrt(nEqMBTot));

      //Fill graph
      Double_t binCenter = (binArray[nofResult+1]-binArray[nofResult])/2 + binArray[nofResult];
      graph->SetPoint(nofResult,binCenter,yieldInt);
      graph->SetPointError(nofResult,r->WidthX()/5,yieldInt);
      // graph_sysUncorr->SetPoint(nofResult,binCenter,num[4]);
      // graph_sysUncorr->SetPointError(nofResult,r->WidthX()/5,num[7]);
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

  // delete graph;
  // delete graph_sysUncorr;
  delete bins;
  delete binArray;

 return graph ;

}

//_____________________________________________________________________________
void AliAnalysisMuMuSpectraCapsulePbP::DrawResults( const char* particle,const char* subresults) const
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

  //Pointer for functions
    TF1* f1 = 0x0;
    TF1* f2 = 0x0;
    TF1* f3 = 0x0;
    TF1* f4 = 0x0;

    const TString sres(subresults);


  // Arrays
  TObjArray* histos = new TObjArray(0x0);// Array to store histograms
  TObjArray* bins=GetSpectra()->Binning()->CreateBinObjArray();// Array to store bins

  if (!bins)
  {
    AliError(Form("Cannot find bins"));
    return;
  }

  // Settings for histo
  Double_t xmin(-1);
  Double_t xmax(-1);
  if ( fSpectraName.Contains(particle))
      {
      xmin = 2;
      xmax = 6;
      }

  //Iterator for bin
  TIter nextBin(bins);

  // Loop on bins
  //==============================================================================
  while ((r = static_cast<AliAnalysisMuMuBinning::Range*>(nextBin())))
  {
    // Make bin a MuMuResult
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
      if (!subresult)
      {
      AliError(Form("Cannot find subresult "));
      return;
      }
      AliDebug(1,Form("subresult(%s) = %p",sr->GetName(),subresult));

      if(!sres.IsNull() && !sres.Contains(sr->GetName())) continue;

      // Get histo
      if ( subresult ) h = (TH1*)subresult->Histo();
      AliDebug(1, Form(" - Histo(%s) = %p ",h->GetTitle(),h));

      // Store it
      if(h) {histos->Add(h);}
      else
      {
        AliError(Form("Cannot set histo result "));
        return;
      }
    }
  }
  //Configure canvas
  Int_t nx(1);
  Int_t ny(1);
  Int_t nofResult = histos->GetEntries(); // # of histo
  if ( nofResult == 2 )
      {
      nx = 2;
      ny=0;
      }
  else if ( nofResult > 2 )
      {
      ny = TMath::Nint(TMath::Sqrt(nofResult));
      nx = TMath::Nint((nofResult/ny) +0.6);
      }
  TCanvas *c = new TCanvas;
  c->Draw();
  c->Divide(nx,ny);
  c->SetTitle(Form("%s",fSpectraName.Data()));
  gStyle->SetOptFit(1112);
  AliDebug(1, Form(" Canvas divided in %dx%d",nx,ny));

  //Iterator on histos + counter
  TIter nextHisto(histos);
  TH1 * h2;
  Int_t n=0;
  // Loop on Pad
  while ((h2 = static_cast<TH1*>(nextHisto())))
  {
    AliDebug(1,Form(" - subcanvas = %d",n));
    h = static_cast<TH1*>(histos->At(n));
    if (h)
    {
      ++n;
      c->cd(n);// got to pad
      if (xmin>0)
        {
          // Loop to configure the pad as you like
          h->GetXaxis()->SetRangeUser(xmin,xmax);
          h->SetTitleSize(10);
        }
      h->DrawCopy("histes");

      //Get fitting functions and draw them
      f1 = h->GetFunction("signal+bck");
      f2 = h->GetFunction("signalJPsi");
      f3 = h->GetFunction("signalPsiP");
      f4 = h->GetFunction("bck");
      if(f1) f1->DrawCopy("same");
      if(f2) f2->DrawCopy("same");
      if(f3) f3->DrawCopy("same");
      if(f4)
          {
            f4->SetLineColor(kBlue);
            f4->SetLineStyle(16);
            f4->DrawCopy("same");
          }
      gPad->Modified();
      gPad->Update();
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
TGraphErrors * AliAnalysisMuMuSpectraCapsulePbP::RpAAsGraphic(Double_t MUL) const
{
  /**
   *
   * Run over each bin, calculate R_pA according to fBinType throught GetValuesFromExternFiles() :
   * Return a graph to be deleted by owner.
   *
   */

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

  AliError("To be implemented !");
  return 0x0;


 //  // Pointers to handle results and subresults and binning
 //  AliAnalysisMuMuResult    * result;
 //  AliAnalysisMuMuBinning   ::Range* r;

 // // Array to store bins for the while loop
 //  TObjArray * bins=GetSpectra()->Binning()->CreateBinObjArray();// (intrinseque 'new')
 //  if (!bins)
 //  {
 //    AliError(Form("Cannot find bins"));
 //    return 0x0;
 //  }
 //  // Array for listed quantities
 //  Double_t num[8]={0.};
 //  //  num[0]   ,  num[1]   ,   num[2]   ,   num[3]  ,  num[4] ,  num[5]  ,   num[6]   ,   num[7]
 //  //  NofJpsi     JPsiStat     JPsiSyst     NormTot    RAA       StatErr     SystCorrErr  SystUnCorrErr
 //  // --------------------------

 //  // Here we define some pointers
 //  TGraphErrors*graph(0x0);
 //  // TGraphErrors*graph_sysUncorr(0x0);

 //  Double_t    * binArray(0x0) ;// (intrinseque 'new')
 //  Int_t binsX = 0;

 //  //________Define histo according to bin type
 //  if (fSpectraName.Contains("-INTEGRATED"))
 //  {
 //    graph           = new TGraphErrors(1);
 //    // graph_sysUncorr = new TGraphErrors(1);
 //    graph->SetTitle(histoName.Data());
 //    graph->SetMinimum(0.);
 //    graph->SetMaximum(1.2);
 //    // graph_sysUncorr->SetFillColorAlpha(5,0.05);
 //  }
 //  else if (fSpectraName.Contains("-PT")|| fSpectraName.Contains("-Y"))
 //  {
 //    binArray =GetSpectra()->Binning()->CreateBinArray();
 //    binsX    = GetSpectra()->Binning()->GetNBinsX();

 //    if (!binArray)
 //    {
 //      AliError(Form("Cannot set binArray"));
 //      return 0x0;
 //    }
 //    if (binsX==0)
 //    {
 //      AliError(Form("Cannot set binsX"));
 //      return 0x0;
 //    }

 //    graph           = new TGraphErrors(binsX);
 //    // graph_sysUncorr = new TGraphErrors(binsX);
 //    graph->SetTitle(histoName.Data());
 //    graph->SetMinimum(0.);
 //    graph->SetMaximum(1.2);
 //    // graph_sysUncorr->SetFillColorAlpha(5,0.05);
 //  }
 //  else
 //  {
 //    cout << "Unknowned Bin type !" << endl;
 //    return 0x0;
 //  }
 //  //________

 //  //________Counters and Iterator for bin
 //  Int_t nofResult = 0;
 //  TIter nextBin(bins);
 //  nextBin.Reset();
 //  //________

 //  // Loop on bins
 //  //==============================================================================
 //  while ((r = static_cast<AliAnalysisMuMuBinning::Range*>(nextBin())))
 //  {
 //    //________Make bin a MuMuResult
 //    result = GetSpectra()->GetResultForBin(*r);
 //    if (!result)
 //    {
 //      AliError(Form("Cannot find result "));
 //      return 0x0;
 //    }
 //    AliDebug(1, Form("result(%s) = %p ",result->GetName(),result));
 //    //________

 //    // Get a string with bin name
 //    TString binAsString = r->AsString();

 //    // Store quantities
 //    num[0] = result->GetValue("NofJPsi");
 //    num[1] = result->GetErrorStat("NofJPsi");
 //    num[2] = result->GetRMS("NofJPsi");

 //    GetValuesFromExternFile(binAsString,&num[0],MUL);

 //    //________Compute R_AA in case of fully integrated spectra
 //    if(fSpectraName.Contains("-INTEGRATED"))
 //    {
 //      //Output messages
 //      cout << Form("") << endl;
 //      cout << Form("  |    %s    |   %.0f  %.0f  %.0f  |  %.3f  %.3f  %.3f  |  %.0f   %.0f  |"  ,binAsString.Data(),num[0],num[1],num[2],num[4],num[5],num[7],fConstArray[0],fConstArray[1]) << endl;

 //      // Add results to TGraphs
 //      graph->SetPoint(nofResult,fConstArray[0],num[4]);
 //      graph->SetPointError(nofResult,fConstArray[1],num[5]);
 //      // graph_sysUncorr->SetPoint(nofResult,0,num[4]);
 //      // graph_sysUncorr->SetPointError(nofResult,0.2,num[7]);
 //    }
 //    else
 //    {
 //      num[4]=num[4]/(r->WidthX());
 //      cout << Form("") << endl;
 //      cout << Form("  | %s |   %.0f  %.0f  %.0f  |  %.3f  %.3f  %.3f  |  %.0f   %.0f  |"  ,binAsString.Data(),num[0],num[1],num[2],num[4],num[5],num[7],fConstArray[0],fConstArray[1]) << endl;
 //      //Fill graph
 //      Double_t binCenter = (binArray[nofResult+1]-binArray[nofResult])/2 + binArray[nofResult] ;
 //      graph->SetPoint(nofResult,binCenter,num[4]);
 //      graph->SetPointError(nofResult,r->WidthX()/5,num[5]);
 //      // graph_sysUncorr->SetPoint(nofResult,binCenter,num[4]);
 //      // graph_sysUncorr->SetPointError(nofResult,r->WidthX()/5,num[7]);
 //    }
 //    //________

 //    nofResult++;
 //  }

 //  // Config. graphics
 //  if(fSpectraName.Contains("INTEGRATED"))graph->GetXaxis()->SetTitle(Form("<NPart>"));
 //  else if (fSpectraName.Contains("-PT"))graph->GetXaxis()->SetTitle(Form("PT"));
 //  else if (fSpectraName.Contains("-Y"))graph->GetXaxis()->SetTitle(Form("Y"));
 //  graph->GetYaxis()->SetTitle("R_{pA}");
 //  graph->SetMarkerColor(4);
 //  graph->SetMarkerStyle(21);

 //  // delete graph;
 //  // delete graph_sysUncorr;
 //  delete bins;
 //  delete binArray;

 // return graph ;
}


//_____________________________________________________________________________
void AliAnalysisMuMuSpectraCapsulePbP::GetValuesFromExternFile(TString sbin, Double_t numArray[], Double_t MUL) const
{
  /**
   *
   * Checks bin type and read files (or not) accordingly. Then computes and stores several results in numArray.
   *
   */

  AliWarning("INNER NORMALIZATION FACTOR, YOU MIGHT CHECK THE CODE !!");

  AliError("To be implemented !");
  return;

  // //________PT and Y case
  // if (fSpectraName.Contains("-PT") || fSpectraName.Contains("-Y"))
  // {
  //   char status;

  //    //________Arrays to store quantities from externFile
  //   float intervalArray[2];
  //   // intervalLow , intervalHight
  //   //      0      ,      1
  //   float valueArray[10];
  //   //    sigmapp   dsigmapp   dsigmappCorr   dsigmappUncorr  AccEff   dAccEff  sysMC   TrajEffError  TriggerError   PairError
  //   //________ﬁ

  //   //________Open file
  //   FILE*  infile;
  //   infile = fopen(fExternFile.Data(),"rb") ;

  //   if (infile != NULL)
  //   {
  //     AliDebug(1, " ==== opening file ==== ");
  //     // Loop until end of file is reached
  //     while(!feof(infile))
  //     {
  //       // Reminder :
  //       // intervalLow  intervalHight  sigma_pp  dsigma_pp  dsigma_pp_Correl  dsigma_pp_Uncorrel  AccEff  dAccEff sysMC TrajEffError  TriggerError   PairError;
  //       // Store value in array
  //       fscanf(infile,"%s %f_%f %f %f %f %f %f %f %f %f %f %f",&status,&intervalArray[0],&intervalArray[1],
  //                                                  &valueArray[0],&valueArray[1],&valueArray[2],&valueArray[3],&valueArray[4],
  //                                                  &valueArray[5],&valueArray[6],&valueArray[7],&valueArray[8],&valueArray[9]);
  //       if(status == 'F') continue; // F = false, T =true
  //       // Make intervalArray a string
  //       TString intervalLow  = TString::Format("%.2f",intervalArray[0]);
  //       TString intervalHigh = TString::Format("%.2f",intervalArray[1]);

  //       // Select the good interval. Since interval is written in <binAsString>, just need them to match
  //       if(sbin.Contains(Form("%s",intervalLow.Data())) && sbin.Contains(Form("%s",intervalHigh.Data())))
  //       {
  //         // Check Point
  //         AliDebug(1,Form(" -- Selected line :"));
  //         AliDebug(1,Form(" -- intervalLow  intervalHight  sigma_pp  dsigma_pp  dsigma_pp_Correl  dsigma_pp_Uncorrel  AccEff  dAccEff  sysMC TrajEffError  TriggerError   PairErrorﬁ"));
  //         AliDebug(1,Form(" --  %.2f  %.2f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f ",
  //         intervalArray[0],intervalArray[1],valueArray[0],valueArray[1],valueArray[2],valueArray[3],valueArray[4],valueArray[5],valueArray[6],valueArray[7],valueArray[8],valueArray[9]));

  //         //Normalization according to centrality bin
  //         if (fSpectraName.Contains("V0M_00.00_90.00"))
  //         {
  //           numArray[3] = fConstArray[2]*BR*MUL*Fnorm*(valueArray[0]/1000.)*(valueArray[4]);
  //         }
  //         else
  //         {
  //           numArray[3] = (1./9.)*fConstArray[2]*BR*MUL*Fnorm*(valueArray[0]/1000.)*(valueArray[4]);
  //         }

  //         numArray[4] = numArray[0]/numArray[3];

  //       }
  //       else
  //       {
  //         AliDebug(1,Form("Not the good interval, so continue ...."));
  //         continue;
  //       }
  //     }
  //     fclose(infile);
  //     AliDebug(1, " ==== Closing file ==== ");
  //   }
  //   else
  //   {
  //     cout << Form("Cannot open configuration file %s ",fExternFile.Data()) << endl;
  //     return;
  //   }
  //   //________

  //   if (fSpectraName.Contains("-PT") )
  //   {
  //     // Normalization factor due to how PP cross-section are calculated
  //     numArray[4] = numArray[4]/1.5;

  //     //Corr
  //     numArray[6] = numArray[4] * TMath::Sqrt(TrajPT                       *TrajPT                       + // Traj. reconstruction Eff.
  //                                             TriggPT                      *TriggPT                      + // Trigg Eff.
  //                                             FnormSyst/Fnorm              *FnormSyst/Fnorm              + // Fnorm Syst
  //                                             fConstArray[3]/fConstArray[2]*fConstArray[3]/fConstArray[2]+ // TAA syst.
  //                                             valueArray[2]/valueArray[0]  *valueArray[2]/valueArray[0]);  // dsigma_pp_Corr/sigma_pp


  //   }
  //   else if (fSpectraName.Contains("-Y"))
  //   {
  //     //Corr
  //     numArray[6] = numArray[4] * TMath::Sqrt(TrajY                        *TrajY                         + // Traj. reconstruction Eff.
  //                                             TriggY                       *TriggY                        + // Trigg Eff.
  //                                             FnormSyst/Fnorm              *FnormSyst/Fnorm               + // Fnorm Syst
  //                                             fConstArray[3]/fConstArray[2]*fConstArray[3]/fConstArray[2] + // TAA syst.
  //                                             valueArray[2]/valueArray[0]  *valueArray[2]/valueArray[0]);   // dsigma_pp_Corr/sigma_pp

  //   }
  //   else
  //   {
  //     AliError("Unowned bin type... I Told you !");
  //     return;
  //   }

  //   //Stat
  //   numArray[5] = numArray[4] * TMath::Sqrt(numArray[1]/numArray[0]    *numArray[1]/numArray[0]     + // Jpsi extraction
  //                                           valueArray[1]/valueArray[0]*valueArray[1]/valueArray[0] ); // dsigma_pp/sigma_pp

  //   //UnCorr
  //   numArray[7] = numArray[4] * TMath::Sqrt(valueArray[6]             *valueArray[6]                   + // MC param.
  //                                           numArray[2]/numArray[0]   *numArray[2]/numArray[0]         + // Signal extraction
  //                                           valueArray[7]             *valueArray[7]                   + // Traj. Eff.
  //                                           valueArray[8]             *valueArray[8]                   + // Trigg. Eff.
  //                                           valueArray[9]             *valueArray[9]                   + // Pair. Eff.
  //                                           valueArray[3]/valueArray[0]*valueArray[3]/valueArray[0]);     // dsigma_pp_Uncorr/sigma_pp
  // }

  // //________Compute R_AA in case of integrated spectra in PT and Y
  // else if(fSpectraName.Contains("-INTEGRATED"))
  // {
  //   //Get quantities

  //   //Normalization according to centrality bin
  //   if (!fSpectraName.Contains("V0M_00.00_90.00")) numArray[3] = BR*fConstArray[2]*(Fnorm*MUL/9)*(sigmaPP/1000)*(fConstArray[8]);
  //   else                                           numArray[3] = BR*fConstArray[2]*Fnorm*MUL*(sigmaPP/1000)*(fConstArray[8]);
  //   numArray[4] = numArray[0]/numArray[3];

  //   //Stat error
  //   numArray[5] = numArray[4] * TMath::Sqrt((numArray[1]/numArray[0])*(numArray[1]/numArray[0]) + // Jpsi extraction
  //                                                    dsigmaPP/sigmaPP*dsigmaPP/sigmaPP          ); // dsigma_pp/sigma_pp
  //   //Corr error
  //   numArray[6] = numArray[4] * TMath::Sqrt(MCParamError        *MCParamError         + // MCParamError
  //                                           fConstArray[4]      *fConstArray[4]       + // Traj. reconstruction Eff.
  //                                           fConstArray[5]      *fConstArray[5]       + // Trig. Eff.
  //                                           fConstArray[6]      *fConstArray[6]       + // Pair Reconst. Eff.
  //                                           dsigmaPPCorr/sigmaPP*dsigmaPPCorr/sigmaPP + // dsigma_pp_Corr/sigma_pp
  //                                           FnormSyst/Fnorm     *FnormSyst/Fnorm);      // Fnorm Syst
  //   //Uncorr error
  //   numArray[7] = numArray[4] * TMath::Sqrt(numArray[2]/numArray[0]      *numArray[2]/numArray[0]      + // Signal extraction
  //                                           fConstArray[1]/fConstArray[0]*fConstArray[1]/fConstArray[0]);// TAA syst.
  // }
  // else
  // {
  //   AliError("Unowned bin type... I Told you !");
  //   return;
  // }
}

//_____________________________________________________________________________
void AliAnalysisMuMuSpectraCapsulePbP::Print(Option_t* opt) const
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
void AliAnalysisMuMuSpectraCapsulePbP::PrintConst() const
{
  /**
   *
   * Print member constants
   *
   */

  //Check point
  if(!GetSpectra()) return ;
  else
  {
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
