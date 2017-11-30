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

#include "AliAnalysisMuMuSpectraProcessor.h"

ClassImp(AliAnalysisMuMuSpectraProcessor)


#include "AliLog.h"
#include "TObject.h"
#include <TString.h>
#include <iostream>
#include <string>
#include "AliAnalysisMuMuResult.h"
#include "AliAnalysisMuMuJpsiResult.h"
#include "AliAnalysisMuMuBinning.h"
#include "TObjArray.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TH1.h"
#include "TH2.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TLatex.h"
#include "TList.h"
#include <fstream>
#include <string>

using std::cout;
using std::endl;
using std::ifstream;


//_____________________________________________________________________________
AliAnalysisMuMuSpectraProcessor::AliAnalysisMuMuSpectraProcessor(
const AliAnalysisMuMuSpectra*  spectra,
const TString                 spectraPath,
const char                  * externFile,
const char                  * externFile2) :
TObject(),
fSpectra(spectra),
fSpectraName(spectraPath),
fExternFile(externFile),
fExternFile2(externFile2),
fPrintFlag(kFALSE)
{
  /// Default ctor
  //Check point
  if (!fSpectra)
  {
    AliError(Form("Cannot find spectra wih name %s Please check the name",fSpectra->GetName()));
    return;
  }
  AliDebug(1, Form(" - spectra(%s) = %p ",fSpectra->GetName(),fSpectra));


  if (fSpectraName.IsNull())
  {
    AliWarning(Form("No spectra name ! "));
    return;
  }
}

//_____________________________________________________________________________
AliAnalysisMuMuSpectraProcessor::~AliAnalysisMuMuSpectraProcessor()
{
  // dtor
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuSpectraProcessor::SetConstantFromExternFile(const char* file, Double_t* constantArray, const TString* spectraName)
{
  /// Set member constants depending on centrality bin from an ewternfile.
  /// If values are empty and can be obtained from a graph provided by the AliAnalysisMuMu Framework, this value is set by default
  /// For the PP capsule ine could be :
  ///
  /// #centrality Low  High  lumi.    lumi (stat)  lumi (syst. %)  Trigg  Trigg (local board)   Traj.err.(%)  MC Input (%)  Matching(%)  AccEff  dAccEff  NofJpsi Stat.Jpsi SystJpsi
  /// PP          PP   PP    0.         0.0         0.0            0.0    0.0                    04           03            01           0.0     0.0      0.0     0.0        0.0
  ///
  /// Note that for the PP case, the centrality limits are irrelevant
  ///
  /// For PbPb capsule :
  /// #centrality Low  High  <npart>    d<npart>    TAA           dTAA      sys.AP(%)           Traj.err.(%)   Trigg.err.(%) Matching(%)   AccEff   dAccEff  NofJpsi Stat.Jpsi SystJpsi
  /// V0M         00   10    359        31.2        23.4          0.351     2.00                04             03             01            0.1297   0.00040  105159  1693      488

    // Reset on fConstant
    for (int i = 0; i < 13; ++i) constantArray[i]=0.;
    Bool_t ok= kFALSE;
    AliDebug(1,Form("Reading from file %s",file));

    //________Open file
    ifstream infile(file,std::ios::in);
    TString line;
    TObjArray* lineArray;

    if (infile){
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
            TString centrality   =  static_cast<TObjString*>(lineArray->At(0))->String().Data();
            TString intervalLow  =  TString::Format("%.2f",static_cast<TObjString*>(lineArray->At(1))->String().Atof());
            TString intervalHigh =  TString::Format("%.2f",static_cast<TObjString*>(lineArray->At(2))->String().Atof());
            AliDebug(1,Form(" --__--__-- interval low = %s",intervalLow.Data()));
            AliDebug(1,Form(" --__--__-- interval high = %s",intervalHigh.Data()));
            if (intervalLow.EqualTo("0.00")) intervalLow ="00.00";

            // Select the good interval for PbPb case. Since interval is written in <binAsString>, just need them to match
            if(spectraName->Contains(Form("%s",centrality.Data()))&& spectraName->Contains(Form("%s_%s",intervalLow.Data(),intervalHigh.Data())) && spectraName->Contains(Form("%s_%s",centrality.Data(),intervalLow.Data()))){
                AliDebug(1,Form(" spectraName = %s",spectraName->Data()));
                AliDebug(1,Form(" -- line selected -- "));
                ok = kTRUE;
                break;
            }
            // PP case
            else if(centrality.Contains("PP")){
                AliDebug(1,Form(" spectraName = %s",spectraName->Data()));
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
            constantArray[i]= static_cast<TObjString*>(lineArray->At(i+3))->String().Atof();
        }
        return ok;
    }
    else return ok;
}


//_____________________________________________________________________________

void AliAnalysisMuMuSpectraProcessor::PrintFitParam(const char* subresult, const char* param) const
{
  /**
   *   [AliAnalysisMuMuSpectraProcessor::PrintFitParam description]
   *   @brief   [bief description]
   *   @details [long description]
   *
   *   @param   subresult [description]
   *   @param   param     [description]
   */

  TCanvas * c = new TCanvas;
  TH1*       h          = 0x0;
  TH1*       hcent      = 0x0;
  TH1*       href       = 0x0;

  // Get param as string
  TString par(param);
  TObjArray * paramArray = par.Tokenize(",");
  int nparam = paramArray->GetEntries();

  // Get subresults as string
  TString ssubres(subresult);
  TObjArray * subresultArray = new TObjArray;
  // Get the name of the subresults if chains empty
  if(ssubres.IsNull()){
    printf("toto\n");
    subresultArray = static_cast<AliAnalysisMuMuResult*>(fSpectra->BinContentArray()->At(0))->SubResults();
    // subresultArray->Print();
    // return;
  } else subresultArray = ssubres.Tokenize(",");
  TObjString* ssubresult;
  TIter nextSubResult(subresultArray);

  c->Divide(2,nparam);

  // Loop on param
  int k = 1;
  TObjString* sparam;
  TIter nextparam(paramArray);
  nextparam.Reset();
  while ( ( sparam = static_cast<TObjString*>(nextparam()) ) )
  {
    AliDebug(1,Form("param %s",sparam->String().Data()));

    if(c)
    {
      // --- First canvas ---
      c->cd(k);

      TLegend*leg = new TLegend(0.1,0.7,0.48,0.9);
      leg->SetHeader(Form("Fit Parameters %s ",sparam->String().Data()));
      leg->SetTextSize(0.03);

      printf("going to subcanvas %d\n",k );
      int i = 1;
      nextSubResult.Reset();
      while ( ( ssubresult = static_cast<TObjString*>(nextSubResult()) ) )
      {
        //Loop over subresults
        AliDebug(1,Form("-----SubResults %s",ssubresult->String().Data()));
        if  ( !fSpectraName.Contains("VS") )
          h = fSpectra->Plot(sparam->String().Data(),ssubresult->String().Data(),kFALSE);
        else if ( fSpectraName.Contains("YVSPT") )
          h = static_cast<TH2*>(fSpectra->Plot(sparam->String().Data(),ssubresult->String().Data(),kFALSE))->ProjectionX();
        else if ( fSpectraName.Contains("PTVSY") )
          h = static_cast<TH2*>(fSpectra->Plot(sparam->String().Data(),ssubresult->String().Data(),kFALSE))->ProjectionX();

        if(!h) {
          AliError(Form("Cannot find histo for SubResults %s",ssubresult->String().Data()));
          return;
        }

        // beautifull histo
        if( i!=3 && i!=5 && i!=10 && i!=11 && i!=12 && i!=13 && i!=14 ) h->SetMarkerColor(i); //nobody likes green and yellow
        else               h->SetMarkerColor(i+5);

        h->SetMarkerSize(1.);
        h->SetMarkerStyle(20+i);
        if(i==11)h->SetMarkerStyle(20+i+3);
        if(i==1)
        {
          h->GetYaxis()->SetTitleSize(0.05);
          h->GetYaxis()->SetLabelSize(0.05);
          h->GetXaxis()->SetLabelSize(0.05);
          h->GetXaxis()->SetTitleSize(0.05);
          h->SetTitle(Form(" %s for bin %s",sparam->String().Data(),fSpectraName.Data()));
          h->GetYaxis()->SetTitle(sparam->String().Data());
          if (fSpectraName.Contains("YVSPT") )
            h->GetXaxis()->SetTitle("PT");
          else if (fSpectraName.Contains("PTVSY") )
            h->GetXaxis()->SetTitle("Y");
        }

        if(! sparam->String().Contains("FitChi2PerNDF"))
        {
          if(i==1)h->DrawCopy();
          else    h->DrawCopy("same");
        }
        else
        {
          if(i==1)h->DrawCopy("p");
          else    h->DrawCopy("samep");
        }

        leg->AddEntry(h,Form("%s with %s",sparam->String().Data(),ssubresult->String().Data()),"p");
        i++;
      }

      leg->Draw("same");

      // --- ratio ---
      c->cd(++k);

      TLegend * leg2 = new TLegend(0.1,0.7,0.48,0.9);
      leg2->SetHeader(Form("Ratio"));
      leg2->SetTextSize(0.03);

      nextSubResult.Reset();
      int j= 1;
      TString refName;
      while ( ( ssubresult = static_cast<TObjString*>(nextSubResult()) ) )
      {
        AliDebug(1,Form("-----SubResults %s",ssubresult->String().Data()));

        if(j==1)
        {
          href    = fSpectra->Plot(sparam->String().Data(),ssubresult->String().Data(),kFALSE);
          if  ( !fSpectraName.Contains("VS") ) href = fSpectra->Plot(sparam->String().Data(),ssubresult->String().Data(),kFALSE);
          else if ( fSpectraName.Contains("YVSPT") )
            href = static_cast<TH2*>(fSpectra->Plot(sparam->String().Data(),ssubresult->String().Data(),kFALSE))->ProjectionX();
          else if ( fSpectraName.Contains("PTVSY") )
            href = static_cast<TH2*>(fSpectra->Plot(sparam->String().Data(),ssubresult->String().Data(),kFALSE))->ProjectionX();

          refName = href->GetName();
          j++;
          continue;
        }

        if  ( !fSpectraName.Contains("VS") )
          h = fSpectra->Plot(sparam->String().Data(),ssubresult->String().Data(),kFALSE);
        else if ( fSpectraName.Contains("YVSPT") )
          h = static_cast<TH2*>(fSpectra->Plot(sparam->String().Data(),ssubresult->String().Data(),kFALSE))->ProjectionX();
        else if ( fSpectraName.Contains("PTVSY") )
          h = static_cast<TH2*>(fSpectra->Plot(sparam->String().Data(),ssubresult->String().Data(),kFALSE))->ProjectionX();

        if(!h || !href )
        {
          AliError(Form("Cannot find histos for SubResults  ratio "));
          return;
        }

        if( j!=3 && j!=5 && j!= 10 && j!=11 && j!=12 && j!=13 && j!=14 ) h->SetMarkerColor(j); //nobody likes green and yellow
        else               h->SetMarkerColor(j+5);
        h->SetMarkerSize(1.);
        h->SetMarkerStyle(20+j);
        if(j==11)h->SetMarkerStyle(20+j+3);
        if(j==2)
        {
          h->GetYaxis()->SetTitleSize(0.05);
          h->GetYaxis()->SetLabelSize(0.05);
          h->GetXaxis()->SetLabelSize(0.05);
          h->GetXaxis()->SetTitleSize(0.05);
          h->SetTitle(Form(" %s Ratio over %s for %s",sparam->String().Data(),refName.Data(),fSpectraName.Data()));
          if (fSpectraName.Contains("YVSPT") ) h->GetXaxis()->SetTitle("PT");
          else if (fSpectraName.Contains("PTVSY") ) h->GetXaxis()->SetTitle("Y");
        }
        h->Divide(href);
        if(! sparam->String().Contains("FitChi2PerNDF"))
        {
          if(j==2)h->DrawCopy();
          else    h->DrawCopy("same");
        }
        else
        {
          if(j==2)h->DrawCopy("p");
          else    h->DrawCopy("samep");
        }
        leg2->AddEntry(h,Form("Results %d over Result 1",j),"pe");

        j++;
      }

      leg2->Draw("same");
      ++k;
    } else {
      nextSubResult.Reset();
      int i= 1;
      while ( ( ssubresult = static_cast<TObjString*>(nextSubResult()) ) ){
        //Loop over subresults
        AliDebug(1,Form("-----SubResults %s",ssubresult->String().Data()));
        h = fSpectra->Plot(sparam->String().Data(),ssubresult->String().Data(),kFALSE);
        if(!h) {
          AliError(Form("Cannot find histo for SubResults %s",ssubresult->String().Data()));
          return;
        }
      }
    }
  }

  delete paramArray;
  delete subresultArray;


}

//_____________________________________________________________________________
void AliAnalysisMuMuSpectraProcessor::PrintNofWhat(const char* what) const
{
  /// Print what number for each results on terminal.

  //Check point
  if(!GetSpectra() || strcmp(what,"")==1 )
    {
      AliError("No Spectra or no arguments given !");
      return ;
    }

  // Pointers to handle results and subresults and binning
  AliAnalysisMuMuResult    * result;
  AliAnalysisMuMuJpsiResult* subresult;
  AliAnalysisMuMuResult    * sr;
  AliAnalysisMuMuBinning   ::Range* r;

  // Array to store bins for the while loop
  TObjArray * bins=GetSpectra()->Binning()->CreateBinObjArray();// (intrinseque 'new')
  if (!bins)
  {
    AliError(Form("Cannot find bins"));
    return;
  }

  //Counters and Iterator for bin
  Int_t nofResult = 0;
  TIter nextBin(bins);
  nextBin.Reset();

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

    Int_t nofSubResult = 0; // Counter for subresult
    TIter nextSubResult(result->SubResults());// Iterator for subresults
    nextSubResult.Reset();

    //Some variables
    TString  binAsString(r->AsString());// Usefull for the coming loop

    // To store subresults values
    Double_t subNofWhat[result->SubResults()->GetEntries()];
    Double_t subNofWhatStatError[result->SubResults()->GetEntries()];
    const char * srName[result->SubResults()->GetEntries()];

    cout << Form(" -_-_-_-_- %s --- %s -_-_-_-_- ",binAsString.Data(),fSpectraName.Data()) << endl;

    int excludedResults =0;

    // --- Loop on subresults ---
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

      //Get quantities
      Double_t NofWhat                  = subresult->GetValue(what);
      Double_t NofWhatErrorStat         = subresult->GetErrorStat(what);

      subNofWhat[nofSubResult]          = NofWhat;
      subNofWhatStatError[nofSubResult] = NofWhatErrorStat;
      srName[nofSubResult]              = sr->GetName();

      //Output messages
      cout << Form(" -------- ") << endl;
      cout << Form(" -- subresult %s :  %.3f +/- %.3f ",sr->GetName(),NofWhat,NofWhatErrorStat) << endl;

      // To check the status
      Int_t fitStatus = subresult->HasValue("FitResult") ? subresult->GetValue("FitResult") : 0;
      Int_t covStatus = subresult->HasValue("CovMatrixStatus") ? subresult->GetValue("CovMatrixStatus") : 3;
      Int_t chi2      = subresult->HasValue("FitChi2PerNDF") ? subresult->GetValue("FitChi2PerNDF") : 1;
      if ( (fitStatus!=0 && fitStatus!=4000) || chi2 > 2.5 /*|| covStatus!=3*/ ){
          printf("Fit most likely excluded, you can check in AliAnalysisMuMuResults (FitResult = %d | Cov. Mat. = %d | chi2 = %d)\n",fitStatus,covStatus,chi2);
          ++excludedResults;
      }

      nofSubResult++;
    }

    // Plot the histograms
    TH1F * h_test = new TH1F(Form("%s_%s",what,r->AsString().Data()),Form("%s_%s",what,r->AsString().Data()),result->SubResults()->GetEntries(),0,result->SubResults()->GetEntries());
    for (int i = 0; i < result->SubResults()->GetEntries(); ++i)
    {
        // the histo we plot
        h_test->SetBinContent(i+1,subNofWhat[i]);
        h_test->SetBinError(i+1,subNofWhatStatError[i]);

        // Here we change the label names
        h_test->GetXaxis()->SetBinLabel(i+1,Form("%s",srName[i]));
    }

    // --- Here we draw ---
    TCanvas *se = new TCanvas;
    h_test->DrawCopy();

    TLine *line1 = new TLine(0,result->GetValue(what),result->SubResults()->GetEntries(),result->GetValue(what));
    line1->SetLineColor(kBlue);
    line1->SetLineWidth(3);

    TLine *line2 = new TLine(0,result->GetValue(what)-result->GetRMS(what),result->SubResults()->GetEntries(),result->GetValue(what)-result->GetRMS(what));
    line2->SetLineColor(kBlue);
    line2->SetLineWidth(3);
    line2->SetLineStyle(3);

    TLine *line3 = new TLine(0,result->GetValue(what)+result->GetRMS(what),result->SubResults()->GetEntries(),result->GetValue(what)+result->GetRMS(what));
    line3->SetLineColor(kBlue);
    line3->SetLineWidth(3);
    line3->SetLineStyle(3);
    line1->Draw("same");
    line2->Draw("same");
    line3->Draw("same");

    TLegend* leg = new TLegend(0.2772713,0.6872277,0.5569771,0.8642429,NULL,"brNDC");
    leg->AddEntry("NULL",Form("number of excluded tests : %d",excludedResults),"P");
    leg->Draw("same");



    cout << Form(" -------- ") << endl;
    cout << Form(" ------ Mean :  %.0f +/- %.0f (%.1f %%) +/- %.0f (%.1f %%) ------ ",
      result->GetValue(what),result->GetErrorStat(what),100*result->GetErrorStat(what)/result->GetValue(what),result->GetRMS(what),100*result->GetRMS(what)/result->GetValue(what))
    << " with " <<excludedResults << " excluded results" << endl;
    cout << "" << endl;
    nofResult++;
  }
}
