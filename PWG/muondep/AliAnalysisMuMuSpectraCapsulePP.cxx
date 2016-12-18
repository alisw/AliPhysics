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

///
/// Class to deal with Spectra after the fitting procedure.
///
/// author: Benjamin Audurier (Subatech)
///

#include "AliAnalysisMuMuSpectraCapsulePP.h"

ClassImp(AliAnalysisMuMuSpectraCapsulePP)

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
#include "AliAnalysisMuMuSpectraCapsulePP.h"
#include <fstream>
#include <string>

using std::cout;
using std::endl;
using std::ifstream;


namespace
{
  const Double_t BR           = 5.96/100; // Branching ratio
  const Double_t BRerr        = 0.03/5.96;   // Branching ratio
  //luminosity
  const Double_t lumi         = 106.28; // nb^-1
  const Double_t lumiStat     = 0.09;   // nb^-1
  const Double_t lumiSys      = 2.1/100; // (%)
  const Double_t Trigg        = 1.5/100; // (%)

}


//_____________________________________________________________________________
 AliAnalysisMuMuSpectraCapsulePP::AliAnalysisMuMuSpectraCapsulePP(
const AliAnalysisMuMuSpectra*  spectra,
const TString                 spectraPath,
const char                  * externFile,
const char                  * externFile2)
:
  AliAnalysisMuMuSpectraCapsule(),
  fSpectra(spectra),
  fSpectraName(spectraPath),
  fExternFile(externFile),
  fExternFile2(externFile2),
  fPrintFlag(kFALSE)
{
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

  if(!AliAnalysisMuMuSpectraCapsule::SetConstantFromExternFile(fExternFile2,&fConstArray[0],&fSpectraName))
  {
    AliWarning(Form("No extern file readed"));
  }
}

//_____________________________________________________________________________
AliAnalysisMuMuSpectraCapsulePP::~AliAnalysisMuMuSpectraCapsulePP()
{
  // dtor
}

TList* AliAnalysisMuMuSpectraCapsulePP::ComputePPCrossSection(const char* what) const
{
  /// Compute the PP cross section
  /// Warning : the cross-section is normalized to bin width only if bin width > 2


  AliAnalysisMuMuBinning* binning = fSpectra->Binning();
  TObjArray* bins = binning->CreateBinObjArray();
  TIter nextBin(bins);
  AliAnalysisMuMuBinning::Range* bin;
  Int_t i(0);
  AliAnalysisMuMuResult* r;

  const Double_t * binArrayX = binning->CreateBinArrayX();
  Int_t nBinX = binning->GetNBinsX();

  TGraphErrors * gCrossSection = new TGraphErrors(nBinX);
  TGraphErrors * gSys = new TGraphErrors(nBinX);
  TString sbin;

  int j= 1;
  while ( ( bin = static_cast<AliAnalysisMuMuBinning::Range*>(nextBin()) ) )
  {
      r = static_cast<AliAnalysisMuMuResult*>(fSpectra->GetResultForBin(bin->AsString()));
      if(j==1 && bin->AsString().Contains("PT"))      sbin ="PT";
      else if(j==1 && bin->AsString().Contains("Y"))  sbin ="Y";

      // PT and Y bin
      if(sbin.Contains("PT") || sbin.Contains("Y")){

        // read exterfile and get the correct value
        float valueArray[4];
        //  valueArray[0],   valueArray[1],   valueArray[2],   valueArray[3]
        //  sysMC(%)         TrajEffError(%)  TriggerError(%)  matchingError(%)
        if(ReadFromFile(bin->AsString(),&valueArray[0])==kFALSE) return 0x0;
        AliDebug(1, " Values correctly read from extern file");
        //Define quantities
        Double_t sigma          =0.;
        Double_t sigmaerrorstat =0.;

        Double_t CorrNofJPsi    =0.;
        Double_t NofJPsi        =0.;
        Double_t NofJPsiError   =0.;
        Double_t AccEff         =0.;
        // Get Jpsi numbers
        CorrNofJPsi  =r->GetValue(what);
        NofJPsi      =r->GetValue("NofJPsi");
        NofJPsiError =r->GetErrorStat("NofJPsi");

        if(fPrintFlag){
          printf("\n");
          printf("%s                        = %f +/- %f\n",what,r->GetValue(what),r->GetErrorStat(what));
          printf("NofJPsi                   = %f +/- %f (%f %%)+/- %f (%f %%)\n ",NofJPsi,NofJPsiError,100*NofJPsiError/NofJPsi,r->GetRMS("NofJPsi"),100*r->GetRMS("NofJPsi")/NofJPsi);
          printf("Systematic MC             = %f \n ",valueArray[0]);
          printf("Tracking Error            = %f \n ",valueArray[1]);
          printf("Trigger  Error            = %f \n ",valueArray[2]);
          printf("matching Error            = %f \n ",valueArray[3]);
          printf("Trigger Error local board = %f %%\n ",100*Trigg);
          printf("Lumiosity                 = %f +/- %f  +/- %f %% \n ",lumi,lumiStat,lumiSys);
          printf("BR                        = %f +/- %f (%%)\n ",BR,100*BRerr);
          printf("\n");
        }


        if(CorrNofJPsi==0. || NofJPsi ==0.){
          printf(" cannot found Corrected NofJpsi or NofJpsi, did you compute AccEff ? Abording...");
          continue;
        }

        // Select Delta y according to bin
        Double_t deltaY =0.;
        if(sbin.Contains("PT") && bin->WidthX() <= 2.0 )deltaY = 1.5; // For pT_0_1,pT_1_2...Y_4_3.75...
        else  deltaY =1.;

        // Compute cross section
        if(bin->WidthX() <= 2.0)sigma = (r->GetValue(what))/(lumi*BR*1000.*bin->WidthX()); // For pT_0_1,pT_1_2...Y_4_3.75...
        else                    sigma = (r->GetValue(what))/(lumi*BR*1000.);               // For pT_0_8,pT_0_12...Y_4_2.5...


        // Compute stat. error on cross section
        sigmaerrorstat = sigma*TMath::Sqrt(
          r->GetErrorStat("NofJPsi")*r->GetErrorStat("NofJPsi")/r->GetValue("NofJPsi")/r->GetValue("NofJPsi"));
        //                                    Signal stat.

         // Get X
        Double_t xmin = bin->Xmin();
        Double_t xmax = bin->Xmax();
        Double_t x    = xmin + (xmax-xmin)/2;

        // Uncorrelated error squared (%)
        Double_t UncError2 =
          (r->GetRMS("NofJPsi")/r->GetValue("NofJPsi"))*(r->GetRMS("NofJPsi")/r->GetValue("NofJPsi")) // Signal
          + pow(valueArray[0]/100.,2)  // systematic MC(%)
          + pow(valueArray[1]/100.,2)  // Tracking Error(%)
          + pow(valueArray[2]/100.,2)  // Trigger  Error(%) (trigger response)
          + pow(valueArray[3]/100.,2)  // matching Error(%)
          + Trigg*Trigg;               // Trigger  Error (plateau)


        // Correlated squared (%)
        Double_t CorrError2 = lumiSys*lumiSys + lumiStat*lumiStat/lumi/lumi +BRerr*BRerr ;
        //                        lumi sys. (%)    lumi stat.                 BR (%)

        // In case of fully integrated results == large pT bins
        Double_t CorrErrorFull2 =  CorrError2 + UncError2;

        printf("  -------- cross section for bin %s = %f +/- %f (stat. %f %%) +/- %f (uncorr. %f %%) +/- %f (corr. %f %%)  (#Delta y = %f) -------- \n"
          ,bin->AsString().Data(),
          sigma/deltaY,
          sigmaerrorstat/deltaY,
          100*sigmaerrorstat/sigma,
          sigma*TMath::Sqrt(UncError2)/deltaY,
          100*TMath::Sqrt(UncError2),
          sigma*TMath::Sqrt(CorrError2)/deltaY,
          100*TMath::Sqrt(CorrError2),
          deltaY);

        printf("  --------  if fully integrated dsigma_corr. = %f (corrFull. %f %%)  -------- \n\n"
          ,sigma*TMath::Sqrt(CorrErrorFull2)/deltaY,
          100*TMath::Sqrt(CorrErrorFull2));

        // Fill graphs
        if(sbin.Contains("Y"))gCrossSection->SetPoint(j,-x,sigma/deltaY);
        else gCrossSection->SetPoint(j,x,sigma/deltaY);
        gCrossSection->SetPointError(j,bin->WidthX()/2,sigmaerrorstat/deltaY);

        if(sbin.Contains("Y"))gSys->SetPoint(j,-x,sigma/deltaY);
        else gSys->SetPoint(j,x,sigma/deltaY);
        gSys->SetPointError(j,bin->WidthX()/4,sigma*TMath::Sqrt(UncError2)/deltaY);


        j++;
      }
      else { // Integrated

        Double_t sigma        =0.;
        Double_t sigmaerror   =0.;
        Double_t sigmasys2    =0.;

        Double_t CorrNofJPsi  =0.;
        Double_t NofJPsi      =0.;
        Double_t NofJPsiError =0.;
        Double_t NofJPsiSys   =0.;

        NofJPsi               = fSpectra->GetResultForBin("INTEGRATED")->GetValue("NofJPsi");
        CorrNofJPsi           = fSpectra->GetResultForBin("INTEGRATED")->GetValue(what);

        NofJPsiSys            = fSpectra->GetResultForBin("INTEGRATED")->GetRMS("NofJPsi");
        NofJPsiError          = fSpectra->GetResultForBin("INTEGRATED")->GetErrorStat("NofJPsi");

        AliDebug(1,Form(""));
        AliDebug(1,Form("%s      = %f +/- %f\n",what,r->GetValue(what),r->GetErrorStat(what)));
        AliDebug(1,Form("NofJPsi = %f +/- %f +/- %f \n ",NofJPsi,NofJPsiError,r->GetRMS("NofJPsi")));
        AliDebug(1,Form(""));

        if(CorrNofJPsi==0. || NofJPsi ==0.|| NofJPsiSys ==0.|| NofJPsiError ==0.){
          printf(" cannot found Corrected NofJpsi or NofJpsi, did you compute AccEff ? Abording...");
          continue;
        }

        sigma      = CorrNofJPsi/(fConstArray[0]*BR*1000.);
        sigmaerror = sigma*TMath::Sqrt(
          r->GetErrorStat("NofJPsi")*r->GetErrorStat("NofJPsi")/r->GetValue("NofJPsi")/r->GetValue("NofJPsi"));
        //                                    Signal stat.
        sigmasys2 =
        NofJPsiSys/NofJPsi*NofJPsiSys/NofJPsi // Signal extraction
        +pow(lumiStat/lumi,2)                 // Lumi stat.
        +pow(fConstArray[2]/100.,2)           // Lumi syst.
        +pow(fConstArray[3]/100.,2)           // Trigger
        +pow(fConstArray[4]/100.,2)           // Trigger Local board
        +pow(fConstArray[5]/100.,2)           // Tracking
        +pow(fConstArray[6]/100.,2)           // MC Input
        +pow(fConstArray[7]/100.,2)           // Matching
        +BRerr*BRerr ;                        // BR


        if(fPrintFlag){
          printf("\n");
          printf("%s                        = %f +/- %f\n",what,r->GetValue(what),r->GetErrorStat(what));
          printf("NofJPsi                   = %f +/- %f +/- %f \n ",NofJPsi,NofJPsiError,r->GetRMS("NofJPsi"));
          printf("Systematic MC             = %f \n ",100*fConstArray[6]);
          printf("Tracking Error            = %f \n ",100*fConstArray[5]);
          printf("Trigger  Error            = %f \n ",100*fConstArray[4]);
          printf("matching Error            = %f \n ",100*fConstArray[7]);
          printf("Trigger Error local board = %f %%\n ",100*fConstArray[4]);
          printf("Lumiosity                 = %f +/- %f  +/- %f %% \n ",lumi,lumiStat,100*fConstArray[2]);
          printf("BR                        = %f +/- %f (%%)\n ",BR,100*BRerr);
          printf("\n");
        }

        printf("integrated cross section for  %s = %f +/- %f +/- %f #mubarn\n",fSpectra->GetName(),sigma,sigmaerror,sigma*TMath::Sqrt(sigmasys2) );
        return 0x0;
      }
  }


  //Add and merge all Graph
  TList* l = new TList();
  l->SetOwner(kTRUE);
  l->Add(gCrossSection);
  l->Add(gSys);

  delete binArrayX;

  return l ;
}

//_____________________________________________________________________________
TGraphErrors* AliAnalysisMuMuSpectraCapsulePP::ComputeYield( const char* what, const TH1* histo, const char* sResName)
{
  /// Compute Yield.
  /// Arguments :
  ///   - what : the yield nominator, i.e NofJPsi, meanPT etc. (null by default)
  ///   - histo : histogramme of Equivalent MinBias

 printf("Not implemented yet \n");
 return 0x0;

}

//_____________________________________________________________________________
void AliAnalysisMuMuSpectraCapsulePP::DrawResults( const char* particle,const char* subresults) const
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
  c->Divide(nx,ny,0,0);
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
      gPad->SetLogy();
      if (xmin>0)
        {
          // Loop to configure the pad as you like
          h->GetXaxis()->SetRangeUser(xmin,xmax);
          h->SetTitleSize(10);
        }
      h->GetXaxis()->SetRangeUser(1.,5.);
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
      f1 = h->GetFunction("signal");
      if(f1) f1->DrawCopy("same");
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
void AliAnalysisMuMuSpectraCapsulePP::Print(Option_t* opt) const
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
void AliAnalysisMuMuSpectraCapsulePP::PrintConst() const
{
    ///
    /// Print member constants on the terminal
    ///

  //Check point
  if(!GetSpectra()) return ;
  else{
      cout <<      " ================================================================ " << endl;
      cout << Form("      Constants for Spectra %s",fSpectraName.Data()) << endl;
      cout <<      " ================================================================ " << endl;
      cout <<  "   Branching ratio = " <<  5.96/100 << endl; //

      cout << Form(" -- Value of Lumi                       = %f",fConstArray[0]) << endl;
      cout << Form(" -- Value of Lumi stat.                 = %f",fConstArray[1]) << endl;
      cout << Form(" -- Value of Lumi syst (%%)              = %f",fConstArray[2]) << endl;
      cout << Form(" -- Value of Trigg syst (%%)             = %f",fConstArray[3]) << endl;
      cout << Form(" -- Value of Trigg Local board (%%)      = %f",fConstArray[4]) << endl;
      cout << Form(" -- Value of Traj. err.                 = %f",fConstArray[5]) << endl;
      cout << Form(" -- Value of MC Input. err.             = %f",fConstArray[6]) << endl;
      cout << Form(" -- Value of Matching. err.             = %f",fConstArray[7]) << endl;
      cout << Form(" -- Value of AccEff                     = %f",fConstArray[8]) << endl;
      cout << Form(" -- Value of dAccEff                    = %f",fConstArray[9]) << endl;
      cout << Form(" -- Value of NofJpsi     from exterfile = %f",fConstArray[10]) << endl;
      cout << Form(" -- Value of StatNofJpsi from exterfile = %f",fConstArray[11]) << endl;
      cout << Form(" -- Value of SystNofJpsi from exterfile = %f",fConstArray[12]) << endl;
  }
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuSpectraCapsulePP::ReadFromFile(TString sbin, float valueArray[]) const
{
    ///
    /// Read extern file lines and store associated values. Exemple of line :
    /// #intervalLow intervalHight  sysMC   TrajEffError(%)  TriggerError(%)  matchingError(%)
    /// 00           01             1.600    4.0              4.6             1.0
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
        AliDebug(1,Form("intervalLow = %s\n", intervalLow.Data()));
        AliDebug(1,Form("intervalHigh = %s\n", intervalHigh.Data()));
        AliDebug(1,Form("sbin = %s\n", sbin.Data()));

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
        for (int i =0 ; i<4 ; i++) {
            valueArray[i]= static_cast<TObjString*>(lineArray->At(i+2))->String().Atof();
        }
        return ok;
    }
    else return ok;
}
