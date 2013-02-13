/* $Id$ */

#include "TFile.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TNamed.h"
#include "TF1.h"
#include "TAttLine.h"
#include "TAttMarker.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TMath.h"
#include "TList.h"
#include <TString.h>
#include <TObject.h>
#include <TROOT.h>
#include <THashList.h>
#include <Rtypes.h>
#include <TPRegexp.h>
#include "TFitResult.h"
#include <TMap.h>
#include <TObjString.h>
#include "TPad.h"
#include "TAxis.h"

namespace RawProduction {

  bool ignoreErrors = false;

  // Fit range
  Double_t rangeMin=0.05 ;
  Double_t rangeMax=0.3 ;

  // Fit limits
  double lowerMass = 0.110; double upperMass = 0.145;
  double lowerWidth = 0.001; double upperWidth = 0.012;

  Double_t PeakPosition(Double_t pt){
    //Fit to the measured peak position
    return 4.99292e-003*exp(-pt*9.32300e-001)+1.34944e-001 ;
  }
  //-----------------------------------------------------------------------------
  Double_t PeakWidth(Double_t pt){
    //fit to the measured peak width
    Double_t a=0.0068, b=0.0025, c=0.000319 ;
    return TMath::Sqrt(a*a+b*b/pt/pt+c*c*pt*pt) ;
  }

  // Pt bin parameters
  Int_t nPtBins=0;
  Double_t ptBinEdges[1000] = {0};
  int fool;
  double GetPtBin(int bin);
  void MakePtBins();

  // various parameters
  const Double_t kMean=0.136 ; //Approximate peak position to facilitate error estimate
  const char format[] = ".pdf"; // say if you want .pdf

  class Input;

  TH1* GetHistogram_cent(Input& input, const TString& name, int centrality);
  TH1* MergeHistogram_cent(Input& input, const TString& name, int newCentIndex, int fromCentIndex, int toCentIndex);
  int kNCentralityBins = 3;

  Double_t CGausPol1(Double_t * x, Double_t * par);
  Double_t CGausPol2(Double_t * x, Double_t * par);
  Double_t CGausPol0(Double_t * x, Double_t * par);
  Double_t CPol1(Double_t * x, Double_t * par);
  Double_t CPol2(Double_t * x, Double_t * par);


  // Output Bin
  class TriggerBin {
  public:
    TriggerBin(const TString& trigger = "kCentral" );
    TriggerBin(const char* trigger);
    const TString& Key() const {return fKey;}
    const TString& Trigger() const {return fTrigger;}
  protected:
    TString fTrigger;
    TString fKey;
  };
  // Object describing the input for the macros
  class Input {
  public:
    Input(const TString& fileName, const TriggerBin& inputBin, TString listPath = "");
    TH1 * GetHistogram(const char* name="");
    const TriggerBin& Bin() const {return fBin;}
  private:
    static TFile* fFile;
    TList* fList;
    TriggerBin fBin;
  };

  // Output Bin
  class TriCenPidBin : public TriggerBin {
  public:
    TriCenPidBin(Int_t centrality, const TString& pid, const TString& trigger);
    Int_t Centrality() const {return fCentrality;}
    const TString& PID() const {return fPID;}
  private:
    Int_t fCentrality;
    TString fPID;
  };
  // Object describing the output of the macros
  class Output {
  public:
    Output(const TString& fileName = "RawProduction.root", const char* options = "UPDATE");
    TH1* GetHistogram(const TString& name, const TriggerBin& inBin);
    void SetDir(const TriggerBin& inBin);
    void Write();
  private:
    TFile* fFile;
  };

  void MakePi0Fit(Input& input, const TriCenPidBin& outBin, Output& output)
  {
    MakePtBins();
    Printf("\nMakePi0Fit(%s)", outBin.Key().Data());
    output.SetDir(outBin);

    TH1F * hTotSelEvents          = (TH1F*) input.GetHistogram("hTotSelEvents");
    TH2F * hCentrality  = (TH2F*) input.GetHistogram("hCenPHOSCells");
    TH1D * hCentralityX = hCentrality->ProjectionX();
    TH2F *hPi0 =    (TH2F*)GetHistogram_cent(input, Form("hPi0%s", outBin.PID().Data()), outBin.Centrality());
    TH2F *hPi0Mix = (TH2F*)GetHistogram_cent(input, Form("hMiPi0%s", outBin.PID().Data()), outBin.Centrality());

    printf("TotSelEvents (4): %.0f \n", hTotSelEvents->GetBinContent(4)) ;
    printf("Centrality:   %.0f \n",     hCentralityX->Integral()) ;

    if( !hPi0 || !hPi0Mix ) {
      Printf(Form("no histogram(0x%p, 0x%p), returning", hPi0, hPi0Mix));
      return;
    }

    // for temp convas for drawing/monitoring
    TCanvas* canvas = new TCanvas("cMakePi0Fit", Form("MakePi0Fit Canvas, %s", outBin.Key().Data()),10,10,1200,800);


    // Peak Parameterization
    //  1. Linear Bg
    TF1 * funcRatioFit1 = new TF1("funcRatioFit1",CGausPol1,0.,1.,5) ;
    funcRatioFit1->SetParNames("A", "m_{0}", "#sigma", "a_{0}", "a_{1}");
    funcRatioFit1->SetLineWidth(2) ;
    funcRatioFit1->SetLineColor(2) ;
    //  2. Quadratic Bg
    TF1 * funcRatioFit2 = new TF1("funcRatioFit2",CGausPol2,0.,1.,6) ;
    funcRatioFit2->SetParNames("A", "m_{0}", "#sigma", "a_{0}", "a_{1}", "a_{2}");
    funcRatioFit2->SetLineWidth(2) ;
    funcRatioFit2->SetLineColor(4) ;
    funcRatioFit2->SetLineStyle(2) ;
    //     Other
    TF1 * fgs = new TF1("gs",CGausPol0,0.,1.,4) ;
    fgs->SetParNames("A", "m_{0}", "#sigma", "B") ;
    fgs->SetLineColor(2) ;
    fgs->SetLineWidth(1) ;
    TF1 * fbg1 = new TF1("bg1",CPol1,0.,1.,2) ;
    TF1 * fbg2 = new TF1("bg2",CPol2,0.,1.,3) ;

    // Adding Histograms:
    //  1. Linear Bg
    TStringToken names("mr1;mr1r;sr1;sr1r;ar1;br1;yr1;yr1int", ";");
    TStringToken titles("Mass;Mass, Ratio Fit;Width;Width, Ratio Fit;a;b;Raw Yield; Raw Yield, integrated", ";");
    while( names.NextToken() && titles.NextToken() ) {
      new TH1D(names.Data(), titles.Data(), nPtBins,ptBinEdges);
      new TH1D(Form("%s_error", names.Data()), titles.Data(), nPtBins,ptBinEdges);
    }
    //  2. Quadratic Bg
    TStringToken names2("mr2;mr2r;sr2;sr2r;ar2;br2;cr2;yr2;yr2int", ";");
    TStringToken titles2("Mass;Mass, Ratio Fit;Width;Width, Ratio Fit;a;b;c;Raw Yield; Raw Yield, integrated", ";");
    while( names2.NextToken() && titles2.NextToken() ) {
      new TH1D(names2.Data(), titles2.Data(), nPtBins,ptBinEdges);
      new TH1D(Form("%s_error", names2.Data()), titles2.Data(), nPtBins,ptBinEdges);
    }

    TH1D* hMixRawRatio =new TH1D("hMixRawRatio", "ratio of statistics in RAW and Mixed", nPtBins, ptBinEdges);

    // Pt slice loop
    for(Int_t ptBin=1; ptBin<=nPtBins; ptBin++){
      //gSystem->Sleep(2000);
      canvas->Clear();
      canvas->Divide(2,2);
      canvas->cd(1);
      printf("pt bin %i, %3.1f to %3.1f, ", ptBin, ptBinEdges[ptBin-1], ptBinEdges[ptBin]);

      //TList* ptBinOutputList = outputList;
      TAxis * ptAxis=hPi0->GetYaxis() ;
      Int_t imin=ptAxis->FindBin(ptBinEdges[ptBin-1]+0.0001);
      Int_t imax=ptAxis->FindBin(ptBinEdges[ptBin]+0.0001) -1;
      if( TMath::Abs( ptAxis->GetBinLowEdge(imin) - ptBinEdges[ptBin -1] ) >0.0001 ) {
	Printf("\nError: hPi0 lower bin edge (%f) different from ptBinEdges lower (%f)", ptAxis->GetBinLowEdge(imin), ptBinEdges[ptBin -1]);
	continue;
      }
      if( TMath::Abs( ptAxis->GetBinLowEdge(imax+1) - ptBinEdges[ptBin] ) >0.0001 ) {
	Printf("\nError: hPi0 upper bin edge (%f) different from ptBinEdges upper (%f)", ptAxis->GetBinLowEdge(imax+1), ptBinEdges[ptBin]);
	continue;
      }

      Double_t pt=(ptBinEdges[ptBin]+ptBinEdges[ptBin-1])/2. ;
      Double_t dpt= ptBinEdges[ptBin] - ptBinEdges[ptBin-1];

      TH1D * hPi0Proj = hPi0->ProjectionX(Form("pt%03i_hPi0",ptBin), imin, imax);
      hPi0Proj->Sumw2();
      TH1D * hPi0ProjMix= hPi0Mix->ProjectionX(Form("pt%03i_hPi0Mix",ptBin), imin, imax) ;
      hPi0ProjMix->Sumw2();

      const Int_t pi0Entries = hPi0Proj->Integral(hPi0Proj->FindBin(lowerMass), hPi0Proj->FindBin(upperMass));
      const Int_t mixEntries = hPi0ProjMix->Integral(hPi0Proj->FindBin(lowerMass), hPi0Proj->FindBin(upperMass));
      if( pi0Entries >0)  {
	hMixRawRatio->SetBinContent(ptBin, mixEntries/pi0Entries);
	hMixRawRatio->SetBinError(ptBin, TMath::Sqrt(mixEntries/pi0Entries/pi0Entries + mixEntries*mixEntries/pi0Entries/pi0Entries/pi0Entries));
      }
      printf("statistics in bin is %i, mixed %i, ", pi0Entries, mixEntries);
      if( pi0Entries < 10 ) {
	Printf("to few entries");
	continue;
      }

      const bool rebin = pi0Entries < 1000;
      if( rebin && pi0Entries >= 500 ) {
	printf("rebin by factor 2, ");
	hPi0Proj->Rebin(2);
	hPi0ProjMix->Rebin(2);
      } else if( rebin && pi0Entries >= 200) {
	printf("rebin by factor 4, ");
	hPi0Proj->Rebin(4);
	hPi0ProjMix->Rebin(4);
      } else if( rebin && pi0Entries >= 10) {
	printf("rebin by factor 5, ");
	hPi0Proj->Rebin(5);
	hPi0ProjMix->Rebin(5);
      }
      std::cout << std::endl;

      hPi0Proj->SetTitle(Form("M_{#gamma#gamma}, %.1f<p_{T}<%.1f, PID=%s, %s", ptBinEdges[ptBin-1], ptBinEdges[ptBin], outBin.PID().Data(), outBin.Trigger().Data()));
      hPi0Proj->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
      hPi0ProjMix->SetTitle(Form("M_{#gamma#gamma}^{Mix}, %.1f<p_{T}<%.1f, PID=%s, %s",ptBinEdges[ptBin-1],ptBinEdges[ptBin],outBin.PID().Data(), outBin.Trigger().Data()));
      hPi0ProjMix->SetXTitle("M_{#gamma#gamma}^{Mix} (GeV/c^{2})");


      // Error Fix
      for(Int_t ib=1; ib<=hPi0Proj->GetNbinsX();ib++){if(hPi0Proj ->GetBinContent(ib)==0)hPi0Proj ->SetBinError(ib,1.);}
      for(Int_t ib=1; ib<=hPi0Proj->GetNbinsX();ib++){if(hPi0ProjMix->GetBinContent(ib)==0)hPi0ProjMix->SetBinError(ib,1.);}

      // Signal-Mix Ratio
      canvas->cd(2);
      TH1D * hPi0Ratio = (TH1D*)hPi0Proj->Clone( Form("pt%03i_hPi0Ratio",ptBin) ) ;
      hPi0Ratio->SetTitle(Form("#frac{M_{#gamma#gamma}}{M_{#gamma#gamma}^{Mix}}, %.1f<p_{T}<%.1f GeV/c", ptBinEdges[ptBin-1], ptBinEdges[ptBin]));
      hPi0Ratio->Divide(hPi0ProjMix) ;
      hPi0Ratio->SetMarkerStyle(20) ;
      hPi0Ratio->SetMarkerSize(0.7) ;

      //     if(outBin.Centrality()==0) rangeMax=0.4 ;
      //     if(ptBin==1){
      //       rangeMin=0.06 ;
      //       rangeMax=0.25 ;
      //     }



      // ================================================
      // Fit Pol1 ratio
      // ================================================
      printf("Pol1 ratio Fit, ");
      canvas->cd(2);

      funcRatioFit1->SetParameters(0.001,0.136,0.0055,0.0002,-0.002) ;
      funcRatioFit1->SetParLimits(0,0.000,1.000) ;
      funcRatioFit1->SetParLimits(1,lowerMass,upperMass) ;
      funcRatioFit1->SetParLimits(2,lowerWidth,upperWidth) ;


      TFitResultPtr ratioFitResultPtr1 = hPi0Ratio->Fit(funcRatioFit1,"MSQ" ,"",rangeMin,rangeMax) ;
      if( int(ratioFitResultPtr1) % 4000 ) // "More" error is acceptable
	ratioFitResultPtr1 = hPi0Ratio->Fit(funcRatioFit1,"MSQ" ,"",rangeMin,rangeMax) ;

      Int_t ratioFitError1 = ratioFitResultPtr1;
      ratioFitError1 = ratioFitError1 % 4000; // "More" error is acceptable
      ratioFitError1 = ratioFitError1 || TMath::Abs( funcRatioFit1->GetParameter(1) - lowerMass) < 0.0001 || TMath::Abs( funcRatioFit1->GetParameter(1) - upperMass) < 0.0001;//  "center" mass converged to limit
      ratioFitError1 = ratioFitError1 || funcRatioFit1->GetParError(1) > (upperMass - lowerMass)/2; // to large of an error
      ratioFitError1 = ratioFitError1 || TMath::Abs( funcRatioFit1->GetParameter(2) - lowerWidth) < 0.0001 || TMath::Abs( funcRatioFit1->GetParameter(2) - upperWidth) < 0.0001;//  st. error converged to limit
      ratioFitError1 = ratioFitError1 || funcRatioFit1->GetParError(2) > (upperWidth - lowerWidth)/2; // to large of an error

      if( ratioFitError1) {
	Printf("in ERROR, %i", ratioFitError1);
	((TH1D*) output.GetHistogram("mr1r_error", outBin))->SetBinContent(ptBin,fgs->GetParameter(1)) ;
	((TH1D*) output.GetHistogram("mr1r_error", outBin))->SetBinError  (ptBin,fgs->GetParError(1) ) ;
	((TH1D*) output.GetHistogram("sr1r_error", outBin))->SetBinContent(ptBin,TMath::Abs(fgs->GetParameter(2))) ;
	((TH1D*) output.GetHistogram("sr1r_error", outBin))->SetBinError  (ptBin,fgs->GetParError(2) ) ;
	((TH1D*) output.GetHistogram("ar1_error", outBin))->SetBinContent(ptBin,funcRatioFit1->GetParameter(3)) ;
	((TH1D*) output.GetHistogram("ar1_error", outBin))->SetBinError  (ptBin,funcRatioFit1->GetParError(3)) ;
	((TH1D*) output.GetHistogram("br1_error", outBin))->SetBinContent(ptBin,funcRatioFit1->GetParameter(4)) ;
	((TH1D*) output.GetHistogram("br1_error", outBin))->SetBinError  (ptBin,funcRatioFit1->GetParError(4)) ;
      }
      if( !ratioFitError1 || ignoreErrors ) {
	Printf("converged, status:%i, covMatrixStatus: %i", ratioFitResultPtr1->Status(), ratioFitResultPtr1->CovMatrixStatus());
	((TH1D*) output.GetHistogram("mr1r", outBin))->SetBinContent(ptBin,fgs->GetParameter(1)) ;
	((TH1D*) output.GetHistogram("mr1r", outBin))->SetBinError  (ptBin,fgs->GetParError(1) ) ;
	((TH1D*) output.GetHistogram("sr1r", outBin))->SetBinContent(ptBin,TMath::Abs(fgs->GetParameter(2))) ;
	((TH1D*) output.GetHistogram("sr1r", outBin))->SetBinError  (ptBin,fgs->GetParError(2) ) ;
	((TH1D*) output.GetHistogram("ar1", outBin))->SetBinContent(ptBin,funcRatioFit1->GetParameter(3)) ;
	((TH1D*) output.GetHistogram("ar1", outBin))->SetBinError  (ptBin,funcRatioFit1->GetParError(3)) ;
	((TH1D*) output.GetHistogram("br1", outBin))->SetBinContent(ptBin,funcRatioFit1->GetParameter(4)) ;
	((TH1D*) output.GetHistogram("br1", outBin))->SetBinError  (ptBin,funcRatioFit1->GetParError(4)) ;
      }



      // ================================================
      // Fit Pol2 ratio
      // ================================================
      printf("Pol2 ratio Fit, ");
      if( ratioFitError1 ) {
	funcRatioFit2->SetParameters(0.001,0.136,0.0055,0.0002,-0.002, 0) ;
	funcRatioFit2->SetParLimits(0,0.000,1.000) ;
	funcRatioFit2->SetParLimits(1,lowerMass,upperMass) ;
	funcRatioFit2->SetParLimits(2,lowerWidth,upperWidth) ;
      } else {
	funcRatioFit2->SetParameters(funcRatioFit1->GetParameters()) ;
	funcRatioFit2->SetParameter(5, 0);
	funcRatioFit2->SetParLimits(0,0.000,1.000) ;
	funcRatioFit2->SetParLimits(1,lowerMass,upperMass) ;
	funcRatioFit2->SetParLimits(2,lowerWidth,upperWidth) ;
      }
      TFitResultPtr ratioFitResultPtr2 = hPi0Ratio->Fit(funcRatioFit2,"+MSQ" ,"",rangeMin,rangeMax) ;
      if( int(ratioFitResultPtr2) != 4000 ) // if error, "More" error is acceptable
	ratioFitResultPtr2  = hPi0Ratio->Fit(funcRatioFit2,"MSQ" ,"",rangeMin,rangeMax) ;

      Int_t ratioFitError2 = ratioFitResultPtr2;
      ratioFitError2 = ratioFitError2 % 4000; // "More" error is acceptable
      ratioFitError2 = ratioFitError2 || TMath::Abs( funcRatioFit2->GetParameter(1) - lowerMass) < 0.0001 || TMath::Abs( funcRatioFit2->GetParameter(1) - upperMass) < 0.0001;  // "center" mass converged to limit
      //ratioFitError2 = ratioFitError2 || funcRatioFit2->GetParError(1) > (upperMass - lowerMass)/2; // to large of an error
      ratioFitError2 = ratioFitError2 || TMath::Abs( funcRatioFit2->GetParameter(2) - lowerWidth) < 0.0001 || TMath::Abs( funcRatioFit2->GetParameter(2) - upperWidth) < 0.0001;  // st. error converged to limit
      //ratioFitError2 = ratioFitError2 || funcRatioFit2->GetParError(2) > (upperWidth - lowerWidth)/2; // to large of an error

      if( ratioFitError2) {
	Printf("in ERROR, %i", ratioFitError2);
	((TH1D*) output.GetHistogram("mr2r_error", outBin))->SetBinContent(ptBin,fgs->GetParameter(1)) ;
	((TH1D*) output.GetHistogram("mr2r_error", outBin))->SetBinError  (ptBin,fgs->GetParError(1) ) ;
	((TH1D*) output.GetHistogram("sr2r_error", outBin))->SetBinContent(ptBin,TMath::Abs(fgs->GetParameter(2))) ;
	((TH1D*) output.GetHistogram("sr2r_error", outBin))->SetBinError  (ptBin,fgs->GetParError(2) ) ;
	((TH1D*) output.GetHistogram("ar2_error", outBin))->SetBinContent(ptBin,funcRatioFit2->GetParameter(3)) ;
	((TH1D*) output.GetHistogram("ar2_error", outBin))->SetBinError  (ptBin,funcRatioFit2->GetParError(3)) ;
	((TH1D*) output.GetHistogram("br2_error", outBin))->SetBinContent(ptBin,funcRatioFit2->GetParameter(4)) ;
	((TH1D*) output.GetHistogram("br2_error", outBin))->SetBinError  (ptBin,funcRatioFit2->GetParError(4)) ;
	((TH1D*) output.GetHistogram("cr2_error", outBin))->SetBinContent(ptBin,funcRatioFit2->GetParameter(5)) ;
	((TH1D*) output.GetHistogram("cr2_error", outBin))->SetBinError  (ptBin,funcRatioFit2->GetParError(5)) ;
      }
      if( !ratioFitError2 || ignoreErrors ) {
	Printf("converged, status:%i, covMatrixStatus: %i", ratioFitResultPtr2->Status(), ratioFitResultPtr2->CovMatrixStatus());
	((TH1D*) output.GetHistogram("mr2r", outBin))->SetBinContent(ptBin,fgs->GetParameter(1)) ;
	((TH1D*) output.GetHistogram("mr2r", outBin))->SetBinError  (ptBin,fgs->GetParError(1) ) ;
	((TH1D*) output.GetHistogram("sr2r", outBin))->SetBinContent(ptBin,TMath::Abs(fgs->GetParameter(2))) ;
	((TH1D*) output.GetHistogram("sr2r", outBin))->SetBinError  (ptBin,fgs->GetParError(2) ) ;
	((TH1D*) output.GetHistogram("ar2", outBin))->SetBinContent(ptBin,funcRatioFit2->GetParameter(3)) ;
	((TH1D*) output.GetHistogram("ar2", outBin))->SetBinError  (ptBin,funcRatioFit2->GetParError(3)) ;
	((TH1D*) output.GetHistogram("br2", outBin))->SetBinContent(ptBin,funcRatioFit2->GetParameter(4)) ;
	((TH1D*) output.GetHistogram("br2", outBin))->SetBinError  (ptBin,funcRatioFit2->GetParError(4)) ;
	((TH1D*) output.GetHistogram("cr2", outBin))->SetBinContent(ptBin,funcRatioFit2->GetParameter(5)) ;
	((TH1D*) output.GetHistogram("cr2", outBin))->SetBinError  (ptBin,funcRatioFit2->GetParError(5)) ;
      }



      // ================================================
      // Plot Ratio Fits
      // ================================================
      hPi0Ratio->GetXaxis()->SetRangeUser(rangeMin, rangeMax);
      hPi0Ratio->Draw();
      canvas->Update();



      // ================================================
      // Pol1 Scaled Background Subtraction
      // ================================================
      canvas->cd(3);
      Double_t intRangeMin = PeakPosition(pt)-3.*PeakWidth(pt) ;
      Double_t intRangeMax = PeakPosition(pt)+3.*PeakWidth(pt) ;
      Int_t    intBinMin   = hPi0Proj->GetXaxis()->FindBin(intRangeMin) ;
      Int_t    intBinMax   = hPi0Proj->GetXaxis()->FindBin(intRangeMax) ;
      Double_t mixInt     = hPi0ProjMix->Integral(intBinMin,intBinMax);

      if( ! ratioFitError1 || true) {
	printf("Pol1 BS Fit, ");
	TH1D * hPi0MixScaledPol1 = (TH1D*)hPi0ProjMix->Clone(Form("pt%03i_hPi0MixScaledPol1", ptBin)) ;
	TH1D * hPi0BSPol1 = (TH1D*)hPi0Proj->Clone(Form("pt%03i_hPi0BSPol1", ptBin)) ;

	// Scale Mix by linear part of ratio, yielding approx background
	fbg1->SetParameters(funcRatioFit1->GetParameter(3), funcRatioFit1->GetParameter(4));
	hPi0MixScaledPol1 ->Multiply(fbg1) ;
	hPi0BSPol1->Add(hPi0MixScaledPol1 ,-1.) ;

	Int_t binPi0 = hPi0BSPol1->FindBin(funcRatioFit1->GetParameter(1));
	Int_t nWidPi0 = 2 * (Int_t) (funcRatioFit1->GetParameter(2)/hPi0BSPol1->GetBinWidth(1));
	Int_t integral = TMath::Abs(hPi0BSPol1->Integral(binPi0-nWidPi0,binPi0+nWidPi0));
	fgs->SetParameters(integral/5., funcRatioFit1->GetParameter(1), funcRatioFit1->GetParameter(2)) ;
	fgs->SetParLimits(0,0.,pi0Entries) ;
	fgs->SetParLimits(1,lowerMass,upperMass) ;
	fgs->SetParLimits(2,lowerWidth,upperWidth) ;

	// Fit
	TFitResultPtr bs1FitResultPtr = hPi0BSPol1->Fit(fgs,"MSQ","",rangeMin,rangeMax) ;
	if( int(bs1FitResultPtr) != 4000 ) // if error, "More" error is acceptable
	  bs1FitResultPtr = hPi0BSPol1->Fit(fgs,"MSQ","",rangeMin,rangeMax) ;

	Int_t bs1FitError = bs1FitResultPtr;
	bs1FitError = bs1FitError % 4000; // "More" error is acceptable
	bs1FitError = bs1FitError || TMath::Abs( fgs->GetParameter(1) - lowerMass) < 0.0001 || TMath::Abs( fgs->GetParameter(1) - upperMass) < 0.0001;  // "center" mass converged to limit
	//bs1FitError = bs1FitError || fgs->GetParError(1) > (upperMass - lowerMass)/2; // to large of an error
	bs1FitError = bs1FitError || TMath::Abs( fgs->GetParameter(2) - lowerWidth) < 0.0001 || TMath::Abs( fgs->GetParameter(2) - upperWidth) < 0.0001;  // st. error converged to limit
	//bs1FitError = bs1FitError || fgs->GetParError(2) > (upperWidth - lowerWidth)/2; // to large of an error

	  Double_t y=fgs->GetParameter(0)/hPi0BSPol1->GetXaxis()->GetBinWidth(1) ;
	  Double_t ey=fgs->GetParError(0)/hPi0BSPol1->GetXaxis()->GetBinWidth(1) ;
	  Double_t npiInt = hPi0BSPol1->Integral(intBinMin,intBinMax) ;
	  Double_t norm   = fbg1->GetParameter(0) ;
	  Double_t normErr= fbg1->GetParError(0) ;
	if( bs1FitError) {
	  Printf("in ERROR, %i", bs1FitError);
	  ((TH1D*) output.GetHistogram("mr1_error", outBin))->SetBinContent(ptBin,fgs->GetParameter(1)) ;
	  ((TH1D*) output.GetHistogram("mr1_error", outBin))->SetBinError  (ptBin,fgs->GetParError(1) ) ;
	  ((TH1D*) output.GetHistogram("sr1_error", outBin))->SetBinContent(ptBin,TMath::Abs(fgs->GetParameter(2))) ;
	  ((TH1D*) output.GetHistogram("sr1_error", outBin))->SetBinError  (ptBin,fgs->GetParError(2) ) ;

	  ((TH1D*) output.GetHistogram("yr1_error", outBin))->SetBinContent(ptBin,y/dpt) ;
	  ((TH1D*) output.GetHistogram("yr1_error", outBin))->SetBinError(ptBin,ey/dpt) ;
	  if(npiInt>0.){
	    ((TH1D*) output.GetHistogram("yr1int_error", outBin))->SetBinContent(ptBin,npiInt/dpt) ;
	    ((TH1D*) output.GetHistogram("yr1int_error", outBin))->SetBinError(ptBin,TMath::Sqrt(npiInt + norm*mixInt + normErr*normErr*mixInt*mixInt + norm*norm*mixInt)/dpt) ;
	  }
	}
	if( !bs1FitError || ignoreErrors ) {
	  Printf("converged, status:%i, covMatrixStatus: %i", bs1FitResultPtr->Status(), bs1FitResultPtr->CovMatrixStatus());
	  ((TH1D*) output.GetHistogram("mr1", outBin))->SetBinContent(ptBin,fgs->GetParameter(1)) ;
	  ((TH1D*) output.GetHistogram("mr1", outBin))->SetBinError  (ptBin,fgs->GetParError(1) ) ;
	  ((TH1D*) output.GetHistogram("sr1", outBin))->SetBinContent(ptBin,TMath::Abs(fgs->GetParameter(2))) ;
	  ((TH1D*) output.GetHistogram("sr1", outBin))->SetBinError  (ptBin,fgs->GetParError(2) ) ;

	  ((TH1D*) output.GetHistogram("yr1", outBin))->SetBinContent(ptBin,y/dpt) ;
	  ((TH1D*) output.GetHistogram("yr1", outBin))->SetBinError(ptBin,ey/dpt) ;
	  if(npiInt>0.){
	    ((TH1D*) output.GetHistogram("yr1int", outBin))->SetBinContent(ptBin,npiInt/dpt) ;
	    ((TH1D*) output.GetHistogram("yr1int", outBin))->SetBinError(ptBin,TMath::Sqrt(npiInt + norm*mixInt + normErr*normErr*mixInt*mixInt + norm*norm*mixInt)/dpt) ;
	    // maybe we should use TH1::IntegralAndError
	  }
	}
	// Ploting
	hPi0BSPol1->SetAxisRange(rangeMin, rangeMax);
	hPi0BSPol1->SetMaximum(hPi0BSPol1->GetMaximum()*1.5) ;
	hPi0BSPol1->SetMinimum(hPi0BSPol1->GetMinimum()*0.9) ;
	hPi0BSPol1->SetMarkerStyle(20) ;
	hPi0BSPol1->SetMarkerSize(0.7) ;
	hPi0BSPol1->Draw();
	canvas->Update();
	//printf("end pt");
      }


      // ================================================
      // Pol2 Scaled Background Subtraction
      // ================================================
      canvas->cd(4);
      fbg2->SetParameters(funcRatioFit2->GetParameter(3),funcRatioFit2->GetParameter(4),funcRatioFit2->GetParameter(5));
      if( ! ratioFitError2 || true) {
	printf("Pol1 Scaled Background Subtraction, ");
	TH1D * hPi0MixScaledPol2 = (TH1D*)hPi0ProjMix->Clone(Form("pt%03i_hPi0MixScaledPol2", ptBin)) ;
	TH1D * hPi0BSPol2     = (TH1D*)hPi0Proj    ->Clone(Form("pt%03i_hPi0BSPol2", ptBin)) ;

	hPi0MixScaledPol2->Multiply(fbg2) ;
	hPi0BSPol2 ->Add(hPi0MixScaledPol2,-1.) ;
	hPi0BSPol2->SetOption();

	Int_t binPi0 = hPi0BSPol2->FindBin(funcRatioFit2->GetParameter(1));
	Int_t nWidPi0 = 2 * (Int_t) (funcRatioFit2->GetParameter(2)/hPi0BSPol2->GetBinWidth(1));
	Int_t integral = TMath::Abs(hPi0BSPol2->Integral(binPi0-nWidPi0,binPi0+nWidPi0));
	fgs->SetParameters(integral/5., funcRatioFit2->GetParameter(1), funcRatioFit2->GetParameter(2)) ;
	fgs->SetParLimits(0,0.,pi0Entries) ;
	fgs->SetParLimits(1,lowerMass,upperMass) ;
	fgs->SetParLimits(2,lowerWidth,upperWidth) ;

      	// Fit
	TFitResultPtr bs2FitResultPtr = hPi0BSPol2->Fit(fgs,"MSQ","",rangeMin,rangeMax) ;
	if( int(bs2FitResultPtr) != 4000 ) // if error, "More" error is acceptable
	  bs2FitResultPtr = hPi0BSPol2->Fit(fgs,"MSQ","",rangeMin,rangeMax) ;

	Int_t bs2FitError = bs2FitResultPtr;
	bs2FitError = bs2FitError % 4000; // "More" error is acceptable
	bs2FitError = bs2FitError || TMath::Abs( fgs->GetParameter(1) - lowerMass) < 0.0001 || TMath::Abs( fgs->GetParameter(1) - upperMass) < 0.0001;  // "center" mass converged to limit
// 	bs2FitError = bs2FitError || fgs->GetParError(1) > (upperMass - lowerMass)/2; // to large of an error
	bs2FitError = bs2FitError || TMath::Abs( fgs->GetParameter(2) - lowerWidth) < 0.0001 || TMath::Abs( fgs->GetParameter(2) - upperWidth) < 0.0001;  // st. error converged to limit
// 	bs2FitError = bs2FitError || fgs->GetParError(2) > (upperWidth - lowerWidth)/2; // to large of an error

	  Double_t y=fgs->GetParameter(0)/hPi0BSPol2->GetXaxis()->GetBinWidth(1) ;
	  Double_t ey=fgs->GetParError(0)/hPi0BSPol2->GetXaxis()->GetBinWidth(1) ;
	  Double_t npiInt = hPi0BSPol2->Integral(intBinMin,intBinMax) ;
	  Double_t norm   = fbg2->GetParameter(0) ;
	  Double_t normErr= fbg2->GetParError(0) ;
	if( bs2FitError ) {
	  Printf("in ERROR, %i", bs2FitError);
	  ((TH1D*) output.GetHistogram("mr2_error", outBin))->SetBinContent(ptBin,fgs->GetParameter(1)) ;
	  ((TH1D*) output.GetHistogram("mr2_error", outBin))->SetBinError  (ptBin,fgs->GetParError(1) ) ;
	  ((TH1D*) output.GetHistogram("sr2_error", outBin))->SetBinContent(ptBin,TMath::Abs(fgs->GetParameter(2))) ;
	  ((TH1D*) output.GetHistogram("sr2_error", outBin))->SetBinError  (ptBin,fgs->GetParError(2) ) ;

	  ((TH1D*) output.GetHistogram("yr2_error", outBin))->SetBinContent(ptBin,y/dpt) ;
	  ((TH1D*) output.GetHistogram("yr2_error", outBin))->SetBinError(ptBin,ey/dpt) ;
	  if(npiInt>0.){
	    ((TH1D*) output.GetHistogram("yr2int_error", outBin))->SetBinContent(ptBin,npiInt/dpt) ;
	    ((TH1D*) output.GetHistogram("yr2int_error", outBin))->SetBinError(ptBin,TMath::Sqrt(npiInt + norm*mixInt + normErr*normErr*mixInt*mixInt + norm*norm*mixInt)/dpt) ;
	    // maybe we should use TH1::IntegralAndError
	    }
	}
	if( !bs2FitError || ignoreErrors ) {
	  Printf("converged, status:%i, covMatrixStatus: %i", bs2FitResultPtr->Status(), bs2FitResultPtr->CovMatrixStatus());
	  ((TH1D*) output.GetHistogram("mr2", outBin))->SetBinContent(ptBin,fgs->GetParameter(1)) ;
	  ((TH1D*) output.GetHistogram("mr2", outBin))->SetBinError  (ptBin,fgs->GetParError(1) ) ;
	  ((TH1D*) output.GetHistogram("sr2", outBin))->SetBinContent(ptBin,TMath::Abs(fgs->GetParameter(2))) ;
	  ((TH1D*) output.GetHistogram("sr2", outBin))->SetBinError  (ptBin,fgs->GetParError(2) ) ;

	  ((TH1D*) output.GetHistogram("yr2", outBin))->SetBinContent(ptBin,y/dpt) ;
	  ((TH1D*) output.GetHistogram("yr2", outBin))->SetBinError(ptBin,ey/dpt) ;
	  if(npiInt>0.){
	    ((TH1D*) output.GetHistogram("yr2int", outBin))->SetBinContent(ptBin,npiInt/dpt) ;
	    ((TH1D*) output.GetHistogram("yr2int", outBin))->SetBinError(ptBin,TMath::Sqrt(npiInt + norm*mixInt + normErr*normErr*mixInt*mixInt + norm*norm*mixInt)/dpt) ;
	    // maybe we should use TH1::IntegralAndError
	    }
	}
	// Plotting
	hPi0BSPol2->SetAxisRange(rangeMin, rangeMax);
	hPi0BSPol2->SetMaximum(hPi0BSPol2->GetMaximum()*1.5) ;
	hPi0BSPol2->SetMinimum(hPi0BSPol2->GetMinimum()*0.9) ;
	hPi0BSPol2->SetMarkerStyle(20) ;
	hPi0BSPol2->SetMarkerSize(0.7) ;
	hPi0BSPol2->Draw();
	canvas->Update();

      }

      canvas->cd(1);
      TH1D * hPi0MixScaled = (TH1D*)hPi0ProjMix ->Clone(Form("%sScaled", hPi0Proj->GetName())) ;
      //hPi0MixScaled->Scale(double(pi0Entries)/mixEntries);
      hPi0MixScaled->Scale(fbg2->Eval(0.134));
      hPi0MixScaled->SetLineColor(kRed);
      hPi0MixScaled->SetTitle(Form("%s, Scaled", hPi0Proj->GetName()));
      hPi0Proj->SetAxisRange(rangeMin, 1.);
      //hPi0Proj->SetMaximum(TMath::Max(hPi0Proj->GetMaximum(), hPi0ProjMix->GetMaximum())*1.2);
      hPi0Proj->SetMinimum(0);
      hPi0Proj->Draw();
      hPi0MixScaled->Draw("same");

      canvas->Update();

      canvas->Print(Form("imgs/%s_ptBin:%03i.pdf", outBin.Key().Data(), ptBin));
      canvas->Print(Form("imgs/%s_ptBin:%03i.png", outBin.Key().Data(), ptBin));

      std::cout << std::endl;
    }// end of Pt slice loop


    //Normalize by the number of events
    Int_t cMin=0, cMax=0;
    if( input.Bin().Trigger().EqualTo("kCentral") )
      switch(outBin.Centrality()) {
      case 0: cMin = 1; cMax = 5; break;
      case 1: cMin = 6; cMax = 10; break;
      case -1: cMin = 1; cMax = 10; break;
      default: Printf("ERROR: cent bin not defined for trigger");
      }
    else if( input.Bin().Trigger().EqualTo("kSemiCentral") )
      switch(outBin.Centrality()) {
      case 0: cMin = 11; cMax = 20; break;
      case 1: cMin = 21; cMax = 30; break;
      case 2: cMin = 31; cMax = 40; break;
      case 3: cMin = 41; cMax = 50; break;
      case -2: cMin = 11; cMax = 20; break;
      case -3: cMin = 21; cMax = 30; break;
      case -4: cMin = 31; cMax = 40; break;
      case -5: cMin = 41; cMax = 50; break;
      default: Printf("ERROR: cent bin not defined for trigger");
      }
    else if ( input.Bin().Trigger().EqualTo("kMB") || input.Bin().Trigger().EqualTo("kPHOSPb") )
      switch(outBin.Centrality()) {
      case 0: cMin = 1; cMax = 5; break;
      case 1: cMin = 6; cMax = 10; break;
      case 2: cMin = 11; cMax = 20; break;
      case 3: cMin = 21; cMax = 30; break;
      case 4: cMin = 31; cMax = 40; break;
      case 5: cMin = 41; cMax = 50; break;
      case 6: cMin = 51; cMax = 80; break;
      case -10: cMin=1; cMax = 80; break;
      case -1: cMin = 1; cMax = 10; break;
      case -2: cMin = 11; cMax = 20; break;
      case -3: cMin = 21; cMax = 30; break;
      case -4: cMin = 31; cMax = 40; break;
      case -5: cMin = 41; cMax = 50; break;
      case -6: cMin = 51; cMax = 80; break;
      default: Printf("ERROR: cent bin not defined for trigger");
      }
    else
      Printf("ERROR: cent bins not defined for trigger, %s", input.Bin().Trigger().Data());

    Double_t nevents = hCentralityX->Integral(cMin,cMax);
    if ( nevents > 0.9 ) {
      ((TH1D*) output.GetHistogram("yr1", outBin)) ->Scale(1./nevents) ;
      ((TH1D*) output.GetHistogram("yr1int", outBin)) ->Scale(1./nevents) ;
      ((TH1D*) output.GetHistogram("yr2", outBin)) ->Scale(1./nevents) ;
      ((TH1D*) output.GetHistogram("yr2int", outBin)) ->Scale(1./nevents) ;

      ((TH1D*) output.GetHistogram("yr1_error", outBin)) ->Scale(1./nevents) ;
      ((TH1D*) output.GetHistogram("yr1int_error", outBin)) ->Scale(1./nevents) ;
      ((TH1D*) output.GetHistogram("yr2_error", outBin)) ->Scale(1./nevents) ;
      ((TH1D*) output.GetHistogram("yr2int_error", outBin)) ->Scale(1./nevents) ;
    } else {
      Printf("WARNING: non positive nEvents in centrality range, cMin:%d, cMax:%d, nEvents:%f", cMin, cMax, nevents );

    }


    // store to file
    //  TFile* outputFile = TFile::Open(saveToFileName.Data(), "UPDATE");
    //outputList->Write(input.KeySuggestion(), TObject::kSingleKey);
    //outputList->Write();
//     outputFile->Write();
//     outputFile->Close();
//     delete outputFile;
    // delete list and content from memory
    //outputList->SetOwner(kTRUE);
    //outputList->Delete("slow");
    //delete outputList;
    //output.Write();
    delete canvas;
  }


  //-----------------------------------------------------------------------------
  double GetPtBin(int bin){
    // recusive function used by MakePtBins

    if( bin==0 )
      return 1.;

    // return GetPtBin(bin-1) * 1.1;

    // if ( bin % 2 )
    //   return GetPtBin(bin-1) + 0.4;
    // else
    //   return GetPtBin(bin-1) + 0.2;

    double previousBin = GetPtBin(bin-1);
    double linInc = 1.;
    double threshold = 5.;
    double logFact = 1 + linInc/threshold;
    if ( previousBin < threshold ) // linear
      return double(int((previousBin + linInc)*10))/10;
    else { // logarithmic
      return double(int((previousBin * logFact)*10))/10;
    }
  }

  //-----------------------------------------------------------------------------
  void MakePtBins() {
    // function for setting Pt bins

    int bin = -1;
    do {
      ++bin;
      ptBinEdges[bin] = GetPtBin(bin);
    } while(ptBinEdges[bin] < 40.);
    nPtBins = bin -2;

    printf("Making Pt Bins:\n");
    for(int b=0; b < nPtBins+1; ++b)
      printf("%.1f, ", ptBinEdges[b]);
    printf("\n N. Bins: %d\n", nPtBins);


    // for(int bin = 0; bin <= nPtBins; ++bin){
    //   ptBinEdges[bin] = GetPtBin(bin);
    //   printf("%.1f, ", ptBinEdges[bin]);
    // }
    // printf("\n");
  }

  //-----------------------------------------------------------------------------
  Double_t CGausPol1(Double_t * x, Double_t * par){
    //Parameterization of Real/Mixed ratio
    Double_t m=par[1] ;
    Double_t s=par[2] ;
    Double_t dx=(x[0]-m)/s ;
    return par[0]*exp(-dx*dx/2.)+par[3]+par[4]*(x[0]-kMean) ;
  }

  //-----------------------------------------------------------------------------
  Double_t CGausPol2(Double_t * x, Double_t * par){
    //Another parameterization of Real/Mixed ratio7TeV/README
    Double_t m=par[1] ;
    Double_t s=par[2] ;
    Double_t dx=(x[0]-m)/s ;
    return par[0]*exp(-dx*dx/2.)+par[3]+par[4]*(x[0]-kMean)+par[5]*(x[0]-kMean)*(x[0]-kMean) ;
  }
  //-----------------------------------------------------------------------------
  Double_t CGausPol0(Double_t * x, Double_t * par){
    //Parameterizatin of signal
    Double_t m=par[1] ;
    Double_t s=par[2] ;
    Double_t dx=(x[0]-m)/s ;
    return par[0]*exp(-dx*dx/2.)/TMath::Sqrt(TMath::TwoPi())/s+par[3] ;
  }
  //-----------------------------------------------------------------------------
  Double_t CPol1(Double_t * x, Double_t * par){
    //Normalizatino of Mixed
    return par[0]+par[1]*(x[0]-kMean) ;
  }
  //-----------------------------------------------------------------------------
  Double_t CPol2(Double_t * x, Double_t * par){
    //Another normalization of  Mixed
    return par[0]+par[1]*(x[0]-kMean)+par[2]*(x[0]-kMean)*(x[0]-kMean) ;
  }


  // Input Definitions
  TFile* Input::fFile = 0;
  Input::Input(const TString& fileName, const RawProduction::TriggerBin& inputBin, TString listPath)
  : fList(0x0), fBin(inputBin.Trigger())
  {
    // File
    if(fFile && !fileName.EqualTo(fFile->GetName())){
      fFile->Close();
      delete fFile;
      fFile = 0x0;
    } else if(! fFile) {
      Printf("Opening %s", fileName.Data());
     fFile = TFile::Open(fileName.Data(), "READ");
    }

    if( listPath.EqualTo("") ) {
      char cstr[256] = "";
      sprintf(cstr, "PHOSPi0Flow_%s/PHOSPi0Flow_%sCoutput1", fBin.Trigger().Data(), fBin.Trigger().Data());
      listPath = cstr;
    }

    Printf("Getting list, %s", listPath.Data());
    fFile->GetObject(listPath.Data(), fList);
    if( !fList )
      Printf("ERROR: list not found");
  }

  TH1* Input::GetHistogram(const char* name){
    TObject* obj = fList->FindObject(name);
    TH1* hist = dynamic_cast<TH1*> (obj);
    if( ! hist)
      Printf("MakePi0FitInput::GetHistogram: Error, could not find object of name: %s or cast to hist", name);
    return hist;
  }

  //OutputBin Definitions
  TriggerBin::TriggerBin(const TString& trigger)
  : fTrigger(trigger), fKey(trigger)
  { }

  TriggerBin::TriggerBin(const char* trigger)
  : fTrigger(trigger), fKey(trigger)
  { }

  TriCenPidBin::TriCenPidBin(Int_t centrality, const TString& pid, const TString& trigger)
  : TriggerBin(trigger), fCentrality(centrality), fPID(pid)
  {
    fKey.Form("c%03i_%s_%s", centrality, pid.Data(), trigger.Data());
  }

  Output::Output(const TString& fileName, const char* options)
  : fFile(0x0)
  {
    fFile = TFile::Open(fileName.Data(), options);
  }

  void Output::SetDir(const TriggerBin& inBin)
  {
    Bool_t success = fFile->cd(inBin.Key().Data());
    if( ! success ) {
      TDirectory* newDir = fFile->mkdir(inBin.Key().Data());
      newDir->cd();
    }
  }

  TH1* Output::GetHistogram(const TString& name, const RawProduction::TriggerBin& inBin)
  {
    TDirectory* dir = fFile->GetDirectory(inBin.Key().Data(), true);
    TH1* hist = dynamic_cast<TH1*>( dir->Get(name.Data()) );
    if( hist )
      return hist;
    else {
      Printf("ERROR: Output::GetHistogram: %s could not be found", name.Data());
      return 0x0;
    }
  }



  void Output::Write()
  {
    fFile->Write();
  }




  TH1* GetHistogram_cent(Input& input, const TString& name, int centrality)
  {
    // Getter (and merger) for histograms following the naming patern of %s_cen%i
    //
    // For certain negeative values, function is defined to merge centrality bins.
    // -1   0-10%
    // -2  10-20%
    // -3  20-30%
    // -4  30-40%
    // -5  40-50%
    // -6  50-80%
    // -10  0-80%
    if( centrality >= 0 ) {
      char cname[256] = "";
      sprintf(cname, "%s_cen%i", name.Data(), centrality);
      input.GetHistogram(cname);
    }

    TH1* hist = 0x0;
    if( input.Bin().Trigger().EqualTo("kMB") || input.Bin().Trigger().EqualTo("kPHOSPb") ) {
      switch(centrality) {
      case -10: hist = MergeHistogram_cent(input, name, centrality, 0, 7); break;
      case -1:  hist = MergeHistogram_cent(input, name, centrality, 0, 2); break;
      case -2:  hist = MergeHistogram_cent(input, name, centrality, 2, 3); break;
      case -3:  hist = MergeHistogram_cent(input, name, centrality, 3, 4); break;
      case -4:  hist = MergeHistogram_cent(input, name, centrality, 4, 5); break;
      case -5:  hist = MergeHistogram_cent(input, name, centrality, 5, 6); break;
      case -6:  hist = MergeHistogram_cent(input, name, centrality, 6, 7); break;
      }
    } else if ( input.Bin().Trigger().EqualTo("kCentral") ) {
      switch( centrality ) {
      case -1: return MergeHistogram_cent(input, name, centrality, 0, 2); break;
      }
    } else if ( input.Bin().Trigger().EqualTo("kSemiCentral") ) {
      switch( centrality ) {
      case -2: return MergeHistogram_cent(input, name, centrality, 0, 1); break;
      case -3: return MergeHistogram_cent(input, name, centrality, 1, 2); break;
      case -4: return MergeHistogram_cent(input, name, centrality, 2, 3); break;
      case -5: return MergeHistogram_cent(input, name, centrality, 3, 4); break;
      }
    }
    // in case not defined above
    if( ! hist ) {
      Printf("ERROR:GetHistogram_cent: %i not possible for %s trigger", centrality, input.Bin().Trigger().Data());
      return 0x0;
    }

    switch(centrality) {
      case -10: hist->SetTitle( Form("%s, 0-80%% centrality", hist->GetTitle())); break;
      case -1:  hist->SetTitle( Form("%s, 0-10%% centrality", hist->GetTitle())); break;
      case -2:  hist->SetTitle(Form("%s, 10-20%% centrality", hist->GetTitle())); break;
      case -3:  hist->SetTitle(Form("%s, 20-30%% centrality", hist->GetTitle())); break;
      case -4:  hist->SetTitle(Form("%s, 30-40%% centrality", hist->GetTitle())); break;
      case -5:  hist->SetTitle(Form("%s, 40-50%% centrality", hist->GetTitle())); break;
      case -6:  hist->SetTitle(Form("%s, 50-80%% centrality", hist->GetTitle())); break;
    }
    return hist;
  }

  TH1* MergeHistogram_cent(Input& input, const TString& name, int newCentIndex, int fromCentIndex, int toCentIndex)
  {
    // Merger (All cent) for histograms following the naming patern of %s_cen%i, from including to excluding
    //
    // Merges centralites bins into one histogram, from including to excluding, and names the histogram using the patern above.
    // If an histogram with that name Allready exist in the current directory (gDirectory), then no merge
    // occurs and this hist. is simply returned.

    char mergeHistName[256] = "";
    sprintf(mergeHistName, "%s_cen%i", name.Data(), newCentIndex);

    // Check if histogram allready exists in current directory.
    TH1* mergeHist = dynamic_cast<TH1*>( gDirectory->FindObject(mergeHistName) );
    if ( mergeHist )
      return mergeHist;

    // If not so; get the first hist, clone, and add the others
    char cname[256] = "";
    sprintf(cname, "%s_cen%i", name.Data(), fromCentIndex);
    TH1* hist0 = input.GetHistogram(cname);
    sprintf(cname, "%s_cen%i", name.Data(), newCentIndex);
    TH1 * histMerged = (TH1*) hist0->Clone(cname);

    for(int cent=fromCentIndex+1; cent < toCentIndex; ++cent) {
      sprintf(cname, "%s_cen%i", name.Data(), cent);
      TH1* nextHist = input.GetHistogram(cname);
      if( ! nextHist ) {Printf("ERROR: Merge of histograms failed"); delete histMerged; return 0x0; }
      histMerged->Add(nextHist);
    }
    return histMerged;
  }



  //-----------------------------------------------------------------------------
  void PPRstyle()
  {

    //////////////////////////////////////////////////////////////////////
    //
    // ROOT style macro for the TRD TDR
    //
    //////////////////////////////////////////////////////////////////////

    gStyle->SetPalette(1);
    gStyle->SetCanvasBorderMode(-1);
    gStyle->SetCanvasBorderSize(1);
    gStyle->SetCanvasColor(10);

    gStyle->SetFrameFillColor(10);
    gStyle->SetFrameBorderSize(1);
    gStyle->SetFrameBorderMode(-1);
    gStyle->SetFrameLineWidth(1.2);
    gStyle->SetFrameLineColor(1);

    gStyle->SetHistFillColor(0);
    gStyle->SetHistLineWidth(1);
    gStyle->SetHistLineColor(1);

    gStyle->SetPadColor(10);
    gStyle->SetPadBorderSize(1);
    gStyle->SetPadBorderMode(-1);

    gStyle->SetStatColor(10);
    gStyle->SetTitleColor(kBlack,"X");
    gStyle->SetTitleColor(kBlack,"Y");

    gStyle->SetLabelSize(0.04,"X");
    gStyle->SetLabelSize(0.04,"Y");
    gStyle->SetLabelSize(0.04,"Z");
    gStyle->SetTitleSize(0.04,"X");
    gStyle->SetTitleSize(0.04,"Y");
    gStyle->SetTitleSize(0.04,"Z");
    gStyle->SetTitleFont(42,"X");
    gStyle->SetTitleFont(42,"Y");
    gStyle->SetTitleFont(42,"X");
    gStyle->SetLabelFont(42,"X");
    gStyle->SetLabelFont(42,"Y");
    gStyle->SetLabelFont(42,"Z");
    gStyle->SetStatFont(42);

    gStyle->SetTitleOffset(1.0,"X");
    gStyle->SetTitleOffset(1.4,"Y");

    gStyle->SetFillColor(kWhite);
    gStyle->SetTitleFillColor(kWhite);

    gStyle->SetOptDate(0);
    gStyle->SetOptTitle(1);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(111);

  }
}


void MakeRawProduction()
{
  RawProduction::Output output;


  //TStringToken triggers("kMB kCentral kSemiCentral kPHOSPb", " ");
  TStringToken triggers("kMB kPHOSPb", " ");
  while(triggers.NextToken()) {
    RawProduction::TriggerBin inBin(triggers);
    RawProduction::Input input("AnalysisResults.root", inBin);
    TStringToken pids("All Allcore Allwou Disp Disp2 Dispcore Dispwou CPV CPVcore CPV2 Both Bothcore", " ");
    //TStringToken pids("Bothcore", " ");
    while(pids.NextToken()) {
      RawProduction::TriCenPidBin outBin(-10, pids, inBin.Trigger());
      RawProduction::MakePi0Fit(input, outBin, output);
    }
  }
  output.Write();

}

void MakeRawProductionAll()
{
  //gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  RawProduction::Output output;

  TStringToken triggers("kMB kCentral kSemiCentral kPHOSPb", " ");
  //TStringToken triggers("kCentral ", " ");
  while(triggers.NextToken()) {
    RawProduction::TriggerBin triggerBin(triggers);
    RawProduction::Input input("AnalysisResults.root", triggerBin);

    RawProduction::MakePi0Fit(input, triggerBin, output);

    //TStringToken pids("All Allcore Allwou Disp Disp2 Dispcore Dispwou CPV CPVcore CPV2 Both Bothcore", " ");
    //TStringToken pids("All Allcore Allwou Disp Disp2 Dispcore Dispwou", " ");
    TStringToken pids("All", " ");
    while(pids.NextToken()) {
      for(int cent = -1; cent > -7; cent--) {
	if(triggers.EqualTo("kCentral") && cent != -1) continue;
	if(triggers.EqualTo("kSemiCentral") && !(-1 > cent && cent > -6 )) continue;

	RawProduction::TriCenPidBin tcpBin(cent, pids, triggerBin.Trigger());
	RawProduction::MakePi0FitTCP(input, tcpBin, output);
      }

      if( triggers.EqualTo("kCentral") || triggers.EqualTo("kSemiCentral") ) continue;
      RawProduction::TriCenPidBin tcpBinc(-10, pids, triggerBin.Trigger());
      RawProduction::MakePi0FitTCP(input, tcpBinc, output);
    }
  }
  output.Write();
}
