///
/// \file InvMassFit.C
/// \ingroup CaloTrackCorrMacrosPlotting
/// \brief Fit invariant mass distributions
///
/// Macro using as input the 2D histograms mass vs pT of AliAnaPi0
/// For a given set of pT bins invariant mass plots are fitted and mass vs pT 
/// and width vs pT and neutral meson spectra plots are obtained.
/// Also pure MC histograms checking the origin of the clusters and generated spectra
/// is used to plot efficiencies and acceptance
///
/// Based on old macros
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)
///

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TString.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TH3F.h>
#include <TH1D.h>
#include <TF1.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TPad.h>
#include <TFile.h>
#include <TLegend.h>
#include <TObject.h>
#include <TDirectoryFile.h>
#include <TGraphErrors.h>
#include <TList.h>
#include <TArrayD.h>
#include <TGaxis.h>
#include <TROOT.h>
#include "PlotUtils.C"

#endif

// Settings
Bool_t  kMix        = kFALSE;              /// use mixed event to constrain combinatorial background
Int_t   kPolN       = 1;                   /// polinomyal type for residual background under the peak
Bool_t  kSumw2      = kTRUE;               /// Apply Root method Sumw2(), off for pT hard prod where already done before
Float_t kNPairCut   = 20;                  /// Minimum number of cluster pairs in pi0 or eta window 
Float_t kFirstTRDSM = -1;                  /// First Calo SM covered by a TRD SM, 6 in 2011 and 4 in 2012-13
TString kHistoStartName = "AnaPi0";        /// Common starting name in histograms
TString kHistoList  = "default";           /// TList name, usually each list is a different trigger
TString kProdName   = "LHC18c3_NystromOn"; /// Input production directory name where input file is located
TString kFileName   = "AnalysisResults";   /// Name of file with input histograms
TString kPlotFormat = "eps";               /// Automatic plots format
TString kCalorimeter= "EMCAL";             /// Calorimeter, EMCAL, DCAL (PHOS)
TString kParticle   = "Pi0";               /// Particle searched: "Pi0", "Eta"
Bool_t  kTrunMixFunc= kTRUE;               /// Use a truncated function to get the mixed event scale
Float_t kChi2NDFMax = 1000;                /// Maximum value of chi2/ndf to define a fit as good (set lower than 1e6)
Int_t   kRebin      = 1;                   /// Rebin invariant mass histograms default binning with this value

// Initializations
Double_t nEvt = 0;
TFile * fil   = 0;
TFile * fout  = 0;
TList * lis   = 0;
TF1   * fitfun= 0;
TDirectoryFile * direc =0;
// Mixing
TF1 *tPolPi0 = 0;
TF1 *tPolEta = 0;

Int_t modColorIndex[]={1 , 1, 2, 2, 3, 3, 4, 4, 7, 7, 6, 6, 2, 3, 4, 7, 6, 2, 2, 3, 3, 4, 4, 6, 6};
Int_t modStyleIndex[]={24,25,25,24,25,24,25,24,25,24,25,21,21,21,21,21,22,26,22,26,22,26,22,26};

///
/// Open the file and the list and the number of analyzed events
/// 
//-----------------------------------------------------------------------------
Bool_t GetFileAndEvents( TString dirName , TString listName /*, TString calo*/ )
{
  fil = new TFile(Form("%s/%s.root",kProdName.Data(),kFileName.Data()),"read");
  
  printf("Open Input File: %s/%s.root\n",kProdName.Data(),kFileName.Data());
  
  if ( !fil ) return kFALSE;
  
  direc  = (TDirectoryFile*) fil->Get(dirName);
  
  if ( !direc && dirName != "" ) return kFALSE;
  
  //printf("dir %p, list %s\n",dir,listName.Data());
  
  if(direc)
    lis = (TList*) direc ->Get(listName);
  else
    lis = (TList*) fil->Get(listName);
//  lis = (TList*) fil->Get(Form("CTC_%s_Trig_%s_SPDPileUp",calo.Data(),listName.Data()));
  
  if ( !lis && listName != "") return kFALSE;
  
  if(!lis)
    nEvt = ((TH1F*) fil->Get("hNEvents"))->GetBinContent(1);
  else
    nEvt = ((TH1F*) lis->FindObject("hNEvents"))->GetBinContent(1);
  
  printf("nEvt = %2.3e\n",nEvt);
  //nEvt = 1;
  
  return kTRUE;
}

///
/// Initialize the fitting function
/// The function definitions are in PlotUtils.C
///
//-----------------------------------------------------------------------------
void SetFitFun()
{
  //Fitting function
  //if(mix) kPolN = 0;
  
  if(kParticle == "Pi0")
  {
    if ( kPolN == 0)
      fitfun = new TF1("fitfun",FunctionGaussPol0,0.100,0.250,4);
    else if      (kPolN == 1)
      fitfun = new TF1("fitfun",FunctionGaussPol1,0.100,0.250,5);
    else if (kPolN == 2)
      fitfun = new TF1("fitfun",FunctionGaussPol2,0.100,0.250,6);
    else if (kPolN == 3)
      fitfun = new TF1("fitfun",FunctionGaussPol3,0.100,0.250,7);
    else
    {
      printf("*** <<< Set Crystal Ball!!! >>> ***\n");
      fitfun = new TF1("fitfun",FunctionCrystalBall,0.100,0.250,6);
    }
    
    if(kPolN < 4)
    {
      fitfun->SetParLimits(0,  kNPairCut/5,kNPairCut*1.e4);
      fitfun->SetParLimits(1,  0.105,0.185);
      fitfun->SetParLimits(2,  0.001,0.040);
    }
    else
    {
      fitfun->SetParLimits(0,  kNPairCut/5,kNPairCut*1.e6);
      fitfun->SetParLimits(2,  0.105,0.185);
      fitfun->SetParLimits(1,  0.001,0.040);
      fitfun->SetParLimits(3,  0.001,0.040);
      fitfun->SetParLimits(4,  0,10);
      fitfun->SetParLimits(5,  1,1.e+6);
    }
    
  } // Pi0
  else  // Eta
  {
    if ( kPolN == 0)
      fitfun = new TF1("fitfun",FunctionGaussPol0,0.400,0.650,4);
    else if      (kPolN == 1)
      fitfun = new TF1("fitfun",FunctionGaussPol1,0.400,0.650,5);
    else if (kPolN == 2)
      fitfun = new TF1("fitfun",FunctionGaussPol2,0.400,0.650,6);
    else if (kPolN == 3)
      fitfun = new TF1("fitfun",FunctionGaussPol3,0.400,0.650,7);
    if(kPolN < 4)
    {
      fitfun->SetParLimits(0,  kNPairCut/10,1.e+6);
      fitfun->SetParLimits(1,  0.20,0.80);
      fitfun->SetParLimits(2,  0.001,0.06);
    }
  }
  
  fitfun->SetLineColor(kRed);
  fitfun->SetLineWidth(2);
  
  fitfun->SetParName(0,"A");
  fitfun->SetParName(1,"m_{0}");
  fitfun->SetParName(2,"#sigma");
  //  fitfun->SetParName(3,"a_{0}");
  //  fitfun->SetParName(4,"a_{1}");
  
  if (kPolN > 1)
  {
    fitfun->SetParName(5,"a_{2}");
    if (kPolN > 2) fitfun->SetParName(6,"a_{3}"); 
  }  
  
  //Mix fit func
  tPolPi0    = new TF1("truncatedPolPi0"   , FunctionTruncatedPolPi0   , 0.02, 0.25, 2);
  tPolEta    = new TF1("truncatedPolEta"   , FunctionTruncatedPolEta   , 0.2, 1, 2);
}

//-----------------------
/// Efficiency calculation called in ProjectAndFit()
//-----------------------
void Efficiency
(Int_t nPt, 
 TArrayD xPt      , TArrayD exPt      ,
 TArrayD mesonPt  , TArrayD emesonPt  , 
 TArrayD mesonPtBC, TArrayD emesonPtBC,
 TString hname)
{
  gStyle->SetOptTitle(1);
  gStyle->SetTitleOffset(2,"Y");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(000000);
  
  TGraphErrors * gPrim    = (TGraphErrors*) fout->Get("Primary");
  TGraphErrors * gPrimAcc = (TGraphErrors*) fout->Get("PrimaryInAcceptance");
  if ( !gPrim || !gPrimAcc ) return ;
    
  TArrayD  effPt        ;  effPt        .Set(nPt); 
  TArrayD  effBCPt      ;  effBCPt      .Set(nPt); 
  TArrayD  effPtAcc     ;  effPtAcc     .Set(nPt); 
  TArrayD  effBCPtAcc   ;  effBCPtAcc   .Set(nPt);  
  TArrayD  effPtErr     ;  effPtErr     .Set(nPt); 
  TArrayD  effBCPtErr   ;  effBCPtErr   .Set(nPt); 
  TArrayD  effPtAccErr  ;  effPtAccErr  .Set(nPt); 
  TArrayD  effBCPtAccErr;  effBCPtAccErr.Set(nPt);  
  
  for(Int_t ibin = 0; ibin < nPt; ibin++ )
  {
    //printf("Bin %d \n",ibin);
    //printf("- Fit Reco %2.3e, Prim %2.3e, PrimAcc %2.3e\n",
    //       mesonPt[ibin],gPrim->GetY()[ibin],gPrimAcc->GetY()[ibin]);
    
    if ( gPrim->GetY()[ibin] > 0 && mesonPt[ibin] > 0 )
    {
      effPt   [ibin] = 100 * mesonPt[ibin] / gPrim->GetY()[ibin];
      effPtErr[ibin] = 100 * GetFractionError( mesonPt[ibin], gPrim->GetY ()[ibin], 
                                              emesonPt[ibin], gPrim->GetEY()[ibin]); // PlotUtils.C
      
      //printf("\t EffxAcc %f, err %f\n",effPt[ibin],effPtErr[ibin]);
    }
    else { effPt[ibin] = 0; effPtErr[ibin] = 0; }
    
    if ( gPrimAcc->GetY()[ibin] > 0 && mesonPt[ibin] > 0 )
    {
      effPtAcc   [ibin] = 100 * mesonPt[ibin] / gPrimAcc->GetY()[ibin];
      effPtAccErr[ibin] = 100 * GetFractionError( mesonPt[ibin], gPrimAcc->GetY ()[ibin], 
                                                 emesonPt[ibin], gPrimAcc->GetEY()[ibin]); // PlotUtils.C
      //printf("\t Eff %f, err %f\n",effPtAcc[ibin],effPtAccErr[ibin]);
    }
    else { effPtAcc[ibin] = 0; effPtAccErr[ibin] = 0; }
    
    if ( kMix || hname.Contains("MCpT") )
    {
      //printf("- BC  Reco %2.3e, Prim %2.3e, PrimAcc %2.3e\n",
      //       mesonPtBC[ibin],gPrim->GetY()[ibin],gPrimAcc->GetY()[ibin]);
      
      if ( gPrim->GetY()[ibin] > 0 && mesonPtBC[ibin] > 0 )
      {
        effBCPt   [ibin] = 100 * mesonPtBC[ibin] / gPrim->GetY()[ibin];
        effBCPtErr[ibin] = 100 * GetFractionError( mesonPtBC[ibin], gPrim->GetY ()[ibin], 
                                                  emesonPtBC[ibin], gPrim->GetEY()[ibin]); // PlotUtils.C
       //printf("\t EffxAcc %f, err %f\n",effBCPt[ibin],effBCPtErr[ibin]);
      }
      else { effBCPt[ibin] = 0; effBCPtErr[ibin] = 0; }
      
      if ( gPrimAcc->GetY()[ibin] > 0 && mesonPtBC[ibin] > 0 )
      {
        effBCPtAcc   [ibin] = 100 * mesonPtBC[ibin] / gPrimAcc->GetY()[ibin];
        effBCPtAccErr[ibin] = 100 * GetFractionError( mesonPtBC[ibin], gPrimAcc->GetY ()[ibin], 
                                                     emesonPtBC[ibin], gPrimAcc->GetEY()[ibin]); // PlotUtils.C
        //printf("\t EffxAcc %f, err %f\n",effBCPtAcc[ibin],effBCPtAccErr[ibin]);
      }
      else { effBCPtAcc[ibin] = 0; effBCPtAccErr[ibin] = 0; }
      
    } // Bin counting
    
  } // pt bin loop
  
  TGraphErrors *gEff = new TGraphErrors(nPt,xPt.GetArray(),effPt.GetArray(),exPt.GetArray(),effPtErr.GetArray());
  gEff->SetName(Form("EfficiencyxAcceptance_%s",hname.Data()));
  gEff->GetHistogram()->SetYTitle("#epsilon_{Reco #times PID #times Acc} (%)");
  gEff->GetHistogram()->SetTitleOffset(1.4,"Y"); 
  gEff->GetHistogram()->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  gEff->GetHistogram()->SetTitleOffset(1.2,"X");
  gEff->GetHistogram()->SetTitle("Reconstruction efficiency #times acceptance");
  gEff->SetMarkerColor(1);
  gEff->SetLineColor(1);
  gEff->SetMarkerStyle(20);
  gEff->Write();
  
  TGraphErrors *gEffAcc = new TGraphErrors(nPt,xPt.GetArray(),effPtAcc.GetArray(),exPt.GetArray(),effPtAccErr.GetArray());
  gEffAcc->SetName(Form("Efficiency_%s",hname.Data()));
  gEffAcc->GetHistogram()->SetYTitle("#epsilon_{Reco #times PID} (%)");
  gEffAcc->GetHistogram()->SetTitleOffset(1.5,"Y"); 
  gEffAcc->GetHistogram()->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  gEffAcc->GetHistogram()->SetTitleOffset(1.2,"X");
  gEffAcc->SetMarkerColor(1);
  gEffAcc->SetLineColor(1);
  gEffAcc->SetMarkerStyle(20);
  gEffAcc->GetHistogram()->SetTitle("Reconstruction efficiency");

  gEffAcc->Write();    
  
  TGraphErrors *gBCEff    = 0;
  TGraphErrors *gBCEffAcc = 0;
  if(kMix)
  {
    gBCEff = new TGraphErrors(nPt,xPt.GetArray(),effBCPt.GetArray(),exPt.GetArray(),effBCPtErr.GetArray());
    gBCEff->SetName(Form("EfficiencyxAcceptance_BC_%s",hname.Data()));
    gBCEff->GetHistogram()->SetYTitle("#epsilon_{Reco #times PID #times Acc}");
    gBCEff->GetHistogram()->SetTitleOffset(1.5,"Y"); 
    gBCEff->GetHistogram()->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    gBCEff->GetHistogram()->SetTitleOffset(1.2,"X");
    gBCEff->SetMarkerColor(4);
    gBCEff->SetLineColor(4);
    gBCEff->SetMarkerStyle(24);
    gBCEff->Write();
    
    gBCEffAcc = new TGraphErrors(nPt,xPt.GetArray(),effBCPtAcc.GetArray(),exPt.GetArray(),effBCPtAccErr.GetArray());
    gBCEffAcc->SetName(Form("Efficiency_BC_%s",hname.Data()));
    gBCEffAcc->GetHistogram()->SetYTitle("#epsilon_{Reco #times PID}");
    gBCEffAcc->GetHistogram()->SetTitleOffset(1.4,"Y"); 
    gBCEffAcc->GetHistogram()->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    gBCEffAcc->GetHistogram()->SetTitleOffset(1.2,"X");
    gBCEffAcc->SetMarkerColor(4);
    gBCEffAcc->SetLineColor(4);
    gBCEffAcc->SetMarkerStyle(24);
    gBCEffAcc->Write();   
  }
  
  // Plot efficiencies
  //
  TCanvas *cEff = new TCanvas(Form("cEff_%s",hname.Data()),
                              Form("Efficiency Graphs for %s",hname.Data()),2*600,600);
  cEff->Divide(2, 1);
  
  cEff->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  //gPad->SetLogy();
  
  //gEff->SetMaximum(1);
  //gEff->SetMinimum(1e-8);
  
  gEff->Draw("AP");
  
  if(kMix)
  {
    gBCEff ->Draw("P");
    
    TLegend * legend =  new TLegend(0.5,0.7,0.9,0.9);
    legend->AddEntry(gEff  ,"From fit","P");
    legend->AddEntry(gBCEff,"From bin counting","P");
    legend->Draw();
  }
  
  cEff->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  //gPad->SetLogy();
  
  //gEffAcc->SetMaximum(1);
  //gEffAcc->SetMinimum(1e-8);
  
  gEffAcc->Draw("AP");
  
  if(kMix)
  {
    gBCEffAcc->Draw("P");
    
    TLegend * legend =  new TLegend(0.5,0.7,0.9,0.9);
    legend->AddEntry(gEffAcc  ,"From fit","P");
    legend->AddEntry(gBCEffAcc,"From bin counting","P");
    legend->Draw();
  }
  
  cEff->Print(Form("%s/%s/%s/%s/Efficiency_%s.%s",
                   kProdName.Data(),kCalorimeter.Data(),kHistoList.Data(),kParticle.Data(),hname.Data(),
                   /*kFileName.Data(),*/ kPlotFormat.Data()));   
}

//-------------------------------------
/// Do energy projections, fits and plotting
/// of invariant mass distributions.
//-------------------------------------
void ProjectAndFit
(
 const Int_t nPt, TString comment, TString leg, TString hname,   
 TArrayD xPtLimits, TArrayD xPt, TArrayD exPt, 
 TH2F * hRe, TH2F * hMi
 )
{
  if ( !hRe ) 
  {
    printf("ProjectAndFit() - Null inv mass 2D histo for %s, skip!\n",hname.Data());
    return;
  }
  
  gStyle->SetOptTitle(1);
  gStyle->SetTitleOffset(2,"Y");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(000000);
  
  Double_t mmin = 0;
  Double_t mmax = 1;
  TString  particleN = "";

  if(kParticle=="Pi0")
  {
    mmin = 0.04;
    mmax = 0.44;
    particleN = " #pi^{0}";
  }
  else // eta
  {
    mmin = 0.25;
    mmax = 0.95;
    particleN = " #eta";
  }
  
  TLegend * pLegendIM[nPt];
  TH1D* hIM           [nPt];
  TH1D* hMix          [nPt];
  TH1D* hMixCorrected [nPt];
  TH1D* hRatio        [nPt];
  TH1D* hSignal       [nPt];
  for(Int_t ipt = 0; ipt< nPt; ipt++)
  {
    hIM          [ipt] = 0 ;
    hMix         [ipt] = 0 ;
    hMixCorrected[ipt] = 0 ;
    hRatio       [ipt] = 0 ;
    hSignal      [ipt] = 0 ;
  }

  TArrayD  mesonPt   ;  mesonPt   .Set(nPt); 
  TArrayD  mesonPtBC ;  mesonPtBC .Set(nPt); 
  TArrayD  mesonMass ;  mesonMass .Set(nPt); 
  TArrayD  mesonWidth;  mesonWidth.Set(nPt);  
  TArrayD emesonPt   ; emesonPt   .Set(nPt); 
  TArrayD emesonPtBC ; emesonPtBC .Set(nPt); 
  TArrayD emesonMass ; emesonMass .Set(nPt); 
  TArrayD emesonWidth; emesonWidth.Set(nPt); 
  
  Int_t rebin2 = 2*kRebin;
  
  Int_t col = TMath::Ceil(TMath::Sqrt(nPt));

  //hRe->Scale(1./nEvt);
  
  // Mix correction factor
  TF1 *line3[nPt];// = new TF1("line3", "[0]+[1]*x+[2]*x*x+[3]*x*x*x", 0.0, 1);
  
  TF1 * fitFunction = 0;
  Bool_t ok = kFALSE;
  
  TCanvas  * cIMModi = new TCanvas(Form("c%s", hname.Data()), Form("%s", leg.Data()), 1200, 1200) ;
  cIMModi->Divide(col, col);
  
  for(Int_t i = 0; i < nPt; i++)
  {    
    cIMModi->cd(i+1) ; 
    //gPad->SetLogy();
    
    Double_t ptMin = xPtLimits[i];
    Double_t ptMax = xPtLimits[i+1];
    //printf("Bin %d (%2.1f, %2.1f)\n",i, ptMin,ptMax);
    
    hIM[i] = hRe->ProjectionY(Form("IM_%s_PtBin%d",hname.Data(),i),
                                            hRe->GetXaxis()->FindBin(ptMin),
                                            hRe->GetXaxis()->FindBin(ptMax));
    hIM[i]->SetTitle(Form("%2.1f < #it{p}_{T, #gamma#gamma} < %2.1f GeV/#it{c}",ptMin,ptMax));
    if(i < nPt-3)
      hIM[i]->Rebin(kRebin);
    else
      hIM[i]->Rebin(rebin2);
    
    if ( kSumw2 )
      hIM[i]->Sumw2();
    
    hIM[i]->SetXTitle("M_{#gamma,#gamma} (GeV/#it{c}^{2})");
    //printf("mmin %f, mmax %f\n",mmin,mmax);
    hIM[i]->SetAxisRange(mmin,mmax,"X");
    hIM[i]->SetLineWidth(2);
    hIM[i]->SetLineColor(4);
    
    Double_t mStep = hIM[i]->GetBinWidth(1);
    
    hSignal[i] = (TH1D*) hIM[i]->Clone();
        
    //--------------------------------------------------
    //   Mix: Project, scale and subract to real pairs
    //--------------------------------------------------
    if ( kMix && hMi )
    {
      hMix[i] = (TH1D*) hMi->ProjectionY(Form("MiMass_PtBin%d_%s",i,hname.Data()),
                                                       hMi->GetXaxis()->FindBin(ptMin), 
                                                       hMi->GetXaxis()->FindBin(ptMax));
      hMix[i]->SetLineColor(1);  
      hMix[i]->SetTitle(Form("%2.2f < #it{p}_{T, #gamma#gamma} < %2.2f GeV/#it{c}",ptMin,ptMax));
      if(i < nPt-3)
        hMix[i]->Rebin(kRebin);
      else
        hMix[i]->Rebin(rebin2);
      
      //         Double_t sptmin = 0.2;
      //         Double_t sptmax = 0.45;
      
      //         Double_t scalereal = hIM->Integral(hIM->FindBin(sptmin),
      //                                                             hIM->FindBin(sptmax));
      //         Double_t scalemix  = hMix->Integral(hMix->FindBin(sptmin),
      //                                                             hMix->FindBin(sptmax));
      //         if(scalemix > 0)
      //           hMix->Scale(scalereal/scalemix);
      //         else{ 
      //           //printf("could not scale:  re %f  mi %f\n",scalereal,scalemix);
      //           continue;
      //         }
      
      //      hMix[i]->SetLineWidth(1);  
      
      if ( kSumw2 ) 
        hMix[i]->Sumw2();
      
      // --- Ratio ---  
      hRatio[i] = (TH1D*)hIM[i]->Clone(Form("RatioRealMix_PtBin%d_%s",i,hname.Data()));
      hRatio[i]->SetAxisRange(mmin,mmax,"X");
      hRatio[i]->Divide(hMix[i]);
      Double_t p0 =0;
      Double_t p1 =0;
      Double_t p2 =0;
      Double_t p3 =0;
      
      // --- Subtract from signal ---
      hSignal[i] ->SetName(Form("Signal_PtBin%d_%s",i,hname.Data()));

      if ( hMix[i]->GetEntries() > kNPairCut*3 )
      {
        hSignal[i] ->SetLineColor(kViolet);

        if ( kTrunMixFunc )
        {
          if ( kParticle == "Pi0" )
          {
            hRatio[i]->Fit("truncatedPolPi0", "NQRL","",0.11,0.4);
            p0=  tPolPi0->GetParameter(0);           
            p1=  tPolPi0->GetParameter(1);
            //p2=  tPolPi0->GetParameter(2);
            //p3=  tPolPi0->GetParameter(3);
          } // pi0
          else  // eta
          {
            hRatio[i]->SetAxisRange(mmin,mmax);
            hRatio[i]->Fit("truncatedPolEta", "NQRL","",0.35,0.95);
            p0 =  tPolEta->GetParameter(0);      
            p1 =  tPolEta->GetParameter(1);
            //p2 =  tPolEta->GetParameter(2);
            //p3 =  tPolEta->GetParameter(3);
          }
          
          //line3[i] = new TF1("line3", "[0]+[1]*x+[2]*x*x+[3]*x*x*x", mmin, mmax);
          line3[i] = new TF1("line3", "[0]+[1]*x", mmin, mmax);
          line3[i]->SetParameters(p0,p1,p2,p3);
          
          //printf("Correct\n");
          hMixCorrected[i] = (TH1D*) hMix[i]->Clone(Form("MixCorrected_PtBin%d_%s",i,hname.Data()));
          for(Int_t j = 0; j< hMix[i]->GetNbinsX(); j++)
          {
            Double_t x          = hMix[i]->GetBinCenter(j);
            Double_t correction = line3[i]->Eval(x);
            Double_t corrected  = hMix[i]->GetBinContent(j)*correction;
            Double_t ecorrected = hMix[i]->GetBinError(j)  *correction;
            //printf("bin %d, center %f, mix %f, corrected %f\n",i,x,hMix[i]->GetBinContent(j),corrected);
            hMixCorrected[i]->SetBinContent(j, corrected);
            hMixCorrected[i]->SetBinError  (j,ecorrected);
          }

          // Subtract
          hSignal[i] ->Add(hMixCorrected[i],-1);
        }
        else
        {
          Float_t scale = hRatio[i]->GetBinContent(hRatio[i]->FindBin(0.7));
          hMix[i]->Scale(scale);
          // Subtract
          hSignal[i] ->Add(hMix[i],-1);
        }
      } // enough mixed entries
    } // mixed event histogtam exists
    
    //----------------------------
    // ---- Fit subtracted signal
    //----------------------------
    Double_t nMax = 0;
    if(kParticle=="Pi0")
    {
      nMax= hSignal[i]->Integral(hIM[i]->FindBin(0.1),
                                 hIM[i]->FindBin(0.2));
    }
    else
    {
      nMax = hSignal[i]->Integral(hIM[i]->FindBin(0.5),
                                  hIM[i]->FindBin(0.7));
    }
    
    if(nMax > kNPairCut)
    {
      fitfun->SetParLimits(0,nMax/100,nMax*100);
      
      if(kParticle=="Pi0")
      {
        if(kPolN < 4 )fitfun->SetParameters(nMax/5,0.135,20,0);
        else         fitfun->SetParameters(nMax/5,20,0.135,20,20,nMax/5);
        
        if(i < nPt-4)
          hSignal[i]->Fit("fitfun","QR","",0.11,0.3);
        else
          hSignal[i]->Fit("fitfun","QR","",0.11,0.3);
      }// pi0
      else // eta
      {
        
        if ( kPolN < 4 ) fitfun->SetParameters(nMax/5,0.54,40,0);
        else             fitfun->SetParameters(nMax/5,40,0.54,40,40,nMax/5);
        //if(i<4)
        hSignal[i]->Fit("fitfun","QR","",0.4,0.7);
        //else if(i < 15)
        //  hSignal[i]->Fit("fitfun","QRL","",0.3,0.8);
        //else 
        //hSignal[i]->Fit("fitfun","QRL","",0.3,0.8);
      }
    }
    else
      printf("Skip bin %d: n max %2.3e < cut %2.3e\n",i, nMax, kNPairCut);
    
    
    //----------------------------
    // ---- Put fit results in arrays
    //----------------------------
    // Init results arrays
     mesonPt   .SetAt(-1,i);
    emesonPt   .SetAt(-1,i);    
     mesonPtBC .SetAt(-1,i);
    emesonPtBC .SetAt(-1,i);
     mesonMass .SetAt(-1,i);
    emesonMass .SetAt(-1,i);
     mesonWidth.SetAt(-1,i);
    emesonWidth.SetAt(-1,i);
    
    fitFunction = (TF1*) hSignal[i]->GetFunction("fitfun");
    
    Float_t chi2ndf = 1e6;
    if ( fitFunction )
    {
      //            printf("ipt %d: Chi/NDF %f, Chi %f, NDF %d\n",i,
      //                   fitFunction->GetChisquare()/fitFunction->GetNDF(),
      //                   fitFunction->GetChisquare(),
      //                   fitFunction->GetNDF());
      
      Float_t chi2 = fitFunction->GetChisquare();
      Int_t   ndf  = fitFunction->GetNDF();
      if( ndf > 0 ) chi2ndf = chi2 / ndf;
      
      if ( chi2ndf < kChi2NDFMax )
      {
        Double_t A    = fitFunction->GetParameter(0);
        Double_t mean = fitFunction->GetParameter(1);
        Double_t sigm = fitFunction->GetParameter(2);
        //              Double_t a0   = fitFunction->GetParameter(3);
        //              Double_t a1   = fitFunction->GetParameter(4);
        //              Double_t a2   = 0;
        //              if (kPolN == 2)
        //                a2 = fitFunction->GetParameter(5);
        
        Double_t eA    = fitFunction->GetParError(0);
        Double_t emean = fitFunction->GetParError(1);
        Double_t esigm = fitFunction->GetParError(2);
        //              Double_t ea0   = fitFunction->GetParError(3);
        //              Double_t ea1   = fitFunction->GetParError(4);
        //              Double_t ea2   = 0;
        //              if (kPolN == 2)
        //                ea2 = fitFunction->GetParError(5);
        pLegendIM[i] = new TLegend(0.48,0.65,0.95,0.93);
        pLegendIM[i]->SetTextSize(0.035);
        pLegendIM[i]->SetFillColor(10);
        pLegendIM[i]->SetBorderSize(1);
        pLegendIM[i]->SetHeader(Form(" %s - %s",particleN.Data(), leg.Data()));
        pLegendIM[i]->AddEntry("",Form("A = %2.1e#pm%2.1e ",A,eA),"");
        pLegendIM[i]->AddEntry("",Form("#mu = %3.1f #pm %3.1f GeV/#it{c}^{2}",mean*1000,emean*1000),"");
        pLegendIM[i]->AddEntry("",Form("#sigma = %3.1f #pm %3.1f GeV/#it{c}^{2}",sigm*1000,esigm*1000),"");
        
        //pLegendIM[i]->AddEntry("",Form("p_{0} = %2.1f #pm %2.1f ",a0,ea0),"");
        //pLegendIM[i]->AddEntry("",Form("p_{1} = %2.1f #pm %2.1f ",a1,ea1),"");
        //if(kPolN==2)pLegendIM[i]->AddEntry("",Form("p_{2} = %2.1f#pm %2.1f ",a2,ea2),"");
        
        Double_t  counts = A*sigm / mStep * TMath::Sqrt(TMath::TwoPi());
        Double_t eCounts =
        TMath::Power(eA/A,2) +
        TMath::Power(esigm/sigm,2);
        eCounts = TMath::Sqrt(eCounts) * counts;
        //eCounts = TMath::Min(eCounts, TMath::Sqrt(counts));
        
         counts/=(nEvt*(exPt.At(i)*2));
        eCounts/=(nEvt*(exPt.At(i)*2));
        
         mesonPt.SetAt( counts,i);
        emesonPt.SetAt(eCounts,i);
        
        //printf("A %e Aerr %e; mu %f muErr %f; sig %f sigErr %f",
        //       A,eA,mean,emean,sigm,esigm);
        //cout << "N(pi0) fit  = " <<  counts << "+-"<<  eCounts << " mstep " << mStep<<endl;
        
         mesonMass .SetAt( mean*1000.,i);
        emesonMass .SetAt(emean*1000.,i);
        
         mesonWidth.SetAt( sigm*1000.,i);
        emesonWidth.SetAt(esigm*1000.,i);
      } //Good fit
      else 
      {
        printf("Bin %d, Bad fit, Chi2 %f ndf %d, ratio %f!\n",
             i,chi2,ndf,chi2ndf);
      }
      
    }
    else printf("Bin %d, NO fit available!\n",i);
    
    // Set the integration window, depending on mass mean and width
    // 2 sigma
    Double_t mass  = mesonMass [i]/1000.;
    Double_t width = mesonWidth[i]/1000.;
    Double_t mMinBin =   mass - 2 * width;
    Double_t mMaxBin =   mass + 2 * width;
    if(mass <= 0 || width <= 0)
    {
      if(kParticle=="Pi0")
      {
        mMinBin = 0.115; mMaxBin = 0.3 ;
      }
      else // eta
      {
        mMinBin = 0.4;  mMaxBin = 0.9 ; 
      }
    }
    
    //
    // Bin counting instead of fit 
    //
    //printf("Signal %p, Mix %p, hname %s\n",hSignal[i], hMix[i], hname.Data());
    if ( hSignal[i] && ( hMix[i] || hname.Contains("MCpT") ) )
    {
      Double_t  countsBin = 0.;//hSignal[i]->Integral(hSignal[i]->FindBin(mMinBin),  hSignal[i]->FindBin(mMaxBin)) ;
      Double_t eCountsBin = 0.;//TMath::Sqrt(countsBin);
     
      // PlotUtils.C
      GetRangeIntegralAndError(hSignal[i], hSignal[i]->FindBin(mMinBin), hSignal[i]->FindBin(mMaxBin), countsBin, eCountsBin );

       countsBin/=(nEvt*(exPt.At(i)*2));
      eCountsBin/=(nEvt*(exPt.At(i)*2));
      
       mesonPtBC.SetAt( countsBin,i);
      emesonPtBC.SetAt(eCountsBin,i);
      
      //printf("pt bin %d: [%2.1f,%2.1f] - Mass window [%2.2f, %2.2f] - N(pi0) BC = %2.3e +- %2.4e\n",
      //       i, xPtLimits[i], xPtLimits[i+1], mMinBin, mMaxBin, countsBin, eCountsBin);      
    }
        
    if ( nMax > kNPairCut && chi2ndf < kChi2NDFMax )
    {
      if( kMix && hMix[i] )
      {
        //    if ( !kMix ) hIM[i]->SetMinimum( 0.1);
        //    else         hIM[i]->SetMinimum(-0.1);
        
        hIM[i]->SetMaximum(hIM[i]->GetMaximum()*1.2);
        hIM[i]->SetMinimum(hSignal[i]->GetMinimum()*0.2);
        
        hIM[i]->Draw("HE");
        pLegendIM[i]->AddEntry(hIM[i],"Raw pairs"   ,"L");

        if ( hMix[i]->GetEntries() > kNPairCut*3 )
        {
          if ( kTrunMixFunc )
          {
            hMixCorrected[i]->Draw("same");
            pLegendIM[i]->AddEntry(hMixCorrected[i],"Mixed pairs" ,"L");
          }
          else            
          {
            hMix[i]->Draw("same");
            pLegendIM[i]->AddEntry(hMix[i],"Mixed pairs" ,"L");
          }
        }

        ok = kTRUE;

        hSignal[i] -> Draw("same");
        pLegendIM[i]->AddEntry(hSignal[i],"Signal pairs","L");
      }
      else 
      {
        ok = kTRUE;
        hSignal[i]->SetMaximum(hSignal[i]->GetMaximum()*1.2);
        hSignal[i]->SetMinimum(hSignal[i]->GetMinimum()*0.8);
        
        pLegendIM[i]->AddEntry(hSignal[i],"Raw pairs","L");
        hSignal[i] -> Draw("HE");
      }

      // Plot mass from pairs originated from pi0/eta
      if ( hname.Contains("All") )
      {
        TH1F* hMesonPairOrigin = (TH1F*) fout->Get(Form("IM_MCpTReco_PtBin%d",i));
        if ( hMesonPairOrigin )
        {
          hMesonPairOrigin->SetLineColor(8);
          hMesonPairOrigin->Draw("same");
          pLegendIM[i]->AddEntry(hMesonPairOrigin,Form("%s origin",particleN.Data()),"L");
        }
      } // MC origin
      
//      if ( fitFunction ) 
//        pLegendIM[i]->AddEntry(fitFunction,"Fit","L");

      pLegendIM[i]->Draw();
    }
    else 
    {
      // Plot just raw pairs
      hIM[i]->Draw("HE");
    }
  } //pT bins
  
  if ( ok ) 
    cIMModi->Print(Form("%s/%s/%s/%s/Mgg_%s.%s",
                        kProdName.Data(),kCalorimeter.Data(),kHistoList.Data(),kParticle.Data(),/*kFileName.Data(),*/
                        hname.Data(),kPlotFormat.Data()));
  
  //xxxx Real / Mixxxxx
  if(kMix)
  {
    //printf("Do real/mix\n");
    Bool_t okR = kFALSE;
    TCanvas  * cRat = new TCanvas(Form("Ratio_%s\n",hname.Data()), Form("Ratio %s\n", leg.Data()), 1200, 1200) ;
    cRat->Divide(col, col);
    
    for(Int_t i = 0; i < nPt; i++)
    {
      cRat->cd(i+1) ;
      //gPad->SetGridy();
      //gPad->SetLog();
      if(!hRatio[i]) 
      {
        //printf("No ratio in pt bin %d continue\n",i);
        continue;
      }
      
      Double_t nMax =0;
      if(kParticle=="Pi0")
        nMax = hRatio[i]->Integral(hRatio[i]->FindBin(0.005),
                                   hRatio[i]->FindBin(0.20));
      else
        nMax = hRatio[i]->Integral(hRatio[i]->FindBin(0.4),
                                   hRatio[i]->FindBin(0.6));
      
      if ( nMax <= 0 ) continue; 

      //printf("Ratio nMax = %e \n",nMax);
            
      hRatio[i]->SetMaximum(nMax/4);
      hRatio[i]->SetMinimum(1e-6);
      hRatio[i]->SetAxisRange(mmin,mmax,"X");
      
      okR = kTRUE;
      
      hRatio[i]->SetYTitle("Real pairs / Mixed pairs");
      hRatio[i]->SetTitleOffset(1.5,"Y");
      
      hRatio[i]->Draw();
      
      //if(line3[i]) line3[i]->Draw("same");
    }
    
    if ( okR )
      cRat->Print(Form("%s/%s/%s/%s/MggRatio_%s.%s",
                       kProdName.Data(),kCalorimeter.Data(),kHistoList.Data(),kParticle.Data(),/*kFileName.Data(),*/ 
                       hname.Data(),kPlotFormat.Data()));
  }   // Mix  
  
  if( !ok ) return ;
  
  //------------------------------
  // Fit parameters 
  //------------------------------
  // Put fit results in TGraphErrors
  
  TGraphErrors* gPt   =  new TGraphErrors(nPt,xPt.GetArray(),mesonPt.GetArray(),exPt.GetArray(),emesonPt.GetArray());
  gPt->SetName(Form("gPt_%s",hname.Data()));
  gPt->GetHistogram()->SetTitle(Form("#it{p}_{T} of reconstructed %s, %s ",particleN.Data(), comment.Data()));
  gPt->GetHistogram()->SetXTitle("#it{p}_{T} (GeV/#it{c})    ");
  gPt->GetHistogram()->SetYTitle(Form("dN_{%s}/d#it{p}_{T} (GeV/#it{c})^{-1}  / N_{events}  ",particleN.Data()));
  //gPt->GetHistogram()->SetAxisRange(0.,30);
  gPt->GetHistogram()->SetTitleOffset(1.5,"Y");
  gPt->SetMarkerStyle(20);
  //gPt->SetMarkerSize(1);
  gPt->SetMarkerColor(1);
  //   gPt->GetHistogram()->SetMaximum(1e8);
  //   gPt->GetHistogram()->SetMinimum(1e-8);
  
  TGraphErrors* gPtBC =  new TGraphErrors(nPt,xPt.GetArray(),mesonPtBC.GetArray(),exPt.GetArray(),emesonPtBC.GetArray());
  gPtBC->SetName(Form("gPtBC_%s",hname.Data()));
  gPtBC->GetHistogram()->SetTitle(Form("#it{p}_{T} of reconstructed %s, %s ",particleN.Data(), comment.Data()));
  gPtBC->GetHistogram()->SetXTitle("#it{p}_{T} (GeV/#it{c})    ");
  gPtBC->GetHistogram()->SetYTitle(Form("d#it{N}_{%s}/d#it{p}_{T} (GeV/#it{c})^{-1}  / #it{N}_{events}  ",particleN.Data()));
  //gPtBC->GetHistogram()->SetAxisRange(0.,30);
  gPtBC->GetHistogram()->SetTitleOffset(1.5,"Y");
  gPtBC->SetMarkerStyle(24);
  //gPtBC->SetMarkerSize(1);
  gPtBC->SetMarkerColor(4);
  //   gPtBC->GetHistogram()->SetMaximum(1e8);
  //   gPtBC->GetHistogram()->SetMinimum(1e-8);
  
  TGraphErrors* gMass =  new TGraphErrors(nPt,xPt.GetArray(),mesonMass.GetArray(),exPt.GetArray(),emesonMass.GetArray());
  gMass->SetName(Form("gMass_%s",hname.Data()));
  gMass->GetHistogram()->SetTitle(Form("mass vs #it{p}_{T} of reconstructed %s, %s ",particleN.Data(), comment.Data()));
  gMass->GetHistogram()->SetXTitle("#it{p}_{T} (GeV/#it{c})    ");
  gMass->GetHistogram()->SetYTitle(Form("mass_{%s} (MeV/#it{c}^{2})    ",particleN.Data()));
  //gMass->GetHistogram()->SetAxisRange(0.,30);
  gMass->GetHistogram()->SetTitleOffset(1.5,"Y");
  gMass->SetMarkerStyle(20);
  //gMass->SetMarkerSize(1);
  gMass->SetMarkerColor(1);
  
  TGraphErrors* gWidth=  new TGraphErrors(nPt,xPt.GetArray(),mesonWidth.GetArray(),exPt.GetArray(),emesonWidth.GetArray());
  gWidth->SetName(Form("gWidth_%s",hname.Data()));
  gWidth->GetHistogram()->SetTitle(Form("Width vs #it{p}_{T} of reconstructed %s, %s ",particleN.Data(), comment.Data()));
  gWidth->GetHistogram()->SetXTitle("#it{p}_{T} (GeV/#it{c})    ");
  gWidth->GetHistogram()->SetYTitle(Form("#sigma_{%s} (MeV/#it{c}^{2})    ",particleN.Data()));
  //gWidth->GetHistogram()->SetAxisRange(0.,30);
  gWidth->GetHistogram()->SetTitleOffset(1.5,"Y");
  gWidth->SetMarkerStyle(20);
  //gWidth->SetMarkerSize(1);
  gWidth->SetMarkerColor(1);
  
  // Plot the fitted results
  
  TCanvas *cFitGraph = new TCanvas(Form("cFitGraph_%s",hname.Data()),
                                   Form("Fit Graphs for %s",hname.Data()),600,600);
  cFitGraph->Divide(2, 2);

  // Mass  
  cFitGraph->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  
  if ( kParticle=="Pi0" )
  {
    gMass->SetMaximum(160);
    gMass->SetMinimum(100);
  }
  else
  {
    gMass->SetMaximum(660);
    gMass->SetMinimum(420);
  }
  
  gMass->Draw("AP");
  
  cFitGraph->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  
  if(kParticle=="Pi0")
  {
    gWidth->SetMaximum(16);
    gWidth->SetMinimum(6);
  }
  else
  {
    gWidth->SetMaximum(60);
    gWidth->SetMinimum(0);
  }
  
  gWidth->Draw("AP");
 
  cFitGraph->cd(3);
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogy();
  
  gPt->Draw("AP");
  
  Double_t maximum = gPt->GetHistogram()->GetMaximum();
  Double_t minimum = gPt->GetHistogram()->GetMinimum();
  //printf("A-maximum %e, minimum %e\n",maximum,minimum);
  
//  Double_t maximum = gPrim->GetHistogram()->GetMaximum();
//  Double_t minimum = gPrim->GetHistogram()->GetMinimum();
//  
//  if ( gPrimAcc->GetMaximum() > maximum ) maximum = gPrimAcc->GetMaximum() ;
//  if ( gPrimAcc->GetMinimum() > minimum ) minimum = gPrimAcc->GetMinimum() ;
//  
//  gPrim->SetMaximum(maximum*10);
//  gPrim->SetMinimum(minimum/10);
  
  if(kMix)
  {
    gPtBC ->Draw("P");
    
    if(maximum < gPtBC->GetHistogram()->GetMaximum()) maximum = gPtBC->GetHistogram()->GetMaximum();
    if(minimum > gPtBC->GetHistogram()->GetMinimum()) minimum = gPtBC->GetHistogram()->GetMaximum();
    //printf("B-maximum %e, minimum %e\n",maximum,minimum);
    
    if(minimum < 0 ) minimum = 2.e-9 ;
    gPtBC->SetMaximum(maximum*2.);
    gPtBC->SetMinimum(minimum/2.);
    
    TLegend * legend =  new TLegend(0.5,0.7,0.9,0.9);
    legend->AddEntry(gPt  ,"From fit","P");
    legend->AddEntry(gPtBC,"From bin counting","P");
    legend->Draw();
  }
  
  gPt->SetMaximum(maximum*2.);
  gPt->SetMinimum(minimum/2.);
  
  cFitGraph->Print(Form("%s/%s/%s/%s/MassWidthPtSpectra_%s.%s",
                        kProdName.Data(),kCalorimeter.Data(),kHistoList.Data(),kParticle.Data(),hname.Data(),
                        /*kFileName.Data(),*/ kPlotFormat.Data())); 
  
  //-----------------------
  // Write results to file
  //-----------------------
  hRe->Write();
  if ( hMi ) hMi->Write();
  
  for(Int_t ipt = 0; ipt < nPt; ipt++)
  {
    if ( hIM[ipt] ) hIM[ipt]->Write(); 

    if ( kMix )
    {
      if (hMix         [ipt] ) hMix         [ipt]->Write();
      if (hMixCorrected[ipt] ) hMixCorrected[ipt]->Write();
      if (hRatio       [ipt] ) hRatio       [ipt]->Write();
      if (hSignal      [ipt] ) hSignal      [ipt]->Write();
    }
  }
  
  gPt   ->Write();
  gMass ->Write();
  gWidth->Write();  
  if ( kMix ) gPtBC->Write();

  // Do some final efficiency calculations if MC   
  Efficiency(nPt, xPt, exPt, mesonPt, emesonPt, mesonPtBC, emesonPtBC, hname);
}

///
/// Plot combinations of graphs with fit results
/// SM/Sector/Side/SM group by SM/Sector/Side/SM group comparisons
///
//-----------------------------------------------------------------------------
void PlotGraphs(TString opt, Int_t first, Int_t last)
{
  gStyle->SetOptTitle(1);
  gStyle->SetTitleOffset(2,"Y");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(000000);
  
  const Int_t nCombi = last-first+1;
  Int_t nCombiActive = 0;
  
  Float_t xmin = 0.7;
  Float_t ymin = 0.3;
  Float_t xmax = 0.9;
  Float_t ymax = 0.9;
  
  if(opt.Contains("Group"))
  {
    xmin = 0.3;
    ymin = 0.7;
  }
 
  if(opt.Contains("Sector"))
  {
    xmin = 0.65;
    ymin = 0.5;
  }
  
  if(opt.Contains("Side"))
  {
    xmin = 0.65;
  }
  
  TLegend * legend =  new TLegend(xmin,ymin,xmax,ymax);
  legend->SetTextSize(0.05);
  
  TString particleN = " #pi^{0}";
  if(kParticle=="Eta") particleN = " #eta";
  
  // recover the graphs from output file
  //printf("%s\n",opt.Data());
  
  
  // Recover SUM
  TString opt2 = opt;
  if(opt.Contains("SMGroup")) opt2 = "SM";
  
  TGraphErrors* gSumPt    = (TGraphErrors*) fout->Get(Form("gPt_Same%s"   ,opt2.Data()));
  TGraphErrors* gSumPtBC  = (TGraphErrors*) fout->Get(Form("gPtBC_Same%s" ,opt2.Data()));
  TGraphErrors* gSumMass  = (TGraphErrors*) fout->Get(Form("gMass_Same%s" ,opt2.Data()));
  TGraphErrors* gSumWidth = (TGraphErrors*) fout->Get(Form("gWidth_Same%s",opt2.Data()));
  //printf("\t Sum %p %p %p %p\n",gSumPt,gSumPtBC,gSumMass,gSumWidth);
  
  if ( !gSumMass ) return ;
  
  gSumPt   ->SetMarkerStyle(20);
  gSumPt   ->SetMarkerColor(1);
  gSumPt   ->SetLineColor  (1);    
  
  gSumMass->SetMarkerStyle(20);
  gSumMass->SetMarkerColor(1);
  gSumMass->SetLineColor  (1);   
  
  gSumWidth->SetMarkerStyle(20);
  gSumWidth->SetMarkerColor(1);
  gSumWidth->SetLineColor  (1);
  
  if ( kMix )
  {
    gSumPtBC ->SetMarkerStyle(20);
    gSumPtBC ->SetMarkerColor(1);
    gSumPtBC ->SetLineColor  (1);   
  }
  
  legend->AddEntry(gSumPt,"Sum","P");

  // Normalize to total number of SMs/Sectors/Sides
  if ( opt != "SMGroup" )
  {
    for (Int_t ibin = 0; ibin < gSumPt->GetN(); ibin++) 
    {
      gSumPt->GetY ()[ibin] /= nCombi;  
      gSumPt->GetEY()[ibin] /= nCombi;  
      if ( kMix )
      {
        gSumPtBC->GetY ()[ibin] /= nCombi; 
        gSumPtBC->GetEY()[ibin] /= nCombi; 
      }
    }
  }
  
  TGraphErrors* gSumTRDnotPt    = 0;
  TGraphErrors* gSumTRDnotPtBC  = 0;
  TGraphErrors* gSumTRDnotMass  = 0;
  TGraphErrors* gSumTRDnotWidth = 0;  
  TGraphErrors* gSumTRDyesPt    = 0;
  TGraphErrors* gSumTRDyesPtBC  = 0;
  TGraphErrors* gSumTRDyesMass  = 0;
  TGraphErrors* gSumTRDyesWidth = 0;
  
  //
  //// TRD covered
  // 
  if ( kFirstTRDSM > 0 && opt == "SM" ) // it could be extended to sector/side
  {
    gSumTRDnotPt    = (TGraphErrors*) fout->Get(Form("gPt_Same%sTRDNot"   ,opt2.Data()));
    gSumTRDnotPtBC  = (TGraphErrors*) fout->Get(Form("gPtBC_Same%sTRDNot" ,opt2.Data()));
    gSumTRDnotMass  = (TGraphErrors*) fout->Get(Form("gMass_Same%sTRDNot" ,opt2.Data()));
    gSumTRDnotWidth = (TGraphErrors*) fout->Get(Form("gWidth_Same%sTRDNot",opt2.Data()));
    printf("\t Sum Not TRD %p %p %p %p\n",gSumTRDnotPt,gSumTRDnotPtBC,gSumTRDnotMass,gSumTRDnotWidth);
        
    gSumTRDnotPt   ->SetMarkerStyle(20);
    gSumTRDnotPt   ->SetMarkerColor(kGray);
    gSumTRDnotPt   ->SetLineColor  (kGray);    
    
    gSumTRDnotMass->SetMarkerStyle(20);
    gSumTRDnotMass->SetMarkerColor(kGray);
    gSumTRDnotMass->SetLineColor  (kGray);   
    
    gSumTRDnotWidth->SetMarkerStyle(20);
    gSumTRDnotWidth->SetMarkerColor(kGray);
    gSumTRDnotWidth->SetLineColor  (kGray);
    
    if ( kMix )
    {
      gSumTRDnotPtBC ->SetMarkerStyle(20);
      gSumTRDnotPtBC ->SetMarkerColor(kGray);
      gSumTRDnotPtBC ->SetLineColor  (kGray);   
    }
    
    legend->AddEntry(gSumTRDnotPt,"w/o TRD","P");
   
    for (Int_t ibin = 0; ibin < gSumPt->GetN(); ibin++) 
    {
      gSumTRDnotPt->GetY ()[ibin] /= kFirstTRDSM;  
      gSumTRDnotPt->GetEY()[ibin] /= kFirstTRDSM;  
      if ( kMix )
      {
        gSumTRDnotPtBC->GetY ()[ibin] /= kFirstTRDSM; 
        gSumTRDnotPtBC->GetEY()[ibin] /= kFirstTRDSM; 
      }
    }
    
    gSumTRDyesPt    = (TGraphErrors*) fout->Get(Form("gPt_Same%sTRDYes"   ,opt2.Data()));
    gSumTRDyesPtBC  = (TGraphErrors*) fout->Get(Form("gPtBC_Same%sTRDYes" ,opt2.Data()));
    gSumTRDyesMass  = (TGraphErrors*) fout->Get(Form("gMass_Same%sTRDYes" ,opt2.Data()));
    gSumTRDyesWidth = (TGraphErrors*) fout->Get(Form("gWidth_Same%sTRDYes",opt2.Data()));
    printf("\t Sum Yes TRD %p %p %p %p\n",gSumTRDyesPt,gSumTRDyesPtBC,gSumTRDyesMass,gSumTRDyesWidth);
    
    gSumTRDyesPt   ->SetMarkerStyle(20);
    gSumTRDyesPt   ->SetMarkerColor(kCyan);
    gSumTRDyesPt   ->SetLineColor  (kCyan);    
    
    gSumTRDyesMass->SetMarkerStyle(20);
    gSumTRDyesMass->SetMarkerColor(kCyan);
    gSumTRDyesMass->SetLineColor  (kCyan);   
    
    gSumTRDyesWidth->SetMarkerStyle(20);
    gSumTRDyesWidth->SetMarkerColor(kCyan);
    gSumTRDyesWidth->SetLineColor  (kCyan);
    
    if ( kMix )
    {
      gSumTRDyesPtBC ->SetMarkerStyle(20);
      gSumTRDyesPtBC ->SetMarkerColor(kCyan);
      gSumTRDyesPtBC ->SetLineColor  (kCyan);   
    }
    
    legend->AddEntry(gSumTRDyesPt,"w/ TRD","P");
    
    Int_t nCovered = 10-kFirstTRDSM;
    for (Int_t ibin = 0; ibin < gSumPt->GetN(); ibin++) 
    {
      gSumTRDyesPt->GetY ()[ibin] /= nCovered;  
      gSumTRDyesPt->GetEY()[ibin] /= nCovered;  
      if ( kMix )
      {
        gSumTRDyesPtBC->GetY ()[ibin] /= nCovered; 
        gSumTRDyesPtBC->GetEY()[ibin] /= nCovered; 
      }
    }
    
  } // TRD
  
  //
  // Get different combinations 
  //
  TGraphErrors* gPt   [nCombi];
  TGraphErrors* gPtBC [nCombi];
  TGraphErrors* gMass [nCombi];
  TGraphErrors* gWidth[nCombi];
  
  for(Int_t icomb = 0; icomb < nCombi; icomb ++)
  {
    gPt   [icomb] = (TGraphErrors*) fout->Get(Form("gPt_%s%d"   ,opt.Data(),icomb+first));
    gPtBC [icomb] = (TGraphErrors*) fout->Get(Form("gPtBC_%s%d" ,opt.Data(),icomb+first));
    gMass [icomb] = (TGraphErrors*) fout->Get(Form("gMass_%s%d" ,opt.Data(),icomb+first));
    gWidth[icomb] = (TGraphErrors*) fout->Get(Form("gWidth_%s%d",opt.Data(),icomb+first));
    //printf("\t %d %p %p %p %p\n",icomb,gPt[icomb],gPtBC[icomb],gMass[icomb],gWidth[icomb]);
    
    if ( !gPt[icomb] ) continue ;
    
    nCombiActive++;
    
    gPt   [icomb]->SetMarkerStyle(modStyleIndex[icomb]);
    gPt   [icomb]->SetMarkerColor(modColorIndex[icomb]);
    gPt   [icomb]->SetLineColor  (modColorIndex[icomb]);    
    
    gMass[icomb]->SetMarkerStyle(modStyleIndex[icomb]);
    gMass[icomb]->SetMarkerColor(modColorIndex[icomb]);
    gMass[icomb]->SetLineColor  (modColorIndex[icomb]);   
    
    gWidth[icomb]->SetMarkerStyle(modStyleIndex[icomb]);
    gWidth[icomb]->SetMarkerColor(modColorIndex[icomb]);
    gWidth[icomb]->SetLineColor  (modColorIndex[icomb]);
        
    gPt   [icomb]->GetHistogram()->SetTitle(Form("#it{p}_{T} of %s from fit, %s "   ,particleN.Data(), opt.Data()));
    gWidth[icomb]->GetHistogram()->SetTitle(Form("%s mass #sigma vs #it{p}_{T}, %s ",particleN.Data(), opt.Data()));
    gMass [icomb]->GetHistogram()->SetTitle(Form("%s mass #mu vs #it{p}_{T}, %s "   ,particleN.Data(), opt.Data()));
    
    if ( kMix )
    {
      gPtBC [icomb]->SetMarkerStyle(modStyleIndex[icomb]);
      gPtBC [icomb]->SetMarkerColor(modColorIndex[icomb]);
      gPtBC [icomb]->SetLineColor  (modColorIndex[icomb]);   
      gPtBC [icomb]->GetHistogram()->SetTitle(Form("#it{p}_{T} of %s bin counted, %s ",particleN.Data(), opt.Data()));
    }
    
    if ( !opt.Contains("Group") )
    {
//      if(opt == "SM" && kCalorimeter == "DCAL")
//        legend->AddEntry(gPt[icomb],Form("%s %d",opt.Data(),icomb+12),"P");
//      else
        legend->AddEntry(gPt[icomb],Form("%s %d",opt.Data(),icomb+first),"P");
    }
    else
    {
      if(icomb == 0) legend->AddEntry(gPt[icomb],"SM0+4+5+6+7+8+9","P");
      if(icomb == 1) legend->AddEntry(gPt[icomb],"SM1+2","P");
      if(icomb == 2) legend->AddEntry(gPt[icomb],"SM3+7","P");
    }
  }
  
  if ( !gMass[0] ) return ;
  
  //
  // PLOT
  //
  TCanvas *cFitGraph = new TCanvas(Form("cFitGraph_%sCombinations",opt.Data()),
                                   Form("Fit Graphs for %s combinations",opt.Data()),1000,1000);
  cFitGraph->Divide(2, 2);
  
  // Mass  
  cFitGraph->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  
  if ( kParticle == "Pi0" )
  {
    gMass[0]->SetMaximum(141);
    gMass[0]->SetMinimum(116);
  }
 else
 {
   gMass[0]->SetMaximum(600);
   gMass[0]->SetMinimum(450);
 }
  gMass[0]->Draw("AP");
  for(Int_t icomb = 1; icomb < nCombi; icomb ++)
  {
    if ( gMass[icomb] ) gMass[icomb]->Draw("P");
  }
  
  gSumMass->Draw("P");
  if ( kFirstTRDSM > 0 && opt == "SM" )
  {
    gSumTRDnotMass->Draw("P");
    gSumTRDyesMass->Draw("P"); 
  }
  
  legend->Draw();

  cFitGraph->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  
  if ( kParticle == "Pi0" )
  {
    gWidth[0]->SetMaximum(17);
    gWidth[0]->SetMinimum(7);
    if(opt=="Side") gWidth[0]->SetMinimum(2);
  }
  else
  {
    gWidth[0]->SetMaximum(50);
    gWidth[0]->SetMinimum(10);
  }
  
  gWidth[0]->Draw("AP");
  for(Int_t icomb = 1; icomb < nCombi; icomb ++)
  {
   if ( gWidth[icomb] ) gWidth[icomb]->Draw("P");
  }
  
  gSumWidth->Draw("P");
  if ( kFirstTRDSM > 0 && opt == "SM" )
  {
    gSumTRDnotWidth->Draw("P");
    gSumTRDyesWidth->Draw("P"); 
  }
  
  legend->Draw();

  cFitGraph->cd(3);

  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogy();
  
  gPt[0]->SetMaximum(1);
  gPt[0]->SetMinimum(1e-8);
  
  gPt[0]->Draw("AP");
  for(Int_t icomb = 1; icomb < nCombi; icomb ++)
  {
    if ( gPt[icomb] ) gPt[icomb]->Draw("P");
  }
  
  gSumPt->Draw("P");
  if ( kFirstTRDSM > 0 && opt == "SM" )
  {
    gSumTRDnotPt->Draw("P");
    gSumTRDyesPt->Draw("P"); 
  }
  
  legend->Draw();

  if ( kMix )
  {
    cFitGraph->cd(4);

    gPad->SetGridx();
    gPad->SetGridy();
    gPad->SetLogy();
    
    gPtBC[0]->SetMaximum(1);
    gPtBC[0]->SetMinimum(1e-8);
    
    gPtBC[0]->Draw("AP");
    
    for(Int_t icomb = 1; icomb < nCombi; icomb ++)
    {
      if ( gPtBC[icomb] ) gPtBC[icomb]->Draw("P");
    }
    
    gSumPtBC->Draw("P");
    if ( kFirstTRDSM > 0 && opt == "SM" )
    {
      gSumTRDnotPtBC->Draw("P");
      gSumTRDyesPtBC->Draw("P"); 
    }
    legend->Draw();
  }
  
  cFitGraph->Print(Form("%s/%s/%s/%s/MassWidthPtSpectra_%sCombinations.%s",
                        kProdName.Data(),kCalorimeter.Data(),kHistoList.Data(),kParticle.Data(),opt.Data(),
                        /*kFileName.Data(),*/kPlotFormat.Data()));
  //
  // Calculate the ratio to the sum
  //
  TGraphErrors* gRatPt   [nCombi];
  TGraphErrors* gRatPtBC [nCombi];
  TGraphErrors* gRatMass [nCombi];
  TGraphErrors* gRatWidth[nCombi];
  
  TLegend * legendR =  new TLegend(xmin,ymin,xmax,ymax);
  legendR->SetTextSize(0.05);
  
  for(Int_t icomb = 0; icomb < nCombi; icomb ++)
  {
    //printf("icomb %d\n",icomb);
    //printf("\t pt %p %p\n",gPt[icomb],gSumPt);
    gRatPt   [icomb] = DivideGraphs(gPt   [icomb],gSumPt    ); // PlotUtils.C
    //printf("\t ptBC %p %p\n",gPtBC[icomb],gSumPtBC);
    gRatPtBC [icomb] = DivideGraphs(gPtBC [icomb],gSumPtBC  ); // PlotUtils.C
    //printf("\t pt %p %p\n",gMass[icomb],gSumMass);
    gRatMass [icomb] = DivideGraphs(gMass [icomb],gSumMass  ); // PlotUtils.C
    //printf("\t pt %p %p\n",gWidth[icomb],gSumWidth);
    gRatWidth[icomb] = DivideGraphs(gWidth[icomb],gSumWidth);  // PlotUtils.C
    
    if(!gRatMass[icomb]) continue;
    
    gRatPt   [icomb]->SetMarkerStyle(modStyleIndex[icomb]);
    gRatPt   [icomb]->SetMarkerColor(modColorIndex[icomb]);
    gRatPt   [icomb]->SetLineColor  (modColorIndex[icomb]);    
    
    gRatMass[icomb]->SetMarkerStyle(modStyleIndex[icomb]);
    gRatMass[icomb]->SetMarkerColor(modColorIndex[icomb]);
    gRatMass[icomb]->SetLineColor  (modColorIndex[icomb]);   
    
    gRatWidth[icomb]->SetMarkerStyle(modStyleIndex[icomb]);
    gRatWidth[icomb]->SetMarkerColor(modColorIndex[icomb]);
    gRatWidth[icomb]->SetLineColor  (modColorIndex[icomb]);
    
    gRatPt   [icomb]->GetHistogram()->SetTitle(Form("#it{p}_{T} of %s from fit, %s "   ,particleN.Data(), opt.Data()));
    gRatWidth[icomb]->GetHistogram()->SetTitle(Form("%s mass #sigma vs #it{p}_{T}, %s ",particleN.Data(), opt.Data()));
    gRatMass [icomb]->GetHistogram()->SetTitle(Form("%s mass #mu vs #it{p}_{T}, %s "   ,particleN.Data(), opt.Data()));
    
    gRatPt   [icomb]->GetHistogram()->SetYTitle("Ratio to sum");    
    gRatWidth[icomb]->GetHistogram()->SetYTitle("Ratio to sum");    
    gRatMass [icomb]->GetHistogram()->SetYTitle("Ratio to sum");    
    
    gRatPt   [icomb]->GetHistogram()->SetXTitle("#it{p}_{T} (GeV/#it{c})");    
    gRatWidth[icomb]->GetHistogram()->SetXTitle("#it{p}_{T} (GeV/#it{c})");    
    gRatMass [icomb]->GetHistogram()->SetXTitle("#it{p}_{T} (GeV/#it{c})");    

    gRatPt   [icomb]->GetHistogram()->SetTitleOffset(1.4,"Y"); 
    gRatWidth[icomb]->GetHistogram()->SetTitleOffset(1.4,"Y"); 
    gRatMass [icomb]->GetHistogram()->SetTitleOffset(1.4,"Y"); 
    
    if ( kMix )
    {
      gRatPtBC [icomb]->SetMarkerStyle(modStyleIndex[icomb]);
      gRatPtBC [icomb]->SetMarkerColor(modColorIndex[icomb]);
      gRatPtBC [icomb]->SetLineColor  (modColorIndex[icomb]);   
      gRatPtBC [icomb]->GetHistogram()->SetTitle(Form("#it{p}_{T} of %s bin counted, %s ",particleN.Data(), opt.Data()));
      gRatPtBC [icomb]->GetHistogram()->SetYTitle("Ratio to sum");   
      gRatPtBC [icomb]->GetHistogram()->SetXTitle("#it{p}_{T} (GeV/#it{c})");    
      gRatPtBC [icomb]->GetHistogram()->SetTitleOffset(1.5,"Y"); 
    }
    
    if ( !opt.Contains("Group") )
    {
//      if(opt == "SM" && kCalorimeter == "DCAL")
//        legendR->AddEntry(gPt[icomb],Form("%s %d",opt.Data(),icomb+12),"P");
//      else
        legendR->AddEntry(gPt[icomb],Form("%s %d",opt.Data(),icomb+first),"P");
    }
    else
    {
      if(icomb == 0) legendR->AddEntry(gPt[icomb],"SM0+4+5+6+7+8+9","P");
      if(icomb == 1) legendR->AddEntry(gPt[icomb],"SM1+2","P");
      if(icomb == 2) legendR->AddEntry(gPt[icomb],"SM3+7","P");
    }
  } // combi loop

  //
  // PLOT
  //
  TCanvas *cFitGraphRatio = new TCanvas(Form("cFitGraphSumRatio_%sCombinations",opt.Data()),
                                        Form("Fit Graphs ratio to sum for %s combinations",opt.Data()),1000,1000);
  cFitGraphRatio->Divide(2, 2);
  
  // Mass  
  cFitGraphRatio->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  
  gRatMass[0]->SetMaximum(1.03);
  gRatMass[0]->SetMinimum(0.97);
  
  gRatMass[0]->Draw("AP");
  for(Int_t icomb = 1; icomb < nCombi; icomb ++)
  {
    if ( gRatMass[icomb] ) gRatMass[icomb]->Draw("P");
  }
  
  legendR->Draw();
  
  cFitGraphRatio->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  
  gRatWidth[0]->SetMaximum(1.5);
  gRatWidth[0]->SetMinimum(0.7);
  
  gRatWidth[0]->Draw("AP");
  for(Int_t icomb = 1; icomb < nCombi; icomb ++)
  {
    if ( gRatWidth[icomb] ) gRatWidth[icomb]->Draw("P");
  }
  
  legendR->Draw();
  
  cFitGraphRatio->cd(3);
  
  gPad->SetGridx();
  gPad->SetGridy();
//  gPad->SetLogy();
  
  gRatPt[0]->SetMaximum(2);
  gRatPt[0]->SetMinimum(0);
  
  gRatPt[0]->Draw("AP");
  for(Int_t icomb = 1; icomb < nCombi; icomb ++)
  {
    if ( gRatPt[icomb] ) gRatPt[icomb]->Draw("P");
  }
  
  legendR->Draw();
  
  if ( kMix )
  {
    cFitGraphRatio->cd(4);
    
    gPad->SetGridx();
    gPad->SetGridy();
//    gPad->SetLogy();
    
    gRatPtBC[0]->SetMaximum(2);
    gRatPtBC[0]->SetMinimum(0);
    
    gRatPtBC[0]->Draw("AP");
    for(Int_t icomb = 1; icomb < nCombi; icomb ++)
    {
      if ( gRatPtBC[icomb] ) gRatPtBC[icomb]->Draw("P");
    }
    
    legendR->Draw();
  }
  
  cFitGraphRatio->Print(Form("%s/%s/%s/%s/MassWidthPtSpectra_%sCombinations_RatioToSum.%s",
                             kProdName.Data(),kCalorimeter.Data(),kHistoList.Data(),kParticle.Data(),opt.Data(),
                             /*kFileName.Data(),*/kPlotFormat.Data()));
 
  
  // Ratio to corresponding TRD case
  if ( kFirstTRDSM > 0 && opt == "SM" ) // it could be extended for Sector/Side
  {
    //
    // Calculate the ratio to the sum
    //
    TGraphErrors* gRatTRDPt   [nCombi];
    TGraphErrors* gRatTRDPtBC [nCombi];
    TGraphErrors* gRatTRDMass [nCombi];
    TGraphErrors* gRatTRDWidth[nCombi];
    
    TLegend * legendR =  new TLegend(xmin,ymin,xmax,ymax);
    legendR->SetTextSize(0.05);
    
    for(Int_t icomb = 0; icomb < nCombi; icomb ++)
    {
      //printf("icomb %d\n",icomb);
      //printf("\t pt %p %p\n",gPt[icomb],gSumPt);
  
      
      if(icomb < kFirstTRDSM)
      {
        gRatTRDPt   [icomb] = DivideGraphs(gPt   [icomb],gSumTRDnotPt    ); // PlotUtils.C
        gRatTRDPtBC [icomb] = DivideGraphs(gPtBC [icomb],gSumTRDnotPtBC  ); // PlotUtils.C
        gRatTRDMass [icomb] = DivideGraphs(gMass [icomb],gSumTRDnotMass  ); // PlotUtils.C
        gRatTRDWidth[icomb] = DivideGraphs(gWidth[icomb],gSumTRDnotWidth);  // PlotUtils.C
      }
      else
      {
        gRatTRDPt   [icomb] = DivideGraphs(gPt   [icomb],gSumTRDyesPt    ); // PlotUtils.C
        gRatTRDPtBC [icomb] = DivideGraphs(gPtBC [icomb],gSumTRDyesPtBC  ); // PlotUtils.C
        gRatTRDMass [icomb] = DivideGraphs(gMass [icomb],gSumTRDyesMass  ); // PlotUtils.C
        gRatTRDWidth[icomb] = DivideGraphs(gWidth[icomb],gSumTRDyesWidth);  // PlotUtils.C
      }
        
      
      if(!gRatTRDMass[icomb]) continue;
      
      gRatTRDPt   [icomb]->SetMarkerStyle(modStyleIndex[icomb]);
      gRatTRDPt   [icomb]->SetMarkerColor(modColorIndex[icomb]);
      gRatTRDPt   [icomb]->SetLineColor  (modColorIndex[icomb]);    
      
      gRatTRDMass[icomb]->SetMarkerStyle(modStyleIndex[icomb]);
      gRatTRDMass[icomb]->SetMarkerColor(modColorIndex[icomb]);
      gRatTRDMass[icomb]->SetLineColor  (modColorIndex[icomb]);   
      
      gRatTRDWidth[icomb]->SetMarkerStyle(modStyleIndex[icomb]);
      gRatTRDWidth[icomb]->SetMarkerColor(modColorIndex[icomb]);
      gRatTRDWidth[icomb]->SetLineColor  (modColorIndex[icomb]);
      
      gRatTRDPt   [icomb]->GetHistogram()->SetTitle(Form("#it{p}_{T} of %s from fit, %s "   ,particleN.Data(), opt.Data()));
      gRatTRDWidth[icomb]->GetHistogram()->SetTitle(Form("%s mass #sigma vs #it{p}_{T}, %s ",particleN.Data(), opt.Data()));
      gRatTRDMass [icomb]->GetHistogram()->SetTitle(Form("%s mass #mu vs #it{p}_{T}, %s "   ,particleN.Data(), opt.Data()));
      
      gRatTRDPt   [icomb]->GetHistogram()->SetYTitle("Ratio to sum");    
      gRatTRDWidth[icomb]->GetHistogram()->SetYTitle("Ratio to sum");    
      gRatTRDMass [icomb]->GetHistogram()->SetYTitle("Ratio to sum");    
      
      gRatTRDPt   [icomb]->GetHistogram()->SetXTitle("#it{p}_{T} (GeV/#it{c})");    
      gRatTRDWidth[icomb]->GetHistogram()->SetXTitle("#it{p}_{T} (GeV/#it{c})");    
      gRatTRDMass [icomb]->GetHistogram()->SetXTitle("#it{p}_{T} (GeV/#it{c})");    
      
      gRatTRDPt   [icomb]->GetHistogram()->SetTitleOffset(1.4,"Y"); 
      gRatTRDWidth[icomb]->GetHistogram()->SetTitleOffset(1.4,"Y"); 
      gRatTRDMass [icomb]->GetHistogram()->SetTitleOffset(1.4,"Y"); 
      
      if ( kMix )
      {
        gRatTRDPtBC [icomb]->SetMarkerStyle(modStyleIndex[icomb]);
        gRatTRDPtBC [icomb]->SetMarkerColor(modColorIndex[icomb]);
        gRatTRDPtBC [icomb]->SetLineColor  (modColorIndex[icomb]);   
        gRatTRDPtBC [icomb]->GetHistogram()->SetTitle(Form("#it{p}_{T} of %s bin counted, %s ",particleN.Data(), opt.Data()));
        gRatTRDPtBC [icomb]->GetHistogram()->SetYTitle("Ratio to sum");   
        gRatTRDPtBC [icomb]->GetHistogram()->SetXTitle("#it{p}_{T} (GeV/#it{c})");    
        gRatTRDPtBC [icomb]->GetHistogram()->SetTitleOffset(1.5,"Y"); 
      }
      
      if ( !opt.Contains("Group") )
        legendR->AddEntry(gPt[icomb],Form("%s %d",opt.Data(),icomb),"P");
      else
      {
        if(icomb == 0) legendR->AddEntry(gPt[icomb],"SM0+4+5+6+7+8+9","P");
        if(icomb == 1) legendR->AddEntry(gPt[icomb],"SM1+2","P");
        if(icomb == 2) legendR->AddEntry(gPt[icomb],"SM3+7","P");
      }
    } // combi loop
    
    //
    // PLOT
    //
    TCanvas *cFitGraphRatioTRD = new TCanvas(Form("cFitGraphSumRatioTRD_%sCombinations",opt.Data()),
                                             Form("Fit Graphs ratio to sum for %s combinations for TRD cases",opt.Data()),1000,1000);
    cFitGraphRatioTRD->Divide(2, 2);
    
    // Mass  
    cFitGraphRatioTRD->cd(1);
    gPad->SetGridx();
    gPad->SetGridy();
    
    gRatTRDMass[0]->SetMaximum(1.03);
    gRatTRDMass[0]->SetMinimum(0.97);
    
    gRatTRDMass[0]->Draw("AP");
    for(Int_t icomb = 1; icomb < nCombi; icomb ++)
    {
      if ( gRatTRDMass[icomb] ) gRatTRDMass[icomb]->Draw("P");
    }
    
    legendR->Draw();
    
    cFitGraphRatioTRD->cd(2);
    gPad->SetGridx();
    gPad->SetGridy();
    
    gRatTRDWidth[0]->SetMaximum(1.5);
    gRatTRDWidth[0]->SetMinimum(0.7);
    
    gRatTRDWidth[0]->Draw("AP");
    for(Int_t icomb = 1; icomb < nCombi; icomb ++)
    {
      if ( gRatTRDWidth[icomb] ) gRatTRDWidth[icomb]->Draw("P");
    }
    
    legendR->Draw();
    
    cFitGraphRatioTRD->cd(3);
    
    gPad->SetGridx();
    gPad->SetGridy();
    //  gPad->SetLogy();
    
    gRatTRDPt[0]->SetMaximum(2);
    gRatTRDPt[0]->SetMinimum(0);
    
    gRatTRDPt[0]->Draw("AP");
    for(Int_t icomb = 1; icomb < nCombi; icomb ++)
    {
      if ( gRatTRDPt[icomb] ) gRatTRDPt[icomb]->Draw("P");
    }
    
    legendR->Draw();
    
    if ( kMix )
    {
      cFitGraphRatioTRD->cd(4);
      
      gPad->SetGridx();
      gPad->SetGridy();
      //    gPad->SetLogy();
      
      gRatTRDPtBC[0]->SetMaximum(2);
      gRatTRDPtBC[0]->SetMinimum(0);
      
      gRatTRDPtBC[0]->Draw("AP");
      for(Int_t icomb = 1; icomb < nCombi; icomb ++)
      {
        if ( gRatTRDPtBC[icomb] ) gRatTRDPtBC[icomb]->Draw("P");
      }
      
      legendR->Draw();
    }
    
    cFitGraphRatioTRD->Print(Form("%s/%s/%s/%s/MassWidthPtSpectra_%sCombinations_RatioToSumTRD.%s",
                               kProdName.Data(),kCalorimeter.Data(),kHistoList.Data(),kParticle.Data(),opt.Data(),
                               /*kFileName.Data(),*/kPlotFormat.Data()));
    
    
    TLegend * legendSums =  new TLegend(xmin,ymin,xmax,ymax);
 
    
    //
    // PLOT
    //
    TCanvas *cFitGraphSums = new TCanvas(Form("cFitGraph_%sSums",opt.Data()),
                                     Form("Fit Graphs for %s Sums",opt.Data()),1000,1000);
    cFitGraphSums->Divide(2, 2);
    
    // Mass  
    cFitGraphSums->cd(1);
    gPad->SetGridx();
    gPad->SetGridy();
    
    if ( kParticle == "Pi0" )
    {
      gSumMass->SetMaximum(135);
      gSumMass->SetMinimum(125);
    }
    else
    {
      gSumMass->SetMaximum(560);
      gSumMass->SetMinimum(520);
    }
    
    gSumMass->Draw("AP");
    
    gSumTRDnotMass->Draw("P");
    gSumTRDyesMass->Draw("P"); 
    
    legendSums->AddEntry(gSumMass,"Sum","P");
    legendSums->AddEntry(gSumTRDnotMass,"w/o TRD","P");
    legendSums->AddEntry(gSumTRDyesMass,"w/ TRD","P");

    legendSums->Draw();
    
    cFitGraphSums->cd(2);
    gPad->SetGridx();
    gPad->SetGridy();
    
    if ( kParticle == "Pi0" )
    {
      gSumWidth->SetMaximum(17);
      gSumWidth->SetMinimum(7);
      if(opt=="Side") gSumWidth->SetMinimum(2);
    }
    else
    {
      gSumWidth->SetMaximum(50);
      gSumWidth->SetMinimum(10);
    }
        
    gSumWidth->Draw("AP");
    gSumTRDnotWidth->Draw("P");
    gSumTRDyesWidth->Draw("P"); 
    
    legendSums->Draw();
    
    cFitGraphSums->cd(3);
    
    gPad->SetGridx();
    gPad->SetGridy();
    gPad->SetLogy();
    
    gSumPt->SetMaximum(1);
    gSumPt->SetMinimum(1e-8);
    
    gSumPt->Draw("AP");
    gSumTRDnotPt->Draw("P");
    gSumTRDyesPt->Draw("P"); 
  
    legendSums->Draw();
    
    if ( kMix )
    {
      cFitGraphSums->cd(4);
      
      gPad->SetGridx();
      gPad->SetGridy();
      gPad->SetLogy();
      
      gSumPtBC->SetMaximum(1);
      gSumPtBC->SetMinimum(1e-8);
      
      gSumPtBC->Draw("AP");
      gSumTRDnotPtBC->Draw("P");
      gSumTRDyesPtBC->Draw("P"); 
      
      legend->Draw();
    }
    
    cFitGraphSums->Print(Form("%s/%s/%s/%s/MassWidthPtSpectra_%sSums.%s",
                          kProdName.Data(),kCalorimeter.Data(),kHistoList.Data(),kParticle.Data(),opt.Data(),
                          /*kFileName.Data(),*/kPlotFormat.Data()));
  }
}

///
/// Get primary generated spectra and detector acceptance
/// Plot them
///
//-----------------------------------------------------------------------------
void PrimarySpectra(Int_t nPt, TArrayD xPtLimits, TArrayD xPt, TArrayD exPt)
{
  gStyle->SetOptTitle(1);
  gStyle->SetTitleOffset(2,"Y");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(000000);
  
  TArrayD   primMesonPt   ;   primMesonPt   .Set(nPt); 
  TArrayD  eprimMesonPt   ;  eprimMesonPt   .Set(nPt); 
  TArrayD   primMesonPtAcc;   primMesonPtAcc.Set(nPt); 
  TArrayD  eprimMesonPtAcc;  eprimMesonPtAcc.Set(nPt); 
  TArrayD   primMesonAcc  ;   primMesonAcc.Set(nPt); 
  TArrayD  eprimMesonAcc  ;  eprimMesonAcc.Set(nPt); 
  
  TH2D* h2PrimMesonPtY    = (TH2D*) fil->Get( Form("%s_hPrim%sRapidity"   ,kHistoStartName.Data(), kParticle.Data() ) );
  TH2D* h2PrimMesonPtYAcc = (TH2D*) fil->Get( Form("%s_hPrim%sAccRapidity",kHistoStartName.Data(), kParticle.Data() ) );
  
  if ( !h2PrimMesonPtY ) return ;
  
  if(kSumw2)
  {
    h2PrimMesonPtY   ->Sumw2();
    h2PrimMesonPtYAcc->Sumw2();
  }
  
//  if ( nEvt < 1 )
//  {
//    h2PrimMesonPtY   ->Scale(1e10);
//    h2PrimMesonPtYAcc->Scale(1e10);
//  }
//  else
  {
    h2PrimMesonPtY   ->Scale(1./nEvt);
    h2PrimMesonPtYAcc->Scale(1./nEvt);
  }
  
  h2PrimMesonPtY   ->Write();
  h2PrimMesonPtYAcc->Write();
  
  TH1D* hY = h2PrimMesonPtY->ProjectionY("Rapidity",-1,-1);
  Int_t binYmin = hY->FindBin(-0.65);
  Int_t binYmax = hY->FindBin(+0.65);
  //printf("Y bin min %d, max %d\n",binYmin,binYmax);
  TH1D* hPrimMesonPt = (TH1D*) h2PrimMesonPtY->ProjectionX("PrimaryPt" ,binYmin   ,binYmax   );
  hPrimMesonPt->Write();
  
  TH1D* hYAcc = h2PrimMesonPtYAcc->ProjectionY("RapidityAcc",-1,-1);
  binYmin = hYAcc->FindBin(-0.65);
  binYmax = hYAcc->FindBin(+0.65);
  //printf("Y bin min %d, max %d\n",binYmin,binYmax);
  TH1D* hPrimMesonPtAcc = (TH1D*) h2PrimMesonPtYAcc->ProjectionX("PrimaryPtAccepted" ,binYmin   ,binYmax   );
  hPrimMesonPtAcc->Write();

  // PlotUtils.C
  ScaleBinBySize(hPrimMesonPt);
  ScaleBinBySize(hPrimMesonPtAcc);

  Double_t integralA    = 0;
  Double_t integralAErr = 0;
  Double_t integralB    = 0;
  Double_t integralBErr = 0;  
  Double_t ptMin = -1;
  Double_t ptMax = -1;
  Int_t    binMin= -1;
  Int_t    binMax= -1;
  for(Int_t ibin = 0; ibin < nPt; ibin++ )
  {
    ptMin  = xPtLimits[ibin];
    ptMax  = xPtLimits[ibin+1];

    binMin = hPrimMesonPt->FindBin(ptMin);
    binMax = hPrimMesonPt->FindBin(ptMax)-1;
    
    // PlotUtils.C
    GetRangeIntegralAndError(hPrimMesonPt   , binMin, binMax, integralA, integralAErr );
    GetRangeIntegralAndError(hPrimMesonPtAcc, binMin, binMax, integralB, integralBErr );
    
     primMesonPt   [ibin] = integralA;
    eprimMesonPt   [ibin] = integralAErr;
     primMesonPtAcc[ibin] = integralB;
    eprimMesonPtAcc[ibin] = integralBErr;
    
    if ( integralA > 0 && integralB > 0 )
    {
       primMesonAcc[ibin] = integralB / integralA ;
      eprimMesonAcc[ibin] = GetFractionError(integralB,integralA,integralBErr,integralAErr); // PlotUtils.C
    }
    else
    {
       primMesonAcc[ibin] = 0;
      eprimMesonAcc[ibin] = 0;
    }
    
//    printf("Bin %d, [%1.1f,%1.1f] bin size %1.1f, content num %2.2e, content den %2.2e: frac %2.3f+%2.3f\n",
//           ibin,xPtLimits[ibin],xPtLimits[ibin+1],2*exPt[ibin],primMesonPtAcc[ibin],primMesonPt[ibin],primMesonAcc[ibin],eprimMesonAcc[ibin] );
    
    primMesonPt   [ibin]/=((exPt[ibin]*4)); // It should be 2 not 4???
    primMesonPtAcc[ibin]/=((exPt[ibin]*4));
    
    //        // Scale biased production
    //        if(kFileName.Contains("pp_7TeV_Pi0"))
    //        {
    //          primMeson[ibin]*=1.e6;
    //          primMeson[ibin]*=360./100.;
    //          
    //          //eprimMeson[ibin]*=1.e6*1.e6;
    //          eprimMeson[ibin]*=360./100.;
    //          
    //          mesonPtBC [0][ibin] /= (0.62 +xPt[ibin]*0.0017 );
    //          emesonPtBC[0][ibin] /= (0.62 +xPt[ibin]*0.0017 );
    //          
    //          primMeson [ibin] /= (0.56 +xPt[ibin]*0.0096 );
    //          eprimMeson[ibin] /= (0.56 +xPt[ibin]*0.0096 );
    //        }
  } // pT bin loop
  
  TGraphErrors * gPrim =  new TGraphErrors(nPt,xPt.GetArray(),primMesonPt.GetArray(),exPt.GetArray(),eprimMesonPt.GetArray());
  gPrim->SetName("Primary");
  gPrim->GetHistogram()->SetYTitle("d #it{N}/d #it{p}_{T}");
  gPrim->GetHistogram()->SetTitleOffset(1.4,"Y"); 
  gPrim->GetHistogram()->SetXTitle("#it{p}_{T,gen} (GeV/#it{c})");
  gPrim->GetHistogram()->SetTitleOffset(1.2,"X");
  gPrim->GetHistogram()->SetTitle("#it{p}_{T} spectra for |#it{Y}| < 0.65");
  gPrim->Write();
  
  TGraphErrors * gPrimAcc =  new TGraphErrors(nPt,xPt.GetArray(),primMesonPtAcc.GetArray(),exPt.GetArray(),eprimMesonPtAcc.GetArray());
  gPrimAcc->SetName("PrimaryInAcceptance");
  gPrimAcc->GetHistogram()->SetYTitle("d#it{N}/d #it{p}_{T}");
  gPrimAcc->GetHistogram()->SetTitle(Form("|#it{Y}| < 0.65 in %s",kCalorimeter.Data()));
  gPrimAcc->Write();    

  // Acceptance
  TH1D* hAcc = (TH1D*) hPrimMesonPtAcc->Clone("hAcceptance");
  hAcc->Divide(hPrimMesonPt);
  
  TGraphErrors * gAcc =  new TGraphErrors(nPt,xPt.GetArray(),primMesonAcc.GetArray(),exPt.GetArray(),eprimMesonAcc.GetArray());
  gAcc->SetName("Acceptance");
  gAcc->GetHistogram()->SetYTitle("Acceptance");
  gAcc->GetHistogram()->SetTitleOffset(1.4,"Y"); 
  gAcc->GetHistogram()->SetXTitle("#it{p}_{T,gen} (GeV/#it{c})");
  gAcc->GetHistogram()->SetTitleOffset(1.2,"X");
  gAcc->GetHistogram()->SetTitle(Form("Acceptance for |#it{Y}| < 0.65 in %s",kCalorimeter.Data()));
  gAcc->Write();    
 
  // Plot spectra and acceptance
  //
  TCanvas *cAcc = new TCanvas(Form("cAcceptance"),
                              Form("Primary generated p_{T} spectra and acceptance"),2*600,600);
  cAcc->Divide(2, 1);
  
  cAcc->cd(1);
  //gPad->SetGridx();
  //gPad->SetGridy();
  gPad->SetLogy();
  
  Double_t maximum = gPrim->GetHistogram()->GetMaximum();
  Double_t minimum = gPrim->GetHistogram()->GetMinimum();

  if ( gPrimAcc->GetMaximum() > maximum ) maximum = gPrimAcc->GetMaximum() ;
  if ( gPrimAcc->GetMinimum() > minimum ) minimum = gPrimAcc->GetMinimum() ;
  
  gPrim->SetMaximum(maximum*10);
  gPrim->SetMinimum(minimum/10);
  
  gPrim   ->Draw("AP");
  gPrimAcc->Draw("P");
  
  gPrim   ->SetMarkerColor(1);
  gPrimAcc->SetMarkerColor(4);  
  gPrim   ->SetLineColor(1);
  gPrimAcc->SetLineColor(4); 
  gPrim   ->SetMarkerStyle(24);
  gPrimAcc->SetMarkerStyle(24);
  
  hPrimMesonPt   ->Draw("same");
  hPrimMesonPtAcc->Draw("same");
 
  hPrimMesonPt   ->Draw("Hsame");
  hPrimMesonPtAcc->Draw("Hsame");
  hPrimMesonPt   ->SetLineColor(1);
  hPrimMesonPtAcc->SetLineColor(4); 
  
  TLegend * legendS =  new TLegend(0.4,0.7,0.9,0.9);
  legendS->AddEntry(gPrim,"|Y| < 0.65","P");
  legendS->AddEntry(gPrimAcc,Form("Both in %s",kCalorimeter.Data()),"P"); 
  
  legendS->AddEntry(hPrimMesonPt,"Histo |Y| < 0.65","L");
  legendS->AddEntry(hPrimMesonPtAcc,Form("Histo Both in %s",kCalorimeter.Data()),"L");
  legendS->Draw();
  
  cAcc->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  //gPad->SetLogy();
  
  //gAcc->SetMaximum(1);
  gAcc->SetMaximum(gAcc->GetHistogram()->GetMaximum()*1.3);
  gAcc->SetMinimum(0);
  
  gAcc->Draw("AP");
  gAcc->SetMarkerColor(1);
  gAcc->SetMarkerStyle(24);
  
  hAcc->Draw("Hsame");
  hAcc->SetLineColor(4);
  
  TLegend * legendA =  new TLegend(0.7,0.75,0.9,0.9);
  legendA->AddEntry(gAcc,"Graph","P");
  legendA->AddEntry(hAcc,"Histo","L");

  legendA->Draw();
  
  cAcc->Print(Form("%s/%s/%s/%s/PrimarySpectraAcceptance.%s",
                   kProdName.Data(),kCalorimeter.Data(),kHistoList.Data(),kParticle.Data(),
                   /*kFileName.Data(),*/ kPlotFormat.Data()));   
}
  
///
/// Main method
///
/// \param prodname   : name of directory with histogram file
/// \param filename   : histogram file name
/// \param histoDir   : TDirectoryFile folder name
/// \param histoList  : TList folder name
/// \param histoStart : begining of histograms name
/// \param calorimeter: "EMCAL","DCAL"
/// \param particle   : "Pi0","Eta", define fitting and plotting ranges for particle
/// \param nPairMin   : minimum number of entries un the pi0 or eta peak integral to do the fits and plotting. Careful in MC scaled productions.
/// \param sumw2      : bool, apply the histograms Sumw2 method to get proper errors when not done before
/// \param rebin      : int, rebin default Inv Mass histograms binning with this value
/// \param firstTRD   : int, First EMCal SM covered by a TRD SM, in 2011 and 4 in 2012-13
/// \param mixed      : bool, use mixed event to constrain combinatorial background
/// \param truncmix   : bool, 1: use truncated function to constrain backgroun; 0: use factor from fixed bin
/// \param pol        : int, polinomyal type for residual background under the peak
/// \param drawAll    : bool, activate plot per pair in any SM
/// \param drawMCAll  : bool, activate plot per pair in any SM and MC
/// \param drawPerSM  : bool, activate plot per pair in same SM
/// \param drawPerSMGr: bool, activate plot per pair in same SM adding SM: 3 groups [3,7], [1,2], [0,4,5,6,8,9], needs drawPerSM=kTRUE
/// \param drawPerSector: bool, activate plot per pair in different continuous SM sectors
/// \param drawPerSide: bool, activate plot per pair in different continuous SM sides 
/// \param plotFormat : define the type of figures: eps, pdf, etc.
///
//-----------------------------------------------------------------------------
void InvMassFit
(TString prodname     = "LHC18c3_NystromOn",//"LHC17l3b_fast", 
 TString filename     = "AnalysisResults",
 TString histoDir     = "Pi0IM_GammaTrackCorr_EMCAL", 
 TString histoList    = "default", 
 TString histoStart   = "",
 TString calorimeter  = "EMCAL",
 TString particle     = "Pi0",
 Float_t nPairMin     = 20, 
 Bool_t  sumw2        = kTRUE,
 Int_t   rebin        = 1,
 Int_t   firstTRD     = 0,
 Bool_t  mixed        = kTRUE,
 Bool_t  truncmix     = kFALSE,
 Int_t   pol          = 1,
 Bool_t  drawAll      = kTRUE,
 Bool_t  drawMCAll    = kFALSE,
 Bool_t  drawPerSM    = kTRUE, 
 Bool_t  drawPerSMGr  = kFALSE, 
 Bool_t  drawPerSector= kFALSE, 
 Bool_t  drawPerSide  = kFALSE, 
 TString plotFormat   = "eps")
{
  // Load plotting utilities
  //gROOT->Macro("$ALICE_PHYSICS/PWGGA/CaloTrackCorrelations/macros/PlotUtils.C"); // Already in include
  //gROOT->Macro("$MACROS/style.C");//Set different root style parameters

  kSumw2      = sumw2;
  kPolN       = pol;
  kMix        = mixed;
  kTrunMixFunc= truncmix;
  kParticle   = particle;
  kProdName   = prodname;
  kFileName   = filename;
  kPlotFormat = plotFormat;
  kCalorimeter= calorimeter;
  kNPairCut   = nPairMin;
  kRebin      = rebin;
  kFirstTRDSM = firstTRD;
  kHistoList  = histoList;
  
  if(histoStart !="")
    kHistoStartName = histoStart;
  else
  {
    kHistoStartName = "AnaPi0_Calo0";
    if(kCalorimeter=="DCAL") 
      kHistoStartName = "AnaPi0_Calo1";
  }
  
  if(drawPerSMGr) drawPerSM = kTRUE;
  
  //---------------------------------------
  // Get input file
  //---------------------------------------
  Int_t ok = GetFileAndEvents(histoDir, histoList/*, calorimeter*/);
  
  // Create output directories
  TString processline = Form(".! mkdir -p %s/%s/%s/%s",
                             prodname.Data(),calorimeter.Data(),histoList.Data(),particle.Data()) ;
  gROOT->ProcessLine(processline.Data());
  
//  kProdName.ReplaceAll("module/","");
//  kProdName.ReplaceAll("TCardChannel","");
//  kProdName.ReplaceAll("/","_");
//  kProdName.ReplaceAll("2_","_");
//  kProdName.ReplaceAll("__","_");
  if( kFileName !="AnalysisResults" && !kFileName.Contains("Scale") )
    kProdName+=kFileName;
  
  printf("Settings: prodname <%s>, filename <%s>, histoDir <%s>, histoList <%s>,\n"
         " \t histo start name <%s>, particle <%s>, calorimeter <%s>, \n"
         " \t mix %d, kPolN %d, sumw2 %d, n pairs cut %2.2e\n",
         kProdName.Data(), kFileName.Data(),histoDir.Data(),histoList.Data(), kHistoStartName.Data(), 
         kParticle.Data(), kCalorimeter.Data(), kMix, kPolN, kSumw2, kNPairCut);
  
  if ( !ok )
  {
    printf("Could not recover file <%p> or dir <%p> or list <%p>\n",fil,direc,lis);
    return ; 
  }
  
  //  kNPairCut/=nEvt;

  //---------------------------------------
  // Set-up output file with processed histograms
  //---------------------------------------
  
  // Open output file
  fout = new TFile(Form("%s/%s/%s/%s/MassWidthPtHistograms.root",
                        kProdName.Data(),kCalorimeter.Data(),kHistoList.Data(),kParticle.Data() /*,kFileName.Data()*/), "recreate");
  printf("Open Output File: %s/%s/%s/%s/MassWidthPtHistograms.root\n",
         kProdName.Data(),kCalorimeter.Data(),kHistoList.Data(),kParticle.Data());//,kFileName.Data());
  
  //---------------------------------------
  // Set the pt bins and total range
  //---------------------------------------

  const Int_t nPtEta = 13;
  Double_t xPtLimitsEta[] = {2.,3.,4.,5.,6.,8.,10.,14.,18.,24.,28.,34.,40.,50.};
  
//  const Int_t nPtPi0 = 7;
//  Double_t xPtLimitsPi0[] = {1.5,2,2.5,3.,3.5,4.,5.,6.};
  
  const Int_t nPtPi0 = 23;
  Double_t xPtLimitsPi0[] = {1.6,1.8,2.,2.2,2.4,2.6,2.8,3.,3.4,3.8,4.2,4.6,5.,5.5,6.,6.5,7.,7.5,8.,9.,10.,12.,15.,17.,20.};
  
  Int_t nPt = nPtPi0;
  if(kParticle == "Eta") nPt = nPtEta;
  
  TArrayD xPtLimits; xPtLimits.Set(nPt+1); 
  TArrayD xPt      ; xPt      .Set(nPt);
  TArrayD exPt     ; exPt     .Set(nPt);
  
  for(Int_t i = 0; i < nPt; i++)
  {
    if(kParticle == "Pi0")
    {
      xPtLimits.SetAt(xPtLimitsPi0[i],i);
      exPt     .SetAt((xPtLimitsPi0[i+1]-xPtLimitsPi0[i])/2,i);
      xPt      .SetAt(xPtLimitsPi0[i]+exPt[i],i);
    }
    else
    {
      xPtLimits.SetAt(xPtLimitsEta[i],i);
      exPt     .SetAt((xPtLimitsEta[i+1]-xPtLimitsEta[i])/2,i);
      xPt      .SetAt(xPtLimitsEta[i]+exPt[i],i);
    }
    
    //printf("%s pT bin %d, pT %2.2f, dpT %2.2f, limit %2.2f\n",kParticle.Data(),i,xPt[i],exPt[i],xPtLimits[i]);
  }
  // Extra entry
  if ( kParticle == "Pi0" ) xPtLimits.SetAt(xPtLimitsPi0[nPt],nPt);
  else                      xPtLimits.SetAt(xPtLimitsEta[nPt],nPt);
  
  TString comment = "";
  TString leg     = "";
  TString hname   = "";   
 
  //---------------------------------------
  // Fitting function and histograms initialization
  //---------------------------------------
  SetFitFun();
  
  TH2F  * hRe ;
  TH2F  * hMi ;
  
  //=======================================
  // Get histograms, project per pT bin and fit
  // Do it for different kind of histograms: 
  //     Total pairs, MC pairs, pairs per SM combinations
  //=======================================
  
  //---------------------------------------
  // MC input
  //---------------------------------------
  
  if ( drawMCAll )
  {
    // Get the generated primary pi0 spectra, for efficiency calculation 
    PrimarySpectra(nPt,xPtLimits, xPt, exPt);

    // Reconstructed pT of real meson pairs
    //hRe = (TH2F *) fil->Get( Form("%s_hMCOrgMass_2", kHistoStartName.Data()) ); // Pi0
    //hRe = (TH2F *) fil->Get( Form("%s_hMCOrgMass_3", kHistoStartName.Data()) ); // Eta
    hRe = (TH2F *) fil->Get( Form("%s_hMC%sMassPtRec", kHistoStartName.Data(), kParticle.Data()) );
    comment = "MC pT reconstructed";
    leg     = "MC pT reconstructed";
    hname   = "MCpTReco" ; 
    ProjectAndFit(nPt,comment,leg,hname,xPtLimits, xPt, exPt, hRe, 0x0 );    
  
    
    // Reconstructed pT of real eta pairs
    hRe = (TH2F *) fil->Get( Form("%s_hMC%sMassPtTrue", kHistoStartName.Data(), kParticle.Data()) );
    comment = "MC pT generated";
    leg     = "MC pT generated";
    hname   = "MCpTGener" ; 
    ProjectAndFit(nPt,comment,leg,hname,xPtLimits, xPt, exPt, hRe, 0x0 );    
  } // MC treatment
    
  //---------------------------------------
  // All SM
  //---------------------------------------
  
  if ( drawAll )
  {
    if ( !lis )
    {
      hRe = (TH2F *) fil->Get(Form("%s_hRe_cen0_pidbit0_asy0_dist1", kHistoStartName.Data()));
      hMi = (TH2F *) fil->Get(Form("%s_hMi_cen0_pidbit0_asy0_dist1", kHistoStartName.Data()));
    }
    else
    {
      //printf("Find: %s_hRe_cen0_pidbit0_asy0_dist1", kHistoStartName.Data());
      hRe = (TH2F *) lis->FindObject(Form("%s_hRe_cen0_pidbit0_asy0_dist1", kHistoStartName.Data()));
      hMi = (TH2F *) lis->FindObject(Form("%s_hMi_cen0_pidbit0_asy0_dist1", kHistoStartName.Data()));
    }
    
    printf("histo Re %p - Mix %p All SM\n",hRe,hMi);
    
    if( !hMi ) kMix = kFALSE; 
    
    comment = "All SM";
    leg     = "All SM";
    hname   = "AllSM";
    ProjectAndFit(nPt,comment,leg,hname,xPtLimits, xPt, exPt, hRe, hMi );    
  }
  
  //-------------------------------------------
  // Histograms for cluster pairs in different 
  // combinations of SM
  //-------------------------------------------
  
  //
  // Pairs in same SM
  // 
  if ( drawPerSM )
  {
    Int_t imod;
    Int_t firstMod = 0;
    Int_t lastMod  = 11;   

    if ( kCalorimeter=="DCAL" )
    {
      firstMod  = 12;
      lastMod   = 19;
    }
    
    // Initialize SM grouping histograms
    // 0 : SM 0+4+5+6+8+9
    // 1 : SM 1+2 
    // 2 : SM 3+7 
    TH2F  * hReSMGroups[3] ;
    TH2F  * hMiSMGroups[3] ;
    TH2F  * hReSM = 0;
    TH2F  * hMiSM = 0;   
    TH2F  * hReSMTRDYes = 0;
    TH2F  * hMiSMTRDYes = 0;   
    TH2F  * hReSMTRDNot = 0;
    TH2F  * hMiSMTRDNot = 0;
    
    for(Int_t imodgroup = 0; imodgroup<=2; imodgroup++)
    {
      hReSMGroups[imodgroup] = hMiSMGroups[imodgroup] = 0x0;
    }
    
    for(imod = firstMod; imod<=lastMod; imod++)
    {
      if(!lis)
      {
        hRe = (TH2F *) fil->Get(Form("%s_hReMod_%d", kHistoStartName.Data(), imod));
        hMi = (TH2F *) fil->Get(Form("%s_hMiMod_%d", kHistoStartName.Data(), imod));
      }
      else
      {
        hRe = (TH2F *) lis->FindObject(Form("%s_hReMod_%d", kHistoStartName.Data(), imod));
        hMi = (TH2F *) lis->FindObject(Form("%s_hMiMod_%d", kHistoStartName.Data(), imod));
      }
      
      //printf("histo mod %d Re %p - Mix %p\n",imod, hRe,hMi);
      
      comment = Form("both clusters in SM %d",imod);
      leg     = Form("SM %d",imod);
      hname   = Form("SM%d" ,imod);
      
      ProjectAndFit(nPt,comment,leg,hname,xPtLimits, xPt, exPt, hRe, hMi );    
      
      // Add all SM
      if( imod == firstMod )
      {
        hReSM = (TH2F*) hRe->Clone("h2D_Re_IM_SM");
        if(kMix) hMiSM = (TH2F*) hMi->Clone("h2D_Mi_IM_SM");
      }
      else
      {
        hReSM->Add( hRe );
        if(kMix) hMiSM ->Add( hMi ); 
      }
      
      // Add TRD covered or not SM
      if( kFirstTRDSM > 0 )
      {
        if ( imod >= kFirstTRDSM && imod < 10 )
        {
          if( imod == kFirstTRDSM )
          {
            hReSMTRDYes = (TH2F*) hRe->Clone("h2D_Re_IM_SM_TRD_Yes");
            if(kMix) hMiSMTRDYes = (TH2F*) hMi->Clone("h2D_Mi_IM_SM_TRD_Yes");
          }
          else
          {
            hReSMTRDYes->Add( hRe );
            if(kMix) hMiSMTRDYes ->Add( hMi ); 
          }
        }
      }
      
      // Group SMs in with particular behaviors
      if ( drawPerSMGr )
      {
        if      ( imod == 0 )
        {
          hReSMGroups[0] = (TH2F*) hRe->Clone("h2D_Re_IM_SM045689");
          if(kMix) 
            hMiSMGroups[0] = (TH2F*) hMi->Clone("h2D_Mi_IM_SM045689");
        }
        else if ( imod == 1 )
        {
          hReSMGroups[1] = (TH2F*) hRe->Clone("h2D_Re_IM_SM12");
          if(kMix) hMiSMGroups[1] = (TH2F*) hMi->Clone("h2D_Mi_IM_SM12");
        }
        else if ( imod == 3 )
        {
          hReSMGroups[2] = (TH2F*) hRe->Clone("h2D_Re_IM_SM37");
          if(kMix) hMiSMGroups[2] = (TH2F*) hMi->Clone("h2D_Mi_IM_SM37");
        }     
        else if ( imod == 2 )
        {
          hReSMGroups[1]->Add(hRe);
          if(kMix) hMiSMGroups[1]->Add(hMi);
        }     
        else if ( imod == 7 )
        {
          hReSMGroups[2]->Add(hRe);
          if(kMix) hMiSMGroups[2]->Add(hMi);
        }
        else if ( imod < 10 )
        {
          hReSMGroups[0]->Add(hRe);
          if(kMix) hMiSMGroups[0]->Add(hMi);
        }
      } // do SM groups addition
    } // SM loop
    
    // Sum of pairs in same SM   
    comment = "both clusters in same SM";
    leg     = "Same SM";
    hname   = "SameSM";
    ProjectAndFit(nPt,comment,leg,hname,xPtLimits, xPt, exPt, hReSM, hMiSM );    

    if ( kFirstTRDSM > 0 )
    {
      // Sum of pairs in same SM  with TRD
      comment = "both clusters in same SM, TRD covered";
      leg     = "Same SM, TRD Yes";
      hname   = "SameSMTRDYes";
      ProjectAndFit(nPt,comment,leg,hname,xPtLimits, xPt, exPt, hReSMTRDYes, hMiSMTRDYes );    
      
      // Sum of pairs in same SM without TRD
      hReSMTRDNot = (TH2F*) hReSM->Clone("h2D_Re_IM_SM_TRD_Not");
      hReSMTRDNot->Add(hReSMTRDYes,-1);
      if(kMix)
      {
        hMiSMTRDNot = (TH2F*) hMiSM->Clone("h2D_Mi_IM_SM_TRD_Not");
        hMiSMTRDNot->Add(hMiSMTRDYes,-1);
      }
      comment = "both clusters in same SM, NOT TRD covered";
      leg     = "Same SM, TRD No";
      hname   = "SameSMTRDNot";
      ProjectAndFit(nPt,comment,leg,hname,xPtLimits, xPt, exPt, hReSMTRDNot, hMiSMTRDNot );    
    }
    
    // Fit parameters final plotting    
    PlotGraphs("SM"     ,firstMod   ,lastMod   );
    
    // Group SMs in with particular behaviors
    if( drawPerSMGr )
    {
      for(Int_t imodgroup = 0; imodgroup<=2; imodgroup++)
      {
        //printf("histo modgroup %d Re %p - Mix %p\n",
        //       imod, hReSMGroups[imodgroup], hMiSMGroups[imodgroup]);
        
        comment = Form("both clusters in same SM, group %d",imodgroup);
        leg     = Form("SMGroup %d",imodgroup);
        hname   = Form("SMGroup%d" ,imodgroup);
        
        if ( imodgroup != 0 )
        {
          hReSMGroups[imodgroup]->Scale(1./2.);
          if ( kMix ) hMiSMGroups[imodgroup]->Scale(1./2.);
        }
        else
        {
          hReSMGroups[imodgroup]->Scale(1./6.);
          if ( kMix ) hMiSMGroups[imodgroup]->Scale(1./6.);
        }
        
        ProjectAndFit(nPt, comment, leg, hname, xPtLimits, xPt, exPt, 
                      hReSMGroups[imodgroup], hMiSMGroups[imodgroup] );    
      }
      
      // Fit parameters final plotting    
      PlotGraphs("SMGroup",0          ,2         );
    } // Per SM groups
  
  } // Per SM
  
  //
  // Pairs in same sector
  //
  if(drawPerSector)
  {
    Int_t firstSector = 0;
    Int_t lastSector  = 5;
    if ( kCalorimeter=="DCAL" )
    {
      firstSector = 6;
      lastSector  = 9;
    }
    
    TH2F  * hReSector = 0;
    TH2F  * hMiSector = 0;
    Int_t isector;
    for(isector = firstSector; isector<=lastSector; isector++)
    { 
      if(!lis)
      {
        hRe = (TH2F *) fil->Get(Form("%s_hReSameSectorEMCALMod_%d", kHistoStartName.Data(), isector));
        hMi = (TH2F *) fil->Get(Form("%s_hMiSameSectorEMCALMod_%d", kHistoStartName.Data(), isector));
      }     
      else 
      {
        hRe = (TH2F *) lis->FindObject(Form("%s_hReSameSectorEMCALMod_%d", kHistoStartName.Data(), isector));
        hMi = (TH2F *) lis->FindObject(Form("%s_hMiSameSectorEMCALMod_%d", kHistoStartName.Data(), isector));
      }
      
      // Add all Sectors
      if( isector == firstSector )
      {
        hReSector = (TH2F*) hRe->Clone("h2D_Re_IM_Sector");
        if(kMix) hMiSector = (TH2F*) hMi->Clone("h2D_Mi_IM_Sector");
      }
      else
      {
        hReSector->Add( hRe );
        if(kMix) hMiSector ->Add( hMi ); 
      }
      
      //printf("histo sector %d Re %p - Mix %p\n",isector, hRe,hMi);
      
      comment = Form("both clusters in Sector %d",isector);
      leg     = Form("Sector %d",isector);
      hname   = Form("Sector%d" ,isector);
      
      ProjectAndFit(nPt,comment,leg,hname,xPtLimits, xPt, exPt, hRe, hMi );    
    }
    
    // Sum of pairs in same Sector
    comment = "both clusters in same Sector";
    leg     = "Same Sector";
    hname   = "SameSector";
    ProjectAndFit(nPt,comment,leg,hname,xPtLimits, xPt, exPt, hReSector, hMiSector );  
    
    // Fit parameters final plotting
    PlotGraphs("Sector" ,firstSector,lastSector);
  } // Per sector
  
  
  // Pairs in same side
  //
  if(drawPerSide)
  {
    Int_t firstSide = 0;
    Int_t lastSide  = 9;   
    
    if ( kCalorimeter=="DCAL" )
    {
      firstSide = 10;
      lastSide  = 15;   
    }  
    
    TH2F  * hReSide = 0;
    TH2F  * hMiSide = 0;
    Int_t iside;
    for(iside = firstSide; iside <= lastSide; iside++)
    { 
      if(!lis)
      {
        hRe = (TH2F *) fil->Get(Form("%s_hReSameSideEMCALMod_%d", kHistoStartName.Data(), iside));
        hMi = (TH2F *) fil->Get(Form("%s_hMiSameSideEMCALMod_%d", kHistoStartName.Data(), iside));
      }
      else
      {
        hRe = (TH2F *) lis->FindObject(Form("%s_hReSameSideEMCALMod_%d", kHistoStartName.Data(), iside));
        hMi = (TH2F *) lis->FindObject(Form("%s_hMiSameSideEMCALMod_%d", kHistoStartName.Data(), iside));
      }
      
      // Add all Sides
      if( iside == firstSide )
      {
        hReSide = (TH2F*) hRe->Clone("h2D_Re_IM_Side");
        if(kMix) hMiSide = (TH2F*) hMi->Clone("h2D_Mi_IM_Side");
      }
      else
      {
        hReSide->Add( hRe );
        if(kMix) hMiSide ->Add( hMi ); 
      }
      
      //printf("histo side %d Re %p - Mix %p\n",iside, hRe,hMi);
      
      comment = Form("both clusters in Side %d",iside);
      leg     = Form("Side %d",iside);
      hname   = Form("Side%d" ,iside);      
      
      ProjectAndFit(nPt,comment,leg,hname,xPtLimits, xPt, exPt, hRe, hMi );    
    }
    
    // Sum of pairs in same Side   
    comment = "both clusters in same Side";
    leg     = "Same Side";
    hname   = "SameSide";
    ProjectAndFit(nPt,comment,leg,hname,xPtLimits, xPt, exPt, hReSide, hMiSide );    
    
    // Fit parameters final plotting
    PlotGraphs("Side"   ,firstSide  ,lastSide  );
    
  } // per Side
  
  fout->Close();  
}
