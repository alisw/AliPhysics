///
/// \file MakePurityCalculationAndComparisons.C
/// \ingroup CaloTrackCorrMacrosPlotting
/// \brief Calculate isolated photon purity
///
/// Example macro to calculate isolated photon purity for different data/MC productions and 
/// different shower shape selection regions with the ABCD method.
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)
///

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TFile.h>
#include <TDirectoryFile.h>
#include <TList.h>
#include <TString.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TObject.h>
#include <TAxis.h>
#include <TGaxis.h>
#include <TLine.h>
#include <TF1.h>
#include <TMath.h>
#include <TGraphErrors.h>
#include "PlotUtils.C"
#endif

Int_t color    [] = {1,4,2,6,8,9,kYellow-6,kCyan,kOrange+2,kViolet,kOrange-2,4,2,6,8,9,kYellow-6,kCyan,kOrange+2,kViolet,kOrange-2};
Int_t lineStyle[] = {1,1,1,1,1,1,        1,    1,        1,      1,        1,2,2,2,2,2,        2,    2,        2,      2,        2,2,2,2,2,2,2,2,2,2,2,2};
Int_t marker   [] = {24,21,20,25,24,24,24,24,24,24,24};
TString format    = ".eps";

// Main input data file location, set them in MakePurityCalculationAndComparison method
TString opt         = "100MeV_New";
TString dataFileDir = "data3";

//------------------------------------------------------------------------------
/// Get the histograms defining the different purity regions ABDC, 
/// and do the double ratios needed for purity calculation in MakePurity()
/// where data and MC double ratios are combined. This method is used the same
/// for that and MC.
///
/// \param titleProd   : Name of analized production
/// \param gjProd      : Name of gamma-jet analized production
/// \param nCases      : in case of MC, the different cases studied, in case of data it is 1
/// \param mcCasesPath : path to MC/Data case input file
/// \param mcCases     : simplified name of MC/Data sample
/// \param mcLef       : name of MC/Data sample for legends in plots
/// \param addGJ       : add the GJ sample (P3) or not (P2)
/// \param gjFactor    : in case of addition, scale the GJ sample by this factor
/// \param usePi0      : Use as bkg region the pi0 shower shape band
/// \param rebin       : if 0 use the energy binning defined globally, if not, use simple binning number
/// \param m02LowCutMin: Signal region, minimum shower shape value
/// \param m02LowCutMax: Signal region, maximum shower shape value
/// \param m02HighCutMin: Bkg region, minimum shower shape value
/// \param m02HighCutMax: Bkg region, maximum shower shape value
/// \param debug       : Bool to activate printfs.
//------------------------------------------------------------------------------
void MakeDoubleRatios
(
   TString titleProd, TString gjProd
 , Int_t nCases, TString * mcCasesPath 
 , TString * mcCases, TString * mcLeg
 , Bool_t  addGJ        = 1 
 , Float_t gjFactor     = 1 
 , Bool_t  usePi0       = 0
 , Int_t   rebin        = 0
 , Float_t m02LowCutMin = 0.10
 , Float_t m02LowCutMax = 0.27//0.27
 , Float_t m02HigCutMin = 0.4//0.40
 , Float_t m02HigCutMax = 3.00
 , Int_t   debug        = kFALSE
)
{  
  // Do some checks and modify stuff
  //
  
  if ( addGJ && titleProd.Contains("LHC") ) 
  {
    if ( debug )
    {
      printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
      printf("In case of data, not possible to add GJ simulation, do nothing, deactivate this option\n");
      printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
    }
    return;
  }
  
  // Create histograms arrays
  //
  const Int_t nMCCases = nCases;
    
  TFile * fileJJ[nMCCases];
  TFile * fileGJ[nMCCases];
  Int_t firstMC = 0;
  
  TH2D* h2M02Iso  [nMCCases];
  TH2D* h2M02Not  [nMCCases];
  TH2D* h2M02IsoGJ[nMCCases];

  TH1D* hM02LowIs [nMCCases];
  TH1D* hM02LowNo [nMCCases];  
  TH1D* hM02HigIs [nMCCases];
  TH1D* hM02HigNo [nMCCases];

  TH1D* hM02LowIsRat[nMCCases];
  TH1D* hM02LowNoRat[nMCCases];  
  TH1D* hM02HigIsRat[nMCCases];
  TH1D* hM02HigNoRat[nMCCases];
 
  TH1D* hDoubleRatioCases[nMCCases];
  TH1D* hDoubleRatio[nMCCases];
 
  // Init the arrays, just in case
  for(Int_t icase = 0; icase < nMCCases; icase++)
  {
    fileJJ[icase] = 0;
    fileGJ[icase] = 0;
    
    hM02LowIs [icase]=0;
    hM02LowNo [icase]=0;  
    hM02HigIs [icase]=0;
    hM02HigNo [icase]=0;
    
    hM02LowIsRat[icase]=0;
    hM02LowNoRat[icase]=0;  
    hM02HigIsRat[icase]=0;
    hM02HigNoRat[icase]=0;
    
    hDoubleRatio[icase]=0;
    hDoubleRatioCases[icase]=0;
  }
  
  // Define energy binning
  //
  // Rebin (>0)
  const Int_t nBins = 14;
  Double_t binE[]={8,9,10,11,12,13,14,16,18,20,25,30,40,50,60};
  // Histogram limits
  Float_t emin = 8;
  Float_t emax = 60;
  if(titleProd =="Low")
  {
    emin = 8;
    emax = 30;  
  }
  if(titleProd =="High")
  {
    emin = 18;
    emax = 60;  
  }  
  
  // String used later to name the output files that depends on the parameters
  // provided
  TString cuts = "";
  if ( usePi0 ) 
    cuts = Form("Low_%2.2f-%2.2f"                 ,m02LowCutMin, m02LowCutMax);
  else
    cuts = Form("Low_%2.2f-%2.2f_High_%2.2f-%2.2f",m02LowCutMin, m02LowCutMax, m02HigCutMin, m02HigCutMax);
  
  TString prodTitleOutput = titleProd;
  if ( addGJ )  prodTitleOutput+=Form("_GJx%2.2f",gjFactor);
  if ( usePi0 ) prodTitleOutput+="_MergedPi0";
  
  //---------------------------------------
  // Open the files and get the histograms
  //---------------------------------------
  
  if ( debug )
    printf("Get files and histograms \n");
  // Simu
  for(Int_t icase = firstMC; icase < nMCCases; icase++)
  {
     if ( debug )
       printf("\t Case %d - %s\n",icase,  mcCases[icase].Data());
    if ( titleProd.Contains("LHC") ) 
    {
      fileJJ[icase] = TFile::Open(Form("%s/%s/%s.root",opt.Data(),dataFileDir.Data(),titleProd.Data()));
    }
    else
      fileJJ[icase] = TFile::Open(Form("%s/%s/%s/ScaledMerged_%s.root",opt.Data(),mcCasesPath[icase].Data(),titleProd.Data(), mcCases[icase].Data()));
    
    if ( addGJ ) 
      fileGJ[icase] = TFile::Open(Form("%s/%s/%s/ScaledMerged_%s.root",opt.Data(),mcCasesPath[icase].Data(),gjProd.Data(), mcCases[icase].Data()));
    
    if ( debug )
    {
      printf("Data/MC file %s %p\n",titleProd.Data(),fileJJ[icase]);
      printf("GJ file %s %p\n",gjProd.Data(),fileGJ[icase]);
    }
    
    if ( !fileJJ[icase] ) 
    {
      h2M02Iso[icase] = 0;
      h2M02Not[icase] = 0;
      
      continue;
    }
    
    //--------------------------------------------------------------------------
    // Get M02 vs E histograms
    //--------------------------------------------------------------------------
    
    h2M02Not[icase] = (TH2D*) fileJJ[icase]->Get("AnaIsolPhoton_hPtLambda0NoIso");
    if ( addGJ ) 
      h2M02Not[icase]->Add((TH2D*) fileGJ[icase]->Get("AnaIsolPhoton_hPtLambda0NoIso"),gjFactor);
    
    h2M02Iso[icase] = (TH2D*) fileJJ[icase]->Get("AnaIsolPhoton_hPtLambda0Iso");
    // P3 does not need GJ in photon region, only in pi0 region
    //if ( addGJ ) 
    //  h2M02Iso[icase]->Add((TH2D*) fileGJ[icase]->Get("AnaIsolPhoton_hPtLambda0Iso"),gjFactor);

    if ( addGJ ) 
      h2M02IsoGJ[icase] = (TH2D*) fileGJ[icase]->Get("AnaIsolPhoton_hPtLambda0Iso");
    
//    h2M02Not[icase] = (TH2D*) fileJJ[icase]->Get("AnaIsolPhoton_TM1_hPtLambda0NoIso");
//    if ( addGJ ) 
//      h2M02Not[icase]->Add((TH2D*) fileGJ[icase]->Get("AnaIsolPhoton_TM1_hPtLambda0NoIso"),gjFactor);
//    h2M02Iso[icase] = (TH2D*) fileJJ[icase]->Get("AnaIsolPhoton_TM1_hPtLambda0Iso");
//    // P3 does not need GJ in photon region, only in pi0 region
//    //if ( addGJ ) 
//    //  h2M02Iso[icase]->Add((TH2D*) fileGJ[icase]->Get("AnaIsolPhoton_TM1_hPtLambda0Iso"),gjFactor);
//       
//    if ( addGJ )
//      h2M02IsoGJ[icase] = (TH2D*) fileGJ[icase]->Get("AnaIsolPhoton_TM1_hPtLambda0Iso");
    
    //--------------------------------------------------------------------------
    // Project M02 vs E plot in signal region
    //--------------------------------------------------------------------------

    Int_t binm02LowCutMin = h2M02Iso[icase]->GetYaxis()->FindBin(m02LowCutMin);
    Int_t binm02LowCutMax = h2M02Iso[icase]->GetYaxis()->FindBin(m02LowCutMax)-1;    
    
    if ( debug )
    {
      printf("Signal region: M02 [%2.2f,%2.2f] - Bins [%d,%d] \n",
             m02LowCutMin,m02LowCutMax,binm02LowCutMin,binm02LowCutMax);
    }
    
    hM02LowNo[icase] = 
    (TH1D*) h2M02Not[icase]->ProjectionX(Form("LowM02_NoIso_%s_%s",titleProd.Data(),cuts.Data()),
                                         binm02LowCutMin,binm02LowCutMax);
   // if(hM02LowNo[icase]) 
   //   printf("Entries %e Integral %e\n",hM02LowNo[icase]->GetEntries(),hM02LowNo[icase]->Integral());
    
    hM02LowIs[icase] = 
    (TH1D*) h2M02Iso[icase]->ProjectionX(Form("LowM02_Iso_%s_%s",titleProd.Data(),cuts.Data()),
                                         binm02LowCutMin,binm02LowCutMax);
    
   // if(hM02LowIs[icase]) 
   //   printf("Entries %e Integral %e\n",hM02LowIs[icase]->GetEntries(),hM02LowIs[icase]->Integral());
    
    //--------------------------------------------------------------------------
    // Project M02 vs E plot in bkg region, 
    // if usePi0=kTRUE, recover the corresponding histogram
    //--------------------------------------------------------------------------

    if ( usePi0 )
    {
      hM02HigNo[icase]  = (TH1D*) fileJJ[icase]->Get("AnaIsolPi0SS_hPtNoIso");
      if(addGJ)hM02HigNo[icase]->Add((TH2D*) fileGJ[icase]->Get("AnaIsolPi0SS_hPtNoIso"),gjFactor);
      
      //if(hM02HigNo[icase]) 
      //  printf("Entries %e Integral %e\n",hM02HigNo[icase]->GetEntries(),hM02HigNo[icase]->Integral());
      
      hM02HigIs[icase]  = (TH1D*) fileJJ[icase]->Get("AnaIsolPi0SS_hE");
      if(addGJ)hM02HigIs[icase]->Add((TH2D*) fileGJ[icase]->Get("AnaIsolPi0SS_hE"),gjFactor);
    
     // if(hM02HigIs[icase]) 
      //  printf("Entries %e Integral %e\n",hM02HigIs[icase]->GetEntries(),hM02HigIs[icase]->Integral());
      
      //hM02HigNo[icase]  = (TH1D*) fileJJ[icase]->Get("AnaIsolPi0SS_TM1_hPtNoIso");
      //if(addGJ)hM02HigNo[icase]->Add((TH2D*) fileGJ[icase]->Get("AnaIsolPi0SS_TM1_hPtNoIso"),gjFactor);
      //hM02HigIs[icase]  = (TH1D*) fileJJ[icase]->Get("AnaIsolPi0SS_TM1_hE");
      //if(addGJ)hM02HigIs[icase]->Add((TH2D*) fileGJ[icase]->Get("AnaIsolPi0SS_TM1_hE"),gjFactor);
    }
    else
    {
      // Project M02 vs E plot
      
      Int_t binm02HigCutMin = h2M02Iso[icase]->GetYaxis()->FindBin(m02HigCutMin);
      Int_t binm02HigCutMax = h2M02Iso[icase]->GetYaxis()->FindBin(m02HigCutMax)-1;
      
      if ( debug )
      {
        printf("Bkg region: M02 [%2.2f,%2.2f] - Bins [%d,%d] \n",
               m02HigCutMin,m02HigCutMax,binm02HigCutMin,binm02HigCutMax);
      }
      
      hM02HigNo[icase] = 
      (TH1D*) h2M02Not[icase]->ProjectionX(Form("HigM02_NoIso_%s_%s",titleProd.Data(),cuts.Data()),
                                           binm02HigCutMin,binm02HigCutMax);
      //if(hM02HigNo[icase]) 
      //  printf("Entries %e Integral %e\n",hM02HigNo[icase]->GetEntries(),hM02HigNo[icase]->Integral());
      
      hM02HigIs[icase] = 
      (TH1D*) h2M02Iso[icase]->ProjectionX(Form("HigM02_Iso_%s_%s",titleProd.Data(),cuts.Data()),
                                           binm02HigCutMin,binm02HigCutMax);
      
      // P3 needs GJ on iso bkg photon region  not in signal
      if ( addGJ )
      {
        hM02HigIs[icase] ->Add 
        ( (TH1D*) h2M02IsoGJ[icase]->ProjectionX(Form("HigM02_GJ_Iso_%s_%s",titleProd.Data(),cuts.Data()),
                                                 binm02HigCutMin,binm02HigCutMax)
         ,gjFactor);
      }
      //if(hM02HigIs[icase]) 
      //  printf("Entries %e Integral %e\n",hM02HigIs[icase]->GetEntries(),hM02HigIs[icase]->Integral());
    }
    
    // Set histogram ranges and titles
    //
    hM02HigIs[icase]->SetMarkerStyle(marker[1]);
    hM02HigIs[icase]->SetMarkerColor(color[icase]);
    hM02HigIs[icase]->SetLineColor  (color[icase]);   
    hM02HigIs[icase]->SetAxisRange(emin,emax,"X");
    hM02HigIs[icase]->SetYTitle("d #sigma / d #it{p}_{T} (mb)");
    
    hM02HigNo[icase]->SetMarkerStyle(marker[3]);
    hM02HigNo[icase]->SetMarkerColor(color[icase]);
    hM02HigNo[icase]->SetLineColor  (color[icase]);   
    hM02HigNo[icase]->SetAxisRange(emin,emax,"X");
    hM02HigNo[icase]->SetYTitle("d #sigma / d #it{p}_{T} (mb)");

    hM02LowIs[icase]->SetMarkerStyle(marker[0]);
    hM02LowIs[icase]->SetMarkerColor(color[icase]);
    hM02LowIs[icase]->SetLineColor  (color[icase]);   
    hM02LowIs[icase]->SetAxisRange(emin,emax,"X");
    hM02LowIs[icase]->SetYTitle("d #sigma / d #it{p}_{T} (mb)");

    hM02LowNo[icase]->SetMarkerStyle(marker[2]);
    hM02LowNo[icase]->SetMarkerColor(color[icase]);
    hM02LowNo[icase]->SetLineColor  (color[icase]);
    hM02LowNo[icase]->SetAxisRange(emin,emax,"X");
    hM02LowNo[icase]->SetYTitle("d #sigma / d #it{p}_{T} (mb)");

    hM02LowNo[icase]->SetTitle(Form("C: %2.2f < #sigma_{long} < %2.2f - #Sigma #it{p}_{T} > 1 GeV/c",
                                m02LowCutMin,m02LowCutMax));
    hM02LowNo[icase]->SetTitleOffset(1.6,"Y");
    
    hM02LowIs[icase]->SetTitle(Form("A: %2.2f < #sigma_{long} < %2.2f - #Sigma #it{p}_{T} < 1 GeV/c",
                                    m02LowCutMin,m02LowCutMax));
    hM02LowIs[icase]->SetTitleOffset(1.6,"Y");
   
    if(usePi0)
    {
      hM02HigIs[icase]->SetTitle("B: Merged #pi^{0} - #Sigma #it{p}_{T} < 1 GeV/c");
      hM02HigIs[icase]->SetTitleOffset(1.6,"Y");
      
      hM02HigNo[icase]->SetTitle("D: Merged #pi^{0} - #Sigma #it{p}_{T} > 1 GeV/c");
      hM02HigNo[icase]->SetTitleOffset(1.6,"Y");
    }
    else
    {
      hM02HigNo[icase]->SetTitle(Form("D: %2.2f < #sigma_{long} < %2.2f - #Sigma #it{p}_{T} > 1 GeV/c",
                                      m02HigCutMin,m02HigCutMax));
      hM02HigNo[icase]->SetTitleOffset(1.6,"Y");
      
      hM02HigIs[icase]->SetTitle(Form("B: %2.2f < #sigma_{long} < %2.2f - #Sigma #it{p}_{T} < 1 GeV/c",
                                      m02HigCutMin,m02HigCutMax));
      hM02HigIs[icase]->SetTitleOffset(1.6,"Y");
    }
    
    if ( titleProd.Contains("LHC") )
    {
      hM02LowIs[icase]->Sumw2();
      hM02HigIs[icase]->Sumw2();
      hM02LowNo[icase]->Sumw2();
      hM02HigNo[icase]->Sumw2();
    }
    
    // Rebin
    if(rebin > 0)
    {
      hM02LowIs[icase]->Rebin(rebin);
      hM02HigIs[icase]->Rebin(rebin);
      hM02LowNo[icase]->Rebin(rebin);
      hM02HigNo[icase]->Rebin(rebin);
    }
    else
    {
      //printf("Rebin with n bins %d\n",nBins);
      hM02LowIs[icase] = (TH1D*) hM02LowIs[icase]->Rebin(nBins,Form("%s_Rebin",hM02LowIs[icase]->GetName()),binE);
      hM02HigIs[icase] = (TH1D*) hM02HigIs[icase]->Rebin(nBins,Form("%s_Rebin",hM02HigIs[icase]->GetName()),binE);
      hM02LowNo[icase] = (TH1D*) hM02LowNo[icase]->Rebin(nBins,Form("%s_Rebin",hM02LowNo[icase]->GetName()),binE);
      hM02HigNo[icase] = (TH1D*) hM02HigNo[icase]->Rebin(nBins,Form("%s_Rebin",hM02HigNo[icase]->GetName()),binE);
      
      if ( debug )
      {
        printf("Rebinned ABCD names: %s, %s,\n %s, %s\n",
               hM02LowIs[icase]->GetName(),hM02HigIs[icase]->GetName(),
               hM02LowNo[icase]->GetName(),hM02HigNo[icase]->GetName());
      }
      
      // From PlotUtils.C
      ScaleBinBySize(hM02LowIs[icase]);
      ScaleBinBySize(hM02LowNo[icase]);
      ScaleBinBySize(hM02HigIs[icase]);
      ScaleBinBySize(hM02HigNo[icase]);
    }
  }
  
  // 
  // Plot
  //
  
  if ( debug )
    printf("PLOT Double Ratios \n");

  TString fileName ;
  
  gStyle->SetPadRightMargin(0.01);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadTopMargin(0.08);

  gStyle->SetTitleFontSize(0.06);
  //  gStyle->SetStatFontSize(0.06);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  
  /// ABCD ///
  TCanvas * cABCD = new TCanvas(Form("cABCD_%s_%s" ,cuts.Data(),prodTitleOutput.Data()),
                                Form("ABCD, %s, %s",cuts.Data(),prodTitleOutput.Data()),1000,1000);
  
  cABCD->Divide(2,2);
  
  TLegend *lABCD = new TLegend(0.3,0.6,0.89,0.89);
  lABCD->SetFillColor(0);
  lABCD->SetFillStyle(0);
  lABCD->SetLineColor(0);
  lABCD->SetBorderSize(0);
  lABCD->SetTextSize(0.05);
  
  /////// A ////////////
  cABCD->cd(1);
      
  gPad->SetLogy();
  
  for(Int_t icase = 0; icase < nMCCases; icase++)
  {
    //printf("icase %d %p\n",icase,hM02LowIs[icase]);
    if(!hM02LowIs[icase]) continue; 
      
    if ( icase==0 ) hM02LowIs[icase]->Draw("H");
    else            hM02LowIs[icase]->Draw("H same");
    
//    if(hM02LowIs[icase]->GetMaximum() > hM02LowIs[0]->GetMaximum())
//      hM02LowIs[0]->SetMaximum(hM02LowIs[icase]->GetMaximum()*1.2);
    
    lABCD->AddEntry(hM02LowIs[icase],Form("%s",mcLeg[icase].Data()),"L");
  }
  
  if ( nCases > 1 )
    lABCD->Draw();
  
  /////// B ////////////

  cABCD->cd(2);
  
  gPad->SetLogy();
  
  for(Int_t icase = 0; icase < nMCCases; icase++)
  {
    if(!hM02HigIs[icase]) continue; 
    
    if ( icase==0 ) hM02HigIs[icase]->Draw("H");
    else            hM02HigIs[icase]->Draw("H same");
 
  }
  
  if ( nCases > 1 )
    lABCD->Draw();
  
  /////// C ////////////
  cABCD->cd(3);
  
  gPad->SetLogy();
  
  for(Int_t icase = 0; icase < nMCCases; icase++)
  {
    if(!hM02LowNo[icase]) continue; 
    
    if ( icase==0 ) hM02LowNo[icase]->Draw("H");
    else            hM02LowNo[icase]->Draw("H same");
    
    //    if(hM02LowNo[icase]->GetMaximum() > hM02LowNo[0]->GetMaximum())
    //      hM02LowNo[0]->SetMaximum(hM02LowNo[icase]->GetMaximum()*1.2);
  }
  
  if ( nCases > 1 )
    lABCD->Draw();
  
  /////// D ////////////
  
  cABCD->cd(4);
  
  gPad->SetLogy();
  
  for(Int_t icase = 0; icase < nMCCases; icase++)
  {
    if(!hM02HigNo[icase]) continue; 
    
    if ( icase==0 ) hM02HigNo[icase]->Draw("H");
    else            hM02HigNo[icase]->Draw("H same");
    
    //    if(hM02HigNo[icase]->GetMaximum() > hM02HigNo[0]->GetMaximum())
    //      hM02HigNo[0]->SetMaximum(hM02HigNo[icase]->GetMaximum()*1.2);
  }

  if ( nCases > 1 )
    lABCD->Draw();
  
  fileName = Form("%s/figures/ABCD_%s_pT_%s",opt.Data(),cuts.Data(),prodTitleOutput.Data());
  fileName+=".eps";
  cABCD->Print(fileName);

  /// ABCD Ratio of the different MC productions///
  if ( nMCCases > 1 )
  {
    TCanvas * cABCDRatio = new TCanvas(Form("cABCDRatio_%s_%s"        ,cuts.Data(),prodTitleOutput.Data()),
                                       Form("ABCD cases Ratio, %s, %s",cuts.Data(),prodTitleOutput.Data()),
                                       1000,1000);
    
    cABCDRatio->Divide(2,2);
    
    TLegend *lABCDRat = new TLegend(0.3,0.6,0.89,0.89);
    lABCDRat->SetFillColor(0);
    lABCDRat->SetFillStyle(0);
    lABCDRat->SetLineColor(0);
    lABCDRat->SetBorderSize(0);
    lABCDRat->SetTextSize(0.05);
    //lABCDRat->SetHeader(Form("Den.: %s",mcLeg[nMCCases-1].Data()));
    lABCDRat->SetHeader(Form("Den.: %s",mcLeg[0].Data()));
    
    /////// A ////////////
    cABCDRatio->cd(1);
    gPad->SetGridy();
    
    //for(Int_t icase = 0; icase < nMCCases-1; icase++)
    for(Int_t icase = 1; icase < nMCCases; icase++)
    {
      if(!hM02LowIs[icase]) continue; 
      
      TString cloneName = Form("M02LowIsCasesRat%s_%s_%s",mcCases[icase].Data(),mcCasesPath[icase].Data(),cuts.Data());
      // Careful, in case of data, take the data name and not the MC case name
      if ( titleProd.Contains("LHC") ) cloneName = Form("M02LowIsCasesRat%s_%s_%s",titleProd.Data(),"",cuts.Data());
      
      //printf(">>>>>>>>>>>>>> Clone name %s\n",cloneName.Data());
      
      hM02LowIsRat[icase] = (TH1D*) hM02LowIs[icase]->Clone(cloneName);
      //hM02LowIsRat[icase]->Divide(hM02LowIs[nMCCases-1]);
      hM02LowIsRat[icase]->Divide(hM02LowIs[0]);
      
      hM02LowIsRat[icase]->SetYTitle("Ratio");
      hM02LowIsRat[icase]->SetMaximum(2.5);
      hM02LowIsRat[icase]->SetMinimum(0.5);
      
      if(icase == 0) hM02LowIsRat[icase]->Draw("H");
      else           hM02LowIsRat[icase]->Draw("H same");
      
      lABCDRat->AddEntry(hM02LowIsRat[icase],Form("Num.:%s",mcLeg[icase].Data()),"L");
    }
    
    lABCDRat->Draw();
    
    /////// B ////////////
    cABCDRatio->cd(2);
    gPad->SetGridy();
    
    for(Int_t icase = 1; icase < nMCCases; icase++)
    //for(Int_t icase = 0; icase < nMCCases-1; icase++)
    {
      if(!hM02HigIs[icase]) continue; 
      
      TString cloneName = Form("M02HighIsCasesRat%s_%s_%s",mcCases[icase].Data(),mcCasesPath[icase].Data(),cuts.Data());
      // Careful, in case of data, take the data name and not the MC case name
      if ( titleProd.Contains("LHC") ) cloneName = Form("M02HighIsCasesRat%s_%s_%s",titleProd.Data(),"",cuts.Data());
      
      hM02HigIsRat[icase] = (TH1D*) hM02HigIs[icase]->Clone(cloneName);
      //hM02HigIsRat[icase]->Divide(hM02HigIs[nMCCases-1]);
      hM02HigIsRat[icase]->Divide(hM02HigIs[0]);
      
      hM02HigIsRat[icase]->SetYTitle("Ratio");
      hM02HigIsRat[icase]->SetMaximum(2.5);
      hM02HigIsRat[icase]->SetMinimum(0.5);
      
      if(icase == 0) hM02HigIsRat[icase]->Draw("H");
      else           hM02HigIsRat[icase]->Draw("H same");
    }
    
    lABCDRat->Draw();
    
    /////// C ////////////
    cABCDRatio->cd(3);
    gPad->SetGridy();
    
    for(Int_t icase = 1; icase < nMCCases; icase++)
    //for(Int_t icase = 0; icase < nMCCases-1; icase++)
    {
      if(!hM02LowNo[icase]) continue; 
      
      TString cloneName = Form("M02LowNoCasesRat%s_%s_%s",mcCases[icase].Data(),mcCasesPath[icase].Data(),cuts.Data());
      // Careful, in case of data, take the data name and not the MC case name
      if ( titleProd.Contains("LHC") ) cloneName = Form("M02LowNoCasesRat%s_%s_%s",titleProd.Data(),"",cuts.Data());
      
      hM02LowNoRat[icase] = (TH1D*) hM02LowNo[icase]->Clone(cloneName);
      hM02LowNoRat[icase]->Divide(hM02LowNo[0]);
      //hM02LowNoRat[icase]->Divide(hM02LowNo[nMCCases-1]);
      
      hM02LowNoRat[icase]->SetYTitle("Ratio");
      hM02LowNoRat[icase]->SetMaximum(2.5);
      hM02LowNoRat[icase]->SetMinimum(0.5);
      
      if(icase == 0) hM02LowNoRat[icase]->Draw("H");
      else           hM02LowNoRat[icase]->Draw("H same");
    }
    
    lABCDRat->Draw();
    
    /////// D ////////////
    cABCDRatio->cd(4);
    gPad->SetGridy();
    
    for(Int_t icase = 1; icase < nMCCases; icase++)
    //for(Int_t icase = 0; icase < nMCCases-1; icase++)
    {
      if(!hM02HigNo[icase]) continue; 
      
      TString cloneName = Form("M02HighNoCasesRat%s_%s_%s",mcCases[icase].Data(),mcCasesPath[icase].Data(),cuts.Data());
      // Careful, in case of data, take the data name and not the MC case name
      if ( titleProd.Contains("LHC") ) cloneName = Form("M02HighNoCasesRat%s_%s_%s",titleProd.Data(),"",cuts.Data());
      
      hM02HigNoRat[icase] = (TH1D*) hM02HigNo[icase]->Clone(cloneName);
      //hM02HigNoRat[icase]->Divide(hM02HigNo[nMCCases-1]);
      hM02HigNoRat[icase]->Divide(hM02HigNo[0]);
      
      hM02HigNoRat[icase]->SetYTitle("Ratio");
      hM02HigNoRat[icase]->SetMaximum(2.5);
      hM02HigNoRat[icase]->SetMinimum(0.5);
      
      if(icase == 0) hM02HigNoRat[icase]->Draw("H");
      else           hM02HigNoRat[icase]->Draw("H same");
    }
    
    lABCDRat->Draw();
    
    
    fileName = Form("%s/figures/ABCD_%s_pT_CasesRatio_%s",opt.Data(),cuts.Data(),prodTitleOutput.Data());
    fileName+=format;
    cABCDRatio->Print(fileName);
  }
  
  /// ABCD Double Ratio ///
  TCanvas * cABCD2Ratio = new TCanvas(Form("cABCD2Ratio_%s_%s"        ,cuts.Data(),prodTitleOutput.Data()),
                                      Form("ABCD Double Ratio, %s, %s",cuts.Data(),prodTitleOutput.Data()),
                                      1000,1000);
  
  TLegend *lABCD2Rat = new TLegend(0.15,0.6,0.5,0.89);
  lABCD2Rat->SetFillColor(0);
  lABCD2Rat->SetFillStyle(0);
  lABCD2Rat->SetLineColor(0);
  lABCD2Rat->SetBorderSize(0);
  lABCD2Rat->SetTextSize(0.05);
  
  gPad->SetGridy();
  
  for(Int_t icase = 0; icase < nMCCases; icase++)
  {
    if(!hM02LowIs[icase]) continue; 
    
    TString cloneName = Form("DoubleRat%s_%s_%s",mcCases[icase].Data(),mcCasesPath[icase].Data(),cuts.Data());
    // Careful, in case of data, take the data name and not the MC case name
    if ( titleProd.Contains("LHC") ) cloneName = Form("DoubleRat%s_%s",titleProd.Data(),cuts.Data());
    
    //printf(">>>>> CloneName %s for %s\n",cloneName.Data(),titleProd.Data());
    
    hDoubleRatio[icase] = (TH1D*) hM02LowIs[icase]->Clone(cloneName);
    hDoubleRatio[icase]->Divide  (hM02LowNo[icase]);
    hDoubleRatio[icase]->Divide  (hM02HigIs[icase]);
    hDoubleRatio[icase]->Multiply(hM02HigNo[icase]);
    
    hDoubleRatio[icase]->SetYTitle("(A*D)/(C*B)");
    hDoubleRatio[icase]->SetTitle(Form("#Sigma #it{p}_{T} > 1 GeV/c - AC: %2.2f < #sigma_{long} < %2.2f - BD:  %2.2f < #sigma_{long} < %2.2f",
                                       m02LowCutMin,m02LowCutMax,m02HigCutMin,m02HigCutMax));
    
    if(usePi0) hDoubleRatio[icase]->SetTitle(Form("#Sigma #it{p}_{T} > 1 GeV/c - AC: %2.2f < #sigma_{long} < %2.2f - BD:Merged #pi^{0} band",
                                                  m02LowCutMin,m02LowCutMax));
    
    hDoubleRatio[icase]->SetMaximum(1.4);
    hDoubleRatio[icase]->SetMinimum(0.7);
    
    if(titleProd.Contains("LHC"))
    {
      hDoubleRatio[icase]->SetMaximum(8);
      hDoubleRatio[icase]->SetMinimum(0.8);      
    }
    
    if(icase == 0) hDoubleRatio[icase]->Draw("H");
    else           hDoubleRatio[icase]->Draw("H same");
    
    lABCD2Rat->AddEntry(hDoubleRatio[icase],Form("%s",mcLeg[icase].Data()),"L");
  }
  
  lABCD2Rat->Draw();
  
  
  fileName = Form("%s/figures/ABCD_%s_DoubleRatio_%s",opt.Data(),cuts.Data(),prodTitleOutput.Data());
  fileName+=format;
  cABCD2Ratio->Print(fileName);
  
  if(nMCCases > 1)
  {
    /// ABCD Double Ratio Cases ///
    TCanvas * cABCD2RatioCases = new TCanvas(Form("cABCD2RatioCases_%s_%s"         ,cuts.Data(),prodTitleOutput.Data()),
                                             Form("ABCD Double Ratio Cases, %s, %s",cuts.Data(),prodTitleOutput.Data()),
                                             1000,1000);
    
    TLegend *lABCD2RatCa = new TLegend(0.15,0.6,0.5,0.89);
    lABCD2RatCa->SetFillColor(0);
    lABCD2RatCa->SetFillStyle(0);
    lABCD2RatCa->SetLineColor(0);
    lABCD2RatCa->SetBorderSize(0);
    lABCD2RatCa->SetTextSize(0.05);
    //lABCD2RatCa->SetHeader(Form("Den.: %s",mcLeg[nMCCases-1].Data()));
    lABCD2RatCa->SetHeader(Form("Den.: %s",mcLeg[0].Data()));
    
    gPad->SetGridy();
    
    //for(Int_t icase = 0; icase < nMCCases-1; icase++)
    for(Int_t icase = 1; icase < nMCCases; icase++)
    {
      if(!hM02LowIs[icase]) continue; 
      
      TString cloneName = Form("DoubleRatCases%s_%s_%s",mcCases[icase].Data(),mcCasesPath[icase].Data(),cuts.Data());
      // Careful, in case of data, take the data name and not the MC case name
      if ( titleProd.Contains("LHC") ) cloneName = Form("DoubleRatCases%s_%s_%s",titleProd.Data(),"",cuts.Data());
      
      hDoubleRatioCases[icase] = (TH1D*) hDoubleRatio[icase]->Clone(cloneName);
      //hDoubleRatioCases[icase]->Divide(hDoubleRatio[nMCCases-1]);
      hDoubleRatioCases[icase]->Divide(hDoubleRatio[0]);
      
      hDoubleRatioCases[icase]->SetTitle(Form("#Sigma #it{p}_{T} > 1 GeV/c - AC: %2.2f < #sigma_{long} < %2.2f - BD:  %2.2f < #sigma_{long} < %2.2f",
                                              m02LowCutMin,m02LowCutMax,m02HigCutMin,m02HigCutMax));    
      
      if(usePi0) hDoubleRatioCases[icase]->SetTitle(Form("#Sigma #it{p}_{T} > 1 GeV/c - AC: %2.2f < #sigma_{long} < %2.2f - BD: Merged #pi^{0} band",
                                                    m02LowCutMin,m02LowCutMax));
      hDoubleRatioCases[icase]->SetYTitle("Ratio");
      hDoubleRatioCases[icase]->SetMaximum(1.4);
      hDoubleRatioCases[icase]->SetMinimum(0.7);
      
      if(icase == 0) hDoubleRatioCases[icase]->Draw("");
      else           hDoubleRatioCases[icase]->Draw("same");
      
      lABCD2RatCa->AddEntry(hDoubleRatioCases[icase],Form("Num.: %s",mcLeg[icase].Data()),"L");
    }
    
    lABCD2RatCa->Draw();
    
    
    fileName = Form("%s/figures/ABCD_%s_DoubleRatio_CaseRatio_%s",opt.Data(),cuts.Data(),prodTitleOutput.Data());
    fileName+=format;
    cABCD2RatioCases->Print(fileName);
  }
  
  //-----------------------------------------------
  // Save the ratios in external file
  //-----------------------------------------------
  fileName = Form("%s/figures/ABCD_%s_DoubleRatios_%s.root",opt.Data(),cuts.Data(),prodTitleOutput.Data());

  TFile * fout = new TFile(fileName,"recreate");

  for(Int_t icase = 0; icase < nMCCases; icase++)
  {
    //printf("icase %d\n",icase);
    if(hM02LowIs[icase])hM02LowIs [icase]->Write();
    if(hM02LowNo[icase])hM02LowNo [icase]->Write();  
    if(hM02HigIs[icase])hM02HigIs [icase]->Write();
    if(hM02HigNo[icase])hM02HigNo [icase]->Write();
   
    if(hM02LowIsRat[icase])hM02LowIsRat[icase]->Write();
    if(hM02LowNoRat[icase])hM02LowNoRat[icase]->Write();  
    if(hM02HigIsRat[icase])hM02HigIsRat[icase]->Write();
    if(hM02HigNoRat[icase])hM02HigNoRat[icase]->Write();
    
    if(hDoubleRatio     [icase])hDoubleRatio     [icase]->Write();
    if(hDoubleRatioCases[icase])hDoubleRatioCases[icase]->Write();
  }

  fout->Close();
}

//------------------------------------------------------------------------------
/// Get the histograms defining the different purity regions ABDC, 
/// and do the double ratios needed for purity calculation in MakePurity()
/// where data and MC double ratios are combined. This method is used the same
/// for that and MC.
///
/// \param dataName    : Name analyzed data sample
/// \param simuName    : Name analyzed MC sample
/// \param nCases      : in case of MC, the different cases studied, in case of data it is 1
/// \param mcCasesPath : path to MC/Data case input file
/// \param mcCases     : simplified name of MC/Data sample
/// \param mcLef       : name of MC/Data sample for legends in plots
/// \param gj          : string to select the addition and scaling of the GJ sample in MC - "","_GJx0.50","_GJx1.00","_GJx2.00"
/// \param usePi0      : Use as bkg region the pi0 shower shape band
/// \param m02LowCutMin: Signal region, minimum shower shape value
/// \param m02LowCutMax: Signal region, maximum shower shape value
/// \param m02HighCutMin: Bkg region, minimum shower shape value
/// \param m02HighCutMax: Bkg region, maximum shower shape value
/// \param debug       : Bool to activate printfs.
//------------------------------------------------------------------------------
void MakePurity
(
 TString dataName, TString simuName
 , Int_t nCases, TString * mcCasesPath, TString * mcCases, TString * mcLeg
 , TString gj       = ""//"_GJx0.50","_GJx1.00","_GJx2.00"
 , Bool_t  usePi0    = kTRUE
 , Float_t m02LowCutMin = 0.10
 , Float_t m02LowCutMax = 0.27
 , Float_t m02HigCutMin = 0.40
 , Float_t m02HigCutMax = 3.00
 , Bool_t  debug        = kFALSE
 )
{  
  // Init histogram arrays and other things
  // 
  const Int_t nMCCases = nCases;
  
  TH1D* hDoubleRatioData = 0;
  TH1D* hDoubleRatioSimu [nMCCases];
  TH1D* hPurityCases     [nMCCases];  
  TH1D* hPurityCasesRatio[nMCCases];  
  
  for(Int_t icase = 0; icase < nMCCases; icase++)
  {
    hDoubleRatioSimu [icase]=0;
    hPurityCases     [icase]=0;
    hPurityCasesRatio[icase]=0;
  }
  
  TString cuts   = "";
  TString merged = "";
  if ( usePi0 )
  {
    merged = "_MergedPi0";
    cuts = Form("Low_%2.2f-%2.2f"                 ,m02LowCutMin,m02LowCutMax);
  }
  else 
    cuts = Form("Low_%2.2f-%2.2f_High_%2.2f-%2.2f",m02LowCutMin,m02LowCutMax,m02HigCutMin,m02HigCutMax);
  
  // Histogram limits
  Float_t emin = 8;
  Float_t emax = 29;
  if ( simuName.Contains("High") )
  {
    emin = 16;
    emax = 45;
  }
  
  // ----------
  // Open files
  // ----------
  TFile*  fileData = TFile::Open(Form("%s/figures/ABCD_%s_DoubleRatios_%s%s.root"  ,opt.Data(),cuts.Data(),dataName.Data(),merged.Data()));
  TFile*  fileSimu = TFile::Open(Form("%s/figures/ABCD_%s_DoubleRatios_%s%s%s.root",opt.Data(),cuts.Data(),simuName.Data(),gj.Data(),merged.Data()));
  
  if ( debug )
  {
    printf("Get files (%p, %p) and histograms \n",fileData,fileSimu);
    printf("\t %s/figures/ABCD_%s_DoubleRatios_%s%s.root\n"  ,opt.Data(),cuts.Data(),dataName.Data(),merged.Data());
    printf("\t %s/figures/ABCD_%s_DoubleRatios_%s%s%s.root\n",opt.Data(),cuts.Data(),simuName.Data(),gj.Data(),merged.Data());
  }
  
  if ( !fileData || !fileSimu ) return;
    
  // Data
  hDoubleRatioData = (TH1D*) fileData->Get(Form("DoubleRat%s_%s",dataName.Data(),cuts.Data())); 
  
  if ( debug )
  {
    printf("\t data histo %p\n",hDoubleRatioData);
    printf("\t \t DoubleRat%s_\n",dataName.Data());
  }                                              
  
  if ( !hDoubleRatioData ) return;
  
  // Simu
  for(Int_t icase = 0; icase < nMCCases; icase++)
  {
    hDoubleRatioSimu[icase] = (TH1D*) fileSimu->Get(Form("DoubleRat%s_%s_%s",mcCases[icase].Data(),mcCasesPath[icase].Data(),cuts.Data()));
    if ( debug )
    {
      printf("\t simu case %d histo %p\n",icase, hDoubleRatioSimu[icase]);
      printf("\t \t DoubleRat%s_%s\n",mcCases[icase].Data(),mcCasesPath[icase].Data());
    }
    
    if ( !hDoubleRatioSimu[icase] ) return;
    
    hPurityCases[icase] = (TH1D*) hDoubleRatioSimu[icase]->Clone(Form("Purity_Data_%s_MC_%s_%s_%s%s",
                                                                      dataName.Data(),mcCases[icase].Data(),mcCasesPath[icase].Data(),cuts.Data(),gj.Data()));
    hPurityCases[icase]->Divide(hPurityCases[icase],hDoubleRatioData,1,1,"");
    
    for(Int_t ibin = 1; ibin < hPurityCases[icase]->GetNbinsX(); ibin++)
    { 
      //      Double_t content = 1.-hDoubleRatioSimu[icase]->GetBinContent(ibin)/hDoubleRatioData->GetBinContent(ibin);
      //      Double_t error = content* TMath::Sqrt( hDoubleRatioSimu[icase]->GetBinError  (ibin)*hDoubleRatioSimu[icase]->GetBinError  (ibin)/
      //                                            (hDoubleRatioSimu[icase]->GetBinContent(ibin)*hDoubleRatioSimu[icase]->GetBinContent(ibin)) +
      //                                             hDoubleRatioData       ->GetBinError  (ibin)*hDoubleRatioData       ->GetBinError  (ibin)/
      //                                            (hDoubleRatioData       ->GetBinContent(ibin)*hDoubleRatioData       ->GetBinContent(ibin)));
      
      Double_t content2 = 1.-hPurityCases[icase]->GetBinContent(ibin);
      Double_t error2   = hPurityCases[icase]->GetBinError(ibin);
      //      printf("ibin %d purity = 1-%2.2f(+-%2.2f)*%2.2f(+-%2.2f) = %2.2f +- %2.2f; p2 = %2.2f +- %2.2f \n",
      //             ibin,hDoubleRatioSimu[icase]->GetBinContent(ibin),hDoubleRatioSimu[icase]->GetBinError(ibin),
      //                  hDoubleRatioData       ->GetBinContent(ibin),hDoubleRatioData       ->GetBinError(ibin),
      //             content,error,content2,error2);
      
      hPurityCases[icase]->SetBinContent(ibin,content2);
      hPurityCases[icase]->SetBinError  (ibin,error2  );
    }
    
    hPurityCases[icase]->SetMarkerStyle(marker[1]);
    hPurityCases[icase]->SetMarkerColor(color[icase]);
    hPurityCases[icase]->SetLineColor  (color[icase]);   
    hPurityCases[icase]->SetAxisRange(emin,emax,"X");
    hPurityCases[icase]->SetYTitle("Purity");
  } // simu case loop
  
  if ( debug )
    printf("** Make purity plots ** \n");
  
  // Ratio to last case
  for(Int_t icase = 1; icase < nMCCases; icase++)
    //for(Int_t icase = 0; icase < nMCCases-1; icase++)
  {
    hPurityCasesRatio[icase] = (TH1D*) hPurityCases[icase]->Clone(Form("Purity_%s_%s_%s%s_Ratio",
                                                                       mcCases[icase].Data(),mcCasesPath[icase].Data(),cuts.Data(),gj.Data())); 
    hPurityCasesRatio[icase]->Divide(hPurityCases[0] );
    //hPurityCasesRatio[icase]->Divide(hPurityCases[nMCCases-1] );
    //hPurityCases[icase]->SetYTitle(Form("Purity X / Purity %s",mcLeg[nMCCases-1].Data()));
    hPurityCases[icase]->SetYTitle(Form("Purity X / Purity %s",mcLeg[0].Data()));
  }
  
  gStyle->SetPadRightMargin(0.01);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadTopMargin(0.08);
  
  gStyle->SetTitleFontSize(0.06);
  //  gStyle->SetStatFontSize(0.06);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  
  /// Purity ///
  TCanvas * cPurity = new TCanvas(Form("cPurity_%s_Data_%s_MC_%s%s%s",cuts.Data(),dataName.Data(),simuName.Data(),gj.Data(),merged.Data()),
                                  Form("Purity: %s_Data_%s_MC_%s%s%s",cuts.Data(),dataName.Data(),simuName.Data(),gj.Data(),merged.Data())
                                  ,1000,1000);
  
  TLegend *lPu = new TLegend(0.3,0.15,0.98,0.4);
  lPu->SetFillColor(0);
  lPu->SetFillStyle(0);
  lPu->SetLineColor(0);
  lPu->SetBorderSize(0);
  lPu->SetTextSize(0.05);
  
  gPad->SetGridy();
  
  for(Int_t icase = 0; icase < nMCCases; icase++)
  {
    hPurityCases[icase]->SetTitle(Form("#Sigma #it{p}_{T} > 1 GeV/c - AC: %2.2f < #sigma_{long} < %2.2f - BD:  %2.2f < #sigma_{long} < %2.2f ",
                                       m02LowCutMin,m02LowCutMax,m02HigCutMin,m02HigCutMax));
    
    if(usePi0) hPurityCases[icase]->SetTitle(Form("#Sigma #it{p}_{T} > 1 GeV/c - AC: %2.2f < #sigma_{long} < %2.2f - BD: Merged #pi^{0} band",
                                                  m02LowCutMin,m02LowCutMax));  
    
    hPurityCases[icase]->SetMaximum(0.85);
    hPurityCases[icase]->SetMinimum(0.0);
    
    if(icase == 0) hPurityCases[icase]->Draw("H");
    else           hPurityCases[icase]->Draw("H same");
    
    lPu->AddEntry(hPurityCases[icase],Form("%s",mcLeg[icase].Data()),"L");
  }
  
  lPu->Draw();
  
  
  TString fileName = Form("%s/figures/Purity_%s_Data_%s_MC_%s%s%s",opt.Data(),cuts.Data(),dataName.Data(),simuName.Data(),gj.Data(),merged.Data());
  
  fileName+=format;
  cPurity->Print(fileName);
  
  if(nMCCases > 1)
  {
    /// Purity Ratio Cases ///
    TCanvas * cPurityRat = new TCanvas(Form("cPurityRat_%s_Data_%s_MC_%s%s%s",cuts.Data(),dataName.Data(),simuName.Data(),gj.Data(),merged.Data()),
                                       Form("Purity Ratio Cases: %s_Data_%s_MC_%s%s%s",cuts.Data(),dataName.Data(),simuName.Data(),gj.Data(),merged.Data()),
                                       1000,1000);
    
    TLegend *lPuCa = new TLegend(0.3,0.7,0.98,0.89);
    lPuCa->SetFillColor(0);
    lPuCa->SetFillStyle(0);
    lPuCa->SetLineColor(0);
    lPuCa->SetBorderSize(0);
    lPuCa->SetTextSize(0.04);
    //lPuCa->SetHeader(Form("Den.: %s",mcLeg[nMCCases-1].Data()));
    lPuCa->SetHeader(Form("Den.: %s",mcLeg[0].Data()));
    
    gPad->SetGridy();
    
    for(Int_t icase = 1; icase < nMCCases; icase++)
      //for(Int_t icase = 0; icase < nMCCases-1; icase++)
    {
      hPurityCasesRatio[icase]->SetTitle(Form("#Sigma #it{p}_{T} > 1 GeV/c - AC: %2.2f < #sigma_{long} < %2.2f - BD:  %2.2f < #sigma_{long} < %2.2f",
                                              m02LowCutMin,m02LowCutMax,m02HigCutMin,m02HigCutMax));  
      if(usePi0) hPurityCasesRatio[icase]->SetTitle(Form("#Sigma #it{p}_{T} > 1 GeV/c - AC: %2.2f < #sigma_{long} < %2.2f - BD: Merged #pi^{0} band",
                                                         m02LowCutMin,m02LowCutMax));  
      
      hPurityCasesRatio[icase]->SetYTitle("Ratio");
      
      hPurityCasesRatio[icase]->SetMaximum(1.5);
      hPurityCasesRatio[icase]->SetMinimum(0.7);
      
      if(icase == 0) hPurityCasesRatio[icase]->Draw("");
      else           hPurityCasesRatio[icase]->Draw("same");
      
      lPuCa->AddEntry(hPurityCasesRatio[icase],Form("Num.: %s",mcLeg[icase].Data()),"L");
    }
    
    lPuCa->Draw();
    
    
    fileName = Form("%s/figures/Purity_%s_CaseRatio_Data_%s_MC_%s%s%s",opt.Data(),cuts.Data(),dataName.Data(),simuName.Data(),gj.Data(),merged.Data());
    fileName+=format;
    cPurityRat->Print(fileName);
  }
  
}

//------------------------------------------------------------------------------
/// Main steering method.
/// Execute here the different selection options
///
/// \param doDoRat   : Execute method MakeDoubleRatios()
/// \param doPurity  : Execute method MakePurity()
/// \param debug     : Bool to activate printfs.
//------------------------------------------------------------------------------
void MakePurityCalculationAndComparisons
(
 Bool_t  doDoRat  = kTRUE,
 Bool_t  doPurity = kTRUE,
 Bool_t  debug    = kFALSE
)
{
  // Declare the input histograms, output of the analysis
  // First the data analysis, after all the MC analysis from different jet-jet productions
  // Do not merge the Gamma-jet to the jet-jet production, only internally.
  //
  Int_t nData = 3;
  TString titleDataMC[] = {"LHC11cd_EMC7","pp_7TeV_Dec_Low","pp_7TeV_Dec_High"};
  TString gjProd        = "pp_7TeV_GJ";

  // In case gamma-jet is added internally, define the extra scaling factors for comparison 
  // in a later stage.
  Int_t nfactors = 3;
  Float_t gjFactor       [] = {        1,         2,         0.5  }; // Used MakeDoubleRatio
  TString gjFactorString [] = {"","_GJx1.00","_GJx2.00","_GJx0.50"}; // Used in MakePurity
  
  // Declare here the MC analysis with different analysis cuts/emulation/settings
  // and their corresponding title for the legends and output files
  //
  Int_t nCases = 2;
  
  TString mcCasesPath[] = 
  {
      "mimic0_Scaled2_v3"
    , "mimic10c_EcellCut_Scaled2_v3"  
  };
  
  TString mcCases[] = 
  {
      "Default"
    , "TCardEmulation"
  };
  
  TString mcLeg[] = 
  {
      "MC default"
    , "MC mimic"
  };
  
  // Declare here the shoer shape cuts for the ABCD regions
  //
  Int_t nShShCutsMaxLow  = 3;
  Int_t nShShCutsMinHigh = 3;
  Float_t m02LowCutMin = 0.10;
  Float_t m02LowCutMax [] = {0.27, 0.30, 0.40};
  Float_t m02HigCutMin [] = {0.30, 0.40, 0.50};
  Float_t m02HigCutMax = 2.00;

  if ( debug ) 
    printf("N data %d, cases %d nShShCuts MaxLow %d MinHigh %d\n",
           nData, nCases, nShShCutsMaxLow, nShShCutsMinHigh);
  
  Int_t   rebin    = 0; // Use defined E binning
  Int_t npi0       = 2;
  Int_t ngj        = 2;
  
//  nfactors         = 1;
//  nShShCutsMaxLow  = 1;
//  nShShCutsMinHigh = 1;
  //nCases           = 1;

  // Calculate and plot double ratios 
  //
  if ( doDoRat )
  {     
    if ( debug )
      printf("=== Make double ratios === \n");
    
    for(Int_t igj = 0; igj < ngj; igj++)
    {
      for(Int_t ifactor = 0; ifactor < nfactors; ifactor++)
      {
        if ( igj == 0 && ifactor > 0 ) continue;
        
        for(Int_t idata = 0; idata < nData; idata++)
        {
          if( idata == 0 && igj > 0 ) continue;
          
          for(Int_t ish = 0; ish < nShShCutsMaxLow; ish++)
          {
            for(Int_t jsh = 0; jsh < nShShCutsMinHigh; jsh++)
            {
              if ( m02LowCutMax[ish] > m02HigCutMin[jsh] ) continue;
              
              for(Int_t ipi0 = 0; ipi0 < npi0; ipi0++)
              {
                if ( ipi0 && (ish > 1 || jsh > 0)  ) 
                  continue;
                
                if ( debug )
                  printf("\t idata %d %s, gj %s, usepi0 %d, addGJ %d, GJ factor %d-%2.1f, ish %d %2.2f, jsh %d %2.2f\n", 
                         idata, titleDataMC[idata].Data(), gjProd.Data(),
                         ipi0,igj,ifactor, gjFactor[ifactor],
                         ish, m02LowCutMax[ish],
                         jsh, m02HigCutMin[jsh]);
                
                Int_t nCases2 = nCases;
                
                // For data, there is only one case in this example
                if ( titleDataMC[idata].Contains("LHC") ) 
                  nCases2 = 1;
                
                MakeDoubleRatios
                (titleDataMC[idata], gjProd,
                 nCases2, mcCasesPath, mcCases, mcLeg, 
                 igj, gjFactor[ifactor], ipi0, rebin,
                 m02LowCutMin     , m02LowCutMax[ish],
                 m02HigCutMin[jsh], m02HigCutMax,
                 debug );
              } // use Pi0
            } // jsh
          } // ish
        } // idata
      } // GJ scale factor
    } // add GJ
  } // make double ratios
  
  
  // Calculate and plot the purity
  //
  //nfactors = 0;
  if ( doPurity )
  {
    if ( debug )
      printf("=== Make purity === \n");
    
    // Data and JJ low
    for(Int_t ifactor = 0; ifactor < nfactors+1; ifactor++)
    {
      for(Int_t ish = 0; ish < nShShCutsMaxLow; ish++)
      {
        for(Int_t jsh = 0; jsh < nShShCutsMinHigh; jsh++)
        {
          if ( m02LowCutMax[ish] > m02HigCutMin[jsh] ) continue;

          for(Int_t ipi0 = 0; ipi0 < npi0; ipi0++)
          {
            if ( ipi0 && (ish > 1 || jsh > 0)  ) 
              continue;
            
            if ( debug )
              printf("\t usepi0 %d, GJ factor %d-%2.1f, ish %d %2.2f, jsh %d %2.2f\n", 
                     ipi0,ifactor, gjFactor[ifactor],
                     ish, m02LowCutMax[ish],
                     jsh, m02HigCutMin[jsh]);
            
            // JJ low
            MakePurity(titleDataMC[0],titleDataMC[1],  
                       nCases, mcCasesPath, mcCases, mcLeg,
                       gjFactorString[ifactor], ipi0,
                       m02LowCutMin     , m02LowCutMax[ish],
                       m02HigCutMin[jsh], m02HigCutMax, debug);
            
            // JJ High
            MakePurity(titleDataMC[0],titleDataMC[2], 
                       nCases, mcCasesPath, mcCases, mcLeg,
                       gjFactorString[ifactor], ipi0,
                       m02LowCutMin     , m02LowCutMax[ish],
                       m02HigCutMin[jsh], m02HigCutMax, debug);
          } // ipi0
        } // jsh
      }// ish
    } // GJ factor
  } // Make purity
}


