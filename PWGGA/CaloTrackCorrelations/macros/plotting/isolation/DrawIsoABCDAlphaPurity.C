///
/// \file DrawIsoABCDAlphaPurity.C
/// \ingroup CaloTrackCorrMacrosPlotting
/// \brief Calculate isolated photon purity
///
/// Example macro to calculate isolated photon purity for different data/MC productions.
/// Input histograms from data, gamma-jet MC and multiple jet-jet MC productions, defined in "period" array.
/// Order matters, data must always be first, then gamma-jet and then the different jet-jet.
/// Different ways depending the input histograms selected, 
///   * 1D: pt with shower shape and isolation cuts already applied, 
///   * 2D: pt vs shower shape, isolation cut already applied
///   * 3D: pt vs shower shape vs sum pT in cone (isolation)
///     
///   Output, plots with ABCD spectrum and ratios, alpha MC correction and purity and corresponding ratios
///   Also purity, alpha and ABCD spectra stored in new root file 
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)
///

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TFile.h>
#include <TString.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#endif


/// Common general definitions
///

/// Name of directory with productions, order matters
TString period [] = {"Data","GJ","JJ","JJlow","JJhigh","MB"};
TString perString[] = 
{
  "Data",
  "GJ",
  "JJ",
  "JJ #it{p}_{T}^{#gamma}>3.5 GeV/#it{c}",
  "JJ #it{p}_{T}^{#gamma}>7 GeV/#it{c}",
  "MB"
};
const Int_t nPeriod  = 5;

TFile* f[nPeriod];

// Histograms points Style
Int_t color[] = {kBlack,kBlue,kGreen-2,kRed,kViolet};
Int_t marker[] = {20,24,25,26,27,28};
Int_t lineWidth = 2;
Bool_t makeRelStatErrPlots = kFALSE;
Bool_t makePurityRatioFit  = kFALSE;

//-----------------------------------------------------------------------------
///  Method to get relative statistical error histogram
//-----------------------------------------------------------------------------
TH1F* GetRelErrorHistogram(TH1F* h,TString tag)
{
  TH1F * hRelErr = (TH1F*) h->Clone(Form("%sRelErr%s",h->GetName(),tag.Data()));
  
  for(Int_t ibin = 1; ibin < h->GetNbinsX();ibin++)
  {
    Double_t content = hRelErr->GetBinContent(ibin);
    Double_t error   = hRelErr->GetBinError(ibin);
    if ( content > 0) 
    {
      //printf("bin %d, content %e, error %e, error/content %e\n",ibin,content,error,error/content);
      hRelErr->SetBinContent(ibin,error/content);
    }
    else
    {
      //printf("bin %d, content %e, error %e, skip!\n",ibin,content,error);
      hRelErr->SetBinContent(ibin,0);
    }
    hRelErr->SetBinError  (ibin,0);
  }
  
  hRelErr->SetYTitle("Stat error / value");
  
  return hRelErr;
}

//-----------------------------------------------------------------------------
/// Method to scale each bin by it size, needed for spectra plotting and different bin size in spectrum
//-----------------------------------------------------------------------------
static void ScaleBinBySize(TH1F* h)
{
  for(Int_t ibin = 1; ibin < h->GetNbinsX();ibin++)
  {
    Double_t width   = h->GetBinWidth(ibin);
    Double_t content = h->GetBinContent(ibin);
    Double_t error   = h->GetBinError(ibin);
    //printf("bin %d, width %f, content %e\n",ibin,width,content);
    //if(content/width > 0.5 * error/width)
    {
      h->SetBinContent(ibin,content/width);
      h->SetBinError  (ibin,error/width);
    }
//    else 
//    {
//      h->SetBinContent(ibin,0);
//      h->SetBinError  (ibin,0);
//    }
      
  }
}

//-----------------------------------------------------------------------------
/// Main execution method doing the plots. It also produces root files with main produced histograms, 
/// alpha MC correlation factor and purity.
/// Pass, 
///   * opt: 1,2,3 - Dimension of the original histograms 1D, 2D, 3D
///   * The different cuts defining the A (signal) - BCD  (background) regions:
///      * Isolation: sumSigMin, sumSigMax ,sumBkgMin, sumBkgMax,
///      * Shower shape: shSigMin, shSigMax ,shBkgMin, shBkgMax ,
///    Note that for op1, they do not have any use, for opt2, only the shower shape cuts select regions
///    and for opt3 all are used.
///   * addGJ : Take the gamma-jet production and add it to the background BCD histograms 
///   * fragInBkg : Consider the fragmentation photons only in the BCD region
///   * fragInBkg :  Remove fragmentation photons from any Jet-Jet distribution
///   * scaleBiasedJJ: constant scale factor to be applied to biased jet-jet MC only (just for testing, not to use)
///   * scaleGJ: constant scale factor to signal gamma-jet simulations, for systematic checks. 
//-----------------------------------------------------------------------------
void Exec( Int_t opt = 2, 
          Bool_t addGJ               = kTRUE, 
          Bool_t fragInBkg           = kFALSE,
          Bool_t removeFragmentation = kTRUE,
          Float_t scaleBiasedJJ = 0.,
          Float_t scaleGJ = 0.,
          Float_t ptmin     = 8.0,  Float_t ptmax     = 75, 
          Float_t sumSigMin = -200, Float_t sumSigMax = 1.5, 
          Float_t sumBkgMin = 2.5 , Float_t sumBkgMax = 20,
          Float_t shSigMin  = 0.1,  Float_t shSigMax  = 0.27, 
          Float_t shBkgMin  = 0.35, Float_t shBkgMax  = 5)
{ 
  
  // Set up the legends of the plots
  //
  TLegend *leg = new TLegend(0.55,0.6,0.95,0.89);
  leg->SetTextSize(0.045);
  leg->SetFillColor(kWhite);
  leg->SetLineColor(0);
  //leg->SetBorderColor(0);
  
  TLegend *leg2 = new TLegend(0.15,0.6,0.45,0.89);
  leg2->SetTextSize(0.045);
  leg2->SetFillColor(kWhite);
  leg2->SetLineColor(0);
  //leg2->SetBorderColor(0);
  
  TLegend *legR = new TLegend(0.2,0.7,0.4,0.89);
  legR->SetTextSize(0.045);
  legR->SetFillColor(kWhite);
  legR->SetLineColor(0);
  
  // Define here the binning of the pT histograms
  const Int_t nBins = 20;
  Double_t edges[] = {10,12,14,16,18,20,25,30,35,40,50,60,70,80,100,120,140,160,180,200,220};
 
//  const Int_t nBins = 17;
//  Double_t edges[] = {10,12,14,16,18,20,25,30,35,40,50,60,70,80,100,140,180,220};
  
  // Define histograms for ABCD distributions
  //
  // opt=3
  TH3F* h3[nPeriod];
  TH3F* h3Fragment[nPeriod];
  
  // opt=2
  TH2F* h2AB[nPeriod];
  TH2F* h2CD[nPeriod];
  TH2F* h2ABFragment[nPeriod];
  TH2F* h2CDFragment[nPeriod];
  
  // opt = 1,2,3
  TH1F* hA[nPeriod];
  TH1F* hB[nPeriod];
  TH1F* hC[nPeriod];
  TH1F* hD[nPeriod];
  
  TH1F* hAFragment[nPeriod];
  TH1F* hBFragment[nPeriod];
  TH1F* hCFragment[nPeriod];
  TH1F* hDFragment[nPeriod];
  
  TH1F* hAGJ[nPeriod];
  
  TH1F* hARatio[nPeriod];
  TH1F* hBRatio[nPeriod];
  TH1F* hCRatio[nPeriod];
  TH1F* hDRatio[nPeriod];
  
  TH1F* hARelErr[nPeriod];
  TH1F* hBRelErr[nPeriod];
  TH1F* hCRelErr[nPeriod];
  TH1F* hDRelErr[nPeriod];
  
  // Prepare strings to be attached to final figures
  //
  TString tagName = "";
  if ( addGJ ) 
  {
    tagName+="_GammaJetInBCD";
    if ( scaleGJ > 0 ) 
      tagName+=Form("Scaledx%2.1f",scaleGJ);
  }
  
  if ( removeFragmentation ) 
    tagName+="_NoFragPhoton";
  else if( fragInBkg )
    tagName+="_FragPhotonInBkg";
  
  // Loop on the different productions, data and MC used in the 
  // calculations
  //
  for(Int_t iprod = 0; iprod < nPeriod; iprod++)
  {
    //TString tag = Form("MC_%s%s",period[iprod].Data(),tagName.Data());
    TString tag = Form("%s",period[iprod].Data());
    //printf("%s\n",tag.Data());
    
    //----------------------------------------------------------------
    // Input is a TH3 histogram pT vs sum pT in cone vs shower shape
    // Not tested
    //
    if ( opt == 3 ) 
    {
      h3[iprod] = (TH3F*) f[iprod]->Get("AnaIsolPhoton_hPtM02SumPtCone");
      //printf("%p\n",h3[iprod]);
      if ( !h3[iprod] )
      {
        printf("No histo: iprod %d\n",iprod);
        hA[iprod] = 0;
        hB[iprod] = 0;
        hC[iprod] = 0;
        hD[iprod] = 0;
        continue;
      }
      
      // Remove fragmentation photons from MCs
      //
      if ( period[iprod]!="Data" )
      {
        h3Fragment[iprod] = (TH3F*) f[iprod]->Get("AnaIsolPhoton_hPtM02SumPtCone_MCPhotonFrag");
        
        if ( removeFragmentation )
          h3[iprod]->Add(h3Fragment[iprod],-1);
      }
      else h3Fragment[iprod] = 0;
      
      if (period[iprod]=="Data" || period[iprod]=="MB")
        h3[iprod]->Sumw2();
      
      // A
      h3[iprod]->SetAxisRange(ptmin    , ptmax    ,"X"); // pT trigger
      h3[iprod]->SetAxisRange(sumSigMin, sumSigMax,"Z"); // Sum pt Cone
      h3[iprod]->SetAxisRange(shSigMin , shSigMax ,"Y"); // M02
      
      hA[iprod] = (TH1F*) h3[iprod]->Project3D(Form("xA%s",tag.Data()));
      
      //B
      h3[iprod]->SetAxisRange(ptmin    , ptmax    ,"X"); // pT trigger
      h3[iprod]->SetAxisRange(sumSigMin, sumSigMax,"Z"); // Sum pt Cone
      h3[iprod]->SetAxisRange(shBkgMin , shBkgMax ,"Y"); // M02
      
      hB[iprod] = (TH1F*) h3[iprod]->Project3D(Form("xB%s",tag.Data()));
      
      //C
      h3[iprod]->SetAxisRange(ptmin    , ptmax    ,"X"); // pT trigger
      h3[iprod]->SetAxisRange(sumBkgMin, sumBkgMax,"Z"); // Sum pt Cone
      h3[iprod]->SetAxisRange(shSigMin , shSigMax ,"Y"); // M02
      
      hC[iprod] = (TH1F*) h3[iprod]->Project3D(Form("xC%s",tag.Data()));
      
      // D
      h3[iprod]->SetAxisRange(ptmin    , ptmax    ,"X"); // pT trigger
      h3[iprod]->SetAxisRange(sumBkgMin, sumBkgMax,"Z"); // Sum pt Cone
      h3[iprod]->SetAxisRange(shBkgMin , shBkgMax ,"Y"); // M02
      
      hD[iprod] = (TH1F*) h3[iprod]->Project3D(Form("xD%s",tag.Data())); 
      
      // Check what pT iso cut and what pT iso gap value was used for Not Isolated case
      // Old analysis, no gap
      hA[iprod]->SetTitle(Form("A: %2.2f<#sigma^{2}_{long}<%2.2f - #it{p}_{T}^{iso} < %2.2f GeV/#it{c}",
                               shSigMin,shSigMax,sumSigMax));
      
      hB[iprod]->SetTitle(Form("B: %2.2f<#sigma^{2}_{long}<%2.2f - #it{p}_{T}^{iso} < %2.2f GeV/#it{c}",
                               shBkgMin,shBkgMax,sumSigMax));    
      
      hC[iprod]->SetTitle(Form("C: %2.2f<#sigma^{2}_{long}<%2.2f - %2.2f < #it{p}_{T}^{iso} < %2.2f GeV/#it{c}",
                               shSigMin,shSigMax, sumBkgMin,sumBkgMax));
      
      hD[iprod]->SetTitle(Form("D: %2.2f<#sigma^{2}_{long}<%2.2f - %2.2f < #it{p}_{T}^{iso} < %2.2f GeV/#it{c}",
                               shBkgMin,shBkgMax, sumBkgMin,sumBkgMax));  
      
      
      if ( period[iprod]!="Data" && removeFragmentation && fragInBkg )
      {
        // A
        h3Fragment[iprod]->SetAxisRange(ptmin    , ptmax    ,"X"); // pT trigger
        h3Fragment[iprod]->SetAxisRange(sumSigMin, sumSigMax,"Z"); // Sum pt Cone
        h3Fragment[iprod]->SetAxisRange(shSigMin , shSigMax ,"Y"); // M02
        
        hAFragment[iprod] = (TH1F*) h3Fragment[iprod]->Project3D(Form("xAFragment%s",tag.Data()));
        
        //B
        h3Fragment[iprod]->SetAxisRange(ptmin    , ptmax    ,"X"); // pT trigger
        h3Fragment[iprod]->SetAxisRange(sumSigMin, sumSigMax,"Z"); // Sum pt Cone
        h3Fragment[iprod]->SetAxisRange(shBkgMin , shBkgMax ,"Y"); // M02
        
        hBFragment[iprod] = (TH1F*) h3Fragment[iprod]->Project3D(Form("xBFragment%s",tag.Data()));
        
        //C
        h3Fragment[iprod]->SetAxisRange(ptmin    , ptmax    ,"X"); // pT trigger
        h3Fragment[iprod]->SetAxisRange(sumBkgMin, sumBkgMax,"Z"); // Sum pt Cone
        h3Fragment[iprod]->SetAxisRange(shSigMin , shSigMax ,"Y"); // M02
        
        hCFragment[iprod] = (TH1F*) h3Fragment[iprod]->Project3D(Form("xCFragment%s",tag.Data()));
        
        // D
        h3Fragment[iprod]->SetAxisRange(ptmin    , ptmax    ,"X"); // pT trigger
        h3Fragment[iprod]->SetAxisRange(sumBkgMin, sumBkgMax,"Z"); // Sum pt Cone
        h3Fragment[iprod]->SetAxisRange(shBkgMin , shBkgMax ,"Y"); // M02
        
        hDFragment[iprod] = (TH1F*) h3Fragment[iprod]->Project3D(Form("xDFragment%s",tag.Data())); 
      }
      else
      {
        hAFragment[iprod] = 0;
        hBFragment[iprod] = 0;
        hCFragment[iprod] = 0;
        hDFragment[iprod] = 0;
      }
    } // opt = 3
    else if ( opt == 1 ) // Not tested!
    {
      hA[iprod] = (TH1F*) f[iprod]->Get("AnaIsolPhoton_hPtIsoNarrow");
      hB[iprod] = (TH1F*) f[iprod]->Get("AnaIsolPhoton_hPtIsoWide");
      hC[iprod] = (TH1F*) f[iprod]->Get("AnaIsolPhoton_hPtNoIsoNarrow");
      hD[iprod] = (TH1F*) f[iprod]->Get("AnaIsolPhoton_hPtNoIsoWide");
      
      // Check what fixed values were used
      hA[iprod]->SetTitle("A: 0.1<#sigma^{2}_{long}<0.3 - #it{p}_{T}^{iso} < 2 GeV/#it{c}");
      hB[iprod]->SetTitle("B: 0.4<#sigma^{2}_{long}<2.0 - #it{p}_{T}^{iso} < 2 GeV/#it{c}");    
      hC[iprod]->SetTitle("C: 0.1<#sigma^{2}_{long}<0.3 - #it{p}_{T}^{iso} > 2.5 GeV/#it{c}");
      hD[iprod]->SetTitle("D: 0.4<#sigma^{2}_{long}<2.0 - #it{p}_{T}^{iso} > 2.5 GeV/#it{c}");   
      
      if ( period[iprod]!="Data" )
      {
        hAFragment[iprod] = (TH1F*) f[iprod]->Get("AnaIsolPhoton_hPtIsoNarrow_MCPhotonFrag");
        hBFragment[iprod] = (TH1F*) f[iprod]->Get("AnaIsolPhoton_hPtIsoWide_MCPhotonFrag");
        hCFragment[iprod] = (TH1F*) f[iprod]->Get("AnaIsolPhoton_hPtNoIsoNarrow_MCPhotonFrag");
        hDFragment[iprod] = (TH1F*) f[iprod]->Get("AnaIsolPhoton_hPtNoIsoWide_MCPhotonFrag");
        
        if ( removeFragmentation )
        {
          hA[iprod]->Add(hAFragment[iprod],-1);
          hB[iprod]->Add(hBFragment[iprod],-1);
          hC[iprod]->Add(hCFragment[iprod],-1);
          hD[iprod]->Add(hDFragment[iprod],-1);
        }
      }
      else
      {
        hAFragment[iprod] = 0;
        hBFragment[iprod] = 0;
        hCFragment[iprod] = 0;
        hDFragment[iprod] = 0;
      }
    } // opt = 1
    else if(opt == 2)
    {
      h2AB[iprod] = (TH2F*) f[iprod]->Get("AnaIsolPhoton_hPtLambda0Iso");
      h2CD[iprod] = (TH2F*) f[iprod]->Get("AnaIsolPhoton_hPtLambda0NoIso");
      
      // Fragmentation photon component
      // remove here if requested
      //
      if ( period[iprod]!="Data" )
      {
        //printf("Get fragment %d\n",iprod);
        h2ABFragment[iprod] = (TH2F*) f[iprod]->Get("AnaIsolPhoton_hPtLambda0Iso_MCPhotonFrag");
        h2CDFragment[iprod] = (TH2F*) f[iprod]->Get("AnaIsolPhoton_hPtLambda0NoIso_MCPhotonFrag");
        
        if ( removeFragmentation )
        {
          h2AB[iprod] ->Add(h2ABFragment[iprod],-1);
          h2CD[iprod] ->Add(h2CDFragment[iprod],-1);
          
          if ( period[iprod]=="MB" ) // Cannot tag fragmentation on MB???
          {
            h2AB[iprod] = (TH2F*) f[iprod]->Get("AnaIsolPhoton_hPtLambda0Iso_MCPi0Decay");
            h2CD[iprod] = (TH2F*) f[iprod]->Get("AnaIsolPhoton_hPtLambda0NoIso_MCPi0Decay");
            h2AB[iprod] ->Add( (TH2F*) f[iprod]->Get("AnaIsolPhoton_hPtLambda0Iso_MCEtaDecay"  ),1);
            h2CD[iprod] ->Add( (TH2F*) f[iprod]->Get("AnaIsolPhoton_hPtLambda0NoIso_MCEtaDecay"),1);
            h2AB[iprod] ->Add( (TH2F*) f[iprod]->Get("AnaIsolPhoton_hPtLambda0Iso_MCEta"  ),1);
            h2CD[iprod] ->Add( (TH2F*) f[iprod]->Get("AnaIsolPhoton_hPtLambda0NoIso_MCEta"),1);
            h2AB[iprod] ->Add( (TH2F*) f[iprod]->Get("AnaIsolPhoton_hPtLambda0Iso_MCPi0"  ),1);
            h2CD[iprod] ->Add( (TH2F*) f[iprod]->Get("AnaIsolPhoton_hPtLambda0NoIso_MCPi0"),1);
            //          h2AB[iprod] ->Add( (TH2F*) f[iprod]->Get("AnaIsolPhoton_hPtLambda0Iso_MCOtherDecay"  ),1); // Direct here?
            //          h2CD[iprod] ->Add( (TH2F*) f[iprod]->Get("AnaIsolPhoton_hPtLambda0NoIso_MCOtherDecay"),1);
            h2AB[iprod] ->Add( (TH2F*) f[iprod]->Get("AnaIsolPhoton_hPtLambda0Iso_MCElectron"  ),1);
            h2CD[iprod] ->Add( (TH2F*) f[iprod]->Get("AnaIsolPhoton_hPtLambda0NoIso_MCElectron"),1);
            h2AB[iprod] ->Add( (TH2F*) f[iprod]->Get("AnaIsolPhoton_hPtLambda0Iso_MCHadron"  ),1);
            h2CD[iprod] ->Add( (TH2F*) f[iprod]->Get("AnaIsolPhoton_hPtLambda0NoIso_MCHadron"),1);
          }
        }
        else
        {
          h2ABFragment[iprod] = 0;
          h2CDFragment[iprod] = 0;
        }
      } // not data
      
//      printf("iprod %d - 2AB %p 2CD %p - Frag 2AB %p 2CD %p \n", iprod,
//             h2AB[iprod],h2CD[iprod],h2ABFragment[iprod],h2CDFragment[iprod]);
      
      if ( !h2AB[iprod] || !h2CD[iprod] ) 
      {
        printf("opt2, no 2D files for %s\n",tag.Data());
        return;
      }
      
      Int_t shSigMinBin = h2AB[iprod]->GetYaxis()->FindBin(shSigMin);
      Int_t shSigMaxBin = h2AB[iprod]->GetYaxis()->FindBin(shSigMax);  
      Int_t shBkgMinBin = h2AB[iprod]->GetYaxis()->FindBin(shBkgMin);
      Int_t shBkgMaxBin = h2AB[iprod]->GetYaxis()->FindBin(shBkgMax);
      
      hA[iprod] = (TH1F*) h2AB[iprod]->ProjectionX(Form("A_%s",tag.Data()),shSigMinBin,shSigMaxBin);
      hB[iprod] = (TH1F*) h2AB[iprod]->ProjectionX(Form("B_%s",tag.Data()),shBkgMinBin,shBkgMaxBin);
      hC[iprod] = (TH1F*) h2CD[iprod]->ProjectionX(Form("C_%s",tag.Data()),shSigMinBin,shSigMaxBin);
      hD[iprod] = (TH1F*) h2CD[iprod]->ProjectionX(Form("D_%s",tag.Data()),shBkgMinBin,shBkgMaxBin);
      
      hA[iprod]->SetTitle(Form("A: %2.2f<#sigma^{2}_{long}<%2.2f - #it{p}_{T}^{iso} < 2 GeV/#it{c}",shSigMin,shSigMax));
      hB[iprod]->SetTitle(Form("B: %2.2f<#sigma^{2}_{long}<%2.2f - #it{p}_{T}^{iso} < 2 GeV/#it{c}",shBkgMin,shBkgMax));    
      hC[iprod]->SetTitle(Form("C: %2.2f<#sigma^{2}_{long}<%2.2f - #it{p}_{T}^{iso} > 2 GeV/#it{c}",shSigMin,shSigMax));
      hD[iprod]->SetTitle(Form("D: %2.2f<#sigma^{2}_{long}<%2.2f - #it{p}_{T}^{iso} > 2 GeV/#it{c}",shBkgMin,shBkgMax));    
      
      // Fragmentation photon component projection
      //
      if ( period[iprod]!="Data" )
      {
        //printf("iprod %d, %p\n",iprod,h2ABFragment[iprod]);
        hAFragment[iprod] = (TH1F*) h2ABFragment[iprod]->ProjectionX(Form("AFragment%s",tag.Data()),shSigMinBin,shSigMaxBin);
        hBFragment[iprod] = (TH1F*) h2ABFragment[iprod]->ProjectionX(Form("BFragment%s",tag.Data()),shBkgMinBin,shBkgMaxBin);
        hCFragment[iprod] = (TH1F*) h2CDFragment[iprod]->ProjectionX(Form("CFragment%s",tag.Data()),shSigMinBin,shSigMaxBin);
        hDFragment[iprod] = (TH1F*) h2CDFragment[iprod]->ProjectionX(Form("DFragment%s",tag.Data()),shBkgMinBin,shBkgMaxBin);
      }
      else
      {
        hAFragment[iprod] = 0;
        hBFragment[iprod] = 0;
        hCFragment[iprod] = 0;
        hDFragment[iprod] = 0;
      }
//      printf("iprod %d - A %p B %p C %p D %p - Frag A %p B %p  C %p D %p \n", iprod,
//             hA[iprod],hB[iprod],hC[iprod],hD[iprod],
//             hAFragment[iprod],hBFragment[iprod],hCFragment[iprod],hDFragment[iprod]);
    } // opt = 2
    else return;
        
    if ( !hA[iprod] ) continue;
    
    
    //----------------------------------------------------------------
    // Apply Sumw2 if needed and scale data and MB distributions
    //
    //printf("Scale\n");

    if ( period[iprod] == "MB" || period[iprod] == "Data" )
    {
      hA[iprod]->Sumw2();
      hB[iprod]->Sumw2();
      hC[iprod]->Sumw2();
      hD[iprod]->Sumw2();
      
      Double_t nEvents = ((TH1F*) f[iprod]->Get("hNEvents"))->GetBinContent(1);
      //printf("MB events %e\n",nEvents);
      
      // Set scale factor for MB or data productions to better compare spectra shape
      Float_t scale = 70./nEvents;
      
      // Arbitrary factor to match data and MC
      if(period[iprod] == "Data") 
        scale = 1./nEvents/200.;
      
      hA[iprod]->Scale(scale);
      hB[iprod]->Scale(scale);
      hC[iprod]->Scale(scale);
      hD[iprod]->Scale(scale);
    }
    
    // Apply arbitrary scaling factors to JJ or GJ if requested
    if ( scaleBiasedJJ > 0 && iprod > 2 )
    {
      //printf("per %d Scale %2.2f\n",iprod,scaleBiasedJJ);
      hA[iprod]->Scale(1./scaleBiasedJJ);
      hB[iprod]->Scale(1./scaleBiasedJJ);
      hC[iprod]->Scale(1./scaleBiasedJJ);
      hD[iprod]->Scale(1./scaleBiasedJJ);
    }
    
    if ( scaleGJ > 0 && iprod == 1)
    {
      //printf("scale GJ\n");
      hA[iprod]->Scale(scaleGJ);
      hB[iprod]->Scale(scaleGJ);
      hC[iprod]->Scale(scaleGJ);
      hD[iprod]->Scale(scaleGJ);
    }
        
    //----------------------------------------------------------------
    // Rebin pT spectrum.
    // It might not be needed for opt=3, since there already non constant binninng used
    //
    //printf("Rebin\n");

    hA[iprod] = (TH1F*) hA[iprod]->Rebin(nBins,Form("A_%s_rb",tag.Data()),edges);
    hB[iprod] = (TH1F*) hB[iprod]->Rebin(nBins,Form("B_%s_rb",tag.Data()),edges);
    hC[iprod] = (TH1F*) hC[iprod]->Rebin(nBins,Form("C_%s_rb",tag.Data()),edges);
    hD[iprod] = (TH1F*) hD[iprod]->Rebin(nBins,Form("D_%s_rb",tag.Data()),edges);
    
    ScaleBinBySize(hA[iprod]);
    ScaleBinBySize(hB[iprod]);
    ScaleBinBySize(hC[iprod]);
    ScaleBinBySize(hD[iprod]);
    
    if ( period[iprod]!="Data" )
    {
      hAFragment[iprod] = (TH1F*) hAFragment[iprod]->Rebin(nBins,Form("AFragment_%s_rb",tag.Data()),edges);
      hBFragment[iprod] = (TH1F*) hBFragment[iprod]->Rebin(nBins,Form("BFragment_%s_rb",tag.Data()),edges);
      hCFragment[iprod] = (TH1F*) hCFragment[iprod]->Rebin(nBins,Form("CFragment_%s_rb",tag.Data()),edges);
      hDFragment[iprod] = (TH1F*) hDFragment[iprod]->Rebin(nBins,Form("DFragment_%s_rb",tag.Data()),edges);
      
      ScaleBinBySize(hAFragment[iprod]);
      ScaleBinBySize(hBFragment[iprod]);
      ScaleBinBySize(hCFragment[iprod]);
      ScaleBinBySize(hDFragment[iprod]);
    }
    
    //    Int_t rb = 2;
    //    hA[iprod]->Rebin(rb);
    //    hB[iprod]->Rebin(rb);
    //    hC[iprod]->Rebin(rb);
    //    hD[iprod]->Rebin(rb);
    
    //----------------------------------------------------------------
    // In case of fragmentation photons only in bkg regions,  
    // since not removed before, remove only from signal region
    //

    if ( !removeFragmentation && fragInBkg && period[iprod] != "Data" )
    {
      //printf("Remove frag from A\n");
      hA[iprod]->Add(hAFragment[iprod],-1); 
    }

    // Add Gamma-Jet to BCD regions in Jet-Jet MCs for Purity_DD calculation
    // make sure data and GJ samples are recovered before
    
    // Just needed later for purity DD comparisons
    hAGJ[iprod] = (TH1F*) hA[iprod]->Clone(Form("%s_AddedGJ",hA[iprod]->GetName()));
    
    if ( addGJ && iprod > 1 )
    {
      //printf("Add GJ %d\n",iprod);
      // Use only for Purity_DD in MC
      hAGJ[iprod]->Add(hA[1]);
      
      hB[iprod]->Add(hB[1]);
      hC[iprod]->Add(hC[1]);
      hD[iprod]->Add(hD[1]);
    }
    
    
    //----------------------------------------------------------------
    // Histogram style settings
    //
    //printf("Set style\n");

    hAGJ[iprod]->SetLineWidth  (lineWidth);
    hAGJ[iprod]->SetLineColor  (color [iprod]);
    hAGJ[iprod]->SetMarkerColor(color [iprod]);
    hAGJ[iprod]->SetMarkerStyle(marker[iprod]);

    hA[iprod]->SetLineWidth  (lineWidth);
    hA[iprod]->SetLineColor  (color [iprod]);
    hA[iprod]->SetMarkerColor(color [iprod]);
    hA[iprod]->SetMarkerStyle(marker[iprod]);
    
    hB[iprod]->SetLineWidth  (lineWidth);
    hB[iprod]->SetLineColor  (color [iprod]);
    hB[iprod]->SetMarkerColor(color [iprod]);
    hB[iprod]->SetMarkerStyle(marker[iprod]);
    
    hC[iprod]->SetLineWidth  (lineWidth);
    hC[iprod]->SetLineColor  (color [iprod]);
    hC[iprod]->SetMarkerColor(color [iprod]);
    hC[iprod]->SetMarkerStyle(marker[iprod]);
    
    hD[iprod]->SetLineWidth  (lineWidth);
    hD[iprod]->SetLineColor  (color [iprod]);
    hD[iprod]->SetMarkerColor(color [iprod]);
    hD[iprod]->SetMarkerStyle(marker[iprod]);
    
    hA[iprod]->SetAxisRange(ptmin,ptmax);
    hB[iprod]->SetAxisRange(ptmin,ptmax);
    hC[iprod]->SetAxisRange(ptmin,ptmax);
    hD[iprod]->SetAxisRange(ptmin,ptmax);
    
    if(period[iprod].Contains("MB"))
    {
      hA[iprod]->SetAxisRange(ptmin,50);
      hB[iprod]->SetAxisRange(ptmin,50);
      hC[iprod]->SetAxisRange(ptmin,50);
      hD[iprod]->SetAxisRange(ptmin,50);
    }
    else if(period[iprod].Contains("low"))
    {
      hA[iprod]->SetAxisRange(ptmin,80);
      hB[iprod]->SetAxisRange(ptmin,80);
      hC[iprod]->SetAxisRange(ptmin,80);
      hD[iprod]->SetAxisRange(ptmin,80);
    }
    else if(period[iprod].Contains("high"))
    {
      hA[iprod]->SetAxisRange(16,ptmax);
      hB[iprod]->SetAxisRange(16,ptmax);
      hC[iprod]->SetAxisRange(16,ptmax);
      hD[iprod]->SetAxisRange(16,ptmax);
    }
    
    leg ->AddEntry(hA[iprod],Form("%s",perString[iprod].Data()),"LP");
    leg2->AddEntry(hA[iprod],Form("%s",perString[iprod].Data()),"LP");
    
    //----------------------------------------------------------------
    // Make relative error histograms and ratios for comparisons biased to unbiased MC
    //
    //printf("Relative error and ratios \n");

    hARelErr[iprod] = GetRelErrorHistogram(hA[iprod],Form("hARelErr_%s",tag.Data()));
    hBRelErr[iprod] = GetRelErrorHistogram(hB[iprod],Form("hBRelErr_%s",tag.Data()));
    hCRelErr[iprod] = GetRelErrorHistogram(hC[iprod],Form("hCRelErr_%s",tag.Data()));
    hDRelErr[iprod] = GetRelErrorHistogram(hD[iprod],Form("hDRelErr_%s",tag.Data()));
    
    // Ratio to unbiased
    if ( iprod > 2 ) // neither GJ nor JJ
    {
      hARatio[iprod] = (TH1F*) hA[iprod]->Clone(Form("hARatioTo%d_%s",iprod,tag.Data()));
      hBRatio[iprod] = (TH1F*) hB[iprod]->Clone(Form("hBRatioTo%d_%s",iprod,tag.Data()));
      hCRatio[iprod] = (TH1F*) hC[iprod]->Clone(Form("hCRatioTo%d_%s",iprod,tag.Data()));
      hDRatio[iprod] = (TH1F*) hD[iprod]->Clone(Form("hDRatioTo%d_%s",iprod,tag.Data()));
      
      hARatio[iprod]->Divide(hARatio[iprod],hA[2],1,1,"B");
      hBRatio[iprod]->Divide(hBRatio[iprod],hB[2],1,1,"B");
      hCRatio[iprod]->Divide(hCRatio[iprod],hC[2],1,1,"B");
      hDRatio[iprod]->Divide(hDRatio[iprod],hD[2],1,1,"B");
      legR->AddEntry(hARatio[iprod],Form("%s",perString[iprod].Data()),"LP");
    }
  }

  //----------------------------------------------------------------
  // From here make ABCD plots
  //----------------------------------------------------------------
  //printf("PLOT ABCD Spectra and Ratios\n");
  
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(000000);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetPadTopMargin(0.1);
  
  TCanvas * cABCD = new TCanvas
  (Form("cABCD%s",tagName.Data()),
   Form("ABCD%s" ,tagName.Data()),2*800,2*600);
  cABCD->Divide(2,2);
  
  cABCD->cd(1);
  gPad->SetLogy();
  //gPad->SetLogx();
  //gPad->SetGridx();
  
  hA[0]->Draw("");

  for(Int_t iprod = 0; iprod < nPeriod; iprod++)
  {
    if ( !hA[iprod] ) continue;  
  
    hA[iprod]->SetMinimum(1e-13);
    hA[iprod]->GetXaxis()->SetMoreLogLabels(kTRUE);
    hA[iprod]->SetYTitle("Cross section (mb)");
    
    hA[iprod]->Draw("same");
  }
  
  leg->Draw("same");

  cABCD->cd(2);
  gPad->SetLogy();
  //gPad->SetLogx();
  //gPad->SetGridx();
  
  hB[0]->Draw("");

  for(Int_t iprod = 0; iprod < nPeriod; iprod++)
  {
    if ( !hB[iprod] ) continue;  
  
    hB[iprod]->SetMinimum(1e-13);
    hB[iprod]->GetXaxis()->SetMoreLogLabels(kTRUE);
    hB[iprod]->SetYTitle("Cross section (mb)");
    
    hB[iprod]->Draw("same");
  }
  
  cABCD->cd(3);
  gPad->SetLogy();
  //gPad->SetLogx();
  //gPad->SetGridx();
  
  hC[0]->Draw("");

  for(Int_t iprod = 0; iprod < nPeriod; iprod++)
  {
    if ( !hC[iprod] ) continue;  
    
    hC[iprod]->SetMinimum(1e-13);
    hC[iprod]->GetXaxis()->SetMoreLogLabels(kTRUE);
    hC[iprod]->SetYTitle("Cross section (mb)");
    
    hC[iprod]->Draw("same");
  }
  
  cABCD->cd(4);
  gPad->SetLogy();
  //gPad->SetLogx();
  //gPad->SetGridx();
  
  hD[0]->Draw("");

  for(Int_t iprod = 0; iprod < nPeriod; iprod++)
  {
    if ( !hD[iprod] ) continue;  
  
    hD[iprod]->SetMinimum(1e-13);
    hD[iprod]->GetXaxis()->SetMoreLogLabels(kTRUE);
    hD[iprod]->SetYTitle("Cross section (mb)");
    
    hD[iprod]->Draw("same");
  }
  
  cABCD->Print(Form("figures/ABCD%s.eps",tagName.Data()));  
  
  //--------------------------------------------------------
  
  if ( makeRelStatErrPlots )
  {
    TCanvas * cABCDRelErr = new TCanvas
    (Form("cABCDRelErr%s",tagName.Data()),
     Form("ABCDRelErr%s" ,tagName.Data()),2*800,2*600);
    cABCDRelErr->Divide(2,2);
    
    cABCDRelErr->cd(1);
    //gPad->SetLogy();
    //gPad->SetLogx();
    //gPad->SetGridx();
    
    hARelErr[0]->Draw("");
    
    for(Int_t iprod = 0; iprod < nPeriod; iprod++)
    {
      if ( !hARelErr[iprod] ) continue;
      
      //hARelErr[iprod]->SetMinimum(0.001);
      hARelErr[iprod]->SetYTitle("Stat err / content");
      hARelErr[iprod]->GetXaxis()->SetMoreLogLabels(kTRUE);
      
      if(hARelErr[iprod]->GetMaximum() > hARelErr[0]->GetMaximum()) 
        hARelErr[0]->SetMaximum(hARelErr[iprod]->GetMaximum()*1.2);
      
      hARelErr[iprod]->Draw("same");
    }
    
    leg2->Draw("same");
    
    cABCDRelErr->cd(2);
    //gPad->SetLogy();
    //gPad->SetLogx();
    //gPad->SetGridx();
    
    hBRelErr[0]->Draw("");
    
    for(Int_t iprod = 0; iprod < nPeriod; iprod++)
    {
      if ( !hBRelErr[iprod] ) continue;
      
      //hBRelErr[iprod]->SetMinimum(0.001);
      hBRelErr[iprod]->SetYTitle("Stat err / content");
      hBRelErr[iprod]->GetXaxis()->SetMoreLogLabels(kTRUE);
      
      if(hBRelErr[iprod]->GetMaximum() > hBRelErr[0]->GetMaximum()) 
        hBRelErr[0]->SetMaximum(hBRelErr[iprod]->GetMaximum()*1.2);
      
      hBRelErr[iprod]->Draw("same");
    }
    leg2->Draw("same");
    
    cABCDRelErr->cd(3);
    //gPad->SetLogy();
    //gPad->SetLogx();
    //gPad->SetGridx();
    
    hCRelErr[0]->Draw("");
    
    for(Int_t iprod = 0; iprod < nPeriod; iprod++)
    {
      if ( !hCRelErr[iprod] ) continue;
      
      //hCRelErr[iprod]->SetMinimum(0.001);
      hCRelErr[iprod]->SetYTitle("Stat err / content");
      hCRelErr[iprod]->GetXaxis()->SetMoreLogLabels(kTRUE);
      
      if(hCRelErr[iprod]->GetMaximum() > hCRelErr[0]->GetMaximum()) 
        hCRelErr[0]->SetMaximum(hCRelErr[iprod]->GetMaximum()*1.2);
      
      hCRelErr[iprod]->Draw("same");
    }
    leg2->Draw("same");
    
    cABCDRelErr->cd(4);
    //gPad->SetLogy();
    //gPad->SetLogx();
    //gPad->SetGridx();
    
    hDRelErr[0]->Draw("");
    
    for(Int_t iprod = 0; iprod < nPeriod; iprod++)
    {
      if ( !hDRelErr[iprod] ) continue;
      
      //hDRelErr[iprod]->SetMinimum(0.001);
      hDRelErr[iprod]->SetYTitle("Stat err / content");
      hDRelErr[iprod]->GetXaxis()->SetMoreLogLabels(kTRUE);
      
      if(hDRelErr[iprod]->GetMaximum() > hDRelErr[0]->GetMaximum()) 
        hDRelErr[0]->SetMaximum(hDRelErr[iprod]->GetMaximum()*1.2);
      
      hDRelErr[iprod]->Draw("same");
    }
    leg2->Draw("same");
    
    cABCDRelErr->Print(Form("figures/ABCDRelErr%s.eps",tagName.Data()));  
  }
  
  //--------------------------------------------------------

  TCanvas * cABCDRatio = new TCanvas
  (Form("cABCDRatio%s",tagName.Data()),
   Form("ABCDRatio%s" ,tagName.Data()),
   2*800,2*600);
  
  cABCDRatio->Divide(2,2);
  
  cABCDRatio->cd(1);
  //gPad->SetLogy();
  gPad->SetLogx();
  gPad->SetGridy();
  
  hARatio[3]->Draw("");
  
  for(Int_t iprod = 3; iprod < nPeriod; iprod++)
  {
    if ( !hARatio[iprod] ) continue;

    hARatio[iprod]->SetMinimum(0.5);
    hARatio[iprod]->SetMaximum(1.5);
    hARatio[iprod]->SetYTitle("JJ biased / unbiased");
    hARatio[iprod]->GetXaxis()->SetMoreLogLabels(kTRUE);
    
    hARatio[iprod]->Draw("same");
  }
  legR->Draw("same");

  cABCDRatio->cd(2);
  //gPad->SetLogy();
  gPad->SetLogx();
  gPad->SetGridy();
  
  hBRatio[3]->Draw("");

  for(Int_t iprod = 3; iprod < nPeriod; iprod++)
  {
    if ( !hBRatio[iprod] ) continue;

    hBRatio[iprod]->SetYTitle("JJ biased / unbiased");
    hBRatio[iprod]->GetXaxis()->SetMoreLogLabels(kTRUE);
    hBRatio[iprod]->SetMinimum(0.5);
    hBRatio[iprod]->SetMaximum(1.5);
    
    hBRatio[iprod]->Draw("same");
  }
  legR->Draw("same");

  cABCDRatio->cd(3);
  //gPad->SetLogy();
  gPad->SetLogx();
  gPad->SetGridy();
  
  hCRatio[3]->Draw("");
  
  for(Int_t iprod = 3; iprod < nPeriod; iprod++)
  {
    if ( !hCRatio[iprod] ) continue;
      
    hCRatio[iprod]->SetMinimum(0.5);
    hCRatio[iprod]->SetMaximum(1.5);
    hCRatio[iprod]->SetYTitle("JJ biased / unbiased");
    hCRatio[iprod]->GetXaxis()->SetMoreLogLabels(kTRUE);
    
    hCRatio[iprod]->Draw("same");
  }
  legR->Draw("same");

  cABCDRatio->cd(4);
  //gPad->SetLogy();
  gPad->SetLogx();
  gPad->SetGridy();
  
  hDRatio[3]->Draw("");  
  for(Int_t iprod = 3; iprod < nPeriod; iprod++)
  {
    if ( !hDRatio[iprod] ) continue;
    
    hDRatio[iprod]->SetYTitle("JJ biased / unbiased");
    hDRatio[iprod]->SetMinimum(0.5);
    hDRatio[iprod]->SetMaximum(1.5);
    hDRatio[iprod]->GetXaxis()->SetMoreLogLabels(kTRUE);
    
    hDRatio[iprod]->Draw("same");
  }
  legR->Draw("same");

  cABCDRatio->Print(Form("figures/ABCDRatio%s.eps",tagName.Data()));  
  
  ////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////
  ///  Alpha
  ////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////
  //printf("Alpha MC calculation");

  gStyle->SetOptTitle(0);
  gStyle->SetPadTopMargin(0.01);

  TH1F* hAlpha      [nPeriod];
  TH1F* hAlphaRelErr[nPeriod];
  TH1F* hAlphaRatio [nPeriod];
  
  TH1F* hACRatio[nPeriod];
  TH1F* hBDRatio[nPeriod];
  TH1F* hCARatio[nPeriod];
  TH1F* hDBRatio[nPeriod];
  
  TLegend *legA = new TLegend(0.15,0.7,0.4,0.95);
  legA->SetTextSize(0.045);
  legA->SetFillColor(kWhite);
  legA->SetLineColor(0);
  //legA->SetBorderColor(0);
  
  // Loop on Jet-Jet productions to get the Alpha MC factor
  //
  for(Int_t iprod = 2; iprod < nPeriod; iprod++)
  {
    if ( !hA[iprod] ) continue;

    //TString tag = Form("MC_%s%s",period[iprod].Data(),tagName.Data());
    TString tag = Form("%s",period[iprod].Data());

    hACRatio[iprod] = (TH1F*) hA[iprod]->Clone(Form("hACRatio%s",tag.Data()));
    hBDRatio[iprod] = (TH1F*) hB[iprod]->Clone(Form("hACRatio%s",tag.Data()));
    
    hACRatio[iprod]->Divide(hACRatio[iprod],hC[iprod],1,1,"B");
    hBDRatio[iprod]->Divide(hBDRatio[iprod],hD[iprod],1,1,"B");
    
    hAlpha[iprod] = (TH1F*) hACRatio[iprod]->Clone(Form("Alpha_%s",tag.Data()));
    hAlpha[iprod]->Divide(hAlpha[iprod],hBDRatio[iprod],1,1,"B");
    
    legA->AddEntry(hAlpha[iprod],Form("%s",perString[iprod].Data()),"LP");
   
    hAlphaRelErr[iprod] = GetRelErrorHistogram(hA[iprod],Form("hARelErr_%s",tag.Data()));
    
    if ( iprod > 2 ) // neither GJ nor JJ
    {
      hAlphaRatio[iprod] = (TH1F*) hAlpha[iprod]->Clone(Form("hAlpha_%s_Ratio",tag.Data()));
            
      hAlphaRatio[iprod]->Divide(hAlphaRatio[iprod],hAlpha[2],1,1,"B");
    }
  }
  
  //--------------------------------------------------------
  //printf("Alpha MC plots");
  
  TCanvas * cAlpha = new TCanvas
  (Form("cAlpha%s",tagName.Data()),
   Form("Alpha MC factor %s",tagName.Data()),
   1*800,1*600);
    
  //gPad->SetLogy();
  gPad->SetLogx();
  //gPad->SetGridx();
  
  hAlpha[2]->Draw("");
  hAlpha[2]->SetTitle("");
  hAlpha[2]->SetMaximum(3.);
  if (!removeFragmentation ) 
    hAlpha[2]->SetMaximum(5.);
  hAlpha[2]->SetMinimum(0.2);
  
  for(Int_t iprod = 2; iprod < nPeriod; iprod++)
  {
    if ( !hAlpha[iprod] ) continue;

    hAlpha[iprod]->SetYTitle("#alpha_{MC}");
   
    hAlpha[iprod]->GetXaxis()->SetMoreLogLabels(kTRUE);
   
    hAlpha[iprod]->Draw("same");
  }
  legA->Draw("same");

  cAlpha->Print(Form("figures/AlphaMC%s.eps",tagName.Data()));
 
  //--------------------------------------------------------

  if ( makeRelStatErrPlots )
  {
    TCanvas * cAlphaRelErr = new TCanvas
    (Form("cAlphaRelErr%s",tagName.Data()),
     Form("Alpha MC factor rel. err %s",tagName.Data()),
     1*800,1*600);
    
    //gPad->SetLogy();
    gPad->SetLogx();
    //gPad->SetGridx();
    
    hAlphaRelErr[2]->Draw("");
    hAlphaRelErr[2]->SetTitle("");
    hAlphaRelErr[2]->SetMaximum(1.2);
    for(Int_t iprod = 2; iprod < nPeriod; iprod++)
    {
      if ( !hAlphaRelErr[iprod] ) continue;
      
      hAlphaRelErr[iprod]->SetYTitle("Statistical error / value");
      
      hAlphaRelErr[iprod]->GetXaxis()->SetMoreLogLabels(kTRUE);
      
      hAlphaRelErr[iprod]->Draw("same");
    }
    legA->Draw("same");
    
    cAlphaRelErr->Print(Form("figures/AlphaMCRelErr%s.eps",tagName.Data()));
  }
  
  //--------------------------------------------------------

  TCanvas * cAlphaRatio = new TCanvas
  (Form("cAlphaRatio%s",tagName.Data()),
   Form("Alpha MC factor rel. err %s",tagName.Data()),
   1*800,1*600);
  
  cAlphaRatio->Divide(1,1);
   
   //c->cd(icen+1);
   //gPad->SetLogy();
   gPad->SetLogx();
   gPad->SetGridy();
   
   hAlphaRatio[4]->Draw("");
   hAlphaRatio[4]->SetTitle("");
   hAlphaRatio[4]->SetMaximum(4);
   hAlphaRatio[4]->SetMinimum(0.);
   for(Int_t iprod = 3; iprod < nPeriod; iprod++)
   {
     if ( !hAlphaRatio[iprod] ) continue;

     hAlphaRatio[iprod]->SetYTitle("#alpha_{MC} ratio JJ biased / unbiased");
     
     hAlphaRatio[iprod]->GetXaxis()->SetMoreLogLabels(kTRUE);
     
     hAlphaRatio[iprod]->Draw("same");
   }
   legR->Draw("same");

   cAlphaRatio->Print(Form("figures/AlphaMCRatio%s.eps",tagName.Data()));
  
  ////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////
  /// Purity
  ////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////
  //printf("Purity calculation");

  TH1F* hPurityDD[nPeriod];
  TH1F* hPurity  [nPeriod];
  TH1F* hPurityMC[nPeriod];
  TH1F* hPurityInput[nPeriod];
  TH1F* hPurityRatio[nPeriod];
  for(Int_t iprod = 0; iprod < nPeriod; iprod++)
  {
    if ( !hA[iprod] ) continue;

    if ( iprod == 1 ) continue; // skip GJ
    
    //TString tag = Form("MC_%s%s",period[iprod].Data(),tagName.Data());
    TString tag = Form("%s",period[iprod].Data());
    
    hCARatio[iprod] = (TH1F*) hC[iprod]->Clone(Form("hCARatio%s",tag.Data()));
    hDBRatio[iprod] = (TH1F*) hD[iprod]->Clone(Form("hDBRatio%s",tag.Data()));

    hCARatio[iprod]->Divide(hCARatio[iprod],hAGJ[iprod],1,1,"B");
    hDBRatio[iprod]->Divide(hDBRatio[iprod],hB[iprod],1,1,"B");
    
    hPurityDD[iprod] = (TH1F*) hCARatio[iprod]->Clone(Form("PurityDD%s",tag.Data()));
    hPurityDD[iprod]->Divide(hPurityDD[iprod],hDBRatio[iprod],1,1,"B");
    
    // Clone PurityDD in data and mutiply by alpha per each MC JJ 
    //
    if ( iprod == 0 )
    {
      for(Int_t iprod2 = 2; iprod2 < nPeriod; iprod2++)
      {
        hPurity[iprod2] = (TH1F*) hPurityDD[0]->Clone(Form("Purity_%s",period[iprod2].Data()));
        hPurity[iprod2]->Multiply(hAlpha[iprod2]);
      }
    }
    else if ( iprod > 1 )
    {
      hPurityMC[iprod] = (TH1F*) hPurityDD[iprod]->Clone(Form("PurityMC_%s",period[iprod].Data()));
      hPurityMC[iprod]->Multiply(hAlpha[iprod]);
    }
    
    // So far data driven contamination, get purity 1-C
    for(Int_t ix = 0; ix < hPurityDD[iprod]->GetNbinsX(); ix++)
    {
      Float_t content = hPurityDD[iprod]->GetBinContent(ix);
      if ( content > 0 ) hPurityDD[iprod]->SetBinContent(ix,1-content);
    }
    
    // Input purity, A region GJ / JJ+GJ
     if ( iprod > 1 )
     {
       hPurityInput[iprod] = (TH1F*) hA[1]->Clone(Form("hPurityInput_%s",period[iprod].Data()));
       hPurityInput[iprod]->Divide(hPurityInput[iprod],hAGJ[iprod],1,1,"B");
     }
    
    // Now purity as 1-C
    //
    if ( iprod < 2 ) continue;
    
    for(Int_t ix = 0; ix < hPurityDD[iprod]->GetNbinsX(); ix++)
    {
      Float_t content = hPurity[iprod]->GetBinContent(ix);
      if ( content > 0 ) hPurity[iprod]->SetBinContent(ix,1-content);
      
      content = hPurityMC[iprod]->GetBinContent(ix);
      if ( content > 0 ) hPurityMC[iprod]->SetBinContent(ix,1-content);
    }

    if ( iprod < 3 ) continue;
    
    hPurityRatio[iprod] = (TH1F*) hPurity[iprod]->Clone(Form("Purity_%s_Ratio",tag.Data()));
    hPurityRatio[iprod]->Divide(hPurityRatio[iprod],hPurity[2],1,1,"B");
    
  } // iprod

  //printf("Purity MC plots");

  //--------------------------------------------------------
  
  TLegend *legPdd = new TLegend(0.15,0.7,0.45,0.97);
  legPdd->SetTextSize(0.045);
  legPdd->SetFillColor(kWhite);
  legPdd->SetLineColor(0);
  //legPdd->SetBorderColor(0);
    
  TCanvas * cPurityDD = new TCanvas
  (Form("cPurityDD%s",tagName.Data()),
   Form("Purity data-driven %s",tagName.Data()),
   1*800,1*600);
  
  cPurityDD->Divide(1,1);
  
  //c->cd(icen+1);
  //gPad->SetLogy();
  gPad->SetLogx();
  //gPad->SetGridx();
  
  hPurityDD[0]->Draw("");
  hPurityDD[0]->SetTitle("");
  hPurityDD[0]->SetYTitle("P_{dd}");
  hPurityDD[0]->SetMaximum(1.2);
  hPurityDD[0]->SetMinimum(0.1);
  legPdd->AddEntry(hPurityDD[0],perString[0].Data(),"PL");
  for(Int_t iprod = 2; iprod < nPeriod; iprod++)
  {
    if ( !hA[iprod] ) continue;

    hPurityDD[iprod]->SetAxisRange(ptmin,ptmax,"X");
    if(period[iprod].Contains("low"))  hPurityDD[iprod]->SetAxisRange(ptmin,60,"X");

    hPurityDD[iprod]->GetXaxis()->SetMoreLogLabels(kTRUE);

    legPdd->AddEntry(hPurityDD[iprod],perString[iprod].Data(),"PL");

    hPurityDD[iprod]->Draw("same");
  }
  legPdd->Draw("same");
  
  hPurityDD[0]->Draw("same");
  
  cPurityDD->Print(Form("figures/PurityDD%s.eps",tagName.Data()));
  
  //--------------------------------------------------------

  TLegend *legP = new TLegend(0.15,0.7,0.45,0.97);
  legP->SetTextSize(0.045);
  legP->SetFillColor(kWhite);
  legP->SetLineColor(0);
  //legP->SetBorderColor(0);
  
  TCanvas * cPurity = new TCanvas
  (Form("cPurity%s",tagName.Data()),
   Form("Purity%s" ,tagName.Data()),
   1*800,1*600);
  cPurity->Divide(1,1);
  
  //c->cd(icen+1);
  //gPad->SetLogy();
  gPad->SetLogx();
  gPad->SetGridy();
  
  hPurity[2]->Draw("");
 
  for(Int_t iprod = 2; iprod < nPeriod; iprod++)
  {
    if ( !hPurity[iprod] ) continue;

    hPurity[iprod]->SetYTitle("Purity");
    hPurity[iprod]->SetMaximum(1.2);
    hPurity[iprod]->SetMinimum(0.1);
    
    hPurity[iprod]->SetMarkerStyle(marker[iprod]);
    hPurity[iprod]->SetMarkerColor(color[iprod]);
    hPurity[iprod]->SetLineColor  (color[iprod]);
    
    hPurity[iprod]->GetXaxis()->SetMoreLogLabels(kTRUE);

    legP->AddEntry(hPurity[iprod],perString[iprod].Data(),"PL");

    hPurity[iprod]->Draw("same");
  }
  legP->Draw("same");
  
  cPurity->Print(Form("figures/Purity%s.eps",tagName.Data()));
  
  //--------------------------------------------------------

  TCanvas * cPurityRatio = new TCanvas
  (Form("cPurityRatio%s",tagName.Data()),
   Form("PurityRatio%s" ,tagName.Data()),
   1*800,1*600);
  cPurityRatio->Divide(1,1);
  
  //c->cd(icen+1);
  //gPad->SetLogy();
  gPad->SetLogx();
  gPad->SetGridy();
  
  hPurityRatio[3]->Draw("");
 
  for(Int_t iprod = 3; iprod < nPeriod; iprod++)
  {
    if ( !hA[iprod] ) continue;

    hPurityRatio[iprod]->SetMarkerStyle(marker[iprod]);
    hPurityRatio[iprod]->SetMarkerColor(color[iprod]);
    hPurityRatio[iprod]->SetLineColor  (color[iprod]);
    
    hPurityRatio[iprod]->SetYTitle("Purity ratio JJ biased / unbiased");
    hPurityRatio[iprod]->SetMaximum(2);
    hPurityRatio[iprod]->SetMinimum(0.5);
    
    hPurityRatio[iprod]->GetXaxis()->SetMoreLogLabels(kTRUE);
      
    hPurityRatio[iprod]->Draw("same");
    
    if ( makePurityRatioFit )
    {
      if(period[iprod].Contains("low")) 
        hPurityRatio[iprod]->Fit("pol0","","NRO",20,60);
      else                            
        hPurityRatio[iprod]->Fit("pol0","","NRO",20,200);
      //hPurityRatio[iprod]->Fit("pol0","","NRO",20,100);
    }
  }
  legR->Draw("same");
  
  cPurityRatio->Print(Form("figures/PurityRatio%s.eps",tagName.Data()));

  //--------------------------------------------------------
//
//  TCanvas * cPurityMCInput = new TCanvas
//  (Form("cPurityMCInput%s",tagName.Data()),
//   Form("Purity MC and input %s",tagName.Data()),
//   1*800,1*600);
//  
//  cPurityMCInput->Divide(1,1);
//  
//  //c->cd(icen+1);
//  //gPad->SetLogy();
//  gPad->SetLogx();
//  //gPad->SetGridx();
//  
//  hPurityMC[2]->Draw("");
//  hPurityMC[2]->SetTitle("");
//  hPurityMC[2]->SetYTitle("Purity MC");
//  hPurityMC[2]->SetMaximum(1.1);
//  hPurityMC[2]->SetMinimum(0.);
//  hPurityMC[2]->GetXaxis()->SetMoreLogLabels(kTRUE);
//  
//  for(Int_t iprod = 2; iprod < nPeriod; iprod++)
//  {
//    hPurityInput[iprod]->SetMarkerStyle(marker[iprod]);
//    hPurityInput[iprod]->SetMarkerColor(kCyan);
//    hPurityInput[iprod]->SetLineColor  (kCyan);
//    hPurityInput[iprod]->SetLineStyle  (1);
//    
//    hPurityInput[iprod]->Draw("sameH");
//  }
//  
//  for(Int_t iprod = 2; iprod < nPeriod; iprod++)
//  {
//    if ( !hPurityMC[iprod] ) continue;
//    hPurityMC[iprod]->Draw("same");
//  }
//  legP->Draw("same");
//  
//  cPurityMCInput->Print(Form("figures/PurityMCvsInput%s.eps",tagName.Data()));
//  
//  //--------------------------------------------------------
//
//  TCanvas * cPurityMCInputRatio = new TCanvas
//  (Form("cPurityMCInputRatio%s",tagName.Data()),
//   Form("Purity MC and input %s",tagName.Data()),
//   1*800,1*600);
//  
//  cPurityMCInputRatio->Divide(1,1);
//  
//  //c->cd(icen+1);
//  //gPad->SetLogy();
//  gPad->SetLogx();
//  //gPad->SetGridx();
//  
//  for(Int_t iprod = 2; iprod < nPeriod; iprod++)
//  {
//    if ( !hPurityMC[iprod] ) continue;
//    
//    hPurityMC[iprod]->Divide(hPurityInput[iprod]);
//    hPurityMC[iprod]->SetYTitle("Purity Ratio MC / Input");
//    hPurityMC[iprod]->SetMaximum(1.2);
//    hPurityMC[iprod]->SetMinimum(0.8);
//    
//    if ( iprod==2 ) hPurityMC[iprod]->Draw("H");
//    else            hPurityMC[iprod]->Draw("sameH");
//  }
//  //legP->Draw("same");
//  
//  cPurityMCInputRatio->Print(Form("figures/PurityMCInputRatio%s.eps",tagName.Data()));
  
  //--------------------------------------------------------
  // Write to file
  //--------------------------------------------------------

  TFile * fout = new TFile(Form("AlphaPurity%s.root",tagName.Data()),"recreate");
  
  for(Int_t iprod = 0; iprod < nPeriod; iprod++)
  {
    hA[iprod]->Write();
    hB[iprod]->Write();
    hC[iprod]->Write();
    hD[iprod]->Write();
  }
  
  for(Int_t iprod = 2; iprod < nPeriod; iprod++)
  {
    hAlpha [iprod]->Write();
    hPurity[iprod]->Write();
  }
  
  fout->Close();
}

//-----------------------------------------------------------------------------
/// Main method
/// It recovers and opens the files. It can in principle execute the plot creation for different optiosn and 
/// cut setting but not tested, better execute once per set of cuts/options.
//-----------------------------------------------------------------------------
void DrawIsoABCDAlphaPurity()
{  
  for(Int_t iprod = 0; iprod < nPeriod; iprod++)
  {
    f[iprod] = TFile::Open(Form("%s/Scaled.root",
                               period[iprod].Data())); 
    if ( !f[iprod] )
    {
      printf("File for %s not found\n",period[iprod].Data());
      return;
    }
  }
  
  Int_t opt = 2; /// Use 1 and 3 only for recent aliroot from 03/2020
 
  Float_t ptmin     = 10.0;  Float_t ptmax     = 200;
  Float_t sumSigMin = -200;  Float_t sumSigMax = 2;
  Float_t sumBkgMin = 3 ;    Float_t sumBkgMax = 100;
  Float_t shSigMin  = 0.1 ;  Float_t shSigMax  = 0.3; 
  Float_t shBkgMin  = 0.4 ;  Float_t shBkgMax  = 1.6;
  //Float_t shBkgMin  = 0.5 ;  Float_t shBkgMax  = 2;
 
  Float_t scaleBiasedJJ = 1.; // Do not change
  Float_t scaleGJ       = 1.; // Change to 0.5 or 2
  
  Bool_t  fragInBkg           = kFALSE; // true only meaningful if removeFragmentation is false
  Bool_t  addGJ               = kTRUE;
  Bool_t  removeFragmentation = kTRUE;
  
  makeRelStatErrPlots = kFALSE;
  makePurityRatioFit  = kFALSE;

  Exec(opt, addGJ, fragInBkg,
       removeFragmentation, 
       scaleBiasedJJ, scaleGJ,
       ptmin    , ptmax,
       sumSigMin, sumSigMax, sumBkgMin, sumBkgMax,
       shSigMin , shSigMax , shBkgMin , shBkgMax  );
}
