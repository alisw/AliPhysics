///
/// \file CalculateParamChi2MCxTalkDataPerSM.C
/// \ingroup CaloTrackCorrMacrosPlotting
/// \brief Calculate Chi2 distributions for different cluster parameters 
///
/// Take the output distributions of CompareDataAndNMCProdPerSM.C and calculate the
/// Chi2 distribution provided a "data" distribution reference compared to multiple 
/// MC distributions. See how to execute it in MakeDataMCComparisonPerSMClusterEbinAndChi2.C.
/// Calculated Chi2 histograms stored in output file for further processing and calculations.
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
#include <TH3F.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TObject.h>
#include <TAxis.h>
#include <TGaxis.h>
#include <TLine.h>
#include <TF1.h>
#include <TMath.h>
#include <TArrayD.h>
#include <TGraphErrors.h>
#endif

void Chi2( TH1D* h, TH1D* hD,TH1D*& hChi2, TH1D*& hDiff, Bool_t debug = kFALSE );

TString fileFormat = ".eps";

Int_t firstEbin = 0;

//------------------------------------------------------------------------------
/// Main method
///
/// \param histoName : Name of histogram with the parameter to study = "SMM02NoCut", "SMM02NoCut","SMM02",
///                    "SMM20LowM02NoCut","SMM20LowM02","SMM20HighM02NoCut","SMM20HighM02","SMNCell"
/// \param inputFileName : String with part of input file name 
/// \param dataProd  : Index of reference data production file
/// \param firstSM   : First SM to be studied
/// \param lastSM    : Last SM to be studied
/// \param binE      : TArrayD with list of energy bin limits
/// \param prodParam : TArrayD with list of parameter that changed in the different MC
/// \param rebin     : Optional additional histogram rebinning
/// \param debug     : Bool to activate printfs.
//------------------------------------------------------------------------------
void CalculateParamChi2MCxTalkDataPerSM
(
 TString histoName  = "SMM02NoCut",
 TString inputFileName = "", 
 Int_t   dataProd  = 0,
 Int_t   firstSM   = 0,
 Int_t   lastSM    = 9,
 TArrayD binE      = 0,
 TArrayD prodParam = 0,
 Int_t   rebin     = 1,
 Bool_t  debug     = kFALSE
 )
{  
  Int_t color[] = {1,4,2,kYellow-2,8,kCyan,kOrange+2,kViolet,kOrange-2,4,2,6,8,9,kYellow-6,kCyan,kOrange+2,kViolet,kOrange-2};
  
  // Define Chi2 minimum finding window
  Float_t xmin = 0.1;
  Float_t xmax = 1.0;
  if ( histoName.Contains("SMM02") )
  {
    xmin = 0.1;
    xmax = 1.0;    
  }
  if ( histoName == "SMNCell" && !histoName.Contains("Module") )        
  { 
    xmin = 2   ; 
    xmax = 12   ; 
  }
  if  ( histoName.Contains("SMM20Low") )
  { 
    xmin = 0.05; xmax = 0.25; 
  }
  if ( histoName.Contains("SMM20Hig") )
  { 
    xmin = 0.05; xmax = 0.30; 
  }
  if ( histoName.Contains("Module") )
  { 
    xmin = 0.0; xmax = 100000; 
  }
  
  // Open the file with projected distributions
  TFile * file = TFile::Open(Form("figures/%s/Projections_%s.root",histoName.Data(),inputFileName.Data()));;
  if(!file) return;
  
  // Energy bins 
  const Int_t nEBins = binE.GetSize()-1;
  Double_t binEErr[nEBins-1];
  for(Int_t ie = 0; ie < nEBins; ie++)  
    binEErr[ie] = (binE.At(ie+1)-binE.At(ie))/2; 
  
  // MC production parameter
  const Int_t nProd = prodParam.GetSize();
  Double_t prodParamE[nProd];
  for(Int_t iprod = 0; iprod < nProd-1; iprod++)  
  {
    prodParamE[iprod] = (prodParam.At(iprod+1)-prodParam.At(iprod))/2; 
    if ( debug ) printf("mc iprod %d: param %2.5f, next %2.5f, err %2.5f\n",
                        iprod,prodParam.At(iprod),prodParam.At(iprod+1),prodParamE[iprod]);
  }
  prodParamE[nProd] = prodParamE[nProd-1];
  
  // Total number of SM inspected
  const Int_t nSM = lastSM-firstSM+1; 
  if ( debug )
    printf("N MC prod %d N e bins %d, nSM %d\n", nProd, nEBins, nSM);
  
  //
  // Recover the histograms from file
  //
  TH1D* hData       [nEBins][nSM];
  TH1D* hSimu[nProd][nEBins][nSM];
  
  for(Int_t iebin = 0; iebin < nEBins; iebin++)
  {
    for(Int_t ism = firstSM; ism <= lastSM; ism++)
    {     
      hData[iebin][ism-firstSM] = (TH1D *) file->Get(Form("Prod%d_Ebin%d_Param%d",dataProd,iebin,ism));
      if ( rebin > 1 ) 
        hData[iebin][ism-firstSM]->Rebin(rebin);
      if ( debug )
        printf("data ie %d, ism %d, %p\n",iebin,ism,hData[iebin][ism-firstSM]);
      
      hSimu[dataProd][iebin][ism-firstSM] = 0;
      for(Int_t iprod = 0; iprod < nProd; iprod++)
      {
        hSimu[iprod][iebin][ism-firstSM] = (TH1D *) file->Get(Form("Prod%d_Ebin%d_Param%d",iprod+dataProd+1,iebin,ism));
        if ( rebin > 1 ) 
          hSimu[iprod][iebin][ism-firstSM]->Rebin(rebin);
        
        if ( debug )
          printf("mc %d ie %d, ism %d, %p\n",iprod,iebin,ism,hSimu[iprod][iebin][ism-firstSM]);
        
      } // prod
    } // sm 
  } // ebin
  
  ////////////////////////////////
  // Make Chi2 and difference distributions
  ////////////////////////////////
  
  TH1D* hChi2 [nProd][nEBins][nSM];
  TH1D* hDiff [nProd][nEBins][nSM];
  
  if ( debug )
    printf("Make Chi2 distributions");
  
  for(Int_t iebin = 0; iebin < nEBins; iebin++)
  {
    for(Int_t iprod = 0; iprod < nProd; iprod++)
    {
      for(Int_t ism = firstSM; ism <= lastSM; ism++)
      {        
        hChi2[iprod][iebin][ism-firstSM] =0;   
        hDiff[iprod][iebin][ism-firstSM] = 0;
        
        if ( debug )
          printf("prod %d, ebin %d, sm %d\n",iprod,iebin,ism);
        
        Chi2(hSimu[iprod][iebin][ism-firstSM],hData[iebin][ism-firstSM],
             hChi2[iprod][iebin][ism-firstSM],hDiff[iprod][iebin][ism-firstSM], debug);
        
        if ( debug )
          printf("\t %p %p\n",hChi2[iprod][iebin][ism-firstSM],hDiff[iprod][iebin][ism-firstSM]);
        
        if(hChi2[iprod][iebin][ism-firstSM])hChi2[iprod][iebin][ism-firstSM]->SetTitle(Form("SM = %d",ism));
        if(hDiff[iprod][iebin][ism-firstSM])hDiff[iprod][iebin][ism-firstSM]->SetTitle(Form("SM = %d",ism));
      }// sm
      
    }// prod
    
  } // energy bin
  
  ////////////////////////////////
  // Make Chi2 and Difference distributions
  // Each file an energy bin, 
  // each file contains pads per SM
  ////////////////////////////////
  if ( debug ) printf("Plot Chi2-Diff\n");
  
  Int_t ncol = 4;
  Int_t nrow = 3;
  GetCanvasColRowNumber(nSM, ncol, nrow); // PlotUtils.C
  
  TString fileName = "";
  
  for(Int_t iebin = 0; iebin < nEBins; iebin++)
  {
    if ( debug ) printf("iebin %d",iebin);
    
    gStyle->SetOptTitle(1);
    gStyle->SetPadTopMargin(0.1);
    
    // Chi2
    //
    TCanvas * cchi2 = new TCanvas
    (Form("cchi2_ebin%d_%s",iebin,histoName.Data()),
     Form("Chi2 %2.1f < E < %2.1f, %s",binE[iebin],binE[iebin+1],histoName.Data()),
     ncol*2000,nrow*2000);
    
    cchi2->Divide(ncol,nrow);
    
    TLegend *l = new TLegend(-0.04,0,1,1);
    l->SetFillColor(0);
    l->SetFillStyle(0);
    l->SetLineColor(0);
    l->SetBorderSize(0);
    l->SetTextSize(0.07);
    if((histoName.Contains("MM02") || histoName.Contains("MM20")) && !histoName.Contains("NoCut") )
      l->SetHeader(Form("    %2.1f < #it{E} < %2.1f GeV, #it{n}^{w}_{cells} > 4",binE[iebin],binE[iebin+1]));
    else
      l->SetHeader(Form("    %2.1f < #it{E} < %2.1f GeV",binE[iebin],binE[iebin+1]));
    
    for(Int_t ism = firstSM; ism <= lastSM; ism++)
    {
      cchi2->cd(ism+1);
      
      gPad->SetLogy();
      gPad->SetGridy();
      
      if ( debug ) printf("iE %d ism %d\n",iebin,ism);
      if ( !hData[iebin][ism-firstSM] ) continue;
      
      for(Int_t iprod = 0; iprod < nProd; iprod++)
      {
        if ( debug ) 
          printf("\t Chi2 prod %d %p %p\n",iprod,hChi2[iprod][iebin][ism-firstSM],hDiff[iprod][iebin][ism-firstSM]);
        
        if(!hChi2[iprod][iebin][ism-firstSM]) continue;
        
        hChi2[iprod][iebin][ism-firstSM]->SetTitleOffset(1.8,"Y");        
        hChi2[iprod][iebin][ism-firstSM]->SetYTitle("(x_{simu}-x_{Data})^{2}/(#sigma_{x,simu}^{2}+#sigma_{x,data}^{2})");
        // hChi2[iprod][iebin][ism-firstSM]->SetAxisRange(0,0.8,"X");
        
        //printf("iebin %d, ism %d 0 %d %p\n",iebin,ism,0,hData[iebin][ism-firstSM]);
        
        if(iprod==1) hChi2[iprod][iebin][ism-firstSM]->Draw("H");
        else         hChi2[iprod][iebin][ism-firstSM]->Draw("H same");
        
        //hChi2[iprod][iebin][ism-firstSM]->SetMaximum(100);
        hChi2[iprod][iebin][ism-firstSM]->SetMinimum(1e-2);
        if(hChi2[1][iebin][ism-firstSM] && hChi2[iprod][iebin][ism-firstSM]->GetMaximum() > hChi2[1][iebin][ism-firstSM]->GetMaximum())
          hChi2[1][iebin][ism-firstSM]->SetMaximum(hChi2[iprod][iebin][ism-firstSM]->GetMaximum()*1.2);
        
        if(ism==0)
          l->AddEntry(hChi2[iprod][iebin][ism-firstSM],Form("%f",prodParam[iprod]),"PL");
      } // iprod
    } // param
    
    cchi2->cd(ncol*nrow);
    l->Draw();
    
    fileName = Form("figures/%s/Comparison_Chi2_Ebin%d_%s",
                    histoName.Data(),iebin,inputFileName.Data());
    fileName+=fileFormat;
    cchi2->Print(fileName);
    
    // Difference
    //
    TCanvas * cdiff = new TCanvas(Form("cdiff_ebin%d_%s",
                                       iebin,histoName.Data()),
                                  Form("MC-Data %2.1f < E < %2.1f, %s",
                                       binE[iebin],binE[iebin+1],
                                       histoName.Data()),
                                  ncol*2000,nrow*2000);
    cdiff->Divide(ncol,nrow);
    
    //        TLegend *l = new TLegend(-0.04,0,1,1);
    //        l->SetFillColor(0);
    //        l->SetFillStyle(0);
    //        l->SetLineColor(0);
    //        l->SetBorderSize(0);
    //        l->SetTextSize(0.07);
    //        l->SetHeader(Form("    %2.1f < #it{E} < %2.1f GeV",binE[iebin],binE[iebin+1]));
    for(Int_t ism = firstSM; ism <= lastSM; ism++)
    {
      cdiff->cd(ism+1);
      
      //gPad->SetLogy();
      gPad->SetGridy();
      
      if( debug ) printf("iE %d ism %d\n",iebin,ism);
      if ( !hData[iebin][ism-firstSM] ) continue;
      
      for(Int_t iprod = 0; iprod < nProd; iprod++)
      {
        if( debug )
          printf("\t Diff A prod %d %p\n",iprod,hDiff[iprod][iebin][ism-firstSM]);
        
        if(!hDiff[iprod][iebin][ism-firstSM]) continue;
        
        hDiff[iprod][iebin][ism-firstSM]->SetTitleOffset(1.8,"Y");
        //h[iebin][0]->SetLineColor(1);
        //h[iebin][0]->SetAxisRange(0.1,2.5,"X");
        
        hDiff[iprod][iebin][ism-firstSM]->SetYTitle("x_{simu}-x_{Data}");
        
        // hDiff[iprod][iebin][ism-firstSM]->SetAxisRange(0,0.8,"X");
        
        //printf("iebin %d, ism %d 0 %d %p\n",iebin,ism,0,hData[iebin][ism-firstSM]);
        
        if(iprod==1) hDiff[iprod][iebin][ism-firstSM]->Draw("H");
        else         hDiff[iprod][iebin][ism-firstSM]->Draw("H same");
        
        //hDiff[iprod][iebin][ism-firstSM]->SetMaximum(100);
        //hDiff[iprod][iebin][ism-firstSM]->SetMinimum(1e-2);
        if(hDiff[1][iebin][ism-firstSM] && hDiff[iprod][iebin][ism-firstSM]->GetMaximum() > hDiff[1][iebin][ism-firstSM]->GetMaximum())
          hDiff[1][iebin][ism-firstSM]->SetMaximum(hDiff[iprod][iebin][ism-firstSM]->GetMaximum()*1.2);
        
        //                if(ism==0)
        //                    l->AddEntry(hDiff[iprod][iebin][ism-firstSM],Form("%s",prodLeg[iprod].Data()),"PL");
      } // iprod
    } // sm
    
    cdiff->cd(ncol*nrow);
    l->Draw();
    
    fileName = Form("figures/%s/Comparison_DiffMCData_Ebin%d_%s",
                    histoName.Data(),iebin,inputFileName.Data());
    fileName+=fileFormat;
    cdiff->Print(fileName);
  } // ebin
  
  ////////////////////////////////
  // For each MC different parameter
  // assing a Chi2 value, do plots with
  // this distribution and the mininimum
  // of this distribution
  ////////////////////////////////
  if ( debug ) printf ("Compare param\n");
  
  const Int_t nX = hData[0][0]->GetNbinsX();
  const Int_t nXbins = nX;
  //Double_t * xaxis = hData[0][0]->GetXaxis()->GetXbins()->GetArray();
  Double_t xaxis[nXbins];
  Double_t xaxisErr[nXbins];
  Double_t yaxisErr[nXbins];
  Double_t yaxis[nXbins];
  Double_t yaxisA[nXbins];
  
  TGraphErrors* gChi2     [nEBins][nSM];
  TGraphErrors* gChi2Min  [nSM];
  TGraphErrors* gDiffMin  [nEBins][nSM];
  
  Double_t chi2axis       [nProd];
  Double_t chi2axisErr    [nProd];
  
  Double_t chi2axisMin    [nEBins];
  Double_t chi2axisMinErr [nEBins];
  
  // Get x array
  for(Int_t ix = 0; ix < nXbins; ix++)
  {
    yaxisErr[ix] = 0.05;
    xaxisErr[ix] = hData[0][0]->GetBinWidth(ix+1)/2;
    xaxis   [ix] = hData[0][0]->GetBinCenter(ix+1);
    //printf("bin %d, x %2.4f err %2.4f\n",ix,xaxis[ix],xaxisErr[ix]);
  }
  
  //
  // Fill y array, get min MC-Data of all cases or Chi2 in selected range
  //
  /// Per SM/Param
  for(Int_t ism = firstSM; ism <= lastSM; ism++)
  {
    for(Int_t iebin = 0; iebin < nEBins; iebin++)
    {
      for(Int_t ix = 0; ix < nXbins; ix++)
      {
        yaxis[ix] = -1;
        Double_t min = 1e6;
        
        for(Int_t iprod = 0; iprod < nProd; iprod++)
        {
          //if( debug )
          //  printf("\t Diff B prod %d %p\n",iprod,hDiff[iprod][iebin][ism-firstSM]);
          
          if ( !hDiff[iprod][iebin][ism-firstSM] ) continue;
          
          Double_t content = TMath::Abs(hDiff[iprod][iebin][ism-firstSM]->GetBinContent(ix+1));
          if(min > content && iprod > 0)
          {
            min = content;
            yaxis[ix] = prodParam[iprod];
          }
        } // prod
          //if(ism==3 && iebin==3) printf("ie %d, ip %d, ix %d, x %2.2f, y %2.2f\n",iebin,ism,ix,xaxis[ix],yaxis[ix]);
      } // x
      
      gDiffMin[iebin][ism-firstSM] = new TGraphErrors(nXbins,xaxis,yaxis,xaxisErr,yaxisErr);
      
      //gDiffMin[iebin][ism-firstSM]->GetHistogram()->SetAxisRange(xmin,xmax,"X");
      if ( hDiff[1][iebin][ism-firstSM] )
        gDiffMin[iebin][ism-firstSM]->GetHistogram()->SetXTitle(hDiff[1][iebin][ism-firstSM]->GetXaxis()->GetTitle());
      
      //
      // Chi2
      //
      Double_t minChi2SM = 10000000;
      for(Int_t iprod = 1; iprod < nProd-1; iprod++)
      {
        if ( !hChi2[iprod][iebin][ism-firstSM]  ) continue;
        
        Int_t minChi2Bin = hChi2[iprod][iebin][ism-firstSM]->FindBin(xmin);
        Int_t maxChi2Bin = hChi2[iprod][iebin][ism-firstSM]->FindBin(xmax);
        Double_t error = 0;
        if ( maxChi2Bin-minChi2Bin > 0 )
        {
          chi2axis   [iprod-1] = hChi2[iprod][iebin][ism-firstSM]->
          IntegralAndError(minChi2Bin,maxChi2Bin,error)/(maxChi2Bin-minChi2Bin);
          chi2axisErr[iprod-1] = error / (maxChi2Bin-minChi2Bin);
          //            if(ism==3 && iebin > 2 && iebin < 5)
          //              printf("sm %d, ie %d, iprod %d, mu %1.2f, chi2 %2.2f, chierr %2.2f\n",
          //                     ism,iebin,iprod,prodParam[iprod-1],chi2axis[iprod],chi2axisErr[iprod]);
        }
        else
        {
          chi2axis   [iprod-1] = 0;
          chi2axisErr[iprod-1] = 0;
        }
        
        if ( minChi2SM > chi2axis[iprod-1] )
        {
          minChi2SM = chi2axis[iprod-1];
          chi2axisMin   [iebin] = prodParam[iprod-1];
          
          chi2axisMinErr[iebin] = 0.1;//0.05;
        }
      } // prod
      
      gChi2[iebin][ism-firstSM] = new TGraphErrors
      (nProd,prodParam.GetArray(),chi2axis,prodParamE,chi2axisErr);
      gChi2[iebin][ism-firstSM]->GetHistogram()->SetXTitle("#mu_{1} (%)");
      
    } // energy bin
    
    gChi2Min[ism-firstSM] = new TGraphErrors(nEBins,binE.GetArray(),chi2axisMin,binEErr,chi2axisMinErr);
    gChi2Min[ism-firstSM]->GetHistogram()->SetXTitle("E_{cluster} (GeV)");
    
  }// sm
  
  ////////////////////////////////
  // Plot difference graphs
  ////////////////////////////////  
  if ( debug )
    printf("Plot Diff Graphs!\n");
  
  // Diff
  //
  TCanvas * cDiffgraph = new TCanvas(Form("cDiff_graph_%s",
                                          histoName.Data()),
                                     Form("Diff, %s",
                                          histoName.Data()),
                                     ncol*2000,nrow*2000);
  cDiffgraph->Divide(ncol,nrow);
  
  TLegend *lg = new TLegend(-0.04,0,1,1);
  lg->SetFillColor(0);
  lg->SetFillStyle(0);
  lg->SetLineColor(0);
  lg->SetBorderSize(0);
  lg->SetTextSize(0.07);
  //l->SetHeader(Form("    %2.1f < #it{E} < %2.1f GeV",binE[iebin],binE[iebin+1]));
  for(Int_t ism = firstSM; ism <= lastSM; ism++)
  {
    cDiffgraph->cd(ism+1);
    
    //gPad->SetLogy();
    gPad->SetGridy();
    //gPad->SetGridx();
    
    //printf("iE %d ism %d\n",iebin,ism);
    if(!gDiffMin[firstEbin][ism-firstSM]) continue;
    
    for(Int_t iebin = firstEbin; iebin < nEBins; iebin++)
    {
      //printf("\t prod %d %p\n",iprod,gDiffMin[iebin][ism-firstSM]);
      if(!gDiffMin[iebin][ism-firstSM]) continue;
      
      gDiffMin[iebin][ism-firstSM]->SetMarkerColor(color[iebin]);
      gDiffMin[iebin][ism-firstSM]->SetLineColor  (color[iebin]);
      gDiffMin[iebin][ism-firstSM]->SetMarkerStyle(24);
      gDiffMin[iebin][ism-firstSM]->SetMarkerSize(3);
      
      gDiffMin[iebin][ism-firstSM]->GetHistogram()->SetTitleOffset(1.8,"Y");
      //h[iebin][0]->SetLineColor(1);
      //h[iebin][0]->SetAxisRange(0.1,2.5,"X");
      
      //hDiff[iprod][iebin][ism-firstSM]->SetYTitle("(x_{simu}-x_{Data})^{2}/(#sigma_{x,simu}^{2}+#sigma_{x,data}^{2})");
      gDiffMin[iebin][ism-firstSM]->GetHistogram()->SetYTitle("#mu_{1} %");
      gDiffMin[iebin][ism-firstSM]->SetTitle(Form("SM %d",ism));
      
      // gDiffMin[iebin][ism-firstSM]->SetAxisRange(0,0.8,"X");
      
      //printf("iebin %d, ism %d 0 %d %p\n",iebin,ism,0,hData[iebin][ism-firstSM]);
      
      if(iebin==firstEbin)
        gDiffMin[iebin][ism-firstSM]->Draw("APL");
      else
        gDiffMin[iebin][ism-firstSM]->Draw("PL");
      
      gDiffMin[iebin][ism-firstSM]->SetMaximum(2);
      gDiffMin[iebin][ism-firstSM]->SetMinimum(-0.1);
      //        if(hDiff[iprod][iebin][ism-firstSM]->GetMaximum() > hData[iebin][ism-firstSM]->GetMaximum())
      //          hData[iebin][ism-firstSM]->SetMaximum(hDiff[iprod][iebin][ism-firstSM]->GetMaximum()*1.2);
      
      if(ism==0)
      {
        lg->AddEntry(gDiffMin[iebin][ism-firstSM],
                     Form("%2.1f < #it{E} < %2.1f GeV", binE[iebin], binE[iebin+1]), "PL");
      }
      
    } // iebin
    
    //printf("\t end\n");
  } // param
  
  cDiffgraph->cd(ncol*nrow);
  lg->Draw();
  
  fileName = Form("figures/%s/Comparison_MinDiff_PerSM_Ebin_%s",
                  histoName.Data(),inputFileName.Data());
  fileName+=fileFormat;
  cDiffgraph->Print(fileName);
  
  // Min diff
  //
  TCanvas * cDiffgraph2 = new TCanvas(Form("cMinDiff_graph2_%s",
                                           histoName.Data()),
                                      Form("MinDiff, %s",
                                           histoName.Data()),
                                      3*2000,2*2000);
  cDiffgraph2->Divide(3,2);
  
  TLegend *lg2 = new TLegend(-0.04,0,1,1);
  lg2->SetFillColor(0);
  lg2->SetFillStyle(0);
  lg2->SetLineColor(0);
  lg2->SetBorderSize(0);
  lg2->SetTextSize(0.07);
  for(Int_t iebin = firstEbin; iebin < nEBins; iebin++)
  {
    cDiffgraph2->cd(iebin-2);
    
    //gPad->SetLogy();
    gPad->SetGridy();
    
    //printf("iE %d ism %d\n",iebin,ism);
    if(!gDiffMin[iebin][0]) continue;
    
    for(Int_t ism = firstSM; ism <= lastSM; ism++)
    {
      //printf("\t prod %d %p\n",iprod,gDiffMin[iebin][ism-firstSM]);
      if(!gDiffMin[iebin][ism-firstSM]) continue;
      
      gDiffMin[iebin][ism-firstSM]->SetMarkerColor(color[ism-firstSM]);
      gDiffMin[iebin][ism-firstSM]->SetLineColor  (color[ism-firstSM]);
      gDiffMin[iebin][ism-firstSM]->SetMarkerStyle(24);
      gDiffMin[iebin][ism-firstSM]->SetMarkerSize(3);
      
      gDiffMin[iebin][ism-firstSM]->GetHistogram()->SetTitleOffset(1.8,"Y");
      //h[iebin][0]->SetLineColor(1);
      //h[iebin][0]->SetAxisRange(0.1,2.5,"X");
      
      //hDiff[iprod][iebin][ism-firstSM]->SetYTitle("(x_{simu}-x_{Data})^{2}/(#sigma_{x,simu}^{2}+#sigma_{x,data}^{2})");
      gDiffMin[iebin][ism-firstSM]->GetHistogram()->SetYTitle("#mu_{1} %");
      
      gDiffMin[iebin][ism-firstSM]->SetTitle(Form("%2.1f < #it{E} < %2.1f GeV",binE[iebin],binE[iebin+1]));
      // gDiffMin[iebin][ism-firstSM]->SetAxisRange(0,0.8,"X");
      
      //printf("iebin %d, ism %d 0 %d %p\n",iebin,ism,0,hData[iebin][ism-firstSM]);
      
      if(ism==0)
        gDiffMin[iebin][ism-firstSM]->Draw("APL");
      else
        gDiffMin[iebin][ism-firstSM]->Draw("PL");
      
      gDiffMin[iebin][ism-firstSM]->SetMaximum(2);
      gDiffMin[iebin][ism-firstSM]->SetMinimum(-0.1);
      //        if(hDiff[iprod][iebin][ism-firstSM]->GetMaximum() > hData[iebin][ism-firstSM]->GetMaximum())
      //          hData[iebin][ism-firstSM]->SetMaximum(hDiff[iprod][iebin][ism-firstSM]->GetMaximum()*1.2);
      
      if(iebin==firstEbin)
      {
        lg2->AddEntry(gDiffMin[iebin][ism-firstSM],(Form("SM %d",ism)),"PL");
      }
      
    } // iebin
    
    //printf("\t end\n");
  } // param
  
  cDiffgraph2->cd(6);
  lg2->Draw();
  
  fileName = Form("figures/%s/Comparison_MinDiff_PerEBin_SM_%s",
                  histoName.Data(),inputFileName.Data());
  fileName+=fileFormat;
  cDiffgraph2->Print(fileName);
  
  ////////////////////////////////
  // Plot  Chi2 graphs
  ////////////////////////////////  
  if ( debug )
    printf("Plot Chi2 Graphs!\n");
  
  // Chi2
  //
  TCanvas * gChi2graph = new TCanvas(Form("cChi2_graph_%s",
                                          histoName.Data()),
                                     Form("Chi2, %s",
                                          histoName.Data()),
                                     ncol*2000,nrow*2000);
  gChi2graph->Divide(ncol,nrow);
  
  TLegend *lgchi2 = new TLegend(-0.04,0,1,1);
  lgchi2->SetFillColor(0);
  lgchi2->SetFillStyle(0);
  lgchi2->SetLineColor(0);
  lgchi2->SetBorderSize(0);
  lgchi2->SetTextSize(0.07);
  //l->SetHeader(Form("    %2.1f < #it{E} < %2.1f GeV",binE[iebin],binE[iebin+1]));
  for(Int_t ism = firstSM; ism <= lastSM; ism++)
  {
    gChi2graph->cd(ism+1);
    
    gPad->SetLogy();
    gPad->SetGridy();
    gPad->SetGridx();
    
    //printf("iE %d ism %d\n",iebin,ism);
    if(!gChi2[firstEbin][ism-firstSM]) continue;
    
    for(Int_t iebin = firstEbin; iebin < nEBins; iebin++)
    {
      //printf("\t prod %d %p\n",iprod,gChi2[iebin][ism-firstSM]);
      if(!gChi2[iebin][ism-firstSM]) continue;
      
      gChi2[iebin][ism-firstSM]->SetMarkerColor(color[iebin]);
      gChi2[iebin][ism-firstSM]->SetLineColor  (color[iebin]);
      gChi2[iebin][ism-firstSM]->SetMarkerStyle(24);
      gChi2[iebin][ism-firstSM]->SetMarkerSize(3);
      
      gChi2[iebin][ism-firstSM]->GetHistogram()->SetTitleOffset(1.6,"Y");
      //h[iebin][0]->SetLineColor(1);
      //h[iebin][0]->SetAxisRange(0.1,2.5,"X");
      
      //hDiff[iprod][iebin][ism-firstSM]->SetYTitle("(x_{simu}-x_{Data})^{2}/(#sigma_{x,simu}^{2}+#sigma_{x,data}^{2})");
      gChi2[iebin][ism-firstSM]->GetHistogram()->SetYTitle("#chi^{2} / #nu");
      gChi2[iebin][ism-firstSM]->SetTitle(Form("SM %d",ism));
      
      // gChi2[iebin][ism-firstSM]->SetAxisRange(0,0.8,"X");
      
      //printf("iebin %d, ism %d 0 %d %p\n",iebin,ism,0,hData[iebin][ism-firstSM]);
      
      if(iebin==firstEbin)
        gChi2[iebin][ism-firstSM]->Draw("APL");
      else
        gChi2[iebin][ism-firstSM]->Draw("PL");
      
      gChi2[iebin][ism-firstSM]->SetMaximum(200);
      gChi2[iebin][ism-firstSM]->SetMinimum(0.25);
      //              if(hDiff[iprod][iebin][ism-firstSM]->GetMaximum() > hData[iebin][ism-firstSM]->GetMaximum())
      //              hData[iebin][ism-firstSM]->SetMaximum(hDiff[iprod][iebin][ism-firstSM]->GetMaximum()*1.2);
      
      if(ism==0)
      {
        lgchi2->AddEntry(gChi2[iebin][ism-firstSM],
                         Form("%2.1f < #it{E} < %2.1f GeV",binE[iebin],binE[iebin+1]),"PL");
      }
      
    } // iebin
    
    //printf("\t end\n");
  } // param
  
  gChi2graph->cd(ncol*nrow);
  lgchi2->Draw();
  
  fileName = Form("figures/%s/Comparison_Chi2_PerSM_Ebin_%s",
                  histoName.Data(),inputFileName.Data());
  fileName+=fileFormat;
  gChi2graph->Print(fileName);
  
  // Min Chi2
  //
  TCanvas * gChi2Mingraph = new TCanvas(Form("cChi2Min_graph_%s",
                                             histoName.Data()),
                                        Form("Diff, %s",
                                             histoName.Data()),
                                        ncol*2000,nrow*2000);
  gChi2Mingraph->Divide(ncol,nrow);
  
  TLegend *lgChi2Min = new TLegend(-0.04,0,1,1);
  lgChi2Min->SetFillColor(0);
  lgChi2Min->SetFillStyle(0);
  lgChi2Min->SetLineColor(0);
  lgChi2Min->SetBorderSize(0);
  lgChi2Min->SetTextSize(0.07);
  
  for(Int_t ism = firstSM; ism <= lastSM; ism++)
  {
    gChi2Mingraph->cd(ism+1);
    
    //gPad->SetLogy();
    //gPad->SetGridy();
    
    //printf("iE %d ism %d\n",iebin,ism);
    if(!gChi2Min[ism-firstSM]) continue;
    
    //printf("\t prod %d %p\n",iprod,gChi2Min[ism-firstSM]);
    if(!gChi2Min[ism-firstSM]) continue;
    
    gChi2Min[ism-firstSM]->SetMarkerColor(color[ism-firstSM]);
    gChi2Min[ism-firstSM]->SetLineColor  (color[ism-firstSM]);
    gChi2Min[ism-firstSM]->SetMarkerStyle(24);
    gChi2Min[ism-firstSM]->SetMarkerSize(3);
    
    gChi2Min[ism-firstSM]->GetHistogram()->SetTitleOffset(1.8,"Y");
    //hData->SetLineColor(1);
    //hData->SetAxisRange(0.1,2.5,"X");
    
    //hDiff[iprod][ism-firstSM]->SetYTitle("(x_{simu}-x_{Data})^{2}/(#sigma_{x,simu}^{2}+#sigma_{x,data}^{2})");
    gChi2Min[ism-firstSM]->GetHistogram()->SetYTitle("#mu_{1} (%)");
    gChi2Min[ism-firstSM]->SetTitle(Form("SM %d",ism));
    
    //      if(firstEbin == 3) gChi2Min[ism-firstSM]->GetHistogram()->SetAxisRange(7,20,"X");
    //      else               gChi2Min[ism-firstSM]->GetHistogram()->SetAxisRange(15,60,"X");
    
    //printf("iebin %d, ism %d 0 %d %p\n",iebin,ism,0,hData[ism-firstSM]);
    
    gChi2Min[ism-firstSM]->Draw("APL");
    
    gChi2Min[ism-firstSM]->SetMaximum(1.7);
    gChi2Min[ism-firstSM]->SetMinimum(-0.1);
    
  } // param
  
  gChi2Mingraph->cd(ncol*nrow);
  lgChi2Min->Draw();
  
  fileName = Form("figures/%s/Comparison_Chi2Min_PerSM_Ebin_%s",
                  histoName.Data(),inputFileName.Data());
  fileName+=fileFormat;
  gChi2Mingraph->Print(fileName);
  
  // Write to file
  TString inputFileName2 = inputFileName;
  inputFileName2.ReplaceAll("ForChi2","");
  
  TFile* fout = new TFile(Form("figures/%s/Chi2Distributions_%s_PerSM_Ebin.root",
                               histoName.Data(),inputFileName2.Data()),"recreate");
  for(Int_t ism = firstSM; ism <= lastSM; ism++)
  {
    gChi2Min[ism-firstSM]->SetName(Form("Chi2Min_SM%d",ism));
    gChi2Min[ism-firstSM]->Write();
    for(Int_t iebin = firstEbin; iebin < nEBins; iebin++)
    {
      gChi2[iebin][ism-firstSM]->SetName(Form("Chi2_SM%d_EBin%d",ism,iebin));
      gChi2[iebin][ism-firstSM]->SetTitle(Form("%2.1f < #it{E} < %2.1f GeV",binE[iebin],binE[iebin+1]));
      gChi2[iebin][ism-firstSM]->Write();
    }
  }
  fout->Close();                      
    
}

//------------------------------------------------------------------------------
/// This method calculates the difference and Chi2 of 2 histograms
///
/// \param h     : Input MC histogram
/// \param hD    : Input Data reference histogram
/// \param hChi2 : Output Chi2 histogram
/// \param hDiff : Output Difference
/// \param debug : Print debugging info
//------------------------------------------------------------------------------
void Chi2( TH1D* h, TH1D* hD, TH1D*& hChi2, TH1D*&  hDiff, Bool_t debug)
{
  if ( !h || !hD ) return ;
  
  hChi2 = (TH1D*) h->Clone(Form("%s_Chi2",h->GetName()));
  hDiff = (TH1D*) h->Clone(Form("%s_Diff",h->GetName()));
  
  //printf("%p %p, N bins %d\n",h,hD,h->GetNbinsX());
  
  for(Int_t ibin = 1; ibin < h->GetNbinsX(); ibin++)
  {
    //printf("bin %d\n",ibin);
    Double_t content  = h ->GetBinContent(ibin);
    Double_t contentD = hD->GetBinContent(ibin);
    Double_t value2   = content-contentD;
    value2*=value2;
    
    Double_t econtent2  = h->GetBinError(ibin);
    Double_t econtent2D = hD->GetBinError(ibin);
    
    Double_t esum = econtent2D+econtent2;
    esum *=esum;
    
    econtent2 *=econtent2;
    econtent2D*=econtent2D;
    
    Double_t esum2 = econtent2D+econtent2;
    
    if ( esum2 > 1e-19 )
    {
      if(debug) printf("ibin %d, value2 %e, esum %e\n",ibin,value2,esum2);
      
      hDiff->SetBinContent(ibin,content-contentD);
      hDiff->SetBinError  (ibin,h->GetBinError(ibin)+hD->GetBinError(ibin));
      
      value2/=esum2;
      esum/=esum2;
      
      if ( debug )
        printf("\t ibin %d, x %2.2f, vS %2.4f, vD %2.4f, |vS-vD| %2.4e |vS-vD|^2 %2.4e, Err|vS-vD|^2 %2.4e\n",
               ibin, hD->GetBinCenter(ibin), content, contentD, TMath::Abs(content-contentD),value2,esum);
      
      hChi2->SetBinContent(ibin,value2); 
      hChi2->SetBinError(ibin,esum); 
    }
    else
    {
      hChi2->SetBinContent(ibin,0);
      hChi2->SetBinError(ibin,0);
      hDiff->SetBinContent(ibin,0);
      hDiff->SetBinError(ibin,0);
    }
    
  }
  if ( debug ) 
    printf("\t func %p %p\n",hChi2,hDiff);
}


