///
/// \file CompareShowerShapeLongPerAnaPerSM.C
/// \ingroup CaloTrackCorrMacrosPlotting
/// \brief Compare a shower shape long axis between different productions 
///
/// Compare different shape long axis in different data and MC productions
/// depending on the SM per cluster energy bis and different cluster selections.
/// Treatment of the output of the class AliAnaPhoton and AliAnaParticleIsolation
/// Input are TH2 histograms where x=energy, y= shower shape long axis. There is a histogram per each SM.
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
#include <TGraphErrors.h>
#include "PlotUtils.C"
#endif

//----------------------------------------------------------------------------
// Global variables
//----------------------------------------------------------------------------

Int_t color    [] = {1,4,2,kYellow-2,8,kCyan,kYellow-6,kCyan,kOrange+2,kViolet,
  kOrange-2,4,2,6,8,9,kYellow-6,kCyan,kOrange+2,kViolet,kOrange-2};
Int_t lineStyle[] = {1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,        
  2,2,2,2,2,2,2,2,2,2,2,2};
Int_t marker   [] = {20,20,20,21,24,24,24,24,24,24,24};

TString fileFormat = ".eps";

//*****
// Set the energy bins and productions to analyze
//*****
const Int_t nEBins = 5; 
const Int_t nProd  = 3;
TFile * file[nProd];

//////////////////////////////////////////////////////////
// Low JJ
Double_t binE   [] = { 8,10,12,14,16,18,20,25};

TString titleName = "LHC11cd_EMC7_MCGJ-JJLow";
TString filePath[] = 
{
  "data/module/TCardChannel3/LHC11cd_EMC7"
  , "simu/module/pp_7TeV_JJ_Dec_GJ/TCardChannel_Mimic0_Scaled2_v3/ScaledMerged"
  , "simu/module/pp_7TeV_JJ_Dec_GJ/TCardChannel_Mimic10c_EcellCut_Scaled2_v3/ScaledMerged"
}  ;

TString dataType[] = 
{
  "Data, pp@ 7 TeV, LHC11c+d"
  , "MC default, GJ+JJ_{p^{EMCal}_{T,#gamma} > 3.5 GeV/#it{c}}"
  , "MC   xTalk, GJ+JJ_{p^{EMCal}_{T,#gamma} > 3.5 GeV/#it{c}}"
}  ;
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//  // High JJ
//  Double_t binE   [] = { 18,20, 25, 30, 35,40,50};
//  
//  TString titleName = "LHC11cd_EMC7_MCGJ-JJHigh";
//  TString filePath[] = 
//  {
//  "data/module/TCardChannel3/LHC11cd_EMC7"
//  , "simu/module/pp_7TeV_JJ_Dec_High_GJ/TCardChannel_Mimic0_Scaled2_v3/ScaledMerged"
//  , "simu/module/pp_7TeV_JJ_Dec_High_GJ/TCardChannel_Mimic10c_EcellCut_Scaled2_v3/ScaledMerged"
//  }  ;
//  
//  TString dataType[] = 
//  {
//    "Data, pp@ 7 TeV, LHC11c+d"
//    , "MC default, GJ+JJ_{p^{EMCal}_{T,#gamma} > 7 GeV/#it{c}}"
//    , "MC   xTalk, GJ+JJ_{p^{EMCal}_{T,#gamma} > 7 GeV/#it{c}}"
//  }  ;
//////////////////////////////////////////////////////////
//----------------------------------------------------------------------------


///---------------------------------------------------------------------------
/// Main method called by CompareShowerShapeLongPerAnaPerSM()
/// \param iana   : input analysis type: 0-Not isolated; 1-Isolated; 2-Inclusive from AliAnaPhoton
/// \param cut    : plot shower shape long axis 0-without n cell cut, 1- with n cell cut 
/// \param firstSM: first SM number to be inspected
/// \param lastSM : last  SM number to be inspected
/// \param plotRat: make ratio plots.
/// \param bAllSM : just recover the histogram filled for all SM if it exists and plot only those.
/// \param debug  : Bool to activate printfs.
///---------------------------------------------------------------------------
void DoIt ( Int_t  iana    = 0, 
            Bool_t cut     = kFALSE, 
            Int_t  firstSM = 0,
            Int_t  lastSM  = 9, 
            Bool_t plotRat = kFALSE,
            Bool_t bAllSM  = kFALSE,
            Bool_t debug   = kFALSE
          )
{    
  Double_t lowE[nEBins];
  Double_t higE[nEBins];
  
  // SM
  const Int_t nSM = lastSM-firstSM+2; // Last entry will be sum of all.
  
  if ( debug )
    printf("Execute ana %d, cut %d, allSM %d, first SM %d, last SM %d, total SM+1 %d\n",
           iana,cut,bAllSM,firstSM,lastSM, nSM);
  
  // histograms
  //printf("N SM = %d\n",nSM);
  TH2F* h2[nProd][nSM];
  TH1D* h [nProd][nEBins][nSM];
  TH1D* hA[nProd][nEBins];
  
  // Merge good/bad SMs
  TH1D* hGroup[nProd][nEBins][3];
  
  Int_t    rebin = 2;
  Double_t xmin = 0.1;
  Double_t xmax = 1.1;
  
  TString isoName [] = {"NoIso","Iso","Photon"};
  TString isoTitle[] = {"NOT Isolated"  ,"Isolated","Inclusive"};
  
  TString cutName [] = {"","_NCellCut"};
  TString cutTitle[] = {"no #it{n}^{w} cut"  ,"#it{n}^{w} > 4"};
 
//  for(Int_t iprod = 0; iprod < nProd; iprod++)
//  {    
//    for(Int_t ism = 0; ism < nSM; ism++)
//    {
//      h2[iprod][ism] = 0;
//    }
//  }
  
  for(Int_t iprod = 0; iprod < nProd; iprod++)
  {
    if ( !file[iprod] ) continue;
    
    for(Int_t ism = 0; ism < nSM-1; ism++)
    {
      h2[iprod][ism] = 0;
      
      if ( iana == 2 && cut == 1 ) 
        continue;
      
      if(iana < 2)
      {
        h2[iprod][ism] = 
        (TH2F*) file[iprod]->Get(Form("AnaIsolPhoton_hPtLambda0%s_SM%d_%s",
                               cutName[cut].Data(),ism+firstSM,isoName[iana].Data()));
      }
      else
      {
        h2[iprod][ism] = 
        (TH2F*) file[iprod]->Get(Form("AnaPhoton_hLam0_SM%d",ism+firstSM));
      }
      
      if ( !h2[iprod][ism] )
      {
        printf("Histogram for ana %d, cut %d, sm %d not found\n",iana,cut,ism);
        continue;
      }
      
//      if ( debug )
//      {
//        printf("iana %d, sm %d, prod %d, %p\n",iana, ism, iprod, h2[iprod][ism]);
//        printf("\t %s Entries %e Integral %e\n",
//               h2[iprod][ism]->GetName(),h2[iprod][ism]->GetEntries(),h2[iprod][ism]->Integral());
//      }      
      
      // Merge all SM
      if ( ism == 0 )
      {
        h2[iprod][nSM-1] = (TH2F*) h2[iprod][ism]->Clone(Form("AllSM_Prod%d_Ana%d_Cut%d",iprod,iana,cut));
        //printf(" Clone ism %d %p\n",ism, h2[iprod][nSM]);
      }
      else
      {
        h2[iprod][nSM-1] -> Add (h2[iprod][ism]);
        //printf(" Add ism %d to %p\n",ism, h2[iprod][nSM]);
      }
    } // ism
    
    // Get the histo for all SM from file
    // or merge all SM, done before
    if ( bAllSM && cut == 0 )
    {
      if ( debug )
      printf("Recover histograms for ALL SM!!!!\n");
      
      if(iana < 2)
        h2[iprod][nSM-1] = 
        (TH2F*) file[iprod]->Get(Form("AnaIsolPhoton_hPtLambda0%s",
                                      isoName[iana].Data())); // hELambda0
      else
        h2[iprod][nSM-1] = (TH2F*) file[iprod]->Get("AnaPhoton_hLam0Pt"); // ("AnaPhoton_hLam0E");
      
      if ( !h2[iprod][nSM-1] )
      {
        printf("Histogram for ana %d, cut %d, All SM not found\n",iana,cut);
      }
 
    } // Recover histo all SM from file
  
  } // prod loop
  
  // Apply Sumw2
  //
  for(Int_t iprod = 0; iprod < nProd; iprod++)
  {
    for(Int_t ism = 0; ism < nSM; ism++)
    {
      if ( !h2[iprod][ism] )
      {
        printf("Sumw2: Histogram for ana %d, cut %d, sm %d not found\n",iana,cut,ism);
        continue;
      }
      
      if ( filePath[iprod].Contains("data") || 
          (filePath[iprod].Contains("simu") && filePath[iprod].Contains("MB"))) 
        h2[iprod][ism]->Sumw2();
    }
  }
  
  //
  // Project
  //
  
  Double_t width = 0;
  for(Int_t ism = 0; ism < nSM; ism++)
  {
    for(Int_t iprod = 0; iprod < nProd; iprod++)
    {      
      for(Int_t iebin = 0; iebin < nEBins; iebin++)
      {        
        if ( !h2[iprod][ism] )
        {
          printf("Projection Histogram for iprod %d, ana %d, cut %d, sm %d not found\n",
                 iprod, iana, cut, ism);
          continue;
        }
        
        Int_t ebinMin = h2[iprod][ism]->GetXaxis()->FindBin(binE[iebin  ]);
        Int_t ebinMax = h2[iprod][ism]->GetXaxis()->FindBin(binE[iebin+1])-1;
        //if(ebinMin==ebinMax) ebinMax = 1000;
        lowE[iebin]  = h2[iprod][ism]->GetXaxis()->GetBinLowEdge(ebinMin);
        width        = h2[iprod][ism]->GetXaxis()->GetBinWidth  (ebinMax);
        higE[iebin]  = h2[iprod][ism]->GetXaxis()->GetBinLowEdge(ebinMax)+width;
                
        h[iprod][iebin][ism] = 
        (TH1D*) h2[iprod][ism]->ProjectionY(Form("Histo%d_BinE%d_SM%d",iprod,iebin,ism),ebinMin,ebinMax);
        
//        if ( debug )
//        {
//          printf("project ism %d, iprod %d, ie %d [%2.1f,%2.1f] GeV\n",
//                 ism, iprod, iebin,lowE[iebin],higE[iebin]);
//          printf("\t %p \n",h[iprod][iebin][ism]);
//        }
        
        if ( !h[iprod][iebin][ism] ) continue;
        
        h[iprod][iebin][ism]->SetLineColor(color[iprod]);
        
        h[iprod][iebin][ism]->SetLineWidth(2);
        
        h[iprod][iebin][ism]->SetLineStyle(lineStyle[iprod]);
        
        h[iprod][iebin][ism]->SetMarkerStyle(marker[iprod]);
        
        h[iprod][iebin][ism]->SetMarkerColor(color[iprod]);
        
        h[iprod][iebin][ism]->SetMarkerSize(0.5);
        
        h[iprod][iebin][ism]->SetXTitle("#sigma_{long}^{2}");
        
        if(rebin > 1) 
          h[iprod][iebin][ism]->Rebin(rebin);
        
        h[iprod][iebin][ism]->SetAxisRange(xmin,xmax,"X");
        
        // Normalize to integral above 0.6 (n cell cut)
        if ( cut )
        {
          Double_t integral = h[iprod][iebin][ism]->Integral(h[iprod][iebin][ism]->FindBin(0.6),
                                                             h[iprod][iebin][ism]->FindBin(2));
          if(integral > 1e-9)
          {
            h[iprod][iebin][ism]->Scale(1./integral);
            h[iprod][iebin][ism]->SetYTitle("Norm. integral 0.6<#lambda^{2}_{0}<2 ");
          }
          else
            h[iprod][iebin][ism] = 0;
        }
        // Normalize to max at 0.25
        else
        {
          Double_t scale = h[iprod][iebin][ism]->GetBinContent(h[iprod][iebin][ism]->FindBin(0.245));
          if(scale > 1e-9)
          {
            h[iprod][iebin][ism]->Scale(1./scale);
            h[iprod][iebin][ism]->SetYTitle("Norm. to max at 0.25");
          }
          else h[iprod][iebin][ism] = 0;
        }
      } // ie
    } // iprod
  } // ism
  
  
  //
  // Merge SM Groups
  //
  for(Int_t iprod = 0; iprod < nProd; iprod++)
  {      
    for(Int_t iebin = 0; iebin < nEBins; iebin++)
    {            
      if ( h[iprod][iebin][0] )
      {
        hGroup[iprod][iebin][0] = (TH1D*) h[iprod][iebin][0]->Clone(Form("SMGroup0_Prod%d_Ana%d_Cut%d",
                                                                         iprod,iana,cut));
        
        hGroup[iprod][iebin][0]->Add( h[iprod][iebin][4] );
        hGroup[iprod][iebin][0]->Add( h[iprod][iebin][5] );
        hGroup[iprod][iebin][0]->Add( h[iprod][iebin][6] );
        hGroup[iprod][iebin][0]->Add( h[iprod][iebin][8] );
        hGroup[iprod][iebin][0]->Add( h[iprod][iebin][9] );
        
        hGroup[iprod][iebin][0]->Scale(1./6.);
      }
      else
      {
        hGroup[iprod][iebin][0] = 0;
      }
      
      if(h[iprod][iebin][1]) 
      {
        hGroup[iprod][iebin][1] = (TH1D*) h[iprod][iebin][1]->Clone(Form("SMGroup1_Prod%d_Ana%d_Cut%d",
                                                                         iprod,iana,cut));
        
         hGroup[iprod][iebin][1]->Add( h[iprod][iebin][1] );
        
         hGroup[iprod][iebin][1]->Scale(1./2.);
      }
      else
      {
         hGroup[iprod][iebin][1] = 0;
      }
    
      if(h[iprod][iebin][3]) 
      {
        hGroup[iprod][iebin][2] = (TH1D*) h[iprod][iebin][3]->Clone(Form("SMGroup2_Prod%d_Ana%d_Cut%d",
                                                                         iprod,iana,cut));
        
        hGroup[iprod][iebin][2]->Add( h[iprod][iebin][7] );
        
        hGroup[iprod][iebin][2]->Scale(1./2.);
      }
      else
      {
        hGroup[iprod][iebin][2] = 0;
      }      
    } // ebins
  } // prod bins
  
    
  // 
  // Plots
  //
  
  if ( debug )
    printf("PLOT\n");

  TString fileName ;
  
  gStyle->SetPadRightMargin(0.01);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetTitleFontSize(0.06);
  //  gStyle->SetStatFontSize(0.06);
  gStyle->SetOptStat(0);
  
  // Plot in file different frames per E bin.
  //
  Int_t ncolE = 4;
  Int_t nrowE = 3;
  GetCanvasColRowNumber(nEBins,ncolE,nrowE); // PlotUtils.C
  
  TCanvas * cA = new TCanvas(Form("cAllSM_iana%d_cut%d",iana,cut),
                             Form("All SM, iana %d, cut %d",iana,cut),
                             ncolE*2000,nrowE*2000);
  
  cA->Divide(ncolE,nrowE);
  
  TLegend *l = 0;
  if(ncolE*nrowE!=nSM)
  {
    l= new TLegend(-0.04,0.,1,1);
    l->SetTextSize(0.06);
  }
  else
  {
    l= new TLegend(0.5,0.6,0.9,0.9);
    l->SetTextSize(0.035);
  }
  l->SetFillColor(0);
  l->SetFillStyle(0);
  l->SetLineColor(0);
  l->SetBorderSize(0);
  l->AddEntry("",Form(isoTitle[iana].Data(), cutTitle[cut].Data()),"");
  
  for(Int_t iebin = 0; iebin < nEBins; iebin++)
  {
    cA->cd(iebin+1);
    
    //printf("ebin %d %p\n",iebin,h[0][iebin][nSM-1]);
    if(!h[0][iebin][nSM-1]) continue;
    
    //gPad->SetLogy();
    
    gStyle->SetOptTitle(1);
    gStyle->SetPadTopMargin(0.08);
    
    h[0][iebin][nSM-1]->Draw("H");
    
    for(Int_t iprod = 0; iprod < nProd; iprod++)
    {
      //printf("iprod %d %p\n",iprod,h[iprod][iebin][nSM-1]);
      if(!h[iprod][iebin][nSM-1]) continue; 
      
      h[iprod][iebin][nSM-1]->SetTitleOffset(1.8,"Y");
      
      //printf("\t low %f high %f\n",lowE[iebin],higE[iebin]);
      h[iprod][iebin][nSM-1]->SetTitle(Form("%2.1f < #it{E}_{T}^{clus} < %2.1f GeV",lowE[iebin],higE[iebin]));
      
      h[iprod][iebin][nSM-1]->Draw("H same");
      
      if(h[iprod][iebin][nSM-1]->GetMaximum() > h[0][iebin][nSM-1]->GetMaximum())
        h[0][iebin][nSM-1]->SetMaximum(h[iprod][iebin][nSM-1]->GetMaximum()*1.2);
      
      if(iebin==0)
        l->AddEntry(h[iprod][iebin][nSM-1],Form("%s",dataType[iprod].Data()),"PL");      
    }  
  }
  
  cA->cd(ncolE*nrowE);
  l->Draw("same");
  
  fileName = Form("figures/Comparison_M02_AllSM_%s%s_%s",
                  isoName[iana].Data(),cutName[cut].Data(),titleName.Data());
  fileName+=fileFormat;
  cA->Print(fileName);
  
  if ( plotRat )
  {
    TCanvas * cAR = new TCanvas(Form("cR_AllSM_iana%d_cut%d",iana, cut),
                                Form("All SM, iana %d, cut %d",iana, cut),
                                ncolE*2000,nrowE*2000);
    
    cAR->Divide(ncolE,nrowE);
    
    for(Int_t iebin = 0; iebin < nEBins; iebin++)
    {
      cAR->cd(iebin+1);
      
      if(!h[0][iebin][nSM-1]) continue;
      
      //if(quantity.Contains("ECell"))gPad->SetLogy();
      
      gStyle->SetOptTitle(1);
      gStyle->SetPadTopMargin(0.08);
      
      h[0][iebin][nSM-1]->Draw("H");
      
      for(Int_t iprod = 1; iprod < nProd; iprod++)
      {
        if(!h[iprod][iebin][nSM-1]) continue; 
        
        TH1F* hRat = (TH1F*) h[iprod][iebin][nSM-1]->Clone(Form("%s_Ratio",h[iprod][iebin][nSM-1]->GetName()));
        hRat->Divide(h[0][iebin][nSM-1]);
        
        if(iprod==1) hRat->Draw("H");
        else         hRat->Draw("H same");
        
        hRat->SetYTitle("Ratio MC to Data");
        hRat->SetMaximum(1.2);
        hRat->SetMinimum(0);
      }
    }
    
    cAR->cd(ncolE*nrowE);
    l->Draw("same");
    
    fileName = Form("figures/Comparison_Ratio_M02_AllSM_%s%s_%s",
                    isoName[iana].Data(),cutName[cut].Data(),titleName.Data());
    fileName+=fileFormat;
    cAR->Print(fileName);
  }
  
  if ( bAllSM ) return;
  
  // One file per E bin
  // Each file all SM and the sum
  //
  
  Int_t ncolSM = 4;
  Int_t nrowSM = 3;
  GetCanvasColRowNumber(nSM,ncolSM,nrowSM); // PlotUtils.C
  
  for(Int_t iebin = 0; iebin < nEBins; iebin++)
  {
    gStyle->SetOptTitle(1);
    gStyle->SetPadTopMargin(0.07);
    
    TCanvas * c = new TCanvas(Form("c_ebin%d_iso%d_cut%d",
                                   iebin,iana,cut),
                              Form("iana %d, cut %d, %2.1f < E < %2.1f",
                                   iana, cut,lowE[iebin],higE[iebin]),
                              ncolSM*2000,nrowSM*2000);
    c->Divide(ncolSM,nrowSM);
    
    TLegend *lSM = 0;
    if(ncolSM*ncolSM!=nSM)
    {
      lSM= new TLegend(-0.04,0,1,1);
      lSM->SetTextSize(0.06);
    }
    else
    {
      lSM= new TLegend(0.5,0.6,0.9,0.9);
      lSM->SetTextSize(0.035);
    }
  
    lSM->SetFillColor(0);
    lSM->SetFillStyle(0);
    lSM->SetLineColor(0);
    lSM->SetBorderSize(0);
    lSM->SetHeader(Form("    %2.1f < #it{E} < %2.1f GeV", lowE[iebin],higE[iebin]));
    lSM->AddEntry("",Form(isoTitle[iana].Data(), cutTitle[cut].Data()),"");
    for(Int_t ism = 0; ism < nSM; ism++)
    {
      c->cd(ism+1);
      
      //gPad->SetLogy();
      
      //printf("iE %d ism %d\n",iebin,ism);
      if(!h[0][iebin][ism]) continue;
      
      for(Int_t iprod = 0; iprod < nProd; iprod++)
      {
        //printf("\t prod %d %p\n",iprod,h[iprod][iebin][ism]);
        if(!h[iprod][iebin][ism]) continue;
        
        h[iprod][iebin][ism]->SetTitle(Form("SM %d",ism));
        
        h[iprod][iebin][ism]->SetTitleOffset(1.8,"Y");
        //h[iebin][0]->SetLineColor(1);
        //h[iebin][0]->SetAxisRange(0.1,2.5,"X");
        
        if ( iprod == 0 ) h[iprod][iebin][ism]->Draw("H");
        else              h[iprod][iebin][ism]->Draw("H same");
        
        if(h[iprod][iebin][ism]->GetMaximum() > h[0][iebin][ism]->GetMaximum())
          h[0][iebin][ism]->SetMaximum(h[iprod][iebin][ism]->GetMaximum()*1.2);
        
        if(ism==0)
          lSM->AddEntry(h[iprod][iebin][ism],Form("%s",dataType[iprod].Data()),"PL");
      } // iprod

      //printf("\t end\n");
    } // param
    
    c->cd(ncolSM*nrowSM);
    lSM->Draw();
    
    fileName = Form("figures/Comparison_M02_Ebin%d_%s%s_%s",
                    iebin,isoName[iana].Data(),cutName[cut].Data(),titleName.Data());
    fileName+=fileFormat;
    c->Print(fileName);
        
    if ( !plotRat ) continue;
    
    TCanvas * cR = new TCanvas(Form("cR_ebin%d_iana%d_cut%d",
                                    iebin,iana,cut),
                               Form("Ratio iana %d, %2.1f < E < %2.1f",
                                    iana, lowE[iebin],higE[iebin]),
                               ncolSM*2000,nrowSM*2000);
    cR->Divide(ncolSM,nrowSM);
    
    for(Int_t ism = 0; ism < nSM; ism++)
    {
      cR->cd(ism+1);
      
      //gPad->SetLogy();
      
      //printf("iE %d ism %d\n",iebin,ism);
      if(!h[1][iebin][ism]) continue;
      
      for(Int_t iprod = 1; iprod < nProd; iprod++)
      {
        //printf("\t prod %d %p\n",iprod,h[iprod][iebin][ism]);
        if(!h[iprod][iebin][ism]) continue;
        
        TH1F* hRat = (TH1F*) h[iprod][iebin][ism]->Clone(Form("%s_Ratio",h[iprod][iebin][ism]->GetName()));
        hRat->Divide(h[0][iebin][ism]);
        
        if(iprod==0) hRat->Draw("H");
        else          hRat->Draw("H same");
        
        hRat->SetYTitle("Ratio MC to Data");
        hRat->SetMaximum(1.2);
        hRat->SetMinimum(0);
        
      } // iprod
      
    } // param
    
    cR->cd(ncolSM*nrowSM);
    lSM->Draw();
    
    fileName = Form("figures/Comparison_Ratio_M02_Ebin%d_%s%s_%s",
                    iebin,isoName[iana].Data(),cutName[cut].Data(),titleName.Data());
    
    fileName+=fileFormat;
    
    cR->Print(fileName);
   
  } // e bin
  
  // Each file is per SM
  // each pad has an E bin
  // last SM is the sum of all, already plotted before.
  for(Int_t ism = 0; ism < nSM-1; ism++)
  {
    gStyle->SetOptTitle(1);
    gStyle->SetPadTopMargin(0.07);
    
    TCanvas * c = new TCanvas(Form("c_sm%d_iana%d_cut%d",
                                   ism,iana,cut),
                              Form("iana %d, cut %d, SM %d",
                                   iana, cut, ism),
                              ncolE*2000,ncolE*2000);
    c->Divide(ncolE,nrowE);
//    
//    TLegend *l = new TLegend(-0.04,0,1,1);
//    l->SetFillColor(0);
//    l->SetFillStyle(0);
//    l->SetLineColor(0);
//    l->SetBorderSize(0);
//    l->SetTextSize(0.06);
//    l->SetHeader(Form("    SM %d", ism));
//    l->AddEntry("",Form(isoTitle[iana].Data(), cutTitle[cut].Data()),"");
//    //if(higE[iebin] > 90) l->SetHeader(Form("#it{E} > %2.1f GeV",lowE[iebin]));
    for(Int_t iebin = 0; iebin < nEBins; iebin++)
    {
      c->cd(iebin+1);
      
      //gPad->SetLogy();
      
      //printf("iE %d ism %d\n",iebin,ism);
      if(!h[0][iebin][ism]) continue;
      
      for(Int_t iprod = 0; iprod < nProd; iprod++)
      {
        //printf("\t prod %d %p\n",iprod,h[iprod][iebin][ism]);
        if(!h[iprod][iebin][ism]) continue;
        
        h[iprod][iebin][ism]->SetTitle(Form("%2.1f < #it{E}^{clus} < %2.1f GeV",lowE[iebin],higE[iebin]));

        h[iprod][iebin][ism]->SetTitleOffset(1.8,"Y");
        //h[iebin][0]->SetLineColor(1);
        //h[iebin][0]->SetAxisRange(0.1,2.5,"X");
        
        if ( iprod == 0 ) h[iprod][iebin][ism]->Draw("H");
        else               h[iprod][iebin][ism]->Draw("H same");
        
        if(h[iprod][iebin][ism]->GetMaximum() > h[0][iebin][ism]->GetMaximum())
          h[0][iebin][ism]->SetMaximum(h[iprod][iebin][ism]->GetMaximum()*1.2);
      } // iprod
      
      //printf("\t end\n");
    } // param
    
    c->cd(ncolE*nrowE);
    l->Draw();
    
    fileName = Form("figures/Comparison_M02_SM%d_%s%s_%s",
                    ism,isoName[iana].Data(),cutName[cut].Data(),titleName.Data());
    fileName+=fileFormat;
    c->Print(fileName);
  } // sm bin
  
}

///---------------------------------------------------------------------------
/// Execute all possible combinations
///
/// \param firstSM: first SM number to be inspected.
/// \param lastSM : last  SM number to be inspected.
/// \param plotRat: make ratio plots.
/// \param bAllSM : just recover the histogram filled for all SM if it exists and plot only those.
/// \param debug  : Bool to activate printfs.
///---------------------------------------------------------------------------
void CompareShowerShapeLongPerAnaPerSM
(Int_t  firstSM = 0,
 Int_t  lastSM  = 9,  
 Bool_t plotRat = kFALSE,
 Bool_t bAllSM  = kFALSE, 
 Bool_t debug   = kFALSE )
{
  for(Int_t iprod = 0; iprod < nProd; iprod++)
  {
    file[iprod] = TFile::Open(Form("%s.root",filePath[iprod].Data()));
    if ( debug )
      printf("Read: %s.root, %p\n",filePath[iprod].Data(),file);
  }
  
  // Execute the analysis
  //
  DoIt(0, 0, firstSM, lastSM, plotRat, bAllSM, debug); // Not iso  
//  DoIt(1, 0, firstSM, lastSM, plotRat, bAllSM, debug); // Iso
//  DoIt(2, 0, firstSM, lastSM, plotRat, bAllSM, debug); // Inclusive: AliAnaPhoton
//  
//  // With cut on Ncells
//  DoIt(0, 1, firstSM, lastSM, plotRat, bAllSM, debug); // Not iso 
//  DoIt(1, 1, firstSM, lastSM, plotRat, bAllSM, debug); // Iso 
}

