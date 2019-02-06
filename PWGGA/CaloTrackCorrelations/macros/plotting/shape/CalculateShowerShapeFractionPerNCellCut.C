///
/// \file CalculateShowerShapeFractionPerNCellCut.C
/// \ingroup CaloTrackCorrMacrosPlotting
/// \brief Get fraction of clusters with photon shape and large number of cells
///
/// Plot the shower shape distribution for isolated, non isolated and inclusive clusters
/// depending if the cluster has or not more than 4 cells with weight.
/// Calculate the fraction of such clusters over the total when they have a shower
/// shape like the one of photons, for different selection cuts 0.1<M02<0.27,0.3,0.4
/// Put the fractions vs E in output files to be processed later for comparison
/// in method  CompareProductions()
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)
///

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TFile.h>
#include <TDirectoryFile.h>
#include <TList.h>
#include <TString.h>
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
#include <TROOT.h>
#include <TStyle.h>
#include "PlotUtils.C"

#endif

//--------------------Global variables------------------------------------------

// Cluster types
TString anaName [] = {"NoIso","Iso","All"};
TString anaTitle[] = {"NOT Isolated","Isolated","Inclusive"};

// Cut cases, 0-no cut, 1-n cell > 4, 2-ncell <=4
const Int_t nHisto = 3;
TString histoLeg[] = {"#it{n}^{w}_{cell}>1", "#it{n}^{w}_{cell} > 4","#it{n}^{w}_{cell} #leq 4"};

// Energy bins
const Int_t nEBins = 18;
Double_t binE   [] = {  2,  3,  4,  5,  6,  7,  8,10,12,14,16,18,20,22, 25,  30,  35,  40,50};
Double_t binEErr[] = {0.5,0.5,0.5,0.5,0.5,0.5,0.5, 1, 1, 1, 1, 1, 1, 1,1.5, 2.5, 2.5, 2.5, 5};

// Shower shape max cut array
const Int_t nShShMaxCut = 3;
Float_t maxShShCut[] = {0.27,0.3,0.4};

//const Int_t nShShMaxCut = 1;
//Float_t maxShShCut[] = {0.3};

Int_t shshRef = 1; // take 0.3 for some plots 

// Histogram settings
Int_t color    [] = {1,4,2,6,8,9,kYellow-6,kCyan,kOrange+2,kViolet,kOrange-2,
  4,2,6,8,9,kYellow-6,kCyan,kOrange+2,kViolet,kOrange-2};
Int_t lineStyle[] = {1,1,1,1,1,1,        1,    1,        1,      
  1,        1,2,2,2,2,2,        2,    2,        2,      2,        
  2,2,2,2,2,2,2,2,2,2,2,2};
Int_t marker   [] = {20,20,20,21,24,24,24,24,24,24,24};

// Plot file format
TString fileFormat = ".eps";

//------------------------------------------------------------------------------
/// Get SM total number, first SM, and n col/rows depending data string name
///
/// \param titleName : String with simplified production name, it must contain data period and trigger for data and calorimeter
/// \param totalSM   : Number of SM
/// \param firstSM   : 0 for EMCal, 12 for DCal
/// \param ncol      : Number of canvas columns
/// \param nrow      : Number of canvas rows
//------------------------------------------------------------------------------
void GetSMNumber(TString titleName, 
                 Int_t & totalSM, Int_t & firstSM,
                 Int_t & ncol   , Int_t & nrow)
{
  totalSM = 10;
  firstSM = 0;
  ncol = 4;
  nrow = 3;
  GetCanvasColRowNumber(totalSM,ncol,nrow); // PlotUtils.C
  
  if(titleName.Contains("15") || titleName.Contains("16") || 
     titleName.Contains("17") || titleName.Contains("18"))
  {
    totalSM = 12;
    GetCanvasColRowNumber(totalSM,ncol,nrow); // PlotUtils.C
    if(titleName.Contains("DCAL") || titleName.Contains("DMC") || 
       titleName.Contains("DG")   || titleName.Contains("DMG"))
    {
      totalSM = 20;   
      firstSM = 12;
      GetCanvasColRowNumber(8,ncol,nrow); // PlotUtils.C
    }
  }
}

//------------------------------------------------------------------------------
/// Open the input histogram file and recover the histogram corresponding to the cluster type selected
/// do the plots and calculate the fraction of clusters and put them in output file.
///
/// \param titleName : String with simplified production name, it must contain data period and trigger for data and calorimeter
/// \param filePath  : String input file path
/// \param iana      : Selected cluster type: 0-not isolated, 1-isolated, 2-inclusive (AliAnaClusterShapeStudies)
/// \param tm        : String with track matching option used
/// \param plotRat   : Bool to activate ratio plots.
/// \param debug     : Bool to activate printfs.
//------------------------------------------------------------------------------
void CalculateAndPlot
(TString titleName, TString filePath, 
 Int_t   iana    = 0, 
 TString tm      = "_TMDep", 
 Bool_t  plotRat = kFALSE,
 Bool_t  debug   = kFALSE)
{  
  // ------------------------------------------------
  // Declare histogram arrays and different settings
  // ------------------------------------------------
  
  // Histogram settings
  Int_t    rebin = 2;
  
  Double_t xmin = 0.1;
  Double_t xmax = 1.1;
  
  // Energy bin settings
  // bins defined as global
  Double_t lowE[nEBins];
  Double_t higE[nEBins];
  Double_t muE [nEBins];
  Double_t errE[nEBins];
  Int_t ncolE = 5;
  Int_t nrowE = 4;
  GetCanvasColRowNumber(nEBins,ncolE,nrowE); // PlotUtils.C
  
  // SM
  //
  Int_t totalSM = 10, firstSM = 0;
  Int_t ncol = 4, nrow = 3;
  GetSMNumber(titleName.Data(),totalSM,firstSM,ncol,nrow);
  const Int_t nSM = totalSM;
  
  // histograms
  TH1D* h [nHisto][nEBins][nSM];
  TH1D* hA[nHisto][nEBins];

  // Shower shape max cut
  Double_t  frac[nEBins][nSM+1][nShShMaxCut];
  Double_t efrac[nEBins][nSM+1][nShShMaxCut];
  
  // Merge good/bad SMs
  const Int_t nGroup = 3;
  TH1D* hGroup[nHisto][nEBins][nGroup];
  
  // ------------------------------------------------
  // Open input file and fill histogram arrays
  // ------------------------------------------------
  
  // Recover data file
  TFile* file = TFile::Open(Form("%s.root",filePath.Data()));
  if ( debug )
    printf("Read: %s.root, %p\n",filePath.Data(),file);
  
  if ( !file ) return;
  
  //titleName.ReplaceAll("/","_");
  
  //
  // Histograms from AliAnaEMCalClusterShape
  // TH3 histogram, treat differently to AliAnaParticleIsolation
  //
  if ( iana == 2 )
  {
    TH3F* h3[nHisto];
    
    TString start = "Shape"; // Careful, it can change.
    h3[0] = (TH3F*) file->Get(Form("%s%s_hSMM02NoCut_Neutral",start.Data(),tm.Data()));
    h3[1] = (TH3F*) file->Get(Form("%s%s_hSMM02_Neutral"     ,start.Data(),tm.Data()));
    if ( debug )
      printf("AliAnaClusterShapeStudies TH3F %p %p\n",h3[0],h3[1]);
    
    if(filePath.Contains("data"))
    {
      h3[0]->Sumw2();
      h3[1]->Sumw2();
    }
    
    if(nHisto > 2) h3[2] = (TH3F*) h3[0]->Clone("Histo2_Iso2");
    h3[2]->Add(h3[1],-1);
    
    //
    // Project
    //
    Double_t width = 0;
    Int_t ebinMin  =-1;
    Int_t ebinMax  =-1;
    Int_t smbin    =-1;
    for(Int_t ihisto = 0; ihisto < nHisto; ihisto++)
    {      
      for(Int_t iebin = 0; iebin < nEBins; iebin++)
      {
        ebinMin     = h3[ihisto]->GetXaxis()->FindBin(binE[iebin  ]);
        ebinMax     = h3[ihisto]->GetXaxis()->FindBin(binE[iebin+1])-1;
        lowE[iebin] = h3[ihisto]->GetXaxis()->GetBinLowEdge(ebinMin);
        width       = h3[ihisto]->GetXaxis()->GetBinWidth  (ebinMax);
        higE[iebin] = h3[ihisto]->GetXaxis()->GetBinLowEdge(ebinMax)+width;
        muE [iebin] = (lowE[iebin]+higE[iebin])/2.;
        errE[iebin] = (higE[iebin]-lowE[iebin])/2.;
        
        for(Int_t ism = firstSM; ism < nSM; ism++)
        {
          smbin = h3[ihisto]->GetYaxis()->FindBin(ism);
          h[ihisto][iebin][ism] = 
          (TH1D*) h3[ihisto]->ProjectionZ(Form("Histo%d_BinE%d_SM%d",ihisto,iebin,ism),ebinMin,ebinMax,smbin,smbin);
      
          if ( !h[ihisto][iebin][ism] ) continue;
          
//          printf("%s, histo %d, iebin %d, [%2.1f,%2.1f] GeV, muE %2.1f errE %2.1f,ism %d, smbin%d integral %e\n",
//                 h[ihisto][iebin][ism]->GetName(),ihisto, 
//                 iebin,lowE[iebin],higE[iebin],muE[iebin],errE[iebin],
//                 ism,smbin,h[ihisto][iebin][ism]->Integral());
          
          h[ihisto][iebin][ism]->SetLineColor(color[ihisto]);
          
          h[ihisto][iebin][ism]->SetLineWidth(2);
          
          h[ihisto][iebin][ism]->SetLineStyle(lineStyle[ihisto]);
          
          h[ihisto][iebin][ism]->SetMarkerStyle(marker[ihisto]);
          
          h[ihisto][iebin][ism]->SetMarkerColor(color[ihisto]);
          
          h[ihisto][iebin][ism]->SetMarkerSize(0.5);
          
          if(rebin > 1) 
            h[ihisto][iebin][ism]->Rebin(rebin);
          
          h[ihisto][iebin][ism]->SetAxisRange(xmin,xmax,"X");
        } // SM histo
      } // E bin 
    } // cut histo
  } // END Histograms from AliAnaEMCalClusterShape
  //
  // Histograms from AliAnaParticleIsolation
  //
  else
  {
    TH2F* h2[nHisto][nSM];
    for(Int_t ism = firstSM; ism < nSM; ism++)
    {
      for(Int_t ihisto = 0; ihisto < nHisto; ihisto++)
      {
        if ( ihisto==0 ) h2[ihisto][ism] = (TH2F*) file->Get(Form("AnaIsolPhoton_hPtLambda0_SM%d_%s"         ,ism,anaName[iana].Data()));
        if ( ihisto==1 ) h2[ihisto][ism] = (TH2F*) file->Get(Form("AnaIsolPhoton_hPtLambda0_NCellCut_SM%d_%s",ism,anaName[iana].Data()));
        
        if ( ihisto==2 ) 
        {
          h2[ihisto][ism] = (TH2F*) h2[0][ism]->Clone(Form("Histo%d_Iso%d_SM%d",ihisto,iana,ism));
          h2[ihisto][ism]->Add(h2[1][ism],-1);
        }
        
        if ( ihisto < 2 )
        {
        if(filePath.Contains("data") ) h2[ihisto][ism]->Sumw2();
        if(filePath.Contains("simu") && filePath.Contains("MB") ) h2[ihisto][ism]->Sumw2();
        }
        
      } // cut histo
    }  //  SM histo
    
    //
    // Project
    //
    Double_t width = 0;
    Int_t ebinMin  =-1;
    Int_t ebinMax  =-1;
    for(Int_t ism = firstSM; ism < nSM; ism++)
    {
      for(Int_t ihisto = 0; ihisto < nHisto; ihisto++)
      {      
        for(Int_t iebin = 0; iebin < nEBins; iebin++)
        {
          ebinMin     = h2[ihisto][ism]->GetXaxis()->FindBin(binE[iebin  ]);
          ebinMax     = h2[ihisto][ism]->GetXaxis()->FindBin(binE[iebin+1])-1;
          lowE[iebin] = h2[ihisto][ism]->GetXaxis()->GetBinLowEdge(ebinMin);
          width       = h2[ihisto][ism]->GetXaxis()->GetBinWidth  (ebinMax);
          higE[iebin] = h2[ihisto][ism]->GetXaxis()->GetBinLowEdge(ebinMax)+width;
          muE [iebin] = (lowE[iebin]+higE[iebin])/2.;
          errE[iebin] = (higE[iebin]-lowE[iebin])/2.;
          
          h[ihisto][iebin][ism] = 
          (TH1D*) h2[ihisto][ism]->ProjectionY(Form("Histo%d_BinE%d_SM%d",ihisto,iebin,ism),ebinMin,ebinMax);
          
          //printf("\t %p \n",h[ihisto][iebin][ism]);
          
          if ( !h[ihisto][iebin][ism] ) continue;
          
//          printf("%s, histo %d, iebin %d, [%2.1f,%2.1f] GeV, muE %2.1f errE %2.1f,ism %d, integral %e\n",
//                 h[ihisto][iebin][ism]->GetName(),ihisto, 
//                 iebin,lowE[iebin],higE[iebin],muE[iebin],errE[iebin],
//                 ism,h[ihisto][iebin][ism]->Integral());
          
          h[ihisto][iebin][ism]->SetLineColor(color[ihisto]);
          
          h[ihisto][iebin][ism]->SetLineWidth(2);
          
          h[ihisto][iebin][ism]->SetLineStyle(lineStyle[ihisto]);
          
          h[ihisto][iebin][ism]->SetMarkerStyle(marker[ihisto]);
          
          h[ihisto][iebin][ism]->SetMarkerColor(color[ihisto]);
          
          h[ihisto][iebin][ism]->SetMarkerSize(0.5);
          
          if(rebin > 1) 
            h[ihisto][iebin][ism]->Rebin(rebin);
          
          h[ihisto][iebin][ism]->SetAxisRange(xmin,xmax,"X");
        } // ie
      } // ihisto
    } // ism
  } // END Histograms from AliAnaParticleIsolation

    
  //
  // Merged SM
  //
  for(Int_t ihisto = 0; ihisto < nHisto; ihisto++)
  {      
    for(Int_t iebin = 0; iebin < nEBins; iebin++)
    {            
      if ( h[ihisto][iebin][firstSM] )
      {
        hA [ihisto][iebin] = (TH1D*) h[ihisto][iebin][firstSM]->Clone(Form("SMAll_Cut%d_Ana%d" ,ihisto,iana));

        for(Int_t ism = firstSM+1; ism < nSM; ism++)
          hA[ihisto][iebin]->Add( h[ihisto][iebin][ism] );
      }
  
      
      if( !titleName.Contains("DCAL") && !titleName.Contains("DMC") && 
          !titleName.Contains("DG")   && !titleName.Contains("DMG"))     
      {
        hGroup[ihisto][iebin][0] = (TH1D*) h[ihisto][iebin][0]->Clone(Form("SM_Group0_Cut%d_Ana%d",ihisto, iana));
        hGroup[ihisto][iebin][1] = (TH1D*) h[ihisto][iebin][1]->Clone(Form("SM_Group1_Cut%d_Ana%d",ihisto, iana));
        hGroup[ihisto][iebin][2] = (TH1D*) h[ihisto][iebin][3]->Clone(Form("SM_Group2_Cut%d_Ana%d",ihisto, iana));
        
        hGroup[ihisto][iebin][0]->Add( h[ihisto][iebin][4] );
        hGroup[ihisto][iebin][0]->Add( h[ihisto][iebin][5] );
        hGroup[ihisto][iebin][0]->Add( h[ihisto][iebin][6] );
        hGroup[ihisto][iebin][0]->Add( h[ihisto][iebin][8] );
        hGroup[ihisto][iebin][0]->Add( h[ihisto][iebin][9] );
        
        hGroup[ihisto][iebin][1]->Add( h[ihisto][iebin][2] );

        hGroup[ihisto][iebin][2]->Add( h[ihisto][iebin][7] );
        
        hGroup[ihisto][iebin][0]->Scale(1./6.);
        hGroup[ihisto][iebin][1]->Scale(1./2.);
        hGroup[ihisto][iebin][2]->Scale(1./2.);
      }
      
    } // ebins
  } // prod bins

  // ---------------------------------------
  // Plot 
  // ---------------------------------------
  if ( debug )
    printf("PLOT\n");
  
  TString fileName ;
  
  gStyle->SetPadRightMargin(0.01);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetTitleFontSize(0.06);
  //  gStyle->SetStatFontSize(0.06);
  gStyle->SetOptStat(0);
  
  // 
  // Plot comparison of 3 cut cases and calculate fractions per E bin
  //
  Int_t shshRef = 1; //Select a max Cut to be shown in the legend.
  
  // Plots per E bin, each par per SM
  for(Int_t iebin = 0; iebin < nEBins; iebin++)
  {
    gStyle->SetOptTitle(1);
    gStyle->SetPadTopMargin(0.1);
    
    TCanvas * c = new TCanvas(Form("c_ebin%d_iana%d_%s",
                                   iebin,iana, titleName.Data()),
                              Form("iana %d, %2.1f < E < %2.1f, %s",
                                   iana, lowE[iebin],higE[iebin], titleName.Data()),
                              ncol*2000,nrow*2000);
    c->Divide(ncol,nrow);
    
    TLegend *l = new TLegend(0,0,1,1);
    l->SetFillColor(0);
    l->SetFillStyle(0);
    l->SetLineColor(0);
    l->SetBorderSize(0);
    l->SetTextSize(0.07);
    l->SetHeader(Form("%2.1f < #it{E} < %2.1f GeV", lowE[iebin],higE[iebin]));
    l->AddEntry("",anaTitle[iana].Data(),"");

    for(Int_t ism = firstSM; ism < nSM; ism++)
    {
      c->cd(ism+1);
      
      //gPad->SetLogy();
      if ( debug )
        printf("iana %d iE %d ism %d: %p %p %p\n",
               iana,iebin,ism,h[0][iebin][ism],h[1][iebin][ism],h[2][iebin][ism]);
      
      if(!h[0][iebin][ism]) continue;
      
      for(Int_t ihisto = 0; ihisto < nHisto; ihisto++)
      {
        //printf("\t prod %d %p\n",ihisto,h[ihisto][iebin][ism]);
        if(!h[ihisto][iebin][ism]) continue;
        
        h[ihisto][iebin][ism]->SetTitle(Form("SM %d",ism));
        
        h[ihisto][iebin][ism]->SetTitleOffset(1.8,"Y");
        //h[iebin][0]->SetLineColor(1);
        //h[iebin][0]->SetAxisRange(0.1,2.5,"X");
        
        if ( ihisto == 0 ) h[ihisto][iebin][ism]->Draw("H");
        else               h[ihisto][iebin][ism]->Draw("H same");
        
        if(h[ihisto][iebin][ism]->GetMaximum() > h[0][iebin][ism]->GetMaximum())
          h[0][iebin][ism]->SetMaximum(h[ihisto][iebin][ism]->GetMaximum()*1.2);
        
        if(ism==firstSM)
          l->AddEntry(h[ihisto][iebin][ism],Form("%s",histoLeg[ihisto].Data()),"PL");
      } // ihisto
      
      // *****
      // Calculate the fraction of clusters with n cell > 4
      // *****
      Double_t integral0    = 0;
      Double_t integral1    = 0; 
      Double_t integral0Err = 0;
      Double_t integral1Err = 0;
      Int_t    binM02Min    = h[0][iebin][ism]->FindBin(0.1);
      for(Int_t imax = 0; imax < nShShMaxCut; imax++)
      {
        Int_t binM02Max = h[0][iebin][ism]->FindBin(maxShShCut[imax]);
        
        GetRangeIntegralAndError(h[0][iebin][ism], binM02Min, binM02Max, integral0, integral0Err );
        GetRangeIntegralAndError(h[1][iebin][ism], binM02Min, binM02Max, integral1, integral1Err );
        
        Double_t ratio = 0;
        if ( integral0 > 0 ) ratio = integral1 / integral0 ;
        Double_t eratio  = GetFractionError(integral1,integral0,integral1Err,integral0Err);
        
        frac [iebin][ism][imax] =  ratio*100;
        efrac[iebin][ism][imax] = eratio*100;
        
        if( debug )
          printf("\t max %2.2f, frac %2.2f, err %2.2f\n",
                 maxShShCut[imax],frac [iebin][ism][imax],efrac[iebin][ism][imax]);
      }
      
      TLegend *lsm = new TLegend(0.15,0.7,0.95,0.9);
      lsm->SetFillColor(0);
      lsm->SetFillStyle(0);
      lsm->SetLineColor(0);
      lsm->SetBorderSize(0);
      lsm->SetTextSize(0.07);
      lsm->SetHeader(Form("        For 0.10<#sigma_{long}^{2}<%2.2f",maxShShCut[shshRef]));
      lsm->AddEntry("",Form("#it{R}_{#sigma_{long}}^{#it{n}^{w}_{cell}} = %2.1f #pm %2.1f %%",frac[iebin][ism][shshRef],efrac[iebin][ism][shshRef]),"") ;
      lsm->Draw("same");
      //printf("\t end\n");
    } // param
    
    c->cd(ncol*nrow);
    l->Draw();
    
    fileName = Form("figures/Comparison_M02_NCellCut_Ebin%d_%s_%s",
                    iebin,anaName[iana].Data(),titleName.Data());
    fileName+=fileFormat;
    c->Print(fileName);
    
    if ( !plotRat ) continue;
    
    TCanvas * cR = new TCanvas(Form("cR_ebin%d_iana%d_%s",
                                    iebin,iana, titleName.Data()),
                               Form("Ratio iana %d, %2.1f < E < %2.1f, %s",
                                    iana, lowE[iebin],higE[iebin], titleName.Data()),
                               ncol*2000,nrow*2000);
    cR->Divide(ncol,nrow);
    
    for(Int_t ism = 0; ism < nSM; ism++)
    {
      cR->cd(ism+1);
      
      //gPad->SetLogy();
      
      //printf("iE %d ism %d\n",iebin,ism);
      if(!h[1][iebin][ism]) continue;
      
      for(Int_t ihisto = 1; ihisto < nHisto; ihisto++)
      {
        //printf("\t prod %d %p\n",ihisto,h[ihisto][iebin][ism]);
        if(!h[ihisto][iebin][ism]) continue;
        
        TH1F* hRat = (TH1F*) h[ihisto][iebin][ism]->Clone(Form("%s_Ratio",h[ihisto][iebin][ism]->GetName()));
        hRat->Divide(h[0][iebin][ism]);
        
        if(ihisto==0) hRat->Draw("H");
        else          hRat->Draw("H same");
        
        hRat->SetYTitle("Ratio MC to Data");
        hRat->SetMaximum(1.2);
        hRat->SetMinimum(0);
        
      } // ihisto
      
    } // param
    
    cR->cd(ncol*nrow);
    l->Draw();
    
    fileName = Form("figures/Comparison_Ratio_M02_NCellCut_Ebin%d_%s_%s",
                    iebin,anaName[iana].Data(),titleName.Data());
    fileName+=fileFormat;
    
    cR->Print(fileName);
    
  } // e bin
  
  // Sum of all SM, each pad one energy
  //
  TCanvas * cA = new TCanvas(Form("cAllSM_iana%d_%s",iana, titleName.Data()),
                            Form("All SM, iana %d, %s",iana, titleName.Data()),
                             ncolE*2000,nrowE*2000);
  
  cA->Divide(ncolE,nrowE);
  
  TLegend *l2 = new TLegend(0,0.,1,1);
  l2->SetFillColor(0);
  l2->SetFillStyle(0);
  l2->SetLineColor(0);
  l2->SetBorderSize(0);
  l2->SetTextSize(0.07);
  l2->AddEntry("",Form("All SM, %s",anaTitle[iana].Data()),"");

  for(Int_t iebin = 0; iebin < nEBins; iebin++)
  {
    cA->cd(iebin+1);
    
    if ( debug )
      printf("All SM, iana %d, ebin %d: %p, %p, %p\n",
             iana,iebin,hA[0][iebin],hA[1][iebin],hA[2][iebin]);
    
    if(!hA[0][iebin]) continue;
    
    //gPad->SetGridy();
    //gPad->SetGridx();
    //gPad->SetLogy();
    
    gStyle->SetOptTitle(1);
    gStyle->SetPadTopMargin(0.1);
    
    hA[0][iebin]->Draw("H");
    
    for(Int_t ihisto = 0; ihisto < nHisto; ihisto++)
    {
      if(!hA[ihisto][iebin]) continue; 
      
      hA[ihisto][iebin]->SetTitleOffset(1.8,"Y");
      //h[iebin][0]->SetAxisRange(0.1,2.5,"X");
      //hA[ihisto][iebin]->SetMinimum(1e-2);
      //printf("iebin %d, ism %d 0 %d %p\n",iebin,ism,0,h[0][iebin]);
      
      hA[ihisto][iebin]->SetTitle(Form("%2.1f < #it{E}^{clus} < %2.1f GeV",lowE[iebin],higE[iebin]));
      if(higE[iebin] > 90) hA[ihisto][iebin]->SetTitle(Form("#it{E} > %2.1f GeV",lowE[iebin]));
      
      hA[ihisto][iebin]->Draw("H same");
      
      if(hA[ihisto][iebin]->GetMaximum() > hA[0][iebin]->GetMaximum())
        hA[0][iebin]->SetMaximum(hA[ihisto][iebin]->GetMaximum()*1.2);
      
      if(iebin==0)
        l2->AddEntry(hA[ihisto][iebin],Form("%s",histoLeg[ihisto].Data()),"PL");      
    }  
    
    // *****
    // Calculate the fraction of clusters with n cell > 4
    // *****
    Double_t integral0    = 0;
    Double_t integral1    = 0; 
    Double_t integral0Err = 0;
    Double_t integral1Err = 0;
    Int_t    binM02Min    = hA[0][iebin]->FindBin(0.1);
    for(Int_t imax = 0; imax < nShShMaxCut; imax++)
    {
      Int_t binM02Max = hA[0][iebin]->FindBin(maxShShCut[imax]);
      
      GetRangeIntegralAndError(hA[0][iebin], binM02Min, binM02Max, integral0, integral0Err );
      GetRangeIntegralAndError(hA[1][iebin], binM02Min, binM02Max, integral1, integral1Err );
      
      Double_t ratio = 0;
      if ( integral0 > 0 ) ratio = integral1 / integral0 ;
      Double_t eratio  = GetFractionError(integral1,integral0,integral1Err,integral0Err);
      
      //    printf("A ratio %2.3f, eratio %2.3f;"
      //           " int0 %2.3e, int1 %2.3e;"
      //           " Eint0 %2.3e, Eint1 %2.3e;"
      //           " sqrt0 %2.3f sqrt1 %2.3f \n",
      //           ratio,eratio,
      //           integral0,integral1,
      //           integral0Err,integral1Err,
      //           TMath::Sqrt(integral0),TMath::Sqrt(integral1));
      
      frac [iebin][nSM][imax] =  ratio*100;
      efrac[iebin][nSM][imax] = eratio*100;
      
      if( debug )
        printf("\t max %2.2f, frac %2.2f, err %2.2f\n",
               maxShShCut[imax],frac [iebin][nSM][imax],efrac[iebin][nSM][imax]);
    }
    
    TLegend *lE = new TLegend(0.15,0.7,0.95,0.9);
    lE->SetFillColor(0);
    lE->SetFillStyle(0);
    lE->SetLineColor(0);
    lE->SetBorderSize(0);
    lE->SetTextSize(0.07);
    lE->SetHeader(Form("        For 0.10<#sigma_{long}^{2}<%2.2f",maxShShCut[shshRef]));
    lE->AddEntry("",Form("#it{R}_{#sigma_{long}}^{#it{n}^{w}_{cell}} = %2.1f #pm %2.1f %%",frac[iebin][nSM][shshRef],efrac[iebin][nSM][shshRef]),"") ;
    lE->Draw("same");
  }
  
  cA->cd(ncolE*nrowE);
  l2->Draw("same");
  
  fileName = Form("figures/Comparison_M02_NCellCut_AllSM_%s_%s",
                  anaName[iana].Data(),titleName.Data());
  fileName+=fileFormat;
  cA->Print(fileName);

  if ( plotRat )
  {
    TCanvas * cAR = new TCanvas(Form("cR_AllSM_iana%d_%s",iana, titleName.Data()),
                                Form("All SM, iana %d, %s",iana, titleName.Data()),
                                ncolE*2000,nrowE*2000);
    
    cAR->Divide(ncolE,nrowE);
    
    for(Int_t iebin = 0; iebin < nEBins; iebin++)
    {
      cAR->cd(iebin+1);
      
      if(!hA[0][iebin]) continue;
      
      //if(quantity.Contains("ECell"))gPad->SetLogy();
      
      gStyle->SetOptTitle(1);
      gStyle->SetPadTopMargin(0.1);
      
      hA[0][iebin]->Draw("H");
      
      for(Int_t ihisto = 1; ihisto < nHisto; ihisto++)
      {
        if(!hA[ihisto][iebin]) continue; 
        
        TH1F* hRat = (TH1F*) hA[ihisto][iebin]->Clone(Form("%s_Ratio",hA[ihisto][iebin]->GetName()));
        hRat->Divide(hA[0][iebin]);
        
        if(ihisto==1) hRat->Draw("H");
        else         hRat->Draw("H same");
        
        hRat->SetYTitle("Ratio MC to Data");
        hRat->SetMaximum(1.2);
        hRat->SetMinimum(0);
      }
    }
    
    cAR->cd(ncolE*nrowE);
    l2->Draw("same");
    
    fileName = Form("figures/Comparison_Ratio_M02_NCellCut_AllSM_%s_%s",
                    anaName[iana].Data(),titleName.Data());
    
    fileName+=fileFormat;
    cAR->Print(fileName);
  }
  
  // Now plot each file per SM, per E bin in each pad
  //
  for(Int_t ism = firstSM; ism < nSM; ism++)
  {
    TCanvas * cSM = new TCanvas(Form("cSM%d_iana%d_%s",ism,iana,titleName.Data()),
                                Form("SM %d, iana %d, %s",ism,iana,titleName.Data()),
                                ncolE*2000,nrowE*2000);
    
    cSM->Divide(ncolE,nrowE);
    
    TLegend *l2 = new TLegend(0,0.,1,1);
    l2->SetFillColor(0);
    l2->SetFillStyle(0);
    l2->SetLineColor(0);
    l2->SetBorderSize(0);
    l2->SetTextSize(0.07);
    l2->AddEntry("",Form("SM %d, %s",ism,anaTitle[iana].Data()),"");
    
    for(Int_t iebin = 0; iebin < nEBins; iebin++)
    {
      cSM->cd(iebin+1);
      
      if(!h[0][iebin][ism]) continue;
      
      //gPad->SetLogy();
      
      gStyle->SetOptTitle(1);
      gStyle->SetPadTopMargin(0.1);
      
      h[0][iebin][ism]->Draw("H");
      
      for(Int_t ihisto = 0; ihisto < nHisto; ihisto++)
      {
        if(!h[ihisto][iebin][ism]) continue;         
        
        h[ihisto][iebin][ism]->SetTitle(Form("%2.1f < #it{E}^{clus} < %2.1f GeV",lowE[iebin],higE[iebin]));
        
        h[ihisto][iebin][ism]->Draw("H same");
        
        if(h[ihisto][iebin][ism]->GetMaximum() > h[0][iebin][ism]->GetMaximum())
          h[0][iebin][ism]->SetMaximum(h[ihisto][iebin][ism]->GetMaximum()*1.2);
        
        if(iebin==0)
          l2->AddEntry(h[ihisto][iebin][ism],Form("%s",histoLeg[ihisto].Data()),"PL");      
      }  
      
      TLegend *lE = new TLegend(0.15,0.7,0.95,0.9);
      lE->SetFillColor(0);
      lE->SetFillStyle(0);
      lE->SetLineColor(0);
      lE->SetBorderSize(0);
      lE->SetTextSize(0.07);
      lE->SetHeader(Form("        For 0.10<#sigma_{long}^{2}<%2.2f",maxShShCut[shshRef]));
      lE->AddEntry("",Form("#it{R}_{#sigma_{long}}^{#it{n}^{w}_{cell}} = %2.1f #pm %2.1f %%",frac[iebin][ism][shshRef],efrac[iebin][ism][shshRef]),"") ;
      lE->Draw("same");
    }
    
    cSM->cd(ncolE*nrowE);
    l2->Draw("same");
    
    fileName = Form("figures/Comparison_M02_NCellCut_SM%d_%s_%s",
                    ism,anaName[iana].Data(),titleName.Data());
    fileName+=fileFormat;
    cSM->Print(fileName);
  } // sm
  
  //=================================================================
  // Put the calculated fractions in GRAPHS
  //=================================================================
  TGraphErrors* gFracPerEn[nSM+1] [nShShMaxCut];
  TGraphErrors* gFracPerSM[nEBins][nShShMaxCut];

  Double_t  fracEn[nEBins];
  Double_t efracEn[nEBins];
  Double_t  fracSM[nSM+1];
  Double_t efracSM[nSM+1];
  
  Double_t binSM    [] = {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22};
  Double_t binSMErr [] = {0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5};
  for(Int_t imax = 0; imax < nShShMaxCut; imax++)
  {
    // Graphs frac vs SM per E bin
    //
    for(Int_t ie = 0; ie < nEBins; ie++)
    {
      for(Int_t ism = firstSM; ism < nSM+1; ism++)
      {
         fracSM[ism] =  frac[ie][ism][imax];
        efracSM[ism] = efrac[ie][ism][imax];      
      }
      
      gFracPerSM[ie][imax] = new TGraphErrors(nSM+1,binSM,fracSM,binSMErr,efracSM);
    } // iE bin
    
    // Graphs frac vs E per SM
    //
    for(Int_t ism = firstSM; ism < nSM+1; ism++)
    {
      for(Int_t ie = 0; ie < nEBins; ie++)
      {
        
        fracEn[ie] =  frac[ie][ism][imax];
        efracEn[ie] = efrac[ie][ism][imax];
      }
      
      gFracPerEn[ism][imax] = new TGraphErrors(nEBins,muE,fracEn,errE,efracEn);
      
      // Set E ranges depending the data type
      //
      if( filePath.Contains("simu") )
      {
        if(filePath.Contains("High"))
        {
          gFracPerEn[ism][imax]->GetHistogram()->SetAxisRange(16,45,"X");
        }
        else if(filePath.Contains("Low"))
        {
          gFracPerEn[ism][imax]->GetHistogram()->SetAxisRange(7,20,"X");
        }      
        else if(filePath.Contains("MB"))
        {
          gFracPerEn[ism][imax]->GetHistogram()->SetAxisRange(2,10,"X");
        }
      } // simu
      else
      {
        if( filePath.Contains("EMC") ||  filePath.Contains("DMC") ||  
           filePath.Contains("EG")   ||  filePath.Contains("DG") )
        {
          gFracPerEn[ism][imax]->GetHistogram()->SetAxisRange(6,50,"X");
        }  
        else if ( filePath.Contains("MB") || filePath.Contains("INT") )
        {
          gFracPerEn[ism][imax]->GetHistogram()->SetAxisRange(2,10,"X");
          if(filePath.Contains("LHC17"))
          {
            gFracPerEn[ism][imax]->GetHistogram()->SetAxisRange(2,14,"X");
          }
        }
      } // simu
    } // ism
  
  } // imax 

  // --------------------------------------------------
  // Plot fractions graphs, each pad one max Shower shape cut
  // --------------------------------------------------
  
//  Float_t colorSM[] = {2,kBlue-3,kBlue+3,kRed-3,8,kYellow-2,kOrange-2,kRed+3,kCyan, kViolet,1};
//  Float_t styleSM[] = {21,24,24,22,21,21,21,22,21,21,20};
  
  Int_t colorSM[]={1 , 1, 2, 2, 3, 3, 4, 4, 7, 7, 6, 6, 2, 3, 4, 7, 6, 2, 2, 3, 3, 4, 4, 6, 6};
  Int_t styleSM[]={24,25,25,24,25,24,25,24,25,24,25,21,21,21,21,21,22,26,22,26,22,26,22,26};
  
  TCanvas * cGraphE = new TCanvas(Form("cGraphE_iana%d_%s",iana, titleName.Data()),
                                  Form("Graph E bin, iana %d, %s",iana, titleName.Data()),
                                  2*2000,2*2000);
  
  cGraphE->Divide(2,2);
  
  TLegend *lSM = new TLegend(0.1,0.1,1,1);
  lSM->SetFillColor(0);
  lSM->SetFillStyle(0);
  lSM->SetLineColor(0);
  lSM->SetBorderSize(0);
  lSM->SetTextSize(0.07);  
  
  for(Int_t imax = 0; imax < nShShMaxCut; imax++)
  {
    cGraphE->cd(imax+2);
    
    //gPad->SetGridy();
    //gPad->SetGridx();
    
    gFracPerEn[nSM][imax]->SetTitle(Form("0.10<#sigma_{long}^{2}<%2.2f, %s",maxShShCut[imax],anaTitle[iana].Data()));
    gFracPerEn[nSM][imax]->Draw("AP");
    
    for(Int_t ism = firstSM; ism < nSM+1; ism++)
    {    
      gFracPerEn[ism][imax]->Draw("P");
      gFracPerEn[ism][imax]->GetYaxis()->SetTitle("#it{R}_{#sigma_{long}}^{#it{n}^{w}_{cell}} (%)");
      if ( iana < 2 ) gFracPerEn[ism][imax]->GetXaxis()->SetTitle("#it{E}_{T} (GeV)");
      else           gFracPerEn[ism][imax]->GetXaxis()->SetTitle("#it{E} (GeV)");
      gFracPerEn[ism][imax]->SetMaximum(50);
      gFracPerEn[ism][imax]->SetMinimum(0);
      
      if(filePath.Contains("Default"))
      {
        gFracPerEn[ism][imax]->SetMaximum(5);
        gFracPerEn[ism][imax]->SetMinimum(0);      
      }
      
      gFracPerEn[ism][imax]->SetMarkerColor(colorSM[ism]);
      gFracPerEn[ism][imax]->SetLineColor  (colorSM[ism]);
      gFracPerEn[ism][imax]->SetMarkerStyle(styleSM[ism]);
      gFracPerEn[ism][imax]->SetMarkerSize (5);
      if(ism!=nSM)
      {
        if ( imax == 0 ) lSM->AddEntry(gFracPerEn[ism][imax],Form("SM %d",ism),"PL");
      }
      else    
      {
        gFracPerEn[ism][imax]->SetMarkerColor(1);
        gFracPerEn[ism][imax]->SetLineColor  (1);
        gFracPerEn[ism][imax]->SetMarkerStyle(20);
        if ( imax == 0 ) lSM->AddEntry(gFracPerEn[ism][imax],"All SM","PL");
      }
    }
  }
  cGraphE->cd(1);
  
  lSM->Draw();
  
  fileName = Form("figures/Comparison_Fraction_M02PhotonToAll_NCellCut_YaxisE_PerSM_%s_%s",
                  anaName[iana].Data(),titleName.Data());
  
  fileName+=fileFormat;
  cGraphE->Print(fileName);
  
  
  gStyle->SetOptTitle(1);

  // Plot each max cut individually 
  //  
  for(Int_t imax = 0; imax < nShShMaxCut; imax++)
  {
    TCanvas * cGraphEMax = new TCanvas(Form("cGraphE_Sh%d_iana%d_%s",imax,iana, titleName.Data()),
                                       Form("Graph E bin, Sh < %2.2f, iana %d, %s", maxShShCut[imax], iana, titleName.Data()),
                                       1*2000,1*2000);
    
    //gPad->SetGridy();
    //gPad->SetGridx();
    
    //gFracPerEn[nSM][imax]->SetMaximum(25);
    //gFracPerEn[nSM][imax]->SetMinimum(0);
    
    TLegend *lSM = new TLegend(0.85,0.2,0.98,0.85);
    lSM->SetFillColor(1);
    lSM->SetFillStyle(1);
    lSM->SetLineColor(0);
    lSM->SetBorderSize(0);
    lSM->SetTextSize(0.03);  
    
    TH1F * hAxis = new TH1F(Form("hAxis_Max%d",imax),
                            Form("0.1<#sigma_{long}^{2}<%2.2f, %s",maxShShCut[imax], anaTitle[iana].Data()),60,5,65);
    hAxis->SetYTitle("#it{R}_{#sigma_{long}}^{#it{n}^{w}_{cell}} (%)");
    if ( iana < 2 )  hAxis->SetXTitle("#it{E}_{T} (GeV)");
    else             hAxis->SetXTitle("#it{E} (GeV)");  
    //hAxis->SetTitle(Form("0.1<#sigma_{long}^{2}<0.3, %s",anaTitle[iana].Data()));
    hAxis->SetMaximum(50);
    hAxis->SetMinimum(0);
    hAxis->Draw("");
    
    gFracPerEn[nSM][imax]->Draw("P");
    lSM->AddEntry(gFracPerEn[nSM][imax],"All SM","PL");
    
    for(Int_t ism = firstSM; ism < nSM+1; ism++)
    {    
      gFracPerEn[ism][imax]->Draw("P");
      if(ism < nSM) lSM->AddEntry(gFracPerEn[ism][imax],Form("SM %d",ism),"PL");   
    }
    
    lSM->Draw();
    
    fileName = Form("figures/Comparison_Fraction_M02PhotonToAll_NCellCut_YaxisE_PerSM_ShMax%d_%s_%s",
                    imax, anaName[iana].Data(),titleName.Data());
    
    fileName+=fileFormat;
    cGraphEMax->Print(fileName);
  }
  
  // --------------------------------------------------
  // Plot fractions graphs vs SM, each pad one max Shower shape cut
  // --------------------------------------------------
  TCanvas * cGraphSM = new TCanvas(Form("cGraphSM_iana%d_%s",iana, titleName.Data()),
                                   Form("Graph E bin, iana %d, %s",iana, titleName.Data()),
                                   2*2000,2*2000);
  
  cGraphSM->Divide(2,2);
  
  TLegend *lE = new TLegend(0.1,0.1,1,1);
  lE->SetFillColor(0);
  lE->SetFillStyle(0);
  lE->SetLineColor(0);
  lE->SetBorderSize(0);
  lE->SetTextSize(0.07);  
  
  for(Int_t imax = 0; imax < nShShMaxCut; imax++)
  {
    cGraphSM->cd(imax+2);
    
    gFracPerSM[0][imax]->SetTitle(Form("0.10<#sigma_{long}^{2}<%2.2f, %s",maxShShCut[imax],anaTitle[iana].Data()));
    gFracPerSM[0][imax]->Draw("AP");
    
    for(Int_t ie = 0; ie < nEBins; ie++)
    {    
      gFracPerSM[ie][imax]->Draw("P");
      gFracPerSM[ie][imax]->GetYaxis()->SetTitle("Fraction %");
      gFracPerSM[ie][imax]->GetXaxis()->SetTitle("SM");
      gFracPerSM[ie][imax]->SetMaximum(50);
      gFracPerSM[ie][imax]->SetMinimum(0);
      if ( filePath.Contains("Default") )
      {
        gFracPerSM[ie][imax]->SetMaximum(5);
        gFracPerSM[ie][imax]->SetMinimum(0);      
      }
      gFracPerSM[ie][imax]->SetMarkerColor(color[ie]);
      gFracPerSM[ie][imax]->SetLineColor  (color[ie]);
      gFracPerSM[ie][imax]->SetMarkerStyle(24);
      gFracPerSM[ie][imax]->SetMarkerSize (5);
      if ( imax == 0 ) lE->AddEntry(gFracPerSM[ie][imax],Form("%2.f <E %2.1f GeV",lowE[ie],higE[ie]),"PL");
    }
  }  
  cGraphSM->cd(1);
  lE->Draw();
  
  fileName = Form("figures/Comparison_Fraction_M02PhotonToAll_NCellCut_YaxisSM_PerE_%s_%s",
                  anaName[iana].Data(),titleName.Data());
  
  fileName+=fileFormat;
  cGraphSM->Print(fileName);
  
  //****************************************************
  // Save fractions in histogram
  //****************************************************
  TFile *fout = new TFile(Form("figures/Fraction_M02PhotonToAll_NCellCut_%s_%s.root",
                               anaName[iana].Data(),titleName.Data()),"recreate");
  for(Int_t imax = 0; imax < nShShMaxCut; imax++)
  {
    for(Int_t ism = firstSM; ism < nSM+1; ism++)
    {    
      gFracPerEn[ism][imax]->SetName(Form("gFractionFrom%1.2f_SM%d_PerE",maxShShCut[imax], ism));
      gFracPerEn[ism][imax]->Write();
    }
    
    for(Int_t ie = 0; ie < nEBins; ie++)
    {
      gFracPerSM[ie][imax]->SetName(Form("gFractionFrom%1.2f_Ebin%d_PerSM",maxShShCut[imax],ie));
      gFracPerSM[ie][imax]->Write();
    }
  }
  fout->Print();
  fout->Close();
  
  file->Close();
  delete file;
}


//------------------------------------------------------------------------------
/// Open the files with output graphs in different productions from CalculateAndPlot 
/// and recover the graph corresponding to the cluster type selected and plot 
/// the production comparison
///
/// \param nProd     : Number of productions
/// \param fileName  : String input file path, file with graphs
/// \param legName   : String input production title
/// \param iana      : Selected cluster type: 0-not isolated, 1-isolated, 2-inclusive (AliAnaClusterShapeStudies)
/// \param debug     : Bool to activate printfs.
//------------------------------------------------------------------------------
void CompareProductions
( Int_t nProd, 
  TString * fileName,  
  TString * legName,
  Int_t  iana  = 0,
  Bool_t debug = kFALSE)
{        
  // SM
  //
  Int_t totalSM = 10, firstSM = 0;
  Int_t ncolSM = 4, nrowSM = 3;
  GetSMNumber(fileName[0].Data(),totalSM,firstSM,ncolSM,nrowSM);
  const Int_t nSM = totalSM+1; // Add last SM containing the sum of all
  
  if ( debug )
    printf("Compare %d productions for ana %d, nSM %d+1, first SM %d, ncol %d, nrow %d\n",
           nProd,iana,totalSM,firstSM,ncolSM,nrowSM);
  
  // Init arrays
  TFile        * file[nProd];
  TGraphErrors * gFrac[nSM][nShShMaxCut][nProd];
  TH1F         * hAxis[nShShMaxCut][nSM];
  
  // Get the files  with graphs
  //
  for(Int_t iprod = 0; iprod < nProd; iprod++)
  {
    file[iprod] = new TFile(Form("figures/Fraction_M02PhotonToAll_NCellCut_%s_%s.root",
                                 anaName[iana].Data(),fileName[iprod].Data()),"read");
    
    if ( debug )
    {
      printf("figures/Fraction_M02PhotonToAll_NCellCut_%s_%s.root \n",
             anaName[iana].Data(),fileName[iprod].Data());
      printf("iprod %d, file %p\n",iprod,file[iprod]);
    }
    
    if ( !file[iprod] ) continue;
    
    for(Int_t ism = firstSM; ism < nSM; ism++)
    {    
      for(Int_t im02 = 0; im02 < nShShMaxCut; im02++)
      {    
        gFrac[ism][im02][iprod] = (TGraphErrors*) file[iprod]->Get(Form("gFractionFrom%1.2f_SM%d_PerE",maxShShCut[im02],ism));
        
        if ( debug )
        {
          printf("\t gFractionFrom%1.2f_SM%d_PerE \n",maxShShCut[im02],ism);
          printf("\t ism %d, im02 %d graph %p \n",ism,im02,gFrac[ism][im02]);
        }
        
        if ( !gFrac[ism][im02][iprod] ) continue;
        
        // Set the plotting ranges, removing the unwanted points
        // Depending on the production name
        if(fileName[iprod].Contains("MB") || fileName[iprod].Contains("INT7"))
          RemovePointsOutOfRangeOrLargeErrorFromGraph(gFrac[ism][im02][iprod], 2., 10.); // PlotUtils.C
        else if( fileName[iprod].Contains("LHC11") && fileName[iprod].Contains("EMC7") )
          RemovePointsOutOfRangeOrLargeErrorFromGraph(gFrac[ism][im02][iprod], 6., 60.); // PlotUtils.C
        
        if(fileName[iprod].Contains("JJ"))
        {
          if(fileName[iprod].Contains("High"))
            RemovePointsOutOfRangeOrLargeErrorFromGraph(gFrac[ism][im02][iprod], 14., 60.); // PlotUtils.C
          else if(fileName[iprod].Contains("Low"))
            RemovePointsOutOfRangeOrLargeErrorFromGraph(gFrac[ism][im02][iprod],  7., 20.); // PlotUtils.C
        }
        
      } // im02
      
    } // ism
    
  } // iprod
  
  
  //
  // PLOTS
  //
  Float_t colorProd[] = { 1, 1,kBlue, kBlue,kBlue,kRed,kRed,kRed,  8,  8,  8};
  Float_t styleProd[] = {20,22,   26,    24,   25,   26,  24, 25, 26, 24, 25};  
  
  gStyle->SetPadRightMargin(0.01);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetTitleFontSize(0.06);
  //  gStyle->SetStatFontSize(0.06);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  gStyle->SetPadTopMargin(0.1);
  
  //
  // Plot: Each frame is a M02 cut, each file is per one SM, SM10 is the sum of all.
  //
  Int_t ncol = 2;
  Int_t nrow = 2;
  
  if(nShShMaxCut==1)
  {
    ncol=1;
    nrow=1;
  }
  
  for(Int_t ism = firstSM; ism < nSM; ism++)
  {    
    TCanvas * cGraphSM = new TCanvas(Form("cGraphSM_iana%d_SM%d",iana,ism),
                                     Form("iana %d, SM%d",iana, ism),
                                     ncol*2000,nrow*2000);
    
    cGraphSM->Divide(ncol,nrow);
    
    gPad->GetLogx();
    gPad->GetLogy();
    
    TLegend *lE = 0;
    if(nShShMaxCut == 3) 
    {
      lE = new TLegend(0.1,0.1,1,1);
      lE->SetTextSize(0.05);  
    }
    else         
    {
      lE = new TLegend(0.4,0.55,0.95,0.9);
      lE->SetTextSize(0.035);  
    }
    lE->SetFillColor(0);
    lE->SetFillStyle(0);
    lE->SetLineColor(0);
    lE->SetBorderSize(0);
    //lE->SetHeader(Form("0.10<#sigma^{2}_{long}<%1.2f, %s", maxShShCut[im02], anaTitle[iana].Data()));  
    
    if ( nShShMaxCut > 1 )
    {
      if (ism==10) lE->SetHeader(Form(" All SM, %s", anaTitle[iana].Data()));
      else         lE->SetHeader(Form("    SM %d, %s",ism, anaTitle[iana].Data()));
    }
    
    for(Int_t im02 = 0; im02 < nShShMaxCut; im02++)
    {    
      cGraphSM->cd(im02+1);
      
      gPad->GetLogx();
      gPad->GetLogy();
      
      TString title = Form("0.10<#sigma^{2}_{long}<%1.2f", maxShShCut[im02]);
      if ( nShShMaxCut == 1 )
      {
        if (ism == 10 ) title = Form("0.10<#sigma^{2}_{long}<%1.2f, all SM, %s", maxShShCut[im02],anaTitle[iana].Data());
        else            title = Form("0.10<#sigma^{2}_{long}<%1.2f, SM %d,  %s", maxShShCut[im02],ism, anaTitle[iana].Data());
      }
     
      hAxis[im02][ism] = new TH1F(Form("hAxisBis_SM%d_M02%d",ism,im02),title.Data(),1000,2,50);
      
      hAxis[im02][ism]->GetYaxis()->SetTitle("#it{R}_{#sigma_{long}}^{#it{n}^{w}_{cell}} (%)");
      if ( iana < 2 ) hAxis[im02][ism]->GetXaxis()->SetTitle("#it{E} (GeV)");
      else            hAxis[im02][ism]->GetXaxis()->SetTitle("#it{E}_{T} (GeV/#it{c})");
      
      if(ism!=3 && ism !=7 ) hAxis[im02][ism]->SetMaximum(45);
      else                   hAxis[im02][ism]->SetMaximum(65);
      
      hAxis[im02][ism]->SetMinimum(0);
      if( (ism == 3 || ism ==7) && im02 == 2 ) hAxis[im02][ism]->SetMaximum(45);
      
      hAxis[im02][ism]->Draw("");
      
      for(Int_t iprod = 0; iprod < nProd; iprod++)
      {
        if(!file[iprod]) continue;
        if(! gFrac[ism][im02][iprod] ) continue;
        
        gFrac[ism][im02][iprod]->SetMarkerColor(colorProd[iprod]);
        gFrac[ism][im02][iprod]->SetLineColor  (colorProd[iprod]);
        gFrac[ism][im02][iprod]->SetMarkerStyle(styleProd[iprod]);
        gFrac[ism][im02][iprod]->SetMarkerSize (5);
        if(im02==0) lE->AddEntry(gFrac[ism][im02][iprod],legName[iprod],"PL");
        
        gFrac[ism][im02][iprod]->Draw("P");
      }
      
      gFrac[ism][im02][0]->Draw("P");
    }
    
    if(nShShMaxCut == 3) cGraphSM->cd(4);
    lE->Draw();
    
    TString fileName = Form("figures/Comparison_Fraction_M02PhotonToAll_NCellCut_XaxisE_PerM02Cut_PerProd_%s_",
                            anaName[iana].Data());
    
    if(ism == 10 ) fileName+="AllSM";
    else           fileName+=Form("SM%d",ism);
    
    fileName+=fileFormat;
    cGraphSM->Print(fileName);
  }
  
  //
  // Plot: Each frame is a SM, each file is per one M02 cut
  //
  for(Int_t im02 = 0; im02 < nShShMaxCut; im02++)
  {    
    TCanvas * cGraphE = new TCanvas(Form("cGraphE_iana%d_%d",iana,im02),
                                    Form("iana %d, %2.2f",iana, maxShShCut[im02]),
                                    ncolSM*2000,nrowSM*2000);
    
    gPad->GetLogx();
    gPad->GetLogy();
    
    cGraphE->Divide(ncolSM,nrowSM);
    
    TLegend *lSM = new TLegend(-0.04,0.1,1,1);
    lSM->SetFillColor(0);
    lSM->SetFillStyle(0);
    lSM->SetLineColor(0);
    lSM->SetBorderSize(0);
    lSM->SetTextSize(0.06);  
    lSM->SetHeader(Form("    0.10<#sigma^{2}_{long}<%1.2f, %s", maxShShCut[im02], anaTitle[iana].Data()));  
    
    for(Int_t ism = 0; ism < nSM; ism++)
    {    
      cGraphE->cd(ism+1);
      
      gPad->GetLogx();
      gPad->GetLogy();
      
      TString title = Form("SM %d",ism);
      if(ism == 10) 
        title = "All SM";
      
      hAxis[im02][ism] = new TH1F(Form("hAxis_SM%d_M02%d",ism,im02),title.Data(),1000,2,50);
      
      hAxis[im02][ism]->GetYaxis()->SetTitle("#it{R}_{#sigma_{long}}^{#it{n}^{w}_{cell}} (%)");
      if ( iana < 2 ) hAxis[im02][ism]->GetXaxis()->SetTitle("#it{E} (GeV)");
      else           hAxis[im02][ism]->GetXaxis()->SetTitle("#it{E}_{T} (GeV/#it{c})");
      
      hAxis[im02][ism]->SetMaximum(35);
      hAxis[im02][ism]->SetMinimum(0);
      
      hAxis[im02][ism]->Draw("");
      
      for(Int_t iprod = 0; iprod < nProd; iprod++)
      {
        if(!file[iprod]) continue;
        if(! gFrac[ism][im02][iprod] ) continue;
        
        gFrac[ism][im02][iprod]->SetMarkerColor(colorProd[iprod]);
        gFrac[ism][im02][iprod]->SetLineColor  (colorProd[iprod]);
        gFrac[ism][im02][iprod]->SetMarkerStyle(styleProd[iprod]);
        gFrac[ism][im02][iprod]->SetMarkerSize (5);
        if(ism==0) lSM->AddEntry(gFrac[ism][im02][iprod],legName[iprod],"PL");
        
        gFrac[ism][im02][iprod]->Draw("P");
      }
      
      gFrac[ism][im02][0]->Draw("P");
    }
    
    cGraphE->cd(nrowSM*ncolSM);
    lSM->Draw();
    
    TString fileName = Form("figures/Comparison_Fraction_M02PhotonToAll_NCellCut_XaxisE_PerSM_PerProd_%s_M02Cut%d",
                            anaName[iana].Data(),im02);
    
    fileName+=fileFormat;
    cGraphE->Print(fileName);
  }
}

//------------------------------------------------------------------------------
/// Main steering method.
/// Execute here the different cluster selection options for not isolated, 
/// isolated or inclusive. Set also here the productions to be analyzed and 
/// call the methods doing the calculation of the fractions and previous plotting
/// and/or the comparison of the fractions of different productions.
///
/// \param doCalc    : Execute method CalculateAndPlot()
/// \param doComp    : Execute method CompareProductions()
/// \param tm        : String with track matching option used
/// \param plotRat   : Bool to activate ratio plots.
/// \param debug     : Bool to activate printfs.
//------------------------------------------------------------------------------
void CalculateShowerShapeFractionPerNCellCut
(
 Bool_t  doCalc  = kTRUE,
 Bool_t  doComp  = kTRUE,
 TString tm      = "_TMDep",
 Bool_t  plotRat = kFALSE,
 Bool_t  debug   = kFALSE
 )
{
  const Int_t nProd = 8; // Make sure number is equal to entries below, if not it crashes
  TString filePath[nProd];
  TString title   [nProd]; 
  TString legend  [nProd]; 
  
  if ( debug ) printf("N prod %d\n",nProd);
  
  Int_t nproditer = 0;
  
  filePath[nproditer++] = "data/module/TCardChannel3/LHC11cd_EMC7";
  title   [nproditer-1] = "LHC11cd_EMC7";     
  legend  [nproditer-1] = "data, LHC11c+d, EMC7";     
  filePath[nproditer++] = "data/module/TCardChannel3/LHC11cd_INT7";
  title   [nproditer-1] = "LHC11cd_INT7";
  legend  [nproditer-1] = "data, LHC11c+d, INT7";     

  filePath[nproditer++] = "simu/module/pp_7TeV_MB/TCardChannel_Mimic0_Scaled2_v3/ScaledMerged"; 
  title   [nproditer-1] = "default_MC_MB";
  legend  [nproditer-1] = "MC default, MB";     
  filePath[nproditer++] = "simu/module/pp_7TeV_MB/TCardChannel_Mimic10c_EcellCut_Scaled2_v3/ScaledMerged";
  legend  [nproditer-1] = "MC mimic, MB";     
  
  filePath[nproditer++] = "simu/module/pp_7TeV_JJ_Dec_GJ/TCardChannel_Mimic0_Scaled2_v3/ScaledMerged"; 
  title   [nproditer-1] = "Default_MC_JJ_DecLow";
  legend  [nproditer-1] = "MC default, #gammaJ-JJ, #it{p}^{EMC}_{T}>3.5 GeV/#it{c}";     
  filePath[nproditer++] = "simu/module/pp_7TeV_JJ_Dec_GJ/TCardChannel_Mimic10c_EcellCut_Scaled2_v3/ScaledMerged";
  title   [nproditer-1] = "Mimic_MC_JJ_DecLow";
  legend  [nproditer-1] = "MC mimic, #gammaJ-JJ, #it{p}^{EMC}_{T}>3.5 GeV/#it{c}";     

  filePath[nproditer++] = "simu/module/pp_7TeV_JJ_Dec_High_GJ/TCardChannel_Mimic0_Scaled2_v3/ScaledMerged"; 
  title   [nproditer-1] = "Default_MC_JJ_DecHigh";
  legend  [nproditer-1] = "MC default, #gammaJ-JJ, #it{p}^{EMC}_{T}>7 GeV/#it{c}";     
  filePath[nproditer++] = "simu/module/pp_7TeV_JJ_Dec_High_GJ/TCardChannel_Mimic10c_EcellCut_Scaled2_v3/ScaledMerged";
  title   [nproditer-1] = "Mimic_MC_JJ_DecHigh";
  legend  [nproditer-1] = "MC mimic, #gammaJ-JJ, #it{p}^{EMC}_{T}>7 GeV/#it{c}";     

  if ( doCalc )
  {
    for(Int_t iprod = 0; iprod < nProd; iprod++)
    {
      CalculateAndPlot(title[iprod],filePath[iprod],0,tm, plotRat, debug); // Not isolated
      CalculateAndPlot(title[iprod],filePath[iprod],1,tm, plotRat, debug); // Isolated
      CalculateAndPlot(title[iprod],filePath[iprod],2,tm, plotRat, debug); // Inclusive
    }
  }
  
  if(doComp)
  {
    CompareProductions(nProd,title,legend,0,debug);
    CompareProductions(nProd,title,legend,1,debug);
    CompareProductions(nProd,title,legend,2,debug);
  }
}
