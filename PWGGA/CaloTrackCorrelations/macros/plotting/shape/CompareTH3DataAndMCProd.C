///
/// \file CompareTH3DataAndMCProd.C
/// \ingroup CaloTrackCorrMacrosPlotting
/// \brief Compare a cluster parameter (shower shape) between different productions 
///
/// Compare different shape parameter histograms in data and MC productions
/// depending on the SM per cluster energy bis.
/// Treatment of the output of the class AliAnaClusterShapeCorrelStudies
/// Input are TH3 histograms where x=energy, y= SM number or T-Card index and z is one cluster parameter
/// See the different options in MakeDataMCComparisonPerSMClusterEbin.C
/// Here plots are made or projected distributions are stored on another file for further processing
/// like in CalculateParamChi2MCxTalkDataPerSM.C
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
#include "PlotUtils.C"
#endif

void DefineHistogramSettings
(TString histoName, 
 Float_t & xmin    , Float_t & xmax, 
 Float_t & xNormMin, Float_t & xNormMax,
 Bool_t & normtomax, Int_t   & rebin     );

Bool_t SetHistogramRangesAndTitles
( TH1D* histo, TString histoName, Int_t iprod,
 Int_t rebin, Float_t xmin, Float_t xmax, 
 Bool_t normtomax, Float_t xNormMin, Float_t xNormMax, 
 Float_t lowEBin, Float_t highEBin, Bool_t debug );


//------------------------------------------------------------------------------
/// Main method
///
/// \param outputFileName : String with part of output file name to be recovered and used by another macro, see  MakeDataMCComparisonPerSMClusterEbinAndChi2.C 
/// \param histoName : Name of histogram with the parameter to study = "SMM02NoCut", "SMM02NoCut","SMM02",
///                    "SMM20LowM02NoCut","SMM20LowM02","SMM20HighM02NoCut","SMM20HighM02","SMNCell" and other see  MakeDataMCComparisonPerSMClusterEbin.C
/// \param nProd     : Number of analyzed productions 
/// \param prod      : Strings array with input files location
/// \param prodLeg   : Strings array with title name of analized productions
/// \param clusterization : String with clusterizer used, in case it is included in the input data file name
/// \param tm        : String with track matching option used
/// \param pid       : String with PID cut used, see AliAnaClusterShapeAnalysis
/// \param titleMC   : Simplified string acronym of input MC
/// \param mcLeg     : Full title string of input MC
/// \param titleData : Simplified string acronym of input data
/// \param daLeg     : Full title string of input data
/// \param firstP    : First SM or T-Card index to be studied
/// \param lastP     : Last SM or T-Card index to be studied
/// \param binE      : TArrayD with list of energy bin limits
/// \param firstMC   : prod should contain first data and then MC arrays, this int indicates which one is the first MC
/// \param plotRatio : Optional plotting of different ratio histograms
/// \param saveHisto : Optional storing of projected histograms into output file, if active, no plots are done
/// \param opt       : Additional string to differenciate output histograms for different analysis 
/// \param debug     : Bool to activate printfs.
//------------------------------------------------------------------------------
void CompareTH3DataAndMCProd
(
 TString & outputFileName, 
 TString histoName = "SMM02NoCut",
 
 const Int_t nProd = 3,
 TString * prod    = 0x0,
 TString * prodLeg = 0x0,
 
 TString clusterization = "",//"_V1_Ecell100_Eseed500",
 TString tm        = "_TMDep",//"TMFix", 
 TString pid       = "_Neutral",//"_Electron", "_Hadron","Neutral"
 
 TString titleMC   = "JJDecLow", // "JJDecHigh", "MB", "JJDecLow"
 TString mcLeg     = "MC: #gammaJ+JJ({p^{EMCal}_{T,#gamma}>3.5 GeV/#it{c})",
 TString titleData = "LHC11cd_EMC7",
 TString daLeg     = "pp@7 TeV, LHC11c+d EMC7",
 
 Int_t   firstP    = 0,
 Int_t   lastP     = 9,
 TArrayD binE      = 0,
 Int_t   firstMC   = 1,

 Bool_t  plotRatio = kFALSE,
 Bool_t  saveHisto = kFALSE,
 TString opt       = "",
 Bool_t  debug     = kFALSE
)
{
  TString fileFormat = ".eps";
  
  if(histoName.Contains("CellModule"))        pid = "";
  if(histoName.Contains("SMEMaxEClusterRat")) pid = "";
  if(histoName.Contains("TCardChannel"))      pid = "";

  Float_t xmin      = -1;
  Float_t xmax      = 10000; 
  Float_t xNormMin  = -1;
  Float_t xNormMax  = 10000;
  Int_t   rebin     = 1;
  Bool_t  normtomax = kFALSE;
  DefineHistogramSettings(histoName, xmin, xmax,  xNormMin, xNormMax, normtomax, rebin);
  
  if ( debug ) printf("Histo: %s, range [%2.2f,%2.2f], Norm %d, range [%2.2f, %2.2f], rebin %d, param range [%d,%d]\n",
                       histoName.Data(), xmin,xmax, normtomax, xNormMin,xNormMax,rebin, firstP,lastP);
  
  // Energy bins and productions
  //
  const Int_t nEBins = binE.GetSize()-1;
  Int_t ncolE = 5;
  Int_t nrowE = 5;
  GetCanvasColRowNumber(nEBins,ncolE,nrowE); // PlotUtils.C
  
  Double_t lowE[nEBins];
  Double_t higE[nEBins];

  // SM or TC bins
  //
  const Int_t nparam = lastP-firstP+1;
  Int_t ncolPar = 5;
  Int_t nrowPar = 5;
  GetCanvasColRowNumber(nparam,ncolPar,nrowPar); // PlotUtils.C

  // Merge good/bad SMs or TC
  //
  const Int_t nGroup = 4;
  Int_t ngroupmax = 3;
  if(histoName.Contains("TCardChannel"))
    ngroupmax = 4;
 // printf("N GROUP MAX %d\n",ngroupmax);
  
  TString groupSM[] = {"0,4,5,6,8,9","1,2","3,7"};
  TString groupTC[] = {"Border","Border+1","Center+1","Center"};
 
  
  // Histogram/file arrays
  //
  TH3F  * h3    [nProd];
  TH1D  * h     [nProd][nEBins][nparam];
  TH1D  * hA    [nProd][nEBins];
  TH1D  * hGroup[nProd][nEBins][nGroup];
 
  TFile * file[nProd];

  //
  //Open the histogram files and get the TH3 histogram
  //
  // Data
  //
  for(Int_t iprod = 0; iprod < firstMC; iprod++)
  {
    file[iprod] = TFile::Open(Form("%s%s.root",prod[iprod].Data(),clusterization.Data()));
    if ( debug ) printf("Read: %s%s.root, %p\n",prod[iprod].Data(),clusterization.Data(),file[iprod]);
    
    if(!file[iprod])
    {
      printf("Missing file %d\n",iprod);
      continue;
    }
    
    h3[iprod] = (TH3F*) file[iprod]->Get(Form("Shape%s_h%s%s"  ,tm.Data(),histoName.Data(),pid.Data()));
    if ( debug ) printf("Data, Shape%s_h%s%s %p\n",tm.Data(),histoName.Data(),pid.Data(),h3[iprod]);
    if ( !h3[iprod] ) return;
    
    h3[iprod]->Sumw2();
    
    //prod[iprod].ReplaceAll("/","_");
   
    if ( debug ) printf("Entries %e Integral %e\n",h3[iprod]->GetEntries(),h3[iprod]->Integral());
  }
  
  // Simu
  //
  for(Int_t iprod = firstMC; iprod < nProd; iprod++)
  {
    if ( !prod[iprod].Contains("root") ) 
      file[iprod] = TFile::Open(Form("%s/ScaledMerged%s.root",prod[iprod].Data(),clusterization.Data()));
    else                                 
      file[iprod] = TFile::Open(Form("%s",prod[iprod].Data()));
    
    if ( debug ) printf("Simu %d file %s, Shape_%s_h%s%s, %p\n",
                        iprod, prod[iprod].Data(), tm.Data(),histoName.Data(),pid.Data(), file[iprod]);

    if(!file[iprod]) 
    {
      h3[iprod] = 0;
      continue;
    }
    
    h3[iprod] = (TH3F*) file[iprod]->Get(Form("Shape%s_h%s%s"  ,tm.Data(),histoName.Data(), pid.Data()));
    if ( debug ) printf("\t histo %p\n",h3[iprod]);
    
    if ( !h3[iprod] ) return;
    
    if ( debug ) printf("Entries %e Integral %e\n",h3[iprod]->GetEntries(),h3[iprod]->Integral());

//    prod[iprod].ReplaceAll("/","_");
//    prod[iprod].ReplaceAll("_sum","");    
  }
    
  //
  // Project
  //
  if ( debug ) printf("Project per E bin, SM/TC bin\n");

  Double_t width = 0;
  for(Int_t iprod = 0; iprod < nProd; iprod++)
  {
    if ( !file[iprod] ) 
      continue;
    
    if ( debug ) printf("iprod %d\n",iprod);
        
    for(Int_t iebin = 0; iebin < nEBins; iebin++)
    {
      if ( debug ) printf("\t ebin %d\n",iebin);
      h3[iprod]->SetAxisRange(-1000,1000,"Y");

      Int_t ebinMin = h3[iprod]->GetXaxis()->FindBin(binE[iebin  ]);
      Int_t ebinMax = h3[iprod]->GetXaxis()->FindBin(binE[iebin+1])-1;
      //if(ebinMin==ebinMax) ebinMax = 1000;
      lowE[iebin]  = h3[iprod]->GetXaxis()->GetBinLowEdge(ebinMin);
      width        = h3[iprod]->GetXaxis()->GetBinWidth  (ebinMax);
      higE[iebin]  = h3[iprod]->GetXaxis()->GetBinLowEdge(ebinMax)+width;
      //if(lowE[iebin] == higE[iebin]) higE[iebin] = 100;
      if ( debug && iprod==0 )
      {
        printf("\t iE %d, ebinMin %d, emin %2.1f, low edge %2.1f, ebinMax %d, emax %2.1f, high edge %2.1f, width %2.2f\n",
               iebin, ebinMin, binE[iebin], lowE[iebin], ebinMax, binE[iebin+1], higE[iebin], width);
      }
            
      for(Int_t ipbin = firstP; ipbin <= lastP; ipbin++)
      {        
        Int_t pbinMax = h3[iprod]->GetYaxis()->FindBin(ipbin);
        Int_t pbinMin = h3[iprod]->GetYaxis()->FindBin(ipbin);
  
        if ( debug ) printf("\t \t ipbin %d - first %d =  %d, last %d, total %d\n",ipbin,firstP,ipbin-firstP, lastP, nparam);
        
        h[iprod][iebin][ipbin-firstP] = 
        (TH1D*) h3[iprod]->ProjectionZ(Form("%sProd%d_BinE%d_BinP%d",histoName.Data(),iprod,iebin,ipbin),ebinMin,ebinMax, pbinMin,pbinMax);
        
        Bool_t ok = SetHistogramRangesAndTitles
        (h[iprod][iebin][ipbin-firstP],histoName, iprod,
         rebin, xmin,xmax, 
         normtomax, xNormMin, xNormMax, lowE[iebin], higE[iebin],debug);
        
        if ( !ok )
          continue;
      }
      
      //
      // All param
      //
      if ( debug ) printf("\t \t Merge all SM\n");
      {
        hA[iprod][iebin] = 
        (TH1D*) h3[iprod]->ProjectionZ(Form("%sProd%d_BinE%d_All",histoName.Data(),iprod,iebin),ebinMin,ebinMax, -1,10000);

        Bool_t ok = SetHistogramRangesAndTitles
        (hA[iprod][iebin],histoName, iprod,
         rebin, xmin,xmax, 
         normtomax, xNormMin, xNormMax,  lowE[iebin],higE[iebin],debug);
        
        if ( !ok )
          printf("Histogram with sum of all not available\n");
      }
      
      if ( debug ) printf("\t \t Merge groups\n");

      //
      // Merged SM
      //  
      Float_t colorGroup[] = {kBlue,kViolet,kRed,kBlack,6,8,3,5}; 

      if ( histoName.Contains("SM") )
      {        
//        printf("\t \t ie %d, iprod %d\n",iebin,iprod);
//        printf("\t \t 0-%p, 1-%p, 2-%p, 3-%p, 4-%p\n",h[iprod][iebin][0],h[iprod][iebin][1], h[iprod][iebin][2],  h[iprod][iebin][3], h[iprod][iebin][4]);
//        printf("\t \t 5-%p, 6-%p, 7-%p, 8-%p, 9-%p\n",h[iprod][iebin][5],h[iprod][iebin][6], h[iprod][iebin][7],  h[iprod][iebin][8], h[iprod][iebin][9]);
        
        hGroup[iprod][iebin][0] = 0;
        hGroup[iprod][iebin][1] = 0;
        hGroup[iprod][iebin][2] = 0;
        
        if ( h[iprod][iebin][0] ) hGroup[iprod][iebin][0] = (TH1D*) h[iprod][iebin][0]->Clone(Form("SM_Group0_Ebin%d_Prod%d",iebin,iprod));
        if ( h[iprod][iebin][1] ) hGroup[iprod][iebin][1] = (TH1D*) h[iprod][iebin][1]->Clone(Form("SM_Group1_Ebin%d_Prod%d",iebin,iprod));
        if ( h[iprod][iebin][3] ) hGroup[iprod][iebin][2] = (TH1D*) h[iprod][iebin][3]->Clone(Form("SM_Group2_Ebin%d_Prod%d",iebin,iprod));
                
        if( hGroup[iprod][iebin][0] )
        {
          hGroup[iprod][iebin][0]->Add( h[iprod][iebin][4] );
          hGroup[iprod][iebin][0]->Add( h[iprod][iebin][5] );
          hGroup[iprod][iebin][0]->Add( h[iprod][iebin][6] );
          hGroup[iprod][iebin][0]->Add( h[iprod][iebin][8] );
          hGroup[iprod][iebin][0]->Add( h[iprod][iebin][9] );
          
          hGroup[iprod][iebin][0]->Scale(1./6.);
        }
        
        if( hGroup[iprod][iebin][1] )
        {  
          hGroup[iprod][iebin][1]->Add( h[iprod][iebin][2] );
          
          hGroup[iprod][iebin][1]->Scale(1./2.);            
        }
        
        if( hGroup[iprod][iebin][2] )
        {  
          hGroup[iprod][iebin][2]->Add( h[iprod][iebin][7] );
          
          hGroup[iprod][iebin][2]->Scale(1./2.);
        }
      } // SM
      
      if ( histoName.Contains("TCardChannel") )
      {        
        hGroup[iprod][iebin][0] = 0;
        hGroup[iprod][iebin][1] = 0;
        hGroup[iprod][iebin][2] = 0;
        hGroup[iprod][iebin][3] = 0;
        
        if ( h[iprod][iebin][0] ) hGroup[iprod][iebin][0] = (TH1D*) h[iprod][iebin][0]->Clone(Form("TC_Group0_Ebin%d_Prod%d",iebin,iprod));
        if ( h[iprod][iebin][1] ) hGroup[iprod][iebin][1] = (TH1D*) h[iprod][iebin][1]->Clone(Form("TC_Group1_Ebin%d_Prod%d",iebin,iprod));
        if ( h[iprod][iebin][2] ) hGroup[iprod][iebin][2] = (TH1D*) h[iprod][iebin][2]->Clone(Form("TC_Group2_Ebin%d_Prod%d",iebin,iprod));
        if ( h[iprod][iebin][3] ) hGroup[iprod][iebin][3] = (TH1D*) h[iprod][iebin][3]->Clone(Form("TC_Group3_Ebin%d_Prod%d",iebin,iprod));
        
        if ( debug )
        {
          printf("ie %d, iprod %d\n",iebin,iprod);
          printf(" 0-%p,  1-%p,  2-%p,  3-%p,  4-%p\n",h[iprod][iebin][ 0],h[iprod][iebin][ 1], h[iprod][iebin][ 2],  h[iprod][iebin][ 3], h[iprod][iebin][ 4]);
          printf(" 5-%p,  6-%p,  7-%p,  8-%p,  9-%p\n",h[iprod][iebin][ 5],h[iprod][iebin][ 6], h[iprod][iebin][ 7],  h[iprod][iebin][ 8], h[iprod][iebin][ 9]);
          printf("10-%p, 11-%p, 12-%p, 13-%p, 14-%p\n",h[iprod][iebin][10],h[iprod][iebin][11], h[iprod][iebin][12],  h[iprod][iebin][13], h[iprod][iebin][14]);
          printf("15-%p, 16-%p \n"                    ,h[iprod][iebin][15],h[iprod][iebin][16]);
        }
        
        if( hGroup[iprod][iebin][0] )
        {
          hGroup[iprod][iebin][0]->Add( h[iprod][iebin][ 7] );
          hGroup[iprod][iebin][0]->Add( h[iprod][iebin][ 8] );
          hGroup[iprod][iebin][0]->Add( h[iprod][iebin][15] );
          
          hGroup[iprod][iebin][0]->Scale(1./4.);
        }
        
        if( hGroup[iprod][iebin][1] )
        {  
          hGroup[iprod][iebin][1]->Add( h[iprod][iebin][ 6] );
          hGroup[iprod][iebin][1]->Add( h[iprod][iebin][ 9] );
          hGroup[iprod][iebin][1]->Add( h[iprod][iebin][14] );
          
          hGroup[iprod][iebin][1]->Scale(1./4.);            
        }

        if( hGroup[iprod][iebin][2] )
        {  
          hGroup[iprod][iebin][2]->Add( h[iprod][iebin][ 5] );
          hGroup[iprod][iebin][2]->Add( h[iprod][iebin][10] );
          hGroup[iprod][iebin][2]->Add( h[iprod][iebin][13] );
          
          hGroup[iprod][iebin][2]->Scale(1./4.);
        }

        if( hGroup[iprod][iebin][3] )
        {  
          hGroup[iprod][iebin][3]->Add( h[iprod][iebin][ 4] );
          hGroup[iprod][iebin][3]->Add( h[iprod][iebin][11] );
          hGroup[iprod][iebin][3]->Add( h[iprod][iebin][12] );
          
          hGroup[iprod][iebin][3]->Scale(1./4.);
        }

      } // TCard

      for(Int_t igroup = 0; igroup < ngroupmax; igroup++)
      {
        if( !hGroup[iprod][iebin][igroup] ) continue;
        hGroup[iprod][iebin][igroup]->SetLineColor  (colorGroup[igroup]);
        hGroup[iprod][iebin][igroup]->SetMarkerColor(colorGroup[igroup]);
      }
      
    } // ebins
  } // prod bins
  
  //
  // Save the projected histograms in file for further processing
  // In such case, do not do plotting.
  //
  TString outputFileName2 = "";//Form("%s%s%s",clusterization.Data(),pid.Data(),tm.Data());
  //printf("Opt name %s\n",opt.Data());
  outputFileName= Form("%s_%s%s%s",titleData.Data(),titleMC.Data(),outputFileName2.Data(),opt.Data());

  if ( debug )
    printf("Output string: <%s>",outputFileName.Data());
  
  if ( saveHisto )
  {
    TFile* fout = new TFile(Form("figures/%s/Projections_%s.root",
                                 histoName.Data(),outputFileName.Data()),"recreate");
    
    for(Int_t iprod = 0; iprod < nProd; iprod++)
    {
      for(Int_t iebin= 0; iebin < nEBins; iebin++)
      {
        for(Int_t ipbin = firstP; ipbin <= lastP; ipbin++)
        {
          if( !h[iprod][iebin][ipbin-firstP] ) continue;
          h[iprod][iebin][ipbin-firstP]->SetName(Form("Prod%d_Ebin%d_Param%d",iprod,iebin,ipbin));
          h[iprod][iebin][ipbin-firstP]->Write();
        }
        
        for(Int_t igroup = 0; igroup < ngroupmax; igroup++)
        {
          if( !hGroup[iprod][iebin][igroup] ) continue;
          hGroup[iprod][iebin][igroup]->SetName(Form("Prod%d_Ebin%d_ParamGroup%d",iprod,iebin,igroup));
          hGroup[iprod][iebin][igroup]->Write();
        }
        
        if ( hA[iprod][iebin] )
        {
          hA[iprod][iebin]->SetName(Form("Prod%d_Ebin%d_AllParam",iprod,iebin));
          hA[iprod][iebin]->Write();
        }
        
      } // iebin
    } // iprod
    fout->Close();
    printf("Create: %s\n",fout->GetName());
    // do no more
    return;
  }
  
  // 
  // Plot per E bin per Param: SM, TC, ...
  //
  
  TString fileName ;
  titleData.ReplaceAll("/","_");
  
  gStyle->SetPadRightMargin(0.01);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetTitleFontSize(0.06);
  //  gStyle->SetStatFontSize(0.06);
  gStyle->SetOptStat(0);

  //---------------------------------------  
  // One file per E bin, each canvas one SM/TC
  //---------------------------------------
  for(Int_t iebin = 0; iebin < nEBins; iebin++)
  {
    gStyle->SetOptTitle(1);
    gStyle->SetPadTopMargin(0.1);
    
    TCanvas * c = new TCanvas(Form("c_ebin%d_%s",
                                   iebin,histoName.Data()),
                              Form("%2.1f < E < %2.1f, %s",
                                   lowE[iebin],higE[iebin], 
                                   histoName.Data()),
                              ncolPar*2000,nrowPar*2000);
    c->Divide(ncolPar,nrowPar);
    
    TLegend *l = new TLegend(-0.04,0,1,1);
    if(ncolPar!=nrowPar) 
    {
      l = new TLegend(-0.04,0.,1,1);
      l->SetTextSize(0.06);
    }
    else  
    {
      l = new TLegend(0.35,0.6,0.9,0.9);
      l->SetTextSize(0.045);
    }
    l->SetFillColor(0);
    l->SetFillStyle(0);
    l->SetLineColor(0);
    l->SetBorderSize(0);
    l->AddEntry("",daLeg,"");
    l->AddEntry("",mcLeg,"");
    if((histoName.Contains("MM02") || histoName.Contains("MM20")) && !histoName.Contains("NoCut") )
      l->SetHeader(Form("   %2.1f < #it{E} < %2.1f GeV, #it{n}^{w}_{cells} > 4",lowE[iebin],higE[iebin]));
    else
      l->SetHeader(Form("   %2.1f < #it{E} < %2.1f GeV",lowE[iebin],higE[iebin]));
    
    for(Int_t ipbin = firstP; ipbin <= lastP; ipbin++)
    {
      c->cd(ipbin+1);
      
      if ( (histoName.Contains("ECell") && !histoName.Contains("Tot")) || histoName.Contains("NLocMax"))gPad->SetLogy();
      gPad->SetGridy();
      
      if(!h[0][iebin][ipbin-firstP]) continue;
      
      for(Int_t iprod = 0; iprod < nProd; iprod++)
      {
        if(!h[iprod][iebin][ipbin-firstP]) continue;
        
        h[iprod][iebin][ipbin-firstP]->SetTitle(Form("SM %d",ipbin));
        
        if(iprod==0) h[iprod][iebin][ipbin-firstP]->Draw("H");
        else         h[iprod][iebin][ipbin-firstP]->Draw("H same");
        
        if(h[iprod][iebin][ipbin-firstP]->GetMaximum() > h[0][iebin][ipbin-firstP]->GetMaximum())
          h[0][iebin][ipbin-firstP]->SetMaximum(h[iprod][iebin][ipbin-firstP]->GetMaximum()*1.2);
        
        if(ipbin==0)
          l->AddEntry(h[iprod][iebin][ipbin-firstP],Form("%s",prodLeg[iprod].Data()),"PL");
      } // iprod
      
      h[0][iebin][ipbin-firstP]->Draw("H same");

    } // param
    
    if(ncolPar!=nrowPar && ncolPar*nrowPar!=nparam)
      c->cd(nparam+1);
    l->Draw();
    
    fileName = Form("figures/%s/ProductionComparisonPerSM_%s_Ebin%d",
                    histoName.Data(),outputFileName.Data(),iebin);
    fileName+=fileFormat;
    c->Print(fileName);
    
    if ( !plotRatio ) continue;
    
    TCanvas * cR = new TCanvas(Form("cRat_ebin_%d_%s",
                                    iebin,histoName.Data()),
                               Form("%2.1f < E < %2.1f, %s",
                                    lowE[iebin],higE[iebin], 
                                    histoName.Data()),
                               ncolPar*2000,nrowPar*2000);
    cR->Divide(ncolPar,nrowPar);
    
    
    for(Int_t ipbin = firstP; ipbin <= lastP; ipbin++)
    {
      cR->cd(ipbin+1);
      
      gPad->SetGridy();
      
      if(!h[1][iebin][ipbin-firstP]) continue;
      
      for(Int_t iprod = 1; iprod < nProd; iprod++)
      {
        if(!h[iprod][iebin][ipbin-firstP]) continue;
        
        TH1F* hRat = (TH1F*) h[iprod][iebin][ipbin-firstP]->Clone(Form("%s_Ratio",h[iprod][iebin][ipbin-firstP]->GetName()));
        hRat->Divide(h[0][iebin][ipbin-firstP]);
        
        if(iprod==0) hRat->Draw("H");
        else         hRat->Draw("H same");
        
        hRat->SetYTitle("Ratio MC to Data");
        hRat->SetMaximum(2);
        hRat->SetMinimum(0);
        
      } // iprod
      
    } // param
    
    cR->cd(ncolPar*nrowPar);
    l->Draw();
    
    fileName = Form("figures/%s/ProductionComparisonPerSM_%s_Ebin%d_Ratios",
                    histoName.Data(),outputFileName.Data(),iebin);
    fileName+=fileFormat;
    cR->Print(fileName);
    
  } // e bin
  
  //---------------------------------------  
  // All SM, each canvas one energy bin
  //---------------------------------------
  TLegend *l2 = 0;
  if(nEBins > 2 && ncolE!=nrowE) 
  {
    l2 = new TLegend(-0.04,0.,1,1);
    l2->SetTextSize(0.06);
  }
  else  
  {
    l2 = new TLegend(0.35,0.7,0.9,0.9);
    l2->SetTextSize(0.025);
  }
  l2->SetFillColor(0);
  l2->SetFillStyle(0);
  l2->SetLineColor(0);
  l2->SetBorderSize(0);
  l2->AddEntry("",daLeg,"");
  l2->AddEntry("",mcLeg,"");
  
  if(ncolE==nrowE) l2->SetTextSize(0.04);
  
  TCanvas * cA = new TCanvas(Form("cAllSM_%s",histoName.Data()),
                             Form("All SM, %s",histoName.Data()),
                             ncolE*2000,nrowE*2000);
  
  cA->Divide(ncolE,nrowE);
  
  for(Int_t iebin = 0; iebin < nEBins; iebin++)
  {
    cA->cd(iebin+1);
    
    if(!hA[0][iebin]) continue;
    
    if((histoName.Contains("ECell") && !histoName.Contains("Tot")) || histoName.Contains("NLocMax"))
      gPad->SetLogy();
    
    gStyle->SetOptTitle(1);
    gStyle->SetPadTopMargin(0.1);
    
    hA[0][iebin]->Draw("H");
    
    for(Int_t iprod = 0; iprod < nProd; iprod++)
    {
      if(!hA[iprod][iebin]) continue; 
      
      hA[iprod][iebin]->SetTitle(Form("%2.1f < #it{E}^{clus} < %2.1f GeV",lowE[iebin],higE[iebin]));
      
      if((histoName.Contains("MM02") || histoName.Contains("MM20")) && !histoName.Contains("NoCut") )
      {
        hA[iprod][iebin]->SetTitle(Form("%2.1f < #it{E}^{clus} < %2.1f GeV, #it{n}^{w}_{cells} > 4",lowE[iebin],higE[iebin]));
      }
      
      hA[iprod][iebin]->Draw("H same");
      
      if(hA[iprod][iebin]->GetMaximum() > hA[0][iebin]->GetMaximum())
        hA[0][iebin]->SetMaximum(hA[iprod][iebin]->GetMaximum()*1.2);
      
      if(hA[iprod][iebin] && iebin==0)  l2->AddEntry(hA[iprod][iebin],Form("%s",prodLeg[iprod].Data()),"L");      
    }
    
    hA[0][iebin]->Draw("H same");
  }
  
  if ( nEBins > 1 && ncolE!=nrowE) 
    cA->cd(ncolE*nrowE);
  
  l2->Draw("same");
  
  fileName = Form("figures/%s/ProductionComparisonPerEbin_%s_AllSM",
                  histoName.Data(),outputFileName.Data());
  fileName+=fileFormat;
  cA->Print(fileName);
  
  if(plotRatio)
  {
    TCanvas * cAR = new TCanvas(Form("cAllSM_Rat_%s",histoName.Data()),
                                Form("All SM, %s",histoName.Data()),
                                ncolE*2000,nrowE*2000);
    
    cAR->Divide(ncolE,nrowE);
    
    for(Int_t iebin = 0; iebin < nEBins; iebin++)
    {
      cAR->cd(iebin+1);
      
      if(!hA[0][iebin]) continue;
      
      //if(histoName.Contains("ECell"))gPad->SetLogy();
      
      gStyle->SetOptTitle(1);
      gStyle->SetPadTopMargin(0.1);
      
      hA[0][iebin]->Draw("H");
      
      for(Int_t iprod = 1; iprod < nProd; iprod++)
      {
        if(!hA[iprod][iebin]) continue; 
        
        TH1F* hRat = (TH1F*) hA[iprod][iebin]->Clone(Form("%s_Ratio",hA[iprod][iebin]->GetName()));
        hRat->Divide(hA[0][iebin]);
        
        if(iprod==1) hRat->Draw("H");
        else         hRat->Draw("H same");
        
        hRat->SetYTitle("Ratio MC to Data");
        hRat->SetMaximum(2);
        hRat->SetMinimum(0);
      }
    }
    
    if ( nEBins > 1 && ncolE!=nrowE) 
      cAR->cd(ncolE*nrowE);
    
    l2->Draw("same");
    
    fileName = Form("figures/%s/ProductionComparisonPerEbin_%s_AllSM_Ratios",
                    histoName.Data(),outputFileName.Data());
    fileName+=fileFormat;
    cAR->Print(fileName);
  }
  
  //---------------------------------------  
  // Each file one SM, each canvas one energy bin
  //---------------------------------------
  for(Int_t ipbin = firstP; ipbin <= lastP; ipbin++)
  {
    gStyle->SetOptTitle(1);
    gStyle->SetPadTopMargin(0.1);
    
    TCanvas * c = new TCanvas(Form("c_par%d_%s",
                                   ipbin,histoName.Data()),
                              Form("param %d, %s",
                                   ipbin, 
                                   histoName.Data()),
                              3*2000,2*2000);
    c->Divide(3,2);
    
    TLegend *l = new TLegend(-0.04,0,1,1);
    l->SetFillColor(0);
    l->SetFillStyle(0);
    l->SetLineColor(0);
    l->SetBorderSize(0);
    l->SetTextSize(0.07);
    l->SetHeader(Form("    SM %d",ipbin));
    for(Int_t iebin = 0; iebin < nEBins; iebin++)
    {
      c->cd(iebin+1);
      
      if(histoName.Contains("ECell") && !histoName.Contains("Tot"))gPad->SetLogy();
      gPad->SetGridy();
      
      if(!h[0][iebin][ipbin-firstP]) continue;
      
      for(Int_t iprod = 0; iprod < nProd; iprod++)
      {
        if(!h[iprod][iebin][ipbin-firstP]) continue;
        
        h[iprod][iebin][ipbin-firstP]->SetTitle(Form("%2.1f < #it{E}^{clus} < %2.1f GeV",lowE[iebin],higE[iebin]));
        if((histoName.Contains("MM02") || histoName.Contains("MM20")) && !histoName.Contains("NoCut") )
          h[iprod][iebin][ipbin-firstP]->SetTitle(Form("%2.1f < #it{E}^{clus} < %2.1f GeV, #it{n}^{w}_{cells} > 4",lowE[iebin],higE[iebin]));
        
        if(iprod==0) h[iprod][iebin][ipbin-firstP]->Draw("H");
        else         h[iprod][iebin][ipbin-firstP]->Draw("H same");
        
        if(h[iprod][iebin][ipbin-firstP]->GetMaximum() > h[0][iebin][ipbin-firstP]->GetMaximum())
          h[0][iebin][ipbin-firstP]->SetMaximum(h[iprod][iebin][ipbin-firstP]->GetMaximum()*1.2);
        
        if(iebin==0)
          l->AddEntry(h[iprod][iebin][ipbin-firstP],Form("%s",prodLeg[iprod].Data()),"PL");
      } // iprod
      
      h[0][iebin][ipbin-firstP]->Draw("H same");
    } // param
    
    c->cd(6);
    l->Draw();
    
    fileName = Form("figures/%s/ProductionComparisonPerEbin_%s_SM%d",
                    histoName.Data(),outputFileName.Data(),ipbin);
    fileName+=fileFormat;
    c->Print(fileName);
    
  } // param
  
  Int_t colorTC    [] = {1,2,6,4,4,6,2,1,1,2,6,4,4,6,2,1};
  Int_t lineStyleTC[] = {1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4};
   
  Int_t colorSM    []= { 1, 1, 2, 2, 3, 3, 4, 4, 7, 7, 6, 6, 8, 8, kOrange, kOrange, kRose, kRose, kMagenta, kMagenta};
  Int_t lineStyleSM[]= { 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2,       1,       2,      1,     2,       1,        2};

  //Int_t markerStyleSM[]={24,25,25,24,25,24,25,24,25,24,25,21,21,21,21,21,22,26,22,26,22,26,22,26};
  
  //---------------------------------------  
  // Each file per production Data/MC, 
  // each canvas one energy bin, each pad several param
  //---------------------------------------
  for(Int_t iprod = 0; iprod < nProd; iprod++)
  {
    gStyle->SetOptTitle(1);
    gStyle->SetPadTopMargin(0.1);
    
    TCanvas * c = new TCanvas(Form("c_iprod%d_%s",iprod,histoName.Data()),
                              Form("%s, %s"      ,prod[iprod].Data(), histoName.Data()),
                              ncolE*2000,nrowE*2000);
    c->Divide(ncolE,nrowE);
    
    TLegend *l = new TLegend(-0.04,0,1,1);
    l->SetFillColor(0);
    l->SetFillStyle(0);
    l->SetLineColor(0);
    l->SetBorderSize(0);
    l->SetTextSize(0.07);
    l->SetHeader(Form("    %s",prodLeg[iprod].Data()));
    for(Int_t iebin = 0; iebin < nEBins; iebin++)
    {
      c->cd(iebin+1);
      gPad->SetGridy();
      
      if(!h[iprod][iebin][0]) continue;
      
      for(Int_t iparam = firstP; iparam <= lastP; iparam++)
      {
        if(!h[iprod][iebin][iparam]) continue;
        
        h[iprod][iebin][iparam]->SetTitle(Form("%2.1f < #it{E} < %2.1f GeV",lowE[iebin],higE[iebin]));
        if((histoName.Contains("MM02") || histoName.Contains("MM20")) && !histoName.Contains("NoCut") )
          h[iprod][iebin][iparam]->SetTitle(Form("%2.1f < #it{E}^{clus} < %2.1f GeV, #it{n}^{w}_{cells} > 4",lowE[iebin],higE[iebin]));
                
        if(histoName.Contains("TCardChan"))
        {
          h[iprod][iebin][iparam]->SetLineColor(colorTC[iparam]);
          h[iprod][iebin][iparam]->SetLineStyle(lineStyleTC[iparam]);
        }
        else
        {
          h[iprod][iebin][iparam]->SetLineColor(colorSM[iparam]);
          h[iprod][iebin][iparam]->SetLineStyle(lineStyleSM[iparam]); 
        }
        
        if(iparam==0) h[iprod][iebin][iparam]->Draw("H");
        else          h[iprod][iebin][iparam]->Draw("H same");
        
        if(h[0][iebin][iparam] && h[iprod][iebin][iparam]->GetMaximum() > h[0][iebin][iparam]->GetMaximum())
          h[0][iebin][iparam]->SetMaximum(h[iprod][iebin][iparam]->GetMaximum()*1.2);
        
        if(iebin==0)
        {
          if(histoName.Contains("SM"))
            l->AddEntry(h[iprod][iebin][iparam],Form("SM = %d",iparam),"PL");
          else             
            l->AddEntry(h[iprod][iebin][iparam],Form("T-Card i = %d",iparam),"PL");
        }
      } // param
    }//iebin  
    c->cd(ncolE*nrowE);
    l->Draw();
    
//    fileName = Form("figures/%s/Comparison_PerParam_%s%s",
//                    histoName.Data(),prod[iprod].Data(),outputFileName2.Data());
    
    if ( prod[iprod].Contains("data") )
      fileName = Form("figures/%s/SMComparison_PerEbin_%s%s",
                    histoName.Data(),titleData.Data(),outputFileName2.Data());
    else
    {
      fileName = Form("figures/%s/SMComparison_PerEbin_%s_%d%s",
                      histoName.Data(),titleMC.Data(),iprod,outputFileName2.Data());
    }
    fileName+=fileFormat;
    c->Print(fileName);
    
    //---------------------------------------  
    // Each file, one production MC/Data
    // Grouped plots, each canvas one energy bin, each pad different groups
    //---------------------------------------
    {
      gStyle->SetOptTitle(1);
      gStyle->SetPadTopMargin(0.1);
      TCanvas * cGroup = new TCanvas(Form("cGroups_%s_%s",histoName.Data(),prod[iprod].Data()),
                                     Form("Groups, %s %s",histoName.Data(),prod[iprod].Data()),
                                     ncolE*2000,nrowE*2000);
      
      cGroup->Divide(ncolE,nrowE);
      
      TLegend *lG = new TLegend(0,0,1,1);
      lG->SetFillColor(0);
      lG->SetFillStyle(0);
      lG->SetLineColor(0);
      lG->SetBorderSize(0);
      lG->SetTextSize(0.07);
      lG->SetHeader(Form("    %s",prodLeg[iprod].Data()));
      

      for(Int_t iebin = 0; iebin < nEBins; iebin++)
      {
        cGroup->cd(iebin+1);
        for(Int_t igroup = 0; igroup < ngroupmax; igroup++)
        {
          if(!hGroup[iprod][iebin][igroup]) continue;
      
          hGroup[iprod][iebin][igroup]->SetTitle(Form("%2.1f < #it{E} < %2.1f GeV",lowE[iebin],higE[iebin]));
          if((histoName.Contains("MM02") || histoName.Contains("MM20")) && !histoName.Contains("NoCut") )
            hGroup[iprod][iebin][igroup]->SetTitle(Form("%2.1f < #it{E}^{clus} < %2.1f GeV, #it{n}^{w}_{cells} > 4",lowE[iebin],higE[iebin]));
          
          if(igroup==0) hGroup[iprod][iebin][igroup]->Draw("H");
          else          hGroup[iprod][iebin][igroup]->Draw("H same");
          
          if(hGroup[0][iebin][igroup] && hGroup[iprod][iebin][igroup]->GetMaximum() > hGroup[0][iebin][igroup]->GetMaximum())
            hGroup[0][iebin][igroup]->SetMaximum(hGroup[iprod][iebin][igroup]->GetMaximum()*1.2);
          
          if(iebin==0)
          {
            if(histoName.Contains("SM"))lG->AddEntry(hGroup[iprod][iebin][igroup],Form("SM = %s",groupSM[igroup].Data()),"PL");
            if(histoName.Contains("TCardChannel"))lG->AddEntry(hGroup[iprod][iebin][igroup],Form("T-Card cell in %s",groupTC[igroup].Data()),"PL");
          }
          
        } // param
      }//iebin  
      cGroup->cd(ncolE*nrowE);
      lG->Draw();
      
//      fileName = Form("figures/%s/Comparison_Groups_PerParam_%s%s",
//                      histoName.Data(),prod[iprod].Data(),outputFileName2.Data());
      
      if ( prod[iprod].Contains("data") )
        fileName = Form("figures/%s/GroupComparison_PerEbin_%s%s",
                        histoName.Data(),titleData.Data(),outputFileName2.Data());
      else
      {
        fileName = Form("figures/%s/GroupComparison_PerEbin_%s_%d%s",
                        histoName.Data(),titleMC.Data(),iprod,outputFileName2.Data());
      }
      fileName+=fileFormat;
      cGroup->Print(fileName);
    } // GROUPS
    
    if(!plotRatio) continue;
    
    TCanvas * cR = new TCanvas(Form("cRat_prod_%d_%s",iprod,histoName.Data()),
                               Form("%s, %s",prod[iprod].Data(),histoName.Data()),
                               ncolE*2000,nrowE*2000);
    cR->Divide(ncolE,nrowE);
    
    for(Int_t iebin = 0; iebin < nEBins; iebin++)
    {
      cR->cd(iebin+1);
      
      //if(histoName.Contains("ECell"))gPad->SetLogy();
      gPad->SetGridy();
      
      if(!h[iprod][iebin][0]) continue;
      
      for(Int_t iparam = 1; iparam < nparam; iparam++)
      {
        if(!h[iprod][iebin][iparam]) continue;
        
        TH1F* hRat = (TH1F*) h[iprod][iebin][iparam]->Clone(Form("%s_Ratio_PerParam_%s",h[iprod][iebin][iparam]->GetName(),prod[iprod].Data()));
        hRat->Divide(hA[iprod][iebin]);
        
        if(iparam==0) hRat->Draw("H");
        else          hRat->Draw("H same");
        
        hRat->SetYTitle("Ratio to no selection");
        hRat->SetMaximum(2);
        hRat->SetMinimum(0);
        
      } // iparam
      
    } // ebin
    
    cR->cd(ncolE*nrowE);
    l->Draw();
    
//    fileName = Form("figures/%s/Comparison_PerParam_%s%s_Ratios",
//                    histoName.Data(),prod[iprod].Data(),outputFileName2.Data());
    
    if ( prod[iprod].Contains("data") )
      fileName = Form("figures/%s/SMComparison_PerEbin_%s%s_Ratios",
                      histoName.Data(),titleData.Data(),outputFileName2.Data());
    else
    {
      fileName = Form("figures/%s/SMComparison_PerEbin_%s_%d%s_Ratios",
                      histoName.Data(),titleMC.Data(),iprod,outputFileName2.Data());
    }
    fileName+=fileFormat;
    cR->Print(fileName);
  }
}

//-------------------------------------------------------------------------
/// Define the histogram x axis range, normalization and rebin, 
/// depending the cluster parameter
//-------------------------------------------------------------------------
void DefineHistogramSettings(TString histoName, 
                             Float_t & xmin    , Float_t & xmax, 
                             Float_t & xNormMin, Float_t & xNormMax,
                             Bool_t & normtomax, Int_t   & rebin     )
{
  xmin      = -1;
  xmax      = 10000;
  xNormMin  = -1; 
  xNormMax  = 10000;
  rebin     = 1;
  normtomax = kFALSE;
  if ( histoName == "SMM02" )          
  { 
    xmin = 0.1; xmax = 1.10; 
    xNormMin = 0.6; xNormMax = 2;
    rebin     = 4;
    normtomax = 0;
  }
  if ( histoName == "SMM02NoCut" )     
  { 
    xmin = 0.1; xmax = 1.0; 
    xNormMin = 0.1; xNormMax = 0.4;
    rebin     = 2;
    normtomax = 1;
  }
  if ( histoName == "SMNCell" && !histoName.Contains("Module") )        
  { 
    xmin = 2   ; xmax = 12   ; 
    rebin     = 1;
    normtomax = 0;
  }
  if  ( histoName.Contains("SMM20Low") )
  { 
    xmin = 0.05; xmax = 0.25; 
    //xNormMin = 0.6; xNormMax = 2;
    rebin     = 1;
    normtomax = 0;
  }
  if ( histoName.Contains("SMM20Hig") )
  { 
    xmin = 0.05; xmax = 0.30; 
    //xNormMin = 0.6; xNormMax = 2;
    rebin     = 1;
    normtomax = 0;
  }
  if ( histoName.Contains("Module") )
  { 
    xmin = 0.0; xmax = 100000; 
    rebin     = 1;
    normtomax = 0;
    
    if ( histoName.Contains("Rat") )
    {
      rebin     = 5;
    }    
    else if ( histoName.Contains("Tot") )
    {
      xmin = 0.0; xmax = 20; 
      rebin     = 1;
      
    }
    else if ( histoName.Contains("ECellModuleMax") || 
              histoName.Contains("ECellModuleOut") )
    {
      xmin = 0.0; xmax = 10; 
      rebin     = 5;
    }
  }
  
  if ( histoName.Contains("TCardChannelM02") )
  { 
    xmin = 0.1; xmax = 0.5; 
    rebin     = 2;
    normtomax = 0;
    if( histoName.Contains("NoCut") ) normtomax = 1;
  }  
  
  if ( histoName.Contains ( "SMEMaxEClusterRat") )        
  { 
    xmin = 0   ; xmax = 1   ; 
    rebin     = 5;
    normtomax = 0;
    if ( histoName.Contains("LowM02")  ) xmin = 0.4   ;
    if ( histoName.Contains("HighM02") ) xmax = 0.8   ;    
  }
}

//-------------------------------------------------------------------------
/// Prepare the histogram for plotting, set the axis ranges, titles,
/// line and markers style colors, etc.
//-------------------------------------------------------------------------
Bool_t SetHistogramRangesAndTitles
(
 TH1D* histo, TString histoName, Int_t iprod, 
 Int_t rebin, Float_t xmin, Float_t xmax, 
 Bool_t normtomax, Float_t xNormMin, Float_t xNormMax, 
 Float_t lowEBin, Float_t highEBin, Bool_t debug
)
{
  if ( !histo )
    return kFALSE;
 
  Int_t color    [] = {1,4,2,kYellow-2,8,kCyan,kOrange+2,kViolet,kOrange-2,4,2,6,8,9,kYellow-6,kCyan,kOrange+2,kViolet,kOrange-2};
  Int_t lineStyle[] = {1,1,1,1,1,1,        1,    1,        1,      1,        1,2,2,2,2,2,        2,    2,        2,      2,        2,2,2,2,2,2,2,2,2,2,2,2};
  Int_t marker   [] = {20,20,20,21,24,24,24,24,24,24,24};
  
  histo->SetLineColor(color[iprod]);
  
  histo->SetLineWidth(2);
  
  histo->SetLineStyle(lineStyle[iprod]);
  
  histo->SetMarkerStyle(marker[iprod]);
  
  histo->SetMarkerColor(color[iprod]);
  
  histo->SetMarkerSize(0.5);
  
  histo->SetTitleOffset(1.8,"Y");
  
  histo->SetTitle(Form("%2.1f < #it{E} < %2.1f GeV",lowEBin,highEBin));
  
  if(rebin > 1) 
  {
    histo->Rebin(rebin);
  }
  
  histo->SetAxisRange(xmin,xmax,"X");
  
  if ( !normtomax )
  {          
    Int_t   xNormMinBin = histo->FindBin(xNormMin);
    Int_t   xNormMaxBin = histo->FindBin(xNormMax);
    
    Double_t integral = histo->Integral(xNormMinBin, xNormMaxBin);
    
    if ( debug ) 
      printf("\t \t Integral %e bins [%d, %d], range [%2.2f, %2.2f]\n", 
                        integral,xNormMinBin,xNormMaxBin,xNormMin,xNormMax);
    
    if ( integral < 1e-10 ) 
      return kFALSE;
    
    histo->Scale(1./integral);
    
    if(!histoName.Contains("Module"))
    {
      if(histoName.Contains("MEtaPhi"))
        histo->SetYTitle("1/ integral d #it{N} / d #sigma_{#eta#varphi}^{2}");
      else if(histoName.Contains("MEta"))
        histo->SetYTitle("1/ integral d #it{N} / d #sigma_{#eta}^{2}");
      else if(histoName.Contains("MPhi"))
        histo->SetYTitle("1/ integral d #it{N} / d #sigma_{#varphi}^{2}");
      else if(histoName.Contains("SMEMaxEClusterRat"))
        histo->SetYTitle("1/ integral d #it{N} / d (#it{E}_{cell}^{max}/#it{E}_{cluster})");
      else if(histoName.Contains("MNCell"))
        histo->SetYTitle("1/ integral d #it{N} / d #it{n}^{w}_{cell}");
      else if(histoName.Contains("MNLocMax"))
        histo->SetYTitle("1/ integral d #it{N} / d #it{n}_{LM}");            
      else if(histoName.Contains("M20"))
        histo->SetYTitle("1/ integral d #it{N} / d #sigma_{short}^{2}");
      else if(histoName.Contains("M02"))
        histo->SetYTitle("1/ integral d #it{N} / d #sigma_{long}^{2}");
    }
    else histo->SetYTitle("1/ integral d #it{N} ");
    
  } // normalize to integral
  else
  {
    Double_t scale = histo->GetBinContent(histo->FindBin(0.245));
    if ( debug ) printf("\t \t \t Norm to max scale %2.2e\n",scale);
    
    if ( scale < 1e-10 )
      return kFALSE;
    
    histo->Scale(1./scale);
    histo->SetYTitle("Norm. to max at 0.25");
    
  } // normalize to maximum at photon peak
  
  if(!histoName.Contains("Module"))
  {
    if(histoName.Contains("MEtaPhi"))
      histo->SetXTitle("#sigma_{#eta#varphi}^{2}");
    else if(histoName.Contains("MEta"))
      histo->SetXTitle("#sigma_{#eta}^{2}");
    else if(histoName.Contains("MPhi"))
      histo->SetXTitle("#sigma_{#varphi}^{2}");
    else if(histoName.Contains("MNCell"))
      histo->SetXTitle("#it{n}_{cell}^{w}");
    else if(histoName.Contains("MNLocMax"))
      histo->SetXTitle("#it{n}_{LM}");      
    else if(histoName.Contains("M20"))
      histo->SetXTitle("#sigma_{#short}^{2}");
    else if(histoName.Contains("M02"))
      histo->SetXTitle("#sigma_{long}^{2}");    
  }
  else if(histoName.Contains("ECell"))
  {
    histo->SetAxisRange(0,highEBin*1.05,"X");
  }
  return kTRUE;
}
