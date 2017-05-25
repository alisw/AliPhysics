///
/// \file DrawProductionComparison.C
/// \ingroup CaloTrackCorrMacrosQA
/// \brief Plot N productions comparison analysis of QA histograms from EMCal PWG-GA wagon
///
/// Macro to plot comparison of different
/// distributions (spectra, correlations)
/// produced in QA trains but different data or cuts on the same data
/// Based on the plots provided by DrawAnaCaloTrackQA.C. 
///
///  *** To be completed with more cases ***
///
/// To execute: root -q -b -l DrawProductionComparison.C'("Pi0IM_GammaTrackCorr_EMCAL","AnalysisResults.root")'
///
/// The trigger name might change depending on the wagon / data type
/// In simulations only the "default" case is available
/// On data, there can be different triggers, depending on the period
/// * default : min bias like triggers, kINT7, kINT1, kMB
/// * EMCAL_L0: kEMC7 L0 EMCal
/// * DCAL_L0 : kEMC7 L0 DCal
/// * EMCAL_L1: kEMCEGA L1 EG1 EMCal
/// * DCAL_L1 : kEMCEGA L1 EG1 DCal 
/// * EMCAL_L2: kEMCEGA L1 EG2 EMCal
/// * DCAL_L2 : kEMCEGA L1 EG2 DCal
/// A plot will be produced for each of the triggers, if they existed in the data.

///
/// The input files must be placed in different directories,
/// each one defined in the string array "prod"m for ex.:
///     TString prod[] = {"AOD142","AOD115","ESD"};
/// and provide a name to be used in the legeds
///      TString prodLeg[] = {"AOD 142","AOD 115","ESD"};
/// Both have to be modified inside the macro.
///
/// The number of productions has to be specified
///     const Int_t nProd = 3;
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)
///

//---------------------------------------------------------
// Set includes and declare methods for compilation

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TFile.h>
#include <TDirectoryFile.h>
#include <TList.h>
#include <TString.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TH1F.h>
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

#endif

void ProcessTrigger
(TString trigName = "default", 
 TString fileName = "AnalysisResults.root", 
 TString listName = "Pi0IM_GammaTrackCorr_EMCAL");

void CaloQA    (Int_t icalo);
void CorrelQA  (Int_t icalo);

void TrackQA ();
void VertexQA();
void CentralityQA();

TObject * GetHisto      (TString histoName, Int_t iprod);
Bool_t    GetFileAndList(TString fileName, TString listName, TString trigName);
//
//---------------------------------------------------------

//-----------------------
// Some global variables

/// productions to be compared, directory name
TString      prod   [] = {"DCAoffPIDoff","DCAonPIDoff","DCAoffPIDon","DCAonPIDon"}; 
/// productions name used in legends
TString      prodLeg[] = {"DCA off - PID off","DCA on - PID off","DCA off - PID on","DCA on - PID on"}; 
const Int_t nProd = 4;       /// total number of productions to compare

TDirectoryFile *dir[nProd];  /// TDirectory file where lists per trigger are stored in train ouput
TList  *list   [nProd];      /// TList with histograms for a given trigger
TFile  *file   [nProd];      /// Input train file
Float_t nEvents[nProd];      /// number of events container

TString histoTag = "";       /// file names tag, basically the trigger and calorimeter combination 
TString format   = "eps";    /// output plots format: eps, pdf, etc.

/// pre-defined colors list 
Int_t color[]={kBlack,kRed,kBlue,kOrange+1,kYellow+1,kGreen+2,kCyan+1,kViolet,kMagenta+2,kGray};
//-----------------------

//_______________________________________________________________________
/// Main method, produce the plots with the comparisons
///
/// \param listName: Name of list with histograms in file (same for all productions)
/// \param fileName: File name (same for all productions)
/// \param fileFormat: define the type of figures: eps, pdf, etc.
//_______________________________________________________________________
void DrawProductionComparison
(
 TString listName   = "Pi0IM_GammaTrackCorr_EMCAL",
 TString fileName   = "AnalysisResults.root",
 TString fileFormat = "eps"
 )
{
  format  = fileFormat;
  
  printf("Open <%s>; Get List : <%s>; format %s\n",fileName.Data(),listName.Data(),format.Data());
  
  // Process each of the triggers
  //
  ProcessTrigger("default" ,fileName,listName);
  ProcessTrigger("EMCAL_L0",fileName,listName);
  ProcessTrigger("EMCAL_L1",fileName,listName);
  ProcessTrigger("EMCAL_L2",fileName,listName);
  ProcessTrigger("DCAL_L0" ,fileName,listName);
  ProcessTrigger("DCAL_L1" ,fileName,listName);
  ProcessTrigger("DCAL_L2" ,fileName,listName);
}

/// 
/// Produce the plots per trigger, options are:
/// * default : min bias like triggers, kINT7, kINT1, kMB
/// * EMCAL_L0: kEMC7 L0 EMCal
/// * DCAL_L0 : kEMC7 L0 DCal
/// * EMCAL_L1: kEMCEGA L1 EG1 EMCal
/// * DCAL_L1 : kEMCEGA L1 EG1 DCal 
/// * EMCAL_L2: kEMCEGA L1 EG2 EMCal
/// * DCAL_L2 : kEMCEGA L1 EG2 DCal
///
/// Input:
/// \param trigName: trigger case name
/// \param fileName: File name (same for all productions)
/// \param listName: Name of list with histograms in file
//_______________________________________________________________________
void ProcessTrigger(TString trigName, TString fileName, TString listName)
{  
  // Access the file and list of histograms, global variables
  // Check first that the requested trigger exist
  Bool_t ok = GetFileAndList(fileName, listName, trigName);
  if ( !ok ) return;
  
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(000000);
  gStyle->SetPadRightMargin(0.02);
  //gStyle->SetPadTopMargin(0.02);
  //gStyle->SetPadLeftMargin(0.15);
  gStyle->SetTitleFontSize(0.05);
  
  Int_t nCalo = 2;
  Int_t calo  = 0;
  if     (trigName.Contains("EMCAL")) { calo = 0 ; nCalo = 1 ; }
  else if(trigName.Contains("DCAL" )) { calo = 1 ; nCalo = 2 ; }
  
  TString caloString [] = {"EMCAL","DCAL"}; 
  
  histoTag = trigName;
  
  // Plot basic Track QA
  TrackQA();
  
  VertexQA();
  
  CentralityQA();
  
  for(Int_t icalo = calo; icalo < nCalo; icalo++)
  {
    if(trigName.Contains("default")) histoTag=Form("%s_%s",caloString[icalo].Data(),trigName.Data());
    
    // Plot basic QA
    CaloQA(icalo);
    
    // Plot basic correlation QA
    CorrelQA(icalo);
  }
}

///
/// Plot basic calorimeter QA histograms.
/// 2 canvases with 2-4 pads
/// * cluster spectra after matching and PID cuts, and their ratio to different productions
/// * cluster track matching residuals in eta and phi directions, and their ratio to different productions
/// 
/// To be updated
///
/// \param icalo: 0 EMCal, 1 DCal
//______________________________________
void CaloQA(Int_t icalo)
{
  // Declare the different histograms, arrays input is production
  TH1F* hRaw [nProd];
  TH1F* hCorr[nProd];
  TH1F* hTM  [nProd];
  TH1F* hShSh[nProd];
  
  TH1F* hRatRaw [nProd-1];
  TH1F* hRatCorr[nProd-1];
  TH1F* hRatTM  [nProd-1];
  TH1F* hRatShSh[nProd-1];
  
  TH2F* h2TrackMatchResEtaNeg[nProd];
  TH2F* h2TrackMatchResEtaPos[nProd];
  TH2F* h2TrackMatchResPhiNeg[nProd];
  TH2F* h2TrackMatchResPhiPos[nProd];
  
  TH1F* hTrackMatchResEtaNeg[nProd];
  TH1F* hTrackMatchResEtaPos[nProd];
  TH1F* hTrackMatchResPhiNeg[nProd];
  TH1F* hTrackMatchResPhiPos[nProd];
  
  TH1F* hRatTrackMatchResEtaNeg[nProd-1];
  TH1F* hRatTrackMatchResEtaPos[nProd-1];
  TH1F* hRatTrackMatchResPhiNeg[nProd-1];
  TH1F* hRatTrackMatchResPhiPos[nProd-1];
  
  TH2F* h2TrackMatchResEtaNegTrackPt[nProd];
  TH2F* h2TrackMatchResEtaPosTrackPt[nProd];
  TH2F* h2TrackMatchResPhiNegTrackPt[nProd];
  TH2F* h2TrackMatchResPhiPosTrackPt[nProd];
  
  TH1F* hTrackMatchResEtaNegTrackPt[nProd];
  TH1F* hTrackMatchResEtaPosTrackPt[nProd];
  TH1F* hTrackMatchResPhiNegTrackPt[nProd];
  TH1F* hTrackMatchResPhiPosTrackPt[nProd];
  
  TH1F* hRatTrackMatchResEtaNegTrackPt[nProd-1];
  TH1F* hRatTrackMatchResEtaPosTrackPt[nProd-1];
  TH1F* hRatTrackMatchResPhiNegTrackPt[nProd-1];
  TH1F* hRatTrackMatchResPhiPosTrackPt[nProd-1];
  
  //Legend for productions
  TLegend lprod(0.6,0.475,0.95,0.675);
  lprod.SetTextSize(0.04);
  lprod.SetBorderSize(0);
  lprod.SetFillColor(0);
  
  for(Int_t iprod = 0; iprod <  nProd; iprod++)
  {      
    hRaw [iprod] = (TH1F*) GetHisto(Form("AnaPhoton_Calo%d_hPt_Cut_0_Open"    ,icalo),iprod);
    hCorr[iprod] = (TH1F*) GetHisto(Form("AnaPhoton_Calo%d_hPt_Cut_4_NCells"  ,icalo),iprod);
    hTM  [iprod] = (TH1F*) GetHisto(Form("AnaPhoton_Calo%d_hPt_Cut_7_Matching",icalo),iprod);
    hShSh[iprod] = (TH1F*) GetHisto(Form("AnaPhoton_Calo%d_hPt_Cut_9_PID"     ,icalo),iprod);
    
    if(!hRaw[iprod]) return;
    
    hRaw [iprod]->Sumw2();
    hCorr[iprod]->Sumw2();
    hTM  [iprod]->Sumw2();
    hShSh[iprod]->Sumw2();
    
    hRaw [iprod]->Scale(1./nEvents[iprod]);
    hCorr[iprod]->Scale(1./nEvents[iprod]);
    hTM  [iprod]->Scale(1./nEvents[iprod]);
    hShSh[iprod]->Scale(1./nEvents[iprod]);
    
    hRaw[iprod]->SetMarkerColor(color[iprod]);
    hRaw[iprod]->SetMarkerStyle(24);
    
    hCorr[iprod]->SetTitle("Cluster spectra with/out cuts");
    hCorr[iprod]->SetYTitle("1/N_{events} dN/dp_{T}");
    
    hCorr[iprod]->SetTitleOffset(1.5,"Y");
    hCorr[iprod]->SetMarkerColor(color[iprod]);
    hCorr[iprod]->SetMarkerStyle(20);
    hCorr[iprod]->SetAxisRange(0.,30.,"X");
    //hCorr[iprod]->SetMaximum(1.1);
    //hCorr[iprod]->SetMinimum(0);
    
    hTM  [iprod]->SetMarkerColor(color[iprod]);
    hTM  [iprod]->SetMarkerStyle(21);
    
    hShSh[iprod]->SetMarkerColor(color[iprod]);
    hShSh[iprod]->SetMarkerStyle(22);
    
    hRaw [iprod]->SetTitleOffset(1.5,"Y");
    hTM  [iprod]->SetTitleOffset(1.5,"Y");
    hShSh[iprod]->SetTitleOffset(1.5,"Y");
    hCorr[iprod]->SetTitleOffset(1.5,"Y");
    
    if(iprod > 0)
    {
      hRatRaw [iprod-1] = (TH1F*)hRaw [iprod]->Clone(Form("hRatRaw%s_%s" ,prod[iprod].Data(),histoTag.Data()));
      hRatCorr[iprod-1] = (TH1F*)hCorr[iprod]->Clone(Form("hRatCorr%s_%s",prod[iprod].Data(),histoTag.Data()));
      hRatTM  [iprod-1] = (TH1F*)hTM  [iprod]->Clone(Form("hRatTM%s_%s"  ,prod[iprod].Data(),histoTag.Data()));
      hRatShSh[iprod-1] = (TH1F*)hShSh[iprod]->Clone(Form("hRatShSh%s_%s",prod[iprod].Data(),histoTag.Data()));
      
      hRatRaw [iprod-1]->Divide(hRatRaw [iprod-1],hRaw [0],1.000,1,"B");
      hRatCorr[iprod-1]->Divide(hRatCorr[iprod-1],hCorr[0],0.975,1,"B");
      hRatTM  [iprod-1]->Divide(hRatTM  [iprod-1],hTM  [0],0.950,1,"B");
      hRatShSh[iprod-1]->Divide(hRatShSh[iprod-1],hShSh[0],0.925,1,"B");
    }
    
    // Cluster-Track Matching Residuals
    
    // E cluster bin
    h2TrackMatchResEtaNeg[iprod] = (TH2F*) GetHisto(Form("AnaPhoton_Calo%d_hTrackMatchedDEtaNegNoCut",icalo),iprod);
    h2TrackMatchResEtaPos[iprod] = (TH2F*) GetHisto(Form("AnaPhoton_Calo%d_hTrackMatchedDEtaPosNoCut",icalo),iprod);
    h2TrackMatchResPhiNeg[iprod] = (TH2F*) GetHisto(Form("AnaPhoton_Calo%d_hTrackMatchedDPhiNegNoCut",icalo),iprod);
    h2TrackMatchResPhiPos[iprod] = (TH2F*) GetHisto(Form("AnaPhoton_Calo%d_hTrackMatchedDPhiPosNoCut",icalo),iprod);
    
    Float_t emin = 0.5;
    Float_t emax = 2.5;
    if ( histoTag.Contains("L0") )
    {
      emin = 2.5;
      emax = 5;
    }
    else if ( histoTag.Contains("L1") )
    {
      emin = 5;
      emax = 15;
    }
    else if ( histoTag.Contains("L2") )
    {
      emin = 5;
      emax = 15;
    }
    
    Float_t binMin = h2TrackMatchResEtaNeg[iprod]->GetXaxis()->FindBin(emin);
    Float_t binMax = h2TrackMatchResEtaNeg[iprod]->GetXaxis()->FindBin(emax);
    
    hTrackMatchResEtaNeg[iprod] = 
    (TH1F*) h2TrackMatchResEtaNeg[iprod]->ProjectionY(Form("TMProjEtaNeg%s_%s",
                                                           prod[iprod].Data(),histoTag.Data()),binMin, binMax);
    hTrackMatchResEtaPos[iprod] = 
    (TH1F*) h2TrackMatchResEtaPos[iprod]->ProjectionY(Form("TMProjEtaPos%s_%s",
                                                           prod[iprod].Data(),histoTag.Data()),binMin, binMax);
    hTrackMatchResPhiNeg[iprod] = 
    (TH1F*) h2TrackMatchResPhiNeg[iprod]->ProjectionY(Form("TMProjPhiNeg%s_%s",
                                                           prod[iprod].Data(),histoTag.Data()),binMin, binMax);
    hTrackMatchResPhiPos[iprod] = 
    (TH1F*) h2TrackMatchResPhiPos[iprod]->ProjectionY(Form("TMProjPhiPos%s_%s",
                                                           prod[iprod].Data(),histoTag.Data()),binMin, binMax);
    
    hTrackMatchResEtaNeg[iprod]->SetXTitle("#Delta #eta");
    hTrackMatchResEtaNeg[iprod]->SetYTitle("Entries / N events");
    hTrackMatchResEtaNeg[iprod]->SetTitle(Form("Track-cluster #eta residuals, %2.1f < #it{E}^{cluster} < %2.1f GeV",emin,emax));
    hTrackMatchResEtaNeg[iprod]->SetAxisRange(-0.025,0.025,"X");
    hTrackMatchResEtaNeg[iprod]->Sumw2();
    hTrackMatchResEtaNeg[iprod]->SetMarkerStyle(24);
    hTrackMatchResEtaNeg[iprod]->SetMarkerColor(color[iprod]);
    
    hTrackMatchResEtaPos[iprod]->Sumw2();
    hTrackMatchResEtaPos[iprod]->SetYTitle("Entries / N events");
    hTrackMatchResEtaPos[iprod]->SetXTitle("#Delta #eta");
    hTrackMatchResEtaPos[iprod]->SetTitle(Form("Track-cluster #eta residuals, %2.1f < #it{E}^{cluster} < %2.1f GeV",emin,emax));
    hTrackMatchResEtaPos[iprod]->SetAxisRange(-0.025,0.025,"X");
    hTrackMatchResEtaPos[iprod]->SetMarkerStyle(25);
    hTrackMatchResEtaPos[iprod]->SetMarkerColor(color[iprod]);
    
    hTrackMatchResPhiNeg[iprod]->SetXTitle("#Delta #varphi");
    hTrackMatchResPhiNeg[iprod]->SetTitle(Form("Track-cluster #varphi residuals, %2.1f < #it{E}^{cluster} < %2.1f GeV",emin,emax));
    hTrackMatchResPhiNeg[iprod]->SetYTitle("entries / N events");
    hTrackMatchResPhiNeg[iprod]->SetAxisRange(-0.025,0.025,"X");
    hTrackMatchResPhiNeg[iprod]->Sumw2();
    hTrackMatchResPhiNeg[iprod]->SetMarkerStyle(24);
    hTrackMatchResPhiNeg[iprod]->SetMarkerColor(color[iprod]);
    
    hTrackMatchResPhiPos[iprod]->Sumw2();
    hTrackMatchResPhiPos[iprod]->SetYTitle("Entries / N events");
    hTrackMatchResPhiPos[iprod]->SetXTitle("#Delta #varphi");
    hTrackMatchResPhiPos[iprod]->SetTitle(Form("Track-cluster #varphi residuals, %2.1f < #it{E}^{cluster} < %2.1f GeV",emin,emax));
    hTrackMatchResPhiPos[iprod]->SetAxisRange(-0.025,0.025,"X");
    hTrackMatchResPhiPos[iprod]->SetMarkerStyle(25);
    hTrackMatchResPhiPos[iprod]->SetMarkerColor(color[iprod]);
    
    hTrackMatchResEtaNeg[iprod]->Scale(1./nEvents[iprod]);
    hTrackMatchResEtaPos[iprod]->Scale(1./nEvents[iprod]);
    hTrackMatchResPhiNeg[iprod]->Scale(1./nEvents[iprod]);
    hTrackMatchResPhiPos[iprod]->Scale(1./nEvents[iprod]);
    
    hTrackMatchResEtaNeg[iprod]->SetTitleOffset(1.5,"Y");
    hTrackMatchResEtaPos[iprod]->SetTitleOffset(1.5,"Y");
    hTrackMatchResPhiNeg[iprod]->SetTitleOffset(1.5,"Y");
    hTrackMatchResPhiPos[iprod]->SetTitleOffset(1.5,"Y");
    
    if(iprod > 0)
    {
      hRatTrackMatchResPhiPos[iprod-1] = 
      (TH1F*)hTrackMatchResPhiPos[iprod]->Clone(Form("hRatPhiPos%s_%s",prod[iprod].Data(),histoTag.Data()));
      hRatTrackMatchResPhiNeg[iprod-1] = 
      (TH1F*)hTrackMatchResPhiNeg[iprod]->Clone(Form("hRatPhiNeg%s_%s",prod[iprod].Data(),histoTag.Data()));
      hRatTrackMatchResEtaPos[iprod-1] = 
      (TH1F*)hTrackMatchResEtaPos[iprod]->Clone(Form("hRatEtaPos%s_%s",prod[iprod].Data(),histoTag.Data()));
      hRatTrackMatchResEtaNeg[iprod-1] = 
      (TH1F*)hTrackMatchResEtaNeg[iprod]->Clone(Form("hRatEtaNeg%s_%s",prod[iprod].Data(),histoTag.Data()));
      
      hRatTrackMatchResPhiPos[iprod-1]->Divide(hRatTrackMatchResPhiPos[iprod-1],hTrackMatchResPhiPos[0],1.000,1,"B");
      hRatTrackMatchResPhiNeg[iprod-1]->Divide(hRatTrackMatchResPhiNeg[iprod-1],hTrackMatchResPhiNeg[0],1.000,1,"B");
      hRatTrackMatchResEtaPos[iprod-1]->Divide(hRatTrackMatchResEtaPos[iprod-1],hTrackMatchResEtaPos[0],1.000,1,"B");
      hRatTrackMatchResEtaNeg[iprod-1]->Divide(hRatTrackMatchResEtaNeg[iprod-1],hTrackMatchResEtaNeg[0],1.000,1,"B");
    }
    
    // pt track bin
    h2TrackMatchResEtaNegTrackPt[iprod] = (TH2F*) GetHisto(Form("AnaPhoton_Calo%d_hTrackMatchedDEtaNegTrackPtNoCut",icalo),iprod);
    h2TrackMatchResEtaPosTrackPt[iprod] = (TH2F*) GetHisto(Form("AnaPhoton_Calo%d_hTrackMatchedDEtaPosTrackPtNoCut",icalo),iprod);
    h2TrackMatchResPhiNegTrackPt[iprod] = (TH2F*) GetHisto(Form("AnaPhoton_Calo%d_hTrackMatchedDPhiNegTrackPtNoCut",icalo),iprod);
    h2TrackMatchResPhiPosTrackPt[iprod] = (TH2F*) GetHisto(Form("AnaPhoton_Calo%d_hTrackMatchedDPhiPosTrackPtNoCut",icalo),iprod);
  
    if ( !h2TrackMatchResEtaNegTrackPt[iprod] ) continue;
    
    binMin = h2TrackMatchResEtaNegTrackPt[iprod]->GetXaxis()->FindBin(emin);
    binMax = h2TrackMatchResEtaNegTrackPt[iprod]->GetXaxis()->FindBin(emax);
    
    hTrackMatchResEtaNegTrackPt[iprod] = 
    (TH1F*) h2TrackMatchResEtaNegTrackPt[iprod]->ProjectionY(Form("TMProjEtaNegTrackPt%s_%s",
                                                           prod[iprod].Data(),histoTag.Data()),binMin, binMax);
    hTrackMatchResEtaPosTrackPt[iprod] = 
    (TH1F*) h2TrackMatchResEtaPosTrackPt[iprod]->ProjectionY(Form("TMProjEtaPosTrackPt%s_%s",
                                                           prod[iprod].Data(),histoTag.Data()),binMin, binMax);
    hTrackMatchResPhiNegTrackPt[iprod] = 
    (TH1F*) h2TrackMatchResPhiNegTrackPt[iprod]->ProjectionY(Form("TMProjPhiNegTrackPt%s_%s",
                                                           prod[iprod].Data(),histoTag.Data()),binMin, binMax);
    hTrackMatchResPhiPosTrackPt[iprod] = 
    (TH1F*) h2TrackMatchResPhiPosTrackPt[iprod]->ProjectionY(Form("TMProjPhiPosTrackPt%s_%s",
                                                           prod[iprod].Data(),histoTag.Data()),binMin, binMax);
    
    hTrackMatchResEtaNegTrackPt[iprod]->SetXTitle("#Delta #eta");
    hTrackMatchResEtaNegTrackPt[iprod]->SetYTitle("Entries / N events");
    hTrackMatchResEtaNegTrackPt[iprod]->SetTitle(Form("Track-cluster #eta residuals, %2.1f < #it{p}_{T}^{track} < %2.1f GeV",emin,emax));
    hTrackMatchResEtaNegTrackPt[iprod]->SetAxisRange(-0.025,0.025,"X");
    hTrackMatchResEtaNegTrackPt[iprod]->Sumw2();
    hTrackMatchResEtaNegTrackPt[iprod]->SetMarkerStyle(24);
    hTrackMatchResEtaNegTrackPt[iprod]->SetMarkerColor(color[iprod]);
    
    hTrackMatchResEtaPosTrackPt[iprod]->Sumw2();
    hTrackMatchResEtaPosTrackPt[iprod]->SetYTitle("Entries / N events");
    hTrackMatchResEtaPosTrackPt[iprod]->SetXTitle("#Delta #eta");
    hTrackMatchResEtaPosTrackPt[iprod]->SetTitle(Form("Track-cluster #eta residuals, %2.1f < #it{p}_{T}^{track} < %2.1f GeV",emin,emax));
    hTrackMatchResEtaPosTrackPt[iprod]->SetAxisRange(-0.025,0.025,"X");
    hTrackMatchResEtaPosTrackPt[iprod]->SetMarkerStyle(25);
    hTrackMatchResEtaPosTrackPt[iprod]->SetMarkerColor(color[iprod]);
    
    hTrackMatchResPhiNegTrackPt[iprod]->SetXTitle("#Delta #varphi");
    hTrackMatchResPhiNegTrackPt[iprod]->SetTitle(Form("Track-cluster #varphi residuals, %2.1f < #it{p}_{T}^{track} < %2.1f GeV",emin,emax));
    hTrackMatchResPhiNegTrackPt[iprod]->SetYTitle("entries / N events");
    hTrackMatchResPhiNegTrackPt[iprod]->SetAxisRange(-0.025,0.025,"X");
    hTrackMatchResPhiNegTrackPt[iprod]->Sumw2();
    hTrackMatchResPhiNegTrackPt[iprod]->SetMarkerStyle(24);
    hTrackMatchResPhiNegTrackPt[iprod]->SetMarkerColor(color[iprod]);
    
    hTrackMatchResPhiPosTrackPt[iprod]->Sumw2();
    hTrackMatchResPhiPosTrackPt[iprod]->SetYTitle("Entries / N events");
    hTrackMatchResPhiPosTrackPt[iprod]->SetXTitle("#Delta #varphi");
    hTrackMatchResPhiPosTrackPt[iprod]->SetTitle(Form("Track-cluster #varphi residuals, %2.1f < #it{p}_{T}^{track} < %2.1f GeV",emin,emax));
    hTrackMatchResPhiPosTrackPt[iprod]->SetAxisRange(-0.025,0.025,"X");
    hTrackMatchResPhiPosTrackPt[iprod]->SetMarkerStyle(25);
    hTrackMatchResPhiPosTrackPt[iprod]->SetMarkerColor(color[iprod]);
    
    hTrackMatchResEtaNegTrackPt[iprod]->Scale(1./nEvents[iprod]);
    hTrackMatchResEtaPosTrackPt[iprod]->Scale(1./nEvents[iprod]);
    hTrackMatchResPhiNegTrackPt[iprod]->Scale(1./nEvents[iprod]);
    hTrackMatchResPhiPosTrackPt[iprod]->Scale(1./nEvents[iprod]);
    
    hTrackMatchResEtaNegTrackPt[iprod]->SetTitleOffset(1.5,"Y");
    hTrackMatchResEtaPosTrackPt[iprod]->SetTitleOffset(1.5,"Y");
    hTrackMatchResPhiNegTrackPt[iprod]->SetTitleOffset(1.5,"Y");
    hTrackMatchResPhiPosTrackPt[iprod]->SetTitleOffset(1.5,"Y");
    
    if(iprod > 0)
    {
      hRatTrackMatchResPhiPosTrackPt[iprod-1] = 
      (TH1F*)hTrackMatchResPhiPosTrackPt[iprod]->Clone(Form("hRatPhiPos%s_%s",prod[iprod].Data(),histoTag.Data()));
      hRatTrackMatchResPhiNegTrackPt[iprod-1] = 
      (TH1F*)hTrackMatchResPhiNegTrackPt[iprod]->Clone(Form("hRatPhiNeg%s_%s",prod[iprod].Data(),histoTag.Data()));
      hRatTrackMatchResEtaPosTrackPt[iprod-1] = 
      (TH1F*)hTrackMatchResEtaPosTrackPt[iprod]->Clone(Form("hRatEtaPos%s_%s",prod[iprod].Data(),histoTag.Data()));
      hRatTrackMatchResEtaNegTrackPt[iprod-1] = 
      (TH1F*)hTrackMatchResEtaNegTrackPt[iprod]->Clone(Form("hRatEtaNeg%s_%s",prod[iprod].Data(),histoTag.Data()));
      
      hRatTrackMatchResPhiPosTrackPt[iprod-1]->Divide(hRatTrackMatchResPhiPosTrackPt[iprod-1],hTrackMatchResPhiPosTrackPt[0],1.000,1,"B");
      hRatTrackMatchResPhiNegTrackPt[iprod-1]->Divide(hRatTrackMatchResPhiNegTrackPt[iprod-1],hTrackMatchResPhiNegTrackPt[0],1.000,1,"B");
      hRatTrackMatchResEtaPosTrackPt[iprod-1]->Divide(hRatTrackMatchResEtaPosTrackPt[iprod-1],hTrackMatchResEtaPosTrackPt[0],1.000,1,"B");
      hRatTrackMatchResEtaNegTrackPt[iprod-1]->Divide(hRatTrackMatchResEtaNegTrackPt[iprod-1],hTrackMatchResEtaNegTrackPt[0],1.000,1,"B");
    }
    
  } // prod loop

  
  /////////////////
  // Make the plots
  /////////////////
  {
    TCanvas * ccalo = new TCanvas(Form("Cluster_%s",histoTag.Data()),"",1000,500);
    ccalo->Divide(2,1);
    
    ccalo->cd(1);
    gPad->SetLogy();
    //gPad->SetLogx();
    
    hCorr[0]->Draw();
    for(Int_t iprod = 0; iprod <  nProd; iprod++)
    {
      hRaw [iprod]->Draw("same");
      hCorr[iprod]->Draw("same");
      hTM  [iprod]->Draw("same");
      hShSh[iprod]->Draw("same");
      
      lprod.AddEntry(hRaw[iprod],prodLeg[iprod],"P");
    }
    
    lprod.Draw();
    
    TLegend lcl(0.55,0.7,0.95,0.89);
    lcl.SetTextSize(0.04);
    lcl.SetBorderSize(0);
    lcl.SetFillColor(0);
    lcl.AddEntry(hRaw [0],"Raw","P");
    lcl.AddEntry(hCorr[0],"No Exotics + non lin.","P");
    lcl.AddEntry(hTM  [0],  "+ Track matching","P");
    lcl.AddEntry(hShSh[0],"+ #lambda^{2}_{0} < 0.4","P");
    lcl.Draw();
    
    ccalo->cd(2);
    //gPad->SetLogy();
    //gPad->SetLogx();
    
    hRatCorr[0]->SetTitle("Cluster spectra ratio");
    hRatCorr[0]->SetYTitle(Form("Ratio data X / %s",prodLeg[0].Data()));
    hRatCorr[0]->SetMinimum(0.850);
    hRatCorr[0]->SetMaximum(1.025);
    hRatCorr[0]->Draw("");
    
    for(Int_t iprod = 0; iprod <  nProd-1; iprod++)
    {
      hRatRaw [iprod]->Draw("same");
      hRatCorr[iprod]->Draw("same");
      hRatTM  [iprod]->Draw("same");
      hRatShSh[iprod]->Draw("same");
    }
    
    TLine l1(0,1,30,1);
    TLine l2(0,0.975,30,0.975);
    TLine l3(0,0.95,30,0.95);
    TLine l4(0,0.925,30,0.925);
    
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    
    ccalo->Print(Form("%s_ClusterSpectraComp.%s",histoTag.Data(),format.Data()));
  }
  
  // Cluster-Track Matching Residual
  {
    TGaxis::SetMaxDigits(3);

    TLine l0Eta(0,hTrackMatchResEtaNeg[0]->GetMinimum(),0,hTrackMatchResEtaNeg[0]->GetMaximum());
    
    TLegend lres(0.6,0.75,0.84,0.89);
    lres.SetTextSize(0.04);
    //lres.SetBorderSize(0);
    lres.SetFillColor(0);
    lres.AddEntry(hTrackMatchResEtaNeg[0],"Negative","P");
    lres.AddEntry(hTrackMatchResEtaPos[0],"Positive","P");
    lres.Draw();
    
    TCanvas * ccalo2 = new TCanvas(Form("MatchingResiduals_%s",histoTag.Data()),"",1000,1000);
    ccalo2->Divide(2,2);
    
    ccalo2->cd(1);
    //gPad->SetLogy();
    
    hTrackMatchResEtaPos[0]->SetMaximum(hTrackMatchResEtaPos[0]->GetMaximum()*1.3);
    
    hTrackMatchResEtaPos[0]->Draw("");
    for(Int_t iprod = 0; iprod < nProd; iprod++)
    {
      hTrackMatchResEtaNeg[iprod]->Draw("same");
      hTrackMatchResEtaPos[iprod]->Draw("same");
    }
    
    l0Eta.SetLineColor(2);
    l0Eta.SetLineWidth(2);
    l0Eta.Draw("same");
    
    lres.Draw();
    lprod.Draw();
    ccalo2->cd(2);
    
    hTrackMatchResPhiPos[0]->Draw("");
    for(Int_t iprod = 0; iprod < nProd; iprod++)
    {
      hTrackMatchResPhiNeg[iprod]->Draw("same");
      hTrackMatchResPhiPos[iprod]->Draw("same");
    }
    
    TLine l0Phi(0,hTrackMatchResPhiNeg[0]->GetMinimum(),0,hTrackMatchResPhiNeg[0]->GetMaximum());
    l0Phi.SetLineColor(2);
    l0Phi.SetLineWidth(2);
    l0Phi.Draw("same");
    
    ccalo2->cd(3);
    //gPad->SetLogy();
    
    hRatTrackMatchResEtaPos[0]->SetMaximum(1.2);
    hRatTrackMatchResEtaPos[0]->SetMinimum(0.8);
    hRatTrackMatchResEtaPos[0]->Draw("");
    hRatTrackMatchResEtaPos[0]->SetYTitle(Form("Ratio data X / %s",prodLeg[0].Data()));
    for(Int_t iprod = 0; iprod < nProd-1; iprod++)
    {
      hRatTrackMatchResEtaNeg[iprod]->Draw("same");
      hRatTrackMatchResEtaPos[iprod]->Draw("same");
    }
    
    //l0.Draw("same");
    
    ccalo2->cd(4);
    
    hRatTrackMatchResPhiPos[0]->SetMaximum(1.2);
    hRatTrackMatchResPhiPos[0]->SetMinimum(0.8);
    hRatTrackMatchResPhiPos[0]->Draw("");
    hRatTrackMatchResPhiPos[0]->SetYTitle(Form("Ratio data X / %s",prodLeg[0].Data()));
    for(Int_t iprod = 0; iprod < nProd-1; iprod++)
    {
      hRatTrackMatchResPhiNeg[iprod]->Draw("same");
      hRatTrackMatchResPhiPos[iprod]->Draw("same");
    }
    
    ccalo2->Print(Form("%s_MatchingResidualsComp.%s",histoTag.Data(),format.Data()));
  }
  
  // Cluster-Track Matching Residual, pT track bin
  if ( h2TrackMatchResEtaNegTrackPt[0] )
  {
    TGaxis::SetMaxDigits(3);
    
    TLine l0Eta(0,hTrackMatchResEtaNeg[0]->GetMinimum(),0,hTrackMatchResEtaNeg[0]->GetMaximum());
    
    TLegend lres(0.6,0.75,0.84,0.89);
    lres.SetTextSize(0.04);
    //lres.SetBorderSize(0);
    lres.SetFillColor(0);
    lres.AddEntry(hTrackMatchResEtaNeg[0],"Negative","P");
    lres.AddEntry(hTrackMatchResEtaPos[0],"Positive","P");
    lres.Draw();
    
    TCanvas * ccalo3 = new TCanvas(Form("MatchingResiduals_TrackPtBin_%s",histoTag.Data()),"",1000,1000);
    ccalo3->Divide(2,2);
    
    ccalo3->cd(1);
    //gPad->SetLogy();
    
    hTrackMatchResEtaPosTrackPt[0]->SetMaximum(hTrackMatchResEtaPosTrackPt[0]->GetMaximum()*1.3);
    
    hTrackMatchResEtaPosTrackPt[0]->Draw("");
    for(Int_t iprod = 0; iprod < nProd; iprod++)
    {
      hTrackMatchResEtaNegTrackPt[iprod]->Draw("same");
      hTrackMatchResEtaPosTrackPt[iprod]->Draw("same");
    }
    
    l0Eta.SetLineColor(2);
    l0Eta.SetLineWidth(2);
    l0Eta.Draw("same");
    
    lres.Draw();
    lprod.Draw();
    ccalo3->cd(2);
    
    hTrackMatchResPhiPosTrackPt[0]->Draw("");
    for(Int_t iprod = 0; iprod < nProd; iprod++)
    {
      hTrackMatchResPhiNegTrackPt[iprod]->Draw("same");
      hTrackMatchResPhiPosTrackPt[iprod]->Draw("same");
    }
    
    TLine l0Phi(0,hTrackMatchResPhiNegTrackPt[0]->GetMinimum(),0,hTrackMatchResPhiNegTrackPt[0]->GetMaximum());
    l0Phi.SetLineColor(2);
    l0Phi.SetLineWidth(2);
    l0Phi.Draw("same");
    
    ccalo3->cd(3);
    //gPad->SetLogy();
    
    hRatTrackMatchResEtaPosTrackPt[0]->SetMaximum(1.2);
    hRatTrackMatchResEtaPosTrackPt[0]->SetMinimum(0.8);
    hRatTrackMatchResEtaPosTrackPt[0]->Draw("");
    hRatTrackMatchResEtaPosTrackPt[0]->SetYTitle(Form("Ratio data X / %s",prodLeg[0].Data()));
    for(Int_t iprod = 0; iprod < nProd-1; iprod++)
    {
      hRatTrackMatchResEtaNegTrackPt[iprod]->Draw("same");
      hRatTrackMatchResEtaPosTrackPt[iprod]->Draw("same");
    }
    
    //l0.Draw("same");
    
    ccalo3->cd(4);
    
    hRatTrackMatchResPhiPosTrackPt[0]->SetMaximum(1.2);
    hRatTrackMatchResPhiPosTrackPt[0]->SetMinimum(0.8);
    hRatTrackMatchResPhiPosTrackPt[0]->Draw("");
    hRatTrackMatchResPhiPosTrackPt[0]->SetYTitle(Form("Ratio data X / %s",prodLeg[0].Data()));
    for(Int_t iprod = 0; iprod < nProd-1; iprod++)
    {
      hRatTrackMatchResPhiNegTrackPt[iprod]->Draw("same");
      hRatTrackMatchResPhiPosTrackPt[iprod]->Draw("same");
    }
    
    ccalo3->Print(Form("%s_MatchingResidualsTrackPtBinComp.%s",histoTag.Data(),format.Data()));
  }
}

///
/// Hybrid Tracks distributions
/// To be updated
//______________________________________
void TrackQA()
{ 
  TH1F * hTrackPt[nProd] ;
  TH1F * hTrackPtSPD[nProd] ;
  TH1F * hTrackPtNoSPD[nProd] ;
  TH1F * hRatTrackPt[nProd-1] ;
  TH1F * hRatTrackPtSPD[nProd-1] ;
  TH1F * hRatTrackPtNoSPD[nProd-1] ;
  
  TH2F * hTrackEtaPhi[nProd] ;
  TH2F * hTrackEtaPhiSPD[nProd] ;
  TH2F * hTrackEtaPhiNoSPD[nProd] ;
  TH1F * hTrackPhi[nProd] ;
  TH1F * hTrackPhiSPD[nProd] ;
  TH1F * hTrackPhiNoSPD[nProd] ;
  TH1F * hRatTrackPhi[nProd-1] ;
  TH1F * hRatTrackPhiSPD[nProd-1] ;
  TH1F * hRatTrackPhiNoSPD[nProd-1] ;
  
  //Legend for productions
  TLegend lprod(0.3,0.475,0.84,0.675);
  lprod.SetTextSize(0.04);
  lprod.SetBorderSize(0);
  lprod.SetFillColor(0);
  
  for(Int_t iprod = 0; iprod <  nProd; iprod++)
  {  
    hTrackPt         [iprod] = (TH1F*) GetHisto("AnaHadrons_hPt"                  ,iprod);
    if ( !hTrackPt[iprod] ) return;
    
    hTrackPtSPD      [iprod] = (TH1F*) GetHisto("AnaHadrons_hPtSPDRefit"          ,iprod);
    hTrackPtNoSPD    [iprod] = (TH1F*) GetHisto("AnaHadrons_hPtNoSPDRefit"        ,iprod);
    hTrackEtaPhiSPD  [iprod] = (TH2F*) GetHisto("AnaHadrons_hEtaPhiSPDRefitPt02"  ,iprod);
    hTrackEtaPhiNoSPD[iprod] = (TH2F*) GetHisto("AnaHadrons_hEtaPhiNoSPDRefitPt02",iprod);
    hTrackEtaPhi     [iprod] = (TH2F*) GetHisto("AnaHadrons_hEtaPhiPositive"      ,iprod);
    hTrackEtaPhi     [iprod]->Add((TH2F*) GetHisto("AnaHadrons_hEtaPhiNegative"   ,iprod));
    
    hTrackPhiSPD     [iprod] = (TH1F*)hTrackEtaPhiSPD  [iprod]->ProjectionY(Form("hTrackPhiSPD%s"  ,prod[iprod].Data()),0,1000);
    hTrackPhiNoSPD   [iprod] = (TH1F*)hTrackEtaPhiNoSPD[iprod]->ProjectionY(Form("hTrackPhiNoSPD%s",prod[iprod].Data()),0,1000);
    hTrackPhi        [iprod] = (TH1F*)hTrackEtaPhi     [iprod]->ProjectionY(Form("hTrackPhi%s"     ,prod[iprod].Data()),0,1000);
    
    hTrackPt     [iprod]->Sumw2();
    hTrackPtSPD  [iprod]->Sumw2();
    hTrackPtNoSPD[iprod]->Sumw2();
    
    hTrackPt     [iprod]->Scale(1./nEvents[iprod]);
    hTrackPtSPD  [iprod]->Scale(1./nEvents[iprod]);
    hTrackPtNoSPD[iprod]->Scale(1./nEvents[iprod]);
    
    hTrackPhi     [iprod]->Sumw2();
    hTrackPhiSPD  [iprod]->Sumw2();
    hTrackPhiNoSPD[iprod]->Sumw2();
    
    hTrackPhi     [iprod]->Scale(1./nEvents[iprod]);
    hTrackPhiSPD  [iprod]->Scale(1./nEvents[iprod]);
    hTrackPhiNoSPD[iprod]->Scale(1./nEvents[iprod]);
    
    hTrackPt[iprod]->SetTitle("Track spectra with/out SPD");
    hTrackPt[iprod]->SetYTitle("1/N_{events} dN/dp_{T}");
    hTrackPt[iprod]->SetTitleOffset(1.5,"Y");
    hTrackPt[iprod]->SetMarkerColor(color[iprod]);
    hTrackPt[iprod]->SetMarkerStyle(20);
    hTrackPt[iprod]->SetAxisRange(0.,30.,"X");
    //hTrackPt[iprod]->SetMaximum(1.1);
    //hTrackPt[iprod]->SetMinimum(0);
    
    hTrackPtSPD[iprod]->SetMarkerColor(color[iprod]);
    hTrackPtSPD[iprod]->SetMarkerStyle(26);
    
    hTrackPtNoSPD[iprod]->SetMarkerColor(color[iprod]);
    hTrackPtNoSPD[iprod]->SetMarkerStyle(25);
    
    hTrackPhi[iprod]->SetTitle("Track #varphi with/out SPD");
    hTrackPhi[iprod]->SetYTitle("1/N_{events} dN/d#varphi");
    hTrackPhi[iprod]->SetTitleOffset(1.5,"Y");
    hTrackPhi[iprod]->SetMarkerColor(color[iprod]);
    hTrackPhi[iprod]->SetMarkerStyle(20);
    hTrackPhi[iprod]->SetAxisRange(0.,30.,"X");
    //hTrackPhi[iprod]->SetMaximum(1.1);
    //hTrackPhi[iprod]->SetMinimum(0);
    
    hTrackPhiSPD[iprod]->SetMarkerColor(color[iprod]);
    hTrackPhiSPD[iprod]->SetMarkerStyle(26);
    
    hTrackPhiNoSPD[iprod]->SetMarkerColor(color[iprod]);
    hTrackPhiNoSPD[iprod]->SetMarkerStyle(25);
    
    if(iprod > 0)
    {
      hRatTrackPhi     [iprod-1] = (TH1F*)hTrackPhi     [iprod]->Clone(Form("hRatTrackPhi%s"     ,prod[iprod].Data()));
      hRatTrackPhiNoSPD[iprod-1] = (TH1F*)hTrackPhiNoSPD[iprod]->Clone(Form("hRatTrackPhiNoSPD%s",prod[iprod].Data()));
      hRatTrackPhiSPD  [iprod-1] = (TH1F*)hTrackPhiSPD  [iprod]->Clone(Form("hRatTrackPhiSPD%s"  ,prod[iprod].Data()));
      
      hRatTrackPhi     [iprod-1]->Divide(hRatTrackPhi     [iprod-1],hTrackPhi     [0],1.000,1,"B");
      hRatTrackPhiSPD  [iprod-1]->Divide(hRatTrackPhiSPD  [iprod-1],hTrackPhiSPD  [0],1.000,1,"B");
      hRatTrackPhiNoSPD[iprod-1]->Divide(hRatTrackPhiNoSPD[iprod-1],hTrackPhiNoSPD[0],1.000,1,"B");
      
      hRatTrackPt     [iprod-1] = (TH1F*)hTrackPt     [iprod]->Clone(Form("hRatTrackPt%s"     ,prod[iprod].Data()));
      hRatTrackPtNoSPD[iprod-1] = (TH1F*)hTrackPtNoSPD[iprod]->Clone(Form("hRatTrackPtNoSPD%s",prod[iprod].Data()));
      hRatTrackPtSPD  [iprod-1] = (TH1F*)hTrackPtSPD  [iprod]->Clone(Form("hRatTrackPtSPD%s"  ,prod[iprod].Data()));
      
      hRatTrackPt     [iprod-1]->Divide(hRatTrackPt     [iprod-1],hTrackPt     [0],1.000,1,"B");
      hRatTrackPtSPD  [iprod-1]->Divide(hRatTrackPtSPD  [iprod-1],hTrackPtSPD  [0],1.000,1,"B");
      hRatTrackPtNoSPD[iprod-1]->Divide(hRatTrackPtNoSPD[iprod-1],hTrackPtNoSPD[0],1.000,1,"B");
    }
  }
  
  
  /////////////////
  // Make the plots
  /////////////////
  
  // Hybrid tracks
  {
    TLegend ltrack(0.6,0.75,0.84,0.89);
    ltrack.SetTextSize(0.04);
    //ltrack.SetBorderSize(0);
    ltrack.SetFillColor(0);
    ltrack.AddEntry(hTrackPt     [0],"All","P");
    ltrack.AddEntry(hTrackPtSPD  [0],"SPD","P");
    ltrack.AddEntry(hTrackPtNoSPD[0],"No SPD","P");
    
    TCanvas * ctrack = new TCanvas(Form("TrackHisto_%s",histoTag.Data()),"",1500,1500);
    ctrack->Divide(2,2);
    ctrack->cd(1);
    gPad->SetLogy();
    hTrackPt[0]->Draw("");
    for(Int_t iprod = 0; iprod < nProd; iprod++)
    {
      hTrackPt     [iprod]->Draw("same");
      hTrackPtSPD  [iprod]->Draw("same");
      hTrackPtNoSPD[iprod]->Draw("same");
    }
    
    ltrack.Draw();
    lprod.Draw();
    
    ctrack->cd(2);
    
    hRatTrackPt[0]->SetMaximum(1.05);
    hRatTrackPt[0]->SetMinimum(0.95);
    hRatTrackPt[0]->Draw("");
    hRatTrackPt[0]->SetYTitle(Form("Ratio data X / %s",prodLeg[0].Data()));
    for(Int_t iprod = 0; iprod < nProd-1; iprod++)
    {
      hRatTrackPt     [iprod]->Draw("same");
      hRatTrackPtSPD  [iprod]->Draw("same");
      hRatTrackPtNoSPD[iprod]->Draw("same");
    }
    
    ctrack->cd(3);
    hTrackPhi[0]->SetMaximum(3.);
    hTrackPhi[0]->SetMinimum(0.);
    hTrackPhi[0]->Draw("");
    for(Int_t iprod = 0; iprod < nProd; iprod++)
    {
      hTrackPhi     [iprod]->Draw("same");
      hTrackPhiSPD  [iprod]->Draw("same");
      hTrackPhiNoSPD[iprod]->Draw("same");
    }
    
    ctrack->cd(4);
    //gPad->SetLogy();
    
    hRatTrackPhi[0]->SetMaximum(1.05);
    hRatTrackPhi[0]->SetMinimum(0.95);
    hRatTrackPhi[0]->Draw("");
    hRatTrackPhi[0]->SetYTitle(Form("Ratio data X / %s",prodLeg[0].Data()));
    for(Int_t iprod = 0; iprod < nProd-1; iprod++)
    {
      hRatTrackPhi     [iprod]->Draw("same");
      hRatTrackPhiSPD  [iprod]->Draw("same");
      hRatTrackPhiNoSPD[iprod]->Draw("same");
    }
    
    ctrack->Print(Form("%s_TrackComp.%s",histoTag.Data(),format.Data()));
  }
}

///
/// cluster-track correlation
/// 
/// To be updated
///
/// \param icalo: 0 EMCal, 1 DCal
//______________________________________
void CorrelQA(Int_t icalo)
{
  TH2F* h2XE[nProd];
  TH2F* h2XEUE[nProd];
  TH1F* hXE[nProd];
  TH1F* hXEUE[nProd];
  TH1F* hRatXE[nProd-1];
  TH1F* hRatXEUE[nProd-1];
  
  //Legend for productions
  TLegend lprod(0.3,0.475,0.84,0.675);
  lprod.SetTextSize(0.04);
  lprod.SetBorderSize(0);
  lprod.SetFillColor(0);
  
  for(Int_t iprod = 0; iprod <  nProd; iprod++)
  {
    h2XE   [iprod]= (TH2F*) GetHisto(Form("AnaPhotonHadronCorr_Calo%d_hXECharged"  ,icalo),iprod);
    h2XEUE [iprod]= (TH2F*) GetHisto(Form("AnaPhotonHadronCorr_Calo%d_hXEUeCharged",icalo),iprod);
    
    if ( !h2XE[iprod] ) return;
    
    Float_t minClusterE = 8;
    TH1F * hTrigger = (TH1F*) GetHisto("AnaPhotonHadronCorr_hPtTrigger",iprod);
    Int_t minClusterEBin = hTrigger->FindBin(minClusterE);
    Float_t nTrig = hTrigger->Integral(minClusterE,100000);
    
    hXE  [iprod] = (TH1F*)h2XE  [iprod]->ProjectionY(Form("hXE%s_%s"  ,prod[iprod].Data(),histoTag.Data()),minClusterEBin,1000);
    hXEUE[iprod] = (TH1F*)h2XEUE[iprod]->ProjectionY(Form("hXEUE%s_%s",prod[iprod].Data(),histoTag.Data()),minClusterEBin,1000);
    
    hXE  [iprod]->Sumw2();
    hXEUE[iprod]->Sumw2();
    
    hXE  [iprod]->Scale(1./nTrig);
    hXEUE[iprod]->Scale(1./nTrig);
    
    hXE[iprod]->SetTitle(Form("#gamma-hadron x_{E}, p_{T,Trig}>%2.1f GeV/c",minClusterE));
    hXE[iprod]->SetYTitle("1/N_{trigger} dN/dx_{E}");
    hXE[iprod]->SetTitleOffset(1.5,"Y");
    hXE[iprod]->SetMarkerColor(color[iprod]);
    hXE[iprod]->SetMarkerStyle(20);
    hXE[iprod]->SetAxisRange(0.,1.,"X");
    //hXE[iprod]->SetMaximum(1.1);
    //hXE[iprod]->SetMinimum(0);
    
    hXEUE[iprod]->SetMarkerColor(color[iprod]);
    hXEUE[iprod]->SetMarkerStyle(25);
    
    if(iprod > 0)
    {
      hRatXE  [iprod-1] = (TH1F*)hXE  [iprod]->Clone(Form("hRatXE%s_%s"  ,prod[iprod].Data(),histoTag.Data()));
      hRatXEUE[iprod-1] = (TH1F*)hXEUE[iprod]->Clone(Form("hRatXEUE%s_%s",prod[iprod].Data(),histoTag.Data()));
      
      hRatXE  [iprod-1]->Divide(hRatXE  [iprod-1],hXE  [0],1.000,1,"B");
      hRatXEUE[iprod-1]->Divide(hRatXEUE[iprod-1],hXEUE[0],1.000,1,"B");
    }
  }
  
  /////////////////
  // Make the plots
  /////////////////
  
  // XE
  {
    TLegend lxe(0.6,0.75,0.84,0.89);
    lxe.SetTextSize(0.04);
    //lxe.SetBorderSize(0);
    lxe.SetFillColor(0);
    lxe.AddEntry(hXE  [0],"Signal+bkg","P");
    lxe.AddEntry(hXEUE[0],"Und. Event","P");
    
    TCanvas * cxe = new TCanvas(Form("XEHisto_%s_%s",histoTag.Data(),histoTag.Data()),"",1000,500);
    cxe->Divide(2,1);
    cxe->cd(1);
    gPad->SetLogy();
    hXE[0]->Draw("");
    for(Int_t iprod = 0; iprod < nProd; iprod++)
    {
      hXE  [iprod]->Draw("same");
      hXEUE[iprod]->Draw("same");
    }
    
    lxe.Draw();
    lprod.Draw();
    
    cxe->cd(2);
    
    hRatXE[0]->SetMaximum(1.05);
    hRatXE[0]->SetMinimum(0.95);
    hRatXE[0]->Draw("");
    hRatXE[0]->SetYTitle(Form("Ratio data X / %s",prodLeg[0].Data()));
    for(Int_t iprod = 0; iprod < nProd-1; iprod++)
    {
      hRatXE  [iprod]->Draw("same");
      hRatXEUE[iprod]->Draw("same");
    }
    
    cxe->Print(Form("%s_XEComp.%s",histoTag.Data(),format.Data()));
  }
  
}

/// Centrality
/// To be updated
//______________________________________
void CentralityQA()    
{
  TH1F* hCen   [nProd];
  TH1F* hRatCen[nProd-1];
  
  //Legend for productions
  TLegend lprod(0.3,0.475,0.84,0.675);
  lprod.SetTextSize(0.04);
  lprod.SetBorderSize(0);
  lprod.SetFillColor(0);
  
  for(Int_t iprod = 0; iprod <  nProd; iprod++)
  {
    hCen[iprod] = (TH1F*) GetHisto("hCentrality",iprod);
    
    if ( !hCen[iprod] ) return;
    
    
    hCen[iprod]->Sumw2();
    
    hCen[iprod]->Scale(1./nEvents[iprod]);
    
    hCen[iprod]->SetTitle("Centrality");
    hCen[iprod]->SetYTitle("1/N_{events} dN/d centrality");
    hCen[iprod]->SetTitleOffset(1.5,"Y");
    hCen[iprod]->SetMarkerColor(color[iprod]);
    hCen[iprod]->SetMarkerStyle(20);
    //hCen[iprod]->SetAxisRange(0.,30.,"X");
    //hCen[iprod]->SetMaximum(1.1);
    //hCen[iprod]->SetMinimum(0);
    
    if(iprod > 0)
    {
      hRatCen[iprod-1] = (TH1F*)hCen[iprod]->Clone(Form("hRatCen%s" ,prod[iprod].Data()));
      
      hRatCen[iprod-1]->Divide(hRatCen[iprod-1],hCen [0],1.000,1,"B");
    }
  }
  
  /////////////////
  // Make the plots
  /////////////////
  // Centrality
  {
    TCanvas * ccen = new TCanvas(Form("Centrality_%s",histoTag.Data()),"",1000,500);
    ccen->Divide(2,1);
    
    ccen->cd(1);
    //gPad->SetLogy();
    
    hCen[0]->Draw();
    for(Int_t iprod = 0; iprod <  nProd; iprod++)
    {
      hCen[iprod]->Draw("same");
    }
    
    lprod.Draw();
    
    ccen->cd(2);
    //gPad->SetLogy();
    
    hRatCen[0]->SetTitle("Centrality");
    hRatCen[0]->SetYTitle(Form("Ratio data X / %s",prodLeg[0].Data()));
    hRatCen[0]->SetMinimum(0.95);
    hRatCen[0]->SetMaximum(1.05);
    hRatCen[0]->Draw("");
    
    for(Int_t iprod = 0; iprod <  nProd-1; iprod++)
    {
      hRatCen [iprod]->Draw("same");
    }
    
    TLine l1(0,1,100,1);
    
    l1.Draw("same");
    
    ccen->Print(Form("%s_CentralityComp.%s",histoTag.Data(),format.Data()));
  }
  
}

///
/// x y z vertex distribution and ratios to different productions.
///
//______________________________________
void VertexQA()
{
  TH1F* hVertex   [3][nProd];
  TH1F* hRatVertex[3][nProd-1];
  
  //Legend for productions
  TLegend lprod(0.6,0.6,0.97,0.85);
  lprod.SetTextSize(0.04);
  lprod.SetBorderSize(0);
  lprod.SetFillColor(0);
  
  for(Int_t iprod = 0; iprod <  nProd; iprod++)
  {
    hVertex[0][iprod] = (TH1F*) GetHisto("hZVertex",iprod);
    hVertex[1][iprod] = (TH1F*) GetHisto("hYVertex",iprod);
    hVertex[2][iprod] = (TH1F*) GetHisto("hXVertex",iprod);
    
    if ( !hVertex[0][iprod] ) return;
    
    hVertex[0][iprod]->Sumw2();
    hVertex[1][iprod]->Sumw2();
    hVertex[2][iprod]->Sumw2();
    
    lprod.AddEntry(hVertex[0][iprod],prodLeg[iprod],"P");
    
    for(Int_t ivertex = 0; ivertex < 3; ivertex++)
    {
      //hVertex[ivertex][iprod]->Sumw2();
      
      hVertex[ivertex][iprod]->Scale(1./nEvents[iprod]);
      
      //hVertex[ivertex][iprod]->SetTitle("Centrality");
      hVertex[ivertex][iprod]->SetYTitle("1/N_{events} dN/ d vertex");
      hVertex[ivertex][iprod]->SetTitleOffset(1.5,"Y");
      hVertex[ivertex][iprod]->SetMarkerColor(color[iprod]);
      hVertex[ivertex][iprod]->SetLineColor(color[iprod]);
      hVertex[ivertex][iprod]->SetMarkerStyle(20);
      if(ivertex==0)hVertex[ivertex][iprod]->SetAxisRange(-10,10.,"X");
      else          hVertex[ivertex][iprod]->SetAxisRange(-1.5,1.5,"X");
      //hVertex[ivertex][iprod]->SetMaximum(1.1);
      //hVertex[ivertex][iprod]->SetMinimum(0);
      
      if(iprod > 0)
      {
        hRatVertex[ivertex][iprod-1] = (TH1F*)hVertex[ivertex][iprod]->Clone(Form("hRatVertex%s_%d" ,prod[iprod].Data(),ivertex));
        
        hRatVertex[ivertex][iprod-1]->Divide(hRatVertex[ivertex][iprod-1],hVertex[ivertex][0],1.000,1,"B");
      }
    }
  }
  
  /////////////////
  // Make the plots
  /////////////////
  // Vertex
  {
    TCanvas * cvertex = new TCanvas(Form("Vertex_%s",histoTag.Data()),"",3*500,2*500);
    cvertex->Divide(3,2);
    Int_t npannel = 1;
    for(Int_t ivertex = 0; ivertex < 3; ivertex++)
    {
      cvertex->cd(npannel);
      //gPad->SetLogy();
      
      //hVertex[ivertex][0]->Draw();
      for(Int_t iprod = 0; iprod <  nProd; iprod++)
      {
        hVertex[ivertex][iprod]->Draw("same");
      }
      
      //lprod.Draw();
      
      cvertex->cd(npannel+3);
      gPad->SetGridy();
      
      //hRatVertex[ivertex][0]->SetTitle("");
      hRatVertex[ivertex][0]->SetYTitle(Form("Ratio data X / %s",prodLeg[0].Data()));
      hRatVertex[ivertex][0]->SetMinimum(0.90);
      hRatVertex[ivertex][0]->SetMaximum(1.10);
      hRatVertex[ivertex][0]->Draw("");
      
      for(Int_t iprod = 0; iprod <  nProd-1; iprod++)
      {
        hRatVertex[ivertex][iprod]->Draw("same");
      }
      lprod.Draw();
      
      npannel++;
    }
    cvertex->Print(Form("%s_VertexComp.%s",histoTag.Data(),format.Data()));
  }
}

///
/// Access the file and list with histograms and number
/// of analyzed events per each production.
///
/// \param fileName: File name (same for all productions)
/// \param listName: Name of list with histograms in file
/// \param trigName: trigger case name
//_____________________________________________________
Bool_t GetFileAndList(TString fileName, TString listName, TString trigName)
{
  Bool_t ok = kTRUE;
  for(Int_t iprod = 0; iprod < nProd; iprod++)
  {
    // First, init to 0
    file   [iprod] = 0;
    dir    [iprod] = 0;
    list   [iprod] = 0;
    nEvents[iprod] = 0;
    
    // now get them
    file[iprod]  = new TFile(Form("%s/%s",prod[iprod].Data(),fileName.Data()),"read");
    if(file[iprod]->Get("hNEvents"))
    {
      nEvents[iprod] = ((TH1F*)file[iprod]->Get("hNEvents"))->GetEntries();
      printf("%s: nEvents %e\n",prod[iprod].Data(),nEvents[iprod]);
    }
    
    dir[iprod] = (TDirectoryFile*) file[iprod]->Get(listName);
    if(dir[iprod])
    {
      list[iprod] = (TList*) dir[iprod]->Get(trigName);
      
      if ( !list[iprod]                    ) continue;
      if (  list[iprod]->GetEntries() <= 0 ) continue;
      
      nEvents[iprod] = ((TH1F*)list[iprod]->FindObject("hNEvents"))->GetEntries();
      printf("%s: nEvents %e\n",prod[iprod].Data(),nEvents[iprod]);
    }
    else
    {
      list[iprod] = (TList*) file[iprod]->Get(trigName);
      
      if ( !list[iprod]                    ) continue;
      if (  list[iprod]->GetEntries() <= 0 ) continue;
      
      nEvents[iprod] = ((TH1F*)list[iprod]->FindObject("hNEvents"))->GetEntries();
      printf("%s: nEvents %e\n",prod[iprod].Data(),nEvents[iprod]);
    }
    
    if(!list[iprod])
    {
      printf("list %s not found\n",trigName.Data());
      ok = kFALSE;
    }
  }
  
  return ok;
}

///
/// Check if the list is available,
/// if not get the histo directly from file
///
/// \return the histogram with the provided name
///
/// \param histoName: histogram name
/// \param iprod: production index
//___________________________________
TObject * GetHisto(TString histoName, Int_t iprod)
{
  if(list[iprod]) return list[iprod]->FindObject(histoName);
  else            return file[iprod]->Get       (histoName);
}
