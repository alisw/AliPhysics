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
/// It is also possible to run comparisons done for pT hard binned productions. 
/// Since the scaling is usually already done, change the bool scaled to "true" 
/// to avoir overnormalization.
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

void Cluster     (Int_t icalo);
void ClusterCells(Int_t icalo);
void Cell        ();
void Correl      (Int_t icalo);
void Isol        (Int_t icalo);
void ShowerShape (Int_t icalo);
void InvMass     (Int_t icalo, TString particle, TString fileName);
void MCParticle  (Int_t icalo, TString particle);
void Track       ();
void Vertex      ();
void Centrality  ();
void PtHard      ();

TObject * GetHisto      (TString histoName, Int_t iprod);
Bool_t    GetFileAndList(TString fileName, TString listName, TString trigName);
//
//---------------------------------------------------------

//-----------------------
// Some global variables

/// productions to be compared, directory name
//TString      prod   [] = {"DCAoffPIDoff","DCAonPIDoff","DCAoffPIDon","DCAonPIDon"}; 
/// productions name used in legends
//TString      prodLeg[] = {"DCA off - PID off","DCA on - PID off","DCA off - PID on","DCA on - PID on"};

TString      prod   [] = {"LHC17l3b_fast_TM","LHC17l4b_fast_TM"}; 
TString      prodLeg[] = {"Geant3, b FAST","Geant4, b FAST"};

const Int_t nProd = 2;       /// total number of productions to compare

TDirectoryFile *dir[nProd];  /// TDirectory file where lists per trigger are stored in train ouput
TList  *histArr[nProd];      /// TList with histograms for a given trigger
TFile  *file   [nProd];      /// Input train file
Float_t nEvents[nProd];      /// number of events container

TString histoTag = "";       /// file names tag, basically the trigger and calorimeter combination 
TString format   = "eps";    /// output plots format: eps, pdf, etc.
TString errType  = "B";//""; /// Error type in histogram division
Bool_t scaled    = kFALSE;   /// In case of pT hard binned productions where the scaling and cross section is already done, do not normalize by N events or get the Sumw2 of histogram

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
 TString fileFormat = "pdf"
 )
{
  format  = fileFormat;
  
  printf("Open <%s>; Get List : <%s>; format %s\n",fileName.Data(),listName.Data(),format.Data());
  
  // Process each of the triggers
  //
  ProcessTrigger("default" ,fileName,listName); // Data min bias, or only one for MC
  
  // EMC/DMC triggered data (not for MC)
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
  Track();

  // Plot basic Global event QAs
  Vertex();
  
  Centrality();
  
  Cell();

  //PtHard(); // Uncomment for pT hard binned MC productions
  
  for(Int_t icalo = calo; icalo < nCalo; icalo++)
  {
    if(trigName.Contains("default")) histoTag=Form("%s_%s",caloString[icalo].Data(),trigName.Data());
    
    // Plot basic calorimeter,  cluster, associated cells,
    // track-matching, shower shape info
    Cluster(icalo);
    ClusterCells(icalo);
    ShowerShape(icalo);

    // Plot clusters Origin, only MC
//    MCParticle(icalo,"Photon");
//    MCParticle(icalo,"PhotonPi0Decay");
//    MCParticle(icalo,"Pi0");
//    MCParticle(icalo,"Electron");

    // Run before InvMassFit.C
    //InvMass(icalo,"Pi0",fileName);
    //InvMass(icalo,"Eta",fileName);
    
    // Plot basic isolation energy 
    Isol(icalo);
    
    // Plot basic correlation 
    Correl(icalo);
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
void Cluster(Int_t icalo)
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
  TLegend lprod(0.65,0.475,0.95,0.675);
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
    
    hRaw [iprod]->Rebin(5);
    hCorr[iprod]->Rebin(5);    
    hTM  [iprod]->Rebin(5);
    hShSh[iprod]->Rebin(5);
    
    if( !scaled )
    {
      hRaw [iprod]->Sumw2();
      hCorr[iprod]->Sumw2();
      hTM  [iprod]->Sumw2();
      hShSh[iprod]->Sumw2();
      
      hRaw [iprod]->Scale(1./nEvents[iprod]);
      hCorr[iprod]->Scale(1./nEvents[iprod]);
      hTM  [iprod]->Scale(1./nEvents[iprod]);
      hShSh[iprod]->Scale(1./nEvents[iprod]);
    }
    
    hRaw[iprod]->SetMarkerColor(color[iprod]);
    hRaw[iprod]->SetMarkerStyle(24);
    
    hCorr[iprod]->SetTitle("Cluster spectra with/out cuts");
    hCorr[iprod]->SetYTitle("1/N_{events} dN/dp_{T}");
    
    hCorr[iprod]->SetTitleOffset(1.5,"Y");
    hCorr[iprod]->SetMarkerColor(color[iprod]);
    hCorr[iprod]->SetMarkerStyle(20);
    hCorr[iprod]->SetAxisRange(10,70.,"X");
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
      
      //      hRatRaw [iprod-1]->Divide(hRatRaw [iprod-1],hRaw [0],1.000,1,errType);
      //      hRatCorr[iprod-1]->Divide(hRatCorr[iprod-1],hCorr[0],0.975,1,errType);
      //      hRatTM  [iprod-1]->Divide(hRatTM  [iprod-1],hTM  [0],0.950,1,errType);
      //      hRatShSh[iprod-1]->Divide(hRatShSh[iprod-1],hShSh[0],0.925,1,errType);
      
      hRatRaw [iprod-1]->Divide(hRatRaw [iprod-1],hRaw [0],1.000,1,errType);
      hRatCorr[iprod-1]->Divide(hRatCorr[iprod-1],hCorr[0],1.000,1,errType);
      hRatTM  [iprod-1]->Divide(hRatTM  [iprod-1],hTM  [0],1.000,1,errType);
      hRatShSh[iprod-1]->Divide(hRatShSh[iprod-1],hShSh[0],1.000,1,errType);
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
    
    if( !scaled )
    {
      hTrackMatchResEtaNeg[iprod]->Sumw2();
      hTrackMatchResEtaPos[iprod]->Sumw2();
      hTrackMatchResPhiNeg[iprod]->Sumw2();
      hTrackMatchResPhiPos[iprod]->Sumw2();
      
      hTrackMatchResEtaNeg[iprod]->Scale(1./nEvents[iprod]);
      hTrackMatchResEtaPos[iprod]->Scale(1./nEvents[iprod]);
      hTrackMatchResPhiNeg[iprod]->Scale(1./nEvents[iprod]);
      hTrackMatchResPhiPos[iprod]->Scale(1./nEvents[iprod]);
    }
    
    hTrackMatchResEtaNeg[iprod]->SetXTitle("#Delta #eta");
    hTrackMatchResEtaNeg[iprod]->SetYTitle("Entries / N events");
    hTrackMatchResEtaNeg[iprod]->SetTitle(Form("Track-cluster #eta residuals, %2.1f < #it{E}^{cluster} < %2.1f GeV",emin,emax));
    hTrackMatchResEtaNeg[iprod]->SetAxisRange(-0.025,0.025,"X");
    hTrackMatchResEtaNeg[iprod]->SetMarkerStyle(24);
    hTrackMatchResEtaNeg[iprod]->SetMarkerColor(color[iprod]);
    
    hTrackMatchResEtaPos[iprod]->SetYTitle("Entries / N events");
    hTrackMatchResEtaPos[iprod]->SetXTitle("#Delta #eta");
    hTrackMatchResEtaPos[iprod]->SetTitle(Form("Track-cluster #eta residuals, %2.1f < #it{E}^{cluster} < %2.1f GeV",emin,emax));
    hTrackMatchResEtaPos[iprod]->SetAxisRange(-0.025,0.025,"X");
    hTrackMatchResEtaPos[iprod]->SetMarkerStyle(25);
    hTrackMatchResEtaPos[iprod]->SetMarkerColor(color[iprod]);
    
    hTrackMatchResPhiNeg[iprod]->SetXTitle("#Delta #varphi");
    hTrackMatchResPhiNeg[iprod]->SetTitle(Form("Track-cluster #varphi residuals, %2.1f < #it{E}^{cluster} < %2.1f GeV",emin,emax));
    hTrackMatchResPhiNeg[iprod]->SetYTitle("entries / N events");
    hTrackMatchResPhiNeg[iprod]->SetAxisRange(-0.035,0.035,"X");
    hTrackMatchResPhiNeg[iprod]->SetMarkerStyle(24);
    hTrackMatchResPhiNeg[iprod]->SetMarkerColor(color[iprod]);
    
    hTrackMatchResPhiPos[iprod]->SetYTitle("Entries / N events");
    hTrackMatchResPhiPos[iprod]->SetXTitle("#Delta #varphi");
    hTrackMatchResPhiPos[iprod]->SetTitle(Form("Track-cluster #varphi residuals, %2.1f < #it{E}^{cluster} < %2.1f GeV",emin,emax));
    hTrackMatchResPhiPos[iprod]->SetAxisRange(-0.035,0.035,"X");
    hTrackMatchResPhiPos[iprod]->SetMarkerStyle(25);
    hTrackMatchResPhiPos[iprod]->SetMarkerColor(color[iprod]);
    
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
      
      hRatTrackMatchResPhiPos[iprod-1]->Divide(hRatTrackMatchResPhiPos[iprod-1],hTrackMatchResPhiPos[0],1.000,1,errType);
      hRatTrackMatchResPhiNeg[iprod-1]->Divide(hRatTrackMatchResPhiNeg[iprod-1],hTrackMatchResPhiNeg[0],1.000,1,errType);
      hRatTrackMatchResEtaPos[iprod-1]->Divide(hRatTrackMatchResEtaPos[iprod-1],hTrackMatchResEtaPos[0],1.000,1,errType);
      hRatTrackMatchResEtaNeg[iprod-1]->Divide(hRatTrackMatchResEtaNeg[iprod-1],hTrackMatchResEtaNeg[0],1.000,1,errType);
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
    
    if( !scaled )
    {
      hTrackMatchResEtaNegTrackPt[iprod]->Sumw2();
      hTrackMatchResEtaPosTrackPt[iprod]->Sumw2();
      hTrackMatchResPhiNegTrackPt[iprod]->Sumw2();
      hTrackMatchResPhiPosTrackPt[iprod]->Sumw2();
      
      hTrackMatchResEtaNegTrackPt[iprod]->Scale(1./nEvents[iprod]);
      hTrackMatchResEtaPosTrackPt[iprod]->Scale(1./nEvents[iprod]);
      hTrackMatchResPhiNegTrackPt[iprod]->Scale(1./nEvents[iprod]);
      hTrackMatchResPhiPosTrackPt[iprod]->Scale(1./nEvents[iprod]);
    }
    
    hTrackMatchResEtaNegTrackPt[iprod]->SetXTitle("#Delta #eta");
    hTrackMatchResEtaNegTrackPt[iprod]->SetYTitle("Entries / N events");
    hTrackMatchResEtaNegTrackPt[iprod]->SetTitle(Form("Track-cluster #eta residuals, %2.1f < #it{p}_{T}^{track} < %2.1f GeV",emin,emax));
    hTrackMatchResEtaNegTrackPt[iprod]->SetAxisRange(-0.025,0.025,"X");
    hTrackMatchResEtaNegTrackPt[iprod]->SetMarkerStyle(24);
    hTrackMatchResEtaNegTrackPt[iprod]->SetMarkerColor(color[iprod]);
    
    hTrackMatchResEtaPosTrackPt[iprod]->SetYTitle("Entries / N events");
    hTrackMatchResEtaPosTrackPt[iprod]->SetXTitle("#Delta #eta");
    hTrackMatchResEtaPosTrackPt[iprod]->SetTitle(Form("Track-cluster #eta residuals, %2.1f < #it{p}_{T}^{track} < %2.1f GeV",emin,emax));
    hTrackMatchResEtaPosTrackPt[iprod]->SetAxisRange(-0.025,0.025,"X");
    hTrackMatchResEtaPosTrackPt[iprod]->SetMarkerStyle(25);
    hTrackMatchResEtaPosTrackPt[iprod]->SetMarkerColor(color[iprod]);
    
    hTrackMatchResPhiNegTrackPt[iprod]->SetXTitle("#Delta #varphi");
    hTrackMatchResPhiNegTrackPt[iprod]->SetTitle(Form("Track-cluster #varphi residuals, %2.1f < #it{p}_{T}^{track} < %2.1f GeV",emin,emax));
    hTrackMatchResPhiNegTrackPt[iprod]->SetYTitle("entries / N events");
    hTrackMatchResPhiNegTrackPt[iprod]->SetAxisRange(-0.035,0.035,"X");
    hTrackMatchResPhiNegTrackPt[iprod]->SetMarkerStyle(24);
    hTrackMatchResPhiNegTrackPt[iprod]->SetMarkerColor(color[iprod]);
    
    hTrackMatchResPhiPosTrackPt[iprod]->SetYTitle("Entries / N events");
    hTrackMatchResPhiPosTrackPt[iprod]->SetXTitle("#Delta #varphi");
    hTrackMatchResPhiPosTrackPt[iprod]->SetTitle(Form("Track-cluster #varphi residuals, %2.1f < #it{p}_{T}^{track} < %2.1f GeV",emin,emax));
    hTrackMatchResPhiPosTrackPt[iprod]->SetAxisRange(-0.035,0.035,"X");
    hTrackMatchResPhiPosTrackPt[iprod]->SetMarkerStyle(25);
    hTrackMatchResPhiPosTrackPt[iprod]->SetMarkerColor(color[iprod]);
    
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
      
      hRatTrackMatchResPhiPosTrackPt[iprod-1]->Divide(hRatTrackMatchResPhiPosTrackPt[iprod-1],hTrackMatchResPhiPosTrackPt[0],1.000,1,errType);
      hRatTrackMatchResPhiNegTrackPt[iprod-1]->Divide(hRatTrackMatchResPhiNegTrackPt[iprod-1],hTrackMatchResPhiNegTrackPt[0],1.000,1,errType);
      hRatTrackMatchResEtaPosTrackPt[iprod-1]->Divide(hRatTrackMatchResEtaPosTrackPt[iprod-1],hTrackMatchResEtaPosTrackPt[0],1.000,1,errType);
      hRatTrackMatchResEtaNegTrackPt[iprod-1]->Divide(hRatTrackMatchResEtaNegTrackPt[iprod-1],hTrackMatchResEtaNegTrackPt[0],1.000,1,errType);
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
    gPad->SetGridy();
    //gPad->SetLogy();
    //gPad->SetLogx();
    
    hRatCorr[0]->SetTitle("Cluster spectra ratio");
    hRatCorr[0]->SetYTitle(Form("Ratio data X / %s",prodLeg[0].Data()));
    //    hRatCorr[0]->SetMinimum(0.850);
    //    hRatCorr[0]->SetMaximum(1.025);
    hRatCorr[0]->SetMinimum(0.95);
    hRatCorr[0]->SetMaximum(1.05);
    hRatCorr[0]->Draw("");
    
    for(Int_t iprod = 0; iprod <  nProd-1; iprod++)
    {
      hRatRaw [iprod]->Draw("same");
      hRatCorr[iprod]->Draw("same");
      hRatTM  [iprod]->Draw("same");
      hRatShSh[iprod]->Draw("same");
    }
    
    //    TLine l1(0,1,30,1);
    //    TLine l2(0,0.975,30,0.975);
    //    TLine l3(0,0.95,30,0.95);
    //    TLine l4(0,0.925,30,0.925);
    //    
    //    l1.Draw("same");
    //    l2.Draw("same");
    //    l3.Draw("same");
    //    l4.Draw("same");
    
    ccalo->Print(Form("%s_ClusterSpectraComp.%s",histoTag.Data(),format.Data()));
  }
  
  // Cluster-Track Matching Residual
  {
    TGaxis::SetMaxDigits(3);
    
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
    
    Double_t max = hTrackMatchResEtaPos[0]->GetMaximum();
    if(max < hTrackMatchResEtaNeg[0]->GetMaximum()) max = hTrackMatchResEtaNeg[0]->GetMaximum();
    
    hTrackMatchResEtaPos[0]->SetMaximum(max*1.2);
    
    hTrackMatchResEtaPos[0]->Draw("");
    for(Int_t iprod = 0; iprod < nProd; iprod++)
    {
      hTrackMatchResEtaNeg[iprod]->Draw("same");
      hTrackMatchResEtaPos[iprod]->Draw("same");
    }
    
    TLine l0Eta(0,hTrackMatchResEtaPos[0]->GetMinimum(),0,hTrackMatchResEtaPos[0]->GetMaximum());
    l0Eta.SetLineColor(2);
    l0Eta.SetLineWidth(2);
    l0Eta.Draw("same");
    
    lres.Draw();
    lprod.Draw();
    ccalo2->cd(2);
    
    max = hTrackMatchResPhiPos[0]->GetMaximum();
    if(max < hTrackMatchResPhiNeg[0]->GetMaximum()) max = hTrackMatchResPhiNeg[0]->GetMaximum();
    
    hTrackMatchResPhiPos[0]->SetMaximum(max*1.2);
    
    hTrackMatchResPhiPos[0]->Draw("");
    for(Int_t iprod = 0; iprod < nProd; iprod++)
    {
      hTrackMatchResPhiNeg[iprod]->Draw("same");
      hTrackMatchResPhiPos[iprod]->Draw("same");
    }
    
    TLine l0Phi(0,hTrackMatchResPhiPos[0]->GetMinimum(),0,hTrackMatchResPhiPos[0]->GetMaximum());
    l0Phi.SetLineColor(2);
    l0Phi.SetLineWidth(2);
    l0Phi.Draw("same");
    
    ccalo2->cd(3);
    //gPad->SetLogy();
    
    hRatTrackMatchResEtaPos[0]->SetMaximum(1.3);
    hRatTrackMatchResEtaPos[0]->SetMinimum(0.9);
    hRatTrackMatchResEtaPos[0]->Draw("");
    hRatTrackMatchResEtaPos[0]->SetYTitle(Form("Ratio data X / %s",prodLeg[0].Data()));
    for(Int_t iprod = 0; iprod < nProd-1; iprod++)
    {
      hRatTrackMatchResEtaNeg[iprod]->Draw("same");
      hRatTrackMatchResEtaPos[iprod]->Draw("same");
    }
    
    //l0.Draw("same");
    
    ccalo2->cd(4);
    
    hRatTrackMatchResPhiPos[0]->SetMaximum(1.3);
    hRatTrackMatchResPhiPos[0]->SetMinimum(0.9);
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
/// Plot basic calorimeter Cell QA histograms.
/// 2 canvases with 2-4 pads
/// * cell spectra 
/// * number of cells distribution
///
//______________________________________
void Cell()
{
  // Declare the different histograms, arrays input is production
  TH1F* hNC1[nProd];
  TH1F* hNC2[nProd];
  TH1F* hAmp[nProd];

  TH1F* hRatNC1[nProd-1];
  TH1F* hRatNC2[nProd-1];
  TH1F* hRatAmp[nProd-1];
  
  //Legend for productions
  TLegend lprod(0.65,0.525,0.95,0.675);
  lprod.SetTextSize(0.04);
  lprod.SetBorderSize(0);
  lprod.SetFillColor(0);
  
  for(Int_t iprod = 0; iprod <  nProd; iprod++)
  {      
    hAmp[iprod] = (TH1F*) GetHisto("QA_Cell_hAmplitude",iprod);
    hNC1[iprod] = (TH1F*) GetHisto("QA_Cell_hNCells",iprod);
    hNC2[iprod] = (TH1F*) GetHisto("QA_Cell_hNCellsCutAmpMin",iprod);
    
    if(!hAmp[iprod]) return;
    
    hAmp[iprod]->Rebin(2);
    hNC1[iprod]->Rebin(2);    
    hNC2[iprod]->Rebin(2);
    
    if( !scaled )
    {
      hAmp[iprod]->Sumw2();
      hNC1[iprod]->Sumw2();
      hNC2[iprod]->Sumw2();
      
      hAmp[iprod]->Scale(1./nEvents[iprod]);
      hNC1[iprod]->Scale(1./nEvents[iprod]);
      hNC2[iprod]->Scale(1./nEvents[iprod]);
    }
    
    hAmp[iprod]->SetMarkerColor(color[iprod]);
    hAmp[iprod]->SetMarkerStyle(20);
    hAmp[iprod]->SetTitle("Cell energy spectra");
    hAmp[iprod]->SetYTitle("1/N_{events} dN/dp_{T}");
    hAmp[iprod]->SetTitleOffset(1.5,"Y");
    hAmp[iprod]->SetAxisRange(0.,10.,"X");
//    hAmp[iprod]->SetMaximum(5);
//    hAmp[iprod]->SetMinimum(1e-10);
    
    
    hNC1[iprod]->SetTitle("Number of EMCal/DCal cells");
    hNC1[iprod]->SetYTitle("1/N_{events} dN/dp_{T}");
    
    hNC1[iprod]->SetTitleOffset(1.5,"Y");
    hNC1[iprod]->SetMarkerColor(color[iprod]);
    hNC1[iprod]->SetMarkerStyle(20);
    hNC1[iprod]->SetAxisRange(0.,100.,"X");
    hNC1[iprod]->SetMaximum(5);
    hNC1[iprod]->SetMinimum(1e-10);
    
    hNC2[iprod]->SetMarkerColor(color[iprod]);
    hNC2[iprod]->SetMarkerStyle(21);
    
    
    hAmp[iprod]->SetTitleOffset(1.5,"Y");
    hNC2[iprod]->SetTitleOffset(1.5,"Y");
    hNC1[iprod]->SetTitleOffset(1.5,"Y");
    
    if(iprod > 0)
    {
      hRatAmp[iprod-1] = (TH1F*)hAmp[iprod]->Clone(Form("hRatAmp%s_%s",prod[iprod].Data(),histoTag.Data()));
      hRatNC1[iprod-1] = (TH1F*)hNC1[iprod]->Clone(Form("hRatNC1%s_%s",prod[iprod].Data(),histoTag.Data()));
      hRatNC2[iprod-1] = (TH1F*)hNC2[iprod]->Clone(Form("hRatNC2%s_%s",prod[iprod].Data(),histoTag.Data()));

      //      hRatAmp[iprod-1]->Divide(hRatAmp [iprod-1],hAmp [0],1.000,1,errType);
      //      hRatNC1[iprod-1]->Divide(hRatNC1[iprod-1],hNC1[0],0.975,1,errType);
      //      hRatNC2[iprod-1]->Divide(hRatNC2  [iprod-1],hNC2  [0],0.950,1,errType);
      
      hRatAmp[iprod-1]->Divide(hRatAmp[iprod-1],hAmp[0],1.000,1,errType);
      hRatNC1[iprod-1]->Divide(hRatNC1[iprod-1],hNC1[0],1.000,1,errType);
      hRatNC2[iprod-1]->Divide(hRatNC2[iprod-1],hNC2[0],1.000,1,errType);
    }
  
  } // prod loop

  
  /////////////////
  // Make the plots
  /////////////////
  {
    TCanvas * cn = new TCanvas(Form("NCell_%s",histoTag.Data()),"",1000,500);
    cn->Divide(2,1);
    
    cn->cd(1);
    gPad->SetLogy();
    //gPad->SetLogx();
    
    hNC1[0]->Draw();
    for(Int_t iprod = 0; iprod <  nProd; iprod++)
    {
      hNC1[iprod]->Draw("same");
      hNC2[iprod]->Draw("same");
      
      lprod.AddEntry(hNC1[iprod],prodLeg[iprod],"P");
    }
    
    lprod.Draw();
    
    TLegend lcl(0.55,0.7,0.95,0.89);
    lcl.SetTextSize(0.04);
    lcl.SetBorderSize(0);
    lcl.SetFillColor(0);
    lcl.AddEntry(hNC1[0],"E >  50 MeV","P");
    lcl.AddEntry(hNC2[0],"E > 200 MeV","P");
    lcl.Draw();
    
    cn->cd(2);
    gPad->SetGridy();
    //gPad->SetLogy();
    //gPad->SetLogx();
    
    hRatNC1[0]->SetTitle("N cell ratio");
    hRatNC1[0]->SetYTitle(Form("Ratio data X / %s",prodLeg[0].Data()));
//    hRatNC1[0]->SetMinimum(0.850);
//    hRatNC1[0]->SetMaximum(1.025);
    hRatNC1[0]->SetMinimum(0.7);
    hRatNC1[0]->SetMaximum(1.3);
    hRatNC1[0]->Draw("");
    
    for(Int_t iprod = 0; iprod <  nProd-1; iprod++)
    {
      hRatNC1[iprod]->Draw("same");
      hRatNC2[iprod]->Draw("same");
    }
    
    cn->Print(Form("%s_NCellComp.%s",histoTag.Data(),format.Data()));
  }

  {
    TCanvas * cE = new TCanvas(Form("Cell_Energy_%s",histoTag.Data()),"",1000,500);
    cE->Divide(2,1);
    
    cE->cd(1);
    gPad->SetLogy();
    //gPad->SetLogx();
    
    hAmp[0]->Draw();
    for(Int_t iprod = 0; iprod <  nProd; iprod++)
    {
      hAmp[iprod]->Draw("same");      
    }
    
    lprod.Draw();
    
   
    cE->cd(2);
    gPad->SetGridy();
    //gPad->SetLogy();
    //gPad->SetLogx();
    
    hRatAmp[0]->SetTitle("Cell E spectra ratio");
    hRatAmp[0]->SetYTitle(Form("Ratio data X / %s",prodLeg[0].Data()));
    //    hRatAmp[0]->SetMinimum(0.850);
    //    hRatAmp[0]->SetMaximum(1.025);
    hRatAmp[0]->SetMinimum(0.7);
    hRatAmp[0]->SetMaximum(1.3);
    hRatAmp[0]->Draw("");
    
    for(Int_t iprod = 0; iprod <  nProd-1; iprod++)
    {
      hRatAmp[iprod]->Draw("same");
    }
    
    cE->Print(Form("%s_CellEnergyComp.%s",histoTag.Data(),format.Data()));
  }

}

///
/// Plot basic shower shape
///
/// \param icalo: 0 EMCal, 1 DCal
//______________________________________
void ClusterCells(Int_t icalo)
{
  // Declare the different histograms, arrays input is production
  
  const Int_t nEbins = 4;
  Float_t ebins [] = {2,4,6,8,10,12,16,20};
  
  TH2F* h2NCell  [nProd];
  TH1F* hNCell   [nProd][nEbins]; 
  TH1F* hNCellRat[nProd][nEbins]; 
  TH2F* h2ECell  [nProd];
  TH1F* hECell   [nProd][nEbins];
  TH1F* hECellRat[nProd][nEbins];
  
  //Legend for productions
  TLegend lprod(0.6,0.7,0.95,0.89);
  lprod.SetTextSize(0.04);
  lprod.SetBorderSize(0);
  lprod.SetFillColor(0);
  
  for(Int_t iprod = 0; iprod <  nProd; iprod++)
  {      
    // E cluster bin
    h2NCell[iprod] = (TH2F*) GetHisto(Form("AnaPhoton_Calo%d_hNCellsE",icalo),iprod);
    h2ECell[iprod] = (TH2F*) GetHisto(Form("AnaPhoton_Calo%d_hCellsE" ,icalo),iprod);
    
    if(!h2NCell[iprod] || !h2ECell[iprod]) return;
    
    for(Int_t ie = 0; ie < nEbins; ie++)
    {
      Float_t binMin = h2NCell[iprod]->GetXaxis()->FindBin(ebins[ie]);
      Float_t binMax = h2NCell[iprod]->GetXaxis()->FindBin(ebins[ie+1])-1;
      
      hNCell[iprod][ie] = 
      (TH1F*) h2NCell[iprod]->ProjectionY(Form("DeltaProjNCell_%s_MC%s_ie%d",
                                               prod[iprod].Data(),histoTag.Data(),ie),
                                          binMin, binMax);
      
      hECell[iprod][ie] = 
      (TH1F*) h2ECell[iprod]->ProjectionY(Form("DeltaProjECell_%s_MC%s_ie%d",
                                               prod[iprod].Data(),histoTag.Data(),ie),
                                          binMin, binMax);
      
      //hNCell[iprod][ie]->SetXTitle("#it{E}_{reco}-#it{E}_{gen} (GeV)");
      hNCell[iprod][ie]->SetYTitle("Entries / N events");
      hNCell[iprod][ie]->SetTitle(Form("%2.1f < #it{E}^{cluster} < %2.1f GeV",ebins[ie],ebins[ie+1]));
      hNCell[iprod][ie]->SetAxisRange(0.,15,"X");
      //hNCell[iprod][ie]->SetMarkerStyle(24);
      hNCell[iprod][ie]->SetMarkerColor(color[iprod]);
      hNCell[iprod][ie]->SetLineColor  (color[iprod]);
      
      if( !scaled )
      {
        hNCell[iprod][ie]->Sumw2();
        hNCell[iprod][ie]->Scale(1./nEvents[iprod]);
      }
      
      hNCell[iprod][ie]->SetTitleOffset(1.5,"Y");
      
      //hECell[iprod][ie]->SetXTitle("#it{E}_{reco}-#it{E}_{gen} (GeV)");
      hECell[iprod][ie]->SetYTitle("Entries / N events");
      hECell[iprod][ie]->SetTitle(Form("%2.1f < #it{E}^{cluster} < %2.1f GeV",ebins[ie],ebins[ie+1]));
      hECell[iprod][ie]->SetAxisRange(0.0,ebins[ie+1],"X");
      //hECell[iprod][ie]->SetMarkerStyle(24);
      hECell[iprod][ie]->SetMarkerColor(color[iprod]);
      hECell[iprod][ie]->SetLineColor  (color[iprod]);
      
      if( !scaled )
      {
        hECell[iprod][ie]->Sumw2();
        hECell[iprod][ie]->Scale(1./nEvents[iprod]);
      }
      
      hECell[iprod][ie]->SetTitleOffset(1.5,"Y");
      
      if(iprod > 0)
      {
        hECellRat[iprod][ie] = (TH1F*) hECell[iprod][ie]->Clone(Form("Ratio_%s",hECell[iprod][ie]->GetName()));
        hNCellRat[iprod][ie] = (TH1F*) hNCell[iprod][ie]->Clone(Form("Ratio_%s",hNCell[iprod][ie]->GetName()));
        
        hECellRat[iprod][ie]->Divide(hECell[iprod][ie],hECell[0][ie],1,1,errType);
        hNCellRat[iprod][ie]->Divide(hNCell[iprod][ie],hNCell[0][ie],1,1,errType);
      }
      else
      {
        hECellRat[iprod][ie] = 0;  
        hNCellRat[iprod][ie] = 0;  
      }
    }
    
    lprod.AddEntry(hNCell[iprod][0],prodLeg[iprod],"LP");
    
  } // prod loop
  
  
  /////////////////
  // Make the plots
  /////////////////
  
  {
    TGaxis::SetMaxDigits(3);
    
    TCanvas * cNCell = new TCanvas(Form("NClusterCell_%s",histoTag.Data()),"",1000,1000);
    cNCell->Divide(2,2);
    
    for(Int_t ie = 0; ie < nEbins; ie++)
    {
      cNCell->cd(ie+1);
      gPad->SetLogy();
      
      hNCell[0][ie]->Draw("H");
      for(Int_t iprod = 0; iprod < nProd; iprod++)
      {
        hNCell[iprod][ie]->Draw("Hsame");
      }
      
      lprod.Draw();
    }
    
    cNCell->Print(Form("%s_NCell.%s",histoTag.Data(),format.Data()));
    
    TCanvas * cECell = new TCanvas(Form("ECell_%s",histoTag.Data()),"",1000,1000);
    cECell->Divide(2,2);
    
    for(Int_t ie = 0; ie < nEbins; ie++)
    {
      cECell->cd(ie+1);
      gPad->SetLogy();
      
      hECell[0][ie]->Draw("H");
      for(Int_t iprod = 0; iprod < nProd; iprod++)
      {
        hECell[iprod][ie]->Draw("Hsame");
      }
      
      lprod.Draw();
    }
    
    cECell->Print(Form("%s_ECell.%s",histoTag.Data(),format.Data()));

    // RATIOS
    
    TCanvas * cNCellR = new TCanvas(Form("RatioNCell_%s",histoTag.Data()),"",1000,1000);
    cNCellR->Divide(2,2);
    
    for(Int_t ie = 0; ie < nEbins; ie++)
    {
      cNCellR->cd(ie+1);
      //gPad->SetLogy();
      
      hNCellRat[1][ie]->Draw("H");
      hNCellRat[1][ie]->SetYTitle(Form("Ratio data X / %s",prodLeg[0].Data()));
      //hNCellRat[1][ie]->SetMaximum(1.2);
      //hNCellRat[1][ie]->SetMinimum(0.6);
      for(Int_t iprod = 1; iprod < nProd; iprod++)
      {
        hNCellRat[iprod][ie]->Draw("Hsame");
      }
      
      lprod.Draw();
    }
    
    cNCellR->Print(Form("%s_NCell_Ratio.%s",histoTag.Data(),format.Data()));
    
    TCanvas * cECellR = new TCanvas(Form("RatioECell_%s",histoTag.Data()),"",1000,1000);
    cECellR->Divide(2,2);
    
    for(Int_t ie = 0; ie < nEbins; ie++)
    {
      cECellR->cd(ie+1);
      //gPad->SetLogy();
      
      hECellRat[1][ie]->Draw("H");
      hECellRat[1][ie]->SetYTitle(Form("Ratio data X / %s",prodLeg[0].Data()));
      hECellRat[1][ie]->SetMaximum(1.2);
      hECellRat[1][ie]->SetMinimum(0.6);
      for(Int_t iprod = 1; iprod < nProd; iprod++)
      {
        hECellRat[iprod][ie]->Draw("Hsame");
      }
      
      lprod.Draw();
    }
    
    cECellR->Print(Form("%s_ECell_Ratio.%s",histoTag.Data(),format.Data()));
    
  }
}


///
/// Plot basic calorimeter QA histograms depending on MC origin of cluster
/// * cluster spectra from particle
/// * reconstructed minus generated energy of cluster particle
///
/// \param icalo: 0 EMCal, 1 DCal
/// \param particle: MC origin of cluster: Photon, PhotonPi0Decay, Pi0 (Merged), Electron, ...
//______________________________________
void MCParticle(Int_t icalo, TString particle)
{
  // Declare the different histograms, arrays input is production
  
  const Int_t nEbins = 4;
  Float_t ebins [] = {1,2,4,6,10};
  
  TH1F* hParticleRecoE[nProd];
  TH1F* hRatParticleRecoE [nProd-1];
  
  TH1F* hParticleGeneE[nProd];
  TH1F* hRatParticleGeneE [nProd-1];

  TH2F* h2ParticleRecoEDelta[nProd];
  TH1F* hParticleRecoEDelta[nProd][nEbins];
  
  //Legend for productions
  TLegend lprod(0.6,0.475,0.95,0.675);
  lprod.SetTextSize(0.04);
  lprod.SetBorderSize(0);
  lprod.SetFillColor(0);

  //Legend for object
  TLegend legPa(0.6,0.7,0.95,0.89);
  legPa.SetTextSize(0.04);
  legPa.SetBorderSize(0);
  legPa.SetFillColor(0);
  
  for(Int_t iprod = 0; iprod <  nProd; iprod++)
  {      
    hParticleRecoE [iprod] = (TH1F*) GetHisto(Form("AnaPhoton_Calo%d_hE_MC%s",icalo,particle.Data()),iprod);
    hParticleGeneE [iprod] = (TH1F*) GetHisto(Form("AnaPhoton_Calo%d_hPtPrim_MC%s",icalo,particle.Data()),iprod);
    
    if(!hParticleRecoE[iprod]) return;
    
    hParticleRecoE[iprod]->Rebin(5);
    hParticleGeneE[iprod]->Rebin(5);
    if( !scaled )
    {
      hParticleRecoE[iprod]->Sumw2();
      hParticleRecoE[iprod]->Scale(1./nEvents[iprod]);      
      
      hParticleGeneE[iprod]->Sumw2();
      hParticleGeneE[iprod]->Scale(1./nEvents[iprod]);
    }
    
    hParticleRecoE[iprod]->SetTitle(Form("Generated and reconstructed %s spectra",particle.Data()));    
    hParticleRecoE[iprod]->SetYTitle("1/N_{events} dN/dp_{T}");
    hParticleRecoE[iprod]->SetTitleOffset(1.5,"Y");
    hParticleRecoE[iprod]->SetAxisRange(10.,70.,"X");
    //hParticleRecoE[iprod]->SetMaximum(1.1);
    //hParticleRecoE[iprod]->SetMinimum(0);
    
    hParticleRecoE[iprod]->SetMarkerColor(color[iprod]);
    hParticleRecoE[iprod]->SetMarkerStyle(20);
    
    hParticleGeneE[iprod]->SetMarkerColor(color[iprod]);
    hParticleGeneE[iprod]->SetMarkerStyle(24);
   
    if(iprod > 0)
    {
      hRatParticleRecoE [iprod-1] = (TH1F*)hParticleRecoE [iprod]->Clone(Form("hRatParticleRecoE%s_MC%s_%s" ,prod[iprod].Data(),particle.Data(),histoTag.Data()));
      
      hRatParticleRecoE [iprod-1]->Divide(hRatParticleRecoE [iprod-1],hParticleRecoE [0],1,1,errType);
      
      hRatParticleGeneE [iprod-1] = (TH1F*)hParticleGeneE [iprod]->Clone(Form("hRatParticleGeneE%s_MC%s_%s" ,prod[iprod].Data(),particle.Data(),histoTag.Data()));
      
      hRatParticleGeneE [iprod-1]->Divide(hRatParticleGeneE [iprod-1],hParticleGeneE [0],1,1,errType);
    }
    
    // Cluster-Track Matching Residuals
    
    // E cluster bin
    h2ParticleRecoEDelta[iprod] = (TH2F*) GetHisto(Form("AnaPhoton_Calo%d_hDeltaE_MC%s",icalo,particle.Data()),iprod);
    
    for(Int_t ie = 0; ie < nEbins; ie++)
    {
      Float_t binMin = h2ParticleRecoEDelta[iprod]->GetXaxis()->FindBin(ebins[ie]);
      Float_t binMax = h2ParticleRecoEDelta[iprod]->GetXaxis()->FindBin(ebins[ie+1])-1;
      
      hParticleRecoEDelta[iprod][ie] = 
      (TH1F*) h2ParticleRecoEDelta[iprod]->ProjectionY(Form("DeltaProj_%s_%s_MC%s_ie%d",
                                                            prod[iprod].Data(),histoTag.Data(),particle.Data(),ie),
                                                       binMin, binMax);
      
      
      hParticleRecoEDelta[iprod][ie]->SetXTitle("#it{E}_{reco}-#it{E}_{gen} (GeV)");
      hParticleRecoEDelta[iprod][ie]->SetYTitle("Entries / N events");
      hParticleRecoEDelta[iprod][ie]->SetTitle(Form("Reco MC %s, %2.1f < #it{E}^{cluster} < %2.1f GeV",particle.Data(),ebins[ie],ebins[ie+1]));
      hParticleRecoEDelta[iprod][ie]->SetAxisRange(-3,3,"X");
      hParticleRecoEDelta[iprod][ie]->SetMarkerStyle(24);
      hParticleRecoEDelta[iprod][ie]->SetMarkerColor(color[iprod]);
      
      if( !scaled )
      {     
        hParticleRecoEDelta[iprod][ie]->Sumw2();
        hParticleRecoEDelta[iprod][ie]->Scale(1./nEvents[iprod]);
      }
      hParticleRecoEDelta[iprod][ie]->SetTitleOffset(1.5,"Y");
    }
    
  } // prod loop
  
  
  /////////////////
  // Make the plots
  /////////////////
  {
    TCanvas * ccalo = new TCanvas(Form("Cluster_MC%s_%s",particle.Data(),histoTag.Data()),"",1000,500);
    ccalo->Divide(2,1);
    
    ccalo->cd(1);
    gPad->SetLogy();
    //gPad->SetLogx();
    
    hParticleRecoE[0]->Draw();
    for(Int_t iprod = 0; iprod <  nProd; iprod++)
    {
      hParticleRecoE [iprod]->Draw("same");
      hParticleGeneE [iprod]->Draw("same");
    
      lprod.AddEntry(hParticleRecoE[iprod],prodLeg[iprod],"P");
      if(iprod == 0)
      {
        legPa.AddEntry(hParticleRecoE[iprod],"Reconstructed","P");
        legPa.AddEntry(hParticleGeneE[iprod],"Generated","P");
      }
    }
    
    lprod.Draw();
    legPa.Draw();
    
    ccalo->cd(2);
    gPad->SetGridy();
    //gPad->SetLogy();
    //gPad->SetLogx();
    
    hRatParticleRecoE[0]->SetTitle(Form("MC %s spectra ratio",particle.Data()));
    hRatParticleRecoE[0]->SetYTitle(Form("Ratio data X / %s",prodLeg[0].Data()));
    //    hRatParticleRecoE[0]->SetMinimum(0.850);
    //    hRatParticleRecoE[0]->SetMaximum(1.025);
    hRatParticleRecoE[0]->SetMinimum(0.95);
    hRatParticleRecoE[0]->SetMaximum(1.05);
    hRatParticleRecoE[0]->Draw("");
    
    for(Int_t iprod = 0; iprod <  nProd-1; iprod++)
    {
      hRatParticleRecoE[iprod]->Draw("same");
      hRatParticleGeneE[iprod]->Draw("same");
    }
    
    ccalo->Print(Form("%s_ClusterSpectra_MC%s.%s",histoTag.Data(),particle.Data(),format.Data()));
  }
  
  // Cluster energy residual
  {
    TGaxis::SetMaxDigits(3);
    
    TCanvas * ccalo2 = new TCanvas(Form("Resolution_MC%s_%s",particle.Data(),histoTag.Data()),"",1000,1000);
    ccalo2->Divide(2,2);
    
    for(Int_t ie = 0; ie < nEbins; ie++)
    {
      ccalo2->cd(ie+1);
      //gPad->SetLogy();
          
      hParticleRecoEDelta[0][ie]->Draw("");
      for(Int_t iprod = 0; iprod < nProd; iprod++)
      {
        hParticleRecoEDelta[iprod][ie]->Draw("same");
      }
      
      lprod.Draw();
    }
    
    ccalo2->Print(Form("%s_ClusterERecoGenDiff_MC%s.%s",histoTag.Data(),particle.Data(),format.Data()));
  }
  
}

///
/// Plot basic shower shape
///
/// \param icalo: 0 EMCal, 1 DCal
//______________________________________
void ShowerShape(Int_t icalo)
{
  // Declare the different histograms, arrays input is production
  
  const Int_t nEbins = 7;
  Float_t ebins [] = {2,4,6,8,10,12,16,20};
  
  TH2F* h2M02[nProd];
  TH1F* hM02 [nProd][nEbins]; 
  TH2F* h2M20[nProd];
  TH1F* hM20 [nProd][nEbins];
  
  //Legend for productions
  TLegend lprod(0.6,0.475,0.95,0.675);
  lprod.SetTextSize(0.04);
  lprod.SetBorderSize(0);
  lprod.SetFillColor(0);
  
  for(Int_t iprod = 0; iprod <  nProd; iprod++)
  {      
    // E cluster bin
    h2M02[iprod] = (TH2F*) GetHisto(Form("AnaPhoton_Calo%d_hLam0E",icalo),iprod);
    h2M20[iprod] = (TH2F*) GetHisto(Form("AnaPhoton_Calo%d_hLam1E",icalo),iprod);
    
    if(!h2M02[iprod] || !h2M20[iprod]) return;
    
    for(Int_t ie = 0; ie < nEbins; ie++)
    {
      Float_t binMin = h2M02[iprod]->GetXaxis()->FindBin(ebins[ie]);
      Float_t binMax = h2M02[iprod]->GetXaxis()->FindBin(ebins[ie+1])-1;
      
      hM02[iprod][ie] = 
      (TH1F*) h2M02[iprod]->ProjectionY(Form("DeltaProjM02_%s_MC%s_ie%d",
                                                            prod[iprod].Data(),histoTag.Data(),ie),
                                                       binMin, binMax);
 
      hM20[iprod][ie] = 
      (TH1F*) h2M20[iprod]->ProjectionY(Form("DeltaProjM20_%s_MC%s_ie%d",
                                             prod[iprod].Data(),histoTag.Data(),ie),
                                        binMin, binMax);
      
      //hM02[iprod][ie]->SetXTitle("#it{E}_{reco}-#it{E}_{gen} (GeV)");
      hM02[iprod][ie]->SetYTitle("Entries / N events");
      hM02[iprod][ie]->SetTitle(Form("%2.1f < #it{E}^{cluster} < %2.1f GeV",ebins[ie],ebins[ie+1]));
      hM02[iprod][ie]->SetAxisRange(0.1,0.8,"X");
      //hM02[iprod][ie]->SetMarkerStyle(24);
      hM02[iprod][ie]->SetMarkerColor(color[iprod]);
      hM02[iprod][ie]->SetLineColor  (color[iprod]);
      
      if( !scaled )
      {
        hM02[iprod][ie]->Sumw2();
        hM02[iprod][ie]->Scale(1./nEvents[iprod]);
      }
      
      hM02[iprod][ie]->SetTitleOffset(1.5,"Y");
      
      
      //hM20[iprod][ie]->SetXTitle("#it{E}_{reco}-#it{E}_{gen} (GeV)");
      hM20[iprod][ie]->SetYTitle("Entries / N events");
      hM20[iprod][ie]->SetTitle(Form("%2.1f < #it{E}^{cluster} < %2.1f GeV",ebins[ie],ebins[ie+1]));
      hM20[iprod][ie]->SetAxisRange(0.0,0.5,"X");
      //hM20[iprod][ie]->SetMarkerStyle(24);
      hM20[iprod][ie]->SetMarkerColor(color[iprod]);
      hM20[iprod][ie]->SetLineColor  (color[iprod]);
      
      if( !scaled )
      {
        hM20[iprod][ie]->Sumw2();
        hM20[iprod][ie]->Scale(1./nEvents[iprod]);
      }
      
      hM20[iprod][ie]->SetTitleOffset(1.5,"Y");
    }
    
    lprod.AddEntry(hM02[iprod][0],prodLeg[iprod],"LP");
    
  } // prod loop
  
  
  /////////////////
  // Make the plots
  /////////////////
 
  {
    TGaxis::SetMaxDigits(3);
    
    TCanvas * cM02 = new TCanvas(Form("M02_%s",histoTag.Data()),"",1000,1000);
    cM02->Divide(2,2);
    
    for(Int_t ie = 0; ie < nEbins; ie++)
    {
      cM02->cd(ie+1);
      //gPad->SetLogy();
      
      hM02[0][ie]->Draw("H");
      for(Int_t iprod = 0; iprod < nProd; iprod++)
      {
        hM02[iprod][ie]->Draw("Hsame");
      }
      
      lprod.Draw();
    }
    
    cM02->Print(Form("%s_M02.%s",histoTag.Data(),format.Data()));
 
    TCanvas * cM20 = new TCanvas(Form("M20_%s",histoTag.Data()),"",1000,1000);
    cM20->Divide(2,2);
    
    for(Int_t ie = 0; ie < nEbins; ie++)
    {
      cM20->cd(ie+1);
      //gPad->SetLogy();
      
      hM20[0][ie]->Draw("H");
      for(Int_t iprod = 0; iprod < nProd; iprod++)
      {
        hM20[iprod][ie]->Draw("Hsame");
      }
      
      lprod.Draw();
    }
    
    cM20->Print(Form("%s_M20.%s",histoTag.Data(),format.Data()));
    
  }
  
}



///
/// Hybrid Tracks distributions
/// To be updated
//______________________________________
void Track()
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
  TLegend lprod(0.5,0.475,0.84,0.675);
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
    
    hTrackPt     [iprod]->Rebin(5);
    hTrackPtSPD  [iprod]->Rebin(5);
    hTrackPtNoSPD[iprod]->Rebin(5);
    
    hTrackPhi     [iprod]->Rebin(5);
    hTrackPhiSPD  [iprod]->Rebin(5);
    hTrackPhiNoSPD[iprod]->Rebin(5);
    
    if( !scaled )
    {
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
    }
    
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
    
    lprod.AddEntry(hTrackPt[iprod],prodLeg[iprod],"P");
    
    if(iprod > 0)
    {
      hRatTrackPhi     [iprod-1] = (TH1F*)hTrackPhi     [iprod]->Clone(Form("hRatTrackPhi%s"     ,prod[iprod].Data()));
      hRatTrackPhiNoSPD[iprod-1] = (TH1F*)hTrackPhiNoSPD[iprod]->Clone(Form("hRatTrackPhiNoSPD%s",prod[iprod].Data()));
      hRatTrackPhiSPD  [iprod-1] = (TH1F*)hTrackPhiSPD  [iprod]->Clone(Form("hRatTrackPhiSPD%s"  ,prod[iprod].Data()));
      
      hRatTrackPhi     [iprod-1]->Divide(hRatTrackPhi     [iprod-1],hTrackPhi     [0],1.000,1,errType);
      hRatTrackPhiSPD  [iprod-1]->Divide(hRatTrackPhiSPD  [iprod-1],hTrackPhiSPD  [0],1.000,1,errType);
      hRatTrackPhiNoSPD[iprod-1]->Divide(hRatTrackPhiNoSPD[iprod-1],hTrackPhiNoSPD[0],1.000,1,errType);
      
      hRatTrackPt     [iprod-1] = (TH1F*)hTrackPt     [iprod]->Clone(Form("hRatTrackPt%s"     ,prod[iprod].Data()));
      hRatTrackPtNoSPD[iprod-1] = (TH1F*)hTrackPtNoSPD[iprod]->Clone(Form("hRatTrackPtNoSPD%s",prod[iprod].Data()));
      hRatTrackPtSPD  [iprod-1] = (TH1F*)hTrackPtSPD  [iprod]->Clone(Form("hRatTrackPtSPD%s"  ,prod[iprod].Data()));
      
      hRatTrackPt     [iprod-1]->Divide(hRatTrackPt     [iprod-1],hTrackPt     [0],1.000,1,errType);
      hRatTrackPtSPD  [iprod-1]->Divide(hRatTrackPtSPD  [iprod-1],hTrackPtSPD  [0],1.000,1,errType);
      hRatTrackPtNoSPD[iprod-1]->Divide(hRatTrackPtNoSPD[iprod-1],hTrackPtNoSPD[0],1.000,1,errType);
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
    //hTrackPt[0]->SetMinimum(1e-8);
    for(Int_t iprod = 0; iprod < nProd; iprod++)
    {
      hTrackPt     [iprod]->Draw("same");
      hTrackPtSPD  [iprod]->Draw("same");
      hTrackPtNoSPD[iprod]->Draw("same");
    }
    
    ltrack.Draw();
    lprod.Draw();
    
    ctrack->cd(2);
    gPad->SetGridy();

    hRatTrackPt[0]->SetMaximum(1.2);
    hRatTrackPt[0]->SetMinimum(0.8);
    hRatTrackPt[0]->Draw("");
    hRatTrackPt[0]->SetYTitle(Form("Ratio data X / %s",prodLeg[0].Data()));
    for(Int_t iprod = 0; iprod < nProd-1; iprod++)
    {
      hRatTrackPt     [iprod]->Draw("same");
      hRatTrackPtSPD  [iprod]->Draw("same");
      hRatTrackPtNoSPD[iprod]->Draw("same");
    }
    
    ctrack->cd(3);
    //hTrackPhi[0]->SetMaximum(0.2);
    //hTrackPhi[0]->SetMinimum(0.);
    hTrackPhi[0]->Draw("");
    for(Int_t iprod = 0; iprod < nProd; iprod++)
    {
      hTrackPhi     [iprod]->Draw("same");
      hTrackPhiSPD  [iprod]->Draw("same");
      hTrackPhiNoSPD[iprod]->Draw("same");
    }
    
    ctrack->cd(4);
    gPad->SetGridy();
    
    hRatTrackPhi[0]->SetMaximum(1.02);
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
/// Pt Hard distributions
//______________________________________
void PtHard()
{ 
  TH1F * hHardPt[nProd] ;
  TH1F * hRatHardPt[nProd-1] ;
  
  //Legend for productions
  TLegend lprod(0.5,0.475,0.84,0.675);
  lprod.SetTextSize(0.04);
  lprod.SetBorderSize(0);
  lprod.SetFillColor(0);
  
  for(Int_t iprod = 0; iprod <  nProd; iprod++)
  {  
    hHardPt         [iprod] = (TH1F*) GetHisto("hPtHard",iprod);
    if ( !hHardPt[iprod] ) return;
    
    hHardPt     [iprod]->Rebin(2);
          
    //hHardPt[iprod]->SetTitle("pT hard");
    hHardPt[iprod]->SetYTitle("1/N_{events} dN/dp_{T}");
    hHardPt[iprod]->SetTitleOffset(1.5,"Y");
    hHardPt[iprod]->SetMarkerColor(color[iprod]);
    hHardPt[iprod]->SetMarkerStyle(20);
    hHardPt[iprod]->SetAxisRange(5.,100.,"X");
    //hHardPt[iprod]->SetMaximum(1.1);
    //hHardPt[iprod]->SetMinimum(0);
    
    lprod.AddEntry(hHardPt[iprod],prodLeg[iprod],"P");
    
    if(iprod > 0)
    {
      hRatHardPt     [iprod-1] = (TH1F*)hHardPt[iprod]->Clone(Form("hRatHardPt%s"     ,prod[iprod].Data()));
      hRatHardPt     [iprod-1]->Divide(hRatHardPt[iprod-1],hHardPt     [0],1.000,1,errType);
    }
  }
  
  /////////////////
  // Make the plots
  /////////////////

  TLegend lhard(0.6,0.75,0.84,0.89);
  lhard.SetTextSize(0.04);
  //lhard.SetBorderSize(0);
  lhard.SetFillColor(0);
  
  TCanvas * chard = new TCanvas(Form("PtHardHisto_%s",histoTag.Data()),"",1000,500);
  chard->Divide(2,1);
  chard->cd(1);
  gPad->SetLogy();
  hHardPt[0]->Draw("");
  //hHardPt[0]->SetMinimum(1e-8);
  for(Int_t iprod = 0; iprod < nProd; iprod++)
  {
    hHardPt[iprod]->Draw("same");
  }
  
  lprod.Draw();
  
  chard->cd(2);
  gPad->SetGridy();
  
  hRatHardPt[0]->SetMaximum(1.05);
  hRatHardPt[0]->SetMinimum(0.95);
  hRatHardPt[0]->Draw("");
  hRatHardPt[0]->SetYTitle(Form("Ratio data X / %s",prodLeg[0].Data()));
  for(Int_t iprod = 0; iprod < nProd-1; iprod++)
  {
    hRatHardPt[iprod]->Draw("same");
  }

  chard->Print(Form("%s_PtHardComp.%s",histoTag.Data(),format.Data()));
  
}


///
/// cluster-track correlation
/// 
/// To be updated
///
/// \param icalo: 0 EMCal, 1 DCal
//______________________________________
void Correl(Int_t icalo)
{
  TH2F* h2XE[nProd];
  TH2F* h2XEUE[nProd];
  TH1F* hXE[nProd];
  TH1F* hXEUE[nProd];
  TH1F* hRatXE[nProd-1];
  TH1F* hRatXEUE[nProd-1];
  
  //Legend for productions
  TLegend lprod(0.5,0.475,0.84,0.675);
  lprod.SetTextSize(0.04);
  lprod.SetBorderSize(0);
  lprod.SetFillColor(0);
  
  for(Int_t iprod = 0; iprod <  nProd; iprod++)
  {
    h2XE   [iprod]= (TH2F*) GetHisto(Form("AnaPhotonHadronCorr_Calo%d_hXECharged"  ,icalo),iprod);
    h2XEUE [iprod]= (TH2F*) GetHisto(Form("AnaPhotonHadronCorr_Calo%d_hXEUeCharged",icalo),iprod);
    
    if ( !h2XE[iprod] ) return;
    
    Float_t minClusterE = 5;
    TH1F * hTrigger = (TH1F*) GetHisto(Form("AnaPhotonHadronCorr_Calo%d_hPtTrigger",icalo),iprod);
    if(!hTrigger)
    {
      printf("Null trigger histo for prod %d\n",iprod);
      return;
    }
    Int_t minClusterEBin = hTrigger->FindBin(minClusterE);
    Float_t nTrig = hTrigger->Integral(minClusterE,100000);
    
    hXE  [iprod] = (TH1F*)h2XE  [iprod]->ProjectionY(Form("hXE%s_%s"  ,prod[iprod].Data(),histoTag.Data()),minClusterEBin,1000);
    hXEUE[iprod] = (TH1F*)h2XEUE[iprod]->ProjectionY(Form("hXEUE%s_%s",prod[iprod].Data(),histoTag.Data()),minClusterEBin,1000);
    
    hXE  [iprod]->Rebin(5);
    hXEUE[iprod]->Rebin(5);
    
    if( !scaled )
    {    
      hXE  [iprod]->Sumw2();
      hXEUE[iprod]->Sumw2();
    }
    
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
    
    lprod.AddEntry(hXE[iprod],prodLeg[iprod],"P");
    
    hXEUE[iprod]->SetMarkerColor(color[iprod]);
    hXEUE[iprod]->SetMarkerStyle(25);
    
    if(iprod > 0)
    {
      hRatXE  [iprod-1] = (TH1F*)hXE  [iprod]->Clone(Form("hRatXE%s_%s"  ,prod[iprod].Data(),histoTag.Data()));
      hRatXEUE[iprod-1] = (TH1F*)hXEUE[iprod]->Clone(Form("hRatXEUE%s_%s",prod[iprod].Data(),histoTag.Data()));
      
      hRatXE  [iprod-1]->Divide(hRatXE  [iprod-1],hXE  [0],1.000,1,errType);
      hRatXEUE[iprod-1]->Divide(hRatXEUE[iprod-1],hXEUE[0],1.000,1,errType);
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
    
    hRatXE[0]->SetMaximum(1.5);
    hRatXE[0]->SetMinimum(0.5);
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


///
/// isolation cone
/// 
/// To be updated
///
/// \param icalo: 0 EMCal, 1 DCal
//______________________________________
void Isol(Int_t icalo)
{
  TH2F* h2Ne[nProd];
  TH2F* h2Ch[nProd];
  TH1F* hNe[nProd];
  TH1F* hCh[nProd];
  TH1F* hRatNe[nProd-1];
  TH1F* hRatCh[nProd-1];
  
  //Legend for productions
  TLegend lprod(0.15,0.15,0.35,0.35);
  lprod.SetTextSize(0.04);
  lprod.SetBorderSize(0);
  lprod.SetFillColor(0);
  
  for(Int_t iprod = 0; iprod <  nProd; iprod++)
  {
    h2Ne [iprod]= (TH2F*) GetHisto(Form("AnaIsolPhoton_Calo%d_hConePtSumCluster",icalo),iprod);
    h2Ch [iprod]= (TH2F*) GetHisto(Form("AnaIsolPhoton_Calo%d_hConePtSumTrack"  ,icalo),iprod);
    
    if ( !h2Ne[iprod] ) return;
    
    Float_t minClusterE = 5;
    TH1F * hTrigger = (TH1F*) GetHisto(Form("AnaIsolPhoton_Calo%d_hE",icalo),iprod);
    hTrigger    ->Add((TH1F*) GetHisto(Form("AnaIsolPhoton_Calo%d_hENoIso",icalo),iprod));
    if(!hTrigger)
    {
      printf("Null trigger histo for prod %d\n",iprod);
      return;
    }
    
    Int_t minClusterEBin = hTrigger->FindBin(minClusterE);
    Float_t nTrig = hTrigger->Integral(minClusterE,100000);
    
    hNe[iprod] = (TH1F*)h2Ne[iprod]->ProjectionY(Form("hNe%s_%s",prod[iprod].Data(),histoTag.Data()),minClusterEBin,1000);
    hCh[iprod] = (TH1F*)h2Ch[iprod]->ProjectionY(Form("hCh%s_%s",prod[iprod].Data(),histoTag.Data()),minClusterEBin,1000);
        
//    hNe[iprod]->Rebin(5);
//    hCh[iprod]->Rebin(5);
    
    if( !scaled )
    {
      hNe[iprod]->Sumw2();
      hCh[iprod]->Sumw2();
    }
    
    hNe[iprod]->Scale(1./nTrig);
    hCh[iprod]->Scale(1./nTrig);
    
    hNe[iprod]->SetTitle(Form("#Sigma p_{T} in cone R=0.4, p_{T,Trig}>%2.1f GeV/c",minClusterE));
    hNe[iprod]->SetYTitle("1/N_{trigger} dN/d#Sigma p_{T}");
    hNe[iprod]->SetTitleOffset(1.5,"Y");
    hNe[iprod]->SetMarkerColor(color[iprod]);
    hNe[iprod]->SetMarkerStyle(20);
    hNe[iprod]->SetAxisRange(0.,30.,"X");
    //hNe[iprod]->SetMaximum(1.1);
    //hNe[iprod]->SetMinimum(0);
    
    lprod.AddEntry(hNe[iprod],prodLeg[iprod],"P");
    
    hCh[iprod]->SetMarkerColor(color[iprod]);
    hCh[iprod]->SetMarkerStyle(25);
    
    if(iprod > 0)
    {
      hRatNe[iprod-1] = (TH1F*)hNe[iprod]->Clone(Form("hRatNe%s_%s",prod[iprod].Data(),histoTag.Data()));
      hRatCh[iprod-1] = (TH1F*)hCh[iprod]->Clone(Form("hRatCh%s_%s",prod[iprod].Data(),histoTag.Data()));
      
      hRatNe[iprod-1]->Divide(hRatNe[iprod-1],hNe[0],1.000,1,errType);
      hRatCh[iprod-1]->Divide(hRatCh[iprod-1],hCh[0],1.000,1,errType);
    }
  }
  
  /////////////////
  // Make the plots
  /////////////////
  
  // Ne
  {
    TLegend lNe(0.6,0.75,0.84,0.89);
    lNe.SetTextSize(0.04);
    //lNe.SetBorderSize(0);
    lNe.SetFillColor(0);
    lNe.AddEntry(hNe[0],"Clusters","P");
    lNe.AddEntry(hCh[0],"Tracks","P");
    
    TCanvas * cNe = new TCanvas(Form("SumPtHisto_%s_%s",histoTag.Data(),histoTag.Data()),"",1000,500);
    cNe->Divide(2,1);
    cNe->cd(1);
    gPad->SetLogy();
    hNe[0]->Draw("");
    for(Int_t iprod = 0; iprod < nProd; iprod++)
    {
      hNe[iprod]->Draw("same");
      hCh[iprod]->Draw("same");
    }
    
    lNe.Draw();
    lprod.Draw();
    
    cNe->cd(2);
    
    hRatNe[0]->SetMaximum(1.5);
    hRatNe[0]->SetMinimum(0.5);
    hRatNe[0]->Draw("");
    hRatNe[0]->SetYTitle(Form("Ratio data X / %s",prodLeg[0].Data()));
    for(Int_t iprod = 0; iprod < nProd-1; iprod++)
    {
      hRatNe  [iprod]->Draw("same");
      hRatCh[iprod]->Draw("same");
    }
    
    cNe->Print(Form("%s_ConeSumPt.%s",histoTag.Data(),format.Data()));
  }
  
}


///
/// invariant mass plots. Before InvMassFit must run.
/// 
/// To be updated
///
/// \param icalo: 0 EMCal, 1 DCal
/// \param particle: Pi0 or Eta
//______________________________________
void InvMass(Int_t icalo, TString particle, TString fileName)
{
  const Int_t nEbins = 12;
  TH1F        * hIM   [nProd][nEbins];
  TGraphErrors* gMass [nProd];
  TGraphErrors* gWidth[nProd];
  TGraphErrors* gPt   [nProd];
  
  TGraphErrors* gRatMass [nProd];
  TGraphErrors* gRatWidth[nProd];
  TGraphErrors* gRatPt   [nProd];
  
  // Open files per production and trigger, get histograms/graphs
  
  //Legend for productions
  TLegend lprod(0.6,0.7,0.84,0.89);
  lprod.SetTextSize(0.04);
  lprod.SetBorderSize(0);
  lprod.SetFillColor(0);
  
  fileName.ReplaceAll(".root","");
  
  TString calorimeter = "EMCAL";
  if(icalo == 1) calorimeter = "DCAL";
  
  for(Int_t iprod = 0; iprod <  nProd; iprod++)
  {
    TFile * filIM = TFile::Open(Form("IMfigures/%s_%s_MassWidthPtHistograms_%s_%s_AllSM.root",
                                     prod[iprod].Data(),calorimeter.Data(),particle.Data(),fileName.Data()));
//    printf("IMfigures/%s_%s_MassWidthPtHistograms_%s_%s_AllSM.root %p\n",
//           prod[iprod].Data(),calorimeter.Data(),particle.Data(),fileName.Data(),filIM);
    
    for(Int_t ie = 0; ie < nEbins; ie++)
    {
      hIM[iprod][ie] = (TH1F*) filIM->Get(Form("IM_Comb0_PtBin%d",ie));
      //    
      // hIM[iprod][ie] ->Rebin(5);
      //    
      // hIM[iprod][ie] ->Scale(1./nTrig);
      hIM[iprod][ie]->SetTitleOffset(1.5,"Y");
      hIM[iprod][ie]->SetYTitle("1/N_{events} dN/dM");
      hIM[iprod][ie]->SetMarkerColor(color[iprod]);
      hIM[iprod][ie]->SetLineColor(color[iprod]);
      // hIM[iprod][ie]->SetMarkerStyle(20);
      // hIM[iprod][ie]->SetAxisRange(0.,1.,"X");
    }
    
    gMass [iprod] = (TGraphErrors*) filIM->Get("gMass_AllSM");
    gWidth[iprod] = (TGraphErrors*) filIM->Get("gWidth_AllSM");
    gPt   [iprod] = (TGraphErrors*) filIM->Get("gPt_AllSM");
    
    gMass [iprod]->SetMarkerColor(color[iprod]);
    gWidth[iprod]->SetMarkerColor(color[iprod]);
    gPt   [iprod]->SetMarkerColor(color[iprod]);

    gMass [iprod]->SetMarkerStyle(20);
    gWidth[iprod]->SetMarkerStyle(20);
    gPt   [iprod]->SetMarkerStyle(20);
    
    lprod.AddEntry(gMass[iprod],prodLeg[iprod],"PL");
    
  } // prod
  
  // Ratio
  const Int_t     nPoints = gMass [0]->GetN();
  Double_t *      x = gMass [0]->GetX();
  Double_t *     ex = gMass [0]->GetEX();
  
  Double_t     massR[nPoints];
  Double_t    emassR[nPoints]; 
  Double_t    widthR[nPoints];
  Double_t   ewidthR[nPoints];
  Double_t       ptR[nPoints];
  Double_t      eptR[nPoints];
  
  Double_t *   massD = gMass [0]->GetY();
  Double_t *  emassD = gMass [0]->GetEY();  
  
  Double_t *  widthD = gWidth[0]->GetY();
  Double_t * ewidthD = gWidth[0]->GetEY();
  
  Double_t *     ptD = gPt   [0]->GetY();
  Double_t *    eptD = gPt   [0]->GetEY();
  
  for(Int_t iprod = 1; iprod <  nProd; iprod++)
  {
    Double_t *   massN = gMass [iprod]->GetY() ;
    Double_t *  emassN = gMass [iprod]->GetEY();  
    
    Double_t *  widthN = gWidth[iprod]->GetY();
    Double_t * ewidthN = gWidth[iprod]->GetEY();
    
    Double_t *     ptN = gPt   [iprod]->GetY();
    Double_t *    eptN = gPt   [iprod]->GetEY();
    
    for(Int_t ie = 0; ie < nPoints ; ie++)
    {
      if ( massD[ie] > 0 && massN [ie] > 0 )
      {
         massR[ie] = massN[ie]/massD[ie];
        emassR[ie] = massR[ie] * TMath::Sqrt((emassN[ie]/massN[ie])*(emassN[ie]/massN[ie]) + 
                                             (emassD[ie]/massD[ie])*(emassD[ie]/massD[ie]));
      }
      
      if ( widthD[ie] > 0 && widthN [ie] > 0 )
      {
         widthR[ie] = widthN[ie]/widthD[ie];
        ewidthR[ie] = widthR[ie] * TMath::Sqrt((ewidthN[ie]/widthN[ie])*(ewidthN[ie]/widthN[ie]) + 
                                               (ewidthD[ie]/widthD[ie])*(ewidthD[ie]/widthD[ie]));
      }
      
      if ( ptD[ie] > 0 && ptN [ie] > 0 )
      {
         ptR[ie] = ptN[ie]/ptD[ie];
        eptR[ie] = ptR[ie] * TMath::Sqrt((eptN[ie]/ptN[ie])*(eptN[ie]/ptN[ie]) + 
                                         (eptD[ie]/ptD[ie])*(eptD[ie]/ptD[ie]));
      }
      
    } // points
    
    gRatMass [iprod] = new TGraphErrors(nPoints,x,massR ,ex,emassR);
    gRatWidth[iprod] = new TGraphErrors(nPoints,x,widthR,ex,ewidthR);
    gRatPt   [iprod] = new TGraphErrors(nPoints,x,ptR   ,ex,eptR);
    
    gRatMass [iprod]->SetMarkerColor(color[iprod]);
    gRatMass [iprod]->SetLineColor(color[iprod]);    
    gRatWidth[iprod]->SetMarkerColor(color[iprod]);
    gRatWidth[iprod]->SetLineColor(color[iprod]);    
    gRatPt   [iprod]->SetMarkerColor(color[iprod]);
    gRatPt   [iprod]->SetLineColor(color[iprod]);
    
  } // prod
  
  /////////////////
  // Make the plots
  /////////////////
  
  // Mass
  {
    TCanvas * cMass = new TCanvas(Form("Mass_%s_%s",particle.Data(), histoTag.Data()),"",1000,500);
    cMass->Divide(2,1);
    
    cMass->cd(1);
    //gPad->SetLogy();
    
    gMass[1]->SetMinimum(120);
    gMass[1]->SetMaximum(170);
    gMass[0]->Draw("AP");
    for(Int_t iprod = 0; iprod <  nProd; iprod++)
    {
      gMass[iprod]->Draw("P");
    }
    
    lprod.Draw();
    
    cMass->cd(2);
    gPad->SetGridy();
    
    gRatMass[1]->SetTitle("Mass Ratio");
    gRatMass[1]->GetHistogram()->SetXTitle("p_{T} (GeV/c)");
    gRatMass[1]->GetHistogram()->SetYTitle(Form("Ratio data X / %s",prodLeg[0].Data()));
    gRatMass[1]->GetHistogram()->SetTitleOffset(1.4,"Y");
    gRatMass[1]->SetMinimum(0.97);
    gRatMass[1]->SetMaximum(1.05);
    gRatMass[1]->Draw("AP");
    
    for(Int_t iprod = 1; iprod <  nProd; iprod++)
    {
      gRatMass [iprod]->Draw("P");
    }
    
      
    cMass->Print(Form("%s_Mass_%s.%s",histoTag.Data(),particle.Data(),format.Data()));
  }
  // Width
  {
    TCanvas * cWidth = new TCanvas(Form("Width_%s_%s",particle.Data(), histoTag.Data()),"",1000,500);
    cWidth->Divide(2,1);
    
    cWidth->cd(1);
    //gPad->SetLogy();
    
    gWidth[0]->SetMinimum(8);
    gWidth[0]->SetMaximum(20);
    gWidth[0]->Draw("AP");
    for(Int_t iprod = 0; iprod <  nProd; iprod++)
    {
      gWidth[iprod]->Draw("P");
    }
    
    lprod.Draw();
    
    cWidth->cd(2);
    gPad->SetGridy();
    
    gRatWidth[1]->SetTitle("Width Ratio");
    gRatWidth[1]->GetHistogram()->SetXTitle("p_{T} (GeV/c)");
    gRatWidth[1]->GetHistogram()->SetYTitle(Form("Ratio data X / %s",prodLeg[0].Data()));
    gRatWidth[1]->GetHistogram()->SetTitleOffset(1.4,"Y");
    gRatWidth[1]->SetMinimum(0.9);
    gRatWidth[1]->SetMaximum(1.1);
    gRatWidth[1]->Draw("AP");
    
    for(Int_t iprod = 1; iprod <  nProd; iprod++)
    {
      gRatWidth [iprod]->Draw("P");
    }
    
    
    cWidth->Print(Form("%s_Width_%s.%s",histoTag.Data(),particle.Data(),format.Data()));
  }
  
  
  // Pt
  {
    TCanvas * cPt = new TCanvas(Form("Pt_%s_%s",particle.Data(), histoTag.Data()),"",1000,500);
    cPt->Divide(2,1);
    
    cPt->cd(1);
    gPad->SetLogy();
    
    gPt[0]->Draw("AP");
    for(Int_t iprod = 0; iprod <  nProd; iprod++)
    {
      gPt[iprod]->Draw("P");
    }
    
    lprod.Draw();
    
    cPt->cd(2);
    gPad->SetGridy();
    
    gRatPt[1]->SetTitle("Pt Ratio");
    gRatPt[1]->GetHistogram()->SetXTitle("p_{T} (GeV/c)");
    gRatPt[1]->GetHistogram()->SetYTitle(Form("Ratio data X / %s",prodLeg[0].Data()));
    gRatPt[1]->GetHistogram()->SetTitleOffset(1.4,"Y");
    gRatPt[1]->SetMinimum(0.7);
    gRatPt[1]->SetMaximum(1.3);
    gRatPt[1]->Draw("AP");
    
    for(Int_t iprod = 1; iprod <  nProd; iprod++)
    {
      gRatPt [iprod]->Draw("P");
    }
    
    
    cPt->Print(Form("%s_Pt_%s.%s",histoTag.Data(),particle.Data(),format.Data()));
  }
  
  // Invariant Mass
  {
    TGaxis::SetMaxDigits(3);
    
    TCanvas * cIM = new TCanvas(Form("IM_%s_%s",particle.Data(),histoTag.Data()),"",1000,1000);
    cIM->Divide(4,3);
    
    for(Int_t ie = 0; ie < nEbins; ie++)
    {
      cIM->cd(ie+1);
      //gPad->SetLogy();
      
      hIM[0][ie]->Draw("");
      for(Int_t iprod = 0; iprod < nProd; iprod++)
      {
        hIM[iprod][ie]->Draw("same");
      }
      
      lprod.Draw();
    }
    
    cIM->Print(Form("%s_IM_%s.%s",histoTag.Data(),particle.Data(),format.Data()));
  }
  
}

/// Centrality
/// To be updated
//______________________________________
void Centrality()    
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
    
    if( !scaled )
    {
      hCen[iprod]->Sumw2();
      hCen[iprod]->Scale(1./nEvents[iprod]);
    }
    
    hCen[iprod]->SetTitle("Centrality");
    hCen[iprod]->SetYTitle("1/N_{events} dN/d centrality");
    hCen[iprod]->SetTitleOffset(1.5,"Y");
    hCen[iprod]->SetMarkerColor(color[iprod]);
    hCen[iprod]->SetMarkerStyle(20);
    //hCen[iprod]->SetAxisRange(0.,30.,"X");
    //hCen[iprod]->SetMaximum(1.1);
    //hCen[iprod]->SetMinimum(0);
    lprod.AddEntry(hCen[iprod],prodLeg[iprod],"P");

    if(iprod > 0)
    {
      hRatCen[iprod-1] = (TH1F*)hCen[iprod]->Clone(Form("hRatCen%s" ,prod[iprod].Data()));
      
      hRatCen[iprod-1]->Divide(hRatCen[iprod-1],hCen [0],1.000,1,errType);
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
void Vertex()
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
    
    lprod.AddEntry(hVertex[0][iprod],prodLeg[iprod],"P");
    
    for(Int_t ivertex = 0; ivertex < 3; ivertex++)
    {
      if( !scaled )
      {
        hVertex[ivertex][iprod]->Sumw2();
        hVertex[ivertex][iprod]->Scale(1./nEvents[iprod]);
      }
      
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
        
        hRatVertex[ivertex][iprod-1]->Divide(hRatVertex[ivertex][iprod-1],hVertex[ivertex][0],1.000,1,errType);
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
    histArr[iprod] = 0;
    nEvents[iprod] = 0;
    
    // now get them
    file[iprod]  = new TFile(Form("%s/%s",prod[iprod].Data(),fileName.Data()),"read");
    if(file[iprod]->Get("hNEvents"))
    {
      nEvents[iprod] = ((TH1F*)file[iprod]->Get("hNEvents"))->GetEntries();
      printf("%s: nEvents %e\n",prod[iprod].Data(),nEvents[iprod]);
      //printf("iprod %d, File %p\n",iprod,file[iprod]);
      continue;
    }
    
    dir[iprod] = (TDirectoryFile*) file[iprod]->Get(listName);
    if(dir[iprod])
    {
      histArr[iprod] = (TList*) dir[iprod]->Get(trigName);
      
      if ( !histArr[iprod]                    ) continue;
      if (  histArr[iprod]->GetEntries() <= 0 ) continue;
      
      nEvents[iprod] = ((TH1F*)histArr[iprod]->FindObject("hNEvents"))->GetEntries();
      printf("%s: nEvents %e\n",prod[iprod].Data(),nEvents[iprod]);
    }
    else
    {
      histArr[iprod] = (TList*) file[iprod]->Get(trigName);
      
      if ( !histArr[iprod]                    ) continue;
      if (  histArr[iprod]->GetEntries() <= 0 ) continue;
      
      nEvents[iprod] = ((TH1F*)histArr[iprod]->FindObject("hNEvents"))->GetEntries();
      printf("%s: nEvents %e\n",prod[iprod].Data(),nEvents[iprod]);
    }
    
    if(!histArr[iprod])
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
  if ( histArr[iprod] ) 
    return histArr[iprod]->FindObject(histoName);
  else  if ( file[iprod] ) 
    return file[iprod]->Get(histoName);
  
  printf("Null file %p or list %p for prod %d and histo %s\n",
         file[iprod],histArr[iprod],iprod, histoName.Data());
  
  return 0x0;
}


