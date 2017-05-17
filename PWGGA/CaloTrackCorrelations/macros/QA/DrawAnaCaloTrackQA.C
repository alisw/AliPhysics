/// \file DrawAnaCaloTrackQA.C
/// \ingroup CaloTrackCorrMacrosQA
/// \brief Plot analysis QA histograms from EMCal PWG-GA wagon
///
/// Macro to plot few selected histograms
/// to QA data productions at 0th order
/// Analysis performed with the wagon
/// AddTaskPi0IMGammaCorrQA.C
/// It generates 7 plots, each containing 2 to 4 pads
///
/// To execute: root -q -b -l DrawAnaCaloTrackQA.C'("Pi0IM_GammaTrackCorr_EMCAL","AnalysisResults.root")'
///
/// The trigger name might change depending on the wagon / data type
/// In simulations only the "default" case is available
/// On data, there can be different triggers, depending on the period
/// * default: min bias like triggers, kINT7, kINT1, kMB
/// * EMCAL_L0: kEMC7 L0 EMCal
/// * DCAL_L0: kEMC7 L0 DCal
/// * EMCAL_L1: kEMCEGA L1 EG1 EMCal
/// * DCAL_L1: kEMCEGA L1 EG1 DCal 
/// * EMCAL_L2: kEMCEGA L1 EG2 EMCal
/// * DCAL_L2: kEMCEGA L1 EG2 DCal
///
/// In case output file is too large, possiblity to dump the list content in a sepate file:  exportToFile = kTRUE
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

void ProcessTrigger(TString trigName    = "default", 
                    Bool_t exportToFile = kFALSE, 
                    Bool_t checkList    = kTRUE);
void CaloQA    (Int_t icalo);
void TrackQA   ();
void Pi0QA     (Int_t icalo);
void IsolQA    (Int_t icalo);
void CorrelQA  (Int_t icalo);
void MCQA      (Int_t icalo);
void ScaleAxis (TAxis *a, Double_t scale);
void ScaleXaxis(TH1   *h, Double_t scale);

TObject * GetHisto(TString histoName);
Bool_t    GetList (TString trigName, Bool_t  exportToFile);
//
//---------------------------------------------------------

//-----------------------
// Some global variables
TDirectoryFile *dir = 0;
TList  *list = 0;
TFile  *file = 0;
TString histoTag = "";
TString format = "eps";
Int_t   color[]={kBlack,kRed,kOrange+1,kYellow+1,kGreen+2,kBlue,kCyan+1,kViolet,kMagenta+2,kGray,kCyan-2,kViolet-2};
//
//-----------------------

///
/// Main method, produce the plots for the 7 different types of triggers:
///
/// * Calorimeter QA: cluster/cell spectra, acceptance, track matching residuals in CaloQA method
/// * Track QA in TrackQA method
/// * Invariant mass plots in Pi0QA method
/// * Cluster isolation, tracks and clusters pT and sum pT in cone, in IsolQA method
/// * Cluster-track correlation plots, azimuthal and xE, in CorrelQA method
/// * Dedicated generated particles QA in MCQA method
///
/// Input:
/// \param listName: Name of list with histograms in file
/// \param fileName: File name
/// \param exportToFile: export list with histograms to separate file, intereting in case of big output file.
/// \param fileFormat: define the type of figures: eps, pdf, etc.
//_______________________________________________________________________
void DrawAnaCaloTrackQA
(
 TString listName     = "Pi0IM_GammaTrackCorr_EMCAL",
 TString fileName     = "AnalysisResults.root",
 Bool_t  exportToFile = kFALSE,
 TString fileFormat   = "eps"
)
{
  format = fileFormat;
  printf("Open <%s>; Get Trigger List : <%s>; Export list? <%d>; format %s\n",
         fileName.Data(),listName.Data(),exportToFile, format.Data());
 
  file  = new TFile(fileName,"read");
  if ( !file ) 
  { 
    printf("File not found, do nothing\n");
    return; 
  }
  
  dir = (TDirectoryFile*) file->Get(listName);
  if ( !dir ) 
  { 
    printf("DirectoryFile not found, do nothing\n");
    return; 
  }
  
  ProcessTrigger("default" ,exportToFile);
  ProcessTrigger("EMCAL_L0",exportToFile);
  ProcessTrigger("EMCAL_L1",exportToFile);
  ProcessTrigger("EMCAL_L2",exportToFile);
  ProcessTrigger("DCAL_L0" ,exportToFile);
  ProcessTrigger("DCAL_L1" ,exportToFile);
  ProcessTrigger("DCAL_L2" ,exportToFile);
  
  file->Close();
}

/// 
/// Produce the plots per trigger, options are:
/// * EMCAL_L0: kEMC7 L0 EMCal
/// * DCAL_L0: kEMC7 L0 DCal
/// * EMCAL_L1: kEMCEGA L1 EG1 EMCal
/// * DCAL_L1: kEMCEGA L1 EG1 DCal 
/// * EMCAL_L2: kEMCEGA L1 EG2 EMCal
/// * DCAL_L2: kEMCEGA L1 EG2 DCal
///
/// Input:
/// \param trigName: File name
/// \param checklist: get the list from file, in case not exported
///
/// This method can be executed directly instead od DrawAnaCaloTrackQA if the
/// list with histograms were exported previously into a separate file 
/// and checkList is set to false.
//_______________________________________________________________________
void ProcessTrigger( TString trigName, Bool_t exportToFile, Bool_t checkList)
{
  // Access the list of histograms, global variables
  if(checkList)
  {
    Bool_t ok = GetList(trigName, exportToFile);
    printf("\t -- Process trigger %s, ok %d\n",trigName.Data(), ok);
    
    if ( !ok ) return;
  }

  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(000000);
  gStyle->SetPadRightMargin(0.15);
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
  
  for(Int_t icalo = calo; icalo < nCalo; icalo++)
  {
    if(trigName.Contains("default")) histoTag=Form("%s_%s",caloString[icalo].Data(),trigName.Data());
    
    // Plot basic QA
    CaloQA(icalo);
    
    // Plot basic Pi0 QA
    Pi0QA(icalo);
    
    // Plot basic isolation QA
    IsolQA(icalo);
    
    // Plot basic correlation QA
    CorrelQA(icalo);
    
    // MC basic QA plots, cluster origins (only if it run on MC)
    MCQA(icalo);
  }
}

///
/// Plot basic calorimeter QA histograms.
/// 
/// \param icalo: 0 EMCal, 1 DCal
///
//______________________________________
void CaloQA(Int_t icalo)
{ 
  TCanvas * ccalo = new TCanvas(Form("CaloHisto_%s",histoTag.Data()),
                                "Cluster spectra and track match residuals",1000,1000);
  ccalo->Divide(2,2);
  
  ccalo->cd(1);
  gPad->SetLogy();
  gPad->SetLogx();
      
  TH1F* hClusterEnergy = (TH1F*) GetHisto(Form("AnaPhoton_Calo%d_hPt_Cut_6_Fidutial",icalo));
  if(!hClusterEnergy) return;
  hClusterEnergy->SetYTitle("entries");
  hClusterEnergy->SetTitle("Cluster-cell energy spectra");
  hClusterEnergy->SetTitleOffset(1.5,"Y");
  hClusterEnergy->Sumw2();
  hClusterEnergy->SetMarkerColor(1);
  hClusterEnergy->SetMarkerStyle(20);
  hClusterEnergy->SetAxisRange(0.,50.,"X");
  hClusterEnergy->Draw();
  
  TH2F* h2CellAmplitude = (TH2F*) GetHisto("QA_Cell_hAmp_Mod");
  TH1F* hCellAmplitude  = 0;
  if(histoTag.Contains("default"))
  {
    if ( icalo == 0 ) hCellAmplitude = (TH1F*) h2CellAmplitude->ProjectionX(Form("CellAmp%s",histoTag.Data()), 1,12);
    else              hCellAmplitude = (TH1F*) h2CellAmplitude->ProjectionX(Form("CellAmp%s",histoTag.Data()),12,20);
  }
  else                hCellAmplitude = (TH1F*) h2CellAmplitude->ProjectionX(Form("CellAmp%s",histoTag.Data()),0,100);
      
  hCellAmplitude->Sumw2();
  hCellAmplitude->SetMarkerColor(4);
  hCellAmplitude->SetMarkerStyle(25);
  hCellAmplitude->Draw("same");
  
  TLegend l(0.4,0.79,0.83,0.89);
  l.SetTextSize(0.04);
  l.AddEntry(hClusterEnergy,"Good Cluster","P");
  l.AddEntry(hCellAmplitude,"Cell","P");
  l.SetBorderSize(0);
  l.SetFillColor(0);
  l.Draw();
  
  ccalo->cd(2);
  //gPad->SetLogy();
  gPad->SetLogx();

  TH1F* hRaw  = (TH1F*) GetHisto(Form("AnaPhoton_Calo%d_hPt_Cut_6_Fidutial",icalo));
  //TH1F* hCorr = (TH1F*) GetHisto(Form("AnaPhoton_Calo%d_hPt_Cut_4_NCells"  ,icalo));
  TH1F* hTM   = (TH1F*) GetHisto(Form("AnaPhoton_Calo%d_hPt_Cut_7_Matching",icalo));
  TH1F* hShSh = (TH1F*) GetHisto(Form("AnaPhoton_Calo%d_hPt_Cut_9_PID"     ,icalo));
  
  //hRaw->Sumw2();
  
//  hCorr->SetTitle("Ratio after cluster cuts application");
//  hCorr->SetYTitle("Selected clusters / Raw clusters");
//  hCorr->SetTitleOffset(1.5,"Y");
//  hCorr->Sumw2();
//  hCorr->SetMarkerColor(1);
//  hCorr->SetMarkerStyle(20);
//  hCorr->Divide(hRaw);
//  hCorr->SetAxisRange(0.,30.,"X");
//  hCorr->SetMaximum(1.1);
//  hCorr->SetMinimum(0);
//  hCorr->Draw();
  
  hTM  ->SetTitle("Ratio after cluster cuts application");
  hTM  ->SetYTitle("Selected clusters / Good clusters");
  hTM  ->SetTitleOffset(1.5,"Y");
  hTM  ->Sumw2();
  hTM  ->SetAxisRange(0.,50.,"X");
  hTM  ->SetMarkerColor(2);
  hTM  ->SetMarkerStyle(21);
  hTM  ->SetMaximum(1.1);
  hTM  ->SetMinimum(0);
  hTM  ->Divide(hRaw);
  hTM  ->Draw("");
  
  hShSh->Sumw2();
  hShSh->SetMarkerColor(4);
  hShSh->SetMarkerStyle(22);
  hShSh->Divide(hRaw);
  hShSh->Draw("same");
  
  TLegend l2(0.45,0.8,0.95,0.93);
  l2.SetTextSize(0.04);
  //l2.AddEntry(hCorr,"No Exotics + non lin.","P");
  l2.AddEntry(hTM,  "+ Track matching","P");
  l2.AddEntry(hShSh,"+ #lambda^{2}_{0} < 0.4","P");
  l2.SetBorderSize(0);
  l2.SetFillColor(0);
  l2.Draw();
  
  
  // Plot track-matching residuals
  TH1F* hTrackMatchResEtaNeg;
  TH1F* hTrackMatchResEtaPos;
  TH1F* hTrackMatchResPhiNeg;
  TH1F* hTrackMatchResPhiPos;

  // first test did not have this histogram, add protection
  TH2F* hTrackMatchResEtaPhi = (TH2F*) GetHisto(Form("AnaPhoton_Calo%d_hTrackMatchedDEtaDPhiPosNoCut",icalo));
  if(hTrackMatchResEtaPhi)
  {
    hTrackMatchResEtaPhi ->Add( (TH2F*) GetHisto(Form("AnaPhoton_Calo%d_hTrackMatchedDEtaDPhiNegNoCut",icalo) ));
    
    ccalo->cd(3);
    gPad->SetLogz();
    
    hTrackMatchResEtaPhi->SetAxisRange(-0.025,0.025,"X");
    hTrackMatchResEtaPhi->SetAxisRange(-0.025,0.025,"Y");
    hTrackMatchResEtaPhi->SetTitleOffset(1.5,"Y");
    hTrackMatchResEtaPhi->SetTitle("Track-cluster residual #Delta#varphi vs #Delta#eta, #it{E}>0.5 GeV");
    hTrackMatchResEtaPhi->SetXTitle("#Delta #eta");
    hTrackMatchResEtaPhi->SetYTitle("#Delta #varphi");
    hTrackMatchResEtaPhi->SetZTitle("entries");
    hTrackMatchResEtaPhi->Draw("colz");
    
    ccalo->cd(4);
    gPad->SetLogy();
    
    TH2F* h2TrackMatchResEtaNeg = (TH2F*) GetHisto(Form("AnaPhoton_Calo%d_hTrackMatchedDEtaNeg",icalo));
    TH2F* h2TrackMatchResEtaPos = (TH2F*) GetHisto(Form("AnaPhoton_Calo%d_hTrackMatchedDEtaPos",icalo));
    TH2F* h2TrackMatchResPhiNeg = (TH2F*) GetHisto(Form("AnaPhoton_Calo%d_hTrackMatchedDPhiNeg",icalo));
    TH2F* h2TrackMatchResPhiPos = (TH2F*) GetHisto(Form("AnaPhoton_Calo%d_hTrackMatchedDPhiPos",icalo));
    
    Float_t binMin = hRaw->FindBin(0.5);
    hTrackMatchResEtaNeg = (TH1F*) h2TrackMatchResEtaNeg->ProjectionY(Form("TMProjEtaNeg%s",histoTag.Data()),binMin, 1000);
    hTrackMatchResEtaPos = (TH1F*) h2TrackMatchResEtaPos->ProjectionY(Form("TMProjEtaPos%s",histoTag.Data()),binMin, 1000);
    hTrackMatchResPhiNeg = (TH1F*) h2TrackMatchResPhiNeg->ProjectionY(Form("TMProjPhiNeg%s",histoTag.Data()),binMin, 1000);
    hTrackMatchResPhiPos = (TH1F*) h2TrackMatchResPhiPos->ProjectionY(Form("TMProjPhiPos%s",histoTag.Data()),binMin, 1000);
        
    hTrackMatchResEtaNeg->SetXTitle("#Delta #eta, #Delta #varphi");
    hTrackMatchResEtaNeg->SetYTitle("entries");
    hTrackMatchResEtaNeg->SetTitle("Track-cluster residuals, #it{E} > 1 GeV");
    hTrackMatchResEtaNeg->SetAxisRange(-0.05,0.05,"X");
    hTrackMatchResEtaNeg->Sumw2();
    hTrackMatchResEtaNeg->SetMarkerStyle(25);
    hTrackMatchResEtaNeg->SetMarkerColor(2);
    hTrackMatchResEtaNeg->Draw("");
    
    hTrackMatchResEtaPos->Sumw2();
    hTrackMatchResEtaPos->SetMarkerStyle(25);
    hTrackMatchResEtaPos->SetMarkerColor(4);
    hTrackMatchResEtaPos->Draw("same");
    
    hTrackMatchResPhiNeg->Sumw2();
    hTrackMatchResPhiNeg->SetMarkerStyle(24);
    hTrackMatchResPhiNeg->SetMarkerColor(2);
    hTrackMatchResPhiNeg->Draw("same");
    
    hTrackMatchResPhiPos->Sumw2();
    hTrackMatchResPhiPos->SetMarkerStyle(24);
    hTrackMatchResPhiPos->SetMarkerColor(4);
    hTrackMatchResPhiPos->Draw("same");
    
    TLine l0(0,hTrackMatchResEtaNeg->GetMinimum(),0,hTrackMatchResEtaNeg->GetMaximum()*1.2);
    l0.Draw("same");
    
    TLegend l3(0.55,0.7,0.83,0.85);
    l3.SetTextSize(0.04);
    l3.AddEntry(hTrackMatchResEtaNeg,"#Delta #eta, Negative","P");
    l3.AddEntry(hTrackMatchResEtaPos,"#Delta #eta, Positive","P");
    l3.AddEntry(hTrackMatchResPhiNeg,"#Delta #varphi, Negative","P");
    l3.AddEntry(hTrackMatchResPhiPos,"#Delta #varphi, Positive","P");
    l3.SetBorderSize(0);
    l3.SetFillColor(0);
    l3.Draw();
  }
  
  ccalo->Print(Form("%s_CaloHisto.%s",histoTag.Data(),format.Data()));
  
  // cleanup
  //
  delete hCellAmplitude;
  delete hTrackMatchResEtaNeg;
  delete hTrackMatchResEtaPos;
  delete hTrackMatchResPhiNeg;
  delete hTrackMatchResPhiPos;

  delete ccalo;
  
  //
  //
  //
  
  TCanvas * ccalo2 = new TCanvas(Form("CaloHisto2_%s",histoTag.Data()),"",500,500);
  ccalo2->Divide(2,2);
  
  ccalo2->cd(3);
//  gPad->SetLogz();
//  TH2F* hCellAmpId   = (TH2F*) GetHisto("QA_hAmpId");
//  hCellAmpId->SetTitle("Cell Id vs energy");
//  hCellAmpId->SetYTitle("Cell Id");
//  //hCellAmpId->SetAxisRange(300.,900.,"Y");
//  hCellAmpId->SetAxisRange(0.,30.,"X");
//  hCellAmpId->SetTitleOffset(1.5,"Y");
//  hCellAmpId->Draw("colz");

  gPad->SetLogz();
  
  TH2F* hClusterTime   = (TH2F*) GetHisto(Form("AnaPhoton_Calo%d_hTimePt",icalo));
  hClusterTime->SetTitle("Cluster energy vs time");
  hClusterTime->SetYTitle("time (ns)");
  //hClusterTime->SetAxisRange(300.,900.,"Y");
  hClusterTime->SetAxisRange(0.,30.,"X");
  hClusterTime->SetTitleOffset(1.5,"Y");
  hClusterTime->Draw("colz");
  
  ccalo2->cd(1);
  gPad->SetLogz();

  TH2F* hCellActivity  = (TH2F*) GetHisto("QA_Cell_hGridCells");
  if(icalo == 0)  hCellActivity->SetAxisRange(  0,127,"Y");
  else            hCellActivity->SetAxisRange(128,220,"Y");
  hCellActivity->SetTitle("Hits per cell (#it{E} > 0.2 GeV)");
  hCellActivity->SetTitleOffset(1.5,"Y");
  hCellActivity->Draw("colz");
  
  ccalo2->cd(2);
  
  TH2F* hCellActivityE = (TH2F*) GetHisto("QA_Cell_hGridCellsE");
  
  if(icalo == 0)  hCellActivityE->SetAxisRange(  0,127,"Y");
  else            hCellActivityE->SetAxisRange(128,220,"Y");
  
  hCellActivityE->SetTitle("Mean energy per cell (#it{E} > 0.2 GeV)");
  
  if(icalo != 1 && !histoTag.Contains("default")) // ratio already done for calo=0
    hCellActivityE->Divide(hCellActivity);
  
  hCellActivityE->SetTitleOffset(1.5,"Y");
  hCellActivityE->Draw("colz");
  
  ccalo2->cd(4);
  gPad->SetLogz();
  
  TH2F* hClusterActivity  = (TH2F*) GetHisto(Form("AnaPhoton_Calo%d_hEBin1_Cluster_ColRow_PID",icalo));
  if(histoTag.Contains("default"))
        hClusterActivity  = (TH2F*) GetHisto(Form("AnaPhoton_Calo%d_hEBin0_Cluster_ColRow_PID",icalo));
  
  if(histoTag.Contains("default")) hClusterActivity->SetTitle("Clusters per col-row 0.5<#it{E}<3 GeV");
  else if(histoTag.Contains("L0")) hClusterActivity->SetTitle("Clusters per col-row #it{E} > 5 GeV");
  else                             hClusterActivity->SetTitle("Clusters per col-row #it{E} > 12 GeV");
  
  hClusterActivity->SetTitleOffset(1.5,"Y");
  hClusterActivity->Draw("colz");
  
  ccalo2->Print(Form("%s_CaloHisto2.%s",histoTag.Data(),format.Data()));

  // cleanup
  //
  delete ccalo2;
}

///
/// Plot basic hybrid tracks histograms in 4 pads:
/// * tracks eta vs phi
/// * track phi distribution per hybrid track component
/// * track TOF
/// * track pT per hybrid track component
///
//______________________________________
void TrackQA()
{
  TCanvas * ctrack = new TCanvas(Form("TrackHisto_%s",histoTag.Data()),
                                 "Hybrid tracks",500,500);
  ctrack->Divide(2,2);
  
  ctrack->cd(1);
  //gPad->SetLogz();
  TH2F * hTrackEtaPhi = (TH2F*) GetHisto("AnaHadrons_hEtaPhiNegative");
  hTrackEtaPhi ->Add((TH2F*) GetHisto("AnaHadrons_hEtaPhiNegative"));
  hTrackEtaPhi ->SetAxisRange(-0.9,0.9,"X");
  hTrackEtaPhi ->SetTitleOffset(1.5,"Y");
  hTrackEtaPhi ->SetTitle("Hybrid tracks #eta vs #varphi #it{p}_{T} > 0.2 GeV/#it{c}");
  hTrackEtaPhi ->Draw("colz");
  
  ctrack->cd(2);
  //gPad->SetLogy();
  TH2F * hTrackEtaPhiSPD   = (TH2F*) GetHisto("AnaHadrons_hEtaPhiSPDRefitPt02");
  TH2F * hTrackEtaPhiNoSPD = (TH2F*) GetHisto("AnaHadrons_hEtaPhiNoSPDRefitPt02");
  
  TH1F* hPhiSPD   = (TH1F*)hTrackEtaPhiSPD  ->ProjectionY("hTrackPhiSPD"  ,0,1000);
  TH1F* hPhiNoSPD = (TH1F*)hTrackEtaPhiNoSPD->ProjectionY("hTrackPhiNoSPD",0,1000);
  //TH1F* hPhi      = (TH1F*)hTrackEtaPhi     ->ProjectionY("hTrackPhi"     ,0,1000);
  TH1F* hPhi      = (TH1F*)hPhiSPD->Clone("hTrackPhi");
  hPhi->Add(hPhiNoSPD);
  hPhi     ->SetTitle("Hybrid track type #varphi, 0.2<#it{p}_{T}<2 GeV/#it{c}");
  hPhi     ->SetLineColor(1);
  hPhiSPD  ->SetLineColor(2);
  hPhiNoSPD->SetLineColor(4);
  
  hPhi     ->SetMinimum(1);
  hPhi     ->SetMaximum(hPhi->GetMaximum()*1.3);
  hPhi     ->SetTitleOffset(1.5,"Y");
  
  hPhi     ->Draw("H");
  hPhiSPD  ->Draw("Hsame");
  hPhiNoSPD->Draw("Hsame");
  
  TLegend l(0.2,0.75,0.4,0.89);
  l.SetTextSize(0.04);
  l.AddEntry(hPhi,"Sum","L");
  l.AddEntry(hPhiSPD  ,"SPD+Refit","L");
  l.AddEntry(hPhiNoSPD,"No SPD+Refit","L");
  l.SetBorderSize(0);
  l.SetFillColor(0);
  l.Draw();
  
  ctrack->cd(3);
  gPad->SetLogy();
  
  TH1F* hTOF = (TH1F*) GetHisto("AnaHadrons_hTOFSignalPtCut");
  hTOF->Draw("colz");
  hTOF->SetTitleOffset(1.5,"Y");

  ctrack->cd(4);
  gPad->SetLogy();
  gPad->SetLogx();
  
  TH1F* hPt      = (TH1F*) GetHisto("AnaHadrons_hPt");
  TH1F* hPtSPD   = (TH1F*) GetHisto("AnaHadrons_hPtSPDRefit");
  TH1F* hPtNoSPD = (TH1F*) GetHisto("AnaHadrons_hPtNoSPDRefit");
  hPt     ->SetLineColor(1);
  hPtSPD  ->SetLineColor(2);
  hPtNoSPD->SetLineColor(4);
  
  hPt     ->SetTitle("Hybrid track type #it{p}_{T}");
  hPt     ->SetTitleOffset(1.5,"Y");

  hPt     ->Draw("");
  hPtSPD  ->Draw("same");
  hPtNoSPD->Draw("same");
  
//  ctrack->cd(3);
//  gPad->SetLogz();
//  
//  TH2F* hPtDCAxy = (TH2F*) GetHisto("AnaHadrons_hPtDCAxy");
//  hPtDCAxy->SetAxisRange(-1,1,"Y");
//  hPtDCAxy->SetAxisRange(0,30,"X");
//  hPtDCAxy->Draw("colz");
//  
//  ctrack->cd(4);
//  gPad->SetLogz();
//  
//  TH2F* hPtDCAz = (TH2F*) GetHisto("AnaHadrons_hPtDCAz");
//  hPtDCAz->SetAxisRange(-1,1,"Y");
//  hPtDCAz->SetAxisRange(0,30,"X");
//  hPtDCAz->Draw("colz");
  
  ctrack->Print(Form("%s_TrackHisto.%s",histoTag.Data(),format.Data()));
  
  // cleanup
  //
  delete hPhi;
  delete hPhiSPD;
  delete hPhiNoSPD;
  
  delete ctrack;
}

///
/// Plot basic invariant mass QA in 4 pads:
/// * Invariant mass vs pT
/// * Invariant mass real/mixed pairs, in pi0 region
/// * Invariant mass real/mixed pairs, in pi0 region per super module
/// * Invariant mass real/mixed pairs, in eta region
/// 
/// \param icalo: 0 EMCal, 1 DCal
///
//_____________________________
void Pi0QA(Int_t icalo)
{
  TCanvas * cpi0 = new TCanvas(Form("Pi0Histo_%s",histoTag.Data()),
                               "Neutral mesons",500,500);
  cpi0->Divide(2,2);
  
  TH2F* hMassE[10];
  TH2F* hMixMassE[10];
  for(Int_t icen = 0; icen < 10; icen++)
  {
    hMassE   [icen] = (TH2F*) GetHisto(Form("AnaPi0_Calo%d_hRe_cen%d_pidbit0_asy0_dist1",icalo,icen));
    hMixMassE[icen] = (TH2F*) GetHisto(Form("AnaPi0_Calo%d_hMi_cen%d_pidbit0_asy0_dist1",icalo,icen));
  }
  if(!hMassE[0]) return;

  // 2D Invariant mass vs E, in PbPb from 60 to 100 %, all in pp
  cpi0->cd(1);
  gPad->SetLogz();
  TH2F* h2DMass;
  
  if(hMassE[1]) // Plot centrality from 60 to 100%
  {
    h2DMass = (TH2F*) hMassE[6]->Clone("h2DMass");
    for(Int_t icen = 7; icen < 10; icen++) h2DMass->Add(hMassE[icen]);
    h2DMass->SetTitle("Inv. mass vs #it{p}_{T,pair}, Cen: 60-100%");
  }
  else
  {
    h2DMass = (TH2F*) hMassE[0]->Clone("hMassProj");
    h2DMass->SetTitle("Inv. mass vs #it{p}_{T,pair}");
  }
  
  h2DMass->SetTitleOffset(1.6,"Y");
  h2DMass->SetAxisRange(0.0,0.7,"Y");
  h2DMass->SetAxisRange(0,30,"X");
  h2DMass->Draw("colz");
  
  // Pi0 Invariant mass projection, in PbPb 6 centrality bins from 0 to 50%, all in pp
  cpi0->cd(2);
  TH1F* hMass   [10];
  TH1F* hMix    [10];
  TH1F* hMassEta[10];
  TH1F* hMassPi0[10];
  
  //Init to 0
  for(Int_t icen=0; icen<10; icen++ )
  {
    hMass   [icen] = 0;
    hMix    [icen] = 0;
    hMassEta[icen] = 0;
    hMassPi0[icen] = 0;
  }
    
  TH1F * hX = (TH1F*) hMassE[0]->ProjectionX("hEPairCen0",0,10000);
  Int_t binmin = hX->FindBin(2);  // Project histo from 2 GeV pairs
  Int_t binmax = hX->FindBin(10); // Project histo up to 10 GeV pairs
  if(histoTag.Contains("L0"))
  {
    binmin = hX->FindBin(5);  // Project histo from 5 GeV pairs
    binmax = hX->FindBin(10); // Project histo up to 10 GeV pairs
  }
  else if(histoTag.Contains("L2"))
  {
    binmin = hX->FindBin(8);  // Project histo from 8 GeV pairs
    binmax = hX->FindBin(12); // Project histo up to 12 GeV pairs
  }
  else if(histoTag.Contains("L1"))
  {
    binmin = hX->FindBin(10);  // Project histo from 10 GeV pairs
    binmax = hX->FindBin(15); // Project histo up to 15 GeV pairs
  }
  
  Float_t maxPi0 = 0;
  Float_t maxEta = 0;
  for(Int_t icen = 0; icen < 6; icen++)
  {
    if(!hMassE[icen]) continue;

    hMass[icen] = (TH1F*) hMassE   [icen]->ProjectionY(Form("hMassCen%d%s",icen,histoTag.Data()),binmin,binmax);
    hMix [icen] = (TH1F*) hMixMassE[icen]->ProjectionY(Form("hMixCen%d%s" ,icen,histoTag.Data()),binmin,binmax);
    hMass[icen]->Sumw2();
    hMix [icen]->Sumw2();
    
    hMassPi0[icen] = (TH1F*) hMass[icen]->Clone(Form("hMassPi0Cen%d%s",icen,histoTag.Data()));
    hMassEta[icen] = (TH1F*) hMass[icen]->Clone(Form("hMassEtaCen%d%s",icen,histoTag.Data()));
    
    hMassPi0[icen]->Divide(hMix[icen]);
    hMassPi0[icen]->Fit("pol0","Q","",0.25,0.35);
    Float_t scale = 1;
    if(hMassPi0[icen]->GetFunction("pol0")) scale = hMassPi0[icen]->GetFunction("pol0")->GetParameter(0);
    //printf("Scale factor %f for cen %d\n",scale,icen);
    hMassPi0[icen]->Scale(1./scale);
    hMassPi0[icen]->SetMarkerStyle(24);
    hMassPi0[icen]->SetMarkerColor(color[icen]);
    hMassPi0[icen]->SetLineColor(color[icen]);
    hMassPi0[icen]->SetAxisRange(0.04,0.24);
    hMassPi0[icen]->SetMarkerSize(0.5);
    
    hMassEta[icen]->Rebin(4);
    hMix    [icen]->Rebin(4);
    hMassEta[icen]->Divide(hMix[icen]);
    hMassEta[icen]->SetMarkerStyle(25);
    hMassEta[icen]->SetMarkerColor(color[icen]);
    hMassEta[icen]->SetLineColor(color[icen]);
    hMassEta[icen]->SetAxisRange(0.4,0.9);
    hMassEta[icen]->SetMarkerSize(0.5);
    hMassEta[icen]->Scale(1./scale);
    
    if(maxEta < hMassEta[icen]->GetMaximum()) maxEta = hMassEta[icen]->GetMaximum();
    if(maxPi0 < hMassPi0[icen]->GetMaximum()) maxPi0 = hMassPi0[icen]->GetMaximum();
  }

  //gPad->SetLogy();
  //gPad->SetGridy();
  hMassPi0[0]->SetMinimum(0.8);
  hMassPi0[0]->SetTitleOffset(1.6,"Y");
  hMassPi0[0]->SetYTitle("Real / Mixed");
  hMassPi0[0]->SetTitle("#pi^{0} peak, 2 < #it{E}_{pair}< 10 GeV");
  if     (histoTag.Contains("L0")) hMassPi0[0]->SetTitle("#pi^{0} peak, 5 < #it{E}_{pair}< 10 GeV");
  else if(histoTag.Contains("L2")) hMassPi0[0]->SetTitle("#pi^{0} peak, 8 < #it{E}_{pair}< 12 GeV");
  else if(histoTag.Contains("L1")) hMassPi0[0]->SetTitle("#pi^{0} peak, 10 < #it{E}_{pair}< 15 GeV");

  hMassPi0[0]->Draw();
  
  if(hMass[1]) // PbPb
  {
    hMassPi0[0]->SetMaximum(maxPi0*1.2);
    hMassPi0[5]->Draw("Hsame");
    hMassPi0[4]->Draw("Hsame");
    hMassPi0[3]->Draw("Hsame");
    hMassPi0[2]->Draw("Hsame");
    hMassPi0[1]->Draw("Hsame");
    hMassPi0[0]->Draw("Hsame");
    //hMass[6]->Draw("Hsame");
    //hMass[7]->Draw("same");
    //hMass[8]->Draw("same");
    //hMass[9]->Draw("same");
    
    TLegend l(0.12,0.6,0.4,0.85);
    l.SetTextSize(0.04);
    l.AddEntry(hMassPi0[0],"0-10%","P");
    l.AddEntry(hMassPi0[1],"10-20%","P");
    l.AddEntry(hMassPi0[2],"20-30%","P");
    l.AddEntry(hMassPi0[3],"30-40%","P");
    l.AddEntry(hMassPi0[4],"40-70%","P");
    l.AddEntry(hMassPi0[5],"50-60%","P");
    l.SetBorderSize(0);
    l.SetFillColor(0);
    l.Draw();
  }

  TLine l1(0.04,1,0.24,1);
  l1.Draw("same");
  
  // Pi0 invariant mass per EMCal super module
  cpi0->cd(3);
  
  TH1F* hSM   [20];
  TH1F* hMixSM[20];
  
  //Init to 0
  for(Int_t ism = 0; ism < 20; ism++)
  {
    hSM   [ism] = 0;
    hMixSM[ism] = 0;
  }
  
  binmin = hX->FindBin(4);  // Project histo from 3 GeV pairs
  binmax = hX->FindBin(20); // Project histo up to 20 GeV pairs
  Float_t maxSM = 0;
  
  Int_t first = 0;
  if(icalo==1) first = 12;
  
  for(Int_t ism = first; ism < 20; ism++)
  {
    TH2F* hTmpSM = (TH2F*) GetHisto(Form("AnaPi0_Calo%d_hReMod_%d",icalo,ism));    
    if(!hTmpSM) continue;
    
    hSM[ism] = (TH1F*) hTmpSM->ProjectionY(Form("hMassSM%d%s",ism,histoTag.Data()),binmin,binmax);
    hSM[ism]->Sumw2();
    hSM[ism]->SetMarkerStyle(26);
    hSM[ism]->Rebin(2);
    //hSM[ism]->Scale(1./hSM[ism]->Integral(0,10000));
    hSM[ism]->SetMarkerColor(color[ism-first]);
    hSM[ism]->SetLineColor(color[ism-first]);
    hSM[ism]->SetMarkerSize(0.5);

    TH2F* hTmpMixSM = (TH2F*) GetHisto(Form("AnaPi0_hMiMod_%d",ism));
    if(hTmpMixSM)
    {
      hMixSM[ism] = (TH1F*) hTmpMixSM->ProjectionY(Form("hMassMixSM%d%s",ism,histoTag.Data()),binmin,binmax);
      hMixSM[ism]->Sumw2();
      hMixSM[ism]->Rebin(2);
      hSM[ism]->Divide(hMixSM[ism]);
      hSM[ism]->Fit("pol0","Q","",0.25,0.35);
      Float_t scale = 1;
      if(hSM[ism]->GetFunction("pol0")) scale = hSM[ism]->GetFunction("pol0")->GetParameter(0);
      //printf("Scale factor %f for cen %d\n",scale,icen);
      hSM[ism]->Scale(1./scale);
      hSM[ism]->SetYTitle("Real / Mixed");
    }
    
    if(maxSM < hSM[ism]->GetMaximum()) maxSM = hSM[ism]->GetMaximum();
  }
  
  hSM[first]->SetTitle("#pi^{0} peak in SM, 4 < #it{E}_{pair}< 10 GeV");
  hSM[first]->SetTitleOffset(1.6,"Y");
  hSM[first]->SetAxisRange(0.04,0.24);
  hSM[first]->SetMaximum(maxSM*1.2);
  hSM[first]->SetMinimum(0.8);

  hSM[first]->Draw("H");
  TLegend lsm(0.12,0.5,0.35,0.85);
  lsm.SetTextSize(0.04);
  lsm.AddEntry(hSM[0],Form("Mod %d",0),"P");
  
  for(Int_t ism = first+1; ism < 20; ism++)
  {
    if(!hSM[ism]) continue;
    
    hSM[ism]->Draw("Hsame");
    lsm.AddEntry(hSM[ism],Form("Mod %d",ism),"P");
  }
  
  lsm.SetBorderSize(0);
  lsm.SetFillColor(0);
  lsm.Draw();
  
  l1.Draw("same");
  
  // Pi0 Invariant mass projection, in PbPb 6 centrality bins from 0 to 50%, all in pp
  cpi0->cd(4);
  
  //gPad->SetLogy();
  //gPad->SetGridy();
  hMassEta[0]->SetMinimum(0.8);
  hMassEta[0]->SetTitleOffset(1.6,"Y");
  hMassEta[0]->SetYTitle("Real / Mixed");
  hMassEta[0]->SetTitle("#eta peak, 2 <#it{E}_{pair}< 10 GeV");
  if     (histoTag.Contains("L0")) hMassEta[0]->SetTitle("#eta peak, 5 < #it{E}_{pair}< 10 GeV");
  else if(histoTag.Contains("L2")) hMassEta[0]->SetTitle("#eta peak, 8 < #it{E}_{pair}< 12 GeV");
  else if(histoTag.Contains("L1")) hMassEta[0]->SetTitle("#eta peak, 10 < #it{E}_{pair}< 15 GeV");

  hMassEta[0]->Draw("H");
  
  if(hMass[1]) // PbPb
  {
    hMassEta[0]->SetMaximum(maxEta*1.2);
    hMassEta[5]->Draw("Hsame");
    hMassEta[4]->Draw("Hsame");
    hMassEta[3]->Draw("Hsame");
    hMassEta[2]->Draw("Hsame");
    hMassEta[1]->Draw("Hsame");
    hMassEta[0]->Draw("Hsame");
    
    
    TLegend l2(0.12,0.6,0.4,0.85);
    l2.SetTextSize(0.04);
    l2.AddEntry(hMassEta[0],"0-10%","P");
    l2.AddEntry(hMassEta[1],"10-20%","P");
    l2.AddEntry(hMassEta[2],"20-30%","P");
    l2.AddEntry(hMassEta[3],"30-40%","P");
    l2.AddEntry(hMassEta[4],"40-70%","P");
    l2.AddEntry(hMassEta[5],"50-60%","P");
    l2.SetBorderSize(0);
    l2.SetFillColor(0);
    l2.Draw();
  }

  cpi0->Print(Form("%s_Pi0Histo.%s",histoTag.Data(),format.Data()));
  
  // cleanup
  //
  delete h2DMass;
  delete hX;
  
  for(Int_t icen=0; icen<10; icen++ )
  {
    if ( hMass   [icen] ) delete hMass   [icen];
    if ( hMix    [icen] ) delete hMix    [icen];
    if ( hMassPi0[icen] ) delete hMassPi0[icen];
    if ( hMassEta[icen] ) delete hMassEta[icen];
  }
  
  for(Int_t ism = first; ism < 20; ism++)
  {
    if ( hSM   [ism] ) delete hSM   [ism];
    if ( hMixSM[ism] ) delete hMixSM[ism];
  }
  
  delete cpi0;
}

///
/// Plot basic candidate cluster isolation histograms in 2 pads:
/// * Cluster spectra, isolated and not isolated
/// * pT distribution of tracks or clusters or both in the isolation cone, or perpendicular cones or eta-band out of cone
/// * sum pT distribution of tracks or clusters or both in the isolation cone, or perpendicular cones or eta-band out of cone
/// * sum pT distribution of tracks or clusters or both in the isolation cone whith subtracted eta-band sum pT
///
//__________________________________________________
void IsolQA(Int_t icalo)
{
  gStyle->SetPadRightMargin(0.02);

  TCanvas * cIsolation = new TCanvas(Form("IsolationHisto_%s",histoTag.Data()),"Isolation cone",1000,1000);
  cIsolation->Divide(2,2);
  
  Float_t minClusterE = 5;
  if      ( histoTag.Contains("L0") ) minClusterE =  5;
  else if ( histoTag.Contains("L2") ) minClusterE = 10;
  else if ( histoTag.Contains("L1") ) minClusterE = 12;
  
  TH1F * hIsolated    = (TH1F*) GetHisto(Form("AnaIsolPhoton_Calo%d_hPt"     ,icalo));
  TH1F * hNotIsolated = (TH1F*) GetHisto(Form("AnaIsolPhoton_Calo%d_hPtNoIso",icalo));
  
  if(!hIsolated) return;
  
  Int_t minClusterEBin = hIsolated->FindBin(minClusterE);
  Float_t nTrig = hIsolated->Integral(minClusterEBin,100000)+hNotIsolated->Integral(minClusterE,100000);
  
  if ( nTrig <=0 ) return ;
    
  //
  // Pt in cone
  //
  cIsolation->cd(1);

  hIsolated   ->Sumw2();
  hNotIsolated->Sumw2();
  hIsolated   ->SetMarkerColor(4);
  hNotIsolated->SetMarkerColor(2);
  hIsolated   ->SetLineColor  (4);
  hNotIsolated->SetLineColor  (2);
  hIsolated   ->SetMarkerStyle(24);
  hNotIsolated->SetMarkerStyle(20);

  hNotIsolated->SetTitle("(non) isolated cluster spectra, #it{R}=0.4, #Sigma #it{p}_{T}<2 GeV/#it{c}");

  hNotIsolated->Draw();  
  hIsolated   ->Draw("same");
  
  TLegend lI(0.4,0.5,0.8,0.88);
  lI.SetTextSize(0.04);
  lI.SetBorderSize(0);
  lI.SetFillColor(0);
  lI.AddEntry(hIsolated   ,"Isolated candidates","P");
  lI.AddEntry(hNotIsolated,"NOT Isolated candidates","P");
  lI.Draw("same");
  
  //
  // Pt in cone
  //
  cIsolation->cd(2);
  gPad->SetLogy();
  
  TLegend l(0.4,0.5,0.8,0.88);
  l.SetTextSize(0.04);
  l.SetBorderSize(0);
  l.SetFillColor(0);
    
  TH2F* h2PtInCone           = (TH2F*) GetHisto(Form("AnaIsolPhoton_Calo%d_hPtInCone"        ,icalo));
  TH2F* h2PtInConeCluster    = (TH2F*) GetHisto(Form("AnaIsolPhoton_Calo%d_hPtClusterInCone" ,icalo));
  TH2F* h2PtInConeTrack      = (TH2F*) GetHisto(Form("AnaIsolPhoton_Calo%d_hPtTrackInCone"   ,icalo));
  TH2F* h2PtInConeTrackPerp  = (TH2F*) GetHisto(Form("AnaIsolPhoton_Calo%d_hPtInPerpCone"    ,icalo));
  TH2F* h2PtInEtaBandTrack   = (TH2F*) GetHisto(Form("AnaIsolPhoton_Calo%d_hEtaBandClusterPt",icalo));
  TH2F* h2PtInEtaBandCluster = (TH2F*) GetHisto(Form("AnaIsolPhoton_Calo%d_hEtaBandClusterPt",icalo));
  
  TH1F* hPtInCone            = (TH1F*) h2PtInCone          ->ProjectionY(Form("PtInCone%2.1fGeV%s"  ,minClusterE,histoTag.Data()),minClusterEBin,10000);
  TH1F* hPtInConeCluster     = (TH1F*) h2PtInConeCluster   ->ProjectionY(Form("PtInConeCl%2.1fGeV%s",minClusterE,histoTag.Data()),minClusterEBin,10000);
  TH1F* hPtInConeTrack       = (TH1F*) h2PtInConeTrack     ->ProjectionY(Form("PtInConeTr%2.1fGeV%s",minClusterE,histoTag.Data()),minClusterEBin,10000);
  TH1F* hPtInConeTrackPerp   = (TH1F*) h2PtInConeTrackPerp ->ProjectionY(Form("PtInConePe%2.1fGeV%s",minClusterE,histoTag.Data()),minClusterEBin,10000);
  TH1F* hPtInEtaBandTrack    = (TH1F*) h2PtInEtaBandTrack  ->ProjectionY(Form("PtInConeET%2.1fGeV%s",minClusterE,histoTag.Data()),minClusterEBin,10000);
  TH1F* hPtInEtaBandCluster  = (TH1F*) h2PtInEtaBandCluster->ProjectionY(Form("PtInConeEC%2.1fGeV%s",minClusterE,histoTag.Data()),minClusterEBin,10000);

  hPtInCone          ->Sumw2();
  hPtInConeCluster   ->Sumw2();
  hPtInConeTrack     ->Sumw2();
  hPtInConeTrackPerp ->Sumw2();
  hPtInEtaBandCluster->Sumw2();
  hPtInEtaBandTrack  ->Sumw2();

  Int_t rb = 1;
  hPtInCone          ->Rebin(rb);
  hPtInConeCluster   ->Rebin(rb);
  hPtInConeTrack     ->Rebin(rb);
  hPtInConeTrackPerp ->Rebin(rb);
  hPtInEtaBandCluster->Rebin(rb);
  hPtInEtaBandTrack  ->Rebin(rb);

  hPtInCone          ->Scale(1./nTrig);
  hPtInConeCluster   ->Scale(1./nTrig);
  hPtInConeTrack     ->Scale(1./nTrig);
  hPtInConeTrackPerp ->Scale(1./nTrig);
  hPtInEtaBandCluster->Scale(1./nTrig);
  hPtInEtaBandTrack  ->Scale(1./nTrig);

  hPtInCone          ->SetAxisRange(0,10);
  hPtInConeCluster   ->SetAxisRange(0,10);
  hPtInConeTrack     ->SetAxisRange(0,10);
  hPtInConeTrackPerp ->SetAxisRange(0,10);
  hPtInEtaBandCluster->SetAxisRange(0,10);
  hPtInEtaBandTrack  ->SetAxisRange(0,10);
  
  hPtInCone          ->SetMarkerStyle(24);
  hPtInConeCluster   ->SetMarkerStyle(20);
  hPtInConeTrack     ->SetMarkerStyle(20);
  hPtInConeTrackPerp ->SetMarkerStyle(27);
  hPtInEtaBandCluster->SetMarkerStyle(21);
  hPtInEtaBandTrack  ->SetMarkerStyle(21);

  hPtInCone          ->SetMarkerColor(1);
  hPtInConeCluster   ->SetMarkerColor(2);
  hPtInConeTrack     ->SetMarkerColor(4);
  hPtInConeTrackPerp ->SetMarkerColor(4);
  hPtInEtaBandCluster->SetMarkerColor(2);
  hPtInEtaBandTrack  ->SetMarkerColor(4);

  hPtInCone          ->SetLineColor(1);
  hPtInConeCluster   ->SetLineColor(2);
  hPtInConeTrack     ->SetLineColor(4);
  hPtInConeTrackPerp ->SetLineColor(4);
  hPtInEtaBandCluster->SetLineColor(2);
  hPtInEtaBandTrack  ->SetLineColor(4);
  
  hPtInCone->SetTitleOffset(1.5,"Y");
  hPtInCone->SetYTitle("Entries / #it{N}_{candidates}");
  hPtInCone->SetTitle(Form("Track/cluster spectra in cone p_{T,cand}>%2.0f GeV/#it{c}, #it{R}=0.4",minClusterE));

  Float_t max = hPtInCone->GetMaximum();
  if(max < hPtInConeTrack     ->GetMaximum()) max = hPtInConeTrack     ->GetMaximum();
  if(max < hPtInConeCluster   ->GetMaximum()) max = hPtInConeCluster   ->GetMaximum();
  if(max < hPtInConeTrackPerp ->GetMaximum()) max = hPtInConeTrackPerp ->GetMaximum();
  if(max < hPtInEtaBandCluster->GetMaximum()) max = hPtInEtaBandCluster->GetMaximum();
  if(max < hPtInEtaBandTrack  ->GetMaximum()) max = hPtInEtaBandTrack  ->GetMaximum();
  hPtInCone->SetMaximum(max*2);
  
  hPtInCone          ->Draw("");
  hPtInConeCluster   ->Draw("same");
  hPtInConeTrack     ->Draw("same");
  hPtInConeTrackPerp ->Draw("same");
  hPtInEtaBandCluster->Draw("same");
  hPtInEtaBandTrack  ->Draw("same");
  
  l.AddEntry(hPtInCone          ,"Tracks+Clusters","P");
  l.AddEntry(hPtInConeCluster   ,"Clusters inside cone","P");
  l.AddEntry(hPtInConeTrack     ,"Tracks inside cone","P");
  l.AddEntry(hPtInConeTrackPerp ,"Tracks inside #perp cones","P");
  l.AddEntry(hPtInEtaBandTrack  ,"Tracks #eta band","P");
  l.AddEntry(hPtInEtaBandCluster,"Clusters #eta band","P");

  l.Draw("same");

  //
  // Sum Pt in cone
  //
  cIsolation->cd(3);
  gPad->SetLogy();
  
  TLegend l2(0.4,0.5,0.8,0.88);
  l2.SetTextSize(0.04);
  l2.SetBorderSize(0);
  l2.SetFillColor(0);
  
  TH2F* h2SumPtCone           = (TH2F*) GetHisto(Form("AnaIsolPhoton_Calo%d_hConePtSum"        ,icalo));
  TH2F* h2SumPtConeCluster    = (TH2F*) GetHisto(Form("AnaIsolPhoton_Calo%d_hConePtSumCluster" ,icalo));
  TH2F* h2SumPtConeTrack      = (TH2F*) GetHisto(Form("AnaIsolPhoton_Calo%d_hConePtSumTrack"   ,icalo));
  TH2F* h2SumPtConeTrackPerp  = (TH2F*) GetHisto(Form("AnaIsolPhoton_Calo%d_hPerpConePtSum"    ,icalo));
  TH2F* h2SumPtEtaBandTrack   = (TH2F*) GetHisto(Form("AnaIsolPhoton_Calo%d_hConeSumPtEtaUENormTrack"  ,icalo));
  TH2F* h2SumPtEtaBandCluster = (TH2F*) GetHisto(Form("AnaIsolPhoton_Calo%d_hConeSumPtEtaUENormCluster",icalo));
  TH2F* h2SumPtConeSubTrack   = (TH2F*) GetHisto(Form("AnaIsolPhoton_Calo%d_hConeSumPtEtaUESubTrack"   ,icalo));
  TH2F* h2SumPtConeSubCluster = (TH2F*) GetHisto(Form("AnaIsolPhoton_Calo%d_hConeSumPtEtaUESubCluster" ,icalo));
  TH2F* h2SumPtConeSub        = (TH2F*) GetHisto(Form("AnaIsolPhoton_Calo%d_hConeSumPtEtaUESub"        ,icalo));
  
  TH1F* hSumPtCone            = (TH1F*) h2SumPtCone          ->ProjectionY(Form("SumPtCone%2.1fGeV%s"   ,minClusterE,histoTag.Data()),minClusterEBin,10000);
  TH1F* hSumPtConeCluster     = (TH1F*) h2SumPtConeCluster   ->ProjectionY(Form("SumPtConeCl%2.1fGeV%s" ,minClusterE,histoTag.Data()),minClusterEBin,10000);
  TH1F* hSumPtConeTrack       = (TH1F*) h2SumPtConeTrack     ->ProjectionY(Form("SumPtConeTr%2.1fGeV%s" ,minClusterE,histoTag.Data()),minClusterEBin,10000);
  TH1F* hSumPtConeTrackPerp   = (TH1F*) h2SumPtConeTrackPerp ->ProjectionY(Form("SumPtConePe%2.1fGeV%s" ,minClusterE,histoTag.Data()),minClusterEBin,10000);
  TH1F* hSumPtEtaBandTrack    = (TH1F*) h2SumPtEtaBandTrack  ->ProjectionY(Form("SumPtConeET%2.1fGeV%s" ,minClusterE,histoTag.Data()),minClusterEBin,10000);
  TH1F* hSumPtEtaBandCluster  = (TH1F*) h2SumPtEtaBandCluster->ProjectionY(Form("SumPtConeEC%2.1fGeV%s" ,minClusterE,histoTag.Data()),minClusterEBin,10000);
  TH1F* hSumPtConeSub         = (TH1F*) h2SumPtConeSub       ->ProjectionY(Form("SumPtConeS%2.1fGeV%s"  ,minClusterE,histoTag.Data()),minClusterEBin,10000);
  TH1F* hSumPtConeSubCluster  = (TH1F*) h2SumPtConeSubCluster->ProjectionY(Form("SumPtConeSCl%2.1fGeV%s",minClusterE,histoTag.Data()),minClusterEBin,10000);
  TH1F* hSumPtConeSubTrack    = (TH1F*) h2SumPtConeSubTrack  ->ProjectionY(Form("SumPtConeSTr%2.1fGeV%s",minClusterE,histoTag.Data()),minClusterEBin,10000);
  
  hSumPtCone          ->Sumw2();
  hSumPtConeCluster   ->Sumw2();
  hSumPtConeTrack     ->Sumw2();
  hSumPtConeSub       ->Sumw2();
  hSumPtConeSubCluster->Sumw2();
  hSumPtConeSubTrack  ->Sumw2();
  hSumPtConeTrackPerp ->Sumw2();
  hSumPtEtaBandCluster->Sumw2();
  hSumPtEtaBandTrack  ->Sumw2();
  
  rb = 1;
  hSumPtCone          ->Rebin(rb);
  hSumPtConeCluster   ->Rebin(rb);
  hSumPtConeTrack     ->Rebin(rb);
  hSumPtConeSub       ->Rebin(rb);
  hSumPtConeSubCluster->Rebin(rb);
  hSumPtConeSubTrack  ->Rebin(rb);
  hSumPtConeTrackPerp ->Rebin(rb);
  hSumPtEtaBandCluster->Rebin(rb);
  hSumPtEtaBandTrack  ->Rebin(rb);
  
  hSumPtCone          ->Scale(1./nTrig);
  hSumPtConeCluster   ->Scale(1./nTrig);
  hSumPtConeTrack     ->Scale(1./nTrig);
  hSumPtConeSub       ->Scale(1./nTrig);
  hSumPtConeSubCluster->Scale(1./nTrig);
  hSumPtConeSubTrack  ->Scale(1./nTrig);
  hSumPtConeTrackPerp ->Scale(1./nTrig);
  hSumPtEtaBandCluster->Scale(1./nTrig);
  hSumPtEtaBandTrack  ->Scale(1./nTrig);
  
  hSumPtCone          ->SetAxisRange(0,30);
  hSumPtConeCluster   ->SetAxisRange(0,30);
  hSumPtConeTrack     ->SetAxisRange(0,30);
  hSumPtConeSub       ->SetAxisRange(-5,30);
  hSumPtConeSubCluster->SetAxisRange(-5,30);
  hSumPtConeSubTrack  ->SetAxisRange(-5,30);
  hSumPtConeTrackPerp ->SetAxisRange(0,30);
  hSumPtEtaBandCluster->SetAxisRange(0,30);
  hSumPtEtaBandTrack  ->SetAxisRange(0,30);
  
  hSumPtCone          ->SetMarkerStyle(24);
  hSumPtConeCluster   ->SetMarkerStyle(20);
  hSumPtConeTrack     ->SetMarkerStyle(20);
  hSumPtConeSub       ->SetMarkerStyle(25);
  hSumPtConeSubCluster->SetMarkerStyle(25);
  hSumPtConeSubTrack  ->SetMarkerStyle(25);
  hSumPtConeTrackPerp ->SetMarkerStyle(27);
  hSumPtEtaBandCluster->SetMarkerStyle(21);
  hSumPtEtaBandTrack  ->SetMarkerStyle(21);
  
  hSumPtCone          ->SetMarkerColor(1);
  hSumPtConeCluster   ->SetMarkerColor(2);
  hSumPtConeTrack     ->SetMarkerColor(4); 
  hSumPtConeSub       ->SetMarkerColor(1);
  hSumPtConeSubCluster->SetMarkerColor(2);
  hSumPtConeSubTrack  ->SetMarkerColor(4);
  hSumPtConeTrackPerp ->SetMarkerColor(4);
  hSumPtEtaBandCluster->SetMarkerColor(2);
  hSumPtEtaBandTrack  ->SetMarkerColor(4);
  
  hSumPtCone          ->SetLineColor(1);
  hSumPtConeCluster   ->SetLineColor(2);
  hSumPtConeTrack     ->SetLineColor(4);
  hSumPtConeSub       ->SetLineColor(1);
  hSumPtConeSubCluster->SetLineColor(2);
  hSumPtConeSubTrack  ->SetLineColor(4);
  hSumPtConeTrackPerp ->SetLineColor(4);
  hSumPtEtaBandCluster->SetLineColor(2);
  hSumPtEtaBandTrack  ->SetLineColor(4);
  
  hSumPtCone->SetTitleOffset(1.5,"Y");
  hSumPtCone->SetYTitle("Entries / #it{N}_{candidates}");
  hSumPtCone->SetTitle(Form("Track/cluster #Sigma #it{p}_{T}, p_{T,cand}>%2.0f GeV/#it{c}, #it{R}=0.4, ",minClusterE));
  
  max = hSumPtCone->GetMaximum();
  if(max < hSumPtConeTrack     ->GetMaximum()) max = hSumPtConeTrack     ->GetMaximum();
  if(max < hSumPtConeCluster   ->GetMaximum()) max = hSumPtConeCluster   ->GetMaximum();
  if(max < hSumPtConeSub       ->GetMaximum()) max = hSumPtConeSub       ->GetMaximum();
  if(max < hSumPtConeSubTrack  ->GetMaximum()) max = hSumPtConeSubTrack  ->GetMaximum();
  if(max < hSumPtConeSubCluster->GetMaximum()) max = hSumPtConeSubCluster->GetMaximum();
  if(max < hSumPtConeTrackPerp ->GetMaximum()) max = hSumPtConeTrackPerp ->GetMaximum();
  if(max < hSumPtEtaBandCluster->GetMaximum()) max = hSumPtEtaBandCluster->GetMaximum();
  if(max < hSumPtEtaBandTrack  ->GetMaximum()) max = hSumPtEtaBandTrack  ->GetMaximum();
  hSumPtCone->SetMaximum(max*2);

  hSumPtCone          ->Draw("");
  hSumPtConeCluster   ->Draw("same");
  hSumPtConeTrack     ->Draw("same");
//  hSumPtConeSub       ->Draw("same");
//  hSumPtConeSubCluster->Draw("same");
//  hSumPtConeSubTrack  ->Draw("same");
  hSumPtConeTrackPerp ->Draw("same");
  hSumPtEtaBandCluster->Draw("same");
  hSumPtEtaBandTrack  ->Draw("same");
  
  l2.AddEntry(hSumPtCone          ,"Tracks+Clusters","P");
  l2.AddEntry(hSumPtConeCluster   ,"Clusters inside cone","P");
  l2.AddEntry(hSumPtConeTrack     ,"Tracks inside cone","P");
  l2.AddEntry(hSumPtConeTrackPerp ,"Tracks inside #perp cones","P");
  l2.AddEntry(hSumPtEtaBandTrack  ,"Tracks #eta band","P");
  l2.AddEntry(hSumPtEtaBandCluster,"Clusters #eta band","P");
//  l2.AddEntry(hSumPtConeSub       ,"Tracks+Clusters-#eta band","P");
//  l2.AddEntry(hSumPtConeSubCluster,"Clusters inside cone-#eta band","P");
//  l2.AddEntry(hSumPtConeSubTrack  ,"Tracks inside cone-#eta band","P");
  
  l2.Draw("same");
  
  //
  // Sum Pt in cone, UE subtracted
  //
  cIsolation->cd(4);
  gPad->SetLogy();
  
  TLegend l3(0.4,0.7,0.8,0.88);
  l3.SetTextSize(0.04);
  l3.SetBorderSize(0);
  l3.SetFillColor(0);

  hSumPtConeSub->SetTitle(Form("Track/Cluster #Sigma #it{p}_{T}-#Sigma #eta band, p_{T,cand}>%2.0f GeV/#it{c}, #it{R}=0.4",minClusterE));

  hSumPtConeSub       ->Draw("");
  hSumPtConeSubCluster->Draw("same");
  hSumPtConeSubTrack  ->Draw("same");
  
  l3.AddEntry(hSumPtConeSub       ,"Tracks+Clusters-#eta band","P");
  l3.AddEntry(hSumPtConeSubCluster,"Clusters inside cone-#eta band","P");
  l3.AddEntry(hSumPtConeSubTrack  ,"Tracks inside cone-#eta band","P");

  l3.Draw("same");

  cIsolation->Print(Form("%s_IsolationHisto.%s",histoTag.Data(),format.Data()));
  
  // cleanup
  //
  delete hPtInCone            ;
  delete hPtInConeCluster     ;
  delete hPtInConeTrack       ;
  delete hPtInConeTrackPerp   ;
  delete hPtInEtaBandTrack    ;
  delete hPtInEtaBandCluster  ;

  delete hSumPtCone           ;
  delete hSumPtConeCluster    ;
  delete hSumPtConeTrack      ;
  delete hSumPtConeTrackPerp  ;
  delete hSumPtEtaBandTrack   ;
  delete hSumPtEtaBandCluster ;

  delete hSumPtConeSub        ;
  delete hSumPtConeSubCluster ;
  delete hSumPtConeSubTrack   ;
  
  delete cIsolation           ;
}

///
/// Plot basic cluster-track correlation histograms in 2 pads:
/// * Azimuthal correlation of high pT trigger cluster and tracks in 4 associated pT bins
/// * xE distribution of tracks correlated to a high pT cluster, in the opposite side or in perpendicular region (UE)
/// 
/// \param icalo: 0 EMCal, 1 DCal
///
//__________________________________________________
void CorrelQA(Int_t icalo)
{
  gStyle->SetPadRightMargin(0.02);

  TCanvas * cCorrelation = new TCanvas(Form("CorrelationHisto_%s",histoTag.Data()),"Trigger cluster - associated track correlation ",1000,500);
  cCorrelation->Divide(2,1);
  
  Float_t minClusterE = 5;
  if      ( histoTag.Contains("L0") ) minClusterE =  8;
  else if ( histoTag.Contains("L2") ) minClusterE = 10;
  else if ( histoTag.Contains("L1") ) minClusterE = 12;
     
  Float_t assocBins[] = {0.5,2.,5.,10.,20.};
  Int_t nAssocBins = 4;
  
  TH1F * hTrigger = (TH1F*) GetHisto(Form("AnaPhotonHadronCorr_Calo%d_hPtTrigger",icalo));
  
  if(!hTrigger) return;
  
  Int_t minClusterEBin = hTrigger->FindBin(minClusterE);
  Float_t nTrig = hTrigger->Integral(minClusterEBin,100000);
  
  if ( nTrig <=0 ) return ;
  
  //Azimuthal correlation
  cCorrelation->cd(1);
  gPad->SetLogy();
  TH1F* hDeltaPhi[4];
  for(Int_t i = 0; i < 4; i++) hDeltaPhi[i] = 0; 
  
  TLegend l(0.35,0.6,0.83,0.85);
  l.SetHeader(Form("p_{T,T} > %2.1f GeV/c",minClusterE));
  l.SetTextSize(0.04);
  l.SetBorderSize(0);
  l.SetFillColor(0);

  for(Int_t ibin = 0; ibin < nAssocBins; ibin++ )
  {
    TH2F* hDeltaPhiE = 
    (TH2F*) GetHisto(Form("AnaPhotonHadronCorr_Calo%d_hDeltaPhiPtAssocPt%2.1f_%2.1f",icalo,assocBins[ibin],assocBins[ibin+1]));
    hDeltaPhi[ibin]  = 
    (TH1F*) hDeltaPhiE->ProjectionY(Form("DeltaPhi%2.1f%s",assocBins[ibin],histoTag.Data()),minClusterEBin,10000);
    hDeltaPhi[ibin]->Sumw2();
    hDeltaPhi[ibin]->Rebin(2);
    hDeltaPhi[ibin]->Scale(1./nTrig);
    
    hDeltaPhi[ibin]->Fit("pol0","Q","",1,2);
    
    Float_t scale = 1;
    if(hDeltaPhi[ibin]->GetFunction("pol0"))
    {
      scale = hDeltaPhi[ibin]->GetFunction("pol0")->GetParameter(0);
      hDeltaPhi[ibin]->GetFunction("pol0")->SetRange(6,7); // move from plot
    }
    hDeltaPhi[ibin]->Scale(1./scale);
    //printf("ibin %d, scale %f\n",ibin,scale);
    
    hDeltaPhi[ibin]->SetAxisRange(-1.6,4.7);
    
    hDeltaPhi[ibin]->SetMarkerStyle(24);
    hDeltaPhi[ibin]->SetMarkerColor(color[ibin]);
    hDeltaPhi[ibin]->SetLineColor(color[ibin]);
    hDeltaPhi[ibin]->SetTitleOffset(1.5,"Y");
    hDeltaPhi[ibin]->SetYTitle("#it{N}_{pairs} / #it{N}_{trig} / ZYAM");
    hDeltaPhi[ibin]->SetTitle("#gamma (#lambda_{0}^{2} < 0.4, neutral cluster) trigger");
    
    l.AddEntry(hDeltaPhi[ibin],Form("%2.1f<#it{p}_{T,A}< %2.0f GeV/c",assocBins[ibin],assocBins[ibin+1]),"P");
  }
  
  hDeltaPhi[2]->SetMaximum(hDeltaPhi[2]->GetMaximum()*10);
  hDeltaPhi[2]->SetMinimum(0.8);
  
  hDeltaPhi[2]->Draw("H");
  hDeltaPhi[1]->Draw("Hsame");
  hDeltaPhi[3]->Draw("Hsame");
  hDeltaPhi[0]->Draw("Hsame");
  
  l.Draw("same");
  
  // xE correlation
  cCorrelation->cd(2);
  gPad->SetLogy();
  
  TLegend l2(0.35,0.6,0.83,0.85);
  l2.SetHeader(Form("p_{T,T} > %2.0f GeV/c",minClusterE));
  l2.SetTextSize(0.04);
  l2.SetBorderSize(0);
  l2.SetFillColor(0);
  
  TH2F* hEXE   = (TH2F*) GetHisto(Form("AnaPhotonHadronCorr_Calo%d_hXECharged"  ,icalo));
  TH2F* hEXEUE = (TH2F*) GetHisto(Form("AnaPhotonHadronCorr_Calo%d_hXEUeCharged",icalo));
  
  TH1F* hXE  = (TH1F*) hEXE->ProjectionY(Form("XE%2.1fGeV%s",minClusterE,histoTag.Data()),minClusterEBin,10000);
  hXE->Sumw2();
  hXE->Rebin(2);
  hXE->Scale(1./nTrig);
  hXE->SetAxisRange(0,1);
  hXE->SetMarkerStyle(24);
  hXE->SetMarkerColor(1);
  hXE->SetLineColor(1);
  hXE->SetTitleOffset(1.5,"Y");
  hXE->SetYTitle("#it{N}_{pairs} / #it{N}_{trig}");
  hXE->SetTitle("#gamma (#lambda_{0}^{2} < 0.4, neutral cluster) trigger");
  l2.AddEntry(hXE,"raw x_{E}","P");
  hXE->Draw();

  TH1F* hXEUE  = (TH1F*) hEXEUE->ProjectionY(Form("XEUE%2.1fGeV%s",minClusterE,histoTag.Data()),minClusterEBin,10000);
  hXEUE->Sumw2();
  hXEUE->Rebin(2);
  hXEUE->Scale(1./nTrig);
  hXEUE->SetAxisRange(0,1);
  hXEUE->SetMarkerStyle(25);
  hXEUE->SetMarkerColor(2);
  hXEUE->SetLineColor(2);
  l2.AddEntry(hXEUE,"raw Und. Event x_{E}","P");
  hXEUE->Draw("same");
  
  l2.Draw("same");

  cCorrelation->Print(Form("%s_CorrelationHisto.%s",histoTag.Data(),format.Data()));
  
  // cleanup
  //  
  for(Int_t i = 0; i < 4; i++) delete hDeltaPhi[i]; 
  
  delete hXE  ;
  delete hXEUE;
  
  delete cCorrelation;
}

///
/// Plot basic generated particle distribution histograms in 4 pads:
/// * pT spectra of generated and reconstructed clusters from photon, photon from pi0, photon from eta, eta and pi0  
/// * ratio of pT spectra of reconstructed over generated
/// * azimuthal distribution of generated particles
/// * pseudorapidity distribution of generated particles
///
/// \param icalo: 0 EMCal, 1 DCal
///
//________________________________________________________
void MCQA(Int_t icalo)
{
  TH1F* hClusterPho    = (TH1F*) GetHisto(Form("AnaPhoton_Calo%d_hPt_MCPhoton",icalo));    
  TH1F* hClusterPi0    = (TH1F*) GetHisto(Form("AnaPhoton_Calo%d_hPt_MCPi0"   ,icalo));       
  TH1F* hClusterEta    = (TH1F*) GetHisto(Form("AnaPhoton_Calo%d_hPt_MCEta"   ,icalo));       
  TH1F* hClusterPhoPi0 = (TH1F*) GetHisto(Form("AnaPhoton_Calo%d_hPt_MCPhotonPi0Decay",icalo));       
  TH1F* hClusterPhoEta = (TH1F*) GetHisto(Form("AnaPhoton_Calo%d_hPt_MCPhotonEtaDecay",icalo)); 
  
  if(!hClusterPho) return;
  
  TH1F* hPrimPho    = (TH1F*) GetHisto(Form("AnaPhoton_Calo%d_hPtPrim_MCPhoton"        ,icalo));
  TH1F* hPrimPhoPi0 = (TH1F*) GetHisto(Form("AnaPhoton_Calo%d_hPtPrim_MCPhotonPi0Decay",icalo));
  TH1F* hPrimPhoEta = (TH1F*) GetHisto(Form("AnaPhoton_Calo%d_hPtPrim_MCPhotonEtaDecay",icalo));
  TH1F* hPrimPi0    = (TH1F*) GetHisto(Form("AnaPi0_Calo%d_hPrimPi0Pt",icalo));
  TH1F* hPrimEta    = (TH1F*) GetHisto(Form("AnaPi0_Calo%d_hPrimEtaPt",icalo));
  
  TCanvas * cmc = new TCanvas(Form("MCHisto_%s",histoTag.Data()),"Cluster origin",1000,1000);
  cmc->Divide(2,2);
  
  cmc->cd(1);
  gPad->SetLogy();
  
  hClusterPho->SetTitle("Cluster origin spectra, primary spectra in Calo acceptance");
  hClusterPho->Sumw2();
  hClusterPho->SetMarkerColor(1);
  hClusterPho->SetMarkerStyle(20);
  hClusterPho->SetAxisRange(0.,50.,"X");
  //hClusterPho->SetXTitle("E_{rec,gen} (GeV)");
  hClusterPho->SetXTitle("E_{rec}, p_{T,gen} (GeV)");
  hClusterPho->Draw("");

  hClusterPhoPi0->Sumw2();
  hClusterPhoPi0->SetMarkerColor(4);
  hClusterPhoPi0->SetMarkerStyle(20);
  hClusterPhoPi0->Draw("same");

  hClusterPhoEta->Sumw2();
  hClusterPhoEta->SetMarkerColor(2);
  hClusterPhoEta->SetMarkerStyle(20);
  hClusterPhoEta->Draw("same");
  
  hClusterPi0->Sumw2();
  hClusterPi0->SetMarkerColor(4);
  hClusterPi0->SetMarkerStyle(21);
  hClusterPi0->Draw("same");

  hClusterEta->Sumw2();
  hClusterEta->SetMarkerColor(2);
  hClusterEta->SetMarkerStyle(22);
  hClusterEta->Draw("same");
  
  hPrimPho->Sumw2();
  hPrimPho->SetMarkerColor(1);
  hPrimPho->SetMarkerStyle(24);
  hPrimPho->Draw("same");

  hPrimPhoPi0->Sumw2();
  hPrimPhoPi0->SetMarkerColor(4);
  hPrimPhoPi0->SetMarkerStyle(24);
  hPrimPhoPi0->Draw("same");
  
  hPrimPhoEta->Sumw2();
  hPrimPhoEta->SetMarkerColor(2);
  hPrimPhoEta->SetMarkerStyle(24);
  hPrimPhoEta->Draw("same");
  
  hPrimPi0->Sumw2();
  hPrimPi0->SetMarkerColor(4);
  hPrimPi0->SetMarkerStyle(25);
  hPrimPi0->Draw("same");
  
  hPrimEta->Sumw2();
  hPrimEta->SetMarkerColor(2);
  hPrimEta->SetMarkerStyle(26);
  hPrimEta->Draw("same");
  
  TLegend lR(0.5,0.5,0.7,0.89);
  lR.SetHeader("reco");
  lR.SetTextSize(0.04);
  lR.AddEntry(hClusterPho,"#gamma","P");
  lR.AddEntry(hClusterPhoPi0,"#gamma_{#pi^{0}}","P");
  lR.AddEntry(hClusterPhoEta,"#gamma_{#eta}","P");
  lR.AddEntry(hClusterPi0,"#pi^{0}","P");
  lR.AddEntry(hClusterEta,"#eta","P");
  lR.SetBorderSize(0);
  lR.SetFillColor(0);
  lR.Draw();
  
  TLegend lG(0.7,0.5,0.83,0.89);
  lG.SetHeader("gener");
  lG.SetTextSize(0.04);
  lG.AddEntry(hPrimPho,"#gamma","P");
  lG.AddEntry(hPrimPhoPi0,"#gamma_{#pi^{0}}","P");
  lG.AddEntry(hPrimPhoEta,"#gamma_{#eta}","P");
  lG.AddEntry(hPrimPi0,"#pi^{0}","P");
  lG.AddEntry(hPrimEta,"#eta","P");
  lG.SetBorderSize(0);
  lG.SetFillColor(0);
  lG.Draw();

  cmc->cd(2);
  gPad->SetLogy();
  TH1F* hRatPho    = (TH1F*) hClusterPho   ->Clone(Form("hGenRecoPho%s"   ,histoTag.Data()));
  TH1F* hRatPi0    = (TH1F*) hClusterPi0   ->Clone(Form("hGenRecoPi0%s"   ,histoTag.Data()));
  TH1F* hRatEta    = (TH1F*) hClusterEta   ->Clone(Form("hGenRecoEta%s"   ,histoTag.Data()));
  TH1F* hRatPhoPi0 = (TH1F*) hClusterPhoPi0->Clone(Form("hGenRecoPhoPi0%s",histoTag.Data()));
  TH1F* hRatPhoEta = (TH1F*) hClusterPhoEta->Clone(Form("hGenRecoPhoEta%s",histoTag.Data()));

  hRatPho   ->Divide(hPrimPho);
  hRatPhoPi0->Divide(hPrimPhoPi0);
  hRatPhoEta->Divide(hPrimPhoEta);
  hRatPi0   ->Divide(hPrimPi0);
  hRatEta   ->Divide(hPrimEta);
  
  hRatPho->SetTitle("Reconstructed cluster / Generated particle in Calo acc.");
  hRatPho->SetYTitle("Ratio");
  hRatPho->SetXTitle("E(GeV)");
  hRatPho->SetMinimum(1e-3);
  hRatPho->SetMaximum(20);
  hRatPho->Draw("");
  hRatPhoPi0->Draw("same");
  hRatPhoEta->Draw("same");  
  hRatPi0->Draw("same");
  hRatEta->Draw("same");

  TLegend l2(0.15,0.62,0.3,0.89);
  l2.SetTextSize(0.04);
  l2.AddEntry(hRatPho,"#gamma","P");
  l2.AddEntry(hRatPhoPi0,"#gamma_{#pi^{0}}","P");
  l2.AddEntry(hRatPhoEta,"#gamma_{#eta}","P"); 
  l2.AddEntry(hRatPi0,"#pi^{0}","P");
  l2.AddEntry(hRatEta,"#eta","P");
  l2.SetBorderSize(0);
  l2.SetFillColor(0);
  l2.Draw();

  cmc->cd(3);
  //gPad->SetLogy();

  TH2F* h2PrimPhoPhi = (TH2F*) GetHisto(Form("AnaPhoton_Calo%d_hPhiPrim_MCPhoton",icalo));
  TH2F* h2PrimPi0Phi = (TH2F*) GetHisto(Form("AnaPi0_Calo%d_hPrimPi0Phi"         ,icalo));
  TH2F* h2PrimEtaPhi = (TH2F*) GetHisto(Form("AnaPi0_Calo%d_hPrimEtaPhi"         ,icalo));
  
  Int_t binMin = hPrimPho->FindBin(3);
  
  TH1F* hPrimPhoPhi = (TH1F*) h2PrimPhoPhi->ProjectionY(Form("PrimPhoPhi%s",histoTag.Data()),binMin,1000);
  TH1F* hPrimPi0Phi = (TH1F*) h2PrimPi0Phi->ProjectionY(Form("PrimPi0Phi%s",histoTag.Data()),binMin,1000);
  TH1F* hPrimEtaPhi = (TH1F*) h2PrimEtaPhi->ProjectionY(Form("PrimEtaPhi%s",histoTag.Data()),binMin,1000);

  hPrimPhoPhi->Sumw2();
  hPrimPi0Phi->Sumw2();
  hPrimEtaPhi->Sumw2();

  hPrimPhoPhi->Scale(1./hPrimPhoPhi->Integral(0,1000));
  hPrimPi0Phi->Scale(1./hPrimPi0Phi->Integral(0,1000));
  hPrimEtaPhi->Scale(1./hPrimEtaPhi->Integral(0,1000));

  Float_t maxPhi = hPrimPhoPhi->GetMaximum();
  if(maxPhi < hPrimPi0Phi->GetMaximum()) maxPhi =  hPrimPi0Phi->GetMaximum();
  if(maxPhi < hPrimEtaPhi->GetMaximum()) maxPhi =  hPrimEtaPhi->GetMaximum();

  Float_t minPhi = hPrimPhoPhi->GetMinimum();
  if(minPhi > hPrimPi0Phi->GetMinimum()) minPhi =  hPrimPi0Phi->GetMinimum();
  if(minPhi > hPrimEtaPhi->GetMinimum()) minPhi =  hPrimEtaPhi->GetMinimum();

  hPrimPi0Phi->SetMaximum(maxPhi*1.1);
  hPrimPi0Phi->SetMinimum(minPhi);
  TGaxis::SetMaxDigits(3);

  hPrimPi0Phi->SetYTitle("1/total entries dN/d#varphi");
  hPrimPi0Phi->SetTitle("Generated particles #varphi for E > 3 GeV");
  hPrimPi0Phi->SetTitleOffset(1.6,"Y");
  hPrimPi0Phi->SetMarkerColor(4);
  hPrimPi0Phi->SetMarkerStyle(21);
  hPrimPi0Phi->Draw("");
  
  hPrimPhoPhi->SetMarkerColor(1);
  hPrimPhoPhi->SetMarkerStyle(20);
  Float_t scale = TMath::RadToDeg();
  ScaleXaxis(hPrimPhoPhi, TMath::RadToDeg());
  hPrimPhoPhi->Draw("same");

  hPrimEtaPhi->SetMarkerColor(2);
  hPrimEtaPhi->SetMarkerStyle(22);
  hPrimEtaPhi->Draw("same");

  cmc->cd(4);
  //gPad->SetLogy();
  
  TH2F* h2PrimPhoEtaP = (TH2F*) GetHisto(Form("AnaPhoton_Calo%d_hYPrim_MCPhoton",icalo));
  TH2F* h2PrimPi0EtaP = (TH2F*) GetHisto(Form("AnaPi0_Calo%d_hPrimPi0Rapidity"  ,icalo));
  TH2F* h2PrimEtaEtaP = (TH2F*) GetHisto(Form("AnaPi0_Calo%d_hPrimEtaRapidity"  ,icalo));

  h2PrimPhoEtaP->Sumw2();
  h2PrimEtaEtaP->Sumw2();
  h2PrimPi0EtaP->Sumw2();
  
  binMin = hPrimPho->FindBin(3);
  
  TH1F* hPrimPhoEtaP = (TH1F*) h2PrimPhoEtaP->ProjectionY(Form("PrimPhoEtaP%s",histoTag.Data()),binMin,1000);
  TH1F* hPrimPi0EtaP = (TH1F*) h2PrimPi0EtaP->ProjectionY(Form("PrimPi0EtaP%s",histoTag.Data()),binMin,1000);
  TH1F* hPrimEtaEtaP = (TH1F*) h2PrimEtaEtaP->ProjectionY(Form("PrimEtaEtaP%s",histoTag.Data()),binMin,1000);
  
  hPrimPhoEtaP->Scale(1./hPrimPhoEtaP->Integral(0,1000));
  hPrimPi0EtaP->Scale(1./hPrimPi0EtaP->Integral(0,1000));
  hPrimEtaEtaP->Scale(1./hPrimEtaEtaP->Integral(0,1000));
  
  Float_t maxEta = hPrimPhoEtaP->GetMaximum();
  if(maxEta < hPrimPi0EtaP->GetMaximum()) maxEta = hPrimPi0EtaP->GetMaximum();
  if(maxEta < hPrimEtaEtaP->GetMaximum()) maxEta = hPrimEtaEtaP->GetMaximum();
  
  Float_t minEta = hPrimPhoEtaP->GetMinimum();
  if(minEta > hPrimPi0EtaP->GetMinimum()) minEta = hPrimPi0EtaP->GetMinimum();
  if(minEta > hPrimEtaEtaP->GetMinimum()) minEta = hPrimEtaEtaP->GetMinimum();
  
  hPrimPi0EtaP->SetMaximum(maxEta*1.1);
  hPrimPi0EtaP->SetMinimum(minEta);
  TGaxis::SetMaxDigits(3);
  
  hPrimPi0EtaP->SetYTitle("1/total entries dN/d#eta");
  hPrimPi0EtaP->SetTitle("Generated particles #eta for E > 3 GeV");
  hPrimPi0EtaP->SetTitleOffset(1.6,"Y");
  hPrimPi0EtaP->SetMarkerColor(4);
  hPrimPi0EtaP->SetMarkerStyle(21);
  hPrimPi0EtaP->Draw("");
  
  hPrimPhoEtaP->SetMarkerColor(1);
  hPrimPhoEtaP->SetMarkerStyle(20);
  scale = TMath::RadToDeg();
  hPrimPhoEtaP->Draw("same");
  
  hPrimEtaEtaP->SetMarkerColor(2);
  hPrimEtaEtaP->SetMarkerStyle(22);
  hPrimEtaEtaP->Draw("same");
  
  cmc->Print(Form("%s_MCHisto.%s",histoTag.Data(),format.Data()));
  
  // cleanup
  //
  delete hPrimPhoPhi  ;
  delete hPrimPi0Phi  ;
  delete hPrimEtaPhi  ;
  delete hPrimPhoEtaP ;
  delete hPrimPi0EtaP ;
  delete hPrimEtaEtaP ;

  delete cmc ;
}

///
/// Open the file and list containing the histograms
/// 
/// \param trigName: name of list of histograms for a particular trigger
/// \param exportToFile: put the list of histograms in a separate file if true
///
//____________________________________________________________________
Bool_t GetList(TString trigName, Bool_t  exportToFile)
{  
  if(list) delete list;
  
  list = (TList*) dir->Get(trigName);
  
  if ( !list ) 
  { 
    printf("List not found, do nothing\n");
    return kFALSE; 
  }

  if ( list->GetEntries() <= 0 ) 
  { 
    printf("No histograms found <%d>, do nothing\n",list->GetEntries());
    return kFALSE; 
  }

  if ( exportToFile )
  {
    TFile * outputFile = new TFile("AnalysisResultsList.root","RECREATE");
    list->Write();
    outputFile->Close();
  }
  
  return kTRUE ;
}

///
/// Check if the list is available,
/// if not get the histo directly from file
///
/// \return the histogram with the provided name
///
/// \param histoName: histogram name
///
//___________________________________
TObject * GetHisto(TString histoName)
{
  if(list) return list->FindObject(histoName);
  else     return file->Get       (histoName);
}

///
/// Scale axis by a constant factor
/// used just to scale degrees to rad in a single histogram in the MC case
///
//___________________________________________________
void ScaleAxis(TAxis *a, Double_t scale)
{
  if (!a) return; // just a precaution
  if (a->GetXbins()->GetSize())
  {
    // an axis with variable bins
    // note: bins must remain in increasing order, hence the "Scale"
    // function must be strictly (monotonically) increasing
    TArrayD X(*(a->GetXbins()));
    for(Int_t i = 0; i < X.GetSize(); i++) X[i] = scale*X[i];
    a->Set((X.GetSize() - 1), X.GetArray()); // new Xbins
  }
  else
  {
    // an axis with fix bins
    // note: we modify Xmin and Xmax only, hence the "Scale" function
    // must be linear (and Xmax must remain greater than Xmin)
    a->Set(a->GetNbins(),
           -           scale*a->GetXmin(), // new Xmin
           -           scale*a->GetXmax()); // new Xmax
  }
  return;
}

///
/// Scale x axis by a constant factor
/// used just to scale degrees to rad in a single histogram in the MC case
///
//___________________________________________________
void ScaleXaxis(TH1 *h, Double_t scale)
{
  if (!h) return; // just a precaution
  ScaleAxis(h->GetXaxis(), scale);
  return;
}



