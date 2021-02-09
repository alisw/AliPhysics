///
/// \file DrawPromptPhotonEfficiency.C
/// \ingroup CaloTrackCorrMacrosPlotting
/// \brief Calculate isolated photon purity
///
/// Example macro to calculate isolated photon efficiency
/// Reconstructed clusters input from TH2 histogram pt vs shower shape
/// Input histograms from analysis of gamma-jet MC production
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
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#endif

///
/// Main method
///
void DrawPromptPhotonEfficiency()
{
  // Open histograms from Gamma-Jet production
  TFile * f = TFile::Open("GJ/Scaled.root");
  
  if ( !f ) return;
  
  Int_t color[] = {1,4,2,6,8,kCyan};
  Int_t rb = 10;
  Float_t markerSize = 1.5;
  
  //------------------------------------------
  // Generated spectra in analysis acceptance
  //------------------------------------------
  
  TH1F * hPromptPrimIso =  (TH1F*) f->Get("AnaIsolPhoton_hPtPrim_MCisoPhotonPrompt"); // Generated and isolated pT
  TH1F * hPromptPrim    =  (TH1F*) f->Get("AnaIsolPhoton_hPtPrim_MCPhotonPrompt"   ); // Generated,  all 
  
  hPromptPrimIso->Rebin(rb);
  hPromptPrim   ->Rebin(rb);
  
  // Ratio
  TH1F* hPromptPrimIsoEff = (TH1F*) hPromptPrimIso->Clone("hPromptPrimIsoEff");
  hPromptPrimIsoEff->Divide(hPromptPrim);

  //-----------------------------
  // Reconstructed spectra vs M02
  //-----------------------------
  
  TH2F * hPromptM02Iso   =  (TH2F*) f->Get("AnaIsolPhoton_hPtLambda0Iso_MCPhotonPrompt"  ); // Reco M02 vs pT, isolated
  TH2F * hPromptM02NoIso =  (TH2F*) f->Get("AnaIsolPhoton_hPtLambda0NoIso_MCPhotonPrompt"); // Reco M02 vs pT, not isolated
   
  // Open shower shape cut: reconstruction efficiency
  //
  TH1F * hPromptIsoNoCut      = (TH1F*) hPromptM02Iso->ProjectionX("hPromptIsoNoCut",
                                                                   hPromptM02Iso->GetYaxis()->FindBin(0.1),
                                                                   hPromptM02Iso->GetYaxis()->FindBin(10));
  TH1F * hPromptNoIsoNoCut    = (TH1F*) hPromptM02NoIso->ProjectionX("hPromptNoIsoNoCut",
                                                                     hPromptM02NoIso->GetYaxis()->FindBin(0.1),
                                                                     hPromptM02NoIso->GetYaxis()->FindBin(10));
  
  // Fixed shower shape photon cut: PID efficiency
  // Note that in the analysis the shower shape cut is not fixed, apply a pT dependent cut.
  // Here for simplicity just a cut at 0.3
  //  
  TH1F * hPromptIsoPhotonId   = (TH1F*) hPromptM02Iso->ProjectionX("hPromptIsoPhotonId",
                                                                   hPromptM02Iso->GetYaxis()->FindBin(0.1),
                                                                   hPromptM02Iso->GetYaxis()->FindBin(0.3));
  TH1F * hPromptNoIsoPhotonId = (TH1F*) hPromptM02NoIso->ProjectionX("hPromptNoIsoPhotonId",
                                                                     hPromptM02NoIso->GetYaxis()->FindBin(0.1),
                                                                     hPromptM02NoIso->GetYaxis()->FindBin(0.3));
  // Rebin
  hPromptIsoNoCut     ->Rebin(rb);
  hPromptIsoPhotonId  ->Rebin(rb);
  hPromptNoIsoNoCut   ->Rebin(rb);
  hPromptNoIsoPhotonId->Rebin(rb);

  //----------------------------------------
  // Add Isolated and non isolated histograms
  //-----------------------------------------
  TH1F * hPromptNoCut =  (TH1F*) hPromptIsoNoCut->Clone("hPromptSumNoCut");
  hPromptNoCut->Add(hPromptNoIsoNoCut);
  
  TH1F * hPromptPhotonId =  (TH1F*) hPromptIsoPhotonId->Clone("hPromptSumPhotonId");
  hPromptPhotonId->Add(hPromptNoIsoPhotonId);
  
  // Efficiency Ratios
  //
  TH1F* hPromptPIDEff = (TH1F*) hPromptPhotonId->Clone("hPromptPIDEff");
  hPromptPIDEff->Divide(hPromptPrim);
 
  TH1F* hPromptNoCutEff = (TH1F*) hPromptNoCut->Clone("hPromptNoCutEff");
  hPromptNoCutEff->Divide(hPromptPrim);
  
  TH1F* hPromptIsoNoCutEff = (TH1F*) hPromptIsoNoCut->Clone("hPromptIsoNoCutEff");
  hPromptIsoNoCutEff->Divide(hPromptPrim);  
  
  TH1F* hPromptIsoNoCutEff2 = (TH1F*) hPromptIsoNoCut->Clone("hPromptIsoNoCutEff2");
  hPromptIsoNoCutEff2->Divide(hPromptPrimIso);  
   
  TH1F* hPromptIsoPIDEff = (TH1F*) hPromptIsoPhotonId->Clone("hPromptPIDEff");
   hPromptIsoPIDEff->Divide(hPromptPrim);  

   //
   // Efficiency to be applied to data
   //
   TH1F* hPromptIsoPIDEff2 = (TH1F*) hPromptIsoPhotonId->Clone("hPromptPIDEff2");
   hPromptIsoPIDEff2->Divide(hPromptPrimIso);  
  
  
  //------------------------------------------------------------
  //------------------ Plot Efficiency -------------------------
  //------------------------------------------------------------
  
  gStyle->SetOptTitle(0);
  gStyle->SetTitleOffset(2.0,"Y");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(000000);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetPadLeftMargin(0.1);
  gStyle->SetPadTopMargin(0.02);
  //gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  
  TCanvas * cEff = new TCanvas("cEff","Efficiency",2*800,1*600);
  
  cEff->Divide(2,1);
  
  cEff->cd(1);
  //gPad->SetLogy();
  //gPad->SetGridy();
  
  hPromptPIDEff->Draw();
  hPromptPIDEff->SetAxisRange(10,100,"X");
  hPromptPIDEff->SetYTitle("selected reconstructed prompt /  generated prompt");  
  hPromptPIDEff->SetMaximum(1.00);
  hPromptPIDEff->SetMinimum(0.00);

  hPromptPIDEff->SetTitleOffset(1.5,"Y");
  hPromptPIDEff->SetMarkerStyle(28);
  hPromptPIDEff->SetMarkerColor(2);
  hPromptPIDEff->SetLineColor(2);
  hPromptPIDEff->SetMarkerSize(markerSize);

  hPromptIsoPIDEff->Draw("same");
  hPromptIsoPIDEff->SetMarkerStyle(27);
  hPromptIsoPIDEff->SetMarkerColor(4);
  hPromptIsoPIDEff->SetLineColor(4);
  hPromptIsoPIDEff->SetMarkerSize(markerSize);

  hPromptNoCutEff->Draw("same");
  hPromptNoCutEff->SetMarkerStyle(21);
  hPromptNoCutEff->SetMarkerColor(kGreen-2);
  hPromptNoCutEff->SetLineColor(kGreen-2);
  hPromptNoCutEff->SetMarkerSize(markerSize);

  hPromptPrimIsoEff->Draw("same");
  hPromptPrimIsoEff->SetMarkerStyle(24);
  hPromptPrimIsoEff->SetMarkerColor(1);
  hPromptPrimIsoEff->SetLineColor(1);
  hPromptPrimIsoEff->SetMarkerSize(markerSize);

  //hPromptIsoNoCutEff->Draw("same");
  hPromptIsoNoCutEff->SetMarkerStyle(25);
  hPromptIsoNoCutEff->SetMarkerColor(4);
  hPromptIsoNoCutEff->SetLineColor(4);
  hPromptIsoNoCutEff->SetMarkerSize(markerSize);

  //hPromptIsoNoCutEff2->Draw("same");
  hPromptIsoNoCutEff2->SetMarkerStyle(25);
  hPromptIsoNoCutEff2->SetMarkerColor(2);
  hPromptIsoNoCutEff2->SetLineColor(2);
  hPromptIsoNoCutEff2->SetMarkerSize(markerSize);

  //
  // Efficiency to be applied to data
  //
  //   hPromptIsoPIDEff2->Draw("same");
  //   hPromptIsoPIDEff2->SetMarkerStyle(21);
  //   hPromptIsoPIDEff2->SetMarkerColor(2);
  //   hPromptIsoPIDEff2->SetLineColor(2);
  //   hPromptIsoPIDEff2->SetMarkerSize(markerSize);

  TLegend *legEff = new TLegend(0.12,0.12,0.5,0.4);
  legEff->SetTextSize(0.045);
  legEff->SetLineColor(0);
  legEff->SetBorderSize(0);
  legEff->AddEntry(hPromptPrimIsoEff  ,"Isolated at generator level","PL");
  legEff->AddEntry(hPromptNoCutEff    ,"Neutral cluster","PL");
  legEff->AddEntry(hPromptPIDEff      ,"Neutral + 0.1<#sigma^{2}<0.3","PL");
  //legEff->AddEntry(hPromptIsoNoCutEff ,"Neutral cluster + isolated","PL");
  //legEff->AddEntry(hPromptIsoNoCutEff2,"Neutral cluster + isolated, denom. isolated","PL");
  legEff->AddEntry(hPromptIsoPIDEff   ,"Neutral + 0.1<#sigma^{2}<0.3 + isolated","PL");
  //legEff->AddEntry(hPromptIsoPIDEff2  ,"Neutral + 0.1<#sigma^{2}<0.3 + isolated, denom. isolated","PL");
  legEff->Draw();

  cEff->cd(2);
  
  //
  // Efficiency to be applied to data
  //
  hPromptIsoPIDEff2->Draw("same");
  hPromptIsoPIDEff2->SetMarkerStyle(20);
  hPromptIsoPIDEff2->SetMarkerColor(kBlack);
  hPromptIsoPIDEff2->SetLineColor(kBlack);
  hPromptIsoPIDEff2->SetAxisRange(10,100,"X");
  hPromptIsoPIDEff2->SetMaximum(1.00);
  hPromptIsoPIDEff2->SetMinimum(0.00);
  hPromptIsoPIDEff2->SetMarkerSize(markerSize);
  hPromptIsoPIDEff2->SetYTitle("selected reconstructed prompt /  generated and isolated prompt");  

  TLegend *legEff2 = new TLegend(0.12,0.12,0.5,0.4);
  legEff2->SetHeader("Final efficiency");
  legEff2->SetTextSize(0.045);
  legEff2->SetLineColor(0);
  legEff2->SetBorderSize(0);
  legEff2->AddEntry(hPromptIsoPIDEff2  ,"Neutral + 0.1<#sigma^{2}<0.3 + isolated, denom. isolated","PL");
  legEff2->Draw();

  cEff->Print("figures/PromptSpectraEfficiencies.eps");
 
}
