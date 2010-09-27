//
// Macro to use the AliHFPtSpectrum class
//  to compute the feed-down corrections for heavy-flavor
//
//  Z.Conesa, September 2010 (zconesa@in2p3.fr)
//

#include <Riostream.h>
#include "TH1D.h"
#include "TH1.h"
#include "TH2F.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"

#include "AliHFPtSpectrum.h"

//
// Macro execution parameters: 
//  0) filename with the theoretical predictions  (direct & feed-down)
//  1) acceptance and reconstruction efficiencies file name (direct & feed-down)
//  2) reconstructed spectra file name 
//  3) output file name
//  4) Set the feed-down calculation option flag: 0=none, 1=fc only, 2=Nb only
//  5) Set the luminosity
//  6) Set the trigger efficiency
//  7-14) If the efficiency histos do not have the right bin width, set the files & histo-names, they'll be computed, if the efficiencies are in file (6), don't set this parameters
//
void HFPtSpectrum ( const char *mcfilename="FeedDownCorrectionMC.root",
		    const char *efffilename="Efficiencies.root",
		    const char *recofilename="Reconstructed.root", const char *recohistoname="hRawSpectrumD0",
		    const char *outfilename="HFPtSpectrum.root",
		    int option=1, double lumi=1.0, double effTrig=1.0,
		    const char *directsimufilename="", const char *directsimuhistoname="CFHFccontainer0_New_3Prong_SelStep0_proj-pt", 
		    const char *directrecofilename="", const char *directrecohistoname="CFHFccontainer0_New_3Prong_SelStep8_proj-pt", 
		    const char *feeddownsimufilename="", const char *feeddownsimuhistoname="CFHFccontainer0allD_New_3Prong_SelStep0_proj-pt", 
		    const char *feeddownrecofilename="", const char *feeddownrecohistoname="CFHFccontainer0allD_New_3Prong_SelStep8_proj-pt") {

  //  Set if calculation considers asymmetric uncertainties or not 
  bool asym = true;

  // Set the meson and decay
  // (only D0 -> K pi & D+-->K pi pi implemented here)
  bool isD0Kpi = true;
  bool isDplusKpipi = false;
  if (isD0Kpi && isDplusKpipi) {
    cout << "Sorry, can not deal with the two corrections at the same time"<<endl;
    return;
  }

  if (option>2) { 
    cout<< "Bad calculation option, should be <=2"<<endl;
    return;
  }
  if (option==0) asym = false;

  //
  // Get the histograms from the files
  //
  TH1D *hDirectMCpt;            // Input MC c-->D spectra
  TH1D *hFeedDownMCpt;          // Input MC b-->D spectra
  TH1D *hDirectMCptMax;        // Input MC maximum c-->D spectra
  TH1D *hDirectMCptMin;        // Input MC minimum c-->D spectra
  TH1D *hFeedDownMCptMax;      // Input MC maximum b-->D spectra
  TH1D *hFeedDownMCptMin;      // Input MC minimum b-->D spectra
  TH1D *hDirectEffpt;           // c-->D Acceptance and efficiency correction
  TH1D *hFeedDownEffpt;         // b-->D Acceptance and efficiency correction
  TH1D *hRECpt;                 // all reconstructed D
  //
  TH1D *hDirectSimulationpt;       // Simulated c--D spectra (used to re-compute the efficiency)
  TH1D *hDirectReconstructionpt;   // Reconstructed c--D spectra (used to re-compute the efficiency)
  TH1D *hFeedDownSimulationpt;     // Simulated b--D spectra (used to re-compute the efficiency)
  TH1D *hFeedDownReconstructionpt; // Reconstructed b--D spectra (used to re-compute the efficiency)
  //

  //
  // Define/Get the input histograms
  //
  TFile * mcfile = new TFile(mcfilename,"read");
  if (isD0Kpi){
    hDirectMCpt = (TH1D*)mcfile->Get("hD0Kpipred_central");
    hFeedDownMCpt = (TH1D*)mcfile->Get("hD0KpifromBpred_central");
    hDirectMCptMax = (TH1D*)mcfile->Get("hD0Kpipred_max");
    hDirectMCptMin = (TH1D*)mcfile->Get("hD0Kpipred_min");
    hFeedDownMCptMax = (TH1D*)mcfile->Get("hD0KpifromBpred_max");
    hFeedDownMCptMin = (TH1D*)mcfile->Get("hD0KpifromBpred_min");
  }
  else if (isDplusKpipi){
    hDirectMCpt = (TH1D*)mcfile->Get("hDpluskpipipred_central");
    hFeedDownMCpt = (TH1D*)mcfile->Get("hDpluskpipifromBpred_central");
    hDirectMCptMax = (TH1D*)mcfile->Get("hDpluskpipipred_max");
    hDirectMCptMin = (TH1D*)mcfile->Get("hDpluskpipipred_min");
    hFeedDownMCptMax = (TH1D*)mcfile->Get("hDpluskpipifromBpred_max");
    hFeedDownMCptMin = (TH1D*)mcfile->Get("hDpluskpipifromBpred_min");
  }
  //
  hDirectMCpt->SetNameTitle("hDirectMCpt","direct MC spectra");
  hFeedDownMCpt->SetNameTitle("hFeedDownMCpt","feed-down MC spectra");
  hDirectMCptMax->SetNameTitle("hDirectMCptMax","max direct MC spectra");
  hDirectMCptMin->SetNameTitle("hDirectMCptMin","min direct MC spectra");
  hFeedDownMCptMax->SetNameTitle("hFeedDownMCptMax","max feed-down MC spectra");
  hFeedDownMCptMin->SetNameTitle("hFeedDownMCptMin","min feed-down MC spectra");
  //
  //
  if (strcmp(directsimufilename,"")!=0 && strcmp(directrecofilename,"")!=0 &&
      strcmp(feeddownsimufilename,"")!=0 && strcmp(feeddownrecofilename,"")!=0 ) {
    if (strcmp(directsimufilename,"")!=0){
      TFile *directSimufile = new TFile(directsimufilename,"read");
      hDirectSimulationpt = (TH1D*)directSimufile->Get(directsimuhistoname);
      hDirectSimulationpt->SetNameTitle("hDirectSimulationpt","hDirectSimulationpt");
    }
    if (strcmp(directrecofilename,"")!=0){
      TFile *directRecofile = new TFile(directrecofilename,"read");
      hDirectReconstructionpt = (TH1D*)directRecofile->Get(directrecohistoname);
      hDirectReconstructionpt->SetNameTitle("hDirectReconstructionpt","hDirectReconstructionpt");
    }
    if (strcmp(feeddownsimufilename,"")!=0){
      TFile *feeddownSimufile = new TFile(feeddownsimufilename,"read");
      hFeedDownSimulationpt = (TH1D*)feeddownSimufile->Get(feeddownsimuhistoname);
      hFeedDownSimulationpt->SetNameTitle("hFeedDownSimulationpt","hFeedDownSimulationpt");
    }
    if (strcmp(feeddownrecofilename,"")!=0){
      TFile *feeddownRecofile = new TFile(feeddownrecofilename,"read");
      hFeedDownReconstructionpt = (TH1D*)feeddownRecofile->Get(feeddownrecohistoname);
      hFeedDownReconstructionpt->SetNameTitle("hFeedDownReconstructionpt","hFeedDownReconstructionpt");
    }
  }
  else {
    TFile * efffile = new TFile(efffilename,"read");
    hDirectEffpt = (TH1D*)efffile->Get("hDirectEffpt");
    hDirectEffpt->SetNameTitle("hDirectEffpt","direct acc x eff");
    hFeedDownEffpt = (TH1D*)efffile->Get("hFeedDownEffpt");
    hFeedDownEffpt->SetNameTitle("hFeedDownEffpt","feed-down acc x eff");
  }
  //
  //
  TFile * recofile = new TFile(recofilename,"read");
  hRECpt = (TH1D*)recofile->Get(recohistoname);
  hRECpt->SetNameTitle("hRECpt","Reconstructed spectra");

  //
  // Define the output histograms
  //
  TFile *out = new TFile(outfilename,"recreate");
  //
  //
  TH1D *histofc = (TH1D*)hRECpt->Clone();
  histofc->SetNameTitle("histofc","fc correction factor");
  histofc->Reset();
  TH1D *histofcMax = (TH1D*)histofc->Clone();
  TH1D *histofcMin = (TH1D*)histofc->Clone();
  TH1D *histoYieldCorr = (TH1D*)hRECpt->Clone();
  histoYieldCorr->SetNameTitle("histoYieldCorrFc","corrected yield");
  histoYieldCorr->Reset();
  TH1D *histoYieldCorrMax = (TH1D*)histoYieldCorr->Clone();
  histoYieldCorrMax->SetNameTitle("histoYieldCorrMax","max corrected yield");
  TH1D *histoYieldCorrMin = (TH1D*)histoYieldCorr->Clone();
  histoYieldCorrMin->SetNameTitle("histoYieldCorrMin","min corrected yield");
  TH1D *histoSigmaCorr = (TH1D*)hRECpt->Clone();
  histoSigmaCorr->SetNameTitle("histoSigmaCorr","corrected invariant cross-section");
  histoSigmaCorr->Reset();
  TH1D *histoSigmaCorrMax = (TH1D*)histoSigmaCorr->Clone();
  histoSigmaCorrMax->SetNameTitle("histoSigmaCorrMax","max corrected invariant cross-section");
  TH1D *histoSigmaCorrMin = (TH1D*)histoSigmaCorr->Clone();
  histoSigmaCorrMin->SetNameTitle("histoSigmaCorrMin","min corrected invariant cross-section");
  //
  int nbins = hRECpt->GetNbinsX();
  TGraphAsymmErrors * gFc = new TGraphAsymmErrors(nbins);
  TGraphAsymmErrors * gYieldCorr = new TGraphAsymmErrors(nbins);
  TGraphAsymmErrors * gSigmaCorr = new TGraphAsymmErrors(nbins);


  //
  // Main functionalities for the calculation
  //

  // Define and set the basic option flags
  AliHFPtSpectrum * spectra = new AliHFPtSpectrum("AliHFPtSpectrum","AliHFPtSpectrum",option);
  spectra->SetFeedDownCalculationOption(option);
  spectra->SetComputeAsymmetricUncertainties(asym);

  // Feed the input histograms
  //  reconstructed spectra
  cout << " Setting the reconstructed spectrum,";
  spectra->SetReconstructedSpectrum(hRECpt);
  // acceptance and efficiency corrections
  cout << " the files to compute the efficiency,";
  if (strcmp(directsimufilename,"")!=0 && strcmp(directrecofilename,"")!=0 &&
      strcmp(feeddownsimufilename,"")!=0 && strcmp(feeddownrecofilename,"")!=0 ) {
    hDirectEffpt = (TH1D*)spectra->EstimateAndSetDirectEfficiencyRecoBin(hDirectSimulationpt,hDirectReconstructionpt);
    hFeedDownEffpt = (TH1D*)spectra->EstimateAndSetFeedDownEfficiencyRecoBin(hFeedDownSimulationpt,hFeedDownReconstructionpt);
    TCanvas *ceff = new TCanvas("ceff","efficiency drawing");
    ceff->Divide(1,2);
    ceff->cd(1);
    hDirectEffpt->Draw();
    ceff->cd(2);
    hFeedDownEffpt->Draw();
    ceff->Update();
  }
  else {
    spectra->SetAccEffCorrection(hDirectEffpt,hFeedDownEffpt);
  }
  // option specific histos (theory)
  cout << " the theoretical spectra";
  if(option==1){
    spectra->SetMCptSpectra(hDirectMCpt,hFeedDownMCpt);
    if(asym) spectra->SetMCptDistributionsBounds(hDirectMCptMax,hDirectMCptMin,hFeedDownMCptMax,hFeedDownMCptMin);
  }
  else if(option==2){
    spectra->SetFeedDownMCptSpectra(hFeedDownMCpt);
    if(asym) spectra->SetFeedDownMCptDistributionsBounds(hFeedDownMCptMax,hFeedDownMCptMin);
  }

  cout << " and the normalization" <<endl;
  // Set normalization factors (uncertainties set to 0. as example)
  spectra->SetLuminosity(lumi,0.);
  spectra->SetTriggerEfficiency(effTrig,0.);

  // Do the calculations
  cout << " Doing the calculation... "<< endl;
  double deltaY = 1.0;
  double branchingRatioC = 1.0;
  double branchingRatioBintoFinalDecay = 1.0; // this is relative to the input theoretical prediction
  spectra->ComputeHFPtSpectrum(deltaY,branchingRatioC,branchingRatioBintoFinalDecay);
  cout << "   ended the calculation, getting the histograms back " << endl;

  //
  // Get the output histograms
  //
  // the corrected yield and cross-section
  histoYieldCorr = (TH1D*)spectra->GetHistoFeedDownCorrectedSpectrum();
  histoSigmaCorr = (TH1D*)spectra->GetHistoCrossSectionFromYieldSpectrum();
  histoYieldCorrMax = (TH1D*)spectra->GetHistoUpperLimitFeedDownCorrectedSpectrum(); 
  histoYieldCorrMin = (TH1D*)spectra->GetHistoLowerLimitFeedDownCorrectedSpectrum(); 
  histoSigmaCorrMax = (TH1D*)spectra->GetHistoUpperLimitCrossSectionFromYieldSpectrum();
  histoSigmaCorrMin = (TH1D*)spectra->GetHistoLowerLimitCrossSectionFromYieldSpectrum();
  if (asym) {
    gSigmaCorr = spectra->GetCrossSectionFromYieldSpectrum();
    gYieldCorr = spectra->GetFeedDownCorrectedSpectrum(); 
  }

  if (option==0 && asym){
    gYieldCorr->SetNameTitle("gYieldCorr","gYieldCorr (uncorr)");
    gSigmaCorr->SetNameTitle("gSigmaCorr","gSigmaCorr (uncorr)");
  }
  if (option==1){
    // fc histos
    histofc = (TH1D*)spectra->GetHistoFeedDownCorrectionFc();
    histofcMax = (TH1D*)spectra->GetHistoUpperLimitFeedDownCorrectionFc();
    histofcMin = (TH1D*)spectra->GetHistoLowerLimitFeedDownCorrectionFc();
    histofcMax->SetNameTitle("hfcMax","max fc correction factor");
    histofcMin->SetNameTitle("histofcMin","min fc correction factor");
    if (asym) {
      gFc = spectra->GetFeedDownCorrectionFc();
      gFc->SetNameTitle("gFc","gFc");
      gYieldCorr->SetNameTitle("gYieldCorr","gYieldCorr (by fc)");
      gSigmaCorr->SetNameTitle("gSigmaCorr","gSigmaCorr (by fc)");
    }
  }
  if (option==2 && asym) {
    gYieldCorr->SetNameTitle("gYieldCorr","gYieldCorr (by Nb)");
    gSigmaCorr->SetNameTitle("gSigmaCorr","gSigmaCorr (by Nb)");
  }

  //
  // Now, plot the results ! :)
  //

  cout << " Drawing the results ! " << endl;

  // control plots
  TCanvas *cTheoryRebin = new TCanvas("cTheoryRebin","control the theoretical spectra rebin");
  cTheoryRebin->Divide(1,2);
  cTheoryRebin->cd(1);
  hDirectMCpt->Draw("");
  TH1D *hDirectMCptRebin = (TH1D*)spectra->GetDirectTheoreticalSpectrum();
  hDirectMCptRebin->SetNameTitle("hDirectMCptRebin","hDirectMCptRebin");
  hDirectMCptRebin->SetLineColor(2);
  hDirectMCptRebin->Draw("same");
  cTheoryRebin->cd(2);
  hFeedDownMCpt->Draw("");
  TH1D *hFeedDownRebin = (TH1D*)spectra->GetFeedDownTheoreticalSpectrum();
  hFeedDownRebin->SetNameTitle("hFeedDownRebin","hFeedDownRebin");
  hFeedDownRebin->SetLineColor(2);
  hFeedDownRebin->Draw("same");
  cTheoryRebin->Update();

  TCanvas *cTheoryRebinLimits = new TCanvas("cTheoryRebinLimits","control the theoretical spectra limits rebin");
  cTheoryRebinLimits->Divide(1,2);
  cTheoryRebinLimits->cd(1);
  hDirectMCptMax->Draw("");
  TH1D *hDirectMCptMaxRebin = (TH1D*)spectra->GetDirectTheoreticalUpperLimitSpectrum();
  hDirectMCptMaxRebin->SetNameTitle("hDirectMCptMaxRebin","hDirectMCptMaxRebin");
  hDirectMCptMaxRebin->SetLineColor(2);
  hDirectMCptMaxRebin->Draw("same");
  hDirectMCptMin->Draw("same");
  TH1D *hDirectMCptMinRebin = (TH1D*)spectra->GetDirectTheoreticalLowerLimitSpectrum();
  hDirectMCptMinRebin->SetNameTitle("hDirectMCptMinRebin","hDirectMCptMinRebin");
  hDirectMCptMinRebin->SetLineColor(2);
  hDirectMCptMinRebin->Draw("same");
  cTheoryRebinLimits->cd(2);
  hFeedDownMCptMax->Draw("");
  TH1D *hFeedDownMaxRebin = (TH1D*)spectra->GetFeedDownTheoreticalUpperLimitSpectrum();
  hFeedDownMaxRebin->SetNameTitle("hFeedDownMaxRebin","hFeedDownMaxRebin");
  hFeedDownMaxRebin->SetLineColor(2);
  hFeedDownMaxRebin->Draw("same");
  hFeedDownMCptMin->Draw("same");
  TH1D *hFeedDownMinRebin = (TH1D*)spectra->GetFeedDownTheoreticalLowerLimitSpectrum();
  hFeedDownMinRebin->SetNameTitle("hFeedDownMinRebin","hFeedDownMinRebin");
  hFeedDownMinRebin->SetLineColor(2);
  hFeedDownMinRebin->Draw("same");
  cTheoryRebinLimits->Update();
  
  if (option==1) {

    TCanvas * cfc = new TCanvas("cfc","Fc");
    histofcMax->Draw("c");
    histofc->Draw("csame");
    histofcMin->Draw("csame");
    cfc->Update();

    if (asym) {
      TH2F *histofcDraw= new TH2F("histofcDraw","histofc (for drawing)",100,0,33.25,100,0.01,1.25);
      histofcDraw->SetStats(0);
      histofcDraw->GetXaxis()->SetTitle("p_{T}  [GeV]");
      histofcDraw ->GetXaxis()->SetTitleSize(0.05);
      histofcDraw->GetXaxis()->SetTitleOffset(0.95);
      histofcDraw->GetYaxis()->SetTitle(" fc ");
      histofcDraw->GetYaxis()->SetTitleSize(0.05);
      TCanvas *cfcAsym = new TCanvas("cfcAsym","Asymmetric fc (TGraphAsymmErr)");
      gFc->SetFillStyle(3001);
      gFc->SetLineWidth(3);
      gFc->SetMarkerStyle(20);
      gFc->SetFillColor(3);
      histofcDraw->Draw();
      gFc->Draw("3Csame");
      gFc->Draw("Xsame");
      cfcAsym->Update();
    }

  }

  //
  // Drawing the results (the raw-reconstructed, the expected, and the corrected spectra)
  //
  TCanvas * cresult = new TCanvas("cresult","corrected yields & sigma");
  hDirectMCpt->SetMarkerStyle(20);
  hDirectMCpt->SetMarkerColor(4);
  hDirectMCpt->Draw("p");
  histoSigmaCorr->SetMarkerStyle(21);
  histoSigmaCorr->SetMarkerColor(2);
  histoSigmaCorr->Draw("psame");
  histoYieldCorr->SetMarkerStyle(22);
  histoYieldCorr->SetMarkerColor(6);
  histoYieldCorr->Draw("psame");
  hRECpt->SetMarkerStyle(23);
  hRECpt->SetMarkerColor(3);
  hRECpt->Draw("psame");
  cresult->SetLogy();
  cresult->Update();

  TCanvas * cresult2 = new TCanvas("cresult2","corrected yield & sigma");
  histoSigmaCorr->SetMarkerStyle(21);
  histoSigmaCorr->SetMarkerColor(2);
  histoSigmaCorr->Draw("p");
  histoYieldCorr->SetMarkerStyle(22);
  histoYieldCorr->SetMarkerColor(6);
  histoYieldCorr->Draw("psame");
  hRECpt->SetMarkerStyle(23);
  hRECpt->SetMarkerColor(3);
  hRECpt->Draw("psame");
  cresult2->SetLogy();
  cresult2->Update();


  if (asym) { 

      TH2F *histoDraw = new TH2F("histoDraw","histo (for drawing)",100,0,33.25,100,50.,1e7);
      float max = 1.1*gYieldCorr->GetMaximum();
      histoDraw->SetAxisRange(0.1,max,"Y");
      histoDraw->SetStats(0);
      histoDraw->GetXaxis()->SetTitle("p_{T}  [GeV]");
      histoDraw->GetXaxis()->SetTitleSize(0.05);
      histoDraw->GetXaxis()->SetTitleOffset(0.95);
      histoDraw->GetYaxis()->SetTitle("#frac{d#N}{dp_{T}} |_{|y|<1} [L & trigger uncorr]");
      histoDraw->GetYaxis()->SetTitleSize(0.05);
      TCanvas * cyieldAsym = new TCanvas("cyieldAsym","Asymmetric corrected yield (TGraphAsymmErr)");
      gYieldCorr->SetFillStyle(3001);
      gYieldCorr->SetLineWidth(3);
      gYieldCorr->SetMarkerStyle(20);
      gYieldCorr->SetFillColor(3);
      histoDraw->Draw();
      gYieldCorr->Draw("3LPsame");
      gYieldCorr->Draw("Xsame");
      cyieldAsym->SetLogy();
      cyieldAsym->Update();

      TH2F *histo2Draw = new TH2F("histo2Draw","histo2 (for drawing)",100,0,33.25,100,50.,1e9);
      max = 1.1*gSigmaCorr->GetMaximum();
      histo2Draw->SetAxisRange(0.1,max,"Y");
      histo2Draw->SetStats(0);
      histo2Draw->GetXaxis()->SetTitle("p_{T}  [GeV]");
      histo2Draw->GetXaxis()->SetTitleSize(0.05);
      histo2Draw->GetXaxis()->SetTitleOffset(0.95);
      histo2Draw->GetYaxis()->SetTitle("#frac{1}{BR} #times #frac{d#sigma}{dp_{T}} |_{|y|<1}");
      histo2Draw->GetYaxis()->SetTitleSize(0.05);
      TCanvas * csigmaAsym = new TCanvas("csigmaAsym","Asymmetric corrected sigma (TGraphAsymmErr)");
      gSigmaCorr->SetFillStyle(3001);
      gSigmaCorr->SetLineWidth(3);
      gSigmaCorr->SetMarkerStyle(21);
      gSigmaCorr->SetFillColor(4);
      histo2Draw->Draw();
      gSigmaCorr->Draw("3LPsame");
      gSigmaCorr->Draw("Xsame");
      csigmaAsym->SetLogy();
      csigmaAsym->Update();
  }
 


  //
  // Write the histograms to the output file
  //
  cout << " Saving the results ! " << endl<< endl;

  out->cd();
  //
  hDirectMCpt->Write();         hFeedDownMCpt->Write();
  hDirectMCptMax->Write();     hDirectMCptMin->Write();
  hFeedDownMCptMax->Write();   hFeedDownMCptMin->Write();
  hDirectEffpt->Write();        hFeedDownEffpt->Write();
  hRECpt->Write();
  //
  histoYieldCorr->Write();
  histoYieldCorrMax->Write();     histoYieldCorrMin->Write();   
  histoSigmaCorr->Write();
  histoSigmaCorrMax->Write();     histoSigmaCorrMin->Write();

  if(asym){
    gYieldCorr->Write();
    gSigmaCorr->Write();
  }

  if(option==1){
    histofc->Write();
    histofcMax->Write();     histofcMin->Write();   
    if(asym) gFc->Write();
  }

  //  out->Close();

}
