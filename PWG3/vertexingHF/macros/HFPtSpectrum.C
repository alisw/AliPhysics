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
#include "TROOT.h"

#include "AliHFSystErr.h"
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
		    Int_t option=1, Double_t lumi=1.0, Double_t effTrig=1.0,
		    const char *directsimufilename="", const char *directsimuhistoname="CFHFccontainer0_New_3Prong_SelStep0_proj-pt", 
		    const char *directrecofilename="", const char *directrecohistoname="CFHFccontainer0_New_3Prong_SelStep8_proj-pt", 
		    const char *feeddownsimufilename="", const char *feeddownsimuhistoname="CFHFccontainer0allD_New_3Prong_SelStep0_proj-pt", 
		    const char *feeddownrecofilename="", const char *feeddownrecohistoname="CFHFccontainer0allD_New_3Prong_SelStep8_proj-pt") {

  //  Set if calculation considers asymmetric uncertainties or not 
  bool asym = true;

  // Set the meson and decay
  // (only D0 -> K pi, D+--> K pi pi & D* --> D0 pi implemented here)
  bool isD0Kpi = true;
  bool isDplusKpipi = false;
  bool isDstarD0pi = false;
  if (isD0Kpi && isDplusKpipi && isDstarD0pi) {
    cout << "Sorry, can not deal with more than one correction at the same time"<<endl;
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
  TH1D *hDirectMCpt=0;           // Input MC c-->D spectra
  TH1D *hFeedDownMCpt=0;         // Input MC b-->D spectra
  TH1D *hDirectMCptMax=0;        // Input MC maximum c-->D spectra
  TH1D *hDirectMCptMin=0;        // Input MC minimum c-->D spectra
  TH1D *hFeedDownMCptMax=0;      // Input MC maximum b-->D spectra
  TH1D *hFeedDownMCptMin=0;      // Input MC minimum b-->D spectra
  TGraphAsymmErrors *gPrediction=0; // Input MC c-->D spectra
  TH1D *hDirectEffpt=0;          // c-->D Acceptance and efficiency correction
  TH1D *hFeedDownEffpt=0;        // b-->D Acceptance and efficiency correction
  TH1D *hRECpt=0;                // all reconstructed D
  //
  TH1D *hDirectSimulationpt=0;       // Simulated c--D spectra (used to re-compute the efficiency)
  TH1D *hDirectReconstructionpt=0;   // Reconstructed c--D spectra (used to re-compute the efficiency)
  TH1D *hFeedDownSimulationpt=0;     // Simulated b--D spectra (used to re-compute the efficiency)
  TH1D *hFeedDownReconstructionpt=0; // Reconstructed b--D spectra (used to re-compute the efficiency)
  //

  //
  // Define/Get the input histograms
  //
  Int_t decay=0;
  TFile * mcfile = new TFile(mcfilename,"read");
  if (isD0Kpi){
    decay = 1;
    hDirectMCpt = (TH1D*)mcfile->Get("hD0Kpipred_central");
    hFeedDownMCpt = (TH1D*)mcfile->Get("hD0KpifromBpred_central");
    hDirectMCptMax = (TH1D*)mcfile->Get("hD0Kpipred_max");
    hDirectMCptMin = (TH1D*)mcfile->Get("hD0Kpipred_min");
    hFeedDownMCptMax = (TH1D*)mcfile->Get("hD0KpifromBpred_max");
    hFeedDownMCptMin = (TH1D*)mcfile->Get("hD0KpifromBpred_min");
    gPrediction = (TGraphAsymmErrors*)mcfile->Get("D0Kpiprediction");
  }
  else if (isDplusKpipi){
    decay = 2;
    hDirectMCpt = (TH1D*)mcfile->Get("hDpluskpipipred_central");
    hFeedDownMCpt = (TH1D*)mcfile->Get("hDpluskpipifromBpred_central");
    hDirectMCptMax = (TH1D*)mcfile->Get("hDpluskpipipred_max");
    hDirectMCptMin = (TH1D*)mcfile->Get("hDpluskpipipred_min");
    hFeedDownMCptMax = (TH1D*)mcfile->Get("hDpluskpipifromBpred_max");
    hFeedDownMCptMin = (TH1D*)mcfile->Get("hDpluskpipifromBpred_min");
    gPrediction = (TGraphAsymmErrors*)mcfile->Get("Dpluskpipiprediction");
  }
  else if(isDstarD0pi){
    decay = 3;
    hDirectMCpt = (TH1D*)mcfile->Get("hDstarD0pipred_central");
    hFeedDownMCpt = (TH1D*)mcfile->Get("hDstarD0pifromBpred_central");
    hDirectMCptMax = (TH1D*)mcfile->Get("hDstarD0pipred_max");
    hDirectMCptMin = (TH1D*)mcfile->Get("hDstarD0pipred_min");
    hFeedDownMCptMax = (TH1D*)mcfile->Get("hDstarD0pifromBpred_max");
    hFeedDownMCptMin = (TH1D*)mcfile->Get("hDstarD0pifromBpred_min");
    gPrediction = (TGraphAsymmErrors*)mcfile->Get("DstarD0piprediction");
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
    hDirectEffpt = (TH1D*)efffile->Get("hEffD");
    hDirectEffpt->SetNameTitle("hDirectEffpt","direct acc x eff");
    hFeedDownEffpt = (TH1D*)efffile->Get("hEffB");
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
  TH1D *histofc=0;
  TH1D *histofcMax=0;
  TH1D *histofcMin=0;
  TH1D *histoYieldCorr=0;
  TH1D *histoYieldCorrMax=0;
  TH1D *histoYieldCorrMin=0;
  TH1D *histoSigmaCorr=0;
  TH1D *histoSigmaCorrMax=0;
  TH1D *histoSigmaCorrMin=0;
  //
  int nbins = hRECpt->GetNbinsX();
  TGraphAsymmErrors * gYieldCorr = 0;
  TGraphAsymmErrors * gSigmaCorr = 0;
  TGraphAsymmErrors * gFcExtreme = 0;
  TGraphAsymmErrors * gFcConservative = 0;
  TGraphAsymmErrors * gYieldCorrExtreme = 0;
  TGraphAsymmErrors * gSigmaCorrExtreme = 0;
  TGraphAsymmErrors * gYieldCorrConservative = 0;
  TGraphAsymmErrors * gSigmaCorrConservative = 0;


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
    spectra->EstimateAndSetDirectEfficiencyRecoBin(hDirectSimulationpt,hDirectReconstructionpt);
    spectra->EstimateAndSetFeedDownEfficiencyRecoBin(hFeedDownSimulationpt,hFeedDownReconstructionpt);
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
  double lumiUnc = 0.10*lumi; // 10% uncertainty on the luminosity
  spectra->SetLuminosity(lumi,lumiUnc);
  spectra->SetTriggerEfficiency(effTrig,0.);

  // Set the global uncertainties on the efficiencies (in percent)
  double globalEffUnc = 0.15; 
  double globalBCEffRatioUnc = 0.15;
  //  double globalEffUnc = 0.; 
  //  double globalBCEffRatioUnc = 0.;
  spectra->SetAccEffPercentageUncertainty(globalEffUnc,globalBCEffRatioUnc);

  // Do the calculations
  cout << " Doing the calculation... "<< endl;
  double deltaY = 1.0;
  double branchingRatioC = 1.0;
  double branchingRatioBintoFinalDecay = 1.0; // this is relative to the input theoretical prediction
  spectra->ComputeHFPtSpectrum(deltaY,branchingRatioC,branchingRatioBintoFinalDecay);
  cout << "   ended the calculation, getting the histograms back " << endl;

  // Set the systematics externally
  AliHFSystErr *systematics = new AliHFSystErr();
  systematics->Init(decay);
  bool combineFeedDown = true;
  spectra->ComputeSystUncertainties(systematics,combineFeedDown);


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
  histoYieldCorr->SetNameTitle("histoYieldCorr","corrected yield");
  histoYieldCorrMax->SetNameTitle("histoYieldCorrMax","max corrected yield");
  histoYieldCorrMin->SetNameTitle("histoYieldCorrMin","min corrected yield");
  histoSigmaCorr->SetNameTitle("histoSigmaCorr","corrected invariant cross-section");
  histoSigmaCorrMax->SetNameTitle("histoSigmaCorrMax","max corrected invariant cross-section");
  histoSigmaCorrMin->SetNameTitle("histoSigmaCorrMin","min corrected invariant cross-section");
  // the efficiencies
  if(!hDirectEffpt) hDirectEffpt = (TH1D*)spectra->GetDirectAccEffCorrection();
  if(!hFeedDownEffpt) hFeedDownEffpt = (TH1D*)spectra->GetFeedDownAccEffCorrection();

  // Get & Rename the TGraphs
  if (asym) {
    gSigmaCorr = spectra->GetCrossSectionFromYieldSpectrum();
    gYieldCorr = spectra->GetFeedDownCorrectedSpectrum(); 
    gSigmaCorrExtreme = spectra->GetCrossSectionFromYieldSpectrumExtreme();
    gYieldCorrExtreme = spectra->GetFeedDownCorrectedSpectrumExtreme(); 
    gSigmaCorrConservative = spectra->GetCrossSectionFromYieldSpectrumConservative();
    gYieldCorrConservative = spectra->GetFeedDownCorrectedSpectrumConservative(); 
  }

  // Get & Rename the TGraphs
  if (option==0 && asym){
    gYieldCorr->SetNameTitle("gYieldCorr","gYieldCorr (uncorr)");
    gSigmaCorr->SetNameTitle("gSigmaCorr","gSigmaCorr (uncorr)");
  }
  if (option==1){
    // fc histos
    histofc = (TH1D*)spectra->GetHistoFeedDownCorrectionFc();
    histofcMax = (TH1D*)spectra->GetHistoUpperLimitFeedDownCorrectionFc();
    histofcMin = (TH1D*)spectra->GetHistoLowerLimitFeedDownCorrectionFc();
    histofc->SetNameTitle("histofc","fc correction factor");
    histofcMax->SetNameTitle("histofcMax","max fc correction factor");
    histofcMin->SetNameTitle("histofcMin","min fc correction factor");
    if (asym) {
      gYieldCorr->SetNameTitle("gYieldCorr","gYieldCorr (by fc)");
      gSigmaCorr->SetNameTitle("gSigmaCorr","gSigmaCorr (by fc)");
      gFcExtreme = spectra->GetFeedDownCorrectionFcExtreme();
      gFcExtreme->SetNameTitle("gFcExtreme","gFcExtreme");
      gYieldCorrExtreme->SetNameTitle("gYieldCorrExtreme","Extreme gYieldCorr (by fc)");
      gSigmaCorrExtreme->SetNameTitle("gSigmaCorrExtreme","Extreme gSigmaCorr (by fc)");
      gFcConservative = spectra->GetFeedDownCorrectionFcConservative();
      gFcConservative->SetNameTitle("gFcConservative","gFcConservative");
      gYieldCorrConservative->SetNameTitle("gYieldCorrConservative","Conservative gYieldCorr (by fc)");
      gSigmaCorrConservative->SetNameTitle("gSigmaCorrConservative","Conservative gSigmaCorr (by fc)");
    }
  }
  if (option==2 && asym) {
    gYieldCorr->SetNameTitle("gYieldCorr","gYieldCorr (by Nb)");
    gSigmaCorr->SetNameTitle("gSigmaCorr","gSigmaCorr (by Nb)");
    gYieldCorrExtreme->SetNameTitle("gYieldCorrExtreme","Extreme gYieldCorr (by Nb)");
    gSigmaCorrExtreme->SetNameTitle("gSigmaCorrExtreme","Extreme gSigmaCorr (by Nb)");
    gYieldCorrConservative->SetNameTitle("gYieldCorrConservative","Conservative gYieldCorr (by Nb)");
    gSigmaCorrConservative->SetNameTitle("gSigmaCorrConservative","Conservative gSigmaCorr (by Nb)");
    gFcConservative = spectra->GetFeedDownCorrectionFcConservative();
    gFcConservative->SetNameTitle("gFcConservative","gFcConservative");
  }

  //
  // Now, plot the results ! :)
  //

  cout << " Drawing the results ! " << endl;

  // control plots
  if (option==1) {

    TCanvas *ceff = new TCanvas("ceff","efficiency drawing");
    ceff->Divide(1,2);
    ceff->cd(1);
    hDirectEffpt->Draw();
    ceff->cd(2);
    hFeedDownEffpt->Draw();
    ceff->Update();

    TCanvas *cTheoryRebin = new TCanvas("cTheoryRebin","control the theoretical spectra rebin");
    cTheoryRebin->Divide(1,2);
    cTheoryRebin->cd(1);
    hDirectMCpt->Draw("");
    TH1D *hDirectMCptRebin = (TH1D*)spectra->GetDirectTheoreticalSpectrum();
    hDirectMCptRebin->SetLineColor(2);
    hDirectMCptRebin->Draw("same");
    cTheoryRebin->cd(2);
    hFeedDownMCpt->Draw("");
    TH1D *hFeedDownRebin = (TH1D*)spectra->GetFeedDownTheoreticalSpectrum();
    hFeedDownRebin->SetLineColor(2);
    hFeedDownRebin->Draw("same");
    cTheoryRebin->Update();
    
    TCanvas *cTheoryRebinLimits = new TCanvas("cTheoryRebinLimits","control the theoretical spectra limits rebin");
    cTheoryRebinLimits->Divide(1,2);
    cTheoryRebinLimits->cd(1);
    hDirectMCptMax->Draw("");
    TH1D *hDirectMCptMaxRebin = (TH1D*)spectra->GetDirectTheoreticalUpperLimitSpectrum();
    hDirectMCptMaxRebin->SetLineColor(2);
    hDirectMCptMaxRebin->Draw("same");
    hDirectMCptMin->Draw("same");
    TH1D *hDirectMCptMinRebin = (TH1D*)spectra->GetDirectTheoreticalLowerLimitSpectrum();
    hDirectMCptMinRebin->SetLineColor(2);
    hDirectMCptMinRebin->Draw("same");
    cTheoryRebinLimits->cd(2);
    hFeedDownMCptMax->Draw("");
    TH1D *hFeedDownMaxRebin = (TH1D*)spectra->GetFeedDownTheoreticalUpperLimitSpectrum();
    hFeedDownMaxRebin->SetLineColor(2);
    hFeedDownMaxRebin->Draw("same");
    hFeedDownMCptMin->Draw("same");
    TH1D *hFeedDownMinRebin = (TH1D*)spectra->GetFeedDownTheoreticalLowerLimitSpectrum();
    hFeedDownMinRebin->SetLineColor(2);
    hFeedDownMinRebin->Draw("same");
    cTheoryRebinLimits->Update();
  }
  
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

      if (gFcExtreme){

// 	for(Int_t item=0; item<gSigmaCorr->GetN(); item++){
// 	  Double_t center=0., value=0.;
// 	  gFcExtreme->GetPoint(item,center,value);
// 	  Double_t highunc = gFcExtreme->GetErrorYhigh(item) / value ;
// 	  Double_t lowunc = gFcExtreme->GetErrorYlow(item) / value ;
// 	  cout << "Fc extreme: i=" << item << ", center=" << center <<", value=" << value << " high unc=" << highunc*100 << "%, low unc=" << lowunc*100 << "%"<<endl;
// 	}
// 	for(Int_t item=0; item<gSigmaCorr->GetN(); item++){
// 	  Double_t center=0., value=0.;
// 	  gFcConservative->GetPoint(item,center,value);
// 	  Double_t highunc = gFcConservative->GetErrorYhigh(item) / value ;
// 	  Double_t lowunc = gFcConservative->GetErrorYlow(item) / value ;
// 	  cout << "Fc conservative: i=" << item << ", center=" << center <<", value=" << value << " high unc=" << highunc*100 << "%, low unc=" << lowunc*100 << "%"<<endl;
// 	}
	TCanvas *cfcExtreme = new TCanvas("cfcExtreme","Extreme Asymmetric fc (TGraphAsymmErr)");
	gFcExtreme->SetFillStyle(3006);
	gFcExtreme->SetLineWidth(3);
	gFcExtreme->SetMarkerStyle(20);
	gFcExtreme->SetFillColor(2);
	histofcDraw->Draw();
	gFcExtreme->Draw("3same");

	if(gFcConservative){
	  gFcConservative->SetFillStyle(3007);
	  gFcConservative->SetFillColor(4);
	  gFcConservative->Draw("3same");
	}

	cfcExtreme->Update();
      }
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

    TCanvas * cyieldExtreme = new TCanvas("cyieldExtreme","Extreme Asymmetric corrected yield (TGraphAsymmErr)");
    histoYieldCorr->Draw();
    gYieldCorrExtreme->SetFillStyle(3002);
    gYieldCorrExtreme->SetLineWidth(3);
    gYieldCorrExtreme->SetMarkerStyle(20);
    gYieldCorrExtreme->SetFillColor(2);
    histoYieldCorr->Draw();
    gYieldCorr->Draw("3same");
    gYieldCorrExtreme->Draw("3same");
    cyieldExtreme->SetLogy();
    cyieldExtreme->Update();
    
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
    gSigmaCorr->SetFillColor(3);
    histo2Draw->Draw();
    gSigmaCorr->Draw("3LPsame");
    gSigmaCorr->Draw("Xsame");
    csigmaAsym->SetLogy();
    csigmaAsym->Update();

//     cout << endl <<" Sytematics (stat approach) " <<endl;
//     for(Int_t item=0; item<gSigmaCorr->GetN(); item++){
//       Double_t center=0., value=0.;
//       gSigmaCorr->GetPoint(item,center,value);
//       Double_t highunc = gSigmaCorr->GetErrorYhigh(item) / value ;
//       Double_t lowunc = gSigmaCorr->GetErrorYlow(item) / value ;
//       cout << "Sigma syst (stat), i=" << item << ", center=" << center <<", value=" << value << " high unc=" << highunc*100 << "%, low unc=" << lowunc*100 << "%"<<endl;
//     }

    TCanvas * csigmaExtreme = new TCanvas("csigmaExtreme","Asymmetric extreme corrected sigma (TGraphAsymmErr)");
    histoSigmaCorr->Draw();
    gSigmaCorr->Draw("3Psame");
    gSigmaCorrExtreme->SetFillStyle(3002);
    gSigmaCorrExtreme->SetLineWidth(3);
    gSigmaCorrExtreme->SetMarkerStyle(21);
    gSigmaCorrExtreme->SetFillColor(2);
    gSigmaCorrExtreme->Draw("3Psame");
    csigmaExtreme->SetLogy();
    csigmaExtreme->Update();
      
//     cout << endl << " Sytematics (Extreme approach)" <<endl;
//     for(Int_t item=0; item<gSigmaCorrExtreme->GetN(); item++){
//       Double_t center=0., value=0.;
//       gSigmaCorrExtreme->GetPoint(item,center,value);
//       Double_t highunc = gSigmaCorrExtreme->GetErrorYhigh(item) / value ;
//       Double_t lowunc = gSigmaCorrExtreme->GetErrorYlow(item) / value ;
//       cout << "Sigma syst (extreme) i=" << item << ", center=" << center <<", value=" << value << " high unc=" << highunc*100 << "%, low unc=" << lowunc*100 << "%"<<endl;
//     }
    
//     cout << endl << " Sytematics (Conservative approach)" <<endl;
//     for(Int_t item=0; item<gSigmaCorrConservative->GetN(); item++){
//       Double_t center=0., value=0.;
//       gSigmaCorrConservative->GetPoint(item,center,value);
//       Double_t highunc = gSigmaCorrConservative->GetErrorYhigh(item) / value ;
//       Double_t lowunc = gSigmaCorrConservative->GetErrorYlow(item) / value ;
//       cout << "Sigma syst (conservative) i=" << item << ", center=" << center <<", value=" << value << " high unc=" << highunc*100 << "%, low unc=" << lowunc*100 << "%"<<endl;
//     }
    
 }
 


  //
  // Write the histograms to the output file
  //
  cout << " Saving the results ! " << endl<< endl;

  out->cd();
  //
  hDirectMCpt->Write();        hFeedDownMCpt->Write();
  hDirectMCptMax->Write();     hDirectMCptMin->Write();
  hFeedDownMCptMax->Write();   hFeedDownMCptMin->Write();
  if(hDirectEffpt) hDirectEffpt->Write();        if(hFeedDownEffpt) hFeedDownEffpt->Write();
  hRECpt->Write();
  //
  histoYieldCorr->Write();
  histoYieldCorrMax->Write();     histoYieldCorrMin->Write();   
  histoSigmaCorr->Write();
  histoSigmaCorrMax->Write();     histoSigmaCorrMin->Write();

  if(asym){
    gYieldCorr->Write();
    gSigmaCorr->Write();
    if(gYieldCorrExtreme) gYieldCorrExtreme->Write();
    if(gSigmaCorrExtreme) gSigmaCorrExtreme->Write();
    if(gYieldCorrConservative) gYieldCorrConservative->Write();
    if(gSigmaCorrConservative) gSigmaCorrConservative->Write();
    if(asym && gFcConservative) gFcConservative->Write();
  }

  if(option==1){
    histofc->Write();
    histofcMax->Write();     histofcMin->Write(); 
    if(asym && gFcExtreme) gFcExtreme->Write();
  }


  // Draw the cross-section 
  //  spectra->DrawSpectrum(gPrediction);

  //  out->Close();

}
