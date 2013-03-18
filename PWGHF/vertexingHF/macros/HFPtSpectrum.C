#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include "TH1D.h"
#include "TH1.h"
#include "TH2F.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLegend.h"
#include "AliHFSystErr.h"
#include "AliHFPtSpectrum.h"
#endif

/* $Id$ */ 

//
// Macro to use the AliHFPtSpectrum class
//  to compute the feed-down corrections for heavy-flavor
//
//  Z.Conesa, September 2010 (zconesa@in2p3.fr)
//



//
// Macro execution parameters: 
//  0) filename with the theoretical predictions  (direct & feed-down)
//  1) acceptance and reconstruction efficiencies file name (direct & feed-down)
//  2) reconstructed spectra file name 
//  3) output file name
//  4) Set the feed-down calculation option flag: knone=none, kfc=fc only, kNb=Nb only
//  5-6) Set the luminosity: the number of events analyzed, and the cross-section of the sample [pb]
//  7) Set whether the yield is for particle + anti-particles or only one of the 'charges'
//  8) Set the centrality class
//  9) Flag to decide if there is need to evaluate the dependence on the energy loss
//

enum centrality{ kpp7, kpp276, k07half, k010, k1020, k020, k2040, k2030, k3040, k4050, k3050, k5060, k4060, k6080, k4080, k80100 };
enum BFDSubtrMethod { knone, kfc, kNb };
enum RaavsEP {kPhiIntegrated, kInPlane, kOutOfPlane};

void HFPtSpectrum ( const char *mcfilename="FeedDownCorrectionMC.root",
		    const char *efffilename="Efficiencies.root",
		    const char *recofilename="Reconstructed.root", const char *recohistoname="hRawSpectrumD0",
		    const char *outfilename="HFPtSpectrum.root",
		    Int_t fdMethod=kNb, Double_t nevents=1.0, Double_t sigma=1.0, // sigma[pb]
		    Bool_t isParticlePlusAntiParticleYield=true, Int_t cc=kpp7, Bool_t PbPbEloss=false, 
		    Int_t isRaavsEP=kPhiIntegrated,const char *epResolfile="") {


  gROOT->Macro("$ALICE_ROOT/PWGHF/vertexingHF/macros/LoadLibraries.C");

  //  Set if calculation considers asymmetric uncertainties or not 
  Bool_t asym = true;

  // Set the meson and decay
  // (only D0 -> K pi, D+--> K pi pi & D* --> D0 pi & D+s -->KKpi implemented here)
  Bool_t isD0Kpi = true;
  Bool_t isDplusKpipi = false;
  Bool_t isDstarD0pi = false;
  Bool_t isDsKKpi = false;
  if (isD0Kpi && isDplusKpipi && isDstarD0pi && isDsKKpi) {
    cout << "Sorry, can not deal with more than one correction at the same time"<<endl;
    return;
  }

  Int_t option=3;
  if (fdMethod==kfc) option=1;
  else if (fdMethod==kNb) option=2;
  else if (fdMethod==knone) { option=0; asym=false; }
  else option=3;

  if (option>2) { 
    cout<< "Bad calculation option, should be <=2"<<endl;
    return;
  }


  //
  // Defining the Tab values for the given centrality class
  // https://twiki.cern.ch/twiki/bin/viewauth/ALICE/CentStudies
  //
  Double_t tab = 1., tabUnc = 0.;
  if ( cc == k07half ) {
    tab = 24.81; tabUnc = 0.8037;
  } else if ( cc == k010 ) {
    tab = 23.48; tabUnc = 0.97;
  } else if ( cc == k1020 ) {
    tab = 14.4318; tabUnc = 0.5733;
  } else if ( cc == k020 ) {
    tab = 18.93; tabUnc = 0.74;
  } else if ( cc == k2040 ) {
    tab = 6.86; tabUnc = 0.28;
  } else if ( cc == k2030 ) {
    tab = 8.73769; tabUnc = 0.370219;
  } else if ( cc == k3040 ) {
    tab = 5.02755; tabUnc = 0.22099;
  } else if ( cc == k4050 ) {
    tab = 2.68327; tabUnc = 0.137073;
  } else if ( cc == k3050 ) {
    tab = 3.87011; tabUnc = 0.183847;
  } else if ( cc == k4060 ) {
    tab = 2.00;  tabUnc= 0.11;
  } else if ( cc == k4080 ) {
    tab = 1.20451; tabUnc = 0.071843;
  } else if ( cc == k5060 ) {
    tab = 1.32884; tabUnc = 0.0929536;
  } else if ( cc == k6080 ) {
    tab = 0.419; tabUnc = 0.033;
  } else if ( cc == k80100 ){
    tab = 0.0690; tabUnc = 0.0062;
  }
  tab *= 1e-9; // to pass from mb^{-1} to pb^{-1}
  tabUnc *= 1e-9;



  //
  // Get the histograms from the files
  //
  TH1D *hDirectMCpt=0;           // Input MC c-->D spectra
  TH1D *hFeedDownMCpt=0;         // Input MC b-->D spectra
  TH1D *hDirectMCptMax=0;        // Input MC maximum c-->D spectra
  TH1D *hDirectMCptMin=0;        // Input MC minimum c-->D spectra
  TH1D *hFeedDownMCptMax=0;      // Input MC maximum b-->D spectra
  TH1D *hFeedDownMCptMin=0;      // Input MC minimum b-->D spectra
  //  TGraphAsymmErrors *gPrediction=0; // Input MC c-->D spectra
  TH1D *hDirectEffpt=0;          // c-->D Acceptance and efficiency correction
  TH1D *hFeedDownEffpt=0;        // b-->D Acceptance and efficiency correction
  TH1D *hRECpt=0;                // all reconstructed D

  //
  // Define/Get the input histograms
  //
  Int_t decay=0;
  TFile * mcfile = new TFile(mcfilename,"read");
  if (isD0Kpi){
    decay = 1;
    hDirectMCpt = (TH1D*)mcfile->Get("hD0Kpipred_central");
    hFeedDownMCpt = (TH1D*)mcfile->Get("hD0KpifromBpred_central_corr");
    hDirectMCptMax = (TH1D*)mcfile->Get("hD0Kpipred_max");
    hDirectMCptMin = (TH1D*)mcfile->Get("hD0Kpipred_min");
    hFeedDownMCptMax = (TH1D*)mcfile->Get("hD0KpifromBpred_max_corr");
    hFeedDownMCptMin = (TH1D*)mcfile->Get("hD0KpifromBpred_min_corr");
    //    gPrediction = (TGraphAsymmErrors*)mcfile->Get("D0Kpiprediction");
  }
  else if (isDplusKpipi){
    decay = 2;
    hDirectMCpt = (TH1D*)mcfile->Get("hDpluskpipipred_central");
    hFeedDownMCpt = (TH1D*)mcfile->Get("hDpluskpipifromBpred_central_corr");
    hDirectMCptMax = (TH1D*)mcfile->Get("hDpluskpipipred_max");
    hDirectMCptMin = (TH1D*)mcfile->Get("hDpluskpipipred_min");
    hFeedDownMCptMax = (TH1D*)mcfile->Get("hDpluskpipifromBpred_max_corr");
    hFeedDownMCptMin = (TH1D*)mcfile->Get("hDpluskpipifromBpred_min_corr");
    //    gPrediction = (TGraphAsymmErrors*)mcfile->Get("Dpluskpipiprediction");
  }
  else if(isDstarD0pi){
    decay = 3;
    hDirectMCpt = (TH1D*)mcfile->Get("hDstarD0pipred_central");
    hFeedDownMCpt = (TH1D*)mcfile->Get("hDstarD0pifromBpred_central_corr");
    hDirectMCptMax = (TH1D*)mcfile->Get("hDstarD0pipred_max");
    hDirectMCptMin = (TH1D*)mcfile->Get("hDstarD0pipred_min");
    hFeedDownMCptMax = (TH1D*)mcfile->Get("hDstarD0pifromBpred_max_corr");
    hFeedDownMCptMin = (TH1D*)mcfile->Get("hDstarD0pifromBpred_min_corr");
    //    gPrediction = (TGraphAsymmErrors*)mcfile->Get("DstarD0piprediction");
  }
  else if (isDsKKpi){
    decay = 4;
    hDirectMCpt = (TH1D*)mcfile->Get("hDsPhipitoKkpipred_central");
    hFeedDownMCpt = (TH1D*)mcfile->Get("hDsPhipitoKkpifromBpred_central_corr");
    hDirectMCptMax = (TH1D*)mcfile->Get("hDsPhipitoKkpipred_max");
    hDirectMCptMin = (TH1D*)mcfile->Get("hDsPhipitoKkpipred_min");
    hFeedDownMCptMax = (TH1D*)mcfile->Get("hDsPhipitoKkpifromBpred_max_corr");
    hFeedDownMCptMin = (TH1D*)mcfile->Get("hDsPhipitoKkpifromBpred_min_corr");
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
  TFile * efffile = new TFile(efffilename,"read");
  hDirectEffpt = (TH1D*)efffile->Get("hDirectEffpt");
  hDirectEffpt->SetNameTitle("hDirectEffpt","direct acc x eff");
  hFeedDownEffpt = (TH1D*)efffile->Get("hFeedDownEffpt");
  hFeedDownEffpt->SetNameTitle("hFeedDownEffpt","feed-down acc x eff");
  //
  //
  TFile * recofile = new TFile(recofilename,"read");
  hRECpt = (TH1D*)recofile->Get(recohistoname);
  hRECpt->SetNameTitle("hRECpt","Reconstructed spectra");

  //
  // Read the file of the EP resolution correction
  TFile *EPf;
  TH1D *hEPresolCorr;
  if(isRaavsEP>0.){
    EPf = new TFile(epResolfile,"read");
    if(isRaavsEP==kInPlane) hEPresolCorr = (TH1D*)EPf->Get("hCorrEPresol_InPlane");
    else if(isRaavsEP==kOutOfPlane) hEPresolCorr = (TH1D*)EPf->Get("hCorrEPresol_OutOfPlane");
    for(Int_t i=1; i<=hRECpt->GetNbinsX(); i++) {
      Double_t value = hRECpt->GetBinContent(i);
      Double_t error = hRECpt->GetBinError(i);
      Double_t pt = hRECpt->GetBinCenter(i);
      Int_t epbin = hEPresolCorr->FindBin( pt );
      Double_t epcorr = hEPresolCorr->GetBinContent( epbin );
      value = value*epcorr;
      error = error*epcorr;
      hRECpt->SetBinContent(i,value);
      hRECpt->SetBinError(i,error);
    }
  }

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
  TH2D *histofcRcb=0;
  TH1D *histofcRcb_px=0;
  TH2D *histoYieldCorrRcb=0;
  TH2D *histoSigmaCorrRcb=0;
  //
  TGraphAsymmErrors * gYieldCorr = 0;
  TGraphAsymmErrors * gSigmaCorr = 0;
  TGraphAsymmErrors * gFcExtreme = 0;
  TGraphAsymmErrors * gFcConservative = 0;
  TGraphAsymmErrors * gYieldCorrExtreme = 0;
  TGraphAsymmErrors * gSigmaCorrExtreme = 0;
  TGraphAsymmErrors * gYieldCorrConservative = 0;
  TGraphAsymmErrors * gSigmaCorrConservative = 0;
  //
  TNtuple * nSigma = 0;


  //
  // Main functionalities for the calculation
  //

  // Define and set the basic option flags
  AliHFPtSpectrum * spectra = new AliHFPtSpectrum("AliHFPtSpectrum","AliHFPtSpectrum",option);
  spectra->SetFeedDownCalculationOption(option);
  spectra->SetComputeAsymmetricUncertainties(asym);
  // Set flag on whether to additional PbPb Eloss hypothesis have to be computed
  spectra->SetComputeElossHypothesis(PbPbEloss);

  // Feed the input histograms
  //  reconstructed spectra
  cout << " Setting the reconstructed spectrum,";
  spectra->SetReconstructedSpectrum(hRECpt);
  // acceptance and efficiency corrections
  cout << " the efficiency,";
  spectra->SetAccEffCorrection(hDirectEffpt,hFeedDownEffpt);
  //    spectra->SetfIsStatUncEff(false);
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
  spectra->SetNormalization(nevents,sigma);
  Double_t lumi = nevents / sigma ;
  Double_t lumiUnc = 0.04*lumi; // 10% uncertainty on the luminosity
  spectra->SetLuminosity(lumi,lumiUnc);
  Double_t effTrig = 1.0;
  spectra->SetTriggerEfficiency(effTrig,0.);
  if(isRaavsEP>0.) spectra->SetIsEventPlaneAnalysis(kTRUE);

  // Set the global uncertainties on the efficiencies (in percent)
  Double_t globalEffUnc = 0.15; 
  Double_t globalBCEffRatioUnc = 0.15;
  spectra->SetAccEffPercentageUncertainty(globalEffUnc,globalBCEffRatioUnc);

  // Set if the yield is for particle+anti-particle or only one type
  spectra->SetIsParticlePlusAntiParticleYield(isParticlePlusAntiParticleYield);

  // Set the Tab parameter and uncertainties
  if ( (cc != kpp7) && (cc != kpp276) ) {
    spectra->SetTabParameter(tab,tabUnc);
  }

  // Do the calculations
  cout << " Doing the calculation... "<< endl;
  Double_t deltaY = 1.0;
  Double_t branchingRatioC = 1.0;
  Double_t branchingRatioBintoFinalDecay = 1.0; // this is relative to the input theoretical prediction
  spectra->ComputeHFPtSpectrum(deltaY,branchingRatioC,branchingRatioBintoFinalDecay);
  cout << "   ended the calculation, getting the histograms back " << endl;

  // Set the systematics externally
  Bool_t combineFeedDown = true;
  AliHFSystErr *systematics = new AliHFSystErr();
  if( cc==kpp276 ) {
    systematics->SetIsLowEnergy(true);
  } else if( cc!=kpp7 )  {
    systematics->SetCollisionType(1);
    if ( cc == k07half ) systematics->SetCentrality("07half");
    else if ( cc == k010 )  systematics->SetCentrality("010");
    else if ( cc == k1020 )  systematics->SetCentrality("1020");
    else if ( cc == k020 )  systematics->SetCentrality("020");
    else if ( cc == k2040 || cc == k2030 || cc == k3040 ) {
      systematics->SetCentrality("2040");
      systematics->SetIsPbPb2010EnergyScan(true);
    }
    else if ( cc == k3050 ) {
      if (isRaavsEP == kPhiIntegrated) systematics->SetCentrality("4080");
      else if (isRaavsEP == kInPlane) systematics->SetCentrality("3050InPlane");
      else if (isRaavsEP == kOutOfPlane) systematics->SetCentrality("3050OutOfPlane");
    }
    else if ( cc == k4060 || cc == k4050 || cc == k5060 )  systematics->SetCentrality("4060");
    else if ( cc == k6080 )  systematics->SetCentrality("6080");
    else if ( cc == k4080 ) systematics->SetCentrality("4080");
    else {
      cout << " Systematics not yet implemented " << endl;
      return;
    }
  } else { systematics->SetCollisionType(0); }
  systematics->Init(decay);
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
  // Get the PbPb Eloss hypothesis histograms
  if(PbPbEloss){
    histofcRcb = spectra->GetHistoFeedDownCorrectionFcVsEloss();
    histoYieldCorrRcb = spectra->GetHistoFeedDownCorrectedSpectrumVsEloss();
    histoSigmaCorrRcb = spectra->GetHistoCrossSectionFromYieldSpectrumVsEloss();
    histofcRcb->SetName("histofcRcb");
    histoYieldCorrRcb->SetName("histoYieldCorrRcb");
    histoSigmaCorrRcb->SetName("histoSigmaCorrRcb");
  }

  // Get & Rename the TGraphs
  gSigmaCorr = spectra->GetCrossSectionFromYieldSpectrum();
  gYieldCorr = spectra->GetFeedDownCorrectedSpectrum();
  if (asym) {
    gSigmaCorrExtreme = spectra->GetCrossSectionFromYieldSpectrumExtreme();
    gYieldCorrExtreme = spectra->GetFeedDownCorrectedSpectrumExtreme();
    gSigmaCorrConservative = spectra->GetCrossSectionFromYieldSpectrumConservative();
    gYieldCorrConservative = spectra->GetFeedDownCorrectedSpectrumConservative();
  }

  // Get & Rename the TGraphs
  if (option==0){
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

  if(PbPbEloss){
    nSigma = spectra->GetNtupleCrossSectionVsEloss();
  }

  //
  // Now, plot the results ! :)
  //

  gROOT->SetStyle("Plain");

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
 
  // Draw the PbPb Eloss hypothesis histograms
  if(PbPbEloss){
    AliHFPtSpectrum *CalcBins;
    gStyle->SetPalette(1);
    TCanvas *canvasfcRcb = new TCanvas("canvasfcRcb","fc vs pt vs Rcb");
    //    histofcRcb->Draw("cont4z");
    histofcRcb->Draw("colz");
    canvasfcRcb->Update();
    canvasfcRcb->cd(2);
    TCanvas *canvasfcRcb1 = new TCanvas("canvasfcRcb1","fc vs pt vs Rcb=1");
    histofcRcb_px = (TH1D*)histofcRcb->ProjectionX("histofcRcb_px",40,40);
    histofcRcb_px->SetLineColor(2);
    if (option==1) {
      histofc->Draw();
      histofcRcb_px->Draw("same");
    } else histofcRcb_px->Draw("");
    canvasfcRcb1->Update();
    TCanvas *canvasfcRcb2 = new TCanvas("canvasfcRcb2","fc vs pt vs Rcb fixed Rcb");
    Int_t bin0 = CalcBins->FindTH2YBin(histofcRcb,0.25);
    Int_t bin1 = CalcBins->FindTH2YBin(histofcRcb,0.5);
    Int_t bin2 = CalcBins->FindTH2YBin(histofcRcb,1.0);
    Int_t bin3 = CalcBins->FindTH2YBin(histofcRcb,1.5);
    Int_t bin4 = CalcBins->FindTH2YBin(histofcRcb,2.0);
    Int_t bin5 = CalcBins->FindTH2YBin(histofcRcb,3.0);
    Int_t bin6 = CalcBins->FindTH2YBin(histofcRcb,4.0);
    TH1D * histofcRcb_px0a = (TH1D*)histofcRcb->ProjectionX("histofcRcb_px0a",bin0,bin0);
    TH1D * histofcRcb_px0 = (TH1D*)histofcRcb->ProjectionX("histofcRcb_px0",bin1,bin1);
    TH1D * histofcRcb_px1 = (TH1D*)histofcRcb->ProjectionX("histofcRcb_px1",bin2,bin2);
    TH1D * histofcRcb_px2 = (TH1D*)histofcRcb->ProjectionX("histofcRcb_px2",bin3,bin3);
    TH1D * histofcRcb_px3 = (TH1D*)histofcRcb->ProjectionX("histofcRcb_px3",bin4,bin4);
    TH1D * histofcRcb_px4 = (TH1D*)histofcRcb->ProjectionX("histofcRcb_px4",bin5,bin5);
    TH1D * histofcRcb_px5 = (TH1D*)histofcRcb->ProjectionX("histofcRcb_px5",bin6,bin6);
    if (option==1) {
      histofc->Draw();
      //      histofcRcb_px->Draw("same");
    } else {
      //      histofcRcb_px->Draw("");
      histofcRcb_px0a->SetLineColor(2);
      histofcRcb_px0a->Draw("");
    }
    histofcRcb_px0a->SetLineColor(2);
    histofcRcb_px0a->Draw("same");
    histofcRcb_px0->SetLineColor(4);
    histofcRcb_px0->Draw("same");
    histofcRcb_px1->SetLineColor(3);
    histofcRcb_px1->Draw("same");
    histofcRcb_px2->SetLineColor(kCyan);
    histofcRcb_px2->Draw("same");
    histofcRcb_px3->SetLineColor(kMagenta+1);
    histofcRcb_px3->Draw("same");
    histofcRcb_px4->SetLineColor(kOrange+7);
    histofcRcb_px4->Draw("same");
    histofcRcb_px5->SetLineColor(kGreen+3);
    histofcRcb_px5->Draw("same");
    TLegend *legrcc = new TLegend(0.8,0.8,0.95,0.9);
    legrcc->SetFillColor(0);
    if (option==1) {
      legrcc->AddEntry(histofcRcb_px0a,"Rc/b=0.25","l");
      legrcc->AddEntry(histofcRcb_px0,"Rc/b=0.5","l");
      legrcc->AddEntry(histofcRcb_px1,"Rc/b=1.0","l");
      legrcc->AddEntry(histofcRcb_px2,"Rc/b=1.5","l");
      legrcc->AddEntry(histofcRcb_px3,"Rc/b=2.0","l");
      legrcc->AddEntry(histofcRcb_px4,"Rc/b=3.0","l");
      legrcc->AddEntry(histofcRcb_px5,"Rc/b=4.0","l");
    }else{
      legrcc->AddEntry(histofcRcb_px0a,"Rb=0.25","l");
      legrcc->AddEntry(histofcRcb_px0,"Rb=0.5","l");
      legrcc->AddEntry(histofcRcb_px1,"Rb=1.0","l");
      legrcc->AddEntry(histofcRcb_px2,"Rb=1.5","l");
      legrcc->AddEntry(histofcRcb_px3,"Rb=2.0","l");
      legrcc->AddEntry(histofcRcb_px4,"Rb=3.0","l");
      legrcc->AddEntry(histofcRcb_px5,"Rb=4.0","l");
    }
    legrcc->Draw();
    canvasfcRcb2->Update();
    TCanvas *canvasYRcb = new TCanvas("canvasYRcb","corrected yield vs pt vs Rcb");
    histoYieldCorrRcb->Draw("cont4z");
    canvasYRcb->Update();
    TCanvas *canvasSRcb = new TCanvas("canvasSRcb","sigma vs pt vs Rcb");
    histoSigmaCorrRcb->Draw("cont4z");
    canvasSRcb->Update();
    TCanvas *canvasSRcb1 = new TCanvas("canvasSRcb1","sigma vs pt vs Rcb fixed Rcb");
    TH1D * histoSigmaCorrRcb_px0a = (TH1D*)histoSigmaCorrRcb->ProjectionX("histoSigmaCorrRcb_px0a",bin0,bin0);
    TH1D * histoSigmaCorrRcb_px0 = (TH1D*)histoSigmaCorrRcb->ProjectionX("histoSigmaCorrRcb_px0",bin1,bin1);
    TH1D * histoSigmaCorrRcb_px1 = (TH1D*)histoSigmaCorrRcb->ProjectionX("histoSigmaCorrRcb_px1",bin2,bin2);
    TH1D * histoSigmaCorrRcb_px2 = (TH1D*)histoSigmaCorrRcb->ProjectionX("histoSigmaCorrRcb_px2",bin3,bin3);
    TH1D * histoSigmaCorrRcb_px3 = (TH1D*)histoSigmaCorrRcb->ProjectionX("histoSigmaCorrRcb_px3",bin4,bin4);
    TH1D * histoSigmaCorrRcb_px4 = (TH1D*)histoSigmaCorrRcb->ProjectionX("histoSigmaCorrRcb_px4",bin5,bin5);
    TH1D * histoSigmaCorrRcb_px5 = (TH1D*)histoSigmaCorrRcb->ProjectionX("histoSigmaCorrRcb_px5",bin6,bin6);
    histoSigmaCorr->Draw();
    histoSigmaCorrRcb_px0a->SetLineColor(2);
    histoSigmaCorrRcb_px0a->Draw("hsame");
    histoSigmaCorrRcb_px0->SetLineColor(4);
    histoSigmaCorrRcb_px0->Draw("hsame");
    histoSigmaCorrRcb_px1->SetLineColor(3);
    histoSigmaCorrRcb_px1->Draw("hsame");
    histoSigmaCorrRcb_px2->SetLineColor(kCyan);
    histoSigmaCorrRcb_px2->Draw("hsame");
    histoSigmaCorrRcb_px3->SetLineColor(kMagenta+1);
    histoSigmaCorrRcb_px3->Draw("hsame");
    histoSigmaCorrRcb_px4->SetLineColor(kOrange+7);
    histoSigmaCorrRcb_px4->Draw("same");
    histoSigmaCorrRcb_px5->SetLineColor(kGreen+3);
    histoSigmaCorrRcb_px5->Draw("same");
    TLegend *legrcb = new TLegend(0.8,0.8,0.95,0.9);
    legrcb->SetFillColor(0);
    if (option==1) {
      legrcb->AddEntry(histoSigmaCorrRcb_px0a,"Rc/b=0.25","l");
      legrcb->AddEntry(histoSigmaCorrRcb_px0,"Rc/b=0.5","l");
      legrcb->AddEntry(histoSigmaCorrRcb_px1,"Rc/b=1.0","l");
      legrcb->AddEntry(histoSigmaCorrRcb_px2,"Rc/b=1.5","l");
      legrcb->AddEntry(histoSigmaCorrRcb_px3,"Rc/b=2.0","l");
      legrcb->AddEntry(histoSigmaCorrRcb_px4,"Rc/b=3.0","l");
      legrcb->AddEntry(histoSigmaCorrRcb_px5,"Rc/b=4.0","l");
    }else{
      legrcb->AddEntry(histoSigmaCorrRcb_px0a,"Rb=0.25","l");
      legrcb->AddEntry(histoSigmaCorrRcb_px0,"Rb=0.5","l");
      legrcb->AddEntry(histoSigmaCorrRcb_px1,"Rb=1.0","l");
      legrcb->AddEntry(histoSigmaCorrRcb_px2,"Rb=1.5","l");
      legrcb->AddEntry(histoSigmaCorrRcb_px3,"Rb=2.0","l");
      legrcb->AddEntry(histoSigmaCorrRcb_px4,"Rb=3.0","l");
      legrcb->AddEntry(histoSigmaCorrRcb_px5,"Rb=4.0","l");
    }
    legrcb->Draw();
    canvasSRcb1->Update();
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

  if(PbPbEloss){
    histofcRcb->Write();    histofcRcb_px->Write();
    histoYieldCorrRcb->Write();
    histoSigmaCorrRcb->Write();
    nSigma->Write();
  }

  gYieldCorr->Write();
  gSigmaCorr->Write();
  if(asym){
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


  TH1D * hStatUncEffcSigma = spectra->GetDirectStatEffUncOnSigma();
  TH1D * hStatUncEffbSigma = spectra->GetFeedDownStatEffUncOnSigma();
  hStatUncEffcSigma->Write();
  hStatUncEffbSigma->Write();
  if(option!=0){
    TH1D * hStatUncEffcFD = spectra->GetDirectStatEffUncOnFc();
    TH1D * hStatUncEffbFD = spectra->GetFeedDownStatEffUncOnFc();
    hStatUncEffcFD->Write();
    hStatUncEffbFD->Write();
  }

  // Draw the cross-section 
  //  spectra->DrawSpectrum(gPrediction);

  //  out->Close();

}
