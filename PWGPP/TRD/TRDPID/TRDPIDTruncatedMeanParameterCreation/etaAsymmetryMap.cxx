//====================================================================================================================================
//
// root macro to parametrise dE/dx vs. pseudorapidity asymmetry   (for >= 4 tracklets)
//
//====================================================================================================================================

#include <iostream>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TMinuit.h>
#include <TList.h>
#include <TF1.h>
//#include "V0Tracks.cxx"
#include "constants.h"
#include "utils.h"


void etaAsymmetryMap() {
    // pdg codes for particle species
    // if no charge seperation is wanted, set tracks.pId() to abs
    Int_t pdgCodes[3] = {11, 211, 2212};

    // get and initialize tree
    printf("etaAsymmetryMap(): get tree %s\n", input.Data());
    V0Tracks tracks;
    tracks.setFile(input.Data());
    tracks.getTree();
    
    // get number of entries from tree
    const Int_t nEntries = tracks.getNumberOfEntries();
    
    // declare histograms
    TH3D * dEdxVsEtaVsBG    = new TH3D("dEdxVsEtaVsBG", "dE/dx vs. #eta vs. #beta#gamma", 180, -0.9, 0.9, 100, 0.3, 1e4, 200, 0, 10);
    TH3D * dEdxVsTPCtglVsBG = new TH3D("dEdxVsTPCtglVsBG", "dE/dx vs. TPCtgl vs. #beta#gamma", 180, -1., 1., 100, 0.3, 1e4, 200, 0, 10);
    TH2D * map              = new TH2D("map",           "corr. map",                      180, -0.9, 0.9, 100, 0.3, 1e4);
    TH2D * mapTPCtgl        = new TH2D("mapTPCtgl",    "corr. map for TPCtgl",            180, -1., 1., 100, 0.3, 1e4);
    TH2D * TPCtglVsBG       = new TH2D("TPCtglVsBG", "TPCtglVsBG", 180, -1., 1., 100, 0.3, 1e4);
    TH2D * Testmap          = new TH2D("Testmap", "Testmap", 180, -1., 1., 100, 0.3, 1e4);
    TH2D * BGVsdEdxMean     = new TH2D("BGVsdEdxMean", "BGVsdEdxMean", 100, 0.3, 1e4, 200, 0, 10);
    TH1D * BGMean           = new TH1D("BGMean", "BGMean", 100, 0.3, 1e4);
    TH1D * BGMeanEta         = new TH1D("BGMeanEta", "BGMeanEta", 100, 0.3, 1e4);
    TH1D * BGMeanNCls[3];

    TH3D * dEdxVsNclsVsBGFull= new TH3D("dEdxVsNclsVsBG", "dEdxVsNclsVsBG", 120, 40, 160, 100, 0.3, 1e4, 200, 0, 10);
    TH3D * dEdxVsNclsVsBG[3];    
    TH2D * NclsMap[3];
    TH2D * NclsVsBG[3];

    TH3D * dEdxVsNclsVsEtaFull= new TH3D("dEdxVsNclsVsBG", "dEdxVsNclsVsBG", 120, 40, 160, 180, -1, 1, 200, 0, 10);
    TH3D * dEdxVsNclsVsEta[3];    
    TH2D * EtaNclsMap[3];

    Double_t BGBins[101];
    Double_t dEdxBins[201];
    BGBins[0]=0.3;
    BGBins[100]=1e4;
    dEdxBins[0]=0;
    dEdxBins[200]=10;

    for (int i=0; i<100; i++) {
      BGBins[i+1]=BGBins[i]+(BGBins[100]-BGBins[0])/100;
      //   cout << i+1 << "  " << BGBins[i+1] << endl;
    }
    for (int i=0; i<200; i++) {
      dEdxBins[i+1]=dEdxBins[i]+(dEdxBins[200]-dEdxBins[0])/200;
      // cout << i+1 << "  " << dEdxBins[i+1] << endl;
    }
    
    TH3D * dEdxVsCentVsBG = new TH3D("dEdxVsCentVsBG", "dEdxVsCentVsBG", NCentBins,  CentBins, 100, BGBins, 200, dEdxBins);
    TH2D * CentVsBG = new TH2D("CentVsBG", "CentVsBG", NCentBins, CentBins, 100, BGBins);
    TH2D * CentMap = new TH2D("CentMap", "CentMap", NCentBins, CentBins, 100, BGBins);
    TH1D * CentMapRef = new TH1D("CentMapRef", "CentMapRef", 100, BGBins);


    // declare array of histograms and set log axis
    for (int i=0; i<3; i++) { // for each #tracklets
      dEdxVsNclsVsBG[i] = new TH3D(Form("dEdxVsNclsVsBG_nCh%i", i+4), Form("dEdxVsNclsVsBG_nCh%i", i+4),  120, 40, 160, 100, 0.3, 1e4, 200, 0, 10);
      NclsMap[i] = new TH2D(Form("NclsMap_nCh%i", i+4), Form( "NclsMap_nCh%i", i+4), 120, 40, 160, 100, 0.3, 1e4);
      NclsVsBG[i] = new TH2D(Form("NclsVsBG_nCh%i", i+4), Form( "NclsVsBG_nCh%i", i+4), 120, 40, 160, 100, 0.3, 1e4);
      BinLogX(NclsMap[i]->GetYaxis());
      BinLogX(NclsVsBG[i]->GetYaxis());
      BinLogX(dEdxVsNclsVsBG[i]->GetYaxis());

      dEdxVsNclsVsEta[i] = new TH3D(Form("dEdxVsNclsVsEta_nCh%i", i+4), Form("dEdxVsNclsVsEta_nCh%i", i+4),  120, 40, 160,180, -1., 1., 200, 0, 10);
      EtaNclsMap[i] = new TH2D(Form("EtaNclsMap_nCh%i", i+4), Form( "EtaNclsMap_nCh%i", i+4), 120, 40, 160, 180, -1, 1.);

      BGMeanNCls[i]    = new TH1D(Form("BGMeanNCls%i", i+4), Form("BGMeanNCls%i", i+4), 100, 0.3, 1e4);
      BinLogX(BGMeanNCls[i]->GetXaxis());
    }

    // set log axis
    BinLogX(dEdxVsEtaVsBG->GetYaxis());
    BinLogX(map->GetYaxis());
    BinLogX(dEdxVsTPCtglVsBG->GetYaxis());
    BinLogX(mapTPCtgl->GetYaxis());
    BinLogX(TPCtglVsBG->GetYaxis());
    BinLogX(dEdxVsNclsVsBGFull->GetYaxis());
    BinLogX(dEdxVsNclsVsEtaFull->GetYaxis());
    BinLogX(BGVsdEdxMean->GetXaxis());
    BinLogX(BGMean->GetXaxis());
    BinLogX(BGMeanEta->GetXaxis());
    BinLogX(dEdxVsCentVsBG->GetYaxis());
    BinLogX(CentMap->GetYaxis());
    BinLogX(CentVsBG->GetYaxis());
    BinLogX(CentMapRef->GetXaxis());
 
    // mean parametrizations (only use nch = 6) - this comment is old: Current code uses  nch >=4.
    printf("etaAsymmetryMap(): getting mean param\n");
    TString EtaSpec, NclsCorr, CentCorr;
    if(!fEtaCorrection) EtaSpec = ""; 
    else EtaSpec=EtaSpecStr;
    if(!NclsCorrection) NclsCorr = "";
    else NclsCorr=NclsCorrStr;
    if (!CentCorrection) CentCorr="";
    else CentCorr=CentCorrStr;


    TFile * meanParamFile = TFile::Open(Form("output/%s_meanParam_pid%i_%s%s%s%s.root", inputSpec.Data(),  pidOpt, FileSpec.Data(),CentCorr.Data(), EtaSpec.Data(), NclsCorr.Data()));
    TDirectory * meanParamFileBasic = (TDirectory*)meanParamFile->Get("Basic");
    TH1D * meanParam = 0;
    meanParam = (TH1D*)meanParamFileBasic->Get("meanFitHist");
    if (!meanParam) {
      cout<< "Error: could not open meanParam" << endl;
      return;
    }

    // loading correction maps. slightly confusing structure at the moment due to iterative structure
    TH2D * etaMap2;
    if (fEtaCorrection) {
      //TString tmp(NclsCorr);
      //NclsCorr="";
      TFile * etaMapFile = TFile::Open(Form("output/%s_etaMap_%s_%s%s.root", inputSpec.Data(),  FileSpec.Data(),CentCorr.Data(), NclsCorr.Data()), "READ");
      //etaMap = (TH2D*)etaMapFile->Get("map");
      etaMap2 = (TH2D*)etaMapFile->Get("mapTPCtgl");
      //NclsCorr=tmp;
    }
   
    TH2D * NclsMaps2[3];
    if (NclsCorrection) {
      TString tmp(EtaSpec);
      EtaSpec="";
      TFile * ClusterMapFile = TFile::Open(Form("output/%s_ClusterMap_%s_%s%s.root", inputSpec.Data(),  FileSpec.Data(), CentCorr.Data(), EtaSpec.Data()), "READ");
      for (int i=0; i<3; i++) {
	NclsMaps2[i] = (TH2D*)ClusterMapFile->Get(Form("NclsMap_nCh%i", i+4));
      }
      EtaSpec=tmp;
    }
    
    TH2D * CentMap2;
    if (CentCorrection) {
      TString tmp(EtaSpec);
      TString tmp2(NclsCorr);
      EtaSpec="";
      NclsCorr="";
      TFile * CentMapFile = TFile::Open(Form("output/%s_CentMap_%s_%s%s.root", inputSpec.Data(), FileSpec.Data(), EtaSpec.Data(), NclsCorr.Data()), "READ");
      CentMap2 = (TH2D*)CentMapFile->Get("CentMap");
      EtaSpec=tmp;
      NclsCorr=tmp2;
    }

    //
    // loop over tree to fill histograms
    //
    printf("etaAsymmetryMap(): looping over tree\n");
    Int_t nb    = 0;
    Int_t iter  = 0;
    Int_t nchTmp = nchMin;
    nchMin=4;
    Int_t run;
    cout << "nEntries of Tree: " <<  nEntries << endl;
    Int_t cut[6]={0,0,0,0,0,0};
    for (Int_t i = 0; i < nEntries; i++) {
        if (i%1000000==0) cout << i << endl;
        nb = tracks.getEntry(i);
        if (nb <= 0) continue;

	if (runCut(tracks.run())) continue;

	if (TMath::Abs(tracks.trdTPCtgl())>1.) continue; // no correction for higher eta
	// similar cuts for cluster or centrality?
        
        // apply general quality cuts on tracks - equal for all in this case (ref Pion)
	if (trackAndLayerCuts(tracks.trdNCh(), tracks.trdNCls(), 211)) continue; // TMath::Abs(tracks.pId()))) continue;
	//if (trackAndLayerCuts(tracks.trdNCh(), tracks.trdNCls(), TMath::Abs(tracks.pId()))) continue;
        
        // apply additional particle identification cuts on tracks
	if (nSigmaCut(tracks.nSigmaTPC(), (Double_t*)cutNSigmaTPC, TMath::Abs(tracks.pId()))) {cut[3]++; continue;};
	if (nSigmaCut(tracks.nSigmaTOF(), (Double_t*)cutNSigmaTOF, TMath::Abs(tracks.pId()))) {cut[4]++; continue;};
        
	// declare correction factor
        Double_t CorrectionFactor = 1;

        // electrons
        if (TMath::Abs(tracks.pId()) == pdgCodes[0]) {
	  // if (betaGamma(momentumMean(tracks.trdMom()), mElectron) < bgMin || betaGamma(momentumMean(tracks.trdMom()), mElectron) >= bgMax) continue;
	  // if (electronMomentumCut(tracks.trdMom(), tracks.nSigmaTPC()[0])) continue;
	  
	  if (fEtaCorrection) CorrectionFactor *= EtaCorrFactorTPCtgl(betaGamma(momentumMean(tracks.trdMom()), mElectron), tracks.trdTPCtgl(), etaMap2);
	  if (NclsCorrection) CorrectionFactor *= NclsCorrFactor(betaGamma(momentumMean(tracks.trdMom()), mElectron), tracks.trdNCls(), NclsMaps2[tracks.trdNCh()-4]);
	  if (CentCorrection) CorrectionFactor *= CentCorrFactor(betaGamma(momentumMean(tracks.trdMom()), mElectron), tracks.centrality(), CentMap2);
	  if (CorrectionFactor<1E-5) continue;	 

	  BGVsdEdxMean->Fill(betaGamma(momentumMean(tracks.trdMom()), mElectron), CorrectionFactor*tracks.trdSig());
	  //if (checkSliceContent(dEdxVsEtaVsBG,tracks.trdEtaLocal()[0], betaGamma(momentumMean(tracks.trdMom()), mElectron)))  
	  //  dEdxVsEtaVsBG->Fill(tracks.trdEtaLocal()[0], betaGamma(momentumMean(tracks.trdMom()), mElectron),  CorrectionFactor*tracks.trdSig());
	  TPCtglVsBG->Fill(tracks.trdTPCtgl(), betaGamma(momentumMean(tracks.trdMom()), mElectron));
	  if (checkSliceContent(TPCtglVsBG, tracks.trdTPCtgl(), betaGamma(momentumMean(tracks.trdMom()), mElectron)))  
	    dEdxVsTPCtglVsBG->Fill(tracks.trdTPCtgl(), betaGamma(momentumMean(tracks.trdMom()), mElectron),  CorrectionFactor*tracks.trdSig());
	
	  //if (checkSliceContent( dEdxVsNclsVsBGFull, tracks.trdNCls(), betaGamma(momentumMean(tracks.trdMom()), mElectron)))  
	  // dEdxVsNclsVsBGFull->Fill(tracks.trdNCls(), betaGamma(momentumMean(tracks.trdMom()), mElectron), CorrectionFactor*tracks.trdSig());
	  NclsVsBG[tracks.trdNCh()-4]->Fill(tracks.trdNCls(), betaGamma(momentumMean(tracks.trdMom()), mElectron)); 
	  if (checkSliceContent( NclsVsBG[tracks.trdNCh()-4], tracks.trdNCls(), betaGamma(momentumMean(tracks.trdMom()), mElectron)))  
	    dEdxVsNclsVsBG[tracks.trdNCh()-4]->Fill(tracks.trdNCls(), betaGamma(momentumMean(tracks.trdMom()), mElectron), CorrectionFactor*tracks.trdSig());
	  // if (checkSliceContent(dEdxVsNclsVsEtaFull, tracks.trdNCls(), tracks.trdTPCtgl()))
	  // dEdxVsNclsVsEtaFull->Fill(tracks.trdNCls(), tracks.trdTPCtgl(), CorrectionFactor*tracks.trdSig());
	  //if (checkSliceContent( dEdxVsNclsVsEta[tracks.trdNCh()-4], tracks.trdNCls(), tracks.trdTPCtgl()))  
	  // dEdxVsNclsVsEta[tracks.trdNCh()-4]->Fill(tracks.trdNCls(),tracks.trdTPCtgl(), CorrectionFactor*tracks.trdSig());
	  CentVsBG->Fill(tracks.centrality(), betaGamma(momentumMean(tracks.trdMom()), mElectron));
	  if (checkSliceContent( CentVsBG, tracks.centrality(), betaGamma(momentumMean(tracks.trdMom()), mElectron)))  
	    dEdxVsCentVsBG->Fill(tracks.centrality(), betaGamma(momentumMean(tracks.trdMom()), mElectron), CorrectionFactor*tracks.trdSig());
        }
        
               
        // pions
        if (TMath::Abs(tracks.pId()) == pdgCodes[1]) {
	  // if (betaGamma(momentumMean(tracks.trdMom()), mPion) < bgMin || betaGamma(momentumMean(tracks.trdMom()), mPion) >= bgMax) continue;
	  // if (momentumMean(tracks.trdMom()) > 0.8 && momentumMean(tracks.trdMom()) < 3) continue;
	  // if (pionMomentumCut(tracks.trdMom())) continue;
	  
	  if (fEtaCorrection) CorrectionFactor *= EtaCorrFactorTPCtgl(betaGamma(momentumMean(tracks.trdMom()), mPion), tracks.trdTPCtgl(), etaMap2);
	  if (NclsCorrection) CorrectionFactor *= NclsCorrFactor(betaGamma(momentumMean(tracks.trdMom()), mPion), tracks.trdNCls(), NclsMaps2[tracks.trdNCh()-4]);
	  if (CentCorrection) CorrectionFactor *= CentCorrFactor(betaGamma(momentumMean(tracks.trdMom()), mPion), tracks.centrality(), CentMap2);
	  if (CorrectionFactor==0) continue;
	  
	  BGVsdEdxMean->Fill(betaGamma(momentumMean(tracks.trdMom()), mPion), CorrectionFactor*tracks.trdSig());
	  //if (checkSliceContent(dEdxVsEtaVsBG,tracks.trdEtaLocal()[0], betaGamma(momentumMean(tracks.trdMom()), mPion)))  
	  // dEdxVsEtaVsBG->Fill(tracks.trdEtaLocal()[0], betaGamma(momentumMean(tracks.trdMom()), mPion), CorrectionFactor*tracks.trdSig());
	  TPCtglVsBG->Fill(tracks.trdTPCtgl(), betaGamma(momentumMean(tracks.trdMom()), mPion));
	  if (checkSliceContent(TPCtglVsBG, tracks.trdTPCtgl(), betaGamma(momentumMean(tracks.trdMom()), mPion)))  
	    dEdxVsTPCtglVsBG->Fill(tracks.trdTPCtgl(), betaGamma(momentumMean(tracks.trdMom()), mPion), CorrectionFactor*tracks.trdSig());
	  // if (checkSliceContent( dEdxVsNclsVsBGFull, tracks.trdNCls(), betaGamma(momentumMean(tracks.trdMom()), mPion)))  
	  // dEdxVsNclsVsBGFull->Fill(tracks.trdNCls(), betaGamma(momentumMean(tracks.trdMom()), mPion), CorrectionFactor*tracks.trdSig());
	  NclsVsBG[tracks.trdNCh()-4]->Fill( tracks.trdNCls(), betaGamma(momentumMean(tracks.trdMom()), mPion));
	  if (checkSliceContent( NclsVsBG[tracks.trdNCh()-4], tracks.trdNCls(), betaGamma(momentumMean(tracks.trdMom()), mPion)))  
	    dEdxVsNclsVsBG[tracks.trdNCh()-4]->Fill(tracks.trdNCls(), betaGamma(momentumMean(tracks.trdMom()), mPion), CorrectionFactor*tracks.trdSig());
	  // if (checkSliceContent(dEdxVsNclsVsEtaFull, tracks.trdNCls(), tracks.trdTPCtgl()))
	  //  dEdxVsNclsVsEtaFull->Fill(tracks.trdNCls(), tracks.trdTPCtgl(), CorrectionFactor*tracks.trdSig());
	  // if (checkSliceContent( dEdxVsNclsVsEta[tracks.trdNCh()-4], tracks.trdNCls(), tracks.trdTPCtgl()))  
	  //dEdxVsNclsVsEta[tracks.trdNCh()-4]->Fill(tracks.trdNCls(),tracks.trdTPCtgl(), CorrectionFactor*tracks.trdSig());
	  CentVsBG->Fill(tracks.centrality(), betaGamma(momentumMean(tracks.trdMom()), mPion));  
       	  if (checkSliceContent( CentVsBG, tracks.centrality(), betaGamma(momentumMean(tracks.trdMom()), mPion)))  
	    dEdxVsCentVsBG->Fill(tracks.centrality(), betaGamma(momentumMean(tracks.trdMom()), mPion), CorrectionFactor*tracks.trdSig());
	}
        
        // proton
        if (TMath::Abs(tracks.pId()) == pdgCodes[2]) {
	  //if (betaGamma(momentumMean(tracks.trdMom()), mProton) < bgMin || bealicetaGamma(momentumMean(tracks.trdMom()), mProton) >= bgMax) continue;
	  //if (betaGamma(momentumMean(tracks.trdMom()), mProton) > 1) continue;
	  // if (protonMomentumCut(tracks.trdMom())) continue;

	  if (fEtaCorrection) CorrectionFactor *= EtaCorrFactorTPCtgl(betaGamma(momentumMean(tracks.trdMom()), mProton), tracks.trdTPCtgl(), etaMap2);
	  if (NclsCorrection) CorrectionFactor *= NclsCorrFactor(betaGamma(momentumMean(tracks.trdMom()), mProton), tracks.trdNCls(), NclsMaps2[tracks.trdNCh()-4]);
	  if (CentCorrection) CorrectionFactor *= CentCorrFactor(betaGamma(momentumMean(tracks.trdMom()), mProton), tracks.centrality(), CentMap2);	 
	  if (CorrectionFactor==0) continue; 

	  BGVsdEdxMean->Fill(betaGamma(momentumMean(tracks.trdMom()), mProton), CorrectionFactor*tracks.trdSig());
	  // if (checkSliceContent(dEdxVsEtaVsBG,tracks.trdEtaLocal()[0], betaGamma(momentumMean(tracks.trdMom()), mProton)))  
	  // dEdxVsEtaVsBG->Fill(tracks.trdEtaLocal()[0], betaGamma(momentumMean(tracks.trdMom()), mProton), CorrectionFactor*tracks.trdSig());
	  TPCtglVsBG->Fill(tracks.trdTPCtgl(), betaGamma(momentumMean(tracks.trdMom()), mProton));
	  if (checkSliceContent(TPCtglVsBG, tracks.trdTPCtgl(), betaGamma(momentumMean(tracks.trdMom()), mProton)))  
	    dEdxVsTPCtglVsBG->Fill(tracks.trdTPCtgl(), betaGamma(momentumMean(tracks.trdMom()), mProton), CorrectionFactor*tracks.trdSig());         
	  //if (checkSliceContent( dEdxVsNclsVsBGFull, tracks.trdNCls(), betaGamma(momentumMean(tracks.trdMom()), mProton)))  
	  // dEdxVsNclsVsBGFull->Fill(tracks.trdNCls(), betaGamma(momentumMean(tracks.trdMom()), mProton), CorrectionFactor*tracks.trdSig());
	  NclsVsBG[tracks.trdNCh()-4]->Fill( tracks.trdNCls(), betaGamma(momentumMean(tracks.trdMom()), mProton));
	  if (checkSliceContent( NclsVsBG[tracks.trdNCh()-4], tracks.trdNCls(), betaGamma(momentumMean(tracks.trdMom()), mProton)))  
	    dEdxVsNclsVsBG[tracks.trdNCh()-4]->Fill(tracks.trdNCls(), betaGamma(momentumMean(tracks.trdMom()), mProton), CorrectionFactor*tracks.trdSig());
	  // if (checkSliceContent(dEdxVsNclsVsEtaFull, tracks.trdNCls(), tracks.trdTPCtgl()))
	  //dEdxVsNclsVsEtaFull->Fill(tracks.trdNCls(), tracks.trdTPCtgl(), CorrectionFactor*tracks.trdSig());
	  //if (checkSliceContent( dEdxVsNclsVsEta[tracks.trdNCh()-4], tracks.trdNCls(), tracks.trdTPCtgl()))  
	  // dEdxVsNclsVsEta[tracks.trdNCh()-4]->Fill(tracks.trdNCls(),tracks.trdTPCtgl(), CorrectionFactor*tracks.trdSig());
	  CentVsBG->Fill(tracks.centrality(), betaGamma(momentumMean(tracks.trdMom()), mProton));
       	  if (checkSliceContent( CentVsBG, tracks.centrality(), betaGamma(momentumMean(tracks.trdMom()), mProton)))  
	    dEdxVsCentVsBG->Fill(tracks.centrality(), betaGamma(momentumMean(tracks.trdMom()), mProton), CorrectionFactor*tracks.trdSig());
        }

    }
    nchMin=nchTmp;

    // declare list
    TList list;
    TList listTPCtgl;
    TList listNcls;
    TList meanList;
    TList centList;
    
    // declare temporary histograms
    TH1D   * tmpHist;
    //TF1     * tmpFit;

    // declare temporary mean-val.
    Double_t tmpHistMean;
    //Double_t tmpHistError;
    //Double_t tmpMax;
    
    //Double_t fitRangeLow;
    //Double_t fitRangeHigh


    // calculation overall signal mean for correction (instead of correction to fit)
  
    // Mean Values
    for (Int_t xBin = BGVsdEdxMean->GetXaxis()->GetFirst(); xBin <= BGVsdEdxMean->GetNbinsX(); xBin++) {
      tmpHist=BGVsdEdxMean->ProjectionY("tmpHist", xBin, xBin, "");
      tmpHist->SetTitle(Form("%s Mean, (%i,%i)", BGVsdEdxMean->GetName(), xBin, xBin));
      tmpHist->SetName(Form("%s_%i_%i", BGVsdEdxMean->GetName(), xBin, xBin));

      // discard low statistics
      if (tmpHist->Integral() < 10) continue;     
      tmpHistMean = tmpHist->GetMean();
      
      // set mean value
      if (tmpHistMean != 0) {
	BGMean->SetBinContent(xBin, tmpHistMean);
      }
      
      meanList.Add(tmpHist);
    }
    meanList.Add(BGMean);
   
    
    // calculate slice mean for NCls Correction
    Int_t Average, AverageChamber;
    Double_t TmpAverage, TmpAverageChamber;
   
    // for (int NCh=0; NCh<3; NCh++) {
    for (Int_t yBin = dEdxVsNclsVsBG[0]->GetYaxis()->GetFirst(); yBin <= dEdxVsNclsVsBG[0]->GetNbinsY(); yBin++) { // for each #cls
      Average=0;
      TmpAverage=0;
      for (int NCh=0; NCh<3; NCh++) { // for each cahmber
	TmpAverageChamber=0;
	AverageChamber=0;
	//for (Int_t xBin = dEdxVsNclsVsBG[NCh]->GetXaxis()->GetFirst()+30+NCh*20; xBin <= dEdxVsNclsVsBG[NCh]->GetXaxis()->GetFirst()+40+NCh*20; xBin++) {
	for (Int_t xBin = dEdxVsNclsVsBG[NCh]->GetXaxis()->GetFirst(); xBin <= dEdxVsNclsVsBG[NCh]->GetNbinsX(); xBin++) { // for each beta*gamma slices
	  //	  for (Int_t xBin = dEdxVsNclsVsBG[NCh]->GetXaxis()->GetFirst(); xBin <= dEdxVsNclsVsBG[NCh]->GetNbinsX(); xBin++) {
	  // project along z
	  tmpHist = dEdxVsNclsVsBG[NCh]->ProjectionZ("tmpHist", xBin, xBin, yBin, yBin);
	  tmpHist->SetTitle(Form("%s projection, (%i,%i)", dEdxVsNclsVsBG[NCh]->GetName(), xBin, yBin));
	    
	  // discard low statistics
	  if (tmpHist->Integral() < 10) continue;
	  
	  TmpAverageChamber += tmpHist->GetMean();
	  AverageChamber++;
	    
	}
	if  (AverageChamber>20) {
	  TmpAverage+=TmpAverageChamber/Double_t(AverageChamber);
	  Average++;
	}
	  
      }
      if  (Average==nchMax-nchMin+1) {
	TmpAverage /= Double_t(Average);
	//cout << yBin <<  " " << TmpAverage << endl;
	BGMeanNCls[0]->SetBinContent(yBin, TmpAverage);
      }
    }
    meanList.Add(BGMeanNCls[0]);

    Int_t AverageSlices;
    // TEST Improved ETAAsymmetryMap with average over slices
    for (Int_t yBin = dEdxVsTPCtglVsBG->GetYaxis()->GetFirst(); yBin <= dEdxVsTPCtglVsBG->GetNbinsY(); yBin++) {
      AverageSlices=0;
      TmpAverage=0;
      for (Int_t xBin = dEdxVsTPCtglVsBG->GetXaxis()->GetFirst()+8; xBin <= dEdxVsTPCtglVsBG->GetNbinsX()-8; xBin++) {
	// project along z
	tmpHist = dEdxVsTPCtglVsBG->ProjectionZ("tmpHist", xBin, xBin, yBin, yBin);
	tmpHist->SetTitle(Form("%s projection, (%i,%i)", dEdxVsTPCtglVsBG->GetName(), xBin, yBin));
            
	// discard low statistics
	if (tmpHist->Integral() < 10) continue;

	TmpAverage += tmpHist->GetMean();
        AverageSlices++;
	
      }
      if  (AverageSlices>120) {
	TmpAverage /= AverageSlices;
	//cout << yBin <<  " " << TmpAverage << endl;
	BGMeanEta->SetBinContent(yBin, TmpAverage);
      }
    }
    meanList.Add(BGMeanEta);

  
    

    printf("ClusterMap(): producing map for NCls vs  TPC tgl variable \n");
    for (Int_t xBin = dEdxVsNclsVsEtaFull->GetXaxis()->GetFirst(); xBin <= dEdxVsNclsVsEtaFull->GetNbinsX(); xBin++) {
      for (Int_t yBin = dEdxVsNclsVsEtaFull->GetYaxis()->GetFirst(); yBin <= dEdxVsNclsVsEtaFull->GetNbinsY(); yBin++) {
	// project along z
	for (Int_t nCh=0; nCh <3; nCh++) {
	  tmpHist = dEdxVsNclsVsEta[nCh]->ProjectionZ("tmpHist", xBin, xBin, yBin, yBin);
	  tmpHist->SetTitle(Form("%s projection, (%i,%i)", dEdxVsNclsVsEta[nCh]->GetName(), xBin, yBin));
          
	  // discard low statistics
	  if (tmpHist->Integral() < 10) continue;
	  tmpHistMean = tmpHist->GetMean();
	  // set mean value
	  if (tmpHistMean != 0) {
	    EtaNclsMap[nCh]->SetBinContent(xBin, yBin, tmpHistMean);
	    //map->SetBinError(   xBin, yBin, tmpHistError);
	  }       
	  // add temporary histogram to list+
	  tmpHist->SetName(Form("%s_%i_%i", dEdxVsNclsVsEta[nCh]->GetName(), xBin, yBin));
	  listNcls.Add(tmpHist);
	}
      }
    }
    
    printf("ClusterMap(): producing corr. map for new NCls variable \n");
    for (Int_t xBin = dEdxVsNclsVsBGFull->GetXaxis()->GetFirst(); xBin <= dEdxVsNclsVsBGFull->GetNbinsX(); xBin++) {
      for (Int_t yBin = dEdxVsNclsVsBGFull->GetYaxis()->GetFirst(); yBin <= dEdxVsNclsVsBGFull->GetNbinsY(); yBin++) {
	// project along z
	for (Int_t nCh=0; nCh <3; nCh++) {
	  
	  tmpHist = dEdxVsNclsVsBG[nCh]->ProjectionZ("tmpHist", xBin, xBin, yBin, yBin);
	  tmpHist->SetTitle(Form("%s projection, (%i,%i)", dEdxVsNclsVsBG[nCh]->GetName(), xBin, yBin));
          
	  // discard low statistics
	  if (tmpHist->Integral() < 10) continue;
	  
	  tmpHistMean = tmpHist->GetMean();
          
	  // set correction factor
	  if (tmpHistMean != 0) {
	    if (CorrectToMean) {
	      NclsMap[nCh]->SetBinContent(xBin, yBin, expMPV(BGMean, NclsMap[nCh]->GetYaxis()->GetBinCenter(yBin))/tmpHistMean); // correction to ClusterSliceMean
	    }
	    else {
	      NclsMap[nCh]->SetBinContent(xBin, yBin, expMPV(meanParam, NclsMap[nCh]->GetYaxis()->GetBinCenter(yBin))/tmpHistMean); // correction to fit
	    }
	    //map->SetBinError(   xBin, yBin, tmpHistError);
	  }
          
	  // add temporary histogram to list+
	  tmpHist->SetName(Form("%s_%i_%i", dEdxVsNclsVsBG[nCh]->GetName(), xBin, yBin));
	  listNcls.Add(tmpHist);
	}
      }
    }
    

    printf("CentralityMap(): producing corr. map for centrality dependence  \n");
    //Correct to cent 50-100 bin
    for (Int_t yBin = dEdxVsCentVsBG->GetYaxis()->GetFirst(); yBin <= dEdxVsCentVsBG->GetNbinsY(); yBin++) {
      // project along z
      tmpHist = dEdxVsCentVsBG->ProjectionZ("tmpHist", NCentBins-1, NCentBins-1, yBin, yBin);
      tmpHist->SetTitle(Form("%s projection, (%i,%i)", dEdxVsCentVsBG->GetName(), NCentBins, yBin));
      // discard low statistics
      if (tmpHist->Integral() < 100) {
	cout << "Error: to low statistic in reference bin" << endl;
	continue;
      }
      tmpHistMean=tmpHist->GetMean();
      CentMapRef->SetBinContent(yBin, tmpHistMean);
    }


    for (Int_t xBin = dEdxVsCentVsBG->GetXaxis()->GetFirst(); xBin <= dEdxVsCentVsBG->GetNbinsX(); xBin++) {
      for (Int_t yBin = dEdxVsCentVsBG->GetYaxis()->GetFirst(); yBin <= dEdxVsCentVsBG->GetNbinsY(); yBin++) {
	// project along z
	tmpHist = dEdxVsCentVsBG->ProjectionZ("tmpHist", xBin, xBin, yBin, yBin);
	tmpHist->SetTitle(Form("%s projection, (%i,%i)", dEdxVsCentVsBG->GetName(), xBin, yBin));
        
	// discard low statistics
	if (tmpHist->Integral() < 10) continue;
	
	tmpHistMean = tmpHist->GetMean();
	
	// set correction factor
	if (tmpHistMean != 0) {
	  // if (CorrectToMean) {
	  CentMap->SetBinContent(xBin, yBin, expMPV(CentMapRef , CentMap->GetYaxis()->GetBinCenter(yBin))/tmpHistMean);
	  //      }
	  // else {
	  //	CentMap->SetBinContent(xBin, yBin, expMPV(meanParam, mapTPCtgl->GetYaxis()->GetBinCenter(yBin))/tmpHistMean); // correct to fit
	      // }
	}
        
	// add temporary histogram to list+
	tmpHist->SetName(Form("%s_%i_%i", dEdxVsCentVsBG->GetName(), xBin, yBin));
	centList.Add(tmpHist);
      }
    }
    centList.Add(CentMap);
    centList.Add(CentMapRef);
    centList.Add(CentVsBG);
    centList.Add(dEdxVsCentVsBG);



    // only for iterative procedures
    /*  
    if (NclsCorrection)
    {
    for (int i =0; i<3; i++) {
    NclsMap[i]->Multiply(NclsMaps2[i]);
    }
    }
    */
 



    printf("etaAsymmetryMap(): producing corr. map for TPCtgl variable \n");
    for (Int_t xBin = dEdxVsTPCtglVsBG->GetXaxis()->GetFirst(); xBin <= dEdxVsTPCtglVsBG->GetNbinsX(); xBin++) {
      for (Int_t yBin = dEdxVsTPCtglVsBG->GetYaxis()->GetFirst(); yBin <= dEdxVsTPCtglVsBG->GetNbinsY(); yBin++) {
	// project along z
	tmpHist = dEdxVsTPCtglVsBG->ProjectionZ("tmpHist", xBin, xBin, yBin, yBin);
	tmpHist->SetTitle(Form("%s projection, (%i,%i)", dEdxVsTPCtglVsBG->GetName(), xBin, yBin));
        
	// discard low statistics
	if (tmpHist->Integral() < 10) continue;
       	tmpHistMean = tmpHist->GetMean();
            
	// set correction factor
	if (tmpHistMean != 0) {
	  if (CorrectToMean) {
	    mapTPCtgl->SetBinContent(xBin, yBin, expMPV(BGMeanEta  , mapTPCtgl->GetYaxis()->GetBinCenter(yBin))/tmpHistMean); // correct to overall mean oder slice mean (BGMeanEta)
	  }
	  else {
	    mapTPCtgl->SetBinContent(xBin, yBin, expMPV(meanParam, mapTPCtgl->GetYaxis()->GetBinCenter(yBin))/tmpHistMean); // correct to fit
	  }
	  //map->SetBinError(   xBin, yBin, tmpHistError);
	}
        
	// set tmpHist and tmpFit names
	//tmpHist->SetName(   Form("%s_%i_%i", dEdxVsEtaVsBG->GetName(), xBin, yBin));
	//tmpFit->SetName(    Form("%s_%i_%i_fit", dEdxVsEtaVsBG->GetName(), xBin, yBin));
        
	// add temporary histogram to list+
	tmpHist->SetName(Form("%s_%i_%i", dEdxVsTPCtglVsBG->GetName(), xBin, yBin));
	listTPCtgl.Add(tmpHist);
	//list.Add(tmpFit);

      }
    }


    /*
    printf("etaAsymmetryMap(): producing eta corr. map\n");
    for (Int_t xBin = dEdxVsEtaVsBG->GetXaxis()->GetFirst(); xBin <= dEdxVsEtaVsBG->GetNbinsX(); xBin++) {
        for (Int_t yBin = dEdxVsEtaVsBG->GetYaxis()->GetFirst(); yBin <= dEdxVsEtaVsBG->GetNbinsY(); yBin++) {
            
	  // project along z
	  tmpHist = dEdxVsEtaVsBG->ProjectionZ("tmpHist", xBin, xBin, yBin, yBin);
	  tmpHist->SetTitle(Form("%s projection, (%i,%i)", dEdxVsEtaVsBG->GetName(), xBin, yBin));
            
	  // discard low statistics
	  if (tmpHist->Integral() < 10) continue;
	    
	  /* old parts!!!
	  // get fit range
	  tmpMax = tmpHist->GetMaximum();a
	  for (Int_t ii = tmpHist->GetXaxis()->GetFirst(); ii <= tmpHist->GetXaxis()->GetLast(); ii++) {
	  if (tmpHist->GetBinContent(ii) >= tmpMax*0.25) {
	  fitRangeLow = tmpHist->GetXaxis()->GetBinLowEdge(ii);
	  break;
	  }
	  }
            
	  for (Int_t ii = tmpHist->GetXaxis()->GetLast(); ii >= tmpHist->GetXaxis()->GetFirst(); ii--) {
	  if (tmpHist->GetBinContent(ii) >= tmpMax*0.25) {
	  fitRangeHigh = tmpHist->GetXaxis()->GetBinUpEdge(ii);
	  break;
	  }
	  }
            
	  // define fit function
	  tmpFit = new TF1("tmpFit", "[0]+[1]*TMath::Gaus(x, [2], [3])", fitRangeLow, fitRangeHigh);
            
	  // initialize fit-parameter
	  tmpFit->SetParameter(0, tmpHist->GetMinimum());
	  tmpFit->SetParameter(1, tmpHist->GetMaximum());
	  tmpFit->SetParameter(2, tmpHist->GetBinCenter(tmpHist->GetMaximumBin()));
	  tmpFit->SetParameter(3, tmpHist->GetRMS());
          
	  // fit tmpHist
	  tmpHist->Fit("tmpFit", "MRQ");
          
	  // get mpv and error from fit
	  tmpHistMean     = tmpFit->GetParameter( 2);
	  tmpHistError    = tmpFit->GetParError(  2);
            
	  // discard "bad" fits
	  if (tmpHistMean < tmpHist->GetXaxis()->GetBinLowEdge(tmpHist->GetXaxis()->GetFirst()) || tmpHistMean > tmpHist->GetXaxis()->GetBinUpEdge(tmpHist->GetXaxis()->GetLast())) continue;
	  //if (tmpHistError > 10e-2) continue;
	  
            
	  tmpHistMean = tmpHist->GetMean();
            
	  // set correction factor
	  if (tmpHistMean != 0) {
	    if (CorrectToMean) {
	      map->SetBinContent( xBin, yBin, expMPV(BGMean, map->GetYaxis()->GetBinCenter(yBin))/tmpHistMean);
	    }
	    else {
	      map->SetBinContent( xBin, yBin, expMPV(meanParam, map->GetYaxis()->GetBinCenter(yBin))/tmpHistMean);
	    }
	    //map->SetBinError(   xBin, yBin, tmpHistError);
	  }
          
	  // set tmpHist and tmpFit names
	  //tmpHist->SetName(   Form("%s_%i_%i", dEdxVsEtaVsBG->GetName(), xBin, yBin));
	  //tmpFit->SetName(    Form("%s_%i_%i_fit", dEdxVsEtaVsBG->GetName(), xBin, yBin));
          
	  // add temporary histogram to list+
	  tmpHist->SetName(Form("%s_%i_%i", dEdxVsEtaVsBG->GetName(), xBin, yBin));
	  list.Add(tmpHist);
	  //list.Add(tmpFit);
        }
    }
    
*/


    TString tmpNcls(NclsCorr);
    TString tmpEta(EtaSpec);


    if (Iteration==2) {
      //map->Multiply(etaMap2);
      for (int i=0; i<3; i++) NclsMap[i]=NclsMaps2[i];
      mapTPCtgl->Multiply(etaMap2);
      //NclsCorr="";
      EtaSpec="";
    }
    if (Iteration==1) {
      //map->Multiply(etaMap2);
      for (int i=0; i<3; i++) NclsMap[i]->Multiply(NclsMaps2[i]);
      mapTPCtgl=etaMap2;
      EtaSpec=""; 
      NclsCorr="";
    }




    //
    // save output
    //

    printf("Centrality Map: writing output\n");
    TFile Centfile(Form("output/%s_CentFile_%s_%s%s%s.root", inputSpec.Data(),  FileSpec.Data(),CentCorr.Data(), EtaSpec.Data(), NclsCorr.Data()), "RECREATE");
    centList.Write();
    Centfile.Close();
    Centfile.Save();

    TFile centMap(Form("output/%s_CentMap_%s_%s%s%s.root", inputSpec.Data(),   FileSpec.Data(), CentCorr.Data(),EtaSpec.Data(), NclsCorr.Data()), "RECREATE");
    CentMap->Write();
    centMap.Close();
    centMap.Save();

    /*
    printf("etaAsymmetryMap(): writing output for eta variable\n");
    // TFile file(Form("output/%s_dedxasym_%incls.root", inputSpec.Data(), ncls), "RECREATE");
    TFile file(Form("output/%s_EtaFile_%s_%s%s%s.root", inputSpec.Data(),  FileSpec.Data(), CentCorr.Data(),EtaSpec.Data(), NclsCorr.Data()), "RECREATE");
    dEdxVsEtaVsBG->Write();
    map->Write();
    list.Write();
    file.Close();
    file.Save();
    
    // TFile etaMap(Form("output/%s_etaMap_%incls.root", inputSpec.Data(), ncls), "RECREATE");
    TFile etaMap(Form("output/%s_EtaMap_%s_%s%s%s.root", inputSpec.Data(),  FileSpec.Data(), CentCorr.Data(),EtaSpec.Data(), NclsCorr.Data()), "RECREATE");
    map->Write();
    etaMap.Close();
    etaMap.Save();
    */

    printf("etaAsymmetryMap(): writing output for TPCtgl variable\n");
    // TFile fileTPCtgl(Form("output/%s_dedxasym_%incls_forTPCtgl.root", inputSpec.Data(), ncls), "RECREATE");
    TFile fileTPCtgl(Form("output/%s_TPCtglFile_%s_%s%s%s.root", inputSpec.Data(), FileSpec.Data(),CentCorr.Data(), EtaSpec.Data(), NclsCorr.Data()), "RECREATE");
    dEdxVsTPCtglVsBG->Write();
    TPCtglVsBG->Write();
    mapTPCtgl->Write();
    listTPCtgl.Write();    
    BGVsdEdxMean->Write();
    BGMean->Write();
    meanList.Write();
    fileTPCtgl.Close();
    fileTPCtgl.Save();
    
    // TFile etaMapTPCtgl(Form("output/%s_etaMap_%incls_forTPCtgl.root", inputSpec.Data(), ncls), "RECREATE");
    TFile etaMapTPCtgl(Form("output/%s_TPCtglMap_%s_%s%s%s.root", inputSpec.Data(),  FileSpec.Data(), CentCorr.Data(),EtaSpec.Data(), NclsCorr.Data()), "RECREATE");
    mapTPCtgl->Write();
    meanList.Write();
    etaMapTPCtgl.Close();
    etaMapTPCtgl.Save();

    printf("etaAsymmetryMap(): writing output for new NCls variable\n");
    // TFile fileTPCtgl(Form("output/%s_dedxasym_%incls_forTPCtgl.root", inputSpec.Data(), ncls), "RECREATE");
    TFile fileNcls(Form("output/%s_ClusterFile_%s_%s%s%s.root", inputSpec.Data(),  FileSpec.Data(), CentCorr.Data(),EtaSpec.Data(), NclsCorr.Data()), "RECREATE");
    dEdxVsNclsVsBGFull->Write();
    dEdxVsNclsVsEtaFull->Write();
    for (Int_t i =0; i<3; i++) {
      dEdxVsNclsVsEta[i]->Write();
      EtaNclsMap[i]->Write();
      
      dEdxVsNclsVsBG[i]->Write();
      NclsMap[i]->Write();

      NclsVsBG[i]->Write();
    }
    listNcls.Write();
    fileNcls.Close();
    fileNcls.Save();
    
    // TFile etaMapTPCtgl(Form("output/%s_etaMap_%incls_forTPCtgl.root", inputSpec.Data(), ncls), "RECREATE");
    TFile nclsMap(Form("output/%s_ClusterMap_%s_%s%s%s.root", inputSpec.Data(),  FileSpec.Data(), CentCorr.Data(), EtaSpec.Data(), NclsCorr.Data()), "RECREATE");
    for (Int_t i = 0; i<3; i++) {
      NclsMap[i]->Write();
    }
    nclsMap.Close();
    nclsMap.Save();

    // TFile etaMapTPCtgl(Form("output/%s_etaMap_%incls_forTPCtgl.root", inputSpec.Data(), ncls), "RECREATE");
    //    TFile TestmapFIle(Form("output/Testmap.root"), "RECREATE");
    //    Testmap->Write();
    //    TestmapFIle.Close();
    //    TestmapFIle.Save();

    NclsCorr= tmpNcls;
    EtaSpec = tmpEta;

    delete dEdxVsEtaVsBG, dEdxVsEtaVsBG, mapTPCtgl, map, TPCtglVsBG, Testmap, BGVsdEdxMean, BGMean, BGMeanEta;
    delete dEdxVsNclsVsBGFull, dEdxVsNclsVsEtaFull;
    for (int i=0; i<3; i++)  delete dEdxVsNclsVsBG[i], NclsVsBG[i], NclsMap[i], dEdxVsNclsVsEta[i], dEdxVsNclsVsEta[i], BGMeanNCls[i];

}

        
