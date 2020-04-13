//====================================================================================================================================
//
// root macro to paramtrize the dE/dx resolution for e, pi, p (using nch = 4 - 6)
//  THIS MACRO NEEDS THE SIGNAL PARAMETRIZATION FOR nch = 4 to nch = 6 -> outdated. Only common parametrization for all chambers available
//
//====================================================================================================================================


#include <iostream>
#include <ctime>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TMinuit.h>
#include <TList.h>
#include <TF1.h>
#include <TStyle.h>
//#include "V0Tracks.cxx"
#include "constants.h"
#include "utils.h"


void resolutionIncl() {
    // if no charge seperation is wanted, set tracks.pId() to abs
    Int_t pdgCodes[3] = {11, 211, 2212};
    
    // get and initialize tree
    V0Tracks tracks;
    tracks.setFile(input.Data());
    tracks.getTree();
    
    // get number of entries from tree
    const Int_t nEntries = tracks.getNumberOfEntries();
    
    // create new TList objects for resolution parametrization
    TList * resList         = new TList;
    TList * resSlicesList   = new TList;

    TList * resListBG         = new TList;
    TList * resSlicesListBG   = new TList;

    TList * resListTPCtgl         = new TList;
    TList * resSlicesListTPCtgl   = new TList;
    
    
    // declare histograms for resolution parametrization
    TH2D * resSigVsTRDncls  = new TH2D("resSigVsTRDncls",   "resSigVsTRDncls", 120, 40, 160, 200, 0, 2);
    TH1D * resNor,      * resMPV,       * resWidth,     * resRes,       * resChi;
    TH2D * SigParamVsTRDncls = new TH2D("SigParamVsTRDncls", "SigParamVsTRDncls", 120, 40, 160, 200, 0, 20);
    TH2D * SigVsTRDncls = new TH2D("SigVsTRDncls", "SigVsTRDncls", 120, 40, 160, 200, 0, 20);

    TH2D * resSigVsTRDnch = new TH2D("resSigVsTRDnch", "resSigVsTRDnch", 3,4,7,200,0,2);
    TH1D * resNorVsnch, * resMPVVsnch, * resWidthVsnch, * resResVsnch, * resChiVsnch;

    TH2D * resSigVsBG       = new TH2D("resSigVsBG",        "resSigVsBG",       100, 0.3,   1e4, 200, 0, 2);
    TH1D * resNorVsBG,  * resMPVVsBG,   * resWidthVsBG, * resResVsBG,   * resChiVsBG;

    TH2D * resSigVsTPCtgl       = new TH2D("resSigVsTPCtgl",        "resSigVsTPCtgl",       180, -1,  1, 200, 0, 2);
    TH1D * resNorVsTPCtgl,  * resMPVVsTPCtgl,   * resWidthVsTPCtgl, * resResVsTPCtgl,   * resChiVsTPCtgl;
    /*
    TH2D * resSigVstrdGlobalPhi      = new TH2D("resSigVstrdGlobalPhi",        "resSigVstrdGlobalPhi",       180, 0, 2*TMath::Pi(), 200, 0, 2);
    TH1D * resNorVstrdGlobalPhi,  * resMPVVstrdGlobalPhi,   * resWidthVstrdGlobalPhi, * resResVstrdGlobalPhi,   * resChiVstrdGlobalPhi;
    
    TH2D * resSigVstrdTheta       = new TH2D("resSigVstrdTheta",        "resSigVstrdTheta",       180, -1, 1, 200, 0, 2);
    TH1D * resNorVstrdTheta,  * resMPVVstrdTheta,   * resWidthVstrdTheta, * resResVstrdTheta,   * resChiVstrdTheta;

    TH2D * resSigVstrdY      = new TH2D("resSigVstrdY",        "resSigVstrdY",       100, 0,  100, 200, 0, 2);
    TH1D * resNorVstrdY,  * resMPVVstrdY,   * resWidthVstrdY, * resResVstrdY,   * resChiVstrdY;

    TH2D * resSigVstrdEtaLocal      = new TH2D("resSigVstrdEtaLocal",        "resSigVstrdEtaLocal",       180, -1,  1, 200, 0, 2);
    TH1D * resNorVstrdEtaLocal,  * resMPVVstrdEtaLocal,   * resWidthVstrdEtaLocal, * resResVstrdEtaLocal,   * resChiVstrdEtaLocal;

    TH2D * resSigVstrdPhiLocal      = new TH2D("resSigVstrdPhiLocal",        "resSigVstrdPhiLocal",       180, 0, 2*TMath::Pi(), 200, 0, 2);
    TH1D * resNorVstrdPhiLocal,  * resMPVVstrdPhiLocal,   * resWidthVstrdPhiLocal, * resResVstrdPhiLocal,   * resChiVstrdPhiLocal;

    */
    // further control histograms
    TH2D * nclsVsBG = new TH2D("nclsVsBG", "nclsVsBG", 100, 0.3, 1e4, 120, 40, 160);
    TH2D * TPCtglVsBG = new TH2D("TPCtglVsBG", "TPCtglVsBG", 100, 0.3, 1e4, 100, -1., 1.);
    TH2D * TPCtglVsncls = new TH2D("TPCtglVsncls", "TPCtglVsncls", 120, 40, 160, 100, -1, 1); 
    TH2D * TPCtglVsncls_ele = new TH2D("TPCtglVsncls_ele", "TPCtglVsncls_ele", 120, 40, 160, 100, -1, 1);
    TH2D * TPCtglVsncls_pion = new TH2D("TPCtglVsncls_pion", "TPCtglVsncls_pion", 120, 40, 160, 100, -1, 1);
    TH2D * TPCtglVsncls_proton = new TH2D("TPCtglVsncls_proton", "TPCtglVsncls_proton", 120, 40, 160, 100, -1, 1);
    TH2D * resSigVsTRDncls_ele = new TH2D("resSigVsTRDncls_ele", "resSigVsTRDncls_ele", 120, 40, 160, 200, 0, 2);
    TH2D * resSigVsTRDncls_pion = new TH2D("resSigVsTRDncls_pion", "resSigVsTRDncls_pion",120, 40, 160, 200, 0, 2);
    TH2D * resSigVsTRDncls_proton = new TH2D("resSigVsTRDncls_proton", "resSigVsTRDncls_proton", 120, 40, 160, 200, 0, 2);
    TH2D * resSigVsTRDtpctgl_ele = new TH2D("resSigVsTRDtpctgl_ele", "resSigVsTRDtpctgl_ele", 180, -1, 1, 200, 0, 2);
    TH2D * resSigVsTRDtpctgl_pion = new TH2D("resSigVsTRDtpctgl_pion", "resSigVsTRDtpctgl_pion",180, -1, 1, 200, 0, 2);
    TH2D * resSigVsTRDtpctgl_proton = new TH2D("resSigVsTRDtpctgl_proton", "resSigVsTRDtpctgl_proton", 180, -1, 1, 200, 0, 2);
    TH1D * resNorVsTRDtpctgl_proton,  * resMPVVsTRDtpctgl_proton,   * resWidthVsTRDtpctgl_proton, * resResVsTRDtpctgl_proton,   * resChiVsTRDtpctgl_proton; 
    TH1D * resNorVsTRDtpctgl_pion,  * resMPVVsTRDtpctgl_pion,   * resWidthVsTRDtpctgl_pion, * resResVsTRDtpctgl_pion,   * resChiVsTRDtpctgl_pion;
    TH1D * resNorVsTRDtpctgl_ele,  * resMPVVsTRDtpctgl_ele,   * resWidthVsTRDtpctgl_ele, * resResVsTRDtpctgl_ele,   * resChiVsTRDtpctgl_ele;
    TH1D * resNorVsTRDncls_ele,  * resMPVVsTRDncls_ele,   * resWidthVsTRDncls_ele, * resResVsTRDncls_ele,   * resChiVsTRDncls_ele;  
    TH1D * resNorVsTRDncls_pion,  * resMPVVsTRDncls_pion,   * resWidthVsTRDncls_pion, * resResVsTRDncls_pion,   * resChiVsTRDncls_pion;  
    TH1D * resNorVsTRDncls_proton,  * resMPVVsTRDncls_proton,   * resWidthVsTRDncls_proton, * resResVsTRDncls_proton,   * resChiVsTRDncls_proton;
    TH2D * TRDCentralityWOCuts = new TH2D("TRDCentralityWOCuts", "TRDCentralityWOCuts",NCentBins,  CentBins, 200, 0, 2 );
    TH2D * TRDCentralityWCuts = new TH2D("TRDCentralityWCuts", "TRDCentralityWCuts", NCentBins,  CentBins, 200, 0, 2 );
    TH1D * resNorVsTRDcent,  * resMPVVsTRDcent,   * resWidthVsTRDcent, * resResVsTRDcent,   * resChiVsTRDcent;

    // load corrections maps -- see comments in signal.cxx
  
    // TPCtgl corr. map - NOTICE: Not really eta map - instead closely related TPCtgl variable
    TString EtaSpec, NclsCorr, CentCorr;
    if (!NclsCorrection) NclsCorr="";
    else NclsCorr=NclsCorrStr;
    if (!fEtaCorrection) EtaSpec ="";
    else EtaSpec = EtaSpecStr;
    if (!CentCorrection) CentCorr="";
    else CentCorr=CentCorrStr;


    // loading correction maps. slightly confusing structure at the moment due to iterative structure
    TH2D * etaMap=0;
    if (fEtaCorrection) {
      //TString tmp(NclsCorr);
      //NclsCorr="";
      TFile * etaMapFile = TFile::Open(Form("output/%s_TPCtglMap_%s_%s%s.root", inputSpec.Data(), FileSpec.Data(),CentCorr.Data(), NclsCorr.Data()), "READ");
      //etaMap = (TH2D*)etaMapFile->Get("map");
      etaMap = (TH2D*)etaMapFile->Get("mapTPCtgl");
      //NclsCorr=tmp;
    }
   
    TH2D * NclsMap[3]={0,0,0};
    if (NclsCorrection) {
      TString tmp(EtaSpec);
      EtaSpec="";
      TFile * ClusterMapFile = TFile::Open(Form("output/%s_ClusterMap_%s_%s%s.root", inputSpec.Data(),  FileSpec.Data(), CentCorr.Data(), EtaSpec.Data()), "READ");
      for (int i=0; i<3; i++) {
	NclsMap[i] = (TH2D*)ClusterMapFile->Get(Form("NclsMap_nCh%i", i+4));
      }
      EtaSpec=tmp;
    }
    
    TH2D * CentMap=0;
    if (CentCorrection) {
      TString tmp(EtaSpec);
      TString tmp2(NclsCorr);
      EtaSpec="";
      NclsCorr="";
      TFile * CentMapFile = TFile::Open(Form("output/%s_CentMap_%s_%s%s.root", inputSpec.Data(),  FileSpec.Data(), EtaSpec.Data(), NclsCorr.Data()), "READ");
      CentMap = (TH2D*)CentMapFile->Get("CentMap");
      EtaSpec=tmp;
      NclsCorr=tmp2;
    }

    if (!CentMap) cout << "No CentMap" << endl;
    if (!etaMap) cout << "No EtaMap" << endl;
    if (!NclsMap[0] || !NclsMap[1] || !NclsMap[2]) cout << "ClusterMaps are missing" << endl;
    
    //===========================================================================================
    //
    // mean parametrization
    //
    //===========================================================================================
    
    // mean parametrizations - old version for separation using different NCh
    /*  TFile * meanParam[3];
    for (Int_t i = 2; i < 3; i++) {
        meanParam[i] = TFile::Open(Form("output/%s_meanParam_pid%i_%inch_%incls_%s%s.root", inputSpec.Data(), pidOpt, 4+i, ncls, FileSpec.Data(), EtaSpec.Data()));
    }
    
    TH1D * meanFitHist[3];
    for (Int_t i = 2; i < 3; i++) {
        meanFitHist[i] = (TH1D*)meanParam[i]->Get("meanFitHist");
    }
    */
    
    TFile* meanParam=0;
    meanParam = TFile::Open(Form("output/%s_meanParam_pid%i_%s%s%s%s.root", inputSpec.Data(), pidOpt, FileSpec.Data(), CentCorr.Data(), EtaSpec.Data(), NclsCorr.Data()));
    if (!meanParam) {
      cout << "Error: meanParam file not found" << endl;
      return;
    }

    TDirectoryFile *meanParamBasic =0;
    meanParamBasic = (TDirectoryFile*)meanParam->Get("Basic");
    if (!meanParamBasic) {
      cout << "Error: FileFolder Basic does not exist" << endl;
      return;
    }

    TH1D* meanFitHist=0;
    meanFitHist = (TH1D*)meanParamBasic->Get("meanFitHist");
    if (!meanFitHist) {
      cout << "Error: meanFitHist not found" << endl;
      return;
    }
    
    //===========================================================================================
    //
    // resolution parametrization
    //
    //===========================================================================================
    
    // prepare signal vs beta*gamma histogram - set log scale
    BinLogX(resSigVsBG->GetXaxis());
    BinLogX(nclsVsBG->GetXaxis());
    BinLogX(TPCtglVsBG->GetXaxis());

    
    // loop over tree to fill signal vs beta*gamma and signal vs trd ncls histogram
    Int_t nb = 0;
    Int_t run= 0;
    cout << "nEntries of Tree: " <<nEntries << endl;
    for (Int_t i = 0; i < nEntries; i++) {
        if (i%1000000==0) cout << i << endl;
        
	nb = tracks.getEntry(i);
        if (nb <= 0) continue;

	if (runCut(tracks.run())) continue;

	// (tracks.centrality()<3) continue;

	//PbPbTest of centrality
	if (TMath::Abs(tracks.pId()) == pdgCodes[0]) {
	  TRDCentralityWOCuts->Fill(tracks.centrality(), tracks.trdSig()/expMPV(meanFitHist, betaGamma(momentumMean(tracks.trdMom()), mElectron)));
	}
	if (TMath::Abs(tracks.pId()) == pdgCodes[1]) {
	  TRDCentralityWOCuts->Fill(tracks.centrality(), tracks.trdSig()/expMPV(meanFitHist, betaGamma(momentumMean(tracks.trdMom()), mPion)));
	}
	if (TMath::Abs(tracks.pId()) == pdgCodes[2]) {
	  TRDCentralityWOCuts->Fill(tracks.centrality(), tracks.trdSig()/expMPV(meanFitHist, betaGamma(momentumMean(tracks.trdMom()), mProton)));
	}

	if (TMath::Abs(tracks.trdTPCtgl())>1.) continue; // no correction for higher eta
	// similar cuts for cluster or centrality?
        
	// cut for number of tracklets, cluster
	if (trackAndLayerCuts(tracks.trdNCh(), tracks.trdNCls(),  TMath::Abs(tracks.pId()))) continue;

	// geometric cuts: LocalTheta, YCut - at the moment empty as no access to this data on AOD level
	if (geometricCuts(tracks.trdEtaLocal(), tracks.trdNCh(), tracks.trdY())) continue;
       
        // apply additional particle identification cuts on tracks
	if (nSigmaCut(tracks.nSigmaTPC(), (Double_t*)cutNSigmaTPC, TMath::Abs(tracks.pId()))) continue;
        if (nSigmaCut(tracks.nSigmaTOF(), (Double_t*)cutNSigmaTOF, TMath::Abs(tracks.pId()))) continue;
        
        // declare eta correction factor
        Double_t CorrectionFactor = 1.;
        
        // electrons
        if (TMath::Abs(tracks.pId()) == pdgCodes[0]) {
	  if (electronMomentumCut(tracks.trdMom(), tracks.nSigmaTPC()[0])) continue;
            
	  if (fEtaCorrection) CorrectionFactor *= EtaCorrFactorTPCtgl(betaGamma(momentumMean(tracks.trdMom()), mElectron), tracks.trdTPCtgl(), etaMap);
	  if (NclsCorrection) CorrectionFactor *= NclsCorrFactor(betaGamma(momentumMean(tracks.trdMom()), mElectron), tracks.trdNCls(), NclsMap[tracks.trdNCh()-4]);
	  if (CentCorrection) CorrectionFactor *= CentCorrFactor(betaGamma(momentumMean(tracks.trdMom()), mElectron), tracks.centrality(), CentMap);
	  if (CorrectionFactor<1E-5) continue;

	  //   if (resSigVsBG->GetBinContent(resSigVsBG->FindBin(betaGamma(momentumMean(tracks.trdMom()), mElectron)))>ResolutionMax) continue;
	  if (FixedResStat) { // to assure that same statistics in all bins 
	    if (resSigVsTRDncls->ProjectionY("tmp", tracks.trdNCls()-39, tracks.trdNCls()-39, "oe")->Integral() <  ResolutionMax)  
	      resSigVsTRDncls->Fill(tracks.trdNCls(), CorrectionFactor*tracks.trdSig()/expMPV(meanFitHist, betaGamma(momentumMean(tracks.trdMom()), mElectron)));
	  }
	  else { // general case
	    resSigVsTRDncls->Fill(tracks.trdNCls(), CorrectionFactor*tracks.trdSig()/expMPV(meanFitHist, betaGamma(momentumMean(tracks.trdMom()), mElectron)));
	  }

	  resSigVsBG->Fill(betaGamma(momentumMean(tracks.trdMom()), mElectron), CorrectionFactor*tracks.trdSig()/expMPV(meanFitHist, betaGamma(momentumMean(tracks.trdMom()), mElectron)));
	  resSigVsTPCtgl->Fill(tracks.trdTPCtgl(), CorrectionFactor*tracks.trdSig()/expMPV(meanFitHist, betaGamma(momentumMean(tracks.trdMom()), mElectron)));
	  SigParamVsTRDncls->Fill(tracks.trdNCls(), expMPV(meanFitHist, betaGamma(momentumMean(tracks.trdMom()), mElectron)));
	  SigVsTRDncls->Fill(tracks.trdNCls(), CorrectionFactor*tracks.trdSig());
	  resSigVsTRDnch->Fill(tracks.trdNCh(), CorrectionFactor*tracks.trdSig()/expMPV(meanFitHist, betaGamma(momentumMean(tracks.trdMom()), mElectron)));
	  resSigVsTRDncls_ele->Fill(tracks.trdNCls(), CorrectionFactor*tracks.trdSig()/expMPV(meanFitHist, betaGamma(momentumMean(tracks.trdMom()), mElectron)));
	  resSigVsTRDtpctgl_ele->Fill(tracks.trdTPCtgl(), CorrectionFactor*tracks.trdSig()/expMPV(meanFitHist, betaGamma(momentumMean(tracks.trdMom()), mElectron)));
	  //  resSigVstrdGlobalPhi->Fill(tracks.trdGlobalPhi(),  CorrectionFactor*tracks.trdSig()/expMPV(meanFitHist, betaGamma(momentumMean(tracks.trdMom()), mElectron)));
	  // resSigVstrdTheta->Fill(tracks.trdTheta(),  CorrectionFactor*tracks.trdSig()/expMPV(meanFitHist, betaGamma(momentumMean(tracks.trdMom()), mElectron)));
	  // resSigVstrdY->Fill(tracks.trdY()[0],  CorrectionFactor*tracks.trdSig()/expMPV(meanFitHist, betaGamma(momentumMean(tracks.trdMom()), mElectron)));
	  // resSigVstrdPhiLocal->Fill(tracks.trdPhiLocal()[0],  CorrectionFactor*tracks.trdSig()/expMPV(meanFitHist, betaGamma(momentumMean(tracks.trdMom()), mElectron)));
	  // resSigVstrdEtaLocal->Fill(tracks.trdEtaLocal()[0],  CorrectionFactor*tracks.trdSig()/expMPV(meanFitHist, betaGamma(momentumMean(tracks.trdMom()), mElectron)));
	  nclsVsBG->Fill(betaGamma(momentumMean(tracks.trdMom()), mElectron), tracks.trdNCls());
	  TPCtglVsBG->Fill(betaGamma(momentumMean(tracks.trdMom()), mElectron), tracks.trdTPCtgl());
	  TPCtglVsncls->Fill(tracks.trdNCls(), tracks.trdTPCtgl());
	  TPCtglVsncls_ele->Fill(tracks.trdNCls(), tracks.trdTPCtgl());

 	  TRDCentralityWCuts->Fill(tracks.centrality(), CorrectionFactor*tracks.trdSig()/expMPV(meanFitHist, betaGamma(momentumMean(tracks.trdMom()), mElectron)));
	  
        }
        
        // pions
        if (TMath::Abs(tracks.pId()) == pdgCodes[1]) {
	  if (pionMomentumCut(tracks.trdMom())) continue;

	  if (fEtaCorrection) CorrectionFactor *= EtaCorrFactorTPCtgl(betaGamma(momentumMean(tracks.trdMom()), mPion), tracks.trdTPCtgl(), etaMap);
	  if (NclsCorrection) CorrectionFactor *= NclsCorrFactor(betaGamma(momentumMean(tracks.trdMom()), mPion), tracks.trdNCls(), NclsMap[tracks.trdNCh()-4]);
	  if (CentCorrection) CorrectionFactor *= CentCorrFactor(betaGamma(momentumMean(tracks.trdMom()), mPion), tracks.centrality(), CentMap);       
	  if (CorrectionFactor==0) continue;     

	  //if (resSigVsBG->GetBinContent(resSigVsBG->FindBin(betaGamma(momentumMean(tracks.trdMom()), mPion)))>ResolutionMax) continue;
	  if (FixedResStat) {
	    if (resSigVsTRDncls->ProjectionY("tmp", tracks.trdNCls()-39, tracks.trdNCls()-39, "oe")->Integral() <  ResolutionMax)  
	      resSigVsTRDncls->Fill(tracks.trdNCls(), CorrectionFactor*tracks.trdSig()/expMPV(meanFitHist, betaGamma(momentumMean(tracks.trdMom()), mPion)));
	  }
	  else {
	    resSigVsTRDncls->Fill(tracks.trdNCls(), CorrectionFactor*tracks.trdSig()/expMPV(meanFitHist, betaGamma(momentumMean(tracks.trdMom()), mPion)));
	  }
	  resSigVsBG->Fill(betaGamma(momentumMean(tracks.trdMom()), mPion),CorrectionFactor*tracks.trdSig()/expMPV(meanFitHist, betaGamma(momentumMean(tracks.trdMom()),   mPion)));
	  resSigVsTPCtgl->Fill(tracks.trdTPCtgl(),CorrectionFactor*tracks.trdSig()/expMPV(meanFitHist, betaGamma(momentumMean(tracks.trdMom()),   mPion)));
	  SigParamVsTRDncls->Fill(tracks.trdNCls(), expMPV(meanFitHist, betaGamma(momentumMean(tracks.trdMom()), mPion)));
	  SigVsTRDncls->Fill(tracks.trdNCls(), CorrectionFactor*tracks.trdSig());

	  resSigVsTRDnch->Fill(tracks.trdNCh(), CorrectionFactor*tracks.trdSig()/expMPV(meanFitHist, betaGamma(momentumMean(tracks.trdMom()), mPion)));
	  resSigVsTRDncls_pion->Fill(tracks.trdNCls(), CorrectionFactor*tracks.trdSig()/expMPV(meanFitHist, betaGamma(momentumMean(tracks.trdMom()), mPion)));
	  resSigVsTRDtpctgl_pion->Fill(tracks.trdTPCtgl(), CorrectionFactor*tracks.trdSig()/expMPV(meanFitHist, betaGamma(momentumMean(tracks.trdMom()), mPion)));
	  //resSigVstrdGlobalPhi->Fill(tracks.trdGlobalPhi(),  CorrectionFactor*tracks.trdSig()/expMPV(meanFitHist, betaGamma(momentumMean(tracks.trdMom()), mPion)));
	  //resSigVstrdTheta->Fill(tracks.trdTheta(),  CorrectionFactor*tracks.trdSig()/expMPV(meanFitHist, betaGamma(momentumMean(tracks.trdMom()), mPion)));
	  //resSigVstrdY->Fill(tracks.trdY()[0],  CorrectionFactor*tracks.trdSig()/expMPV(meanFitHist, betaGamma(momentumMean(tracks.trdMom()), mPion)));
	  //resSigVstrdPhiLocal->Fill(tracks.trdPhiLocal()[0],  CorrectionFactor*tracks.trdSig()/expMPV(meanFitHist, betaGamma(momentumMean(tracks.trdMom()), mPion)));
	  //resSigVstrdEtaLocal->Fill(tracks.trdEtaLocal()[0],  CorrectionFactor*tracks.trdSig()/expMPV(meanFitHist, betaGamma(momentumMean(tracks.trdMom()), mPion)));
	  nclsVsBG->Fill(betaGamma(momentumMean(tracks.trdMom()), mPion), tracks.trdNCls());
	  TPCtglVsBG->Fill(betaGamma(momentumMean(tracks.trdMom()), mPion), tracks.trdTPCtgl());
	  TPCtglVsncls->Fill(tracks.trdNCls(), tracks.trdTPCtgl());
	  TPCtglVsncls_pion->Fill(tracks.trdNCls(), tracks.trdTPCtgl());

	  TRDCentralityWCuts->Fill(tracks.centrality(), CorrectionFactor*tracks.trdSig()/expMPV(meanFitHist, betaGamma(momentumMean(tracks.trdMom()), mPion)));
        }
        
        // protons
        if (pidOpt == 3 && TMath::Abs(tracks.pId()) == pdgCodes[2]) {
	  if (protonMomentumCut(tracks.trdMom())) continue;	  

	  if (fEtaCorrection) CorrectionFactor *= EtaCorrFactorTPCtgl(betaGamma(momentumMean(tracks.trdMom()), mProton), tracks.trdTPCtgl(), etaMap);
	  if (NclsCorrection) CorrectionFactor *= NclsCorrFactor(betaGamma(momentumMean(tracks.trdMom()), mProton), tracks.trdNCls(), NclsMap[tracks.trdNCh()-4]);
	  if (CentCorrection) CorrectionFactor *= CentCorrFactor(betaGamma(momentumMean(tracks.trdMom()), mProton), tracks.centrality(), CentMap);
	  if (CorrectionFactor==0) continue;

	  // if (resSigVsBG->GetBinContent(resSigVsBG->FindBin(betaGamma(momentumMean(tracks.trdMom()), mProton)))>ResolutionMax) continue;
	  if (FixedResStat) {
	    if (resSigVsTRDncls->ProjectionY("tmp", tracks.trdNCls()-39, tracks.trdNCls()-39, "oe")->Integral() <  ResolutionMax)  
	      resSigVsTRDncls->Fill(tracks.trdNCls(), CorrectionFactor*tracks.trdSig()/expMPV(meanFitHist, betaGamma(momentumMean(tracks.trdMom()), mProton)));
	  }
	  else {
	    resSigVsTRDncls->Fill(tracks.trdNCls(), CorrectionFactor*tracks.trdSig()/expMPV(meanFitHist, betaGamma(momentumMean(tracks.trdMom()), mProton)));
	  }
	 
	  resSigVsBG->Fill(betaGamma(momentumMean(tracks.trdMom()), mProton), CorrectionFactor*tracks.trdSig()/expMPV(meanFitHist, betaGamma(momentumMean(tracks.trdMom()),mProton)));
	  resSigVsTPCtgl->Fill(tracks.trdTPCtgl(), CorrectionFactor*tracks.trdSig()/expMPV(meanFitHist, betaGamma(momentumMean(tracks.trdMom()),mProton)));
	  SigParamVsTRDncls->Fill(tracks.trdNCls(), expMPV(meanFitHist, betaGamma(momentumMean(tracks.trdMom()), mProton)));
	  SigVsTRDncls->Fill(tracks.trdNCls(), CorrectionFactor*tracks.trdSig());

	  resSigVsTRDnch->Fill(tracks.trdNCh(), CorrectionFactor*tracks.trdSig()/expMPV(meanFitHist, betaGamma(momentumMean(tracks.trdMom()), mProton)));
	  resSigVsTRDncls_proton->Fill(tracks.trdNCls(), CorrectionFactor*tracks.trdSig()/expMPV(meanFitHist, betaGamma(momentumMean(tracks.trdMom()), mProton)));
	  resSigVsTRDtpctgl_proton->Fill(tracks.trdTPCtgl(), CorrectionFactor*tracks.trdSig()/expMPV(meanFitHist, betaGamma(momentumMean(tracks.trdMom()), mProton)));
	  // resSigVstrdGlobalPhi->Fill(tracks.trdGlobalPhi(),  CorrectionFactor*tracks.trdSig()/expMPV(meanFitHist, betaGamma(momentumMean(tracks.trdMom()), mProton)));
	  // resSigVstrdTheta->Fill(tracks.trdTheta(),  CorrectionFactor*tracks.trdSig()/expMPV(meanFitHist, betaGamma(momentumMean(tracks.trdMom()), mProton)));
	  // resSigVstrdY->Fill(tracks.trdY()[0],  CorrectionFactor*tracks.trdSig()/expMPV(meanFitHist, betaGamma(momentumMean(tracks.trdMom()), mProton)));
	  // resSigVstrdPhiLocal->Fill(tracks.trdPhiLocal()[0],  CorrectionFactor*tracks.trdSig()/expMPV(meanFitHist, betaGamma(momentumMean(tracks.trdMom()), mProton)));
	  // resSigVstrdEtaLocal->Fill(tracks.trdEtaLocal()[0],  CorrectionFactor*tracks.trdSig()/expMPV(meanFitHist, betaGamma(momentumMean(tracks.trdMom()), mProton)));
	  nclsVsBG->Fill(betaGamma(momentumMean(tracks.trdMom()), mProton), tracks.trdNCls());
	  TPCtglVsBG->Fill(betaGamma(momentumMean(tracks.trdMom()), mProton), tracks.trdTPCtgl());
	  TPCtglVsncls->Fill(tracks.trdNCls(), tracks.trdTPCtgl());  
	  TPCtglVsncls_proton->Fill(tracks.trdNCls(), tracks.trdTPCtgl());

	  TRDCentralityWCuts->Fill(tracks.centrality(), CorrectionFactor*tracks.trdSig()/expMPV(meanFitHist, betaGamma(momentumMean(tracks.trdMom()), mProton)));
        }
    }

    // add histograms to list
    resList->Add(resSigVsBG);
    resList->Add(resSigVsTPCtgl);
    resList->Add(resSigVsTRDncls);
    resList->Add(SigParamVsTRDncls);
    resList->Add(SigVsTRDncls);
    resList->Add(resSigVsTRDncls_ele);
    resList->Add(resSigVsTRDncls_pion);
    resList->Add(resSigVsTRDncls_proton);
    resList->Add(resSigVsTRDnch);
    resListBG->Add(nclsVsBG);
    resListBG->Add(TPCtglVsBG);
    resListTPCtgl->Add(resSigVsTRDtpctgl_ele);
    resListTPCtgl->Add(resSigVsTRDtpctgl_pion);
    resListTPCtgl->Add(resSigVsTRDtpctgl_proton);
    resListTPCtgl->Add(TPCtglVsncls);
    resListTPCtgl->Add(TPCtglVsncls_ele);
    resListTPCtgl->Add(TPCtglVsncls_pion);
    resListTPCtgl->Add(TPCtglVsncls_proton);
    //resList->Add(resSigVstrdGlobalPhi);
    //resList->Add(resSigVstrdTheta);
    //resList->Add(resSigVstrdY);
    //resList->Add(resSigVstrdPhiLocal);
    //resList->Add(resSigVstrdEtaLocal);

    TRDCentralityWCuts->FitSlicesY(0,0,-1,10);
    TH1D* TRDCentralityWCuts_1 = (TH1D*)gDirectory->Get("TRDCentralityWCuts_1");
    resList->Add(TRDCentralityWCuts);
    resList->Add(TRDCentralityWCuts_1);
    resList->Add(TRDCentralityWOCuts);

    cout <<  "produce resulting histograms from signal vs bg histogram" << endl;
    FitSlicesY(resSigVsBG, resNorVsBG, resMPVVsBG, resWidthVsBG, resResVsBG, resChiVsBG, minNumberOfEntriesRes, maxFitSliceFitRMSRes, resSlicesListBG);   
    resListBG->Add(resNorVsBG);
    resListBG->Add(resMPVVsBG);
    resListBG->Add(resWidthVsBG);
    resListBG->Add(resResVsBG);
    resListBG->Add(resChiVsBG);

    cout << "produce resulting histograms from signal vs TPCtgl" << endl;
    FitSlicesY(resSigVsTPCtgl, resNorVsTPCtgl, resMPVVsTPCtgl, resWidthVsTPCtgl, resResVsTPCtgl, resChiVsTPCtgl, minNumberOfEntriesRes, maxFitSliceFitRMSRes, resSlicesListTPCtgl);
    resListTPCtgl->Add(resNorVsTPCtgl);
    resListTPCtgl->Add(resMPVVsTPCtgl);
    resListTPCtgl->Add(resWidthVsTPCtgl);
    resListTPCtgl->Add(resResVsTPCtgl);
    resListTPCtgl->Add(resChiVsTPCtgl);

    cout <<  "produce resulting histograms from signal vs TPCtgl_ele" << endl;
    FitSlicesY(resSigVsTRDtpctgl_ele, resNorVsTRDtpctgl_ele, resMPVVsTRDtpctgl_ele, resWidthVsTRDtpctgl_ele, resResVsTRDtpctgl_ele, resChiVsTRDtpctgl_ele, minNumberOfEntriesRes/10., maxFitSliceFitRMSRes, resSlicesListTPCtgl);   
    resListTPCtgl->Add(resNorVsTRDtpctgl_ele);
    resListTPCtgl->Add(resMPVVsTRDtpctgl_ele);
    resListTPCtgl->Add(resWidthVsTRDtpctgl_ele);
    resListTPCtgl->Add(resResVsTRDtpctgl_ele);
    resListTPCtgl->Add(resChiVsTRDtpctgl_ele);

    cout << "produce resulting histograms from signal vs TPCtgl pion" << endl;
    FitSlicesY(resSigVsTRDtpctgl_pion, resNorVsTRDtpctgl_pion, resMPVVsTRDtpctgl_pion, resWidthVsTRDtpctgl_pion, resResVsTRDtpctgl_pion, resChiVsTRDtpctgl_pion, minNumberOfEntriesRes/10., maxFitSliceFitRMSRes, resSlicesListTPCtgl);
    resListTPCtgl->Add(resNorVsTRDtpctgl_pion);
    resListTPCtgl->Add(resMPVVsTRDtpctgl_pion);
    resListTPCtgl->Add(resWidthVsTRDtpctgl_pion);
    resListTPCtgl->Add(resResVsTRDtpctgl_pion);
    resListTPCtgl->Add(resChiVsTRDtpctgl_pion);

    cout << " produce resulting histograms from signal vs TPCtgl proton" << endl;
    FitSlicesY(resSigVsTRDtpctgl_proton, resNorVsTRDtpctgl_proton, resMPVVsTRDtpctgl_proton, resWidthVsTRDtpctgl_proton, resResVsTRDtpctgl_proton, resChiVsTRDtpctgl_proton, minNumberOfEntriesRes/10., maxFitSliceFitRMSRes, resSlicesListTPCtgl);
    resListTPCtgl->Add(resNorVsTRDtpctgl_proton);
    resListTPCtgl->Add(resMPVVsTRDtpctgl_proton);
    resListTPCtgl->Add(resWidthVsTRDtpctgl_proton);
    resListTPCtgl->Add(resResVsTRDtpctgl_proton);
    resListTPCtgl->Add(resChiVsTRDtpctgl_proton);

    /*
    // produce resulting histograms from signal vs trdGlobalPhi 
    FitSlicesY(resSigVstrdGlobalPhi, resNorVstrdGlobalPhi, resMPVVstrdGlobalPhi, resWidthVstrdGlobalPhi, resResVstrdGlobalPhi, resChiVstrdGlobalPhi, minNumberOfEntriesRes, maxFitSliceFitRMSRes, resSlicesList);
    // add resulting histograms to list
    resList->Add(resNorVstrdGlobalPhi);
    resList->Add(resMPVVstrdGlobalPhi);
    resList->Add(resWidthVstrdGlobalPhi);
    resList->Add(resResVstrdGlobalPhi);
    resList->Add(resChiVstrdGlobalPhi);
   
    // produce resulting histograms from signal vs trdTheta 
    FitSlicesY(resSigVstrdTheta, resNorVstrdTheta, resMPVVstrdTheta, resWidthVstrdTheta, resResVstrdTheta, resChiVstrdTheta, minNumberOfEntriesRes, maxFitSliceFitRMSRes, resSlicesList);  
    // add resulting histograms to list
    resList->Add(resNorVstrdTheta);
    resList->Add(resMPVVstrdTheta);
    resList->Add(resWidthVstrdTheta);
    resList->Add(resResVstrdTheta);
    resList->Add(resChiVstrdTheta);

    // produce resulting histograms from signal vs trdY 
    FitSlicesY(resSigVstrdY, resNorVstrdY, resMPVVstrdY, resWidthVstrdY, resResVstrdY, resChiVstrdY, minNumberOfEntriesRes, maxFitSliceFitRMSRes, resSlicesList);
    // add resulting histograms to list
    resList->Add(resNorVstrdY);
    resList->Add(resMPVVstrdY);
    resList->Add(resWidthVstrdY);
    resList->Add(resResVstrdY);
    resList->Add(resChiVstrdY);

    // produce resulting histograms from signal vs PhiLocal 
    FitSlicesY(resSigVstrdPhiLocal, resNorVstrdPhiLocal, resMPVVstrdPhiLocal, resWidthVstrdPhiLocal, resResVstrdPhiLocal, resChiVstrdPhiLocal, minNumberOfEntriesRes, maxFitSliceFitRMSRes, resSlicesList);
    // add resulting histograms to list
    resList->Add(resNorVstrdPhiLocal);
    resList->Add(resMPVVstrdPhiLocal);
    resList->Add(resWidthVstrdPhiLocal);
    resList->Add(resResVstrdPhiLocal);
    resList->Add(resChiVstrdPhiLocal);

    // produce resulting histograms from signal vs etaLocal 
    FitSlicesY(resSigVstrdEtaLocal, resNorVstrdEtaLocal, resMPVVstrdEtaLocal, resWidthVstrdEtaLocal, resResVstrdEtaLocal, resChiVstrdEtaLocal, minNumberOfEntriesRes, maxFitSliceFitRMSRes, resSlicesList);
    
    // add resulting histograms to list
    resList->Add(resNorVstrdEtaLocal);
    resList->Add(resMPVVstrdEtaLocal);
    resList->Add(resWidthVstrdEtaLocal);
    resList->Add(resResVstrdEtaLocal);
    resList->Add(resChiVstrdEtaLocal);
    */

    cout << "Signal vs Centrality" << endl;
    FitSlicesY(TRDCentralityWCuts, resNorVsTRDcent, resMPVVsTRDcent, resWidthVsTRDcent, resResVsTRDcent, resChiVsTRDcent,  minNumberOfEntriesRes/10., maxFitSliceFitRMSRes , resSlicesList);   
    resList->Add(resNorVsTRDcent);
    resList->Add(resMPVVsTRDcent);
    resList->Add(resWidthVsTRDcent);
    resList->Add(resResVsTRDcent);
    resList->Add(resChiVsTRDcent);




    // ========== Main Resolution Part ==================
    
    cout << "Main Resolution NCls" << endl;
    // produce resulting histograms from signal vs trd ncls & signal vs bg histogram
    FitSlicesY(resSigVsTRDncls, resNor, resMPV, resWidth, resRes, resChi, minNumberOfEntriesRes, maxFitSliceFitRMSRes, resSlicesList);
    resList->Add(resNor);
    resList->Add(resMPV);
    resList->Add(resWidth);
    resList->Add(resRes);
    resList->Add(resChi);

    cout << "NCls proton" << endl;
   // produce resulting histograms from signal vs TRDncls
    FitSlicesY(resSigVsTRDncls_proton, resNorVsTRDncls_proton, resMPVVsTRDncls_proton, resWidthVsTRDncls_proton, resResVsTRDncls_proton, resChiVsTRDncls_proton, minNumberOfEntriesRes/10., maxFitSliceFitRMSRes, resSlicesList);   
    resList->Add(resNorVsTRDncls_proton);
    resList->Add(resMPVVsTRDncls_proton);
    resList->Add(resWidthVsTRDncls_proton);
    resList->Add(resResVsTRDncls_proton);
    resList->Add(resChiVsTRDncls_proton);

    cout << "NCls pion" << endl;
   // produce resulting histograms from signal vs TRDncls
    FitSlicesY(resSigVsTRDncls_pion, resNorVsTRDncls_pion, resMPVVsTRDncls_pion, resWidthVsTRDncls_pion, resResVsTRDncls_pion, resChiVsTRDncls_pion, minNumberOfEntriesRes/5., maxFitSliceFitRMSRes, resSlicesList);
    resList->Add(resNorVsTRDncls_pion);
    resList->Add(resMPVVsTRDncls_pion);
    resList->Add(resWidthVsTRDncls_pion);
    resList->Add(resResVsTRDncls_pion);
    resList->Add(resChiVsTRDncls_pion);

    cout << "NCls ele" << endl;
   // produce resulting histograms from signal vs TRDncls
    FitSlicesY(resSigVsTRDncls_ele, resNorVsTRDncls_ele, resMPVVsTRDncls_ele, resWidthVsTRDncls_ele, resResVsTRDncls_ele, resChiVsTRDncls_ele, minNumberOfEntriesRes/10., maxFitSliceFitRMSRes, resSlicesList);
    resList->Add(resNorVsTRDncls_ele);
    resList->Add(resMPVVsTRDncls_ele);
    resList->Add(resWidthVsTRDncls_ele);
    resList->Add(resResVsTRDncls_ele);
    resList->Add(resChiVsTRDncls_ele);

    /*
    // produce resulting histograms from signal vs trd nch  & signal vs bg histogram (commented out because of  sometimes problems due to missing statistics.)
    FitSlicesY(resSigVsTRDnch, resNorVsnch, resMPVVsnch, resWidthVsnch, resResVsnch, resChiVsnch, minNumberOfEntriesRes, maxFitSliceFitRMSRes, resSlicesList);
    // add resulting histograms to list
    resList->Add(resNorVsnch);
    resList->Add(resMPVVsnch);
    resList->Add(resWidthVsnch);
    resList->Add(resResVsnch);
    resList->Add(resChiVsnch);
    */
    
    /* Just for fun to evaluate the statistics dependence of the resolution
    TH1D * ResStat =  new TH1D("ResStat", "ResStat", 100, 0, 100000);
    //  ResStat->SetBinContent(2, 1.381);
    ResStat->SetBinContent(10 , 1.278 );
    ResStat->SetBinContent(20, 1.247);
    ResStat->SetBinContent(40, 1.235);
    ResStat->SetBinContent(60, 1.23);
    //    ResStat->SetBinContent(
    
    // resolution fit start parameters
    Double_t resParsE[] = {ResStat->GetMinimum(), 1};
    Double_t resErrsE[10], resChisE[10];
   
    // produce power law fit of resFinal histogram
    const Int_t kresfailE = ChisquareFit(ResStat, ResFunc, 2, resParsE, resErrsE, resChisE, 0, kTRUE);   // resWidth instead of resRes
  
    // get fit histogram and TF1 object
    TH1D * resFitHistE   = GetHfit("resFitHistE", ResFunc, resParsE, 0, 100000, kFALSE);
    TF1 * resFitFuncE    = new TF1("resFitFuncE", ResFunc, 0, 100000, 2);
    
    // adjust TF1 object
    resFitFuncE->SetParameters(resParsE);
    resFitFuncE->SetChisquare(resChisE[0]);
    resFitFuncE->SetNDF(resChisE[1]);
    
    // set resFitHist title
    if (kresfailE) {
        resFitHistE->SetTitle("fit fail\n");
    }
    else {
        resFitHistE->SetTitle(Form("%.3e , %.3e , (%.2f/%.0f) %s %s", TMath::Abs(resParsE[0]), TMath::Abs(resParsE[1]), resChisE[0], resChisE[1], inputFile.Data(), resFitHistE->GetTitle()));
    }co
    
    // add fit histogram and function to list
    resList->Add(ResStat);
    resList->Add(resFitHistE);
    resList->Add(resFitFuncE);
  
    */

    cout << "Fit Sqrt(N)" << endl;
    // resolution fit start parameters
    Double_t resPars[] = {resRes->GetMinimum(), 1};
    Double_t resErrs[10], resChis[10];
   
    // produce power law fit of resFinal histogram
    const Int_t kresfail = ChisquareFit(resWidth, ResFunc, 2, resPars, resErrs, resChis, 0, kTRUE);   
    // modified: resWidth instead of resRes -- if everything works fine, difference is negilgible
  
    // get fit histogram and TF1 object
    TH1D * resFitHist   = GetHfit("resFitHist", ResFunc, resPars, 40, 160, kFALSE);
    TF1  * resFitFunc    = new TF1("resFitFunc", ResFunc, 40, 160, 2);
    
    // adjust TF1 object
    resFitFunc->SetParameters(resPars);
    resFitFunc->SetChisquare(resChis[0]);
    resFitFunc->SetNDF(resChis[1]);
    
    // set resFitHist title
    if (kresfail) {
      resFitHist->SetTitle("fit fail\n");
    }
    else {
      resFitHist->SetTitle(Form("%.3e , %.3e , (%.2f/%.0f) %s %s", TMath::Abs(resPars[0]), TMath::Abs(resPars[1]), resChis[0], resChis[1], inputFile.Data(), resFitHist->GetTitle()));
    }
    // add fit histogram and function to list
    resList->Add(resFitHist);
    resList->Add(resFitFunc);
   
    /*
    // CHECK RESOLUTION FOR CHAMBER
    Double_t resParsNch[] = {resResVsnch->GetMinimum(), 1};
    Double_t resErrsNch[10], resChisNch[10];
    
    // produce power law fit of resFinal histogram
    const Int_t kresfailNch = ChisquareFit(resWidthVsnch, ResFunc, 2, resParsNch, resErrsNch, resChisNch, 0, kTRUE);   // resWidth instead of resRes
    
    // get fit histogram and TF1 object
    TH1D * resFitHistNch   = GetHfit("resFitHist", ResFunc, resParsNch, 4, 7, kFALSE);
    TF1 * resFitFuncNch    = new TF1("resFitFunc", ResFunc, 4, 7, 2);
    
    // adjust TF1 object
    resFitFuncNch->SetParameters(resParsNch);
    resFitFuncNch->SetChisquare(resChisNch[0]);
    resFitFuncNch->SetNDF(resChisNch[1]);
    
    // set resFitHist title
    if (kresfailNch) {
        resFitHistNch->SetTitle("fit fail\n");
    }
    else {
        resFitHistNch->SetTitle(Form("%.3e , %.3e , (%.2f/%.0f) %s %s", TMath::Abs(resParsNch[0]), TMath::Abs(resParsNch[1]), resChisNch[0], resChisNch[1], inputFile.Data(), resFitHistNch->GetTitle()));
    }
    
    // add fit histogram and function to list
    resList->Add(resFitHistNch);
    resList->Add(resFitFuncNch);
    */

  
    // control GaussFrac dependence
    /*    GaussFrac = 0.0; 
    TH1D * resNorC,  * resMPVC,       * resWidthC,     * resResC,       * resChiC;
    TList * resSlicesListC      = new TList;

    FitSlicesY(resSigVsTRDncls, resNorC, resMPVC, resWidthC, resResC, resChiC, minNumberOfEntriesRes, maxFitSliceFitRMSRes, resSlicesListC);
    resList->Add(resNorC);
    resList->Add(resMPVC);
    resList->Add(resWidthC);
    resList->Add(resResC);
    resList->Add(resChiC);
    
    Double_t resParsC[] = {resResC->GetMinimum(), 1};
    Double_t resErrsC[10], resChisC[10];
    
    // produce power law fit of resFinal histogram
    const Int_t kresfailC = ChisquareFit(resWidthC, ResFunc, 2, resParsC, resErrsC, resChisC, 0, kTRUE);  
  
    // get fit histogram and TF1 object
    TH1D * resFitHistC   = GetHfit("resFitHistC", ResFunc, resParsC, 60, 140, kFALSE);
    TF1 * resFitFuncC    = new TF1("resFitFuncC", ResFunc, 60, 140, 2);
    
    // adjust TF1 object
    resFitFuncC->SetParameters(resParsC);
    resFitFuncC->SetChisquare(resChisC[0]);
    resFitFuncC->SetNDF(resChisC[1]);
    
    // set resFitHist title
    if (kresfailC) {
      resFitHistC->SetTitle("fit fail\n");
    }
    else {
      resFitHistC->SetTitle(Form("%.3e , %.3e , (%.2f/%.0f) %s %s", TMath::Abs(resParsC[0]), TMath::Abs(resParsC[1]), resChisC[0], resChisC[1], inputFile.Data(), resFitHistC->GetTitle()));
    }
    
    // add fit histogram and function to list
    resList->Add(resFitHistC);
    resList->Add(resFitFuncC);

   
    TH1D* RatioFit = new TH1D("RatioFit Diff (0.00-0.50)/0.00", "RatioFit Diff (0.00-0.50)/0.00", 1000, 60, 140);
    TH1D* RatioWidth = new TH1D("RatioWidth Diff (0.00-0.50)/0.00", "RatioWidth Diff (0.00-0.50)/0.00", 120, 40, 160);
    
    RatioFit->Add(resFitHistC, resFitHist, 1, -1);
    RatioFit->Divide(resFitHist);
    RatioWidth->Add(resWidthC, resWidth, 1, -1);
    RatioWidth->Divide(resWidth);

    TCanvas * ControlGaussFrac = new TCanvas("ControlGaussFrac", "ControlGaussFrac", 800, 800);
    ControlGaussFrac->Divide(2,2);
    ControlGaussFrac->cd(1);
    gStyle->SetOptStat(0);
    resWidthC->SetLineColor(1);
    resWidthC->Draw("");
    resWidthC->SetMaximum(0.2);
    resWidthC->SetMinimum(0.1);
    resWidth->SetLineColor(2);
    resWidth->Draw("same");

    ControlGaussFrac->cd(2);
    resFitHistC->SetLineColor(1);
    resFitHistC->Draw("");
    resFitHistC->SetMaximum(0.2);
    resFitHistC->SetMinimum(0.1);
    resFitHist->SetLineColor(2);
    resFitHist->Draw("same");

    ControlGaussFrac->cd(3);
    RatioWidth->Draw("");
    RatioWidth->SetMaximum(0.1);
    RatioWidth->SetMinimum(-0.1);

    ControlGaussFrac->cd(4);
    RatioFit->Draw("");
    RatioFit->SetMaximum(0.1);
    RatioFit->SetMinimum(-0.1);

    ControlGaussFrac->SaveAs(Form("output/ControlGaussFrac_Resolution_%i_%s%s%s.pdf",nch, FileSpec.Data(), EtaSpec.Data(), NclsCorr.Data()));

    resList->Add(RatioWidth);
    resList->Add(RatioFit);
    */


    
    //===========================================================================================
    //
    // save output
    //
    //===========================================================================================
    
    // ctime for time + date output to parameter .txt file
    std::time_t time = std::time(NULL);
    
    // histograms
    TFile resFile(Form("output/%s_resParam_Iter%i_pid%i_%s%s%s%s.root", inputSpec.Data(),Iteration, pidOpt, FileSpec.Data(), CentCorr.Data(), EtaSpec.Data(), NclsCorr.Data()),  "RECREATE");
    resList->Write();
    resSlicesList->Write("slices", TObject::kSingleKey);

    TDirectory *BG = resFile.mkdir("BG");
    BG->cd();    // make the "tof" directory the current directory
    resListBG->Write();
    resSlicesListBG->Write("slices", TObject::kSingleKey);
    
    TDirectory *TPCtgl = resFile.mkdir("TPCtgl");
    TPCtgl->cd();    // make the "tof" directory the current directory
    resListTPCtgl->Write();
    resSlicesListTPCtgl->Write("slices", TObject::kSingleKey);


    resFile.Close();
    resFile.Save();
    
    // parameters
    FILE* resParamFile = fopen(Form("output/%s_resParam_Iter%i_pid%i_%s%s%s%s.txt", inputSpec.Data(), Iteration, pidOpt, FileSpec.Data(), CentCorr.Data(), EtaSpec.Data(), NclsCorr.Data()), "a");
    fprintf(resParamFile, "\n%s\n", std::ctime(&time));
    fprintf(resParamFile, "res vs ncls: \n");
    fprintf(resParamFile, "chi^2 = %.1f\nndf = %.1f\nchi^2/ndf = %.1f\n\n", resChis[0], resChis[1], resChis[0]/resChis[1]);
    for (Int_t i = 0; i < 2; i++) {
        fprintf(resParamFile, "par[%i] = %.3e\t\terr[%i] = %.3e\n", i, TMath::Abs(resPars[i]), i, resErrs[i]);
    }
    fclose(resParamFile);


    delete  resList;
    delete  resSlicesList;
    delete  resSigVsTRDncls; 
    delete  SigParamVsTRDncls;
    delete  SigVsTRDncls;
    delete  resSigVsTRDnch;
    delete  resSigVsBG;
    delete  resSigVsTPCtgl;
    //delete  resSigVstrdGlobalPhi;
    //delete  resSigVstrdThTPCtgl;
    //delete  resSigVstrdY;
    //delete  resSigVstrdEtaLocal;
    //delete  resSigVstrdPhiLocal;
    delete  nclsVsBG;
    delete  TPCtglVsBG;
    delete  TPCtglVsncls;
    delete  resSigVsTRDncls_ele;
    delete  resSigVsTRDncls_pion;
    delete  resSigVsTRDncls_proton;
    delete  resSigVsTRDtpctgl_ele;
    delete  resSigVsTRDtpctgl_pion;
    delete  resSigVsTRDtpctgl_proton;
    delete  resFitFunc;
   
}
