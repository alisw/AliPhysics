//====================================================================================================================================
//
// root macro to paramtrize the dE/dx (vs beta*gamma) for e, pi, p
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
#include <TDirectory.h>
//#include "V0Tracks.h"
#include "constants.h"
#include "utils.h"


// Main class for the signal parametrization. The macro is steered by the include of utils.h


void signal() {
  // QAPlots requires tree variables (TRDNcls, 
  Bool_t QAPlots = kTRUE;

  // if no charge seperation is wanted, set tracks.pId() to abs
  Int_t pdgCodes[3] = {11, 211, 2212};
    
  // get and initialize tree
  V0Tracks tracks;
  tracks.setFile(input.Data());
  tracks.getTree();
    
  // get number of entries from tree
  const Int_t nEntries=tracks.getNumberOfEntries();
    
  // create new TList objects for mean and resolution parametrization
  TList * meanList            = new TList; // main output list
  TList * QAList              = new TList; // general QA list
  TList * ClusterQAList       = new TList; // cluster QA list
  TList * EtaQAList           = new TList; // eta QA list
  TList * meanSlicesList      = new TList; // slice fit results
  TList * dummyList           = new TList; // 
    
  // declare histograms for mean and resolution parametrization
  TH2D * meanSigVsBG        = new TH2D("meanSigVsBG", "meanSigVsBG", 100, 0.3, 1e4, 200, 0, 10);
  TH2D * meanSigVsBG_ele    = new TH2D("meanSigVsBG_ele", "meanSigVsBG_ele", 100, 0.3, 1e4, 200, 0, 10);
  TH2D * meanSigVsBG_proton = new TH2D("meanSigVsBG_proton", "meanSigVsBG_proton", 100, 0.3, 1e4, 200, 0, 10);
  TH2D * meanSigVsBG_pion   = new TH2D("meanSigVsBG_pion", "meanSigVsBG_pion", 100, 0.3, 1e4, 200, 0, 10);
  TH2D * meanSigVsBGnorm;
  TH1D * meanNor, * meanMPV, *meanMPV_ele, *meanMPV_pion, *meanMPV_proton, * meanWidth, * meanRes, * meanChi; // fit results

  // only for Sign dependence
  TH2D * meanSigVsBG_plus        = new TH2D("meanSigVsBG_plus", "meanSigVsBG_plus", 100, 0.3, 1e4, 200, 0, 10);
  TH2D * meanSigVsBG_plus_ele    = new TH2D("meanSigVsBG_plus_ele", "meanSigVsBG_plus_ele", 100, 0.3, 1e4, 200, 0, 10);
  TH2D * meanSigVsBG_plus_proton = new TH2D("meanSigVsBG_plus_proton", "meanSigVsBG_plus_proton", 100, 0.3, 1e4, 200, 0, 10);
  TH2D * meanSigVsBG_plus_pion   = new TH2D("meanSigVsBG_plus_pion", "meanSigVsBG_plus_pion", 100, 0.3, 1e4, 200, 0, 10);
  TH1D * meanMPV_plus, *meanMPV_plus_ele, *meanMPV_plus_pion, *meanMPV_plus_proton;

  TH2D * meanSigVsBG_minus        = new TH2D("meanSigVsBG_minus", "meanSigVsBG_minus", 100, 0.3, 1e4, 200, 0, 10);
  TH2D * meanSigVsBG_minus_ele    = new TH2D("meanSigVsBG_minus_ele", "meanSigVsBG_minus_ele", 100, 0.3, 1e4, 200, 0, 10);
  TH2D * meanSigVsBG_minus_proton = new TH2D("meanSigVsBG_minus_proton", "meanSigVsBG_minus_proton", 100, 0.3, 1e4, 200, 0, 10);
  TH2D * meanSigVsBG_minus_pion   = new TH2D("meanSigVsBG_minus_pion", "meanSigVsBG_minus_pion", 100, 0.3, 1e4, 200, 0, 10);
  TH1D * meanMPV_minus, *meanMPV_minus_ele, *meanMPV_minus_pion, *meanMPV_minus_proton;

  // declare histograms for cluster vs chamber number
  TH2D * clusterVsChamber = new TH2D("clusterVsChamber", "clusterVsChamber", 150, 10, 160, 6, 1, 7);
  TH2D * clusterVsChamberWCut = new TH2D("clusterVsChamberWCut", "clusterVsChamberWCut", 150, 10, 160, 6, 1, 7);
  TH2D * clusterVsChamber_ele = new TH2D("clusterVsChamber_ele", "clusterVsChamber_ele", 150, 10, 160, 6, 1, 7);
  TH2D * clusterVsChamberWCut_ele = new TH2D("clusterVsChamberWCut_ele", "clusterVsChamberWCut_ele", 150, 10, 160, 6, 1, 7);
  TH2D * clusterVsChamber_pion = new TH2D("clusterVsChamber_pion", "clusterVsChamber_pion", 150, 10, 160, 6, 1, 7);
  TH2D * clusterVsChamberWCut_pion = new TH2D("clusterVsChamberWCut_pion", "clusterVsChamberWCut_pion", 150, 10, 160, 6, 1, 7);
  TH2D * clusterVsChamber_proton = new TH2D("clusterVsChamber_proton", "clusterVsChamber_proton", 150, 10, 160, 6, 1, 7);
  TH2D * clusterVsChamberWCut_proton = new TH2D("clusterVsChamberWCut_proton", "clusterVsChamberWCut_proton", 150, 10, 160, 6, 1, 7);

  // declare histograms for cluster vs centrality
  TH2D * clusterVsCentrality = new TH2D("clusterVsCentrality", "clusterVsCentrality",  NCentBins, CentBins, 150, 10, 160);
  TH2D * clusterVsCentrality_ele = new TH2D("clusterVsCentrality_ele", "clusterVsCentrality_ele",  NCentBins, CentBins, 150, 10, 160);
  TH2D * clusterVsCentrality_pion = new TH2D("clusterVsCentrality_pion", "clusterVsCentrality_pion",  NCentBins, CentBins, 150, 10, 160);
  TH2D * clusterVsCentrality_proton = new TH2D("clusterVsCentrality_proton", "clusterVsCentrality_proton",  NCentBins, CentBins, 150, 10, 160);
    
  TH1D * ClusterVSChamber_ele[3][NCentBins];
  TH1D * ClusterVSChamber_pion[3][NCentBins];
  TH1D * ClusterVSChamber_proton[3][NCentBins];
  for (int i = 0; i<3; i++) { // ntracklets 4,5,6
    for (int j=0; j<NCentBins; j++) { // nCentbins
      ClusterVSChamber_ele[i][j] = new TH1D(Form("ClusterVsChamber_ele_%iNch_%iCentBin", i+4, j), Form("ClusterVsChamber_ele_%iNch_%iCentBin", i+4, j), 150, 10, 160);
      ClusterVSChamber_pion[i][j] = new TH1D(Form("ClusterVsChamber_pion_%iNch_%iCentBin", i+4, j), Form("ClusterVsChamber_pion_%iNch_%iCentBin", i+4, j), 150, 10, 160);
      ClusterVSChamber_proton[i][j] = new TH1D(Form("ClusterVsChamber_proton_%iNch_%iCentBin", i+4, j), Form("ClusterVsChamber_proton_%iNch_%iCentBin", i+4, j), 150, 10, 160);
    }
  }

  TH1D * ClusterVSChamberEta_ele[3][NEtaBins];
  TH1D * ClusterVSChamberEta_pion[3][NEtaBins];
  TH1D * ClusterVSChamberEta_proton[3][NEtaBins];
  for (int i = 0; i<3; i++) { // ntracklets 4,5,6
    for (int j=0; j<NEtaBins; j++) { // nEtaBins
      ClusterVSChamberEta_ele[i][j] = new TH1D(Form("ClusterVsChamberEta_ele_%iNch_%iEtaBin", i+4, j), Form("ClusterVsChamberEta_ele_%iNch_%iEtaBin", i+4, j), 150, 10, 160);
      ClusterVSChamberEta_pion[i][j] = new TH1D(Form("ClusterVsChamberEta_pion_%iNch_%iEtaBin", i+4, j), Form("ClusterVsChamberEta_pion_%iNch_%iEtaBin", i+4, j), 150, 10, 160);
      ClusterVSChamberEta_proton[i][j] = new TH1D(Form("ClusterVsChamberEta_proton_%iNch_%iEtaBin", i+4, j), Form("ClusterVsChamberEta_proton_%iNch_%iEtaBin", i+4, j), 150, 10, 160);
    }
  }


  TH2D * clusterVsBG = new TH2D("clusterVsBG", "clusterVsBG",100, 0.3, 1e4, 150, 10, 160);
  BinLogX(clusterVsBG->GetXaxis());
  TH2D * clusterVsBG_ele = new TH2D("clusterVsBG_ele", "clusterVsBG_ele",100, 0.3, 1e4, 150, 10, 160);
  BinLogX(clusterVsBG_ele->GetXaxis());
  TH2D * clusterVsBG_pion = new TH2D("clusterVsBG_pion", "clusterVsBG_pion",100, 0.3, 1e4, 150, 10, 160);
  BinLogX(clusterVsBG_pion->GetXaxis());
  TH2D * clusterVsBG_proton = new TH2D("clusterVsBG_proton", "clusterVsBG_proton",100, 0.3, 1e4, 150, 10, 160);
  BinLogX(clusterVsBG_proton->GetXaxis());
  TH1D * MPVclusterVsBG, * MPVclusterVsBG_ele, * MPVclusterVsBG_proton, * MPVclusterVsBG_pion; 
  
  /*
    TH1D * MPVclusterVsBG=new TH1D("MPVclusterVsBG", "MPVclusterVsBG",100, 0.3, 1e4);
    BinLogX(MPVclusterVsBG->GetXaxis());
    TH1D * MPVclusterVsBG_ele=new TH1D("MPVclusterVsBG_ele", "MPVclusterVsBG_ele",100, 0.3, 1e4);
    BinLogX(MPVclusterVsBG_ele->GetXaxis());
    TH1D * MPVclusterVsBG_pion=new TH1D("MPVclusterVsBG_pion", "MPVclusterVsBG_pion",100, 0.3, 1e4);
    BinLogX(MPVclusterVsBG_pion->GetXaxis());
    TH1D * MPVclusterVsBG_proton=new TH1D("MPVclusterVsBG_proton", "MPVclusterVsBG_proton",100, 0.3, 1e4);
    BinLogX(MPVclusterVsBG_proton->GetXaxis());
  */


  TH2D * GenClusVsBG = new TH2D("GenClusVsBG", "GenClusVsBG",100, 0.3, 1e4, 150, 10, 160);
  BinLogX(GenClusVsBG->GetXaxis());
  TH2D * GenClusVsBG_ele = new TH2D("GenClusVsBG_ele", "GenClusVsBG_ele",100, 0.3, 1e4, 150, 10, 160);
  BinLogX(GenClusVsBG_ele->GetXaxis());
  TH2D * GenClusVsBG_pion = new TH2D("GenClusVsBG_pion", "GenClusVsBG_pion",100, 0.3, 1e4, 150, 10, 160);
  BinLogX(GenClusVsBG_pion->GetXaxis());
  TH2D * GenClusVsBG_proton = new TH2D("GenClusVsBG_proton", "GenClusVsBG_proton",100, 0.3, 1e4, 150, 10, 160);
  BinLogX(GenClusVsBG_proton->GetXaxis());
  TH1D * MPVGenClusVsBG, * MPVGenClusVsBG_ele, * MPVGenClusVsBG_pion, * MPVGenClusVsBG_proton;

  /*=new TH1D("MPVGenClusVsBG", "MPVGenClusVsBG",100, 0.3, 1e4);
    BinLogX(MPVGenClusVsBG->GetXaxis());
    TH1D * MPVGenClusVsBG_ele=new TH1D("MPVGenClusVsBG_ele", "MPVGenClusVsBG_ele",100, 0.3, 1e4);
    BinLogX(MPVGenClusVsBG_ele->GetXaxis());
    TH1D * MPVGenClusVsBG_pion=new TH1D("MPVGenClusVsBG_pion", "MPVGenClusVsBG_pion",100, 0.3, 1e4);
    BinLogX(MPVGenClusVsBG_pion->GetXaxis());
    TH1D * MPVGenClusVsBG_proton=new TH1D("MPVGenClusVsBG_proton", "MPVGenClusVsBG_proton",100, 0.3, 1e4);
    BinLogX(MPVGenClusVsBG_proton->GetXaxis());
  */

  TH2D * EtaTPCtglVsBG = new TH2D("EtaTPCtglVsBG", "EtaTPCtglVsBG",100, 0.3, 1e4, 180, -1, 1);
  BinLogX(EtaTPCtglVsBG->GetXaxis());
  TH2D * EtaTPCtglVsBG_ele = new TH2D("EtaTPCtglVsBG_ele", "EtaTPCtglVsBG_ele",100, 0.3, 1e4, 180, -1, 1);
  BinLogX(EtaTPCtglVsBG_ele->GetXaxis());
  TH2D * EtaTPCtglVsBG_pion = new TH2D("EtaTPCtglVsBG_pion", "EtaTPCtglVsBG_pion",100, 0.3, 1e4, 180, -1, 1);
  BinLogX(EtaTPCtglVsBG_pion->GetXaxis());
  TH2D * EtaTPCtglVsBG_proton = new TH2D("EtaTPCtglVsBG_proton", "EtaTPCtglVsBG_proton",100, 0.3, 1e4, 180, -1, 1);
  BinLogX(EtaTPCtglVsBG_proton->GetXaxis());
  TH1D * MPVEtaTPCtglVsBG, * MPVEtaTPCtglVsBG_ele, * MPVEtaTPCtglVsBG_pion, * MPVEtaTPCtglVsBG_proton;
  /*
    =new TH1D("MPVEtaTPCtglVsBG", "MPVEtaTPCtglVsBG",100, 0.3, 1e4);
    BinLogX(MPVEtaTPCtglVsBG->GetXaxis());
    TH1D * MPVEtaTPCtglVsBG_ele=new TH1D("MPVEtaTPCtglVsBG_ele", "MPVEtaTPCtglVsBG_ele",100, 0.3, 1e4);
    BinLogX(MPVEtaTPCtglVsBG_ele->GetXaxis());
    TH1D * MPVEtaTPCtglVsBG_pion=new TH1D("MPVEtaTPCtglVsBG_pion", "MPVEtaTPCtglVsBG_pion",100, 0.3, 1e4);
    BinLogX(MPVEtaTPCtglVsBG_pion->GetXaxis());
    TH1D * MPVEtaTPCtglVsBG_proton=new TH1D("MPVEtaTPCtglVsBG_proton", "MPVEtaTPCtglVsBG_proton",100, 0.3, 1e4);
    BinLogX(MPVEtaTPCtglVsBG_proton->GetXaxis());
  */    

  
  TH2D * clusterVSBG[3][NCentBins];
  TH2D * clusterVSBG_ele[3][NCentBins];
  TH2D * clusterVSBG_pion[3][NCentBins];
  TH2D * clusterVSBG_proton[3][NCentBins];
  TH1D * MPVclusterVSBG[3][NCentBins];
  TH1D * MPVclusterVSBG_ele[3][NCentBins];
  TH1D * MPVclusterVSBG_pion[3][NCentBins];
  TH1D * MPVclusterVSBG_proton[3][NCentBins];
  
  
  TH2D * clusterVSBG_Nch[3];
  TH2D * clusterVSBG_Nch_ele[3];
  TH2D * clusterVSBG_Nch_pion[3];
  TH2D * clusterVSBG_Nch_proton[3];
  TH1D * MPVclusterVSBG_Nch[3];
  TH1D * MPVclusterVSBG_Nch_ele[3];
  TH1D * MPVclusterVSBG_Nch_pion[3];
  TH1D * MPVclusterVSBG_Nch_proton[3];
  

  TH2D * GenClusVSBG[3][NCentBins];
  TH2D * GenClusVSBG_ele[3][NCentBins];
  TH2D * GenClusVSBG_pion[3][NCentBins];
  TH2D * GenClusVSBG_proton[3][NCentBins];
  TH1D * MPVGenClusVSBG[3][NCentBins];
  TH1D * MPVGenClusVSBG_ele[3][NCentBins];
  TH1D * MPVGenClusVSBG_pion[3][NCentBins];
  TH1D * MPVGenClusVSBG_proton[3][NCentBins];
  
  
  TH2D * GenClusVSBG_Nch[3];
  TH2D * GenClusVSBG_Nch_ele[3];
  TH2D * GenClusVSBG_Nch_pion[3];
  TH2D * GenClusVSBG_Nch_proton[3];
  TH1D * MPVGenClusVSBG_Nch[3];
  TH1D * MPVGenClusVSBG_Nch_ele[3];
  TH1D * MPVGenClusVSBG_Nch_pion[3];
  TH1D * MPVGenClusVSBG_Nch_proton[3];
  
  TH2D * meanSigVSBG[3][NCentBins];
  TH2D * meanSigVSBG_ele[3][NCentBins];
  TH2D * meanSigVSBG_pion[3][NCentBins];
  TH2D * meanSigVSBG_proton[3][NCentBins];
  TH1D * MPVmeanSigVSBG[3][NCentBins];
  TH1D * MPVmeanSigVSBG_ele[3][NCentBins];
  TH1D * MPVmeanSigVSBG_pion[3][NCentBins];
  TH1D * MPVmeanSigVSBG_proton[3][NCentBins];
  
  TH2D * meanSigVSBGEta[3][NEtaBins];
  TH2D * meanSigVSBGEta_ele[3][NEtaBins];
  TH2D * meanSigVSBGEta_pion[3][NEtaBins];
  TH2D * meanSigVSBGEta_proton[3][NEtaBins];
  TH1D * MPVmeanSigVSBGEta[3][NEtaBins];
  TH1D * MPVmeanSigVSBGEta_ele[3][NEtaBins];
  TH1D * MPVmeanSigVSBGEta_pion[3][NEtaBins];
  TH1D * MPVmeanSigVSBGEta_proton[3][NEtaBins];
  
  TH2D * meanSigVSBG_Nch[3];
  TH2D * meanSigVSBG_Nch_ele[3];
  TH2D * meanSigVSBG_Nch_pion[3];
  TH2D * meanSigVSBG_Nch_proton[3];
  TH1D * MPVmeanSigVSBG_Nch[3];
  TH1D * MPVmeanSigVSBG_Nch_ele[3];
  TH1D * MPVmeanSigVSBG_Nch_pion[3];
  TH1D * MPVmeanSigVSBG_Nch_proton[3];

  TH2D * EtaTPCtglVSBG[3][NCentBins];
  TH2D * EtaTPCtglVSBG_ele[3][NCentBins];
  TH2D * EtaTPCtglVSBG_pion[3][NCentBins];
  TH2D * EtaTPCtglVSBG_proton[3][NCentBins];
  TH1D * MPVEtaTPCtglVSBG[3][NCentBins];
  TH1D * MPVEtaTPCtglVSBG_ele[3][NCentBins];
  TH1D * MPVEtaTPCtglVSBG_pion[3][NCentBins];
  TH1D * MPVEtaTPCtglVSBG_proton[3][NCentBins];
  
  TH2D * EtaTPCtglVSBG_Nch[3];
  TH2D * EtaTPCtglVSBG_Nch_ele[3];
  TH2D * EtaTPCtglVSBG_Nch_pion[3];
  TH2D * EtaTPCtglVSBG_Nch_proton[3];
  TH1D * MPVEtaTPCtglVSBG_Nch[3];
  TH1D * MPVEtaTPCtglVSBG_Nch_ele[3];
  TH1D * MPVEtaTPCtglVSBG_Nch_pion[3];
  TH1D * MPVEtaTPCtglVSBG_Nch_proton[3];


  for (int i = 0; i<3; i++) {   // for each nChamber 4,5,6
    clusterVSBG_Nch[i]= new TH2D(Form("ClusterVSBG_Nch_%iNch", i+4), Form("ClusterVSBG_Nch_%iNch", i+4), 100, 0.3, 1e4, 150, 10, 160); 
    BinLogX(clusterVSBG_Nch[i]->GetXaxis());

    clusterVSBG_Nch_ele[i]= new TH2D(Form("ClusterVSBG_Nch_ele_%iNch", i+4), Form("ClusterVSBG_Nch_ele_%iNch", i+4), 100, 0.3, 1e4, 150, 10, 160); 
    BinLogX(clusterVSBG_Nch_ele[i]->GetXaxis());
    clusterVSBG_Nch_pion[i]= new TH2D(Form("ClusterVSBG_Nch_pion_%iNch", i+4), Form("ClusterVSBG_Nch_pion_%iNch", i+4), 100, 0.3, 1e4, 150, 10, 160); 
    BinLogX(clusterVSBG_Nch_pion[i]->GetXaxis());
    clusterVSBG_Nch_proton[i]= new TH2D(Form("ClusterVSBG_Nch_proton_%iNch", i+4), Form("ClusterVSBG_Nch_proton_%iNch", i+4), 100, 0.3, 1e4, 150, 10, 160);
    BinLogX(clusterVSBG_Nch_proton[i]->GetXaxis());
	/*
	MPVclusterVSBG_Nch[i]= new TH1D(Form("MPVClusterVSBG_Nch_%iNch", i+4), Form("MPVClusterVSBG_Nch_%iNch", i+4), 100, 0.3, 1e4); 
	BinLogX(MPVclusterVSBG_Nch[i]->GetXaxis());
	MPVclusterVSBG_Nch_ele[i]= new TH1D(Form("MPVClusterVSBG_Nch_ele_%iNch", i+4), Form("MPVClusterVSBG_Nch_ele_%iNch", i+4), 100, 0.3,1e4); 
	BinLogX(MPVclusterVSBG_Nch_ele[i]->GetXaxis());
	MPVclusterVSBG_Nch_pion[i]= new TH1D(Form("MPVClusterVSBG_Nch_pion_%iNch", i+4), Form("MPVClusterVSBG_Nch_pion_%iNch", i+4), 100, 0.3, 1e4); 
	BinLogX(MPVclusterVSBG_Nch_pion[i]->GetXaxis());
	MPVclusterVSBG_Nch_proton[i]= new TH1D(Form("MPVClusterVSBG_Nch_proton_%iNch", i+4), Form("MPVClusterVSBG_Nch_proton_%iNch", i+4), 100, 0.3, 1e4);
	BinLogX(MPVclusterVSBG_Nch_proton[i]->GetXaxis());
	*/

    GenClusVSBG_Nch[i]= new TH2D(Form("GenClusVSBG_Nch_%iNch", i+4), Form("GenClusVSBG_Nch_%iNch", i+4), 100, 0.3, 1e4, 150, 10, 160); 
    BinLogX(GenClusVSBG_Nch[i]->GetXaxis());
    GenClusVSBG_Nch_ele[i]= new TH2D(Form("GenClusVSBG_Nch_ele_%iNch", i+4), Form("GenClusVSBG_Nch_ele_%iNch", i+4), 100, 0.3, 1e4, 150, 10, 160); 
    BinLogX(GenClusVSBG_Nch_ele[i]->GetXaxis());
    GenClusVSBG_Nch_pion[i]= new TH2D(Form("GenClusVSBG_Nch_pion_%iNch", i+4), Form("GenClusVSBG_Nch_pion_%iNch", i+4), 100, 0.3, 1e4, 150, 10, 160); 
    BinLogX(GenClusVSBG_Nch_pion[i]->GetXaxis());
    GenClusVSBG_Nch_proton[i]= new TH2D(Form("GenClusVSBG_Nch_proton_%iNch", i+4), Form("GenClusVSBG_Nch_proton_%iNch", i+4), 100, 0.3, 1e4, 150, 10, 160);
    BinLogX(GenClusVSBG_Nch_proton[i]->GetXaxis());
	
	/*
	MPVGenClusVSBG_Nch[i]= new TH1D(Form("MPVGenClusVSBG_Nch_%iNch", i+4), Form("MPVGenClusVSBG_Nch_%iNch", i+4), 100, 0.3, 1e4); 
	BinLogX(MPVGenClusVSBG_Nch[i]->GetXaxis());
	MPVGenClusVSBG_Nch_ele[i]= new TH1D(Form("MPVGenClusVSBG_Nch_ele_%iNch", i+4), Form("MPVGenClusVSBG_Nch_ele_%iNch", i+4), 100, 0.3,1e4); 
	BinLogX(MPVGenClusVSBG_Nch_ele[i]->GetXaxis());
	MPVGenClusVSBG_Nch_pion[i]= new TH1D(Form("MPVGenClusVSBG_Nch_pion_%iNch", i+4), Form("MPVGenClusVSBG_Nch_pion_%iNch", i+4), 100, 0.3, 1e4); 
	BinLogX(MPVGenClusVSBG_Nch_pion[i]->GetXaxis());
	MPVGenClusVSBG_Nch_proton[i]= new TH1D(Form("MPVGenClusVSBG_Nch_proton_%iNch", i+4), Form("MPVGenClusVSBG_Nch_proton_%iNch", i+4), 100, 0.3, 1e4);
	BinLogX(MPVGenClusVSBG_Nch_proton[i]->GetXaxis());
	*/

    meanSigVSBG_Nch[i]= new TH2D(Form("MeanSigVSBG_Nch_%iNch", i+4), Form("MeanSigVSBG_Nch_%iNch", i+4), 100, 0.3, 1e4, 200, 0, 10); 
    BinLogX(meanSigVSBG_Nch[i]->GetXaxis());
    meanSigVSBG_Nch_ele[i]= new TH2D(Form("MeanSigVSBG_Nch_ele_%iNch", i+4), Form("MeanSigVSBG_Nch_ele_%iNch", i+4), 100, 0.3, 1e4, 200, 0, 10); 
    BinLogX(meanSigVSBG_Nch_ele[i]->GetXaxis());
    meanSigVSBG_Nch_pion[i]= new TH2D(Form("MeanSigVSBG_Nch_pion_%iNch", i+4), Form("MeanSigVSBG_Nch_pion_%iNch", i+4), 100, 0.3, 1e4, 200, 0, 10); 
    BinLogX(meanSigVSBG_Nch_pion[i]->GetXaxis());
    meanSigVSBG_Nch_proton[i]= new TH2D(Form("MeanSigVSBG_Nch_proton_%iNch", i+4), Form("MeanSigVSBG_Nch_proton_%iNch", i+4), 100, 0.3, 1e4, 200, 0, 10);
    BinLogX(meanSigVSBG_Nch_proton[i]->GetXaxis());
	/*
	MPVmeanSigVSBG_Nch[i]= new TH1D(Form("MPVMeanSigVSBG_Nch_%iNch", i+4), Form("MPVMeanSigVSBG_Nch_%iNch", i+4), 100, 0.3, 1e4); 
	BinLogX(MPVmeanSigVSBG_Nch[i]->GetXaxis());
	MPVmeanSigVSBG_Nch_ele[i]= new TH1D(Form("MPVMeanSigVSBG_Nch_ele_%iNch", i+4), Form("MPVMeanSigVSBG_Nch_ele_%iNch", i+4), 100, 0.3,1e4); 
	BinLogX(MPVmeanSigVSBG_Nch_ele[i]->GetXaxis());
	MPVmeanSigVSBG_Nch_pion[i]= new TH1D(Form("MPVMeanSigVSBG_Nch_pion_%iNch", i+4), Form("MPVMeanSigVSBG_Nch_pion_%iNch", i+4), 100, 0.3, 1e4); 
	BinLogX(MPVmeanSigVSBG_Nch_pion[i]->GetXaxis());
	MPVmeanSigVSBG_Nch_proton[i]= new TH1D(Form("MPVMeanSigVSBG_Nch_proton_%iNch", i+4), Form("MPVMeanSigVSBG_Nch_proton_%iNch", i+4), 100, 0.3, 1e4);
	BinLogX(MPVmeanSigVSBG_Nch_proton[i]->GetXaxis());

	*/
    EtaTPCtglVSBG_Nch[i]= new TH2D(Form("EtaTPCtglVSBG_Nch_%iNch", i+4), Form("EtaTPCtglVSBG_Nch_%iNch", i+4), 100, 0.3, 1e4, 180, -1, 1); 
    BinLogX(EtaTPCtglVSBG_Nch[i]->GetXaxis());
    EtaTPCtglVSBG_Nch_ele[i]= new TH2D(Form("EtaTPCtglVSBG_Nch_ele_%iNch", i+4), Form("EtaTPCtglVSBG_Nch_ele_%iNch", i+4), 100, 0.3, 1e4, 180, -1, 1); 
    BinLogX(EtaTPCtglVSBG_Nch_ele[i]->GetXaxis());
    EtaTPCtglVSBG_Nch_pion[i]= new TH2D(Form("EtaTPCtglVSBG_Nch_pion_%iNch", i+4), Form("EtaTPCtglVSBG_Nch_pion_%iNch", i+4), 100, 0.3, 1e4, 180, -1, 1); 
    BinLogX(EtaTPCtglVSBG_Nch_pion[i]->GetXaxis());
    EtaTPCtglVSBG_Nch_proton[i]= new TH2D(Form("EtaTPCtglVSBG_Nch_proton_%iNch", i+4), Form("EtaTPCtglVSBG_Nch_proton_%iNch", i+4), 100, 0.3, 1e4, 180, -1, 1);
    BinLogX(EtaTPCtglVSBG_Nch_proton[i]->GetXaxis());
	/*
	MPVEtaTPCtglVSBG_Nch[i]= new TH1D(Form("MPVEtaTPCtglVSBG_Nch_%iNch", i+4), Form("MPVEtaTPCtglVSBG_Nch_%iNch", i+4), 100, 0.3, 1e4); 
	BinLogX(MPVEtaTPCtglVSBG_Nch[i]->GetXaxis());
	MPVEtaTPCtglVSBG_Nch_ele[i]= new TH1D(Form("MPVEtaTPCtglVSBG_Nch_ele_%iNch", i+4), Form("MPVEtaTPCtglVSBG_Nch_ele_%iNch", i+4), 100, 0.3,1e4); 
	BinLogX(MPVEtaTPCtglVSBG_Nch_ele[i]->GetXaxis());
	MPVEtaTPCtglVSBG_Nch_pion[i]= new TH1D(Form("MPVEtaTPCtglVSBG_Nch_pion_%iNch", i+4), Form("MPVEtaTPCtglVSBG_Nch_pion_%iNch", i+4), 100, 0.3, 1e4); 
	BinLogX(MPVEtaTPCtglVSBG_Nch_pion[i]->GetXaxis());
	MPVEtaTPCtglVSBG_Nch_proton[i]= new TH1D(Form("MPVEtaTPCtglVSBG_Nch_proton_%iNch", i+4), Form("MPVEtaTPCtglVSBG_Nch_proton_%iNch", i+4), 100, 0.3, 1e4);
	BinLogX(MPVEtaTPCtglVSBG_Nch_proton[i]->GetXaxis());
	*/


    for (int j=0; j<NCentBins; j++) {
      clusterVSBG[i][j]= new TH2D(Form("ClusterVSBG_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), Form("ClusterVSBG_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), 100, 0.3, 1e4, 150, 10, 160); 
      BinLogX(clusterVSBG[i][j]->GetXaxis());
      
      clusterVSBG_ele[i][j]= new TH2D(Form("ClusterVSBG_ele_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), Form("ClusterVSBG_ele_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), 100, 0.3, 1e4, 150, 10, 160); 
      BinLogX(clusterVSBG_ele[i][j]->GetXaxis());
      clusterVSBG_pion[i][j]= new TH2D(Form("ClusterVSBG_pion_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), Form("ClusterVSBG_pion_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), 100, 0.3, 1e4, 150, 10, 160); 
      BinLogX(clusterVSBG_pion[i][j]->GetXaxis());
      clusterVSBG_proton[i][j]= new TH2D(Form("ClusterVSBG_proton_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), Form("ClusterVSBG_proton_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), 100, 0.3, 1e4, 150, 10, 160); 
      BinLogX(clusterVSBG_proton[i][j]->GetXaxis());
      
	/*
	MPVclusterVSBG[i][j]= new TH1D(Form("MPVClusterVSBG_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), Form("MPVClusterVSBG_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]),  100, 0.3, 1e4); 
	BinLogX(MPVclusterVSBG[i][j]->GetXaxis());

	MPVclusterVSBG_ele[i][j]= new TH1D(Form("MPVClusterVSBG_ele_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), Form("MPVClusterVSBG_ele_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]),  100, 0.3, 1e4); 
	BinLogX(MPVclusterVSBG_ele[i][j]->GetXaxis());
	MPVclusterVSBG_pion[i][j]= new TH1D(Form("MPVClusterVSBG_pion_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), Form("MPVClusterVSBG_pion_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]),  100, 0.3, 1e4); 
	BinLogX(MPVclusterVSBG_pion[i][j]->GetXaxis());
	MPVclusterVSBG_proton[i][j]= new TH1D(Form("MPVClusterVSBG_proton_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), Form("MPVClusterVSBG_proton_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]),  100, 0.3, 1e4); 
	BinLogX(MPVclusterVSBG_proton[i][j]->GetXaxis());
	*/


	//
      GenClusVSBG[i][j]= new TH2D(Form("GenClusVSBG_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), Form("GenClusVSBG_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), 100, 0.3, 1e4, 150, 10, 160); 
      BinLogX(GenClusVSBG[i][j]->GetXaxis());
      GenClusVSBG_ele[i][j]= new TH2D(Form("GenClusVSBG_ele_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), Form("GenClusVSBG_ele_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), 100, 0.3, 1e4, 150, 10, 160); 
      BinLogX(GenClusVSBG_ele[i][j]->GetXaxis());
      GenClusVSBG_pion[i][j]= new TH2D(Form("GenClusVSBG_pion_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), Form("GenClusVSBG_pion_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), 100, 0.3, 1e4, 150, 10, 160); 
      BinLogX(GenClusVSBG_pion[i][j]->GetXaxis());
      GenClusVSBG_proton[i][j]= new TH2D(Form("GenClusVSBG_proton_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), Form("GenClusVSBG_proton_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), 100, 0.3, 1e4, 150, 10, 160); 
      BinLogX(GenClusVSBG_proton[i][j]->GetXaxis());
	/*
	MPVGenClusVSBG[i][j]= new TH1D(Form("MPVGenClusVSBG_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), Form("MPVGenClusVSBG_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]),  100, 0.3, 1e4); 
	BinLogX(MPVGenClusVSBG[i][j]->GetXaxis());

	MPVGenClusVSBG_ele[i][j]= new TH1D(Form("MPVGenClusVSBG_ele_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), Form("MPVGenClusVSBG_ele_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]),  100, 0.3, 1e4); 
	BinLogX(MPVGenClusVSBG_ele[i][j]->GetXaxis());
	MPVGenClusVSBG_pion[i][j]= new TH1D(Form("MPVGenClusVSBG_pion_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), Form("MPVGenClusVSBG_pion_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]),  100, 0.3, 1e4); 
	BinLogX(MPVGenClusVSBG_pion[i][j]->GetXaxis());
	MPVGenClusVSBG_proton[i][j]= new TH1D(Form("MPVGenClusVSBG_proton_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), Form("MPVGenClusVSBG_proton_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]),  100, 0.3, 1e4); 
	BinLogX(MPVGenClusVSBG_proton[i][j]->GetXaxis());
	*/

      
      meanSigVSBG[i][j]= new TH2D(Form("MeanSigVSBG_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), Form("MeanSigVSBG_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), 100, 0.3, 1e4, 200, 0, 10); 
      BinLogX(meanSigVSBG[i][j]->GetXaxis());
      
      meanSigVSBG_ele[i][j]= new TH2D(Form("MeanSigVSBG_ele_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), Form("MeanSigVSBG_ele_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), 100, 0.3, 1e4, 200, 0, 10); 
      BinLogX(meanSigVSBG_ele[i][j]->GetXaxis());
      meanSigVSBG_pion[i][j]= new TH2D(Form("MeanSigVSBG_pion_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), Form("MeanSigVSBG_pion_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), 100, 0.3, 1e4, 200, 0, 10); 
      BinLogX(meanSigVSBG_pion[i][j]->GetXaxis());
      meanSigVSBG_proton[i][j]= new TH2D(Form("MeanSigVSBG_proton_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), Form("MeanSigVSBG_proton_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), 100, 0.3, 1e4, 200, 0, 10); 
      BinLogX(meanSigVSBG_proton[i][j]->GetXaxis());
      /*
	MPVmeanSigVSBG[i][j]= new TH1D(Form("MPVMeanSigVSBG_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), Form("MPVMeanSigVSBG_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]),  100, 0.3, 1e4); 
	BinLogX(MPVmeanSigVSBG[i][j]->GetXaxis());


	MPVmeanSigVSBG_ele[i][j]= new TH1D(Form("MPVMeanSigVSBG_ele_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), Form("MPVMeanSigVSBG_ele_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]),  100, 0.3, 1e4); 
	BinLogX(MPVmeanSigVSBG_ele[i][j]->GetXaxis());
	MPVmeanSigVSBG_pion[i][j]= new TH1D(Form("MPVMeanSigVSBG_pion_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), Form("MPVMeanSigVSBG_pion_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]),  100, 0.3, 1e4); 
	BinLogX(MPVmeanSigVSBG_pion[i][j]->GetXaxis());
	MPVmeanSigVSBG_proton[i][j]= new TH1D(Form("MPVMeanSigVSBG_proton_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), Form("MPVMeanSigVSBG_proton_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]),  100, 0.3, 1e4); 
	BinLogX(MPVmeanSigVSBG_proton[i][j]->GetXaxis());

	*/

      EtaTPCtglVSBG[i][j]= new TH2D(Form("EtaTPCtglVSBG_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), Form("EtaTPCtglVSBG_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), 100, 0.3, 1e4, 180, -1, 1); 
      BinLogX(EtaTPCtglVSBG[i][j]->GetXaxis());
      
      EtaTPCtglVSBG_ele[i][j]= new TH2D(Form("EtaTPCtglVSBG_ele_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), Form("EtaTPCtglVSBG_ele_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), 100, 0.3, 1e4, 180, -1, 1); 
      BinLogX(EtaTPCtglVSBG_ele[i][j]->GetXaxis());
      EtaTPCtglVSBG_pion[i][j]= new TH2D(Form("EtaTPCtglVSBG_pion_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), Form("EtaTPCtglVSBG_pion_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), 100, 0.3, 1e4, 180, -1, 1); 
      BinLogX(EtaTPCtglVSBG_pion[i][j]->GetXaxis());
      EtaTPCtglVSBG_proton[i][j]= new TH2D(Form("EtaTPCtglVSBG_proton_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), Form("EtaTPCtglVSBG_proton_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), 100, 0.3, 1e4, 180, -1, 1); 
      BinLogX(EtaTPCtglVSBG_proton[i][j]->GetXaxis());

	/*
	MPVEtaTPCtglVSBG[i][j]= new TH1D(Form("MPVEtaTPCtglVSBG_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), Form("MPVEtaTPCtglVSBG_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]),  100, 0.3, 1e4); 
	BinLogX(MPVEtaTPCtglVSBG[i][j]->GetXaxis());


	MPVEtaTPCtglVSBG_ele[i][j]= new TH1D(Form("MPVEtaTPCtglVSBG_ele_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), Form("MPVEtaTPCtglVSBG_ele_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]),  100, 0.3, 1e4); 
	BinLogX(MPVEtaTPCtglVSBG_ele[i][j]->GetXaxis());

	MPVEtaTPCtglVSBG_pion[i][j]= new TH1D(Form("MPVEtaTPCtglVSBG_pion_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), Form("MPVEtaTPCtglVSBG_pion_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]),  100, 0.3, 1e4); 
	BinLogX(MPVEtaTPCtglVSBG_pion[i][j]->GetXaxis());

	MPVEtaTPCtglVSBG_proton[i][j]= new TH1D(Form("MPVEtaTPCtglVSBG_proton_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]), Form("MPVEtaTPCtglVSBG_proton_%iNch_%f-%fCent", i+4,CentBins[j], CentBins[j+1]),  100, 0.3, 1e4); 
	BinLogX(MPVEtaTPCtglVSBG_proton[i][j]->GetXaxis());

	*/
    }
    for (int j=0; j<NEtaBins; j++) {
      meanSigVSBGEta[i][j]= new TH2D(Form("MeanSigVSBGEta_%iNch_%f-%fEta", i+4,EtaBins[j], EtaBins[j+1]), Form("MeanSigVSBGEta_%iNch_%f-%fEta", i+4,EtaBins[j], EtaBins[j+1]), 100, 0.3, 1e4, 200, 0, 10); 
      BinLogX(meanSigVSBGEta[i][j]->GetXaxis());
      meanSigVSBGEta_ele[i][j]= new TH2D(Form("MeanSigVSBGEta_ele_%iNch_%f-%fEta", i+4,EtaBins[j], EtaBins[j+1]), Form("MeanSigVSBGEta_ele_%iNch_%f-%fEta", i+4,EtaBins[j], EtaBins[j+1]), 100, 0.3, 1e4, 200, 0, 10); 
      BinLogX(meanSigVSBGEta_ele[i][j]->GetXaxis());
      meanSigVSBGEta_pion[i][j]= new TH2D(Form("MeanSigVSBGEta_pion_%iNch_%f-%fEta", i+4,EtaBins[j], EtaBins[j+1]), Form("MeanSigVSBGEta_pion_%iNch_%f-%fEta", i+4,EtaBins[j], EtaBins[j+1]), 100, 0.3, 1e4, 200, 0, 10); 
      BinLogX(meanSigVSBGEta_pion[i][j]->GetXaxis());
      meanSigVSBGEta_proton[i][j]= new TH2D(Form("MeanSigVSBGEta_proton_%iNch_%f-%fEta", i+4,EtaBins[j], EtaBins[j+1]), Form("MeanSigVSBGEta_proton_%iNch_%f-%fEta", i+4,EtaBins[j], EtaBins[j+1]), 100, 0.3, 1e4, 200, 0, 10); 
      BinLogX(meanSigVSBGEta_proton[i][j]->GetXaxis());

	/*
	MPVmeanSigVSBGEta[i][j]= new TH1D(Form("MPVMeanSigVSBGEta_%iNch_%f-%fEta", i+4,EtaBins[j], EtaBins[j+1]), Form("MPVMeanSigVSBGEta_%iNch_%f-%fEta", i+4,EtaBins[j], EtaBins[j+1]),  100, 0.3, 1e4); 
	BinLogX(MPVmeanSigVSBGEta[i][j]->GetXaxis());


	MPVmeanSigVSBGEta_ele[i][j]= new TH1D(Form("MPVMeanSigVSBGEta_ele_%iNch_%f-%fEta", i+4,EtaBins[j], EtaBins[j+1]), Form("MPVMeanSigVSBGEta_ele_%iNch_%f-%fEta", i+4,EtaBins[j], EtaBins[j+1]),  100, 0.3, 1e4); 
	BinLogX(MPVmeanSigVSBGEta_ele[i][j]->GetXaxis());
	MPVmeanSigVSBGEta_pion[i][j]= new TH1D(Form("MPVMeanSigVSBGEta_pion_%iNch_%f-%fEta", i+4,EtaBins[j], EtaBins[j+1]), Form("MPVMeanSigVSBGEta_pion_%iNch_%f-%fEta", i+4,EtaBins[j], EtaBins[j+1]),  100, 0.3, 1e4); 
	BinLogX(MPVmeanSigVSBGEta_pion[i][j]->GetXaxis());
	MPVmeanSigVSBGEta_proton[i][j]= new TH1D(Form("MPVMeanSigVSBGEta_proton_%iNch_%f-%fEta", i+4,EtaBins[j], EtaBins[j+1]), Form("MPVMeanSigVSBGEta_proton_%iNch_%f-%fEta", i+4,EtaBins[j], EtaBins[j+1]),  100, 0.3, 1e4); 
	BinLogX(MPVmeanSigVSBGEta_proton[i][j]->GetXaxis());
	*/
    }
  }


  // PbPb Test
  TH2D * TRDCentrality = new TH2D("TRDCentrality", "TRDCentrality", NCentBins, CentBins, 200, 0, 10 );
  TH2D * TPCCentrality = new TH2D("TPCCentrality", "TPCCentrality", NCentBins, CentBins, 200,-5,5);
  TH2D * TOFCentrality = new TH2D("TOFCentrality", "TOFCentrality", NCentBins, CentBins, 200, -5, 5);
  TH2F * TPCEta = new TH2F("TPCEta", "TPCEta",  180, -1.1,  1.1, 200, -5,5);
  TH2F * TOFEta = new TH2F("TOFEta", "TOfEta", 180, -1.1, 1.1, 200, -5,5);
  TH2F * TPCEta_ele = new TH2F("TPCEta_ele", "TPCEta_ele",  180, -1.1,  1.1, 200, -5,5);
  TH2F * TOFEta_ele = new TH2F("TOFEta_ele", "TOfEta_ele", 180, -1.1, 1.1, 200, -5,5);
  TH2F * TPCEta_pion = new TH2F("TPCEta_pion", "TPCEta_pion",  180, -1.1,  1.1, 200, -5,5);
  TH2F * TOFEta_pion = new TH2F("TOFEta_pion", "TOfEta_pion", 180, -1.1, 1.1, 200, -5,5);
  TH2F * TPCEta_proton = new TH2F("TPCEta_proton", "TPCEta_proton",  180, -1.1,  1.1, 200, -5,5);
  TH2F * TOFEta_proton = new TH2F("TOFEta_proton", "TOfEta_proton", 180, -1.1, 1.1, 200, -5,5);
  TH2D * DummyEta = new TH2D("DummyEta", "DummyEta", NEtaBins, EtaBins, 200, 0, 100);
  


  /////////////////////////////////////////////////////////////////////
  
  /////////////// MAIN ANALYSIS //////////////////////////////////////

  ////////////////////////////////////////////////////////////////////


  // Check if corrections are activated
  // eta corr. map - NOTICE: Not really eta map - instead closely related TPCtgl variable is used
  TString EtaSpec, NclsCorr, CentCorr;
  if (!NclsCorrection) NclsCorr="";
  else NclsCorr=NclsCorrStr;
  if (!fEtaCorrection) EtaSpec ="";
  else EtaSpec = EtaSpecStr;
  if (!CentCorrection) CentCorr="";
  else CentCorr=CentCorrStr;

  // loading correction maps. slightly confusing structure at the moment due to iterative structure - make sure you running the code in the right order
  TH2D * etaMap;
  if (fEtaCorrection) {
    //TString tmp(NclsCorr);
    //NclsCorr="";
    TFile * etaMapFile = TFile::Open(Form("output/%s_TPCtglMap_%s_%s%s.root", inputSpec.Data(), FileSpec.Data(),CentCorr.Data(), NclsCorr.Data()), "READ");
    //etaMap = (TH2D*)etaMapFile->Get("map");
    etaMap = (TH2D*)etaMapFile->Get("mapTPCtgl");
    //NclsCorr=tmp;
  }
   
  TH2D * NclsMap[3];
  if (NclsCorrection) {
    TString tmp(EtaSpec);
    EtaSpec="";
    TFile * ClusterMapFile = TFile::Open(Form("output/%s_ClusterMap_%s_%s%s.root", inputSpec.Data(),    FileSpec.Data(), CentCorr.Data(), EtaSpec.Data()), "READ");
    for (int i=0; i<3; i++) {
      NclsMap[i] = (TH2D*)ClusterMapFile->Get(Form("NclsMap_nCh%i", i+4));
    }
      EtaSpec=tmp;
  }
    
  TH2D * CentMap;
  if (CentCorrection) {
    TString tmp(EtaSpec);
    TString tmp2(NclsCorr);
    EtaSpec="";
    NclsCorr="";
    TFile * CentMapFile = TFile::Open(Form("output/%s_CentMap_%s_%s%s.root", inputSpec.Data(),   FileSpec.Data(), EtaSpec.Data(), NclsCorr.Data()), "READ");
    CentMap = (TH2D*)CentMapFile->Get("CentMap");
    EtaSpec=tmp;
    NclsCorr=tmp2;
  }

    
  //===========================================================================================
  //
  // mean parametrization
  //
  //===========================================================================================
    
  // prepare signal vs beta*gamma histogram
  BinLogX(meanSigVsBG->GetXaxis()); // log binning
  BinLogX(meanSigVsBG_ele->GetXaxis()); // log binning
  BinLogX(meanSigVsBG_pion->GetXaxis()); // log binning
  BinLogX(meanSigVsBG_proton->GetXaxis()); // log binning
  meanSigVsBG->SetMinimum(1); // minimum

  // prepare signal vs beta*gamma histogram
  BinLogX(meanSigVsBG_minus->GetXaxis()); // log binning
  BinLogX(meanSigVsBG_minus_ele->GetXaxis()); // log binning
  BinLogX(meanSigVsBG_minus_pion->GetXaxis()); // log binning
  BinLogX(meanSigVsBG_minus_proton->GetXaxis()); // log binning
  meanSigVsBG_minus->SetMinimum(1); // minimum

  // prepare signal vs beta*gamma histogram
  BinLogX(meanSigVsBG_plus->GetXaxis()); // log binning
  BinLogX(meanSigVsBG_plus_ele->GetXaxis()); // log binning
  BinLogX(meanSigVsBG_plus_pion->GetXaxis()); // log binning
  BinLogX(meanSigVsBG_plus_proton->GetXaxis()); // log binning
  meanSigVsBG_plus->SetMinimum(1); // minimum
    
  // loop over tree to fill signal vs beta*gamma histogram
  Int_t nb = 0;
  Int_t run = 0;

  cout << "Number of Entries: " <<  nEntries << endl; 
  for (Int_t i = 0; i < nEntries; i++) { // loop over tree
    nb = tracks.getEntry(i);
    if (nb <= 0) continue;
    if (i%1000000==0) cout << i << " " <<  tracks.run() << endl;
    
    if (runCut(tracks.run())) continue; // skip specified runs
    
    if (TMath::Abs(tracks.trdTPCtgl())>1.) continue; // no correction for higher eta available
    // similar cuts for cluster or centrality?

    // PbPbTest
    // if (tracks.centrality()>=100) continue;

    // PID cuts for e, pi, proton - w/o exclusion cut for protons 
    if (nSigmaCut(tracks.nSigmaTOF(), (Double_t*)cutNSigmaTOF, TMath::Abs(tracks.pId()))) continue;
    if (nSigmaCut(tracks.nSigmaTPC(), (Double_t*)cutNSigmaTPC, TMath::Abs(tracks.pId()))) continue;

    if (QAPlots) {
      if (TMath::Abs(tracks.pId()) == pdgCodes[0]) { // electron
	TRDCentrality->Fill( tracks.centrality(), tracks.trdSig());
	clusterVsCentrality->Fill(tracks.centrality(), tracks.trdNCls());
	clusterVsChamber->Fill(tracks.trdNCls(), tracks.trdNCh()); // histogram without cuts!

	TPCEta->Fill(tracks.trdTPCtgl(), tracks.nSigmaTPC()[0]);
	TOFEta->Fill(tracks.trdTPCtgl(), tracks.nSigmaTOF()[0]);
	TPCEta_ele->Fill(tracks.trdTPCtgl(), tracks.nSigmaTPC()[0]);
	TOFEta_ele->Fill(tracks.trdTPCtgl(), tracks.nSigmaTOF()[0]);
	
	TPCCentrality->Fill(tracks.centrality(), tracks.nSigmaTPC()[0]); //, tracks.centrality());
	TOFCentrality->Fill(tracks.centrality(), tracks.nSigmaTOF()[0]);
	clusterVsChamber_ele->Fill(tracks.trdNCls(), tracks.trdNCh()); // histogram without cuts!
	if (tracks.trdNCh()>3) {
	  ClusterVSChamber_ele[tracks.trdNCh()-4][TOFCentrality->GetXaxis()->FindBin(tracks.centrality())-1]->Fill(tracks.trdNCls());
	  clusterVsCentrality_ele->Fill(tracks.centrality(), tracks.trdNCls());
	}
      }
      if (TMath::Abs(tracks.pId()) == pdgCodes[1]) { // pion
	TRDCentrality->Fill( tracks.centrality(), tracks.trdSig());
	clusterVsCentrality->Fill(tracks.centrality(), tracks.trdNCls());
	clusterVsChamber->Fill(tracks.trdNCls(), tracks.trdNCh()); // histogram without cuts!

	TPCEta->Fill(tracks.trdTPCtgl(), tracks.nSigmaTPC()[1]);
	TOFEta->Fill(tracks.trdTPCtgl(), tracks.nSigmaTOF()[1]);
	TPCEta_pion->Fill(tracks.trdTPCtgl(), tracks.nSigmaTPC()[1]);
	TOFEta_pion->Fill(tracks.trdTPCtgl(), tracks.nSigmaTOF()[1]);
	
	TPCCentrality->Fill(tracks.centrality(), tracks.nSigmaTPC()[1]); //, tracks.centrality());
	TOFCentrality->Fill(tracks.centrality(), tracks.nSigmaTOF()[1]);
	clusterVsChamber_pion->Fill(tracks.trdNCls(), tracks.trdNCh()); // histogram without cuts!
	
	if (tracks.trdNCh()>3) {
	  ClusterVSChamber_pion[tracks.trdNCh()-4][TOFCentrality->GetXaxis()->FindBin(tracks.centrality())-1]->Fill(tracks.trdNCls());
	  clusterVsCentrality_pion->Fill(tracks.centrality(), tracks.trdNCls());
	}
      }
      if (TMath::Abs(tracks.pId()) == pdgCodes[2]) { // proton
	TRDCentrality->Fill( tracks.centrality(), tracks.trdSig());
	clusterVsCentrality->Fill(tracks.centrality(), tracks.trdNCls());
	clusterVsChamber->Fill(tracks.trdNCls(), tracks.trdNCh()); // histogram without cuts!
      
	TPCEta->Fill(tracks.trdTPCtgl(), tracks.nSigmaTPC()[2]);
	TOFEta->Fill(tracks.trdTPCtgl(), tracks.nSigmaTOF()[2]);
	TPCEta_proton->Fill(tracks.trdTPCtgl(), tracks.nSigmaTPC()[2]);
	TOFEta_proton->Fill(tracks.trdTPCtgl(), tracks.nSigmaTOF()[2]);
	TPCCentrality->Fill(tracks.centrality(), tracks.nSigmaTPC()[2]); //, tracks.centrality());
	TOFCentrality->Fill(tracks.centrality(), tracks.nSigmaTOF()[2]);
	clusterVsChamber_proton->Fill(tracks.trdNCls(), tracks.trdNCh()); // histogram without cuts
	if (tracks.trdNCh()>3) {
	  ClusterVSChamber_proton[tracks.trdNCh()-4][TOFCentrality->GetXaxis()->FindBin(tracks.centrality())-1]->Fill(tracks.trdNCls());
	  clusterVsCentrality_proton->Fill(tracks.centrality(), tracks.trdNCls());
	}
      }

     
      if (TMath::Abs(tracks.pId()) == pdgCodes[0]) { // electron
	clusterVsBG->Fill( betaGamma(momentumMean(tracks.trdMom()), mElectron), tracks.trdNCls());
	clusterVsBG_ele->Fill( betaGamma(momentumMean(tracks.trdMom()), mElectron), tracks.trdNCls());
      
	GenClusVsBG->Fill( betaGamma(momentumMean(tracks.trdMom()), mElectron), tracks.trdNClsGeneral());
	GenClusVsBG_ele->Fill( betaGamma(momentumMean(tracks.trdMom()), mElectron), tracks.trdNClsGeneral());

	EtaTPCtglVsBG->Fill( betaGamma(momentumMean(tracks.trdMom()), mElectron), abs(tracks.trdTPCtgl()));
	EtaTPCtglVsBG_ele->Fill( betaGamma(momentumMean(tracks.trdMom()), mElectron), abs(tracks.trdTPCtgl()));

	//  cout << tracks.trdNCh()-4 << TOFCentrality->GetXaxis()->FindBin(tracks.centrality()) << endl;
	if (tracks.trdNCh()>3) {
	  clusterVSBG[tracks.trdNCh()-4][TOFCentrality->GetXaxis()->FindBin(tracks.centrality())-1]->Fill(betaGamma(momentumMean(tracks.trdMom()), mElectron), tracks.trdNCls());
	  clusterVSBG_ele[tracks.trdNCh()-4][TOFCentrality->GetXaxis()->FindBin(tracks.centrality())-1]->Fill(betaGamma(momentumMean(tracks.trdMom()), mElectron), tracks.trdNCls());
	  clusterVSBG_Nch[tracks.trdNCh()-4]->Fill(betaGamma(momentumMean(tracks.trdMom()), mElectron), tracks.trdNCls());
	  clusterVSBG_Nch_ele[tracks.trdNCh()-4]->Fill(betaGamma(momentumMean(tracks.trdMom()), mElectron), tracks.trdNCls());
	
	  GenClusVSBG[tracks.trdNCh()-4][TOFCentrality->GetXaxis()->FindBin(tracks.centrality())-1]->Fill(betaGamma(momentumMean(tracks.trdMom()), mElectron), tracks.trdNClsGeneral());
	  GenClusVSBG_ele[tracks.trdNCh()-4][TOFCentrality->GetXaxis()->FindBin(tracks.centrality())-1]->Fill(betaGamma(momentumMean(tracks.trdMom()), mElectron), tracks.trdNClsGeneral());
	  GenClusVSBG_Nch[tracks.trdNCh()-4]->Fill(betaGamma(momentumMean(tracks.trdMom()), mElectron), tracks.trdNClsGeneral());
	  GenClusVSBG_Nch_ele[tracks.trdNCh()-4]->Fill(betaGamma(momentumMean(tracks.trdMom()), mElectron), tracks.trdNClsGeneral());
	
	  meanSigVSBG[tracks.trdNCh()-4][TOFCentrality->GetXaxis()->FindBin(tracks.centrality())-1]->Fill(betaGamma(momentumMean(tracks.trdMom()), mElectron), tracks.trdSig());
	  meanSigVSBG_ele[tracks.trdNCh()-4][TOFCentrality->GetXaxis()->FindBin(tracks.centrality())-1]->Fill(betaGamma(momentumMean(tracks.trdMom()), mElectron), tracks.trdSig());
	  meanSigVSBGEta[tracks.trdNCh()-4][DummyEta->GetXaxis()->FindBin(tracks.trdTPCtgl())-1]->Fill(betaGamma(momentumMean(tracks.trdMom()), mElectron), tracks.trdSig());
	  meanSigVSBGEta_ele[tracks.trdNCh()-4][DummyEta->GetXaxis()->FindBin(tracks.trdTPCtgl())-1]->Fill(betaGamma(momentumMean(tracks.trdMom()), mElectron), tracks.trdSig());

	  meanSigVSBG_Nch[tracks.trdNCh()-4]->Fill(betaGamma(momentumMean(tracks.trdMom()), mElectron), tracks.trdSig());
	  meanSigVSBG_Nch_ele[tracks.trdNCh()-4]->Fill(betaGamma(momentumMean(tracks.trdMom()), mElectron), tracks.trdSig());

	  EtaTPCtglVSBG[tracks.trdNCh()-4][TOFCentrality->GetXaxis()->FindBin(tracks.centrality())-1]->Fill(betaGamma(momentumMean(tracks.trdMom()), mElectron), abs(tracks.trdTPCtgl()));
	  EtaTPCtglVSBG_ele[tracks.trdNCh()-4][TOFCentrality->GetXaxis()->FindBin(tracks.centrality())-1]->Fill(betaGamma(momentumMean(tracks.trdMom()), mElectron), abs(tracks.trdTPCtgl()));
	  EtaTPCtglVSBG_Nch[tracks.trdNCh()-4]->Fill(betaGamma(momentumMean(tracks.trdMom()), mElectron), abs(tracks.trdTPCtgl()));
	  EtaTPCtglVSBG_Nch_ele[tracks.trdNCh()-4]->Fill(betaGamma(momentumMean(tracks.trdMom()), mElectron), abs(tracks.trdTPCtgl()));
	}
      }
    
      if (TMath::Abs(tracks.pId()) == pdgCodes[1]) {
	clusterVsBG->Fill( betaGamma(momentumMean(tracks.trdMom()), mPion), tracks.trdNCls());
	clusterVsBG_pion->Fill( betaGamma(momentumMean(tracks.trdMom()), mPion), tracks.trdNCls());
      
	GenClusVsBG->Fill( betaGamma(momentumMean(tracks.trdMom()), mPion), tracks.trdNClsGeneral());
	GenClusVsBG_pion->Fill( betaGamma(momentumMean(tracks.trdMom()), mPion), tracks.trdNClsGeneral());
            
	EtaTPCtglVsBG->Fill( betaGamma(momentumMean(tracks.trdMom()), mPion), abs(tracks.trdTPCtgl()));
	EtaTPCtglVsBG_pion->Fill( betaGamma(momentumMean(tracks.trdMom()), mPion), abs(tracks.trdTPCtgl()));
      
	if (tracks.trdNCh()>3) {
	  clusterVSBG[tracks.trdNCh()-4][TOFCentrality->GetXaxis()->FindBin(tracks.centrality())-1]->Fill(betaGamma(momentumMean(tracks.trdMom()), mPion), tracks.trdNCls());
	  clusterVSBG_pion[tracks.trdNCh()-4][TOFCentrality->GetXaxis()->FindBin(tracks.centrality())-1]->Fill(betaGamma(momentumMean(tracks.trdMom()), mPion), tracks.trdNCls());
	  clusterVSBG_Nch[tracks.trdNCh()-4]->Fill(betaGamma(momentumMean(tracks.trdMom()), mPion), tracks.trdNCls());
	  clusterVSBG_Nch_pion[tracks.trdNCh()-4]->Fill(betaGamma(momentumMean(tracks.trdMom()), mPion), tracks.trdNCls());
	
	  GenClusVSBG[tracks.trdNCh()-4][TOFCentrality->GetXaxis()->FindBin(tracks.centrality())-1]->Fill(betaGamma(momentumMean(tracks.trdMom()), mPion), tracks.trdNClsGeneral());
	  GenClusVSBG_pion[tracks.trdNCh()-4][TOFCentrality->GetXaxis()->FindBin(tracks.centrality())-1]->Fill(betaGamma(momentumMean(tracks.trdMom()), mPion), tracks.trdNClsGeneral());
	  GenClusVSBG_Nch[tracks.trdNCh()-4]->Fill(betaGamma(momentumMean(tracks.trdMom()), mPion), tracks.trdNClsGeneral());
	  GenClusVSBG_Nch_pion[tracks.trdNCh()-4]->Fill(betaGamma(momentumMean(tracks.trdMom()), mPion), tracks.trdNClsGeneral());
	

	  meanSigVSBG[tracks.trdNCh()-4][TOFCentrality->GetXaxis()->FindBin(tracks.centrality())-1]->Fill(betaGamma(momentumMean(tracks.trdMom()), mPion), tracks.trdSig());
	  meanSigVSBG_pion[tracks.trdNCh()-4][TOFCentrality->GetXaxis()->FindBin(tracks.centrality())-1]->Fill(betaGamma(momentumMean(tracks.trdMom()), mPion), tracks.trdSig());
	  meanSigVSBGEta[tracks.trdNCh()-4][DummyEta->GetXaxis()->FindBin(tracks.trdTPCtgl())-1]->Fill(betaGamma(momentumMean(tracks.trdMom()), mPion), tracks.trdSig());
	  meanSigVSBGEta_pion[tracks.trdNCh()-4][DummyEta->GetXaxis()->FindBin(tracks.trdTPCtgl())-1]->Fill(betaGamma(momentumMean(tracks.trdMom()), mPion), tracks.trdSig());
	  meanSigVSBG_Nch[tracks.trdNCh()-4]->Fill(betaGamma(momentumMean(tracks.trdMom()), mPion), tracks.trdSig());
	  meanSigVSBG_Nch_pion[tracks.trdNCh()-4]->Fill(betaGamma(momentumMean(tracks.trdMom()), mPion), tracks.trdSig());


	  EtaTPCtglVSBG[tracks.trdNCh()-4][TOFCentrality->GetXaxis()->FindBin(tracks.centrality())-1]->Fill(betaGamma(momentumMean(tracks.trdMom()), mPion), abs(tracks.trdTPCtgl()));
	  EtaTPCtglVSBG_pion[tracks.trdNCh()-4][TOFCentrality->GetXaxis()->FindBin(tracks.centrality())-1]->Fill(betaGamma(momentumMean(tracks.trdMom()), mPion), abs(tracks.trdTPCtgl()));
	  EtaTPCtglVSBG_Nch[tracks.trdNCh()-4]->Fill(betaGamma(momentumMean(tracks.trdMom()), mPion), abs(tracks.trdTPCtgl()));
	  EtaTPCtglVSBG_Nch_pion[tracks.trdNCh()-4]->Fill(betaGamma(momentumMean(tracks.trdMom()), mPion), abs(tracks.trdTPCtgl()));
	}	
      }	
      
      if (TMath::Abs(tracks.pId()) == pdgCodes[2]) {
	clusterVsBG->Fill( betaGamma(momentumMean(tracks.trdMom()), mProton), tracks.trdNCls());
	clusterVsBG_proton->Fill( betaGamma(momentumMean(tracks.trdMom()), mProton), tracks.trdNCls());
      
	GenClusVsBG->Fill( betaGamma(momentumMean(tracks.trdMom()), mProton), tracks.trdNClsGeneral());
	GenClusVsBG_proton->Fill( betaGamma(momentumMean(tracks.trdMom()), mProton), tracks.trdNClsGeneral());
      
	EtaTPCtglVsBG->Fill( betaGamma(momentumMean(tracks.trdMom()), mProton), abs(tracks.trdTPCtgl()));
	EtaTPCtglVsBG_proton->Fill( betaGamma(momentumMean(tracks.trdMom()), mProton), abs(tracks.trdTPCtgl()));
      
	if (tracks.trdNCh()>3) {
	  clusterVSBG[tracks.trdNCh()-4][TOFCentrality->GetXaxis()->FindBin(tracks.centrality())-1]->Fill(betaGamma(momentumMean(tracks.trdMom()), mProton), tracks.trdNCls());
	  clusterVSBG_proton[tracks.trdNCh()-4][TOFCentrality->GetXaxis()->FindBin(tracks.centrality())-1]->Fill(betaGamma(momentumMean(tracks.trdMom()), mProton), tracks.trdNCls());
	  clusterVSBG_Nch[tracks.trdNCh()-4]->Fill(betaGamma(momentumMean(tracks.trdMom()), mProton), tracks.trdNCls());
	  clusterVSBG_Nch_proton[tracks.trdNCh()-4]->Fill(betaGamma(momentumMean(tracks.trdMom()), mProton), tracks.trdNCls());
	
	  GenClusVSBG[tracks.trdNCh()-4][TOFCentrality->GetXaxis()->FindBin(tracks.centrality())-1]->Fill(betaGamma(momentumMean(tracks.trdMom()), mProton), tracks.trdNClsGeneral());
	  GenClusVSBG_proton[tracks.trdNCh()-4][TOFCentrality->GetXaxis()->FindBin(tracks.centrality())-1]->Fill(betaGamma(momentumMean(tracks.trdMom()), mProton), tracks.trdNClsGeneral());
	  GenClusVSBG_Nch[tracks.trdNCh()-4]->Fill(betaGamma(momentumMean(tracks.trdMom()), mProton), tracks.trdNClsGeneral());
	  GenClusVSBG_Nch_proton[tracks.trdNCh()-4]->Fill(betaGamma(momentumMean(tracks.trdMom()), mProton), tracks.trdNClsGeneral());
	
	  EtaTPCtglVSBG[tracks.trdNCh()-4][TOFCentrality->GetXaxis()->FindBin(tracks.centrality())-1]->Fill(betaGamma(momentumMean(tracks.trdMom()), mProton), abs(tracks.trdTPCtgl()));
	  EtaTPCtglVSBG_proton[tracks.trdNCh()-4][TOFCentrality->GetXaxis()->FindBin(tracks.centrality())-1]->Fill(betaGamma(momentumMean(tracks.trdMom()), mProton), abs(tracks.trdTPCtgl()));
	  EtaTPCtglVSBG_Nch[tracks.trdNCh()-4]->Fill(betaGamma(momentumMean(tracks.trdMom()), mProton), abs(tracks.trdTPCtgl()));
	  EtaTPCtglVSBG_Nch_proton[tracks.trdNCh()-4]->Fill(betaGamma(momentumMean(tracks.trdMom()), mProton), abs(tracks.trdTPCtgl()));

	  meanSigVSBG[tracks.trdNCh()-4][TOFCentrality->GetXaxis()->FindBin(tracks.centrality())-1]->Fill(betaGamma(momentumMean(tracks.trdMom()), mProton), tracks.trdSig());
	  meanSigVSBG_proton[tracks.trdNCh()-4][TOFCentrality->GetXaxis()->FindBin(tracks.centrality())-1]->Fill(betaGamma(momentumMean(tracks.trdMom()), mProton), tracks.trdSig());
	  meanSigVSBGEta[tracks.trdNCh()-4][DummyEta->GetXaxis()->FindBin(tracks.trdTPCtgl())-1]->Fill(betaGamma(momentumMean(tracks.trdMom()), mProton), tracks.trdSig());
	  meanSigVSBGEta_proton[tracks.trdNCh()-4][DummyEta->GetXaxis()->FindBin(tracks.trdTPCtgl())-1]->Fill(betaGamma(momentumMean(tracks.trdMom()), mProton), tracks.trdSig());
	  meanSigVSBG_Nch[tracks.trdNCh()-4]->Fill(betaGamma(momentumMean(tracks.trdMom()), mProton), tracks.trdSig());
	  meanSigVSBG_Nch_proton[tracks.trdNCh()-4]->Fill(betaGamma(momentumMean(tracks.trdMom()), mProton), tracks.trdSig());
	}
      }
    }
    //================================================================================================

    // cut for number of tracklets, cluster
    if (trackAndLayerCuts(tracks.trdNCh(), tracks.trdNCls(), TMath::Abs(tracks.pId()))) continue;

    // geometric cuts: LocalTheta, YCut - at the moment emty as no access to this data on ESD level
    if (geometricCuts(tracks.trdEtaLocal(), tracks.trdNCh(), tracks.trdY())) continue;
        
    //if (momentumMean(tracks.trdMom())<0.8) continue; // test without low momentum particles (problems for Aleph parametrization) -- not really helpful
	
    // declare eta&ncls correction factor TH1F* LHC13P2_Cluster=0;
    Double_t CorrectionFactor = 1;
        
    // electrons
    if (TMath::Abs(tracks.pId()) == pdgCodes[0]) {

      if (electronMomentumCut(tracks.trdMom(), tracks.nSigmaTPC()[0])) continue;

      if (fEtaCorrection) CorrectionFactor *= EtaCorrFactorTPCtgl(betaGamma(momentumMean(tracks.trdMom()), mElectron), tracks.trdTPCtgl(), etaMap); 
      // for TPCtgl
      if (NclsCorrection) CorrectionFactor *= NclsCorrFactor(betaGamma(momentumMean(tracks.trdMom()), mElectron), tracks.trdNCls(), NclsMap[tracks.trdNCh()-4]);
      if (CentCorrection) CorrectionFactor *= CentCorrFactor(betaGamma(momentumMean(tracks.trdMom()), mElectron), tracks.centrality(), CentMap);
      if (CorrectionFactor<1E-5) continue;

      if ((tracks.pId()) == pdgCodes[0]) // positron
	{
	  meanSigVsBG_plus->Fill(betaGamma(momentumMean(tracks.trdMom()), mElectron), CorrectionFactor*tracks.trdSig());
	  meanSigVsBG_plus_ele->Fill(betaGamma(momentumMean(tracks.trdMom()), mElectron), CorrectionFactor*tracks.trdSig());
	}
      else // electron
	{
	  meanSigVsBG_minus->Fill(betaGamma(momentumMean(tracks.trdMom()), mElectron), CorrectionFactor*tracks.trdSig());
	  meanSigVsBG_minus_ele->Fill(betaGamma(momentumMean(tracks.trdMom()), mElectron), CorrectionFactor*tracks.trdSig());
	}
      meanSigVsBG->Fill(betaGamma(momentumMean(tracks.trdMom()), mElectron), CorrectionFactor*tracks.trdSig());
      meanSigVsBG_ele->Fill(betaGamma(momentumMean(tracks.trdMom()), mElectron), CorrectionFactor*tracks.trdSig());
      clusterVsChamberWCut_ele->Fill(tracks.trdNCls(), tracks.trdNCh()); // histogram without cuts!
    }
        
    // pions
    if (TMath::Abs(tracks.pId()) == pdgCodes[1]) {
      if (pionMomentumCut(tracks.trdMom())) continue;
      
      if (fEtaCorrection) CorrectionFactor *= EtaCorrFactorTPCtgl(betaGamma(momentumMean(tracks.trdMom()), mPion), tracks.trdTPCtgl(), etaMap);
      if (NclsCorrection) CorrectionFactor *= NclsCorrFactor(betaGamma(momentumMean(tracks.trdMom()), mPion), tracks.trdNCls(), NclsMap[tracks.trdNCh()-4]);
      if (CentCorrection) CorrectionFactor *= CentCorrFactor(betaGamma(momentumMean(tracks.trdMom()), mPion), tracks.centrality(), CentMap);
      if (CorrectionFactor==0) continue;

      if (tracks.pId()==pdgCodes[1]) // pi+
	{
	  meanSigVsBG_plus->Fill(betaGamma(momentumMean(tracks.trdMom()), mPion), CorrectionFactor*tracks.trdSig());
	  meanSigVsBG_plus_pion->Fill(betaGamma(momentumMean(tracks.trdMom()), mPion), CorrectionFactor*tracks.trdSig());
	}
      else // pi-
	{
	  meanSigVsBG_minus->Fill(betaGamma(momentumMean(tracks.trdMom()), mPion), CorrectionFactor*tracks.trdSig());
	  meanSigVsBG_minus_pion->Fill(betaGamma(momentumMean(tracks.trdMom()), mPion), CorrectionFactor*tracks.trdSig());
	}
      meanSigVsBG->Fill(betaGamma(momentumMean(tracks.trdMom()), mPion), CorrectionFactor*tracks.trdSig());
      meanSigVsBG_pion->Fill(betaGamma(momentumMean(tracks.trdMom()), mPion), CorrectionFactor*tracks.trdSig());
      clusterVsChamberWCut_pion->Fill(tracks.trdNCls(), tracks.trdNCh()); // histogram without cuts!
    }
        
    // protons
    if (pidOpt == 3 && TMath::Abs(tracks.pId()) == pdgCodes[2]) {
      if (protonMomentumCut(tracks.trdMom())) continue;
      
      if (fEtaCorrection) CorrectionFactor *= EtaCorrFactorTPCtgl(betaGamma(momentumMean(tracks.trdMom()), mProton), tracks.trdTPCtgl(), etaMap);
      if (NclsCorrection) CorrectionFactor *= NclsCorrFactor(betaGamma(momentumMean(tracks.trdMom()), mProton), tracks.trdNCls(), NclsMap[tracks.trdNCh()-4]);
      if (CentCorrection) CorrectionFactor *= CentCorrFactor(betaGamma(momentumMean(tracks.trdMom()), mProton), tracks.centrality(), CentMap);
      if (CorrectionFactor==0) continue;

      if (tracks.pId()==pdgCodes[2])  //proton
	{
	  meanSigVsBG_plus->Fill(betaGamma(momentumMean(tracks.trdMom()), mProton), CorrectionFactor*tracks.trdSig());
	  meanSigVsBG_plus_proton->Fill(betaGamma(momentumMean(tracks.trdMom()), mProton), CorrectionFactor*tracks.trdSig());
	}
      else
	{
	  meanSigVsBG_minus->Fill(betaGamma(momentumMean(tracks.trdMom()), mProton), CorrectionFactor*tracks.trdSig());
	  meanSigVsBG_minus_proton->Fill(betaGamma(momentumMean(tracks.trdMom()), mProton), CorrectionFactor*tracks.trdSig());
	}
      meanSigVsBG->Fill(betaGamma(momentumMean(tracks.trdMom()), mProton), CorrectionFactor*tracks.trdSig());
      meanSigVsBG_proton->Fill(betaGamma(momentumMean(tracks.trdMom()), mProton), CorrectionFactor*tracks.trdSig());
      clusterVsChamberWCut_proton->Fill(tracks.trdNCls(), tracks.trdNCh()); // histogram without cuts!
    }
    
    clusterVsChamberWCut->Fill(tracks.trdNCls(), tracks.trdNCh()); // clusterVsChamber histogram after cuts!
  }

  // add clusterVsChamberHistograms
  if (QAPlots) {
    ClusterQAList->Add(clusterVsChamber);
    ClusterQAList->Add(clusterVsChamberWCut);
    ClusterQAList->Add(clusterVsChamber_ele);
    ClusterQAList->Add(clusterVsChamberWCut_ele);
    ClusterQAList->Add(clusterVsChamber_pion);
    ClusterQAList->Add(clusterVsChamberWCut_pion);
    ClusterQAList->Add(clusterVsChamber_proton);
    ClusterQAList->Add(clusterVsChamberWCut_proton);
  
    for (int i=0; i<3; i++) {
      for (int j=0; j<NCentBins; j++) {
	ClusterQAList->Add(ClusterVSChamber_ele[i][j]);
	ClusterQAList->Add(ClusterVSChamber_pion[i][j]);
	ClusterQAList->Add(ClusterVSChamber_proton[i][j]);
      }
    }
    ClusterQAList->Add(clusterVsCentrality);
    ClusterQAList->Add(clusterVsCentrality_ele);
    ClusterQAList->Add(clusterVsCentrality_pion);
    ClusterQAList->Add(clusterVsCentrality_proton);
    
    // add ClusterVSBG Histograms
    ClusterQAList->Add(clusterVsBG);
    ClusterQAList->Add(clusterVsBG_ele);
    ClusterQAList->Add(clusterVsBG_pion);
    ClusterQAList->Add(clusterVsBG_proton);
    for (int i=0; i<3; i++) {
      ClusterQAList->Add(clusterVSBG_Nch[i]);
      ClusterQAList->Add(clusterVSBG_Nch_ele[i]);
      ClusterQAList->Add(clusterVSBG_Nch_pion[i]);
      ClusterQAList->Add(clusterVSBG_Nch_proton[i]); 
      for (int j=0; j<NCentBins; j++) {
	ClusterQAList->Add(clusterVSBG[i][j]);
	ClusterQAList->Add(clusterVSBG_ele[i][j]);
	ClusterQAList->Add(clusterVSBG_pion[i][j]);
	ClusterQAList->Add(clusterVSBG_proton[i][j]);
      }
    }
    
    for (int i=0; i<3; i++) {
      FitSlicesY(clusterVSBG_Nch[i], meanNor, MPVclusterVSBG_Nch[i], meanWidth, meanRes, meanChi, 10,20, dummyList);
      FitSlicesY(clusterVSBG_Nch_ele[i], meanNor, MPVclusterVSBG_Nch_ele[i], meanWidth, meanRes, meanChi, 10, 20, dummyList);
      FitSlicesY(clusterVSBG_Nch_pion[i], meanNor, MPVclusterVSBG_Nch_pion[i], meanWidth, meanRes, meanChi, 10, 20, dummyList);
      FitSlicesY(clusterVSBG_Nch_proton[i], meanNor, MPVclusterVSBG_Nch_proton[i], meanWidth, meanRes, meanChi, 10, 20, dummyList);
      ClusterQAList->Add(MPVclusterVSBG_Nch[i]);
      ClusterQAList->Add(MPVclusterVSBG_Nch_ele[i]);
      ClusterQAList->Add(MPVclusterVSBG_Nch_pion[i]);
      ClusterQAList->Add(MPVclusterVSBG_Nch_proton[i]);
      
      for (int j=0; j<NCentBins; j++) {
	FitSlicesY(clusterVSBG[i][j], meanNor, MPVclusterVSBG[i][j], meanWidth, meanRes, meanChi, 10,20, dummyList);
	FitSlicesY(clusterVSBG_ele[i][j], meanNor, MPVclusterVSBG_ele[i][j], meanWidth, meanRes, meanChi, 10, 20, dummyList);
	FitSlicesY(clusterVSBG_pion[i][j], meanNor, MPVclusterVSBG_pion[i][j], meanWidth, meanRes, meanChi, 10, 20, dummyList);
	FitSlicesY(clusterVSBG_proton[i][j], meanNor, MPVclusterVSBG_proton[i][j], meanWidth, meanRes, meanChi, 10, 20, dummyList);
	ClusterQAList->Add(MPVclusterVSBG[i][j]);
	ClusterQAList->Add(MPVclusterVSBG_ele[i][j]);
	ClusterQAList->Add(MPVclusterVSBG_pion[i][j]);
	ClusterQAList->Add(MPVclusterVSBG_proton[i][j]);
      }
    }

    // add GenClusVSBG Histograms
    ClusterQAList->Add(GenClusVsBG);
    ClusterQAList->Add(GenClusVsBG_ele);
    ClusterQAList->Add(GenClusVsBG_pion);
    ClusterQAList->Add(GenClusVsBG_proton);
    for (int i=0; i<3; i++) {
      ClusterQAList->Add(GenClusVSBG_Nch[i]);
      ClusterQAList->Add(GenClusVSBG_Nch_ele[i]);
      ClusterQAList->Add(GenClusVSBG_Nch_pion[i]);
      ClusterQAList->Add(GenClusVSBG_Nch_proton[i]);
      
      for (int j=0; j<NCentBins; j++) {
	ClusterQAList->Add(GenClusVSBG[i][j]);
	ClusterQAList->Add(GenClusVSBG_ele[i][j]);
	ClusterQAList->Add(GenClusVSBG_pion[i][j]);
	ClusterQAList->Add(GenClusVSBG_proton[i][j]);
      }
    }
    
    for (int i=0; i<3; i++) {
      FitSlicesY(GenClusVSBG_Nch[i], meanNor, MPVGenClusVSBG_Nch[i], meanWidth, meanRes, meanChi, 10,20, dummyList);
      FitSlicesY(GenClusVSBG_Nch_ele[i], meanNor, MPVGenClusVSBG_Nch_ele[i], meanWidth, meanRes, meanChi, 10, 20, dummyList);
      FitSlicesY(GenClusVSBG_Nch_pion[i], meanNor, MPVGenClusVSBG_Nch_pion[i], meanWidth, meanRes, meanChi, 10, 20, dummyList);
      FitSlicesY(GenClusVSBG_Nch_proton[i], meanNor, MPVGenClusVSBG_Nch_proton[i], meanWidth, meanRes, meanChi, 10, 20, dummyList);
      ClusterQAList->Add(MPVGenClusVSBG_Nch[i]);
      ClusterQAList->Add(MPVGenClusVSBG_Nch_ele[i]);
      ClusterQAList->Add(MPVGenClusVSBG_Nch_pion[i]);
      ClusterQAList->Add(MPVGenClusVSBG_Nch_proton[i]);
      
      for (int j=0; j<NCentBins; j++) {
	FitSlicesY(GenClusVSBG[i][j], meanNor, MPVGenClusVSBG[i][j], meanWidth, meanRes, meanChi, 10,20, dummyList);
	FitSlicesY(GenClusVSBG_ele[i][j], meanNor, MPVGenClusVSBG_ele[i][j], meanWidth, meanRes, meanChi, 10, 20, dummyList);
	FitSlicesY(GenClusVSBG_pion[i][j], meanNor, MPVGenClusVSBG_pion[i][j], meanWidth, meanRes, meanChi, 10, 20, dummyList);
	FitSlicesY(GenClusVSBG_proton[i][j], meanNor, MPVGenClusVSBG_proton[i][j], meanWidth, meanRes, meanChi, 10, 20, dummyList);
	ClusterQAList->Add(MPVGenClusVSBG[i][j]);
	ClusterQAList->Add(MPVGenClusVSBG_ele[i][j]);
	ClusterQAList->Add(MPVGenClusVSBG_pion[i][j]);
	ClusterQAList->Add(MPVGenClusVSBG_proton[i][j]);
      }
    }
    
    // add EtaTPCtglVSBG Histograms
    EtaQAList->Add(EtaTPCtglVsBG);
    EtaQAList->Add(EtaTPCtglVsBG_ele);
    EtaQAList->Add(EtaTPCtglVsBG_pion);
    EtaQAList->Add(EtaTPCtglVsBG_proton);
    for (int i=0; i<3; i++) {
      EtaQAList->Add(EtaTPCtglVSBG_Nch[i]);
      EtaQAList->Add(EtaTPCtglVSBG_Nch_ele[i]);
      EtaQAList->Add(EtaTPCtglVSBG_Nch_pion[i]);
      EtaQAList->Add(EtaTPCtglVSBG_Nch_proton[i]);
      
      for (int j=0; j<NCentBins; j++) {
	EtaQAList->Add(EtaTPCtglVSBG[i][j]);
	EtaQAList->Add(EtaTPCtglVSBG_ele[i][j]);
	EtaQAList->Add(EtaTPCtglVSBG_pion[i][j]);
	EtaQAList->Add(EtaTPCtglVSBG_proton[i][j]);
      }
    }
    /*
      FitSlicesY(EtaTPCtglVsBG, meanNor, MPVEtaTPCtglVsBG, meanWidth, meanRes, meanChi, 100, 20, dummyList);
      FitSlicesY(EtaTPCtglVsBG_ele, meanNor, MPVEtaTPCtglVsBG_ele, meanWidth, meanRes, meanChi, 100, 20, dummyList);
      FitSlicesY(EtaTPCtglVsBG_pion, meanNor, MPVEtaTPCtglVsBG_pion, meanWidth, meanRes, meanChi, 100, 10, dummyList);
      FitSlicesY(EtaTPCtglVsBG_proton, meanNor, MPVEtaTPCtglVsBG_proton, meanWidth, meanRes, meanChi, 100, 20, dummyList);
    */
    SliceAverage(EtaTPCtglVsBG, MPVEtaTPCtglVsBG, dummyList);
    SliceAverage(EtaTPCtglVsBG_ele, MPVEtaTPCtglVsBG_ele, dummyList);
    SliceAverage(EtaTPCtglVsBG_pion, MPVEtaTPCtglVsBG_pion, dummyList);
    SliceAverage(EtaTPCtglVsBG_proton, MPVEtaTPCtglVsBG_proton, dummyList);
    EtaQAList->Add(MPVEtaTPCtglVsBG);
    EtaQAList->Add(MPVEtaTPCtglVsBG_ele);
    EtaQAList->Add(MPVEtaTPCtglVsBG_pion);
    EtaQAList->Add(MPVEtaTPCtglVsBG_proton);
    for (int i=0; i<3; i++) {
      /*FitSlicesY(EtaTPCtglVSBG_Nch[i], meanNor, MPVEtaTPCtglVSBG_Nch[i], meanWidth, meanRes, meanChi, 10,20, dummyList);
	FitSlicesY(EtaTPCtglVSBG_Nch_ele[i], meanNor, MPVEtaTPCtglVSBG_Nch_ele[i], meanWidth, meanRes, meanChi, 10, 20, dummyList);
	FitSlicesY(EtaTPCtglVSBG_Nch_pion[i], meanNor, MPVEtaTPCtglVSBG_Nch_pion[i], meanWidth, meanRes, meanChi, 10, 20, dummyList);
	FitSlicesY(EtaTPCtglVSBG_Nch_proton[i], meanNor, MPVEtaTPCtglVSBG_Nch_proton[i], meanWidth, meanRes, meanChi, 10, 20, dummyList);*/
      SliceAverage(EtaTPCtglVSBG_Nch[i], MPVEtaTPCtglVSBG_Nch[i], dummyList);
      SliceAverage(EtaTPCtglVSBG_Nch_ele[i], MPVEtaTPCtglVSBG_Nch_ele[i], dummyList);
      SliceAverage(EtaTPCtglVSBG_Nch_pion[i], MPVEtaTPCtglVSBG_Nch_pion[i], dummyList);
      SliceAverage(EtaTPCtglVSBG_Nch_proton[i], MPVEtaTPCtglVSBG_Nch_proton[i], dummyList);
      EtaQAList->Add(MPVEtaTPCtglVSBG_Nch[i]);
      EtaQAList->Add(MPVEtaTPCtglVSBG_Nch_ele[i]);
      EtaQAList->Add(MPVEtaTPCtglVSBG_Nch_pion[i]);
      EtaQAList->Add(MPVEtaTPCtglVSBG_Nch_proton[i]);
      
      for (int j=0; j<NCentBins; j++) {
	/*FitSlicesY(EtaTPCtglVSBG[i][j], meanNor, MPVEtaTPCtglVSBG[i][j], meanWidth, meanRes, meanChi, 10,20, dummyList);
	  FitSlicesY(EtaTPCtglVSBG_ele[i][j], meanNor, MPVEtaTPCtglVSBG_ele[i][j], meanWidth, meanRes, meanChi, 10, 20, dummyList);
	  FitSlicesY(EtaTPCtglVSBG_pion[i][j], meanNor, MPVEtaTPCtglVSBG_pion[i][j], meanWidth, meanRes, meanChi, 10, 20, dummyList);
	  FitSlicesY(EtaTPCtglVSBG_proton[i][j], meanNor, MPVEtaTPCtglVSBG_proton[i][j], meanWidth, meanRes, meanChi, 10, 20, dummyList);*/
	SliceAverage(EtaTPCtglVSBG[i][j], MPVEtaTPCtglVSBG[i][j], dummyList);
	SliceAverage(EtaTPCtglVSBG_ele[i][j], MPVEtaTPCtglVSBG_ele[i][j], dummyList);
	SliceAverage(EtaTPCtglVSBG_pion[i][j], MPVEtaTPCtglVSBG_pion[i][j], dummyList);
	SliceAverage(EtaTPCtglVSBG_proton[i][j], MPVEtaTPCtglVSBG_proton[i][j], dummyList);
	EtaQAList->Add(MPVEtaTPCtglVSBG[i][j]);
	EtaQAList->Add(MPVEtaTPCtglVSBG_ele[i][j]);
	EtaQAList->Add(MPVEtaTPCtglVSBG_pion[i][j]);
	EtaQAList->Add(MPVEtaTPCtglVSBG_proton[i][j]);
      }
    }
  }



  for (int i=0; i<3; i++) {
    meanList->Add(meanSigVSBG_Nch[i]);
    meanList->Add(meanSigVSBG_Nch_ele[i]);
    meanList->Add(meanSigVSBG_Nch_pion[i]);
    meanList->Add(meanSigVSBG_Nch_proton[i]);

    for (int j=0; j<NCentBins; j++) {
      meanList->Add(meanSigVSBG[i][j]);
      meanList->Add(meanSigVSBG_ele[i][j]);
      meanList->Add(meanSigVSBG_pion[i][j]);
      meanList->Add(meanSigVSBG_proton[i][j]);
    }
    for (int j=0; j<NEtaBins; j++) {
      meanList->Add(meanSigVSBGEta[i][j]);
      meanList->Add(meanSigVSBGEta_ele[i][j]);
      meanList->Add(meanSigVSBGEta_pion[i][j]);
      meanList->Add(meanSigVSBGEta_proton[i][j]);
    }
  }
  
  for (int i=0; i<3; i++) {
    FitSlicesY(meanSigVSBG_Nch[i], meanNor, MPVmeanSigVSBG_Nch[i], meanWidth, meanRes, meanChi, 10,20, dummyList);
    FitSlicesY(meanSigVSBG_Nch_ele[i], meanNor, MPVmeanSigVSBG_Nch_ele[i], meanWidth, meanRes, meanChi, 10, 20, dummyList);
    FitSlicesY(meanSigVSBG_Nch_pion[i], meanNor, MPVmeanSigVSBG_Nch_pion[i], meanWidth, meanRes, meanChi, 10, 20, dummyList);
    FitSlicesY(meanSigVSBG_Nch_proton[i], meanNor, MPVmeanSigVSBG_Nch_proton[i], meanWidth, meanRes, meanChi, 10, 20, dummyList);
    meanList->Add(MPVmeanSigVSBG_Nch[i]);
    meanList->Add(MPVmeanSigVSBG_Nch_ele[i]);
    meanList->Add(MPVmeanSigVSBG_Nch_pion[i]);
    meanList->Add(MPVmeanSigVSBG_Nch_proton[i]);
     
    for (int j=0; j<NCentBins; j++) {
      FitSlicesY(meanSigVSBG[i][j], meanNor, MPVmeanSigVSBG[i][j], meanWidth, meanRes, meanChi, 10,20, dummyList);
      FitSlicesY(meanSigVSBG_ele[i][j], meanNor, MPVmeanSigVSBG_ele[i][j], meanWidth, meanRes, meanChi, 10, 20, dummyList);
      FitSlicesY(meanSigVSBG_pion[i][j], meanNor, MPVmeanSigVSBG_pion[i][j], meanWidth, meanRes, meanChi, 10, 20, dummyList);
      FitSlicesY(meanSigVSBG_proton[i][j], meanNor, MPVmeanSigVSBG_proton[i][j], meanWidth, meanRes, meanChi, 10, 20, dummyList);
      meanList->Add(MPVmeanSigVSBG[i][j]);
      meanList->Add(MPVmeanSigVSBG_ele[i][j]);
      meanList->Add(MPVmeanSigVSBG_pion[i][j]);
      meanList->Add(MPVmeanSigVSBG_proton[i][j]);
    }

    for (int j=0; j<NEtaBins; j++) {
      FitSlicesY(meanSigVSBGEta[i][j], meanNor, MPVmeanSigVSBGEta[i][j], meanWidth, meanRes, meanChi, 10,20, dummyList);
      FitSlicesY(meanSigVSBGEta_ele[i][j], meanNor, MPVmeanSigVSBGEta_ele[i][j], meanWidth, meanRes, meanChi, 10, 20, dummyList);
      FitSlicesY(meanSigVSBGEta_pion[i][j], meanNor, MPVmeanSigVSBGEta_pion[i][j], meanWidth, meanRes, meanChi, 10, 20, dummyList);
      FitSlicesY(meanSigVSBGEta_proton[i][j], meanNor, MPVmeanSigVSBGEta_proton[i][j], meanWidth, meanRes, meanChi, 10, 20, dummyList);
      meanList->Add(MPVmeanSigVSBGEta[i][j]);
      meanList->Add(MPVmeanSigVSBGEta_ele[i][j]);
      meanList->Add(MPVmeanSigVSBGEta_pion[i][j]);
      meanList->Add(MPVmeanSigVSBGEta_proton[i][j]);
    }
  }

  // add signal vs beta*gamma histogram to lists
  meanList->Add(meanSigVsBG_plus_ele);
  // produce resulting histograms from signal vs beta*gamma for each beta*gamma slice
  FitSlicesY(meanSigVsBG_plus_ele, meanNor, meanMPV_plus_ele, meanWidth, meanRes, meanChi, minNumberOfEntries, maxFitSliceFitRMS, meanSlicesList);
  meanList->Add(meanMPV_plus_ele);
  
  // add signal vs beta*gamma histogram to lists
  meanList->Add(meanSigVsBG_plus_proton);
  // produce resulting histograms from signal vs beta*gamma for each beta*gamma slice
  FitSlicesY(meanSigVsBG_plus_proton, meanNor, meanMPV_plus_proton, meanWidth, meanRes, meanChi, minNumberOfEntries, maxFitSliceFitRMS, meanSlicesList);
  meanList->Add(meanMPV_plus_proton);
  
  // add signal vs beta*gamma histogram to lists
  meanList->Add(meanSigVsBG_plus_pion);
  // produce resulting histograms from signal vs beta*gamma for each beta*gamma slice
  FitSlicesY(meanSigVsBG_plus_pion, meanNor, meanMPV_plus_pion, meanWidth, meanRes, meanChi, minNumberOfEntries, maxFitSliceFitRMS, meanSlicesList);
  meanList->Add(meanMPV_plus_pion);

  // add signal vs beta*gamma histogram to lists
  meanList->Add(meanSigVsBG_plus);
  // produce resulting histograms from signal vs beta*gamma for each beta*gamma slice
  FitSlicesY(meanSigVsBG_plus, meanNor, meanMPV_plus, meanWidth, meanRes, meanChi, minNumberOfEntries, maxFitSliceFitRMS, meanSlicesList);
  // add resulting histograms to list
  meanList->Add(meanMPV_plus);
  
  // add signal vs beta*gamma histogram to lists
  meanList->Add(meanSigVsBG_minus_ele);
  // produce resulting histograms from signal vs beta*gamma for each beta*gamma slice
  FitSlicesY(meanSigVsBG_minus_ele, meanNor, meanMPV_minus_ele, meanWidth, meanRes, meanChi, minNumberOfEntries, maxFitSliceFitRMS, meanSlicesList);
  meanList->Add(meanMPV_minus_ele);
  
  // add signal vs beta*gamma histogram to lists
  meanList->Add(meanSigVsBG_minus_proton);
  // produce resulting histograms from signal vs beta*gamma for each beta*gamma slice
  FitSlicesY(meanSigVsBG_minus_proton, meanNor, meanMPV_minus_proton, meanWidth, meanRes, meanChi, minNumberOfEntries, maxFitSliceFitRMS, meanSlicesList);
  meanList->Add(meanMPV_minus_proton);
  
  // add signal vs beta*gamma histogram to lists
  meanList->Add(meanSigVsBG_minus_pion);
  // produce resulting histograms from signal vs beta*gamma for each beta*gamma slice
  FitSlicesY(meanSigVsBG_minus_pion, meanNor, meanMPV_minus_pion, meanWidth, meanRes, meanChi, minNumberOfEntries, maxFitSliceFitRMS, meanSlicesList);
  meanList->Add(meanMPV_minus_pion);
  
  // add signal vs beta*gamma histogram to lists
  meanList->Add(meanSigVsBG_minus);
  // produce resulting histograms from signal vs beta*gamma for each beta*gamma slice
  FitSlicesY(meanSigVsBG_minus, meanNor, meanMPV_minus, meanWidth, meanRes, meanChi, minNumberOfEntries, maxFitSliceFitRMS, meanSlicesList);
  // add resulting histograms to list
  meanList->Add(meanMPV_minus);

  // add signal vs beta*gamma histogram to lists
  meanList->Add(meanSigVsBG_ele);
  // produce resulting histograms from signal vs beta*gamma for each beta*gamma slice
  FitSlicesY(meanSigVsBG_ele, meanNor, meanMPV_ele, meanWidth, meanRes, meanChi, minNumberOfEntries, maxFitSliceFitRMS, meanSlicesList);
  meanList->Add(meanMPV_ele);
  
  // add signal vs beta*gamma histogram to lists
  meanList->Add(meanSigVsBG_proton);
  // produce resulting histograms from signal vs beta*gamma for each beta*gamma slice
  FitSlicesY(meanSigVsBG_proton, meanNor, meanMPV_proton, meanWidth, meanRes, meanChi, minNumberOfEntries, maxFitSliceFitRMS, meanSlicesList);
  meanList->Add(meanMPV_proton);
  
  // add signal vs beta*gamma histogram to lists
  meanList->Add(meanSigVsBG_pion);
  // produce resulting histograms from signal vs beta*gamma for each beta*gamma slice
  FitSlicesY(meanSigVsBG_pion, meanNor, meanMPV_pion, meanWidth, meanRes, meanChi, minNumberOfEntries, maxFitSliceFitRMS, meanSlicesList);
  meanList->Add(meanMPV_pion);
  
  // add signal vs beta*gamma histogram to lists
  meanList->Add(meanSigVsBG);
    // produce resulting histograms from signal vs beta*gamma for each beta*gamma slice
  FitSlicesY(meanSigVsBG, meanNor, meanMPV, meanWidth, meanRes, meanChi, minNumberOfEntries, maxFitSliceFitRMS, meanSlicesList);
  
  // add resulting histograms to list
  meanList->Add(meanNor);
  meanList->Add(meanMPV);
  meanList->Add(meanWidth);
  meanList->Add(meanRes);
  meanList->Add(meanChi);
    
  // produce signals vs beta*gamma histogram normalized to content of slices in y
  meanSigVsBGnorm = NormalHist(meanSigVsBG, 50, 1);
  
  // add normalized histogram to list
  meanList->Add(meanSigVsBGnorm);

  // start parameters for dE/dx + TR fits -- small dependence of results on initial parametrization!!!

  // Double_t meanPars[]={ 0.3, 1.622e+00 , 7.646e+00 , 4.214e-01 , 2.453e+00 , 1.748e-01 , 2.000e+00 , 4.679e-01};
  
  Double_t meanPars[] = {6.5e-01 , 2.8e+00 , 7.5e+00 , 2.7e+01 , 1.8e+00 , 6.5e-01 , 1.10e-00 , 1.60e-01 };
  if (inputSpec=="/LHC13bcNew" || inputSpec =="/LHC15nCorrect") {
    //  Double_t meanPars2[] = {9.507e-01 , 1.575e+00 , 7.434e+00 , 4.380e+01 , 1.636e+00 , 8.624e-01 , 5.300e-02 , 3.190e-03 }; // 2018
    // Double_t meanPars2[] = {0.72, 1.33, 7.7, 0.313, 3.05, 0.035, 2.17, 0.66};
    Double_t meanPars2[] = {0.5, 1, 8, 0.5, 3, 0.1, 2, 0.5};
    for (Int_t Par=0; Par<8; Par++) {
      meanPars[Par] = meanPars2[Par];
    }
  }
  // Double_t meanPars[] = { 1., 3.0, 0.1, 2., 0.5}
  // Double_t meanPars[] = { 0.313, 3.05, 0.035, 2.17, 0.66};
  // Double_t meanPars[] = {5.446e-01, 2.297e+00, 7.654e+00,3.200e-01, 2.875e+00,  3.166e-02	,2.314e+00	, 8.559e-01	};
  // Double_t meanPars[]= {2.011e+00 , 7.840e-01 , 8.174e+00 , 1.196e+00 , 2.078e+00 , 4.837e-01 , 7.089e-01 , 2.346e-08 };
  // Double_t meanPars[] = {0.55, 2.1, 7.9, 1.2, 1.9, 0.54, 1.6, 0.49};
  // Double_t meanPars[] = {4.422e-01, 2.544e+00, 7.575e+00, 3.069e-01,2.859e+00, 2.667e-02, 2.388e+00	, 8.688e-01  };
  // Double_t meanPars[] = {0.6607, 1.518, 7.832, 0.3160, 3.013, 0.03837, 2.214, 0.8232};
  // Double_t meanPars[] = {6.607e-01, 1.518e+00, 7.832e+00, 3.160e-01, 3.013e+00, 3.837e-02, 2.214e+00, 8.232e-01};
  Double_t meanErrs[100], meanChis[100];

  // produce dE/dx + TR fit of MPV distribution of fitted slices
  const Int_t kfail = ChisquareFit( meanMPV, MeandEdxTR, 8, meanPars, meanErrs, meanChis, 0, kTRUE); // meanMPV without threshold!!
  
  // get fit histogram and TF1 object
  TH1D *  meanFitHist = GetHfit("meanFitHist", MeandEdxTR, meanPars,  0.3, 1e4, kTRUE);
  TF1 *   meanFitFunc = new TF1("meanFitFunc", MeandEdxTR,            0.3, 1e4, 8);
    
  // adjust TF1 object
  meanFitFunc->SetParameters(meanPars);
  meanFitFunc->SetChisquare(meanChis[0]);
  meanFitFunc->SetNDF(meanChis[1]);
  
  // set meanFitHist title
  if (kfail) {      
    meanFitHist->SetTitle(Form("Fit fail %d!\n", kfail));
  }
  else {
    meanFitHist->SetTitle(Form("%.3e , %.3e , %.3e , %.3e , %.3e , %.3e , %.3e , %.3e , (%.2f/%d) %s %s", meanFitFunc->GetParameter(0), meanFitFunc->GetParameter(1), meanFitFunc->GetParameter(2), meanFitFunc->GetParameter(3), meanFitFunc->GetParameter(4), meanFitFunc->GetParameter(5), meanFitFunc->GetParameter(6), meanFitFunc->GetParameter(7), meanFitFunc->GetChisquare(), meanFitFunc->GetNDF(), inputFile.Data(), meanMPV->GetTitle()));
  }
  
  // add fit histogram and function to list
  meanList->Add(meanFitHist);
  meanList->Add(meanFitFunc);
  
  QAList->Add(TPCEta);
  QAList->Add(TOFEta);
  QAList->Add(TPCEta_ele);
  QAList->Add(TOFEta_ele);
  QAList->Add(TPCEta_pion);
  QAList->Add(TOFEta_pion);
  QAList->Add(TPCEta_proton);
  QAList->Add(TOFEta_proton);

  TRDCentrality->FitSlicesY(0,0,-1,10);
  TH1D* TRDCentrality_1 = (TH1D*)gDirectory->Get("TRDCentrality_1");
  QAList->Add(TRDCentrality);
  QAList->Add(TRDCentrality_1);
  QAList->Add(TPCCentrality);
  QAList->Add(TOFCentrality);
    

  /*
  // control GaussFrac dependence - creates Testplots for dependence of results/fit on gaussian fraction used for fitting the slices
  GaussFrac = 0.0; 
  TH1D * meanNorC, * meanMPVC, * meanWidthC, * meanResC, * meanChiC;
  TList * meanSlicesListC = new TList;
  FitSlicesY(meanSigVsBG, meanNorC, meanMPVC, meanWidthC, meanResC, meanChiC, minNumberOfEntries, maxFitSliceFitRMS, meanSlicesListC);
  
  // add resulting histograms to list
  meanList->Add(meanNorC);
  meanList->Add(meanMPVC);
  meanList->Add(meanWidthC);
  meanList->Add(meanResC);
  meanList->Add(meanChiC);
  
  Double_t meanParsC[] = {0.72, 1.33, 7.7, 0.313, 3.05, 0.035, 2.17, 0.66};
  Double_t meanErrsC[100], meanChisC[100];
  
  const Int_t kfailC = ChisquareFit( meanMPVC, MeandEdxTR, 8, meanParsC, meanErrsC, meanChisC, 0, kTRUE); // meanMPV
  
  TH1D *  meanFitHistC = GetHfit("meanFitHistC", MeandEdxTR, meanParsC,  0.3, 1e4, kTRUE);
  TF1 *   meanFitFuncC = new TF1("meanFitFuncC", MeandEdxTR,            0.3, 1e4, 8);
  
  meanFitFuncC->SetParameters(meanParsC);
  meanFitFuncC->SetChisquare(meanChisC[0]);
  meanFitFuncC->SetNDF(meanChisC[1]);
  
  if (kfailC) {      
  meanFitHistC->SetTitle(Form("Fit fail %d!\n", kfail));
  }
  else {
  meanFitHistC->SetTitle(Form("%.3e , %.3e , %.3e , %.3e , %.3e , %.3e , %.3e , %.3e , (%.2f/%d) %s %s", meanFitFuncC->GetParameter(0), meanFitFuncC->GetParameter(1), meanFitFuncC->GetParameter(2), meanFitFuncC->GetParameter(3), meanFitFuncC->GetParameter(4), meanFitFuncC->GetParameter(5), meanFitFuncC->GetParameter(6), meanFitFuncC->GetParameter(7), meanFitFuncC->GetChisquare(), meanFitFuncC->GetNDF(), inputFile.Data(), meanMPVC->GetTitle()));
    }
    
    // add fit histogram and function to list
    meanList->Add(meanFitHistC);
    meanList->Add(meanFitFuncC);
    
    TH1D* RatioFit = new TH1D("RatioFit Diff (0.00-0.50)/0.00", "RatioFit Diff (0.00-0.50)/0.00", 1000, 0.3, 1e4);
    TH1D* RatioMPV = new TH1D("RatioMPV Diff (0.00-0.50)/0.00", "RatioMPV Diff (0.00-0.50)/0.00", 100, 0.3, 1e4);
    BinLogX(RatioMPV->GetXaxis());
    BinLogX(RatioFit->GetXaxis()); 

    RatioFit->Add(meanFitHistC, meanFitHist, 1, -1);
    RatioFit->Divide(meanFitHist);
    RatioMPV->Add(meanMPVC, meanMPV, 1, -1);
    RatioMPV->Divide(meanMPV);

    TCanvas * ControlGaussFrac = new TCanvas("ControlGaussFrac", "ControlGaussFrac", 800, 800);
    ControlGaussFrac->Divide(2,2);
    ControlGaussFrac->cd(1);
    gStyle->SetOptStat(0);
    ControlGaussFrac->cd(1)->SetLogx();
    meanMPVC->SetLineColor(1);
    meanMPVC->Draw("");
    meanMPVC->SetMaximum(2.5);
    meanMPVC->SetMinimum(0.8);
    meanMPV->SetLineColor(2);
    meanMPV->Draw("same");

    ControlGaussFrac->cd(2);
    ControlGaussFrac->cd(2)->SetLogx();
    meanFitHistC->SetLineColor(1);
    meanFitHistC->Draw("");
    meanFitHistC->SetMaximum(2.5);
    meanFitHistC->SetMinimum(0.8);
    meanFitHist->SetLineColor(2);
    meanFitHist->Draw("same");

    ControlGaussFrac->cd(3);
    ControlGaussFrac->cd(3)->SetLogx();
    RatioMPV->Draw("");
    RatioMPV->SetMaximum(0.12);
    RatioMPV->SetMinimum(-0.12);

    ControlGaussFrac->cd(4);
    ControlGaussFrac->cd(4)->SetLogx();
    RatioFit->Draw("");
    RatioFit->SetMaximum(0.12);
    RatioFit->SetMinimum(-0.12);

    ControlGaussFrac->SaveAs(Form("output/ControlGaussFrac_%i_%s%s%s.pdf",FileSpec.Data(), EtaSpec.Data(), NclsCorr.Data()));

    meanList->Add(RatioMPV);
    meanList->Add(RatioFit);
    */


  //===========================================================================================
  //
  // Save output
  //
  //===========================================================================================
  
  // ctime for time + date output to parameter .txt file
  std::time_t time = std::time(NULL);
  
  // histograms
  TFile meanFile(Form("output/%s_meanParam_pid%i_%s%s%s%s.root", inputSpec.Data(), pidOpt,  FileSpec.Data(),CentCorr.Data(), EtaSpec.Data(), NclsCorr.Data()), "RECREATE");
  TDirectory *Basic = meanFile.mkdir("Basic");
  Basic->cd();    // make the "tof" directory the current directory
  meanList->Write();


  TDirectory *QAFolder = meanFile.mkdir("DetectorQA");
  QAFolder->cd();    // make the "tof" directory the current directory
  QAList->Write();

  TDirectory *ClusterQAFolder = meanFile.mkdir("ClusterQA");
  ClusterQAFolder->cd();    // make the "tof" directory the current directory
  ClusterQAList->Write();

  TDirectory *EtaQAFolder = meanFile.mkdir("EtaQA");
  EtaQAFolder->cd();    // make the "tof" directory the current directory
  EtaQAList->Write();
  
  Basic->cd();
    
  TCanvas *DataAndFitResults = new TCanvas("DataAndFitResults", "DataAndFitResults", 800, 800);
  DataAndFitResults->SetLogx();
  meanFitHist->SetLineColor(1);
  meanFitHist->Draw();
  // meanFitFunc->Draw("same");
  meanMPV->SetLineColor(2);
  meanMPV->Draw("same");
  DataAndFitResults->Write();
  DataAndFitResults->SaveAs(Form("output/FitResults_%s%s%s%s.pdf", FileSpec.Data(), CentCorr.Data(), EtaSpec.Data(), NclsCorr.Data()));
  DataAndFitResults->Close();
  meanSlicesList->Write("slices", TObject::kSingleKey);
  meanFile.Close();

  // parameters
  FILE* meanParamFile = fopen(Form("output/%s_meanParam_pid%i_%s%s%s%s.txt", inputSpec.Data(),pidOpt,  FileSpec.Data(), CentCorr.Data(), EtaSpec.Data(), NclsCorr.Data()), "a");
  cout << std::ctime(&time) << endl;
  fprintf(meanParamFile, "%s", std::ctime(&time));
  fprintf(meanParamFile, "chi^2 = %.1f\nndf = %.1f\nchi^2/ndf = %.1f\n\n", meanChis[0], meanChis[1], meanChis[0]/meanChis[1]);
  for (Int_t i = 0; i < 8; i++) {
    fprintf(meanParamFile, "par[%i] = %.3e\t\terr[%i] = %.3e\n", i, meanPars[i], i, meanErrs[i]);
  }
  fclose(meanParamFile);
  

  delete meanList;
  delete meanSlicesList;
  delete meanSigVsBG;
  delete clusterVsChamber;
  delete clusterVsChamberWCut;
  delete meanFitFunc;
}
