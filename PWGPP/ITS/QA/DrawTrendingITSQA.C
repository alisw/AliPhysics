#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TPad.h>
#include <TGraphErrors.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <Riostream.h>
#include <iostream>
#include <fstream>
#endif


/*
  
  Macro to draw the ITS QA trending plots by accessing the std tree.
  To be mainly used with the automatic scripts to fill the QA repository.
  Launch with 
  aliroot -x -b -q "DrawTrendingITSQA.C" 
  The macro produces one png file for each trending variables
  and a .root file with the histograms
*/


Int_t DrawTrendingITSQA(TString mergedTrendFile = "trending.root", // trending tree file 
			Bool_t displayAll = kFALSE) //set to kTRUE to display trending for expert plots
{
  //
  //reads merged trending.root file and draws trending plots from tree
  //
  if (!mergedTrendFile) {
    Printf("Cannot open merged trend file with ITS QA");
    return 1;
  }
  
  char  outfilename[200]= "ProductionQA.hist.root";
  TString plotDir(".");
 
  TFile * fin = TFile::Open(mergedTrendFile.Data());
  TTree * ttree = (TTree*) fin->Get("trending");
  if (!ttree){
    Printf("Invalid trending tree.");
    return 2;
  }


  //************************************************ SDD ************************************************//
  
  Int_t nrun,nEvents, nEventsTriggered;
  Float_t minDrTime,errminDrTime;
  Float_t  meanDrTime,errmeanDrTime,EmptyModulesSDD;
  Float_t fracTrackWithClu1,fracTrackWithClu2,errfracTrackWithClu1,errfracTrackWithClu2;
  Float_t fracTrackWithClu3,fracTrackWithClu4,errfracTrackWithClu3,errfracTrackWithClu4;
  Float_t fracTrackWithClu5,fracTrackWithClu6,errfracTrackWithClu5,errfracTrackWithClu6;
  Float_t MPVdEdxLay3,errMPVdEdxLay3,MPVdEdxLay4,errMPVdEdxLay4;
  Float_t MPVdEdxTB0,errMPVdEdxTB0,MPVdEdxTB5,errMPVdEdxTB5;
  Float_t fracExtra,errfracExtra;
  
  ttree->SetBranchAddress("nrun",&nrun);
  ttree->SetBranchAddress("nEvents",&nEvents);
  ttree->SetBranchAddress("nEventsTriggered",&nEventsTriggered);
  ttree->SetBranchAddress("minDrTime",&minDrTime);
  ttree->SetBranchAddress("errminDrTime",&errminDrTime);
  ttree->SetBranchAddress("meanDrTime",&meanDrTime);
  ttree->SetBranchAddress("errmeanDrTime",&errmeanDrTime); //mean time
  ttree->SetBranchAddress("fracTrackWithClu1",&fracTrackWithClu1); //fraction of tracks with cluster in layer 1
  ttree->SetBranchAddress("errfracTrackWithClu1",&errfracTrackWithClu1); //error fraction of tracks with cluster in layer 1
  ttree->SetBranchAddress("fracTrackWithClu2",&fracTrackWithClu2); //fraction of tracks with cluster in layer 2
  ttree->SetBranchAddress("errfracTrackWithClu2",&errfracTrackWithClu2); //error fraction of tracks with cluster in layer 2
  ttree->SetBranchAddress("fracTrackWithClu3",&fracTrackWithClu3);//fraction of tracks with cluster in layer 3
  ttree->SetBranchAddress("errfracTrackWithClu3",&errfracTrackWithClu3);//error fraction of tracks with cluster in layer 3
  ttree->SetBranchAddress("fracTrackWithClu4",&fracTrackWithClu4);//fraction of tracks with cluster in layer 4
  ttree->SetBranchAddress("errfracTrackWithClu4",&errfracTrackWithClu4);//error fraction of tracks with cluster in layer 4
  ttree->SetBranchAddress("fracTrackWithClu5",&fracTrackWithClu5);//fraction of tracks with cluster in layer 5
  ttree->SetBranchAddress("errfracTrackWithClu5",&errfracTrackWithClu5);//error fraction of tracks with cluster in layer 5
  ttree->SetBranchAddress("fracTrackWithClu6",&fracTrackWithClu6);//fraction of tracks with cluster in layer 6
  ttree->SetBranchAddress("errfracTrackWithClu6",&errfracTrackWithClu6);//error fraction of tracks with cluster in layer 6
  ttree->SetBranchAddress("EmptyModulesSDD",&EmptyModulesSDD); // Number of empty SSD  modules
  ttree->SetBranchAddress("fracExtra",&fracExtra); // fraction of extra clusters in SDD
  ttree->SetBranchAddress("errfracExtra",&errfracExtra); // fraction of extra clusters in SDD    
  ttree->SetBranchAddress("MPVdEdxLay3",&MPVdEdxLay3); // most probable value of dE/dx distribution of SDD Layer 3
  ttree->SetBranchAddress("errMPVdEdxLay3",&errMPVdEdxLay3); // error  most probable value of dE/dx distribution of SDD Layer 3
  ttree->SetBranchAddress("MPVdEdxLay4",&MPVdEdxLay4); // most probable value of dE/dx distribution of SDD Layer 4
  ttree->SetBranchAddress("errMPVdEdxLay4",&errMPVdEdxLay4); // error  most probable value of dE/dx distribution of SDD Layer 4
  ttree->SetBranchAddress("MPVdEdxTB0",&MPVdEdxTB0); // most probable value of dE/dx distribution of SDD - small drift time
  ttree->SetBranchAddress("errMPVdEdxTB0",&errMPVdEdxTB0); // most probable value of dE/dx distribution of SDD - small drift time
  ttree->SetBranchAddress("MPVdEdxTB5",&MPVdEdxTB5); // most probable value of dE/dx distribution of SDD - large drift time
  ttree->SetBranchAddress("errMPVdEdxTB5",&errMPVdEdxTB5); // most probable value of dE/dx distribution of SDD - large drift time    


    
  
  //************************************************ VERTEX ************************************************//
    
  Float_t meanVtxTRKx,meanVtxTRKy,meanVtxTRKz;
  Float_t meanVtxSPDx,meanVtxSPDy,meanVtxSPDz;
  Float_t sigmaVtxTRKx,sigmaVtxTRKy,sigmaVtxTRKz;
  Float_t sigmaVtxSPDx,sigmaVtxSPDy,sigmaVtxSPDz;
  Float_t meanVtxTRKxErr,meanVtxTRKyErr,meanVtxTRKzErr;
  Float_t meanVtxSPDxErr,meanVtxSPDyErr,meanVtxSPDzErr;
  Float_t sigmaVtxTRKxErr,sigmaVtxTRKyErr,sigmaVtxTRKzErr;
  Float_t sigmaVtxSPDxErr,sigmaVtxSPDyErr,sigmaVtxSPDzErr;

  ttree->SetBranchAddress("meanVtxTRKx",&meanVtxTRKx); // mean of tracks vertex position - x
  ttree->SetBranchAddress("meanVtxTRKy",&meanVtxTRKy); // mean of tracks vertex position - y
  ttree->SetBranchAddress("meanVtxTRKz",&meanVtxTRKz); // mean of tracks vertex position - z
  ttree->SetBranchAddress("meanVtxTRKxErr",&meanVtxTRKxErr); // error mean of tracks vertex position - x
  ttree->SetBranchAddress("meanVtxTRKyErr",&meanVtxTRKyErr); // error mean of tracks vertex position - y
  ttree->SetBranchAddress("meanVtxTRKzErr",&meanVtxTRKzErr); // error mean of tracks vertex position - z
  ttree->SetBranchAddress("meanVtxSPDx",&meanVtxSPDx); // mean of SPD vertex position - x
  ttree->SetBranchAddress("meanVtxSPDy",&meanVtxSPDy); // mean of SPD vertex position - y
  ttree->SetBranchAddress("meanVtxSPDz",&meanVtxSPDz); // mean of SPD vertex position - z
  ttree->SetBranchAddress("meanVtxSPDxErr",&meanVtxSPDxErr); // error mean of SPD vertex position - x
  ttree->SetBranchAddress("meanVtxSPDyErr",&meanVtxSPDyErr); // error mean of SPD vertex position - y
  ttree->SetBranchAddress("meanVtxSPDzErr",&meanVtxSPDzErr); // error mean of SPD vertex position - z
  ttree->SetBranchAddress("sigmaVtxTRKx",&sigmaVtxTRKx); // sigma of tracks vertex position - x
  ttree->SetBranchAddress("sigmaVtxTRKy",&sigmaVtxTRKy); // sigma of tracks vertex position - y
  ttree->SetBranchAddress("sigmaVtxTRKz",&sigmaVtxTRKz); // sigma of tracks vertex position - z
  ttree->SetBranchAddress("sigmaVtxTRKxErr",&sigmaVtxTRKxErr); // error sigma of tracks vertex position - x
  ttree->SetBranchAddress("sigmaVtxTRKyErr",&sigmaVtxTRKyErr); // error sigma of tracks vertex position - y
  ttree->SetBranchAddress("sigmaVtxTRKzErr",&sigmaVtxTRKzErr); // error sigma of tracks vertex position - z
  ttree->SetBranchAddress("sigmaVtxSPDx",&sigmaVtxSPDx); // sigma of tracks vertex position - x
  ttree->SetBranchAddress("sigmaVtxSPDy",&sigmaVtxSPDy); // sigma of tracks vertex position - y
  ttree->SetBranchAddress("sigmaVtxSPDz",&sigmaVtxSPDz); // sigma of tracks vertex position - z
  ttree->SetBranchAddress("sigmaVtxSPDxErr",&sigmaVtxSPDxErr); // error sigma of tracks vertex position - x
  ttree->SetBranchAddress("sigmaVtxSPDyErr",&sigmaVtxSPDyErr); // error sigma of tracks vertex position - y
  ttree->SetBranchAddress("sigmaVtxSPDzErr",&sigmaVtxSPDzErr); // error sigma of tracks vertex position - z
 

  //************************************************ SSD **************************************************//
     
  Float_t MPVL5,MPVErrL5;
  Float_t MPVL6,MPVErrL6;
  Float_t ChargeRatioL5,ChargeRatioErrL5;
  Float_t ChargeRatioL6,ChargeRatioErrL6;
    

  //************************************************ MATCHING *********************************************//
  
  Float_t FracSPD1;
  Float_t errFracSPD1;
  Float_t FracSPD2;
  Float_t errFracSPD2;
  Float_t Eff6Pt02;
  Float_t errEff6Pt02;
  Float_t Eff6Pt1;
  Float_t errEff6Pt1;
  Float_t Eff6Pt10;
  Float_t errEff6Pt10;
  Float_t Eff5Pt02;
  Float_t errEff5Pt02;
  Float_t Eff5Pt1;
  Float_t errEff5Pt1;
  Float_t Eff5Pt10;
  Float_t errEff5Pt10;
  Float_t Eff4Pt02;
  Float_t errEff4Pt02;
  Float_t Eff4Pt1;
  Float_t errEff4Pt1;
  Float_t Eff4Pt10;
  Float_t errEff4Pt10;
  Float_t Eff3Pt02;
  Float_t errEff3Pt02;
  Float_t Eff3Pt1;
  Float_t errEff3Pt1;
  Float_t Eff3Pt10;
  Float_t errEff3Pt10;
  Float_t Eff2Pt02;
  Float_t errEff2Pt02;
  Float_t Eff2Pt1;
  Float_t errEff2Pt1;
  Float_t Eff2Pt10;
  Float_t errEff2Pt10;
  Float_t EffSPDPt02;
  Float_t errEffSPDPt02;
  Float_t EffSPDPt1;
  Float_t errEffSPDPt1;
  Float_t EffSPDPt10;
  Float_t errEffSPDPt10;
  Float_t EffoneSPDPt02;
  Float_t errEffoneSPDPt02;
  Float_t EffoneSPDPt1;
  Float_t errEffoneSPDPt1;
  Float_t EffoneSPDPt10;
  Float_t errEffoneSPDPt10;
  Float_t EffTOTPt02;
  Float_t errEffTOTPt02;
  Float_t EffTOTPt1;
  Float_t errEffTOTPt1;
  Float_t EffTOTPt10;
  Float_t errEffTOTPt10;



 
  ttree->SetBranchAddress("MPVL5",&MPVL5); // Most Probable Value dEdx Layer 5
  ttree->SetBranchAddress("MPVErrL5",&MPVErrL5); // Most Probable Value error dEdx Layer 5
  ttree->SetBranchAddress("MPVL6",&MPVL6); // Most Probable Value dEdx Layer 6
  ttree->SetBranchAddress("MPVErrL6",&MPVErrL6); // Most Probable Value error dEdx Layer 6
  ttree->SetBranchAddress("ChargeRatioL5",&ChargeRatioL5); // Charge ratio (2 sides of SSD) Layer 5
  ttree->SetBranchAddress("ChargeRatioErrL5",&ChargeRatioErrL5); // Charge ratio error (2 sides of SSD) Layer 5
  ttree->SetBranchAddress("ChargeRatioL6",&ChargeRatioL6); // Charge ratio (2 sides of SSD) Layer 6
  ttree->SetBranchAddress("ChargeRatioErrL6",&ChargeRatioErrL6); // Charge ratio error(2 sides of SSD) Layer 6


  ttree->SetBranchAddress("Eff6Pt02",&Eff6Pt02); // matching efficiency low pt 6 clusters
  ttree->SetBranchAddress("errEff6Pt02",&errEff6Pt02); // error matching efficiency low pt 6 clusters
  ttree->SetBranchAddress("Eff5Pt02",&Eff5Pt02); // matching efficiency low pt 5 clusters
  ttree->SetBranchAddress("errEff5Pt02",&errEff5Pt02); // error matching efficiency low pt 5 clusters
  ttree->SetBranchAddress("Eff4Pt02",&Eff4Pt02); // matching efficiency low pt 4 clusters
  ttree->SetBranchAddress("errEff4Pt02",&errEff4Pt02); // error matching efficiency low pt 4 clusters
  ttree->SetBranchAddress("Eff3Pt02",&Eff3Pt02); // matching efficiency low pt 3 clusters
  ttree->SetBranchAddress("errEff3Pt02",&errEff3Pt02); // error matching efficiency low pt 3 clusters
  ttree->SetBranchAddress("Eff2Pt02",&Eff2Pt02); // matching efficiency low pt 2 clusters
  ttree->SetBranchAddress("errEff2Pt02",&errEff2Pt02); // error matching efficiency low pt 2 clusters
  ttree->SetBranchAddress("EffSPDPt02",&EffSPDPt02); // matching efficiency low pt 2 SPD
  ttree->SetBranchAddress("errEffSPDPt02",&errEffSPDPt02); // error matching efficiency low pt 2 SPD
  ttree->SetBranchAddress("EffoneSPDPt02",&EffoneSPDPt02); // matching efficiency low pt 6 one SPD
  ttree->SetBranchAddress("errEffoneSPDPt02",&errEffoneSPDPt02); // error matching efficiency low pt one SPD
  ttree->SetBranchAddress("EffTOTPt02",&EffTOTPt02); // matching efficiency low pt
  ttree->SetBranchAddress("errEffTOTPt02",&errEffTOTPt02); // error matching efficiency low pt
  ttree->SetBranchAddress("Eff6Pt1",&Eff6Pt1); // matching efficiency mid pt 6 clusters
  ttree->SetBranchAddress("errEff6Pt1",&errEff6Pt1); // error matching efficiency mid pt 6 clusters
  ttree->SetBranchAddress("Eff5Pt1",&Eff5Pt1); // matching efficiency mid pt 5 clusters
  ttree->SetBranchAddress("errEff5Pt1",&errEff5Pt1); // error matching efficiency mid pt 5 clusters
  ttree->SetBranchAddress("Eff4Pt1",&Eff4Pt1); // matching efficiency mid pt 4 clusters
  ttree->SetBranchAddress("errEff4Pt1",&errEff4Pt1); // error matching efficiency mid pt 4 clusters
  ttree->SetBranchAddress("Eff3Pt1",&Eff3Pt1); // matching efficiency mid pt 3 clusters
  ttree->SetBranchAddress("errEff3Pt1",&errEff3Pt1); // error matching efficiency mid pt 3 clusters
  ttree->SetBranchAddress("Eff2Pt1",&Eff2Pt1); // matching efficiency mid pt 2 clusters
  ttree->SetBranchAddress("errEff2Pt1",&errEff2Pt1); // error matching efficiency mid pt 2 clusters
  ttree->SetBranchAddress("EffSPDPt1",&EffSPDPt1); // matching efficiency mid pt 2 SPD
  ttree->SetBranchAddress("errEffSPDPt1",&errEffSPDPt1); // error matching efficiency mid pt 2 SPD
  ttree->SetBranchAddress("EffoneSPDPt1",&EffoneSPDPt1); // matching efficiency mid pt 6 one SPD
  ttree->SetBranchAddress("errEffoneSPDPt1",&errEffoneSPDPt1); // error matching efficiency mid pt one SPD
  ttree->SetBranchAddress("EffTOTPt1",&EffTOTPt1); // matching efficiency mid pt
  ttree->SetBranchAddress("errEffTOTPt1",&errEffTOTPt1); // error matching efficiency mid pt
  ttree->SetBranchAddress("Eff6Pt10",&Eff6Pt10); // matching efficiency high pt 6 clusters
  ttree->SetBranchAddress("errEff6Pt10",&errEff6Pt10); // error matching efficiency high pt 6 clusters
  ttree->SetBranchAddress("Eff5Pt10",&Eff5Pt10); // matching efficiency high pt 5 clusters
  ttree->SetBranchAddress("errEff5Pt10",&errEff5Pt10); // error matching efficiency high pt 5 clusters
  ttree->SetBranchAddress("Eff4Pt10",&Eff4Pt10); // matching efficiency high pt 4 clusters
  ttree->SetBranchAddress("errEff4Pt10",&errEff4Pt10); // error matching efficiency high pt 4 clusters
  ttree->SetBranchAddress("Eff3Pt10",&Eff3Pt10); // matching efficiency high pt 3 clusters
  ttree->SetBranchAddress("errEff3Pt10",&errEff3Pt10); // error matching efficiency high pt 3 clusters
  ttree->SetBranchAddress("Eff2Pt10",&Eff2Pt10); // matching efficiency high pt 2 clusters
  ttree->SetBranchAddress("errEff2Pt10",&errEff2Pt10); // error matching efficiency high pt 2 clusters
  ttree->SetBranchAddress("EffSPDPt10",&EffSPDPt10); // matching efficiency high pt 2 SPD
  ttree->SetBranchAddress("errEffSPDPt10",&errEffSPDPt10); // error matching efficiency high pt 2 SPD
  ttree->SetBranchAddress("EffoneSPDPt10",&EffoneSPDPt10); // matching efficiency high pt 6 one SPD
  ttree->SetBranchAddress("errEffoneSPDPt10",&errEffoneSPDPt10); // error matching efficiency high pt one SPD
  ttree->SetBranchAddress("EffTOTPt10",&EffTOTPt10); // matching efficiency high pt
  ttree->SetBranchAddress("errEffTOTPt10",&errEffTOTPt10); // error matching efficiency high pt  
  ttree->SetBranchAddress("FracSPD1",&FracSPD1); // fraction SPD layers active on 1 layer
  ttree->SetBranchAddress("errFracSPD1",&errFracSPD1);
  ttree->SetBranchAddress("FracSPD2",&FracSPD2); // fraction SPD layers active on 1 layer
  ttree->SetBranchAddress("errFracSPD2",&errFracSPD2);
    
  Int_t nRuns=ttree->GetEntries();
  TList lista;


  //****************************** SDD  PLOTS *****************************************************//
  
  TH1F* histonEvents=new TH1F("histonEvents","",nRuns,0.,nRuns);
  TH1F* histonEventsTriggered=new TH1F("histoEventsTriggered","",nRuns,0.,nRuns);
  TH1F* histominTime=new TH1F("histominTime","",nRuns,0.,nRuns);
  TH1F* histomeanTime=new TH1F("histomeanTime","",nRuns,0.,nRuns);
  TH1F* histofracExtra=new TH1F("histofracExtra","",nRuns,0.,nRuns);
  TH1F* histodEdxTB0=new TH1F("histodEdxTB0","",nRuns,0.,nRuns);
  TH1F* histodEdxTB5=new TH1F("histodEdxTB5","",nRuns,0.,nRuns);
  TH1F* histodEdxLay3=new TH1F("histodEdxLay3","",nRuns,0.,nRuns);
  TH1F* histodEdxLay4=new TH1F("histodEdxLay4","",nRuns,0.,nRuns);
  TH1F* histoTrackClu1=new TH1F("histoTrackClu1","",nRuns,0.,nRuns);
  TH1F* histoTrackClu2=new TH1F("histoTrackClu2","",nRuns,0.,nRuns);
  TH1F* histoTrackClu3=new TH1F("histoTrackClu3","",nRuns,0.,nRuns);
  TH1F* histoTrackClu4=new TH1F("histoTrackClu4","",nRuns,0.,nRuns);
  TH1F* histoTrackClu5=new TH1F("histoTrackClu5","",nRuns,0.,nRuns);
  TH1F* histoTrackClu6=new TH1F("histoTrackClu6","",nRuns,0.,nRuns);
  TH1F* histoNmodEmpty=new TH1F("histoNmodEmpty","",nRuns,0.,nRuns);

  //****************************** VERTEX  PLOTS *****************************************************//

  TH1F *hMeanVx = new TH1F("hMeanVx","Track Vertex Vx Distribution",nRuns,0.,nRuns);
  TH1F *hMeanVy = new TH1F("hMeanVy","Track Vertex Vy Distribution",nRuns,0.,nRuns);
  TH1F *hMeanVz = new TH1F("hMeanVz","Track Vertex Vz Distribution",nRuns,0.,nRuns);
  TH1F *hSigmaVx = new TH1F("hSigmaVx","Track Vertex SigmaVx Distribution",nRuns,0.,nRuns);
  TH1F *hSigmaVy = new TH1F("hSigmaVy","Track Vertex SigmaVy Distribution",nRuns,0.,nRuns);
  TH1F *hSigmaVz = new TH1F("hSigmaVz","Track Vertex SigmaVz Distribution",nRuns,0.,nRuns);
  TH1F *hMeanVxSPD = new TH1F("hMeanVxSPD","Track Vertex Vx Distribution",nRuns,0.,nRuns);
  TH1F *hMeanVySPD = new TH1F("hMeanVySPD","Track Vertex Vy Distribution",nRuns,0.,nRuns);
  TH1F *hMeanVzSPD = new TH1F("hMeanVzSPD","Track Vertex Vz Distribution",nRuns,0.,nRuns);
  TH1F *hSigmaVxSPD = new TH1F("hSigmaVxSPD","Track Vertex SigmaVx Distribution",nRuns,0.,nRuns);
  TH1F *hSigmaVySPD = new TH1F("hSigmaVySPD","Track Vertex SigmaVy Distribution",nRuns,0.,nRuns);
  TH1F *hSigmaVzSPD = new TH1F("hSigmaVzSPD","Track Vertex SigmaVz Distribution",nRuns,0.,nRuns);
  

  //****************************** SSD  PLOTS *****************************************************//
  TH1F* histodEdxLay5 = new TH1F("histodEdxLay5","",nRuns,0.,nRuns);
  TH1F* histodEdxLay6 = new TH1F("histodEdxLay6","",nRuns,0.,nRuns);
  TH1F* histoChargeRatioLay5 = new TH1F("histoChargeRatioLay5","",nRuns,0.,nRuns);
  TH1F* histoChargeRatioLay6 = new TH1F("histoChargeRatioLay6","",nRuns,0.,nRuns);
    

 //************************************************ MATCHING PLOTS *********************************************//
  TH1F *hFracSPD1 = new TH1F("hFracSPD1","SPD inner; run number; Fraction of HSs",nRuns,0.,nRuns);
  hFracSPD1->SetLineColor(kGreen+2);
  hFracSPD1->SetMarkerColor(kGreen+2);
  hFracSPD1->SetMarkerStyle(20);
  TH1F *hFracSPD2 = new TH1F("hFracSPD2","SPD outer; run number; Fraction of HSs",nRuns,0.,nRuns);
  hFracSPD2->SetLineColor(kYellow+2);
  hFracSPD2->SetMarkerColor(kYellow+2);
  hFracSPD2->SetMarkerStyle(20);

  TH1F *hEffSPDPt02 = new TH1F("hEffSPDPt02","Efficiency; run number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEffSPDPt02->SetLineWidth(2);
  hEffSPDPt02->SetLineColor(kAzure+1);
  hEffSPDPt02->SetMarkerColor(kAzure+1);
  hEffSPDPt02->SetMarkerStyle(20);
  TH1F *hEffSPDPt1 = new TH1F("hEffSPDPt1","Efficiency; run number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEffSPDPt1->SetLineWidth(2);
  hEffSPDPt1->SetLineColor(kAzure+1);
  hEffSPDPt1->SetMarkerColor(kAzure+1);
  hEffSPDPt1->SetMarkerStyle(20);
  TH1F *hEffSPDPt10 = new TH1F("hEffSPDPt10","Efficiency; run number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEffSPDPt10->SetLineWidth(2);
  hEffSPDPt10->SetLineColor(kAzure+1);
  hEffSPDPt10->SetMarkerColor(kAzure+1);
  hEffSPDPt10->SetMarkerStyle(20);

  TH1F *hEffoneSPDPt02 = new TH1F("hEffoneSPDPt02","Efficiency; run number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEffoneSPDPt02->SetLineWidth(2);
  hEffoneSPDPt02->SetLineColor(kGray);
  hEffoneSPDPt02->SetMarkerColor(kGray);
  hEffoneSPDPt02->SetMarkerStyle(20);
  TH1F *hEffoneSPDPt1 = new TH1F("hEffoneSPDPt1","Efficiency; run number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEffoneSPDPt1->SetLineWidth(2);
  hEffoneSPDPt1->SetLineColor(kGray);
  hEffoneSPDPt1->SetMarkerColor(kGray);
  hEffoneSPDPt1->SetMarkerStyle(20);
  TH1F *hEffoneSPDPt10 = new TH1F("hEffoneSPDPt10","Efficiency; run number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEffoneSPDPt10->SetLineWidth(2);
  hEffoneSPDPt10->SetLineColor(kGray);
  hEffoneSPDPt10->SetMarkerColor(kGray);
  hEffoneSPDPt10->SetMarkerStyle(20);

  TH1F *hEff2Pt02 = new TH1F("hEff2Pt02","Efficiency; run number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEff2Pt02->SetLineWidth(2);
  hEff2Pt02->SetLineColor(kViolet);
  hEff2Pt02->SetMarkerColor(kViolet);
  hEff2Pt02->SetMarkerStyle(20);
  TH1F *hEff2Pt1 = new TH1F("hEff2Pt1","Efficiency; run number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEff2Pt1->SetLineWidth(2);
  hEff2Pt1->SetLineColor(kViolet);
  hEff2Pt1->SetMarkerColor(kViolet);
  hEff2Pt1->SetMarkerStyle(20);
  TH1F *hEff2Pt10 = new TH1F("hEff2Pt10","Efficiency; run number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEff2Pt10->SetLineWidth(2);
  hEff2Pt10->SetLineColor(kViolet);
  hEff2Pt10->SetMarkerColor(kViolet);
  hEff2Pt10->SetMarkerStyle(20);

  TH1F *hEff3Pt02 = new TH1F("hEff3Pt02","Efficiency; run number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEff3Pt02->SetLineWidth(2);
  hEff3Pt02->SetLineColor(6);
  hEff3Pt02->SetMarkerColor(6);
  hEff3Pt02->SetMarkerStyle(20);
  TH1F *hEff3Pt1 = new TH1F("hEff3Pt1","Efficiency; run number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEff3Pt1->SetLineWidth(2);
  hEff3Pt1->SetLineColor(6);
  hEff3Pt1->SetMarkerColor(6);
  hEff3Pt1->SetMarkerStyle(20);
  TH1F *hEff3Pt10 = new TH1F("hEff3Pt10","Efficiency; run number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEff3Pt10->SetLineWidth(2);
  hEff3Pt10->SetLineColor(6);
  hEff3Pt10->SetMarkerColor(6);
  hEff3Pt10->SetMarkerStyle(20);

  TH1F *hEff4Pt02 = new TH1F("hEff4Pt02","Efficiency; run number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEff4Pt02->SetLineWidth(2);
  hEff4Pt02->SetLineColor(4);
  hEff4Pt02->SetMarkerColor(4);
  hEff4Pt02->SetMarkerStyle(20);
  TH1F *hEff4Pt1 = new TH1F("hEff4Pt1","Efficiency; run number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEff4Pt1->SetLineWidth(2);
  hEff4Pt1->SetLineColor(4);
  hEff4Pt1->SetMarkerColor(4);
  hEff4Pt1->SetMarkerStyle(20);
  TH1F *hEff4Pt10 = new TH1F("hEff4Pt10","Efficiency; run number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEff4Pt10->SetLineWidth(2);
  hEff4Pt10->SetLineColor(4);
  hEff4Pt10->SetMarkerColor(4);
  hEff4Pt10->SetMarkerStyle(20);

  TH1F *hEff5Pt02 = new TH1F("hEff5Pt02","Efficiency; run number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEff5Pt02->SetLineWidth(2);
  hEff5Pt02->SetLineColor(3);
  hEff5Pt02->SetMarkerColor(3);
  hEff5Pt02->SetMarkerStyle(20);
  TH1F *hEff5Pt1 = new TH1F("hEff5Pt1","Efficiency; run number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEff5Pt1->SetLineWidth(2);
  hEff5Pt1->SetLineColor(3);
  hEff5Pt1->SetMarkerColor(3);
  hEff5Pt1->SetMarkerStyle(20);
  TH1F *hEff5Pt10 = new TH1F("hEff5Pt10","Efficiency; run number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEff5Pt10->SetLineWidth(3);
  hEff5Pt10->SetLineColor(3);
  hEff5Pt10->SetMarkerColor(3);
  hEff5Pt10->SetMarkerStyle(20);

  TH1F *hEff6Pt02 = new TH1F("hEff6Pt02","Efficiency; run number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEff6Pt02->SetLineWidth(2);
  hEff6Pt02->SetLineColor(2);
  hEff6Pt02->SetMarkerColor(2);
  hEff6Pt02->SetMarkerStyle(20);
  TH1F *hEff6Pt1 = new TH1F("hEff6Pt1","Efficiency; run number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEff6Pt1->SetLineWidth(2);
  hEff6Pt1->SetLineColor(2);
  hEff6Pt1->SetMarkerColor(2);
  hEff6Pt1->SetMarkerStyle(20);
  TH1F *hEff6Pt10 = new TH1F("hEff6Pt10","Efficiency; run number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEff6Pt10->SetLineWidth(2);
  hEff6Pt10->SetLineColor(2);
  hEff6Pt10->SetMarkerColor(2);
  hEff6Pt10->SetMarkerStyle(20);

  TH1F *hEffTOTPt02 = new TH1F("hEffTOTPt02","Efficiency; run number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEffTOTPt02->SetLineWidth(2);
  hEffTOTPt02->SetLineColor(kBlue+2);
  hEffTOTPt02->SetMarkerColor(kBlue+2);
  hEffTOTPt02->SetMarkerStyle(20);
  TH1F *hEffTOTPt1 = new TH1F("hEffTOTPt1","Efficiency; run number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEffTOTPt1->SetLineWidth(2);
  hEffTOTPt1->SetLineColor(kBlue+2);
  hEffTOTPt1->SetMarkerColor(kBlue+2);
  hEffTOTPt1->SetMarkerStyle(20);
  TH1F *hEffTOTPt10 = new TH1F("hEffTOTPt10","Efficiency; run number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEffTOTPt10->SetLineWidth(2);
  hEffTOTPt10->SetLineColor(kBlue+2);
  hEffTOTPt10->SetMarkerColor(kBlue+2);
  hEffTOTPt10->SetMarkerStyle(20);

  
  lista.Add(histonEvents);
  lista.Add(histonEventsTriggered);
  lista.Add(histominTime);
  lista.Add(histomeanTime);
  lista.Add(histofracExtra);
  lista.Add(histodEdxTB0);
  lista.Add(histodEdxTB5);
  lista.Add(histodEdxLay3); 
  lista.Add(histodEdxLay4);
  lista.Add(histoTrackClu1);
  lista.Add(histoTrackClu2);
  lista.Add(histoTrackClu3);
  lista.Add(histoTrackClu4);
  lista.Add(histoTrackClu5);
  lista.Add(histoTrackClu6);
  lista.Add(histoNmodEmpty);

  lista.Add(hMeanVx);
  lista.Add(hMeanVy);
  lista.Add(hMeanVz);
  lista.Add(hSigmaVx);
  lista.Add(hSigmaVy);
  lista.Add(hSigmaVz);
  lista.Add(hMeanVxSPD);
  lista.Add(hMeanVySPD); 
  lista.Add(hMeanVzSPD);
  lista.Add(hSigmaVxSPD);
  lista.Add(hSigmaVySPD);
  lista.Add(hSigmaVzSPD);

  lista.Add(histodEdxLay5);
  lista.Add(histodEdxLay6);
  lista.Add(histoChargeRatioLay5);
  lista.Add(histoChargeRatioLay6);
 
  lista.Add(hFracSPD1);
  lista.Add(hFracSPD2);
  lista.Add(hEffSPDPt02);
  lista.Add(hEffSPDPt1);
  lista.Add(hEffSPDPt10);
  lista.Add(hEffoneSPDPt02);
  lista.Add(hEffoneSPDPt1);
  lista.Add(hEffoneSPDPt10);
  lista.Add(hEff2Pt02);
  lista.Add(hEff2Pt1);
  lista.Add(hEff2Pt10);
  lista.Add(hEff3Pt02);
  lista.Add(hEff3Pt1);
  lista.Add(hEff3Pt10);
  lista.Add(hEff4Pt02);
  lista.Add(hEff4Pt1);
  lista.Add(hEff4Pt10);
  lista.Add(hEff5Pt02);
  lista.Add(hEff5Pt1);
  lista.Add(hEff5Pt10);
  lista.Add(hEff6Pt02);
  lista.Add(hEff6Pt1);
  lista.Add(hEff6Pt10);
  lista.Add(hEffTOTPt02);
  lista.Add(hEffTOTPt1);
  lista.Add(hEffTOTPt10);
  //lista.Add();
  

  for(Int_t i=0; i<nRuns;i++)
    {
      ttree->GetEntry(i);

      //======= SDD ========
      histoTrackClu1->SetBinContent(i+1,fracTrackWithClu1);
      histoTrackClu1->SetBinError(i+1,errfracTrackWithClu1);
      histoTrackClu1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      histoTrackClu2->SetBinContent(i+1,fracTrackWithClu2);
      histoTrackClu2->SetBinError(i+1,errfracTrackWithClu2);
      histoTrackClu2->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      histoTrackClu3->SetBinContent(i+1,fracTrackWithClu3);
      histoTrackClu3->SetBinError(i+1,errfracTrackWithClu3);
      histoTrackClu3->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      histoTrackClu4->SetBinContent(i+1,fracTrackWithClu4);
      histoTrackClu4->SetBinError(i+1,errfracTrackWithClu4);
      histoTrackClu4->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      histoTrackClu5->SetBinContent(i+1,fracTrackWithClu5);
      histoTrackClu5->SetBinError(i+1,errfracTrackWithClu5);
      histoTrackClu5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      histoTrackClu6->SetBinContent(i+1,fracTrackWithClu6);
      histoTrackClu6->SetBinError(i+1,errfracTrackWithClu6);
      histoTrackClu6->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      histonEvents->SetBinContent(i+1,nEvents);
      histonEventsTriggered->SetBinContent(i+1,nEventsTriggered);
      histonEvents->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      histonEventsTriggered->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      histominTime->SetBinContent(i+1,minDrTime);
      histominTime->SetBinError(i+1,errminDrTime);
      histominTime->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      histomeanTime->SetBinContent(i+1,meanDrTime);
      histomeanTime->SetBinError(i+1,errmeanDrTime);
      histomeanTime->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      histofracExtra->SetBinContent(i+1,fracExtra);
      histofracExtra->SetBinError(i+1,errfracExtra);
      histofracExtra->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      histodEdxTB0->SetBinContent(i+1,MPVdEdxTB0);
      histodEdxTB0->SetBinError(i+1,errMPVdEdxTB0);
      histodEdxTB0->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      histodEdxTB5->SetBinContent(i+1,MPVdEdxTB5);
      histodEdxTB5->SetBinError(i+1,errMPVdEdxTB5);
      histodEdxTB5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      histodEdxLay3->SetBinContent(i+1,MPVdEdxLay3);
      histodEdxLay3->SetBinError(i+1,errMPVdEdxLay3);
      histodEdxLay3->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      histodEdxLay4->SetBinContent(i+1,MPVdEdxLay4);
      histodEdxLay4->SetBinError(i+1,errMPVdEdxLay4);
      histodEdxLay4->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      histoNmodEmpty->SetBinContent(i+1,EmptyModulesSDD);
      histoNmodEmpty->SetBinError(i+1,0.000001);
      histoNmodEmpty->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));


      //======= vertex =========
      hMeanVx->SetBinContent(i+1,meanVtxTRKx);
      hMeanVx->SetBinError(i+1,meanVtxTRKxErr);
      hMeanVx->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun)); 
      hMeanVy->SetBinContent(i+1,meanVtxTRKy);
      hMeanVy->SetBinError(i+1,meanVtxTRKyErr);
      hMeanVy->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      hMeanVz->SetBinContent(i+1,meanVtxTRKz);
      hMeanVz->SetBinError(i+1,meanVtxTRKzErr);
      hMeanVz->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));   
      hSigmaVx->SetBinContent(i+1,sigmaVtxTRKx);
      hSigmaVx->SetBinError(i+1,sigmaVtxTRKxErr);
      hSigmaVx->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      hSigmaVy->SetBinContent(i+1,sigmaVtxTRKy);
      hSigmaVy->SetBinError(i+1,sigmaVtxTRKyErr);
      hSigmaVy->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      hSigmaVz->SetBinContent(i+1,sigmaVtxTRKz);
      hSigmaVz->SetBinError(i+1,sigmaVtxTRKzErr);
      hSigmaVz->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      hMeanVxSPD->SetBinContent(i+1,meanVtxSPDx);
      hMeanVxSPD->SetBinError(i+1,meanVtxSPDxErr);
      hMeanVxSPD->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      hMeanVySPD->SetBinContent(i+1,meanVtxSPDy);
      hMeanVySPD->SetBinError(i+1,meanVtxSPDyErr);
      hMeanVySPD->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      hMeanVzSPD->SetBinContent(i+1,meanVtxSPDz);
      hMeanVzSPD->SetBinError(i+1,meanVtxSPDzErr);
      hMeanVzSPD->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      hSigmaVxSPD->SetBinContent(i+1,sigmaVtxSPDx);
      hSigmaVxSPD->SetBinError(i+1,sigmaVtxSPDxErr);
      hSigmaVxSPD->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      hSigmaVySPD->SetBinContent(i+1,sigmaVtxSPDy);
      hSigmaVySPD->SetBinError(i+1,sigmaVtxSPDyErr);
      hSigmaVySPD->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      hSigmaVzSPD->SetBinContent(i+1,sigmaVtxSPDz);
      hSigmaVzSPD->SetBinError(i+1,sigmaVtxSPDzErr);
      hSigmaVzSPD->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));

      //======= SSD =======
      histodEdxLay5->SetBinContent(i+1,MPVL5);
      histodEdxLay5->SetBinError(i+1,MPVErrL5);
      histodEdxLay5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      
      histodEdxLay6->SetBinContent(i+1,MPVL6);
      histodEdxLay6->SetBinError(i+1,MPVErrL6);
      histodEdxLay6->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      
      histoChargeRatioLay5->SetBinContent(i+1,ChargeRatioL5);
      histoChargeRatioLay5->SetBinError(i+1,ChargeRatioErrL5);
      histoChargeRatioLay5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      
      histoChargeRatioLay6->SetBinContent(i+1,ChargeRatioL6);
      histoChargeRatioLay6->SetBinError(i+1,ChargeRatioErrL6);
      histoChargeRatioLay6->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));

      //======= MATCHING =======
      hFracSPD1->SetBinContent(i+1,FracSPD1);
      hFracSPD1->SetBinError(i+1,.01);
      hFracSPD1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      hFracSPD2->SetBinContent(i+1,FracSPD2);
      hFracSPD2->SetBinError(i+1,.01);
      hFracSPD2->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));

      hEff6Pt02->SetBinContent(i+1,Eff6Pt02);
      hEff6Pt02->SetBinError(i+1,errEff6Pt02);
      hEff6Pt02->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      hEff6Pt1->SetBinContent(i+1,Eff6Pt1);
      hEff6Pt1->SetBinError(i+1,errEff6Pt1);
      hEff6Pt1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      hEff6Pt10->SetBinContent(i+1,Eff6Pt10);
      hEff6Pt10->SetBinError(i+1,errEff6Pt10);
      hEff6Pt10->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      
      hEff5Pt02->SetBinContent(i+1,Eff5Pt02);
      hEff5Pt02->SetBinError(i+1,errEff5Pt02);
      hEff5Pt02->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      hEff5Pt1->SetBinContent(i+1,Eff5Pt1);
      hEff5Pt1->SetBinError(i+1,errEff5Pt1);
      hEff5Pt1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      hEff5Pt10->SetBinContent(i+1,Eff5Pt10);
      hEff5Pt10->SetBinError(i+1,errEff5Pt10);
      hEff5Pt10->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      
      hEff4Pt02->SetBinContent(i+1,Eff4Pt02);
      hEff4Pt02->SetBinError(i+1,errEff4Pt02);
      hEff4Pt02->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      hEff4Pt1->SetBinContent(i+1,Eff4Pt1);
      hEff4Pt1->SetBinError(i+1,errEff4Pt1);
      hEff4Pt1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      hEff4Pt10->SetBinContent(i+1,Eff4Pt10);
      hEff4Pt10->SetBinError(i+1,errEff4Pt10);
      hEff4Pt10->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      
      hEff3Pt02->SetBinContent(i+1,Eff3Pt02);
      hEff3Pt02->SetBinError(i+1,errEff3Pt02);
      hEff3Pt02->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      hEff3Pt1->SetBinContent(i+1,Eff3Pt1);
      hEff3Pt1->SetBinError(i+1,errEff3Pt1);
      hEff3Pt1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      hEff3Pt10->SetBinContent(i+1,Eff3Pt10);
      hEff3Pt10->SetBinError(i+1,errEff3Pt10);
      hEff3Pt10->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));      
      
      hEff2Pt02->SetBinContent(i+1,Eff2Pt02);
      hEff2Pt02->SetBinError(i+1,errEff2Pt02);
      hEff2Pt02->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));  
      hEff2Pt1->SetBinContent(i+1,Eff2Pt1);
      hEff2Pt1->SetBinError(i+1,errEff2Pt1);
      hEff2Pt1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));  
      hEff2Pt10->SetBinContent(i+1,Eff2Pt10);
      hEff2Pt10->SetBinError(i+1,errEff2Pt10);
      hEff2Pt10->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));  
        
      hEffoneSPDPt02->SetBinContent(i+1,EffoneSPDPt02);
      hEffoneSPDPt02->SetBinError(i+1,errEffoneSPDPt02);
      hEffoneSPDPt02->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      hEffoneSPDPt1->SetBinContent(i+1,EffoneSPDPt1);
      hEffoneSPDPt1->SetBinError(i+1,errEffoneSPDPt1);
      hEffoneSPDPt1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      hEffoneSPDPt10->SetBinContent(i+1,EffoneSPDPt10);
      hEffoneSPDPt10->SetBinError(i+1,errEffoneSPDPt10);
      hEffoneSPDPt10->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));

      hEffSPDPt02->SetBinContent(i+1,EffSPDPt02);
      hEffSPDPt02->SetBinError(i+1,errEffSPDPt02);
      hEffSPDPt02->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));  
      hEffSPDPt1->SetBinContent(i+1,EffSPDPt1);
      hEffSPDPt1->SetBinError(i+1,errEffSPDPt1);
      hEffSPDPt1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun)); 
      hEffSPDPt10->SetBinContent(i+1,EffSPDPt10);
      hEffSPDPt10->SetBinError(i+1,errEffSPDPt10);
      hEffSPDPt10->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun)); 
      
      hEffTOTPt02->SetBinContent(i+1,EffTOTPt02);
      hEffTOTPt02->SetBinError(i+1,errEffTOTPt02);
      hEffTOTPt02->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      hEffTOTPt1->SetBinContent(i+1,EffTOTPt1);
      hEffTOTPt1->SetBinError(i+1,errEffTOTPt1);
      hEffTOTPt1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      hEffTOTPt10->SetBinContent(i+1,EffTOTPt10);
      hEffTOTPt10->SetBinError(i+1,errEffTOTPt10);
      hEffTOTPt10->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    
    }

  //--------  Draw Vertex histograms ---------
  TCanvas *cVertexDisto=new TCanvas("cVertexDisto","cVertexDisto",1200,800);
  cVertexDisto->Divide(3,2);
  
  cVertexDisto->cd(1);
  hMeanVx->SetLineWidth(2);
  hMeanVx->SetLineColor(kBlue+2);
  hMeanVx->SetMarkerColor(kBlue+2);
  hMeanVx->SetMarkerStyle(20);
  hMeanVx->SetMinimum(-0.1);
  hMeanVx->SetMaximum(0.1);
  hMeanVx->GetYaxis()->SetTitle("Vertex X coordinate");
  hMeanVx->Draw();
  hMeanVxSPD->SetLineWidth(2);
  hMeanVxSPD->SetLineColor(2);
  hMeanVxSPD->SetMarkerColor(2);
  hMeanVxSPD->SetMarkerStyle(20);
  hMeanVxSPD->Draw("same");

  cVertexDisto->cd(2);
  hMeanVy->SetLineWidth(2);
  hMeanVy->SetLineColor(kBlue+2);
  hMeanVy->SetMarkerColor(kBlue+2);
  hMeanVy->SetMarkerStyle(20);
  hMeanVy->SetMinimum(-0.4);
  hMeanVy->SetMaximum(0.4);
  hMeanVy->GetYaxis()->SetTitle("Vertex Y coordinate");
  hMeanVy->Draw();
  hMeanVySPD->SetLineWidth(2);
  hMeanVySPD->SetLineColor(2);
  hMeanVySPD->SetMarkerColor(2);
  hMeanVySPD->SetMarkerStyle(20);
  hMeanVySPD->Draw("same");
  
  cVertexDisto->cd(3);
  hMeanVz->SetLineWidth(2);
  hMeanVz->SetLineColor(kBlue+2);
  hMeanVz->SetMarkerColor(kBlue+2);
  hMeanVz->SetMarkerStyle(20);
     hMeanVz->SetMinimum(-9.);
     hMeanVz->SetMaximum(9.);
  hMeanVz->GetYaxis()->SetTitle("Vertex Z coordinate");
  hMeanVz->Draw();
  hMeanVzSPD->SetLineWidth(2);
  hMeanVzSPD->SetLineColor(2);
  hMeanVzSPD->SetMarkerColor(2);
  hMeanVzSPD->SetMarkerStyle(20);
  hMeanVzSPD->Draw("same");

  cVertexDisto->cd(4);
  hSigmaVx->SetLineWidth(2);
  hSigmaVx->SetLineColor(kBlue+2);
  hSigmaVx->SetMarkerColor(kBlue+2);
  hSigmaVx->SetMarkerStyle(20);
  hSigmaVx->SetMinimum(0.);
  hSigmaVx->SetMaximum(0.2);
  hSigmaVx->GetYaxis()->SetTitle("sigma on x coordinate");
  hSigmaVx->Draw();
  hSigmaVxSPD->SetLineWidth(2);
  hSigmaVxSPD->SetLineColor(2);
  hSigmaVxSPD->SetMarkerColor(2);
  hSigmaVxSPD->SetMarkerStyle(20);
  hSigmaVxSPD->Draw("same");

  cVertexDisto->cd(5);
  hSigmaVy->SetLineWidth(2);
  hSigmaVy->SetLineColor(kBlue+2);
  hSigmaVy->SetMarkerColor(kBlue+2);
  hSigmaVy->SetMarkerStyle(20);
  hSigmaVy->SetMinimum(0.);
  hSigmaVy->SetMaximum(0.2);
  hSigmaVy->GetYaxis()->SetTitle("sigma on y coordinate");
  hSigmaVy->Draw();
  hSigmaVySPD->SetLineWidth(2);
  hSigmaVySPD->SetLineColor(2);
  hSigmaVySPD->SetMarkerColor(2);
  hSigmaVySPD->SetMarkerStyle(20);
  hSigmaVySPD->Draw("same");

  cVertexDisto->cd(6);
  hSigmaVz->SetLineWidth(2);
  hSigmaVz->SetLineColor(kBlue+2);
  hSigmaVz->SetMarkerColor(kBlue+2);
  hSigmaVz->SetMarkerStyle(20);
    hSigmaVz->SetMinimum(-2.);
    hSigmaVz->SetMaximum(10.);
  hSigmaVz->GetYaxis()->SetTitle("sigma on z coordinate");
  hSigmaVz->Draw();
  hSigmaVzSPD->SetLineWidth(2);
  hSigmaVzSPD->SetLineColor(2);
  hSigmaVzSPD->SetMarkerColor(2);
  hSigmaVzSPD->SetMarkerStyle(20);
  hSigmaVzSPD->Draw("same");
 
  cVertexDisto->SaveAs("Vertex_trend.pdf");
  //pdfFileNames+=" Vertex_trend.pdf";
  cVertexDisto->Update();

  
  TCanvas* cFracTrackWithClLay=new TCanvas("cFracTrackWithClLay","Track With points in ITS");
  histoTrackClu1->SetLineColor(kGray+1);
  histoTrackClu1->SetMarkerColor(kGray+1);
  histoTrackClu1->SetMarkerStyle(24);
    histoTrackClu1->SetMinimum(0);
    histoTrackClu1->SetMaximum(1);
  histoTrackClu1->Draw();
  histoTrackClu2->SetLineColor(kGray+2);
  histoTrackClu2->SetMarkerColor(kGray+2);
  histoTrackClu2->SetMarkerStyle(26);
  histoTrackClu2->Draw("same");
  histoTrackClu3->Draw("same");
  histoTrackClu3->SetLineColor(1);
  histoTrackClu3->SetMarkerStyle(20);
  histoTrackClu3->Draw("same");
  histoTrackClu3->SetMinimum(0.);
  histoTrackClu3->SetMaximum(1.05);
  histoTrackClu4->SetLineColor(2);
  histoTrackClu4->SetMarkerColor(2);
  histoTrackClu4->SetMarkerStyle(22);
  histoTrackClu4->Draw("same");
  histoTrackClu5->SetLineColor(4);
  histoTrackClu5->SetMarkerColor(4);
  histoTrackClu5->SetMarkerStyle(29);
  histoTrackClu5->Draw("same");
  histoTrackClu6->SetLineColor(kBlue+1);
  histoTrackClu6->SetMarkerColor(kBlue+1);
  histoTrackClu6->SetMarkerStyle(30);
  histoTrackClu6->SetTitle("Fraction of tracks with points in ITS");
  histoTrackClu6->Draw("same");
  histoTrackClu3->GetYaxis()->SetTitle("Fraction of Tracks with Points in ITS Layers");
    histoTrackClu3->SetMinimum(0);
    histoTrackClu3->SetMaximum(1);
    histoTrackClu6->SetMinimum(0);
    histoTrackClu6->SetMaximum(1);
    histoTrackClu1->GetYaxis()->SetRangeUser(0,1);


  TLegend* leg3=new TLegend(0.7,0.15,0.88,0.35);
  TLegendEntry* ent;
  ent=leg3->AddEntry(histoTrackClu1,"Layer1","PL");
  ent->SetTextColor(histoTrackClu1->GetMarkerColor());
  ent=leg3->AddEntry(histoTrackClu2,"Layer2","PL");
  ent->SetTextColor(histoTrackClu2->GetMarkerColor());
  ent=leg3->AddEntry(histoTrackClu3,"Layer3","PL");
  ent->SetTextColor(histoTrackClu3->GetMarkerColor());
  ent=leg3->AddEntry(histoTrackClu4,"Layer4","PL");
  ent->SetTextColor(histoTrackClu4->GetMarkerColor());
  ent=leg3->AddEntry(histoTrackClu5,"Layer5","PL");
  ent->SetTextColor(histoTrackClu5->GetMarkerColor());
  ent=leg3->AddEntry(histoTrackClu6,"Layer6","PL");
  ent->SetTextColor(histoTrackClu6->GetMarkerColor());

  leg3->SetFillStyle(0);
  leg3->Draw();
   cFracTrackWithClLay->SaveAs("TrackPoints_trend.pdf");
 //  pdfFileNames+=" TrackPoints_trend.pdf";
  cFracTrackWithClLay->Update();



  TCanvas* cNEvents=new TCanvas("cNEvents","Number of Events",1200,800);
  histonEvents->Draw();
  histonEventsTriggered->Draw("same");
   
  TCanvas* cDriftTimeCharge=new TCanvas("cDriftTimeCharge","SDD DriftTime & Charge",1200,800);
  cDriftTimeCharge->Divide(2,2);

  cDriftTimeCharge->cd(1);
  histominTime->Draw();
  histominTime->SetMinimum(450);
  histominTime->SetMaximum(550);
  histominTime->GetYaxis()->SetTitle("Minimum Drift Time (ns)");
  TLatex* td1=new TLatex(0.2,0.85,"SDD minimum drift time (ref=505ns)");
  td1->SetNDC();
  td1->SetTextColor(1);
  td1->Draw();

  cDriftTimeCharge->cd(2);
  histomeanTime->Draw();
  histomeanTime->SetMinimum(3200);
  histomeanTime->SetMaximum(3300);
  histomeanTime->GetYaxis()->SetTitle("Average Drift Time (ns)");
  TLatex* td2=new TLatex(0.2,0.85,"SDD average drift time");
  td2->SetNDC();
  td2->Draw();

  cDriftTimeCharge->cd(3);
  gPad->SetGridy();
  histodEdxTB0->SetLineColor(1);
  histodEdxTB0->SetMarkerStyle(20);
  histodEdxTB0->Draw();
  histodEdxTB0->SetMinimum(70.);
  histodEdxTB0->SetMaximum(100.);
  histodEdxTB5->SetLineColor(4);
  histodEdxTB5->SetMarkerColor(4);
  histodEdxTB5->SetMarkerStyle(23);
  histodEdxTB5->Draw("same");
  histodEdxTB0->GetYaxis()->SetTitle("MPV of dE/dx (keV/300 #mum)");

  TLegend* leg2=new TLegend(0.6,0.15,0.88,0.35);
  ent=leg2->AddEntry(histodEdxTB0,"Small drift time","PL");
  ent=leg2->AddEntry(histodEdxTB5,"Large drift time","PL");
  ent->SetTextColor(histodEdxTB5->GetMarkerColor());
  leg2->SetFillStyle(0);
  leg2->Draw();

  TLatex* tc1=new TLatex(0.2,0.85,"SDD charge in different drift regions");
  tc1->SetNDC();
  tc1->SetTextColor(1);
  tc1->Draw();

  cDriftTimeCharge->cd(4);
  gPad->SetGridy();
  histodEdxLay3->SetLineColor(1);
  histodEdxLay3->SetMarkerStyle(20);
  histodEdxLay3->GetYaxis()->SetTitle("MPV of dE/dx (keV/300 #mum)");
  histodEdxLay3->Draw();
  histodEdxLay3->SetMinimum(70.);
  histodEdxLay3->SetMaximum(100.);
  histodEdxLay4->SetLineColor(4);
  histodEdxLay4->SetMarkerColor(4);
  histodEdxLay4->SetMarkerStyle(23);
  histodEdxLay4->Draw("same");

  histodEdxLay5->SetLineColor(6);
  histodEdxLay5->SetMarkerColor(6);
  histodEdxLay5->SetMarkerStyle(22);
  histodEdxLay5->Draw("same");
  histodEdxLay5->SetMinimum(0.);
  histodEdxLay6->SetLineColor(7);
  histodEdxLay6->SetMarkerColor(7);
  histodEdxLay6->SetMarkerStyle(24);
  histodEdxLay6->Draw("same");
    
  TLegend* leg2b=new TLegend(0.6,0.15,0.88,0.35);
  ent=leg2b->AddEntry(histodEdxLay3,"Layer 3","PL");
  ent=leg2b->AddEntry(histodEdxLay4,"Layer 4","PL");
  ent=leg2b->AddEntry(histodEdxLay5,"Layer 5","PL");
  ent=leg2b->AddEntry(histodEdxLay6,"Layer 6","PL");
  ent->SetTextColor(histodEdxLay4->GetMarkerColor());
  leg2b->SetFillStyle(0);
  leg2b->Draw();
  TLatex* tcDriftTimeCharge=new TLatex(0.2,0.85,"SDD and SSD charge in different layers");
  tcDriftTimeCharge->SetNDC();
  tcDriftTimeCharge->SetTextColor(1);
  tcDriftTimeCharge->Draw();

  cDriftTimeCharge->SaveAs("SDD_SSD_drift_charge_trend.pdf");
//   pdfFileNames+=" SDD_SSD_drift_charge_trend.pdf";
  cDriftTimeCharge->Update();


  TCanvas *cChargeRatio = new TCanvas("cChargeRatio","Charge ratio");
  cChargeRatio->cd();
  histoChargeRatioLay5->SetLineColor(6);
  histoChargeRatioLay5->SetMarkerColor(6);
  histoChargeRatioLay5->SetMarkerStyle(20);
  histoChargeRatioLay5->SetMinimum(-0.01);
  histoChargeRatioLay5->SetMaximum(+0.01);
  histoChargeRatioLay5->GetYaxis()->SetTitle("SSD charge ratio");
  histoChargeRatioLay5->Draw();

  histoChargeRatioLay6->SetLineColor(7);
  histoChargeRatioLay6->SetMarkerColor(7);
  histoChargeRatioLay6->SetMarkerStyle(22);
  histoChargeRatioLay6->Draw("same");

  TLegend* legCR=new TLegend(0.7,0.65,0.88,0.75);
  ent=legCR->AddEntry(histoChargeRatioLay5,"Layer5","PL");
  ent->SetTextColor(histoChargeRatioLay5->GetMarkerColor());
  ent=legCR->AddEntry(histoChargeRatioLay6,"Layer6","PL");
  ent->SetTextColor(histoChargeRatioLay6->GetMarkerColor());
  legCR->SetFillStyle(0);
  legCR->Draw();
  TLatex* tc3=new TLatex(0.2,0.85,"SSD charge ratio in different layers");
  tc3->SetNDC();
  tc3->SetTextColor(1);
  tc3->Draw();
  
  cChargeRatio->SaveAs("SSD_chargeratio_trend.pdf");
 // pdfFileNames+=" SSD_chargeratio_trend.pdf";
  cChargeRatio->Update();

  
  TCanvas *cTPCMatching=new TCanvas("cTPCMatching","TPC-ITS matching efficiency",1200,300);
  cTPCMatching->Divide(3,1);

  cTPCMatching->cd(1);
  hEff6Pt02->SetMinimum(0);
  hEff6Pt02->GetYaxis()->SetRangeUser(0,1);
  hEff6Pt02->Draw();
  hEff5Pt02->Draw("same");
  hEff4Pt02->Draw("same");
  hEff3Pt02->Draw("same");
  hEff2Pt02->Draw("same");
  hEffSPDPt02->Draw("same");  
  hEffoneSPDPt02->Draw("same");
  hEffTOTPt02->Draw("same");

  TLegend* lpt02=new TLegend(0.9,0.8,1,1);
  lpt02->AddEntry(hEff6Pt02,"6 cls","l");
  lpt02->AddEntry(hEff5Pt02,"5 cls","l");
  lpt02->AddEntry(hEff4Pt02,"4 cls","l");
  lpt02->AddEntry(hEff3Pt02,"3 cls","l");
  lpt02->AddEntry(hEff2Pt02,"2 cls","l");
  lpt02->AddEntry(hEffSPDPt02,"2SPD + any","l");
  lpt02->AddEntry(hEffoneSPDPt02,">=1SPD + any","l");
  lpt02->AddEntry(hEffTOTPt02,">=2","l");
  lpt02->Draw("same");
  TLatex* tpc1=new TLatex(0.2,0.85,"TPC-ITS match eff Pt=0.2");
  tpc1->SetNDC();
  tpc1->SetTextColor(1);
  tpc1->Draw();

  cTPCMatching->cd(2);
  hEff6Pt1->GetYaxis()->SetRangeUser(0,1);
  hEff6Pt1->Draw();
  hEff5Pt1->Draw("same");
  hEff4Pt1->Draw("same");
  hEff3Pt1->Draw("same");
  hEff2Pt1->Draw("same");
  hEffSPDPt1->Draw("same");
  hEffoneSPDPt1->Draw("same");
  hEffTOTPt1->Draw("same");
  
  TLegend* lpt1=new TLegend(0.9,0.8,1,1);
  lpt1->AddEntry(hEff6Pt1,"6 cls","l");
  lpt1->AddEntry(hEff5Pt1,"5 cls","l");
  lpt1->AddEntry(hEff4Pt1,"4 cls","l");
  lpt1->AddEntry(hEff3Pt1,"3 cls","l");
  lpt1->AddEntry(hEff2Pt1,"2 cls","l");
  lpt1->AddEntry(hEffSPDPt1,"2SPD + any","l");
  lpt1->AddEntry(hEffoneSPDPt1,">=1SPD + any","l");
  lpt1->AddEntry(hEffTOTPt02,">=2","l");
  lpt1->Draw("same");
  TLatex* tpc2=new TLatex(0.2,0.75,"TPC-ITS match eff Pt=1");
  tpc2->SetNDC();
  tpc2->SetTextColor(1);
  tpc2->Draw();
  
  cTPCMatching->cd(3);
  hEff6Pt10->GetYaxis()->SetRangeUser(0,1);
  hEff6Pt10->Draw();
  hEff5Pt10->Draw("same");
  hEff4Pt10->Draw("same");
  hEff3Pt10->Draw("same");
  hEff2Pt10->Draw("same");
  hEffSPDPt10->Draw("same");
  hEffoneSPDPt10->Draw("same");
  hEffTOTPt10->Draw("same");
  
  TLegend* lpt10=new TLegend(0.9,0.8,1,1);
  lpt10->AddEntry(hEff6Pt10,"6 cls","l");
  lpt10->AddEntry(hEff5Pt10,"5 cls","l");
  lpt10->AddEntry(hEff4Pt10,"4 cls","l");
  lpt10->AddEntry(hEff3Pt10,"3 cls","l");
  lpt10->AddEntry(hEff2Pt10,"2 cls","l");
  lpt10->AddEntry(hEffSPDPt10,"2SPD + any","l");
  lpt10->AddEntry(hEffoneSPDPt10,">=1SPD + any","l");
  lpt10->AddEntry(hEffTOTPt02,">=2","l");
  lpt10->Draw("same");
  TLatex* tpc3=new TLatex(0.2,0.75,"TPC-ITS match eff Pt=10");
  tpc3->SetNDC();
  tpc3->SetTextColor(1);
  tpc3->Draw();

  cTPCMatching->SaveAs("TPC_ITSMatch_trend.pdf");
 // pdfFileNames+=" TPC_ITSMatch_trend.pdf";
  cTPCMatching->Update();
  
  
  TCanvas *cPixel=new TCanvas("cPixel","SPD on");
  hFracSPD1->SetMaximum(1.2);
  hFracSPD1->SetMinimum(0);
  hFracSPD1->Draw("p");
  hFracSPD2->Draw("same,p");

  TLegend* lSPD=new TLegend(0.9,0.8,1,1);
  lSPD->AddEntry(hFracSPD1,"Frac. SPD1 ON","l");
  lSPD->AddEntry(hFracSPD2,"Frac. SPD2 ON","l");
  lSPD->Draw();
  TLatex* tSPD=new TLatex(0.2,0.85,"Fraction of SPD half staves ON");
  tSPD->SetNDC();
  tSPD->SetTextColor(1);
  tSPD->Draw();

  cPixel->SaveAs("Pixel_trend.pdf");
//  pdfFileNames+=" Pixel_trend.pdf";
  cPixel->Update();



  TFile * fout=new TFile(outfilename,"recreate");
  fout->cd();
  lista.Write();
  fout->Close();

  
  return 0;
}
