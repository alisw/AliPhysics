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
TString pdfFileNames=" ";

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


  //************************************************ SDD  (41) ********************************************//
  
    Int_t nrun,nEvents,nEventsTriggered;
    Float_t minDrTime,errminDrTime,meanDrTime,errmeanDrTime;
    Float_t fracTrackWithClu1,fracTrackWithClu2,errfracTrackWithClu1,errfracTrackWithClu2;
    Float_t fracTrackWithClu3,fracTrackWithClu4,errfracTrackWithClu3,errfracTrackWithClu4;
    Float_t fracTrackWithClu5,fracTrackWithClu6,errfracTrackWithClu5,errfracTrackWithClu6;
    Float_t fracExtra,errfracExtra,fracEvWithSDD,errfracEvWithSDD;
    Float_t MPVdEdxLay3,errMPVdEdxLay3,MPVdEdxLay4,errMPVdEdxLay4;
    Float_t MPVdEdxTB0,errMPVdEdxTB0,MPVdEdxTB5,errMPVdEdxTB5;
    Float_t fracDead3,errfracDead3,fracDead4,errfracDead4;
    Float_t FlagSDD1, FlagSDD2; // flag on fraction of SDD anodes ON
    Float_t FlagMinTime, FlagMeanTime, FlagdEdx3, FlagdEdx4; // flag on SDD time and charge parameters
  
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
//  ttree->SetBranchAddress("EmptyModulesSDD",&EmptyModulesSDD); // Number of empty SSD  modules
    ttree->SetBranchAddress("fracExtra",&fracExtra); // fraction of extra clusters in SDD
    ttree->SetBranchAddress("errfracExtra",&errfracExtra); // fraction of extra clusters in SDD
//    ttree->SetBranchAddress("fracEvWithSDD",&fracEvWithSDD); // fraction of tracks with SSD in DAQ
//    ttree->SetBranchAddress("errfracEvWithSDD",&errfracEvWithSDD); // fraction of tracks with SSD in DAQ
  ttree->SetBranchAddress("MPVdEdxLay3",&MPVdEdxLay3); // most probable value of dE/dx distribution of SDD Layer 3
  ttree->SetBranchAddress("errMPVdEdxLay3",&errMPVdEdxLay3); // error  most probable value of dE/dx distribution of SDD Layer 3
  ttree->SetBranchAddress("MPVdEdxLay4",&MPVdEdxLay4); // most probable value of dE/dx distribution of SDD Layer 4
  ttree->SetBranchAddress("errMPVdEdxLay4",&errMPVdEdxLay4); // error  most probable value of dE/dx distribution of SDD Layer 4
  ttree->SetBranchAddress("MPVdEdxTB0",&MPVdEdxTB0); // most probable value of dE/dx distribution of SDD - small drift time
  ttree->SetBranchAddress("errMPVdEdxTB0",&errMPVdEdxTB0); // most probable value of dE/dx distribution of SDD - small drift time
  ttree->SetBranchAddress("MPVdEdxTB5",&MPVdEdxTB5); // most probable value of dE/dx distribution of SDD - large drift time
  ttree->SetBranchAddress("errMPVdEdxTB5",&errMPVdEdxTB5); // most probable value of dE/dx distribution of SDD - large drift time    
    ttree->SetBranchAddress("fracDead3",&fracDead3); // fraction of bad SDD modules layer 3
    ttree->SetBranchAddress("errfracDead3",&errfracDead3); // fraction of bad SDD modules layer 3
    ttree->SetBranchAddress("fracDead4",&fracDead4); // fraction of bad SDD modules layer 4
    ttree->SetBranchAddress("errfracDead4",&errfracDead4); // fraction of bad SDD modules layer 4
    ttree->SetBranchAddress("FlagSDD1",&FlagSDD1); // flag on fraction of SDD1 anodes ON
    ttree->SetBranchAddress("FlagSDD2",&FlagSDD2); // flag on fraction of SDD2 anodes ON
    ttree->SetBranchAddress("FlagMinTime",&FlagMinTime); // flag on min drift time
    ttree->SetBranchAddress("FlagMeanTime",&FlagMeanTime); // flag on mean drift time
    ttree->SetBranchAddress("FlagdEdx3",&FlagdEdx3); // flag on layer 3 dEdx MPV
    ttree->SetBranchAddress("FlagdEdx4",&FlagdEdx4); // flag on layer 4 dEdx MPV


  //************************************************ VERTEX  (26) ************************************************//
    
  Float_t meanVtxTRKx,meanVtxTRKy,meanVtxTRKz;
  Float_t meanVtxSPDx,meanVtxSPDy,meanVtxSPDz;
  Float_t sigmaVtxTRKx,sigmaVtxTRKy,sigmaVtxTRKz;
  Float_t sigmaVtxSPDx,sigmaVtxSPDy,sigmaVtxSPDz;
  Float_t meanVtxTRKxErr,meanVtxTRKyErr,meanVtxTRKzErr;
  Float_t meanVtxSPDxErr,meanVtxSPDyErr,meanVtxSPDzErr;
  Float_t sigmaVtxTRKxErr,sigmaVtxTRKyErr,sigmaVtxTRKzErr;
  Float_t sigmaVtxSPDxErr,sigmaVtxSPDyErr,sigmaVtxSPDzErr;
  Float_t meanVtxOCDBx,meanVtxOCDBy,meanVtxOCDBz;
  Float_t sigmaVtxOCDBx,sigmaVtxOCDBy,sigmaVtxOCDBz;
  Float_t diamondX=-999.,diamondY=-999.,diamondZ=-999.;
  Float_t diamondSigX=-999.,diamondSigY=-999.,diamondSigZ=-999.;
  Float_t pileupSPD, errpileupSPD;

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
  ttree->SetBranchAddress("diamondX",&meanVtxOCDBx); // mean of OCDB vertex position - x
  ttree->SetBranchAddress("diamondY",&meanVtxOCDBy); // mean of OCDB vertex position - y
  ttree->SetBranchAddress("diamondZ",&meanVtxOCDBz); // mean of OCDB vertex position - z
  ttree->SetBranchAddress("diamondSigX",&sigmaVtxOCDBx); // sigma of OCDB vertex position - x
  ttree->SetBranchAddress("diamondSigY",&sigmaVtxOCDBy); // sigma of OCDB vertex position - y
  ttree->SetBranchAddress("diamondSigZ",&sigmaVtxOCDBz); // sigma of CODB vertex position - z

  ttree->SetBranchAddress("pileupSPD",&pileupSPD); // fraction of events with SPD pileup vertex
  ttree->SetBranchAddress("errpileupSPD",&errpileupSPD); // fraction of events with SPD pileup vertex


  //************************************************ SSD (25) **************************************************//
     
    Float_t MPVL5,MPVErrL5;
    Float_t MPVL6,MPVErrL6;
    Float_t ChargeRatioL5,ChargeRatioErrL5;
    Float_t ChargeRatioL6,ChargeRatioErrL6;
    Float_t EmptyModulesSSD;
    Float_t FracBadn5,errFracBadn5,FracBadp5,errFracBadp5,FracBadn6,errFracBadn6,FracBadp6,errFracBadp6;
    Float_t FlagSSD1n,FlagSSD1p,FlagSSD2n,FlagSSD2p; // flag on fraction of SSD bad strips
    Float_t FlagChR5,FlagChR6,FlagdEdx5,FlagdEdx6; // flag on SSD CR and MPV

    ttree->SetBranchAddress("MPVL5",&MPVL5); // Most Probable Value dEdx Layer 5
    ttree->SetBranchAddress("MPVErrL5",&MPVErrL5); // Most Probable Value error dEdx Layer 5
    ttree->SetBranchAddress("MPVL6",&MPVL6); // Most Probable Value dEdx Layer 6
    ttree->SetBranchAddress("MPVErrL6",&MPVErrL6); // Most Probable Value error dEdx Layer 6
    ttree->SetBranchAddress("ChargeRatioL5",&ChargeRatioL5); // Charge ratio (2 sides of SSD) Layer 5
    ttree->SetBranchAddress("ChargeRatioErrL5",&ChargeRatioErrL5); // Charge ratio error (2 sides of SSD) Layer 5
    ttree->SetBranchAddress("ChargeRatioL6",&ChargeRatioL6); // Charge ratio (2 sides of SSD) Layer 6
    ttree->SetBranchAddress("ChargeRatioErrL6",&ChargeRatioErrL6); // Charge ratio error(2 sides of SSD) Layer 6
    ttree->SetBranchAddress("EmptyModulesSSD",&EmptyModulesSSD); // Number of empty SSD  modules
    ttree->SetBranchAddress("FracBadn5",&FracBadn5); // fraction of bad n-strips layer 5
    ttree->SetBranchAddress("errFracBadn5",&errFracBadn5); // fraction of bad n-strips layer 5
    ttree->SetBranchAddress("FracBadp5",&FracBadp5); // fraction of bad p-strips layer 5
    ttree->SetBranchAddress("errFracBadp5",&errFracBadp5); // fraction of bad p-strips layer 5
    ttree->SetBranchAddress("FracBadn6",&FracBadn6); // fraction of bad n-strips layer 6
    ttree->SetBranchAddress("errFracBadn6",&errFracBadn6); // fraction of bad n-strips layer 6
    ttree->SetBranchAddress("FracBadp6",&FracBadp6); // fraction of bad p-strips layer 6
    ttree->SetBranchAddress("errFracBadp6",&errFracBadp6); // fraction of bad p-strips layer 6
    ttree->SetBranchAddress("FlagSSD1n",&FlagSSD1n); // flag on layer 5 bad n-strips fraction
    ttree->SetBranchAddress("FlagSSD1p",&FlagSSD1p); // flag on layer 5 bad p-strips fraction
    ttree->SetBranchAddress("FlagSSD2n",&FlagSSD2n); // flag on layer 6 bad n-strips fraction
    ttree->SetBranchAddress("FlagSSD2p",&FlagSSD2p); // flag on layer 6 bad p-strips fraction
    ttree->SetBranchAddress("FlagChR5",&FlagChR5); // flag on layer 5 charge ratio value
    ttree->SetBranchAddress("FlagChR6",&FlagChR6); // flag on layer 6 charge ratio value
    ttree->SetBranchAddress("FlagdEdx5",&FlagdEdx5); // flag on layer 5 dEdx MPV
    ttree->SetBranchAddress("FlagdEdx6",&FlagdEdx6); // flag on layer 6 dEdx MPV


    //************************************************ SPD (4) *********************************************//
    

    Float_t FracSPD1,errFracSPD1,FracSPD2,errFracSPD2,FlagSPD1,FlagSPD2;

    ttree->SetBranchAddress("FracSPD1",&FracSPD1); // fraction of SPD1 HS ON
    ttree->SetBranchAddress("errFracSPD1",&errFracSPD1); // fraction of SPD1 HS ON
    ttree->SetBranchAddress("FracSPD2",&FracSPD2); // fraction of SPD2 HS ON
    ttree->SetBranchAddress("errFracSPD2",&errFracSPD2); // fraction of SPD2 HS ON
    ttree->SetBranchAddress("FlagSPD1",&FlagSPD1); // flag on fraction of SPD1 HS ON
    ttree->SetBranchAddress("FlagSPD2",&FlagSPD2); // flag on fraction of SPD2 HS ON

    //************************************************ MATCHING *********************************************//
  
    ///// TPC-ITS (52+24+2 variables)
    Float_t Eff6Pt02,errEff6Pt02,Eff6Pt1,errEff6Pt1,Eff6Pt10,errEff6Pt10;
    Float_t Eff5Pt02,errEff5Pt02,Eff5Pt1,errEff5Pt1,Eff5Pt10,errEff5Pt10;
    Float_t Eff4Pt02,errEff4Pt02,Eff4Pt1,errEff4Pt1,Eff4Pt10,errEff4Pt10;
    Float_t Eff3Pt02,errEff3Pt02,Eff3Pt1,errEff3Pt1,Eff3Pt10,errEff3Pt10;
    Float_t Eff2Pt02,errEff2Pt02,Eff2Pt1,errEff2Pt1,Eff2Pt10,errEff2Pt10;
    Float_t EffSPDPt02,errEffSPDPt02,EffSPDPt1,errEffSPDPt1,EffSPDPt10,errEffSPDPt10;
    Float_t EffoneSPDPt02,errEffoneSPDPt02,EffoneSPDPt1,errEffoneSPDPt1,EffoneSPDPt10,errEffoneSPDPt10;
    Float_t EffTOTPt02,errEffTOTPt02,EffTOTPt1,errEffTOTPt1,EffTOTPt10,errEffTOTPt10;
    Float_t FracTrackMI1,errFracTrackMI1,FracTrackMI2,errFracTrackMI2,FracTrackMI3,errFracTrackMI3;
    Float_t FracTrackMI4,errFracTrackMI4,FracTrackMI5,errFracTrackMI5,FracTrackMI6,errFracTrackMI6;
    Float_t FracTrackSA1,errFracTrackSA1,FracTrackSA2,errFracTrackSA2,FracTrackSA3,errFracTrackSA3;
    Float_t FracTrackSA4,errFracTrackSA4,FracTrackSA5,errFracTrackSA5,FracTrackSA6,errFracTrackSA6;
    
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
    
    ttree->SetBranchAddress("FracTrackMI1",&FracTrackMI1); // fraction of global tracks with hit in ITS layer 1
    ttree->SetBranchAddress("errFracTrackMI1",&errFracTrackMI1);
    ttree->SetBranchAddress("FracTrackMI2",&FracTrackMI2); // fraction of global tracks with hit in ITS layer 2
    ttree->SetBranchAddress("errFracTrackMI2",&errFracTrackMI2);
    ttree->SetBranchAddress("FracTrackMI3",&FracTrackMI3); // fraction of global tracks with hit in ITS layer 3
    ttree->SetBranchAddress("errFracTrackMI3",&errFracTrackMI3);
    ttree->SetBranchAddress("FracTrackMI4",&FracTrackMI4); // fraction of global tracks with hit in ITS layer 4
    ttree->SetBranchAddress("errFracTrackMI4",&errFracTrackMI4);
    ttree->SetBranchAddress("FracTrackMI5",&FracTrackMI5); // fraction of global tracks with hit in ITS layer 5
    ttree->SetBranchAddress("errFracTrackMI5",&errFracTrackMI5);
    ttree->SetBranchAddress("FracTrackMI6",&FracTrackMI6); // fraction of global tracks with hit in ITS layer 6
    ttree->SetBranchAddress("errFracTrackMI6",&errFracTrackMI6);
    ttree->SetBranchAddress("FracTrackSA1",&FracTrackSA1); // fraction of SA tracks with hit in ITS layer 1
    ttree->SetBranchAddress("errFracTrackSA1",&errFracTrackSA1);
    ttree->SetBranchAddress("FracTrackSA2",&FracTrackSA2); // fraction of SA tracks with hit in ITS layer 2
    ttree->SetBranchAddress("errFracTrackSA2",&errFracTrackSA2);
    ttree->SetBranchAddress("FracTrackSA3",&FracTrackSA3); // fraction of SA tracks with hit in ITS layer 3
    ttree->SetBranchAddress("errFracTrackSA3",&errFracTrackSA3);
    ttree->SetBranchAddress("FracTrackSA4",&FracTrackSA4); // fraction of SA tracks with hit in ITS layer 4
    ttree->SetBranchAddress("errFracTrackSA4",&errFracTrackSA4);
    ttree->SetBranchAddress("FracTrackSA5",&FracTrackSA5); // fraction of SA tracks with hit in ITS layer 5
    ttree->SetBranchAddress("errFracTrackSA5",&errFracTrackSA5);
    ttree->SetBranchAddress("FracTrackSA6",&FracTrackSA6); // fraction of SA tracks with hit in ITS layer 6
    ttree->SetBranchAddress("errFracTrackSA6",&errFracTrackSA6);

    ///// TOF-ITS (48 variables)
    Float_t Eff6Pt02TOF,errEff6Pt02TOF,Eff6Pt1TOF,errEff6Pt1TOF,Eff6Pt10TOF,errEff6Pt10TOF;
    Float_t Eff5Pt02TOF,errEff5Pt02TOF,Eff5Pt1TOF,errEff5Pt1TOF,Eff5Pt10TOF,errEff5Pt10TOF;
    Float_t Eff4Pt02TOF,errEff4Pt02TOF,Eff4Pt1TOF,errEff4Pt1TOF,Eff4Pt10TOF,errEff4Pt10TOF;
    Float_t Eff3Pt02TOF,errEff3Pt02TOF,Eff3Pt1TOF,errEff3Pt1TOF,Eff3Pt10TOF,errEff3Pt10TOF;
    Float_t Eff2Pt02TOF,errEff2Pt02TOF,Eff2Pt1TOF,errEff2Pt1TOF,Eff2Pt10TOF,errEff2Pt10TOF;
    Float_t EffSPDPt02TOF,errEffSPDPt02TOF,EffSPDPt1TOF,errEffSPDPt1TOF,EffSPDPt10TOF,errEffSPDPt10TOF;
    Float_t EffoneSPDPt02TOF,errEffoneSPDPt02TOF,EffoneSPDPt1TOF,errEffoneSPDPt1TOF,EffoneSPDPt10TOF,errEffoneSPDPt10TOF;
    Float_t EffTOTPt02TOF,errEffTOTPt02TOF,EffTOTPt1TOF,errEffTOTPt1TOF,EffTOTPt10TOF,errEffTOTPt10TOF;

    ttree->SetBranchAddress("Eff6Pt02TOF",&Eff6Pt02TOF); // matching efficiency low pt 6 clusters
    ttree->SetBranchAddress("errEff6Pt02TOF",&errEff6Pt02TOF); // error matching efficiency low pt 6 clusters
    ttree->SetBranchAddress("Eff5Pt02TOF",&Eff5Pt02TOF); // matching efficiency low pt 5 clusters
    ttree->SetBranchAddress("errEff5Pt02TOF",&errEff5Pt02TOF); // error matching efficiency low pt 5 clusters
    ttree->SetBranchAddress("Eff4Pt02TOF",&Eff4Pt02TOF); // matching efficiency low pt 4 clusters
    ttree->SetBranchAddress("errEff4Pt02TOF",&errEff4Pt02TOF); // error matching efficiency low pt 4 clusters
    ttree->SetBranchAddress("Eff3Pt02TOF",&Eff3Pt02TOF); // matching efficiency low pt 3 clusters
    ttree->SetBranchAddress("errEff3Pt02TOF",&errEff3Pt02TOF); // error matching efficiency low pt 3 clusters
    ttree->SetBranchAddress("Eff2Pt02TOF",&Eff2Pt02TOF); // matching efficiency low pt 2 clusters
    ttree->SetBranchAddress("errEff2Pt02TOF",&errEff2Pt02TOF); // error matching efficiency low pt 2 clusters
    ttree->SetBranchAddress("EffSPDPt02TOF",&EffSPDPt02TOF); // matching efficiency low pt 2 SPD
    ttree->SetBranchAddress("errEffSPDPt02TOF",&errEffSPDPt02TOF); // error matching efficiency low pt 2 SPD
    ttree->SetBranchAddress("EffoneSPDPt02TOF",&EffoneSPDPt02TOF); // matching efficiency low pt 6 one SPD
    ttree->SetBranchAddress("errEffoneSPDPt02TOF",&errEffoneSPDPt02TOF); // error matching efficiency low pt one SPD
    ttree->SetBranchAddress("EffTOTPt02TOF",&EffTOTPt02TOF); // matching efficiency low pt
    ttree->SetBranchAddress("errEffTOTPt02TOF",&errEffTOTPt02TOF); // error matching efficiency low pt
    
    ttree->SetBranchAddress("Eff6Pt1TOF",&Eff6Pt1TOF); // matching efficiency mid pt 6 clusters
    ttree->SetBranchAddress("errEff6Pt1TOF",&errEff6Pt1TOF); // error matching efficiency mid pt 6 clusters
    ttree->SetBranchAddress("Eff5Pt1TOF",&Eff5Pt1TOF); // matching efficiency mid pt 5 clusters
    ttree->SetBranchAddress("errEff5Pt1TOF",&errEff5Pt1TOF); // error matching efficiency mid pt 5 clusters
    ttree->SetBranchAddress("Eff4Pt1TOF",&Eff4Pt1TOF); // matching efficiency mid pt 4 clusters
    ttree->SetBranchAddress("errEff4Pt1TOF",&errEff4Pt1TOF); // error matching efficiency mid pt 4 clusters
    ttree->SetBranchAddress("Eff3Pt1TOF",&Eff3Pt1TOF); // matching efficiency mid pt 3 clusters
    ttree->SetBranchAddress("errEff3Pt1TOF",&errEff3Pt1TOF); // error matching efficiency mid pt 3 clusters
    ttree->SetBranchAddress("Eff2Pt1TOF",&Eff2Pt1TOF); // matching efficiency mid pt 2 clusters
    ttree->SetBranchAddress("errEff2Pt1TOF",&errEff2Pt1TOF); // error matching efficiency mid pt 2 clusters
    ttree->SetBranchAddress("EffSPDPt1TOF",&EffSPDPt1TOF); // matching efficiency mid pt 2 SPD
    ttree->SetBranchAddress("errEffSPDPt1TOF",&errEffSPDPt1TOF); // error matching efficiency mid pt 2 SPD
    ttree->SetBranchAddress("EffoneSPDPt1TOF",&EffoneSPDPt1TOF); // matching efficiency mid pt 6 one SPD
    ttree->SetBranchAddress("errEffoneSPDPt1TOF",&errEffoneSPDPt1TOF); // error matching efficiency mid pt one SPD
    ttree->SetBranchAddress("EffTOTPt1TOF",&EffTOTPt1TOF); // matching efficiency mid pt
    ttree->SetBranchAddress("errEffTOTPt1TOF",&errEffTOTPt1TOF); // error matching efficiency mid pt

    ttree->SetBranchAddress("Eff6Pt10TOF",&Eff6Pt10TOF); // matching efficiency high pt 6 clusters
    ttree->SetBranchAddress("errEff6Pt10TOF",&errEff6Pt10TOF); // error matching efficiency high pt 6 clusters
    ttree->SetBranchAddress("Eff5Pt10TOF",&Eff5Pt10TOF); // matching efficiency high pt 5 clusters
    ttree->SetBranchAddress("errEff5Pt10TOF",&errEff5Pt10TOF); // error matching efficiency high pt 5 clusters
    ttree->SetBranchAddress("Eff4Pt10TOF",&Eff4Pt10TOF); // matching efficiency high pt 4 clusters
    ttree->SetBranchAddress("errEff4Pt10TOF",&errEff4Pt10TOF); // error matching efficiency high pt 4 clusters
    ttree->SetBranchAddress("Eff3Pt10TOF",&Eff3Pt10TOF); // matching efficiency high pt 3 clusters
    ttree->SetBranchAddress("errEff3Pt10TOF",&errEff3Pt10TOF); // error matching efficiency high pt 3 clusters
    ttree->SetBranchAddress("Eff2Pt10TOF",&Eff2Pt10TOF); // matching efficiency high pt 2 clusters
    ttree->SetBranchAddress("errEff2Pt10TOF",&errEff2Pt10TOF); // error matching efficiency high pt 2 clusters
    ttree->SetBranchAddress("EffSPDPt10TOF",&EffSPDPt10TOF); // matching efficiency high pt 2 SPD
    ttree->SetBranchAddress("errEffSPDPt10TOF",&errEffSPDPt10TOF); // error matching efficiency high pt 2 SPD
    ttree->SetBranchAddress("EffoneSPDPt10TOF",&EffoneSPDPt10TOF); // matching efficiency high pt 6 one SPD
    ttree->SetBranchAddress("errEffoneSPDPt10TOF",&errEffoneSPDPt10TOF); // error matching efficiency high pt one SPD
    ttree->SetBranchAddress("EffTOTPt10TOF",&EffTOTPt10TOF); // matching efficiency high pt
    ttree->SetBranchAddress("errEffTOTPt10TOF",&errEffTOTPt10TOF); // error matching efficiency high pt
   
    
    //*********************************** ITSsa (32+12 arrays) *********************************************//
    
    Float_t NITSTPCPtBin0,NITSTPCPtBin1,NITSTPCPtBin2,NITSsaPtBin0,NITSsaPtBin1,NITSsaPtBin2,NITSpureSAPtBin0;
    Float_t NITSpureSAPtBin1,NITSpureSAPtBin2,ratioPtBin0,ratioPtBin1,ratioPtBin2,NcluITSpSA,errNcluITSpSA;
    Float_t dedx4_3,errdedx4_3,PtpionpSA,errPtpionpSA,NclupSA0,errNclupSA0,NclupSA1,errNclupSA1,NclupSA2;
    Float_t errNclupSA2,NclupSA3,errNclupSA3,NclupSA4,errNclupSA4,NclupSA5,errNclupSA5;
    Float_t chi2TPCITS,chi2ITSpureSA;
    Float_t occ_eta_1[2],occ_eta_2[2],occ_eta_3[2],occ_eta_4[2],occ_eta_5[2],occ_eta_6[2];
    Float_t occ_phi_1[40],occ_phi_2[40],occ_phi_3[40],occ_phi_4[40],occ_phi_5[40],occ_phi_6[40];
    
    ttree->SetBranchAddress("NITSTPCPtBin0",&NITSTPCPtBin0); // matching efficiency high pt
    ttree->SetBranchAddress("NITSTPCPtBin1",&NITSTPCPtBin1); // matching efficiency high pt
    ttree->SetBranchAddress("NITSTPCPtBin2",&NITSTPCPtBin2); // matching efficiency high pt
    ttree->SetBranchAddress("NITSsaPtBin0",&NITSsaPtBin0); // matching efficiency high pt
    ttree->SetBranchAddress("NITSsaPtBin1",&NITSsaPtBin1); // matching efficiency high pt
    ttree->SetBranchAddress("NITSsaPtBin2",&NITSsaPtBin2); // matching efficiency high pt
    ttree->SetBranchAddress("NITSpureSAPtBin0",&NITSpureSAPtBin0); // matching efficiency high pt
    ttree->SetBranchAddress("NITSpureSAPtBin1",&NITSpureSAPtBin1); // matching efficiency high pt
    ttree->SetBranchAddress("NITSpureSAPtBin2",&NITSpureSAPtBin2); // matching efficiency high pt
    ttree->SetBranchAddress("ratioPtBin0",&ratioPtBin0); // matching efficiency high pt
    ttree->SetBranchAddress("ratioPtBin1",&ratioPtBin1); // matching efficiency high pt
    ttree->SetBranchAddress("ratioPtBin2",&ratioPtBin2); // matching efficiency high pt
    ttree->SetBranchAddress("NcluITSpSA",&NcluITSpSA); // matching efficiency high pt
    ttree->SetBranchAddress("errNcluITSpSA",&errNcluITSpSA); // matching efficiency high pt
    ttree->SetBranchAddress("dedx4_3",&dedx4_3); // matching efficiency high pt
    ttree->SetBranchAddress("errdedx4_3",&errdedx4_3); // matching efficiency high pt
    ttree->SetBranchAddress("PtpionpSA",&PtpionpSA); // matching efficiency high pt
    ttree->SetBranchAddress("errPtpionpSA",&errPtpionpSA); // matching efficiency high pt
    ttree->SetBranchAddress("NclupSA0",&NclupSA0); // matching efficiency high pt
    ttree->SetBranchAddress("errNclupSA0",&errNclupSA0); // matching efficiency high pt
    ttree->SetBranchAddress("NclupSA1",&NclupSA1); // matching efficiency high pt
    ttree->SetBranchAddress("errNclupSA1",&errNclupSA1); // matching efficiency high pt
    ttree->SetBranchAddress("NclupSA2",&NclupSA2); // matching efficiency high pt
    ttree->SetBranchAddress("errNclupSA2",&errNclupSA2); // matching efficiency high pt
    ttree->SetBranchAddress("NclupSA3",&NclupSA3); // matching efficiency high pt
    ttree->SetBranchAddress("errNclupSA3",&errNclupSA3); // matching efficiency high pt
    ttree->SetBranchAddress("NclupSA4",&NclupSA4); // matching efficiency high pt
    ttree->SetBranchAddress("errNclupSA4",&errNclupSA4); // matching efficiency high pt
    ttree->SetBranchAddress("NclupSA5",&NclupSA5); // matching efficiency high pt
    ttree->SetBranchAddress("errNclupSA5",&errNclupSA5); // matching efficiency high pt

//    ttree->SetBranchAddress("chi2TPCITS",&chi2TPCITS); // mean chi2 for TPCITS tracks
//    ttree->SetBranchAddress("chi2ITSpureSA",&chi2ITSpureSA); // mean chi2 for ITSpureSA tracks
    ttree->SetBranchAddress("chi1TPCITS",&chi2TPCITS); // mean chi2 for TPCITS tracks
    ttree->SetBranchAddress("chi1ITSpureSA",&chi2ITSpureSA); // mean chi2 for ITSpureSA tracks

    ttree->SetBranchAddress("occ_eta_1",occ_eta_1); // eta occupancy on layer 1 for TPCITS tracks
    ttree->SetBranchAddress("occ_eta_2",occ_eta_2); // eta occupancy on layer 2 for TPCITS tracks
    ttree->SetBranchAddress("occ_eta_3",occ_eta_3); // eta occupancy on layer 3 for TPCITS tracks
    ttree->SetBranchAddress("occ_eta_4",occ_eta_4); // eta occupancy on layer 4 for TPCITS tracks
    ttree->SetBranchAddress("occ_eta_5",occ_eta_5); // eta occupancy on layer 5 for TPCITS tracks
    ttree->SetBranchAddress("occ_eta_6",occ_eta_6); // eta occupancy on layer 6 for TPCITS tracks
    
    ttree->SetBranchAddress("occ_phi_1",occ_phi_1); // phi occupancy on layer 1 for TPCITS tracks
    ttree->SetBranchAddress("occ_phi_2",occ_phi_2); // phi occupancy on layer 1 for TPCITS tracks
    ttree->SetBranchAddress("occ_phi_3",occ_phi_3); // phi occupancy on layer 1 for TPCITS tracks
    ttree->SetBranchAddress("occ_phi_4",occ_phi_4); // phi occupancy on layer 1 for TPCITS tracks
    ttree->SetBranchAddress("occ_phi_5",occ_phi_5); // phi occupancy on layer 1 for TPCITS tracks
    ttree->SetBranchAddress("occ_phi_6",occ_phi_6); // phi occupancy on layer 1 for TPCITS tracks

   
    //*********************************** pileup SPD (10) *********************************************//
    
    Float_t npilvtx,errnpilvtx,ntrklpil,errntrklpil,ntrklnopil,errntrklnopil,ncl1pil,errncl1pil,ncl1nopil,errncl1nopil;

    ttree->SetBranchAddress("npilvtx",&npilvtx); // matching efficiency high pt
    ttree->SetBranchAddress("errnpilvtx",&errnpilvtx); // matching efficiency high pt
    ttree->SetBranchAddress("ntrklpil",&ntrklpil); // matching efficiency high pt
    ttree->SetBranchAddress("errntrklpil",&errntrklpil); // matching efficiency high pt
    ttree->SetBranchAddress("ntrklnopil",&ntrklnopil); // matching efficiency high pt
    ttree->SetBranchAddress("errntrklnopil",&errntrklnopil); // matching efficiency high pt
    ttree->SetBranchAddress("ncl1pil",&ncl1pil); // matching efficiency high pt
    ttree->SetBranchAddress("errncl1pil",&errncl1pil); // matching efficiency high pt
    ttree->SetBranchAddress("ncl1nopil",&ncl1nopil); // matching efficiency high pt
    ttree->SetBranchAddress("errncl1nopil",&errncl1nopil); // matching efficiency high pt

    ttree->SetBranchAddress("npilvtx",&npilvtx); // number of pileup vtx/event
    ttree->SetBranchAddress("errnpilvtx",&errnpilvtx); // number of pileup vtx/event
    ttree->SetBranchAddress("ntrklpil",&ntrklpil); // number of tracklets of pileup vtx
    ttree->SetBranchAddress("errntrklpil",&errntrklpil); // number of tracklets of pileup vtx
    ttree->SetBranchAddress("ntrklnopil",&ntrklnopil); // number of tracklets of pileup vtx
    ttree->SetBranchAddress("errntrklnopil",&errntrklnopil); // number of tracklets of no-pileup vtx
    ttree->SetBranchAddress("ncl1pil",&ncl1pil); // number of clusters on layer 2 for pileup vtx
    ttree->SetBranchAddress("errncl1pil",&errncl1pil); // number of clusters on layer 2 for pileup vtx
    ttree->SetBranchAddress("ncl1nopil",&ncl1nopil); // number of clusters on layer 2 for pileup vtx
    ttree->SetBranchAddress("errncl1nopil",&errncl1nopil); // number of clusters on layer 2 for no-pileup vtx


    //*********************************** ITS PID TPCITS tracks (18) *********************************************//

    Float_t nsigmapi02,errnsigmapi02,nsigmapi05,errnsigmapi05,nsigmapi1,errnsigmapi1,nsigmapi3,errnsigmapi3;

    ttree->SetBranchAddress("nsigmapi02",&nsigmapi02); // nsigma for pions at 0.2 GeV/c
    ttree->SetBranchAddress("errnsigmapi02",&errnsigmapi02); // nsigma for pions at 0.2 GeV/c
    ttree->SetBranchAddress("nsigmapi05",&nsigmapi05); // nsigma for pions at 0.5 GeV/c
    ttree->SetBranchAddress("errnsigmapi05",&errnsigmapi05); // nsigma for pions at 0.5 GeV/c
    ttree->SetBranchAddress("nsigmapi1",&nsigmapi1); // nsigma for pions at 1.0 GeV/c
    ttree->SetBranchAddress("errnsigmapi1",&errnsigmapi1); // nsigma for pions at 1.0 GeV/c
    ttree->SetBranchAddress("nsigmapi3",&nsigmapi3); // nsigma for pions at 3 GeV/c
    ttree->SetBranchAddress("errnsigmapi3",&errnsigmapi3); // nsigma for pions at 3 GeV/c
    

    //*********************************** DCA TPCITS tracks (32) *********************************************//

    Float_t mdca05,errmdca05,rmsdca05,errrmsdca05,mdca1,errmdca1,rmsdca1,errrmsdca1,mdca5,errmdca5,rmsdca5,errrmsdca5,mdca10,errmdca10,rmsdca10,errrmsdca10;
    Float_t mdcaz05,errmdcaz05,rmsdcaz05,errrmsdcaz05,mdcaz1,errmdcaz1,rmsdcaz1,errrmsdcaz1,mdcaz5,errmdcaz5,rmsdcaz5,errrmsdcaz5,mdcaz10,errmdcaz10,rmsdcaz10,errrmsdcaz10;

    ttree->SetBranchAddress("mdca05",&mdca05); // mean DCA at 0.5 GeV/c
    ttree->SetBranchAddress("errmdca05",&errmdca05); // mean DCA at 0.5 GeV/c
    ttree->SetBranchAddress("rmsdca05",&rmsdca05); // rms DCA at 0.5 GeV/c
    ttree->SetBranchAddress("errrmsdca05",&errrmsdca05); // rms DCA at 0.5 GeV/c
    ttree->SetBranchAddress("mdca1",&mdca1); // mean DCA at 1.0 GeV/c
    ttree->SetBranchAddress("errmdca1",&errmdca1); // mean DCA at 1 GeV/c
    ttree->SetBranchAddress("rmsdca1",&rmsdca1); // rms DCA at 1 GeV/c
    ttree->SetBranchAddress("errrmsdca1",&errrmsdca1); // rms DCA at 1 GeV/c
    ttree->SetBranchAddress("mdca5",&mdca5); // mean DCA at 5 GeV/c
    ttree->SetBranchAddress("errmdca5",&errmdca5); // mean DCA at 5 GeV/c
    ttree->SetBranchAddress("rmsdca5",&rmsdca5); // rms DCA at 5 GeV/c
    ttree->SetBranchAddress("errrmsdca5",&errrmsdca5); // rms DCA at 5 GeV/c
    ttree->SetBranchAddress("mdca10",&mdca10); // mean DCA at 10 GeV/c
    ttree->SetBranchAddress("errmdca10",&errmdca10); // mean DCA at 10 GeV/c
    ttree->SetBranchAddress("rmsdca10",&rmsdca10); // rms DCA at 10 GeV/c
    ttree->SetBranchAddress("errrmsdca10",&errrmsdca10); // rms DCA at 010 GeV/c

    ttree->SetBranchAddress("mdcaz05",&mdcaz05); // mean DCAz at 0.5 GeV/c
    ttree->SetBranchAddress("errmdcaz05",&errmdcaz05); // mean DCAz at 0.5 GeV/c
    ttree->SetBranchAddress("rmsdcaz05",&rmsdcaz05); // rms DCAz at 0.5 GeV/c
    ttree->SetBranchAddress("errrmsdcaz05",&errrmsdcaz05); // rms DCAz at 0.5 GeV/c
    ttree->SetBranchAddress("mdcaz1",&mdcaz1); // mean DCAz at 1.0 GeV/c
    ttree->SetBranchAddress("errmdcaz1",&errmdcaz1); // mean DCAz at 1 GeV/c
    ttree->SetBranchAddress("rmsdcaz1",&rmsdcaz1); // rms DCAz at 1 GeV/c
    ttree->SetBranchAddress("errrmsdcaz1",&errrmsdcaz1); // rms DCAz at 1 GeV/c
    ttree->SetBranchAddress("mdcaz5",&mdcaz5); // mean DCAz at 5 GeV/c
    ttree->SetBranchAddress("errmdcaz5",&errmdcaz5); // mean DCAz at 5 GeV/c
    ttree->SetBranchAddress("rmsdcaz5",&rmsdcaz5); // rms DCAz at 5 GeV/c
    ttree->SetBranchAddress("errrmsdcaz5",&errrmsdcaz5); // rms DCAz at 5 GeV/c
    ttree->SetBranchAddress("mdcaz10",&mdcaz10); // mean DCAz at 10 GeV/c
    ttree->SetBranchAddress("errmdcaz10",&errmdcaz10); // mean DCAz at 10 GeV/c
    ttree->SetBranchAddress("rmsdcaz10",&rmsdcaz10); // rms DCAz at 10 GeV/c
    ttree->SetBranchAddress("errrmsdcaz10",&errrmsdcaz10); // rms DCAz at 010 GeV/c
   
        //********************************************************************************************************//
    
  Int_t nRuns=ttree->GetEntries();
  TList lista;


  //****************************** SDD  PLOTS  (18) *****************************************************//
  
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
    histoTrackClu1->SetLineColor(kGray+1);
    histoTrackClu1->SetMarkerColor(kGray+1);
    histoTrackClu1->SetMarkerStyle(24);
  TH1F* histoTrackClu2=new TH1F("histoTrackClu2","",nRuns,0.,nRuns);
    histoTrackClu2->SetLineColor(kGray+2);
    histoTrackClu2->SetMarkerColor(kGray+2);
    histoTrackClu2->SetMarkerStyle(26);
  TH1F* histoTrackClu3=new TH1F("histoTrackClu3","",nRuns,0.,nRuns);
    histoTrackClu3->SetLineColor(1);
    histoTrackClu3->SetMarkerStyle(20);
    histoTrackClu3->SetMinimum(0.);
    histoTrackClu3->SetMaximum(1.05);
  TH1F* histoTrackClu4=new TH1F("histoTrackClu4","",nRuns,0.,nRuns);
    histoTrackClu4->SetLineColor(2);
    histoTrackClu4->SetMarkerColor(2);
    histoTrackClu4->SetMarkerStyle(22);
  TH1F* histoTrackClu5=new TH1F("histoTrackClu5","",nRuns,0.,nRuns);
    histoTrackClu5->SetLineColor(4);
    histoTrackClu5->SetMarkerColor(4);
    histoTrackClu5->SetMarkerStyle(29);
  TH1F* histoTrackClu6=new TH1F("histoTrackClu6","",nRuns,0.,nRuns);
    histoTrackClu6->SetLineColor(kBlue+1);
    histoTrackClu6->SetMarkerColor(kBlue+1);
    histoTrackClu6->SetMarkerStyle(30);
//  TH1F* histoNmodEmpty=new TH1F("histoNmodEmpty","",nRuns,0.,nRuns);
//
  TH1F* histoFracDead3=new TH1F("histoFracDead3","",nRuns,0.,nRuns);
  TH1F* histoFracDead4=new TH1F("histoFracDead4","",nRuns,0.,nRuns);
    TH1F* histoEvwSDD=new TH1F("histoEvwSDD","",nRuns,0.,nRuns);

    TH1F* histoFlagON3=new TH1F("histoFlagON3","",nRuns,0.,nRuns);
    TH1F* histoFlagON4=new TH1F("histoFlagON4","",nRuns,0.,nRuns);
    TH1F* histoFlagMPV3=new TH1F("histoFlagMPV3","",nRuns,0.,nRuns);
    TH1F* histoFlagMPV4=new TH1F("histoFlagMPV4","",nRuns,0.,nRuns);
    TH1F* histoFlagminT=new TH1F("histoFlagminT","",nRuns,0.,nRuns);
    TH1F* histoFlagmeanT=new TH1F("histoFlagmeanT","",nRuns,0.,nRuns);

  //****************************** VERTEX  PLOTS  (13) *****************************************************//

  TH1F *hMeanVx = new TH1F("hMeanVx","Track Vertex Vx Distribution",nRuns,0.,nRuns);
    hMeanVx->SetLineWidth(2);
    hMeanVx->SetLineColor(kBlue+2);
    hMeanVx->SetMarkerColor(kBlue+2);
    hMeanVx->SetMarkerStyle(20);

  TH1F *hMeanVy = new TH1F("hMeanVy","Track Vertex Vy Distribution",nRuns,0.,nRuns);
    hMeanVy->SetLineWidth(2);
    hMeanVy->SetLineColor(kBlue+2);
    hMeanVy->SetMarkerColor(kBlue+2);
    hMeanVy->SetMarkerStyle(20);

  TH1F *hMeanVz = new TH1F("hMeanVz","Track Vertex Vz Distribution",nRuns,0.,nRuns);
    hMeanVz->SetLineWidth(2);
    hMeanVz->SetLineColor(kBlue+2);
    hMeanVz->SetMarkerColor(kBlue+2);
    hMeanVz->SetMarkerStyle(20);

  TH1F *hSigmaVx = new TH1F("hSigmaVx","Track Vertex SigmaVx Distribution",nRuns,0.,nRuns);
    hSigmaVx->SetLineWidth(2);
    hSigmaVx->SetLineColor(kBlue+2);
    hSigmaVx->SetMarkerColor(kBlue+2);
    hSigmaVx->SetMarkerStyle(20);

    TH1F *hSigmaVy = new TH1F("hSigmaVy","Track Vertex SigmaVy Distribution",nRuns,0.,nRuns);
    hSigmaVy->SetLineWidth(2);
    hSigmaVy->SetLineColor(kBlue+2);
    hSigmaVy->SetMarkerColor(kBlue+2);
    hSigmaVy->SetMarkerStyle(20);

  TH1F *hSigmaVz = new TH1F("hSigmaVz","Track Vertex SigmaVz Distribution",nRuns,0.,nRuns);
    hSigmaVz->SetLineWidth(2);
    hSigmaVz->SetLineColor(kBlue+2);
    hSigmaVz->SetMarkerColor(kBlue+2);
    hSigmaVz->SetMarkerStyle(20);

TH1F *hMeanVxSPD = new TH1F("hMeanVxSPD","Track Vertex Vx Distribution",nRuns,0.,nRuns);
    hMeanVxSPD->SetLineWidth(2);
    hMeanVxSPD->SetLineColor(2);
    hMeanVxSPD->SetMarkerColor(2);
    hMeanVxSPD->SetMarkerStyle(20);
    
  TH1F *hMeanVySPD = new TH1F("hMeanVySPD","Track Vertex Vy Distribution",nRuns,0.,nRuns);
    hMeanVySPD->SetLineWidth(2);
    hMeanVySPD->SetLineColor(2);
    hMeanVySPD->SetMarkerColor(2);
    hMeanVySPD->SetMarkerStyle(20);
    
  TH1F *hMeanVzSPD = new TH1F("hMeanVzSPD","Track Vertex Vz Distribution",nRuns,0.,nRuns);
    hMeanVzSPD->SetLineWidth(2);
    hMeanVzSPD->SetLineColor(2);
    hMeanVzSPD->SetMarkerColor(2);
    hMeanVzSPD->SetMarkerStyle(20);
    
  TH1F *hSigmaVxSPD = new TH1F("hSigmaVxSPD","Track Vertex SigmaVx Distribution",nRuns,0.,nRuns);
    hSigmaVxSPD->SetLineWidth(2);
    hSigmaVxSPD->SetLineColor(2);
    hSigmaVxSPD->SetMarkerColor(2);
    hSigmaVxSPD->SetMarkerStyle(20);
  TH1F *hSigmaVySPD = new TH1F("hSigmaVySPD","Track Vertex SigmaVy Distribution",nRuns,0.,nRuns);
    hSigmaVySPD->SetLineWidth(2);
    hSigmaVySPD->SetLineColor(2);
    hSigmaVySPD->SetMarkerColor(2);
    hSigmaVySPD->SetMarkerStyle(20);
  TH1F *hSigmaVzSPD = new TH1F("hSigmaVzSPD","Track Vertex SigmaVz Distribution",nRuns,0.,nRuns);
    hSigmaVzSPD->SetLineWidth(2);
    hSigmaVzSPD->SetLineColor(2);
    hSigmaVzSPD->SetMarkerColor(2);
    hSigmaVzSPD->SetMarkerStyle(20);

TH1F *hMeanVxOCDB = new TH1F("hMeanVxOCDB","Diamond OCDB Vx",nRuns,0.,nRuns);
    hMeanVxOCDB->SetLineWidth(2);
    hMeanVxOCDB->SetLineColor(kMagenta+2);
    hMeanVxOCDB->SetMarkerColor(kMagenta+2);
    hMeanVxOCDB->SetMarkerStyle(25);
    
  TH1F *hMeanVyOCDB = new TH1F("hMeanVyOCDB","Diamond OCDB Vy",nRuns,0.,nRuns);
    hMeanVyOCDB->SetLineWidth(2);
    hMeanVyOCDB->SetLineColor(kMagenta+2);
    hMeanVyOCDB->SetMarkerColor(kMagenta+2);
    hMeanVyOCDB->SetMarkerStyle(25);
    
  TH1F *hMeanVzOCDB = new TH1F("hMeanVzOCDB","Diamond OCDB Vz",nRuns,0.,nRuns);
    hMeanVzOCDB->SetLineWidth(2);
    hMeanVzOCDB->SetLineColor(kMagenta+2);
    hMeanVzOCDB->SetMarkerColor(kMagenta+2);
    hMeanVzOCDB->SetMarkerStyle(25);
    
  TH1F *hSigmaVxOCDB = new TH1F("hSigmaVxOCDB","Diamond OCDB SigmaVx",nRuns,0.,nRuns);
    hSigmaVxOCDB->SetLineWidth(2);
    hSigmaVxOCDB->SetLineColor(kMagenta+2);
    hSigmaVxOCDB->SetMarkerColor(kMagenta+2);
    hSigmaVxOCDB->SetMarkerStyle(25);
  TH1F *hSigmaVyOCDB = new TH1F("hSigmaVyOCDB","Diamond OCDB SigmaVy",nRuns,0.,nRuns);
    hSigmaVyOCDB->SetLineWidth(2);
    hSigmaVyOCDB->SetLineColor(kMagenta+2);
    hSigmaVyOCDB->SetMarkerColor(kMagenta+2);
    hSigmaVyOCDB->SetMarkerStyle(25);
  TH1F *hSigmaVzOCDB = new TH1F("hSigmaVzOCDB","Diamond OCDB SigmaVz",nRuns,0.,nRuns);
    hSigmaVzOCDB->SetLineWidth(2);
    hSigmaVzOCDB->SetLineColor(kMagenta+2);
    hSigmaVzOCDB->SetMarkerColor(kMagenta+2);
    hSigmaVzOCDB->SetMarkerStyle(25);
 

  TH1F *hpileupSPD = new TH1F("hpileupSPD","Fraction of tracks with SPD pileup",nRuns,0.,nRuns);
        hpileupSPD->SetMarkerStyle(20);
  

  //****************************** SSD  PLOTS (16) *****************************************************//
  TH1F* histodEdxLay5 = new TH1F("histodEdxLay5","",nRuns,0.,nRuns);
  TH1F* histodEdxLay6 = new TH1F("histodEdxLay6","",nRuns,0.,nRuns);
  TH1F* histoChargeRatioLay5 = new TH1F("histoChargeRatioLay5","",nRuns,0.,nRuns);
  TH1F* histoChargeRatioLay6 = new TH1F("histoChargeRatioLay6","",nRuns,0.,nRuns);
    
    TH1F* histoFracBadn5=new TH1F("histoFracBadn5","",nRuns,0.,nRuns);
    TH1F* histoFracBadp5=new TH1F("histoFracBadp5","",nRuns,0.,nRuns);
    TH1F* histoFracBadn6=new TH1F("histoFracBadn6","",nRuns,0.,nRuns);
    TH1F* histoFracBadp6=new TH1F("histoFracBadp6","",nRuns,0.,nRuns);
    
    TH1F* histoFlagON5n=new TH1F("histoFlagON5n","",nRuns,0.,nRuns);
    TH1F* histoFlagON5p=new TH1F("histoFlagON5p","",nRuns,0.,nRuns);
    TH1F* histoFlagON6n=new TH1F("histoFlagON6n","",nRuns,0.,nRuns);
    TH1F* histoFlagON6p=new TH1F("histoFlagON6p","",nRuns,0.,nRuns);
    TH1F* histoFlagMPV5=new TH1F("histoFlagMPV5","",nRuns,0.,nRuns);
    TH1F* histoFlagMPV6=new TH1F("histoFlagMPV6","",nRuns,0.,nRuns);
    TH1F* histoFlagCR5=new TH1F("histoFlagCR5","",nRuns,0.,nRuns);
    TH1F* histoFlagCR6=new TH1F("histoFlagCR6","",nRuns,0.,nRuns);

    //************************************************ SPD PLOTS  *********************************************//
    TH1F* histoFracSPD1=new TH1F("histoFracSPD1","SPD inner; nrun number; Fraction of HSs",nRuns,0.,nRuns);
    histoFracSPD1->SetLineColor(kGreen+2);
    histoFracSPD1->SetMarkerColor(kGreen+2);
    histoFracSPD1->SetMarkerStyle(20);
    TH1F* histoFracSPD2=new TH1F("histoFracSPD2","SPD outer; nrun number; Fraction of HSs",nRuns,0.,nRuns);
    histoFracSPD2->SetLineColor(kYellow+2);
    histoFracSPD2->SetMarkerColor(kYellow+2);
    histoFracSPD2->SetMarkerStyle(20);
    TH1F* histoFlagON1=new TH1F("histoFlagON1","SPD inner; nrun number; Flag on HSs ON",nRuns,0.,nRuns);
    histoFlagON1->SetLineColor(kGreen+2);
    histoFlagON1->SetMarkerColor(kGreen+2);
    histoFlagON1->SetMarkerStyle(20);
    TH1F* histoFlagON2=new TH1F("histoFlagON2","SPD outer; nrun number; Flag on HSs ON",nRuns,0.,nRuns);
    histoFlagON2->SetLineColor(kYellow+2);
    histoFlagON2->SetMarkerColor(kYellow+2);
    histoFlagON2->SetMarkerStyle(20);
    
    //************************************************ MATCHING PLOTS: TPCITS *********************************************//

  TH1F *hEffSPDPt02 = new TH1F("hEffSPDPt02","Efficiency; nrun number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEffSPDPt02->SetLineWidth(2);
  hEffSPDPt02->SetLineColor(kAzure+1);
  hEffSPDPt02->SetMarkerColor(kAzure+1);
  hEffSPDPt02->SetMarkerStyle(20);
  TH1F *hEffSPDPt1 = new TH1F("hEffSPDPt1","Efficiency; nrun number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEffSPDPt1->SetLineWidth(2);
  hEffSPDPt1->SetLineColor(kAzure+1);
  hEffSPDPt1->SetMarkerColor(kAzure+1);
  hEffSPDPt1->SetMarkerStyle(20);
  TH1F *hEffSPDPt10 = new TH1F("hEffSPDPt10","Efficiency; nrun number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEffSPDPt10->SetLineWidth(2);
  hEffSPDPt10->SetLineColor(kAzure+1);
  hEffSPDPt10->SetMarkerColor(kAzure+1);
  hEffSPDPt10->SetMarkerStyle(20);

  TH1F *hEffoneSPDPt02 = new TH1F("hEffoneSPDPt02","Efficiency; nrun number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEffoneSPDPt02->SetLineWidth(2);
  hEffoneSPDPt02->SetLineColor(kGray);
  hEffoneSPDPt02->SetMarkerColor(kGray);
  hEffoneSPDPt02->SetMarkerStyle(20);
  TH1F *hEffoneSPDPt1 = new TH1F("hEffoneSPDPt1","Efficiency; nrun number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEffoneSPDPt1->SetLineWidth(2);
  hEffoneSPDPt1->SetLineColor(kGray);
  hEffoneSPDPt1->SetMarkerColor(kGray);
  hEffoneSPDPt1->SetMarkerStyle(20);
  TH1F *hEffoneSPDPt10 = new TH1F("hEffoneSPDPt10","Efficiency; nrun number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEffoneSPDPt10->SetLineWidth(2);
  hEffoneSPDPt10->SetLineColor(kGray);
  hEffoneSPDPt10->SetMarkerColor(kGray);
  hEffoneSPDPt10->SetMarkerStyle(20);

  TH1F *hEff2Pt02 = new TH1F("hEff2Pt02","Efficiency; nrun number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEff2Pt02->SetLineWidth(2);
  hEff2Pt02->SetLineColor(kViolet);
  hEff2Pt02->SetMarkerColor(kViolet);
  hEff2Pt02->SetMarkerStyle(20);
  TH1F *hEff2Pt1 = new TH1F("hEff2Pt1","Efficiency; nrun number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEff2Pt1->SetLineWidth(2);
  hEff2Pt1->SetLineColor(kViolet);
  hEff2Pt1->SetMarkerColor(kViolet);
  hEff2Pt1->SetMarkerStyle(20);
  TH1F *hEff2Pt10 = new TH1F("hEff2Pt10","Efficiency; nrun number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEff2Pt10->SetLineWidth(2);
  hEff2Pt10->SetLineColor(kViolet);
  hEff2Pt10->SetMarkerColor(kViolet);
  hEff2Pt10->SetMarkerStyle(20);

  TH1F *hEff3Pt02 = new TH1F("hEff3Pt02","Efficiency; nrun number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEff3Pt02->SetLineWidth(2);
  hEff3Pt02->SetLineColor(6);
  hEff3Pt02->SetMarkerColor(6);
  hEff3Pt02->SetMarkerStyle(20);
  TH1F *hEff3Pt1 = new TH1F("hEff3Pt1","Efficiency; nrun number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEff3Pt1->SetLineWidth(2);
  hEff3Pt1->SetLineColor(6);
  hEff3Pt1->SetMarkerColor(6);
  hEff3Pt1->SetMarkerStyle(20);
  TH1F *hEff3Pt10 = new TH1F("hEff3Pt10","Efficiency; nrun number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEff3Pt10->SetLineWidth(2);
  hEff3Pt10->SetLineColor(6);
  hEff3Pt10->SetMarkerColor(6);
  hEff3Pt10->SetMarkerStyle(20);

  TH1F *hEff4Pt02 = new TH1F("hEff4Pt02","Efficiency; nrun number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEff4Pt02->SetLineWidth(2);
  hEff4Pt02->SetLineColor(4);
  hEff4Pt02->SetMarkerColor(4);
  hEff4Pt02->SetMarkerStyle(20);
  TH1F *hEff4Pt1 = new TH1F("hEff4Pt1","Efficiency; nrun number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEff4Pt1->SetLineWidth(2);
  hEff4Pt1->SetLineColor(4);
  hEff4Pt1->SetMarkerColor(4);
  hEff4Pt1->SetMarkerStyle(20);
  TH1F *hEff4Pt10 = new TH1F("hEff4Pt10","Efficiency; nrun number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEff4Pt10->SetLineWidth(2);
  hEff4Pt10->SetLineColor(4);
  hEff4Pt10->SetMarkerColor(4);
  hEff4Pt10->SetMarkerStyle(20);

  TH1F *hEff5Pt02 = new TH1F("hEff5Pt02","Efficiency; nrun number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEff5Pt02->SetLineWidth(2);
  hEff5Pt02->SetLineColor(3);
  hEff5Pt02->SetMarkerColor(3);
  hEff5Pt02->SetMarkerStyle(20);
  TH1F *hEff5Pt1 = new TH1F("hEff5Pt1","Efficiency; nrun number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEff5Pt1->SetLineWidth(2);
  hEff5Pt1->SetLineColor(3);
  hEff5Pt1->SetMarkerColor(3);
  hEff5Pt1->SetMarkerStyle(20);
  TH1F *hEff5Pt10 = new TH1F("hEff5Pt10","Efficiency; nrun number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEff5Pt10->SetLineWidth(3);
  hEff5Pt10->SetLineColor(3);
  hEff5Pt10->SetMarkerColor(3);
  hEff5Pt10->SetMarkerStyle(20);

  TH1F *hEff6Pt02 = new TH1F("hEff6Pt02","Efficiency; nrun number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEff6Pt02->SetLineWidth(2);
  hEff6Pt02->SetLineColor(2);
  hEff6Pt02->SetMarkerColor(2);
  hEff6Pt02->SetMarkerStyle(20);
  TH1F *hEff6Pt1 = new TH1F("hEff6Pt1","Efficiency; nrun number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEff6Pt1->SetLineWidth(2);
  hEff6Pt1->SetLineColor(2);
  hEff6Pt1->SetMarkerColor(2);
  hEff6Pt1->SetMarkerStyle(20);
  TH1F *hEff6Pt10 = new TH1F("hEff6Pt10","Efficiency; nrun number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEff6Pt10->SetLineWidth(2);
  hEff6Pt10->SetLineColor(2);
  hEff6Pt10->SetMarkerColor(2);
  hEff6Pt10->SetMarkerStyle(20);

  TH1F *hEffTOTPt02 = new TH1F("hEffTOTPt02","Efficiency; nrun number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEffTOTPt02->SetLineWidth(2);
  hEffTOTPt02->SetLineColor(kBlue+2);
  hEffTOTPt02->SetMarkerColor(kBlue+2);
  hEffTOTPt02->SetMarkerStyle(20);
  TH1F *hEffTOTPt1 = new TH1F("hEffTOTPt1","Efficiency; nrun number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEffTOTPt1->SetLineWidth(2);
  hEffTOTPt1->SetLineColor(kBlue+2);
  hEffTOTPt1->SetMarkerColor(kBlue+2);
  hEffTOTPt1->SetMarkerStyle(20);
  TH1F *hEffTOTPt10 = new TH1F("hEffTOTPt10","Efficiency; nrun number; TPC+ITS / TPC",nRuns,0.,nRuns);
  hEffTOTPt10->SetLineWidth(2);
  hEffTOTPt10->SetLineColor(kBlue+2);
  hEffTOTPt10->SetMarkerColor(kBlue+2);
  hEffTOTPt10->SetMarkerStyle(20);
                            
        TH1F* histoTrackMI1=new TH1F("histoTrackMI1","",nRuns,0.,nRuns);
        TH1F* histoTrackMI2=new TH1F("histoTrackMI2","",nRuns,0.,nRuns);
        TH1F* histoTrackMI3=new TH1F("histoTrackMI3","",nRuns,0.,nRuns);
        TH1F* histoTrackMI4=new TH1F("histoTrackMI4","",nRuns,0.,nRuns);
        TH1F* histoTrackMI5=new TH1F("histoTrackMI5","",nRuns,0.,nRuns);
        TH1F* histoTrackMI6=new TH1F("histoTrackMI6","",nRuns,0.,nRuns);
        TH1F* histoTrackSA1=new TH1F("histoTrackSA1","",nRuns,0.,nRuns);
        TH1F* histoTrackSA2=new TH1F("histoTrackSA2","",nRuns,0.,nRuns);
        TH1F* histoTrackSA3=new TH1F("histoTrackSA3","",nRuns,0.,nRuns);
        TH1F* histoTrackSA4=new TH1F("histoTrackSA4","",nRuns,0.,nRuns);
        TH1F* histoTrackSA5=new TH1F("histoTrackSA5","",nRuns,0.,nRuns);
        TH1F* histoTrackSA6=new TH1F("histoTrackSA6","",nRuns,0.,nRuns);

//************************************************ MATCHING PLOTS: TOFITS *********************************************//

        TH1F *hEffTOFSPDPt02 = new TH1F("hEffTOFSPDPt02","Efficiency - P_{T} = 0.5; nrun number; TPC+ITS / TPC",nRuns,0.,nRuns);
        hEffTOFSPDPt02->SetLineWidth(2);
        hEffTOFSPDPt02->SetLineColor(kAzure+1);
        hEffTOFSPDPt02->SetMarkerColor(kAzure+1);
        hEffTOFSPDPt02->SetMarkerStyle(20);
        TH1F *hEffTOFSPDPt1 = new TH1F("hEffTOFSPDPt1","Efficiency; nrun number; TOF+ITS / TPC",nRuns,0.,nRuns);
        hEffTOFSPDPt1->SetLineWidth(2);
        hEffTOFSPDPt1->SetLineColor(kAzure+1);
        hEffTOFSPDPt1->SetMarkerColor(kAzure+1);
        hEffTOFSPDPt1->SetMarkerStyle(20);
        TH1F *hEffTOFSPDPt10 = new TH1F("hEffTOFSPDPt10","Efficiency; nrun number; TOF+ITS / TPC",nRuns,0.,nRuns);
        hEffTOFSPDPt10->SetLineWidth(2);
        hEffTOFSPDPt10->SetLineColor(kAzure+1);
        hEffTOFSPDPt10->SetMarkerColor(kAzure+1);
        hEffTOFSPDPt10->SetMarkerStyle(20);
                            
        TH1F *hEffTOFoneSPDPt02 = new TH1F("hEffTOFoneSPDPt02","Efficiency; nrun number; TOF+ITS / TPC",nRuns,0.,nRuns);
        hEffTOFoneSPDPt02->SetLineWidth(2);
        hEffTOFoneSPDPt02->SetLineColor(kGray);
        hEffTOFoneSPDPt02->SetMarkerColor(kGray);
        hEffTOFoneSPDPt02->SetMarkerStyle(20);
        TH1F *hEffTOFoneSPDPt1 = new TH1F("hEffTOFoneSPDPt1","Efficiency; nrun number; TOF+ITS / TPC",nRuns,0.,nRuns);
        hEffTOFoneSPDPt1->SetLineWidth(2);
        hEffTOFoneSPDPt1->SetLineColor(kGray);
        hEffTOFoneSPDPt1->SetMarkerColor(kGray);
        hEffTOFoneSPDPt1->SetMarkerStyle(20);
        TH1F *hEffTOFoneSPDPt10 = new TH1F("hEffTOFoneSPDPt10","Efficiency; nrun number; TOF+ITS / TPC",nRuns,0.,nRuns);
        hEffTOFoneSPDPt10->SetLineWidth(2);
        hEffTOFoneSPDPt10->SetLineColor(kGray);
        hEffTOFoneSPDPt10->SetMarkerColor(kGray);
        hEffTOFoneSPDPt10->SetMarkerStyle(20);
                            
        TH1F *hEffTOF2Pt02 = new TH1F("hEffTOF2Pt02","Efficiency; nrun number; TOF+ITS / TPC",nRuns,0.,nRuns);
        hEffTOF2Pt02->SetLineWidth(2);
        hEffTOF2Pt02->SetLineColor(kViolet);
        hEffTOF2Pt02->SetMarkerColor(kViolet);
        hEffTOF2Pt02->SetMarkerStyle(20);
        TH1F *hEffTOF2Pt1 = new TH1F("hEffTOF2Pt1","Efficiency; nrun number; TOF+ITS / TPC",nRuns,0.,nRuns);
        hEffTOF2Pt1->SetLineWidth(2);
        hEffTOF2Pt1->SetLineColor(kViolet);
        hEffTOF2Pt1->SetMarkerColor(kViolet);
        hEffTOF2Pt1->SetMarkerStyle(20);
        TH1F *hEffTOF2Pt10 = new TH1F("hEffTOF2Pt10","Efficiency; nrun number; TOF+ITS / TPC",nRuns,0.,nRuns);
        hEffTOF2Pt10->SetLineWidth(2);
        hEffTOF2Pt10->SetLineColor(kViolet);
        hEffTOF2Pt10->SetMarkerColor(kViolet);
        hEffTOF2Pt10->SetMarkerStyle(20);
                            
        TH1F *hEffTOF3Pt02 = new TH1F("hEffTOF3Pt02","Efficiency; nrun number; TOF+ITS / TPC",nRuns,0.,nRuns);
        hEffTOF3Pt02->SetLineWidth(2);
        hEffTOF3Pt02->SetLineColor(6);
        hEffTOF3Pt02->SetMarkerColor(6);
        hEffTOF3Pt02->SetMarkerStyle(20);
        TH1F *hEffTOF3Pt1 = new TH1F("hEffTOF3Pt1","Efficiency; nrun number; TOF+ITS / TPC",nRuns,0.,nRuns);
        hEffTOF3Pt1->SetLineWidth(2);
        hEffTOF3Pt1->SetLineColor(6);
        hEffTOF3Pt1->SetMarkerColor(6);
        hEffTOF3Pt1->SetMarkerStyle(20);
        TH1F *hEffTOF3Pt10 = new TH1F("hEffTOF3Pt10","Efficiency; nrun number; TOF+ITS / TPC",nRuns,0.,nRuns);
        hEffTOF3Pt10->SetLineWidth(2);
        hEffTOF3Pt10->SetLineColor(6);
        hEffTOF3Pt10->SetMarkerColor(6);
        hEffTOF3Pt10->SetMarkerStyle(20);
    
        TH1F *hEffTOF4Pt02 = new TH1F("hEffTOF4Pt02","Efficiency; nrun number; TOF+ITS / TPC",nRuns,0.,nRuns);
        hEffTOF4Pt02->SetLineWidth(2);
        hEffTOF4Pt02->SetLineColor(4);
        hEffTOF4Pt02->SetMarkerColor(4);
        hEffTOF4Pt02->SetMarkerStyle(20);
        TH1F *hEffTOF4Pt1 = new TH1F("hEffTOF4Pt1","Efficiency; nrun number; TOF+ITS / TPC",nRuns,0.,nRuns);
        hEffTOF4Pt1->SetLineWidth(2);
        hEffTOF4Pt1->SetLineColor(4);
        hEffTOF4Pt1->SetMarkerColor(4);
        hEffTOF4Pt1->SetMarkerStyle(20);
        TH1F *hEffTOF4Pt10 = new TH1F("hEffTOF4Pt10","Efficiency; nrun number; TOF+ITS / TPC",nRuns,0.,nRuns);
        hEffTOF4Pt10->SetLineWidth(2);
        hEffTOF4Pt10->SetLineColor(4);
        hEffTOF4Pt10->SetMarkerColor(4);
        hEffTOF4Pt10->SetMarkerStyle(20);
                            
        TH1F *hEffTOF5Pt02 = new TH1F("hEffTOF5Pt02","Efficiency; nrun number; TOF+ITS / TPC",nRuns,0.,nRuns);
        hEffTOF5Pt02->SetLineWidth(2);
        hEffTOF5Pt02->SetLineColor(3);
        hEffTOF5Pt02->SetMarkerColor(3);
        hEffTOF5Pt02->SetMarkerStyle(20);
        TH1F *hEffTOF5Pt1 = new TH1F("hEffTOF5Pt1","Efficiency; nrun number; TOF+ITS / TPC",nRuns,0.,nRuns);
        hEffTOF5Pt1->SetLineWidth(2);
        hEffTOF5Pt1->SetLineColor(3);
        hEffTOF5Pt1->SetMarkerColor(3);
        hEffTOF5Pt1->SetMarkerStyle(20);
         TH1F *hEffTOF5Pt10 = new TH1F("hEffTOF5Pt10","Efficiency; nrun number; TOF+ITS / TPC",nRuns,0.,nRuns);
        hEffTOF5Pt10->SetLineWidth(3);
        hEffTOF5Pt10->SetLineColor(3);
        hEffTOF5Pt10->SetMarkerColor(3);
        hEffTOF5Pt10->SetMarkerStyle(20);
                            
        TH1F *hEffTOF6Pt02 = new TH1F("hEffTOF6Pt02","Efficiency; nrun number; TOF+ITS / TPC",nRuns,0.,nRuns);
        hEffTOF6Pt02->SetLineWidth(2);
        hEffTOF6Pt02->SetLineColor(2);
        hEffTOF6Pt02->SetMarkerColor(2);
        hEffTOF6Pt02->SetMarkerStyle(20);
        TH1F *hEffTOF6Pt1 = new TH1F("hEffTOF6Pt1","Efficiency; nrun number; TOF+ITS / TPC",nRuns,0.,nRuns);
        hEffTOF6Pt1->SetLineWidth(2);
        hEffTOF6Pt1->SetLineColor(2);
        hEffTOF6Pt1->SetMarkerColor(2);
        hEffTOF6Pt1->SetMarkerStyle(20);
        TH1F *hEffTOF6Pt10 = new TH1F("hEffTOF6Pt10","Efficiency; nrun number; TOF+ITS / TPC",nRuns,0.,nRuns);
        hEffTOF6Pt10->SetLineWidth(2);
        hEffTOF6Pt10->SetLineColor(2);
        hEffTOF6Pt10->SetMarkerColor(2);
        hEffTOF6Pt10->SetMarkerStyle(20);
                            
        TH1F *hEffTOFTOTPt02 = new TH1F("hEffTOFTOTPt02","Efficiency; nrun number; TOF+ITS / TPC",nRuns,0.,nRuns);
        hEffTOFTOTPt02->SetLineWidth(2);
        hEffTOFTOTPt02->SetLineColor(kBlue+2);
        hEffTOFTOTPt02->SetMarkerColor(kBlue+2);
        hEffTOFTOTPt02->SetMarkerStyle(20);
        TH1F *hEffTOFTOTPt1 = new TH1F("hEffTOFTOTPt1","Efficiency; nrun number; TOF+ITS / TPC",nRuns,0.,nRuns);
        hEffTOFTOTPt1->SetLineWidth(2);
        hEffTOFTOTPt1->SetLineColor(kBlue+2);
        hEffTOFTOTPt1->SetMarkerColor(kBlue+2);
        hEffTOFTOTPt1->SetMarkerStyle(20);
        TH1F *hEffTOFTOTPt10 = new TH1F("hEffTOFTOTPt10","Efficiency; nrun number; TOF+ITS / TPC",nRuns,0.,nRuns);
        hEffTOFTOTPt10->SetLineWidth(2);
        hEffTOFTOTPt10->SetLineColor(kBlue+2);
        hEffTOFTOTPt10->SetMarkerColor(kBlue+2);
        hEffTOFTOTPt10->SetMarkerStyle(20);
                        
//************************************************ pileup PLOTS (5) *********************************************//
    
        TH1F *hNumPilVtx = new TH1F("hNumPilVtx","Mean pileup vertices number",nRuns,0.,nRuns);
        hNumPilVtx->SetLineWidth(2);
        hNumPilVtx->SetLineColor(kBlue+2);
        hNumPilVtx->SetMarkerColor(kBlue+2);
        hNumPilVtx->SetMarkerStyle(20);
        hNumPilVtx->GetYaxis()->SetRangeUser(0.,3.);
        hNumPilVtx->GetXaxis()->SetTitle("nrun number");
        hNumPilVtx->GetYaxis()->SetTitle("pileup vertices mean number");
                            
                            
        TH1F *hNumPilTrkl = new TH1F("hNumPilTrkl","Mean Tracklets number for prim. vertex w/ pileup",nRuns,0.,nRuns);
        hNumPilTrkl->SetLineWidth(2);
        hNumPilTrkl->SetLineColor(2);
        hNumPilTrkl->SetMarkerColor(2);
        hNumPilTrkl->SetMarkerStyle(20);
        hNumPilTrkl->GetYaxis()->SetRangeUser(0.,60.);
        hNumPilTrkl->GetXaxis()->SetTitle("nrun number");
        hNumPilTrkl->GetYaxis()->SetTitle("pileup vertex Tracklet number");
                            
        TH1F *hNumNoPilTrkl = new TH1F("hNumNoPilTrkl","Mean Tracklets number for prim. vertex w/out pileup",nRuns,0.,nRuns);
        hNumNoPilTrkl->SetLineWidth(2);
        hNumNoPilTrkl->SetLineColor(kBlue+2);
        hNumNoPilTrkl->SetMarkerColor(kBlue+2);
        hNumNoPilTrkl->SetMarkerStyle(24);
                            
        TH1F *hNumPilCL1 = new TH1F("hNumPilCL1","Mean SPD1 cluster number for prim. vertex w/ pileup tagging",nRuns,0.,nRuns);
        hNumPilCL1->SetLineWidth(2);
        hNumPilCL1->SetLineColor(2);
        hNumPilCL1->SetMarkerColor(2);
        hNumPilCL1->SetMarkerStyle(20);
        hNumPilCL1->GetYaxis()->SetRangeUser(0.,100.);
        hNumPilCL1->GetXaxis()->SetTitle("nrun number");
        hNumPilCL1->GetYaxis()->SetTitle("pileup vertex SPD1 cluster number");
                            
        TH1F *hNumNoPilCL1 = new TH1F("hNumNoPilCL1","Mean SPD1 cluster number for prim. vertex w/out pileup tagging",nRuns,0.,nRuns);
        hNumNoPilCL1->SetLineWidth(2);
        hNumNoPilCL1->SetLineColor(kBlue+2);
        hNumNoPilCL1->SetMarkerColor(kBlue+2);
        hNumNoPilCL1->SetMarkerStyle(24);
                            
//******************************************** ITSpureSA PLOTS (32) *********************************************//
                
        TH1F *h0=new TH1F("h0","h0",nRuns,-0.5,nRuns-0.5);
        TH1F *h1=new TH1F("h1","h1",nRuns,-0.5,nRuns-0.5);
        TH1F *h2=new TH1F("h2","h2",nRuns,-0.5,nRuns-0.5);
        TH1F *h3=new TH1F("h3","h4",nRuns,-0.5,nRuns-0.5);
        TH1F *h4=new TH1F("h4","h5",nRuns,-0.5,nRuns-0.5);
        TH1F *h5=new TH1F("h5","h5",nRuns,-0.5,nRuns-0.5);
        TH1F *h6=new TH1F("h6","h6",nRuns,-0.5,nRuns-0.5);
        TH1F *h7=new TH1F("h7","h7",nRuns,-0.5,nRuns-0.5);
        TH1F *h8=new TH1F("h8","h8",nRuns,-0.5,nRuns-0.5);
        TH1F *h9=new TH1F("h9","h9",nRuns,-0.5,nRuns-0.5);
        TH1F *h10=new TH1F("h10","h10",nRuns,-0.5,nRuns-0.5);
        TH1F *h11=new TH1F("h11","h11",nRuns,-0.5,nRuns-0.5);
        TH1F *h12=new TH1F("h12","h12",nRuns,-0.5,nRuns-0.5);
        TH1F *h13=new TH1F("h13","h13",nRuns,-0.5,nRuns-0.5);
        TH1F *h14=new TH1F("h14","h14",nRuns,-0.5,nRuns-0.5);
        TH1F *h15=new TH1F("h15","h15",nRuns,-0.5,nRuns-0.5);
        TH1F *h16=new TH1F("h16","h16",nRuns,-0.5,nRuns-0.5);
        TH1F *h17=new TH1F("h17","h17",nRuns,-0.5,nRuns-0.5);
                         
    TH1F *hchi2=new TH1F("hchi2","hchi2",nRuns,-0.5,nRuns-0.5);
    hchi2->SetLineWidth(2);
    //    hchi2->SetLineColor(kRed);
    hchi2->SetMarkerColor(kRed);
    hchi2->SetMarkerStyle(20);
    
    TH1F *hchi2SA=new TH1F("hchi2SA","hchi2SA",nRuns,-0.5,nRuns-0.5);
    hchi2SA->SetLineWidth(2);
    //    hchi2SA->SetLineColor(kBlue);
    hchi2SA->SetMarkerColor(kBlue);
    hchi2SA->SetMarkerStyle(20);
    
    TH2F *hOccEta1=new TH2F("hOccEta1","hOccEta1",nRuns,-0.5,nRuns-0.5,2,0.5,2.5);
    hOccEta1->SetTitle("Layer 1 - #h");
    TH2F *hOccPhi1=new TH2F("hOccPhi1","hOccPhi1",nRuns,-0.5,nRuns-0.5,40,0.5,40.5);
    hOccPhi1->SetTitle("Layer 1 - #phi");
    TH2F *hOccEta2=new TH2F("hOccEta2","hOccEta2",nRuns,-0.5,nRuns-0.5,2,0.5,2.5);
    hOccEta2->SetTitle("Layer 2 - #h");
    TH2F *hOccPhi2=new TH2F("hOccPhi2","hOccPhi2",nRuns,-0.5,nRuns-0.5,40,0.5,40.5);
    hOccPhi2->SetTitle("Layer 2 - #phi");
    TH2F *hOccEta3=new TH2F("hOccEta3","hOccEta3",nRuns,-0.5,nRuns-0.5,2,0.5,2.5);
    hOccEta3->SetTitle("Layer 3 - #h");
    TH2F *hOccPhi3=new TH2F("hOccPhi3","hOccPhi3",nRuns,-0.5,nRuns-0.5,40,0.5,40.5);
    hOccPhi3->SetTitle("Layer 3 - #phi");
    TH2F *hOccEta4=new TH2F("hOccEta4","hOccEta4",nRuns,-0.5,nRuns-0.5,2,0.5,2.5);
    hOccEta4->SetTitle("Layer 4 - #h");
    TH2F *hOccPhi4=new TH2F("hOccPhi4","hOccPhi4",nRuns,-0.5,nRuns-0.5,40,0.5,40.5);
    hOccPhi4->SetTitle("Layer 4 - #phi");
    TH2F *hOccEta5=new TH2F("hOccEta5","hOccEta5",nRuns,-0.5,nRuns-0.5,2,0.5,2.5);
    hOccEta5->SetTitle("Layer 5 - #h");
    TH2F *hOccPhi5=new TH2F("hOccPhi5","hOccPhi5",nRuns,-0.5,nRuns-0.5,40,0.5,40.5);
    hOccPhi5->SetTitle("Layer 5 - #phi");
    TH2F *hOccEta6=new TH2F("hOccEta6","hOccEta6",nRuns,-0.5,nRuns-0.5,2,0.5,2.5);
    hOccEta6->SetTitle("Layer 6 - #h");
    TH2F *hOccPhi6=new TH2F("hOccPhi6","hOccPhi6",nRuns,-0.5,nRuns-0.5,40,0.5,40.5);
    hOccPhi6->SetTitle("Layer 6 - #phi");
    
    //**************************************** PID PLOTS (4) *********************************************//
    
    TH1F *hNsig02 = new TH1F("hNsig02","Pion nsigma pt=0.2 GeV/c",nRuns,0.,nRuns);
    hNsig02->SetLineWidth(2);
    hNsig02->SetLineColor(kBlue);
    hNsig02->SetMarkerColor(kBlue);
    hNsig02->SetMarkerStyle(20);
    
    TH1F *hNsig05 = new TH1F("hNsig05","Pion nsigma pt=0.5 GeV/c",nRuns,0.,nRuns);
    hNsig05->SetLineWidth(2);
    hNsig05->SetLineColor(kRed);
    hNsig05->SetMarkerColor(kRed);
    hNsig05->SetMarkerStyle(21);
    
    TH1F *hNsig1 = new TH1F("hNsig1","Pion nsigma pt=1.0 GeV/c",nRuns,0.,nRuns);
    hNsig1->SetLineWidth(2);
    hNsig1->SetLineColor(kGreen+1);
    hNsig1->SetMarkerColor(kGreen+1);
    hNsig1->SetMarkerStyle(22);
    
    TH1F *hNsig3 = new TH1F("hNsig3","Pion nsigma pt=2.0 GeV/c",nRuns,0.,nRuns);
    hNsig3->SetLineWidth(2);
    hNsig3->SetLineColor(kViolet);
    hNsig3->SetMarkerColor(kViolet);
    hNsig3->SetMarkerStyle(23);

    //**************************************** DCA PLOTS (16) *********************************************//
    
    TH1F *hmDCA05 = new TH1F("hmDCA05","mean DCA for TPCITS tracks, 0.55<pt<0.56 GeV/c",nRuns,0.,nRuns);
    hmDCA05->SetLineWidth(2);
    hmDCA05->SetLineColor(kBlue);
    hmDCA05->SetMarkerColor(kBlue);
    hmDCA05->SetMarkerStyle(20);
    
    TH1F *hrmsDCA05 = new TH1F("hrmsDCA05","RMS DCA for TPCITS tracks, 0.55<pt<0.56 GeV/c",nRuns,0.,nRuns);
    hrmsDCA05->SetLineWidth(2);
    hrmsDCA05->SetLineColor(kBlue);
    hrmsDCA05->SetMarkerColor(kBlue);
    hrmsDCA05->SetMarkerStyle(20);
    
    TH1F *hmDCA1 = new TH1F("hmDCA1","mean DCA for TPCITS tracks, 1.05<pt<1.07 GeV/c",nRuns,0.,nRuns);
    hmDCA1->SetLineWidth(2);
    hmDCA1->SetLineColor(kRed);
    hmDCA1->SetMarkerColor(kRed);
    hmDCA1->SetMarkerStyle(21);
    
    TH1F *hrmsDCA1 = new TH1F("hrmsDCA1","RMS DCA for TPCITS tracks, 1.05<pt<1.07 GeV/c",nRuns,0.,nRuns);
    hrmsDCA1->SetLineWidth(2);
    hrmsDCA1->SetLineColor(kRed);
    hrmsDCA1->SetMarkerColor(kRed);
    hrmsDCA1->SetMarkerStyle(21);
    
    TH1F *hmDCA5 = new TH1F("hmDCA5","mean DCA for TPCITS tracks, 4.1<pt<5.2 GeV/c",nRuns,0.,nRuns);
    hmDCA5->SetLineWidth(2);
    hmDCA5->SetLineColor(kOrange+1);
    hmDCA5->SetMarkerColor(kOrange+1);
    hmDCA5->SetMarkerStyle(22);
    
    TH1F *hrmsDCA5 = new TH1F("hrmsDCA5","RMS DCA for TPCITS tracks, 4.1<pt<5.2 GeV/c",nRuns,0.,nRuns);
    hrmsDCA5->SetLineWidth(2);
    hrmsDCA5->SetLineColor(kOrange+1);
    hrmsDCA5->SetMarkerColor(kOrange+1);
    hrmsDCA5->SetMarkerStyle(22);
    
    TH1F *hmDCA10 = new TH1F("hmDCA10","mean DCA for TPCITS tracks, 7.0<pt<8.8 GeV/c",nRuns,0.,nRuns);
    hmDCA10->SetLineWidth(2);
    hmDCA10->SetLineColor(kMagenta+2);
    hmDCA10->SetMarkerColor(kMagenta+2);
    hmDCA10->SetMarkerStyle(23);
    
    TH1F *hrmsDCA10 = new TH1F("hrmsDCA10","RMS DCA for TPCITS tracks, 7.0<pt<8.8 GeV/c",nRuns,0.,nRuns);
    hrmsDCA10->SetLineWidth(2);
    hrmsDCA10->SetLineColor(kMagenta+2);
    hrmsDCA10->SetMarkerColor(kMagenta+2);
    hrmsDCA10->SetMarkerStyle(23);
    
    TH1F *hmDCAz05 = new TH1F("hmDCAz05","mean DCAz for TPCITS tracks, 0.55<pt<0.56 GeV/c",nRuns,0.,nRuns);
    hmDCAz05->SetLineWidth(2);
    hmDCAz05->SetLineColor(kBlue);
    hmDCAz05->SetMarkerColor(kBlue);
    hmDCAz05->SetMarkerStyle(20);
    
    TH1F *hrmsDCAz05 = new TH1F("hrmsDCAz05","RMS DCAz for TPCITS tracks, 0.55<pt<0.56 GeV/c",nRuns,0.,nRuns);
    hrmsDCAz05->SetLineWidth(2);
    hrmsDCAz05->SetLineColor(kBlue);
    hrmsDCAz05->SetMarkerColor(kBlue);
    hrmsDCAz05->SetMarkerStyle(20);
    
    TH1F *hmDCAz1 = new TH1F("hmDCAz1","mean DCAz for TPCITS tracks, 1.05<pt<1.07 GeV/c",nRuns,0.,nRuns);
    hmDCAz1->SetLineWidth(2);
    hmDCAz1->SetLineColor(kRed);
    hmDCAz1->SetMarkerColor(kRed);
    hmDCAz1->SetMarkerStyle(21);
    
    TH1F *hrmsDCAz1 = new TH1F("hrmsDCAz1","RMS DCAz for TPCITS tracks, 1.05<pt<1.07 GeV/c",nRuns,0.,nRuns);
    hrmsDCAz1->SetLineWidth(2);
    hrmsDCAz1->SetLineColor(kRed);
    hrmsDCAz1->SetMarkerColor(kRed);
    hrmsDCAz1->SetMarkerStyle(21);
    
    TH1F *hmDCAz5 = new TH1F("hmDCAz5","mean DCAz for TPCITS tracks, 4.1<pt<5.2 GeV/c",nRuns,0.,nRuns);
    hmDCAz5->SetLineWidth(2);
    hmDCAz5->SetLineColor(kOrange+1);
    hmDCAz5->SetMarkerColor(kOrange+1);
    hmDCAz5->SetMarkerStyle(22);
    
    TH1F *hrmsDCAz5 = new TH1F("hrmsDCAz5","RMS DCAz for TPCITS tracks, 4.1<pt<5.2 GeV/c",nRuns,0.,nRuns);
    hrmsDCAz5->SetLineWidth(2);
    hrmsDCAz5->SetLineColor(kOrange+1);
    hrmsDCAz5->SetMarkerColor(kOrange+1);
    hrmsDCAz5->SetMarkerStyle(22);
    
    TH1F *hmDCAz10 = new TH1F("hmDCAz10","mean DCAz for TPCITS tracks, 7.0<pt<8.8 GeV/c",nRuns,0.,nRuns);
    hmDCAz10->SetLineWidth(2);
    hmDCAz10->SetLineColor(kMagenta+2);
    hmDCAz10->SetMarkerColor(kMagenta+2);
    hmDCAz10->SetMarkerStyle(23);
    
    TH1F *hrmsDCAz10 = new TH1F("hrmsDCAz10","RMS DCAz for TPCITS tracks, 7.0<pt<8.8 GeV/c",nRuns,0.,nRuns);
    hrmsDCAz10->SetLineWidth(2);
    hrmsDCAz10->SetLineColor(kMagenta+2);
    hrmsDCAz10->SetMarkerColor(kMagenta+2);
    hrmsDCAz10->SetMarkerStyle(23);

    //****************************************************************************************************//


  // lista
  
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
//  lista.Add(histoNmodEmpty);
   lista.Add(histoFracDead3);
    lista.Add(histoFracDead4);
    lista.Add(histoEvwSDD);
    lista.Add(histoFlagON3);
    lista.Add(histoFlagON4);
    lista.Add(histoFlagMPV3);
    lista.Add(histoFlagMPV4);
    lista.Add(histoFlagminT);
    lista.Add(histoFlagmeanT);


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
  lista.Add(hpileupSPD);

  lista.Add(histodEdxLay5);
  lista.Add(histodEdxLay6);
  lista.Add(histoChargeRatioLay5);
  lista.Add(histoChargeRatioLay6);
    lista.Add(histoFracBadn5);
    lista.Add(histoFracBadp5);
    lista.Add(histoFracBadn6);
    lista.Add(histoFracBadp6);
    lista.Add(histoFlagON5n);
    lista.Add(histoFlagON5p);
    lista.Add(histoFlagON6n);
    lista.Add(histoFlagON6p);
    lista.Add(histoFlagMPV5);
    lista.Add(histoFlagMPV6);
    lista.Add(histoFlagCR5);
    lista.Add(histoFlagCR6);

    lista.Add(histoFracSPD1);
    lista.Add(histoFracSPD2);
    lista.Add(histoFlagON1);
    lista.Add(histoFlagON2);
    
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
    lista.Add(histoTrackMI1);
    lista.Add(histoTrackMI2);
    lista.Add(histoTrackMI3);
    lista.Add(histoTrackMI4);
    lista.Add(histoTrackMI5);
    lista.Add(histoTrackMI6);
    lista.Add(histoTrackSA1);
    lista.Add(histoTrackSA2);
    lista.Add(histoTrackSA3);
    lista.Add(histoTrackSA4);
    lista.Add(histoTrackSA5);
    lista.Add(histoTrackSA6);
    
    lista.Add(hEffTOFSPDPt02);
    lista.Add(hEffTOFSPDPt1);
    lista.Add(hEffTOFSPDPt10);
    lista.Add(hEffTOFoneSPDPt02);
    lista.Add(hEffTOFoneSPDPt1);
    lista.Add(hEffTOFoneSPDPt10);
    lista.Add(hEffTOF2Pt02);
    lista.Add(hEffTOF2Pt1);
    lista.Add(hEffTOF2Pt10);
    lista.Add(hEffTOF3Pt02);
    lista.Add(hEffTOF3Pt1);
    lista.Add(hEffTOF3Pt10);
    lista.Add(hEffTOF4Pt02);
    lista.Add(hEffTOF4Pt1);
    lista.Add(hEffTOF4Pt10);
    lista.Add(hEffTOF5Pt02);
    lista.Add(hEffTOF5Pt1);
    lista.Add(hEffTOF5Pt10);
    lista.Add(hEffTOF6Pt02);
    lista.Add(hEffTOF6Pt1);
    lista.Add(hEffTOF6Pt10);
    lista.Add(hEffTOFTOTPt02);
    lista.Add(hEffTOFTOTPt1);
    lista.Add(hEffTOFTOTPt10);
                           
    lista.Add(hNumPilVtx);
    lista.Add(hNumPilTrkl);
    lista.Add(hNumNoPilTrkl);
    lista.Add(hNumPilCL1);
    lista.Add(hNumNoPilCL1);
    
    lista.Add(h0);
    lista.Add(h1);
    lista.Add(h2);
    lista.Add(h3);
    lista.Add(h4);
    lista.Add(h5);
    lista.Add(h6);
    lista.Add(h7);
    lista.Add(h8);
    lista.Add(h9);
    lista.Add(h10);
    lista.Add(h11);
    lista.Add(h12);
    lista.Add(h13);
    lista.Add(h14);
    lista.Add(h15);
    lista.Add(h16);
    lista.Add(h17);
    lista.Add(hchi2);
    lista.Add(hchi2SA);
    lista.Add(hOccEta1);
    lista.Add(hOccPhi1);
    lista.Add(hOccEta2);
    lista.Add(hOccPhi2);
    lista.Add(hOccEta3);
    lista.Add(hOccPhi3);
    lista.Add(hOccEta4);
    lista.Add(hOccPhi5);
    lista.Add(hOccEta5);
    lista.Add(hOccPhi5);
    lista.Add(hOccEta6);
    lista.Add(hOccPhi6);
    
    lista.Add(hNsig02);
    lista.Add(hNsig05);
    lista.Add(hNsig1);
    lista.Add(hNsig3);

    lista.Add(hmDCA05);
    lista.Add(hrmsDCA05);
    lista.Add(hmDCA1);
    lista.Add(hrmsDCA1);
    lista.Add(hmDCA5);
    lista.Add(hrmsDCA5);
    lista.Add(hmDCA10);
    lista.Add(hrmsDCA10);
    lista.Add(hmDCAz05);
    lista.Add(hrmsDCAz05);
    lista.Add(hmDCAz1);
    lista.Add(hrmsDCAz1);
    lista.Add(hmDCAz5);
    lista.Add(hrmsDCAz5);
    lista.Add(hmDCAz10);
    lista.Add(hrmsDCAz10);

  //lista.Add();
  
    //******************************************** FILLING ************************************************//

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
        histoFracDead3->SetBinContent(i+1,fracDead3);
        histoFracDead3->SetBinError(i+1,errfracDead3);
        histoFracDead3->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        histoFracDead4->SetBinContent(i+1,fracDead4);
        histoFracDead4->SetBinError(i+1,errfracDead4);
        histoFracDead4->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));

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
//      histoNmodEmpty->SetBinContent(i+1,EmptyModulesSDD);
//      histoNmodEmpty->SetBinError(i+1,0.000001);
//      histoNmodEmpty->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));

        histoEvwSDD->SetBinContent(i+1,fracEvWithSDD);
        histoEvwSDD->SetBinError(i+1,errfracEvWithSDD);
        histoEvwSDD->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));

        histoFlagON3->SetBinContent(i+1,FlagSDD1);
        histoFlagON3->SetBinError(i+1,0.01);
        histoFlagON3->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        histoFlagON4->SetBinContent(i+1,FlagSDD2);
        histoFlagON4->SetBinError(i+1,0.01);
        histoFlagON4->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        histoFlagMPV3->SetBinContent(i+1,FlagdEdx3);
        histoFlagMPV3->SetBinError(i+1,0.01);
        histoFlagMPV3->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        histoFlagMPV4->SetBinContent(i+1,FlagdEdx4*2.);
        histoFlagMPV4->SetBinError(i+1,0.01);
        histoFlagMPV4->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        histoFlagminT->SetBinContent(i+1,FlagMinTime);
        histoFlagminT->SetBinError(i+1,0.01);
        histoFlagminT->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        histoFlagmeanT->SetBinContent(i+1,FlagMeanTime*2.);
        histoFlagmeanT->SetBinError(i+1,0.01);
        histoFlagmeanT->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));


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
      hMeanVxOCDB->SetBinContent(i+1,meanVtxOCDBx);
      hMeanVxOCDB->SetBinError(i+1,0.00000001);
      hMeanVxOCDB->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      hMeanVyOCDB->SetBinContent(i+1,meanVtxOCDBy);
      hMeanVyOCDB->SetBinError(i+1,0.00000001);
      hMeanVyOCDB->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      hMeanVzOCDB->SetBinContent(i+1,meanVtxOCDBz);
      hMeanVzOCDB->SetBinError(i+1,0.00000001);
      hMeanVzOCDB->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      hSigmaVxOCDB->SetBinContent(i+1,sigmaVtxOCDBx);
      hSigmaVxOCDB->SetBinError(i+1,0.00000001);
      hSigmaVxOCDB->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      hSigmaVyOCDB->SetBinContent(i+1,sigmaVtxOCDBy);
      hSigmaVyOCDB->SetBinError(i+1,0.00000001);
      hSigmaVyOCDB->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      hSigmaVzOCDB->SetBinContent(i+1,sigmaVtxOCDBz);
      hSigmaVzOCDB->SetBinError(i+1,0.00000001);
      hSigmaVzOCDB->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));

      hpileupSPD->SetBinContent(i+1,pileupSPD);
      hpileupSPD->SetBinError(i+1,errpileupSPD);
      hpileupSPD->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));


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
        
        histoFracBadn5->SetBinContent(i+1,FracBadn5);
        histoFracBadn5->SetBinError(i+1,errFracBadn5);
        histoFracBadn5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        
        histoFracBadp5->SetBinContent(i+1,FracBadp5);
        histoFracBadp5->SetBinError(i+1,errFracBadp5);
        histoFracBadp5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        
        histoFracBadn6->SetBinContent(i+1,FracBadn6);
        histoFracBadn6->SetBinError(i+1,errFracBadn6);
        histoFracBadn6->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        
        histoFracBadp6->SetBinContent(i+1,FracBadp6);
        histoFracBadp6->SetBinError(i+1,errFracBadp6);
        histoFracBadp6->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));

        histoFlagON5n->SetBinContent(i+1,FlagSSD1n);
        histoFlagON5n->SetBinError(i+1,0.01);
        histoFlagON5n->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        histoFlagON5p->SetBinContent(i+1,FlagSSD1p*1.2);
        histoFlagON5p->SetBinError(i+1,0.01);
        histoFlagON5p->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));

        histoFlagON6n->SetBinContent(i+1,FlagSSD2n);
        histoFlagON6n->SetBinError(i+1,0.01);
        histoFlagON6n->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        histoFlagON6p->SetBinContent(i+1,FlagSSD2p*1.2);
        histoFlagON6p->SetBinError(i+1,0.01);
        histoFlagON6p->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));

        histoFlagMPV5->SetBinContent(i+1,FlagdEdx5);
        histoFlagMPV5->SetBinError(i+1,0.01);
        histoFlagMPV5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        histoFlagMPV6->SetBinContent(i+1,FlagdEdx6*2.);
        histoFlagMPV6->SetBinError(i+1,0.01);
        histoFlagMPV6->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));

        histoFlagCR5->SetBinContent(i+1,FlagChR5);
        histoFlagCR5->SetBinError(i+1,0.01);
        histoFlagCR5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        histoFlagCR6->SetBinContent(i+1,FlagChR6*2.);
        histoFlagCR6->SetBinError(i+1,0.01);
        histoFlagCR6->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    
      //======= SPD =======
      histoFracSPD1->SetBinContent(i+1,FracSPD1);
      histoFracSPD1->SetBinError(i+1,errFracSPD1);
      histoFracSPD1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      histoFracSPD2->SetBinContent(i+1,FracSPD2);
      histoFracSPD2->SetBinError(i+1,errFracSPD2);
      histoFracSPD2->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        histoFlagON1->SetBinContent(i+1,FlagSPD1);
        histoFlagON1->SetBinError(i+1,0.01);
        histoFlagON1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        histoFlagON2->SetBinContent(i+1,FlagSPD2*2);
        histoFlagON2->SetBinError(i+1,0.01);
        histoFlagON2->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));

        //======= MATCHING =======
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
    
        histoTrackMI1->SetBinContent(i+1,FracTrackMI1);
        histoTrackMI1->SetBinError(i+1,errFracTrackMI1);
        histoTrackMI1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        histoTrackMI2->SetBinContent(i+1,FracTrackMI2);
        histoTrackMI2->SetBinError(i+1,errFracTrackMI2);
        histoTrackMI2->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        histoTrackMI3->SetBinContent(i+1,FracTrackMI3);
        histoTrackMI3->SetBinError(i+1,errFracTrackMI3);
        histoTrackMI3->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        histoTrackMI4->SetBinContent(i+1,FracTrackMI4);
        histoTrackMI4->SetBinError(i+1,errFracTrackMI4);
        histoTrackMI4->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        histoTrackMI5->SetBinContent(i+1,FracTrackMI5);
        histoTrackMI5->SetBinError(i+1,errFracTrackMI5);
        histoTrackMI5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        histoTrackMI6->SetBinContent(i+1,FracTrackMI6);
        histoTrackMI6->SetBinError(i+1,errFracTrackMI6);
        histoTrackMI6->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        histoTrackSA1->SetBinContent(i+1,FracTrackSA1);
        histoTrackSA1->SetBinError(i+1,errFracTrackSA1);
        histoTrackSA1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        histoTrackSA2->SetBinContent(i+1,FracTrackSA2);
        histoTrackSA2->SetBinError(i+1,errFracTrackSA2);
        histoTrackSA2->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        histoTrackSA3->SetBinContent(i+1,FracTrackSA3);
        histoTrackSA3->SetBinError(i+1,errFracTrackSA3);
        histoTrackSA3->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        histoTrackSA4->SetBinContent(i+1,FracTrackSA4);
        histoTrackSA4->SetBinError(i+1,errFracTrackSA4);
        histoTrackSA4->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        histoTrackSA5->SetBinContent(i+1,FracTrackSA5);
        histoTrackSA5->SetBinError(i+1,errFracTrackSA5);
        histoTrackSA5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        histoTrackSA6->SetBinContent(i+1,FracTrackSA6);
        histoTrackSA6->SetBinError(i+1,errFracTrackSA6);
        histoTrackSA6->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));

        hEffTOF6Pt02->SetBinContent(i+1,Eff6Pt02TOF);
        hEffTOF6Pt02->SetBinError(i+1,errEff6Pt02TOF);
        hEffTOF6Pt02->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        hEffTOF6Pt1->SetBinContent(i+1,Eff6Pt1TOF);
        hEffTOF6Pt1->SetBinError(i+1,errEff6Pt1TOF);
        hEffTOF6Pt1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        hEffTOF6Pt10->SetBinContent(i+1,Eff6Pt10TOF);
        hEffTOF6Pt10->SetBinError(i+1,errEff6Pt10TOF);
        hEffTOF6Pt10->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        hEffTOF5Pt02->SetBinContent(i+1,Eff5Pt02TOF);
        hEffTOF5Pt02->SetBinError(i+1,errEff5Pt02TOF);
        hEffTOF5Pt02->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        hEffTOF5Pt1->SetBinContent(i+1,Eff5Pt1TOF);
        hEffTOF5Pt1->SetBinError(i+1,errEff5Pt1TOF);
        hEffTOF5Pt1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        hEffTOF5Pt10->SetBinContent(i+1,Eff5Pt10TOF);
        hEffTOF5Pt10->SetBinError(i+1,errEff5Pt10TOF);
        hEffTOF5Pt10->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        hEffTOF4Pt02->SetBinContent(i+1,Eff4Pt02TOF);
        hEffTOF4Pt02->SetBinError(i+1,errEff4Pt1TOF);
        hEffTOF4Pt02->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        hEffTOF4Pt1->SetBinContent(i+1,Eff4Pt1TOF);
        hEffTOF4Pt1->SetBinError(i+1,errEff4Pt1TOF);
        hEffTOF4Pt1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        hEffTOF4Pt10->SetBinContent(i+1,Eff4Pt10TOF);
        hEffTOF4Pt10->SetBinError(i+1,errEff4Pt10TOF);
        hEffTOF4Pt10->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        hEffTOF3Pt02->SetBinContent(i+1,Eff3Pt02TOF);
        hEffTOF3Pt02->SetBinError(i+1,errEff3Pt02TOF);
        hEffTOF3Pt02->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        hEffTOF3Pt1->SetBinContent(i+1,Eff3Pt1TOF);
        hEffTOF3Pt1->SetBinError(i+1,errEff3Pt1TOF);
        hEffTOF3Pt1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        hEffTOF3Pt10->SetBinContent(i+1,Eff3Pt10TOF);
        hEffTOF3Pt10->SetBinError(i+1,errEff3Pt10TOF);
        hEffTOF3Pt10->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        hEffTOF2Pt02->SetBinContent(i+1,Eff2Pt02TOF);
        hEffTOF2Pt02->SetBinError(i+1,errEff2Pt02TOF);
        hEffTOF2Pt02->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        hEffTOF2Pt1->SetBinContent(i+1,Eff2Pt1TOF);
        hEffTOF2Pt1->SetBinError(i+1,errEff2Pt1TOF);
        hEffTOF2Pt1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        hEffTOF2Pt10->SetBinContent(i+1,Eff2Pt10TOF);
        hEffTOF2Pt10->SetBinError(i+1,errEff2Pt10TOF);
        hEffTOF2Pt10->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        hEffTOFSPDPt02->SetBinContent(i+1,EffSPDPt02TOF);
        hEffTOFSPDPt02->SetBinError(i+1,errEffSPDPt02TOF);
        hEffTOFSPDPt02->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        hEffTOFSPDPt1->SetBinContent(i+1,EffSPDPt1TOF);
        hEffTOFSPDPt1->SetBinError(i+1,errEffSPDPt1TOF);
        hEffTOFSPDPt1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        hEffTOFSPDPt10->SetBinContent(i+1,EffSPDPt10TOF);
        hEffTOFSPDPt10->SetBinError(i+1,errEffSPDPt10TOF);
        hEffTOFSPDPt10->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        hEffTOFoneSPDPt02->SetBinContent(i+1,EffoneSPDPt02TOF);
        hEffTOFoneSPDPt02->SetBinError(i+1,errEffoneSPDPt02TOF);
        hEffTOFoneSPDPt02->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        hEffTOFoneSPDPt1->SetBinContent(i+1,EffoneSPDPt1TOF);
        hEffTOFoneSPDPt1->SetBinError(i+1,errEffoneSPDPt1TOF);
        hEffTOFoneSPDPt1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        hEffTOFoneSPDPt10->SetBinContent(i+1,EffoneSPDPt10TOF);
        hEffTOFoneSPDPt10->SetBinError(i+1,errEffoneSPDPt10TOF);
        hEffTOFoneSPDPt10->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        hEffTOFTOTPt02->SetBinContent(i+1,EffTOTPt02TOF);
        hEffTOFTOTPt02->SetBinError(i+1,errEffTOTPt02TOF);
        hEffTOFTOTPt02->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        hEffTOFTOTPt1->SetBinContent(i+1,EffTOTPt1TOF);
        hEffTOFTOTPt1->SetBinError(i+1,errEffTOTPt1TOF);
        hEffTOFTOTPt1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        hEffTOFTOTPt10->SetBinContent(i+1,EffTOTPt10TOF);
        hEffTOFTOTPt10->SetBinError(i+1,errEffTOTPt10TOF);
        hEffTOFTOTPt10->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        
        //======= PILEUP =======
       hNumPilVtx->SetBinContent(i+1,npilvtx);
        hNumPilVtx->SetBinError(i+1,errnpilvtx);
        hNumPilVtx->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        hNumPilTrkl->SetBinContent(i+1,ntrklpil);
        hNumPilTrkl->SetBinError(i+1,errntrklpil);
        hNumPilTrkl->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        hNumNoPilTrkl->SetBinContent(i+1,ntrklnopil);
        hNumNoPilTrkl->SetBinError(i+1,errntrklnopil);
        hNumNoPilTrkl->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        hNumPilCL1->SetBinContent(i+1,ncl1pil);
        hNumPilCL1->SetBinError(i+1,errncl1pil);
        hNumPilCL1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        hNumNoPilCL1->SetBinContent(i+1,ncl1nopil);
        hNumNoPilCL1->SetBinError(i+1,errncl1nopil);
        hNumNoPilCL1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
  
        //======= ITSpureSA =======
        h0->SetBinContent(i+1,NITSpureSAPtBin0);
        h0->GetXaxis()->SetBinLabel(i+1,Form("%.0d",nrun));
        h0->GetYaxis()->SetTitle("ITSpureSA tracks");
        h0->GetXaxis()->SetTitle("run number");
        h1->SetBinContent(i+1,NITSpureSAPtBin1);
        h1->GetXaxis()->SetBinLabel(i+1,Form("%.0d",nrun));
        h1->GetYaxis()->SetTitle("ITSpureSA tracks");
        h1->GetXaxis()->SetTitle("run number");
        h2->SetBinContent(i+1,NITSpureSAPtBin2);
        h2->GetXaxis()->SetBinLabel(i+1,Form("%.0d",nrun));
        h2->GetYaxis()->SetTitle("ITSpureSA tracks");
        h2->GetXaxis()->SetTitle("run number");
        h3->SetBinContent(i+1,NITSTPCPtBin0);
        h3->GetXaxis()->SetBinLabel(i+1,Form("%.0d",nrun));
        h3->GetYaxis()->SetTitle("ITS+TPC tracks");
        h3->GetXaxis()->SetTitle("run number");
        h4->SetBinContent(i+1,NITSTPCPtBin1);
        h4->GetXaxis()->SetBinLabel(i+1,Form("%.0d",nrun));
        h4->GetYaxis()->SetTitle("ITS+TPC tracks");
        h4->GetXaxis()->SetTitle("run number");
        h5->SetBinContent(i+1,NITSTPCPtBin2);
        h5->GetXaxis()->SetBinLabel(i+1,Form("%.0d",nrun));
        h5->GetYaxis()->SetTitle("ITS+TPC tracks");
        h5->GetXaxis()->SetTitle("run number");
        h6->SetBinContent(i+1,ratioPtBin0);
        h6->GetXaxis()->SetBinLabel(i+1,Form("%.0d",nrun));
        h6->GetYaxis()->SetTitle("(TPC+ITS)+ITSsa/ITSpureSA");
        h6->GetXaxis()->SetTitle("run number");
        h7->SetBinContent(i+1,ratioPtBin1);
        h7->GetXaxis()->SetBinLabel(i+1,Form("%.0d",nrun));
        h7->GetYaxis()->SetTitle("(TPC+ITS)+ITSsa/ITSpureSA");
        h7->GetXaxis()->SetTitle("run number");
        h8->SetBinContent(i+1,ratioPtBin2);
        h8->GetXaxis()->SetBinLabel(i+1,Form("%.0d",nrun));
        h8->GetYaxis()->SetTitle("(TPC+ITS)+ITSsa/ITSpureSA");
        h8->GetXaxis()->SetTitle("run number");
        h9->SetBinContent(i+1,NcluITSpSA);
        h9->SetBinError(i+1,errNcluITSpSA);
        h9->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        h9->GetXaxis()->SetTitle("run number");
        h9->GetYaxis()->SetTitle("mean cluster number (N>3)");
        h10->SetBinContent(i+1,dedx4_3);
        h10->SetBinError(i+1,errdedx4_3);
        h10->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        h10->GetXaxis()->SetTitle("run number");
        h10->GetYaxis()->SetTitle("N(dedx4clu)/N(dedx3clu)");
        h11->SetBinContent(i+1,PtpionpSA);
        h11->SetBinError(i+1,errPtpionpSA);
        h11->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        h11->GetXaxis()->SetTitle("run number");
        h11->GetYaxis()->SetTitle("mean Pt (GeV/c)");
        h12->SetBinContent(i+1,NclupSA0);
        h12->SetBinError(i+1,errNclupSA0);
        h12->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        h12->GetXaxis()->SetTitle("run number");
        h12->GetYaxis()->SetTitle("Fraction of pureSA tracks with cluster in ITS layers");
        h13->SetBinContent(i+1,NclupSA1);
        h13->SetBinError(i+1,errNclupSA1);
        h13->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        h13->GetXaxis()->SetTitle("run number");
        h14->SetBinContent(i+1,NclupSA2);
        h14->SetBinError(i+1,errNclupSA2);
        h14->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        h14->GetXaxis()->SetTitle("run number");
        h15->SetBinContent(i+1,NclupSA3);
        h15->SetBinError(i+1,errNclupSA3);
        h15->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        h15->GetXaxis()->SetTitle("run number");
        h16->SetBinContent(i+1,NclupSA4);
        h16->SetBinError(i+1,errNclupSA4);
        h16->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        h16->GetXaxis()->SetTitle("run number");
        h17->SetBinContent(i+1,NclupSA5);
        h17->SetBinError(i+1,errNclupSA5);
        h17->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        h17->GetXaxis()->SetTitle("run number");
        
        Double_t Lowbin[3]={0.1,0.5,0.9};
        Double_t Upbin[3]={0.2,0.6,1};
        gStyle->SetOptTitle(0);
        gStyle->SetOptStat(0);
        gStyle->SetFillColor(0);
        
        h0->SetTitle(Form("%.1f GeV/c < pt < %.1f GeV/c",Lowbin[0],Upbin[0]));
        h1->SetTitle(Form("%.1f GeV/c < pt < %.1f GeV/c",Lowbin[1],Upbin[1]));
        h2->SetTitle(Form("%.1f GeV/c < pt < %.1f GeV/c",Lowbin[2],Upbin[2]));
        h3->SetTitle(Form("%.1f GeV/c < pt < %.1f GeV/c",Lowbin[0],Upbin[0]));
        h4->SetTitle(Form("%.1f GeV/c < pt < %.1f GeV/c",Lowbin[1],Upbin[1]));
        h5->SetTitle(Form("%.1f GeV/c < pt < %.1f GeV/c",Lowbin[2],Upbin[2]));
        h6->SetTitle(Form("%.1f GeV/c < pt < %.1f GeV/c",Lowbin[0],Upbin[0]));
        h7->SetTitle(Form("%.1f GeV/c < pt < %.1f GeV/c",Lowbin[1],Upbin[1]));
        h8->SetTitle(Form("%.1f GeV/c < pt < %.1f GeV/c",Lowbin[2],Upbin[2]));

        h0->SetMinimum(100);
        //   h0->SetMaximum(2);
        h1->SetMinimum(100);
        //   h1->SetMaximum(2);
        h2->SetMinimum(100);
        //   h2->SetMaximum(2);
        h3->SetMinimum(100);
        h4->SetMinimum(100);
        h5->SetMinimum(100);
        
        h8->SetMinimum(0.);
        h8->SetMaximum(2.5);
        
        h9->SetMinimum(0.);
        h9->SetMaximum(6.5);
        h10->SetMinimum(0.);
        h10->SetMaximum(2.5);
        h11->SetMinimum(0.);
        h11->SetMaximum(1.5);
        h12->SetMinimum(0.);
        h12->SetMaximum(1.2);
        
        h0->SetMarkerStyle(22);
        h1->SetMarkerStyle(23);
        h2->SetMarkerStyle(24);
        h0->SetMarkerColor(2);
        h1->SetMarkerColor(1);
        h2->SetMarkerColor(4);
        h0->SetLineColor(2);
        h1->SetLineColor(1);
        h2->SetLineColor(4);
        
        h3->SetMarkerStyle(22);
        h4->SetMarkerStyle(23);
        h5->SetMarkerStyle(24);
        h3->SetMarkerColor(2);
        h4->SetMarkerColor(1);
        h5->SetMarkerColor(4);
        h3->SetLineColor(2);
        h4->SetLineColor(1);
        h5->SetLineColor(4);
        
        h6->SetMarkerStyle(22);
        h7->SetMarkerStyle(23);
        h8->SetMarkerStyle(24);
        h6->SetMarkerColor(2);
        h7->SetMarkerColor(1);
        h8->SetMarkerColor(4);
        h6->SetLineColor(2);
        h7->SetLineColor(1);
        h8->SetLineColor(4);

    hchi2->SetBinContent(i+1,chi2TPCITS);
    hchi2->SetBinError(i+1,0.01);
    hchi2->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hchi2SA->SetBinContent(i+1,chi2ITSpureSA);
    hchi2SA->SetBinError(i+1,0.01);
    hchi2SA->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
        hchi2->GetXaxis()->SetTitle("run number");
        hchi2->GetYaxis()->SetTitle("#chi^{2} for TPCITS tracks");
        hchi2SA->GetXaxis()->SetTitle("run number");
        hchi2SA->GetYaxis()->SetTitle("#chi^{2} for ITSpureSA tracks");

    for(Int_t j=1;j<3;j++){
        hOccEta1->SetBinContent(i+1,j,occ_eta_1[j-1]);
        hOccEta2->SetBinContent(i+1,j,occ_eta_2[j-1]);
        hOccEta3->SetBinContent(i+1,j,occ_eta_3[j-1]);
        hOccEta4->SetBinContent(i+1,j,occ_eta_4[j-1]);
        hOccEta5->SetBinContent(i+1,j,occ_eta_5[j-1]);
        hOccEta6->SetBinContent(i+1,j,occ_eta_6[j-1]);
        }

    for(Int_t j=1;j<40;j++){
        hOccPhi1->SetBinContent(i+1,j,occ_phi_1[j-1]);
        hOccPhi2->SetBinContent(i+1,j,occ_phi_2[j-1]);
        hOccPhi3->SetBinContent(i+1,j,occ_phi_3[j-1]);
        hOccPhi4->SetBinContent(i+1,j,occ_phi_4[j-1]);
        hOccPhi5->SetBinContent(i+1,j,occ_phi_5[j-1]);
        hOccPhi6->SetBinContent(i+1,j,occ_phi_6[j-1]);
    }
    hOccEta1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hOccEta2->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hOccEta3->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hOccEta4->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hOccEta5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hOccEta6->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hOccPhi1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hOccPhi2->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hOccPhi3->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hOccPhi4->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hOccPhi5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hOccPhi6->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));

    //======= PID =======
    hNsig02->SetBinContent(i+1,nsigmapi02);
    hNsig02->SetBinError(i+1,errnsigmapi02);
    hNsig02->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    
    hNsig05->SetBinContent(i+1,nsigmapi05);
    hNsig05->SetBinError(i+1,errnsigmapi05);
    hNsig05->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    
    hNsig1->SetBinContent(i+1,nsigmapi1);
    hNsig1->SetBinError(i+1,errnsigmapi1);
    hNsig1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    
    hNsig3->SetBinContent(i+1,nsigmapi3);
    hNsig3->SetBinError(i+1,errnsigmapi3);
    hNsig3->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));

    //======= DCA =======
    hmDCA05->SetBinContent(i+1,mdca05);
    hmDCA05->SetBinError(i+1,errmdca05);
    hmDCA05->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    
    hrmsDCA05->SetBinContent(i+1,rmsdca05);
    hrmsDCA05->SetBinError(i+1,errrmsdca05);
    hrmsDCA05->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    
    hmDCA1->SetBinContent(i+1,mdca1);
    hmDCA1->SetBinError(i+1,errmdca1);
    hmDCA1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    
    hrmsDCA1->SetBinContent(i+1,rmsdca1);
    hrmsDCA1->SetBinError(i+1,errrmsdca1);
    hrmsDCA1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    
    hmDCA5->SetBinContent(i+1,mdca5);
    hmDCA5->SetBinError(i+1,errmdca5);
    hmDCA5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    
    hrmsDCA5->SetBinContent(i+1,rmsdca5);
    hrmsDCA5->SetBinError(i+1,errrmsdca5);
    hrmsDCA5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    
    hmDCA10->SetBinContent(i+1,mdca10);
    hmDCA10->SetBinError(i+1,errmdca10);
    hmDCA10->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    
    hrmsDCA10->SetBinContent(i+1,rmsdca10);
    hrmsDCA10->SetBinError(i+1,errrmsdca10);
    hrmsDCA10->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    //
    hmDCAz05->SetBinContent(i+1,mdcaz05);
    hmDCAz05->SetBinError(i+1,errmdcaz05);
    hmDCAz05->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    
    hrmsDCAz05->SetBinContent(i+1,rmsdcaz05);
    hrmsDCAz05->SetBinError(i+1,errrmsdcaz05);
    hrmsDCAz05->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    
    hmDCAz1->SetBinContent(i+1,mdcaz1);
    hmDCAz1->SetBinError(i+1,errmdcaz1);
    hmDCAz1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    
    hrmsDCAz1->SetBinContent(i+1,rmsdcaz1);
    hrmsDCAz1->SetBinError(i+1,errrmsdcaz1);
    hrmsDCAz1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    
    hmDCAz5->SetBinContent(i+1,mdcaz5);
    hmDCAz5->SetBinError(i+1,errmdcaz5);
    hmDCAz5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    
    hrmsDCAz5->SetBinContent(i+1,rmsdcaz5);
    hrmsDCAz5->SetBinError(i+1,errrmsdcaz5);
    hrmsDCAz5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    
    hmDCAz10->SetBinContent(i+1,mdcaz10);
    hmDCAz10->SetBinError(i+1,errmdcaz10);
    hmDCAz10->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    
    hrmsDCAz10->SetBinContent(i+1,rmsdcaz10);
    hrmsDCAz10->SetBinError(i+1,errrmsdcaz10);
    hrmsDCAz10->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    
    } // loop on runs

// ========================== DRAWING HISTOGRAMS ==================================
  //--------  Draw Vertex histograms -------

gStyle->SetOptStat(0);
TCanvas *cVertexDisto;
    Double_t ymin = 0.;
    Double_t ymax = 0.;

if(hMeanVx->GetEntries()>0){
    cVertexDisto=new TCanvas("cVertexDisto","cVertexDisto",1200,800);
    cVertexDisto->Divide(3,2);
    cVertexDisto->cd(1);
    //
    ymin = hMeanVx->GetMinimum();
    ymax = hMeanVx->GetMaximum();
    if(ymin>-10.)hMeanVx->SetMinimum(0.9*ymin);
    else hMeanVx->SetMinimum(0.)
    hMeanVx->SetMaximum(1.1*ymax);
    //
//    hMeanVx->SetMinimum(0.065);
//    hMeanVx->SetMaximum(0.085);
    hMeanVx->GetYaxis()->SetTitle("Vertex X coordinate");
    hMeanVx->GetXaxis()->SetTitle("run number");
    hMeanVx->Draw();
    if(hMeanVxSPD->GetBinContent(1)>-500.)hMeanVxSPD->Draw("same");
    if(hMeanVxOCDB->GetBinContent(1)>-500.)hMeanVxOCDB->Draw("same");
    TLegend* legVtx=new TLegend(0.70,0.83,1.00,0.93);
    legVtx->SetFillColor(kWhite);
    legVtx->SetFillStyle(1001);
    TLegendEntry* entVtx;
    entVtx=legVtx->AddEntry(hMeanVx,"Tracks vertex","PL");
    entVtx->SetTextColor(hMeanVx->GetMarkerColor());
    entVtx=legVtx->AddEntry(hMeanVxSPD,"Tracklets vertex","PL");
    entVtx->SetTextColor(hMeanVxSPD->GetMarkerColor());
    entVtx=legVtx->AddEntry(hMeanVxOCDB,"OCDB diamond","PL");
    entVtx->SetTextColor(hMeanVxOCDB->GetMarkerColor());
    legVtx->Draw();
    
    cVertexDisto->cd(2);
    ymin = hMeanVy->GetMinimum();
    ymax = hMeanVy->GetMaximum();
    if(ymin>-10.)hMeanVy->SetMinimum(0.9*ymin);
    else hMeanVy->SetMinimum(0.)
    hMeanVy->SetMaximum(1.1*ymax);

//        if(hMeanVySPD->GetEntries()>0){  // pp runs

//            hMeanVy->SetMinimum(0.32);
//            hMeanVy->SetMaximum(0.38);
//        }
//        else{                       // PbPb runs
//            hMeanVy->SetMinimum(0.32);
//            hMeanVy->SetMaximum(0.38);
//        }
    hMeanVy->GetYaxis()->SetTitle("Vertex Y coordinate");
    hMeanVy->GetXaxis()->SetTitle("run number");
    hMeanVy->Draw();
    if(hMeanVySPD->GetBinContent(1)>-500.)hMeanVySPD->Draw("same");
    if(hMeanVyOCDB->GetBinContent(1)>-500.)hMeanVyOCDB->Draw("same");
    legVtx->Draw();

    cVertexDisto->cd(3);
    ymin = hMeanVz->GetMinimum();
    ymax = hMeanVz->GetMaximum();
    if(ymin>-10.)hMeanVz->SetMinimum(0.9*ymin);
    else hMeanVz->SetMinimum(0.)
    hMeanVz->SetMaximum(1.1*ymax);
//        hMeanVz->SetMinimum(-8.);
//        hMeanVz->SetMaximum(8.);
        hMeanVz->GetYaxis()->SetTitle("Vertex Z coordinate");
        hMeanVz->GetXaxis()->SetTitle("run number");
        hMeanVz->Draw();
        if(hMeanVzSPD->GetBinContent(1)>-500.)hMeanVzSPD->Draw("same");
	if(hMeanVzOCDB->GetBinContent(1)>-500.)hMeanVzOCDB->Draw("same");
    legVtx->Draw();

    cVertexDisto->cd(4);
        hSigmaVx->SetMinimum(0.);
        hSigmaVx->SetMaximum(0.1);
        hSigmaVx->GetYaxis()->SetTitle("sigma on x coordinate");
        hSigmaVx->GetXaxis()->SetTitle("run number");
        hSigmaVx->Draw();
        if(hSigmaVxSPD->GetBinContent(1)>0)hSigmaVxSPD->Draw("same");
	if(hSigmaVxOCDB->GetBinContent(1)>-500.)hSigmaVxOCDB->Draw("same");
    legVtx->Draw();

    cVertexDisto->cd(5);
        hSigmaVy->SetMinimum(0.);
        hSigmaVy->SetMaximum(0.1);
        hSigmaVy->GetYaxis()->SetTitle("sigma on y coordinate");
        hSigmaVy->GetXaxis()->SetTitle("run number");
        hSigmaVy->Draw();
        if(hSigmaVySPD->GetBinContent(1)>0)hSigmaVySPD->Draw("same");
	if(hSigmaVyOCDB->GetBinContent(1)>-500.)hSigmaVyOCDB->Draw("same");
    legVtx->Draw();

    cVertexDisto->cd(6);
        //   hSigmaVz->SetMinimum(6.);
        //   hSigmaVz->SetMaximum(10.);
        hSigmaVz->GetYaxis()->SetTitle("sigma on z coordinate");
        hSigmaVz->GetXaxis()->SetTitle("run number");
        ymin = hSigmaVz->GetMinimum();
        ymax = hSigmaVz->GetMaximum();
        if(ymin>0.)hSigmaVz->SetMinimum(0.9*ymin);
        else hSigmaVz->SetMinimum(0.)
        hMeanVz->SetMaximum(1.1*ymax);
        hSigmaVz->Draw();
        if(hSigmaVzSPD->GetBinContent(1)>0)hSigmaVzSPD->Draw("same");
	if(hSigmaVzOCDB->GetBinContent(1)>-500.)hSigmaVzOCDB->Draw("same");
    legVtx->Draw();
        cVertexDisto->SaveAs("Vertex_trend.pdf");
        //    pdfFileNames+=" Vertex_trend.pdf";
        }


//
    //--------  Draw histograms of fraction of global tracks w/ hits in layers  -------

TCanvas* cMI;
if(histoTrackMI1->GetEntries()>0){
    cMI=new TCanvas("cMI"," Reconstructed Track With points in ITS");
    cMI->SetGrid();
    histoTrackMI3->Draw();
    histoTrackMI3->SetLineColor(1);
    histoTrackMI3->SetMarkerStyle(20);
    histoTrackMI3->Draw();
    histoTrackMI3->SetMinimum(0.);
    histoTrackMI3->SetMaximum(1.05);
    histoTrackMI4->SetLineColor(2);
    histoTrackMI4->SetMarkerColor(2);
    histoTrackMI4->SetMarkerStyle(22);
    histoTrackMI4->Draw("same");
    histoTrackMI1->SetLineColor(kGray+1);
    histoTrackMI1->SetMarkerColor(kGray+1);
    histoTrackMI1->SetMarkerStyle(24);
    histoTrackMI1->Draw("same");
    histoTrackMI2->SetLineColor(kGray+2);
    histoTrackMI2->SetMarkerColor(kGray+2);
    histoTrackMI2->SetMarkerStyle(26);
    histoTrackMI2->Draw("same");
    histoTrackMI5->SetLineColor(4);
    histoTrackMI5->SetMarkerColor(4);
    histoTrackMI5->SetMarkerStyle(29);
    histoTrackMI5->Draw("same");
    histoTrackMI6->SetLineColor(kBlue+1);
    histoTrackMI6->SetMarkerColor(kBlue+1);
    histoTrackMI6->SetMarkerStyle(30);
    histoTrackMI6->SetTitle("Fraction of rec. tracks with points in ITS");
    histoTrackMI6->Draw("same");
    histoTrackMI3->GetYaxis()->SetTitle("Fraction of Global Tracks with Points in ITS Layers");
    histoTrackMI3->GetXaxis()->SetTitle("run number");
    TLegend* leg3MI=new TLegend(0.7,0.15,0.88,0.35);
    TLegendEntry* entMI;
    entMI=leg3MI->AddEntry(histoTrackMI1,"Layer1","PL");
    entMI->SetTextColor(histoTrackMI1->GetMarkerColor());
    entMI=leg3MI->AddEntry(histoTrackMI2,"Layer2","PL");
    entMI->SetTextColor(histoTrackMI2->GetMarkerColor());
    entMI=leg3MI->AddEntry(histoTrackMI3,"Layer3","PL");
    entMI->SetTextColor(histoTrackMI3->GetMarkerColor());
    entMI=leg3MI->AddEntry(histoTrackMI4,"Layer4","PL");
    entMI->SetTextColor(histoTrackMI4->GetMarkerColor());
    entMI=leg3MI->AddEntry(histoTrackMI5,"Layer5","PL");
    entMI->SetTextColor(histoTrackMI5->GetMarkerColor());
    entMI=leg3MI->AddEntry(histoTrackMI6,"Layer6","PL");
    entMI->SetTextColor(histoTrackMI6->GetMarkerColor());
    
    leg3MI->SetFillStyle(0);
    leg3MI->Draw();
    cMI->SaveAs("TrackPointsMI_trend.pdf");
    //    pdfFileNames+=" TrackPointsMI_trend.pdf";
    cMI->Update();
}
  
    //--------  Draw histograms of fraction of SA tracks w/ hits in layers -------

TCanvas* cSA;
if(histoTrackSA3->GetEntries()>0 || histoTrackSA4->GetEntries()>0){
    cSA=new TCanvas("cSA"," SA Track With points in ITS");
    histoTrackSA3->Draw();
    histoTrackSA3->SetLineColor(1);
    histoTrackSA3->SetMarkerStyle(20);
    histoTrackSA3->Draw();
    histoTrackSA3->SetMinimum(0.);
    histoTrackSA3->SetMaximum(1.05);
    histoTrackSA4->SetLineColor(2);
    histoTrackSA4->SetMarkerColor(2);
    histoTrackSA4->SetMarkerStyle(22);
    histoTrackSA4->Draw("same");
    histoTrackSA1->SetLineColor(kGray+1);
    histoTrackSA1->SetMarkerColor(kGray+1);
    histoTrackSA1->SetMarkerStyle(24);
    histoTrackSA1->Draw("same");
    histoTrackSA2->SetLineColor(kGray+2);
    histoTrackSA2->SetMarkerColor(kGray+2);
    histoTrackSA2->SetMarkerStyle(26);
    histoTrackSA2->Draw("same");
    histoTrackSA5->SetLineColor(4);
    histoTrackSA5->SetMarkerColor(4);
    histoTrackSA5->SetMarkerStyle(29);
    histoTrackSA5->Draw("same");
    histoTrackSA6->SetLineColor(kBlue+1);
    histoTrackSA6->SetMarkerColor(kBlue+1);
    histoTrackSA6->SetMarkerStyle(30);
    histoTrackSA6->SetTitle("Fraction of rec. tracks with points in ITS");
    histoTrackSA6->Draw("same");
    histoTrackSA3->GetYaxis()->SetTitle("SA Tracks with Points in ITS Layers");
    histoTrackSA3->GetXaxis()->SetTitle("run number");
    TLegend* leg3SA=new TLegend(0.7,0.15,0.88,0.35);
    TLegendEntry* entSA;
    entSA=leg3SA->AddEntry(histoTrackSA1,"Layer1","PL");
    entSA->SetTextColor(histoTrackSA1->GetMarkerColor());
    entSA=leg3SA->AddEntry(histoTrackSA2,"Layer2","PL");
    entSA->SetTextColor(histoTrackSA2->GetMarkerColor());
    entSA=leg3SA->AddEntry(histoTrackSA3,"Layer3","PL");
    entSA->SetTextColor(histoTrackSA3->GetMarkerColor());
    entSA=leg3SA->AddEntry(histoTrackSA4,"Layer4","PL");
    entSA->SetTextColor(histoTrackSA4->GetMarkerColor());
    entSA=leg3SA->AddEntry(histoTrackSA5,"Layer5","PL");
    entSA->SetTextColor(histoTrackSA5->GetMarkerColor());
    entSA=leg3SA->AddEntry(histoTrackSA6,"Layer6","PL");
    entSA->SetTextColor(histoTrackSA6->GetMarkerColor());
    
    leg3SA->SetFillStyle(0);
    leg3SA->Draw();
    cSA->SaveAs("TrackPointsSA_trend.pdf");
    //    pdfFileNames+=" TrackPointsSA_trend.pdf";
    cSA->Update();
}

    //--------  Draw histograms of fraction of global tracks w/ point in layers, events w/ SDD in trigger cluster -------
//
TCanvas* c5;
if(histoTrackClu3->GetEntries()>0){
    c5=new TCanvas("c5","Track With points in ITS");
    histoTrackClu3->Draw();
    histoTrackClu3->SetLineColor(1);
    histoTrackClu3->SetMarkerStyle(20);
    histoTrackClu3->Draw();
    histoTrackClu3->SetMinimum(0.);
    histoTrackClu3->SetMaximum(1.05);
    histoTrackClu4->SetLineColor(2);
    histoTrackClu4->SetMarkerColor(2);
    histoTrackClu4->SetMarkerStyle(22);
    histoTrackClu4->Draw("same");
    histoTrackClu1->SetLineColor(kGray+1);
    histoTrackClu1->SetMarkerColor(kGray+1);
    histoTrackClu1->SetMarkerStyle(24);
    histoTrackClu1->Draw("same");
    histoTrackClu2->SetLineColor(kGray+2);
    histoTrackClu2->SetMarkerColor(kGray+2);
    histoTrackClu2->SetMarkerStyle(26);
    histoTrackClu2->Draw("same");
    histoTrackClu5->SetLineColor(4);
    histoTrackClu5->SetMarkerColor(4);
    histoTrackClu5->SetMarkerStyle(29);
    histoTrackClu5->Draw("same");
    histoTrackClu6->SetLineColor(kBlue+1);
    histoTrackClu6->SetMarkerColor(kBlue+1);
    histoTrackClu6->SetMarkerStyle(30);
    histoTrackClu6->SetTitle("Fraction of tracks with points in ITS (SDD in TC)");
    histoTrackClu6->Draw("same");
    histoTrackClu3->GetYaxis()->SetTitle("Fraction of Tracks with Points in ITS Layers");
    histoTrackClu3->GetXaxis()->SetTitle("run number");
    TLatex* txt=new TLatex(0.15,0.85,"Only events with SDD in read out Trigger Cluster");
    txt->SetNDC();
    txt->SetTextColor(1);
    txt->SetTextSize(0.05);
    txt->Draw();
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
    c5->SaveAs("TrackPoints_trend.pdf");
    //    pdfFileNames+=" TrackPoints_trend.pdf";
    c5->Update();
}

    //--------  Draw histograms of number of events per run -------

TCanvas* c22;
if(histonEvents->GetEntries()>0){
    c22=new TCanvas("c22","Number of Events in Run",1000,600);
    histonEvents->SetTitle("Run Events number");
    histonEvents->Draw();
    histonEvents->GetYaxis()->SetTitleOffset(1.2);
    histonEvents->GetYaxis()->SetTitle("number of events");
    histonEvents->GetXaxis()->SetTitle("run number");
    histonEventsTriggered->Draw("same");
    c22->SaveAs("RunEvents_trend.pdf");
    //    pdfFileNames+=" RunEvents_trend.pdf";
    c22->Update();
}
//

    //--------  Draw histograms of fraction of active modules in layers 3 and 4 -------

TCanvas* cfrac;
if(histoFracDead3->GetEntries()>0){
    cfrac=new TCanvas("cfrac","Fraction of SDD good anodes",900,900);
    cfrac->Divide(1,3);
    cfrac->cd(1);
    histoFracDead3->SetMinimum(0.0);
    histoFracDead3->SetMaximum(1.2);
    histoFracDead3->SetMarkerStyle(20);
    histoFracDead3->SetMarkerColor(kOrange+1);
    histoFracDead3->SetLineColor(kOrange+1);
    histoFracDead3->GetYaxis()->SetRangeUser(0.,1.2);
    histoFracDead3->Draw();
    histoFracDead3->GetYaxis()->SetTitle("Fraction of good anodes");
    histoFracDead3->GetXaxis()->SetTitle("run number");
    TLatex* tf3=new TLatex(0.2,0.8,"SDD good anodes - Layer 3 (total: 84*512)");
    tf3->SetNDC();
    tf3->SetTextColor(kOrange+1);
    tf3->Draw();
    histoFracDead4->SetMinimum(0.0);
    histoFracDead4->SetMaximum(1.2);
    histoFracDead4->SetMarkerStyle(20);
    histoFracDead4->SetMarkerColor(kAzure+1);
    histoFracDead4->SetLineColor(kAzure+1);
    histoFracDead4->GetYaxis()->SetRangeUser(0.,1.2);
    histoFracDead4->Draw("same");
    histoFracDead4->GetYaxis()->SetTitle("Fraction of good anodes");
    histoFracDead4->GetXaxis()->SetTitle("run number");
    TLatex* tf4=new TLatex(0.2,0.5,"SDD good anodes - Layer 4 (total: 176*512)");
    tf4->SetNDC();
    tf4->SetTextColor(kAzure+1);
    tf4->Draw();
    TLegend* legSDD=new TLegend(0.90,0.85,0.95,1.00);
    legSDD->SetFillColor(kWhite);
    legSDD->SetFillStyle(1001);
    TLegendEntry* entSDD;
    entSDD=legSDD->AddEntry(histoFracDead3,"Layer 3","PL");
    entSDD->SetTextColor(histoFracDead3->GetMarkerColor());
    entSDD=legSDD->AddEntry(histoFracDead4,"Layer 4","PL");
    entSDD->SetTextColor(histoFracDead4->GetMarkerColor());
    legSDD->Draw();
    
    cfrac->cd(2);
    histoFlagON3->SetMarkerStyle(20);
    histoFlagON3->SetMarkerColor(kOrange+1);
    histoFlagON3->SetLineColor(kOrange+1);
    histoFlagON3->SetMinimum(0.0);
    histoFlagON3->SetMinimum(1.6);
    histoFlagON3->Draw();
    histoFlagON3->GetYaxis()->SetTitle("Status flag");
    histoFlagON3->GetXaxis()->SetTitle("run number");
    TLatex* tf3_1=new TLatex(0.2,0.40,"Status flag Layer 3: ON>0.8? 1=OK, 0=BAD");
    tf3_1->SetNDC();
    tf3_1->SetTextColor(1);
    tf3_1->Draw();
    cfrac->cd(3);
    histoFlagON4->SetMarkerStyle(20);
    histoFlagON4->SetMarkerColor(kAzure+1);
    histoFlagON4->SetLineColor(kAzure+1);
    histoFlagON4->SetMinimum(0.0);
    histoFlagON4->SetMinimum(1.6);
    histoFlagON4->Draw();
    histoFlagON4->GetYaxis()->SetTitle("Status flag");
    histoFlagON4->GetXaxis()->SetTitle("run number");
    TLatex* tf4_1=new TLatex(0.2,0.40,"Status flag Layer 4: ON>0.75? 1=OK, 0=BAD");
    tf4_1->SetNDC();
    tf4_1->SetTextColor(1);
    tf4_1->Draw();
    
    cfrac->SaveAs("SDDanodesON_trend.pdf");
    //    pdfFileNames+=" SDDanodesON_trend.pdf";
    cfrac->Update();
}

    //--------  Draw histograms of SDD drift time and dEdx, SSD dEdx -------

TCanvas* c2;
if(histodEdxLay5->GetEntries()>0){
    c2=new TCanvas("c2","SDD DriftTime & Charge",1200,800);
    c2->Divide(2,2);
    c2->cd(1);
    histominTime->Draw();
    histominTime->SetMinimum(450);
    histominTime->SetMaximum(550);
    histominTime->GetYaxis()->SetTitle("Minimum Drift Time (ns)");
    histominTime->GetXaxis()->SetTitle("run number");
    TLatex* td1=new TLatex(0.2,0.8,"SDD minimum drift time (ref=495ns)");
    td1->SetNDC();
    td1->SetTextColor(1);
    td1->Draw();
    c2->cd(2);
    histomeanTime->Draw();
    histomeanTime->SetMinimum(3150);
    histomeanTime->SetMaximum(3350);
    histomeanTime->GetYaxis()->SetTitle("Average Drift Time (ns)");
    histomeanTime->GetXaxis()->SetTitle("run number");
    TLatex* td2=new TLatex(0.2,0.85,"SDD average drift time");
    td2->SetNDC();
    // td2>SetTextColor(1);
    td2->Draw();
    
    c2->cd(3);
    gPad->SetGridy();
    histodEdxTB0->SetLineColor(1);
    histodEdxTB0->SetMarkerStyle(20);
    histodEdxTB0->Draw();
    histodEdxTB0->SetMinimum(70.);
    histodEdxTB0->SetMaximum(100.);
    histodEdxTB5->SetLineColor(4);
    histodEdxTB5->SetMarkerColor(4);
    histodEdxTB5->SetMarkerStyle(23);
    // histodEdxTB5->SetMinimum(90);
    // histodEdxTB5->SetMaximum(120);
    histodEdxTB5->Draw("same");
    histodEdxTB0->GetYaxis()->SetTitle("MPV of dE/dx (keV/300 #mum)");
    histodEdxTB0->GetXaxis()->SetTitle("run number");
    TLegend* leg2=new TLegend(0.6,0.15,0.88,0.35);
    TLegendEntry* ent=leg2->AddEntry(histodEdxTB0,"Small drift time","PL");
    ent=leg2->AddEntry(histodEdxTB5,"Large drift time","PL");
    ent->SetTextColor(histodEdxTB5->GetMarkerColor());
    leg2->SetFillStyle(0);
    leg2->Draw();
    TLatex* tc1=new TLatex(0.2,0.85,"SDD charge in different drift regions");
    tc1->SetNDC();
    tc1->SetTextColor(1);
    tc1->Draw();
    c2->cd(4);
    gPad->SetGridy();
    histodEdxLay3->SetLineColor(1);
    histodEdxLay3->SetMarkerStyle(20);
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
    histodEdxLay6->SetLineColor(9);
    histodEdxLay6->SetMarkerColor(9);
    histodEdxLay6->SetMarkerStyle(24);
    histodEdxLay6->Draw("same");
    
    histodEdxLay3->GetYaxis()->SetTitle("MPV of dE/dx (keV/300 #mum)");
    histodEdxLay3->GetXaxis()->SetTitle("run number");
    
    
    TLegend* leg2b=new TLegend(0.6,0.15,0.88,0.35);
    ent=leg2b->AddEntry(histodEdxLay3,"Layer 3","PL");
    ent=leg2b->AddEntry(histodEdxLay4,"Layer 4","PL");
    ent=leg2b->AddEntry(histodEdxLay5,"Layer 5","PL");
    ent=leg2b->AddEntry(histodEdxLay6,"Layer 6","PL");
    ent->SetTextColor(histodEdxLay4->GetMarkerColor());
    leg2b->SetFillStyle(0);
    leg2b->Draw();
    TLatex* tc2=new TLatex(0.2,0.85,"SDD and SSD charge in different layers");
    tc2->SetNDC();
    tc2->SetTextColor(1);
    tc2->Draw();
    c2->SaveAs("SDD_SSD_drift_charge_trend.pdf");
    //  pdfFileNames+=" SDD_SSD_drift_charge_trend.pdf";
    c2->Update();
}

    //--------  Draw histograms of status flag for SDD tdrift&dEdx, SSD CR and dEdx-------

TCanvas* c2b;
if(histoFlagCR5->GetEntries()>0){
    c2b=new TCanvas("c2b","SDD & SSD ALARM FLAGS",1200,800);
    c2b->Divide(2,2);
    c2b->cd(1);
    histoFlagminT->SetLineColor(2);
    histoFlagminT->SetMarkerColor(2);
    histoFlagminT->SetMarkerStyle(20);
    histoFlagminT->SetMarkerSize(0.5);
    histoFlagminT->Draw();
    histoFlagminT->SetMinimum(0);
    histoFlagminT->SetMaximum(2.5);
//    histoFlagminT->GetYaxis()->SetTitle("Minimum Drift Time Alarms");
    histoFlagminT->GetXaxis()->SetTitle("run number");
    
    c2b->cd(1);
    histoFlagmeanT->SetLineColor(4);
    histoFlagmeanT->SetMarkerColor(4);
    histoFlagmeanT->SetMarkerSize(0.5);
    histoFlagmeanT->SetMarkerStyle(21);
    histoFlagmeanT->Draw("same");
//    histoFlagmeanT->GetYaxis()->SetTitle("Mean Drift Time Alarms");
    histoFlagmeanT->GetXaxis()->SetTitle("run number");
    TLatex* td1=new TLatex(0.2,0.45,"SDD: Min Time>485 ns: OK=1, ALARM=0");
    td1->SetNDC();
    td1->SetTextColor(1);
    td1->Draw();
    TLatex* td1a=new TLatex(0.2,0.8,"SDD: Mean Time>3150 ns: OK=2, ALARM=0");
    td1a->SetNDC();
    td1a->SetTextColor(1);
    td1a->Draw();
    
    c2b->cd(2);
    histoFlagMPV3->SetMarkerColor(1);
    histoFlagMPV3->SetMarkerStyle(20);
    histoFlagMPV3->SetLineColor(1);
    histoFlagMPV3->SetMarkerSize(0.5);
    histoFlagMPV3->Draw();
    histoFlagMPV3->SetMinimum(0);
    histoFlagMPV3->SetMaximum(2.5);
//    histoFlagMPV3->GetYaxis()->SetTitle("dEdx Layer 3 - MPV alarm ");
    histoFlagMPV3->GetXaxis()->SetTitle("run number");

    histoFlagMPV4->SetMarkerColor(kBlue);
    histoFlagMPV4->SetMarkerStyle(23);
    histoFlagMPV4->SetLineColor(kBlue);
    histoFlagMPV4->SetMarkerSize(0.5);
    histoFlagMPV4->Draw("same");
    histoFlagMPV4->SetMinimum(-0.5);
    histoFlagMPV4->SetMaximum(1.5);
//    histoFlagMPV4->GetYaxis()->SetTitle("dEdx Layer4  - MPV alarm");
    histoFlagMPV4->GetXaxis()->SetTitle("run number");
   TLatex* td1b=new TLatex(0.2,0.45,"SDD Layer 3 80<dEdx<86: OK=1, ALARM=0");
    td1b->SetNDC();
    td1b->SetTextColor(1);
    td1b->Draw();
    TLatex* td1c=new TLatex(0.2,0.8,"SDD Layer 4 80<dEdx4<86: OK=2, ALARM=0");
    td1c->SetNDC();
    td1c->SetTextColor(1);
    td1c->Draw();

    c2b->cd(3);
    histoFlagCR5->SetMarkerColor(kMagenta+2);
    histoFlagCR5->SetMarkerStyle(20);
    histoFlagCR5->SetLineColor(kMagenta+2);
    histoFlagCR5->SetMarkerSize(0.5);
    histoFlagCR5->Draw();
    histoFlagCR5->SetMinimum(0);
    histoFlagCR5->SetMaximum(2.5);
    histoFlagCR5->GetYaxis()->SetTitle("Charge Ratio Layer5 & Layer6 Alarms");
    histoFlagCR5->GetXaxis()->SetTitle("run number");
    histoFlagCR6->SetMarkerColor(9);
    histoFlagCR6->SetMarkerStyle(22);
    histoFlagCR6->SetLineColor(9);
    histoFlagCR6->SetMarkerSize(0.5);
    histoFlagCR6->Draw("same");
    TLatex* td1bb=new TLatex(0.2,0.45,"SSD Layer 5 -0.01<CR<0.01: OK=1, ALARM=0");
    td1bb->SetNDC();
    td1bb->SetTextColor(1);
    td1bb->Draw();
    TLatex* td1cb=new TLatex(0.2,0.8,"SSD Layer 6 -0.01<CR<0.01: OK=2, ALARM=0");
    td1cb->SetNDC();
    td1cb->SetTextColor(1);
    td1cb->Draw();
    
    c2b->cd(4);
    histoFlagMPV5->SetMarkerColor(kMagenta+2);
    histoFlagMPV5->SetMarkerStyle(20);
    histoFlagMPV5->SetLineColor(kMagenta+2);
    histoFlagMPV5->SetMarkerSize(0.5);
    histoFlagMPV5->Draw();
    histoFlagMPV5->SetMinimum(0);
    histoFlagMPV5->SetMaximum(2.5);
    histoFlagMPV5->GetYaxis()->SetTitle("dEdx MPV Layer5 & Layer6 Alarms");
    histoFlagMPV5->GetXaxis()->SetTitle("run number");
    histoFlagMPV6->SetMarkerColor(9);
    histoFlagMPV6->SetMarkerStyle(22);
    histoFlagMPV6->SetLineColor(9);
    histoFlagMPV6->SetMarkerSize(0.5);
    histoFlagMPV6->Draw("same");
    TLatex* td1b2=new TLatex(0.2,0.45,"SSD Layer 5 82<dEdx<83: OK=1, ALARM=0");
    td1b2->SetNDC();
    td1b2->SetTextColor(1);
    td1b2->Draw();
    TLatex* td1c2=new TLatex(0.2,0.8,"SSD Layer 6 82<dEdx<83: OK=2, ALARM=0");
    td1c2->SetNDC();
    td1c2->SetTextColor(1);
    td1c2->Draw();
    
    c2b->SaveAs("SDDSSD_alarm_trend.pdf");
    //  pdfFileNames+=" SDD_SSD_drift_charge_trend.pdf";
    c2b->Update();
}

    //--------  Draw histograms of SSD CR layers 5 and 6 -------

TCanvas *c7;
if(histoChargeRatioLay5->GetEntries()>0){
    c7=new TCanvas("c7","Charge ratio");
    c7->cd();
    histoChargeRatioLay5->SetLineColor(6);
    histoChargeRatioLay5->SetMarkerColor(6);
    histoChargeRatioLay5->SetMarkerStyle(20);
    histoChargeRatioLay5->SetMinimum(-0.01);
    histoChargeRatioLay5->SetMaximum(+0.01);
    histoChargeRatioLay5->Draw();
    histoChargeRatioLay6->SetLineColor(9);
    histoChargeRatioLay6->SetMarkerColor(9);
    histoChargeRatioLay6->SetMarkerStyle(22);
    histoChargeRatioLay5->GetYaxis()->SetTitle("SSD charge ratio");
    histoChargeRatioLay5->GetXaxis()->SetTitle("run number");
    histoChargeRatioLay6->Draw("same");
    TLegend* legCR=new TLegend(0.7,0.65,0.88,0.75);
    TLegendEntry* ent=legCR->AddEntry(histoChargeRatioLay5,"Layer5","PL");
    ent->SetTextColor(histoChargeRatioLay5->GetMarkerColor());
    ent=legCR->AddEntry(histoChargeRatioLay6,"Layer6","PL");
    ent->SetTextColor(histoChargeRatioLay6->GetMarkerColor());
    legCR->SetFillStyle(0);
    legCR->Draw();
    TLatex* tc3=new TLatex(0.2,0.85,"SSD charge ratio in different layers");
    tc3->SetNDC();
    tc3->SetTextColor(1);
    tc3->Draw();
    c7->SaveAs("SSD_chargeratio_trend.pdf");
    //  pdfFileNames+=" SSD_chargeratio_trend.pdf";
}
//

    //--------  Draw histograms of fraction of SSD n/p bad strips layers 5&6-------

TCanvas *c8;
if(histoFracBadn5->GetEntries()>0){
    c8=new TCanvas("c8","Fraction of SSD bad strips",900,900);
    c8->Divide(1,2);
    c8->cd(1);
    histoFracBadn5->SetLineColor(6);
    histoFracBadn5->SetMarkerColor(6);
    histoFracBadn5->SetMarkerStyle(20);
    histoFracBadn5->SetMinimum(0);
    histoFracBadn5->SetMaximum(0.2);
    histoFracBadn5->GetXaxis()->SetTitle("run number");
    histoFracBadn5->Draw();
    histoFracBadp5->SetLineColor(kMagenta+2);
    histoFracBadp5->SetMarkerColor(kMagenta+2);
    histoFracBadp5->SetMarkerStyle(21);
    histoFracBadp5->SetMinimum(0);
    histoFracBadp5->SetMaximum(0.2);
    histoFracBadp5->Draw("same");
    TLegend* legBad=new TLegend(0.7,0.65,0.88,0.75);
    TLegendEntry* ent=legBad->AddEntry(histoFracBadn5,"Layer5 n-side","PL");
    ent->SetTextColor(histoFracBadn5->GetMarkerColor());
    ent=legBad->AddEntry(histoFracBadp5,"Layer5 p-side","PL");
    ent->SetTextColor(histoFracBadp5->GetMarkerColor());
    legBad->SetFillStyle(0);
    legBad->Draw();
    TLatex* tcbad5=new TLatex(0.2,0.85,"SSD bad strips fractions Layer 5");
    tcbad5->SetNDC();
    tcbad5->SetTextColor(1);
    tcbad5->Draw();
    
    c8->cd(2);
    histoFracBadn6->SetLineColor(9);
    histoFracBadn6->SetMarkerColor(9);
    histoFracBadn6->SetMarkerStyle(20);
    histoFracBadn6->SetMinimum(0);
    histoFracBadn6->SetMaximum(0.2);
    histoFracBadn6->GetXaxis()->SetTitle("run number");
    histoFracBadn6->Draw();
    histoFracBadp6->SetLineColor(38);
    histoFracBadp6->SetMarkerColor(38);
    histoFracBadp6->SetMarkerStyle(21);
    histoFracBadp6->SetMinimum(0);
    histoFracBadp6->SetMaximum(0.2);
    histoFracBadp6->Draw("same");
    TLegend* legBad2=new TLegend(0.7,0.65,0.88,0.75);
    //    TLegendEntry* ent;
    ent=legBad2->AddEntry(histoFracBadn6,"Layer6 n-side","PL");
    ent->SetTextColor(histoFracBadn6->GetMarkerColor());
    ent=legBad2->AddEntry(histoFracBadp6,"Layer6 p-side","PL");
    ent->SetTextColor(histoFracBadp6->GetMarkerColor());
    legBad2->SetFillStyle(0);
    legBad2->Draw();
    TLatex* tcbad6=new TLatex(0.2,0.85,"SSD bad strips fractions Layer 6");
    tcbad6->SetNDC();
    tcbad6->SetTextColor(1);
    tcbad6->Draw();
    c8->SaveAs("SSD_BadStripsFrac_trend.pdf");
}

    //--------  Draw histograms of status flag for SSD bad n/p strips fractions -------

    TCanvas *c8b;
if(histoFlagON5n->GetEntries()>0){
    c8b=new TCanvas("c8b","Automatic QA variables for SSD",900,900);
    c8b->Divide(1,2);
    c8b->cd(1);
    
    histoFlagON5n->SetLineColor(6);
    histoFlagON5n->SetMarkerColor(6);
    histoFlagON5n->SetMarkerStyle(20);
    histoFlagON5n->SetMinimum(0.);
    histoFlagON5n->SetMaximum(1.6);
    histoFlagON5n->GetXaxis()->SetTitle("run number");
    histoFlagON5n->Draw();
    histoFlagON5p->SetLineColor(kMagenta+2);
    histoFlagON5p->SetMarkerColor(kMagenta+2);
    histoFlagON5p->SetMarkerStyle(21);
    histoFlagON5p->SetMinimum(0);
    histoFlagON5p->SetMaximum(1.6);
    histoFlagON5p->Draw("same");
//    hVarSSD2n->Draw("same");
//    hVarSSD2p->Draw("same");
    TLegend* legBad3=new TLegend(0.7,0.28,0.88,0.48);
    TLegendEntry* ent=legBad3->AddEntry(histoFlagON5n,"Layer5 n-side","PL");
    ent->SetTextColor(histoFlagON5n->GetMarkerColor());
    ent=legBad3->AddEntry(histoFlagON5p,"Layer5 p-side","PL");
    ent->SetTextColor(histoFlagON5p->GetMarkerColor());
    legBad3->SetFillStyle(0);
    legBad3->Draw();
    TLatex* tcbad5a=new TLatex(0.2,0.85,"Status flag Layer 5: OFF<0.2 1,1.2=OK, 0=BAD");
    tcbad5a->SetNDC();
    tcbad5a->SetTextColor(1);
    tcbad5a->Draw();
    
    c8b->cd(2);
    histoFlagON6n->SetLineColor(9);
    histoFlagON6n->SetMarkerColor(9);
    histoFlagON6n->SetMarkerStyle(20);
    histoFlagON6n->SetMinimum(0);
    histoFlagON6n->SetMaximum(1.6);
    histoFlagON6n->GetXaxis()->SetTitle("run number");
    histoFlagON6n->Draw();
    histoFlagON6p->SetLineColor(kMagenta+2);
    histoFlagON6p->SetLineColor(38);
    histoFlagON6p->SetMarkerColor(38);
    histoFlagON6p->SetMarkerStyle(21);
    histoFlagON6p->SetMinimum(0);
    histoFlagON6p->SetMaximum(1.6);
    histoFlagON6p->Draw("same");
    histoFlagON6p->Draw("same");
    histoFlagON6p->Draw("same");
    TLegend* legBad3b=new TLegend(0.7,0.28,0.88,0.48);
    TLegendEntry* entr=legBad3b->AddEntry(histoFlagON6n,"Layer6 n-side","PL");
//    entr=legBad3b->AddEntry(histoFlagON6n,"Layer6 n-side","PL");
    entr->SetTextColor(histoFlagON6p->GetMarkerColor());
    entr=legBad3b->AddEntry(histoFlagON6p,"Layer6 p-side","PL");
    entr->SetTextColor(histoFlagON6p->GetMarkerColor());
    legBad3b->Draw();
    TLatex* tcbad5b=new TLatex(0.2,0.85,"Status flag Layer 6: OFF<0.2 1,1.2=OK, 0=BAD");
    tcbad5b->SetNDC();
    tcbad5b->SetTextColor(1);
    tcbad5b->Draw();
    
    c8b->SaveAs("SSD_auto_trend.pdf");
    c8b->Update();
}

    //--------  Draw histograms TPCITS and TOFITS matching efficiency -------

TCanvas *cpt02;
if(hEff6Pt02->GetEntries()>0 && hEff6Pt02->GetBinContent(1)>0.){
    cpt02=new TCanvas("cpt02","TPC-ITS matching efficiency",1200,1000);
    cpt02->Divide(3,2);
    cpt02->cd(1);
    hEff6Pt02->SetMinimum(0);
    hEff6Pt02->Draw();
    hEff5Pt02->Draw("same");
    hEff4Pt02->Draw("same");
    hEff3Pt02->Draw("same");
    hEff2Pt02->Draw("same");
    hEffSPDPt02->Draw("same");
    hEffoneSPDPt02->Draw("same");
    hEffTOTPt02->Draw("same");
    hEff6Pt02->GetYaxis()->SetRangeUser(0,1.1);
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
    TLatex* tpc1=new TLatex(0.2,0.85,"TPCITS match eff Pt=0.2");
    tpc1->SetNDC();
    tpc1->SetTextColor(1);
    tpc1->Draw();
    
    cpt02->cd(2);
    hEff6Pt1->Draw();
    hEff5Pt1->Draw("same");
    hEff4Pt1->Draw("same");
    hEff3Pt1->Draw("same");
    hEff2Pt1->Draw("same");
    hEffSPDPt1->Draw("same");
    hEffoneSPDPt1->Draw("same");
    hEffTOTPt1->Draw("same");
    hEff6Pt1->GetYaxis()->SetRangeUser(0,1.1);
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
    TLatex* tpc2=new TLatex(0.2,0.85,"TPCITS match eff Pt=1");
    tpc2->SetNDC();
    tpc2->SetTextColor(1);
    tpc2->Draw();
    
    cpt02->cd(3);
    hEff6Pt10->Draw();
    hEff5Pt10->Draw("same");
    hEff4Pt10->Draw("same");
    hEff3Pt10->Draw("same");
    hEff2Pt10->Draw("same");
    hEffSPDPt10->Draw("same");
    hEffoneSPDPt10->Draw("same");
    hEffTOTPt10->Draw("same");
    hEff6Pt10->GetYaxis()->SetRangeUser(0,1.1);
    
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
    
    TLatex* tpc3=new TLatex(0.2,0.85,"TPCITS match eff Pt=10");
    tpc3->SetNDC();
    tpc3->SetTextColor(1);
    tpc3->Draw();
    
    // Plot TOF
    cpt02->cd(4);
    hEffTOF6Pt02->SetMinimum(0);
    hEffTOF6Pt02->Draw();
    hEffTOF5Pt02->Draw("same");
    hEffTOF4Pt02->Draw("same");
    hEffTOF3Pt02->Draw("same");
    hEffTOF2Pt02->Draw("same");
    hEffTOFSPDPt02->Draw("same");
    hEffTOFoneSPDPt02->Draw("same");
    hEffTOFTOTPt02->Draw("same");
    hEffTOF6Pt02->GetYaxis()->SetRangeUser(0,1.1);
    TLegend* lptTOF02=new TLegend(0.9,0.8,1,1);
    lptTOF02->AddEntry(hEffTOF6Pt02,"6 cls","l");
    lptTOF02->AddEntry(hEffTOF5Pt02,"5 cls","l");
    lptTOF02->AddEntry(hEffTOF4Pt02,"4 cls","l");
    lptTOF02->AddEntry(hEffTOF3Pt02,"3 cls","l");
    lptTOF02->AddEntry(hEffTOF2Pt02,"2 cls","l");
    lptTOF02->AddEntry(hEffTOFSPDPt02,"2SPD + any","l");
    lptTOF02->AddEntry(hEffTOFoneSPDPt02,">=1SPD + any","l");
    lptTOF02->AddEntry(hEffTOFTOTPt02,">=2","l");
    lptTOF02->Draw("same");
    TLatex* tof1=new TLatex(0.2,0.85,"TOFITS match eff Pt=0.5");
    tof1->SetNDC();
    tof1->SetTextColor(1);
    tof1->Draw();
    
    cpt02->cd(5);
    hEffTOF6Pt1->Draw();
    hEffTOF5Pt1->Draw("same");
    hEffTOF4Pt1->Draw("same");
    hEffTOF3Pt1->Draw("same");
    hEffTOF2Pt1->Draw("same");
    hEffTOFSPDPt1->Draw("same");
    hEffTOFoneSPDPt1->Draw("same");
    hEffTOFTOTPt1->Draw("same");
    hEffTOF6Pt1->GetYaxis()->SetRangeUser(0,1.1);
    
    TLegend* lptTOF1=new TLegend(0.9,0.8,1,1);
    lptTOF1->AddEntry(hEffTOF6Pt1,"6 cls","l");
    lptTOF1->AddEntry(hEffTOF5Pt1,"5 cls","l");
    lptTOF1->AddEntry(hEffTOF4Pt1,"4 cls","l");
    lptTOF1->AddEntry(hEffTOF3Pt1,"3 cls","l");
    lptTOF1->AddEntry(hEffTOF2Pt1,"2 cls","l");
    lptTOF1->AddEntry(hEffTOFSPDPt1,"2SPD + any","l");
    lptTOF1->AddEntry(hEffTOFoneSPDPt1,">=1SPD + any","l");
    lptTOF1->AddEntry(hEffTOFTOTPt02,">=2","l");
    lptTOF1->Draw("same");
    TLatex* tof2=new TLatex(0.2,0.85,"TOFITS match eff Pt=1");
    tof2->SetNDC();
    tof2->SetTextColor(1);
    tof2->Draw();
    
    cpt02->cd(6);
    hEffTOF6Pt10->Draw();
    hEffTOF5Pt10->Draw("same");
    hEffTOF4Pt10->Draw("same");
    hEffTOF3Pt10->Draw("same");
    hEffTOF2Pt10->Draw("same");
    hEffTOFSPDPt10->Draw("same");
    hEffTOFoneSPDPt10->Draw("same");
    hEffTOFTOTPt10->Draw("same");
    hEffTOF6Pt10->GetYaxis()->SetRangeUser(0,1.1);
    
    TLegend* lptTOF10=new TLegend(0.9,0.8,1,1);
    lptTOF10->AddEntry(hEffTOF6Pt10,"6 cls","l");
    lptTOF10->AddEntry(hEffTOF5Pt10,"5 cls","l");
    lptTOF10->AddEntry(hEffTOF4Pt10,"4 cls","l");
    lptTOF10->AddEntry(hEffTOF3Pt10,"3 cls","l");
    lptTOF10->AddEntry(hEffTOF2Pt10,"2 cls","l");
    lptTOF10->AddEntry(hEffTOFSPDPt10,"2SPD + any","l");
    lptTOF10->AddEntry(hEffTOFoneSPDPt10,">=1SPD + any","l");
    lptTOF10->AddEntry(hEffTOFTOTPt02,">=2","l");
    lptTOF10->Draw("same");
    
    TLatex* tof3=new TLatex(0.2,0.85,"TOFITS match eff Pt=10");
    tof3->SetNDC();
    tof3->SetTextColor(1);
    tof3->Draw();
    
    cpt02->SaveAs("TPCTOFMatch_trend.pdf");
    //    pdfFileNames+=" TPCTOFMatch_trend.pdf";
}

    //--------  Draw histogram of fraction of SPD pileup events -------

TCanvas *cPileUp;
if(hpileupSPD->GetEntries()>0){
    cPileUp=new TCanvas("cPileUp","cPileUp",1200,800);
    cPileUp->Divide(1,2);
    cPileUp->cd(1);
    hpileupSPD->SetMinimum(0.0);
    hpileupSPD->SetMaximum(1.0);
    hpileupSPD->GetYaxis()->SetTitle("Events with pileup / Events with vertex");
    hpileupSPD->GetXaxis()->SetTitle("run number");
    hpileupSPD->Draw();
    TLatex* tpil=new TLatex(0.2,0.80,"Fraction of pile-up events");
    tpil->SetNDC();
    tpil->SetTextColor(1);
    tpil->Draw();
    cPileUp->cd(2);
    hEff6Pt1->Draw();
    hEff5Pt1->Draw("same");
    hEff4Pt1->Draw("same");
    hEff3Pt1->Draw("same");
    hEff2Pt1->Draw("same");
    hEffSPDPt1->Draw("same");
    hEffoneSPDPt1->Draw("same");
    hEffTOTPt1->Draw("same");
    hEff6Pt1->GetYaxis()->SetRangeUser(0,1.1);
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
    TLatex* tpc2=new TLatex(0.2,0.75,"TPCITS match eff Pt=1");
    tpc2->SetNDC();
    tpc2->SetTextColor(1);
    tpc2->Draw();
    cPileUp->SaveAs("Pileup_trend.pdf");
    //    pdfFileNames+=" Pileup_trend.pdf";
}

    gStyle->SetOptStat(0);
    
    //--------  Draw histograms on pileup vertices and SPD1 cluster number  -------

TCanvas* cpu;
    if(hNumPilVtx->GetEntries()>0){
        cpu=new TCanvas("cPileUpVtx","Pileup vertices from SPD",1200,800);
        cpu->Divide(1,2);
        cpu->cd(1);
        hNumPilVtx->Draw();
        cpu->cd(2);
        hNumPilCL1->SetMaximum(150.);
        hNumPilCL1->Draw();
        hNumNoPilCL1->Draw("same");
        TLegend* legpu=new TLegend(0.5,0.70,0.88,0.85);
        //    TLegend* legpu=new TLegend(0.,0.,0.,0.);
        TLegendEntry* entpu;
        entpu=legpu->AddEntry(hNumPilCL1,"pileup tagged prim. vtx SPD1 clusters","PL");
        entpu->SetTextColor(hNumPilCL1->GetMarkerColor());
        entpu=legpu->AddEntry(hNumNoPilCL1,"no pileup tagged prim. vtx SPD1 clusters","PL");
        entpu->SetTextColor(hNumNoPilCL1->GetMarkerColor());
        
        legpu->SetFillStyle(0);
        legpu->Draw();
        
        cpu->SaveAs("PileupVtx_trend.pdf");
        //    pdfFileNames+=" PileupVtx_trend.pdf";
    }

    
    //--------  Draw histograms of fraction of SPD HS ON layers 1 and 2 -------
    TCanvas *cPixel;
    if(histoFracSPD1->GetEntries()>0){
        cPixel=new TCanvas("cPixel","SPD on",800,900);
        cPixel->Divide(1,2);
        cPixel->cd(1);
        histoFracSPD1->SetMaximum(1.2);
        histoFracSPD1->SetMinimum(0);
        histoFracSPD1->Draw("p");
        histoFracSPD2->Draw("same,p");
        
        TLegend* lSPD=new TLegend(0.8,0.8,1,1);
        lSPD->AddEntry(histoFracSPD1,"Frac. SPD1 ON","Pl");
        lSPD->AddEntry(histoFracSPD2,"Frac. SPD2 ON","Pl");
        lSPD->Draw();
        TLatex* tSPD=new TLatex(0.2,0.3,"Fraction of SPD half staves ON");
        tSPD->SetNDC();
        tSPD->SetTextColor(1);
        tSPD->Draw();
        TLatex* tf1a=new TLatex(0.2,0.58,"SPD HS ON - Layer 1 (total: 40)");
        tf1a->SetNDC();
        tf1a->SetTextColor(1);
        tf1a->Draw();
        TLatex* tf4a=new TLatex(0.2,0.78,"SPD HS ON - Layer 2 (total: 80)");
        tf4a->SetNDC();
        tf4a->SetTextColor(1);
        tf4a->Draw();
        
        
        cPixel->cd(2);
        histoFlagON1->SetMaximum(2.5);
        histoFlagON1->SetMinimum(0.0);
        histoFlagON1->Draw("p");
        TLatex* tf4a_3=new TLatex(0.2,0.30,"Status flag Layer 1: 1=OK, 0=ALARM");
        tf4a_3->SetNDC();
        tf4a_3->SetTextColor(1);
        tf4a_3->Draw();
        histoFlagON2->Draw("same,p");
        TLatex* tf4a_4=new TLatex(0.2,0.80,"Status flag Layer 2: 2=OK, 0=ALARM");
        tf4a_4->SetNDC();
        tf4a_4->SetTextColor(1);
        tf4a_4->Draw();
        TLegend* lSPD2=new TLegend(0.8,0.8,1,1);
        lSPD2->AddEntry(histoFlagON1,"Status SPD1 ON","Pl");
        lSPD2->AddEntry(histoFlagON2,"Status SPD2 ON","Pl");
        lSPD2->Draw();

        
        cPixel->SaveAs("Pixel_trend.pdf");
        //  pdfFileNames+=" Pixel_trend.pdf";
        cPixel->Update();
    }
    
    //--------  Draw histograms of number of ITSpureSA tracks/event -------

    TCanvas *c;
    if(h0->GetEntries()>0){
        c=new TCanvas("ITS pure SA tracks","ITS pure SA tracks");
        c->SetLogy();
        c->SetGridy();
        h0->Draw("p");
        h1->Draw("psame");
        h2->Draw("psame");
        c->BuildLegend(0.11,0.15,0.45,0.30);
        TLatex* ti1=new TLatex(0.11,0.40,"ITS pure SA tracks (normalized to number of events)");
        ti1->SetNDC();
        ti1->SetTextColor(1);
        ti1->Draw();
        c->SaveAs("ITSsa_trend.pdf");
        //   pdfFileNames+=" ITSsa_trend.pdf";
    }
    
    //--------  Draw histograms of number of TPCITS tracks/event -------

    TCanvas *c2a;
    if(h3->GetEntries()>0){
        c2a=new TCanvas("ITS+TPC tracks","ITS+TPC tracks");
        c2a->SetLogy();
        c2a->SetGridy();
        h3->Draw("p");
        h4->Draw("psame");
        h5->Draw("psame");
        c2a->BuildLegend(0.11,0.15,0.45,0.30);
        TLatex* ti2=new TLatex(0.11,0.40,"ITS+TPC tracks (normalized to number of events)");
        ti2->SetNDC();
        ti2->SetTextColor(1);
        ti2->Draw();
        c2a->SaveAs("ITSTPC_trend.pdf");
        //    pdfFileNames+=" ITSTPC_trend.pdf";
    }
    
    //--------  Draw histograms of (TPCITS+ITSsa)/ITSpureSA ratio -------

    TCanvas *c3;
    if(h8->GetEntries()>0){
        c3=new TCanvas("(ITSTPC+ITSsa)/ITSpureSA","(ITSTPC+ITSsa)/ITSpureSA");
        c3->SetGridy();
        h8->Draw("p");
        h7->Draw("psame");
        h6->Draw("psame");
        c3->BuildLegend();
        TLatex* ti3=new TLatex(0.10,0.8,"(ITSTPC+ITSsa)/ITSpureSA ");
        ti3->SetNDC();
        ti3->Draw();
        c3->SaveAs("tracks_ratio_trend.pdf");
        //    pdfFileNames+=" tracks_ratio_trend.pdf";
    }
    
    //--------  Draw histograms of pureSA tracks mean cluster number (N>=3) -------

    TCanvas* cSAb;
    if(h9->GetEntries()>0 && h9->GetBinContent(1)>0.){
        cSAb=new TCanvas("cSAb","Mean cluster number (N>3) ITSsa");
        h9->SetMarkerStyle(20);
        h9->SetMarkerColor(2);
        TLatex* tl = new TLatex(0.2,0.85,"pureSA tracks");
        tl->SetNDC();
        tl->SetTextColor(1);
        tl->SetTextSize(0.05);
        h9->Draw();
        tl->Draw();
        cSAb->SaveAs("meanclu_SA_trend.pdf");
        //    pdfFileNames+=" meanclu_SA_trend.pdf";
    }

    //--------  Draw histograms of pureSA dEdx tracks(4cls)/dEdx tracks(3cls) ratio -------

    TCanvas* cSA2;
    if(h10->GetEntries()>0 && h10->GetBinContent(1)>0.){
        cSA2=new TCanvas("cSA2","dedx(4clu)/dedx(3clu) ITSsa");
        cSA2->SetGrid();
        h10->SetMarkerStyle(20);
        h10->SetMarkerColor(3);
        h10->Draw();
        TLatex* tl = new TLatex(0.2,0.85,"pureSA tracks");
        tl->Draw();
        cSA2->SaveAs("dedx4_3_SA_trend.pdf");
        //    pdfFileNames+=" dedx4_3_SA_trend.pdf";
    }

    //--------  Draw histograms of pureSA tracks mean pt (all particles or pions only) -------

    TCanvas* cSA3;
    if(h11->GetEntries()>0 && h11->GetBinContent(1)>0.){
        cSA3=new TCanvas("cSA3","mean Pt ITSsa");
        h11->SetMarkerStyle(20);
        h11->SetMarkerColor(4);
        h11->Draw();
        TLatex* tl = new TLatex(0.2,0.85,"pureSA tracks");
        tl->Draw();
        cSA3->SaveAs("piPt_SA_trend.pdf");
        //    pdfFileNames+=" piPt_SA_trend.pdf";
    }

    //--------  Draw histograms fractions of pureSA tracks with point in layers -------

    TCanvas* cSA4;
    if(h12->GetEntries()>0 && h12->GetBinContent(1)>0.){
        cSA4=new TCanvas("cSA4","Fraction of tracks with clusters in layers");
        h12->GetYaxis()->SetRange(0.,1.);
        h12->SetMarkerStyle(24);
        h12->SetMarkerColor(kGray+1);
        h12->SetLineColor(kGray+1);
        h12->Draw();
        //    tl->Draw();
        h13->SetMarkerStyle(26);
        h13->SetMarkerColor(kGray+2);
        h13->SetLineColor(kGray+2);
        h13->Draw("same");
        h14->SetMarkerStyle(20);
        h14->SetMarkerColor(1);
        h14->SetLineColor(1);
        h14->Draw("same");
        h15->SetMarkerStyle(22);
        h15->SetMarkerColor(2);
        h15->SetLineColor(2);
        h15->Draw("same");
        h16->SetMarkerStyle(29);
        h16->SetMarkerColor(4);
        h16->SetLineColor(4);
        h16->Draw("same");
        h17->SetMarkerStyle(30);
        h17->SetMarkerColor(kBlue+1);
        h17->SetLineColor(kBlue+1);
        h17->Draw("same");
        TLegend* legpSA=new TLegend(0.7,0.15,0.88,0.35);
        TLegendEntry* entpSA;
        entpSA=legpSA->AddEntry(h12,"Layer1","PL");
        entpSA->SetTextColor(h12->GetMarkerColor());
        entpSA=legpSA->AddEntry(h13,"Layer2","PL");
        entpSA->SetTextColor(h13->GetMarkerColor());
        entpSA=legpSA->AddEntry(h14,"Layer3","PL");
        entpSA->SetTextColor(h14->GetMarkerColor());
        entpSA=legpSA->AddEntry(h15,"Layer4","PL");
        entpSA->SetTextColor(h15->GetMarkerColor());
        entpSA=legpSA->AddEntry(h16,"Layer5","PL");
        entpSA->SetTextColor(h16->GetMarkerColor());
        entpSA=legpSA->AddEntry(h17,"Layer6","PL");
        entpSA->SetTextColor(h17->GetMarkerColor());
        
        legpSA->SetFillStyle(0);
        legpSA->Draw();
        
        
        cSA4->SaveAs("Frac_track_SA_trend.pdf");
        //    pdfFileNames+=" Frac_track_SA_trend.pdf";
    }

    // ------------------------  Draw chi2 plots -----------------------------------
    
    TCanvas* cChi2;
    if(hchi2->GetEntries()>0 && hchi2->GetBinContent(1)>0.){
        cChi2=new TCanvas("cChi2","#chi^{2} for TPCITS and ITSpureSA tracks",800, 1200);
        cChi2->Divide(1,2);
        cChi2->cd(1);
        hchi2->SetMinimum(0.0);
        hchi2->SetMaximum(5.0);
        hchi2->Draw("P");
        cChi2->cd(2);
        hchi2SA->SetMinimum(0.0);
        hchi2SA->SetMaximum(2.0);
        hchi2SA->Draw("P");
        
        cChi2->SaveAs("Chi2_tracks_trend.pdf");
        //    pdfFileNames+=" Chi2_tracks_trend.pdf";
    }

    // -------------------- Draw eta-phi TPCITS distributions per layer
/*
    TCanvas* cEtaPhi1;
    if(hOccEta1->GetEntries()>0&&hOccEta1->GetMean()>0){
        cEtaPhi1=new TCanvas("cEtaPhi1","#eta and #phi cluster occupancy, Layer 1",800, 1200);
        cEtaPhi1->Divide(1,2);
        cEtaPhi1->cd(1);
        //    hOccEta1->SetMarkerStyle(20);
        hOccEta1->GetYaxis()->SetNdivisions(2);
        hOccEta1->SetLabelOffset(0.01,"Y");
        hOccEta1->Draw("colz");
        cEtaPhi1->cd(2);
        //    hOccPhi1->SetMarkerStyle(20);
        hOccPhi1->Draw("colz");
        
        cEtaPhi1->SaveAs("Layer1_eta_phi_trend.pdf");
        //    pdfFileNames+=" Chi2_tracks_trend.pdf";
    }
    
    TCanvas* cEtaPhi2;
    if(hOccEta2->GetEntries()>0&&hOccEta2->GetMean()>0){
        cEtaPhi2=new TCanvas("cEtaPhi2","#eta and #phi cluster occupancy, Layer 2",800, 1200);
        cEtaPhi2->Divide(1,2);
        cEtaPhi2->cd(1);
        //    hOccEta2->SetMarkerStyle(20);
        hOccEta2->GetYaxis()->SetNdivisions(2);
        hOccEta2->SetLabelOffset(0.01,"Y");
        hOccEta2->Draw("colz");
        cEtaPhi2->cd(2);
        //    hOccPhi2->SetMarkerStyle(20);
        hOccPhi2->Draw("colz");
        
        cEtaPhi2->SaveAs("Layer2_eta_phi_trend.pdf");
        //    pdfFileNames+=" Chi2_tracks_trend.pdf";
    }
    
    TCanvas* cEtaPhi3;
    if(hOccEta3->GetEntries()>0&&hOccEta3->GetMean()>0){
        cEtaPhi3=new TCanvas("cEtaPhi3","#eta and #phi cluster occupancy, Layer 3",800, 1200);
        cEtaPhi3->Divide(1,2);
        cEtaPhi3->cd(1);
        //    hOccEta3->SetMarkerStyle(20);
        hOccEta3->GetYaxis()->SetNdivisions(2);
        hOccEta3->SetLabelOffset(0.01,"Y");
        hOccEta3->Draw("colz");
        cEtaPhi3->cd(2);
        //    hOccPhi3->SetMarkerStyle(20);
        hOccPhi3->Draw("colz");
        
        cEtaPhi3->SaveAs("Layer3_eta_phi_trend.pdf");
        //    pdfFileNames+=" Chi2_tracks_trend.pdf";
    }
    
    TCanvas* cEtaPhi4;
    if(hOccEta4->GetEntries()>0&&hOccEta4->GetMean()>0){
        cEtaPhi4=new TCanvas("cEtaPhi4","#eta and #phi cluster occupancy, Layer 4",800, 1200);
        cEtaPhi4->Divide(1,2);
        cEtaPhi4->cd(1);
        //    hOccEta4->SetMarkerStyle(20);
        hOccEta4->GetYaxis()->SetNdivisions(2);
        hOccEta4->SetLabelOffset(0.01,"Y");
        hOccEta4->Draw("colz");
        cEtaPhi4->cd(2);
        //    hOccPhi4->SetMarkerStyle(20);
        hOccPhi4->Draw("colz");
        
        cEtaPhi4->SaveAs("Layer4_eta_phi_trend.pdf");
        //    pdfFileNames+=" Chi2_tracks_trend.pdf";
    }
    
    TCanvas* cEtaPhi5;
    if(hOccEta5->GetEntries()>0&&hOccEta5->GetMean()>0){
        cEtaPhi5=new TCanvas("cEtaPhi5","#eta and #phi cluster occupancy, Layer 5",800, 1200);
        cEtaPhi5->Divide(1,2);
        cEtaPhi5->cd(1);
        //    hOccEta5->SetMarkerStyle(20);
        hOccEta5->SetLabelOffset(0.01,"Y");
        hOccEta5->GetYaxis()->SetNdivisions(2);
        hOccEta5->Draw("colz");
        cEtaPhi5->cd(2);
        //    hOccPhi5->SetMarkerStyle(20);
        hOccPhi5->Draw("colz");
        
        cEtaPhi5->SaveAs("Layer5_eta_phi_trend.pdf");
        //    pdfFileNames+=" Chi2_tracks_trend.pdf";
    }
    
    
    TCanvas* cEtaPhi6;
    if(hOccEta6->GetEntries()>0&&hOccEta6->GetMean()>0){
        cEtaPhi6=new TCanvas("cEtaPhi6","#eta and #phi cluster occupancy, Layer 6",800, 1200);
        cEtaPhi6->Divide(1,2);
        cEtaPhi6->cd(1);
        //    hOccEta6->SetMarkerStyle(20);
        hOccEta6->GetYaxis()->SetNdivisions(2);
        hOccEta6->Draw("colz");
        cEtaPhi6->cd(2);
        //    hOccPhi6->SetMarkerStyle(20);
        hOccPhi6->Draw("colz");
        
        cEtaPhi6->SaveAs("Layer6_eta_phi_trend.pdf");
        //    pdfFileNames+=" Chi2_tracks_trend.pdf";
    }
*/
    // ------------------------- Draw ITS PID histos
    
    TCanvas* cPID;
    if(hNsig02->GetEntries()>0 && hNsig02->GetBinContent(1)<100.){
        cPID=new TCanvas("cPID"," ITS PID for TPCITS pion tracks");
        hNsig02->SetMinimum(-2.);
        hNsig02->SetMaximum(1.0);
        hNsig02->SetTitle("ITS n#sigma pions - TPCITS tracks");
        hNsig02->GetXaxis()->SetTitle("run number");
        hNsig02->GetYaxis()->SetTitle("nsigma for #pi");
        
        hNsig02->Draw();
        hNsig05->Draw("same");
        hNsig1->Draw("same");
        hNsig3->Draw("same");
        
        TLegend* legPID=new TLegend(0.7,0.65,0.88,0.85);
        TLegendEntry* entPID;
        entPID=legPID->AddEntry(hNsig02,"0.2 GeV/c","PL");
        entPID->SetTextColor(hNsig02->GetMarkerColor());
        entPID=legPID->AddEntry(hNsig05,"0.5 GeV/c","PL");
        entPID->SetTextColor(hNsig05->GetMarkerColor());
        entPID=legPID->AddEntry(hNsig1,"1.0 GeV/c","PL");
        entPID->SetTextColor(hNsig1->GetMarkerColor());
        entPID=legPID->AddEntry(hNsig3,"2.0 GeV/c","PL");
        entPID->SetTextColor(hNsig3->GetMarkerColor());
        
        legPID->SetFillStyle(0);
        legPID->Draw();
        TLatex* tpid = new TLatex(0.2,0.85,"ITS PID - TPCITS pion tracks");
        tpid->SetNDC();
        tpid->SetTextColor(1);
        tpid->SetTextSize(0.05);
        tpid->Draw();
        
        cPID->SaveAs("ITSPID_trend.pdf");
        //    pdfFileNames+=" TrackPointsMI_trend.pdf";
        cPID->Update();
    }
    
    // ---------------------- Draw TPCITS tracks DCA histos
    
    TCanvas* cDCA;
    if(hmDCA05->GetEntries()>0){
        cDCA=new TCanvas("cDCA"," DCA for TPCITS tracks w/ 6 cls in ITS",800,1200);
        cDCA->Divide(1,2);
        cDCA->cd(1);
        
        hmDCA05->SetMinimum(-50.);
        hmDCA05->SetMaximum(50.);
        hmDCA05->SetTitle("DCA mean, pt 0.55-0.65 GeV/c - TPCITS tracks");
        hmDCA05->GetXaxis()->SetTitle("run number");
        hmDCA05->GetYaxis()->SetTitle("DCA mean (#mum)");
        
        hmDCA05->Draw();
        hmDCA1->Draw("same");
        hmDCA5->Draw("same");
        hmDCA10->Draw("same");
        
        TLegend* legDCA=new TLegend(0.7,0.65,0.88,0.85);
        TLegendEntry* entDCA;
        entDCA=legDCA->AddEntry(hmDCA05,"0.5 GeV/c","PL");
        entDCA->SetTextColor(hmDCA05->GetMarkerColor());
        entDCA=legDCA->AddEntry(hmDCA1,"1.0 GeV/c","PL");
        entDCA->SetTextColor(hmDCA1->GetMarkerColor());
        entDCA=legDCA->AddEntry(hmDCA5,"4.5 GeV/c","PL");
        entDCA->SetTextColor(hmDCA5->GetMarkerColor());
        entDCA=legDCA->AddEntry(hmDCA10,"10.0 GeV/c","PL");
        entDCA->SetTextColor(hmDCA10->GetMarkerColor());
        
        legDCA->SetFillStyle(0);
        legDCA->Draw();
        
        TLatex* tdca = new TLatex(0.2,0.85,"DCA mean - TPCITS tracks w/6 cls in ITS");
        tdca->SetNDC();
        tdca->SetTextColor(1);
        tdca->SetTextSize(0.05);
        tdca->Draw();
        
        cDCA->cd(2);
        hrmsDCA05->SetMinimum(0.);
        hrmsDCA05->SetMaximum(300.);
        hrmsDCA05->SetTitle("DCA RMS, pt 0.55-0.65 GeV/c - TPCITS tracks");
        hrmsDCA05->GetXaxis()->SetTitle("run number");
        hrmsDCA05->GetYaxis()->SetTitle("DCA RMS (#mum)");
        
        hrmsDCA05->Draw();
        hrmsDCA1->Draw("same");
        hrmsDCA5->Draw("same");
        hrmsDCA10->Draw("same");
        
        TLegend* legDCA2=new TLegend(0.7,0.65,0.88,0.85);
        TLegendEntry* entDCA2;
        entDCA2=legDCA2->AddEntry(hrmsDCA05,"0.5 GeV/c","PL");
        entDCA2->SetTextColor(hrmsDCA05->GetMarkerColor());
        entDCA2=legDCA2->AddEntry(hrmsDCA1,"1.0 GeV/c","PL");
        entDCA2->SetTextColor(hrmsDCA1->GetMarkerColor());
        entDCA2=legDCA2->AddEntry(hrmsDCA5,"4.5 GeV/c","PL");
        entDCA2->SetTextColor(hrmsDCA5->GetMarkerColor());
        entDCA2=legDCA2->AddEntry(hrmsDCA10,"10.0 GeV/c","PL");
        entDCA2->SetTextColor(hrmsDCA10->GetMarkerColor());
        
        legDCA2->SetFillStyle(0);
        legDCA2->Draw();
        
        TLatex* tdca2 = new TLatex(0.2,0.85,"DCA RMS - TPCITS tracks w/6 cls in ITS");
        tdca2->SetNDC();
        tdca2->SetTextColor(1);
        tdca2->SetTextSize(0.05);
        tdca2->Draw();
        
        cDCA->SaveAs("DCAtracks_trend.pdf");
        //    pdfFileNames+=" DCAtracks_trend.pdf";
        cDCA->Update();
    }

    // ---------------------- Draw TPCITS tracks DCAz histos
    
    TCanvas* cDCAz;
    if(hmDCAz05->GetEntries()>0){
        cDCAz=new TCanvas("cDCAz"," DCAz for TPCITS tracks w/ 6 cls in ITS",800,1200);
        cDCAz->Divide(1,2);
        cDCAz->cd(1);
        hmDCAz05->SetMinimum(-80.);
        hmDCAz05->SetMaximum(80.);
        hmDCAz05->SetTitle("DCAz mean, pt 0.55-0.65 GeV/c - TPCITS tracks");
        hmDCAz05->GetXaxis()->SetTitle("run number");
        hmDCAz05->GetYaxis()->SetTitle("DCAz mean (#mum)");
        
        hmDCAz05->Draw();
        hmDCAz1->Draw("same");
        hmDCAz5->Draw("same");
        hmDCAz10->Draw("same");
        
        TLegend* legDCAz=new TLegend(0.7,0.65,0.88,0.85);
        TLegendEntry* entDCAz;
        entDCAz=legDCAz->AddEntry(hmDCAz05,"0.5 GeV/c","PL");
        entDCAz->SetTextColor(hmDCAz05->GetMarkerColor());
        entDCAz=legDCAz->AddEntry(hmDCAz1,"1.0 GeV/c","PL");
        entDCAz->SetTextColor(hmDCAz1->GetMarkerColor());
        entDCAz=legDCAz->AddEntry(hmDCAz5,"4.5 GeV/c","PL");
        entDCAz->SetTextColor(hmDCAz5->GetMarkerColor());
        entDCAz=legDCAz->AddEntry(hmDCAz10,"10.0 GeV/c","PL");
        entDCAz->SetTextColor(hmDCAz10->GetMarkerColor());
        
        legDCAz->SetFillStyle(0);
        legDCAz->Draw();
        
        TLatex* tdcaz = new TLatex(0.2,0.85,"DCAz mean - TPCITS tracks w/6 cls in ITS");
        tdcaz->SetNDC();
        tdcaz->SetTextColor(1);
        tdcaz->SetTextSize(0.05);
        tdcaz->Draw();
        
        cDCAz->cd(2);
        hrmsDCAz05->SetMinimum(0.);
        hrmsDCAz05->SetMaximum(400.);
        hrmsDCAz05->SetTitle("DCAz RMS, pt 0.55-0.65 GeV/c - TPCITS tracks");
        hrmsDCAz05->GetXaxis()->SetTitle("run number");
        hrmsDCAz05->GetYaxis()->SetTitle("DCAz RMS (#mum)");
        
        hrmsDCAz05->Draw();
        hrmsDCAz1->Draw("same");
        hrmsDCAz5->Draw("same");
        hrmsDCAz10->Draw("same");
        
        TLegend* legDCAz2=new TLegend(0.7,0.65,0.88,0.85);
        TLegendEntry* entDCAz2;
        entDCAz2=legDCAz2->AddEntry(hrmsDCAz05,"0.5 GeV/c","PL");
        entDCAz2->SetTextColor(hrmsDCAz05->GetMarkerColor());
        entDCAz2=legDCAz2->AddEntry(hrmsDCAz1,"1.0 GeV/c","PL");
        entDCAz2->SetTextColor(hrmsDCAz1->GetMarkerColor());
        entDCAz2=legDCAz2->AddEntry(hrmsDCAz5,"4.5 GeV/c","PL");
        entDCAz2->SetTextColor(hrmsDCAz5->GetMarkerColor());
        entDCAz2=legDCAz2->AddEntry(hrmsDCAz10,"10.0 GeV/c","PL");
        entDCAz2->SetTextColor(hrmsDCAz10->GetMarkerColor());
        
        legDCAz2->SetFillStyle(0);
        legDCAz2->Draw();
        
        TLatex* tdcaz2 = new TLatex(0.2,0.85,"DCAz RMS - TPCITS tracks w/6 cls in ITS");
        tdcaz2->SetNDC();
        tdcaz2->SetTextColor(1);
        tdcaz2->SetTextSize(0.05);
        tdcaz2->Draw();
        
        cDCAz->SaveAs("DCAztracks_trend.pdf");
        //    pdfFileNames+=" DCAtracks_trend.pdf";
        cDCAz->Update();
    }

    /// ------------------------------------ merging .pdf outputs
    
    if(hMeanVx->GetEntries()>0) pdfFileNames+=" Vertex_trend.pdf";
    if(histoTrackMI3->GetEntries()>0 && histoTrackMI3->GetBinContent(1)>0.) pdfFileNames+=" TrackPointsMI_trend.pdf";
    if(h12->GetEntries()>0 && h12->GetBinContent(1)>0.) pdfFileNames+=" Frac_track_SA_trend.pdf";
//    if(histoEvwSDD->GetEntries()>0 && histoEvwSDD->GetBinContent(1)>0.) pdfFileNames+=" NoFast_trend.pdf";
    if(histonEvents->GetEntries()>0 && histonEvents->GetBinContent(1)>0.) pdfFileNames+=" RunEvents_trend.pdf";
    if(histoFracDead3->GetEntries()>0) pdfFileNames+=" SDDanodesON_trend.pdf";
    if(histodEdxLay5->GetEntries()>0) pdfFileNames+=" SDD_SSD_drift_charge_trend.pdf";
    if(histoFlagCR5->GetEntries()>0) pdfFileNames+=" SDDSSD_alarm_trend.pdf";
    if(histoFracBadn5->GetEntries()>0) pdfFileNames+=" SSD_BadStripsFrac_trend.pdf";
    if(histoFlagON5n->GetEntries()>0) pdfFileNames+=" SSD_auto_trend.pdf";
    if(histoChargeRatioLay5->GetEntries()>0) pdfFileNames+=" SSD_chargeratio_trend.pdf";
    if(hNsig02->GetEntries()>0 && hNsig02->GetBinContent(1)<100.) pdfFileNames+=" ITSPID_trend.pdf";
    if(hEff6Pt02->GetEntries()>0 && hEff6Pt02->GetBinContent(1)>0.) pdfFileNames+=" TPCTOFMatch_trend.pdf";
    if(hpileupSPD->GetEntries()>0) pdfFileNames+=" Pileup_trend.pdf";
    if(hNumPilVtx->GetEntries()>0) pdfFileNames+=" PileupVtx_trend.pdf";
    if(histoFracSPD1->GetEntries()>0) pdfFileNames+=" Pixel_trend.pdf";
    if(h0->GetEntries()>0) pdfFileNames+=" ITSsa_trend.pdf";
    if(h3->GetEntries()>0) pdfFileNames+=" ITSTPC_trend.pdf";
    if(h8->GetEntries()>0) pdfFileNames+=" tracks_ratio_trend.pdf";
    if(h9->GetEntries()>0 && h9->GetBinContent(1)>0.) pdfFileNames+=" meanclu_SA_trend.pdf";
    if(h10->GetEntries()>0 && h10->GetBinContent(1)>0.) pdfFileNames+=" dedx4_3_SA_trend.pdf";
    if(h11->GetEntries()>0 && h11->GetBinContent(1)>0.) pdfFileNames+=" piPt_SA_trend.pdf";
    if(hchi2->GetEntries()>0 && hchi2->GetBinContent(1)>0.) pdfFileNames+=" Chi2_tracks_trend.pdf";
    if(hmDCA05->GetEntries()>0) pdfFileNames+=" DCAtracks_trend.pdf";
    if(hmDCAz05->GetEntries()>0) pdfFileNames+=" DCAztracks_trend.pdf";
//    if(hOccEta1->GetEntries()>0 && hOccEta1->GetMean()>0) pdfFileNames+=" Layer1_eta_phi_trend.pdf";
//    if(hOccEta2->GetEntries()>0 && hOccEta2->GetMean()>0) pdfFileNames+=" Layer2_eta_phi_trend.pdf";
//    if(hOccEta3->GetEntries()>0 && hOccEta3->GetMean()>0) pdfFileNames+=" Layer3_eta_phi_trend.pdf";
//    if(hOccEta4->GetEntries()>0 && hOccEta4->GetMean()>0) pdfFileNames+=" Layer4_eta_phi_trend.pdf";
//    if(hOccEta5->GetEntries()>0 && hOccEta5->GetMean()>0) pdfFileNames+=" Layer5_eta_phi_trend.pdf";
//    if(hOccEta6->GetEntries()>0 && hOccEta6->GetMean()>0) pdfFileNames+=" Layer6_eta_phi_trend.pdf";
    
        // ------------------------------- merge the pdf files to mergedITS_trend.pdf
    TString command("gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=merged");
    command=command+"ITS_trend.pdf "+pdfFileNames;
    gSystem->Exec(command.Data());
    printf(" Merging the pdf file:  %s \n",command.Data());

    
    // ---------------------------- merge TCanvases to file: mergedITS_trend.root
    TFile *f2=new TFile("mergedITS_trend.root","RECREATE");
    
    
    if(hMeanVx->GetEntries()>0) cVertexDisto->Write();
    if(histoTrackMI3->GetEntries()>0 && histoTrackMI3->GetBinContent(1)>0.) cMI->Write();
    if(histoTrackSA3->GetEntries()>0) cSA->Write();
    if(histoTrackClu3->GetEntries()>0)c5->Write();
//    if(histoEvwSDD->GetEntries()>0 && histoEvwSDD->GetBinContent(1)>0.) cev->Write();
    if(histonEvents->GetEntries()>0 && histonEvents->GetBinContent(1)>0.) c22->Write();
    if(histoFracDead3->GetEntries()>0) cfrac->Write();
    if(histodEdxLay5->GetEntries()>0) c2->Write();
    if(histoFlagCR5->GetEntries()>0) c2b->Write();
    if(histoFracBadn5->GetEntries()>0) c8->Write();
    if(histoFlagON5n->GetEntries()>0) c8b->Write();
    if(histoChargeRatioLay5->GetEntries()>0) c7->Write();
    if(hEff6Pt02->GetEntries()>0 && hEff6Pt02->GetBinContent(1)>0.) cpt02->Write();
    if(hpileupSPD->GetEntries()>0) cPileUp->Write();
    if(hNumPilVtx->GetEntries()>0) cpu->Write();
    if(histoFracSPD1->GetEntries()>0) cPixel->Write();
    if(h0->GetEntries()>0) c->Write();
    if(h3->GetEntries()>0) c2a->Write();
    if(h8->GetEntries()>0) c3->Write();
    if(h9->GetEntries()>0 && h9->GetBinContent(1)>0.) cSAb->Write();
    if(h10->GetEntries()>0 && h10->GetBinContent(1)>0.) cSA2->Write();
    if(h11->GetEntries()>0 && h11->GetBinContent(1)>0.) cSA3->Write();
    if(h12->GetEntries()>0 && h12->GetBinContent(1)>0.) cSA4->Write();
    if(hchi2->GetEntries()>0 && hchi2->GetBinContent(1)>0.) cChi2->Write();
    if(hmDCA05->GetEntries()>0) cDCA->Write();
    if(hmDCAz05->GetEntries()>0) cDCAz->Write();
//    if(hOccEta1->GetEntries()>0 && hOccEta1->GetMean()>0) cEtaPhi1->Write();
//    if(hOccEta2->GetEntries()>0 && hOccEta2->GetMean()>0) cEtaPhi2->Write();
//    if(hOccEta3->GetEntries()>0 && hOccEta3->GetMean()>0) cEtaPhi3->Write();
//    if(hOccEta4->GetEntries()>0 && hOccEta4->GetMean()>0) cEtaPhi4->Write();
//    if(hOccEta5->GetEntries()>0 && hOccEta5->GetMean()>0) cEtaPhi5->Write();
//    if(hOccEta6->GetEntries()>0 && hOccEta6->GetMean()>0) cEtaPhi6->Write();
    if(hNsig02->GetEntries()>0 && hNsig02->GetBinContent(1)<100.) cPID->Write();
    
    f2->Close();

    
    
  return 0;
}
