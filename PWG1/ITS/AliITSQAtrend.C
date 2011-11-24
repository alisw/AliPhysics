#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TPad.h>
#include <TGraphErrors.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TGrid.h>
#include <TGridResult.h>
#include <TMath.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <Riostream.h>
#include "AliITSgeomTGeo.h"
#endif

TString pdfFileNames="";
void MakePlot(TString ntupleFileName="TrendingITS.root");
void PlotITSSA(TFile *fil,Int_t *myIndex);
void FillITSSAntuple(TFile* f,TNtuple* nt, Int_t nrun);
void AliITSQAtrend(TString runListFile="LHC11hNo.txt",TString ntupleFileName="TrendingITS.root");

////////////////////////////////////////////////////////////////
//   Please, read this comment before using this macro 
//
//   INPUT FILE: a text file (by default LHC11hNo.txt) which contains
//   a list of the complete path+file name of the QAresults root files 
//   (without the alien:// prefix).
//   One file per line. The order is irrelevant.
//
//   USAGE:
//
//   Function AliITSQAtrend():  
//   it looks for a local root file named TrendingITS.root,
//   where the ntuples used to build the trending plots are
//   stored. This file is used to generate incrementally the ntuple contents:
//   when you add a new entry in the LHC11hNo.txt file and you have a
//   local copy of the  TrendingITS.root file, only the additional run will
//   be processed. The whole list is processed only the first time you use the
//   macro. Please, bear in mind that this macro is RAM-intensive: all the
//   ntuples are kept in memory. It is better to add few runs to the list at 
//   each time, according to the RAM of your computer. 
//   The function AliITSQAtrend does not produce any plot.
//
//   Function MakePlot():
//   it produces the plots. For each canvas a PDF file is created.
//   A PDF file with all the canvases merged is also produced
////////////////////////////////////////////////////////////////

/* $Id$ */

void MakePlot(TString ntupleFileName){
  //TrendingITS.root
  TFile* fil=new TFile(ntupleFileName.Data(),"read");
  if(!fil){
    printf("File with ntuple does not exist\n");
    return;
  }
  TNtuple* ntsdd=(TNtuple*)fil->Get("ntsdd");

  Float_t nrun;
  Float_t meanTrPts3,errmeanTrPts3,meanTrPts4,errmeanTrPts4;
  Float_t minDrTime,errminDrTime,meanDrTime,errmeanDrTime;
  Float_t fracTrackWithClu1,fracTrackWithClu2,errfracTrackWithClu1,errfracTrackWithClu2;
  Float_t fracTrackWithClu3,fracTrackWithClu4,errfracTrackWithClu3,errfracTrackWithClu4;
  Float_t fracTrackWithClu5,fracTrackWithClu6,errfracTrackWithClu5,errfracTrackWithClu6;
  Float_t fracExtra,errfracExtra;
  Float_t meandEdxLay3,errmeandEdxLay3,meandEdxLay4,errmeandEdxLay4;
  Float_t meandEdxTB0,errmeandEdxTB0,meandEdxTB5,errmeandEdxTB5;
  Float_t nMod95,nMod80,nMod60,nModEmpty;
  
  ntsdd->SetBranchAddress("nrun",&nrun);
  ntsdd->SetBranchAddress("fracTrackWithClu1",&fracTrackWithClu1);
  ntsdd->SetBranchAddress("errfracTrackWithClu1",&errfracTrackWithClu1);
  ntsdd->SetBranchAddress("fracTrackWithClu2",&fracTrackWithClu2);
  ntsdd->SetBranchAddress("errfracTrackWithClu2",&errfracTrackWithClu2);
  ntsdd->SetBranchAddress("fracTrackWithClu3",&fracTrackWithClu3);
  ntsdd->SetBranchAddress("errfracTrackWithClu3",&errfracTrackWithClu3);
  ntsdd->SetBranchAddress("fracTrackWithClu4",&fracTrackWithClu4);
  ntsdd->SetBranchAddress("errfracTrackWithClu4",&errfracTrackWithClu4);
  ntsdd->SetBranchAddress("fracTrackWithClu5",&fracTrackWithClu5);
  ntsdd->SetBranchAddress("errfracTrackWithClu5",&errfracTrackWithClu5);
  ntsdd->SetBranchAddress("fracTrackWithClu6",&fracTrackWithClu6);
  ntsdd->SetBranchAddress("errfracTrackWithClu6",&errfracTrackWithClu6);
  ntsdd->SetBranchAddress("nMod95",&nMod95);
  ntsdd->SetBranchAddress("nMod80",&nMod80);
  ntsdd->SetBranchAddress("nMod60",&nMod60);
  ntsdd->SetBranchAddress("nModEmpty",&nModEmpty);

  ntsdd->SetBranchAddress("meanTrPts3",&meanTrPts3);
  ntsdd->SetBranchAddress("errmeanTrPts3",&errmeanTrPts3);
  ntsdd->SetBranchAddress("meanTrPts4",&meanTrPts4);
  ntsdd->SetBranchAddress("errmeanTrPts4",&errmeanTrPts4);
  ntsdd->SetBranchAddress("minDrTime",&minDrTime);
  ntsdd->SetBranchAddress("errminDrTime",&errminDrTime);
  ntsdd->SetBranchAddress("meanDrTime",&meanDrTime);
  ntsdd->SetBranchAddress("errmeanDrTime",&errmeanDrTime);
  ntsdd->SetBranchAddress("fracExtra",&fracExtra);
  ntsdd->SetBranchAddress("errfracExtra",&errfracExtra);
  ntsdd->SetBranchAddress("meandEdxTB0",&meandEdxTB0);
  ntsdd->SetBranchAddress("errmeandEdxTB0",&errmeandEdxTB0);
  ntsdd->SetBranchAddress("meandEdxTB5",&meandEdxTB5);
  ntsdd->SetBranchAddress("errmeandEdxTB5",&errmeandEdxTB5);
  ntsdd->SetBranchAddress("meandEdxLay3",&meandEdxLay3);
  ntsdd->SetBranchAddress("errmeandEdxLay3",&errmeandEdxLay3);
  ntsdd->SetBranchAddress("meandEdxLay4",&meandEdxLay4);
  ntsdd->SetBranchAddress("errmeandEdxLay4",&errmeandEdxLay4);

  TH1F* histotrp3=new TH1F("histotrp3","",(Int_t)ntsdd->GetEntries(),0.,ntsdd->GetEntries());
  TH1F* histotrp4=new TH1F("histotrp4","",(Int_t)ntsdd->GetEntries(),0.,ntsdd->GetEntries());
  TH1F* histominTime=new TH1F("histominTime","",(Int_t)ntsdd->GetEntries(),0.,ntsdd->GetEntries());
  TH1F* histomeanTime=new TH1F("histomeanTime","",(Int_t)ntsdd->GetEntries(),0.,ntsdd->GetEntries());
  TH1F* histofracExtra=new TH1F("histofracExtra","",(Int_t)ntsdd->GetEntries(),0.,ntsdd->GetEntries());
  TH1F* histodEdxTB0=new TH1F("histodEdxTB0","",(Int_t)ntsdd->GetEntries(),0.,ntsdd->GetEntries());
  TH1F* histodEdxTB5=new TH1F("histodEdxTB5","",(Int_t)ntsdd->GetEntries(),0.,ntsdd->GetEntries());
  TH1F* histodEdxLay3=new TH1F("histodEdxLay3","",(Int_t)ntsdd->GetEntries(),0.,ntsdd->GetEntries());
  TH1F* histodEdxLay4=new TH1F("histodEdxLay4","",(Int_t)ntsdd->GetEntries(),0.,ntsdd->GetEntries());
  TH1F* histoTrackClu1=new TH1F("histoTrackClu1","",(Int_t)ntsdd->GetEntries(),0.,ntsdd->GetEntries());
  TH1F* histoTrackClu2=new TH1F("histoTrackClu2","",(Int_t)ntsdd->GetEntries(),0.,ntsdd->GetEntries());
  TH1F* histoTrackClu3=new TH1F("histoTrackClu3","",(Int_t)ntsdd->GetEntries(),0.,ntsdd->GetEntries());
  TH1F* histoTrackClu4=new TH1F("histoTrackClu4","",(Int_t)ntsdd->GetEntries(),0.,ntsdd->GetEntries());
  TH1F* histoTrackClu5=new TH1F("histoTrackClu5","",(Int_t)ntsdd->GetEntries(),0.,ntsdd->GetEntries());
  TH1F* histoTrackClu6=new TH1F("histoTrackClu6","",(Int_t)ntsdd->GetEntries(),0.,ntsdd->GetEntries());

  TH1F* histoNmodEffBelow95=new TH1F("histoNmodEffBelow95","",(Int_t)ntsdd->GetEntries(),0.,ntsdd->GetEntries());
  TH1F* histoNmodEffBelow80=new TH1F("histoNmodEffBelow80","",(Int_t)ntsdd->GetEntries(),0.,ntsdd->GetEntries());
  TH1F* histoNmodEffBelow60=new TH1F("histoNmodEffBelow60","",(Int_t)ntsdd->GetEntries(),0.,ntsdd->GetEntries());
  TH1F* histoNmodEmpty=new TH1F("histoNmodEmpty","",(Int_t)ntsdd->GetEntries(),0.,ntsdd->GetEntries());
  // Sort entries according to run number
  // Same order is assumed for all the subsequent ntuples
  Int_t nr=ntsdd->GetEntries();
  Int_t *myIndex = new Int_t [nr];
  Int_t *noRuns = new Int_t [nr];
  for(Int_t i=0; i<nr;i++){
    ntsdd->GetEvent(i);
    noRuns[i]=nrun;
  }
  TMath::Sort(nr,noRuns,myIndex,kFALSE);
  for(Int_t i=0; i<ntsdd->GetEntries();i++){
    ntsdd->GetEvent(myIndex[i]);
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
    histominTime->SetBinContent(i+1,minDrTime);
    histominTime->SetBinError(i+1,errminDrTime);
    histominTime->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    histomeanTime->SetBinContent(i+1,meanDrTime);
    histomeanTime->SetBinError(i+1,errmeanDrTime);
    histomeanTime->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    histotrp3->SetBinContent(i+1,meanTrPts3);
    histotrp3->SetBinError(i+1,errmeanTrPts3);
    histotrp3->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    histotrp4->SetBinContent(i+1,meanTrPts4);
    histotrp4->SetBinError(i+1,errmeanTrPts3);
    histotrp4->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    histofracExtra->SetBinContent(i+1,fracExtra);
    histofracExtra->SetBinError(i+1,errfracExtra);
    histofracExtra->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    histodEdxTB0->SetBinContent(i+1,meandEdxTB0);
    histodEdxTB0->SetBinError(i+1,errmeandEdxTB0);
    histodEdxTB0->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    histodEdxTB5->SetBinContent(i+1,meandEdxTB5);
    histodEdxTB5->SetBinError(i+1,errmeandEdxTB5);
    histodEdxTB5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    histodEdxLay3->SetBinContent(i+1,meandEdxLay3);
    histodEdxLay3->SetBinError(i+1,errmeandEdxLay3);
    histodEdxLay3->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    histodEdxLay4->SetBinContent(i+1,meandEdxLay4);
    histodEdxLay4->SetBinError(i+1,errmeandEdxLay4);
    histodEdxLay4->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));

    histoNmodEffBelow95->SetBinContent(i+1,nMod95);
    histoNmodEffBelow95->SetBinError(i+1,0.0000001);
    histoNmodEffBelow95->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    histoNmodEffBelow80->SetBinContent(i+1,nMod80);
    histoNmodEffBelow80->SetBinError(i+1,0.0000001);
    histoNmodEffBelow80->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    histoNmodEffBelow60->SetBinContent(i+1,nMod60);
    histoNmodEffBelow60->SetBinError(i+1,0.0000001);
    histoNmodEffBelow60->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    histoNmodEmpty->SetBinContent(i+1,nModEmpty);
    histoNmodEmpty->SetBinError(i+1,0.000001);
    histoNmodEmpty->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
  }


  //---SSD

  TNtuple* ntssd=(TNtuple*)fil->Get("ntssd");

  Float_t  nrunSSD, meandEdxLay5,errmeandEdxLay5,meandEdxLay6,errmeandEdxLay6;
  Float_t ChargeRatioL5,errChargeratioL5,ChargeRatioL6, errChargeratioL6, moduleOff;

  ntssd->SetBranchAddress("nrun",&nrunSSD);
  ntssd->SetBranchAddress("meandEdxLay5",&meandEdxLay5);
  ntssd->SetBranchAddress("errmeandEdxLay5",&errmeandEdxLay5);
  ntssd->SetBranchAddress("meandEdxLay6",&meandEdxLay6);
  ntssd->SetBranchAddress("errmeandEdxLay6",&errmeandEdxLay6);
  ntssd->SetBranchAddress("ChargeRatioL5",&ChargeRatioL5);
  ntssd->SetBranchAddress("errChargeratioL5",&errChargeratioL5);
  ntssd->SetBranchAddress("ChargeRatioL6",&ChargeRatioL6);
  ntssd->SetBranchAddress("errChargeratioL6",&errChargeratioL6);
  ntssd->SetBranchAddress("moduleOff",&moduleOff);

  TH1F* histodEdxLay5=new TH1F("histodEdxLay5","",(Int_t)ntssd->GetEntries(),0.,ntssd->GetEntries());
  TH1F* histodEdxLay6=new TH1F("histodEdxLay6","",(Int_t)ntssd->GetEntries(),0.,ntssd->GetEntries());

  TH1F* histoChargeRatioLay5=new TH1F("histoChargeRatioLay5","",(Int_t)ntssd->GetEntries(),0.,ntssd->GetEntries());
  TH1F* histoChargeRatioLay6=new TH1F("histoChargeRatioLay6","",(Int_t)ntssd->GetEntries(),0.,ntssd->GetEntries());
 
  TH1F* histoEmpty=new TH1F("histoEmpty","",(Int_t)ntssd->GetEntries(),0.,ntssd->GetEntries());
  
  for(Int_t i=0; i<ntssd->GetEntries();i++){

    ntssd->GetEvent(myIndex[i]);

    histodEdxLay5->SetBinContent(i+1,meandEdxLay5);
    histodEdxLay5->SetBinError(i+1,errmeandEdxLay5);
    histodEdxLay5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunSSD));

    histodEdxLay6->SetBinContent(i+1,meandEdxLay6);
    histodEdxLay6->SetBinError(i+1,errmeandEdxLay6);
    histodEdxLay6->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunSSD));

    histoChargeRatioLay5->SetBinContent(i+1,ChargeRatioL5);
    histoChargeRatioLay5->SetBinError(i+1,errChargeratioL5);
    histoChargeRatioLay5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunSSD));
    
    histoChargeRatioLay6->SetBinContent(i+1,ChargeRatioL6);
    histoChargeRatioLay6->SetBinError(i+1,errChargeratioL6);
    histoChargeRatioLay6->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunSSD));

    histoEmpty->SetBinContent(i+1,moduleOff);
    histoEmpty->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunSSD));

  }

  //Matching

  TNtuple* ntmatching=(TNtuple*)fil->Get("ntmatching");

  Float_t nrunMatch;
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

  ntmatching->SetBranchAddress("nrun",&nrunMatch);
  //  ntmatching->SetBranchAddress("nrunMatch",&nrunMatch);
  ntmatching->SetBranchAddress("FracSPD1",&FracSPD1);
  ntmatching->SetBranchAddress("errFracSPD1",&errFracSPD1);
  ntmatching->SetBranchAddress("FracSPD2",&FracSPD2);
  ntmatching->SetBranchAddress("errFracSPD2",&errFracSPD2);
  ntmatching->SetBranchAddress("Eff6Pt02",&Eff6Pt02);
  ntmatching->SetBranchAddress("errEff6Pt02",&errEff6Pt02);
  ntmatching->SetBranchAddress("Eff6Pt1",&Eff6Pt1);
  ntmatching->SetBranchAddress("errEff6Pt1",&errEff6Pt1);
  ntmatching->SetBranchAddress("Eff6Pt10",&Eff6Pt10);
  ntmatching->SetBranchAddress("errEff6Pt10",&errEff6Pt10);
  ntmatching->SetBranchAddress("Eff5Pt02",&Eff5Pt02);
  ntmatching->SetBranchAddress("errEff5Pt02",&errEff5Pt02);
  ntmatching->SetBranchAddress("Eff5Pt1",&Eff5Pt1);
  ntmatching->SetBranchAddress("errEff5Pt1",&errEff5Pt1);
  ntmatching->SetBranchAddress("Eff5Pt10",&Eff5Pt10);
  ntmatching->SetBranchAddress("errEff5Pt10",&errEff5Pt10);
  ntmatching->SetBranchAddress("Eff4Pt02",&Eff4Pt02);
  ntmatching->SetBranchAddress("errEff4Pt02",&errEff4Pt02);
  ntmatching->SetBranchAddress("Eff4Pt1",&Eff4Pt1);
  ntmatching->SetBranchAddress("errEff4Pt1",&errEff4Pt1);
  ntmatching->SetBranchAddress("Eff4Pt10",&Eff4Pt10);
  ntmatching->SetBranchAddress("errEff4Pt10",&errEff4Pt10);
  ntmatching->SetBranchAddress("Eff3Pt02",&Eff3Pt02);
  ntmatching->SetBranchAddress("errEff3Pt02",&errEff3Pt02);
  ntmatching->SetBranchAddress("Eff3Pt1",&Eff3Pt1);
  ntmatching->SetBranchAddress("errEff3Pt1",&errEff3Pt1);
  ntmatching->SetBranchAddress("Eff3Pt10",&Eff3Pt10);
  ntmatching->SetBranchAddress("errEff3Pt10",&errEff3Pt10);
  ntmatching->SetBranchAddress("Eff2Pt02",&Eff2Pt02);
  ntmatching->SetBranchAddress("errEff2Pt02",&errEff2Pt02);
  ntmatching->SetBranchAddress("Eff2Pt1",&Eff2Pt1);
  ntmatching->SetBranchAddress("errEff2Pt1",&errEff2Pt1);
  ntmatching->SetBranchAddress("Eff2Pt10",&Eff2Pt10);
  ntmatching->SetBranchAddress("errEff2Pt10",&errEff2Pt10);
  ntmatching->SetBranchAddress("EffSPDPt02",&EffSPDPt02);
  ntmatching->SetBranchAddress("errEffSPDPt02",&errEffSPDPt02);  
  ntmatching->SetBranchAddress("EffSPDPt1",&EffSPDPt1);
  ntmatching->SetBranchAddress("errEffSPDPt1",&errEffSPDPt1);  
  ntmatching->SetBranchAddress("EffSPDPt10",&EffSPDPt10);
  ntmatching->SetBranchAddress("errEffSPDPt10",&errEffSPDPt10);  
  ntmatching->SetBranchAddress("EffoneSPDPt02",&EffoneSPDPt02);
  ntmatching->SetBranchAddress("errEffoneSPDPt02",&errEffoneSPDPt02);
  ntmatching->SetBranchAddress("EffoneSPDPt1",&EffoneSPDPt1);
  ntmatching->SetBranchAddress("errEffoneSPDPt1",&errEffoneSPDPt1);
  ntmatching->SetBranchAddress("EffoneSPDPt10",&EffoneSPDPt10);
  ntmatching->SetBranchAddress("errEffoneSPDPt10",&errEffoneSPDPt10);
  ntmatching->SetBranchAddress("EffTOTPt02",&EffTOTPt02);
  ntmatching->SetBranchAddress("errEffTOTPt02",&errEffTOTPt02);
  ntmatching->SetBranchAddress("EffTOTPt1",&EffTOTPt1);
  ntmatching->SetBranchAddress("errEffTOTPt1",&errEffTOTPt1);
  ntmatching->SetBranchAddress("EffTOTPt10",&EffTOTPt10);
  ntmatching->SetBranchAddress("errEffTOTPt10",&errEffTOTPt10);


  TH1F *hFracSPD1 = new TH1F("hFracSPD1","SPD inner; run number; Fraction of HSs",(Int_t)ntmatching->GetEntries(),0.,ntmatching->GetEntries());
  hFracSPD1->SetLineColor(3);
  hFracSPD1->SetMarkerColor(3);
  hFracSPD1->SetMarkerStyle(20);

  TH1F *hFracSPD2 = new TH1F("hFracSPD2","SPD outer; run number; Fraction of HSs",(Int_t)ntmatching->GetEntries(),0.,ntmatching->GetEntries());
  hFracSPD2->SetLineColor(8);
  hFracSPD2->SetMarkerColor(8);
  hFracSPD2->SetMarkerStyle(20);

  TH1F *hEffSPDPt02 = new TH1F("hEffSPDPt02","Efficiency - P_{T} = 0.2; run number; TPC+ITS / TPC",(Int_t)ntmatching->GetEntries(),0.,ntmatching->GetEntries());
  hEffSPDPt02->SetLineWidth(2);
  hEffSPDPt02->SetLineColor(kAzure+1);
  hEffSPDPt02->SetMarkerColor(kAzure+1);
  hEffSPDPt02->SetMarkerStyle(20);

  TH1F *hEffSPDPt1 = new TH1F("hEffSPDPt1","Efficiency; run number; TPC+ITS / TPC",(Int_t)ntmatching->GetEntries(),0.,ntmatching->GetEntries());
  hEffSPDPt1->SetLineWidth(2);
  hEffSPDPt1->SetLineColor(kAzure+1);
  hEffSPDPt1->SetMarkerColor(kAzure+1);
  hEffSPDPt1->SetMarkerStyle(20);

  TH1F *hEffSPDPt10 = new TH1F("hEffSPDPt10","Efficiency; run number; TPC+ITS / TPC",(Int_t)ntmatching->GetEntries(),0.,ntmatching->GetEntries());
  hEffSPDPt10->SetLineWidth(2);
  hEffSPDPt10->SetLineColor(kAzure+1);
  hEffSPDPt10->SetMarkerColor(kAzure+1);
  hEffSPDPt10->SetMarkerStyle(20);

  TH1F *hEffoneSPDPt02 = new TH1F("hEffoneSPDPt02","Efficiency; run number; TPC+ITS / TPC",(Int_t)ntmatching->GetEntries(),0.,ntmatching->GetEntries());
  hEffoneSPDPt02->SetLineWidth(2);
  hEffoneSPDPt02->SetLineColor(kGray);
  hEffoneSPDPt02->SetMarkerColor(kGray);
  hEffoneSPDPt02->SetMarkerStyle(20);
  TH1F *hEffoneSPDPt1 = new TH1F("hEffoneSPDPt1","Efficiency; run number; TPC+ITS / TPC",(Int_t)ntmatching->GetEntries(),0.,ntmatching->GetEntries());
  hEffoneSPDPt1->SetLineWidth(2);
  hEffoneSPDPt1->SetLineColor(kGray);
  hEffoneSPDPt1->SetMarkerColor(kGray);
  hEffoneSPDPt1->SetMarkerStyle(20);
  TH1F *hEffoneSPDPt10 = new TH1F("hEffoneSPDPt10","Efficiency; run number; TPC+ITS / TPC",(Int_t)ntmatching->GetEntries(),0.,ntmatching->GetEntries());
  hEffoneSPDPt10->SetLineWidth(2);
  hEffoneSPDPt10->SetLineColor(kGray);
  hEffoneSPDPt10->SetMarkerColor(kGray);
  hEffoneSPDPt10->SetMarkerStyle(20);

  TH1F *hEff2Pt02 = new TH1F("hEff2Pt02","Efficiency; run number; TPC+ITS / TPC",(Int_t)ntmatching->GetEntries(),0.,ntmatching->GetEntries());
  hEff2Pt02->SetLineWidth(2);
  hEff2Pt02->SetLineColor(kViolet);
  hEff2Pt02->SetMarkerColor(kViolet);
  hEff2Pt02->SetMarkerStyle(20);
  TH1F *hEff2Pt1 = new TH1F("hEff2Pt1","Efficiency; run number; TPC+ITS / TPC",(Int_t)ntmatching->GetEntries(),0.,ntmatching->GetEntries());
  hEff2Pt1->SetLineWidth(2);
  hEff2Pt1->SetLineColor(kViolet);
  hEff2Pt1->SetMarkerColor(kViolet);
  hEff2Pt1->SetMarkerStyle(20);
  TH1F *hEff2Pt10 = new TH1F("hEff2Pt10","Efficiency; run number; TPC+ITS / TPC",(Int_t)ntmatching->GetEntries(),0.,ntmatching->GetEntries());
  hEff2Pt10->SetLineWidth(2);
  hEff2Pt10->SetLineColor(kViolet);
  hEff2Pt10->SetMarkerColor(kViolet);
  hEff2Pt10->SetMarkerStyle(20);

  TH1F *hEff3Pt02 = new TH1F("hEff3Pt02","Efficiency; run number; TPC+ITS / TPC",(Int_t)ntmatching->GetEntries(),0.,ntmatching->GetEntries());
  hEff3Pt02->SetLineWidth(2);
  hEff3Pt02->SetLineColor(6);
  hEff3Pt02->SetMarkerColor(6);
  hEff3Pt02->SetMarkerStyle(20);
  TH1F *hEff3Pt1 = new TH1F("hEff3Pt1","Efficiency; run number; TPC+ITS / TPC",(Int_t)ntmatching->GetEntries(),0.,ntmatching->GetEntries());
  hEff3Pt1->SetLineWidth(2);
  hEff3Pt1->SetLineColor(6);
  hEff3Pt1->SetMarkerColor(6);
  hEff3Pt1->SetMarkerStyle(20);
  TH1F *hEff3Pt10 = new TH1F("hEff3Pt10","Efficiency; run number; TPC+ITS / TPC",(Int_t)ntmatching->GetEntries(),0.,ntmatching->GetEntries());
  hEff3Pt10->SetLineWidth(2);
  hEff3Pt10->SetLineColor(6);
  hEff3Pt10->SetMarkerColor(6);
  hEff3Pt10->SetMarkerStyle(20);

  TH1F *hEff4Pt02 = new TH1F("hEff4Pt02","Efficiency; run number; TPC+ITS / TPC",(Int_t)ntmatching->GetEntries(),0.,ntmatching->GetEntries());
  hEff4Pt02->SetLineWidth(2);
  hEff4Pt02->SetLineColor(4);
  hEff4Pt02->SetMarkerColor(4);
  hEff4Pt02->SetMarkerStyle(20);
  TH1F *hEff4Pt1 = new TH1F("hEff4Pt1","Efficiency; run number; TPC+ITS / TPC",(Int_t)ntmatching->GetEntries(),0.,ntmatching->GetEntries());
  hEff4Pt1->SetLineWidth(2);
  hEff4Pt1->SetLineColor(4);
  hEff4Pt1->SetMarkerColor(4);
  hEff4Pt1->SetMarkerStyle(20);
  TH1F *hEff4Pt10 = new TH1F("hEff4Pt10","Efficiency; run number; TPC+ITS / TPC",(Int_t)ntmatching->GetEntries(),0.,ntmatching->GetEntries());
  hEff4Pt10->SetLineWidth(2);
  hEff4Pt10->SetLineColor(4);
  hEff4Pt10->SetMarkerColor(4);
  hEff4Pt10->SetMarkerStyle(20);

  TH1F *hEff5Pt02 = new TH1F("hEff5Pt02","Efficiency; run number; TPC+ITS / TPC",(Int_t)ntmatching->GetEntries(),0.,ntmatching->GetEntries());
  hEff5Pt02->SetLineWidth(2);
  hEff5Pt02->SetLineColor(3);
  hEff5Pt02->SetMarkerColor(3);
  hEff5Pt02->SetMarkerStyle(20);
  TH1F *hEff5Pt1 = new TH1F("hEff5Pt1","Efficiency; run number; TPC+ITS / TPC",(Int_t)ntmatching->GetEntries(),0.,ntmatching->GetEntries());
  hEff5Pt1->SetLineWidth(2);
  hEff5Pt1->SetLineColor(3);
  hEff5Pt1->SetMarkerColor(3);
  hEff5Pt1->SetMarkerStyle(20);
  TH1F *hEff5Pt10 = new TH1F("hEff5Pt10","Efficiency; run number; TPC+ITS / TPC",(Int_t)ntmatching->GetEntries(),0.,ntmatching->GetEntries());
  hEff5Pt10->SetLineWidth(3);
  hEff5Pt10->SetLineColor(3);
  hEff5Pt10->SetMarkerColor(3);
  hEff5Pt10->SetMarkerStyle(20);

  TH1F *hEff6Pt02 = new TH1F("hEff6Pt02","Efficiency; run number; TPC+ITS / TPC",(Int_t)ntmatching->GetEntries(),0.,ntmatching->GetEntries());
  hEff6Pt02->SetLineWidth(2);
  hEff6Pt02->SetLineColor(2);
  hEff6Pt02->SetMarkerColor(2);
  hEff6Pt02->SetMarkerStyle(20);
  TH1F *hEff6Pt1 = new TH1F("hEff6Pt1","Efficiency; run number; TPC+ITS / TPC",(Int_t)ntmatching->GetEntries(),0.,ntmatching->GetEntries());
  hEff6Pt1->SetLineWidth(2);
  hEff6Pt1->SetLineColor(2);
  hEff6Pt1->SetMarkerColor(2);
  hEff6Pt1->SetMarkerStyle(20);
  TH1F *hEff6Pt10 = new TH1F("hEff6Pt10","Efficiency; run number; TPC+ITS / TPC",(Int_t)ntmatching->GetEntries(),0.,ntmatching->GetEntries());
  hEff6Pt10->SetLineWidth(2);
  hEff6Pt10->SetLineColor(2);
  hEff6Pt10->SetMarkerColor(2);
  hEff6Pt10->SetMarkerStyle(20);


  TH1F *hEffTOTPt02 = new TH1F("hEffTOTPt02","Efficiency; run number; TPC+ITS / TPC",(Int_t)ntmatching->GetEntries(),0.,ntmatching->GetEntries());
  hEffTOTPt02->SetLineWidth(2);
  hEffTOTPt02->SetLineColor(kBlue+2);
  hEffTOTPt02->SetMarkerColor(kBlue+2);
  hEffTOTPt02->SetMarkerStyle(20);
  TH1F *hEffTOTPt1 = new TH1F("hEffTOTPt1","Efficiency; run number; TPC+ITS / TPC",(Int_t)ntmatching->GetEntries(),0.,ntmatching->GetEntries());
  hEffTOTPt1->SetLineWidth(2);
  hEffTOTPt1->SetLineColor(kBlue+2);
  hEffTOTPt1->SetMarkerColor(kBlue+2);
  hEffTOTPt1->SetMarkerStyle(20);
  TH1F *hEffTOTPt10 = new TH1F("hEffTOTPt10","Efficiency; run number; TPC+ITS / TPC",(Int_t)ntmatching->GetEntries(),0.,ntmatching->GetEntries());
  hEffTOTPt10->SetLineWidth(2);
  hEffTOTPt10->SetLineColor(kBlue+2);
  hEffTOTPt10->SetMarkerColor(kBlue+2);
  hEffTOTPt10->SetMarkerStyle(20);

  Int_t nEntriesMatch=ntmatching->GetEntries();
  
  for(Int_t i=0;i<nEntriesMatch;i++){
   
    ntmatching->GetEvent(myIndex[i]);
    //    Int_t bin=nrunMatch;

    // fill histos
        
    hFracSPD1->SetBinContent(i,FracSPD1);
    hFracSPD1->SetBinError(i,.01);
    hFracSPD1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));

    //    cout<<FracSPD1<<endl;

    hFracSPD2->SetBinContent(i+1,FracSPD2);
    hFracSPD2->SetBinError(i+1,.01);
    hFracSPD2->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));
    
    //    cout<<FracSPD2<<endl;

    //-------------------------
      
    hEff6Pt02->SetBinContent(i+1,Eff6Pt02);
    hEff6Pt02->SetBinError(i+1,errEff6Pt02);
    hEff6Pt02->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));
      
    hEff6Pt1->SetBinContent(i+1,Eff6Pt1);
    hEff6Pt1->SetBinError(i+1,errEff6Pt1);
    hEff6Pt1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));
      
    hEff6Pt10->SetBinContent(i+1,Eff6Pt10);
    hEff6Pt10->SetBinError(i+1,errEff6Pt10);
    hEff6Pt10->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));

    hEff5Pt02->SetBinContent(i+1,Eff5Pt02);
    hEff5Pt02->SetBinError(i+1,errEff5Pt02);
    hEff5Pt02->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));
	    
    hEff5Pt1->SetBinContent(i+1,Eff5Pt1);
    hEff5Pt1->SetBinError(i+1,errEff5Pt1);
    hEff5Pt1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));
	    
    hEff5Pt10->SetBinContent(i+1,Eff5Pt10);
    hEff5Pt10->SetBinError(i+1,errEff5Pt10);
    hEff5Pt10->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));
      
    hEff4Pt02->SetBinContent(i+1,Eff4Pt02);
    hEff4Pt02->SetBinError(i+1,errEff4Pt1);
    hEff4Pt02->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));

    hEff4Pt1->SetBinContent(i+1,Eff4Pt1);
    hEff4Pt1->SetBinError(i+1,errEff4Pt1);
    hEff4Pt1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));

    hEff4Pt10->SetBinContent(i+1,Eff4Pt10);
    hEff4Pt10->SetBinError(i+1,errEff4Pt10);
    hEff4Pt10->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));
      
    hEff3Pt02->SetBinContent(i+1,Eff3Pt02);
    hEff3Pt02->SetBinError(i+1,errEff3Pt02);
    hEff3Pt02->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));
	    
    hEff3Pt1->SetBinContent(i+1,Eff3Pt1);
    hEff3Pt1->SetBinError(i+1,errEff3Pt1);
    hEff3Pt1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));

    hEff3Pt10->SetBinContent(i+1,Eff3Pt10);
    hEff3Pt10->SetBinError(i+1,errEff3Pt10);
    hEff3Pt10->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));      
      
    hEff2Pt02->SetBinContent(i+1,Eff2Pt02);
    hEff2Pt02->SetBinError(i+1,errEff2Pt02);
    hEff2Pt02->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));  

    hEff2Pt1->SetBinContent(i+1,Eff2Pt1);
    hEff2Pt1->SetBinError(i+1,errEff2Pt1);
    hEff2Pt1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));  

    hEff2Pt10->SetBinContent(i+1,Eff2Pt10);
    hEff2Pt10->SetBinError(i+1,errEff2Pt10);
    hEff2Pt10->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));  


    hEffSPDPt02->SetBinContent(i+1,EffSPDPt02);
    hEffSPDPt02->SetBinError(i+1,errEffSPDPt02);
    hEffSPDPt02->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));  
      
    hEffSPDPt1->SetBinContent(i+1,EffSPDPt1);
    hEffSPDPt1->SetBinError(i+1,errEffSPDPt1);
    hEffSPDPt1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch)); 
      
    hEffSPDPt10->SetBinContent(i+1,EffSPDPt10);
    hEffSPDPt10->SetBinError(i+1,errEffSPDPt10);
    hEffSPDPt10->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch)); 
      
    hEffoneSPDPt02->SetBinContent(i+1,EffoneSPDPt02);
    hEffoneSPDPt02->SetBinError(i+1,errEffoneSPDPt02);
    hEffoneSPDPt02->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));

    hEffoneSPDPt1->SetBinContent(i+1,EffoneSPDPt1);
    hEffoneSPDPt1->SetBinError(i+1,errEffoneSPDPt1);
    hEffoneSPDPt1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));
      
    hEffoneSPDPt10->SetBinContent(i+1,EffoneSPDPt10);
    hEffoneSPDPt10->SetBinError(i+1,errEffoneSPDPt10);
    hEffoneSPDPt10->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));

    hEffTOTPt02->SetBinContent(i+1,EffTOTPt02);
    hEffTOTPt02->SetBinError(i+1,errEffTOTPt02);
    hEffTOTPt02->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));
      
    hEffTOTPt1->SetBinContent(i+1,EffTOTPt1);
    hEffTOTPt1->SetBinError(i+1,errEffTOTPt1);
    hEffTOTPt1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));
      
    hEffTOTPt10->SetBinContent(i+1,EffTOTPt10);
    hEffTOTPt10->SetBinError(i+1,errEffTOTPt10);
    hEffTOTPt10->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));  

  }
  
  //-----------------------------------

  //---vertex

  TNtuple* ntvertex=(TNtuple*)fil->Get("ntvertex");

  Float_t nrunVertex,Vx,errVx,sigmaVx,errsigmaVx,Vy,errVy,sigmaVy,errsigmaVy,Vz,errVz,sigmaVz,errsigmaVz;

  ntvertex->SetBranchAddress("nrun",&nrunVertex);
  ntvertex->SetBranchAddress("Vx",&Vx);
  ntvertex->SetBranchAddress("errVx",&errVx);
  ntvertex->SetBranchAddress("sigmaVx",&sigmaVx);
  ntvertex->SetBranchAddress("errsigmaVx",&errsigmaVx);
  ntvertex->SetBranchAddress("Vy",&Vy);
  ntvertex->SetBranchAddress("errVy",&errVy);
  ntvertex->SetBranchAddress("sigmaVy",&sigmaVy);
  ntvertex->SetBranchAddress("errsigmaVy",&errsigmaVy);
  ntvertex->SetBranchAddress("Vz",&Vz);
  ntvertex->SetBranchAddress("errVz",&errVz);
  ntvertex->SetBranchAddress("sigmaVz",&sigmaVz);
  ntvertex->SetBranchAddress("errsigmaVz",&errsigmaVz);


  TH1F *hVx = new TH1F("hVx","Track Vertex Vx Distribution",(Int_t)ntvertex->GetEntries(),0.,ntvertex->GetEntries());
  hVx->SetLineWidth(2);
  hVx->SetLineColor(kBlue+2);
  hVx->SetMarkerColor(kBlue+2);
  hVx->SetMarkerStyle(20);

 TH1F *hVy = new TH1F("hVy","Track Vertex Vy Distribution",(Int_t)ntvertex->GetEntries(),0.,ntvertex->GetEntries());
  hVy->SetLineWidth(2);
  hVy->SetLineColor(kBlue+2);
  hVy->SetMarkerColor(kBlue+2);
  hVy->SetMarkerStyle(20);

 TH1F *hVz = new TH1F("hVz","Track Vertex Vz Distribution",(Int_t)ntvertex->GetEntries(),0.,ntvertex->GetEntries());
  hVz->SetLineWidth(2);
  hVz->SetLineColor(kBlue+2);
  hVz->SetMarkerColor(kBlue+2);
  hVz->SetMarkerStyle(20);

  TH1F *hSigmaVx = new TH1F("hSigmaVx","Track Vertex SigmaVx Distribution",(Int_t)ntvertex->GetEntries(),0.,ntvertex->GetEntries());
  hSigmaVx->SetLineWidth(2);
  hSigmaVx->SetLineColor(kBlue+2);
  hSigmaVx->SetMarkerColor(kBlue+2);
  hSigmaVx->SetMarkerStyle(20);

 TH1F *hSigmaVy = new TH1F("hSigmaVy","Track Vertex SigmaVy Distribution",(Int_t)ntvertex->GetEntries(),0.,ntvertex->GetEntries());
  hSigmaVy->SetLineWidth(2);
  hSigmaVy->SetLineColor(kBlue+2);
  hSigmaVy->SetMarkerColor(kBlue+2);
  hSigmaVy->SetMarkerStyle(20);

 TH1F *hSigmaVz = new TH1F("hSigmaVz","Track Vertex SigmaVz Distribution",(Int_t)ntvertex->GetEntries(),0.,ntvertex->GetEntries());
  hSigmaVz->SetLineWidth(2);
  hSigmaVz->SetLineColor(kBlue+2);
  hSigmaVz->SetMarkerColor(kBlue+2);
  hSigmaVz->SetMarkerStyle(20);


  Int_t nEntriesVertex=ntvertex->GetEntries();
  
  for(Int_t i=0;i<nEntriesVertex;i++){

    ntvertex->GetEvent(myIndex[i]);
    
    //    cout<<Vx<<endl;

    hVx->SetBinContent(i+1,Vx);
    hVx->SetBinError(i+1,errVx);
    hVx->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunVertex));
  
    hVy->SetBinContent(i+1,Vy);
    hVy->SetBinError(i+1,errVy);
    hVy->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunVertex));
    
    hVz->SetBinContent(i+1,Vz);
    hVz->SetBinError(i+1,errVz);
    hVz->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunVertex));
  
    
    hSigmaVx->SetBinContent(i+1,sigmaVx);
    hSigmaVx->SetBinError(i+1,errsigmaVx);
    hSigmaVx->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunVertex));

    hSigmaVy->SetBinContent(i+1,sigmaVy);
    hSigmaVy->SetBinError(i+1,errsigmaVy);
    hSigmaVy->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunVertex));

    hSigmaVz->SetBinContent(i+1,sigmaVz);
    hSigmaVz->SetBinError(i+1,errsigmaVz);
    hSigmaVz->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunVertex));
    
  }
 //--------  Draw Vertex histograms ---------
  TCanvas *cVertexDisto=new TCanvas("cVertexDisto","cVertexDisto",1200,800);
  cVertexDisto->Divide(3,2);
  cVertexDisto->cd(1);
  hVx->SetMinimum(0.03);
  hVx->SetMaximum(0.08);
  hVx->GetYaxis()->SetTitle("Vertex X coordinate");
  hVx->Draw();
  cVertexDisto->cd(2);
  hVy->SetMinimum(0.25);
  hVy->SetMaximum(0.30);
  hVy->GetYaxis()->SetTitle("Vertex Y coordinate");
  hVy->Draw();
  cVertexDisto->cd(3);
  hVz->SetMinimum(-1.);
  hVz->SetMaximum(1.);
  hVz->GetYaxis()->SetTitle("Vertex Z coordinate");
  hVz->Draw();
  cVertexDisto->cd(4);
  hSigmaVx->SetMinimum(0.);
  hSigmaVx->SetMaximum(0.01);
  hSigmaVx->GetYaxis()->SetTitle("sigma on x coordinate");
  hSigmaVx->Draw();
  cVertexDisto->cd(5);
  hSigmaVy->SetMinimum(0.);
  hSigmaVy->SetMaximum(0.01);
  hSigmaVy->GetYaxis()->SetTitle("sigma on y coordinate");
  hSigmaVy->Draw();
  cVertexDisto->cd(6);
  hSigmaVz->SetMinimum(6.);
  hSigmaVz->SetMaximum(10.);
  hSigmaVz->GetYaxis()->SetTitle("sigma on z coordinate");
  hSigmaVz->Draw();
  cVertexDisto->SaveAs("Vertex_trend.pdf");
    pdfFileNames+=" Vertex_trend.pdf";
 //-----------------------------------

  gStyle->SetOptStat(0);

  TCanvas* c5=new TCanvas("c5","Track With points in ITS");
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
  histoTrackClu6->SetTitle("Fraction of tracks with points in ITS");
  histoTrackClu6->Draw("same");
  histoTrackClu3->GetYaxis()->SetTitle("Fraction of Tracks with Points in ITS Layers");
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
    pdfFileNames+=" TrackPoints_trend.pdf";
  c5->Update();



  TCanvas* c2=new TCanvas("c2","SDD DriftTime & Charge",1200,800);
  c2->Divide(2,2);
  c2->cd(1);
  histominTime->Draw();
  histominTime->SetMinimum(450);
  histominTime->SetMaximum(550);
  histominTime->GetYaxis()->SetTitle("Minimum Drift Time (ns)");
  TLatex* td1=new TLatex(0.2,0.85,"SDD minimum drift time (ref=505ns)");
  td1->SetNDC();
  td1->SetTextColor(1);
  td1->Draw();
  c2->cd(2);
  histomeanTime->Draw();
  histomeanTime->SetMinimum(3200);
  histomeanTime->SetMaximum(3300);
  histomeanTime->GetYaxis()->SetTitle("Average Drift Time (ns)");
  TLatex* td2=new TLatex(0.2,0.85,"SDD average drift time");
  td2->SetNDC();
  // td2>SetTextColor(1);
  td2->Draw();

  // TCanvas* c4=new TCanvas("c4","Charge");
  c2->cd(3);
  histodEdxTB0->SetLineColor(1);
  histodEdxTB0->SetMarkerStyle(20);
  histodEdxTB0->Draw();
  histodEdxTB0->SetMinimum(90.);
  histodEdxTB0->SetMaximum(120.);
  histodEdxTB5->SetLineColor(4);
  histodEdxTB5->SetMarkerColor(4);
  histodEdxTB5->SetMarkerStyle(23);
  // histodEdxTB5->SetMinimum(90);
  // histodEdxTB5->SetMaximum(120);
  histodEdxTB5->Draw("same");
  histodEdxTB0->GetYaxis()->SetTitle("<dE/dx> (keV/300 #mum)");
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
  // TCanvas* c4b=new TCanvas("c4b","Charge per Layer");
  c2->cd(4);
  histodEdxLay3->SetLineColor(1);
  histodEdxLay3->SetMarkerStyle(20);
  histodEdxLay3->Draw();
  histodEdxLay3->SetMinimum(90.);
  histodEdxLay3->SetMaximum(120.);
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
    
  histodEdxLay3->GetYaxis()->SetTitle("<dE/dx> (keV/300 #mum)");
  
  TLegend* leg2b=new TLegend(0.6,0.15,0.88,0.35);
  ent=leg2b->AddEntry(histodEdxLay3,"Layer 3","PL");
  ent=leg2b->AddEntry(histodEdxLay4,"Layer 4","PL");
  ent=leg2b->AddEntry(histodEdxLay5,"Layer 5","PL");
  ent=leg2b->AddEntry(histodEdxLay6,"Layer 6","PL");
  ent->SetTextColor(histodEdxLay4->GetMarkerColor());
  leg2b->SetFillStyle(0);
  leg2b->Draw();
  // c4b->Update();
  TLatex* tc2=new TLatex(0.2,0.85,"SDD and SSD charge in different layers");
  tc2->SetNDC();
  tc2->SetTextColor(1);
  tc2->Draw();
  c2->SaveAs("SDD_SSD_drift_charge_trend.pdf");  
  pdfFileNames+=" SDD_SSD_drift_charge_trend.pdf";
  c2->Update();

  TCanvas *c7=new TCanvas("c7","Charge ratio");
  c7->cd();
  histoChargeRatioLay5->SetLineColor(6);
  histoChargeRatioLay5->SetMarkerColor(6);
  histoChargeRatioLay5->SetMarkerStyle(20);
  histoChargeRatioLay5->SetMinimum(-0.01);
  histoChargeRatioLay5->SetMaximum(+0.01);
  histoChargeRatioLay5->Draw();
  histoChargeRatioLay6->SetLineColor(7);
  histoChargeRatioLay6->SetMarkerColor(7);
  histoChargeRatioLay6->SetMarkerStyle(22);
  histoChargeRatioLay6->GetYaxis()->SetTitle("SSD charge ratio");
  histoChargeRatioLay6->Draw("same");
  TLegend* legCR=new TLegend(0.7,0.75,0.88,0.85);
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
  c7->SaveAs("SSD_chargeratio_trend.pdf");
      pdfFileNames+=" SSD_chargeratio_trend.pdf";
  // TCanvas *c8=new TCanvas("c8","Masked modules");
  // c8->cd();
  // histoEmpty->Draw();

  TCanvas *cpt02=new TCanvas("cpt02","TPC-ITS matching efficiency",1200,1000);
  cpt02->Divide(2,2);
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
  hEff6Pt02->GetYaxis()->SetRangeUser(0,1);
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

  // TCanvas *cpt1=new TCanvas("cpt1","cpt1");
  // cpt1->cd(1);
  cpt02->cd(2);

  hEff6Pt1->Draw();
  hEff5Pt1->Draw("same");
  hEff4Pt1->Draw("same");
  hEff3Pt1->Draw("same");
  hEff2Pt1->Draw("same");
  hEffSPDPt1->Draw("same");
  hEffoneSPDPt1->Draw("same");
  hEffTOTPt1->Draw("same");
  hEff6Pt1->GetYaxis()->SetRangeUser(0,1);

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


  // TCanvas *cpt10=new TCanvas("cpt10","cpt10");
  cpt02->cd(3);

  hEff6Pt10->Draw();
  hEff5Pt10->Draw("same");
  hEff4Pt10->Draw("same");
  hEff3Pt10->Draw("same");
  hEff2Pt10->Draw("same");
  hEffSPDPt10->Draw("same");
  hEffoneSPDPt10->Draw("same");
  hEffTOTPt10->Draw("same");
  hEff6Pt10->GetYaxis()->SetRangeUser(0,1);

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

  TLatex* tpc3=new TLatex(0.2,0.75,"TPCITS match eff Pt=10");
  tpc3->SetNDC();
  tpc3->SetTextColor(1);
  tpc3->Draw();
 
  cpt02->cd(4);


  //  TCanvas *cSPD = new TCanvas("cSPD","cSPD",0,0,1000,300);
  // cSPD->SetGridy();
  hFracSPD1->SetMaximum(1.2);
  hFracSPD1->SetMaximum(0);
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
  cpt02->SaveAs("TPCITS_trend.pdf");
      pdfFileNames+=" TPCITS_trend.pdf";
  cpt02->Update();

  //-----------  ITS SA -----------
  PlotITSSA(fil,myIndex);
 // merge the pdf files
  TString command("gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=merged");
  command=command+"ITS_trend.pdf "+pdfFileNames;
  gSystem->Exec(command.Data());
  printf(" Merging the pdf file:  %s \n",command.Data());
  delete [] myIndex;
  delete [] noRuns;
}


void PlotITSSA(TFile *fil,Int_t *myIndex){
  Double_t Lowbin[3]={0.1,0.5,0.9};
  Double_t Upbin[3]={0.2,0.6,1};
  gROOT->SetStyle("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetFillColor(0);
  gStyle->SetTextFont(32);

  TNtuple* nt = (TNtuple*)fil->Get("ntITSsa");
 Int_t nruns=nt->GetEntries();
  printf("Events = %d\n",nruns);
  Float_t ITSA[3];
  Float_t TPIT[3];
  Float_t RAT[3];
  Float_t run;
  nt->SetBranchAddress("run",&run);
  nt->SetBranchAddress("NITSsaPtBin0",&ITSA[0]);
  nt->SetBranchAddress("NITSsaPtBin1",&ITSA[1]);
  nt->SetBranchAddress("NITSsaPtBin2",&ITSA[2]);
  nt->SetBranchAddress("NITSTPCPtBin0",&TPIT[0]);
  nt->SetBranchAddress("NITSTPCPtBin1",&TPIT[1]);
  nt->SetBranchAddress("NITSTPCPtBin2",&TPIT[2]);
  nt->SetBranchAddress("ratioPtBin0",&RAT[0]);
  nt->SetBranchAddress("ratioPtBin1",&RAT[1]);
  nt->SetBranchAddress("ratioPtBin2",&RAT[2]);
  TH1F *h0=new TH1F("h0","h0",nruns,-0.5,nruns-0.5);
  TH1F *h1=new TH1F("h1","h1",nruns,-0.5,nruns-0.5);
  TH1F *h2=new TH1F("h2","h2",nruns,-0.5,nruns-0.5);
  TH1F *h3=new TH1F("h3","h4",nruns,-0.5,nruns-0.5);
  TH1F *h4=new TH1F("h4","h5",nruns,-0.5,nruns-0.5);
  TH1F *h5=new TH1F("h5","h5",nruns,-0.5,nruns-0.5);
  TH1F *h6=new TH1F("h6","h6",nruns,-0.5,nruns-0.5);
  TH1F *h7=new TH1F("h7","h7",nruns,-0.5,nruns-0.5);
  TH1F *h8=new TH1F("h8","h8",nruns,-0.5,nruns-0.5);
  for(Int_t iev=0;iev<nruns;iev++){
    nt->GetEvent(myIndex[iev]);
    //cout<<"Numeri TPCITS "<<TPIT[0]<<" "<<TPIT[1]<<" "<<TPIT[2]<<endl;
    h0->Fill(iev,ITSA[0]);
    h1->Fill(iev,ITSA[1]);
    h2->Fill(iev,ITSA[2]);
    h3->Fill(iev,TPIT[0]);
    h4->Fill(iev,TPIT[1]);
    h5->Fill(iev,TPIT[2]);
    h6->Fill(iev,RAT[0]);
    h7->Fill(iev,RAT[1]);
    h8->Fill(iev,RAT[2]);
    h0->GetXaxis()->SetBinLabel(iev+1,Form("%.0f",run));
    h1->GetXaxis()->SetBinLabel(iev+1,Form("%.0f",run));
    h2->GetXaxis()->SetBinLabel(iev+1,Form("%.0f",run));
    h3->GetXaxis()->SetBinLabel(iev+1,Form("%.0f",run));
    h4->GetXaxis()->SetBinLabel(iev+1,Form("%.0f",run));
    h5->GetXaxis()->SetBinLabel(iev+1,Form("%.0f",run));
    h6->GetXaxis()->SetBinLabel(iev+1,Form("%.0f",run));
    h7->GetXaxis()->SetBinLabel(iev+1,Form("%.0f",run));
    h8->GetXaxis()->SetBinLabel(iev+1,Form("%.0f",run));
    Printf("%f   %f   %f",ITSA[0],ITSA[1],ITSA[2]);
  }
  h0->Print("all");
  h0->GetYaxis()->SetTitle("ITSsa tracks");
  h0->GetXaxis()->SetTitle("run");
  h1->GetYaxis()->SetTitle("ITSsa tracks");
  h1->GetXaxis()->SetTitle("run");
  h2->GetYaxis()->SetTitle("ITSsa tracks");
  h2->GetXaxis()->SetTitle("run");
  h3->GetYaxis()->SetTitle("ITS+TPC tracks");
  h3->GetXaxis()->SetTitle("run");
  h4->GetYaxis()->SetTitle("ITS+TPC tracks");
  h4->GetXaxis()->SetTitle("run");
  h5->GetYaxis()->SetTitle("ITS+TPC tracks");
  h5->GetXaxis()->SetTitle("run");
  h6->GetYaxis()->SetTitle("(TPC+ITS)/ITS");
  h6->GetXaxis()->SetTitle("run");
  h7->GetYaxis()->SetTitle("(TPC+ITS)/ITS");
  h7->GetXaxis()->SetTitle("run");
  h8->GetYaxis()->SetTitle("(TPC+ITS)/ITS");
  h8->GetXaxis()->SetTitle("run");
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
  h8->SetMaximum(20.);

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

  TCanvas *c=new TCanvas();
  c->SetLogy();
  c->SetGridy();
  h0->Draw("p");
  h1->Draw("psame");
  h2->Draw("psame");
c->BuildLegend(0.11,0.15,0.45,0.30);
  TLatex* ti1=new TLatex(0.11,0.40,"ITS standalone tracks (normalized to number of events)");
  ti1->SetNDC();
  ti1->SetTextColor(1);
  ti1->Draw();
  c->SaveAs("ITSsa_trend.pdf");
      pdfFileNames+=" ITSsa_trend.pdf";
  TCanvas *c2=new TCanvas();
  c2->SetLogy();
  c2->SetGridy();
  h3->Draw("p");
  h4->Draw("psame");
  h5->Draw("psame");
 c2->BuildLegend(0.11,0.15,0.45,0.30);
  TLatex* ti2=new TLatex(0.11,0.40,"ITS+TPC tracks (normalized to number of events)");
  ti2->SetNDC();
  ti2->SetTextColor(1);
  ti2->Draw();
   c2->SaveAs("ITSTPC_trend.pdf");
      pdfFileNames+=" ITSTPC_trend.pdf";
 
  TCanvas *c3=new TCanvas();
  c3->SetGridy();
  h8->Draw("p");
  h7->Draw("psame");
  h6->Draw("psame");
   c3->BuildLegend();
  TLatex* ti3=new TLatex(0.15,0.8,"ratio (ITS+TPC)/ITS ");
  ti3->SetNDC();
  ti3->Draw();
   c3->SaveAs("tracks_ratio_trend.pdf");
      pdfFileNames+=" tracks_ratio_trend.pdf";

}

void AliITSQAtrend(TString runListFile,TString ntupleFileName){

  TGrid::Connect("alien://");
  

  
  //-----------SDD
  
  const Int_t nVariables=35;
  TNtuple* ntsdd=new TNtuple("ntsdd","SDD trending","nrun:fracTrackWithClu1:errfracTrackWithClu1:fracTrackWithClu2:errfracTrackWithClu2:fracTrackWithClu3:errfracTrackWithClu3:fracTrackWithClu4:errfracTrackWithClu4:fracTrackWithClu5:errfracTrackWithClu5:fracTrackWithClu6:errfracTrackWithClu6:meanTrPts3:errmeanTrPts3:meanTrPts4:errmeanTrPts4:minDrTime:errminDrTime:meanDrTime:errmeanDrTime:fracExtra:errfracExtra:meandEdxLay3:errmeandEdxLay3:meandEdxLay4:errmeandEdxLay4:meandEdxTB0:errmeandEdxTB0:meandEdxTB5:errmeandEdxTB5:nMod95:nMod80:nMod60:nModEmpty");
  Float_t xnt[nVariables];
  
  //--------------SSD
    
  const Int_t nVariablesSSD=10;
  TNtuple* ntssd=new TNtuple("ntssd","SSD trending","nrun:meandEdxLay5:errmeandEdxLay5:meandEdxLay6:errmeandEdxLay6:ChargeRatioL5:errChargeratioL5:ChargeRatioL6:errChargeratioL6:moduleOff");
  Float_t xntSSD[nVariablesSSD];
  
  //----Matching

  const Int_t nVariablesMatching=60;
  TNtuple* ntmatching=new TNtuple("ntmatching","Matching Efficiency","nrun:FracSPD1:errFracSPD1:FracSPD2:errFracSPD2:Eff6Pt02:errEff6Pt02:Eff6Pt1:errEff6Pt1:Eff6Pt10:errEff6Pt10:Eff5Pt02:errEff5Pt02:Eff5Pt1:errEff5Pt1:Eff5Pt10:errEff5Pt10:Eff4Pt02:errEff4Pt02:Eff4Pt1:errEff4Pt1:Eff4Pt10:errEff4Pt10:Eff3Pt02:errEff3Pt02:Eff3Pt1:errEff3Pt1:Eff3Pt10:errEff3Pt10:Eff2Pt02:errEff2Pt02:Eff2Pt1:errEff2Pt1:Eff2Pt10:errEff2Pt10:EffSPDPt02:errEffSPDPt02:EffSPDPt1:errEffSPDPt1:EffSPDPt10:errEffSPDPt10:EffoneSPDPt02:errEffoneSPDPt02:EffoneSPDPt1:errEffoneSPDPt1:EffoneSPDPt10:errEffoneSPDPt10:EffTOTPt02:errEffTOTPt02:EffTOTPt1:errEffTOTPt1:EffTOTPt10:errEffTOTPt10");

  Float_t xntMatching[nVariablesMatching];
  
  //--------------------------------

  //----QA Vertex

      const Int_t nVariablesVertex=15;
      TNtuple* ntvertex=new TNtuple("ntvertex","QA Vertex","nrun:Vx:errVx:sigmaVx:errsigmaVx:Vy:errVy:sigmaVy:errsigmaVy:Vz:errVz:sigmaVz:errsigmaVz");
      
      Float_t xntVertex[nVariablesVertex];
  
  //--------------------------------

  //-----  Tracking ITS SA
  const Int_t nVariablesSA=13;
  Float_t xntSA[nVariablesSA];
  TNtuple *ntSA=new TNtuple("ntITSsa","ntITSsa","run:NITSTPCPtBin0:NITSTPCPtBin1:NITSTPCPtBin2:NITSsaPtBin0:NITSsaPtBin1:NITSsaPtBin2:NITSpureSAPtBin0:NITSpureSAPtBin1:NITSpureSAPtBin2:ratioPtBin0:ratioPtBin1:ratioPtBin2");

  //---------------------------------


  
  TBits* readRun=new TBits(999999);
  readRun->ResetAllBits();
  //    if(!useExternalList){
  if(!gSystem->Exec(Form("ls -l %s > /dev/null 2>&1",ntupleFileName.Data()))){
    TFile* oldfil=new TFile(ntupleFileName.Data());
    
    TNtuple* ntmp=(TNtuple*)oldfil->Get("ntsdd");
    
    TNtuple* ntmpSSD=(TNtuple*)oldfil->Get("ntssd");

    TNtuple* ntmpMatching=(TNtuple*)oldfil->Get("ntmatching");

    TNtuple* ntmpVertex=(TNtuple*)oldfil->Get("ntvertex");

    TNtuple* ntmpSA=(TNtuple*)oldfil->Get("ntITSsa");

    
    //-------SDD
    
    Bool_t isOK=kFALSE;
    if(ntmp){
      if(ntmp->GetNvar()==ntsdd->GetNvar()){
	isOK=kTRUE;
	TObjArray* arr1=(TObjArray*)ntsdd->GetListOfBranches();
	TObjArray* arr2=(TObjArray*)ntmp->GetListOfBranches();
	  for(Int_t iV=0; iV<ntmp->GetNvar(); iV++){
	    TString vnam1=arr1->At(iV)->GetName();
	    TString vnam2=arr2->At(iV)->GetName();
	    if(vnam1!=vnam2) isOK=kFALSE;
	    ntmp->SetBranchAddress(vnam2.Data(),&xnt[iV]);
	  }
	  if(isOK){
	    for(Int_t nE=0; nE<ntmp->GetEntries(); nE++){
	      ntmp->GetEvent(nE);
	      Int_t theRun=(Int_t)(xnt[0]+0.0001);
	      readRun->SetBitNumber(theRun);
	      ntsdd->Fill(xnt);
	    }
	  }
	}
      }
      if(!isOK){
	printf("Ntuple in local file not OK -> will be recreated\n");
      }
	
      //----------SSD----------
	
      Bool_t isOKSSD=kFALSE;
      if(ntmpSSD){
	if(ntmpSSD->GetNvar()==ntssd->GetNvar()){
	  isOKSSD=kTRUE;
	  TObjArray* arr1ssd=(TObjArray*)ntssd->GetListOfBranches();
	  TObjArray* arr2ssd=(TObjArray*)ntmpSSD->GetListOfBranches();
	  for(Int_t iV=0; iV<ntmpSSD->GetNvar(); iV++){
	    TString vnam1=arr1ssd->At(iV)->GetName();
	    TString vnam2=arr2ssd->At(iV)->GetName();
	    if(vnam1!=vnam2) isOKSSD=kFALSE;
	    ntmpSSD->SetBranchAddress(vnam2.Data(),&xntSSD[iV]);
	  }
	  if(isOKSSD){
	    for(Int_t nE=0; nE<ntmpSSD->GetEntries(); nE++){
	      ntmpSSD->GetEvent(nE);
	      Int_t theRun=(Int_t)(xntSSD[0]+0.0001);
	      readRun->SetBitNumber(theRun);
	      ntssd->Fill(xntSSD);
	    }
	  }
	}
      }
      if(!isOKSSD){
	printf("\n\nNtuple SSD in local file not OK -> will be recreated\n\n");
      }

      //---------Matching---------
	
      Bool_t isOKMatching=kFALSE;
      if(ntmpMatching){
	if(ntmpMatching->GetNvar()==ntmatching->GetNvar()){
	  isOKMatching=kTRUE;
	  TObjArray* arr1matching=(TObjArray*)ntmatching->GetListOfBranches();
	  TObjArray* arr2matching=(TObjArray*)ntmpMatching->GetListOfBranches();
	  for(Int_t iV=0; iV<ntmpMatching->GetNvar(); iV++){
	    TString vnam1=arr1matching->At(iV)->GetName();
	    TString vnam2=arr2matching->At(iV)->GetName();
	    if(vnam1!=vnam2) isOKMatching=kFALSE;
	    ntmpMatching->SetBranchAddress(vnam2.Data(),&xntMatching[iV]);
	  }
	  if(isOKMatching){
	    for(Int_t nE=0; nE<ntmpMatching->GetEntries(); nE++){
	      ntmpMatching->GetEvent(nE);
	      Int_t theRun=(Int_t)(xntMatching[0]+0.0001);
	      readRun->SetBitNumber(theRun);
	      ntmatching->Fill(xntMatching);
	    }
	  }
	}
      }
      if(!isOKMatching){
	printf("\n\nNtuple Matching in local file not OK -> will be recreated\n\n");
      }
      //-----------------------

      //---------Vertex QA---------
	
      Bool_t isOKVertex=kFALSE;
      if(ntmpVertex){
	if(ntmpVertex->GetNvar()==ntvertex->GetNvar()){
	  isOKVertex=kTRUE;
	  TObjArray* arr1vertex=(TObjArray*)ntvertex->GetListOfBranches();
	  TObjArray* arr2vertex=(TObjArray*)ntmpVertex->GetListOfBranches();
	  for(Int_t iV=0; iV<ntmpVertex->GetNvar(); iV++){
	    TString vnam1=arr1vertex->At(iV)->GetName();
	    TString vnam2=arr2vertex->At(iV)->GetName();
	    if(vnam1!=vnam2) isOKVertex=kFALSE;
	    ntmpVertex->SetBranchAddress(vnam2.Data(),&xntVertex[iV]);
	  }
	  if(isOKVertex){
	    for(Int_t nE=0; nE<ntmpVertex->GetEntries(); nE++){
	      ntmpVertex->GetEvent(nE);
	      Int_t theRun=(Int_t)(xntVertex[0]+0.0001);
	      readRun->SetBitNumber(theRun);
	      ntvertex->Fill(xntVertex);
	    }
	  }
	}
      }
      if(!isOKVertex){
	printf("\n\nNtuple Vertex in local file not OK -> will be recreated\n\n");
      }
      //-----------------------

      //-----------Tracking SA ----------------
      Bool_t isOKSA = kFALSE;
      if(ntmpSA){
	if(ntmpSA->GetNvar()==ntSA->GetNvar()){
	  isOKSA = kTRUE;
	  TObjArray* arr1SA=(TObjArray*)ntSA->GetListOfBranches();
	  TObjArray* arr2SA=(TObjArray*)ntmpSA->GetListOfBranches();
	  for(Int_t iV=0; iV<ntmpSA->GetNvar(); iV++){
	    TString vnam1=arr1SA->At(iV)->GetName();
	    TString vnam2=arr2SA->At(iV)->GetName();
	    if(vnam1!=vnam2) isOKSA=kFALSE;
	    ntmpSA->SetBranchAddress(vnam2.Data(),&xntSA[iV]);
	  }
	  if(isOKSA){
	    for(Int_t nE=0; nE<ntmpSA->GetEntries(); nE++){
	      ntmpSA->GetEvent(nE);
	      Int_t theRun=(Int_t)(xntSA[0]+0.0001);
	      readRun->SetBitNumber(theRun);
	      ntSA->Fill(xntSA);
	    }
	  }
	}
      }

      oldfil->Close();
      delete oldfil;
  }

#define MAX_LINES 200
#define MAX_LINE_LEN 255
  
  char strings[MAX_LINES][MAX_LINE_LEN];
  ifstream in(runListFile.Data());
  int j = 0;
  Int_t nrun=0;
  Int_t runNumb[MAX_LINES];
  while ( in ) {
    in.getline(strings[j], MAX_LINE_LEN);
    TString aux(strings[j]);
    if(aux.Length()<27)continue;
    aux=aux.Remove(0,27);
    aux=aux.Remove(6,aux.Length());  
    runNumb[j]=atoi(aux.Data());
    printf("%d ) - path %s \n",runNumb[j],strings[j]);
    j++;
    nrun++;
  }

  printf("\n *******************   Loop on runs *********** \n");
  for(Int_t jru=0;jru<nrun;jru++) {
    printf("jru=%d - run number= %d \n",jru,runNumb[jru]);
    Int_t iRun=runNumb[jru];
    if(readRun->TestBitNumber(iRun))printf("Run %d - already processed\n",iRun);
    if(readRun->TestBitNumber(iRun))continue;
    //cout << "Value from file is " <<t << endl;
    
    printf("%s\n",strings[jru]);
  
  
    if(!gGrid||!gGrid->IsConnected()) {
      printf("gGrid not found! exit macro\n");
      return;
    }
    
    TFile *f=TFile::Open(Form("alien://%s",strings[jru])); 
    
    TDirectoryFile* df=(TDirectoryFile*)f->Get("SDD_Performance");
    if(!df){
      printf("Run %d SDD_Performance MISSING -> Exit\n",iRun);
      continue;
    }
    
    TList* l=(TList*)df->Get("coutputRP");
    if(!df){
      printf("Run %d coutputRP TList MISSING -> Exit\n",iRun);
      continue;
    }  
    
    //-------------------

    
    //------------SDD

    TH1F* hcllay=(TH1F*)l->FindObject("hCluInLay");
      Float_t fracT[6]={0.,0.,0.,0.,0.,0.};
      Float_t efracT[6]={0.,0.,0.,0.,0.,0.};
      if(hcllay->GetBinContent(1)>0){
	for(Int_t iLay=0; iLay<6; iLay++){
	  fracT[iLay]=hcllay->GetBinContent(iLay+2)/hcllay->GetBinContent(1);
	  efracT[iLay]=TMath::Sqrt(fracT[iLay]*(1-fracT[iLay])/hcllay->GetBinContent(1));
	}
      }
      TH1F* hmodT=(TH1F*)l->FindObject("hTPMod");
      TH1F* hgamod=(TH1F*)l->FindObject("hGAMod");
      Int_t bestMod=0;
      for(Int_t iMod=0; iMod<260;iMod++){
	Int_t gda=(Int_t)hgamod->GetBinContent(iMod+1);
	if(gda>bestMod) bestMod=gda;
      }
      Int_t nChunks=1;
      if(bestMod>512){
	nChunks=(Int_t)(bestMod/512.+0.5);
      }
      hgamod->Scale(1./nChunks);

      TH1F* hev=(TH1F*)l->FindObject("hNEvents");
      Int_t nTotEvents=hev->GetBinContent(2);
      Int_t nTrigEvents=hev->GetBinContent(3);
      Int_t nEvents=nTotEvents;
      printf("Run %d Number of Events = %d Triggered=%d\n",iRun,nTotEvents,nTrigEvents);
      if(nTrigEvents>0){ 
	nEvents=nTrigEvents;
      }
      if(nTotEvents==0) continue;
      Int_t nModGood3=0;
      Int_t nModGood4=0;
      Int_t nModBadAn=0;
      Float_t sumtp3=0;
      Float_t sumtp4=0;
      Float_t sumEtp3=0;
      Float_t sumEtp4=0;
      for(Int_t iMod=0; iMod<260; iMod++){
	Float_t tps=hmodT->GetBinContent(iMod+1);
	Float_t ga=hgamod->GetBinContent(iMod+1);
	if(ga<500) nModBadAn++;
	Float_t tpsN=0.;
	Float_t etpsN=0.;
	if(ga>0){
	  tpsN=tps/ga/(Float_t)nEvents;
	  etpsN=TMath::Sqrt(tps)/ga/(Float_t)nEvents;
	  if(iMod<84){
	    sumtp3+=tpsN;
	    sumEtp3+=(etpsN*etpsN);
	    nModGood3++;
	  }
	  else{
	    sumtp4+=tpsN;
	    sumEtp4+=(etpsN*etpsN);
	    nModGood4++;
	  }
	}
      }

      TH1F* hapmod=(TH1F*)l->FindObject("hAllPmod");
      TH1F* hgpmod=(TH1F*)l->FindObject("hGoodPmod");
      //     TH1F* hmpmod=(TH1F*)l->FindObject("hMissPmod");
      TH1F* hbrmod=(TH1F*)l->FindObject("hBadRegmod");
      TH1F* hskmod=(TH1F*)l->FindObject("hSkippedmod");
      TH1F* hoamod=(TH1F*)l->FindObject("hOutAccmod");
      TH1F* hnrmod=(TH1F*)l->FindObject("hNoRefitmod");
      Int_t nBelow95=0;
      Int_t nBelow80=0;
      Int_t nBelow60=0;
      Int_t nZeroP=0;
      for(Int_t imod=0; imod<260;imod++){
	Float_t numer=hgpmod->GetBinContent(imod+1)+hbrmod->GetBinContent(imod+1)+hoamod->GetBinContent(imod+1)+hnrmod->GetBinContent(imod+1)+hskmod->GetBinContent(imod+1);
	Float_t denom=hapmod->GetBinContent(imod+1);
	if(denom>0){
	  Float_t eff=numer/denom;
	  if(eff<0.95) nBelow95++;
	  if(eff<0.80) nBelow80++;
	  if(eff<0.60) nBelow60++;
	}
	if(hmodT->GetBinContent(imod+1)<1.){
	  nZeroP++;
	}	
      }

      TH1F* htimT=(TH1F*)l->FindObject("hDrTimTPAll");
      TH1F* htimTe=(TH1F*)l->FindObject("hDrTimTPExtra");
      
      Double_t fracExtra=0.;
      Double_t errFracExtra=0.;
      if(htimT->GetEntries()>0){
	fracExtra=htimTe->GetEntries()/htimT->GetEntries();
	errFracExtra=TMath::Sqrt(htimTe->GetEntries())/htimT->GetEntries();
      }
      Double_t averPoints=0.;
      Double_t cntBins=0.;
      for(Int_t iBin=1; iBin<=htimT->GetNbinsX(); iBin++){
	Float_t tim=htimT->GetBinCenter(iBin);
	if(tim>2000. && tim<4000.){
	  averPoints+=htimT->GetBinContent(iBin);
	  cntBins+=1;
	}
      }
      Double_t minTime=-999.;
      Double_t errMinTime=0.;
      if(cntBins>0){ 
	averPoints/=cntBins;
	for(Int_t iBin=1; iBin<=htimT->GetNbinsX(); iBin++){
	  if(htimT->GetBinContent(iBin)>0.5*averPoints){
	    minTime=htimT->GetBinCenter(iBin);
	    errMinTime=0.5*htimT->GetBinWidth(iBin);
	    break;
	  }
	}
      }

      TH2F* hdedxmod=(TH2F*)l->FindObject("hdEdxVsMod");
      TH1D* hdedxLay3=hdedxmod->ProjectionY("hdedxLay3",1,84);
      TH1D* hdedxLay4=hdedxmod->ProjectionY("hdedxLay4",85,260);
      
      TH1F* hSigTim0=(TH1F*)l->FindObject("hSigTimeInt0");
      TH1F* hSigTim5=(TH1F*)l->FindObject("hSigTimeInt5");

      Int_t index=0;
      xnt[index++]=iRun;
      xnt[index++]=fracT[0];
      xnt[index++]=efracT[0];
      xnt[index++]=fracT[1];
      xnt[index++]=efracT[1];
      xnt[index++]=fracT[2];
      xnt[index++]=efracT[2];
      xnt[index++]=fracT[3];
      xnt[index++]=efracT[3];
      xnt[index++]=fracT[4];
      xnt[index++]=efracT[4];
      xnt[index++]=fracT[5];
      xnt[index++]=efracT[5];
      xnt[index++]=sumtp3/nModGood3;
      xnt[index++]=TMath::Sqrt(sumEtp3)/nModGood3;
      xnt[index++]=sumtp4/nModGood4;
      xnt[index++]=TMath::Sqrt(sumEtp4)/nModGood4;
      xnt[index++]=minTime;
      xnt[index++]=errMinTime;
      xnt[index++]=htimT->GetMean();
      xnt[index++]=htimT->GetMeanError();
      xnt[index++]=fracExtra;
      xnt[index++]=errFracExtra;
      xnt[index++]=hdedxLay3->GetMean();
      xnt[index++]=hdedxLay3->GetMeanError();
      xnt[index++]=hdedxLay4->GetMean();
      xnt[index++]=hdedxLay4->GetMeanError();
      xnt[index++]=hSigTim0->GetMean();
      xnt[index++]=hSigTim0->GetMeanError();
      xnt[index++]=hSigTim5->GetMean();
      xnt[index++]=hSigTim5->GetMeanError();
      xnt[index++]=(Float_t)nBelow95;
      xnt[index++]=(Float_t)nBelow80;
      xnt[index++]=(Float_t)nBelow60;
      xnt[index++]=(Float_t)nZeroP;
      ntsdd->Fill(xnt);

      cout<<"\n\nirun sDD"<<iRun<<endl<<endl;

      //-------SSD

      //TFile* fSSD=TFile::Open(fileNameLong.Data());  
    
      //TDirectoryFile* dfSSD=(TDirectoryFile*)fSSD->Get("PWG1dEdxSSDQA");
      TDirectoryFile* dfSSD=(TDirectoryFile*)f->Get("PWG1dEdxSSDQA");
      if(!dfSSD){
	printf("Run %d SDD_Performance MISSING -> Exit\n",iRun);
	continue;
      }
      
      TList* lSSD=(TList*)dfSSD->Get("SSDdEdxQA");
      if(!dfSSD){
	printf("Run %d coutputRP TList MISSING -> Exit\n",iRun);
	continue;
      }
      //
      
      TH2F* QAchargeRatio=(TH2F*)lSSD->FindObject("QAChargeRatio");
      TH2F* QAcharge=(TH2F*)lSSD->FindObject("QACharge");
    
      Int_t biny = QAcharge->GetXaxis()->FindBin(747);
      Int_t maxy = QAcharge->GetXaxis()->GetXmax();
      //     Int_t miny = QAcharge->GetXaxis()->GetXmin();
    
      Int_t  contEmpty=0;
      Int_t  contFull=0;
    
      TH1D *hChargeL5=QAcharge->ProjectionY("hChargeL5",0,biny);
      TH1D *hChargeL6=QAcharge->ProjectionY("hChargeL6",biny,maxy);
    
      cout<<  hChargeL5->GetMean()<< " " <<  hChargeL5->GetRMS()<<endl;
      cout<<  hChargeL6->GetMean()<< " " <<  hChargeL6->GetRMS()<<endl;
    
      TH1D *hChargeRatioL5=QAchargeRatio->ProjectionY("hChargeRatioL5",0,biny);
      TH1D *hChargeRatioL6=QAchargeRatio->ProjectionY("hChargeRatioL6",biny,maxy);
    
      cout<<  hChargeRatioL5->GetMean()<< " " <<  hChargeRatioL5->GetRMS()<<endl;
      cout<<  hChargeRatioL6->GetMean()<< " " <<  hChargeRatioL6->GetRMS()<<endl;
    
      if(QAcharge->GetEntries()< 45000)
	contEmpty=1;
      
      else{
	for(Int_t i =0;i<1698;i++){
	  
	  TString tmpQ("Q");
	  tmpQ+=i;
	  
	  TH1D* fHist1DQ= QAcharge->ProjectionY(tmpQ,i+1,i+1);
	  Double_t mean=fHist1DQ->GetMean();

	  if(TMath::Abs(mean)<1.0 ||fHist1DQ->GetEntries()<10)
	    contEmpty++;
	  
	  else 
	    contFull++;
	  
	}
      }

      cout<<"contFull: " <<contFull<<" contEmpty: "<<contEmpty<<endl;
      cout<<hChargeL5->GetMean()<<endl;

      Int_t indexSSD=0;
      xntSSD[indexSSD++]=iRun;
      xntSSD[indexSSD++]=(Float_t)hChargeL5->GetMean();
      xntSSD[indexSSD++]=(Float_t)hChargeL5->GetMeanError();
      xntSSD[indexSSD++]=(Float_t)hChargeL6->GetMean();
      xntSSD[indexSSD++]=(Float_t)hChargeL6->GetMeanError();
      xntSSD[indexSSD++]=(Float_t)hChargeRatioL5->GetMean();
      xntSSD[indexSSD++]=(Float_t)hChargeRatioL5->GetMeanError();
      xntSSD[indexSSD++]=(Float_t)hChargeRatioL6->GetMean();
      xntSSD[indexSSD++]=(Float_t)hChargeRatioL6->GetMeanError();
      xntSSD[indexSSD++]=(Float_t)contEmpty;
      ntssd->Fill(xntSSD);						

      cout<<xnt<<endl<<endl;

      cout<<"\n\nirun ssd "<<iRun<<endl<<endl;

      cout<< iRun<<" "<<hChargeL5->GetMean()<<" "<< hChargeL5->GetRMS()<<" "<<hChargeL6->GetMean()<<" "<<hChargeRatioL5->GetMean()<<" "<<" "<<hChargeRatioL5->GetRMS()<<" "<<hChargeRatioL6->GetMean()<<" "<<hChargeRatioL6->GetRMS()<<" "<<contEmpty<<endl;

      cout<<xntSSD[0]<<" "<<xntSSD[1]<<" "<<xntSSD[2]<<" "<<xntSSD[3]<<endl;;
      
      //--------------matching

      TDirectoryFile *dirMatch=(TDirectoryFile*)f->GetDirectory("ITS_Performance");
      TList *list=NULL;
      TList *listSPD=NULL;
      
      if(dirMatch) list = (TList*)dirMatch->Get("cOutputITS_3500_10000");
      dirMatch=(TDirectoryFile*)f->GetDirectory("SPD_Performance");
      if(dirMatch) listSPD = (TList*)dirMatch->Get("coutput1");
  
      // if(!list) return kFALSE;
      
      Float_t ioValues[30];
      Float_t ioErrors[30];
      for(Int_t jj=0;jj<30;jj++){
	ioValues[jj]=0.;
	ioErrors[jj]=0.;
      }
      
      Float_t ptbin=0;
      
      TH1F *hFiredChip = (TH1F*)listSPD->FindObject("hFiredChip");
      Int_t nHSsInner=0,nHSsOuter=0;
      for(Int_t i=0;i<400;i++) if(hFiredChip->GetBinContent(i)>0) nHSsInner++;
      for(Int_t i=400;i<1200;i++) if(hFiredChip->GetBinContent(i)>0) nHSsOuter++;
      nHSsInner = (Int_t)(nHSsInner/10);
      nHSsOuter = (Int_t)(nHSsOuter/10);
      // hnHSsSPD->SetBinContent(1,nHSsInner);
  // hnHSsSPD->SetBinContent(2,nHSsOuter);
      
      ioValues[0]=(Float_t)nHSsInner/40.;
      ioValues[1]=(Float_t)nHSsOuter/80.;

      TH1F *fHistPtTPCInAcc = (TH1F*)list->FindObject("fHistPtTPCInAcc");
      
      TH1F *fHistPtITSMI6InAcc = (TH1F*)list->FindObject("fHistPtITSMI6InAcc");
      TH1F *fHistPtITSMI5InAcc = (TH1F*)list->FindObject("fHistPtITSMI5InAcc");
      TH1F *fHistPtITSMI4InAcc = (TH1F*)list->FindObject("fHistPtITSMI4InAcc");
      TH1F *fHistPtITSMI3InAcc = (TH1F*)list->FindObject("fHistPtITSMI3InAcc");
      TH1F *fHistPtITSMI2InAcc = (TH1F*)list->FindObject("fHistPtITSMI2InAcc");
      TH1F *fHistPtITSMISPDInAcc = (TH1F*)list->FindObject("fHistPtITSMISPDInAcc");
      TH1F *fHistPtITSMIoneSPDInAcc = (TH1F*)list->FindObject("fHistPtITSMIoneSPDInAcc");
      

      TH1F *fHistPtITSMIge2InAcc = (TH1F*)fHistPtITSMI6InAcc->Clone("fHistPtITSMIge2InAcc");
      fHistPtITSMIge2InAcc->Add(fHistPtITSMI5InAcc);
      fHistPtITSMIge2InAcc->Add(fHistPtITSMI4InAcc);
      fHistPtITSMIge2InAcc->Add(fHistPtITSMI3InAcc);
      fHistPtITSMIge2InAcc->Add(fHistPtITSMI2InAcc);

      

      fHistPtITSMI6InAcc->Divide(fHistPtITSMI6InAcc,fHistPtTPCInAcc,1,1,"B");
      ptbin=fHistPtITSMI6InAcc->FindBin(0.201);
      ioValues[2]=fHistPtITSMI6InAcc->GetBinContent(ptbin);
      ioErrors[2]=fHistPtITSMI6InAcc->GetBinError(ptbin);
      ptbin=fHistPtITSMI6InAcc->FindBin(1.001);
      ioValues[3]=fHistPtITSMI6InAcc->GetBinContent(ptbin);
      ioErrors[3]=fHistPtITSMI6InAcc->GetBinError(ptbin);
      ptbin=fHistPtITSMI6InAcc->FindBin(10.001);
      ioValues[4]=fHistPtITSMI6InAcc->GetBinContent(ptbin);
      ioErrors[4]=fHistPtITSMI6InAcc->GetBinError(ptbin);

      fHistPtITSMI5InAcc->Divide(fHistPtITSMI5InAcc,fHistPtTPCInAcc,1,1,"B");
      
      ptbin=fHistPtITSMI5InAcc->FindBin(0.201);
      ioValues[5]=fHistPtITSMI5InAcc->GetBinContent(ptbin);
      ioErrors[5]=fHistPtITSMI5InAcc->GetBinError(ptbin);
      ptbin=fHistPtITSMI5InAcc->FindBin(1.001);
      ioValues[6]=fHistPtITSMI5InAcc->GetBinContent(ptbin);
      ioErrors[6]=fHistPtITSMI5InAcc->GetBinError(ptbin);
      ptbin=fHistPtITSMI5InAcc->FindBin(10.001);
      ioValues[7]=fHistPtITSMI5InAcc->GetBinContent(ptbin);
      ioErrors[7]=fHistPtITSMI5InAcc->GetBinError(ptbin);
      
      fHistPtITSMI4InAcc->Divide(fHistPtITSMI4InAcc,fHistPtTPCInAcc,1,1,"B");
      
      ptbin=fHistPtITSMI4InAcc->FindBin(0.201);
      ioValues[9]=fHistPtITSMI4InAcc->GetBinContent(ptbin);
  ioErrors[9]=fHistPtITSMI4InAcc->GetBinError(ptbin);
  ptbin=fHistPtITSMI4InAcc->FindBin(1.001);
  ioValues[10]=fHistPtITSMI4InAcc->GetBinContent(ptbin);
  ioErrors[10]=fHistPtITSMI4InAcc->GetBinError(ptbin);
  ptbin=fHistPtITSMI4InAcc->FindBin(10.001);
  ioValues[11]=fHistPtITSMI4InAcc->GetBinContent(ptbin);
  ioErrors[11]=fHistPtITSMI4InAcc->GetBinError(ptbin);

  fHistPtITSMI3InAcc->Divide(fHistPtITSMI3InAcc,fHistPtTPCInAcc,1,1,"B");

  ptbin=fHistPtITSMI3InAcc->FindBin(0.201);
  ioValues[12]=fHistPtITSMI3InAcc->GetBinContent(ptbin);
  ioErrors[12]=fHistPtITSMI3InAcc->GetBinError(ptbin);
  ptbin=fHistPtITSMI3InAcc->FindBin(1.001);
  ioValues[13]=fHistPtITSMI3InAcc->GetBinContent(ptbin);
  ioErrors[13]=fHistPtITSMI3InAcc->GetBinError(ptbin);
  ptbin=fHistPtITSMI3InAcc->FindBin(10.001);
  ioValues[14]=fHistPtITSMI3InAcc->GetBinContent(ptbin);
  ioErrors[14]=fHistPtITSMI3InAcc->GetBinError(ptbin);

  fHistPtITSMI2InAcc->Divide(fHistPtITSMI2InAcc,fHistPtTPCInAcc,1,1,"B");

  ptbin=fHistPtITSMI2InAcc->FindBin(0.201);
  ioValues[15]=fHistPtITSMI2InAcc->GetBinContent(ptbin);
  ioErrors[15]=fHistPtITSMI2InAcc->GetBinError(ptbin);
  ptbin=fHistPtITSMI2InAcc->FindBin(1.001);
  ioValues[16]=fHistPtITSMI2InAcc->GetBinContent(ptbin);
  ioErrors[16]=fHistPtITSMI2InAcc->GetBinError(ptbin);
  ptbin=fHistPtITSMI2InAcc->FindBin(10.001);
  ioValues[17]=fHistPtITSMI2InAcc->GetBinContent(ptbin);
  ioErrors[17]=fHistPtITSMI2InAcc->GetBinError(ptbin);

  fHistPtITSMISPDInAcc->Divide(fHistPtITSMISPDInAcc,fHistPtTPCInAcc,1,1,"B");
  ptbin=fHistPtITSMISPDInAcc->FindBin(0.201);
  ioValues[18]=fHistPtITSMISPDInAcc->GetBinContent(ptbin);
  ioErrors[18]=fHistPtITSMISPDInAcc->GetBinError(ptbin);
  ptbin=fHistPtITSMISPDInAcc->FindBin(1.001);
  ioValues[19]=fHistPtITSMISPDInAcc->GetBinContent(ptbin);
  ioErrors[19]=fHistPtITSMISPDInAcc->GetBinError(ptbin);
  ptbin=fHistPtITSMISPDInAcc->FindBin(10.001);
  ioValues[20]=fHistPtITSMISPDInAcc->GetBinContent(ptbin);
  ioErrors[20]=fHistPtITSMISPDInAcc->GetBinError(ptbin);

  fHistPtITSMIoneSPDInAcc->Divide(fHistPtITSMIoneSPDInAcc,fHistPtTPCInAcc,1,1,"B");

  ptbin=fHistPtITSMIoneSPDInAcc->FindBin(0.201);
  ioValues[21]=fHistPtITSMIoneSPDInAcc->GetBinContent(ptbin);
  ioErrors[21]=fHistPtITSMIoneSPDInAcc->GetBinError(ptbin);
  ptbin=fHistPtITSMIoneSPDInAcc->FindBin(1.001);
  ioValues[22]=fHistPtITSMIoneSPDInAcc->GetBinContent(ptbin);
  ioErrors[22]=fHistPtITSMIoneSPDInAcc->GetBinError(ptbin);
  ptbin=fHistPtITSMIoneSPDInAcc->FindBin(10.001);
  ioValues[23]=fHistPtITSMIoneSPDInAcc->GetBinContent(ptbin);
  ioErrors[23]=fHistPtITSMIoneSPDInAcc->GetBinError(ptbin);

  
  fHistPtITSMIge2InAcc->Divide(fHistPtITSMIge2InAcc,fHistPtTPCInAcc,1,1,"B");
  ptbin=fHistPtITSMIge2InAcc->FindBin(0.201);
  ioValues[24]=fHistPtITSMIge2InAcc->GetBinContent(ptbin);
  ioErrors[24]=fHistPtITSMIge2InAcc->GetBinError(ptbin);
  ptbin=fHistPtITSMIge2InAcc->FindBin(1.001);
  ioValues[25]=fHistPtITSMIge2InAcc->GetBinContent(ptbin);
  ioErrors[25]=fHistPtITSMIge2InAcc->GetBinError(ptbin);
  ptbin=fHistPtITSMIge2InAcc->FindBin(10.001);
  ioValues[26]=fHistPtITSMIge2InAcc->GetBinContent(ptbin);
  ioErrors[26]=fHistPtITSMIge2InAcc->GetBinError(ptbin);

  Int_t indexMatching=0;
  xntMatching[indexMatching++]=iRun;
  xntMatching[indexMatching++]=ioValues[0];
  xntMatching[indexMatching++]=ioErrors[0];
  xntMatching[indexMatching++]=ioValues[1];
  xntMatching[indexMatching++]=ioErrors[1];
  xntMatching[indexMatching++]=ioValues[2];
  xntMatching[indexMatching++]=ioErrors[2];
  xntMatching[indexMatching++]=ioValues[3];
  xntMatching[indexMatching++]=ioErrors[3];
  xntMatching[indexMatching++]=ioValues[4];
  xntMatching[indexMatching++]=ioErrors[4];
  xntMatching[indexMatching++]=ioValues[5];
  xntMatching[indexMatching++]=ioErrors[5];
  xntMatching[indexMatching++]=ioValues[6];
  xntMatching[indexMatching++]=ioErrors[6];
  xntMatching[indexMatching++]=ioValues[7];
  xntMatching[indexMatching++]=ioErrors[7];
  xntMatching[indexMatching++]=ioValues[8];
  xntMatching[indexMatching++]=ioErrors[8];
  xntMatching[indexMatching++]=ioValues[9];
  xntMatching[indexMatching++]=ioErrors[9];
  xntMatching[indexMatching++]=ioValues[10];
  xntMatching[indexMatching++]=ioErrors[10];
  xntMatching[indexMatching++]=ioValues[11];
  xntMatching[indexMatching++]=ioErrors[11];
  xntMatching[indexMatching++]=ioValues[12];
  xntMatching[indexMatching++]=ioErrors[12];
  xntMatching[indexMatching++]=ioValues[13];
  xntMatching[indexMatching++]=ioErrors[13];
  xntMatching[indexMatching++]=ioValues[14];
  xntMatching[indexMatching++]=ioErrors[14];
  xntMatching[indexMatching++]=ioValues[15];
  xntMatching[indexMatching++]=ioErrors[15];
  xntMatching[indexMatching++]=ioValues[16];
  xntMatching[indexMatching++]=ioErrors[16];
  xntMatching[indexMatching++]=ioValues[17];
  xntMatching[indexMatching++]=ioErrors[17];
  xntMatching[indexMatching++]=ioValues[18];
  xntMatching[indexMatching++]=ioErrors[18];
  xntMatching[indexMatching++]=ioValues[19];
  xntMatching[indexMatching++]=ioErrors[19];
  xntMatching[indexMatching++]=ioValues[20];
  xntMatching[indexMatching++]=ioErrors[20];
  xntMatching[indexMatching++]=ioValues[21];
  xntMatching[indexMatching++]=ioErrors[21];
  xntMatching[indexMatching++]=ioValues[22];
  xntMatching[indexMatching++]=ioErrors[22];
  xntMatching[indexMatching++]=ioValues[23];
  xntMatching[indexMatching++]=ioErrors[23];
  xntMatching[indexMatching++]=ioValues[24];
  xntMatching[indexMatching++]=ioErrors[24];
  xntMatching[indexMatching++]=ioValues[25];
  xntMatching[indexMatching++]=ioErrors[25];
  xntMatching[indexMatching++]=ioValues[26];
  xntMatching[indexMatching++]=ioErrors[26];
  
  ntmatching->Fill(xntMatching);						

      //---------------------------

  TDirectoryFile *dirVertex = (TDirectoryFile*)f->Get("Vertex_Performance");
  if(!dirVertex){
    Printf("Vertex directory not found... check!");
  }
  
  TList *lt = (TList*)dirVertex->Get("cOutputVtxESD");
	
  TH1F *xVtxTRK = (TH1F*)lt->FindObject("fhTRKVertexX");
  TH1F *yVtxTRK = (TH1F*)lt->FindObject("fhTRKVertexY");
  TH1F *zVtxTRK = (TH1F*)lt->FindObject("fhTRKVertexZ");

  TF1 *fxTRK = new TF1("gausx", "gaus", -1, 1);
  xVtxTRK->Fit("gausx", "M");

  TF1 *fyTRK = new TF1("gausy", "gaus", -1, 1);
  yVtxTRK->Fit("gausy","M");
  cout<<fyTRK->GetParameter(1)<<endl;
  cout<<fyTRK->GetParError(1)<<endl;
  cout<<fyTRK->GetParameter(2)<<endl;
  cout<<fyTRK->GetParError(2)<<endl;

  TF1 *fzTRK = new TF1("gausz", "gaus", -1, 1);
  zVtxTRK->Fit("gausz","M");
  cout<<fzTRK->GetParameter(1)<<endl;
  cout<<fzTRK->GetParError(1)<<endl;
  cout<<fzTRK->GetParameter(2)<<endl;
  cout<<fzTRK->GetParError(2)<<endl;


  Int_t indexVertex=0;
  xntVertex[indexVertex++]=iRun;
  xntVertex[indexVertex++]=(Float_t)fxTRK->GetParameter(1);
  xntVertex[indexVertex++]=(Float_t)fxTRK->GetParError(1);
  xntVertex[indexVertex++]=(Float_t)fxTRK->GetParameter(2);
  xntVertex[indexVertex++]=(Float_t)fxTRK->GetParError(2);
  xntVertex[indexVertex++]=(Float_t)fyTRK->GetParameter(1);
  xntVertex[indexVertex++]=(Float_t)fyTRK->GetParError(1);
  xntVertex[indexVertex++]=(Float_t)fyTRK->GetParameter(2);
  xntVertex[indexVertex++]=(Float_t)fyTRK->GetParError(2);
  xntVertex[indexVertex++]=(Float_t)fzTRK->GetParameter(1);
  xntVertex[indexVertex++]=(Float_t)fzTRK->GetParError(1);
  xntVertex[indexVertex++]=(Float_t)fzTRK->GetParameter(2);
  xntVertex[indexVertex++]=(Float_t)fzTRK->GetParError(2);
  ntvertex->Fill(xntVertex);	

					
  
  //---------------------------

  //--------------   ITS SA ---------------

  FillITSSAntuple(f,ntSA,iRun);

  } // loop on runs

  TFile* outfil=new TFile(ntupleFileName.Data(),"recreate");
  outfil->cd();
  ntsdd->Write();
  ntssd->Write();
  ntmatching->Write();
  ntvertex->Write();
  ntSA->Write();
  outfil->Close();
  delete outfil;
  delete ntsdd;
  delete ntssd;
  delete ntmatching;
  delete ntvertex;
  delete ntSA;

}

//____________________________________________________________________________
void FillITSSAntuple(TFile* f,TNtuple* nt, Int_t nrun){
  static const Int_t nVariables=13;
  static Float_t xnt[nVariables];
  TH1F *hPtTPCITS=0x0;
  TH1F *hPtITSsa=0x0;
  TH1F *hPtITSpureSA=0x0;
  Double_t Lowbin[3]={0.1,0.5,0.9};
  Double_t Upbin[3]={0.2,0.6,1};
  Double_t NTPCITS[3];
  Double_t NITSsa[3];
  Double_t NITSpureSA[3];
  Double_t Ratio[3];
  TDirectory *dirFile=(TDirectory*)f->Get("ITSsaTracks");
  TList *cOutput = (TList*)dirFile->Get("clistITSsaTracks"); 
  // histogram with number of events: in the first cell there is the number of the events
  // in the second, the number of events with SPD vertex.
  // normalization will be done to the second number
  TH1F *hnev =(TH1F*)cOutput->FindObject("hNEvents");
  Double_t noEvents = hnev->GetBinContent(1);
  if(noEvents<1.)noEvents=1.;   // protection to avoid division by zero
  hPtTPCITS=(TH1F*)cOutput->FindObject("hPtTPCITS");
  hPtITSsa=(TH1F*)cOutput->FindObject("hPtITSsa");
  hPtITSpureSA=(TH1F*)cOutput->FindObject("hPtITSpureSA");

  for(Int_t ibin=0;ibin<=2;ibin++){
    NTPCITS[ibin]=hPtTPCITS->Integral(hPtTPCITS->FindBin(Lowbin[ibin]),hPtTPCITS->FindBin(Upbin[ibin]))/noEvents;
    NITSsa[ibin]=hPtITSsa->Integral(hPtITSsa->FindBin(Lowbin[ibin]),hPtITSsa->FindBin(Upbin[ibin]))/noEvents;
    NITSpureSA[ibin]=hPtITSpureSA->Integral(hPtITSpureSA->FindBin(Lowbin[ibin]),hPtITSpureSA->FindBin(Upbin[ibin]))/noEvents;
    if(NTPCITS[ibin]!=0 && NITSsa[ibin]!=0)Ratio[ibin]=NTPCITS[ibin]/NITSsa[ibin];
    else Ratio[ibin]=0;
  }
 
  Int_t index=0;
  xnt[index++]=(Float_t)nrun;
  xnt[index++]=NTPCITS[0];
  xnt[index++]=NTPCITS[1];
  xnt[index++]=NTPCITS[2];
  xnt[index++]=NITSsa[0];
  xnt[index++]=NITSsa[1];
  xnt[index++]=NITSsa[2];
  xnt[index++]=NITSpureSA[0];
  xnt[index++]=NITSpureSA[1];
  xnt[index++]=NITSpureSA[2];
  xnt[index++]=Ratio[0];
  xnt[index++]=Ratio[1];
  xnt[index++]=Ratio[2];
  nt->Fill(xnt);
}
