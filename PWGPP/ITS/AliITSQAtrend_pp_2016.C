#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TH1F.h>
#include <TH1D.h>
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
#include <iostream>
#include <fstream>
//#include <AliITSgeomTGeo.h>
#endif
  
TString pdfFileNames=" ";
void MakePlot(Int_t run1=-1,Int_t run2=999999,TString ntupleFileName="TrendingITS_pp_2016.root");
//void PlotITSSA(TFile *fil,Int_t run1, Int_t run2);
void FillITSSAntuple(TFile* f,TNtuple* nt, Int_t nrun);
void FillSDDntuple(TFile* f,TNtuple* nt, Int_t iRun, Float_t *xnt);
void FillSSDntuple(TFile* f,TNtuple* ntssd, Int_t iRun, Float_t *xntSSD);
void FillMatchntuple(TFile* f,TNtuple* ntmatching, Int_t iRun, Float_t *xntMatching);
void FillMatchntupleTOF(TFile* f,TNtuple* ntmatchingTOF, Int_t iRun, Float_t *xntMatchingTOF);
void FillVTXntuple(TFile* f,TNtuple* ntvertex, Int_t iRun, Float_t *xntVertex);
void FillPileUpntuple(TFile* f,TNtuple* ntpileup, Int_t iRun, Float_t *xntPileup);
void AliITSQAtrend_pp_2016(TString runListFile="lista_pp2_2016.txt",TString ntupleFileName="TrendingITS_pp_2016.root");
Double_t LangausFun(Double_t *x, Double_t *par);
Bool_t WriteInputTextFileFromMonalisaListOfRuns(TString outtxtfilename,Int_t* listofrunsfromMonalisa,Int_t nruns,TString pathbeforRunN="alice/data/2012/LHC12a/",TString pathafterRunN="pass1/QAresults.root");
Int_t RunsToPlot(Int_t run1, Int_t run2,Int_t nr,Int_t *noRuns, Int_t *myIndex);
ofstream myfile("outfile1_pp_fast.txt",ios::app);

////////////////////////////////////////////////////////////////
//    THIS VERSION OPTIMIZED FOR PP
//    TO BE MODIFIED WHEN PASSING FROM cpass1 to pass1  
//
//   Please, read this comment before using this macro 
//
//   INPUT FILE: a text file (by default LHC12a.txt) which contains
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
//   when you add a new entry in the LHC12a.txt file (see below) and you have a
//   local copy of the  TrendingITS.root file, only the additional run will
//   be processed. The whole list is processed only the first time you use the
//   macro. Please, bear in mind that this macro is RAM-intensive: all the
//   ntuples are kept in memory. It is better to add few runs to the list at 
//   each time, according to the RAM of your computer. 
//   The function AliITSQAtrend does not produce any plot.
//
//   The LHC12a.txt file contains the QA input file name, with its path and
//   the number of the LHC fill.
//   Example of LHC12a.txt (input file complete path + FILL number):
/*
/alice/data/2012/LHC12a/000176661/cpass1/QAresults_barrel.root  2469
/alice/data/2012/LHC12a/000176701/cpass1/QAresults_barrel.root  2470
*/
//
//   Function MakePlot(run1,run2):
//   it produces the plots. For each canvas a PDF file is created.
//   A PDF file with all the canvases merged is also produced
//   The first two argument define a range for the runs to be displayed
//   These two arguments are optional: by default, all the runs
//   found in the ntuples are displayed
////////////////////////////////////////////////////////////////

/* $Id: AliITSQAtrend_pp.C 57226 2012-06-18 08:06:23Z masera $ */

void MakePlot(Int_t run1,Int_t run2,TString ntupleFileName){
  // Check run range
  if(run1>=run2){
    printf("******   ERROR: invalid run range %d - %d\n",run1,run2);
    return;
  }

  TFile* fil=new TFile(ntupleFileName.Data(),"read");
  if(!fil){
    printf("File with ntuple does not exist\n");
    return;
  }
  TNtuple* ntsdd=(TNtuple*)fil->Get("ntsdd");  // aggiungere variabile da hGAmod per frazione di moduli ON

  Float_t nrun,nEvents, nEventsTriggered;
  Float_t meanTrPts3,errmeanTrPts3,meanTrPts4,errmeanTrPts4;
  Float_t minDrTime,errminDrTime,meanDrTime,errmeanDrTime;
  Float_t fracTrackWithClu1,fracTrackWithClu2,errfracTrackWithClu1,errfracTrackWithClu2;
  Float_t fracTrackWithClu3,fracTrackWithClu4,errfracTrackWithClu3,errfracTrackWithClu4;
  Float_t fracTrackWithClu5,fracTrackWithClu6,errfracTrackWithClu5,errfracTrackWithClu6;
  Float_t fracExtra,errfracExtra;
  Float_t meandEdxLay3,errmeandEdxLay3,meandEdxLay4,errmeandEdxLay4;
  Float_t meandEdxTB0,errmeandEdxTB0,meandEdxTB5,errmeandEdxTB5;
  Float_t MPVdEdxLay3,errMPVdEdxLay3,MPVdEdxLay4,errMPVdEdxLay4;
  Float_t MPVdEdxTB0,errMPVdEdxTB0,MPVdEdxTB5,errMPVdEdxTB5;
  Float_t nMod95,nMod80,nMod60,nModEmpty;
  Float_t fracEvWithSDD, errfracEvWithSDD;
    Float_t fracDead3, errfracDead3, fracDead4, errfracDead4;
  
  ntsdd->SetBranchAddress("nrun",&nrun);
  ntsdd->SetBranchAddress("nEvents",&nEvents);
  ntsdd->SetBranchAddress("nEventsTriggered",&nEventsTriggered);
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
    ntsdd->SetBranchAddress("fracEvWithSDD",&fracEvWithSDD);
    ntsdd->SetBranchAddress("errfracEvWithSDD",&errfracEvWithSDD);
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
  ntsdd->SetBranchAddress("MPVdEdxTB0",&MPVdEdxTB0);
  ntsdd->SetBranchAddress("errMPVdEdxTB0",&errMPVdEdxTB0);
  ntsdd->SetBranchAddress("MPVdEdxTB5",&MPVdEdxTB5);
  ntsdd->SetBranchAddress("errMPVdEdxTB5",&errMPVdEdxTB5);
  ntsdd->SetBranchAddress("MPVdEdxLay3",&MPVdEdxLay3);
  ntsdd->SetBranchAddress("errMPVdEdxLay3",&errMPVdEdxLay3);
  ntsdd->SetBranchAddress("MPVdEdxLay4",&MPVdEdxLay4);
  ntsdd->SetBranchAddress("errMPVdEdxLay4",&errMPVdEdxLay4);
    ntsdd->SetBranchAddress("fracDead3",&fracDead3);
    ntsdd->SetBranchAddress("errfracDead3",&errfracDead3);
    ntsdd->SetBranchAddress("fracDead4",&fracDead4);
    ntsdd->SetBranchAddress("errfracDead4",&errfracDead4);

  // Sort entries according to run number in the chosen range
  Int_t nr=ntsdd->GetEntries();
  Int_t *myIndex = new Int_t [nr];
  Int_t *noRuns = new Int_t [nr];
  for(Int_t i=0; i<nr;i++){
    ntsdd->GetEvent(i);
    Int_t intrun = static_cast<Int_t>(nrun+0.01);
    noRuns[i]=intrun;
  }
  printf("\n ======== PROCESSING SDD NTUPLE \n");
  Int_t kRunsToPlot = RunsToPlot(run1,run2,nr,noRuns,myIndex);

  TH1F* histotrp3=new TH1F("histotrp3","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* histonEvents=new TH1F("histonEvents","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* histonEventsTriggered=new TH1F("histoEventsTriggered","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* histotrp4=new TH1F("histotrp4","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* histominTime=new TH1F("histominTime","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* histomeanTime=new TH1F("histomeanTime","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* histofracExtra=new TH1F("histofracExtra","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* histodEdxTB0=new TH1F("histodEdxTB0","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* histodEdxTB5=new TH1F("histodEdxTB5","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* histodEdxLay3=new TH1F("histodEdxLay3","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* histodEdxLay4=new TH1F("histodEdxLay4","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* histoTrackClu1=new TH1F("histoTrackClu1","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* histoTrackClu2=new TH1F("histoTrackClu2","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* histoTrackClu3=new TH1F("histoTrackClu3","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* histoTrackClu4=new TH1F("histoTrackClu4","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* histoTrackClu5=new TH1F("histoTrackClu5","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* histoTrackClu6=new TH1F("histoTrackClu6","",kRunsToPlot,0.,kRunsToPlot);
    TH1F* histoEvwSDD=new TH1F("histoEvwSDD","",kRunsToPlot,0.,kRunsToPlot);
    TH1F* histoFracDead3=new TH1F("histoFracDead3","",kRunsToPlot,0.,kRunsToPlot);
    TH1F* histoFracDead4=new TH1F("histoFracDead4","",kRunsToPlot,0.,kRunsToPlot);

  TH1F* histoNmodEffBelow95=new TH1F("histoNmodEffBelow95","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* histoNmodEffBelow80=new TH1F("histoNmodEffBelow80","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* histoNmodEffBelow60=new TH1F("histoNmodEffBelow60","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* histoNmodEmpty=new TH1F("histoNmodEmpty","",kRunsToPlot,0.,kRunsToPlot);

  for(Int_t i=0; i<kRunsToPlot;i++){
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
//
      histoEvwSDD->SetBinContent(i+1,fracEvWithSDD);
      histoEvwSDD->SetBinError(i+1,errfracEvWithSDD);
      histoEvwSDD->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));

      histoFracDead3->SetBinContent(i+1,fracDead3);
      histoFracDead3->SetBinError(i+1,errfracDead3);
      histoFracDead3->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));

      histoFracDead4->SetBinContent(i+1,fracDead4);
      histoFracDead4->SetBinError(i+1,errfracDead4);
      histoFracDead4->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
//
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
    histotrp3->SetBinContent(i+1,meanTrPts3);
    histotrp3->SetBinError(i+1,errmeanTrPts3);
    histotrp3->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    histotrp4->SetBinContent(i+1,meanTrPts4);
    histotrp4->SetBinError(i+1,errmeanTrPts3);
    histotrp4->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
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
  Float_t  MPVdEdxLay5,errMPVdEdxLay5,MPVdEdxLay6,errMPVdEdxLay6;
  Float_t ChargeRatioL5,errChargeratioL5,ChargeRatioL6, errChargeratioL6, moduleOff;
    Float_t FracBadn5,errFracBadn5,FracBadp5,errFracBadp5,FracBadn6,errFracBadn6,FracBadp6,errFracBadp6;

  ntssd->SetBranchAddress("nrun",&nrunSSD);
  ntssd->SetBranchAddress("meandEdxLay5",&meandEdxLay5);
  ntssd->SetBranchAddress("errmeandEdxLay5",&errmeandEdxLay5);
  ntssd->SetBranchAddress("meandEdxLay6",&meandEdxLay6);
  ntssd->SetBranchAddress("errmeandEdxLay6",&errmeandEdxLay6);
  ntssd->SetBranchAddress("MPVdEdxLay5",&MPVdEdxLay5);
  ntssd->SetBranchAddress("errMPVdEdxLay5",&errMPVdEdxLay5);
  ntssd->SetBranchAddress("MPVdEdxLay6",&MPVdEdxLay6);
  ntssd->SetBranchAddress("errMPVdEdxLay6",&errMPVdEdxLay6);
  ntssd->SetBranchAddress("ChargeRatioL5",&ChargeRatioL5);
  ntssd->SetBranchAddress("errChargeratioL5",&errChargeratioL5);
  ntssd->SetBranchAddress("ChargeRatioL6",&ChargeRatioL6);
  ntssd->SetBranchAddress("errChargeratioL6",&errChargeratioL6);
  ntssd->SetBranchAddress("moduleOff",&moduleOff);
    ntssd->SetBranchAddress("FracBadn5",&FracBadn5);
    ntssd->SetBranchAddress("errFracBadn5",&errFracBadn5);
    ntssd->SetBranchAddress("FracBadp5",&FracBadp5);
    ntssd->SetBranchAddress("errFracBadp5",&errFracBadp5);
    ntssd->SetBranchAddress("FracBadn6",&FracBadn6);
    ntssd->SetBranchAddress("errFracBadn6",&errFracBadn6);
    ntssd->SetBranchAddress("FracBadp6",&FracBadp6);
    ntssd->SetBranchAddress("errFracBadp6",&errFracBadp6);

  nr=ntssd->GetEntries();
  delete []myIndex;
  delete []noRuns;
  myIndex = new Int_t [nr];
  noRuns = new Int_t [nr];
  for(Int_t i=0; i<nr;i++){
    ntssd->GetEvent(i);
    Int_t intrun = static_cast<Int_t>(nrunSSD+0.01);
    noRuns[i]=intrun;
  }
  printf("\n ======== PROCESSING SSD NTUPLE \n");
  kRunsToPlot = RunsToPlot(run1,run2,nr,noRuns,myIndex);

  TH1F* histodEdxLay5=new TH1F("histodEdxLay5","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* histodEdxLay6=new TH1F("histodEdxLay6","",kRunsToPlot,0.,kRunsToPlot);

  TH1F* histoChargeRatioLay5=new TH1F("histoChargeRatioLay5","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* histoChargeRatioLay6=new TH1F("histoChargeRatioLay6","",kRunsToPlot,0.,kRunsToPlot);
 
  TH1F* histoEmpty=new TH1F("histoEmpty","",kRunsToPlot,0.,kRunsToPlot);
  
  TH1F* histoFracBadn5=new TH1F("histoFracBadn5","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* histoFracBadp5=new TH1F("histoFracBadp5","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* histoFracBadn6=new TH1F("histoFracBadn6","",kRunsToPlot,0.,kRunsToPlot);
  TH1F* histoFracBadp6=new TH1F("histoFracBadp6","",kRunsToPlot,0.,kRunsToPlot);
    
    for(Int_t i=0; i<kRunsToPlot;i++){

    ntssd->GetEvent(myIndex[i]);

    histodEdxLay5->SetBinContent(i+1,MPVdEdxLay5);
    histodEdxLay5->SetBinError(i+1,errMPVdEdxLay5);
    histodEdxLay5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunSSD));

    histodEdxLay6->SetBinContent(i+1,MPVdEdxLay6);
    histodEdxLay6->SetBinError(i+1,errMPVdEdxLay6);
    histodEdxLay6->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunSSD));

    histoChargeRatioLay5->SetBinContent(i+1,ChargeRatioL5);
    histoChargeRatioLay5->SetBinError(i+1,errChargeratioL5);
    histoChargeRatioLay5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunSSD));
    
    histoChargeRatioLay6->SetBinContent(i+1,ChargeRatioL6);
    histoChargeRatioLay6->SetBinError(i+1,errChargeratioL6);
    histoChargeRatioLay6->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunSSD));

    histoEmpty->SetBinContent(i+1,moduleOff);
    histoEmpty->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunSSD));
        
    histoFracBadn5->SetBinContent(i+1,FracBadn5);
        cout << " run " << i << " FracBadn5 = " << FracBadn5 << endl;
    histoFracBadn5->SetBinError(i+1,errFracBadn5);
    histoFracBadn5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunSSD));

    histoFracBadp5->SetBinContent(i+1,FracBadp5);
        cout << " run " << i << " FracBadp5 = " << FracBadp5 << endl;
    histoFracBadp5->SetBinError(i+1,errFracBadp5);
    histoFracBadp5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunSSD));

    histoFracBadn6->SetBinContent(i+1,FracBadn6);
        cout << " run " << i << " FracBadn6 = " << FracBadn6 << endl;
    histoFracBadn6->SetBinError(i+1,errFracBadn6);
    histoFracBadn6->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunSSD));
        
    histoFracBadp6->SetBinContent(i+1,FracBadp6);
        cout << " run " << i << " FracBadp6 = " << FracBadp6 << endl;
    histoFracBadp6->SetBinError(i+1,errFracBadp6);
    histoFracBadp6->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunSSD));
        
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
    Float_t FracTrackMI1;
    Float_t errFracTrackMI1;
    Float_t FracTrackMI2;
    Float_t errFracTrackMI2;
    Float_t FracTrackMI3;
    Float_t errFracTrackMI3;
    Float_t FracTrackMI4;
    Float_t errFracTrackMI4;
    Float_t FracTrackMI5;
    Float_t errFracTrackMI5;
    Float_t FracTrackMI6;
    Float_t errFracTrackMI6;
    Float_t FracTrackSA1;
    Float_t errFracTrackSA1;
    Float_t FracTrackSA2;
    Float_t errFracTrackSA2;
    Float_t FracTrackSA3;
    Float_t errFracTrackSA3;
    Float_t FracTrackSA4;
    Float_t errFracTrackSA4;
    Float_t FracTrackSA5;
    Float_t errFracTrackSA5;
    Float_t FracTrackSA6;
    Float_t errFracTrackSA6;
    

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
    ntmatching->SetBranchAddress("FracTrackMI1",&FracTrackMI1);
    ntmatching->SetBranchAddress("errFracTrackMI1",&errFracTrackMI1);
    ntmatching->SetBranchAddress("FracTrackMI2",&FracTrackMI2);
    ntmatching->SetBranchAddress("errFracTrackMI2",&errFracTrackMI2);
    ntmatching->SetBranchAddress("FracTrackMI3",&FracTrackMI3);
    ntmatching->SetBranchAddress("errFracTrackMI3",&errFracTrackMI3);
    ntmatching->SetBranchAddress("FracTrackMI4",&FracTrackMI4);
    ntmatching->SetBranchAddress("errFracTrackMI4",&errFracTrackMI4);
    ntmatching->SetBranchAddress("FracTrackMI5",&FracTrackMI5);
    ntmatching->SetBranchAddress("errFracTrackMI5",&errFracTrackMI5);
    ntmatching->SetBranchAddress("FracTrackMI6",&FracTrackMI6);
    ntmatching->SetBranchAddress("errFracTrackMI6",&errFracTrackMI6);
    ntmatching->SetBranchAddress("FracTrackSA1",&FracTrackSA1);
    ntmatching->SetBranchAddress("errFracTrackSA1",&errFracTrackSA1);
    ntmatching->SetBranchAddress("FracTrackSA2",&FracTrackSA2);
    ntmatching->SetBranchAddress("errFracTrackSA2",&errFracTrackSA2);
    ntmatching->SetBranchAddress("FracTrackSA3",&FracTrackSA3);
    ntmatching->SetBranchAddress("errFracTrackSA3",&errFracTrackSA3);
    ntmatching->SetBranchAddress("FracTrackSA4",&FracTrackSA4);
    ntmatching->SetBranchAddress("errFracTrackSA4",&errFracTrackSA4);
    ntmatching->SetBranchAddress("FracTrackSA5",&FracTrackSA5);
    ntmatching->SetBranchAddress("errFracTrackSA5",&errFracTrackSA5);
    ntmatching->SetBranchAddress("FracTrackSA6",&FracTrackSA6);
    ntmatching->SetBranchAddress("errFracTrackSA6",&errFracTrackSA6);

  nr=ntmatching->GetEntries();
  delete []myIndex;
  delete []noRuns;
  myIndex = new Int_t [nr];
  noRuns = new Int_t [nr];
  for(Int_t i=0; i<nr;i++){
    ntmatching->GetEvent(i);
    Int_t intrun = static_cast<Int_t>(nrunMatch+0.01);
    noRuns[i]=intrun;
  }
  printf("\n ======== PROCESSING NTMATCHING NTUPLE \n");
  kRunsToPlot = RunsToPlot(run1,run2,nr,noRuns,myIndex);

  TH1F *hFracSPD1 = new TH1F("hFracSPD1","SPD inner; run number; Fraction of HSs",kRunsToPlot,0.,kRunsToPlot);
  hFracSPD1->SetLineColor(kGreen+2);
  hFracSPD1->SetMarkerColor(kGreen+2);
  hFracSPD1->SetMarkerStyle(20);

  TH1F *hFracSPD2 = new TH1F("hFracSPD2","SPD outer; run number; Fraction of HSs",kRunsToPlot,0.,kRunsToPlot);
  hFracSPD2->SetLineColor(kYellow+2);
  hFracSPD2->SetMarkerColor(kYellow+2);
  hFracSPD2->SetMarkerStyle(20);

  TH1F *hEffSPDPt02 = new TH1F("hEffSPDPt02","Efficiency - P_{T} = 0.2; run number; TPC+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEffSPDPt02->SetLineWidth(2);
  hEffSPDPt02->SetLineColor(kAzure+1);
  hEffSPDPt02->SetMarkerColor(kAzure+1);
  hEffSPDPt02->SetMarkerStyle(20);

  TH1F *hEffSPDPt1 = new TH1F("hEffSPDPt1","Efficiency; run number; TPC+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEffSPDPt1->SetLineWidth(2);
  hEffSPDPt1->SetLineColor(kAzure+1);
  hEffSPDPt1->SetMarkerColor(kAzure+1);
  hEffSPDPt1->SetMarkerStyle(20);

  TH1F *hEffSPDPt10 = new TH1F("hEffSPDPt10","Efficiency; run number; TPC+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEffSPDPt10->SetLineWidth(2);
  hEffSPDPt10->SetLineColor(kAzure+1);
  hEffSPDPt10->SetMarkerColor(kAzure+1);
  hEffSPDPt10->SetMarkerStyle(20);

  TH1F *hEffoneSPDPt02 = new TH1F("hEffoneSPDPt02","Efficiency; run number; TPC+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEffoneSPDPt02->SetLineWidth(2);
  hEffoneSPDPt02->SetLineColor(kGray);
  hEffoneSPDPt02->SetMarkerColor(kGray);
  hEffoneSPDPt02->SetMarkerStyle(20);
  TH1F *hEffoneSPDPt1 = new TH1F("hEffoneSPDPt1","Efficiency; run number; TPC+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEffoneSPDPt1->SetLineWidth(2);
  hEffoneSPDPt1->SetLineColor(kGray);
  hEffoneSPDPt1->SetMarkerColor(kGray);
  hEffoneSPDPt1->SetMarkerStyle(20);
  TH1F *hEffoneSPDPt10 = new TH1F("hEffoneSPDPt10","Efficiency; run number; TPC+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEffoneSPDPt10->SetLineWidth(2);
  hEffoneSPDPt10->SetLineColor(kGray);
  hEffoneSPDPt10->SetMarkerColor(kGray);
  hEffoneSPDPt10->SetMarkerStyle(20);

  TH1F *hEff2Pt02 = new TH1F("hEff2Pt02","Efficiency; run number; TPC+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEff2Pt02->SetLineWidth(2);
  hEff2Pt02->SetLineColor(kViolet);
  hEff2Pt02->SetMarkerColor(kViolet);
  hEff2Pt02->SetMarkerStyle(20);
  TH1F *hEff2Pt1 = new TH1F("hEff2Pt1","Efficiency; run number; TPC+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEff2Pt1->SetLineWidth(2);
  hEff2Pt1->SetLineColor(kViolet);
  hEff2Pt1->SetMarkerColor(kViolet);
  hEff2Pt1->SetMarkerStyle(20);
  TH1F *hEff2Pt10 = new TH1F("hEff2Pt10","Efficiency; run number; TPC+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEff2Pt10->SetLineWidth(2);
  hEff2Pt10->SetLineColor(kViolet);
  hEff2Pt10->SetMarkerColor(kViolet);
  hEff2Pt10->SetMarkerStyle(20);

  TH1F *hEff3Pt02 = new TH1F("hEff3Pt02","Efficiency; run number; TPC+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEff3Pt02->SetLineWidth(2);
  hEff3Pt02->SetLineColor(6);
  hEff3Pt02->SetMarkerColor(6);
  hEff3Pt02->SetMarkerStyle(20);
  TH1F *hEff3Pt1 = new TH1F("hEff3Pt1","Efficiency; run number; TPC+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEff3Pt1->SetLineWidth(2);
  hEff3Pt1->SetLineColor(6);
  hEff3Pt1->SetMarkerColor(6);
  hEff3Pt1->SetMarkerStyle(20);
  TH1F *hEff3Pt10 = new TH1F("hEff3Pt10","Efficiency; run number; TPC+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEff3Pt10->SetLineWidth(2);
  hEff3Pt10->SetLineColor(6);
  hEff3Pt10->SetMarkerColor(6);
  hEff3Pt10->SetMarkerStyle(20);

  TH1F *hEff4Pt02 = new TH1F("hEff4Pt02","Efficiency; run number; TPC+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEff4Pt02->SetLineWidth(2);
  hEff4Pt02->SetLineColor(4);
  hEff4Pt02->SetMarkerColor(4);
  hEff4Pt02->SetMarkerStyle(20);
  TH1F *hEff4Pt1 = new TH1F("hEff4Pt1","Efficiency; run number; TPC+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEff4Pt1->SetLineWidth(2);
  hEff4Pt1->SetLineColor(4);
  hEff4Pt1->SetMarkerColor(4);
  hEff4Pt1->SetMarkerStyle(20);
  TH1F *hEff4Pt10 = new TH1F("hEff4Pt10","Efficiency; run number; TPC+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEff4Pt10->SetLineWidth(2);
  hEff4Pt10->SetLineColor(4);
  hEff4Pt10->SetMarkerColor(4);
  hEff4Pt10->SetMarkerStyle(20);

  TH1F *hEff5Pt02 = new TH1F("hEff5Pt02","Efficiency; run number; TPC+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEff5Pt02->SetLineWidth(2);
  hEff5Pt02->SetLineColor(3);
  hEff5Pt02->SetMarkerColor(3);
  hEff5Pt02->SetMarkerStyle(20);
  TH1F *hEff5Pt1 = new TH1F("hEff5Pt1","Efficiency; run number; TPC+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEff5Pt1->SetLineWidth(2);
  hEff5Pt1->SetLineColor(3);
  hEff5Pt1->SetMarkerColor(3);
  hEff5Pt1->SetMarkerStyle(20);
  TH1F *hEff5Pt10 = new TH1F("hEff5Pt10","Efficiency; run number; TPC+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEff5Pt10->SetLineWidth(3);
  hEff5Pt10->SetLineColor(3);
  hEff5Pt10->SetMarkerColor(3);
  hEff5Pt10->SetMarkerStyle(20);

  TH1F *hEff6Pt02 = new TH1F("hEff6Pt02","Efficiency; run number; TPC+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEff6Pt02->SetLineWidth(2);
  hEff6Pt02->SetLineColor(2);
  hEff6Pt02->SetMarkerColor(2);
  hEff6Pt02->SetMarkerStyle(20);
  TH1F *hEff6Pt1 = new TH1F("hEff6Pt1","Efficiency; run number; TPC+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEff6Pt1->SetLineWidth(2);
  hEff6Pt1->SetLineColor(2);
  hEff6Pt1->SetMarkerColor(2);
  hEff6Pt1->SetMarkerStyle(20);
  TH1F *hEff6Pt10 = new TH1F("hEff6Pt10","Efficiency; run number; TPC+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEff6Pt10->SetLineWidth(2);
  hEff6Pt10->SetLineColor(2);
  hEff6Pt10->SetMarkerColor(2);
  hEff6Pt10->SetMarkerStyle(20);


  TH1F *hEffTOTPt02 = new TH1F("hEffTOTPt02","Efficiency; run number; TPC+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEffTOTPt02->SetLineWidth(2);
  hEffTOTPt02->SetLineColor(kBlue+2);
  hEffTOTPt02->SetMarkerColor(kBlue+2);
  hEffTOTPt02->SetMarkerStyle(20);
  TH1F *hEffTOTPt1 = new TH1F("hEffTOTPt1","Efficiency; run number; TPC+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEffTOTPt1->SetLineWidth(2);
  hEffTOTPt1->SetLineColor(kBlue+2);
  hEffTOTPt1->SetMarkerColor(kBlue+2);
  hEffTOTPt1->SetMarkerStyle(20);
  TH1F *hEffTOTPt10 = new TH1F("hEffTOTPt10","Efficiency; run number; TPC+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEffTOTPt10->SetLineWidth(2);
  hEffTOTPt10->SetLineColor(kBlue+2);
  hEffTOTPt10->SetMarkerColor(kBlue+2);
  hEffTOTPt10->SetMarkerStyle(20);
    
    TH1F* histoTrackMI1=new TH1F("histoTrackMI1","",kRunsToPlot,0.,kRunsToPlot);
    TH1F* histoTrackMI2=new TH1F("histoTrackMI2","",kRunsToPlot,0.,kRunsToPlot);
    TH1F* histoTrackMI3=new TH1F("histoTrackMI3","",kRunsToPlot,0.,kRunsToPlot);
    TH1F* histoTrackMI4=new TH1F("histoTrackMI4","",kRunsToPlot,0.,kRunsToPlot);
    TH1F* histoTrackMI5=new TH1F("histoTrackMI5","",kRunsToPlot,0.,kRunsToPlot);
    TH1F* histoTrackMI6=new TH1F("histoTrackMI6","",kRunsToPlot,0.,kRunsToPlot);
    TH1F* histoTrackSA1=new TH1F("histoTrackSA1","",kRunsToPlot,0.,kRunsToPlot);
    TH1F* histoTrackSA2=new TH1F("histoTrackSA2","",kRunsToPlot,0.,kRunsToPlot);
    TH1F* histoTrackSA3=new TH1F("histoTrackSA3","",kRunsToPlot,0.,kRunsToPlot);
    TH1F* histoTrackSA4=new TH1F("histoTrackSA4","",kRunsToPlot,0.,kRunsToPlot);
    TH1F* histoTrackSA5=new TH1F("histoTrackSA5","",kRunsToPlot,0.,kRunsToPlot);
    TH1F* histoTrackSA6=new TH1F("histoTrackSA6","",kRunsToPlot,0.,kRunsToPlot);
    

  //  Int_t nEntriesMatch=ntmatching->GetEntries();
  
  for(Int_t i=0;i<kRunsToPlot;i++){
   
    ntmatching->GetEvent(myIndex[i]);
    //    Int_t bin=nrunMatch;

    // fill histos
    //    cout<<i<<") "<<"Index= "<<myIndex[i]<<" nrun= "<<nrunMatch<<", FracSPD1= "<<FracSPD1<<", FracSPD2"<<FracSPD2<<endl;        
    hFracSPD1->SetBinContent(i+1,FracSPD1);
    hFracSPD1->SetBinError(i+1,.01);
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

      histoTrackMI1->SetBinContent(i+1,FracTrackMI1);
      histoTrackMI1->SetBinError(i+1,errFracTrackMI1);
      histoTrackMI1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));

      histoTrackMI2->SetBinContent(i+1,FracTrackMI2);
      histoTrackMI2->SetBinError(i+1,errFracTrackMI2);
      histoTrackMI2->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));

      histoTrackMI3->SetBinContent(i+1,FracTrackMI3);
      histoTrackMI3->SetBinError(i+1,errFracTrackMI3);
      histoTrackMI3->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));

      histoTrackMI4->SetBinContent(i+1,FracTrackMI4);
      histoTrackMI4->SetBinError(i+1,errFracTrackMI4);
      histoTrackMI4->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));
  
      histoTrackMI5->SetBinContent(i+1,FracTrackMI5);
      histoTrackMI5->SetBinError(i+1,errFracTrackMI5);
      histoTrackMI5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));

      histoTrackMI6->SetBinContent(i+1,FracTrackMI6);
      histoTrackMI6->SetBinError(i+1,errFracTrackMI6);
      histoTrackMI6->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));

      histoTrackSA1->SetBinContent(i+1,FracTrackSA1);
      histoTrackSA1->SetBinError(i+1,errFracTrackSA1);
      histoTrackSA1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));
      
      histoTrackSA2->SetBinContent(i+1,FracTrackSA2);
      histoTrackSA2->SetBinError(i+1,errFracTrackSA2);
      histoTrackSA2->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));
      
      histoTrackSA3->SetBinContent(i+1,FracTrackSA3);
      histoTrackSA3->SetBinError(i+1,errFracTrackSA3);
      histoTrackSA3->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));
      
      histoTrackSA4->SetBinContent(i+1,FracTrackSA4);
      histoTrackSA4->SetBinError(i+1,errFracTrackSA4);
      histoTrackSA4->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));
      
      histoTrackSA5->SetBinContent(i+1,FracTrackSA5);
      histoTrackSA5->SetBinError(i+1,errFracTrackSA5);
      histoTrackSA5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));
      
      histoTrackSA6->SetBinContent(i+1,FracTrackSA6);
      histoTrackSA6->SetBinError(i+1,errFracTrackSA6);
      histoTrackSA6->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));

  }
  
  // Matching TOF 

  TNtuple* ntmatchingTOF=(TNtuple*)fil->Get("ntmatchingTOF");

  if(!ntmatchingTOF){
    cout<<"Problem with TOF ntple"<<endl;
    return;
  }

  Float_t nrunMatchTOF;
  Float_t FracSPD1TOF;
  Float_t errFracSPD1TOF;
  Float_t FracSPD2TOF;
  Float_t errFracSPD2TOF;
  Float_t Eff6Pt02TOF;
  Float_t errEff6Pt02TOF;
  Float_t Eff6Pt1TOF;
  Float_t errEff6Pt1TOF;
  Float_t Eff6Pt10TOF;
  Float_t errEff6Pt10TOF;
  Float_t Eff5Pt02TOF;
  Float_t errEff5Pt02TOF;
  Float_t Eff5Pt1TOF;
  Float_t errEff5Pt1TOF;
  Float_t Eff5Pt10TOF;
  Float_t errEff5Pt10TOF;
  Float_t Eff4Pt02TOF;
  Float_t errEff4Pt02TOF;
  Float_t Eff4Pt1TOF;
  Float_t errEff4Pt1TOF;
  Float_t Eff4Pt10TOF;
  Float_t errEff4Pt10TOF;
  Float_t Eff3Pt02TOF;
  Float_t errEff3Pt02TOF;
  Float_t Eff3Pt1TOF;
  Float_t errEff3Pt1TOF;
  Float_t Eff3Pt10TOF;
  Float_t errEff3Pt10TOF;
  Float_t Eff2Pt02TOF;
  Float_t errEff2Pt02TOF;
  Float_t Eff2Pt1TOF;
  Float_t errEff2Pt1TOF;
  Float_t Eff2Pt10TOF;
  Float_t errEff2Pt10TOF;
  Float_t EffSPDPt02TOF;
  Float_t errEffSPDPt02TOF;
  Float_t EffSPDPt1TOF;
  Float_t errEffSPDPt1TOF;
  Float_t EffSPDPt10TOF;
  Float_t errEffSPDPt10TOF;
  Float_t EffoneSPDPt02TOF;
  Float_t errEffoneSPDPt02TOF;
  Float_t EffoneSPDPt1TOF;
  Float_t errEffoneSPDPt1TOF;
  Float_t EffoneSPDPt10TOF;
  Float_t errEffoneSPDPt10TOF;
  Float_t EffTOTPt02TOF;
  Float_t errEffTOTPt02TOF;
  Float_t EffTOTPt1TOF;
  Float_t errEffTOTPt1TOF;
  Float_t EffTOTPt10TOF;
  Float_t errEffTOTPt10TOF;

  ntmatchingTOF->SetBranchAddress("nrun",            &nrunMatchTOF);
  ntmatchingTOF->SetBranchAddress("FracSPD1",        &FracSPD1TOF);
  ntmatchingTOF->SetBranchAddress("errFracSPD1",     &errFracSPD1TOF);
  ntmatchingTOF->SetBranchAddress("FracSPD2",        &FracSPD2TOF);
  ntmatchingTOF->SetBranchAddress("errFracSPD2",     &errFracSPD2TOF);
  ntmatchingTOF->SetBranchAddress("Eff6Pt02",        &Eff6Pt02TOF);
  ntmatchingTOF->SetBranchAddress("errEff6Pt02",     &errEff6Pt02TOF);
  ntmatchingTOF->SetBranchAddress("Eff6Pt1",         &Eff6Pt1TOF);
  ntmatchingTOF->SetBranchAddress("errEff6Pt1",      &errEff6Pt1TOF);
  ntmatchingTOF->SetBranchAddress("Eff6Pt10",        &Eff6Pt10TOF);
  ntmatchingTOF->SetBranchAddress("errEff6Pt10",     &errEff6Pt10TOF);
  ntmatchingTOF->SetBranchAddress("Eff5Pt02",        &Eff5Pt02TOF);
  ntmatchingTOF->SetBranchAddress("errEff5Pt02",     &errEff5Pt02TOF);
  ntmatchingTOF->SetBranchAddress("Eff5Pt1",         &Eff5Pt1TOF);
  ntmatchingTOF->SetBranchAddress("errEff5Pt1",      &errEff5Pt1TOF);
  ntmatchingTOF->SetBranchAddress("Eff5Pt10",        &Eff5Pt10TOF);
  ntmatchingTOF->SetBranchAddress("errEff5Pt10",     &errEff5Pt10TOF);
  ntmatchingTOF->SetBranchAddress("Eff4Pt02",        &Eff4Pt02TOF);
  ntmatchingTOF->SetBranchAddress("errEff4Pt02",     &errEff4Pt02TOF);
  ntmatchingTOF->SetBranchAddress("Eff4Pt1",         &Eff4Pt1TOF);
  ntmatchingTOF->SetBranchAddress("errEff4Pt1",      &errEff4Pt1TOF);
  ntmatchingTOF->SetBranchAddress("Eff4Pt10",        &Eff4Pt10TOF);
  ntmatchingTOF->SetBranchAddress("errEff4Pt10",     &errEff4Pt10TOF);
  ntmatchingTOF->SetBranchAddress("Eff3Pt02",        &Eff3Pt02TOF);
  ntmatchingTOF->SetBranchAddress("errEff3Pt02",     &errEff3Pt02TOF);
  ntmatchingTOF->SetBranchAddress("Eff3Pt1",         &Eff3Pt1TOF);
  ntmatchingTOF->SetBranchAddress("errEff3Pt1",      &errEff3Pt1TOF);
  ntmatchingTOF->SetBranchAddress("Eff3Pt10",        &Eff3Pt10TOF);
  ntmatchingTOF->SetBranchAddress("errEff3Pt10",     &errEff3Pt10TOF);
  ntmatchingTOF->SetBranchAddress("Eff2Pt02",        &Eff2Pt02TOF);
  ntmatchingTOF->SetBranchAddress("errEff2Pt02",     &errEff2Pt02TOF);
  ntmatchingTOF->SetBranchAddress("Eff2Pt1",         &Eff2Pt1TOF);
  ntmatchingTOF->SetBranchAddress("errEff2Pt1",      &errEff2Pt1TOF);
  ntmatchingTOF->SetBranchAddress("Eff2Pt10",        &Eff2Pt10TOF);
  ntmatchingTOF->SetBranchAddress("errEff2Pt10",     &errEff2Pt10TOF);
  ntmatchingTOF->SetBranchAddress("EffSPDPt02",      &EffSPDPt02TOF);
  ntmatchingTOF->SetBranchAddress("errEffSPDPt02",   &errEffSPDPt02TOF);  
  ntmatchingTOF->SetBranchAddress("EffSPDPt1",       &EffSPDPt1TOF);
  ntmatchingTOF->SetBranchAddress("errEffSPDPt1",    &errEffSPDPt1TOF);  
  ntmatchingTOF->SetBranchAddress("EffSPDPt10",      &EffSPDPt10TOF);
  ntmatchingTOF->SetBranchAddress("errEffSPDPt10",   &errEffSPDPt10TOF);  
  ntmatchingTOF->SetBranchAddress("EffoneSPDPt02",   &EffoneSPDPt02TOF);
  ntmatchingTOF->SetBranchAddress("errEffoneSPDPt02",&errEffoneSPDPt02TOF);
  ntmatchingTOF->SetBranchAddress("EffoneSPDPt1",    &EffoneSPDPt1TOF);
  ntmatchingTOF->SetBranchAddress("errEffoneSPDPt1", &errEffoneSPDPt1TOF);
  ntmatchingTOF->SetBranchAddress("EffoneSPDPt10",   &EffoneSPDPt10TOF);
  ntmatchingTOF->SetBranchAddress("errEffoneSPDPt10",&errEffoneSPDPt10TOF);
  ntmatchingTOF->SetBranchAddress("EffTOTPt02",      &EffTOTPt02TOF);
  ntmatchingTOF->SetBranchAddress("errEffTOTPt02",   &errEffTOTPt02TOF);
  ntmatchingTOF->SetBranchAddress("EffTOTPt1",       &EffTOTPt1TOF);
  ntmatchingTOF->SetBranchAddress("errEffTOTPt1",    &errEffTOTPt1TOF);
  ntmatchingTOF->SetBranchAddress("EffTOTPt10",      &EffTOTPt10TOF);
  ntmatchingTOF->SetBranchAddress("errEffTOTPt10",   &errEffTOTPt10TOF);

  // //------------------------------------
  // nr=ntmatching->GetEntries();
  // delete []myIndex;
  // delete []noRuns;
  // myIndex = new Int_t [nr];
  // noRuns = new Int_t [nr];
  // for(Int_t i=0; i<nr;i++){
  //   ntmatching->GetEvent(i);
  //   Int_t intrun = static_cast<Int_t>(nrunMatch+0.01);
  //   noRuns[i]=intrun;
  // }
  // printf("\n ======== PROCESSING NTMATCHING NTUPLE \n");
  // kRunsToPlot = RunsToPlot(run1,run2,nr,noRuns,myIndex);

  // //--------------------------------------

  nr=ntmatchingTOF->GetEntries();
  delete []myIndex;
  delete []noRuns;
  myIndex = new Int_t [nr];
  noRuns = new Int_t [nr];
  for(Int_t i=0; i<nr;i++){
    ntmatchingTOF->GetEvent(i);
    Int_t intrun = static_cast<Int_t>(nrunMatchTOF+0.01);
    noRuns[i]=intrun;
  }
  printf("\n ======== PROCESSING NTMATCHING NTUPLE \n");
  kRunsToPlot = RunsToPlot(run1,run2,nr,noRuns,myIndex);

  TH1F *hEffTOFSPDPt02 = new TH1F("hEffTOFSPDPt02","Efficiency - P_{T} = 0.5; run number; TPC+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEffTOFSPDPt02->SetLineWidth(2);
  hEffTOFSPDPt02->SetLineColor(kAzure+1);
  hEffTOFSPDPt02->SetMarkerColor(kAzure+1);
  hEffTOFSPDPt02->SetMarkerStyle(20);

  TH1F *hEffTOFSPDPt1 = new TH1F("hEffTOFSPDPt1","Efficiency; run number; TOF+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEffTOFSPDPt1->SetLineWidth(2);
  hEffTOFSPDPt1->SetLineColor(kAzure+1);
  hEffTOFSPDPt1->SetMarkerColor(kAzure+1);
  hEffTOFSPDPt1->SetMarkerStyle(20);

  TH1F *hEffTOFSPDPt10 = new TH1F("hEffTOFSPDPt10","Efficiency; run number; TOF+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEffTOFSPDPt10->SetLineWidth(2);
  hEffTOFSPDPt10->SetLineColor(kAzure+1);
  hEffTOFSPDPt10->SetMarkerColor(kAzure+1);
  hEffTOFSPDPt10->SetMarkerStyle(20);

  TH1F *hEffTOFoneSPDPt02 = new TH1F("hEffTOFoneSPDPt02","Efficiency; run number; TOF+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEffTOFoneSPDPt02->SetLineWidth(2);
  hEffTOFoneSPDPt02->SetLineColor(kGray);
  hEffTOFoneSPDPt02->SetMarkerColor(kGray);
  hEffTOFoneSPDPt02->SetMarkerStyle(20);

  TH1F *hEffTOFoneSPDPt1 = new TH1F("hEffTOFoneSPDPt1","Efficiency; run number; TOF+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEffTOFoneSPDPt1->SetLineWidth(2);
  hEffTOFoneSPDPt1->SetLineColor(kGray);
  hEffTOFoneSPDPt1->SetMarkerColor(kGray);
  hEffTOFoneSPDPt1->SetMarkerStyle(20);

  TH1F *hEffTOFoneSPDPt10 = new TH1F("hEffTOFoneSPDPt10","Efficiency; run number; TOF+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEffTOFoneSPDPt10->SetLineWidth(2);
  hEffTOFoneSPDPt10->SetLineColor(kGray);
  hEffTOFoneSPDPt10->SetMarkerColor(kGray);
  hEffTOFoneSPDPt10->SetMarkerStyle(20);

  TH1F *hEffTOF2Pt02 = new TH1F("hEffTOF2Pt02","Efficiency; run number; TOF+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEffTOF2Pt02->SetLineWidth(2);
  hEffTOF2Pt02->SetLineColor(kViolet);
  hEffTOF2Pt02->SetMarkerColor(kViolet);
  hEffTOF2Pt02->SetMarkerStyle(20);

  TH1F *hEffTOF2Pt1 = new TH1F("hEffTOF2Pt1","Efficiency; run number; TOF+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEffTOF2Pt1->SetLineWidth(2);
  hEffTOF2Pt1->SetLineColor(kViolet);
  hEffTOF2Pt1->SetMarkerColor(kViolet);
  hEffTOF2Pt1->SetMarkerStyle(20);

  TH1F *hEffTOF2Pt10 = new TH1F("hEffTOF2Pt10","Efficiency; run number; TOF+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEffTOF2Pt10->SetLineWidth(2);
  hEffTOF2Pt10->SetLineColor(kViolet);
  hEffTOF2Pt10->SetMarkerColor(kViolet);
  hEffTOF2Pt10->SetMarkerStyle(20);

  TH1F *hEffTOF3Pt02 = new TH1F("hEffTOF3Pt02","Efficiency; run number; TOF+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEffTOF3Pt02->SetLineWidth(2);
  hEffTOF3Pt02->SetLineColor(6);
  hEffTOF3Pt02->SetMarkerColor(6);
  hEffTOF3Pt02->SetMarkerStyle(20);

  TH1F *hEffTOF3Pt1 = new TH1F("hEffTOF3Pt1","Efficiency; run number; TOF+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEffTOF3Pt1->SetLineWidth(2);
  hEffTOF3Pt1->SetLineColor(6);
  hEffTOF3Pt1->SetMarkerColor(6);
  hEffTOF3Pt1->SetMarkerStyle(20);

  TH1F *hEffTOF3Pt10 = new TH1F("hEffTOF3Pt10","Efficiency; run number; TOF+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEffTOF3Pt10->SetLineWidth(2);
  hEffTOF3Pt10->SetLineColor(6);
  hEffTOF3Pt10->SetMarkerColor(6);
  hEffTOF3Pt10->SetMarkerStyle(20);

  TH1F *hEffTOF4Pt02 = new TH1F("hEffTOF4Pt02","Efficiency; run number; TOF+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEffTOF4Pt02->SetLineWidth(2);
  hEffTOF4Pt02->SetLineColor(4);
  hEffTOF4Pt02->SetMarkerColor(4);
  hEffTOF4Pt02->SetMarkerStyle(20);

  TH1F *hEffTOF4Pt1 = new TH1F("hEffTOF4Pt1","Efficiency; run number; TOF+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEffTOF4Pt1->SetLineWidth(2);
  hEffTOF4Pt1->SetLineColor(4);
  hEffTOF4Pt1->SetMarkerColor(4);
  hEffTOF4Pt1->SetMarkerStyle(20);

  TH1F *hEffTOF4Pt10 = new TH1F("hEffTOF4Pt10","Efficiency; run number; TOF+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEffTOF4Pt10->SetLineWidth(2);
  hEffTOF4Pt10->SetLineColor(4);
  hEffTOF4Pt10->SetMarkerColor(4);
  hEffTOF4Pt10->SetMarkerStyle(20);

  TH1F *hEffTOF5Pt02 = new TH1F("hEffTOF5Pt02","Efficiency; run number; TOF+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEffTOF5Pt02->SetLineWidth(2);
  hEffTOF5Pt02->SetLineColor(3);
  hEffTOF5Pt02->SetMarkerColor(3);
  hEffTOF5Pt02->SetMarkerStyle(20);

  TH1F *hEffTOF5Pt1 = new TH1F("hEffTOF5Pt1","Efficiency; run number; TOF+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEffTOF5Pt1->SetLineWidth(2);
  hEffTOF5Pt1->SetLineColor(3);
  hEffTOF5Pt1->SetMarkerColor(3);
  hEffTOF5Pt1->SetMarkerStyle(20);

  TH1F *hEffTOF5Pt10 = new TH1F("hEffTOF5Pt10","Efficiency; run number; TOF+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEffTOF5Pt10->SetLineWidth(3);
  hEffTOF5Pt10->SetLineColor(3);
  hEffTOF5Pt10->SetMarkerColor(3);
  hEffTOF5Pt10->SetMarkerStyle(20);

  TH1F *hEffTOF6Pt02 = new TH1F("hEffTOF6Pt02","Efficiency; run number; TOF+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEffTOF6Pt02->SetLineWidth(2);
  hEffTOF6Pt02->SetLineColor(2);
  hEffTOF6Pt02->SetMarkerColor(2);
  hEffTOF6Pt02->SetMarkerStyle(20);

  TH1F *hEffTOF6Pt1 = new TH1F("hEffTOF6Pt1","Efficiency; run number; TOF+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEffTOF6Pt1->SetLineWidth(2);
  hEffTOF6Pt1->SetLineColor(2);
  hEffTOF6Pt1->SetMarkerColor(2);
  hEffTOF6Pt1->SetMarkerStyle(20);

  TH1F *hEffTOF6Pt10 = new TH1F("hEffTOF6Pt10","Efficiency; run number; TOF+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEffTOF6Pt10->SetLineWidth(2);
  hEffTOF6Pt10->SetLineColor(2);
  hEffTOF6Pt10->SetMarkerColor(2);
  hEffTOF6Pt10->SetMarkerStyle(20);

  TH1F *hEffTOFTOTPt02 = new TH1F("hEffTOFTOTPt02","Efficiency; run number; TOF+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEffTOFTOTPt02->SetLineWidth(2);
  hEffTOFTOTPt02->SetLineColor(kBlue+2);
  hEffTOFTOTPt02->SetMarkerColor(kBlue+2);
  hEffTOFTOTPt02->SetMarkerStyle(20);

  TH1F *hEffTOFTOTPt1 = new TH1F("hEffTOFTOTPt1","Efficiency; run number; TOF+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEffTOFTOTPt1->SetLineWidth(2);
  hEffTOFTOTPt1->SetLineColor(kBlue+2);
  hEffTOFTOTPt1->SetMarkerColor(kBlue+2);
  hEffTOFTOTPt1->SetMarkerStyle(20);

  TH1F *hEffTOFTOTPt10 = new TH1F("hEffTOFTOTPt10","Efficiency; run number; TOF+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
  hEffTOFTOTPt10->SetLineWidth(2);
  hEffTOFTOTPt10->SetLineColor(kBlue+2);
  hEffTOFTOTPt10->SetMarkerColor(kBlue+2);
  hEffTOFTOTPt10->SetMarkerStyle(20);

  //  Int_t nEntriesMatch=ntmatching->GetEntries();

  for(Int_t i=0;i<kRunsToPlot;i++){
   
    ntmatchingTOF->GetEvent(myIndex[i]);
    //    Int_t bin=nrunMatch;

    // fill histos

    //-------------------------

    hEffTOF6Pt02->SetBinContent(i+1,Eff6Pt02TOF);
    hEffTOF6Pt02->SetBinError(i+1,errEff6Pt02TOF);
    hEffTOF6Pt02->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatchTOF));

    hEffTOF6Pt1->SetBinContent(i+1,Eff6Pt1TOF);
    hEffTOF6Pt1->SetBinError(i+1,errEff6Pt1TOF);
    hEffTOF6Pt1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatchTOF));

    hEffTOF6Pt10->SetBinContent(i+1,Eff6Pt10TOF);
    hEffTOF6Pt10->SetBinError(i+1,errEff6Pt10TOF);
    hEffTOF6Pt10->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatchTOF));

    hEffTOF5Pt02->SetBinContent(i+1,Eff5Pt02TOF);
    hEffTOF5Pt02->SetBinError(i+1,errEff5Pt02TOF);
    hEffTOF5Pt02->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatchTOF));
	   
    hEffTOF5Pt1->SetBinContent(i+1,Eff5Pt1TOF);
    hEffTOF5Pt1->SetBinError(i+1,errEff5Pt1TOF);
    hEffTOF5Pt1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatchTOF));
	   
    hEffTOF5Pt10->SetBinContent(i+1,Eff5Pt10TOF);
    hEffTOF5Pt10->SetBinError(i+1,errEff5Pt10TOF);
    hEffTOF5Pt10->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatchTOF));
      
    hEffTOF4Pt02->SetBinContent(i+1,Eff4Pt02TOF);
    hEffTOF4Pt02->SetBinError(i+1,errEff4Pt1TOF);
    hEffTOF4Pt02->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatchTOF));

    hEffTOF4Pt1->SetBinContent(i+1,Eff4Pt1TOF);
    hEffTOF4Pt1->SetBinError(i+1,errEff4Pt1TOF);
    hEffTOF4Pt1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatchTOF));

    hEffTOF4Pt10->SetBinContent(i+1,Eff4Pt10TOF);
    hEffTOF4Pt10->SetBinError(i+1,errEff4Pt10TOF);
    hEffTOF4Pt10->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatchTOF));

    hEffTOF3Pt02->SetBinContent(i+1,Eff3Pt02TOF);
    hEffTOF3Pt02->SetBinError(i+1,errEff3Pt02TOF);
    hEffTOF3Pt02->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatchTOF));
	    
    hEffTOF3Pt1->SetBinContent(i+1,Eff3Pt1TOF);
    hEffTOF3Pt1->SetBinError(i+1,errEff3Pt1TOF);
    hEffTOF3Pt1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatchTOF));

    hEffTOF3Pt10->SetBinContent(i+1,Eff3Pt10TOF);
    hEffTOF3Pt10->SetBinError(i+1,errEff3Pt10TOF);
    hEffTOF3Pt10->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatchTOF));      
      
    hEffTOF2Pt02->SetBinContent(i+1,Eff2Pt02TOF);
    hEffTOF2Pt02->SetBinError(i+1,errEff2Pt02TOF);
    hEffTOF2Pt02->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatchTOF));  

    hEffTOF2Pt1->SetBinContent(i+1,Eff2Pt1TOF);
    hEffTOF2Pt1->SetBinError(i+1,errEff2Pt1TOF);
    hEffTOF2Pt1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatchTOF));  

    hEffTOF2Pt10->SetBinContent(i+1,Eff2Pt10TOF);
    hEffTOF2Pt10->SetBinError(i+1,errEff2Pt10TOF);
    hEffTOF2Pt10->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatchTOF));  

    hEffTOFSPDPt02->SetBinContent(i+1,EffSPDPt02TOF);
    hEffTOFSPDPt02->SetBinError(i+1,errEffSPDPt02TOF);
    hEffTOFSPDPt02->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatchTOF));  
      
    hEffTOFSPDPt1->SetBinContent(i+1,EffSPDPt1TOF);
    hEffTOFSPDPt1->SetBinError(i+1,errEffSPDPt1TOF);
    hEffTOFSPDPt1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatchTOF)); 

    hEffTOFSPDPt10->SetBinContent(i+1,EffSPDPt10TOF);
    hEffTOFSPDPt10->SetBinError(i+1,errEffSPDPt10TOF);
    hEffTOFSPDPt10->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatchTOF)); 
     
    hEffTOFoneSPDPt02->SetBinContent(i+1,EffoneSPDPt02TOF);
    hEffTOFoneSPDPt02->SetBinError(i+1,errEffoneSPDPt02TOF);
    hEffTOFoneSPDPt02->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatchTOF));

    hEffTOFoneSPDPt1->SetBinContent(i+1,EffoneSPDPt1TOF);
    hEffTOFoneSPDPt1->SetBinError(i+1,errEffoneSPDPt1TOF);
    hEffTOFoneSPDPt1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatchTOF));      

    hEffTOFoneSPDPt10->SetBinContent(i+1,EffoneSPDPt10TOF);
    hEffTOFoneSPDPt10->SetBinError(i+1,errEffoneSPDPt10TOF);
    hEffTOFoneSPDPt10->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatchTOF));

    hEffTOFTOTPt02->SetBinContent(i+1,EffTOTPt02TOF);
    hEffTOFTOTPt02->SetBinError(i+1,errEffTOTPt02TOF);
    hEffTOFTOTPt02->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatchTOF));
      
    hEffTOFTOTPt1->SetBinContent(i+1,EffTOTPt1TOF);
    hEffTOFTOTPt1->SetBinError(i+1,errEffTOTPt1TOF);
    hEffTOFTOTPt1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatchTOF));
      
    hEffTOFTOTPt10->SetBinContent(i+1,EffTOTPt10TOF);
    hEffTOFTOTPt10->SetBinError(i+1,errEffTOTPt10TOF);
    hEffTOFTOTPt10->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatchTOF));  
  }


  //-----------------------------------

  //---vertex

  TNtuple* ntvertex=(TNtuple*)fil->Get("ntvertex");

  Float_t nrunVertex,Vx,errVx,sigmaVx,errsigmaVx,Vy,errVy,sigmaVy,errsigmaVy,Vz,errVz,sigmaVz,errsigmaVz,VxSPD,errVxSPD,sigmaVxSPD,errsigmaVxSPD,VySPD,errVySPD,sigmaVySPD,errsigmaVySPD,VzSPD,errVzSPD,sigmaVzSPD,errsigmaVzSPD,pileupSPD,errpileupSPD;

  ntvertex->SetBranchAddress("nrun",&nrunVertex);
  ntvertex->SetBranchAddress("VxTRK",&Vx);
  ntvertex->SetBranchAddress("errVxTRK",&errVx);
  ntvertex->SetBranchAddress("sigmaVxTRK",&sigmaVx);
  ntvertex->SetBranchAddress("errsigmaVxTRK",&errsigmaVx);
  ntvertex->SetBranchAddress("VyTRK",&Vy);
  ntvertex->SetBranchAddress("errVyTRK",&errVy);
  ntvertex->SetBranchAddress("sigmaVyTRK",&sigmaVy);
  ntvertex->SetBranchAddress("errsigmaVyTRK",&errsigmaVy);
  ntvertex->SetBranchAddress("VzTRK",&Vz);
  ntvertex->SetBranchAddress("errVzTRK",&errVz);
  ntvertex->SetBranchAddress("sigmaVzTRK",&sigmaVz);
  ntvertex->SetBranchAddress("errsigmaVzTRK",&errsigmaVz);
  ntvertex->SetBranchAddress("VxSPD",&VxSPD);
  ntvertex->SetBranchAddress("errVxSPD",&errVxSPD);
  ntvertex->SetBranchAddress("sigmaVxSPD",&sigmaVxSPD);
  ntvertex->SetBranchAddress("errsigmaVxSPD",&errsigmaVxSPD);
  ntvertex->SetBranchAddress("VySPD",&VySPD);
  ntvertex->SetBranchAddress("errVySPD",&errVySPD);
  ntvertex->SetBranchAddress("sigmaVySPD",&sigmaVySPD);
  ntvertex->SetBranchAddress("errsigmaVySPD",&errsigmaVySPD);
  ntvertex->SetBranchAddress("VzSPD",&VzSPD);
  ntvertex->SetBranchAddress("errVzSPD",&errVzSPD);
  ntvertex->SetBranchAddress("sigmaVzSPD",&sigmaVzSPD);
  ntvertex->SetBranchAddress("errsigmaVzSPD",&errsigmaVzSPD);
    ntvertex->SetBranchAddress("pileupSPD",&pileupSPD);
    ntvertex->SetBranchAddress("errpileupSPD",&errpileupSPD);

  nr=ntvertex->GetEntries();
  delete []myIndex;
  delete []noRuns;
  myIndex = new Int_t [nr];
  noRuns = new Int_t [nr];
  for(Int_t i=0; i<nr;i++){
    ntvertex->GetEvent(i);
    Int_t intrun = static_cast<Int_t>(nrunVertex+0.01);
    noRuns[i]=intrun;
  }
  printf("\n ======== PROCESSING VERTEX NTUPLE \n");
  kRunsToPlot = RunsToPlot(run1,run2,nr,noRuns,myIndex);


  TH1F *hVx = new TH1F("hVx","Track Vertex Vx Distribution",kRunsToPlot,0.,kRunsToPlot);
  hVx->SetLineWidth(2);
  hVx->SetLineColor(kBlue+2);
  hVx->SetMarkerColor(kBlue+2);
  hVx->SetMarkerStyle(20);

 TH1F *hVy = new TH1F("hVy","Track Vertex Vy Distribution",kRunsToPlot,0.,kRunsToPlot);
  hVy->SetLineWidth(2);
  hVy->SetLineColor(kBlue+2);
  hVy->SetMarkerColor(kBlue+2);
  hVy->SetMarkerStyle(20);

 TH1F *hVz = new TH1F("hVz","Track Vertex Vz Distribution",kRunsToPlot,0.,kRunsToPlot);
  hVz->SetLineWidth(2);
  hVz->SetLineColor(kBlue+2);
  hVz->SetMarkerColor(kBlue+2);
  hVz->SetMarkerStyle(20);

  TH1F *hSigmaVx = new TH1F("hSigmaVx","Track Vertex SigmaVx Distribution",kRunsToPlot,0.,kRunsToPlot);
  hSigmaVx->SetLineWidth(2);
  hSigmaVx->SetLineColor(kBlue+2);
  hSigmaVx->SetMarkerColor(kBlue+2);
  hSigmaVx->SetMarkerStyle(20);

 TH1F *hSigmaVy = new TH1F("hSigmaVy","Track Vertex SigmaVy Distribution",kRunsToPlot,0.,kRunsToPlot);
  hSigmaVy->SetLineWidth(2);
  hSigmaVy->SetLineColor(kBlue+2);
  hSigmaVy->SetMarkerColor(kBlue+2);
  hSigmaVy->SetMarkerStyle(20);

 TH1F *hSigmaVz = new TH1F("hSigmaVz","Track Vertex SigmaVz Distribution",kRunsToPlot,0.,kRunsToPlot);
  hSigmaVz->SetLineWidth(2);
  hSigmaVz->SetLineColor(kBlue+2);
  hSigmaVz->SetMarkerColor(kBlue+2);
  hSigmaVz->SetMarkerStyle(20);

 TH1F *hVxSPD = new TH1F("hVxSPD","Track Vertex Vx Distribution",kRunsToPlot,0.,kRunsToPlot);
  hVxSPD->SetLineWidth(2);
  hVxSPD->SetLineColor(2);
  hVxSPD->SetMarkerColor(2);
  hVxSPD->SetMarkerStyle(20);

 TH1F *hVySPD = new TH1F("hVySPD","Track Vertex Vy Distribution",kRunsToPlot,0.,kRunsToPlot);
  hVySPD->SetLineWidth(2);
  hVySPD->SetLineColor(2);
  hVySPD->SetMarkerColor(2);
  hVySPD->SetMarkerStyle(20);

 TH1F *hVzSPD = new TH1F("hVzSPD","Track Vertex Vz Distribution",kRunsToPlot,0.,kRunsToPlot);
  hVzSPD->SetLineWidth(2);
  hVzSPD->SetLineColor(2);
  hVzSPD->SetMarkerColor(2);
  hVzSPD->SetMarkerStyle(20);

  TH1F *hSigmaVxSPD = new TH1F("hSigmaVxSPD","Track Vertex SigmaVx Distribution",kRunsToPlot,0.,kRunsToPlot);
  hSigmaVxSPD->SetLineWidth(2);
  hSigmaVxSPD->SetLineColor(2);
  hSigmaVxSPD->SetMarkerColor(2);
  hSigmaVxSPD->SetMarkerStyle(20);

 TH1F *hSigmaVySPD = new TH1F("hSigmaVySPD","Track Vertex SigmaVy Distribution",kRunsToPlot,0.,kRunsToPlot);
  hSigmaVySPD->SetLineWidth(2);
  hSigmaVySPD->SetLineColor(2);
  hSigmaVySPD->SetMarkerColor(2);
  hSigmaVySPD->SetMarkerStyle(20);

 TH1F *hSigmaVzSPD = new TH1F("hSigmaVzSPD","Track Vertex SigmaVz Distribution",kRunsToPlot,0.,kRunsToPlot);
  hSigmaVzSPD->SetLineWidth(2);
  hSigmaVzSPD->SetLineColor(2);
  hSigmaVzSPD->SetMarkerColor(2);
  hSigmaVzSPD->SetMarkerStyle(20);

TH1F *hpileupSPD = new TH1F("hpileupSPD","Fraction of tracks with SPD pileup",kRunsToPlot,0.,kRunsToPlot);
   hpileupSPD->SetMarkerStyle(20);

  //  Int_t nEntriesVertex=ntvertex->GetEntries();
  
  for(Int_t i=0;i<kRunsToPlot;i++){

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

    hVxSPD->SetBinContent(i+1,VxSPD);
    hVxSPD->SetBinError(i+1,errVxSPD);
    hVxSPD->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunVertex));
  
    hVySPD->SetBinContent(i+1,VySPD);
    hVySPD->SetBinError(i+1,errVySPD);
    hVySPD->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunVertex));
    
    hVzSPD->SetBinContent(i+1,VzSPD);
    hVzSPD->SetBinError(i+1,errVzSPD);
    hVzSPD->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunVertex));
  
    
    hSigmaVxSPD->SetBinContent(i+1,sigmaVxSPD);
    hSigmaVxSPD->SetBinError(i+1,errsigmaVxSPD);
    hSigmaVxSPD->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunVertex));

    hSigmaVySPD->SetBinContent(i+1,sigmaVySPD);
    hSigmaVySPD->SetBinError(i+1,errsigmaVySPD);
    hSigmaVySPD->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunVertex));

    hSigmaVzSPD->SetBinContent(i+1,sigmaVzSPD);
    hSigmaVzSPD->SetBinError(i+1,errsigmaVzSPD);
    hSigmaVzSPD->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunVertex));

      hpileupSPD->SetBinContent(i+1,pileupSPD);
      hpileupSPD->SetBinError(i+1,errpileupSPD);
      hpileupSPD->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunVertex));

  }
 //--------  Draw Vertex histograms ---------
    
  TCanvas *cVertexDisto=new TCanvas("cVertexDisto","cVertexDisto",1200,800);
  cVertexDisto->Divide(3,2);
  cVertexDisto->cd(1);
  hVx->SetMinimum(0.06);
  hVx->SetMaximum(0.08);
    hVx->GetYaxis()->SetTitle("Vertex X coordinate");
    hVx->GetXaxis()->SetTitle("run number");
  hVx->Draw();
  hVxSPD->Draw("same");
  cVertexDisto->cd(2);
  hVy->SetMinimum(0.36);
  hVy->SetMaximum(0.37);
  hVy->GetYaxis()->SetTitle("Vertex Y coordinate");
    hVy->GetXaxis()->SetTitle("run number");
  hVy->Draw();
  hVySPD->Draw("same");
  cVertexDisto->cd(3);
//   hVz->SetMinimum(-1.);
//   hVz->SetMaximum(1.);
  hVz->GetYaxis()->SetTitle("Vertex Z coordinate");
    hVz->GetXaxis()->SetTitle("run number");
  hVz->Draw();
  hVzSPD->Draw("same");
  cVertexDisto->cd(4);
  hSigmaVx->SetMinimum(0.);
  hSigmaVx->SetMaximum(0.1);
  hSigmaVx->GetYaxis()->SetTitle("sigma on x coordinate");
    hSigmaVx->GetXaxis()->SetTitle("run number");
  hSigmaVx->Draw();
  hSigmaVxSPD->Draw("same");
  cVertexDisto->cd(5);
 hSigmaVy->SetMinimum(0.);
 hSigmaVy->SetMaximum(0.2);
  hSigmaVy->GetYaxis()->SetTitle("sigma on y coordinate");
    hSigmaVy->GetXaxis()->SetTitle("run number");
  hSigmaVy->Draw();
  hSigmaVySPD->Draw("same");
  cVertexDisto->cd(6);
//   hSigmaVz->SetMinimum(6.);
//   hSigmaVz->SetMaximum(10.);
  hSigmaVz->GetYaxis()->SetTitle("sigma on z coordinate");
    hSigmaVz->GetXaxis()->SetTitle("run number");
  hSigmaVz->Draw();
  hSigmaVzSPD->Draw("same");
  cVertexDisto->SaveAs("Vertex_trend.pdf");
//    pdfFileNames+=" Vertex_trend.pdf";
    
    
    //---Pileup


    TNtuple* ntpu=(TNtuple*)fil->Get("ntPileUp");
    
    Float_t nrunpu,npilvtx,errnpilvtx,ntrklpil,errntrklpil,ntrklnopil,errntrklnopil,ncl1pil,errncl1pil,ncl1nopil,errncl1nopil;
    
    ntpu->SetBranchAddress("run",&nrunpu);
    ntpu->SetBranchAddress("npilvtx",&npilvtx);
    ntpu->SetBranchAddress("errnpilvtx",&errnpilvtx);
    ntpu->SetBranchAddress("ntrklpil",&ntrklpil);
    ntpu->SetBranchAddress("errntrklpil",&errntrklpil);
    ntpu->SetBranchAddress("ntrklnopil",&ntrklnopil);
    ntpu->SetBranchAddress("errntrklnopil",&errntrklnopil);
    ntpu->SetBranchAddress("ncl1pil",&ncl1pil);
    ntpu->SetBranchAddress("errncl1pil",&errncl1pil);
    ntpu->SetBranchAddress("ncl1nopil",&ncl1nopil);
    ntpu->SetBranchAddress("errncl1nopil",&errncl1nopil);
    
    nr=ntpu->GetEntries();
    delete []myIndex;
    delete []noRuns;
    myIndex = new Int_t [nr];
    noRuns = new Int_t [nr];
    for(Int_t i=0; i<nr;i++){
        ntpu->GetEvent(i);
        Int_t intrun = static_cast<Int_t>(nrunpu+0.01);
        noRuns[i]=intrun;
    }
    
    printf("\n ======== PROCESSING PILEUP NTUPLE \n");
    kRunsToPlot = RunsToPlot(run1,run2,nr,noRuns,myIndex);
    
    
    TH1F *hNumPilVtx = new TH1F("hNumPilVtx","Mean pileup vertices number",kRunsToPlot,0.,kRunsToPlot);
    hNumPilVtx->SetLineWidth(2);
    hNumPilVtx->SetLineColor(kBlue+2);
    hNumPilVtx->SetMarkerColor(kBlue+2);
    hNumPilVtx->SetMarkerStyle(20);
    hNumPilVtx->GetYaxis()->SetRangeUser(0.,3.);
    hNumPilVtx->GetXaxis()->SetTitle("run number");
    hNumPilVtx->GetYaxis()->SetTitle("pileup vertices mean number");
    

    TH1F *hNumPilTrkl = new TH1F("hNumPilTrkl","Mean Tracklets number for prim. vertex w/ pileup",kRunsToPlot,0.,kRunsToPlot);
    hNumPilTrkl->SetLineWidth(2);
    hNumPilTrkl->SetLineColor(2);
    hNumPilTrkl->SetMarkerColor(2);
    hNumPilTrkl->SetMarkerStyle(20);
    hNumPilTrkl->GetYaxis()->SetRangeUser(0.,60.);
    hNumPilTrkl->GetXaxis()->SetTitle("run number");
    hNumPilTrkl->GetYaxis()->SetTitle("pileup vertex Tracklet number");

    TH1F *hNumNoPilTrkl = new TH1F("hNumNoPilTrkl","Mean Tracklets number for prim. vertex w/out pileup",kRunsToPlot,0.,kRunsToPlot);
    hNumNoPilTrkl->SetLineWidth(2);
    hNumNoPilTrkl->SetLineColor(kBlue+2);
    hNumNoPilTrkl->SetMarkerColor(kBlue+2);
    hNumNoPilTrkl->SetMarkerStyle(24);

    TH1F *hNumPilCL1 = new TH1F("hNumPilCL1","Mean SPD1 cluster number for prim. vertex w/ pileup tagging",kRunsToPlot,0.,kRunsToPlot);
    hNumPilCL1->SetLineWidth(2);
    hNumPilCL1->SetLineColor(2);
    hNumPilCL1->SetMarkerColor(2);
    hNumPilCL1->SetMarkerStyle(20);
    hNumPilCL1->GetYaxis()->SetRangeUser(0.,100.);
    hNumPilCL1->GetXaxis()->SetTitle("run number");
    hNumPilCL1->GetYaxis()->SetTitle("pileup vertex SPD1 cluster number");
    
    TH1F *hNumNoPilCL1 = new TH1F("hNumNoPilCL1","Mean SPD1 cluster number for prim. vertex w/out pileup tagging",kRunsToPlot,0.,kRunsToPlot);
    hNumNoPilCL1->SetLineWidth(2);
    hNumNoPilCL1->SetLineColor(kBlue+2);
    hNumNoPilCL1->SetMarkerColor(kBlue+2);
    hNumNoPilCL1->SetMarkerStyle(24);
    
    
    for(Int_t i=0;i<kRunsToPlot;i++){
        
        ntpu->GetEvent(myIndex[i]);
        
        //    cout<<Vx<<endl;
        
        hNumPilVtx->SetBinContent(i+1,npilvtx);
        hNumPilVtx->SetBinError(i+1,errnpilvtx);
        hNumPilVtx->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunpu));
        printf("number of pileup vtx for run %d = %f \n",(Int_t)nrunpu, npilvtx);
//        cout << "run " << nrunpu << " numero di vertici PU " << npilvtx << endl;

        hNumPilTrkl->SetBinContent(i+1,ntrklpil);
        hNumPilTrkl->SetBinError(i+1,errntrklpil);
        hNumPilTrkl->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunpu));
    
        hNumNoPilTrkl->SetBinContent(i+1,ntrklnopil);
        hNumNoPilTrkl->SetBinError(i+1,errntrklnopil);
        hNumNoPilTrkl->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunpu));
        
        hNumPilCL1->SetBinContent(i+1,ncl1pil);
        hNumPilCL1->SetBinError(i+1,errncl1pil);
        hNumPilCL1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunpu));
        
        hNumNoPilCL1->SetBinContent(i+1,ncl1nopil);
        hNumNoPilCL1->SetBinError(i+1,errncl1nopil);
        hNumNoPilCL1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunpu));
    }
    
 //-----------------------------------

  gStyle->SetOptStat(0);

//
    TCanvas* cMI=new TCanvas("cMI"," Reconstructed Track With points in ITS");
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
    
    //
    TCanvas* cSA=new TCanvas("cSA"," SA Track With points in ITS");
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

    //
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
    
//
    TCanvas* cev=new TCanvas("cev","Fraction of Events with SDD in Trigger Cluster");
    histoEvwSDD->SetLineColor(kBlue);
    histoEvwSDD->SetMarkerColor(kBlue);
    histoEvwSDD->SetMarkerStyle(22);
    histoEvwSDD->Draw();
    histoEvwSDD->GetYaxis()->SetTitle("Fraction of Events with SDD in Trigger Cluster");
    histoEvwSDD->GetYaxis()->SetRangeUser(0.,1.2);
    histoEvwSDD->GetXaxis()->SetTitle("run number");

    cev->SaveAs("NoFast_trend.pdf");
//    pdfFileNames+=" NoFast_trend.pdf";
    cev->Update();
//
    
    TCanvas* c22=new TCanvas("c22","Number of Events in Run",1000,600);
    histonEvents->SetTitle("Run Events number");
    histonEvents->Draw();
    histonEvents->GetYaxis()->SetTitleOffset(1.2);
    histonEvents->GetYaxis()->SetTitle("number of events");
    histonEvents->GetXaxis()->SetTitle("run number");
    histonEventsTriggered->Draw("same");
    c22->SaveAs("RunEvents_trend.pdf");
//    pdfFileNames+=" RunEvents_trend.pdf";
    c22->Update();

//
    
    TCanvas* cfrac=new TCanvas("cfrac","Fraction of SDD modules ON",900,900);
    cfrac->Divide(1,2);
    cfrac->cd(1);
    histoFracDead3->SetMarkerStyle(20);
    histoFracDead3->SetMarkerColor(kOrange+1);
    histoFracDead3->SetLineColor(kOrange+1);
    histoFracDead3->GetYaxis()->SetRangeUser(0.,1.2);
    histoFracDead3->Draw();
    histoFracDead3->GetYaxis()->SetTitle("Fraction of Modules ON");
    histoFracDead3->GetXaxis()->SetTitle("run number");
    TLatex* tf3=new TLatex(0.2,0.8,"SDD modules ON - Layer 3 (total: 84)");
    tf3->SetNDC();
    tf3->SetTextColor(1);
    tf3->Draw();
    cfrac->cd(2);
    histoFracDead4->SetMarkerStyle(20);
    histoFracDead4->SetMarkerColor(kAzure+1);
    histoFracDead4->SetLineColor(kAzure+1);
    histoFracDead4->GetYaxis()->SetRangeUser(0.,1.2);
    histoFracDead4->Draw();
    histoFracDead4->GetYaxis()->SetTitle("Fraction of Modules ON");
    histoFracDead4->GetXaxis()->SetTitle("run number");
    TLatex* tf4=new TLatex(0.2,0.8,"SDD modules ON - Layer 4 (total: 176)");
    tf4->SetNDC();
    tf4->SetTextColor(1);
    tf4->Draw();

    cfrac->SaveAs("SDDmodulesON_trend.pdf");
//    pdfFileNames+=" SDDmodulesON_trend.pdf";
    cfrac->Update();

    
    //

  TCanvas* c2=new TCanvas("c2","SDD DriftTime & Charge",1200,800);
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
  ent=leg2->AddEntry(histodEdxTB0,"Small drift time","PL");
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

  TCanvas *c7=new TCanvas("c7","Charge ratio");
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
//  pdfFileNames+=" SSD_chargeratio_trend.pdf";
  //

    TCanvas *c8=new TCanvas("c8","Fraction of SSD bad strips",900,900);
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
    ent=legBad->AddEntry(histoFracBadn5,"Layer5 n-side","PL");
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
    
    
  TCanvas *cpt02=new TCanvas("cpt02","TPC-ITS matching efficiency",1200,1000);
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
    
    // TCanvas *cpt1=new TCanvas("cpt1","cpt1");
    // cpt1->cd(1);
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
    
    TFile *matching_histo =new TFile("Match.root","RECREATE");
    cpt02->Write();
    matching_histo->Close(); 
    
    cpt02->SaveAs("TPCTOFMatch_trend.pdf");
//    pdfFileNames+=" TPCTOFMatch_trend.pdf";
    

    TCanvas *cPileUp=new TCanvas("cPileUp","cPileUp",1200,800);
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
/*    TLegend* lpt1=new TLegend(0.9,0.8,1,1);
    lpt1->AddEntry(hEff6Pt1,"6 cls","l");
    lpt1->AddEntry(hEff5Pt1,"5 cls","l");
    lpt1->AddEntry(hEff4Pt1,"4 cls","l");
    lpt1->AddEntry(hEff3Pt1,"3 cls","l");
    lpt1->AddEntry(hEff2Pt1,"2 cls","l");
    lpt1->AddEntry(hEffSPDPt1,"2SPD + any","l");
    lpt1->AddEntry(hEffoneSPDPt1,">=1SPD + any","l");
    lpt1->AddEntry(hEffTOTPt02,">=2","l");
 */
    lpt1->Draw("same");
/*    TLatex* tpc2=new TLatex(0.2,0.75,"TPCITS match eff Pt=1");
    tpc2->SetNDC();
    tpc2->SetTextColor(1);
*/
    tpc2->Draw();
    cPileUp->SaveAs("Pileup_trend.pdf");
//    pdfFileNames+=" Pileup_trend.pdf";

    gStyle->SetOptStat(0);
    
    TCanvas* cpu=new TCanvas("cPileUpVtx","Pileup vertices from SPD",1200,800);
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

    
    
    //-----------------------------------
    TCanvas *cPixel=new TCanvas("cPixel","SPD on");
    
  cPixel->cd(1);
  hFracSPD1->SetMaximum(1.2);
  hFracSPD1->SetMinimum(0);
  hFracSPD1->Draw("p");
  hFracSPD2->Draw("same,p");

  TLegend* lSPD=new TLegend(0.9,0.8,1,1);
  lSPD->AddEntry(hFracSPD1,"Frac. SPD1 ON","l");
  lSPD->AddEntry(hFracSPD2,"Frac. SPD2 ON","l");
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
  cPixel->SaveAs("Pixel_trend.pdf");
//  pdfFileNames+=" Pixel_trend.pdf";
  cPixel->Update();
    
   //-----------  ITS SA -----------
//  PlotITSSA(fil,run1,run2);

    Double_t Lowbin[3]={0.1,0.5,0.9};
    Double_t Upbin[3]={0.2,0.6,1};
//    gROOT->SetStyle("Plain");
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetFillColor(0);
//    gStyle->SetTextFont(32);
    
    TNtuple* nt = (TNtuple*)fil->Get("ntITSsa");
    
    Float_t ITSA[3];
    Float_t TPIT[3];
    Float_t RAT[3];
    Float_t run;
    Float_t NcluITSpSA,errNcluITSpSA,dedx4_3,errdedx4_3,PtpionpSA,errPtpionpSA;
    Float_t NclupSA0,errNclupSA0,NclupSA1,errNclupSA1,NclupSA2,errNclupSA2,NclupSA3,errNclupSA3,NclupSA4,errNclupSA4,NclupSA5,errNclupSA5;
    nt->SetBranchAddress("run",&run);
    //  nt->SetBranchAddress("NITSsaPtBin0",&ITSA[0]); //Pb-Pb
    // nt->SetBranchAddress("NITSsaPtBin1",&ITSA[1]); //Pb-Pb
    // nt->SetBranchAddress("NITSsaPtBin2",&ITSA[2]); //Pb-Pb
    nt->SetBranchAddress("NITSpureSAPtBin0",&ITSA[0]);
    nt->SetBranchAddress("NITSpureSAPtBin1",&ITSA[1]);
    nt->SetBranchAddress("NITSpureSAPtBin2",&ITSA[2]);
    nt->SetBranchAddress("NITSTPCPtBin0",&TPIT[0]);
    nt->SetBranchAddress("NITSTPCPtBin1",&TPIT[1]);
    nt->SetBranchAddress("NITSTPCPtBin2",&TPIT[2]);
    nt->SetBranchAddress("ratioPtBin0",&RAT[0]);
    nt->SetBranchAddress("ratioPtBin1",&RAT[1]);
    nt->SetBranchAddress("ratioPtBin2",&RAT[2]);
    nt->SetBranchAddress("NcluITSpSA",&NcluITSpSA);
    nt->SetBranchAddress("errNcluITSpSA",&errNcluITSpSA);
    nt->SetBranchAddress("dedx4_3",&dedx4_3);
    nt->SetBranchAddress("errdedx4_3",&errdedx4_3);
    nt->SetBranchAddress("PtpionpSA",&PtpionpSA);
    nt->SetBranchAddress("errPtpionpSA",&errPtpionpSA);
    nt->SetBranchAddress("NclupSA0",&NclupSA0);
    nt->SetBranchAddress("errNclupSA0",&errNclupSA0);
    nt->SetBranchAddress("NclupSA1",&NclupSA1);
    nt->SetBranchAddress("errNclupSA1",&errNclupSA1);
    nt->SetBranchAddress("NclupSA2",&NclupSA2);
    nt->SetBranchAddress("errNclupSA2",&errNclupSA2);
    nt->SetBranchAddress("NclupSA3",&NclupSA3);
    nt->SetBranchAddress("errNclupSA3",&errNclupSA3);
    nt->SetBranchAddress("NclupSA4",&NclupSA4);
    nt->SetBranchAddress("errNclupSA4",&errNclupSA4);
    nt->SetBranchAddress("NclupSA5",&NclupSA5);
    nt->SetBranchAddress("errNclupSA5",&errNclupSA5);
    cout << "SPD1 " << NclupSA0 << " errore " << errNclupSA0 << endl;
    
    Int_t nru=nt->GetEntries();
    Int_t *myIndex2 = new Int_t [nru];
    Int_t *noRuns2 = new Int_t [nru];
    for(Int_t i=0; i<nru;i++){
        nt->GetEvent(i);
        Int_t intrun = static_cast<Int_t>(run+0.01);
        noRuns2[i]=intrun;
    }
    printf("\n ======== PROCESSING ITS SA NTUPLE \n");
//    Int_t kRunsToPlot = RunsToPlot(run1,run2,nr,noRuns2,myIndex2);
    kRunsToPlot = RunsToPlot(run1,run2,nr,noRuns2,myIndex2);
    TH1F *h0=new TH1F("h0","h0",kRunsToPlot,-0.5,kRunsToPlot-0.5);
    TH1F *h1=new TH1F("h1","h1",kRunsToPlot,-0.5,kRunsToPlot-0.5);
    TH1F *h2=new TH1F("h2","h2",kRunsToPlot,-0.5,kRunsToPlot-0.5);
    TH1F *h3=new TH1F("h3","h4",kRunsToPlot,-0.5,kRunsToPlot-0.5);
    TH1F *h4=new TH1F("h4","h5",kRunsToPlot,-0.5,kRunsToPlot-0.5);
    TH1F *h5=new TH1F("h5","h5",kRunsToPlot,-0.5,kRunsToPlot-0.5);
    TH1F *h6=new TH1F("h6","h6",kRunsToPlot,-0.5,kRunsToPlot-0.5);
    TH1F *h7=new TH1F("h7","h7",kRunsToPlot,-0.5,kRunsToPlot-0.5);
    TH1F *h8=new TH1F("h8","h8",kRunsToPlot,-0.5,kRunsToPlot-0.5);
    TH1F *h9=new TH1F("h9","h9",kRunsToPlot,-0.5,kRunsToPlot-0.5);
    TH1F *h10=new TH1F("h10","h10",kRunsToPlot,-0.5,kRunsToPlot-0.5);
    TH1F *h11=new TH1F("h11","h11",kRunsToPlot,-0.5,kRunsToPlot-0.5);
    TH1F *h12=new TH1F("h12","h12",kRunsToPlot,-0.5,kRunsToPlot-0.5);
    TH1F *h13=new TH1F("h13","h13",kRunsToPlot,-0.5,kRunsToPlot-0.5);
    TH1F *h14=new TH1F("h14","h14",kRunsToPlot,-0.5,kRunsToPlot-0.5);
    TH1F *h15=new TH1F("h15","h15",kRunsToPlot,-0.5,kRunsToPlot-0.5);
    TH1F *h16=new TH1F("h16","h16",kRunsToPlot,-0.5,kRunsToPlot-0.5);
    TH1F *h17=new TH1F("h17","h17",kRunsToPlot,-0.5,kRunsToPlot-0.5);
    
    
    for(Int_t iev=0;iev<kRunsToPlot;iev++){
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
        h9->SetBinContent(iev+1,NcluITSpSA);
        h9->SetBinError(iev+1,errNcluITSpSA);
        h10->SetBinContent(iev+1,dedx4_3);
        h10->SetBinError(iev+1,errdedx4_3);
        h11->SetBinContent(iev+1,PtpionpSA);
        h11->SetBinError(iev+1,errPtpionpSA);
        h12->SetBinContent(iev+1,NclupSA0);
        h12->SetBinError(iev+1,errNclupSA0);
        h13->SetBinContent(iev+1,NclupSA1);
        h13->SetBinError(iev+1,errNclupSA1);
        h14->SetBinContent(iev+1,NclupSA2);
        h14->SetBinError(iev+1,errNclupSA2);
        h15->SetBinContent(iev+1,NclupSA3);
        h15->SetBinError(iev+1,errNclupSA3);
        h16->SetBinContent(iev+1,NclupSA4);
        h16->SetBinError(iev+1,errNclupSA4);
        h17->SetBinContent(iev+1,NclupSA5);
        h17->SetBinError(iev+1,errNclupSA5);
        
        h0->GetXaxis()->SetBinLabel(iev+1,Form("%.0f",run));
        h1->GetXaxis()->SetBinLabel(iev+1,Form("%.0f",run));
        h2->GetXaxis()->SetBinLabel(iev+1,Form("%.0f",run));
        h3->GetXaxis()->SetBinLabel(iev+1,Form("%.0f",run));
        h4->GetXaxis()->SetBinLabel(iev+1,Form("%.0f",run));
        h5->GetXaxis()->SetBinLabel(iev+1,Form("%.0f",run));
        h6->GetXaxis()->SetBinLabel(iev+1,Form("%.0f",run));
        h7->GetXaxis()->SetBinLabel(iev+1,Form("%.0f",run));
        h8->GetXaxis()->SetBinLabel(iev+1,Form("%.0f",run));
        h9->GetXaxis()->SetBinLabel(iev+1,Form("%d",(Int_t)run));
        h10->GetXaxis()->SetBinLabel(iev+1,Form("%d",(Int_t)run));
        h11->GetXaxis()->SetBinLabel(iev+1,Form("%d",(Int_t)run));
        h12->GetXaxis()->SetBinLabel(iev+1,Form("%d",(Int_t)run));
        h13->GetXaxis()->SetBinLabel(iev+1,Form("%d",(Int_t)run));
        h14->GetXaxis()->SetBinLabel(iev+1,Form("%d",(Int_t)run));
        h15->GetXaxis()->SetBinLabel(iev+1,Form("%d",(Int_t)run));
        h16->GetXaxis()->SetBinLabel(iev+1,Form("%d",(Int_t)run));
        h17->GetXaxis()->SetBinLabel(iev+1,Form("%d",(Int_t)run));
        
        //    Printf("%f   %f   %f",ITSA[0],ITSA[1],ITSA[2]);
    }
    h0->Print("all");
    // h0->GetYaxis()->SetTitle("ITSsa tracks");
    // h0->GetXaxis()->SetTitle("run");
    // h1->GetYaxis()->SetTitle("ITSsa tracks");
    // h1->GetXaxis()->SetTitle("run");
    // h2->GetYaxis()->SetTitle("ITSsa tracks");
    // h2->GetXaxis()->SetTitle("run");
    h0->GetYaxis()->SetTitle("ITSpureSA tracks");
    h0->GetXaxis()->SetTitle("run number");
    h1->GetYaxis()->SetTitle("ITSpureSA tracks");
    h1->GetXaxis()->SetTitle("run number");
    h2->GetYaxis()->SetTitle("ITSpureSA tracks");
    h2->GetXaxis()->SetTitle("run number");
    h3->GetYaxis()->SetTitle("ITS+TPC tracks");
    h3->GetXaxis()->SetTitle("run number");
    h4->GetYaxis()->SetTitle("ITS+TPC tracks");
    h4->GetXaxis()->SetTitle("run number");
    h5->GetYaxis()->SetTitle("ITS+TPC tracks");
    h5->GetXaxis()->SetTitle("run number");
    // h6->GetYaxis()->SetTitle("(TPC+ITS)/ITS");
    // h6->GetXaxis()->SetTitle("run");
    // h7->GetYaxis()->SetTitle("(TPC+ITS)/ITS");
    // h7->GetXaxis()->SetTitle("run");
    // h8->GetYaxis()->SetTitle("(TPC+ITS)/ITS");
    // h8->GetXaxis()->SetTitle("run");
    h6->GetYaxis()->SetTitle("(TPC+ITS)+ITSsa/ITSpureSA");
    h6->GetXaxis()->SetTitle("run number");
    h7->GetYaxis()->SetTitle("(TPC+ITS)+ITSsa/ITSpureSA");
    h7->GetXaxis()->SetTitle("run number");
    h8->GetYaxis()->SetTitle("(TPC+ITS)+ITSsa/ITSpureSA");
    h8->GetXaxis()->SetTitle("run number");
    h9->GetXaxis()->SetTitle("run number");
    h9->GetYaxis()->SetTitle("mean cluster number (N>3)");
    h10->GetXaxis()->SetTitle("run number");
    h10->GetYaxis()->SetTitle("N(dedx4clu)/N(dedx3clu)");
    h11->GetXaxis()->SetTitle("run number");
    h11->GetYaxis()->SetTitle("mean Pt (GeV/c)");
    h12->GetXaxis()->SetTitle("run number");
    h12->GetYaxis()->SetTitle("Fraction of pureSA tracks with cluster in ITS layers");
    h13->GetXaxis()->SetTitle("run number");
    h14->GetXaxis()->SetTitle("run number");
    h15->GetXaxis()->SetTitle("run number");
    h16->GetXaxis()->SetTitle("run number");
    h17->GetXaxis()->SetTitle("run number");
    
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
    
    TCanvas *c=new TCanvas("ITS pure SA tracks","ITS pure SA tracks");
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
    
    TCanvas *c2a=new TCanvas("ITS+TPC tracks","ITS+TPC tracks");
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
    
    TCanvas *c3=new TCanvas("(ITSTPC+ITSsa)/ITSpureSA","(ITSTPC+ITSsa)/ITSpureSA");
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
    
    TString name1="File_tracksITS.root";
    
    TFile *f=new TFile(name1,"RECREATE");
    f->cd();
    c->Write();
    c2->Write();
    c3->Write();
    f->Close();
    
    TCanvas* cSAb=new TCanvas("cSAb","Mean cluster number (N>3) ITSsa");
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
    
    TCanvas* cSA2=new TCanvas("cSA2","dedx(4clu)/dedx(3clu) ITSsa");
    h10->SetMarkerStyle(20);
    h10->SetMarkerColor(3);
    h10->Draw();
    tl->Draw();
    cSA2->SaveAs("dedx4_3_SA_trend.pdf");
    //    pdfFileNames+=" dedx4_3_SA_trend.pdf";
    
    TCanvas* cSA3=new TCanvas("cSA3","mean pion Pt ITSsa");
    h11->SetMarkerStyle(20);
    h11->SetMarkerColor(4);
    h11->Draw();
    tl->Draw();
    cSA3->SaveAs("piPt_SA_trend.pdf");
    //    pdfFileNames+=" piPt_SA_trend.pdf";
    
    TCanvas* cSA4=new TCanvas("cSA4","Fraction of tracks with clusters in layers");
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
    
 // create output file order pdf and root canvases
    
    
    // histos nel file .pdf nell'ordine voluto
    pdfFileNames+=" Vertex_trend.pdf";
    pdfFileNames+=" TrackPointsMI_trend.pdf";
    pdfFileNames+=" Frac_track_SA_trend.pdf";
    pdfFileNames+=" NoFast_trend.pdf";
    pdfFileNames+=" RunEvents_trend.pdf";
    pdfFileNames+=" SDDmodulesON_trend.pdf";
    pdfFileNames+=" SDD_SSD_drift_charge_trend.pdf";
    pdfFileNames+=" SSD_BadStripsFrac_trend.pdf";
    pdfFileNames+=" SSD_chargeratio_trend.pdf";
    pdfFileNames+=" TPCTOFMatch_trend.pdf";
    pdfFileNames+=" Pileup_trend.pdf";
    pdfFileNames+=" PileupVtx_trend.pdf";
    pdfFileNames+=" Pixel_trend.pdf";
    pdfFileNames+=" ITSsa_trend.pdf";
    pdfFileNames+=" ITSTPC_trend.pdf";
    pdfFileNames+=" tracks_ratio_trend.pdf";
    pdfFileNames+=" meanclu_SA_trend.pdf";
    pdfFileNames+=" dedx4_3_SA_trend.pdf";
    pdfFileNames+=" piPt_SA_trend.pdf";

    
    // merge the pdf files
  TString command("gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=merged");
  command=command+"ITS_trend_2016.pdf "+pdfFileNames;
  gSystem->Exec(command.Data());
  printf(" Merging the pdf file:  %s \n",command.Data());
  delete [] myIndex;
  delete [] noRuns;

    // TString nameO=ntupleFileName;
    // nameO+="_File.root";
    TFile *f2=new TFile("nameO.root","RECREATE");
    
    cVertexDisto->Write();
    cMI->Write();
    cSA->Write();
    c5->Write();
    cev->Write();
    c22->Write();
    cfrac->Write();
    c2->Write();
    c8->Write();
    c7->Write();
    cpt02->Write();
    cPileUp->Write();
    cpu->Write();
    cPixel->Write();
    c->Write();
    c2a->Write();
    c3->Write();
    cSAb->Write();
    cSA2->Write();
    cSA3->Write();
    cSA4->Write();
    f2->Close();
    
    TString NameFileOut=ntupleFileName;
    NameFileOut+="_Canvas.root";
    
    TString command1("rm -fr");
    command1=command1+" "+NameFileOut;
    cout<<"\n\n\nCommand 1: "<<command1<<endl;
    gSystem->Exec(command1.Data());
    TString command2("hadd ");
    command2=command2+" "+NameFileOut;
    command2=command2+" nameO.root File_tracksITS.root";
    gSystem->Exec(command2.Data());
    printf(" Merging output:  %s \n",command2.Data());

}

//____________________________________________________________________________
/*void PlotITSSA(TFile *fil, Int_t run1, Int_t run2){
  Double_t Lowbin[3]={0.1,0.5,0.9};
  Double_t Upbin[3]={0.2,0.6,1};
  gROOT->SetStyle("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetFillColor(0);
  gStyle->SetTextFont(32);
  
  TNtuple* nt = (TNtuple*)fil->Get("ntITSsa");

  Float_t ITSA[3];
  Float_t TPIT[3];
  Float_t RAT[3];
  Float_t run;
    Float_t NcluITSpSA,errNcluITSpSA,dedx4_3,errdedx4_3,PtpionpSA,errPtpionpSA;
    Float_t NclupSA0,errNclupSA0,NclupSA1,errNclupSA1,NclupSA2,errNclupSA2,NclupSA3,errNclupSA3,NclupSA4,errNclupSA4,NclupSA5,errNclupSA5;
  nt->SetBranchAddress("run",&run);
  //  nt->SetBranchAddress("NITSsaPtBin0",&ITSA[0]); //Pb-Pb
  // nt->SetBranchAddress("NITSsaPtBin1",&ITSA[1]); //Pb-Pb
  // nt->SetBranchAddress("NITSsaPtBin2",&ITSA[2]); //Pb-Pb
  nt->SetBranchAddress("NITSpureSAPtBin0",&ITSA[0]);
  nt->SetBranchAddress("NITSpureSAPtBin1",&ITSA[1]);
  nt->SetBranchAddress("NITSpureSAPtBin2",&ITSA[2]);
  nt->SetBranchAddress("NITSTPCPtBin0",&TPIT[0]);
  nt->SetBranchAddress("NITSTPCPtBin1",&TPIT[1]);
  nt->SetBranchAddress("NITSTPCPtBin2",&TPIT[2]);
  nt->SetBranchAddress("ratioPtBin0",&RAT[0]);
  nt->SetBranchAddress("ratioPtBin1",&RAT[1]);
  nt->SetBranchAddress("ratioPtBin2",&RAT[2]);
    nt->SetBranchAddress("NcluITSpSA",&NcluITSpSA);
    nt->SetBranchAddress("errNcluITSpSA",&errNcluITSpSA);
    nt->SetBranchAddress("dedx4_3",&dedx4_3);
    nt->SetBranchAddress("errdedx4_3",&errdedx4_3);
    nt->SetBranchAddress("PtpionpSA",&PtpionpSA);
    nt->SetBranchAddress("errPtpionpSA",&errPtpionpSA);
    nt->SetBranchAddress("NclupSA0",&NclupSA0);
    nt->SetBranchAddress("errNclupSA0",&errNclupSA0);
    nt->SetBranchAddress("NclupSA1",&NclupSA1);
    nt->SetBranchAddress("errNclupSA1",&errNclupSA1);
    nt->SetBranchAddress("NclupSA2",&NclupSA2);
    nt->SetBranchAddress("errNclupSA2",&errNclupSA2);
    nt->SetBranchAddress("NclupSA3",&NclupSA3);
    nt->SetBranchAddress("errNclupSA3",&errNclupSA3);
    nt->SetBranchAddress("NclupSA4",&NclupSA4);
    nt->SetBranchAddress("errNclupSA4",&errNclupSA4);
    nt->SetBranchAddress("NclupSA5",&NclupSA5);
    nt->SetBranchAddress("errNclupSA5",&errNclupSA5);
    cout << "SPD1 " << NclupSA0 << " errore " << errNclupSA0 << endl;

  Int_t nr=nt->GetEntries();
  Int_t *myIndex = new Int_t [nr];
  Int_t *noRuns = new Int_t [nr];
  for(Int_t i=0; i<nr;i++){
    nt->GetEvent(i);
    Int_t intrun = static_cast<Int_t>(run+0.01);
    noRuns[i]=intrun;
  }
  printf("\n ======== PROCESSING ITS SA NTUPLE \n");
  Int_t kRunsToPlot = RunsToPlot(run1,run2,nr,noRuns,myIndex);
  TH1F *h0=new TH1F("h0","h0",kRunsToPlot,-0.5,kRunsToPlot-0.5);
  TH1F *h1=new TH1F("h1","h1",kRunsToPlot,-0.5,kRunsToPlot-0.5);
  TH1F *h2=new TH1F("h2","h2",kRunsToPlot,-0.5,kRunsToPlot-0.5);
  TH1F *h3=new TH1F("h3","h4",kRunsToPlot,-0.5,kRunsToPlot-0.5);
  TH1F *h4=new TH1F("h4","h5",kRunsToPlot,-0.5,kRunsToPlot-0.5);
  TH1F *h5=new TH1F("h5","h5",kRunsToPlot,-0.5,kRunsToPlot-0.5);
  TH1F *h6=new TH1F("h6","h6",kRunsToPlot,-0.5,kRunsToPlot-0.5);
  TH1F *h7=new TH1F("h7","h7",kRunsToPlot,-0.5,kRunsToPlot-0.5);
  TH1F *h8=new TH1F("h8","h8",kRunsToPlot,-0.5,kRunsToPlot-0.5);
    TH1F *h9=new TH1F("h9","h9",kRunsToPlot,-0.5,kRunsToPlot-0.5);
    TH1F *h10=new TH1F("h10","h10",kRunsToPlot,-0.5,kRunsToPlot-0.5);
    TH1F *h11=new TH1F("h11","h11",kRunsToPlot,-0.5,kRunsToPlot-0.5);
    TH1F *h12=new TH1F("h12","h12",kRunsToPlot,-0.5,kRunsToPlot-0.5);
    TH1F *h13=new TH1F("h13","h13",kRunsToPlot,-0.5,kRunsToPlot-0.5);
    TH1F *h14=new TH1F("h14","h14",kRunsToPlot,-0.5,kRunsToPlot-0.5);
    TH1F *h15=new TH1F("h15","h15",kRunsToPlot,-0.5,kRunsToPlot-0.5);
    TH1F *h16=new TH1F("h16","h16",kRunsToPlot,-0.5,kRunsToPlot-0.5);
    TH1F *h17=new TH1F("h17","h17",kRunsToPlot,-0.5,kRunsToPlot-0.5);

    
  for(Int_t iev=0;iev<kRunsToPlot;iev++){
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
      h9->SetBinContent(iev+1,NcluITSpSA);
      h9->SetBinError(iev+1,errNcluITSpSA);
      h10->SetBinContent(iev+1,dedx4_3);
      h10->SetBinError(iev+1,errdedx4_3);
      h11->SetBinContent(iev+1,PtpionpSA);
      h11->SetBinError(iev+1,errPtpionpSA);
      h12->SetBinContent(iev+1,NclupSA0);
      h12->SetBinError(iev+1,errNclupSA0);
      h13->SetBinContent(iev+1,NclupSA1);
      h13->SetBinError(iev+1,errNclupSA1);
      h14->SetBinContent(iev+1,NclupSA2);
      h14->SetBinError(iev+1,errNclupSA2);
      h15->SetBinContent(iev+1,NclupSA3);
      h15->SetBinError(iev+1,errNclupSA3);
      h16->SetBinContent(iev+1,NclupSA4);
      h16->SetBinError(iev+1,errNclupSA4);
      h17->SetBinContent(iev+1,NclupSA5);
      h17->SetBinError(iev+1,errNclupSA5);

    h0->GetXaxis()->SetBinLabel(iev+1,Form("%.0f",run));
    h1->GetXaxis()->SetBinLabel(iev+1,Form("%.0f",run));
    h2->GetXaxis()->SetBinLabel(iev+1,Form("%.0f",run));
    h3->GetXaxis()->SetBinLabel(iev+1,Form("%.0f",run));
    h4->GetXaxis()->SetBinLabel(iev+1,Form("%.0f",run));
    h5->GetXaxis()->SetBinLabel(iev+1,Form("%.0f",run));
    h6->GetXaxis()->SetBinLabel(iev+1,Form("%.0f",run));
    h7->GetXaxis()->SetBinLabel(iev+1,Form("%.0f",run));
    h8->GetXaxis()->SetBinLabel(iev+1,Form("%.0f",run));
      h9->GetXaxis()->SetBinLabel(iev+1,Form("%d",(Int_t)run));
      h10->GetXaxis()->SetBinLabel(iev+1,Form("%d",(Int_t)run));
      h11->GetXaxis()->SetBinLabel(iev+1,Form("%d",(Int_t)run));
      h12->GetXaxis()->SetBinLabel(iev+1,Form("%d",(Int_t)run));
      h13->GetXaxis()->SetBinLabel(iev+1,Form("%d",(Int_t)run));
      h14->GetXaxis()->SetBinLabel(iev+1,Form("%d",(Int_t)run));
      h15->GetXaxis()->SetBinLabel(iev+1,Form("%d",(Int_t)run));
      h16->GetXaxis()->SetBinLabel(iev+1,Form("%d",(Int_t)run));
      h17->GetXaxis()->SetBinLabel(iev+1,Form("%d",(Int_t)run));

    //    Printf("%f   %f   %f",ITSA[0],ITSA[1],ITSA[2]);
  }
  h0->Print("all");
  // h0->GetYaxis()->SetTitle("ITSsa tracks");
  // h0->GetXaxis()->SetTitle("run");
  // h1->GetYaxis()->SetTitle("ITSsa tracks");
  // h1->GetXaxis()->SetTitle("run");
  // h2->GetYaxis()->SetTitle("ITSsa tracks");
  // h2->GetXaxis()->SetTitle("run");
  h0->GetYaxis()->SetTitle("ITSpureSA tracks");
  h0->GetXaxis()->SetTitle("run");
  h1->GetYaxis()->SetTitle("ITSpureSA tracks");
  h1->GetXaxis()->SetTitle("run");
  h2->GetYaxis()->SetTitle("ITSpureSA tracks");
  h2->GetXaxis()->SetTitle("run");
  h3->GetYaxis()->SetTitle("ITS+TPC tracks");
  h3->GetXaxis()->SetTitle("run");
  h4->GetYaxis()->SetTitle("ITS+TPC tracks");
  h4->GetXaxis()->SetTitle("run");
  h5->GetYaxis()->SetTitle("ITS+TPC tracks");
  h5->GetXaxis()->SetTitle("run");
  // h6->GetYaxis()->SetTitle("(TPC+ITS)/ITS");
  // h6->GetXaxis()->SetTitle("run");
  // h7->GetYaxis()->SetTitle("(TPC+ITS)/ITS");
  // h7->GetXaxis()->SetTitle("run");
  // h8->GetYaxis()->SetTitle("(TPC+ITS)/ITS");
  // h8->GetXaxis()->SetTitle("run");
  h6->GetYaxis()->SetTitle("(TPC+ITS)+ITSsa/ITSpureSA");
  h6->GetXaxis()->SetTitle("run");
  h7->GetYaxis()->SetTitle("(TPC+ITS)+ITSsa/ITSpureSA");
  h7->GetXaxis()->SetTitle("run");
  h8->GetYaxis()->SetTitle("(TPC+ITS)+ITSsa/ITSpureSA");
  h8->GetXaxis()->SetTitle("run");
    h9->GetXaxis()->SetTitle("run");
    h9->GetYaxis()->SetTitle("mean cluster number (N>3)");
    h10->GetXaxis()->SetTitle("run");
    h10->GetYaxis()->SetTitle("N(dedx4clu)/N(dedx3clu)");
    h11->GetXaxis()->SetTitle("run");
    h11->GetYaxis()->SetTitle("mean pion Pt (GeV/c)");
    h12->GetXaxis()->SetTitle("run");
    h12->GetYaxis()->SetTitle("Fraction of pureSA tracks with cluster in ITS layers");
    h13->GetXaxis()->SetTitle("run");
    h14->GetXaxis()->SetTitle("run");
    h15->GetXaxis()->SetTitle("run");
    h16->GetXaxis()->SetTitle("run");
    h17->GetXaxis()->SetTitle("run");
    
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

  TCanvas *c=new TCanvas("ITS pure SA tracks","ITS pure SA tracks");
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
      pdfFileNames+=" ITSsa_trend.pdf";
    
  TCanvas *c2=new TCanvas("ITS+TPC tracks","ITS+TPC tracks");
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
 
  TCanvas *c3=new TCanvas("(ITSTPC+ITSsa)/ITSpureSA","(ITSTPC+ITSsa)/ITSpureSA");
  c3->SetGridy();
  h8->Draw("p");
  h7->Draw("psame");
  h6->Draw("psame");
   c3->BuildLegend();
  TLatex* ti3=new TLatex(0.10,0.8,"(ITSTPC+ITSsa)/ITSpureSA ");
  ti3->SetNDC();
  ti3->Draw();
   c3->SaveAs("tracks_ratio_trend.pdf");
      pdfFileNames+=" tracks_ratio_trend.pdf";

    TString name1="File_tracksITS.root";
    
    TFile *f=new TFile(name1,"RECREATE");
    f->cd();
    c->Write();
    c2->Write();
    c3->Write();
    f->Close();
    
    TCanvas* cSAb=new TCanvas("cSAb","Mean cluster number (N>3) ITSsa");
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

    TCanvas* cSA2=new TCanvas("cSA2","dedx(4clu)/dedx(3clu) ITSsa");
    h10->SetMarkerStyle(20);
    h10->SetMarkerColor(3);
    h10->Draw();
    tl->Draw();
    cSA2->SaveAs("dedx4_3_SA_trend.pdf");
//    pdfFileNames+=" dedx4_3_SA_trend.pdf";
 
    TCanvas* cSA3=new TCanvas("cSA3","mean pion Pt ITSsa");
    h11->SetMarkerStyle(20);
    h11->SetMarkerColor(4);
    h11->Draw();
    tl->Draw();
    cSA3->SaveAs("piPt_SA_trend.pdf");
//    pdfFileNames+=" piPt_SA_trend.pdf";

    TCanvas* cSA4=new TCanvas("cSA4","Fraction of tracks with clusters in layers");
    h12->GetYaxis()->SetRange(0.,1.);
    h12->SetMarkerStyle(24);
    h12->SetMarkerColor(kGray+1);
    h12->SetLineColor(kGray+1);
    h12->Draw();
    tl->Draw();
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
    pdfFileNames+=" Frac_track_SA_trend.pdf";

}

 */


//____________________________________________________________________________
void AliITSQAtrend_pp_2016(TString runListFile,TString ntupleFileName){

  TGrid::Connect("alien://");
// myfile<<"vediamo se funziona"<<endl;

  
  //-----------SDD
  
  const Int_t nVariables=51;
  TNtuple* ntsdd=new TNtuple("ntsdd","SDD trending","nrun:nEvents:nEventsTriggered:fracTrackWithClu1:errfracTrackWithClu1:fracTrackWithClu2:errfracTrackWithClu2:fracTrackWithClu3:errfracTrackWithClu3:fracTrackWithClu4:errfracTrackWithClu4:fracTrackWithClu5:errfracTrackWithClu5:fracTrackWithClu6:errfracTrackWithClu6:fracEvWithSDD:errfracEvWithSDD:meanTrPts3:errmeanTrPts3:meanTrPts4:errmeanTrPts4:minDrTime:errminDrTime:meanDrTime:errmeanDrTime:fracExtra:errfracExtra:meandEdxLay3:errmeandEdxLay3:meandEdxLay4:errmeandEdxLay4:meandEdxTB0:errmeandEdxTB0:meandEdxTB5:errmeandEdxTB5:MPVdEdxLay3:errMPVdEdxLay3:MPVdEdxLay4:errMPVdEdxLay4:MPVdEdxTB0:errMPVdEdxTB0:MPVdEdxTB5:errMPVdEdxTB5:nMod95:nMod80:nMod60:nModEmpty:fracDead3:errfracDead3:fracDead4:errfracDead4");
  Float_t xnt[nVariables];
  
  //--------------SSD
    
  const Int_t nVariablesSSD=14+8;
  TNtuple* ntssd=new TNtuple("ntssd","SSD trending","nrun:meandEdxLay5:errmeandEdxLay5:meandEdxLay6:errmeandEdxLay6:MPVdEdxLay5:errMPVdEdxLay5:MPVdEdxLay6:errMPVdEdxLay6:ChargeRatioL5:errChargeratioL5:ChargeRatioL6:errChargeratioL6:moduleOff:FracBadn5:errFracBadn5:FracBadp5:errFracBadp5:FracBadn6:errFracBadn6:FracBadp6:errFracBadp6");
  Float_t xntSSD[nVariablesSSD];
  
  //----Matching

  const Int_t nVariablesMatching=60+24;
  TNtuple* ntmatching=new TNtuple("ntmatching","Matching Efficiency","nrun:FracSPD1:errFracSPD1:FracSPD2:errFracSPD2:Eff6Pt02:errEff6Pt02:Eff6Pt1:errEff6Pt1:Eff6Pt10:errEff6Pt10:Eff5Pt02:errEff5Pt02:Eff5Pt1:errEff5Pt1:Eff5Pt10:errEff5Pt10:Eff4Pt02:errEff4Pt02:Eff4Pt1:errEff4Pt1:Eff4Pt10:errEff4Pt10:Eff3Pt02:errEff3Pt02:Eff3Pt1:errEff3Pt1:Eff3Pt10:errEff3Pt10:Eff2Pt02:errEff2Pt02:Eff2Pt1:errEff2Pt1:Eff2Pt10:errEff2Pt10:EffSPDPt02:errEffSPDPt02:EffSPDPt1:errEffSPDPt1:EffSPDPt10:errEffSPDPt10:EffoneSPDPt02:errEffoneSPDPt02:EffoneSPDPt1:errEffoneSPDPt1:EffoneSPDPt10:errEffoneSPDPt10:EffTOTPt02:errEffTOTPt02:EffTOTPt1:errEffTOTPt1:EffTOTPt10:errEffTOTPt10:FracTrackMI1:errFracTrackMI1:FracTrackMI2:errFracTrackMI2:FracTrackMI3:errFracTrackMI3:FracTrackMI4:errFracTrackMI4:FracTrackMI5:errFracTrackMI5:FracTrackMI6:errFracTrackMI6:FracTrackSA1:errFracTrackSA1:FracTrackSA2:errFracTrackSA2:FracTrackSA3:errFracTrackSA3:FracTrackSA4:errFracTrackSA4:FracTrackSA5:errFracTrackSA5:FracTrackSA6:errFracTrackSA6");

  Float_t xntMatching[nVariablesMatching];
  
  //--------------------------------
  //----Matching TOF

  const Int_t nVariablesMatchingTOF=60;

  TNtuple* ntmatchingTOF=new TNtuple("ntmatchingTOF","MatchingTOF Efficiency","nrun:FracSPD1:errFracSPD1:FracSPD2:errFracSPD2:Eff6Pt02:errEff6Pt02:Eff6Pt1:errEff6Pt1:Eff6Pt10:errEff6Pt10:Eff5Pt02:errEff5Pt02:Eff5Pt1:errEff5Pt1:Eff5Pt10:errEff5Pt10:Eff4Pt02:errEff4Pt02:Eff4Pt1:errEff4Pt1:Eff4Pt10:errEff4Pt10:Eff3Pt02:errEff3Pt02:Eff3Pt1:errEff3Pt1:Eff3Pt10:errEff3Pt10:Eff2Pt02:errEff2Pt02:Eff2Pt1:errEff2Pt1:Eff2Pt10:errEff2Pt10:EffSPDPt02:errEffSPDPt02:EffSPDPt1:errEffSPDPt1:EffSPDPt10:errEffSPDPt10:EffoneSPDPt02:errEffoneSPDPt02:EffoneSPDPt1:errEffoneSPDPt1:EffoneSPDPt10:errEffoneSPDPt10:EffTOTPt02:errEffTOTPt02:EffTOTPt1:errEffTOTPt1:EffTOTPt10:errEffTOTPt10");

  Float_t xntMatchingTOF[nVariablesMatchingTOF];

  //--------------------------------

  //----QA Vertex

  const Int_t nVariablesVertex=29;
  TNtuple* ntvertex=new TNtuple("ntvertex","QA Vertex","nrun:VxTRK:errVxTRK:sigmaVxTRK:errsigmaVxTRK:VyTRK:errVyTRK:sigmaVyTRK:errsigmaVyTRK:VzTRK:errVzTRK:sigmaVzTRK:errsigmaVzTRK:VxSPD:errVxSPD:sigmaVxSPD:errsigmaVxSPD:VySPD:errVySPD:sigmaVySPD:errsigmaVySPD:VzSPD:errVzSPD:sigmaVzSPD:errsigmaVzSPD:pileupSPD:errpileupSPD");
      
  Float_t xntVertex[nVariablesVertex];
  
  //--------------------------------

  //-----  Tracking ITS SA
  const Int_t nVariablesSA=13+18;
  Float_t xntSA[nVariablesSA];
  TNtuple *ntSA=new TNtuple("ntITSsa","ntITSsa","run:NITSTPCPtBin0:NITSTPCPtBin1:NITSTPCPtBin2:NITSsaPtBin0:NITSsaPtBin1:NITSsaPtBin2:NITSpureSAPtBin0:NITSpureSAPtBin1:NITSpureSAPtBin2:ratioPtBin0:ratioPtBin1:ratioPtBin2:NcluITSpSA:errNcluITSpSA:dedx4_3:errdedx4_3:PtpionpSA:errPtpionpSA:NclupSA0:errNclupSA0:NclupSA1:errNclupSA1:NclupSA2:errNclupSA2:NclupSA3:errNclupSA3:NclupSA4:errNclupSA4:NclupSA5:errNclupSA5");

  //---------------------------------

//-----  PileUp from SPD vertices
    const Int_t nVariablesPU=11;
    Float_t xntPU[nVariablesPU];
    TNtuple *ntPU=new TNtuple("ntPileUp","ntPileUp","run:npilvtx:errnpilvtx:ntrklpil:errntrklpil:ntrklnopil:errntrklnopil:ncl1pil:errncl1pil:ncl1nopil:errncl1nopil");
    
    //---------------------------------

  
  TBits* readRun=new TBits(999999);
  readRun->ResetAllBits();
  //    if(!useExternalList){
  if(!gSystem->Exec(Form("ls -l %s > /dev/null 2>&1",ntupleFileName.Data()))){
    TFile* oldfil=new TFile(ntupleFileName.Data());
    
    TNtuple* ntmp=(TNtuple*)oldfil->Get("ntsdd");
    
    TNtuple* ntmpSSD=(TNtuple*)oldfil->Get("ntssd");

    TNtuple* ntmpMatching=(TNtuple*)oldfil->Get("ntmatching");

    TNtuple* ntmpMatchingTOF=(TNtuple*)oldfil->Get("ntmatchingTOF");

    TNtuple* ntmpVertex=(TNtuple*)oldfil->Get("ntvertex");

      TNtuple* ntmpSA=(TNtuple*)oldfil->Get("ntITSsa");

      TNtuple* ntmpPU=(TNtuple*)oldfil->Get("ntPileUp");

    
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
    //---------MatchingTOF---------

    Bool_t isOKMatchingTOF=kFALSE;
    if(ntmpMatchingTOF){
      if(ntmpMatchingTOF->GetNvar()==ntmatchingTOF->GetNvar()){
	isOKMatchingTOF=kTRUE;
	TObjArray* arr1matchingTOF=(TObjArray*)ntmatchingTOF->GetListOfBranches();
	TObjArray* arr2matchingTOF=(TObjArray*)ntmpMatchingTOF->GetListOfBranches();
	for(Int_t iV=0; iV<ntmpMatchingTOF->GetNvar(); iV++){
	  TString vnam1=arr1matchingTOF->At(iV)->GetName();
	  TString vnam2=arr2matchingTOF->At(iV)->GetName();
	  if(vnam1!=vnam2) isOKMatchingTOF=kFALSE;
	  ntmpMatchingTOF->SetBranchAddress(vnam2.Data(),&xntMatchingTOF[iV]);
	}
	if(isOKMatchingTOF){
	  for(Int_t nE=0; nE<ntmpMatchingTOF->GetEntries(); nE++){
	    ntmpMatchingTOF->GetEvent(nE);
	    Int_t theRun=(Int_t)(xntMatchingTOF[0]+0.0001);
	    readRun->SetBitNumber(theRun);
	    ntmatchingTOF->Fill(xntMatchingTOF);
	  }
	}
      }
    }
    if(!isOKMatchingTOF){

      printf("\n\nNtuple MatchingTOF in local file not OK -> will be recreated\n\n");

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

      //-----------------------
      
      //-----------PileUp ----------------
      Bool_t isOKPU = kFALSE;
      if(ntmpPU){
          if(ntmpPU->GetNvar()==ntPU->GetNvar()){
              isOKPU = kTRUE;
              TObjArray* arr1PU=(TObjArray*)ntPU->GetListOfBranches();
              TObjArray* arr2PU=(TObjArray*)ntmpPU->GetListOfBranches();
              for(Int_t iV=0; iV<ntmpPU->GetNvar(); iV++){
                  TString vnam1=arr1PU->At(iV)->GetName();
                  TString vnam2=arr2PU->At(iV)->GetName();
                  if(vnam1!=vnam2) isOKPU=kFALSE;
                  ntmpPU->SetBranchAddress(vnam2.Data(),&xntPU[iV]);
              }
              if(isOKPU){
                  for(Int_t nE=0; nE<ntmpPU->GetEntries(); nE++){
                      ntmpPU->GetEvent(nE);
                      Int_t theRun=(Int_t)(xntSA[0]+0.0001);
                      readRun->SetBitNumber(theRun);
                      ntPU->Fill(xntPU);
                  }
              }
          }
      }
      
    oldfil->Close();
    delete oldfil;
  }

#define MAX_LINES 300
#define MAX_LINE_LEN 255
  
  char strings[MAX_LINES][MAX_LINE_LEN];
  ifstream in(runListFile.Data());
  int j = 0;
  Int_t nrun=0;
  Int_t runNumb[MAX_LINES];
  Int_t fillNumb[MAX_LINES];
  Bool_t goout = kFALSE;
  while ( in ) {
    in.getline(strings[j], MAX_LINE_LEN);
    TString aux(strings[j]);
    TString auxrun(strings[j]);
    TString auxfill(strings[j]);
    Int_t lentrail=0;
    Int_t lenfill=0;
    if(aux.Contains("LHC11h/")){
      lentrail = 27;
    }
    else if(aux.Contains("LHC11h_2/")){
      lentrail = 29;
    }
    else if(aux.Contains("LHC12a17b/")){
      lentrail = 26;
    }
    else if(aux.Contains("LHC11c/")){
      lentrail = 27;
      lenfill=30;
    }    else if(aux.Contains("LHC12b/")){
      lentrail = 27;
      lenfill=30;
    }
    else if(aux.Contains("LHC12c/")  || aux.Contains("LHC12d/") ){
      lentrail = 27;
    }
    else if(aux.Contains("LHC12a/")){
      lentrail = 27;
      lenfill=30;
    }
  else if(aux.Contains("LHC12i/")){
      lentrail = 27;
      lenfill=30;
    }
    else if(aux.Contains("LHC12f/")){
      lentrail = 27;
      lenfill=30;
    }
    else if(aux.Contains("g/germain/")){
        lentrail = 50;
        lenfill=30;
    }
    else if(aux.Contains("LHC16k/")){
        lentrail = 27;
        lenfill=30;
    }
    else {
      if(!aux.IsNull())printf("Unrecognised path name %s \n",aux.Data());
      goout = kTRUE;
    }
    if(goout)break;
    if(aux.Length()<lentrail)continue;
    auxrun=aux.Remove(0,lentrail);
    auxfill=aux.Remove(0,36);
    // cout<<" auxfill= "<<auxfill<<endl;
    //   aux=aux.Remove(6,aux.Length());
    //    aux=aux.Remove();
    runNumb[j]=atoi(auxrun.Data());
    fillNumb[j]=atoi(auxfill.Data());
    printf("%d ) - fill %d - path %s \n",runNumb[j],fillNumb[j],strings[j]);
    j++;
    nrun++;
  }

  printf("\n *******************   Loop on runs *********** \n");
  Int_t filenotfound=0;
  for(Int_t jru=0;jru<nrun;jru++) {
    printf("jru=%d - run number= %d \n",jru,runNumb[jru]);
    Int_t iRun=runNumb[jru];
    if(readRun->TestBitNumber(iRun))printf("Run %d - already processed\n",iRun);
    if(readRun->TestBitNumber(iRun))continue;
    //cout << "Value from file is " <<t << endl;
    
    //    printf("%s\n",strings[jru]);
  
      myfile<<endl<<"Processing run "<<runNumb[jru]<<endl;

    if(!gGrid||!gGrid->IsConnected()) {
      printf("gGrid not found! exit macro\n");
      return;
    }
    
    TFile *f=TFile::Open(Form("alien://%s",strings[jru])); 
    if(!f) {
      printf("File not found, continue with next one\n");
    myfile<<"File not found, continue with next one\n";

      filenotfound++;
      continue;
    }
      
    //------- check ITS presence in QAresults.root file

      TDirectoryFile* df=(TDirectoryFile*)f->Get("Vertex_Performance");
      if(!df) continue;
      
      TList* l=(TList*)df->Get("cOutputVtxESD");
      if(!l) continue;
      
      TH1F* h=(TH1F*)l->FindObject("fhSPDVertexX"); // number of processed events
      if(h->GetEntries()==0) {
          printf("ITS not present in run %d, continue with next one\n", iRun);
          continue;
      }
      
    //-------SDD
      myfile << "SDD " << endl;
    FillSDDntuple(f,ntsdd,iRun,xnt);
  

    //-------SSD
      myfile << "SSD " << endl;

    FillSSDntuple(f,ntssd,iRun,xntSSD);

 
    //--------------matching
      myfile << "Matching " << endl;

    FillMatchntuple(f,ntmatching,iRun,xntMatching);

    //--------------matching TOF
      myfile << "Matching TOF " << endl;

    FillMatchntupleTOF(f,ntmatchingTOF,iRun,xntMatchingTOF);

    //------------- Vertex
      myfile << "Vertex " << endl;

      FillVTXntuple(f,ntvertex,iRun,xntVertex);
	
  
    //---------------------------

    //--------------   ITS SA ---------------
    cout<<"ITS - SA"<<endl;
    myfile << "ITS - SA " << endl;

    //    cout<<f<<" "<<ntSA<<" "<<iRun<<endl;
    FillITSSAntuple(f,ntSA,iRun);

      //--------------   PileUp ---------------
      //
      myfile << "PileUp " << endl;

      FillPileUpntuple(f,ntPU,iRun,xntPU);

  } // loop on runs

  printf("%d runs skipped because QA file not found\n",filenotfound);
    myfile<<filenotfound<<" skipped because no QA file"<<endl;

  TFile* outfil=new TFile(ntupleFileName.Data(),"recreate");
  outfil->cd();
  ntsdd->Write();
  ntssd->Write();
  ntmatching->Write();
  ntmatchingTOF->Write();
  ntvertex->Write();
  ntSA->Write();
  ntPU->Write();
  outfil->Close();
  delete outfil;
  delete ntsdd;
  delete ntssd;
  delete ntmatching;
  delete ntvertex;
  delete ntSA;
  delete ntPU;
    
}

//____________________________________________________________________________
void FillITSSAntuple(TFile* f,TNtuple* nt, Int_t nrun){
  static const Int_t nVariables=31;
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
  myfile<<"number of events with spd vertex "<<noEvents<<endl;
  if(noEvents<1.)noEvents=1.;   // protection to avoid division by zero
  hPtTPCITS=(TH1F*)cOutput->FindObject("hPtTPCITS");
  myfile<<"number of ITSTPC Tracks "<<hPtTPCITS->GetEntries()<<endl;
  hPtITSsa=(TH1F*)cOutput->FindObject("hPtITSsa");
  
  myfile<<"number of ITSsa Tracks "<<hPtITSsa->GetEntries()<<endl;

  hPtITSpureSA=(TH1F*)cOutput->FindObject("hPtITSpureSA");
  
  myfile<<"number of ITSpureSA Tracks "<<hPtITSpureSA->GetEntries()<<endl;


  for(Int_t ibin=0;ibin<=2;ibin++){
    NTPCITS[ibin]=hPtTPCITS->Integral(hPtTPCITS->FindBin(Lowbin[ibin]),hPtTPCITS->FindBin(Upbin[ibin]))/noEvents;
    NITSsa[ibin]=hPtITSsa->Integral(hPtITSsa->FindBin(Lowbin[ibin]),hPtITSsa->FindBin(Upbin[ibin]))/noEvents;
    NITSpureSA[ibin]=hPtITSpureSA->Integral(hPtITSpureSA->FindBin(Lowbin[ibin]),hPtITSpureSA->FindBin(Upbin[ibin]))/noEvents;
    //    if(NTPCITS[ibin]!=0 && NITSsa[ibin]!=0)Ratio[ibin]=NTPCITS[ibin]/NITSsa[ibin];
    Double_t totaltrks=NTPCITS[ibin]+NITSsa[ibin];
  if(totaltrks!=0 && NITSpureSA[ibin]!=0 )Ratio[ibin]=totaltrks/NITSpureSA[ibin];
    else Ratio[ibin]=0;
  }
    

 // Elena
    TH1F* hNcluITSpSA =(TH1F*)cOutput->FindObject("hNcluITSpureSA"); // numero cluster se N>3
    TH1F* hdedx4 =(TH1F*)cOutput->FindObject("hdedxvsP4clsITSpureSA");
    TH1F* hdedx3 =(TH1F*)cOutput->FindObject("hdedxvsP3clsITSpureSA");
    TH1F* hPtpSA =(TH1F*)cOutput->FindObject("hPtITSpureSA");  // <-- nuovo, controllare che esista
    TH1F* hPtpionpSA =(TH1F*)cOutput->FindObject("hPtITSpureSAPion");  //
    TH2F* hNcluPion2d = (TH2F*)cOutput->FindObject("hCluInLayITSpureSAPion");
    TH2F* hNclu2d = (TH2F*)cOutput->FindObject("hCluInLayVsPtITSpureSA"); // <-- nuovo, controllare che esista
    TH1D* hNcluPion = hNcluPion2d->ProjectionY();
    TH1D* hNclu;
    if(hNclu2d)hNclu = hNclu2d->ProjectionY();
    
    Float_t e4, e3;        // quiqui
    e4 = (Float_t)hdedx4->GetEntries();
    e3 = (Float_t)hdedx3->GetEntries();
    Float_t r43 = 0.;
    Float_t erre43=0.;
    if(hdedx3->GetEntries()>0){
        r43 =e4/e3;
        erre43=TMath::Sqrt(e4/(e3*e3)*(1+r43));
    }

    
    Float_t fracTpi[6]={0.,0.,0.,0.,0.,0.};
    Float_t efracTpi[6]={0.,0.,0.,0.,0.,0.};
    if(hNclu && hNclu->GetBinContent(1)>0){
//    if(hNcluPion->GetBinContent(1)>0){
        for(Int_t iLay=0; iLay<6; iLay++){
//            fracTpi[iLay]=hNcluPion->GetBinContent(iLay+2)/hNcluPion->GetBinContent(1);
//            efracTpi[iLay]=TMath::Sqrt(fracTpi[iLay]*(1-fracTpi[iLay])/hNcluPion->GetBinContent(1));
            fracTpi[iLay]=hNclu->GetBinContent(iLay+2)/hNclu->GetBinContent(1);
            efracTpi[iLay]=TMath::Sqrt(fracTpi[iLay]*(1-fracTpi[iLay])/hNclu->GetBinContent(1));
        }
    }
// Elena
    
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


// Elena
    xnt[index++]=hNcluITSpSA->GetMean();
    xnt[index++]=hNcluITSpSA->GetMeanError();
    xnt[index++]=r43;
    xnt[index++]=erre43;
//    xnt[index++]=hPtpionpSA->GetMean();
//    xnt[index++]=hPtpionpSA->GetMeanError();
    if(hPtpSA) xnt[index++]=hPtpSA->GetMean();
    else xnt[index++]=0.;
    if(hPtpSA)xnt[index++]=hPtpSA->GetMeanError();
    else xnt[index++]=0.;
    xnt[index++]=fracTpi[0];
    xnt[index++]=efracTpi[0];
    xnt[index++]=fracTpi[1];
    xnt[index++]=efracTpi[1];
    xnt[index++]=fracTpi[2];
    xnt[index++]=efracTpi[2];
    xnt[index++]=fracTpi[3];
    xnt[index++]=efracTpi[3];
    xnt[index++]=fracTpi[4];
    xnt[index++]=efracTpi[4];
    xnt[index++]=fracTpi[5];
    xnt[index++]=efracTpi[5];
     // Elena

  nt->Fill(xnt);
}

//_____________________________________________________________________________
void FillSDDntuple(TFile* f,TNtuple* nt, Int_t iRun, Float_t *xnt){
  TDirectoryFile* df=(TDirectoryFile*)f->Get("SDD_Performance");
  if(!df){
    printf("Run %d SDD_Performance MISSING -> Exit\n",iRun);
 myfile<<" SDD_Performance MISSING -> Exit\n";

    return;
  } 
   
    
  TList* l=(TList*)df->Get("coutputRP");
  if(!l){
    printf("Run %d coutputRP TList MISSING -> Exit\n",iRun);
myfile<<" coutputRP TList MISSING -> Exit\n";

    return;
  }  
    


  cout<<"SDD - QA"<<endl;

    TH1F* hnev=(TH1F*)l->FindObject("hNEvents"); // number of processed events
    
    if(hnev->GetEntries()==0){
        
        printf("Run %d hNEvents EMPTY -> Continue\n",iRun);
        myfile<<"Run  hnev EMPTY -> Continue\n";
//        return;
        
    }
    
    Float_t fracE=0.;
    Float_t efracE=0.;
    if(hnev->GetBinContent(2)>0){
             fracE=hnev->GetBinContent(4)/hnev->GetBinContent(2); // bin 4 = with SDD; bin 2 = selected events
             efracE=TMath::Sqrt(fracE*(1-fracE)/hnev->GetBinContent(2));
        }
    
    
    TH1F* hcllay=(TH1F*)l->FindObject("hCluInLay"); // questo histo potrebbe essere sostituito da un altro
                                                  // o da altri due presi da ITS_Performance
    
  if(hcllay->GetEntries()==0){

    printf("Run %d hcllay EMPTY -> Continue\n",iRun);
    myfile<<"Run  hcllay EMPTY -> Continue\n";

// return;
      
  }

  Float_t fracT[6]={0.,0.,0.,0.,0.,0.};
  Float_t efracT[6]={0.,0.,0.,0.,0.,0.};
  if(hcllay->GetBinContent(1)>0){
    for(Int_t iLay=0; iLay<6; iLay++){
      fracT[iLay]=hcllay->GetBinContent(iLay+2)/hcllay->GetBinContent(1);
      efracT[iLay]=TMath::Sqrt(fracT[iLay]*(1-fracT[iLay])/hcllay->GetBinContent(1));
        }
    }
  
    
  TH1F* hmodT=(TH1F*)l->FindObject("hTPMod");
    
  if(hmodT->GetEntries()==0){
    printf("Run %d hmodT EMPTY -> Continue\n",iRun);
    myfile<<"Run hmodT EMPTY -> Continue\n";
      // return;

  }
    
  TH1F* hgamod=(TH1F*)l->FindObject("hGAMod");
    
  if(hgamod->GetEntries()==0){
    printf("Run %d hgamod EMPTY -> Continue\n",iRun);
    myfile<<"Run  hgamod EMPTY -> Continue\n";

      // return;
  }

  Int_t bestMod=0;
  //
  Int_t deadMod3=0, deadMod4=0;
  //

  for(Int_t iMod=0; iMod<260;iMod++){
    Int_t gda=(Int_t)hgamod->GetBinContent(iMod+1);
    if(gda>bestMod) bestMod=gda;
      if(gda == 0) {
        if(iMod<84) deadMod3+=1;
        else deadMod4+=1;
      }
  }
    
  Int_t nChunks=1;
  if(bestMod>512){
    nChunks=(Int_t)(bestMod/512.+0.5);
  }
  hgamod->Scale(1./nChunks);
    Float_t fracDead3 = (Float_t)deadMod3/84;
    Float_t fracDead4 = (Float_t)deadMod4/176;
    Float_t errfracDead3 = TMath::Sqrt(fracDead3*(1-fracDead3)/84);         //
    Float_t errfracDead4 = TMath::Sqrt(fracDead4*(1-fracDead4)/176);
    
  TH1F* hev=(TH1F*)l->FindObject("hNEvents");

  if(hev->GetEntries()==0){
    printf("Run %d hev EMPTY -> Continue\n",iRun);
      myfile<<"Run  hev EMPTY -> Continue\n";

      // return;
  }
    
  Int_t nTotEvents=hev->GetBinContent(1);
  Int_t nTrigEvents=hev->GetBinContent(2);
  Int_t nEvents=nTotEvents;
  printf("Run %d Number of Events = %d Triggered=%d\n",iRun,nTotEvents,nTrigEvents);
myfile<<"Number of events "<<nTotEvents<<" Triggered "<<nTrigEvents<<endl;

  if(nTrigEvents>0){ 
    nEvents=nTrigEvents;
  }
  if(nTotEvents==0) return;
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
      
  if(hapmod->GetEntries()==0){
    printf("Run %d hapmod EMPTY -> Continue\n",iRun);
      // return;
  }


  TH1F* hgpmod=(TH1F*)l->FindObject("hGoodPmod");
  if(hgpmod->GetEntries()==0){
    printf("Run %d hgpmod EMPTY -> Continue\n",iRun);
myfile<<"hgpmod EMPTY"<<endl;

      // return;
  }
      
  //     TH1F* hmpmod=(TH1F*)l->FindObject("hMissPmod");
  TH1F* hbrmod=(TH1F*)l->FindObject("hBadRegmod");
  if(hbrmod->GetEntries()==0){
    printf("Run %d hbrmod EMPTY -> Continue\n",iRun);
    myfile<<"hbrmod EMPTY"<<endl;

      // return;
  }

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

  if(htimT->GetEntries()==0){
    printf("Run %d htimT EMPTY -> Continue\n",iRun);
    myfile<<"htimt empty"<<endl;

      // return;
  }

  TH1F* htimTe=(TH1F*)l->FindObject("hDrTimTPExtra");
      
  if(htimTe->GetEntries()==0){
    printf("Run %d htimTe EMPTY -> Continue\n",iRun);
    myfile<<"htimTe empty"<<endl;

      // return;
  }

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

  if(hdedxmod->GetEntries()==0){
    printf("Run %d hdedxmod EMPTY -> Continue\n",iRun);
    myfile<<"hdedxmod empty"<<endl;

      // return;
  }

  TH1D* hdedxLay3=hdedxmod->ProjectionY("hdedxLay3",1,84);
  TH1D* hdedxLay4=hdedxmod->ProjectionY("hdedxLay4",85,260);
      
  TH1F* hSigTim0=(TH1F*)l->FindObject("hSigTimeInt0");
  if(hSigTim0->GetEntries()==0){
    printf("Run %d hSigTim0 EMPTY -> Continue\n",iRun);
    myfile<<"hSigTim0 empty"<<endl;

      // return;
  }

  TH1F* hSigTim5=(TH1F*)l->FindObject("hSigTimeInt5");
  if(hSigTim5->GetEntries()==0){
    printf("Run %d hSigTim5  EMPTY -> Continue\n",iRun);
    myfile<<"hSigTim5 empty"<<endl;

      // return;
  }
  //Fitting the same distributions in order to have the MPV
  TF1 *lfunLay3 = new TF1("LangausFunLay3",LangausFun,50.,300.,4); 
  lfunLay3->SetParameter(0,5.);
  lfunLay3->SetParameter(1,80.);  // questo resta a 80 se il fit non viene fatto!!
  lfunLay3->SetParameter(2,hdedxLay3->GetEntries()/10.);
  lfunLay3->SetParameter(3,10.);
  lfunLay3->SetParLimits(3,0.,20);
//    hdedxLay3->Fit(lfunLay3,"NQLR");
    hdedxLay3->Fit(lfunLay3,"NQLR");
  TF1 *lfunLay4 = new TF1("LangausFunLay4",LangausFun,50.,300.,4);
  lfunLay4->SetParameter(0,5.);
  lfunLay4->SetParameter(1,80.); // questo resta a 80 se il fit non viene fatto!!
  lfunLay4->SetParameter(2,hdedxLay4->GetEntries()/10.);
  lfunLay4->SetParameter(3,10.);
  lfunLay4->SetParLimits(3,0.,20);
  hdedxLay4->Fit(lfunLay4,"NQLR");
  TF1 *lfunTim0 = new TF1("LangausFunTim0",LangausFun,50.,300.,4); 
  lfunTim0->SetParameter(0,5.);
  lfunTim0->SetParameter(1,80.); // questo resta a 80 se il fit non viene fatto!!
  lfunTim0->SetParameter(2,hSigTim0->GetEntries()/10.);
  lfunTim0->SetParameter(3,10.);
  lfunTim0->SetParLimits(3,0.,20);
  hSigTim0->Fit(lfunTim0,"NQLR");
  TF1 *lfunTim5 = new TF1("LangausFunTim5",LangausFun,50.,300.,4); 
  lfunTim5->SetParameter(0,5.);
  lfunTim5->SetParameter(1,80.); // questo resta a 80 se il fit non viene fatto!!
  lfunTim5->SetParameter(2,hSigTim5->GetEntries()/10.);
  lfunTim5->SetParameter(3,10.);
  lfunTim5->SetParLimits(3,0.,20);
  hSigTim5->Fit(lfunTim5,"NQLR");
     
  Int_t index=0;
  xnt[index++]=iRun;
    xnt[index++]=nTotEvents;
    xnt[index++]=nTrigEvents;
    // inserire qui riempimento numero eventi in ntupla
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
//
    xnt[index++]=fracE;
    xnt[index++]=efracE;
//
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


    if(hdedxLay3->GetEntries()>0) {xnt[index++]=lfunLay3->GetParameter(1); xnt[index++]=lfunLay3->GetParError(1);}
    else {xnt[index++]=0.; xnt[index++]=0.;}
    if(hdedxLay4->GetEntries()>0) {xnt[index++]=lfunLay4->GetParameter(1); xnt[index++]=lfunLay4->GetParError(1);}
    else {xnt[index++]=0.; xnt[index++]=0.;}

    if(hSigTim0->GetEntries()>0) {xnt[index++]=lfunTim0->GetParameter(1); xnt[index++]=lfunTim0->GetParError(1);}
    else {xnt[index++]=0.; xnt[index++]=0.;}
    if(hSigTim5->GetEntries()>0) {xnt[index++]=lfunTim5->GetParameter(1); xnt[index++]=lfunTim5->GetParError(1);}
    else {xnt[index++]=0.; xnt[index++]=0.;}
    
  xnt[index++]=(Float_t)nBelow95;
  xnt[index++]=(Float_t)nBelow80;
  xnt[index++]=(Float_t)nBelow60;
  xnt[index++]=(Float_t)nZeroP;
    xnt[index++]=1.-fracDead3;
    xnt[index++]=errfracDead3;
    xnt[index++]=1.-fracDead4;
    xnt[index++]=errfracDead4;
  nt->Fill(xnt);
  

}

//_____________________________________________________________________________
void FillSSDntuple(TFile* f,TNtuple* ntssd, Int_t iRun, Float_t *xntSSD){

  cout<<"SSD - QA"<<endl;

  TDirectoryFile* dfSSD=(TDirectoryFile*)f->Get("PWGPPdEdxSSDQA");
  if(!dfSSD){
    printf("Run %d SSD_Performance MISSING -> Exit\n",iRun);
    myfile<<"SDD Performance missing"<<endl;

    return;
  }
      
  TList* lSSD=(TList*)dfSSD->Get("SSDdEdxQA");
  if(!dfSSD){
    printf("Run %d coutputRP TList MISSING -> Exit\n",iRun);
    myfile<<"coutputRP empty"<<endl;

    return;
  }
  //
      
  TH2F* QAchargeRatio=(TH2F*)lSSD->FindObject("QAChargeRatioSA");
      
  if(QAchargeRatio->GetEntries()==0){
    printf("Run %d QAchargeRatio EMPTY -> Continue\n",iRun);
    myfile<<"QA charge ratio empty"<<endl;

      // return;
  }

  TH2F* QAcharge=(TH2F*)lSSD->FindObject("QAChargeSA");
      
  if(QAcharge->GetEntries()==0){
    printf("Run %d QAcharge EMPTY -> Continue\n",iRun);
    myfile<<"QA charge  empty"<<endl;

      // return;
  }
 
  Int_t biny = QAcharge->GetXaxis()->FindBin(747);
  Int_t maxy = QAcharge->GetXaxis()->GetXmax();
  //     Int_t miny = QAcharge->GetXaxis()->GetXmin();
    
  Int_t  contEmpty=0;
  Int_t  contFull=0;
    
  TH1D *hChargeL5=QAcharge->ProjectionY("hChargeL5",0,biny);
  TH1D *hChargeL6=QAcharge->ProjectionY("hChargeL6",biny,maxy);
    
  //     cout<<  hChargeL5->GetMean()<< " " <<  hChargeL5->GetRMS()<<endl;
  //     cout<<  hChargeL6->GetMean()<< " " <<  hChargeL6->GetRMS()<<endl;
    
  TH1D *hChargeRatioL5=QAchargeRatio->ProjectionY("hChargeRatioL5",0,biny);
  TH1D *hChargeRatioL6=QAchargeRatio->ProjectionY("hChargeRatioL6",biny,maxy);
    
  //     cout<<  hChargeRatioL5->GetMean()<< " " <<  hChargeRatioL5->GetRMS()<<endl;
  //     cout<<  hChargeRatioL6->GetMean()<< " " <<  hChargeRatioL6->GetRMS()<<endl;
    
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

  //     cout<<"contFull: " <<contFull<<" contEmpty: "<<contEmpty<<endl;
  //     cout<<hChargeL5->GetMean()<<endl;

  //Fitting dE/dx Distr in order to have the MPV
  TF1 *lfunLay5 = new TF1("LangausFunLay5",LangausFun,50.,300.,4); 
  lfunLay5->SetParameter(0,5.);
  lfunLay5->SetParameter(1,80.);
  lfunLay5->SetParameter(2,hChargeL5->GetEntries()/10.);
  lfunLay5->SetParameter(3,10.);
  lfunLay5->SetParLimits(3,0.,20);
  hChargeL5->Fit(lfunLay5,"NQLR");
  TF1 *lfunLay6 = new TF1("LangausFunLay6",LangausFun,50.,300.,4); 
  lfunLay6->SetParameter(0,5.);
  lfunLay6->SetParameter(1,80.);
  lfunLay6->SetParameter(2,hChargeL6->GetEntries()/10.);
  lfunLay6->SetParameter(3,10.);
  lfunLay6->SetParLimits(3,0.,20);
  hChargeL6->Fit(lfunLay6,"NQLR");
    
 // Elena: fraction of bad_n/p_strips for layer 5 and 6
    
    TH1F* bad_p=(TH1F*)lSSD->FindObject("Bad-p-strips");
    TH1F* bad_n=(TH1F*)lSSD->FindObject("Bad-n-strips");
 // find number of merged subjobs
    Int_t max_p=0;
    Float_t n_subjobs=0;
//    Int_t max_n=0; // not used: same maximum for both histos
    for(Int_t i=1;i<1699;i++) if(bad_p->GetBinContent(i)>max_p)max_p=bad_p->GetBinContent(i);
    n_subjobs=max_p/768.;
// find the number of bad n/p strips per layer
    bad_p->Scale(1/n_subjobs);
    bad_n->Scale(1/n_subjobs);
    Int_t bad_n5=0, bad_p5=0, bad_n6=0, bad_p6=0;
    for(Int_t j=1; j<749; j++) {bad_n5=bad_n5+bad_n->GetBinContent(j); bad_p5=bad_p5+bad_p->GetBinContent(j);}
    for(Int_t j=749; j<1699; j++) {bad_n6=bad_n6+bad_n->GetBinContent(j); bad_p6=bad_p6+bad_p->GetBinContent(j);}
//
  Int_t indexSSD=0;
  xntSSD[indexSSD++]=iRun;
  xntSSD[indexSSD++]=(Float_t)hChargeL5->GetMean();
  xntSSD[indexSSD++]=(Float_t)hChargeL5->GetMeanError();
  xntSSD[indexSSD++]=(Float_t)hChargeL6->GetMean();
  xntSSD[indexSSD++]=(Float_t)hChargeL6->GetMeanError();
    
    if(hChargeL5->GetEntries()>0) {xntSSD[indexSSD++]=lfunLay5->GetParameter(1); xntSSD[indexSSD++]=lfunLay5->GetParError(1);}
    else {xntSSD[indexSSD++]=0.; xntSSD[indexSSD++]=0.;}

    if(hChargeL6->GetEntries()>0) {xntSSD[indexSSD++]=lfunLay6->GetParameter(1); xntSSD[indexSSD++]=lfunLay6->GetParError(1);}
    else {xntSSD[indexSSD++]=0.; xntSSD[indexSSD++]=0.;}

    xntSSD[indexSSD++]=(Float_t)hChargeRatioL5->GetMean();
  xntSSD[indexSSD++]=(Float_t)hChargeRatioL5->GetMeanError();
  xntSSD[indexSSD++]=(Float_t)hChargeRatioL6->GetMean();
  xntSSD[indexSSD++]=(Float_t)hChargeRatioL6->GetMeanError();
  xntSSD[indexSSD++]=(Float_t)contEmpty;
// fractions of bad n/p strips per layer
    xntSSD[indexSSD++]=(Float_t)bad_n5/574464.;
    xntSSD[indexSSD++]=TMath::Sqrt((Float_t)bad_n5)/574464.;
    xntSSD[indexSSD++]=(Float_t)bad_p5/574464.;
    xntSSD[indexSSD++]=TMath::Sqrt((Float_t)bad_p5)/574464.;
    xntSSD[indexSSD++]=(Float_t)bad_n6/729600.;
    xntSSD[indexSSD++]=TMath::Sqrt((Float_t)bad_n6)/729600.;
    xntSSD[indexSSD++]=(Float_t)bad_p6/729600.;
    xntSSD[indexSSD++]=TMath::Sqrt((Float_t)bad_p6)/729600.;
    
  ntssd->Fill(xntSSD);
}

//_____________________________________________________________________________
void FillMatchntuple(TFile* f,TNtuple* ntmatching, Int_t iRun, Float_t *xntMatching){
    cout<<"Tracking TPC"<<endl;

    TDirectoryFile *dirMatch=(TDirectoryFile*)f->GetDirectory("ITS_Performance");
    TList *list=NULL;
    TList *listSPD=NULL;
      
    if(dirMatch) {
      //	list = (TList*)dirMatch->Get("cOutputITS_3500_10000"); //LHC11h
      //	if(!list)list = (TList*)dirMatch->Get("cOutputITS"); // LHC11e
      list = (TList*)dirMatch->Get("cOutputITS"); // LHC12e
    }
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

    if(hFiredChip->GetEntries()==0){
      printf("Run %d hFiredChip EMPTY -> Continue\n",iRun);
      myfile<<"hFiredChip  empty"<<endl;
        // return;
    }

    Int_t nHSsInner=0,nHSsOuter=0;
    for(Int_t i=0;i<400;i++) if(hFiredChip->GetBinContent(i)>0) nHSsInner++;
    for(Int_t i=400;i<1200;i++) if(hFiredChip->GetBinContent(i)>0) nHSsOuter++;
    nHSsInner = (Int_t)(nHSsInner/10);
    nHSsOuter = (Int_t)(nHSsOuter/10);
    // hnHSsSPD->SetBinContent(1,nHSsInner);
    // hnHSsSPD->SetBinContent(2,nHSsOuter);
      
    ioValues[0]=(Float_t)nHSsInner/40.;
    ioValues[1]=(Float_t)nHSsOuter/80.;

    //
    TH1F* hclmapMI=(TH1F*)list->FindObject("fHistClusterMapITSMI"); // tracce ricostruite
    TH1F* hclMI=(TH1F*)list->FindObject("fHistNclsITSMI"); //

    TH1F* hclmapSA=(TH1F*)list->FindObject("fHistClusterMapITSSA"); // tracce SA
    TH1F* hclSA=(TH1F*)list->FindObject("fHistNclsITSSA"); //


    if(hclmapMI->GetEntries()==0){
        printf("Run %d hclmapMI EMPTY -> Continue\n",iRun);
        myfile<<"Run  hclmapMI EMPTY -> Continue\n";
        // return;
    }


    if(hclMI->GetEntries()==0){
        printf("Run %d hclMI EMPTY -> Continue\n",iRun);
        myfile<<"Run  hclMI EMPTY -> Continue\n";
        // return;
    }


    if(hclmapSA->GetEntries()==0){
        
        printf("Run %d hclmapSA EMPTY -> Continue\n",iRun);
        myfile<<"Run  hclmapSA EMPTY -> Continue\n";
        // return;
    }


    if(hclSA->GetEntries()==0){
        printf("Run %d hclSA EMPTY -> Continue\n",iRun);
        myfile<<"Run  hclSA EMPTY -> Continue\n";
        // return;
    }

    
    Float_t fracTMI[6]={0.,0.,0.,0.,0.,0.};
    Float_t efracTMI[6]={0.,0.,0.,0.,0.,0.};
    if(hclMI->GetEntries()>0){
        for(Int_t iLay=0; iLay<6; iLay++){
//            fracTMI[iLay]=hclmapMI->GetBinContent(iLay+1)/hclMI->GetEntries();
//            efracTMI[iLay]=TMath::Sqrt(fracTMI[iLay]*(1-fracTMI[iLay])/hclMI->GetEntries());
              fracTMI[iLay]=hclmapMI->GetBinContent(iLay+1)/(hclMI->GetEntries()-hclMI->GetBinContent(1));
              efracTMI[iLay]=TMath::Sqrt(fracTMI[iLay]*(1-fracTMI[iLay])/(hclMI->GetEntries()-hclMI->GetBinContent(1)));
        }
    }

    Float_t fracTSA[6]={0.,0.,0.,0.,0.,0.};
    Float_t efracTSA[6]={0.,0.,0.,0.,0.,0.};
    if(hclSA->GetEntries()>0){
        for(Int_t iLay=0; iLay<6; iLay++){
            fracTSA[iLay]=hclmapSA->GetBinContent(iLay+1)/hclSA->GetEntries();
            efracTSA[iLay]=TMath::Sqrt(fracTSA[iLay]*(1-fracTSA[iLay])/hclSA->GetEntries());
        }
    }

    //
    
    TH1F *fHistPtTPCInAcc = (TH1F*)list->FindObject("fHistPtTPCInAcc");
    Int_t check1=0;

 if(fHistPtTPCInAcc->GetEntries()==0){check1=1;
      printf("Run %dfHistPtTPCInAcc  EMPTY -> Continue\n",iRun);
     // return;
    }

    TH1F *fHistPtITSMI6InAcc = (TH1F*)list->FindObject("fHistPtITSMI6InAcc");

    if(fHistPtITSMI6InAcc->GetEntries()==0){check1=1;
      printf("Run %d fHistPtITSMI6InAcc EMPTY -> Continue\n",iRun);
        // return;
    }

    TH1F *fHistPtITSMI5InAcc = (TH1F*)list->FindObject("fHistPtITSMI5InAcc");
    if(fHistPtITSMI5InAcc->GetEntries()==0){check1=1;
      printf("Run %d fHistPtITSMI5InAcc EMPTY -> Continue\n",iRun);
        // return;
    }
      
    TH1F *fHistPtITSMI4InAcc = (TH1F*)list->FindObject("fHistPtITSMI4InAcc");
    if(fHistPtITSMI5InAcc->GetEntries()==0){check1=1;
      printf("Run %d fHistPtITSMI4InAcc EMPTY -> Continue\n",iRun);
        // return;
    }
      
    TH1F *fHistPtITSMI3InAcc = (TH1F*)list->FindObject("fHistPtITSMI3InAcc");
    if(fHistPtITSMI3InAcc->GetEntries()==0){
      printf("Run %d fHistPtITSMI3InAcc EMPTY -> Continue\n",iRun);
        // return;
    }
    TH1F *fHistPtITSMI2InAcc = (TH1F*)list->FindObject("fHistPtITSMI2InAcc");
      
    if(fHistPtITSMI2InAcc->GetEntries()==0){check1=1;
      printf("Run %d fHistPtITSMI2InAcc EMPTY -> Continue\n",iRun);
        // return;
    }
    TH1F *fHistPtITSMISPDInAcc = (TH1F*)list->FindObject("fHistPtITSMISPDInAcc");
    if(fHistPtITSMISPDInAcc->GetEntries()==0){check1=1;
      printf("Run %d fHistPtITSMISPDInAcc EMPTY -> Continue\n",iRun);
        // return;
    }
      
    TH1F *fHistPtITSMIoneSPDInAcc = (TH1F*)list->FindObject("fHistPtITSMIoneSPDInAcc");
      
    if(fHistPtITSMIoneSPDInAcc->GetEntries()==0){check1=1;
      printf("Run %d fHistPtITSMIoneSPDInAcc  EMPTY -> Continue\n",iRun);
        // return;
    }
    if(check1==1)myfile<<"Something is missing for marching"<<endl;

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
    ioValues[8]=fHistPtITSMI4InAcc->GetBinContent(ptbin);
    ioErrors[8]=fHistPtITSMI4InAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMI4InAcc->FindBin(1.001);
    ioValues[9]=fHistPtITSMI4InAcc->GetBinContent(ptbin);
    ioErrors[9]=fHistPtITSMI4InAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMI4InAcc->FindBin(10.001);
    ioValues[10]=fHistPtITSMI4InAcc->GetBinContent(ptbin);
    ioErrors[10]=fHistPtITSMI4InAcc->GetBinError(ptbin);

    fHistPtITSMI3InAcc->Divide(fHistPtITSMI3InAcc,fHistPtTPCInAcc,1,1,"B");

    ptbin=fHistPtITSMI3InAcc->FindBin(0.201);
    ioValues[11]=fHistPtITSMI3InAcc->GetBinContent(ptbin);
    ioErrors[11]=fHistPtITSMI3InAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMI3InAcc->FindBin(1.001);
    ioValues[12]=fHistPtITSMI3InAcc->GetBinContent(ptbin);
    ioErrors[12]=fHistPtITSMI3InAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMI3InAcc->FindBin(10.001);
    ioValues[13]=fHistPtITSMI3InAcc->GetBinContent(ptbin);
    ioErrors[13]=fHistPtITSMI3InAcc->GetBinError(ptbin);

    fHistPtITSMI2InAcc->Divide(fHistPtITSMI2InAcc,fHistPtTPCInAcc,1,1,"B");

    ptbin=fHistPtITSMI2InAcc->FindBin(0.201);
    ioValues[14]=fHistPtITSMI2InAcc->GetBinContent(ptbin);
    ioErrors[14]=fHistPtITSMI2InAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMI2InAcc->FindBin(1.001);
    ioValues[15]=fHistPtITSMI2InAcc->GetBinContent(ptbin);
    ioErrors[15]=fHistPtITSMI2InAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMI2InAcc->FindBin(10.001);
    ioValues[16]=fHistPtITSMI2InAcc->GetBinContent(ptbin);
    ioErrors[16]=fHistPtITSMI2InAcc->GetBinError(ptbin);

    fHistPtITSMISPDInAcc->Divide(fHistPtITSMISPDInAcc,fHistPtTPCInAcc,1,1,"B");
    ptbin=fHistPtITSMISPDInAcc->FindBin(0.201);
    ioValues[17]=fHistPtITSMISPDInAcc->GetBinContent(ptbin);
    ioErrors[17]=fHistPtITSMISPDInAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMISPDInAcc->FindBin(1.001);
    ioValues[18]=fHistPtITSMISPDInAcc->GetBinContent(ptbin);
    ioErrors[18]=fHistPtITSMISPDInAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMISPDInAcc->FindBin(10.001);
    ioValues[19]=fHistPtITSMISPDInAcc->GetBinContent(ptbin);
    ioErrors[19]=fHistPtITSMISPDInAcc->GetBinError(ptbin);

    fHistPtITSMIoneSPDInAcc->Divide(fHistPtITSMIoneSPDInAcc,fHistPtTPCInAcc,1,1,"B");

    ptbin=fHistPtITSMIoneSPDInAcc->FindBin(0.201);
    ioValues[20]=fHistPtITSMIoneSPDInAcc->GetBinContent(ptbin);
    ioErrors[20]=fHistPtITSMIoneSPDInAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMIoneSPDInAcc->FindBin(1.001);
    ioValues[21]=fHistPtITSMIoneSPDInAcc->GetBinContent(ptbin);
    ioErrors[21]=fHistPtITSMIoneSPDInAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMIoneSPDInAcc->FindBin(10.001);
    ioValues[22]=fHistPtITSMIoneSPDInAcc->GetBinContent(ptbin);
    ioErrors[22]=fHistPtITSMIoneSPDInAcc->GetBinError(ptbin);

  
    fHistPtITSMIge2InAcc->Divide(fHistPtITSMIge2InAcc,fHistPtTPCInAcc,1,1,"B");
    ptbin=fHistPtITSMIge2InAcc->FindBin(0.201);
    ioValues[23]=fHistPtITSMIge2InAcc->GetBinContent(ptbin);
    ioErrors[23]=fHistPtITSMIge2InAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMIge2InAcc->FindBin(1.001);
    ioValues[24]=fHistPtITSMIge2InAcc->GetBinContent(ptbin);
    ioErrors[24]=fHistPtITSMIge2InAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMIge2InAcc->FindBin(10.001);
    ioValues[25]=fHistPtITSMIge2InAcc->GetBinContent(ptbin);
    ioErrors[25]=fHistPtITSMIge2InAcc->GetBinError(ptbin);

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
//    xntMatching[indexMatching++]=ioValues[26];
//    xntMatching[indexMatching++]=ioErrors[26];
    xntMatching[indexMatching++]=fracTMI[0];
    xntMatching[indexMatching++]=efracTMI[0];
    xntMatching[indexMatching++]=fracTMI[1];
    xntMatching[indexMatching++]=efracTMI[1];
    xntMatching[indexMatching++]=fracTMI[2];
    xntMatching[indexMatching++]=efracTMI[2];
    xntMatching[indexMatching++]=fracTMI[3];
    xntMatching[indexMatching++]=efracTMI[3];
    xntMatching[indexMatching++]=fracTMI[4];
    xntMatching[indexMatching++]=efracTMI[4];
    xntMatching[indexMatching++]=fracTMI[5];
    xntMatching[indexMatching++]=efracTMI[5];
    xntMatching[indexMatching++]=fracTSA[0];
    xntMatching[indexMatching++]=efracTSA[0];
    xntMatching[indexMatching++]=fracTSA[1];
    xntMatching[indexMatching++]=efracTSA[1];
    xntMatching[indexMatching++]=fracTSA[2];
    xntMatching[indexMatching++]=efracTSA[2];
    xntMatching[indexMatching++]=fracTSA[3];
    xntMatching[indexMatching++]=efracTSA[3];
    xntMatching[indexMatching++]=fracTSA[4];
    xntMatching[indexMatching++]=efracTSA[4];
    xntMatching[indexMatching++]=fracTSA[5];
    xntMatching[indexMatching++]=efracTSA[5];
    
  
    ntmatching->Fill(xntMatching);						
}
 
//_________________________________________________________________________

// Matching with TOF

void FillMatchntupleTOF(TFile* f,TNtuple* ntmatchingTOF, Int_t iRun, Float_t *xntMatchingTOF){

  cout<<"Tracking TOF"<<endl;

  //  printf("\n\nMATCHING WITH TOF!!!!!!!!!\n",iRun);

  TDirectoryFile *dirMatch=(TDirectoryFile*)f->GetDirectory("ITS_Performance");
  TList *list=NULL;
  TList *listSPD=NULL;

  if(dirMatch) {
   //	list = (TList*)dirMatch->Get("cOutputITS_3500_10000"); //LHC11h
    //	if(!list)list = (TList*)dirMatch->Get("cOutputITS"); // LHC11e
    list = (TList*)dirMatch->Get("cOutputITS"); // LHC12e
  }
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

  TH1F *fHistPtTPCInAcc = (TH1F*)list->FindObject("fHistPtTPCInAccTOFbc0");
  Int_t check1=0;

   if(fHistPtTPCInAcc->GetEntries()==0){
     check1=1;
     printf("Run %dfHistPtTPCInAccTOFbc0  EMPTY -> Continue\n",iRun);
       // return;
   }
    
  TH1F *fHistPtITSMI6InAcc = (TH1F*)list->FindObject("fHistPtITSMI6InAccTOFbc0");

  if(fHistPtITSMI6InAcc->GetEntries()==0){
     check1=1;
     printf("Run %d fHistPtITSMI6InAccTOFbc0 EMPTY -> Continue\n",iRun);
      // return;
   }

  TH1F *fHistPtITSMI5InAcc = (TH1F*)list->FindObject("fHistPtITSMI5InAccTOFbc0");

   if(fHistPtITSMI5InAcc->GetEntries()==0){
     check1=1;
     printf("Run %d fHistPtITSMI5InAccTOFbc0 EMPTY -> Continue\n",iRun);
       // return;
   }

  TH1F *fHistPtITSMI4InAcc = (TH1F*)list->FindObject("fHistPtITSMI4InAccTOFbc0");

   if(fHistPtITSMI4InAcc->GetEntries()==0){
     check1=1;
     printf("Run %d fHistPtITSMI5InAccTOFbc0 EMPTY -> Continue\n",iRun);
       // return;
   }
      
  TH1F *fHistPtITSMI3InAcc = (TH1F*)list->FindObject("fHistPtITSMI3InAccTOFbc0");

   if(fHistPtITSMI3InAcc->GetEntries()==0){
     printf("Run %d fHistPtITSMI3InAccTOFbc0 EMPTY -> Continue\n",iRun);
       // return;
   }
 
  TH1F *fHistPtITSMI2InAcc = (TH1F*)list->FindObject("fHistPtITSMI2InAccTOFbc0");
  
   if(fHistPtITSMI2InAcc->GetEntries()==0){
     check1=1;
     printf("Run %d fHistPtITSMI2InAccTOFbc0 EMPTY -> Continue\n",iRun);
       // return;
   }
  
  TH1F *fHistPtITSMISPDInAcc = (TH1F*)list->FindObject("fHistPtITSMISPDInAccTOFbc0");
 
   if(fHistPtITSMISPDInAcc->GetEntries()==0){
     check1=1;
     printf("Run %d fHistPtITSMISPDInAccTOFbc0 EMPTY -> Continue\n",iRun);
       // return;
   }

   TH1F *fHistPtITSMIoneSPDInAcc = (TH1F*)list->FindObject("fHistPtITSMIoneSPDInAccTOFbc0");

   if(fHistPtITSMIoneSPDInAcc->GetEntries()==0){
     check1=1;
     printf("Run %d fHistPtITSMIoneSPDInAccTOFbc0  EMPTY -> Continue\n",iRun);
       // return;
   }

  if(check1==1)myfile<<"Something is missing for marching"<<endl;

  TH1F *fHistPtITSMIge2InAcc = (TH1F*)fHistPtITSMI6InAcc->Clone("fHistPtITSMIge2InAccTOFbc0");
  fHistPtITSMIge2InAcc->Add(fHistPtITSMI5InAcc);
  fHistPtITSMIge2InAcc->Add(fHistPtITSMI4InAcc);
  fHistPtITSMIge2InAcc->Add(fHistPtITSMI3InAcc);
  fHistPtITSMIge2InAcc->Add(fHistPtITSMI2InAcc);
 
 fHistPtITSMI6InAcc->Divide(fHistPtITSMI6InAcc,fHistPtTPCInAcc,1,1,"B");
  ptbin=fHistPtITSMI6InAcc->FindBin(0.501);
  ioValues[2]=fHistPtITSMI6InAcc->GetBinContent(ptbin);
  ioErrors[2]=fHistPtITSMI6InAcc->GetBinError(ptbin);
  ptbin=fHistPtITSMI6InAcc->FindBin(1.001);
  ioValues[3]=fHistPtITSMI6InAcc->GetBinContent(ptbin);
  ioErrors[3]=fHistPtITSMI6InAcc->GetBinError(ptbin);
  ptbin=fHistPtITSMI6InAcc->FindBin(10.001);
  ioValues[4]=fHistPtITSMI6InAcc->GetBinContent(ptbin);
  ioErrors[4]=fHistPtITSMI6InAcc->GetBinError(ptbin);

  fHistPtITSMI5InAcc->Divide(fHistPtITSMI5InAcc,fHistPtTPCInAcc,1,1,"B");
  ptbin=fHistPtITSMI5InAcc->FindBin(0.501);
  ioValues[5]=fHistPtITSMI5InAcc->GetBinContent(ptbin);
  ioErrors[5]=fHistPtITSMI5InAcc->GetBinError(ptbin);
  ptbin=fHistPtITSMI5InAcc->FindBin(1.001);
  ioValues[6]=fHistPtITSMI5InAcc->GetBinContent(ptbin);
  ioErrors[6]=fHistPtITSMI5InAcc->GetBinError(ptbin);
  ptbin=fHistPtITSMI5InAcc->FindBin(10.001);
  ioValues[7]=fHistPtITSMI5InAcc->GetBinContent(ptbin);
  ioErrors[7]=fHistPtITSMI5InAcc->GetBinError(ptbin);

  fHistPtITSMI4InAcc->Divide(fHistPtITSMI4InAcc,fHistPtTPCInAcc,1,1,"B");
  ptbin=fHistPtITSMI4InAcc->FindBin(0.501);
  ioValues[8]=fHistPtITSMI4InAcc->GetBinContent(ptbin);
  ioErrors[8]=fHistPtITSMI4InAcc->GetBinError(ptbin);
  ptbin=fHistPtITSMI4InAcc->FindBin(1.001);
  ioValues[9]=fHistPtITSMI4InAcc->GetBinContent(ptbin);
  ioErrors[9]=fHistPtITSMI4InAcc->GetBinError(ptbin);
  ptbin=fHistPtITSMI4InAcc->FindBin(10.001);
  ioValues[10]=fHistPtITSMI4InAcc->GetBinContent(ptbin);
  ioErrors[10]=fHistPtITSMI4InAcc->GetBinError(ptbin);

  fHistPtITSMI3InAcc->Divide(fHistPtITSMI3InAcc,fHistPtTPCInAcc,1,1,"B");
  ptbin=fHistPtITSMI3InAcc->FindBin(0.501);
  ioValues[11]=fHistPtITSMI3InAcc->GetBinContent(ptbin);
  ioErrors[11]=fHistPtITSMI3InAcc->GetBinError(ptbin);
  ptbin=fHistPtITSMI3InAcc->FindBin(1.001);
  ioValues[12]=fHistPtITSMI3InAcc->GetBinContent(ptbin);
  ioErrors[12]=fHistPtITSMI3InAcc->GetBinError(ptbin);
  ptbin=fHistPtITSMI3InAcc->FindBin(10.001);
  ioValues[13]=fHistPtITSMI3InAcc->GetBinContent(ptbin);
  ioErrors[13]=fHistPtITSMI3InAcc->GetBinError(ptbin);

  fHistPtITSMI2InAcc->Divide(fHistPtITSMI2InAcc,fHistPtTPCInAcc,1,1,"B");
  ptbin=fHistPtITSMI2InAcc->FindBin(0.501);
  ioValues[14]=fHistPtITSMI2InAcc->GetBinContent(ptbin);
  ioErrors[14]=fHistPtITSMI2InAcc->GetBinError(ptbin);
  ptbin=fHistPtITSMI2InAcc->FindBin(1.001);
  ioValues[15]=fHistPtITSMI2InAcc->GetBinContent(ptbin);
  ioErrors[15]=fHistPtITSMI2InAcc->GetBinError(ptbin);
  ptbin=fHistPtITSMI2InAcc->FindBin(10.001);
  ioValues[16]=fHistPtITSMI2InAcc->GetBinContent(ptbin);
  ioErrors[16]=fHistPtITSMI2InAcc->GetBinError(ptbin);

  fHistPtITSMISPDInAcc->Divide(fHistPtITSMISPDInAcc,fHistPtTPCInAcc,1,1,"B");
  ptbin=fHistPtITSMISPDInAcc->FindBin(0.501);
  ioValues[17]=fHistPtITSMISPDInAcc->GetBinContent(ptbin);
  ioErrors[17]=fHistPtITSMISPDInAcc->GetBinError(ptbin);
  ptbin=fHistPtITSMISPDInAcc->FindBin(1.001);
  ioValues[18]=fHistPtITSMISPDInAcc->GetBinContent(ptbin);
  ioErrors[18]=fHistPtITSMISPDInAcc->GetBinError(ptbin);
  ptbin=fHistPtITSMISPDInAcc->FindBin(10.001);
  ioValues[19]=fHistPtITSMISPDInAcc->GetBinContent(ptbin);
  ioErrors[19]=fHistPtITSMISPDInAcc->GetBinError(ptbin);

  fHistPtITSMIoneSPDInAcc->Divide(fHistPtITSMIoneSPDInAcc,fHistPtTPCInAcc,1,1,"B");
  ptbin=fHistPtITSMIoneSPDInAcc->FindBin(0.501);
  ioValues[20]=fHistPtITSMIoneSPDInAcc->GetBinContent(ptbin);
  ioErrors[20]=fHistPtITSMIoneSPDInAcc->GetBinError(ptbin);
  ptbin=fHistPtITSMIoneSPDInAcc->FindBin(1.001);
  ioValues[21]=fHistPtITSMIoneSPDInAcc->GetBinContent(ptbin);
  ioErrors[21]=fHistPtITSMIoneSPDInAcc->GetBinError(ptbin);
  ptbin=fHistPtITSMIoneSPDInAcc->FindBin(10.001);
  ioValues[22]=fHistPtITSMIoneSPDInAcc->GetBinContent(ptbin);
  ioErrors[22]=fHistPtITSMIoneSPDInAcc->GetBinError(ptbin);

  fHistPtITSMIge2InAcc->Divide(fHistPtITSMIge2InAcc,fHistPtTPCInAcc,1,1,"B");
  ptbin=fHistPtITSMIge2InAcc->FindBin(0.501);
  ioValues[23]=fHistPtITSMIge2InAcc->GetBinContent(ptbin);
  ioErrors[23]=fHistPtITSMIge2InAcc->GetBinError(ptbin);
  ptbin=fHistPtITSMIge2InAcc->FindBin(1.001);
  ioValues[24]=fHistPtITSMIge2InAcc->GetBinContent(ptbin);
  ioErrors[24]=fHistPtITSMIge2InAcc->GetBinError(ptbin);
  ptbin=fHistPtITSMIge2InAcc->FindBin(10.001);
  ioValues[25]=fHistPtITSMIge2InAcc->GetBinContent(ptbin);
  ioErrors[25]=fHistPtITSMIge2InAcc->GetBinError(ptbin);

  Int_t indexMatchingTOF=0;
  xntMatchingTOF[indexMatchingTOF++]=iRun;
  xntMatchingTOF[indexMatchingTOF++]=ioValues[0];
  xntMatchingTOF[indexMatchingTOF++]=ioErrors[0];
  xntMatchingTOF[indexMatchingTOF++]=ioValues[1];
  xntMatchingTOF[indexMatchingTOF++]=ioErrors[1];
  xntMatchingTOF[indexMatchingTOF++]=ioValues[2];
  xntMatchingTOF[indexMatchingTOF++]=ioErrors[2];
  xntMatchingTOF[indexMatchingTOF++]=ioValues[3];
  xntMatchingTOF[indexMatchingTOF++]=ioErrors[3];
  xntMatchingTOF[indexMatchingTOF++]=ioValues[4];
  xntMatchingTOF[indexMatchingTOF++]=ioErrors[4];
  xntMatchingTOF[indexMatchingTOF++]=ioValues[5];
  xntMatchingTOF[indexMatchingTOF++]=ioErrors[5];
  xntMatchingTOF[indexMatchingTOF++]=ioValues[6];
  xntMatchingTOF[indexMatchingTOF++]=ioErrors[6];
  xntMatchingTOF[indexMatchingTOF++]=ioValues[7];
  xntMatchingTOF[indexMatchingTOF++]=ioErrors[7];
  xntMatchingTOF[indexMatchingTOF++]=ioValues[8];
  xntMatchingTOF[indexMatchingTOF++]=ioErrors[8];
  xntMatchingTOF[indexMatchingTOF++]=ioValues[9];
  xntMatchingTOF[indexMatchingTOF++]=ioErrors[9];
  xntMatchingTOF[indexMatchingTOF++]=ioValues[10];
  xntMatchingTOF[indexMatchingTOF++]=ioErrors[10];
  xntMatchingTOF[indexMatchingTOF++]=ioValues[11];
  xntMatchingTOF[indexMatchingTOF++]=ioErrors[11];
  xntMatchingTOF[indexMatchingTOF++]=ioValues[12];
  xntMatchingTOF[indexMatchingTOF++]=ioErrors[12];
  xntMatchingTOF[indexMatchingTOF++]=ioValues[13];
  xntMatchingTOF[indexMatchingTOF++]=ioErrors[13];
  xntMatchingTOF[indexMatchingTOF++]=ioValues[14];
  xntMatchingTOF[indexMatchingTOF++]=ioErrors[14];
  xntMatchingTOF[indexMatchingTOF++]=ioValues[15];
  xntMatchingTOF[indexMatchingTOF++]=ioErrors[15];
  xntMatchingTOF[indexMatchingTOF++]=ioValues[16];
  xntMatchingTOF[indexMatchingTOF++]=ioErrors[16];
  xntMatchingTOF[indexMatchingTOF++]=ioValues[17];
  xntMatchingTOF[indexMatchingTOF++]=ioErrors[17];
  xntMatchingTOF[indexMatchingTOF++]=ioValues[18];
  xntMatchingTOF[indexMatchingTOF++]=ioErrors[18];
  xntMatchingTOF[indexMatchingTOF++]=ioValues[19];
  xntMatchingTOF[indexMatchingTOF++]=ioErrors[19];
  xntMatchingTOF[indexMatchingTOF++]=ioValues[20];
  xntMatchingTOF[indexMatchingTOF++]=ioErrors[20];
  xntMatchingTOF[indexMatchingTOF++]=ioValues[21];
  xntMatchingTOF[indexMatchingTOF++]=ioErrors[21];
  xntMatchingTOF[indexMatchingTOF++]=ioValues[22];
  xntMatchingTOF[indexMatchingTOF++]=ioErrors[22];
  xntMatchingTOF[indexMatchingTOF++]=ioValues[23];
  xntMatchingTOF[indexMatchingTOF++]=ioErrors[23];
  xntMatchingTOF[indexMatchingTOF++]=ioValues[24];
  xntMatchingTOF[indexMatchingTOF++]=ioErrors[24];
  xntMatchingTOF[indexMatchingTOF++]=ioValues[25];
  xntMatchingTOF[indexMatchingTOF++]=ioErrors[25];
//  xntMatchingTOF[indexMatchingTOF++]=ioValues[26];
//  xntMatchingTOF[indexMatchingTOF++]=ioErrors[26];

  ntmatchingTOF->Fill(xntMatchingTOF);						

}
  
//_____________________________________________________________________________
void FillVTXntuple(TFile* f,TNtuple* ntvertex, Int_t iRun, Float_t *xntVertex){
   cout<<"Primary Vertex"<<endl;

    TDirectoryFile *dirVertex = (TDirectoryFile*)f->Get("Vertex_Performance");
    if(!dirVertex){
      Printf("Vertex directory not found... check!");
      myfile<<"Vertex directory not found"<<endl;

    }
  
    TList *lt = (TList*)dirVertex->Get("cOutputVtxESD");
	
     TH1F *xVtxTRK = (TH1F*)lt->FindObject("fhTRKVertexX");

     // VERTEX Tracks
     if(xVtxTRK->GetEntries()==0){
       printf("Run %d xVtxTRK EMPTY -> Continue\n",iRun);
         // return;
     }
      TH1F *yVtxTRK = (TH1F*)lt->FindObject("fhTRKVertexY");

     if(yVtxTRK->GetEntries()==0){
       printf("Run %d yVtxTRK EMPTY -> Continue\n",iRun);
         // return;
     }

     TH1F *zVtxTRK = (TH1F*)lt->FindObject("fhTRKVertexZ");

     if(zVtxTRK->GetEntries()==0){
       printf("Run %d zVtxTRK EMPTY -> Continue\n",iRun);
         // return;
     }

 TH1F *xVtxSPD = (TH1F*)lt->FindObject("fhSPDVertexX");

    if(xVtxSPD->GetEntries()==0){
      printf("Run %d xVtxSOD EMPTY -> Continue\n",iRun);
      myfile<<"xVtxSPD not found"<<endl;

        // return;
   }

    TH1F *yVtxSPD = (TH1F*)lt->FindObject("fhSPDVertexY");

    if(yVtxSPD->GetEntries()==0){
      printf("Run %d yVtxSPD EMPTY -> Continue\n",iRun);
      myfile<<"yVtxSPD not found"<<endl;

        // return;
    }

    TH1F *zVtxSPD = (TH1F*)lt->FindObject("fhSPDVertexZ");

    if(zVtxSPD->GetEntries()==0){
      printf("Run %d zVtxSPD EMPTY -> Continue\n",iRun);
      myfile<<"zVtxSPD not found"<<endl;

        // return;
    }


    TH1F *zVtxSPDpil = (TH1F*)lt->FindObject("fhSPDVertexZPile");

    if(zVtxSPDpil->GetEntries()==0){
        printf("Run %d zVtxSPDpil EMPTY -> Continue\n",iRun);
        myfile<<"VtxSPDpil not found"<<endl;
        
        // return;
    }

    
    TF1 *fxTRK = new TF1("gausx", "gaus", -1, 1);
    xVtxTRK->Fit("gausx", "NQRL");

    TF1 *fyTRK = new TF1("gausy", "gaus", -1, 1);
    yVtxTRK->Fit("gausy","NQLR");
//     cout<<fyTRK->GetParameter(1)<<endl;
//     cout<<fyTRK->GetParError(1)<<endl;
//     cout<<fyTRK->GetParameter(2)<<endl;
//     cout<<fyTRK->GetParError(2)<<endl;

    TF1 *fzTRK = new TF1("gausz", "gaus", -1, 1);
    zVtxTRK->Fit("gausz","NQRL");
    TF1 *fxSPD = new TF1("gausxSPD", "gaus", -1, 1);
    xVtxSPD->Fit("gausxSPD", "NQRL");

    TF1 *fySPD = new TF1("gausySPD", "gaus", -1, 1);
    yVtxSPD->Fit("gausySPD","NQLR");
//     cout<<fyTRK->GetParameter(1)<<endl;
//     cout<<fyTRK->GetParError(1)<<endl;
//     cout<<fyTRK->GetParameter(2)<<endl;
//     cout<<fyTRK->GetParError(2)<<endl;

    TF1 *fzSPD = new TF1("gauszSPD", "gaus", -1, 1);
    zVtxSPD->Fit("gauszSPD","NQRL");
//     cout<<fzTRK->GetParameter(1)<<endl;
//     cout<<fzTRK->GetParError(1)<<endl;
//     cout<<fzTRK->GetParameter(2)<<endl;
//     cout<<fzTRK->GetParError(2)<<endl;


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
    xntVertex[indexVertex++]=(Float_t)fxSPD->GetParameter(1);
    xntVertex[indexVertex++]=(Float_t)fxSPD->GetParError(1);
    xntVertex[indexVertex++]=(Float_t)fxSPD->GetParameter(2);
    xntVertex[indexVertex++]=(Float_t)fxSPD->GetParError(2);
    xntVertex[indexVertex++]=(Float_t)fySPD->GetParameter(1);
    xntVertex[indexVertex++]=(Float_t)fySPD->GetParError(1);
    xntVertex[indexVertex++]=(Float_t)fySPD->GetParameter(2);
    xntVertex[indexVertex++]=(Float_t)fySPD->GetParError(2);
    xntVertex[indexVertex++]=(Float_t)fzSPD->GetParameter(1);
    xntVertex[indexVertex++]=(Float_t)fzSPD->GetParError(1);
    xntVertex[indexVertex++]=(Float_t)fzSPD->GetParameter(2);
    xntVertex[indexVertex++]=(Float_t)fzSPD->GetParError(2);
    
    xntVertex[indexVertex++]=(Float_t)zVtxSPDpil->GetEntries()/(Float_t)zVtxSPD->GetEntries();
    xntVertex[indexVertex++]=TMath::Sqrt((Float_t)zVtxSPDpil->GetEntries())/(Float_t)zVtxSPD->GetEntries();
    
    ntvertex->Fill(xntVertex);	

}


//_____________________________________________________________________________
void FillPileUpntuple(TFile* f,TNtuple* ntPU, Int_t iRun, Float_t *xntPU){
    cout<<"Pileup from SPD vertices"<<endl;
    
    TDirectoryFile *dirPU = (TDirectoryFile*)f->Get("CheckPileupQA");
    if(!dirPU){
        Printf("Pileup directory not found... check!");
        myfile<<"Pileup directory not found"<<endl;
        
    }
    
    TList *lt = (TList*)dirPU->Get("clistPileupSPDQA");
    
    TH1F *hnPilVtx = (TH1F*)lt->FindObject("hNOfPileupVertSPD");
    
    // Pileup vertices from SPD
    if(hnPilVtx->GetEntries()==0){
        printf("Run %d hnPilVtx EMPTY -> Continue\n",iRun);
        // return;
    }
    TH1F *hnTrklPil = (TH1F*)lt->FindObject("hNtracklPilSPD");
    
    if(hnTrklPil->GetEntries()==0){
        printf("Run %d hnTrklPil EMPTY -> Continue\n",iRun);
        // return;
    }
    
    TH1F *hnTrklNoPil = (TH1F*)lt->FindObject("hNtracklNoPilSPD");
    
    if(hnTrklNoPil->GetEntries()==0){
        printf("Run %d hnTrklNoPil EMPTY -> Continue\n",iRun);
        // return;
    }
    
    TH1F *hnCl1Pil = (TH1F*)lt->FindObject("hNCL1PilSPD");
    
    if(hnCl1Pil->GetEntries()==0){
        printf("Run %d hnCl1Pil EMPTY -> Continue\n",iRun);
        myfile<<"hnCl1Pil not found"<<endl;
        
        // return;
    }
    
    TH1F *hnCl1NoPil = (TH1F*)lt->FindObject("hNCL1NoPilSPD");
    
    if(hnCl1NoPil->GetEntries()==0){
        printf("Run %d hnCl1NoPil EMPTY -> Continue\n",iRun);
        myfile<<"hnCl1NoPil not found"<<endl;
        
        // return;
    }
    
    Float_t meanPilVtx = 0.0, meanPilVtx_n=0;
    Float_t errmeanPilVtx = 0.0;
    Float_t errmeanPilVtx_n = 0.0;
    Float_t errmeanPilVtx_d = 0.0;

    if(hnPilVtx->GetEntries()){
        for(Int_t i=2;i<12;i++){
            meanPilVtx_n = meanPilVtx_n + (i-1)*hnPilVtx->GetBinContent(i);
            meanPilVtx = meanPilVtx_n/hnPilVtx->Integral(2,12);
            errmeanPilVtx_n = meanPilVtx + (i-1)*(i-1)*hnPilVtx->GetBinContent(i);
    }
        errmeanPilVtx_n = TMath::Sqrt(errmeanPilVtx_n);
        errmeanPilVtx_d = TMath::Sqrt(hnPilVtx->Integral(2,12));
        errmeanPilVtx = (errmeanPilVtx_n/meanPilVtx_n)*(errmeanPilVtx_n/meanPilVtx_n)+(errmeanPilVtx_d/hnPilVtx->Integral(2,12))*(errmeanPilVtx_d/hnPilVtx->Integral(2,12));
        errmeanPilVtx = TMath::Sqrt(errmeanPilVtx);
    }
    
    Int_t indexPil=0;
    xntPU[indexPil++]=iRun;
    xntPU[indexPil++]=meanPilVtx;
    xntPU[indexPil++]=errmeanPilVtx;
    xntPU[indexPil++]=hnTrklPil->GetMean();
    xntPU[indexPil++]=hnTrklPil->GetMeanError();
    xntPU[indexPil++]=hnTrklNoPil->GetMean();
    xntPU[indexPil++]=hnTrklNoPil->GetMeanError();
    xntPU[indexPil++]=hnCl1Pil->GetMean();
    xntPU[indexPil++]=hnCl1Pil->GetMeanError();
    xntPU[indexPil++]=hnCl1NoPil->GetMean();
    xntPU[indexPil++]=hnCl1NoPil->GetMeanError();

    ntPU->Fill(xntPU);
    
}


//_____________________________________________________________________________
Double_t LangausFun(Double_t *x, Double_t *par) {
  
  //Fit parameters:
  //par[0]=Width (scale) parameter of Landau density
  //par[1]=Most Probable (MP, location) parameter of Landau density
  //par[2]=Total area (integral -inf to inf, normalization constant)
  //par[3]=Width (sigma) of convoluted Gaussian function
  //
  //In the Landau distribution (represented by the CERNLIB approximation), 
  //the maximum is located at x=-0.22278298 with the location parameter=0.
  //This shift is corrected within this function, so that the actual
  //maximum is identical to the MP parameter.
  
  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location
  
  // Control constants
  Double_t np = 100.0;      // number of convolution steps
  Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
  
  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;
  
  
  // MP shift correction
  mpc = par[1] - mpshift * par[0]; 
  
  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];
  
  step = (xupp-xlow) / np;
  
  // Convolution integral of Landau and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
    
    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }
  
  return (par[2] * step * sum * invsq2pi / par[3]);
}

Bool_t WriteInputTextFileFromMonalisaListOfRuns(TString outtxtfilename,Int_t* listofrunsfromMonalisa,Int_t nruns,TString pathbeforRunN,TString pathafterRunN){

  // e.g. Int_t listofrunsfromMonalisa[3]={1111,2222,3333,4444}. Copy&paste the list of run given in Monalisa! http://alimonitor.cern.ch/raw/raw_details.jsp?timesel=0&filter_jobtype=LHC+period+LHC12c+-+CPass1+%28reconstruction%29#

  ofstream outfile;
  outfile.open(outtxtfilename.Data());
  cout<<"Writing..."<<endl;
  for(Int_t i=0;i<nruns;i++){
    outfile<<pathbeforRunN.Data()<<"000"<<listofrunsfromMonalisa[i]<<"/"<<pathafterRunN.Data()<<endl;
    cout<<pathbeforRunN.Data()<<"000"<<listofrunsfromMonalisa[i]<<"/"<<pathafterRunN.Data()<<endl;
  }

  cout<<"Done"<<endl;
  outfile.close();
  return kTRUE;
}

//______________________________________________________________________
Int_t RunsToPlot(Int_t run1, Int_t run2,Int_t nr,Int_t *noRuns, Int_t *myIndex){
  // Sort entries according to run number in the chosen range
  Int_t kRunsToPlot=0;
  printf("Processing runs from %d up to %d\n",run1,run2);
  for(Int_t i=0; i<nr;i++){
    printf("Run %d\n",noRuns[i]);
    if(noRuns[i]>=run1 && noRuns[i]<=run2){
      printf("Accepting run number %d in position %d\n",noRuns[i],kRunsToPlot);
      kRunsToPlot++;
    }
    else { 
      printf("Rejecting run number %d - out of range\n",noRuns[i]);
      noRuns[i]=run2+10; 
    }
  }
  TMath::Sort(nr,noRuns,myIndex,kFALSE);
  printf("Total number of runs accepted for display %d\n",kRunsToPlot);
  if(kRunsToPlot==0)return 0;
  for(Int_t i=0;i<kRunsToPlot;i++)printf("Position %d ) Run: %d\n",i,noRuns[myIndex[i]]);
  return kRunsToPlot;
}