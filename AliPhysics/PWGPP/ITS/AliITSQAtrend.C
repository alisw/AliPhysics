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
void MakePlot(Int_t run1=-1,Int_t run2=999999,TString ntupleFileName="TrendingITS_2016.root");
//void PlotITSSA(TFile *fil,Int_t run1, Int_t run2);
void FillITSSAntuple(TFile* f,TNtuple* nt, Int_t nrun, Float_t *xntSA);
void FillSDDntuple(TFile* f,TNtuple* nt, Int_t iRun, Float_t *xnt);
void FillSSDntuple(TFile* f,TNtuple* ntssd, Int_t iRun, Float_t *xntSSD);
void FillMatchntuple(TFile* f,TNtuple* ntmatching, Int_t iRun, Float_t *xntMatching);
void FillMatchntupleTOF(TFile* f,TNtuple* ntmatchingTOF, Int_t iRun, Float_t *xntMatchingTOF);
void FillVTXntuple(TFile* f,TNtuple* ntvertex, Int_t iRun, Float_t *xntVertex);
void FillPileUpntuple(TFile* f,TNtuple* ntpileup, Int_t iRun, Float_t *xntPileup);
void FillPIDntuple(TFile* f,TNtuple* nt, Int_t nrun, Float_t *xntPID);
void FillDCAntuple(TFile* f,TNtuple* nt, Int_t nrun, Float_t *xntDCA);
void AliITSQAtrend(TString runListFile="lista_2016.txt",TString ntupleFileName="TrendingITS_2016.root");
Double_t LangausFun(Double_t *x, Double_t *par);

Bool_t WriteInputTextFileFromMonalisaListOfRuns(TString outtxtfilename,Int_t* listofrunsfromMonalisa,Int_t nruns,TString pathbeforRunN="alice/data/2012/LHC12a/",TString pathafterRunN="pass1/QAresults.root");
Int_t RunsToPlot(Int_t run1, Int_t run2,Int_t nr,Int_t *noRuns, Int_t *myIndex);
ofstream myfile("outfile1_fast.txt",ios::app);

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

    // histos for automatic QA
    TH1F* hFlagminTime=new TH1F("hFlagminTime","",kRunsToPlot,0.,kRunsToPlot);
    TH1F* hFlagmeanTime=new TH1F("hFlagmeanTime","",kRunsToPlot,0.,kRunsToPlot);
    TH1F* hFlagdEdx3=new TH1F("hFlagdEdx3","",kRunsToPlot,0.,kRunsToPlot);
    TH1F* hFlagdEdx4=new TH1F("hFlagdEdx4","",kRunsToPlot,0.,kRunsToPlot);
    
    TH1F *hVarSDD1 = new TH1F("hVarSDD1","SDD inner; run number; #Delta N/N modules SDD1",kRunsToPlot,0.,kRunsToPlot);
    hVarSDD1->SetLineColor(kGreen+2);
    hVarSDD1->SetMarkerColor(kGreen+2);
    hVarSDD1->SetMarkerStyle(20);
    
    TH1F *hVarSDD2 = new TH1F("hVarSDD2","SDD outer; run number; #Delta N/N modules SDD2",kRunsToPlot,0.,kRunsToPlot);
    hVarSDD2->SetLineColor(kYellow+2);
    hVarSDD2->SetMarkerColor(kYellow+2);
    hVarSDD2->SetMarkerStyle(20);
    
    TH1F *hFlagSDD1 = new TH1F("hFlagSDD1","SDD inner; run number; SDD1 alarm flag",kRunsToPlot,0.,kRunsToPlot);
    hFlagSDD1->SetLineColor(kGreen+2);
    hFlagSDD1->SetMarkerColor(kGreen+2);
    hFlagSDD1->SetMarkerStyle(20);
    
    TH1F *hFlagSDD2 = new TH1F("hFlagSDD2","SDD outer; run number; SDD2 alarm flag",kRunsToPlot,0.,kRunsToPlot);
    hFlagSDD2->SetLineColor(kYellow+2);
    hFlagSDD2->SetMarkerColor(kYellow+2);
    hFlagSDD2->SetMarkerStyle(20);
    
    //
    Float_t minT=1.0, meanT=2.0,dEdx3=1.0,dEdx4=2.0;

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
      minT=1.0;
      if(490.>minDrTime || minDrTime>510.) minT=0.5;
      hFlagminTime->SetBinContent(i+1,minT);
      hFlagminTime->SetBinError(i+1,0.01);
      hFlagminTime->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    histomeanTime->SetBinContent(i+1,meanDrTime);
    histomeanTime->SetBinError(i+1,errmeanDrTime);
    histomeanTime->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      meanT=2.0;
      if(3180.>meanDrTime || meanDrTime>3240.) meanT=1.5;
      hFlagmeanTime->SetBinContent(i+1,meanT);
      hFlagmeanTime->SetBinError(i+1,0.01);
      hFlagmeanTime->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
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
      dEdx3=1.0;
      if(83.>MPVdEdxLay3 || MPVdEdxLay3>85.) dEdx3=0.5;
      hFlagdEdx3->SetBinContent(i+1,dEdx3);
      hFlagdEdx3->SetBinError(i+1,0.01);
      hFlagdEdx3->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    histodEdxLay4->SetBinContent(i+1,MPVdEdxLay4);
    histodEdxLay4->SetBinError(i+1,errMPVdEdxLay4);
    histodEdxLay4->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
      dEdx4=2.0;
      if(83.>MPVdEdxLay4 || MPVdEdxLay4>85.) dEdx4=1.5;
      hFlagdEdx4->SetBinContent(i+1,dEdx4);
      hFlagdEdx4->SetBinError(i+1,0.01);
      hFlagdEdx4->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
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

    // filling histos automatic QA: check variations
    
    Float_t diff1 =0.0, diff2=0.2, flag1=1.0, flag2=2.0;
    
    for(Int_t t=0;t<kRunsToPlot;t++){
        
        if(t==0){
            diff1=0.0;
            diff2=0.2;
        }
        else{
            diff1=(histoFracDead3->GetBinContent(t+1)-histoFracDead3->GetBinContent(t))/histoFracDead3->GetBinContent(t);
            diff2=(histoFracDead4->GetBinContent(t+1)-histoFracDead4->GetBinContent(t))/histoFracDead4->GetBinContent(t) +0.2;
            //            cout << "i = " << t+1 << " Frac[i] = " << hFracSPD1->GetBinContent(t) << " Frac[i-1] = " << hFracSPD1->GetBinContent(t-1) << endl;
        }
        flag1=1.0;
        if(TMath::Abs(diff1) > 0.01) flag1=0.5; // variazione di >= 1 modulo (0.012)
        flag2=2.0;
        if(TMath::Abs(diff2-0.2) > 0.005) flag2=1.5; // variazione di >= 1 HS (0.006)
        
        
        hVarSDD1->SetBinContent(t+1,diff1);
        hVarSDD1->SetBinError(t+1,0.01);
        hVarSDD1->GetXaxis()->SetBinLabel(t+1,Form("%d",(Int_t)noRuns[t]));
        hVarSDD2->SetBinContent(t+1,diff2);
        hVarSDD2->SetBinError(t+1,0.01);
        hVarSDD2->GetXaxis()->SetBinLabel(t+1,Form("%d",(Int_t)noRuns[t]));
        
        hFlagSDD1->SetBinContent(t+1,flag1);
        hFlagSDD1->SetBinError(t+1,0.01);
        hFlagSDD1->GetXaxis()->SetBinLabel(t+1,Form("%d",(Int_t)noRuns[t]));
        hFlagSDD2->SetBinContent(t+1,flag2);
        hFlagSDD2->SetBinError(t+1,0.01);
        hFlagSDD2->GetXaxis()->SetBinLabel(t+1,Form("%d",(Int_t)noRuns[t]));
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
    
    // histos for automatic QA
    TH1F* hFlagChR5=new TH1F("hFlagChR5","",kRunsToPlot,0.,kRunsToPlot);
    TH1F* hFlagChR6=new TH1F("hFlagChR6","",kRunsToPlot,0.,kRunsToPlot);
    TH1F* hFlagdEdx5=new TH1F("hFlagdEdx5","",kRunsToPlot,0.,kRunsToPlot);
    TH1F* hFlagdEdx6=new TH1F("hFlagdEdx6","",kRunsToPlot,0.,kRunsToPlot);
    
    
    TH1F *hVarSSD1n = new TH1F("hVarSSD1n","SSD inner; run number; #Delta N/N n-strips SSD1",kRunsToPlot,0.,kRunsToPlot);
    hVarSSD1n->SetLineColor(6);
    hVarSSD1n->SetMarkerColor(6);
    hVarSSD1n->SetMarkerStyle(20);
    
    TH1F *hVarSSD1p = new TH1F("hVarSSD1p","SSD inner; run number; #Delta N/N p-strips SSD1",kRunsToPlot,0.,kRunsToPlot);
    hVarSSD1p->SetLineColor(kMagenta+2);
    hVarSSD1p->SetMarkerColor(kMagenta+2);
    hVarSSD1p->SetMarkerStyle(21);
    
    TH1F *hVarSSD2n = new TH1F("hVarSSD2n","SDD outer; run number; #Delta N/N p-strips SSD2",kRunsToPlot,0.,kRunsToPlot);
    hVarSSD2n->SetLineColor(9);
    hVarSSD2n->SetMarkerColor(9);
    hVarSSD2n->SetMarkerStyle(20);
    
    TH1F *hVarSSD2p = new TH1F("hVarSSD2p","SDD outer; run number; #Delta N/N p-strips SSD2",kRunsToPlot,0.,kRunsToPlot);
    hVarSSD2p->SetLineColor(38);
    hVarSSD2p->SetMarkerColor(38);
    hVarSSD2p->SetMarkerStyle(21);
    
    TH1F *hFlagSSD1n = new TH1F("hFlagSSD1n","SSD inner; run number; SSD1 n-strips alarm flag",kRunsToPlot,0.,kRunsToPlot);
    hFlagSSD1n->SetLineColor(6);
    hFlagSSD1n->SetMarkerColor(6);
    hFlagSSD1n->SetMarkerStyle(20);
    
    TH1F *hFlagSSD1p = new TH1F("hFlagSSD1p","SSD inner; run number; SSD1 p-strips alarm flag",kRunsToPlot,0.,kRunsToPlot);
    hFlagSSD1p->SetLineColor(kMagenta+2);
    hFlagSSD1p->SetMarkerColor(kMagenta+2);
    hFlagSSD1p->SetMarkerStyle(21);
    
    TH1F *hFlagSSD2n = new TH1F("hFlagSSD2n","SSD inner; run number; SSD2 n-strips alarm flag",kRunsToPlot,0.,kRunsToPlot);
    hFlagSSD2n->SetLineColor(9);
    hFlagSSD2n->SetMarkerColor(9);
    hFlagSSD2n->SetMarkerStyle(20);
    
    TH1F *hFlagSSD2p = new TH1F("hFlagSSD2p","SSD inner; run number; SSD2 p-strips alarm flag",kRunsToPlot,0.,kRunsToPlot);
    hFlagSSD2p->SetLineColor(38);
    hFlagSSD2p->SetMarkerColor(38);
    hFlagSSD2p->SetMarkerStyle(21);
    //
    Float_t dEdx5, dEdx6, ChR5, ChR6;

    for(Int_t i=0; i<kRunsToPlot;i++){

    ntssd->GetEvent(myIndex[i]);

    histodEdxLay5->SetBinContent(i+1,MPVdEdxLay5);
    histodEdxLay5->SetBinError(i+1,errMPVdEdxLay5);
    histodEdxLay5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunSSD));
        dEdx5=1.0;
        if(82.>MPVdEdxLay5 || MPVdEdxLay5>83.) dEdx5=0.5;
        hFlagdEdx5->SetBinContent(i+1,dEdx5);
        hFlagdEdx5->SetBinError(i+1,0.01);
        hFlagdEdx5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunSSD));

    histodEdxLay6->SetBinContent(i+1,MPVdEdxLay6);
    histodEdxLay6->SetBinError(i+1,errMPVdEdxLay6);
    histodEdxLay6->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunSSD));
        dEdx6=2.0;
        if(82.>MPVdEdxLay6 || MPVdEdxLay6>83.) dEdx6=1.5;
        hFlagdEdx6->SetBinContent(i+1,dEdx6);
        hFlagdEdx6->SetBinError(i+1,0.01);
        hFlagdEdx6->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunSSD));

    histoChargeRatioLay5->SetBinContent(i+1,ChargeRatioL5);
    histoChargeRatioLay5->SetBinError(i+1,errChargeratioL5);
    histoChargeRatioLay5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunSSD));
        ChR5=1.0;
        if(-0.01>ChargeRatioL5 || ChargeRatioL5>0.01) ChR5=0.5;
        hFlagChR5->SetBinContent(i+1,ChR5);
        hFlagChR5->SetBinError(i+1,0.01);
        hFlagChR5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunSSD));

    histoChargeRatioLay6->SetBinContent(i+1,ChargeRatioL6);
    histoChargeRatioLay6->SetBinError(i+1,errChargeratioL6);
    histoChargeRatioLay6->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunSSD));
        ChR6=2.0;
        if(-0.01>ChargeRatioL6 || ChargeRatioL6>0.01) ChR6=1.5;
        hFlagChR6->SetBinContent(i+1,ChR6);
        hFlagChR6->SetBinError(i+1,0.01);
        hFlagChR6->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunSSD));

    histoEmpty->SetBinContent(i+1,moduleOff);
    histoEmpty->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunSSD));
        
    histoFracBadn5->SetBinContent(i+1,FracBadn5);
//        cout << " run " << i << " FracBadn5 = " << FracBadn5 << endl;
    histoFracBadn5->SetBinError(i+1,errFracBadn5);
    histoFracBadn5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunSSD));

    histoFracBadp5->SetBinContent(i+1,FracBadp5);
//        cout << " run " << i << " FracBadp5 = " << FracBadp5 << endl;
    histoFracBadp5->SetBinError(i+1,errFracBadp5);
    histoFracBadp5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunSSD));

    histoFracBadn6->SetBinContent(i+1,FracBadn6);
//        cout << " run " << i << " FracBadn6 = " << FracBadn6 << endl;
    histoFracBadn6->SetBinError(i+1,errFracBadn6);
    histoFracBadn6->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunSSD));
        
    histoFracBadp6->SetBinContent(i+1,FracBadp6);
//        cout << " run " << i << " FracBadp6 = " << FracBadp6 << endl;
    histoFracBadp6->SetBinError(i+1,errFracBadp6);
    histoFracBadp6->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunSSD));
        
  }

    // filling histos automatic QA: check variations
    
    Float_t diff1n =0.0, diff1p=0.0, diff2n=0.2, diff2p=0.2, flag1n=1.0, flag1p=1.0,flag2n=2.0, flag2p=2.0;
    
    for(Int_t t=0;t<kRunsToPlot;t++){
        
        if(t==0){
            diff1n=0.0,diff1p=0.1;
            diff2n=0.2,diff2p=0.3;
        }
        else{
            diff1n=(histoFracBadn5->GetBinContent(t+1)-histoFracBadn5->GetBinContent(t))/histoFracBadn5->GetBinContent(t);
            diff1p=(histoFracBadp5->GetBinContent(t+1)-histoFracBadp5->GetBinContent(t))/histoFracBadp5->GetBinContent(t)+ 0.1;
            diff2n=(histoFracBadn6->GetBinContent(t+1)-histoFracBadn6->GetBinContent(t))/histoFracBadn6->GetBinContent(t) +0.2;
            diff2p=(histoFracBadp6->GetBinContent(t+1)-histoFracBadp6->GetBinContent(t))/histoFracBadp6->GetBinContent(t) +0.3;
        }
        
        flag1n=1.0;
        if(TMath::Abs(diff1n) > 0.01) flag1n=0.5; // variazione di >= 0.01 ...??
        flag1p=1.0;
        if(TMath::Abs(diff1p-0.1) > 0.01) flag1p=0.5; // variazione di >= 0.01 ...??
        flag2n=2.0;
        if(TMath::Abs(diff2n-0.2) > 0.01) flag2n=1.5; // variazione di >= 0.01 ... ??
        flag2p=2.0;
        if(TMath::Abs(diff2p-0.3) > 0.01) flag2p=1.5; // variazione di >= 0.01 ... ??
        
        hVarSSD1n->SetBinContent(t+1,diff1n);
        hVarSSD1n->SetBinError(t+1,0.01);
        hVarSSD1n->GetXaxis()->SetBinLabel(t+1,Form("%d",(Int_t)noRuns[t]));
        hVarSSD1p->SetBinContent(t+1,diff1p);
        hVarSSD1p->SetBinError(t+1,0.01);
        hVarSSD1p->GetXaxis()->SetBinLabel(t+1,Form("%d",(Int_t)noRuns[t]));
        hVarSSD2n->SetBinContent(t+1,diff2n);
        hVarSSD2n->SetBinError(t+1,0.01);
        hVarSSD2n->GetXaxis()->SetBinLabel(t+1,Form("%d",(Int_t)noRuns[t]));
        hVarSSD2p->SetBinContent(t+1,diff2p);
        hVarSSD2p->SetBinError(t+1,0.01);
        hVarSSD2p->GetXaxis()->SetBinLabel(t+1,Form("%d",(Int_t)noRuns[t]));
        
        hFlagSSD1n->SetBinContent(t+1,flag1n);
        hFlagSSD1n->SetBinError(t+1,0.01);
        hFlagSSD1n->GetXaxis()->SetBinLabel(t+1,Form("%d",(Int_t)noRuns[t]));
        hFlagSSD1p->SetBinContent(t+1,flag1p);
        hFlagSSD1p->SetBinError(t+1,0.01);
        hFlagSSD1p->GetXaxis()->SetBinLabel(t+1,Form("%d",(Int_t)noRuns[t]));
        hFlagSSD2n->SetBinContent(t+1,flag2n);
        hFlagSSD2n->SetBinError(t+1,0.01);
        hFlagSSD2n->GetXaxis()->SetBinLabel(t+1,Form("%d",(Int_t)noRuns[t]));
        hFlagSSD2p->SetBinContent(t+1,flag2p);
        hFlagSSD2p->SetBinError(t+1,0.01);
        hFlagSSD2p->GetXaxis()->SetBinLabel(t+1,Form("%d",(Int_t)noRuns[t]));
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
  printf("\n ======== PROCESSING NTMATCHING TPC NTUPLE \n");
  kRunsToPlot = RunsToPlot(run1,run2,nr,noRuns,myIndex);

  TH1F *hFracSPD1 = new TH1F("hFracSPD1","SPD inner; run number; Fraction of HSs",kRunsToPlot,0.,kRunsToPlot);
  hFracSPD1->SetLineColor(kGreen+2);
  hFracSPD1->SetMarkerColor(kGreen+2);
  hFracSPD1->SetMarkerStyle(20);

  TH1F *hFracSPD2 = new TH1F("hFracSPD2","SPD outer; run number; Fraction of HSs",kRunsToPlot,0.,kRunsToPlot);
  hFracSPD2->SetLineColor(kYellow+2);
  hFracSPD2->SetMarkerColor(kYellow+2);
  hFracSPD2->SetMarkerStyle(20);

    // histos for automatic QA
    
    TH1F *hVarSPD1 = new TH1F("hVarSPD1","SPD inner; run number; #Delta N/N HS SPD1",kRunsToPlot,0.,kRunsToPlot);
    hVarSPD1->SetLineColor(kGreen+2);
    hVarSPD1->SetMarkerColor(kGreen+2);
    hVarSPD1->SetMarkerStyle(20);
    
    TH1F *hVarSPD2 = new TH1F("hVarSPD2","SPD outer; run number; #Delta N/N HS SPD2",kRunsToPlot,0.,kRunsToPlot);
    hVarSPD2->SetLineColor(kYellow+2);
    hVarSPD2->SetMarkerColor(kYellow+2);
    hVarSPD2->SetMarkerStyle(20);
    
    TH1F *hFlagSPD1 = new TH1F("hFlagSPD1","SPD inner; run number; SPD1 alarm flag",kRunsToPlot,0.,kRunsToPlot);
    hFlagSPD1->SetLineColor(kGreen+2);
    hFlagSPD1->SetMarkerColor(kGreen+2);
    hFlagSPD1->SetMarkerStyle(20);
    
    TH1F *hFlagSPD2 = new TH1F("hFlagSPD2","SPD outer; run number; SPD2 alarm flag",kRunsToPlot,0.,kRunsToPlot);
    hFlagSPD2->SetLineColor(kYellow+2);
    hFlagSPD2->SetMarkerColor(kYellow+2);
    hFlagSPD2->SetMarkerStyle(20);
    
    //

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

/*    // histos for automatic QA
    TH1F *hEffoneSPDPt02norm = new TH1F("hEffoneSPDPt02norm","Efficiency; run number; TPC+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
    hEffoneSPDPt02norm->SetLineWidth(2);
    hEffoneSPDPt02norm->SetLineColor(kGray+2);
    hEffoneSPDPt02norm->SetMarkerColor(kGray+2);
    hEffoneSPDPt02norm->SetMarkerStyle(20);
    TH1F *hEffoneSPDPt1norm = new TH1F("hEffoneSPDPt1norm","Efficiency; run number; TPC+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
    hEffoneSPDPt1norm->SetLineWidth(2);
    hEffoneSPDPt1norm->SetLineColor(kGray+2);
    hEffoneSPDPt1norm->SetMarkerColor(kGray+2);
    hEffoneSPDPt1norm->SetMarkerStyle(20);
    TH1F *hEffoneSPDPt10norm = new TH1F("hEffoneSPDPt10norm","Efficiency; run number; TPC+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
    hEffoneSPDPt10norm->SetLineWidth(2);
    hEffoneSPDPt10norm->SetLineColor(kGray+2);
    hEffoneSPDPt10norm->SetMarkerColor(kGray+2);
    hEffoneSPDPt10norm->SetMarkerStyle(20);
    
    TH1F *hEff456Pt02norm = new TH1F("hEff456Pt02norm","Efficiency; run number; TPC+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
    hEff456Pt02norm->SetLineWidth(2);
    hEff456Pt02norm->SetLineColor(kRed);
    hEff456Pt02norm->SetMarkerColor(kRed);
    hEff456Pt02norm->SetMarkerStyle(20);
    TH1F *hEff456Pt1norm = new TH1F("hEff456Pt1norm","Efficiency; run number; TPC+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
    hEff456Pt1norm->SetLineWidth(2);
    hEff456Pt1norm->SetLineColor(kRed);
    hEff456Pt1norm->SetMarkerColor(kRed);
    hEff456Pt1norm->SetMarkerStyle(20);
    TH1F *hEff456Pt10norm = new TH1F("hEff456Pt10norm","Efficiency; run number; TPC+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
    hEff456Pt10norm->SetLineWidth(2);
    hEff456Pt10norm->SetLineColor(kRed);
    hEff456Pt10norm->SetMarkerColor(kRed);
    hEff456Pt10norm->SetMarkerStyle(20);
 // automatic QA
 */
    
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
    
/*    TH1F* histoTrackMI1norm=new TH1F("histoTrackMI1norm","",kRunsToPlot,0.,kRunsToPlot);
    TH1F* histoTrackMI2norm=new TH1F("histoTrackMI2norm","",kRunsToPlot,0.,kRunsToPlot);
    TH1F* histoTrackMI3norm=new TH1F("histoTrackMI3norm","",kRunsToPlot,0.,kRunsToPlot);
    TH1F* histoTrackMI4norm=new TH1F("histoTrackMI4norm","",kRunsToPlot,0.,kRunsToPlot);
    TH1F* histoTrackMI5norm=new TH1F("histoTrackMI5norm","",kRunsToPlot,0.,kRunsToPlot);
    TH1F* histoTrackMI6norm=new TH1F("histoTrackMI6norm","",kRunsToPlot,0.,kRunsToPlot);
*/
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

/*      // automatic QA
      hEffoneSPDPt02norm->SetBinContent(i+1,EffoneSPDPt02/(0.33*hFracSPD1->GetBinContent(i+1)+0.66*hFracSPD2->GetBinContent(i+1)));
      hEffoneSPDPt02norm->SetBinError(i+1,errEffoneSPDPt02/(0.33*hFracSPD1->GetBinContent(i+1)+0.66*hFracSPD2->GetBinContent(i+1)));
      hEffoneSPDPt02norm->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));
//      cout << "run " << noRuns[i] << ", norm SPD1/2 = " << 0.33*hFracSPD1->GetBinContent(i+1)+0.66*hFracSPD2->GetBinContent(i+1) << endl;
      
      hEffoneSPDPt1norm->SetBinContent(i+1,EffoneSPDPt1/(0.33*hFracSPD1->GetBinContent(i+1)+0.66*hFracSPD2->GetBinContent(i+1)));
      hEffoneSPDPt1norm->SetBinError(i+1,errEffoneSPDPt1/(0.33*hFracSPD1->GetBinContent(i+1)+0.66*hFracSPD2->GetBinContent(i+1)));
      hEffoneSPDPt1norm->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));
      
      hEffoneSPDPt10norm->SetBinContent(i+1,EffoneSPDPt10/(0.33*hFracSPD1->GetBinContent(i+1)+0.66*hFracSPD2->GetBinContent(i+1)));
      hEffoneSPDPt10norm->SetBinError(i+1,errEffoneSPDPt10/(0.33*hFracSPD1->GetBinContent(i+1)+0.66*hFracSPD2->GetBinContent(i+1)));
      hEffoneSPDPt10norm->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));
      
      Float_t norm5 = 1-0.5*(histoFracBadn5->GetBinContent(i+1)+histoFracBadp5->GetBinContent(i+1));
      Float_t norm6 = 1-0.5*(histoFracBadn6->GetBinContent(i+1)+histoFracBadp6->GetBinContent(i+1));
      Float_t normtot = hFracSPD1->GetBinContent(i+1)+hFracSPD2->GetBinContent(i+1)+histoFracDead3->GetBinContent(i+1);
      normtot = normtot+histoFracDead3->GetBinContent(i+1)+norm5+norm6;
      normtot /=6.;
//      cout << "run " << noRuns[i] << ", normtot456 = " << normtot << endl << endl;
      
      hEff456Pt02norm->SetBinContent(i+1,(Eff4Pt02+Eff5Pt02+Eff6Pt02)/normtot);
      hEff456Pt02norm->SetBinError(i+1,(errEff4Pt02+errEff5Pt02+errEff6Pt02)/normtot);
      hEff456Pt02norm->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));
      
      hEff456Pt1norm->SetBinContent(i+1,(Eff4Pt1+Eff5Pt1+Eff6Pt1)/normtot);
      hEff456Pt1norm->SetBinError(i+1,(errEff4Pt1+errEff5Pt1+errEff6Pt1)/normtot);
      hEff456Pt1norm->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));
      
      hEff456Pt10norm->SetBinContent(i+1,(Eff4Pt10+Eff5Pt10+Eff6Pt10)/normtot);
      hEff456Pt10norm->SetBinError(i+1,(errEff4Pt10+errEff5Pt10+errEff6Pt10)/normtot);
      hEff456Pt10norm->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));
// automatic QA
*/
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

/*      histoTrackMI1norm->SetBinContent(i+1,FracTrackMI1/hFracSPD1->GetBinContent(i+1));
      histoTrackMI1norm->SetBinError(i+1,errFracTrackMI1/hFracSPD1->GetBinContent(i+1));
      histoTrackMI1norm->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));
      
      histoTrackMI2norm->SetBinContent(i+1,FracTrackMI2/hFracSPD2->GetBinContent(i+1));
      histoTrackMI2norm->SetBinError(i+1,errFracTrackMI2/hFracSPD2->GetBinContent(i+1));
      histoTrackMI2norm->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));
      
      histoTrackMI3norm->SetBinContent(i+1,FracTrackMI3/histoFracDead3->GetBinContent(i+1));
      histoTrackMI3norm->SetBinError(i+1,errFracTrackMI3/histoFracDead3->GetBinContent(i+1));
      histoTrackMI3norm->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));
      
      histoTrackMI4norm->SetBinContent(i+1,FracTrackMI4/histoFracDead4->GetBinContent(i+1));
      histoTrackMI4norm->SetBinError(i+1,errFracTrackMI4/histoFracDead4->GetBinContent(i+1));
      histoTrackMI4norm->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));
      
      histoTrackMI5norm->SetBinContent(i+1,FracTrackMI5/(1-0.5*(histoFracBadn5->GetBinContent(i+1)+histoFracBadp5->GetBinContent(i+1))));
      histoTrackMI5norm->SetBinError(i+1,errFracTrackMI5/(1-0.5*(histoFracBadn5->GetBinContent(i+1)+histoFracBadp5->GetBinContent(i+1))));
      histoTrackMI5norm->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));
      
      histoTrackMI6norm->SetBinContent(i+1,FracTrackMI6/(1-0.5*(histoFracBadn6->GetBinContent(i+1)+histoFracBadp6->GetBinContent(i+1))));
      histoTrackMI6norm->SetBinError(i+1,errFracTrackMI6/(1-0.5*(histoFracBadn6->GetBinContent(i+1)+histoFracBadp6->GetBinContent(i+1))));
      histoTrackMI6norm->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));
*/
  }

    // filling histos automatic QA: check variations
    
    diff1 =0.0, diff2=0.2, flag1=1.0, flag2=2.0;
    
    for(Int_t t=0;t<kRunsToPlot;t++){
        
        if(t==0){
            diff1=0.0;
            diff2=0.2;
        }
        else{
            diff1=(hFracSPD1->GetBinContent(t+1)-hFracSPD1->GetBinContent(t))/hFracSPD1->GetBinContent(t);
            diff2=(hFracSPD2->GetBinContent(t+1)-hFracSPD2->GetBinContent(t))/hFracSPD2->GetBinContent(t) +0.2;
//            cout << "i = " << t+1 << " Frac[i] = " << hFracSPD1->GetBinContent(t) << " Frac[i-1] = " << hFracSPD1->GetBinContent(t-1) << endl;
        }
        flag1=1.0;
        if(TMath::Abs(diff1) > 0.02) flag1=0.5; // variazione di >= 1 HS (0.025)
        flag2=2.0;
        if(TMath::Abs(diff2-0.2) > 0.01) flag2=1.5; // variazione di >= 1 HS (0.0125)
        
        
        hVarSPD1->SetBinContent(t+1,diff1);
        hVarSPD1->SetBinError(t+1,0.01);
        hVarSPD1->GetXaxis()->SetBinLabel(t+1,Form("%d",(Int_t)noRuns[t]));
        hVarSPD2->SetBinContent(t+1,diff2);
        hVarSPD2->SetBinError(t+1,0.01);
        hVarSPD2->GetXaxis()->SetBinLabel(t+1,Form("%d",(Int_t)noRuns[t]));
        
        hFlagSPD1->SetBinContent(t+1,flag1);
        hFlagSPD1->SetBinError(t+1,0.01);
        hFlagSPD1->GetXaxis()->SetBinLabel(t+1,Form("%d",(Int_t)noRuns[t]));
        hFlagSPD2->SetBinContent(t+1,flag2);
        hFlagSPD2->SetBinError(t+1,0.01);
        hFlagSPD2->GetXaxis()->SetBinLabel(t+1,Form("%d",(Int_t)noRuns[t]));
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
  printf("\n ======== PROCESSING NTMATCHING TOF NTUPLE \n");
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

/*    // automatic QA
    TH1F *hEffTOFoneSPDPt02norm = new TH1F("hEffTOFoneSPDPt02norm","Efficiency; run number; TOF+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
    hEffTOFoneSPDPt02norm->SetLineWidth(2);
    hEffTOFoneSPDPt02norm->SetLineColor(kGray+2);
    hEffTOFoneSPDPt02norm->SetMarkerColor(kGray+2);
    hEffTOFoneSPDPt02norm->SetMarkerStyle(20);
    TH1F *hEffTOFoneSPDPt1norm = new TH1F("hEffTOFoneSPDPt1norm","Efficiency; run number; TOF+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
    hEffTOFoneSPDPt1norm->SetLineWidth(2);
    hEffTOFoneSPDPt1norm->SetLineColor(kGray+2);
    hEffTOFoneSPDPt1norm->SetMarkerColor(kGray+2);
    hEffTOFoneSPDPt1norm->SetMarkerStyle(20);
    TH1F *hEffTOFoneSPDPt10norm = new TH1F("hEffTOFoneSPDPt10norm","Efficiency; run number; TOF+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
    hEffTOFoneSPDPt10norm->SetLineWidth(2);
    hEffTOFoneSPDPt10norm->SetLineColor(kGray+2);
    hEffTOFoneSPDPt10norm->SetMarkerColor(kGray+2);
    hEffTOFoneSPDPt10norm->SetMarkerStyle(20);
    
    TH1F *hEffTOF456Pt02norm = new TH1F("hEffTOF456Pt02norm","Efficiency; run number; TOF+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
    hEffTOF456Pt02norm->SetLineWidth(2);
    hEffTOF456Pt02norm->SetLineColor(kRed);
    hEffTOF456Pt02norm->SetMarkerColor(kRed);
    hEffTOF456Pt02norm->SetMarkerStyle(20);
    TH1F *hEffTOF456Pt1norm = new TH1F("hEffTOF456Pt1norm","Efficiency; run number; TOF+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
    hEffTOF456Pt1norm->SetLineWidth(2);
    hEffTOF456Pt1norm->SetLineColor(kRed);
    hEffTOF456Pt1norm->SetMarkerColor(kRed);
    hEffTOF456Pt1norm->SetMarkerStyle(20);
    TH1F *hEffTOF456Pt10norm = new TH1F("hEffTOF456Pt10norm","Efficiency; run number; TOF+ITS / TPC",kRunsToPlot,0.,kRunsToPlot);
    hEffTOF456Pt10norm->SetLineWidth(2);
    hEffTOF456Pt10norm->SetLineColor(kRed);
    hEffTOF456Pt10norm->SetMarkerColor(kRed);
    hEffTOF456Pt10norm->SetMarkerStyle(20);
// automatic QA
*/
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

/*      // automatic QA
      
      hEffTOFoneSPDPt02norm->SetBinContent(i+1,EffoneSPDPt02TOF/(0.5*(hFracSPD1->GetBinContent(i+1)+hFracSPD2->GetBinContent(i+1))));
      hEffTOFoneSPDPt02norm->SetBinError(i+1,errEffoneSPDPt02TOF/(0.5*(hFracSPD1->GetBinContent(i+1)+hFracSPD2->GetBinContent(i+1))));
      hEffTOFoneSPDPt02norm->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));
      
      hEffTOFoneSPDPt1norm->SetBinContent(i+1,EffoneSPDPt1TOF/(0.5*(hFracSPD1->GetBinContent(i+1)+hFracSPD2->GetBinContent(i+1))));
      hEffTOFoneSPDPt1norm->SetBinError(i+1,errEffoneSPDPt1TOF/(0.5*(hFracSPD1->GetBinContent(i+1)+hFracSPD2->GetBinContent(i+1))));
      hEffTOFoneSPDPt1norm->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));
      
      hEffTOFoneSPDPt10norm->SetBinContent(i+1,EffoneSPDPt10TOF/(0.5*(hFracSPD1->GetBinContent(i+1)+hFracSPD2->GetBinContent(i+1))));
      hEffTOFoneSPDPt10norm->SetBinError(i+1,errEffoneSPDPt10TOF/(0.5*(hFracSPD1->GetBinContent(i+1)+hFracSPD2->GetBinContent(i+1))));
      hEffTOFoneSPDPt10norm->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));
      
      Float_t norm5 = 1-0.5*(histoFracBadn5->GetBinContent(i+1)+histoFracBadp5->GetBinContent(i+1));
      Float_t norm6 = 1-0.5*(histoFracBadn6->GetBinContent(i+1)+histoFracBadp6->GetBinContent(i+1));
      Float_t normtot = hFracSPD1->GetBinContent(i+1)+hFracSPD2->GetBinContent(i+1)+histoFracDead3->GetBinContent(i+1);
      normtot = normtot+histoFracDead3->GetBinContent(i+1)+norm5+norm6;
      normtot /=6.;
      
      hEffTOF456Pt02norm->SetBinContent(i+1,(Eff4Pt02TOF+Eff5Pt02TOF+Eff6Pt02TOF)/normtot);
      hEffTOF456Pt02norm->SetBinError(i+1,(errEff4Pt02TOF+errEff5Pt02TOF+errEff6Pt02TOF)/normtot);
      hEffTOF456Pt02norm->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));
      
      hEffTOF456Pt1norm->SetBinContent(i+1,(Eff4Pt1TOF+Eff5Pt1TOF+Eff6Pt1TOF)/normtot);
      hEffTOF456Pt1norm->SetBinError(i+1,(errEff4Pt1TOF+errEff5Pt1TOF+errEff6Pt1TOF)/normtot);
      hEffTOF456Pt1norm->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));
      
      hEffTOF456Pt10norm->SetBinContent(i+1,(Eff4Pt10TOF+Eff5Pt10TOF+Eff6Pt10TOF)/normtot);
      hEffTOF456Pt10norm->SetBinError(i+1,(errEff4Pt10TOF+errEff5Pt10TOF+errEff6Pt10TOF)/normtot);
      hEffTOF456Pt10norm->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunMatch));
// automatic QA
 */
  }

  //-----------------------------------

    //--- ITS PID
    //    cout << "inizio PID" << endl;
    
    TNtuple* ntPID=(TNtuple*)fil->Get("ntPIDITS");
    //    cout << "ntupla aperta" << endl;
    
    
    Float_t nrunPID,nsigpi02,errnsigpi02,nsigpi05,errnsigpi05,nsigpi1,errnsigpi1,nsigpi3,errnsigpi3;
    
    ntPID->SetBranchAddress("run",&nrunPID);
    ntPID->SetBranchAddress("nsigmapi02",&nsigpi02);
    ntPID->SetBranchAddress("errnsigmapi02",&errnsigpi02);
    ntPID->SetBranchAddress("nsigmapi05",&nsigpi05);
    ntPID->SetBranchAddress("errnsigmapi05",&errnsigpi05);
    ntPID->SetBranchAddress("nsigmapi1",&nsigpi1);
    ntPID->SetBranchAddress("errnsigmapi1",&errnsigpi1);
    ntPID->SetBranchAddress("nsigmapi3",&nsigpi3);
    ntPID->SetBranchAddress("errnsigmapi3",&errnsigpi3);
    
    
    nr=ntPID->GetEntries();
    delete []myIndex;
    delete []noRuns;
    myIndex = new Int_t [nr];
    noRuns = new Int_t [nr];
    for(Int_t i=0; i<nr;i++){
        ntPID->GetEvent(i);
        Int_t intrun = static_cast<Int_t>(nrunPID+0.01);
        noRuns[i]=intrun;
    }
    printf("\n ======== PROCESSING ITS PID NTUPLE \n");
    kRunsToPlot = RunsToPlot(run1,run2,nr,noRuns,myIndex);
    
    
    TH1F *hNsig02 = new TH1F("hNsig02","Pion nsigma pt=0.2 GeV/c",kRunsToPlot,0.,kRunsToPlot);
    hNsig02->SetLineWidth(2);
    hNsig02->SetLineColor(kBlue);
    hNsig02->SetMarkerColor(kBlue);
    hNsig02->SetMarkerStyle(20);
    
    TH1F *hNsig05 = new TH1F("hNsig05","Pion nsigma pt=0.5 GeV/c",kRunsToPlot,0.,kRunsToPlot);
    hNsig05->SetLineWidth(2);
    hNsig05->SetLineColor(kRed);
    hNsig05->SetMarkerColor(kRed);
    hNsig05->SetMarkerStyle(21);
    
    TH1F *hNsig1 = new TH1F("hNsig1","Pion nsigma pt=1.0 GeV/c",kRunsToPlot,0.,kRunsToPlot);
    hNsig1->SetLineWidth(2);
    hNsig1->SetLineColor(kGreen+1);
    hNsig1->SetMarkerColor(kGreen+1);
    hNsig1->SetMarkerStyle(22);
    
    TH1F *hNsig3 = new TH1F("hNsig3","Pion nsigma pt=2.0 GeV/c",kRunsToPlot,0.,kRunsToPlot);
    hNsig3->SetLineWidth(2);
    hNsig3->SetLineColor(kViolet);
    hNsig3->SetMarkerColor(kViolet);
    hNsig3->SetMarkerStyle(23);
    
    //  Int_t nEntriesVertex=ntvertex->GetEntries();
    
    for(Int_t i=0;i<kRunsToPlot;i++){
        
        ntPID->GetEvent(myIndex[i]);
        
        //    cout<<Vx<<endl;
        
        hNsig02->SetBinContent(i+1,nsigpi02);
        //      cout << "nsigma a 0.2 GeV/c = " << nsigpi02 << endl;
        hNsig02->SetBinError(i+1,errnsigpi02);
        hNsig02->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunPID));
        
        hNsig05->SetBinContent(i+1,nsigpi05);
        //      cout << "nsigma a 0.5 GeV/c = " << nsigpi05 << endl;
        hNsig05->SetBinError(i+1,errnsigpi05);
        hNsig05->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunPID));
        
        hNsig1->SetBinContent(i+1,nsigpi1);
        //      cout << "nsigma a 1.0 GeV/c = " << nsigpi1 << endl;
        hNsig1->SetBinError(i+1,errnsigpi1);
        hNsig1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunPID));
        
        hNsig3->SetBinContent(i+1,nsigpi3);
        //      cout << "nsigma a 3.0 GeV/c = " << nsigpi3 << endl;
        hNsig3->SetBinError(i+1,errnsigpi3);
        hNsig3->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunPID));
        
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
    gStyle->SetOptStat(0);
    TCanvas *cVertexDisto;
    if(hVx->GetEntries()>0){
   cVertexDisto=new TCanvas("cVertexDisto","cVertexDisto",1200,800);
  cVertexDisto->Divide(3,2);
  cVertexDisto->cd(1);
  hVx->SetMinimum(0.065);
  hVx->SetMaximum(0.105);
    hVx->GetYaxis()->SetTitle("Vertex X coordinate");
    hVx->GetXaxis()->SetTitle("run number");
  hVx->Draw();
  if(hVxSPD->GetBinContent(1)>0)hVxSPD->Draw("same");
        TLegend* legVtx=new TLegend(0.70,0.83,1.00,0.93);
        legVtx->SetFillColor(kWhite);
        legVtx->SetFillStyle(1001);
        TLegendEntry* entVtx;
        entVtx=legVtx->AddEntry(hVx,"Tracks vertex","PL");
        entVtx->SetTextColor(hVx->GetMarkerColor());
        entVtx=legVtx->AddEntry(hVxSPD,"Tracklets vertex","PL");
        entVtx->SetTextColor(hVxSPD->GetMarkerColor());
        legVtx->Draw();

  cVertexDisto->cd(2);
    if(hVySPD->GetEntries()>0){  // pp runs and pPb runs
  hVy->SetMinimum(0.28);
  hVy->SetMaximum(0.38);
    }
    else{                       // PbPb runs
        hVy->SetMinimum(0.32);
        hVy->SetMaximum(0.38);
    }
  hVy->GetYaxis()->SetTitle("Vertex Y coordinate");
    hVy->GetXaxis()->SetTitle("run number");
  hVy->Draw();
  if(hVySPD->GetBinContent(1)>0)hVySPD->Draw("same");
  legVtx->Draw();

  cVertexDisto->cd(3);
//   hVz->SetMinimum(-1.);
//   hVz->SetMaximum(1.);
  hVz->GetYaxis()->SetTitle("Vertex Z coordinate");
    hVz->GetXaxis()->SetTitle("run number");
  hVz->Draw();
  if(hVzSPD->GetBinContent(1)>0)hVzSPD->Draw("same");
        legVtx->Draw();

  cVertexDisto->cd(4);
  hSigmaVx->SetMinimum(0.);
  hSigmaVx->SetMaximum(0.1);
  hSigmaVx->GetYaxis()->SetTitle("sigma on x coordinate");
    hSigmaVx->GetXaxis()->SetTitle("run number");
  hSigmaVx->Draw();
  if(hSigmaVxSPD->GetBinContent(1)>0)hSigmaVxSPD->Draw("same");
        legVtx->Draw();

  cVertexDisto->cd(5);
 hSigmaVy->SetMinimum(0.);
 hSigmaVy->SetMaximum(0.2);
  hSigmaVy->GetYaxis()->SetTitle("sigma on y coordinate");
    hSigmaVy->GetXaxis()->SetTitle("run number");
  hSigmaVy->Draw();
  if(hSigmaVySPD->GetBinContent(1)>0)hSigmaVySPD->Draw("same");
        legVtx->Draw();
  cVertexDisto->cd(6);
//   hSigmaVz->SetMinimum(6.);
//   hSigmaVz->SetMaximum(10.);
  hSigmaVz->GetYaxis()->SetTitle("sigma on z coordinate");
    hSigmaVz->GetXaxis()->SetTitle("run number");
  hSigmaVz->Draw();
  if(hSigmaVzSPD->GetBinContent(1)>0)hSigmaVzSPD->Draw("same");
        legVtx->Draw();

  cVertexDisto->SaveAs("Vertex_trend.pdf");
//    pdfFileNames+=" Vertex_trend.pdf";
    }
    
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
//        printf("number of pileup vtx for run %d = %f \n",(Int_t)nrunpu, npilvtx);
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
    
    //--- TPCITS tracks DCA
    //    cout << "inizio DCA" << endl;
    
    TNtuple* ntDCA=(TNtuple*)fil->Get("ntDCAtracks");
    //    cout << "ntupla aperta" << endl;
    
    
    Float_t nrunDCA,mdca05,errmdca05,rmsdca05,errrmsdca05,mdca1,errmdca1,rmsdca1,errrmsdca1,mdca5,errmdca5,rmsdca5,errrmsdca5,mdca10,errmdca10,rmsdca10,errrmsdca10,mdcaz05,errmdcaz05,rmsdcaz05,errrmsdcaz05,mdcaz1,errmdcaz1,rmsdcaz1,errrmsdcaz1,mdcaz5,errmdcaz5,rmsdcaz5,errrmsdcaz5,mdcaz10,errmdcaz10,rmsdcaz10,errrmsdcaz10;
    
    ntDCA->SetBranchAddress("run",&nrunDCA);
    ntDCA->SetBranchAddress("dcamean_05",&mdca05);
    ntDCA->SetBranchAddress("errdcamean_05",&errmdca05);
    ntDCA->SetBranchAddress("dcaRMS_05",&rmsdca05);
    ntDCA->SetBranchAddress("errdcaRMS_05",&errrmsdca05);
    ntDCA->SetBranchAddress("dcamean_1",&mdca1);
    ntDCA->SetBranchAddress("errdcamean_1",&errmdca1);
    ntDCA->SetBranchAddress("dcaRMS_1",&rmsdca1);
    ntDCA->SetBranchAddress("errdcaRMS_1",&errrmsdca1);
    ntDCA->SetBranchAddress("dcamean_5",&mdca5);
    ntDCA->SetBranchAddress("errdcamean_5",&errmdca5);
    ntDCA->SetBranchAddress("dcaRMS_5",&rmsdca5);
    ntDCA->SetBranchAddress("errdcaRMS_5",&errrmsdca5);
    ntDCA->SetBranchAddress("dcamean_10",&mdca10);
    ntDCA->SetBranchAddress("errdcamean_10",&errmdca10);
    ntDCA->SetBranchAddress("dcaRMS_10",&rmsdca10);
    ntDCA->SetBranchAddress("errdcaRMS_10",&errrmsdca10);

    ntDCA->SetBranchAddress("dcazmean_05",&mdcaz05);
    ntDCA->SetBranchAddress("errdcazmean_05",&errmdcaz05);
    ntDCA->SetBranchAddress("dcazRMS_05",&rmsdcaz05);
    ntDCA->SetBranchAddress("errdcazRMS_05",&errrmsdcaz05);
    ntDCA->SetBranchAddress("dcazmean_1",&mdcaz1);
    ntDCA->SetBranchAddress("errdcazmean_1",&errmdcaz1);
    ntDCA->SetBranchAddress("dcazRMS_1",&rmsdcaz1);
    ntDCA->SetBranchAddress("errdcazRMS_1",&errrmsdcaz1);
    ntDCA->SetBranchAddress("dcazmean_5",&mdcaz5);
    ntDCA->SetBranchAddress("errdcazmean_5",&errmdcaz5);
    ntDCA->SetBranchAddress("dcazRMS_5",&rmsdcaz5);
    ntDCA->SetBranchAddress("errdcazRMS_5",&errrmsdcaz5);
    ntDCA->SetBranchAddress("dcazmean_10",&mdcaz10);
    ntDCA->SetBranchAddress("errdcazmean_10",&errmdcaz10);
    ntDCA->SetBranchAddress("dcazRMS_10",&rmsdcaz10);
    ntDCA->SetBranchAddress("errdcazRMS_10",&errrmsdcaz10);
    
    
    nr=ntDCA->GetEntries();
    delete []myIndex;
    delete []noRuns;
    myIndex = new Int_t [nr];
    noRuns = new Int_t [nr];
    for(Int_t i=0; i<nr;i++){
        ntDCA->GetEvent(i);
        Int_t intrun = static_cast<Int_t>(nrunDCA+0.01);
        noRuns[i]=intrun;
    }
    printf("\n ======== PROCESSING DCA NTUPLE \n");
    kRunsToPlot = RunsToPlot(run1,run2,nr,noRuns,myIndex);
    
    
    TH1F *hmDCA05 = new TH1F("hmDCA05","mean DCA for TPCITS tracks, 0.55<pt<0.56 GeV/c",kRunsToPlot,0.,kRunsToPlot);
    hmDCA05->SetLineWidth(2);
    hmDCA05->SetLineColor(kBlue);
    hmDCA05->SetMarkerColor(kBlue);
    hmDCA05->SetMarkerStyle(20);
    
    TH1F *hrmsDCA05 = new TH1F("hrmsDCA05","RMS DCA for TPCITS tracks, 0.55<pt<0.56 GeV/c",kRunsToPlot,0.,kRunsToPlot);
    hrmsDCA05->SetLineWidth(2);
    hrmsDCA05->SetLineColor(kBlue);
    hrmsDCA05->SetMarkerColor(kBlue);
    hrmsDCA05->SetMarkerStyle(20);
    
    TH1F *hmDCA1 = new TH1F("hmDCA1","mean DCA for TPCITS tracks, 1.05<pt<1.07 GeV/c",kRunsToPlot,0.,kRunsToPlot);
    hmDCA1->SetLineWidth(2);
    hmDCA1->SetLineColor(kRed);
    hmDCA1->SetMarkerColor(kRed);
    hmDCA1->SetMarkerStyle(21);
    
    TH1F *hrmsDCA1 = new TH1F("hrmsDCA1","RMS DCA for TPCITS tracks, 1.05<pt<1.07 GeV/c",kRunsToPlot,0.,kRunsToPlot);
    hrmsDCA1->SetLineWidth(2);
    hrmsDCA1->SetLineColor(kRed);
    hrmsDCA1->SetMarkerColor(kRed);
    hrmsDCA1->SetMarkerStyle(21);
    
    TH1F *hmDCA5 = new TH1F("hmDCA5","mean DCA for TPCITS tracks, 4.1<pt<5.2 GeV/c",kRunsToPlot,0.,kRunsToPlot);
    hmDCA5->SetLineWidth(2);
    hmDCA5->SetLineColor(kOrange+1);
    hmDCA5->SetMarkerColor(kOrange+1);
    hmDCA5->SetMarkerStyle(22);
    
    TH1F *hrmsDCA5 = new TH1F("hrmsDCA5","RMS DCA for TPCITS tracks, 4.1<pt<5.2 GeV/c",kRunsToPlot,0.,kRunsToPlot);
    hrmsDCA5->SetLineWidth(2);
    hrmsDCA5->SetLineColor(kOrange+1);
    hrmsDCA5->SetMarkerColor(kOrange+1);
    hrmsDCA5->SetMarkerStyle(22);
    
    TH1F *hmDCA10 = new TH1F("hmDCA10","mean DCA for TPCITS tracks, 7.0<pt<8.8 GeV/c",kRunsToPlot,0.,kRunsToPlot);
    hmDCA10->SetLineWidth(2);
    hmDCA10->SetLineColor(kMagenta+2);
    hmDCA10->SetMarkerColor(kMagenta+2);
    hmDCA10->SetMarkerStyle(23);
    
    TH1F *hrmsDCA10 = new TH1F("hrmsDCA10","RMS DCA for TPCITS tracks, 7.0<pt<8.8 GeV/c",kRunsToPlot,0.,kRunsToPlot);
    hrmsDCA10->SetLineWidth(2);
    hrmsDCA10->SetLineColor(kMagenta+2);
    hrmsDCA10->SetMarkerColor(kMagenta+2);
    hrmsDCA10->SetMarkerStyle(23);
   
    TH1F *hmDCAz05 = new TH1F("hmDCAz05","mean DCAz for TPCITS tracks, 0.55<pt<0.56 GeV/c",kRunsToPlot,0.,kRunsToPlot);
    hmDCAz05->SetLineWidth(2);
    hmDCAz05->SetLineColor(kBlue);
    hmDCAz05->SetMarkerColor(kBlue);
    hmDCAz05->SetMarkerStyle(20);
    
    TH1F *hrmsDCAz05 = new TH1F("hrmsDCAz05","RMS DCAz for TPCITS tracks, 0.55<pt<0.56 GeV/c",kRunsToPlot,0.,kRunsToPlot);
    hrmsDCAz05->SetLineWidth(2);
    hrmsDCAz05->SetLineColor(kBlue);
    hrmsDCAz05->SetMarkerColor(kBlue);
    hrmsDCAz05->SetMarkerStyle(20);
    
    TH1F *hmDCAz1 = new TH1F("hmDCAz1","mean DCAz for TPCITS tracks, 1.05<pt<1.07 GeV/c",kRunsToPlot,0.,kRunsToPlot);
    hmDCAz1->SetLineWidth(2);
    hmDCAz1->SetLineColor(kRed);
    hmDCAz1->SetMarkerColor(kRed);
    hmDCAz1->SetMarkerStyle(21);
    
    TH1F *hrmsDCAz1 = new TH1F("hrmsDCAz1","RMS DCAz for TPCITS tracks, 1.05<pt<1.07 GeV/c",kRunsToPlot,0.,kRunsToPlot);
    hrmsDCAz1->SetLineWidth(2);
    hrmsDCAz1->SetLineColor(kRed);
    hrmsDCAz1->SetMarkerColor(kRed);
    hrmsDCAz1->SetMarkerStyle(21);
    
    TH1F *hmDCAz5 = new TH1F("hmDCAz5","mean DCAz for TPCITS tracks, 4.1<pt<5.2 GeV/c",kRunsToPlot,0.,kRunsToPlot);
    hmDCAz5->SetLineWidth(2);
    hmDCAz5->SetLineColor(kOrange+1);
    hmDCAz5->SetMarkerColor(kOrange+1);
    hmDCAz5->SetMarkerStyle(22);
    
    TH1F *hrmsDCAz5 = new TH1F("hrmsDCAz5","RMS DCAz for TPCITS tracks, 4.1<pt<5.2 GeV/c",kRunsToPlot,0.,kRunsToPlot);
    hrmsDCAz5->SetLineWidth(2);
    hrmsDCAz5->SetLineColor(kOrange+1);
    hrmsDCAz5->SetMarkerColor(kOrange+1);
    hrmsDCAz5->SetMarkerStyle(22);
    
    TH1F *hmDCAz10 = new TH1F("hmDCAz10","mean DCAz for TPCITS tracks, 7.0<pt<8.8 GeV/c",kRunsToPlot,0.,kRunsToPlot);
    hmDCAz10->SetLineWidth(2);
    hmDCAz10->SetLineColor(kMagenta+2);
    hmDCAz10->SetMarkerColor(kMagenta+2);
    hmDCAz10->SetMarkerStyle(23);
    
    TH1F *hrmsDCAz10 = new TH1F("hrmsDCAz10","RMS DCAz for TPCITS tracks, 7.0<pt<8.8 GeV/c",kRunsToPlot,0.,kRunsToPlot);
    hrmsDCAz10->SetLineWidth(2);
    hrmsDCAz10->SetLineColor(kMagenta+2);
    hrmsDCAz10->SetMarkerColor(kMagenta+2);
    hrmsDCAz10->SetMarkerStyle(23);

    
    for(Int_t i=0;i<kRunsToPlot;i++){
        
        ntDCA->GetEvent(myIndex[i]);
        
        hmDCA05->SetBinContent(i+1,mdca05);
        hmDCA05->SetBinError(i+1,errmdca05);
        hmDCA05->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunDCA));
        
        hrmsDCA05->SetBinContent(i+1,rmsdca05);
        hrmsDCA05->SetBinError(i+1,errrmsdca05);
        hrmsDCA05->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunDCA));
        
        hmDCA1->SetBinContent(i+1,mdca1);
        hmDCA1->SetBinError(i+1,errmdca1);
        hmDCA1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunDCA));
        
        hrmsDCA1->SetBinContent(i+1,rmsdca1);
        hrmsDCA1->SetBinError(i+1,errrmsdca1);
        hrmsDCA1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunDCA));
        
        hmDCA5->SetBinContent(i+1,mdca5);
        hmDCA5->SetBinError(i+1,errmdca5);
        hmDCA5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunDCA));
        
        hrmsDCA5->SetBinContent(i+1,rmsdca5);
        hrmsDCA5->SetBinError(i+1,errrmsdca5);
        hrmsDCA5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunDCA));

        hmDCA10->SetBinContent(i+1,mdca10);
        hmDCA10->SetBinError(i+1,errmdca10);
        hmDCA10->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunDCA));
        
        hrmsDCA10->SetBinContent(i+1,rmsdca10);
        hrmsDCA10->SetBinError(i+1,errrmsdca10);
        hrmsDCA10->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunDCA));
//
        hmDCAz05->SetBinContent(i+1,mdcaz05);
        hmDCAz05->SetBinError(i+1,errmdcaz05);
        hmDCAz05->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunDCA));
        
        hrmsDCAz05->SetBinContent(i+1,rmsdcaz05);
        hrmsDCAz05->SetBinError(i+1,errrmsdcaz05);
        hrmsDCAz05->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunDCA));
        
        hmDCAz1->SetBinContent(i+1,mdcaz1);
        hmDCAz1->SetBinError(i+1,errmdcaz1);
        hmDCAz1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunDCA));
        
        hrmsDCAz1->SetBinContent(i+1,rmsdcaz1);
        hrmsDCAz1->SetBinError(i+1,errrmsdcaz1);
        hrmsDCAz1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunDCA));
        
        hmDCAz5->SetBinContent(i+1,mdcaz5);
        hmDCAz5->SetBinError(i+1,errmdcaz5);
        hmDCAz5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunDCA));
        
        hrmsDCAz5->SetBinContent(i+1,rmsdcaz5);
        hrmsDCAz5->SetBinError(i+1,errrmsdcaz5);
        hrmsDCAz5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunDCA));
        
        hmDCAz10->SetBinContent(i+1,mdcaz10);
        hmDCAz10->SetBinError(i+1,errmdcaz10);
        hmDCAz10->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunDCA));
        
        hrmsDCAz10->SetBinContent(i+1,rmsdcaz10);
        hrmsDCAz10->SetBinError(i+1,errrmsdcaz10);
        hrmsDCAz10->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrunDCA));

    }
    //-----------------------------------

    
 //-----------------------------------

  gStyle->SetOptStat(0);

//
    TCanvas* cMI;
    if(histoTrackMI3->GetEntries()>0){
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
/*
    TCanvas* cMInorm=new TCanvas("cMInorm"," Reconstructed Track With points in ITS");
    histoTrackMI3norm->Draw();
    histoTrackMI3norm->SetLineColor(1);
    histoTrackMI3norm->SetMarkerStyle(20);
    histoTrackMI3norm->Draw();
    histoTrackMI3norm->SetMinimum(0.);
    histoTrackMI3norm->SetMaximum(1.05);
    histoTrackMI4norm->SetLineColor(2);
    histoTrackMI4norm->SetMarkerColor(2);
    histoTrackMI4norm->SetMarkerStyle(22);
    histoTrackMI4norm->Draw("same");
    histoTrackMI1norm->SetLineColor(kGray+1);
    histoTrackMI1norm->SetMarkerColor(kGray+1);
    histoTrackMI1norm->SetMarkerStyle(24);
    histoTrackMI1norm->Draw("same");
    histoTrackMI2norm->SetLineColor(kGray+2);
    histoTrackMI2norm->SetMarkerColor(kGray+2);
    histoTrackMI2norm->SetMarkerStyle(26);
    histoTrackMI2norm->Draw("same");
    histoTrackMI5norm->SetLineColor(4);
    histoTrackMI5norm->SetMarkerColor(4);
    histoTrackMI5norm->SetMarkerStyle(29);
    histoTrackMI5norm->Draw("same");
    histoTrackMI6norm->SetLineColor(kBlue+1);
    histoTrackMI6norm->SetMarkerColor(kBlue+1);
    histoTrackMI6norm->SetMarkerStyle(30);
    histoTrackMI6norm->SetTitle("Fraction of rec. tracks with points in ITS norm.");
    histoTrackMI6norm->Draw("same");
    histoTrackMI3norm->GetYaxis()->SetTitle("Norm. frac. of Global Tracks w/ Points in ITS Layers");
    histoTrackMI3norm->GetXaxis()->SetTitle("run number");
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
    cMInorm->SaveAs("TrackPointsMI_norm_trend.pdf");
    //    pdfFileNames+=" TrackPointsMI_trend.pdf";
    cMInorm->Update();
*/
    //
    TCanvas* cSA;
    if(histoTrackSA3->GetEntries()>0){
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
    
//
    TCanvas* cev;
    if(histoEvwSDD->GetEntries()>0){
        cev=new TCanvas("cev","Fraction of Events with SDD in Trigger Cluster");
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
    }
        //

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
    
    TCanvas* cfrac;
    if(histoFracDead3->GetEntries()>0){
        cfrac=new TCanvas("cfrac","Fraction of SDD modules ON",900,900);
        cfrac->Divide(1,3);
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
    histoFracDead4->SetMarkerStyle(20);
    histoFracDead4->SetMarkerColor(kAzure+1);
    histoFracDead4->SetLineColor(kAzure+1);
    histoFracDead4->GetYaxis()->SetRangeUser(0.,1.2);
    histoFracDead4->Draw("same");
    histoFracDead4->GetYaxis()->SetTitle("Fraction of Modules ON");
    histoFracDead4->GetXaxis()->SetTitle("run number");
    TLatex* tf4=new TLatex(0.2,0.5,"SDD modules ON - Layer 4 (total: 176)");
    tf4->SetNDC();
    tf4->SetTextColor(1);
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
    hVarSDD1->SetMarkerStyle(20);
    hVarSDD1->SetMarkerColor(kOrange+1);
    hVarSDD1->SetLineColor(kOrange+1);
    hVarSDD1->GetYaxis()->SetRangeUser(-0.2,0.5);
    hVarSDD1->Draw();
    hVarSDD1->GetYaxis()->SetTitle("#Delta N/N Modules ON");
    hVarSDD1->GetXaxis()->SetTitle("run number");
    TLatex* tf3_1=new TLatex(0.2,0.40,"#DeltaN/N -- Layer 3");
    tf3_1->SetNDC();
    tf3_1->SetTextColor(1);
    tf3_1->Draw();
    hVarSDD2->SetMarkerStyle(20);
    hVarSDD2->SetMarkerColor(kAzure+1);
    hVarSDD2->SetLineColor(kAzure+1);
    hVarSDD2->GetYaxis()->SetRangeUser(-0.2,0.5);
    hVarSDD2->Draw("same");
    hVarSDD2->GetYaxis()->SetTitle("#Delta N/N Modules ON");
    hVarSDD2->GetXaxis()->SetTitle("run number");
    TLatex* tf4_1=new TLatex(0.2,0.65,"#DeltaN/N +0.2 -- Layer 4");
    tf4_1->SetNDC();
    tf4_1->SetTextColor(1);
    tf4_1->Draw();
        legSDD->Draw();
    
    
    cfrac->cd(3);
    hFlagSDD1->SetMarkerStyle(20);
    hFlagSDD1->SetMarkerColor(kOrange+1);
    hFlagSDD1->SetLineColor(kOrange+1);
    hFlagSDD1->GetYaxis()->SetRangeUser(0.,2.5);
    hFlagSDD1->Draw();
    hFlagSDD1->GetYaxis()->SetTitle("SDD1 alarm flag");
    hFlagSDD1->GetXaxis()->SetTitle("run number");
    TLatex* tf3_2=new TLatex(0.2,0.30,"Status flag Layer 3: 1=OK, 0.5=ALARM");
    tf3_2->SetNDC();
    tf3_2->SetTextColor(1);
    tf3_2->Draw();
    
    hFlagSDD2->SetMarkerStyle(20);
    hFlagSDD2->SetMarkerColor(kAzure+1);
    hFlagSDD2->SetLineColor(kAzure+1);
    hFlagSDD2->GetYaxis()->SetRangeUser(0.,2.5);
    hFlagSDD2->Draw("same");
    hFlagSDD2->GetYaxis()->SetTitle("SDD2 alarm flag");
    hFlagSDD2->GetXaxis()->SetTitle("run number");
    TLatex* tf4_2=new TLatex(0.2,0.80,"Status flag Layer 4: 2=OK, 1.5=ALARM");
    tf4_2->SetNDC();
    tf4_2->SetTextColor(1);
    tf4_2->Draw();
        legSDD->Draw();

    cfrac->SaveAs("SDDmodulesON_trend.pdf");
//    pdfFileNames+=" SDDmodulesON_trend.pdf";
    cfrac->Update();
    }
    
    //
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
  leg2->SetFillStyle(1001);
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
  leg2b->SetFillStyle(1001);
  leg2b->Draw();
  TLatex* tc2=new TLatex(0.2,0.85,"SDD and SSD charge in different layers");
  tc2->SetNDC();
  tc2->SetTextColor(1);
  tc2->Draw();
  c2->SaveAs("SDD_SSD_drift_charge_trend.pdf");  
//  pdfFileNames+=" SDD_SSD_drift_charge_trend.pdf";
  c2->Update();
    }
    
    TCanvas* c2b;
    if(hFlagChR5->GetEntries()>0){
        c2b=new TCanvas("c2b","SDD & SSD ALARM FLAGS",1200,800);
        c2b->Divide(2,2);
    c2b->cd(1);
    hFlagminTime->SetLineColor(2);
    hFlagminTime->SetMarkerColor(2);
    hFlagminTime->SetMarkerStyle(20);
    hFlagminTime->SetMarkerSize(0.5);
    hFlagminTime->Draw();
    hFlagminTime->SetMinimum(0);
    hFlagminTime->SetMaximum(2.5);
    hFlagminTime->GetYaxis()->SetTitle("Minimum & Mean Drift Time Alarms");
    hFlagminTime->GetXaxis()->SetTitle("run number");
    hFlagmeanTime->SetLineColor(4);
    hFlagmeanTime->SetMarkerColor(4);
    hFlagmeanTime->SetMarkerSize(0.5);
    hFlagmeanTime->SetMarkerStyle(21);
    hFlagmeanTime->Draw("same");
    TLatex* td1=new TLatex(0.12,0.45,"SDD: 490 ns<Min Time<510 ns: OK=1, ALARM=0.5");
    td1->SetNDC();
    td1->SetTextColor(1);
    td1->Draw();
    TLatex* td1a=new TLatex(0.12,0.8,"SDD: 3180 ns<Mean Time<3240 ns: OK=2, ALARM=1.5");
    td1a->SetNDC();
    td1a->SetTextColor(1);
    td1a->Draw();
    c2b->cd(2);
    hFlagdEdx3->SetMarkerColor(1);
    hFlagdEdx3->SetMarkerStyle(20);
    hFlagdEdx3->SetLineColor(1);
    hFlagdEdx3->SetMarkerSize(0.5);
    hFlagdEdx3->Draw();
    hFlagdEdx3->SetMinimum(0);
    hFlagdEdx3->SetMaximum(2.5);
    hFlagdEdx3->GetYaxis()->SetTitle("dEdx MPV Layer3 & Layer4 Alarms");
    hFlagdEdx3->GetXaxis()->SetTitle("run number");
    hFlagdEdx4->SetMarkerColor(kBlue);
    hFlagdEdx4->SetMarkerStyle(23);
    hFlagdEdx4->SetLineColor(kBlue);
    hFlagdEdx4->SetMarkerSize(0.5);
    hFlagdEdx4->Draw("same");
    TLatex* td1b=new TLatex(0.12,0.45,"SDD1: 83<dEdx3<85: OK=1, ALARM=0.5");
    td1b->SetNDC();
    td1b->SetTextColor(1);
    td1b->Draw();
    TLatex* td1c=new TLatex(0.12,0.8,"SDD2: 83<dEdx4<85: OK=2, ALARM=1.5");
    td1c->SetNDC();
    td1c->SetTextColor(1);
    td1c->Draw();
    
    c2b->cd(3);
    hFlagChR5->SetMarkerColor(kMagenta+2);
    hFlagChR5->SetMarkerStyle(20);
    hFlagChR5->SetLineColor(kMagenta+2);
    hFlagChR5->SetMarkerSize(0.5);
    hFlagChR5->Draw();
    hFlagChR5->SetMinimum(0);
    hFlagChR5->SetMaximum(2.5);
    hFlagChR5->GetYaxis()->SetTitle("Charge Ratio Layer5 & Layer6 Alarms");
    hFlagChR5->GetXaxis()->SetTitle("run number");
    hFlagChR6->SetMarkerColor(9);
    hFlagChR6->SetMarkerStyle(22);
    hFlagChR6->SetLineColor(9);
    hFlagChR6->SetMarkerSize(0.5);
    hFlagChR6->Draw("same");
    TLatex* td1bb=new TLatex(0.12,0.45,"SSD: -0.01<Charge Ratio L5<0.01: OK=1, ALARM=0.5");
    td1bb->SetNDC();
    td1bb->SetTextColor(1);
    td1bb->Draw();
    TLatex* td1cb=new TLatex(0.12,0.8,"SSD: -0.01<Charge Ratio L6<0.01: OK=2, ALARM=1.5");
    td1cb->SetNDC();
    td1cb->SetTextColor(1);
    td1cb->Draw();
    
    c2b->cd(4);
    hFlagdEdx5->SetMarkerColor(kMagenta+2);
    hFlagdEdx5->SetMarkerStyle(20);
    hFlagdEdx5->SetLineColor(kMagenta+2);
    hFlagdEdx5->SetMarkerSize(0.5);
    hFlagdEdx5->Draw();
    hFlagdEdx5->SetMinimum(0);
    hFlagdEdx5->SetMaximum(2.5);
    hFlagdEdx5->GetYaxis()->SetTitle("dEdx MPV Layer5 & Layer6 Alarms");
    hFlagdEdx5->GetXaxis()->SetTitle("run number");
    hFlagdEdx6->SetMarkerColor(9);
    hFlagdEdx6->SetMarkerStyle(22);
    hFlagdEdx6->SetLineColor(9);
    hFlagdEdx6->SetMarkerSize(0.5);
    hFlagdEdx6->Draw("same");
    TLatex* td1b2=new TLatex(0.12,0.45,"SSD: 82<dEdx5<83: OK=1, ALARM=0.5");
    td1b2->SetNDC();
    td1b2->SetTextColor(1);
    td1b2->Draw();
    TLatex* td1c2=new TLatex(0.12,0.8,"SSD: 82<dEdx6<83: OK=2, ALARM=1.5");
    td1c2->SetNDC();
    td1c2->SetTextColor(1);
    td1c2->Draw();
    
    c2b->SaveAs("SDDSSD_alarm_trend.pdf");
    //  pdfFileNames+=" SDD_SSD_drift_charge_trend.pdf";
    c2b->Update();
    }
    
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
    
    TCanvas *c8b;
    if(hVarSSD1n->GetEntries()>0){
        c8b=new TCanvas("c8b","Automatic QA variables for SSD",900,900);
        c8b->Divide(1,2);
    c8b->cd(1);
    hVarSSD1n->SetMinimum(-0.1);
    hVarSSD1n->SetMaximum(0.5);
    hVarSSD1n->GetXaxis()->SetTitle("run number");
    hVarSSD1n->Draw();
        hVarSSD1n->SetTitle("");
    hVarSSD1p->SetLineColor(kMagenta+2);
    hVarSSD1p->SetMinimum(0);
    hVarSSD1p->SetMaximum(0.2);
    hVarSSD1p->Draw("same");
    hVarSSD2n->Draw("same");
    hVarSSD2p->Draw("same");
    TLegend* legBad3txt=new TLegend(0.8,0.75,0.98,0.95);
     TLegendEntry* enttxt=legBad3txt->AddEntry(hVarSSD1n,"Layer5 n-side","PL");
    enttxt->SetTextColor(hVarSSD1n->GetMarkerColor());
    enttxt=legBad3txt->AddEntry(hVarSSD1p,"Layer5 p-side + 0.1","PL");
    enttxt->SetTextColor(hVarSSD1p->GetMarkerColor());
    enttxt=legBad3txt->AddEntry(hVarSSD2n,"Layer6 n-side + 0.2","PL");
    enttxt->SetTextColor(hVarSSD2n->GetMarkerColor());
    enttxt=legBad3txt->AddEntry(hVarSSD2p,"Layer6 p-side + 0.3","PL");
    enttxt->SetTextColor(hVarSSD2p->GetMarkerColor());
    legBad3txt->SetFillStyle(1001);
    legBad3txt->Draw();
    TLatex* tcbad5a=new TLatex(0.2,0.85,"#DeltaN/N bad n/p strips Layer 5 & 6");
    tcbad5a->SetNDC();
    tcbad5a->SetTextColor(1);
    tcbad5a->Draw();
    
    c8b->cd(2);
    hFlagSSD1n->SetMinimum(0);
    hFlagSSD1n->SetMaximum(2.5);
    hFlagSSD1n->GetXaxis()->SetTitle("run number");
    hFlagSSD1n->Draw();
        hFlagSSD1n->SetTitle("");
    hFlagSSD1p->SetLineColor(kMagenta+2);
    hFlagSSD1p->SetMinimum(0);
    hFlagSSD1p->SetMaximum(0.2);
    hFlagSSD1p->Draw("same");
    hFlagSSD2n->Draw("same");
    hFlagSSD2p->Draw("same");
        TLegend* legBad3=new TLegend(0.8,0.80,0.98,1.00);
        TLegendEntry* ent=legBad3->AddEntry(hVarSSD1n,"Layer5 n-side","PL");
        ent->SetTextColor(hVarSSD1n->GetMarkerColor());
        ent=legBad3->AddEntry(hVarSSD1p,"Layer5 p-side","PL");
        ent->SetTextColor(hVarSSD1p->GetMarkerColor());
        ent=legBad3->AddEntry(hVarSSD2n,"Layer6 n-side","PL");
        ent->SetTextColor(hVarSSD2n->GetMarkerColor());
        ent=legBad3->AddEntry(hVarSSD2p,"Layer6 p-side","PL");
        ent->SetTextColor(hVarSSD2p->GetMarkerColor());
        legBad3->SetFillStyle(1001);
        legBad3->Draw();
    legBad3->Draw();
    TLatex* tcbad5b=new TLatex(0.2,0.85,"SSD Layer 5 & 6 n/p strips alarm flag: tolerance = 1%");
    tcbad5b->SetNDC();
    tcbad5b->SetTextColor(1);
    tcbad5b->Draw();
    
    c8b->SaveAs("SSD_auto_trend.pdf");
    c8b->Update();
    }
    
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
    
/*
    TCanvas *cpt02norm=new TCanvas("cpt02norm","TPC-ITS matching efficiency norm.",1200,1000);
    cpt02norm->Divide(3,2);
    cpt02norm->cd(1);
    hEffoneSPDPt02norm->SetMinimum(0);
    hEffoneSPDPt02norm->Draw();
    hEff456Pt02norm->Draw("same");
    hEffoneSPDPt02norm->GetYaxis()->SetRangeUser(0,1.1);
    TLegend* lpt02a=new TLegend(0.9,0.8,1,1);
    lpt02a->AddEntry(hEff456Pt02norm,"4+5+6 cls","l");
    lpt02a->AddEntry(hEffoneSPDPt02norm,">=1SPD + any","l");
    lpt02a->Draw("same");
    tpc1->Draw();
    
    cpt02norm->cd(2);
    hEffoneSPDPt1norm->SetMinimum(0);
    hEffoneSPDPt1norm->Draw();
    hEff456Pt1norm->Draw("same");
    hEffoneSPDPt1norm->GetYaxis()->SetRangeUser(0,1.1);
    TLegend* lpt1a=new TLegend(0.9,0.8,1,1);
    lpt1a->AddEntry(hEff456Pt1norm,"4+5+6 cls","l");
    lpt1a->AddEntry(hEffoneSPDPt1norm,">=1SPD + any","l");
    lpt1a->Draw("same");
    tpc2->Draw();
    
    cpt02norm->cd(3);
    hEffoneSPDPt10norm->SetMinimum(0);
    hEffoneSPDPt10norm->Draw();
    hEff456Pt10norm->Draw("same");
    hEffoneSPDPt10norm->GetYaxis()->SetRangeUser(0,1.1);
    
    TLegend* lpt10a=new TLegend(0.9,0.8,1,1);
    lpt10a->AddEntry(hEff456Pt10norm,"4+5+6 cls","l");
    lpt10a->AddEntry(hEffoneSPDPt10norm,">=1SPD + any","l");
    lpt10a->Draw("same");
    tpc3->Draw();
    
    // Plot TOF
    cpt02norm->cd(4);
    hEffTOFoneSPDPt02norm->SetMinimum(0);
    hEffTOFoneSPDPt02norm->Draw();
    hEffTOF456Pt02norm->Draw("same");
    hEffTOFoneSPDPt02norm->GetYaxis()->SetRangeUser(0,1.1);
    
    TLegend* lpt0TOF2a=new TLegend(0.9,0.8,1,1);
    lpt0TOF2a->AddEntry(hEffTOF456Pt02norm,"4+5+6 cls","l");
    lpt0TOF2a->AddEntry(hEffTOFoneSPDPt02norm,">=1SPD + any","l");
    lpt0TOF2a->Draw("same");
    tof1->Draw();
    
    
    cpt02norm->cd(5);
    hEffTOFoneSPDPt1norm->SetMinimum(0);
    hEffTOFoneSPDPt1norm->Draw();
    hEffTOF456Pt1norm->Draw("same");
    hEffTOFoneSPDPt1norm->GetYaxis()->SetRangeUser(0,1.1);
    
    TLegend* lpt0TOF1a=new TLegend(0.9,0.8,1,1);
    lpt0TOF1a->AddEntry(hEffTOF456Pt1norm,"4+5+6 cls","l");
    lpt0TOF1a->AddEntry(hEffTOFoneSPDPt1norm,">=1SPD + any","l");
    lpt0TOF1a->Draw("same");
    tof2->Draw();
    
    cpt02norm->cd(6);
    
    hEffTOFoneSPDPt10norm->SetMinimum(0);
    hEffTOFoneSPDPt10norm->Draw();
    hEffTOF456Pt10norm->Draw("same");
    hEffTOFoneSPDPt10norm->GetYaxis()->SetRangeUser(0,1.1);
    
    TLegend* lpt0TOF10a=new TLegend(0.9,0.8,1,1);
    lpt0TOF10a->AddEntry(hEffTOF456Pt10norm,"4+5+6 cls","l");
    lpt0TOF10a->AddEntry(hEffTOFoneSPDPt10norm,">=1SPD + any","l");
    lpt0TOF10a->Draw("same");
    tof3->Draw();
    
    cpt02norm->SaveAs("TPCTOFMatch_norm_trend.pdf");
    //    pdfFileNames+=" TPCTOFMatch_trend.pdf";
*/
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
    
    
    //-----------------------------------
    TCanvas *cPixel;
    if(hFracSPD1->GetEntries()>0){
        cPixel=new TCanvas("cPixel","SPD on",800,900);
        cPixel->Divide(1,3);
    cPixel->cd(1);
    hFracSPD1->SetMaximum(1.2);
    hFracSPD1->SetMinimum(0);
        hFracSPD1->SetTitle("");
    hFracSPD1->Draw("p");
    hFracSPD2->Draw("same,p");
    
    TLegend* lSPD=new TLegend(0.8,0.8,1,1);
    lSPD->AddEntry(hFracSPD1,"Frac. SPD1 ON","Pl");
    lSPD->AddEntry(hFracSPD2,"Frac. SPD2 ON","Pl");
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
    hVarSPD1->SetMaximum(0.3);
    hVarSPD1->SetMinimum(-0.1);
    hVarSPD1->Draw("p");
        hVarSPD1->SetTitle("");
    TLatex* tf4a_1=new TLatex(0.2,0.35,"#DeltaN/N -- Layer 1");
    tf4a_1->SetNDC();
    tf4a_1->SetTextColor(1);
    tf4a_1->Draw();
    hVarSPD2->Draw("same,p");
    TLatex* tf4a_2=new TLatex(0.2,0.75,"#DeltaN/N +0.2 -- Layer 2");
    tf4a_2->SetNDC();
    tf4a_2->SetTextColor(1);
    tf4a_2->Draw();
        TLegend* lSPD1=new TLegend(0.85,0.8,0.98,1);
        lSPD1->AddEntry(hFracSPD1,"SPD1","Pl");
        lSPD1->AddEntry(hFracSPD2,"SPD2 + 0.2","Pl");
        lSPD1->Draw();
    
    cPixel->cd(3);
    hFlagSPD1->SetMaximum(2.5);
    hFlagSPD1->SetMinimum(-0.0);
    hFlagSPD1->Draw("p");
        hFlagSPD1->SetTitle("");
    TLatex* tf4a_3=new TLatex(0.2,0.30,"Status flag Layer 1: 1=OK, 0.5=ALARM");
    tf4a_3->SetNDC();
    tf4a_3->SetTextColor(1);
    tf4a_3->Draw();
    hFlagSPD2->Draw("same,p");
    TLatex* tf4a_4=new TLatex(0.2,0.80,"Status flag Layer 2: 2=OK, 1.5=ALARM");
    tf4a_4->SetNDC();
    tf4a_4->SetTextColor(1);
    tf4a_4->Draw();
        TLegend* lSPD1a=new TLegend(0.9,0.8,0.98,1);
        lSPD1a->AddEntry(hFracSPD1,"SPD1","Pl");
        lSPD1a->AddEntry(hFracSPD2,"SPD2","Pl");
        lSPD1a->Draw();
    
    cPixel->SaveAs("Pixel_trend.pdf");
    //  pdfFileNames+=" Pixel_trend.pdf";
    cPixel->Update();
    }
    
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
    Float_t chi2TPCITS, chi2ITSpureSA;
    
    Float_t occ_phi_1_1,occ_phi_1_2,occ_phi_1_3,occ_phi_1_4,occ_phi_1_5,occ_phi_1_6,occ_phi_1_7,occ_phi_1_8,occ_phi_1_9,occ_phi_1_10;
    Float_t occ_phi_1_11,occ_phi_1_12,occ_phi_1_13,occ_phi_1_14,occ_phi_1_15,occ_phi_1_16,occ_phi_1_17,occ_phi_1_18,occ_phi_1_19,occ_phi_1_20;
    Float_t occ_phi_1_21,occ_phi_1_22,occ_phi_1_23,occ_phi_1_24,occ_phi_1_25,occ_phi_1_26,occ_phi_1_27,occ_phi_1_28,occ_phi_1_29,occ_phi_1_30;
    Float_t occ_phi_1_31,occ_phi_1_32,occ_phi_1_33,occ_phi_1_34,occ_phi_1_35,occ_phi_1_36,occ_phi_1_37,occ_phi_1_38,occ_phi_1_39,occ_phi_1_40;
    Float_t occ_eta_1_1,occ_eta_1_2;
    
    Float_t occ_phi_2_1,occ_phi_2_2,occ_phi_2_3,occ_phi_2_4,occ_phi_2_5,occ_phi_2_6,occ_phi_2_7,occ_phi_2_8,occ_phi_2_9,occ_phi_2_10;
    Float_t occ_phi_2_11,occ_phi_2_12,occ_phi_2_13,occ_phi_2_14,occ_phi_2_15,occ_phi_2_16,occ_phi_2_17,occ_phi_2_18,occ_phi_2_19,occ_phi_2_20;
    Float_t occ_phi_2_21,occ_phi_2_22,occ_phi_2_23,occ_phi_2_24,occ_phi_2_25,occ_phi_2_26,occ_phi_2_27,occ_phi_2_28,occ_phi_2_29,occ_phi_2_30;
    Float_t occ_phi_2_31,occ_phi_2_32,occ_phi_2_33,occ_phi_2_34,occ_phi_2_35,occ_phi_2_36,occ_phi_2_37,occ_phi_2_38,occ_phi_2_39,occ_phi_2_40;
    Float_t occ_eta_2_1,occ_eta_2_2;
   
    Float_t occ_phi_3_1,occ_phi_3_2,occ_phi_3_3,occ_phi_3_4,occ_phi_3_5,occ_phi_3_6,occ_phi_3_7,occ_phi_3_8,occ_phi_3_9,occ_phi_3_10;
    Float_t occ_phi_3_11,occ_phi_3_12,occ_phi_3_13,occ_phi_3_14,occ_phi_3_15,occ_phi_3_16,occ_phi_3_17,occ_phi_3_18,occ_phi_3_19,occ_phi_3_20;
    Float_t occ_phi_3_21,occ_phi_3_22,occ_phi_3_23,occ_phi_3_24,occ_phi_3_25,occ_phi_3_26,occ_phi_3_27,occ_phi_3_28,occ_phi_3_29,occ_phi_3_30;
    Float_t occ_phi_3_31,occ_phi_3_32,occ_phi_3_33,occ_phi_3_34,occ_phi_3_35,occ_phi_3_36,occ_phi_3_37,occ_phi_3_38,occ_phi_3_39,occ_phi_3_40;
    Float_t occ_eta_3_1,occ_eta_3_2;
    
    Float_t occ_phi_4_1,occ_phi_4_2,occ_phi_4_3,occ_phi_4_4,occ_phi_4_5,occ_phi_4_6,occ_phi_4_7,occ_phi_4_8,occ_phi_4_9,occ_phi_4_10;
    Float_t occ_phi_4_11,occ_phi_4_12,occ_phi_4_13,occ_phi_4_14,occ_phi_4_15,occ_phi_4_16,occ_phi_4_17,occ_phi_4_18,occ_phi_4_19,occ_phi_4_20;
    Float_t occ_phi_4_21,occ_phi_4_22,occ_phi_4_23,occ_phi_4_24,occ_phi_4_25,occ_phi_4_26,occ_phi_4_27,occ_phi_4_28,occ_phi_4_29,occ_phi_4_30;
    Float_t occ_phi_4_31,occ_phi_4_32,occ_phi_4_33,occ_phi_4_34,occ_phi_4_35,occ_phi_4_36,occ_phi_4_37,occ_phi_4_38,occ_phi_4_39,occ_phi_4_40;
    Float_t occ_eta_4_1,occ_eta_4_2;

    Float_t occ_phi_5_1,occ_phi_5_2,occ_phi_5_3,occ_phi_5_4,occ_phi_5_5,occ_phi_5_6,occ_phi_5_7,occ_phi_5_8,occ_phi_5_9,occ_phi_5_10;
    Float_t occ_phi_5_11,occ_phi_5_12,occ_phi_5_13,occ_phi_5_14,occ_phi_5_15,occ_phi_5_16,occ_phi_5_17,occ_phi_5_18,occ_phi_5_19,occ_phi_5_20;
    Float_t occ_phi_5_21,occ_phi_5_22,occ_phi_5_23,occ_phi_5_24,occ_phi_5_25,occ_phi_5_26,occ_phi_5_27,occ_phi_5_28,occ_phi_5_29,occ_phi_5_30;
    Float_t occ_phi_5_31,occ_phi_5_32,occ_phi_5_33,occ_phi_5_34,occ_phi_5_35,occ_phi_5_36,occ_phi_5_37,occ_phi_5_38,occ_phi_5_39,occ_phi_5_40;
    Float_t occ_eta_5_1,occ_eta_5_2;

    Float_t occ_phi_6_1,occ_phi_6_2,occ_phi_6_3,occ_phi_6_4,occ_phi_6_5,occ_phi_6_6,occ_phi_6_7,occ_phi_6_8,occ_phi_6_9,occ_phi_6_10;
    Float_t occ_phi_6_11,occ_phi_6_12,occ_phi_6_13,occ_phi_6_14,occ_phi_6_15,occ_phi_6_16,occ_phi_6_17,occ_phi_6_18,occ_phi_6_19,occ_phi_6_20;
    Float_t occ_phi_6_21,occ_phi_6_22,occ_phi_6_23,occ_phi_6_24,occ_phi_6_25,occ_phi_6_26,occ_phi_6_27,occ_phi_6_28,occ_phi_6_29,occ_phi_6_30;
    Float_t occ_phi_6_31,occ_phi_6_32,occ_phi_6_33,occ_phi_6_34,occ_phi_6_35,occ_phi_6_36,occ_phi_6_37,occ_phi_6_38,occ_phi_6_39,occ_phi_6_40;
    Float_t occ_eta_6_1,occ_eta_6_2;

    
    nt->SetBranchAddress("run",&run);
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
    nt->SetBranchAddress("chi2TPCITS",&chi2TPCITS);
    nt->SetBranchAddress("chi2ITSpureSA",&chi2ITSpureSA);
    nt->SetBranchAddress("occ_eta_1_1",&occ_eta_1_1);
    nt->SetBranchAddress("occ_eta_1_2",&occ_eta_1_2);
    nt->SetBranchAddress("occ_phi_1_1",&occ_phi_1_1);
    nt->SetBranchAddress("occ_phi_1_2",&occ_phi_1_2);
    nt->SetBranchAddress("occ_phi_1_3",&occ_phi_1_3);
    nt->SetBranchAddress("occ_phi_1_4",&occ_phi_1_4);
    nt->SetBranchAddress("occ_phi_1_5",&occ_phi_1_5);
    nt->SetBranchAddress("occ_phi_1_6",&occ_phi_1_6);
    nt->SetBranchAddress("occ_phi_1_7",&occ_phi_1_7);
    nt->SetBranchAddress("occ_phi_1_8",&occ_phi_1_8);
    nt->SetBranchAddress("occ_phi_1_9",&occ_phi_1_9);
    nt->SetBranchAddress("occ_phi_1_10",&occ_phi_1_10);
    nt->SetBranchAddress("occ_phi_1_11",&occ_phi_1_11);
    nt->SetBranchAddress("occ_phi_1_12",&occ_phi_1_12);
    nt->SetBranchAddress("occ_phi_1_13",&occ_phi_1_13);
    nt->SetBranchAddress("occ_phi_1_14",&occ_phi_1_14);
    nt->SetBranchAddress("occ_phi_1_15",&occ_phi_1_15);
    nt->SetBranchAddress("occ_phi_1_16",&occ_phi_1_16);
    nt->SetBranchAddress("occ_phi_1_17",&occ_phi_1_17);
    nt->SetBranchAddress("occ_phi_1_18",&occ_phi_1_18);
    nt->SetBranchAddress("occ_phi_1_19",&occ_phi_1_19);
    nt->SetBranchAddress("occ_phi_1_20",&occ_phi_1_20);
    nt->SetBranchAddress("occ_phi_1_21",&occ_phi_1_21);
    nt->SetBranchAddress("occ_phi_1_22",&occ_phi_1_22);
    nt->SetBranchAddress("occ_phi_1_23",&occ_phi_1_23);
    nt->SetBranchAddress("occ_phi_1_24",&occ_phi_1_24);
    nt->SetBranchAddress("occ_phi_1_25",&occ_phi_1_25);
    nt->SetBranchAddress("occ_phi_1_26",&occ_phi_1_26);
    nt->SetBranchAddress("occ_phi_1_27",&occ_phi_1_27);
    nt->SetBranchAddress("occ_phi_1_28",&occ_phi_1_28);
    nt->SetBranchAddress("occ_phi_1_29",&occ_phi_1_29);
    nt->SetBranchAddress("occ_phi_1_30",&occ_phi_1_30);
    nt->SetBranchAddress("occ_phi_1_31",&occ_phi_1_31);
    nt->SetBranchAddress("occ_phi_1_32",&occ_phi_1_32);
    nt->SetBranchAddress("occ_phi_1_33",&occ_phi_1_33);
    nt->SetBranchAddress("occ_phi_1_34",&occ_phi_1_34);
    nt->SetBranchAddress("occ_phi_1_35",&occ_phi_1_35);
    nt->SetBranchAddress("occ_phi_1_36",&occ_phi_1_36);
    nt->SetBranchAddress("occ_phi_1_37",&occ_phi_1_37);
    nt->SetBranchAddress("occ_phi_1_38",&occ_phi_1_38);
    nt->SetBranchAddress("occ_phi_1_39",&occ_phi_1_39);
    nt->SetBranchAddress("occ_phi_1_40",&occ_phi_1_40);
    
    nt->SetBranchAddress("occ_eta_2_1",&occ_eta_2_1);
    nt->SetBranchAddress("occ_eta_2_2",&occ_eta_2_2);
    nt->SetBranchAddress("occ_phi_2_1",&occ_phi_2_1);
    nt->SetBranchAddress("occ_phi_2_2",&occ_phi_2_2);
    nt->SetBranchAddress("occ_phi_2_3",&occ_phi_2_3);
    nt->SetBranchAddress("occ_phi_2_4",&occ_phi_2_4);
    nt->SetBranchAddress("occ_phi_2_5",&occ_phi_2_5);
    nt->SetBranchAddress("occ_phi_2_6",&occ_phi_2_6);
    nt->SetBranchAddress("occ_phi_2_7",&occ_phi_2_7);
    nt->SetBranchAddress("occ_phi_2_8",&occ_phi_2_8);
    nt->SetBranchAddress("occ_phi_2_9",&occ_phi_2_9);
    nt->SetBranchAddress("occ_phi_2_10",&occ_phi_2_10);
    nt->SetBranchAddress("occ_phi_2_11",&occ_phi_2_11);
    nt->SetBranchAddress("occ_phi_2_12",&occ_phi_2_12);
    nt->SetBranchAddress("occ_phi_2_13",&occ_phi_2_13);
    nt->SetBranchAddress("occ_phi_2_14",&occ_phi_2_14);
    nt->SetBranchAddress("occ_phi_2_15",&occ_phi_2_15);
    nt->SetBranchAddress("occ_phi_2_16",&occ_phi_2_16);
    nt->SetBranchAddress("occ_phi_2_17",&occ_phi_2_17);
    nt->SetBranchAddress("occ_phi_2_18",&occ_phi_2_18);
    nt->SetBranchAddress("occ_phi_2_19",&occ_phi_2_19);
    nt->SetBranchAddress("occ_phi_2_20",&occ_phi_2_20);
    nt->SetBranchAddress("occ_phi_2_21",&occ_phi_2_21);
    nt->SetBranchAddress("occ_phi_2_22",&occ_phi_2_22);
    nt->SetBranchAddress("occ_phi_2_23",&occ_phi_2_23);
    nt->SetBranchAddress("occ_phi_2_24",&occ_phi_2_24);
    nt->SetBranchAddress("occ_phi_2_25",&occ_phi_2_25);
    nt->SetBranchAddress("occ_phi_2_26",&occ_phi_2_26);
    nt->SetBranchAddress("occ_phi_2_27",&occ_phi_2_27);
    nt->SetBranchAddress("occ_phi_2_28",&occ_phi_2_28);
    nt->SetBranchAddress("occ_phi_2_29",&occ_phi_2_29);
    nt->SetBranchAddress("occ_phi_2_30",&occ_phi_2_30);
    nt->SetBranchAddress("occ_phi_2_31",&occ_phi_2_31);
    nt->SetBranchAddress("occ_phi_2_32",&occ_phi_2_32);
    nt->SetBranchAddress("occ_phi_2_33",&occ_phi_2_33);
    nt->SetBranchAddress("occ_phi_2_34",&occ_phi_2_34);
    nt->SetBranchAddress("occ_phi_2_35",&occ_phi_2_35);
    nt->SetBranchAddress("occ_phi_2_36",&occ_phi_2_36);
    nt->SetBranchAddress("occ_phi_2_37",&occ_phi_2_37);
    nt->SetBranchAddress("occ_phi_2_38",&occ_phi_2_38);
    nt->SetBranchAddress("occ_phi_2_39",&occ_phi_2_39);
    nt->SetBranchAddress("occ_phi_2_40",&occ_phi_2_40);

    nt->SetBranchAddress("occ_eta_3_1",&occ_eta_3_1);
    nt->SetBranchAddress("occ_eta_3_2",&occ_eta_3_2);
    nt->SetBranchAddress("occ_phi_3_1",&occ_phi_3_1);
    nt->SetBranchAddress("occ_phi_3_2",&occ_phi_3_2);
    nt->SetBranchAddress("occ_phi_3_3",&occ_phi_3_3);
    nt->SetBranchAddress("occ_phi_3_4",&occ_phi_3_4);
    nt->SetBranchAddress("occ_phi_3_5",&occ_phi_3_5);
    nt->SetBranchAddress("occ_phi_3_6",&occ_phi_3_6);
    nt->SetBranchAddress("occ_phi_3_7",&occ_phi_3_7);
    nt->SetBranchAddress("occ_phi_3_8",&occ_phi_3_8);
    nt->SetBranchAddress("occ_phi_3_9",&occ_phi_3_9);
    nt->SetBranchAddress("occ_phi_3_10",&occ_phi_3_10);
    nt->SetBranchAddress("occ_phi_3_11",&occ_phi_3_11);
    nt->SetBranchAddress("occ_phi_3_12",&occ_phi_3_12);
    nt->SetBranchAddress("occ_phi_3_13",&occ_phi_3_13);
    nt->SetBranchAddress("occ_phi_3_14",&occ_phi_3_14);
    nt->SetBranchAddress("occ_phi_3_15",&occ_phi_3_15);
    nt->SetBranchAddress("occ_phi_3_16",&occ_phi_3_16);
    nt->SetBranchAddress("occ_phi_3_17",&occ_phi_3_17);
    nt->SetBranchAddress("occ_phi_3_18",&occ_phi_3_18);
    nt->SetBranchAddress("occ_phi_3_19",&occ_phi_3_19);
    nt->SetBranchAddress("occ_phi_3_20",&occ_phi_3_20);
    nt->SetBranchAddress("occ_phi_3_21",&occ_phi_3_21);
    nt->SetBranchAddress("occ_phi_3_22",&occ_phi_3_22);
    nt->SetBranchAddress("occ_phi_3_23",&occ_phi_3_23);
    nt->SetBranchAddress("occ_phi_3_24",&occ_phi_3_24);
    nt->SetBranchAddress("occ_phi_3_25",&occ_phi_3_25);
    nt->SetBranchAddress("occ_phi_3_26",&occ_phi_3_26);
    nt->SetBranchAddress("occ_phi_3_27",&occ_phi_3_27);
    nt->SetBranchAddress("occ_phi_3_28",&occ_phi_3_28);
    nt->SetBranchAddress("occ_phi_3_29",&occ_phi_3_29);
    nt->SetBranchAddress("occ_phi_3_30",&occ_phi_3_30);
    nt->SetBranchAddress("occ_phi_3_31",&occ_phi_3_31);
    nt->SetBranchAddress("occ_phi_3_32",&occ_phi_3_32);
    nt->SetBranchAddress("occ_phi_3_33",&occ_phi_3_33);
    nt->SetBranchAddress("occ_phi_3_34",&occ_phi_3_34);
    nt->SetBranchAddress("occ_phi_3_35",&occ_phi_3_35);
    nt->SetBranchAddress("occ_phi_3_36",&occ_phi_3_36);
    nt->SetBranchAddress("occ_phi_3_37",&occ_phi_3_37);
    nt->SetBranchAddress("occ_phi_3_38",&occ_phi_3_38);
    nt->SetBranchAddress("occ_phi_3_39",&occ_phi_3_39);
    nt->SetBranchAddress("occ_phi_3_40",&occ_phi_3_40);

    nt->SetBranchAddress("occ_eta_4_1",&occ_eta_4_1);
    nt->SetBranchAddress("occ_eta_4_2",&occ_eta_4_2);
    nt->SetBranchAddress("occ_phi_4_1",&occ_phi_4_1);
    nt->SetBranchAddress("occ_phi_4_2",&occ_phi_4_2);
    nt->SetBranchAddress("occ_phi_4_3",&occ_phi_4_3);
    nt->SetBranchAddress("occ_phi_4_4",&occ_phi_4_4);
    nt->SetBranchAddress("occ_phi_4_5",&occ_phi_4_5);
    nt->SetBranchAddress("occ_phi_4_6",&occ_phi_4_6);
    nt->SetBranchAddress("occ_phi_4_7",&occ_phi_4_7);
    nt->SetBranchAddress("occ_phi_4_8",&occ_phi_4_8);
    nt->SetBranchAddress("occ_phi_4_9",&occ_phi_4_9);
    nt->SetBranchAddress("occ_phi_4_10",&occ_phi_4_10);
    nt->SetBranchAddress("occ_phi_4_11",&occ_phi_4_11);
    nt->SetBranchAddress("occ_phi_4_12",&occ_phi_4_12);
    nt->SetBranchAddress("occ_phi_4_13",&occ_phi_4_13);
    nt->SetBranchAddress("occ_phi_4_14",&occ_phi_4_14);
    nt->SetBranchAddress("occ_phi_4_15",&occ_phi_4_15);
    nt->SetBranchAddress("occ_phi_4_16",&occ_phi_4_16);
    nt->SetBranchAddress("occ_phi_4_17",&occ_phi_4_17);
    nt->SetBranchAddress("occ_phi_4_18",&occ_phi_4_18);
    nt->SetBranchAddress("occ_phi_4_19",&occ_phi_4_19);
    nt->SetBranchAddress("occ_phi_4_20",&occ_phi_4_20);
    nt->SetBranchAddress("occ_phi_4_21",&occ_phi_4_21);
    nt->SetBranchAddress("occ_phi_4_22",&occ_phi_4_22);
    nt->SetBranchAddress("occ_phi_4_23",&occ_phi_4_23);
    nt->SetBranchAddress("occ_phi_4_24",&occ_phi_4_24);
    nt->SetBranchAddress("occ_phi_4_25",&occ_phi_4_25);
    nt->SetBranchAddress("occ_phi_4_26",&occ_phi_4_26);
    nt->SetBranchAddress("occ_phi_4_27",&occ_phi_4_27);
    nt->SetBranchAddress("occ_phi_4_28",&occ_phi_4_28);
    nt->SetBranchAddress("occ_phi_4_29",&occ_phi_4_29);
    nt->SetBranchAddress("occ_phi_4_30",&occ_phi_4_30);
    nt->SetBranchAddress("occ_phi_4_31",&occ_phi_4_31);
    nt->SetBranchAddress("occ_phi_4_32",&occ_phi_4_32);
    nt->SetBranchAddress("occ_phi_4_33",&occ_phi_4_33);
    nt->SetBranchAddress("occ_phi_4_34",&occ_phi_4_34);
    nt->SetBranchAddress("occ_phi_4_35",&occ_phi_4_35);
    nt->SetBranchAddress("occ_phi_4_36",&occ_phi_4_36);
    nt->SetBranchAddress("occ_phi_4_37",&occ_phi_4_37);
    nt->SetBranchAddress("occ_phi_4_38",&occ_phi_4_38);
    nt->SetBranchAddress("occ_phi_4_39",&occ_phi_4_39);
    nt->SetBranchAddress("occ_phi_4_40",&occ_phi_4_40);

    nt->SetBranchAddress("occ_eta_5_1",&occ_eta_5_1);
    nt->SetBranchAddress("occ_eta_5_2",&occ_eta_5_2);
    nt->SetBranchAddress("occ_phi_5_1",&occ_phi_5_1);
    nt->SetBranchAddress("occ_phi_5_2",&occ_phi_5_2);
    nt->SetBranchAddress("occ_phi_5_3",&occ_phi_5_3);
    nt->SetBranchAddress("occ_phi_5_4",&occ_phi_5_4);
    nt->SetBranchAddress("occ_phi_5_5",&occ_phi_5_5);
    nt->SetBranchAddress("occ_phi_5_6",&occ_phi_5_6);
    nt->SetBranchAddress("occ_phi_5_7",&occ_phi_5_7);
    nt->SetBranchAddress("occ_phi_5_8",&occ_phi_5_8);
    nt->SetBranchAddress("occ_phi_5_9",&occ_phi_5_9);
    nt->SetBranchAddress("occ_phi_5_10",&occ_phi_5_10);
    nt->SetBranchAddress("occ_phi_5_11",&occ_phi_5_11);
    nt->SetBranchAddress("occ_phi_5_12",&occ_phi_5_12);
    nt->SetBranchAddress("occ_phi_5_13",&occ_phi_5_13);
    nt->SetBranchAddress("occ_phi_5_14",&occ_phi_5_14);
    nt->SetBranchAddress("occ_phi_5_15",&occ_phi_5_15);
    nt->SetBranchAddress("occ_phi_5_16",&occ_phi_5_16);
    nt->SetBranchAddress("occ_phi_5_17",&occ_phi_5_17);
    nt->SetBranchAddress("occ_phi_5_18",&occ_phi_5_18);
    nt->SetBranchAddress("occ_phi_5_19",&occ_phi_5_19);
    nt->SetBranchAddress("occ_phi_5_20",&occ_phi_5_20);
    nt->SetBranchAddress("occ_phi_5_21",&occ_phi_5_21);
    nt->SetBranchAddress("occ_phi_5_22",&occ_phi_5_22);
    nt->SetBranchAddress("occ_phi_5_23",&occ_phi_5_23);
    nt->SetBranchAddress("occ_phi_5_24",&occ_phi_5_24);
    nt->SetBranchAddress("occ_phi_5_25",&occ_phi_5_25);
    nt->SetBranchAddress("occ_phi_5_26",&occ_phi_5_26);
    nt->SetBranchAddress("occ_phi_5_27",&occ_phi_5_27);
    nt->SetBranchAddress("occ_phi_5_28",&occ_phi_5_28);
    nt->SetBranchAddress("occ_phi_5_29",&occ_phi_5_29);
    nt->SetBranchAddress("occ_phi_5_30",&occ_phi_5_30);
    nt->SetBranchAddress("occ_phi_5_31",&occ_phi_5_31);
    nt->SetBranchAddress("occ_phi_5_32",&occ_phi_5_32);
    nt->SetBranchAddress("occ_phi_5_33",&occ_phi_5_33);
    nt->SetBranchAddress("occ_phi_5_34",&occ_phi_5_34);
    nt->SetBranchAddress("occ_phi_5_35",&occ_phi_5_35);
    nt->SetBranchAddress("occ_phi_5_36",&occ_phi_5_36);
    nt->SetBranchAddress("occ_phi_5_37",&occ_phi_5_37);
    nt->SetBranchAddress("occ_phi_5_38",&occ_phi_5_38);
    nt->SetBranchAddress("occ_phi_5_39",&occ_phi_5_39);
    nt->SetBranchAddress("occ_phi_5_40",&occ_phi_5_40);

    nt->SetBranchAddress("occ_eta_6_1",&occ_eta_6_1);
    nt->SetBranchAddress("occ_eta_6_2",&occ_eta_6_2);
    nt->SetBranchAddress("occ_phi_6_1",&occ_phi_6_1);
    nt->SetBranchAddress("occ_phi_6_2",&occ_phi_6_2);
    nt->SetBranchAddress("occ_phi_6_3",&occ_phi_6_3);
    nt->SetBranchAddress("occ_phi_6_4",&occ_phi_6_4);
    nt->SetBranchAddress("occ_phi_6_5",&occ_phi_6_5);
    nt->SetBranchAddress("occ_phi_6_6",&occ_phi_6_6);
    nt->SetBranchAddress("occ_phi_6_7",&occ_phi_6_7);
    nt->SetBranchAddress("occ_phi_6_8",&occ_phi_6_8);
    nt->SetBranchAddress("occ_phi_6_9",&occ_phi_6_9);
    nt->SetBranchAddress("occ_phi_6_10",&occ_phi_6_10);
    nt->SetBranchAddress("occ_phi_6_11",&occ_phi_6_11);
    nt->SetBranchAddress("occ_phi_6_12",&occ_phi_6_12);
    nt->SetBranchAddress("occ_phi_6_13",&occ_phi_6_13);
    nt->SetBranchAddress("occ_phi_6_14",&occ_phi_6_14);
    nt->SetBranchAddress("occ_phi_6_15",&occ_phi_6_15);
    nt->SetBranchAddress("occ_phi_6_16",&occ_phi_6_16);
    nt->SetBranchAddress("occ_phi_6_17",&occ_phi_6_17);
    nt->SetBranchAddress("occ_phi_6_18",&occ_phi_6_18);
    nt->SetBranchAddress("occ_phi_6_19",&occ_phi_6_19);
    nt->SetBranchAddress("occ_phi_6_20",&occ_phi_6_20);
    nt->SetBranchAddress("occ_phi_6_21",&occ_phi_6_21);
    nt->SetBranchAddress("occ_phi_6_22",&occ_phi_6_22);
    nt->SetBranchAddress("occ_phi_6_23",&occ_phi_6_23);
    nt->SetBranchAddress("occ_phi_6_24",&occ_phi_6_24);
    nt->SetBranchAddress("occ_phi_6_25",&occ_phi_6_25);
    nt->SetBranchAddress("occ_phi_6_26",&occ_phi_6_26);
    nt->SetBranchAddress("occ_phi_6_27",&occ_phi_6_27);
    nt->SetBranchAddress("occ_phi_6_28",&occ_phi_6_28);
    nt->SetBranchAddress("occ_phi_6_29",&occ_phi_6_29);
    nt->SetBranchAddress("occ_phi_6_30",&occ_phi_6_30);
    nt->SetBranchAddress("occ_phi_6_31",&occ_phi_6_31);
    nt->SetBranchAddress("occ_phi_6_32",&occ_phi_6_32);
    nt->SetBranchAddress("occ_phi_6_33",&occ_phi_6_33);
    nt->SetBranchAddress("occ_phi_6_34",&occ_phi_6_34);
    nt->SetBranchAddress("occ_phi_6_35",&occ_phi_6_35);
    nt->SetBranchAddress("occ_phi_6_36",&occ_phi_6_36);
    nt->SetBranchAddress("occ_phi_6_37",&occ_phi_6_37);
    nt->SetBranchAddress("occ_phi_6_38",&occ_phi_6_38);
    nt->SetBranchAddress("occ_phi_6_39",&occ_phi_6_39);
    nt->SetBranchAddress("occ_phi_6_40",&occ_phi_6_40);

    //    cout << "SPD1 " << NclupSA0 << " errore " << errNclupSA0 << endl;
    
    nr=nt->GetEntries();
    delete []myIndex;
    delete []noRuns;
    myIndex = new Int_t [nr];
    noRuns = new Int_t [nr];
    for(Int_t i=0; i<nr;i++){
        nt->GetEvent(i);
        Int_t intrun = static_cast<Int_t>(run+0.01);
        noRuns[i]=intrun;
    }
    printf("\n ======== PROCESSING ITS SA NTUPLE \n");
//    cout << "nru = " << nru << " nr = " << nr << endl;
//    Int_t kRunsToPlot = RunsToPlot(run1,run2,nr,noRuns2,myIndex2);
    kRunsToPlot = RunsToPlot(run1,run2,nr,noRuns,myIndex);
    
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
    
/*    TH1F *h12norm=new TH1F("h12norm","",kRunsToPlot,-0.5,kRunsToPlot-0.5);
    TH1F *h13norm=new TH1F("h13norm","",kRunsToPlot,-0.5,kRunsToPlot-0.5);
    TH1F *h14norm=new TH1F("h14norm","",kRunsToPlot,-0.5,kRunsToPlot-0.5);
    TH1F *h15norm=new TH1F("h15norm","",kRunsToPlot,-0.5,kRunsToPlot-0.5);
    TH1F *h16norm=new TH1F("h16norm","",kRunsToPlot,-0.5,kRunsToPlot-0.5);
    TH1F *h17norm=new TH1F("h17norm","",kRunsToPlot,-0.5,kRunsToPlot-0.5);
*/
    TH1F *hchi2=new TH1F("hchi2","hchi2",kRunsToPlot,-0.5,kRunsToPlot-0.5);
    hchi2->SetLineWidth(2);
//    hchi2->SetLineColor(kRed);
    hchi2->SetMarkerColor(kRed);
    hchi2->SetMarkerStyle(20);

    TH1F *hchi2SA=new TH1F("hchi2SA","hchi2SA",kRunsToPlot,-0.5,kRunsToPlot-0.5);
    hchi2SA->SetLineWidth(2);
//    hchi2SA->SetLineColor(kBlue);
    hchi2SA->SetMarkerColor(kBlue);
    hchi2SA->SetMarkerStyle(20);

    TH2F *hOccEta1=new TH2F("hOccEta1","hOccEta1",kRunsToPlot,-0.5,kRunsToPlot-0.5,2,0.5,2.5);
    hOccEta1->SetTitle("Layer 1 - #h");
    TH2F *hOccPhi1=new TH2F("hOccPhi1","hOccPhi1",kRunsToPlot,-0.5,kRunsToPlot-0.5,40,0.5,40.5);
    hOccPhi1->SetTitle("Layer 1 - #phi");
    TH2F *hOccEta2=new TH2F("hOccEta2","hOccEta2",kRunsToPlot,-0.5,kRunsToPlot-0.5,2,0.5,2.5);
    hOccEta2->SetTitle("Layer 2 - #h");
    TH2F *hOccPhi2=new TH2F("hOccPhi2","hOccPhi2",kRunsToPlot,-0.5,kRunsToPlot-0.5,40,0.5,40.5);
    hOccPhi2->SetTitle("Layer 2 - #phi");
    TH2F *hOccEta3=new TH2F("hOccEta3","hOccEta3",kRunsToPlot,-0.5,kRunsToPlot-0.5,2,0.5,2.5);
    hOccEta3->SetTitle("Layer 3 - #h");
    TH2F *hOccPhi3=new TH2F("hOccPhi3","hOccPhi3",kRunsToPlot,-0.5,kRunsToPlot-0.5,40,0.5,40.5);
    hOccPhi3->SetTitle("Layer 3 - #phi");
    TH2F *hOccEta4=new TH2F("hOccEta4","hOccEta4",kRunsToPlot,-0.5,kRunsToPlot-0.5,2,0.5,2.5);
    hOccEta4->SetTitle("Layer 4 - #h");
    TH2F *hOccPhi4=new TH2F("hOccPhi4","hOccPhi4",kRunsToPlot,-0.5,kRunsToPlot-0.5,40,0.5,40.5);
    hOccPhi4->SetTitle("Layer 4 - #phi");
    TH2F *hOccEta5=new TH2F("hOccEta5","hOccEta5",kRunsToPlot,-0.5,kRunsToPlot-0.5,2,0.5,2.5);
    hOccEta5->SetTitle("Layer 5 - #h");
    TH2F *hOccPhi5=new TH2F("hOccPhi5","hOccPhi5",kRunsToPlot,-0.5,kRunsToPlot-0.5,40,0.5,40.5);
    hOccPhi5->SetTitle("Layer 5 - #phi");
    TH2F *hOccEta6=new TH2F("hOccEta6","hOccEta6",kRunsToPlot,-0.5,kRunsToPlot-0.5,2,0.5,2.5);
    hOccEta6->SetTitle("Layer 6 - #h");
    TH2F *hOccPhi6=new TH2F("hOccPhi6","hOccPhi6",kRunsToPlot,-0.5,kRunsToPlot-0.5,40,0.5,40.5);
    hOccPhi6->SetTitle("Layer 6 - #phi");

    
    for(Int_t iev=0;iev<kRunsToPlot;iev++){
        nt->GetEvent(myIndex[iev]);
//        cout << "numero ordine run: " << iev << " ITSSA numero di run = " << run << endl;
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
        
/*        h12norm->SetBinContent(iev+1,NclupSA0/hFracSPD1->GetBinContent(iev+1));
        h12norm->SetBinError(iev+1,errNclupSA0/hFracSPD1->GetBinContent(iev+1));
        h13norm->SetBinContent(iev+1,NclupSA1/hFracSPD2->GetBinContent(iev+1));
        h13norm->SetBinError(iev+1,errNclupSA1/hFracSPD2->GetBinContent(iev+1));
        h14norm->SetBinContent(iev+1,NclupSA2/histoFracDead3->GetBinContent(iev+1));
        h14norm->SetBinError(iev+1,errNclupSA2/histoFracDead3->GetBinContent(iev+1));
        h15norm->SetBinContent(iev+1,NclupSA3/histoFracDead4->GetBinContent(iev+1));
        h15norm->SetBinError(iev+1,errNclupSA3/histoFracDead4->GetBinContent(iev+1));
        h16norm->SetBinContent(iev+1,NclupSA4/(1-0.5*(histoFracBadn5->GetBinContent(iev+1)+histoFracBadp5->GetBinContent(iev+1))));
        h16norm->SetBinError(iev+1,errNclupSA4/(1-0.5*(histoFracBadn5->GetBinContent(iev+1)+histoFracBadp5->GetBinContent(iev+1))));
        h17norm->SetBinContent(iev+1,NclupSA5/(1-0.5*(histoFracBadn6->GetBinContent(iev+1)+histoFracBadp6->GetBinContent(iev+1))));
        h17norm->SetBinError(iev+1,errNclupSA5/(1-0.5*(histoFracBadn6->GetBinContent(iev+1)+histoFracBadp6->GetBinContent(iev+1))));
*/
        hchi2->SetBinContent(iev+1,chi2TPCITS);
        hchi2->SetBinError(iev+1,0.01);
        hchi2SA->SetBinContent(iev+1,chi2ITSpureSA);
        hchi2SA->SetBinError(iev+1,0.01);
        
        
        hOccEta1->SetBinContent(iev+1,1,occ_eta_1_1); // first eta bin
        hOccEta1->SetBinContent(iev+1,2,occ_eta_1_2); // second eta bin
        hOccPhi1->SetBinContent(iev+1,1,occ_phi_1_1); // first phi bin
        hOccPhi1->SetBinContent(iev+1,2,occ_phi_1_2); // second phi bin
        hOccPhi1->SetBinContent(iev+1,3,occ_phi_1_3); // third phi bin
        hOccPhi1->SetBinContent(iev+1,4,occ_phi_1_4); // fourth phi bin
        hOccPhi1->SetBinContent(iev+1,5,occ_phi_1_5); // fifth phi bin
        hOccPhi1->SetBinContent(iev+1,6,occ_phi_1_6); // sixth phi bin
        hOccPhi1->SetBinContent(iev+1,7,occ_phi_1_7); // seventh phi bin
        hOccPhi1->SetBinContent(iev+1,8,occ_phi_1_8); // eighth phi bin
        hOccPhi1->SetBinContent(iev+1,9,occ_phi_1_9); // nineth phi bin
        hOccPhi1->SetBinContent(iev+1,10,occ_phi_1_10); // tenth phi bin
        hOccPhi1->SetBinContent(iev+1,11,occ_phi_1_11); // first phi bin
        hOccPhi1->SetBinContent(iev+1,12,occ_phi_1_12); // second phi bin
        hOccPhi1->SetBinContent(iev+1,13,occ_phi_1_13); // third phi bin
        hOccPhi1->SetBinContent(iev+1,14,occ_phi_1_14); // fourth phi bin
        hOccPhi1->SetBinContent(iev+1,15,occ_phi_1_15); // fifth phi bin
        hOccPhi1->SetBinContent(iev+1,16,occ_phi_1_16); // sixth phi bin
        hOccPhi1->SetBinContent(iev+1,17,occ_phi_1_17); // seventh phi bin
        hOccPhi1->SetBinContent(iev+1,18,occ_phi_1_18); // eighth phi bin
        hOccPhi1->SetBinContent(iev+1,19,occ_phi_1_19); // nineth phi bin
        hOccPhi1->SetBinContent(iev+1,20,occ_phi_1_20); // tenth phi bin
        hOccPhi1->SetBinContent(iev+1,21,occ_phi_1_21); // first phi bin
        hOccPhi1->SetBinContent(iev+1,22,occ_phi_1_22); // second phi bin
        hOccPhi1->SetBinContent(iev+1,23,occ_phi_1_23); // third phi bin
        hOccPhi1->SetBinContent(iev+1,24,occ_phi_1_24); // fourth phi bin
        hOccPhi1->SetBinContent(iev+1,25,occ_phi_1_25); // fifth phi bin
        hOccPhi1->SetBinContent(iev+1,26,occ_phi_1_26); // sixth phi bin
        hOccPhi1->SetBinContent(iev+1,27,occ_phi_1_27); // seventh phi bin
        hOccPhi1->SetBinContent(iev+1,28,occ_phi_1_28); // eighth phi bin
        hOccPhi1->SetBinContent(iev+1,29,occ_phi_1_29); // nineth phi bin
        hOccPhi1->SetBinContent(iev+1,30,occ_phi_1_30); // tenth phi bin
        hOccPhi1->SetBinContent(iev+1,31,occ_phi_1_31); // first phi bin
        hOccPhi1->SetBinContent(iev+1,32,occ_phi_1_32); // second phi bin
        hOccPhi1->SetBinContent(iev+1,33,occ_phi_1_33); // third phi bin
        hOccPhi1->SetBinContent(iev+1,34,occ_phi_1_34); // fourth phi bin
        hOccPhi1->SetBinContent(iev+1,35,occ_phi_1_35); // fifth phi bin
        hOccPhi1->SetBinContent(iev+1,36,occ_phi_1_36); // sixth phi bin
        hOccPhi1->SetBinContent(iev+1,37,occ_phi_1_37); // seventh phi bin
        hOccPhi1->SetBinContent(iev+1,38,occ_phi_1_38); // eighth phi bin
        hOccPhi1->SetBinContent(iev+1,39,occ_phi_1_39); // nineth phi bin
        hOccPhi1->SetBinContent(iev+1,40,occ_phi_1_40); // tenth phi bin

        hOccEta2->SetBinContent(iev+1,1,occ_eta_2_1); // first eta bin
        hOccEta2->SetBinContent(iev+1,2,occ_eta_2_2); // second eta bin
        hOccPhi2->SetBinContent(iev+1,1,occ_phi_2_1); // first phi bin
        hOccPhi2->SetBinContent(iev+1,2,occ_phi_2_2); // second phi bin
        hOccPhi2->SetBinContent(iev+1,3,occ_phi_2_3); // third phi bin
        hOccPhi2->SetBinContent(iev+1,4,occ_phi_2_4); // fourth phi bin
        hOccPhi2->SetBinContent(iev+1,5,occ_phi_2_5); // fifth phi bin
        hOccPhi2->SetBinContent(iev+1,6,occ_phi_2_6); // sixth phi bin
        hOccPhi2->SetBinContent(iev+1,7,occ_phi_2_7); // seventh phi bin
        hOccPhi2->SetBinContent(iev+1,8,occ_phi_2_8); // eighth phi bin
        hOccPhi2->SetBinContent(iev+1,9,occ_phi_2_9); // nineth phi bin
        hOccPhi2->SetBinContent(iev+1,10,occ_phi_2_10); // tenth phi bin
        hOccPhi2->SetBinContent(iev+1,11,occ_phi_2_11); // first phi bin
        hOccPhi2->SetBinContent(iev+1,12,occ_phi_2_12); // second phi bin
        hOccPhi2->SetBinContent(iev+1,13,occ_phi_2_13); // third phi bin
        hOccPhi2->SetBinContent(iev+1,14,occ_phi_2_14); // fourth phi bin
        hOccPhi2->SetBinContent(iev+1,15,occ_phi_2_15); // fifth phi bin
        hOccPhi2->SetBinContent(iev+1,16,occ_phi_2_16); // sixth phi bin
        hOccPhi2->SetBinContent(iev+1,17,occ_phi_2_17); // seventh phi bin
        hOccPhi2->SetBinContent(iev+1,18,occ_phi_2_18); // eighth phi bin
        hOccPhi2->SetBinContent(iev+1,19,occ_phi_2_19); // nineth phi bin
        hOccPhi2->SetBinContent(iev+1,20,occ_phi_2_20); // tenth phi bin
        hOccPhi2->SetBinContent(iev+1,21,occ_phi_2_21); // first phi bin
        hOccPhi2->SetBinContent(iev+1,22,occ_phi_2_22); // second phi bin
        hOccPhi2->SetBinContent(iev+1,23,occ_phi_2_23); // third phi bin
        hOccPhi2->SetBinContent(iev+1,24,occ_phi_2_24); // fourth phi bin
        hOccPhi2->SetBinContent(iev+1,25,occ_phi_2_25); // fifth phi bin
        hOccPhi2->SetBinContent(iev+1,26,occ_phi_2_26); // sixth phi bin
        hOccPhi2->SetBinContent(iev+1,27,occ_phi_2_27); // seventh phi bin
        hOccPhi2->SetBinContent(iev+1,28,occ_phi_2_28); // eighth phi bin
        hOccPhi2->SetBinContent(iev+1,29,occ_phi_2_29); // nineth phi bin
        hOccPhi2->SetBinContent(iev+1,30,occ_phi_2_30); // tenth phi bin
        hOccPhi2->SetBinContent(iev+1,31,occ_phi_2_31); // first phi bin
        hOccPhi2->SetBinContent(iev+1,32,occ_phi_2_32); // second phi bin
        hOccPhi2->SetBinContent(iev+1,33,occ_phi_2_33); // third phi bin
        hOccPhi2->SetBinContent(iev+1,34,occ_phi_2_34); // fourth phi bin
        hOccPhi2->SetBinContent(iev+1,35,occ_phi_2_35); // fifth phi bin
        hOccPhi2->SetBinContent(iev+1,36,occ_phi_2_36); // sixth phi bin
        hOccPhi2->SetBinContent(iev+1,37,occ_phi_2_37); // seventh phi bin
        hOccPhi2->SetBinContent(iev+1,38,occ_phi_2_38); // eighth phi bin
        hOccPhi2->SetBinContent(iev+1,39,occ_phi_2_39); // nineth phi bin
        hOccPhi2->SetBinContent(iev+1,40,occ_phi_2_40); // tenth phi bin

        hOccEta3->SetBinContent(iev+1,1,occ_eta_3_1); // first eta bin
        hOccEta3->SetBinContent(iev+1,2,occ_eta_3_2); // second eta bin
        hOccPhi3->SetBinContent(iev+1,1,occ_phi_3_1); // first phi bin
        hOccPhi3->SetBinContent(iev+1,2,occ_phi_3_2); // second phi bin
        hOccPhi3->SetBinContent(iev+1,3,occ_phi_3_3); // third phi bin
        hOccPhi3->SetBinContent(iev+1,4,occ_phi_3_4); // fourth phi bin
        hOccPhi3->SetBinContent(iev+1,5,occ_phi_3_5); // fifth phi bin
        hOccPhi3->SetBinContent(iev+1,6,occ_phi_3_6); // sixth phi bin
        hOccPhi3->SetBinContent(iev+1,7,occ_phi_3_7); // seventh phi bin
        hOccPhi3->SetBinContent(iev+1,8,occ_phi_3_8); // eighth phi bin
        hOccPhi3->SetBinContent(iev+1,9,occ_phi_3_9); // nineth phi bin
        hOccPhi3->SetBinContent(iev+1,10,occ_phi_3_10); // tenth phi bin
        hOccPhi3->SetBinContent(iev+1,11,occ_phi_3_11); // first phi bin
        hOccPhi3->SetBinContent(iev+1,12,occ_phi_3_12); // second phi bin
        hOccPhi3->SetBinContent(iev+1,13,occ_phi_3_13); // third phi bin
        hOccPhi3->SetBinContent(iev+1,14,occ_phi_3_14); // fourth phi bin
        hOccPhi3->SetBinContent(iev+1,15,occ_phi_3_15); // fifth phi bin
        hOccPhi3->SetBinContent(iev+1,16,occ_phi_3_16); // sixth phi bin
        hOccPhi3->SetBinContent(iev+1,17,occ_phi_3_17); // seventh phi bin
        hOccPhi3->SetBinContent(iev+1,18,occ_phi_3_18); // eighth phi bin
        hOccPhi3->SetBinContent(iev+1,19,occ_phi_3_19); // nineth phi bin
        hOccPhi3->SetBinContent(iev+1,20,occ_phi_3_20); // tenth phi bin
        hOccPhi3->SetBinContent(iev+1,21,occ_phi_3_21); // first phi bin
        hOccPhi3->SetBinContent(iev+1,22,occ_phi_3_22); // second phi bin
        hOccPhi3->SetBinContent(iev+1,23,occ_phi_3_23); // third phi bin
        hOccPhi3->SetBinContent(iev+1,24,occ_phi_3_24); // fourth phi bin
        hOccPhi3->SetBinContent(iev+1,25,occ_phi_3_25); // fifth phi bin
        hOccPhi3->SetBinContent(iev+1,26,occ_phi_3_26); // sixth phi bin
        hOccPhi3->SetBinContent(iev+1,27,occ_phi_3_27); // seventh phi bin
        hOccPhi3->SetBinContent(iev+1,28,occ_phi_3_28); // eighth phi bin
        hOccPhi3->SetBinContent(iev+1,29,occ_phi_3_29); // nineth phi bin
        hOccPhi3->SetBinContent(iev+1,30,occ_phi_3_30); // tenth phi bin
        hOccPhi3->SetBinContent(iev+1,31,occ_phi_3_31); // first phi bin
        hOccPhi3->SetBinContent(iev+1,32,occ_phi_3_32); // second phi bin
        hOccPhi3->SetBinContent(iev+1,33,occ_phi_3_33); // third phi bin
        hOccPhi3->SetBinContent(iev+1,34,occ_phi_3_34); // fourth phi bin
        hOccPhi3->SetBinContent(iev+1,35,occ_phi_3_35); // fifth phi bin
        hOccPhi3->SetBinContent(iev+1,36,occ_phi_3_36); // sixth phi bin
        hOccPhi3->SetBinContent(iev+1,37,occ_phi_3_37); // seventh phi bin
        hOccPhi3->SetBinContent(iev+1,38,occ_phi_3_38); // eighth phi bin
        hOccPhi3->SetBinContent(iev+1,39,occ_phi_3_39); // nineth phi bin
        hOccPhi3->SetBinContent(iev+1,40,occ_phi_3_40); // tenth phi bin

        hOccEta4->SetBinContent(iev+1,1,occ_eta_4_1); // first eta bin
        hOccEta4->SetBinContent(iev+1,2,occ_eta_4_2); // second eta bin
        hOccPhi4->SetBinContent(iev+1,1,occ_phi_4_1); // first phi bin
        hOccPhi4->SetBinContent(iev+1,2,occ_phi_4_2); // second phi bin
        hOccPhi4->SetBinContent(iev+1,3,occ_phi_4_3); // third phi bin
        hOccPhi4->SetBinContent(iev+1,4,occ_phi_4_4); // fourth phi bin
        hOccPhi4->SetBinContent(iev+1,5,occ_phi_4_5); // fifth phi bin
        hOccPhi4->SetBinContent(iev+1,6,occ_phi_4_6); // sixth phi bin
        hOccPhi4->SetBinContent(iev+1,7,occ_phi_4_7); // seventh phi bin
        hOccPhi4->SetBinContent(iev+1,8,occ_phi_4_8); // eighth phi bin
        hOccPhi4->SetBinContent(iev+1,9,occ_phi_4_9); // nineth phi bin
        hOccPhi4->SetBinContent(iev+1,10,occ_phi_4_10); // tenth phi bin
        hOccPhi4->SetBinContent(iev+1,11,occ_phi_4_11); // first phi bin
        hOccPhi4->SetBinContent(iev+1,12,occ_phi_4_12); // second phi bin
        hOccPhi4->SetBinContent(iev+1,13,occ_phi_4_13); // third phi bin
        hOccPhi4->SetBinContent(iev+1,14,occ_phi_4_14); // fourth phi bin
        hOccPhi4->SetBinContent(iev+1,15,occ_phi_4_15); // fifth phi bin
        hOccPhi4->SetBinContent(iev+1,16,occ_phi_4_16); // sixth phi bin
        hOccPhi4->SetBinContent(iev+1,17,occ_phi_4_17); // seventh phi bin
        hOccPhi4->SetBinContent(iev+1,18,occ_phi_4_18); // eighth phi bin
        hOccPhi4->SetBinContent(iev+1,19,occ_phi_4_19); // nineth phi bin
        hOccPhi4->SetBinContent(iev+1,20,occ_phi_4_20); // tenth phi bin
        hOccPhi4->SetBinContent(iev+1,21,occ_phi_4_21); // first phi bin
        hOccPhi4->SetBinContent(iev+1,22,occ_phi_4_22); // second phi bin
        hOccPhi4->SetBinContent(iev+1,23,occ_phi_4_23); // third phi bin
        hOccPhi4->SetBinContent(iev+1,24,occ_phi_4_24); // fourth phi bin
        hOccPhi4->SetBinContent(iev+1,25,occ_phi_4_25); // fifth phi bin
        hOccPhi4->SetBinContent(iev+1,26,occ_phi_4_26); // sixth phi bin
        hOccPhi4->SetBinContent(iev+1,27,occ_phi_4_27); // seventh phi bin
        hOccPhi4->SetBinContent(iev+1,28,occ_phi_4_28); // eighth phi bin
        hOccPhi4->SetBinContent(iev+1,29,occ_phi_4_29); // nineth phi bin
        hOccPhi4->SetBinContent(iev+1,30,occ_phi_4_30); // tenth phi bin
        hOccPhi4->SetBinContent(iev+1,31,occ_phi_4_31); // first phi bin
        hOccPhi4->SetBinContent(iev+1,32,occ_phi_4_32); // second phi bin
        hOccPhi4->SetBinContent(iev+1,33,occ_phi_4_33); // third phi bin
        hOccPhi4->SetBinContent(iev+1,34,occ_phi_4_34); // fourth phi bin
        hOccPhi4->SetBinContent(iev+1,35,occ_phi_4_35); // fifth phi bin
        hOccPhi4->SetBinContent(iev+1,36,occ_phi_4_36); // sixth phi bin
        hOccPhi4->SetBinContent(iev+1,37,occ_phi_4_37); // seventh phi bin
        hOccPhi4->SetBinContent(iev+1,38,occ_phi_4_38); // eighth phi bin
        hOccPhi4->SetBinContent(iev+1,39,occ_phi_4_39); // nineth phi bin
        hOccPhi4->SetBinContent(iev+1,40,occ_phi_4_40); // tenth phi bin

        hOccEta5->SetBinContent(iev+1,1,occ_eta_2_1); // first eta bin
        hOccEta5->SetBinContent(iev+1,2,occ_eta_5_2); // second eta bin
        hOccPhi5->SetBinContent(iev+1,1,occ_phi_5_1); // first phi bin
        hOccPhi5->SetBinContent(iev+1,2,occ_phi_5_2); // second phi bin
        hOccPhi5->SetBinContent(iev+1,3,occ_phi_5_3); // third phi bin
        hOccPhi5->SetBinContent(iev+1,4,occ_phi_5_4); // fourth phi bin
        hOccPhi5->SetBinContent(iev+1,5,occ_phi_5_5); // fifth phi bin
        hOccPhi5->SetBinContent(iev+1,6,occ_phi_5_6); // sixth phi bin
        hOccPhi5->SetBinContent(iev+1,7,occ_phi_5_7); // seventh phi bin
        hOccPhi5->SetBinContent(iev+1,8,occ_phi_5_8); // eighth phi bin
        hOccPhi5->SetBinContent(iev+1,9,occ_phi_5_9); // nineth phi bin
        hOccPhi5->SetBinContent(iev+1,10,occ_phi_5_10); // tenth phi bin
        hOccPhi5->SetBinContent(iev+1,11,occ_phi_5_11); // first phi bin
        hOccPhi5->SetBinContent(iev+1,12,occ_phi_5_12); // second phi bin
        hOccPhi5->SetBinContent(iev+1,13,occ_phi_5_13); // third phi bin
        hOccPhi5->SetBinContent(iev+1,14,occ_phi_5_14); // fourth phi bin
        hOccPhi5->SetBinContent(iev+1,15,occ_phi_5_15); // fifth phi bin
        hOccPhi5->SetBinContent(iev+1,16,occ_phi_5_16); // sixth phi bin
        hOccPhi5->SetBinContent(iev+1,17,occ_phi_5_17); // seventh phi bin
        hOccPhi5->SetBinContent(iev+1,18,occ_phi_5_18); // eighth phi bin
        hOccPhi5->SetBinContent(iev+1,19,occ_phi_5_19); // nineth phi bin
        hOccPhi5->SetBinContent(iev+1,20,occ_phi_5_20); // tenth phi bin
        hOccPhi5->SetBinContent(iev+1,21,occ_phi_5_21); // first phi bin
        hOccPhi5->SetBinContent(iev+1,22,occ_phi_5_22); // second phi bin
        hOccPhi5->SetBinContent(iev+1,23,occ_phi_5_23); // third phi bin
        hOccPhi5->SetBinContent(iev+1,24,occ_phi_5_24); // fourth phi bin
        hOccPhi5->SetBinContent(iev+1,25,occ_phi_5_25); // fifth phi bin
        hOccPhi5->SetBinContent(iev+1,26,occ_phi_5_26); // sixth phi bin
        hOccPhi5->SetBinContent(iev+1,27,occ_phi_5_27); // seventh phi bin
        hOccPhi5->SetBinContent(iev+1,28,occ_phi_5_28); // eighth phi bin
        hOccPhi5->SetBinContent(iev+1,29,occ_phi_5_29); // nineth phi bin
        hOccPhi5->SetBinContent(iev+1,30,occ_phi_5_30); // tenth phi bin
        hOccPhi5->SetBinContent(iev+1,31,occ_phi_5_31); // first phi bin
        hOccPhi5->SetBinContent(iev+1,32,occ_phi_5_32); // second phi bin
        hOccPhi5->SetBinContent(iev+1,33,occ_phi_5_33); // third phi bin
        hOccPhi5->SetBinContent(iev+1,34,occ_phi_5_34); // fourth phi bin
        hOccPhi5->SetBinContent(iev+1,35,occ_phi_5_35); // fifth phi bin
        hOccPhi5->SetBinContent(iev+1,36,occ_phi_5_36); // sixth phi bin
        hOccPhi5->SetBinContent(iev+1,37,occ_phi_5_37); // seventh phi bin
        hOccPhi5->SetBinContent(iev+1,38,occ_phi_5_38); // eighth phi bin
        hOccPhi5->SetBinContent(iev+1,39,occ_phi_5_39); // nineth phi bin
        hOccPhi5->SetBinContent(iev+1,40,occ_phi_5_40); // tenth phi bin

        hOccEta6->SetBinContent(iev+1,1,occ_eta_6_1); // first eta bin
        hOccEta6->SetBinContent(iev+1,2,occ_eta_6_2); // second eta bin
        hOccPhi6->SetBinContent(iev+1,1,occ_phi_6_1); // first phi bin
        hOccPhi6->SetBinContent(iev+1,2,occ_phi_6_2); // second phi bin
        hOccPhi6->SetBinContent(iev+1,3,occ_phi_6_3); // third phi bin
        hOccPhi6->SetBinContent(iev+1,4,occ_phi_6_4); // fourth phi bin
        hOccPhi6->SetBinContent(iev+1,5,occ_phi_6_5); // fifth phi bin
        hOccPhi6->SetBinContent(iev+1,6,occ_phi_6_6); // sixth phi bin
        hOccPhi6->SetBinContent(iev+1,7,occ_phi_6_7); // seventh phi bin
        hOccPhi6->SetBinContent(iev+1,8,occ_phi_6_8); // eighth phi bin
        hOccPhi6->SetBinContent(iev+1,9,occ_phi_6_9); // nineth phi bin
        hOccPhi6->SetBinContent(iev+1,10,occ_phi_6_10); // tenth phi bin
        hOccPhi6->SetBinContent(iev+1,11,occ_phi_6_11); // first phi bin
        hOccPhi6->SetBinContent(iev+1,12,occ_phi_6_12); // second phi bin
        hOccPhi6->SetBinContent(iev+1,13,occ_phi_6_13); // third phi bin
        hOccPhi6->SetBinContent(iev+1,14,occ_phi_6_14); // fourth phi bin
        hOccPhi6->SetBinContent(iev+1,15,occ_phi_6_15); // fifth phi bin
        hOccPhi6->SetBinContent(iev+1,16,occ_phi_6_16); // sixth phi bin
        hOccPhi6->SetBinContent(iev+1,17,occ_phi_6_17); // seventh phi bin
        hOccPhi6->SetBinContent(iev+1,18,occ_phi_6_18); // eighth phi bin
        hOccPhi6->SetBinContent(iev+1,19,occ_phi_6_19); // nineth phi bin
        hOccPhi6->SetBinContent(iev+1,20,occ_phi_6_20); // tenth phi bin
        hOccPhi6->SetBinContent(iev+1,21,occ_phi_6_21); // first phi bin
        hOccPhi6->SetBinContent(iev+1,22,occ_phi_6_22); // second phi bin
        hOccPhi6->SetBinContent(iev+1,23,occ_phi_6_23); // third phi bin
        hOccPhi6->SetBinContent(iev+1,24,occ_phi_6_24); // fourth phi bin
        hOccPhi6->SetBinContent(iev+1,25,occ_phi_6_25); // fifth phi bin
        hOccPhi6->SetBinContent(iev+1,26,occ_phi_6_26); // sixth phi bin
        hOccPhi6->SetBinContent(iev+1,27,occ_phi_6_27); // seventh phi bin
        hOccPhi6->SetBinContent(iev+1,28,occ_phi_6_28); // eighth phi bin
        hOccPhi6->SetBinContent(iev+1,29,occ_phi_6_29); // nineth phi bin
        hOccPhi6->SetBinContent(iev+1,30,occ_phi_6_30); // tenth phi bin
        hOccPhi6->SetBinContent(iev+1,31,occ_phi_6_31); // first phi bin
        hOccPhi6->SetBinContent(iev+1,32,occ_phi_6_32); // second phi bin
        hOccPhi6->SetBinContent(iev+1,33,occ_phi_6_33); // third phi bin
        hOccPhi6->SetBinContent(iev+1,34,occ_phi_6_34); // fourth phi bin
        hOccPhi6->SetBinContent(iev+1,35,occ_phi_6_35); // fifth phi bin
        hOccPhi6->SetBinContent(iev+1,36,occ_phi_6_36); // sixth phi bin
        hOccPhi6->SetBinContent(iev+1,37,occ_phi_6_37); // seventh phi bin
        hOccPhi6->SetBinContent(iev+1,38,occ_phi_6_38); // eighth phi bin
        hOccPhi6->SetBinContent(iev+1,39,occ_phi_6_39); // nineth phi bin
        hOccPhi6->SetBinContent(iev+1,40,occ_phi_6_40); // tenth phi bin

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
/*        h12norm->GetXaxis()->SetBinLabel(iev+1,Form("%d",(Int_t)run));
        h13norm->GetXaxis()->SetBinLabel(iev+1,Form("%d",(Int_t)run));
        h14norm->GetXaxis()->SetBinLabel(iev+1,Form("%d",(Int_t)run));
        h15norm->GetXaxis()->SetBinLabel(iev+1,Form("%d",(Int_t)run));
        h16norm->GetXaxis()->SetBinLabel(iev+1,Form("%d",(Int_t)run));
        h17norm->GetXaxis()->SetBinLabel(iev+1,Form("%d",(Int_t)run));
 */
        hchi2->GetXaxis()->SetBinLabel(iev+1,Form("%d",(Int_t)run));
        hchi2SA->GetXaxis()->SetBinLabel(iev+1,Form("%d",(Int_t)run));
        hOccEta1->GetXaxis()->SetBinLabel(iev+1,Form("%d",(Int_t)run));
        hOccEta2->GetXaxis()->SetBinLabel(iev+1,Form("%d",(Int_t)run));
        hOccEta3->GetXaxis()->SetBinLabel(iev+1,Form("%d",(Int_t)run));
        hOccEta4->GetXaxis()->SetBinLabel(iev+1,Form("%d",(Int_t)run));
        hOccEta5->GetXaxis()->SetBinLabel(iev+1,Form("%d",(Int_t)run));
        hOccEta6->GetXaxis()->SetBinLabel(iev+1,Form("%d",(Int_t)run));
        hOccPhi1->GetXaxis()->SetBinLabel(iev+1,Form("%d",(Int_t)run));
        hOccPhi2->GetXaxis()->SetBinLabel(iev+1,Form("%d",(Int_t)run));
        hOccPhi3->GetXaxis()->SetBinLabel(iev+1,Form("%d",(Int_t)run));
        hOccPhi4->GetXaxis()->SetBinLabel(iev+1,Form("%d",(Int_t)run));
        hOccPhi5->GetXaxis()->SetBinLabel(iev+1,Form("%d",(Int_t)run));
        hOccPhi6->GetXaxis()->SetBinLabel(iev+1,Form("%d",(Int_t)run));
       
        //    Printf("%f   %f   %f",ITSA[0],ITSA[1],ITSA[2]);
    }
//    h0->Print("all");
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
/*    h12norm->GetXaxis()->SetTitle("run number");
    h12norm->GetYaxis()->SetTitle("Norm. fraction of pureSA tracks w/ cls in ITS layers");
    h13norm->GetXaxis()->SetTitle("run number");
    h14norm->GetXaxis()->SetTitle("run number");
    h15norm->GetXaxis()->SetTitle("run number");
    h16norm->GetXaxis()->SetTitle("run number");
    h17norm->GetXaxis()->SetTitle("run number");
 */
    hchi2->GetXaxis()->SetTitle("run number");
    hchi2->GetYaxis()->SetTitle("#chi^{2} for TPCITS tracks");
    hchi2SA->GetXaxis()->SetTitle("run number");
    hchi2SA->GetYaxis()->SetTitle("#chi^{2} for ITSpureSA tracks");
    hOccEta1->GetXaxis()->SetTitle("run number");
    hOccEta1->GetYaxis()->SetTitle("#eta: 1=-1<#eta<0; 2=0<#eta<1");
    hOccPhi1->GetXaxis()->SetTitle("run number");
    hOccPhi1->GetYaxis()->SetTitle("#phi, 9 degrees/bin");
    hOccEta2->GetXaxis()->SetTitle("run number");
    hOccEta2->GetYaxis()->SetTitle("#eta: 1=-1<#eta<0; 2=0<#eta<1");
    hOccPhi2->GetXaxis()->SetTitle("run number");
    hOccPhi2->GetYaxis()->SetTitle("#phi, 9 degrees/bin");
    hOccEta3->GetXaxis()->SetTitle("run number");
    hOccEta3->GetYaxis()->SetTitle("#eta: 1=-1<#eta<0; 2=0<#eta<1");
    hOccPhi3->GetXaxis()->SetTitle("run number");
    hOccPhi3->GetYaxis()->SetTitle("#phi, 9 degrees/bin");
    hOccEta4->GetXaxis()->SetTitle("run number");
    hOccEta4->GetYaxis()->SetTitle("#eta: 1=-1<#eta<0; 2=0<#eta<1");
    hOccPhi4->GetXaxis()->SetTitle("run number");
    hOccPhi4->GetYaxis()->SetTitle("#phi, 9 degrees/bin");
    hOccEta5->GetXaxis()->SetTitle("run number");
    hOccEta5->GetYaxis()->SetTitle("#eta: 1=-1<#eta<0; 2=0<#eta<1");
    hOccPhi5->GetXaxis()->SetTitle("run number");
    hOccPhi5->GetYaxis()->SetTitle("#phi, 9 degrees/bin");
    hOccEta6->GetXaxis()->SetTitle("run number");
    hOccEta6->GetYaxis()->SetTitle("#eta: 1=-1<#eta<0; 2=0<#eta<1");
    hOccPhi6->GetXaxis()->SetTitle("run number");
    hOccPhi6->GetYaxis()->SetTitle("#phi, 9 degrees/bin");

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
//    h12norm->SetMinimum(0.);
//    h12norm->SetMaximum(1.2);
    
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

    TString name1="File_tracksITS.root";
    TFile *f=new TFile(name1,"RECREATE");
    f->cd();
    c->Write();
    c2a->Write();
    c3->Write();
    f->Close();
    
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
    
    TCanvas* cSA3;
    if(h11->GetEntries()>0 && h11->GetBinContent(1)>0.){
        cSA3=new TCanvas("cSA3","mean pion Pt ITSsa");
        h11->SetMarkerStyle(20);
    h11->SetMarkerColor(4);
    h11->Draw();
        TLatex* tl = new TLatex(0.2,0.85,"pureSA tracks");
    tl->Draw();
    cSA3->SaveAs("piPt_SA_trend.pdf");
    //    pdfFileNames+=" piPt_SA_trend.pdf";
    }
    
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
/*
    TCanvas* cSA4n=new TCanvas("cSA4n","Normalized fraction of tracks with clusters in layers");
    h12norm->GetYaxis()->SetRange(0.,1.);
    h12norm->SetMarkerStyle(24);
    h12norm->SetMarkerColor(kGray+1);
    h12norm->SetLineColor(kGray+1);
    h12norm->Draw();
    //    tl->Draw();
    h13norm->SetMarkerStyle(26);
    h13norm->SetMarkerColor(kGray+2);
    h13norm->SetLineColor(kGray+2);
    h13norm->Draw("same");
    h14norm->SetMarkerStyle(20);
    h14norm->SetMarkerColor(1);
    h14norm->SetLineColor(1);
    h14norm->Draw("same");
    h15norm->SetMarkerStyle(22);
    h15norm->SetMarkerColor(2);
    h15norm->SetLineColor(2);
    h15norm->Draw("same");
    h16norm->SetMarkerStyle(29);
    h16norm->SetMarkerColor(4);
    h16norm->SetLineColor(4);
    h16norm->Draw("same");
    h17norm->SetMarkerStyle(30);
    h17norm->SetMarkerColor(kBlue+1);
    h17norm->SetLineColor(kBlue+1);
    h17norm->Draw("same");
    
    legpSA->SetFillStyle(0);
    legpSA->Draw();
    
    cSA4n->SaveAs("Frac_track_SA_norm_trend.pdf");
    //    pdfFileNames+=" Frac_track_SA_trend.pdf";
 */

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

    // plot eta-phi cluster occupancy per layer

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
    
    // plot ITS PID histos
    
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
    
    // plot TPCITS tracks DCA histos
    
    TCanvas* cDCA;
    if(hmDCA05->GetEntries()>0 && hmDCA05->GetBinContent(1)>-500.){
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
    
    // plot TPCITS tracks DCAz histos
    
    TCanvas* cDCAz;
    if(hmDCAz05->GetEntries()>0 && hmDCAz05->GetBinContent(1)>-500.){
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

    // create output file order pdf and root canvases
    // histos nel file .pdf nell'ordine voluto
    if(hVx->GetEntries()>0) pdfFileNames+=" Vertex_trend.pdf";
    if(histoTrackMI3->GetEntries()>0 && histoTrackMI3->GetBinContent(1)>0.) pdfFileNames+=" TrackPointsMI_trend.pdf";
//    pdfFileNames+=" TrackPointsMI_norm_trend.pdf";
    if(h12->GetEntries()>0 && h12->GetBinContent(1)>0.) pdfFileNames+=" Frac_track_SA_trend.pdf";
//    pdfFileNames+=" Frac_track_SA_norm_trend.pdf";
    if(histoEvwSDD->GetEntries()>0 && histoEvwSDD->GetBinContent(1)>0.) pdfFileNames+=" NoFast_trend.pdf";
    if(histonEvents->GetEntries()>0 && histonEvents->GetBinContent(1)>0.) pdfFileNames+=" RunEvents_trend.pdf";
    if(histoFracDead3->GetEntries()>0) pdfFileNames+=" SDDmodulesON_trend.pdf";
    if(histodEdxLay5->GetEntries()>0) pdfFileNames+=" SDD_SSD_drift_charge_trend.pdf";
    if(hFlagChR5->GetEntries()>0) pdfFileNames+=" SDDSSD_alarm_trend.pdf";
    if(histoFracBadn5->GetEntries()>0) pdfFileNames+=" SSD_BadStripsFrac_trend.pdf";
    if(hVarSSD1n->GetEntries()>0) pdfFileNames+=" SSD_auto_trend.pdf";
    if(histoChargeRatioLay5->GetEntries()>0) pdfFileNames+=" SSD_chargeratio_trend.pdf";
    if(hNsig02->GetEntries()>0 && hNsig02->GetBinContent(1)<100.) pdfFileNames+=" ITSPID_trend.pdf";
    if(hEff6Pt02->GetEntries()>0 && hEff6Pt02->GetBinContent(1)>0.) pdfFileNames+=" TPCTOFMatch_trend.pdf";
//    pdfFileNames+=" TPCTOFmatch_norm_trend.pdf";
    if(hpileupSPD->GetEntries()>0) pdfFileNames+=" Pileup_trend.pdf";
    if(hNumPilVtx->GetEntries()>0) pdfFileNames+=" PileupVtx_trend.pdf";
    if(hFracSPD1->GetEntries()>0) pdfFileNames+=" Pixel_trend.pdf";
    if(h0->GetEntries()>0) pdfFileNames+=" ITSsa_trend.pdf";
    if(h3->GetEntries()>0) pdfFileNames+=" ITSTPC_trend.pdf";
    if(h8->GetEntries()>0) pdfFileNames+=" tracks_ratio_trend.pdf";
    if(h9->GetEntries()>0 && h9->GetBinContent(1)>0.) pdfFileNames+=" meanclu_SA_trend.pdf";
    if(h10->GetEntries()>0 && h10->GetBinContent(1)>0.) pdfFileNames+=" dedx4_3_SA_trend.pdf";
    if(h11->GetEntries()>0 && h11->GetBinContent(1)>0.) pdfFileNames+=" piPt_SA_trend.pdf";
    if(hchi2->GetEntries()>0 && hchi2->GetBinContent(1)>0.) pdfFileNames+=" Chi2_tracks_trend.pdf";
    if(hmDCA05->GetEntries()>0 && hmDCA05->GetBinContent(1)>-500.) pdfFileNames+=" DCAtracks_trend.pdf";
    if(hmDCAz05->GetEntries()>0 && hmDCAz05->GetBinContent(1)>-500.) pdfFileNames+=" DCAztracks_trend.pdf";
    if(hOccEta1->GetEntries()>0 && hOccEta1->GetMean()>0) pdfFileNames+=" Layer1_eta_phi_trend.pdf";
    if(hOccEta2->GetEntries()>0 && hOccEta2->GetMean()>0) pdfFileNames+=" Layer2_eta_phi_trend.pdf";
    if(hOccEta3->GetEntries()>0 && hOccEta3->GetMean()>0) pdfFileNames+=" Layer3_eta_phi_trend.pdf";
    if(hOccEta4->GetEntries()>0 && hOccEta4->GetMean()>0) pdfFileNames+=" Layer4_eta_phi_trend.pdf";
    if(hOccEta5->GetEntries()>0 && hOccEta5->GetMean()>0) pdfFileNames+=" Layer5_eta_phi_trend.pdf";
    if(hOccEta6->GetEntries()>0 && hOccEta6->GetMean()>0) pdfFileNames+=" Layer6_eta_phi_trend.pdf";

    
    // merge the pdf files
  TString command("gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=merged");
  command=command+"ITS_trend_2016.pdf "+pdfFileNames;
  gSystem->Exec(command.Data());
  printf(" Merging the pdf file:  %s \n",command.Data());
  delete [] myIndex;
  delete [] noRuns;

     TString nameO=ntupleFileName;
     nameO+="_File.root";
    TFile *f2=new TFile("nameO.root","RECREATE");
    

    if(hVx->GetEntries()>0) cVertexDisto->Write();
    if(histoTrackMI3->GetEntries()>0 && histoTrackMI3->GetBinContent(1)>0.) cMI->Write();
//    cMInorm->Write();
    if(histoTrackSA3->GetEntries()>0) cSA->Write();
    if(histoTrackClu3->GetEntries()>0)c5->Write();
    if(histoEvwSDD->GetEntries()>0 && histoEvwSDD->GetBinContent(1)>0.) cev->Write();
    if(histonEvents->GetEntries()>0 && histonEvents->GetBinContent(1)>0.) c22->Write();
    if(histoFracDead3->GetEntries()>0) cfrac->Write();
    if(histodEdxLay5->GetEntries()>0) c2->Write();
    if(hFlagChR5->GetEntries()>0) c2b->Write();
    if(histoFracBadn5->GetEntries()>0) c8->Write();
    if(hVarSSD1n->GetEntries()>0) c8b->Write();
    if(histoChargeRatioLay5->GetEntries()>0) c7->Write();
    if(hEff6Pt02->GetEntries()>0 && hEff6Pt02->GetBinContent(1)>0.) cpt02->Write();
    if(hpileupSPD->GetEntries()>0) cPileUp->Write();
    if(hNumPilVtx->GetEntries()>0) cpu->Write();
    if(hFracSPD1->GetEntries()>0) cPixel->Write();
    if(h0->GetEntries()>0) c->Write();
    if(h3->GetEntries()>0) c2a->Write();
    if(h8->GetEntries()>0) c3->Write();
    if(h9->GetEntries()>0 && h9->GetBinContent(1)>0.) cSAb->Write();
    if(h10->GetEntries()>0 && h10->GetBinContent(1)>0.) cSA2->Write();
    if(h11->GetEntries()>0 && h11->GetBinContent(1)>0.) cSA3->Write();
    if(h12->GetEntries()>0 && h12->GetBinContent(1)>0.) cSA4->Write();
    if(hchi2->GetEntries()>0 && hchi2->GetBinContent(1)>0.) cChi2->Write();
    if(hmDCA05->GetEntries()>0 && hmDCA05->GetBinContent(1)>-500.) cDCA->Write();
    if(hmDCAz05->GetEntries()>0 && hmDCAz05->GetBinContent(1)>-500.) cDCAz->Write();
    if(hOccEta1->GetEntries()>0 && hOccEta1->GetMean()>0) cEtaPhi1->Write();
    if(hOccEta2->GetEntries()>0 && hOccEta2->GetMean()>0) cEtaPhi2->Write();
    if(hOccEta3->GetEntries()>0 && hOccEta3->GetMean()>0) cEtaPhi3->Write();
    if(hOccEta4->GetEntries()>0 && hOccEta4->GetMean()>0) cEtaPhi4->Write();
    if(hOccEta5->GetEntries()>0 && hOccEta5->GetMean()>0) cEtaPhi5->Write();
    if(hOccEta6->GetEntries()>0 && hOccEta6->GetMean()>0) cEtaPhi6->Write();
    if(hNsig02->GetEntries()>0 && hNsig02->GetBinContent(1)<100.) cPID->Write();
 
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
//____________________________________________________________________________
void AliITSQAtrend(TString runListFile,TString ntupleFileName){

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

  const Int_t nVariablesMatchingTOF=60; // 53 actually

  TNtuple* ntmatchingTOF=new TNtuple("ntmatchingTOF","MatchingTOF Efficiency","nrun:FracSPD1:errFracSPD1:FracSPD2:errFracSPD2:Eff6Pt02:errEff6Pt02:Eff6Pt1:errEff6Pt1:Eff6Pt10:errEff6Pt10:Eff5Pt02:errEff5Pt02:Eff5Pt1:errEff5Pt1:Eff5Pt10:errEff5Pt10:Eff4Pt02:errEff4Pt02:Eff4Pt1:errEff4Pt1:Eff4Pt10:errEff4Pt10:Eff3Pt02:errEff3Pt02:Eff3Pt1:errEff3Pt1:Eff3Pt10:errEff3Pt10:Eff2Pt02:errEff2Pt02:Eff2Pt1:errEff2Pt1:Eff2Pt10:errEff2Pt10:EffSPDPt02:errEffSPDPt02:EffSPDPt1:errEffSPDPt1:EffSPDPt10:errEffSPDPt10:EffoneSPDPt02:errEffoneSPDPt02:EffoneSPDPt1:errEffoneSPDPt1:EffoneSPDPt10:errEffoneSPDPt10:EffTOTPt02:errEffTOTPt02:EffTOTPt1:errEffTOTPt1:EffTOTPt10:errEffTOTPt10");

  Float_t xntMatchingTOF[nVariablesMatchingTOF];

  //--------------------------------

  //----QA Vertex

  const Int_t nVariablesVertex=29;
  TNtuple* ntvertex=new TNtuple("ntvertex","QA Vertex","nrun:VxTRK:errVxTRK:sigmaVxTRK:errsigmaVxTRK:VyTRK:errVyTRK:sigmaVyTRK:errsigmaVyTRK:VzTRK:errVzTRK:sigmaVzTRK:errsigmaVzTRK:VxSPD:errVxSPD:sigmaVxSPD:errsigmaVxSPD:VySPD:errVySPD:sigmaVySPD:errsigmaVySPD:VzSPD:errVzSPD:sigmaVzSPD:errsigmaVzSPD:pileupSPD:errpileupSPD");
      
  Float_t xntVertex[nVariablesVertex];
  
  //--------------------------------

  //-----  Tracking ITS SA
/*  const Int_t nVariablesSA=13+18+2+504; // 144 eta-phi cluster occupancy per layer --> 504 for 40 phi bins (2+2)*6 eta + (40+40)*6 phi
  Float_t xntSA[nVariablesSA];
  TNtuple *ntSA=new TNtuple("ntITSsa","ntITSsa","run:NITSTPCPtBin0:NITSTPCPtBin1:NITSTPCPtBin2:NITSsaPtBin0:NITSsaPtBin1:NITSsaPtBin2:NITSpureSAPtBin0:NITSpureSAPtBin1:NITSpureSAPtBin2:ratioPtBin0:ratioPtBin1:ratioPtBin2:NcluITSpSA:errNcluITSpSA:dedx4_3:errdedx4_3:PtpionpSA:errPtpionpSA:NclupSA0:errNclupSA0:NclupSA1:errNclupSA1:NclupSA2:errNclupSA2:NclupSA3:errNclupSA3:NclupSA4:errNclupSA4:NclupSA5:errNclupSA5:chi2TPCITS:chi2ITSpureSA:occ_eta_1_1:erocc_eta_1_1:occ_eta_1_2:erocc_eta_1_2:occ_phi_1_1:erocc_phi_1_1:occ_phi_1_2:errocc_phi_1_2:occ_phi_1_3:erocc_phi_1_3:occ_phi_1_4:erocc_phi_1_4:occ_phi_1_5:erocc_phi_1_5:occ_phi_1_6:erocc_phi_1_6:occ_phi_1_7:erocc_phi_1_7:occ_phi_1_8:erocc_phi_1_8:occ_phi_1_9:erocc_phi_1_9:occ_phi_1_10:erocc_phi_1_10:occ_phi_1_11:erocc_phi_1_11:occ_phi_1_12:errocc_phi_1_12:occ_phi_1_13:erocc_phi_1_13:occ_phi_1_14:erocc_phi_1_14:occ_phi_1_15:erocc_phi_1_15:occ_phi_1_16:erocc_phi_1_16:occ_phi_1_17:erocc_phi_1_17:occ_phi_1_18:erocc_phi_1_18:occ_phi_1_19:erocc_phi_1_19:occ_phi_1_20:erocc_phi_1_20:occ_phi_1_21:erocc_phi_1_21:occ_phi_1_22:errocc_phi_1_22:occ_phi_1_23:erocc_phi_1_23:occ_phi_1_24:erocc_phi_1_24:occ_phi_1_25:erocc_phi_1_25:occ_phi_1_26:erocc_phi_1_26:occ_phi_1_27:erocc_phi_1_27:occ_phi_1_28:erocc_phi_1_28:occ_phi_1_29:erocc_phi_1_29:occ_phi_1_30:erocc_phi_1_30:occ_phi_1_31:erocc_phi_1_31:occ_phi_1_32:errocc_phi_1_32:occ_phi_1_33:erocc_phi_1_33:occ_phi_1_34:erocc_phi_1_34:occ_phi_1_35:erocc_phi_1_35:occ_phi_1_36:erocc_phi_1_36:occ_phi_1_37:erocc_phi_1_37:occ_phi_1_38:erocc_phi_1_38:occ_phi_1_39:erocc_phi_1_39:occ_phi_1_40:erocc_phi_1_40:occ_eta_2_1:erocc_eta_2_1:occ_eta_2_2:erocc_eta_2_2:occ_phi_2_1:erocc_phi_2_1:occ_phi_2_2:errocc_phi_2_2:occ_phi_2_3:erocc_phi_2_3:occ_phi_2_4:erocc_phi_2_4:occ_phi_2_5:erocc_phi_2_5:occ_phi_2_6:erocc_phi_2_6:occ_phi_2_7:erocc_phi_2_7:occ_phi_2_8:erocc_phi_2_8:occ_phi_2_9:erocc_phi_2_9:occ_phi_2_10:erocc_phi_2_10:occ_phi_2_11:erocc_phi_2_11:occ_phi_2_12:errocc_phi_2_12:occ_phi_2_13:erocc_phi_2_13:occ_phi_2_14:erocc_phi_2_14:occ_phi_2_15:erocc_phi_2_15:occ_phi_2_16:erocc_phi_2_16:occ_phi_2_17:erocc_phi_2_17:occ_phi_2_18:erocc_phi_2_18:occ_phi_2_19:erocc_phi_2_19:occ_phi_2_20:erocc_phi_2_20:occ_phi_2_21:erocc_phi_2_21:occ_phi_2_22:errocc_phi_2_22:occ_phi_2_23:erocc_phi_2_23:occ_phi_2_24:erocc_phi_2_24:occ_phi_2_25:erocc_phi_2_25:occ_phi_2_26:erocc_phi_2_26:occ_phi_2_27:erocc_phi_2_27:occ_phi_2_28:erocc_phi_2_28:occ_phi_2_29:erocc_phi_2_29:occ_phi_2_30:erocc_phi_2_30:occ_phi_2_31:erocc_phi_2_31:occ_phi_2_32:errocc_phi_2_32:occ_phi_2_33:erocc_phi_2_33:occ_phi_2_34:erocc_phi_2_34:occ_phi_2_35:erocc_phi_2_35:occ_phi_2_36:erocc_phi_2_36:occ_phi_2_37:erocc_phi_2_37:occ_phi_2_38:erocc_phi_2_38:occ_phi_2_39:erocc_phi_2_39:occ_phi_2_40:erocc_phi_2_40:occ_eta_3_1:erocc_eta_3_1:occ_eta_3_2:erocc_eta_3_2:occ_phi_3_1:erocc_phi_3_1:occ_phi_3_2:errocc_phi_3_2:occ_phi_3_3:erocc_phi_3_3:occ_phi_3_4:erocc_phi_3_4:occ_phi_3_5:erocc_phi_3_5:occ_phi_3_6:erocc_phi_3_6:occ_phi_3_7:erocc_phi_3_7:occ_phi_3_8:erocc_phi_3_8:occ_phi_3_9:erocc_phi_3_9:occ_phi_3_10:erocc_phi_3_10:occ_phi_3_11:erocc_phi_3_11:occ_phi_3_12:errocc_phi_3_12:occ_phi_3_13:erocc_phi_3_13:occ_phi_3_14:erocc_phi_3_14:occ_phi_3_15:erocc_phi_3_15:occ_phi_3_16:erocc_phi_3_16:occ_phi_3_17:erocc_phi_3_17:occ_phi_3_18:erocc_phi_3_18:occ_phi_3_19:erocc_phi_3_19:occ_phi_3_20:erocc_phi_3_20:occ_phi_3_21:erocc_phi_3_21:occ_phi_3_22:errocc_phi_3_22:occ_phi_3_23:erocc_phi_3_23:occ_phi_3_24:erocc_phi_3_24:occ_phi_3_25:erocc_phi_3_25:occ_phi_3_26:erocc_phi_3_26:occ_phi_3_27:erocc_phi_3_27:occ_phi_3_28:erocc_phi_3_28:occ_phi_3_29:erocc_phi_3_29:occ_phi_3_30:erocc_phi_3_30:occ_phi_3_31:erocc_phi_3_31:occ_phi_3_32:errocc_phi_3_32:occ_phi_3_33:erocc_phi_3_33:occ_phi_3_34:erocc_phi_3_34:occ_phi_3_35:erocc_phi_3_35:occ_phi_3_36:erocc_phi_3_36:occ_phi_3_37:erocc_phi_3_37:occ_phi_3_38:erocc_phi_3_38:occ_phi_3_39:erocc_phi_3_39:occ_phi_3_40:erocc_phi_3_40:occ_eta_4_1:erocc_eta_4_1:occ_eta_4_2:erocc_eta_4_2:occ_phi_4_1:erocc_phi_4_1:occ_phi_4_2:errocc_phi_4_2:occ_phi_4_3:erocc_phi_4_3:occ_phi_4_4:erocc_phi_4_4:occ_phi_4_5:erocc_phi_4_5:occ_phi_4_6:erocc_phi_4_6:occ_phi_4_7:erocc_phi_4_7:occ_phi_4_8:erocc_phi_4_8:occ_phi_4_9:erocc_phi_4_9:occ_phi_4_10:erocc_phi_4_10:occ_phi_4_11:erocc_phi_4_11:occ_phi_4_12:errocc_phi_4_12:occ_phi_4_13:erocc_phi_4_13:occ_phi_4_14:erocc_phi_4_14:occ_phi_4_15:erocc_phi_4_15:occ_phi_4_16:erocc_phi_4_16:occ_phi_4_17:erocc_phi_4_17:occ_phi_4_18:erocc_phi_4_18:occ_phi_4_19:erocc_phi_4_19:occ_phi_4_20:erocc_phi_4_20:occ_phi_4_21:erocc_phi_4_21:occ_phi_4_22:errocc_phi_4_22:occ_phi_4_23:erocc_phi_4_23:occ_phi_4_24:erocc_phi_4_24:occ_phi_4_25:erocc_phi_4_25:occ_phi_4_26:erocc_phi_4_26:occ_phi_4_27:erocc_phi_4_27:occ_phi_4_28:erocc_phi_4_28:occ_phi_4_29:erocc_phi_4_29:occ_phi_4_30:erocc_phi_4_30:occ_phi_4_31:erocc_phi_4_31:occ_phi_4_32:errocc_phi_4_32:occ_phi_4_33:erocc_phi_4_33:occ_phi_4_34:erocc_phi_4_34:occ_phi_4_35:erocc_phi_4_35:occ_phi_4_36:erocc_phi_4_36:occ_phi_4_37:erocc_phi_4_37:occ_phi_4_38:erocc_phi_4_38:occ_phi_4_39:erocc_phi_4_39:occ_phi_4_40:erocc_phi_4_40:occ_eta_5_1:erocc_eta_5_1:occ_eta_5_2:erocc_eta_5_2:occ_phi_5_1:erocc_phi_5_1:occ_phi_5_2:errocc_phi_5_2:occ_phi_5_3:erocc_phi_5_3:occ_phi_5_4:erocc_phi_5_4:occ_phi_5_5:erocc_phi_5_5:occ_phi_5_6:erocc_phi_5_6:occ_phi_5_7:erocc_phi_5_7:occ_phi_5_8:erocc_phi_5_8:occ_phi_5_9:erocc_phi_5_9:occ_phi_5_10:erocc_phi_5_10:occ_phi_5_11:erocc_phi_5_11:occ_phi_5_12:errocc_phi_5_12:occ_phi_5_13:erocc_phi_5_13:occ_phi_5_14:erocc_phi_5_14:occ_phi_5_15:erocc_phi_5_15:occ_phi_5_16:erocc_phi_5_16:occ_phi_5_17:erocc_phi_5_17:occ_phi_5_18:erocc_phi_5_18:occ_phi_5_19:erocc_phi_5_19:occ_phi_5_20:erocc_phi_5_20:occ_phi_5_21:erocc_phi_5_21:occ_phi_5_22:errocc_phi_5_22:occ_phi_5_23:erocc_phi_5_23:occ_phi_5_24:erocc_phi_5_24:occ_phi_5_25:erocc_phi_5_25:occ_phi_5_26:erocc_phi_5_26:occ_phi_5_27:erocc_phi_5_27:occ_phi_5_28:erocc_phi_5_28:occ_phi_5_29:erocc_phi_5_29:occ_phi_5_30:erocc_phi_5_30:occ_phi_5_31:erocc_phi_5_31:occ_phi_5_32:errocc_phi_5_32:occ_phi_5_33:erocc_phi_5_33:occ_phi_5_34:erocc_phi_5_34:occ_phi_5_35:erocc_phi_5_35:occ_phi_5_36:erocc_phi_5_36:occ_phi_5_37:erocc_phi_5_37:occ_phi_5_38:erocc_phi_5_38:occ_phi_5_39:erocc_phi_5_39:occ_phi_5_40:erocc_phi_5_40:occ_eta_6_1:erocc_eta_6_1:occ_eta_6_2:erocc_eta_6_2:occ_phi_6_1:erocc_phi_6_1:occ_phi_6_2:errocc_phi_6_2:occ_phi_6_3:erocc_phi_6_3:occ_phi_6_4:erocc_phi_6_4:occ_phi_6_5:erocc_phi_6_5:occ_phi_6_6:erocc_phi_6_6:occ_phi_6_7:erocc_phi_6_7:occ_phi_6_8:erocc_phi_6_8:occ_phi_6_9:erocc_phi_6_9:occ_phi_6_10:erocc_phi_6_10:occ_phi_6_11:erocc_phi_6_11:occ_phi_6_12:errocc_phi_6_12:occ_phi_6_13:erocc_phi_6_13:occ_phi_6_14:erocc_phi_6_14:occ_phi_6_15:erocc_phi_6_15:occ_phi_6_16:erocc_phi_6_16:occ_phi_6_17:erocc_phi_6_17:occ_phi_6_18:erocc_phi_6_18:occ_phi_6_19:erocc_phi_6_19:occ_phi_6_20:erocc_phi_6_20:occ_phi_6_21:erocc_phi_6_21:occ_phi_6_22:errocc_phi_6_22:occ_phi_6_23:erocc_phi_6_23:occ_phi_6_24:erocc_phi_6_24:occ_phi_6_25:erocc_phi_6_25:occ_phi_6_26:erocc_phi_6_26:occ_phi_6_27:erocc_phi_6_27:occ_phi_6_28:erocc_phi_6_28:occ_phi_6_29:erocc_phi_6_29:occ_phi_6_30:erocc_phi_6_30:occ_phi_6_31:erocc_phi_6_31:occ_phi_6_32:errocc_phi_6_32:occ_phi_6_33:erocc_phi_6_33:occ_phi_6_34:erocc_phi_6_34:occ_phi_6_35:erocc_phi_6_35:occ_phi_6_36:erocc_phi_6_36:occ_phi_6_37:erocc_phi_6_37:occ_phi_6_38:erocc_phi_6_38:occ_phi_6_39:erocc_phi_6_39:occ_phi_6_40:erocc_phi_6_40");
*/
    const Int_t nVariablesSA=13+18+2+252; // 144 eta-phi cluster occupancy per layer --> 252 for 40 phi bins 2*6 eta + 40*6 phi (no errori)
    Float_t xntSA[nVariablesSA];
    TNtuple *ntSA=new TNtuple("ntITSsa","ntITSsa","run:NITSTPCPtBin0:NITSTPCPtBin1:NITSTPCPtBin2:NITSsaPtBin0:NITSsaPtBin1:NITSsaPtBin2:NITSpureSAPtBin0:NITSpureSAPtBin1:NITSpureSAPtBin2:ratioPtBin0:ratioPtBin1:ratioPtBin2:NcluITSpSA:errNcluITSpSA:dedx4_3:errdedx4_3:PtpionpSA:errPtpionpSA:NclupSA0:errNclupSA0:NclupSA1:errNclupSA1:NclupSA2:errNclupSA2:NclupSA3:errNclupSA3:NclupSA4:errNclupSA4:NclupSA5:errNclupSA5:chi2TPCITS:chi2ITSpureSA:occ_eta_1_1:occ_eta_1_2:occ_phi_1_1:occ_phi_1_2:occ_phi_1_3:occ_phi_1_4:occ_phi_1_5:occ_phi_1_6:occ_phi_1_7:occ_phi_1_8:occ_phi_1_9:occ_phi_1_10:occ_phi_1_11:occ_phi_1_12:occ_phi_1_13:occ_phi_1_14:occ_phi_1_15:occ_phi_1_16:occ_phi_1_17:occ_phi_1_18:occ_phi_1_19:occ_phi_1_20:occ_phi_1_21:occ_phi_1_22:occ_phi_1_23:occ_phi_1_24:occ_phi_1_25:occ_phi_1_26:occ_phi_1_27:occ_phi_1_28:occ_phi_1_29:occ_phi_1_30:occ_phi_1_31:occ_phi_1_32:occ_phi_1_33:occ_phi_1_34:occ_phi_1_35:occ_phi_1_36:occ_phi_1_37:occ_phi_1_38:occ_phi_1_39:occ_phi_1_40:occ_eta_2_1:occ_eta_2_2:occ_phi_2_1:occ_phi_2_2:occ_phi_2_3:occ_phi_2_4:occ_phi_2_5:occ_phi_2_6:occ_phi_2_7:occ_phi_2_8:occ_phi_2_9:occ_phi_2_10:occ_phi_2_11:occ_phi_2_12:occ_phi_2_13:occ_phi_2_14:occ_phi_2_15:occ_phi_2_16:occ_phi_2_17:occ_phi_2_18:occ_phi_2_19:occ_phi_2_20:occ_phi_2_21:occ_phi_2_22:occ_phi_2_23:occ_phi_2_24:occ_phi_2_25:occ_phi_2_26:occ_phi_2_27:occ_phi_2_28:occ_phi_2_29:occ_phi_2_30:occ_phi_2_31:occ_phi_2_32:occ_phi_2_33:occ_phi_2_34:occ_phi_2_35:occ_phi_2_36:occ_phi_2_37:occ_phi_2_38:occ_phi_2_39:occ_phi_2_40:occ_eta_3_1:occ_eta_3_2:occ_phi_3_1:occ_phi_3_2:occ_phi_3_3:occ_phi_3_4:occ_phi_3_5:occ_phi_3_6:occ_phi_3_7:occ_phi_3_8:occ_phi_3_9:occ_phi_3_10:occ_phi_3_11:occ_phi_3_12:occ_phi_3_13:occ_phi_3_14:occ_phi_3_15:occ_phi_3_16:occ_phi_3_17:occ_phi_3_18:occ_phi_3_19:occ_phi_3_20:occ_phi_3_21:occ_phi_3_22:occ_phi_3_23:occ_phi_3_24:occ_phi_3_25:occ_phi_3_26:occ_phi_3_27:occ_phi_3_28:occ_phi_3_29:occ_phi_3_30:occ_phi_3_31:occ_phi_3_32:occ_phi_3_33:occ_phi_3_34:occ_phi_3_35:occ_phi_3_36:occ_phi_3_37:occ_phi_3_38:occ_phi_3_39:occ_phi_3_40:occ_eta_4_1:occ_eta_4_2:occ_phi_4_1:occ_phi_4_2:occ_phi_4_3:occ_phi_4_4:occ_phi_4_5:occ_phi_4_6:occ_phi_4_7:occ_phi_4_8:occ_phi_4_9:occ_phi_4_10:occ_phi_4_11:occ_phi_4_12:occ_phi_4_13:occ_phi_4_14:occ_phi_4_15:occ_phi_4_16:occ_phi_4_17:occ_phi_4_18:occ_phi_4_19:occ_phi_4_20:occ_phi_4_21:occ_phi_4_22:occ_phi_4_23:occ_phi_4_24:occ_phi_4_25:occ_phi_4_26:occ_phi_4_27:occ_phi_4_28:occ_phi_4_29:occ_phi_4_30:occ_phi_4_31:occ_phi_4_32:occ_phi_4_33:occ_phi_4_34:occ_phi_4_35:occ_phi_4_36:occ_phi_4_37:occ_phi_4_38:occ_phi_4_39:occ_phi_4_40:occ_eta_5_1:occ_eta_5_2:occ_phi_5_1:occ_phi_5_2:occ_phi_5_3:occ_phi_5_4:occ_phi_5_5:occ_phi_5_6:occ_phi_5_7:occ_phi_5_8:occ_phi_5_9:occ_phi_5_10:occ_phi_5_11:occ_phi_5_12:occ_phi_5_13:occ_phi_5_14:occ_phi_5_15:occ_phi_5_16:occ_phi_5_17:occ_phi_5_18:occ_phi_5_19:occ_phi_5_20:occ_phi_5_21:occ_phi_5_22:occ_phi_5_23:occ_phi_5_24:occ_phi_5_25:occ_phi_5_26:occ_phi_5_27:occ_phi_5_28:occ_phi_5_29:occ_phi_5_30:occ_phi_5_31:occ_phi_5_32:occ_phi_5_33:occ_phi_5_34:occ_phi_5_35:occ_phi_5_36:occ_phi_5_37:occ_phi_5_38:occ_phi_5_39:occ_phi_5_40:occ_eta_6_1:occ_eta_6_2:occ_phi_6_1:occ_phi_6_2:occ_phi_6_3:occ_phi_6_4:occ_phi_6_5:occ_phi_6_6:occ_phi_6_7:occ_phi_6_8:occ_phi_6_9:occ_phi_6_10:occ_phi_6_11:occ_phi_6_12:occ_phi_6_13:occ_phi_6_14:occ_phi_6_15:occ_phi_6_16:occ_phi_6_17:occ_phi_6_18:occ_phi_6_19:occ_phi_6_20:occ_phi_6_21:occ_phi_6_22:occ_phi_6_23:occ_phi_6_24:occ_phi_6_25:occ_phi_6_26:occ_phi_6_27:occ_phi_6_28:occ_phi_6_29:occ_phi_6_30:occ_phi_6_31:occ_phi_6_32:occ_phi_6_33:occ_phi_6_34:occ_phi_6_35:occ_phi_6_36:occ_phi_6_37:occ_phi_6_38:occ_phi_6_39:occ_phi_6_40",64000);
    
    //---------------------------------

//-----  PileUp from SPD vertices
    const Int_t nVariablesPU=11;
    Float_t xntPU[nVariablesPU];
    TNtuple *ntPU=new TNtuple("ntPileUp","ntPileUp","run:npilvtx:errnpilvtx:ntrklpil:errntrklpil:ntrklnopil:errntrklnopil:ncl1pil:errncl1pil:ncl1nopil:errncl1nopil");
    
    //---------------------------------
    //-----  ITS PID for TPCITS tracks
    const Int_t nVariablesPID=9;
    Float_t xntPID[nVariablesPID];
    TNtuple *ntPID=new TNtuple("ntPIDITS","ntPIDITS","run:nsigmapi02:errnsigmapi02:nsigmapi05:errnsigmapi05:nsigmapi1:errnsigmapi1:nsigmapi3:errnsigmapi3");
    
    //---------------------------------

    //---------------------------------
    //-----  DCA for TPCITS tracks
    const Int_t nVariablesDCA=2*4*4+1; // 4 pt bins, mean, ermean, RMS, erRMS, meanz ermeanz, RMSz, erRMSz
    Float_t xntDCA[nVariablesDCA];
    TNtuple *ntDCA=new TNtuple("ntDCAtracks","ntDCAtracks","run:dcamean_05:errdcamean_05:dcaRMS_05:errdcaRMS_05:dcazmean_05:errdcazmean_05:dcazRMS_05:errdcazRMS_05:dcamean_1:errdcamean_1:dcaRMS_1:errdcaRMS_1:dcazmean_1:errdcazmean_1:dcazRMS_1:errdcazRMS_1:dcamean_5:errdcamean_5:dcaRMS_5:errdcaRMS_5:dcazmean_5:errdcazmean_5:dcazRMS_5:errdcazRMS_5:dcamean_10:errdcamean_10:dcaRMS_10:errdcaRMS_10:dcazmean_10:errdcazmean_10:dcazRMS_10:errdcazRMS_10");
    
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

      TNtuple* ntmpPID=(TNtuple*)oldfil->Get("ntPIDITS");

      TNtuple* ntmpDCA=(TNtuple*)oldfil->Get("ntDCAtracks");

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
      
      //-----------------------
      
      //-----------PID ITS ----------------
      Bool_t isOKPID = kFALSE;
      if(ntmpPID){
          if(ntmpPID->GetNvar()==ntPID->GetNvar()){
              isOKPID = kTRUE;
              TObjArray* arr1PID=(TObjArray*)ntPID->GetListOfBranches();
              TObjArray* arr2PID=(TObjArray*)ntmpPID->GetListOfBranches();
              for(Int_t iV=0; iV<ntmpPID->GetNvar(); iV++){
                  TString vnam1=arr1PID->At(iV)->GetName();
                  TString vnam2=arr2PID->At(iV)->GetName();
                  if(vnam1!=vnam2) isOKPID=kFALSE;
                  ntmpPID->SetBranchAddress(vnam2.Data(),&xntPID[iV]);
              }
              if(isOKPID){
                  for(Int_t nE=0; nE<ntmpPID->GetEntries(); nE++){
                      ntmpPID->GetEvent(nE);
                      Int_t theRun=(Int_t)(xntPID[0]+0.0001);
                      readRun->SetBitNumber(theRun);
                      ntPID->Fill(xntPID);
                  }
              }
          }
      }
      
      //-----------------------
      
      //----------- Tracks DCA ----------------
      Bool_t isOKDCA = kFALSE;
      if(ntmpDCA){
          if(ntmpDCA->GetNvar()==ntDCA->GetNvar()){
              isOKDCA = kTRUE;
              TObjArray* arr1DCA=(TObjArray*)ntDCA->GetListOfBranches();
              TObjArray* arr2DCA=(TObjArray*)ntmpDCA->GetListOfBranches();
              for(Int_t iV=0; iV<ntmpDCA->GetNvar(); iV++){
                  TString vnam1=arr1DCA->At(iV)->GetName();
                  TString vnam2=arr2DCA->At(iV)->GetName();
                  if(vnam1!=vnam2) isOKDCA=kFALSE;
                  ntmpDCA->SetBranchAddress(vnam2.Data(),&xntDCA[iV]);
              }
              if(isOKDCA){
                  for(Int_t nE=0; nE<ntmpDCA->GetEntries(); nE++){
                      ntmpDCA->GetEvent(nE);
                      Int_t theRun=(Int_t)(xntDCA[0]+0.0001);
                      readRun->SetBitNumber(theRun);
                      ntDCA->Fill(xntDCA);
                  }
              }
          }
      }

      //-----------------------

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
    else if(aux.Contains("LHC13b/")){
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
      
    //------- check ITS presence in QAresults.root file: pp & PbPb runs!

      TDirectoryFile* df=(TDirectoryFile*)f->Get("Vertex_Performance");
      if(!df) continue;
      
      TList* l=(TList*)df->Get("cOutputVtxESD");
      if(!l) continue;
      
      TH1F* h=(TH1F*)l->FindObject("fhSPDVertexZonly"); // number of processed events
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
	
    //--------------   ITS SA ---------------
    cout<<"ITS - SA"<<endl;
    myfile << "ITS - SA " << endl;

    //    cout<<f<<" "<<ntSA<<" "<<iRun<<endl;
    FillITSSAntuple(f,ntSA,iRun,xntSA);

      //--------------   PileUp ---------------
      //
      myfile << "PileUp " << endl;

      FillPileUpntuple(f,ntPU,iRun,xntPU);

      //--------------   PID ITS ---------------
      //
      cout<<"ITS - PID"<<endl;
      myfile << "PID " << endl;
      
      FillPIDntuple(f,ntPID,iRun,xntPID);
      
      //--------------   Tracks DCA ---------------
      //
      cout<<"DCA - TPCITS tracks"<<endl;
      myfile << "DCA " << endl;
      
      FillDCAntuple(f,ntDCA,iRun,xntDCA);
      

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
    ntPID->Write();
    ntDCA->Write();
  outfil->Close();
  delete outfil;
  delete ntsdd;
  delete ntssd;
  delete ntmatching;
  delete ntvertex;
  delete ntSA;
  delete ntPU;
    delete ntPID;
    delete ntDCA;
    
}

//____________________________________________________________________________
void FillITSSAntuple(TFile* f,TNtuple* nt, Int_t nrun, Float_t *xntSA){
//  static const Int_t nVariables=13+18+2+252;
//  static Float_t xnt[nVariables];
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
    if(!dirFile){
        printf("Run %d ITSsaTracks MISSING -> Exit\n",nrun);
        myfile<<" ITSsaTracks MISSING -> Exit\n";
        return;
    }
  TList *cOutput = (TList*)dirFile->Get("clistITSsaTracks");
    if(!cOutput){
        printf("Run %d clistITSsaTracks MISSING -> Exit\n",nrun);
        myfile<<" clistITSsaTracks MISSING -> Exit\n";
        return;
    }

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
    TH1F* hChi2TPCITS=(TH1F*)cOutput->FindObject("hChi2TPCITS");
    TH1F* hChi2ITSpureSA=(TH1F*)cOutput->FindObject("hChi2ITSpureSA");


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
    TH1F* hPtpSA =(TH1F*)cOutput->FindObject("hPtITSpureSA");  //
    TH1F* hPtpionpSA =(TH1F*)cOutput->FindObject("hPtITSpureSAPion");  //
    TH2F* hNcluPion2d = (TH2F*)cOutput->FindObject("hCluInLayITSpureSAPion");
    TH1D* hNcluPion;
    if(hNcluPion2d) hNcluPion = hNcluPion2d->ProjectionY();
    TH2F* hNclu2d = (TH2F*)cOutput->FindObject("hCluInLayVsPtITSpureSA"); // <-- check existence
    TH1D* hNclu;
    Bool_t is2d = kFALSE;
    if(hNclu2d) {hNclu = hNclu2d->ProjectionY(); is2d = kTRUE;}
    
    Float_t e4, e3;
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

    if(is2d){
        if(hNclu->GetBinContent(1)>0){
        for(Int_t iLay=0; iLay<6; iLay++){
//            fracTpi[iLay]=hNcluPion->GetBinContent(iLay+2)/hNcluPion->GetBinContent(1);
//            efracTpi[iLay]=TMath::Sqrt(fracTpi[iLay]*(1-fracTpi[iLay])/hNcluPion->GetBinContent(1));
            fracTpi[iLay]=hNclu->GetBinContent(iLay+2)/hNclu->GetBinContent(1);
            efracTpi[iLay]=TMath::Sqrt(fracTpi[iLay]*(1-fracTpi[iLay])/hNclu->GetBinContent(1));
        }
        }
    }
    else
    {
        for(Int_t iLay=0; iLay<6; iLay++){
            if(hNcluPion->GetBinContent(1)>0){
                fracTpi[iLay]=hNcluPion->GetBinContent(iLay+2)/hNcluPion->GetBinContent(1);
                efracTpi[iLay]=TMath::Sqrt(fracTpi[iLay]*(1-fracTpi[iLay])/hNcluPion->GetBinContent(1));
        }
        }
    }

// Elena
    
    // eta-phi cluster plots in single layers
    TH2F* hEtaPhiTracksLay1TPCITS=(TH2F*)cOutput->FindObject("hEtaPhiTracksLay1TPCITS");
    TH2F* hEtaPhiTracksLay2TPCITS=(TH2F*)cOutput->FindObject("hEtaPhiTracksLay2TPCITS");
    TH2F* hEtaPhiTracksLay3TPCITS=(TH2F*)cOutput->FindObject("hEtaPhiTracksLay3TPCITS");
    TH2F* hEtaPhiTracksLay4TPCITS=(TH2F*)cOutput->FindObject("hEtaPhiTracksLay4TPCITS");
    TH2F* hEtaPhiTracksLay5TPCITS=(TH2F*)cOutput->FindObject("hEtaPhiTracksLay5TPCITS");
    TH2F* hEtaPhiTracksLay6TPCITS=(TH2F*)cOutput->FindObject("hEtaPhiTracksLay6TPCITS");
    TH1D* clus_phi_1, *clus_phi_2, *clus_phi_3, *clus_phi_4, *clus_phi_5, *clus_phi_6;
    TH1D* clus_eta_1, *clus_eta_2, *clus_eta_3, *clus_eta_4, *clus_eta_5, *clus_eta_6;
//    Float_t phi[10] = {1.,2.,3.,4.,5.,6.,7.,8.,9.,10.}; // serve se si scrivono i valori effettivi
//    Float_t eta[2]={1.,2.}; //  ma serve al momento del plot
    Float_t occ_phi_1[40] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
//    Float_t erocc_phi_1[40] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    Float_t occ_eta_1[2] = {0.,0.};
//    Float_t erocc_eta_1[2] = {0.,0.};
    Float_t occ_phi_2[40] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
//    Float_t erocc_phi_2[40] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    Float_t occ_eta_2[2] = {0.,0.};
//    Float_t erocc_eta_2[2] = {0.,0.};
    Float_t occ_phi_3[40] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
//    Float_t erocc_phi_3[40] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    Float_t occ_eta_3[2] = {0.,0.};
//    Float_t erocc_eta_3[2] = {0.,0.};
    Float_t occ_phi_4[40] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
//    Float_t erocc_phi_4[40] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    Float_t occ_eta_4[2] = {0.,0.};
//    Float_t erocc_eta_4[2] = {0.,0.};
    Float_t occ_phi_5[40] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
//    Float_t erocc_phi_5[40] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    Float_t occ_eta_5[2] = {0.,0.};
//    Float_t erocc_eta_5[2] = {0.,0.};
    Float_t occ_phi_6[40] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
//    Float_t erocc_phi_6[40] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    Float_t occ_eta_6[2] = {0.,0.};
//    Float_t erocc_eta_6[2] = {0.,0.};
    
    
    if(hEtaPhiTracksLay1TPCITS){
        clus_phi_1 = hEtaPhiTracksLay1TPCITS->ProjectionY();
        clus_phi_1->Rebin(5); // 40 bins in phi, 9 degrees/bin
        for(Int_t phibin=0;phibin<40;phibin++){
            occ_phi_1[phibin]=clus_phi_1->GetBinContent(phibin+1)/clus_phi_1->GetEntries();
//            erocc_phi_1[phibin]= TMath::Sqrt(occ_phi_1[phibin]*(1-occ_phi_1[phibin])/clus_phi_1->GetEntries());
        }
        clus_eta_1 = hEtaPhiTracksLay1TPCITS->ProjectionX();
        clus_eta_1->Rebin(25); // 2 bins in phi, HS, HL
         for(Int_t etabin=0;etabin<2;etabin++){
            occ_eta_1[etabin]=clus_eta_1->GetBinContent(etabin+1)/clus_eta_1->GetEntries();
//            erocc_eta_1[etabin]= TMath::Sqrt(occ_eta_1[etabin]*(1-occ_eta_1[etabin])/clus_eta_1->GetEntries());
//            cout << "etabin = " << etabin << " bin center = " << clus_eta_1->GetBinCenter(etabin+1) << " occ_eta_1[etabin] = " << occ_eta_1[etabin] << endl;
        }
    }
    
    
    if(hEtaPhiTracksLay2TPCITS){
        clus_phi_2 = hEtaPhiTracksLay2TPCITS->ProjectionY();
        clus_phi_2->Rebin(5); // 40 bins in phi, 36 degrees/bin
        for(Int_t phibin=0;phibin<40;phibin++){
            occ_phi_2[phibin]=clus_phi_2->GetBinContent(phibin+1)/clus_phi_2->GetEntries();
//            erocc_phi_2[phibin]= TMath::Sqrt(occ_phi_2[phibin]*(1-occ_phi_2[phibin])/clus_phi_2->GetEntries());
        }
        clus_eta_2 = hEtaPhiTracksLay2TPCITS->ProjectionX();
        clus_eta_2->Rebin(25); // 2 bins in phi, HS, HL
        for(Int_t etabin=0;etabin<2;etabin++){
            occ_eta_2[etabin]=clus_eta_2->GetBinContent(etabin+1)/clus_eta_2->GetEntries();
//            erocc_eta_2[etabin]= TMath::Sqrt(occ_eta_2[etabin]*(1-occ_eta_2[etabin])/clus_eta_2->GetEntries());
        }
    }

    if(hEtaPhiTracksLay3TPCITS){
        clus_phi_3 = hEtaPhiTracksLay3TPCITS->ProjectionY();
        clus_phi_3->Rebin(5); // 40 bins in phi, 36 degrees/bin
        for(Int_t phibin=0;phibin<40;phibin++){
            occ_phi_3[phibin]=clus_phi_3->GetBinContent(phibin+1)/clus_phi_3->GetEntries();
//            erocc_phi_3[phibin]= TMath::Sqrt(occ_phi_3[phibin]*(1-occ_phi_3[phibin])/clus_phi_3->GetEntries());
        }
        clus_eta_3 = hEtaPhiTracksLay3TPCITS->ProjectionX();
        clus_eta_3->Rebin(25); // 2 bins in phi, HS, HL
        for(Int_t etabin=0;etabin<2;etabin++){
            occ_eta_3[etabin]=clus_eta_3->GetBinContent(etabin+1)/clus_eta_3->GetEntries();
//            erocc_eta_3[etabin]= TMath::Sqrt(occ_eta_3[etabin]*(1-occ_eta_3[etabin])/clus_eta_3->GetEntries());
        }
    }

    if(hEtaPhiTracksLay4TPCITS){
        clus_phi_4 = hEtaPhiTracksLay4TPCITS->ProjectionY();
        clus_phi_4->Rebin(5); // 40 bins in phi, 36 degrees/bin
        for(Int_t phibin=0;phibin<40;phibin++){
            occ_phi_4[phibin]=clus_phi_4->GetBinContent(phibin+1)/clus_phi_4->GetEntries();
//            erocc_phi_4[phibin]= TMath::Sqrt(occ_phi_4[phibin]*(1-occ_phi_4[phibin])/clus_phi_4->GetEntries());
        }
        clus_eta_4 = hEtaPhiTracksLay4TPCITS->ProjectionX();
        clus_eta_4->Rebin(25); // 2 bins in phi, HS, HL
        for(Int_t etabin=0;etabin<2;etabin++){
            occ_eta_4[etabin]=clus_eta_4->GetBinContent(etabin+1)/clus_eta_4->GetEntries();
//            erocc_eta_4[etabin]= TMath::Sqrt(occ_eta_4[etabin]*(1-occ_eta_4[etabin])/clus_eta_4->GetEntries());
        }
    }
    
    if(hEtaPhiTracksLay5TPCITS){
        clus_phi_5 = hEtaPhiTracksLay5TPCITS->ProjectionY();
        clus_phi_5->Rebin(5); // 40 bins in phi, 36 degrees/bin
        for(Int_t phibin=0;phibin<40;phibin++){
            occ_phi_5[phibin]=clus_phi_5->GetBinContent(phibin+1)/clus_phi_5->GetEntries();
//            erocc_phi_5[phibin]= TMath::Sqrt(occ_phi_5[phibin]*(1-occ_phi_5[phibin])/clus_phi_5->GetEntries());
        }
        clus_eta_5 = hEtaPhiTracksLay5TPCITS->ProjectionX();
        clus_eta_5->Rebin(25); // 2 bins in phi, HS, HL
        for(Int_t etabin=0;etabin<2;etabin++){
            occ_eta_5[etabin]=clus_eta_5->GetBinContent(etabin+1)/clus_eta_5->GetEntries();
//            erocc_eta_5[etabin]= TMath::Sqrt(occ_eta_5[etabin]*(1-occ_eta_5[etabin])/clus_eta_5->GetEntries());
        }
    }
    
    if(hEtaPhiTracksLay6TPCITS){
        clus_phi_6 = hEtaPhiTracksLay6TPCITS->ProjectionY();
        clus_phi_6->Rebin(5); // 40 bins in phi, 36 degrees/bin
        for(Int_t phibin=0;phibin<40;phibin++){
            occ_phi_6[phibin]=clus_phi_6->GetBinContent(phibin+1)/clus_phi_6->GetEntries();
//            erocc_phi_6[phibin]= TMath::Sqrt(occ_phi_6[phibin]*(1-occ_phi_6[phibin])/clus_phi_6->GetEntries());
        }
        clus_eta_6 = hEtaPhiTracksLay6TPCITS->ProjectionX();
        clus_eta_6->Rebin(25); // 2 bins in phi, HS, HL
        for(Int_t etabin=0;etabin<2;etabin++){
            occ_eta_6[etabin]=clus_eta_6->GetBinContent(etabin+1)/clus_eta_6->GetEntries();
//            erocc_eta_6[etabin]= TMath::Sqrt(occ_eta_6[etabin]*(1-occ_eta_6[etabin])/clus_eta_6->GetEntries());
        }
    }

    
    Int_t index=0;
  xntSA[index++]=(Float_t)nrun;
  xntSA[index++]=NTPCITS[0];
  xntSA[index++]=NTPCITS[1];
  xntSA[index++]=NTPCITS[2];
  xntSA[index++]=NITSsa[0];
  xntSA[index++]=NITSsa[1];
  xntSA[index++]=NITSsa[2];
  xntSA[index++]=NITSpureSA[0];
  xntSA[index++]=NITSpureSA[1];
  xntSA[index++]=NITSpureSA[2];
  xntSA[index++]=Ratio[0];
  xntSA[index++]=Ratio[1];
  xntSA[index++]=Ratio[2];
// 13

// Elena
    xntSA[index++]=hNcluITSpSA->GetMean();
    xntSA[index++]=hNcluITSpSA->GetMeanError();
    xntSA[index++]=r43;
    xntSA[index++]=erre43;
//    xnt[index++]=hPtpionpSA->GetMean();
//    xnt[index++]=hPtpionpSA->GetMeanError();
    if(hPtpSA) xntSA[index++]=hPtpSA->GetMean();
    else xntSA[index++]=0.;
    if(hPtpSA)xntSA[index++]=hPtpSA->GetMeanError();
    else xntSA[index++]=0.;
    xntSA[index++]=fracTpi[0];
    xntSA[index++]=efracTpi[0];
    xntSA[index++]=fracTpi[1];
    xntSA[index++]=efracTpi[1];
    xntSA[index++]=fracTpi[2];
    xntSA[index++]=efracTpi[2];
    xntSA[index++]=fracTpi[3];
    xntSA[index++]=efracTpi[3];
    xntSA[index++]=fracTpi[4];
    xntSA[index++]=efracTpi[4];
    xntSA[index++]=fracTpi[5];
    xntSA[index++]=efracTpi[5];
    // SA31

    // chi2
    xntSA[index++]=hChi2TPCITS->GetMean();
    xntSA[index++]=hChi2ITSpureSA->GetMean();
    // 33
    
    // eta-phi occupancy (504 variables)
    for(Int_t etabin=0;etabin<2;etabin++){
        xntSA[index++]= occ_eta_1[etabin];
//        xntSA[index++]= erocc_eta_1[etabin];
    }

    for(Int_t phibin=0;phibin<40;phibin++){
        xntSA[index++]= occ_phi_1[phibin];
//        xntSA[index++]= erocc_phi_1[phibin];
    }
    
    for(Int_t etabin=0;etabin<2;etabin++){
        xntSA[index++]= occ_eta_2[etabin];
//        xntSA[index++]= erocc_eta_2[etabin];
    }
    for(Int_t phibin=0;phibin<40;phibin++){
        xntSA[index++]= occ_phi_2[phibin];
//        xntSA[index++]= erocc_phi_2[phibin];
    }
    
    for(Int_t etabin=0;etabin<2;etabin++){
        xntSA[index++]= occ_eta_3[etabin];
//        xntSA[index++]= erocc_eta_3[etabin];
    }
    for(Int_t phibin=0;phibin<40;phibin++){
        xntSA[index++]= occ_phi_3[phibin];
//        xntSA[index++]= erocc_phi_3[phibin];
    }

    for(Int_t etabin=0;etabin<2;etabin++){
        xntSA[index++]= occ_eta_4[etabin];
//        xntSA[index++]= erocc_eta_4[etabin];
    }
    for(Int_t phibin=0;phibin<40;phibin++){
        xntSA[index++]= occ_phi_4[phibin];
//        xntSA[index++]= erocc_phi_4[phibin];
    }

    for(Int_t etabin=0;etabin<2;etabin++){
        xntSA[index++]= occ_eta_5[etabin];
//        xntSA[index++]= erocc_eta_5[etabin];
    }

    for(Int_t phibin=0;phibin<40;phibin++){
        xntSA[index++]= occ_phi_5[phibin];
//        xntSA[index++]= erocc_phi_5[phibin];
//        cout << " layer 5, erocc = " << erocc_phi_5[phibin] << " phibin = " << phibin << endl;
    }

    for(Int_t etabin=0;etabin<2;etabin++){
        xntSA[index++]= occ_eta_6[etabin];
//        xntSA[index++]= erocc_eta_6[etabin];
    }
    for(Int_t phibin=0;phibin<40;phibin++){
        xntSA[index++]= occ_phi_6[phibin];
//        xntSA[index++]= erocc_phi_6[phibin];
    }
    // 285

  nt->Fill(xntSA);
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
    if(!QAchargeRatio) QAchargeRatio=(TH2F*)lSSD->FindObject("QAChargeRatio");
      
  if(QAchargeRatio->GetEntries()==0){
    printf("Run %d QAchargeRatio EMPTY -> Continue\n",iRun);
    myfile<<"QA charge ratio empty"<<endl;

      // return;
  }

  TH2F* QAcharge=(TH2F*)lSSD->FindObject("QAChargeSA");
    if(!QAcharge) QAcharge=(TH2F*)lSSD->FindObject("QACharge");

    
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
    Float_t bad_n5=0, bad_p5=0, bad_n6=0, bad_p6=0;
    Float_t n_subjobs=0;
    if(bad_n && bad_p){
    //    Int_t max_n=0; // not used: same maximum for both histos
    for(Int_t i=1;i<1699;i++) if(bad_p->GetBinContent(i)>max_p)max_p=bad_p->GetBinContent(i);
    n_subjobs=max_p/768.;
// find the number of bad n/p strips per layer
    bad_p->Scale(1/n_subjobs);
    bad_n->Scale(1/n_subjobs);
    for(Int_t j=1; j<749; j++) {bad_n5=bad_n5+bad_n->GetBinContent(j); bad_p5=bad_p5+bad_p->GetBinContent(j);}
    for(Int_t j=749; j<1699; j++) {bad_n6=bad_n6+bad_n->GetBinContent(j); bad_p6=bad_p6+bad_p->GetBinContent(j);}
    }
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
    TDirectoryFile *dirMatchSPD=(TDirectoryFile*)f->GetDirectory("SPD_Performance");
    
    if(!dirMatch && !dirMatchSPD){
        printf("Run %d ITS_Performance and SPD_Performance MISSING -> Exit\n",iRun);
        myfile<<"ITS_Performance and SDD Performance missing"<<endl;
        
        return;
    }
    
    TList *list=NULL;
    TList *listSPD=NULL;
      
    if(dirMatch) list = (TList*)dirMatch->Get("cOutputITS"); // LHC12e
    if(!list){
        printf("Run %d coutputITS TList MISSING -> Exit\n",iRun);
        myfile<<"coutputITS empty"<<endl;
        
        return;
    }


    if(dirMatchSPD) listSPD = (TList*)dirMatchSPD->Get("coutput1");
    if(!listSPD){
        printf("Run %d coutput1 ITSTList MISSING -> Exit\n",iRun);
        myfile<<"coutput1 empty"<<endl;
        
        return;
    }
  
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
    TDirectoryFile *dirMatchSPD=(TDirectoryFile*)f->GetDirectory("SPD_Performance");

    if(!dirMatch && !dirMatchSPD){
        printf("Run %d ITS_Performence and SPD_Performance MISSING -> Exit\n",iRun);
        myfile<<"ITS_Performence and SDD Performance missing"<<endl;
        
        return;
    }
    
    TList *list=NULL;
    TList *listSPD=NULL;
    
    if(dirMatch) list = (TList*)dirMatch->Get("cOutputITS"); // LHC12e
    if(!list){
        printf("Run %d coutputITS TList MISSING -> Exit\n",iRun);
        myfile<<"coutputITS empty"<<endl;
        
        return;
    }
    
    
    if(dirMatchSPD) listSPD = (TList*)dirMatchSPD->Get("coutput1");
    if(!listSPD){
        printf("Run %d coutput1 ITSTList MISSING -> Exit\n",iRun);
        myfile<<"coutput1 empty"<<endl;
        
        return;
    }

    
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
        return;
    }
  
    TList *lt = (TList*)dirVertex->Get("cOutputVtxESD");
    if(!lt){
        Printf("cOutputVtxESD TList not found... check!");
        myfile<<"cOutputVtxESD TList not found"<<endl;
        return;
    }
	
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
      printf("Run %d xVtxSPD EMPTY -> Continue\n",iRun);
      myfile<<"xVtxSPD not found"<<endl;

        // return;
   }

    TH1F *yVtxSPD = (TH1F*)lt->FindObject("fhSPDVertexY");

    if(yVtxSPD->GetEntries()==0){
      printf("Run %d yVtxSPD EMPTY -> Continue\n",iRun);
      myfile<<"yVtxSPD not found"<<endl;

        // return;
    }

    TH1F *zVtxSPD = (TH1F*)lt->FindObject("fhSPDVertexZ"); // pp runs
    if(zVtxSPD->GetEntries()==0) zVtxSPD = (TH1F*)lt->FindObject("fhSPDVertexZonly"); // PbPb runs!!!

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
        return;
    }
    
    TList *lt = (TList*)dirPU->Get("clistPileupSPDQA");
    if(!lt){
        Printf("clistPileupSPDQA TList not found... check!");
        myfile<<"clistPileupSPDQA TList not found"<<endl;
        return;
    }
    
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

    if(hnPilVtx->GetBinContent(2)){
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

/////////////////////////////////////////////////////////////////////////////////////////////////////

void FillPIDntuple(TFile* f,TNtuple* nt, Int_t nrun, Float_t *xntPID){
    
//    static const Int_t nVariables=9;
//    static Float_t xnt[nVariables];
    
    Float_t  ptbin=0.0;
    Float_t ptvalue[4]={0.201, 0.501, 1.001,2.001};
    Int_t bin[4]={25,57,81,105};
    //=fHistPtITSMI6InAcc->FindBin(0.201);
    Double_t NsigmaTPCITS[4]={100.,100.,100.,100.};
    Double_t Nsigma_err_TPCITS[4]={0.,0.,0.,0.};
    
    TDirectoryFile *QAdata = (TDirectoryFile*)f->Get("PIDqa");
    if(!QAdata){
        Printf("PIDqa directory not found... return!");
        return;
    }
    
    TList *fListData = (TList*)QAdata->Get("PIDqa");
    if(!fListData){
        Printf("PIDqa/PIDqa directory not found... return!");
        return;
    }

    TList *fListITS = (TList*)fListData->FindObject("ITS");
    if(!fListITS){
        Printf("fListITS TList not found... return!");
        return;
    }
    
    
    TH2F *nsigma_pi = (TH2F *)fListITS->FindObject("hNsigmaP_ITS_pion");
    //    nsigma_pi->SetTitle(Form("%s %s",nsigma_pi->GetTitle(),fTitleData.Data()));
    
    /*    TH1D *nsigma_pi_1d = (TH1D *) nsigma_pi->ProfileX();
     for(Int_t l=0;l<4;l++){
     ptbin = nsigma_pi_1d->FindBin(ptvalue[l]);
     NsigmaTPCITS[l]=nsigma_pi_1d->GetBinContent(ptbin);
     Nsigma_err_TPCITS[l]=nsigma_pi_1d->GetBinError(ptbin);
     //        cout << "ptvalue = " << ptvalue[l] << " nsigma = " << NsigmaTPCITS[l] << endl;
     }
     */
    Int_t biny=-1;
    Float_t binycont =-100;
    
    if(nsigma_pi->GetEntries()>0){
        for(Int_t l=0;l<4;l++){
            TH1D *nsigma_pi_1d=(TH1D *) nsigma_pi->ProjectionY("pro",bin[l],bin[l]);
            biny = nsigma_pi_1d->GetMaximumBin();
            NsigmaTPCITS[l] = nsigma_pi_1d->GetBinLowEdge(biny)+0.5*nsigma_pi_1d->GetBinWidth(biny);
            Nsigma_err_TPCITS[l] = 0.5*nsigma_pi_1d->GetBinWidth(biny);
        }
    }
    Int_t index=0;
    xntPID[index++]=(Float_t)nrun;
    xntPID[index++]=NsigmaTPCITS[0];
    xntPID[index++]=Nsigma_err_TPCITS[0];
    xntPID[index++]=NsigmaTPCITS[1];
    xntPID[index++]=Nsigma_err_TPCITS[1];
    xntPID[index++]=NsigmaTPCITS[2];
    xntPID[index++]=Nsigma_err_TPCITS[2];
    xntPID[index++]=NsigmaTPCITS[3];
    xntPID[index++]=Nsigma_err_TPCITS[3];
    
    
    nt->Fill(xntPID);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

void FillDCAntuple(TFile* f,TNtuple* nt, Int_t nrun, Float_t *xntDCA){
    
//    static const Int_t nVariables=17+16; /// 16: DCArphi, 16 DCAz
//    static Float_t xnt[nVariables];
    
    TDirectoryFile *DCAdata = (TDirectoryFile*)f->Get("ImpParRes_Performance");
    if(!DCAdata){
        Printf("ImpParRes_Performance directory not found... return!");
        return;
    }
    
    TDirectoryFile *fListDCA = (TDirectoryFile*)DCAdata->Get("coutputd0allPointRec_0_1000000");
    if(!fListDCA){
        Printf("coutputd0allPointRec_0_1000000 directory not found... return!");
        return;
    }
    
    Float_t dcamean[4]={-1000.,-1000.,-1000.,-1000.};
    Float_t err_dcamean[4]={0.01,0.01,0.01,0.01};
    Float_t dcaRMS[4]={0.0,0.0,0.0,0.0};
    Float_t err_dcaRMS[4]={0.0,0.0,0.0,0.0};
    
    Float_t dcazmean[4]={-1000.,-1000.,-1000.,-1000.};
    Float_t err_dcazmean[4]={0.01,0.01,0.01,0.01};
    Float_t dcazRMS[4]={0.0,0.0,0.0,0.0};
    Float_t err_dcazRMS[4]={0.0,0.0,0.0,0.0};
    
    TF1 *h = new TF1("h","gaus",-10000.,10000.);
    Double_t xmin=0.0, xmax=0.0;
    
    TH1F *hDCA_05 = (TH1F *)fListDCA->FindObject("d0allpointrphiRec_5");
    if(hDCA_05->GetEntries()>0){
        xmin=hDCA_05->GetMean()-3*hDCA_05->GetRMS();
        xmax=hDCA_05->GetMean()+3*hDCA_05->GetRMS();
        hDCA_05->Fit("h","NR,Q","",xmin,xmax);
//     dcamean[0]=hDCA_05->GetMean();
//     err_dcamean[0]=hDCA_05->GetMeanError();
//     dcaRMS[0]=hDCA_05->GetRMS();
//     err_dcaRMS[0]=hDCA_05->GetRMSError();
       dcamean[0]=h->GetParameter(1);
       err_dcamean[0]=h->GetParError(1);
       dcaRMS[0]=h->GetParameter(2);
       err_dcaRMS[0]=h->GetParError(2);

    }
                 
    TH1F *hDCA_1 = (TH1F *)fListDCA->FindObject("d0allpointrphiRec_9");
    if(hDCA_1->GetEntries()>0){
        xmin=hDCA_1->GetMean()-3*hDCA_1->GetRMS();
        xmax=hDCA_1->GetMean()+3*hDCA_1->GetRMS();
        hDCA_1->Fit("h","NR,Q","",xmin,xmax);
//     dcamean[1]=hDCA_1->GetMean();
//     err_dcamean[1]=hDCA_1->GetMeanError();
//     dcaRMS[1]=hDCA_1->GetRMS();
//     err_dcaRMS[1]=hDCA_1->GetRMSError();
        dcamean[1]=h->GetParameter(1);
        err_dcamean[1]=h->GetParError(1);
        dcaRMS[1]=h->GetParameter(2);
        err_dcaRMS[1]=h->GetParError(2);
    }
                 
    TH1F *hDCA_5 = (TH1F *)fListDCA->FindObject("d0allpointrphiRec_18");
    if(hDCA_5->GetEntries()>0){
        xmin=hDCA_5->GetMean()-3*hDCA_5->GetRMS();
        xmax=hDCA_5->GetMean()+3*hDCA_5->GetRMS();
        hDCA_5->Fit("h","NR,Q","",xmin,xmax);
//     dcamean[2]=hDCA_5->GetMean();
//     err_dcamean[2]=hDCA_5->GetMeanError();
//     dcaRMS[2]=hDCA_5->GetRMS();
//     err_dcaRMS[2]=hDCA_5->GetRMSError();
        dcamean[2]=h->GetParameter(1);
        err_dcamean[2]=h->GetParError(1);
        dcaRMS[2]=h->GetParameter(2);
        err_dcaRMS[2]=h->GetParError(2);
    }
    
    TH1F *hDCA_10 = (TH1F *)fListDCA->FindObject("d0allpointrphiRec_20");
//    TH1F *hDCA_10b = (TH1F *)fListDCA->FindObject("d0allpointrphiRec_21");
//    TH1F *hDCA_10c = (TH1F *)fListDCA->FindObject("d0allpointrphiRec_22");
//    TH1F *hDCA_10d = (TH1F *)fListDCA->FindObject("d0allpointrphiRec_23");
//    hDCA_10->Add(hDCA_10,hDCA_10b,1.,1.);
//    hDCA_10->Add(hDCA_10,hDCA_10c,1.,1.);
//    hDCA_10->Add(hDCA_10,hDCA_10d,1.,1.);
    if(hDCA_10->GetEntries()>0){
        xmin=hDCA_10->GetMean()-3*hDCA_10->GetRMS();
        xmax=hDCA_10->GetMean()+3*hDCA_10->GetRMS();
        hDCA_10->Fit("h","NR,Q","",xmin,xmax);
//     dcamean[3]=hDCA_10->GetMean();
//     err_dcamean[3]=hDCA_10->GetMeanError();
//     dcaRMS[3]=hDCA_10->GetRMS();
//     err_dcaRMS[3]=hDCA_10->GetRMSError();
        dcamean[3]=h->GetParameter(1);
        err_dcamean[3]=h->GetParError(1);
        dcaRMS[3]=h->GetParameter(2);
        err_dcaRMS[3]=h->GetParError(2);
    }

    //
// DCAz
//
    TH1F *hDCAz_05 = (TH1F *)fListDCA->FindObject("d0allpointzRec_5");
    if(hDCAz_05->GetEntries()>0){
        xmin=hDCAz_05->GetMean()-3*hDCAz_05->GetRMS();
        xmax=hDCAz_05->GetMean()+3*hDCAz_05->GetRMS();
        hDCAz_05->Fit("h","NR,Q","",xmin,xmax);
        //     dcamean[0]=hDCA_05->GetMean();
        //     err_dcamean[0]=hDCA_05->GetMeanError();
        //     dcaRMS[0]=hDCA_05->GetRMS();
        //     err_dcaRMS[0]=hDCA_05->GetRMSError();
        dcazmean[0]=h->GetParameter(1);
        err_dcazmean[0]=h->GetParError(1);
        dcazRMS[0]=h->GetParameter(2);
        err_dcazRMS[0]=h->GetParError(2);
        
    }
    
    TH1F *hDCAz_1 = (TH1F *)fListDCA->FindObject("d0allpointzRec_9");
    if(hDCAz_1->GetEntries()>0){
        xmin=hDCAz_1->GetMean()-3*hDCAz_1->GetRMS();
        xmax=hDCAz_1->GetMean()+3*hDCAz_1->GetRMS();
        hDCAz_1->Fit("h","NR,Q","",xmin,xmax);
        //     dcamean[1]=hDCA_1->GetMean();
        //     err_dcamean[1]=hDCA_1->GetMeanError();
        //     dcaRMS[1]=hDCA_1->GetRMS();
        //     err_dcaRMS[1]=hDCA_1->GetRMSError();
        dcazmean[1]=h->GetParameter(1);
        err_dcazmean[1]=h->GetParError(1);
        dcazRMS[1]=h->GetParameter(2);
        err_dcazRMS[1]=h->GetParError(2);
    }
    
    TH1F *hDCAz_5 = (TH1F *)fListDCA->FindObject("d0allpointzRec_18");
    if(hDCAz_5->GetEntries()>0){
        xmin=hDCAz_5->GetMean()-3*hDCAz_5->GetRMS();
        xmax=hDCAz_5->GetMean()+3*hDCAz_5->GetRMS();
        hDCAz_5->Fit("h","NR,Q","",xmin,xmax);
        //     dcamean[2]=hDCA_5->GetMean();
        //     err_dcamean[2]=hDCA_5->GetMeanError();
        //     dcaRMS[2]=hDCA_5->GetRMS();
        //     err_dcaRMS[2]=hDCA_5->GetRMSError();
        dcazmean[2]=h->GetParameter(1);
        err_dcazmean[2]=h->GetParError(1);
        dcazRMS[2]=h->GetParameter(2);
        err_dcazRMS[2]=h->GetParError(2);
    }
    
    TH1F *hDCAz_10 = (TH1F *)fListDCA->FindObject("d0allpointzRec_20");
    //    TH1F *hDCA_10b = (TH1F *)fListDCA->FindObject("d0allpointrphiRec_21");
    //    TH1F *hDCA_10c = (TH1F *)fListDCA->FindObject("d0allpointrphiRec_22");
    //    TH1F *hDCA_10d = (TH1F *)fListDCA->FindObject("d0allpointrphiRec_23");
    //    hDCA_10->Add(hDCA_10,hDCA_10b,1.,1.);
    //    hDCA_10->Add(hDCA_10,hDCA_10c,1.,1.);
    //    hDCA_10->Add(hDCA_10,hDCA_10d,1.,1.);
    if(hDCAz_10->GetEntries()>0){
        xmin=hDCAz_10->GetMean()-3*hDCAz_10->GetRMS();
        xmax=hDCAz_10->GetMean()+3*hDCAz_10->GetRMS();
        hDCAz_10->Fit("h","NR,Q","",xmin,xmax);
        //     dcamean[3]=hDCA_10->GetMean();
        //     err_dcamean[3]=hDCA_10->GetMeanError();
        //     dcaRMS[3]=hDCA_10->GetRMS();
        //     err_dcaRMS[3]=hDCA_10->GetRMSError();
        dcazmean[3]=h->GetParameter(1);
        err_dcazmean[3]=h->GetParError(1);
        dcazRMS[3]=h->GetParameter(2);
        err_dcazRMS[3]=h->GetParError(2);
    }

    
    
    Int_t index=0;
    xntDCA[index++]=(Float_t)nrun;
    for(Int_t ij=0;ij<4;ij++){
        xntDCA[index++]=dcamean[ij];
//        cout<< "dcamean["<<ij<<"] = " << dcamean[ij] << endl;
        xntDCA[index++]=err_dcamean[ij];
//        cout<< "err_dcamean["<<ij<<"] = " << err_dcamean[ij] << endl;
        xntDCA[index++]=dcaRMS[ij];
//        cout <<"dcaRMS["<<ij<<"] = " << dcaRMS[ij] << endl;
        xntDCA[index++]=err_dcaRMS[ij];
//        cout<< "err_dcaRMS["<<ij<<"] = " << err_dcaRMS[ij] << endl;
        xntDCA[index++]=dcazmean[ij];
//        cout<< "dcazmean["<<ij<<"] = " << dcazmean[ij] << endl;
        xntDCA[index++]=err_dcazmean[ij];
//        cout<< "err_dcazmean["<<ij<<"] = " << err_dcazmean[ij] << endl;
        xntDCA[index++]=dcazRMS[ij];
//        cout <<"dcazRMS["<<ij<<"] = " << dcazRMS[ij] << endl;
        xntDCA[index++]=err_dcazRMS[ij];
//        cout << "err_dcazRMS["<<ij<<"] = " << err_dcazRMS[ij] << endl;
    }
    
    nt->Fill(xntDCA);
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
//    printf("Run %d\n",noRuns[i]);
    if(noRuns[i]>=run1 && noRuns[i]<=run2){
//      printf("Accepting run number %d in position %d\n",noRuns[i],kRunsToPlot);
      kRunsToPlot++;
    }
    else { 
//      printf("Rejecting run number %d - out of range\n",noRuns[i]);
      noRuns[i]=run2+10; 
    }
  }
  TMath::Sort(nr,noRuns,myIndex,kFALSE);
  printf("Total number of runs accepted for display %d\n",kRunsToPlot);
  if(kRunsToPlot==0)return 0;
//  for(Int_t i=0;i<kRunsToPlot;i++)printf("Position %d ) Run: %d\n",i,noRuns[myIndex[i]]);
  return kRunsToPlot;
}
