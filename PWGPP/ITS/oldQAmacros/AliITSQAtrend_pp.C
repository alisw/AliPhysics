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
//#include <AliITSgeomTGeo.h>
#endif

TString pdfFileNames="";
void MakePlot(Int_t run1=-1,Int_t run2=999999,TString ntupleFileName="TrendingITS.root");
void PlotITSSA(TFile *fil,Int_t run1, Int_t run2);
void FillITSSAntuple(TFile* f,TNtuple* nt, Int_t nrun);
void FillSDDntuple(TFile* f,TNtuple* nt, Int_t iRun, Float_t *xnt);
void FillSSDntuple(TFile* f,TNtuple* ntssd, Int_t iRun, Float_t *xntSSD);
void FillMatchntuple(TFile* f,TNtuple* ntmatching, Int_t iRun, Float_t *xntMatching);
void FillVTXntuple(TFile* f,TNtuple* ntvertex, Int_t iRun, Float_t *xntVertex);
void AliITSQAtrend(TString runListFile="LHC12d.txt",TString ntupleFileName="TrendingITS.root");
Double_t LangausFun(Double_t *x, Double_t *par);
Bool_t WriteInputTextFileFromMonalisaListOfRuns(TString outtxtfilename,Int_t* listofrunsfromMonalisa,Int_t nruns,TString pathbeforRunN="alice/data/2012/LHC12a/",TString pathafterRunN="cpass1/QAresults_barrel.root");
Int_t RunsToPlot(Int_t run1, Int_t run2,Int_t nr,Int_t *noRuns, Int_t *myIndex);

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
/alice/data/2012/LHC12d/000184780/ESDs/high_lumi/QAresults.root 2848
*/
//
//   Function MakePlot(run1,run2):
//   it produces the plots. For each canvas a PDF file is created.
//   A PDF file with all the canvases merged is also produced
//   The first two argument define a range for the runs to be displayed
//   These two arguments are optional: by default, all the runs
//   found in the ntuples are displayed
////////////////////////////////////////////////////////////////

/* $Id$ */

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
  Float_t MPVdEdxLay3,errMPVdEdxLay3,MPVdEdxLay4,errMPVdEdxLay4;
  Float_t MPVdEdxTB0,errMPVdEdxTB0,MPVdEdxTB5,errMPVdEdxTB5;
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
  ntsdd->SetBranchAddress("MPVdEdxTB0",&MPVdEdxTB0);
  ntsdd->SetBranchAddress("errMPVdEdxTB0",&errMPVdEdxTB0);
  ntsdd->SetBranchAddress("MPVdEdxTB5",&MPVdEdxTB5);
  ntsdd->SetBranchAddress("errMPVdEdxTB5",&errMPVdEdxTB5);
  ntsdd->SetBranchAddress("MPVdEdxLay3",&MPVdEdxLay3);
  ntsdd->SetBranchAddress("errMPVdEdxLay3",&errMPVdEdxLay3);
  ntsdd->SetBranchAddress("MPVdEdxLay4",&MPVdEdxLay4);
  ntsdd->SetBranchAddress("errMPVdEdxLay4",&errMPVdEdxLay4);

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

  }
  
  //-----------------------------------

  //---vertex

  TNtuple* ntvertex=(TNtuple*)fil->Get("ntvertex");

  Float_t nrunVertex,Vx,errVx,sigmaVx,errsigmaVx,Vy,errVy,sigmaVy,errsigmaVy,Vz,errVz,sigmaVz,errsigmaVz,VxSPD,errVxSPD,sigmaVxSPD,errsigmaVxSPD,VySPD,errVySPD,sigmaVySPD,errsigmaVySPD,VzSPD,errVzSPD,sigmaVzSPD,errsigmaVzSPD;

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
    
  }
 //--------  Draw Vertex histograms ---------
  TCanvas *cVertexDisto=new TCanvas("cVertexDisto","cVertexDisto",1200,800);
  cVertexDisto->Divide(3,2);
  cVertexDisto->cd(1);
 //  hVx->SetMinimum(0.03);
//   hVx->SetMaximum(0.08);
  hVx->GetYaxis()->SetTitle("Vertex X coordinate");
  hVx->Draw();
  hVxSPD->Draw("same");
  cVertexDisto->cd(2);
//   hVy->SetMinimum(0.25);
//   hVy->SetMaximum(0.30);
  hVy->GetYaxis()->SetTitle("Vertex Y coordinate");
  hVy->Draw();
  hVySPD->Draw("same");
  cVertexDisto->cd(3);
//   hVz->SetMinimum(-1.);
//   hVz->SetMaximum(1.);
  hVz->GetYaxis()->SetTitle("Vertex Z coordinate");
  hVz->Draw();
  hVzSPD->Draw("same");
  cVertexDisto->cd(4);
//   hSigmaVx->SetMinimum(0.);
//   hSigmaVx->SetMaximum(0.01);
  hSigmaVx->GetYaxis()->SetTitle("sigma on x coordinate");
  hSigmaVx->Draw();
  hSigmaVxSPD->Draw("same");
  cVertexDisto->cd(5);
//   hSigmaVy->SetMinimum(0.);
//   hSigmaVy->SetMaximum(0.01);
  hSigmaVy->GetYaxis()->SetTitle("sigma on y coordinate");
  hSigmaVy->Draw();
  hSigmaVySPD->Draw("same");
  cVertexDisto->cd(6);
//   hSigmaVz->SetMinimum(6.);
//   hSigmaVz->SetMaximum(10.);
  hSigmaVz->GetYaxis()->SetTitle("sigma on z coordinate");
  hSigmaVz->Draw();
  hSigmaVzSPD->Draw("same");
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
  histodEdxLay6->SetLineColor(7);
  histodEdxLay6->SetMarkerColor(7);
  histodEdxLay6->SetMarkerStyle(24);
  histodEdxLay6->Draw("same");
    
  histodEdxLay3->GetYaxis()->SetTitle("MPV of dE/dx (keV/300 #mum)");
  
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
  cpt02->SaveAs("TPCITS_trend.pdf");
      pdfFileNames+=" TPCITS_trend.pdf";
  cpt02->Update();

  //-----------  ITS SA -----------
  PlotITSSA(fil,run1,run2);
 // merge the pdf files
  TString command("gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=merged");
  command=command+"ITS_trend.pdf "+pdfFileNames;
  gSystem->Exec(command.Data());
  printf(" Merging the pdf file:  %s \n",command.Data());
  delete [] myIndex;
  delete [] noRuns;
}

//____________________________________________________________________________
void PlotITSSA(TFile *fil, Int_t run1, Int_t run2){
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
    h0->GetXaxis()->SetBinLabel(iev+1,Form("%.0f",run));
    h1->GetXaxis()->SetBinLabel(iev+1,Form("%.0f",run));
    h2->GetXaxis()->SetBinLabel(iev+1,Form("%.0f",run));
    h3->GetXaxis()->SetBinLabel(iev+1,Form("%.0f",run));
    h4->GetXaxis()->SetBinLabel(iev+1,Form("%.0f",run));
    h5->GetXaxis()->SetBinLabel(iev+1,Form("%.0f",run));
    h6->GetXaxis()->SetBinLabel(iev+1,Form("%.0f",run));
    h7->GetXaxis()->SetBinLabel(iev+1,Form("%.0f",run));
    h8->GetXaxis()->SetBinLabel(iev+1,Form("%.0f",run));
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
  TLatex* ti1=new TLatex(0.11,0.40,"ITS pure SA tracks (normalized to number of events)");
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
  TLatex* ti3=new TLatex(0.10,0.8,"(ITSTPC+ITSsa)/ITSpureSA ");
  ti3->SetNDC();
  ti3->Draw();
   c3->SaveAs("tracks_ratio_trend.pdf");
      pdfFileNames+=" tracks_ratio_trend.pdf";

}



//____________________________________________________________________________
void AliITSQAtrend(TString runListFile,TString ntupleFileName){

  TGrid::Connect("alien://");
  

  
  //-----------SDD
  
  const Int_t nVariables=43;
  TNtuple* ntsdd=new TNtuple("ntsdd","SDD trending","nrun:fracTrackWithClu1:errfracTrackWithClu1:fracTrackWithClu2:errfracTrackWithClu2:fracTrackWithClu3:errfracTrackWithClu3:fracTrackWithClu4:errfracTrackWithClu4:fracTrackWithClu5:errfracTrackWithClu5:fracTrackWithClu6:errfracTrackWithClu6:meanTrPts3:errmeanTrPts3:meanTrPts4:errmeanTrPts4:minDrTime:errminDrTime:meanDrTime:errmeanDrTime:fracExtra:errfracExtra:meandEdxLay3:errmeandEdxLay3:meandEdxLay4:errmeandEdxLay4:meandEdxTB0:errmeandEdxTB0:meandEdxTB5:errmeandEdxTB5:MPVdEdxLay3:errMPVdEdxLay3:MPVdEdxLay4:errMPVdEdxLay4:MPVdEdxTB0:errMPVdEdxTB0:MPVdEdxTB5:errMPVdEdxTB5:nMod95:nMod80:nMod60:nModEmpty");
  Float_t xnt[nVariables];
  
  //--------------SSD
    
  const Int_t nVariablesSSD=14;
  TNtuple* ntssd=new TNtuple("ntssd","SSD trending","nrun:meandEdxLay5:errmeandEdxLay5:meandEdxLay6:errmeandEdxLay6:MPVdEdxLay5:errMPVdEdxLay5:MPVdEdxLay6:errMPVdEdxLay6:ChargeRatioL5:errChargeratioL5:ChargeRatioL6:errChargeratioL6:moduleOff");
  Float_t xntSSD[nVariablesSSD];
  
  //----Matching

  const Int_t nVariablesMatching=60;
  TNtuple* ntmatching=new TNtuple("ntmatching","Matching Efficiency","nrun:FracSPD1:errFracSPD1:FracSPD2:errFracSPD2:Eff6Pt02:errEff6Pt02:Eff6Pt1:errEff6Pt1:Eff6Pt10:errEff6Pt10:Eff5Pt02:errEff5Pt02:Eff5Pt1:errEff5Pt1:Eff5Pt10:errEff5Pt10:Eff4Pt02:errEff4Pt02:Eff4Pt1:errEff4Pt1:Eff4Pt10:errEff4Pt10:Eff3Pt02:errEff3Pt02:Eff3Pt1:errEff3Pt1:Eff3Pt10:errEff3Pt10:Eff2Pt02:errEff2Pt02:Eff2Pt1:errEff2Pt1:Eff2Pt10:errEff2Pt10:EffSPDPt02:errEffSPDPt02:EffSPDPt1:errEffSPDPt1:EffSPDPt10:errEffSPDPt10:EffoneSPDPt02:errEffoneSPDPt02:EffoneSPDPt1:errEffoneSPDPt1:EffoneSPDPt10:errEffoneSPDPt10:EffTOTPt02:errEffTOTPt02:EffTOTPt1:errEffTOTPt1:EffTOTPt10:errEffTOTPt10");

  Float_t xntMatching[nVariablesMatching];
  
  //--------------------------------

  //----QA Vertex

  const Int_t nVariablesVertex=27;
  TNtuple* ntvertex=new TNtuple("ntvertex","QA Vertex","nrun:VxTRK:errVxTRK:sigmaVxTRK:errsigmaVxTRK:VyTRK:errVyTRK:sigmaVyTRK:errsigmaVyTRK:VzTRK:errVzTRK:sigmaVzTRK:errsigmaVzTRK:VxSPD:errVxSPD:sigmaVxSPD:errsigmaVxSPD:VySPD:errVySPD:sigmaVySPD:errsigmaVySPD:VzSPD:errVzSPD:sigmaVzSPD:errsigmaVzSPD");
      
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
  Int_t fillNumb[MAX_LINES];
  Bool_t goout = kFALSE;
  while ( in ) {
    in.getline(strings[j], MAX_LINE_LEN);
    TString aux(strings[j]);
    TString auxrun(strings[j]);
    TString auxfill(strings[j]);
    TString aux2(strings[j]);
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
    else if(aux.Contains("LHC12b/")){
      lentrail = 27;
      lenfill=36;
    }
    else if(aux.Contains("LHC12c/")  || aux.Contains("LHC12d/") ){
      lentrail = 27;
      lenfill=36;
    }
    else if(aux.Contains("LHC12a/")){
      lentrail = 27;
      lenfill=36;
    }
  
    
    else {
      if(!aux.IsNull())printf("Unrecognised path name %s \n",aux.Data());
      goout = kTRUE;
    }
    if(goout)break;
    if(aux.Length()<lentrail)continue;
    auxrun=aux.Remove(0,lentrail);
    auxfill=aux.Remove(0,lenfill);
    runNumb[j]=atoi(auxrun.Data());
    fillNumb[j]=atoi(auxfill.Data());
    aux2.Remove(aux2.Length()-5);
    printf("%d ) - fill %d - path %s \n",runNumb[j],fillNumb[j],aux2.Data());
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

    if(!gGrid||!gGrid->IsConnected()) {
      printf("gGrid not found! exit macro\n");
      return;
    }
    TString aux(strings[jru]);
    aux.Remove(aux.Length()-5);
    TFile *f=TFile::Open(Form("alien://%s",aux.Data())); 
    if(!f) {
      printf("File not found, continue with next one\n");
      filenotfound++;
      continue;
    }
    
    //-------SDD
    FillSDDntuple(f,ntsdd,iRun,xnt);
  

    //-------SSD

    FillSSDntuple(f,ntssd,iRun,xntSSD);

 
    //--------------matching

    FillMatchntuple(f,ntmatching,iRun,xntMatching);

    //------------- Vertex
    FillVTXntuple(f,ntvertex,iRun,xntVertex);
	
  
    //---------------------------

    //--------------   ITS SA ---------------
    cout<<"ITS - SA"<<endl;
    //    cout<<f<<" "<<ntSA<<" "<<iRun<<endl;
    FillITSSAntuple(f,ntSA,iRun);

  } // loop on runs

  printf("%d runs skipped because QA file not found\n",filenotfound);
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
    //    if(NTPCITS[ibin]!=0 && NITSsa[ibin]!=0)Ratio[ibin]=NTPCITS[ibin]/NITSsa[ibin];
    Double_t totaltrks=NTPCITS[ibin]+NITSsa[ibin];
  if(totaltrks!=0 && NITSpureSA[ibin]!=0 )Ratio[ibin]=totaltrks/NITSpureSA[ibin];
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

//_____________________________________________________________________________
void FillSDDntuple(TFile* f,TNtuple* nt, Int_t iRun, Float_t *xnt){
  TDirectoryFile* df=(TDirectoryFile*)f->Get("SDD_Performance");
  if(!df){
    printf("Run %d SDD_Performance MISSING -> Exit\n",iRun);
    return;
  } 
   
    
  TList* l=(TList*)df->Get("coutputRP");
  if(!l){
    printf("Run %d coutputRP TList MISSING -> Exit\n",iRun);
    return;
  }  
    


  cout<<"SDD - QA"<<endl;

  TH1F* hcllay=(TH1F*)l->FindObject("hCluInLay");
    
  if(hcllay->GetEntries()==0){

    printf("Run %d hcllay EMPTY -> Return\n",iRun);
    return;
      
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
    return;
  }
    
  TH1F* hgamod=(TH1F*)l->FindObject("hGAMod");
    
  if(hgamod->GetEntries()==0){
    printf("Run %d hgamod EMPTY -> Continue\n",iRun);
    return;
  }

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

  if(hev->GetEntries()==0){
    printf("Run %d hev EMPTY -> Continue\n",iRun);
    return;
  }
    
  Int_t nTotEvents=hev->GetBinContent(2);
  Int_t nTrigEvents=hev->GetBinContent(3);
  Int_t nEvents=nTotEvents;
  printf("Run %d Number of Events = %d Triggered=%d\n",iRun,nTotEvents,nTrigEvents);
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
    return;
  }


  TH1F* hgpmod=(TH1F*)l->FindObject("hGoodPmod");
  if(hgpmod->GetEntries()==0){
    printf("Run %d hgpmod EMPTY -> Continue\n",iRun);
    return;
  }
      
  //     TH1F* hmpmod=(TH1F*)l->FindObject("hMissPmod");
  TH1F* hbrmod=(TH1F*)l->FindObject("hBadRegmod");
  if(hbrmod->GetEntries()==0){
    printf("Run %d hbrmod EMPTY -> Continue\n",iRun);
    return;
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
    return;
  }

  TH1F* htimTe=(TH1F*)l->FindObject("hDrTimTPExtra");
      
  if(htimTe->GetEntries()==0){
    printf("Run %d htimTe EMPTY -> Continue\n",iRun);
    return;
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
    return;
  }

  TH1D* hdedxLay3=hdedxmod->ProjectionY("hdedxLay3",1,84);
  TH1D* hdedxLay4=hdedxmod->ProjectionY("hdedxLay4",85,260);
      
  TH1F* hSigTim0=(TH1F*)l->FindObject("hSigTimeInt0");
  if(hSigTim0->GetEntries()==0){
    printf("Run %d hSigTim0 EMPTY -> Continue\n",iRun);
    return;
  }

  TH1F* hSigTim5=(TH1F*)l->FindObject("hSigTimeInt5");
  if(hSigTim5->GetEntries()==0){
    printf("Run %d hSigTim5  EMPTY -> Continue\n",iRun);
    return;
  }
  //Fitting the same distributions in order to have the MPV
  TF1 *lfunLay3 = new TF1("LangausFunLay3",LangausFun,50.,300.,4); 
  lfunLay3->SetParameter(0,5.);
  lfunLay3->SetParameter(1,80.);
  lfunLay3->SetParameter(2,hdedxLay3->GetEntries()/10.);
  lfunLay3->SetParameter(3,10.);
  lfunLay3->SetParLimits(3,0.,20);
  hdedxLay3->Fit(lfunLay3,"NQLR");
  TF1 *lfunLay4 = new TF1("LangausFunLay4",LangausFun,50.,300.,4); 
  lfunLay4->SetParameter(0,5.);
  lfunLay4->SetParameter(1,80.);
  lfunLay4->SetParameter(2,hdedxLay4->GetEntries()/10.);
  lfunLay4->SetParameter(3,10.);
  lfunLay4->SetParLimits(3,0.,20);
  hdedxLay4->Fit(lfunLay4,"NQLR");
  TF1 *lfunTim0 = new TF1("LangausFunTim0",LangausFun,50.,300.,4); 
  lfunTim0->SetParameter(0,5.);
  lfunTim0->SetParameter(1,80.);
  lfunTim0->SetParameter(2,hSigTim0->GetEntries()/10.);
  lfunTim0->SetParameter(3,10.);
  lfunTim0->SetParLimits(3,0.,20);
  hSigTim0->Fit(lfunTim0,"NQLR");
  TF1 *lfunTim5 = new TF1("LangausFunTim5",LangausFun,50.,300.,4); 
  lfunTim5->SetParameter(0,5.);
  lfunTim5->SetParameter(1,80.);
  lfunTim5->SetParameter(2,hSigTim5->GetEntries()/10.);
  lfunTim5->SetParameter(3,10.);
  lfunTim5->SetParLimits(3,0.,20);
  hSigTim5->Fit(lfunTim5,"NQLR");
     
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
  xnt[index++]=lfunLay3->GetParameter(1);
  xnt[index++]=lfunLay3->GetParError(1);
  xnt[index++]=lfunLay4->GetParameter(1);
  xnt[index++]=lfunLay4->GetParError(1);
  xnt[index++]=lfunTim0->GetParameter(1);
  xnt[index++]=lfunTim0->GetParError(1);
  xnt[index++]=lfunTim5->GetParameter(1);
  xnt[index++]=lfunTim5->GetParError(1);
  xnt[index++]=(Float_t)nBelow95;
  xnt[index++]=(Float_t)nBelow80;
  xnt[index++]=(Float_t)nBelow60;
  xnt[index++]=(Float_t)nZeroP;
  nt->Fill(xnt);
  

}

//_____________________________________________________________________________
void FillSSDntuple(TFile* f,TNtuple* ntssd, Int_t iRun, Float_t *xntSSD){

  cout<<"SSD - QA"<<endl;

  TDirectoryFile* dfSSD=(TDirectoryFile*)f->Get("PWGPPdEdxSSDQA");
  if(!dfSSD){
    printf("Run %d SSD_Performance MISSING -> Exit\n",iRun);
    return;
  }
      
  TList* lSSD=(TList*)dfSSD->Get("SSDdEdxQA");
  if(!dfSSD){
    printf("Run %d coutputRP TList MISSING -> Exit\n",iRun);
    return;
  }
  //
      
  TH2F* QAchargeRatio=(TH2F*)lSSD->FindObject("QAChargeRatio");
      
  if(QAchargeRatio->GetEntries()==0){
    printf("Run %d QAchargeRatio EMPTY -> Return\n",iRun);
    return;
  }

  TH2F* QAcharge=(TH2F*)lSSD->FindObject("QACharge");
      
  if(QAcharge->GetEntries()==0){
    printf("Run %d QAcharge EMPTY -> Return\n",iRun);
    return;
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
      
  Int_t indexSSD=0;
  xntSSD[indexSSD++]=iRun;
  xntSSD[indexSSD++]=(Float_t)hChargeL5->GetMean();
  xntSSD[indexSSD++]=(Float_t)hChargeL5->GetMeanError();
  xntSSD[indexSSD++]=(Float_t)hChargeL6->GetMean();
  xntSSD[indexSSD++]=(Float_t)hChargeL6->GetMeanError();
  xntSSD[indexSSD++]=(Float_t)lfunLay5->GetParameter(1);
  xntSSD[indexSSD++]=(Float_t)lfunLay5->GetParError(1);
  xntSSD[indexSSD++]=(Float_t)lfunLay6->GetParameter(1);
  xntSSD[indexSSD++]=(Float_t)lfunLay6->GetParError(1);
  xntSSD[indexSSD++]=(Float_t)hChargeRatioL5->GetMean();
  xntSSD[indexSSD++]=(Float_t)hChargeRatioL5->GetMeanError();
  xntSSD[indexSSD++]=(Float_t)hChargeRatioL6->GetMean();
  xntSSD[indexSSD++]=(Float_t)hChargeRatioL6->GetMeanError();
  xntSSD[indexSSD++]=(Float_t)contEmpty;
  ntssd->Fill(xntSSD);						
}

//_____________________________________________________________________________
void FillMatchntuple(TFile* f,TNtuple* ntmatching, Int_t iRun, Float_t *xntMatching){
    cout<<"Tracking"<<endl;

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
      printf("Run %d hFiredChip EMPTY -> Return\n",iRun);
      return;
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

    TH1F *fHistPtTPCInAcc = (TH1F*)list->FindObject("fHistPtTPCInAcc");

    if(fHistPtTPCInAcc->GetEntries()==0){
      printf("Run %dfHistPtTPCInAcc  EMPTY -> Return\n",iRun);
      return;
    }

    TH1F *fHistPtITSMI6InAcc = (TH1F*)list->FindObject("fHistPtITSMI6InAcc");

    if(fHistPtITSMI6InAcc->GetEntries()==0){
      printf("Run %d fHistPtITSMI6InAcc EMPTY -> Return\n",iRun);
      return;
    }

    TH1F *fHistPtITSMI5InAcc = (TH1F*)list->FindObject("fHistPtITSMI5InAcc");
    if(fHistPtITSMI5InAcc->GetEntries()==0){
      printf("Run %d fHistPtITSMI5InAcc EMPTY -> Return\n",iRun);
      return;
    }
      
    TH1F *fHistPtITSMI4InAcc = (TH1F*)list->FindObject("fHistPtITSMI4InAcc");
    if(fHistPtITSMI5InAcc->GetEntries()==0){
      printf("Run %d fHistPtITSMI5InAcc EMPTY -> Return\n",iRun);
      return;
    }
      
    TH1F *fHistPtITSMI3InAcc = (TH1F*)list->FindObject("fHistPtITSMI3InAcc");
    if(fHistPtITSMI3InAcc->GetEntries()==0){
      printf("Run %d fHistPtITSMI3InAcc EMPTY -> Return\n",iRun);
      return;
    }
    TH1F *fHistPtITSMI2InAcc = (TH1F*)list->FindObject("fHistPtITSMI2InAcc");
      
    if(fHistPtITSMI2InAcc->GetEntries()==0){
      printf("Run %d fHistPtITSMI2InAcc EMPTY -> Return\n",iRun);
      return;
    }
    TH1F *fHistPtITSMISPDInAcc = (TH1F*)list->FindObject("fHistPtITSMISPDInAcc");
    if(fHistPtITSMISPDInAcc->GetEntries()==0){
      printf("Run %d fHistPtITSMISPDInAcc EMPTY -> Return\n",iRun);
      return;
    }
      
    TH1F *fHistPtITSMIoneSPDInAcc = (TH1F*)list->FindObject("fHistPtITSMIoneSPDInAcc");
      
    if(fHistPtITSMIoneSPDInAcc->GetEntries()==0){
      printf("Run %d fHistPtITSMIoneSPDInAcc  EMPTY -> Return\n",iRun);
      return;
    }

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
}

//_____________________________________________________________________________
void FillVTXntuple(TFile* f,TNtuple* ntvertex, Int_t iRun, Float_t *xntVertex){
   cout<<"Primary Vertex"<<endl;

    TDirectoryFile *dirVertex = (TDirectoryFile*)f->Get("Vertex_Performance");
    if(!dirVertex){
      Printf("Vertex directory not found... check!");
    }
  
    TList *lt = (TList*)dirVertex->Get("cOutputVtxESD");
	
    TH1F *xVtxTRK = (TH1F*)lt->FindObject("fhTRKVertexX");

    if(xVtxTRK->GetEntries()==0){
      printf("Run %d xVtxTRK EMPTY -> Return\n",iRun);
      return;
    }

    TH1F *yVtxTRK = (TH1F*)lt->FindObject("fhTRKVertexY");

    if(yVtxTRK->GetEntries()==0){
      printf("Run %d yVtxTRK EMPTY -> Return\n",iRun);
      return;
    }

    TH1F *zVtxTRK = (TH1F*)lt->FindObject("fhTRKVertexZ");

    if(zVtxTRK->GetEntries()==0){
      printf("Run %d zVtxTRK EMPTY -> Return\n",iRun);
      return;
    }

 TH1F *xVtxSPD = (TH1F*)lt->FindObject("fhSPDVertexX");

    if(xVtxSPD->GetEntries()==0){
      printf("Run %d xVtxSOD EMPTY -> Return\n",iRun);
      return;
    }

    TH1F *yVtxSPD = (TH1F*)lt->FindObject("fhSPDVertexY");

    if(yVtxSPD->GetEntries()==0){
      printf("Run %d yVtxSPD EMPTY -> Return\n",iRun);
      return;
    }

    TH1F *zVtxSPD = (TH1F*)lt->FindObject("fhSPDVertexZ");

    if(zVtxSPD->GetEntries()==0){
      printf("Run %d zVtxSPD EMPTY -> Return\n",iRun);
      return;
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
    ntvertex->Fill(xntVertex);	

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
  printf("Total number of runs accepted fot display %d\n",kRunsToPlot);
  if(kRunsToPlot==0)return 0;
  for(Int_t i=0;i<kRunsToPlot;i++)printf("Position %d ) Run: %d\n",i,noRuns[myIndex[i]]);
  return kRunsToPlot;
}

