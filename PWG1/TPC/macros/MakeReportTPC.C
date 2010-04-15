/*
   //gSystem->AddIncludePath("-I$ALICE_ROOT/PWG1/TPC");
   //.L $ALICE_ROOT/PWG1/TPC/macros/MakeReportTPC.C+
   //
   //MakeReportTPC();

   gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
   gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+")
   AliXRDPROOFtoolkit tool;
   TChain * chain = tool.MakeChain("summaryTPCQA.txt","tpcQA",0,200000);
*/


#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TSystem.h"
#include "TFile.h"
#include "TGrid.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "THnSparse.h"
#include "TTreeStream.h"
#include "AliPerformanceTPC.h"
#endif

TObject *pTPCObject=0;
TTreeSRedirector  *pcstream=0;
void Init();
void ReadObjects(const char * path=0);
void MakeReportTPC(Int_t run, const char *path=0 );
void AnalyzeNCL();
void AnalyzeDrift();

void Init(){
  //
  //
  //
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libTENDER.so");
  gSystem->Load("libCORRFW.so");
  gSystem->Load("libPWG0base.so");
  gSystem->Load("libPWG0dep.so");
  gSystem->Load("libPWG0selectors.so");
  gSystem->Load("libPWG1.so");
  gSystem->Load("libTPCcalib"); 
  gSystem->Setenv("alien_CLOSE_SE","ALICE::GSI::SE");
  TGrid::Connect("alien://",0,0,"t"); 
  gSystem->Setenv("alien_CLOSE_SE","ALICE::GSI::SE");    
}


void ReadObjects(const char * path){
  //
  //
  //
  TFile *f =0;
  if (path) f=TFile::Open(Form("%s/TPC.PerformanceQAQA.root",path));
  if (!path) f=TFile::Open("TPC.Performance.root");
  if (!f){
    printf("File %s not available\n", path);
    return;
  }
  TList *list = (TList*)f->Get("TPC");
  if (!list){
    printf("QA %s not available\n", path);
    return;
  }
  pTPCObject= (AliPerformanceTPC*)list->FindObject("AliPerformanceTPC");
}


void MakeReportTPC(Int_t run, const char *path ){
  //
  // make a tpcQA report
  //  typical variables exported to the tree
  //
  Init();
  ReadObjects(path);
  
  pcstream= new TTreeSRedirector("TPC.PerformanceSummary.root");
  (*pcstream)<<"tpcQA"<<"run="<<run;
  AnalyzeNCL();
  AnalyzeDrift();
  (*pcstream)<<"tpcQA"<<"\n";
  delete pcstream;
}


void AnalyzeNCL(){
  //
  // NCL statistic
  //    
  // variables:
  static Double_t meanTPCnclF=0;
  static Double_t rmsTPCnclF=0;
  static Double_t slopeATPCnclF=0;
  static Double_t slopeCTPCnclF=0;
  static Double_t meanTPCncl=0;
  static Double_t rmsTPCncl=0;
  static Double_t slopeATPCncl=0;
  static Double_t slopeCTPCncl=0;
  AliPerformanceTPC * pTPC=  (AliPerformanceTPC *)pTPCObject;

  TH1* his1D=0;
  TProfile* hprof=0;
  static TF1 *fpol1 = new TF1("fpol1","pol1");
  //
  // all clusters
  // eta cut - +-1
  // pt cut  - +-0.250 GeV
  pTPC->GetTPCTrackHisto()->GetAxis(5)->SetRangeUser(-1.,1.);
  pTPC->GetTPCTrackHisto()->GetAxis(7)->SetRangeUser(0.25,10);
  his1D = pTPC->GetTPCTrackHisto()->Projection(0);
  meanTPCncl= his1D->GetMean();
  rmsTPCncl= his1D->GetRMS();
  delete his1D;
  hprof = pTPC->GetTPCTrackHisto()->Projection(0,5)->ProfileX();
  hprof->Fit(fpol1,"QNR","QNR",0,0.8);
  slopeATPCncl= fpol1->GetParameter(1);
  hprof->Fit(fpol1,"QNR","QNR",-0.8,0.0);
  slopeCTPCncl= fpol1->GetParameter(1);
  delete hprof;
  //
  // findable clusters
  //
  pTPC->GetTPCTrackHisto()->GetAxis(2)->SetRangeUser(0.4,1.1);
  his1D = pTPC->GetTPCTrackHisto()->Projection(2);
  meanTPCnclF= his1D->GetMean();
  rmsTPCnclF= his1D->GetRMS();
   delete his1D;
  his1D = pTPC->GetTPCTrackHisto()->Projection(2,5)->ProfileX();
  his1D->Fit(fpol1,"QNR","QNR",0,0.8);
  slopeATPCnclF= fpol1->GetParameter(1);
  his1D->Fit(fpol1,"QNR","QNR",-0.8,0.0);
  slopeCTPCnclF= fpol1->GetParameter(1);
  delete his1D;
  printf("Cluster QA report\n");
  printf("meanTPCnclF=\t%f\n",meanTPCnclF);
  printf("rmsTPCnclF=\t%f\n",rmsTPCnclF);
  printf("slopeATPCnclF=\t%f\n",slopeATPCnclF);
  printf("slopeCTPCnclF=\t%f\n",slopeCTPCnclF);
  printf("meanTPCncl=\t%f\n",meanTPCncl);
  printf("rmsTPCncl=\t%f\n",rmsTPCncl);
  printf("slopeATPCncl=\t%f\n",slopeATPCncl);
  printf("slopeCTPCncl=\t%f\n",slopeCTPCncl);
  //
  // dump results to the tree
  //
  (*pcstream)<<"tpcQA"<<
    "meanTPCnclF="<<meanTPCnclF <<
    "rmsTPCnclF="<<rmsTPCnclF <<
    "slopeATPCnclF="<< slopeATPCnclF<<
    "slopeCTPCnclF="<< slopeCTPCnclF<<
    "meanTPCncl="<<meanTPCncl <<
    "rmsTPCncl="<< rmsTPCncl<<
    "slopeATPCncl="<< slopeATPCncl<<
    "slopeCTPCncl="<< slopeCTPCncl;

    
}

void AnalyzeDrift(){
  //
  // Analyze drift velocity imperfections
  //
  // variables:
  static Double_t offsetdZA=0;
  static Double_t slopedZA=0;
  static Double_t offsetdZC=0;
  static Double_t slopedZC=0;
  AliPerformanceTPC * pTPC=  (AliPerformanceTPC *)pTPCObject;

  TH1* his1D=0;
  TH2* his2D=0;
  static TF1 *fpol1 = new TF1("fpol1","pol1");
  TObjArray arrayFit;
  his2D = pTPC->GetTPCTrackHisto()->Projection(4,5);
  his2D->FitSlicesY(0,0,-1,10,"QNR",&arrayFit);
  delete his2D;
  his1D = (TH1*) arrayFit.At(1);
  his1D->Fit(fpol1,"QNR","QNR",-0.8,-0.1);
  offsetdZC=fpol1->GetParameter(0);
  slopedZC=fpol1->GetParameter(1);
  his1D->Fit(fpol1,"QNR","QNR",0.1,0.8);
  offsetdZA=fpol1->GetParameter(0);
  slopedZA=fpol1->GetParameter(1);
  //
  printf("Drift velocity QA report\n");
  printf("offsetdZA\t%f\n",offsetdZA);
  printf("slopedZA\t%f\n",slopedZA);
  printf("offsetdZC\t%f\n",offsetdZC);
  printf("slopedZC\t%f\n",slopedZC);
  //
  // dump drift QA values
  //
  (*pcstream)<<"tpcQA"<<
    "offsetdZA="<< offsetdZA<<
    "slopedZA="<< slopedZA<<
    "offsetdZC="<< offsetdZC<<
    "slopedZC="<<slopedZC;
  
}
