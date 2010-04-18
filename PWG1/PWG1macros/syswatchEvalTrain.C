/*
  Make a Analysis/summary of syswatch.root files created in Analysis train 

  Summary created:  Memory usage - per train/task
                    CPU usage    - per train/task
                    IO load      - per hostName
 
  Input:    syswatch.txt          - list of syswatch.root files
  Output:   syswatchSummary.root  - file with defaut pictures
                                  - and summary information of forms of trees


 //  
 // input list of syswatch files to create a chain
 // Local creation e.g using:
 //

 .L $ALICE_ROOT/PWG1/PWG1macros/syswatchEvalTrain.C

*/


#if !defined(__CINT__) || defined(__MAKECINT__)
#include <fstream>
#include "TSystem.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TChain.h"
#include "TH1.h"
#include "TCut.h"
#include "THashList.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMap.h"
#include "TTreeStream.h"
#endif



TChain * chainWatch = 0; 
TTreeSRedirector * pcstream = 0; 
TObjArray * picArray    = new TObjArray;
TObjArray * taskSummary = new TObjArray;
Int_t maxEntries   = 100000;   // maximal amount of points to visualize

//
// export variables
//
Double_t ftimeSlope = 0;   // time slope per event
Double_t fmemSlope  = 0;   // mem slope  per event
Double_t fmemOffset = 0;   // mem offset at the initialization
Double_t fmeandVM   = 0;   // mean delta VM per task
Double_t fmeandT    = 0;   // mean delta T  per task
TCut cutStep1="1";   // skipt the points if too many
//
// functions:
//
void AddToChain(TString inputList);
void SetAliasTrain();
void GetSummaryTrain();
void GetSummaryTask();
void GetSummaryHost();
void PrintPairSummary(TObjArray *array);


void syswatchEvalTrain(){
  AddToChain("syswatch.txt");
  SetAliasTrain();
  Int_t allEntries = chainWatch->GetEntries();
  if (allEntries>maxEntries) cutStep1 = Form("((%s%d)==0)","Entry$%",allEntries/maxEntries+1);
  pcstream = new TTreeSRedirector("syswatchSummary.root"); 
  pcstream->GetFile()->cd();
  //
  GetSummaryTrain();
  GetSummaryTask();
  GetSummaryHost();
  //picArray->Write("Pictures");
  //taskSummary->Write("TaskSummary",TObjectArray::kSingleKey);
  delete pcstream;
}

void PrintPairSummary(TObjArray *array){
  
  for (Int_t id1=0; id1<array->GetEntries(); id1++){
    TPair *pair= (TPair*)array->At(id1);
    printf("%s\t%s\n",pair->GetName(),pair->Value()->GetName());    
  }
}



void AddToChain(TString inputList){
  //
  // add the files form inputList into chain
  //
  if (!chainWatch) chainWatch = new TChain("syswatch");
  ifstream in;
  in.open(inputList.Data());
  Int_t counter=0;
  TString currentFile;
  while(in.good()) {
    in >> currentFile;
    if (!currentFile.Contains(".root")) continue;
    chainWatch->Add(currentFile.Data());
    printf("%d\t%s\n",counter,currentFile.Data());
    counter++;
  }
}



void SetAliasTrain(){
  //
  // Time, Memory consumption and host Id  stored for given snapshot
  // Snapshots are identified using 3 ids: id0, id1 and id2
  // Each sanphsot is identeified also using snapshot name
  //
  // Meaning of ids depends on the type of ananlysis
  //
  // Alieses according train conventions
  // id conventions in the ANALYSIS is following:
  //    id1 - task id number
  //        - 1000 reserved for the Handlers_BeginEvent
  chainWatch->SetAlias("EventNumber","id0");
  chainWatch->SetAlias("isTask","id1<200&&id1>=0");
  chainWatch->SetAlias("isReading","id1==1000");
  chainWatch->SetAlias("isEvent","id0>10");  // skip first events
  chainWatch->SetAlias("hostName","hname");
  chainWatch->SetAlias("stampName","sname");
}

void GetSummaryTrain(){
  //
  // Get the Summary info per train
  //  
  TF1 ftime("ftime","[0]*x");
  TF1 fmem("fmem","pol1");
  ftime.SetLineColor(2);
  fmem.SetLineColor(2);
  TGraph *gr = 0;
  Int_t entries=0;
  //
  //
  //
  entries = chainWatch->Draw("T:EventNumber*0.001",cutStep1+"isTask&&isEvent","goff");
  gr=new TGraphErrors(entries, chainWatch->GetV2(),chainWatch->GetV1(),0,0);
  gr->SetName("TvEventNumber");
  gr->SetTitle("T:EventNumber");
  gr->GetXaxis()->SetTitle("Event Nr. (10^3)");
  gr->GetYaxis()->SetTitle("T (s)");
  gr->Fit(&ftime,"ROB=0.7");
  gr->Draw("ap");
  picArray->AddLast(gr->Clone());
  pcstream->GetFile()->cd();
  gr->Write();
  //
  entries = chainWatch->Draw("0.001*VM:0.001*EventNumber",cutStep1+"isTask&&isEvent");
  gr=new TGraphErrors(entries, chainWatch->GetV2(),chainWatch->GetV1(),0,0);
  gr->SetName("VMvEventNumber");
  gr->SetTitle("VM:EventNumber");
  gr->GetXaxis()->SetTitle("Event Nr. (10^3)");
  gr->GetYaxis()->SetTitle("Virtual memory (GBy)");
  gr->Fit(&fmem,"ROB=0.7");
  gr->Draw("ap");
  picArray->AddLast(gr->Clone());
  pcstream->GetFile()->cd();
  gr->Write();

  fmemSlope=fmem.GetParameter(1);
  fmemOffset=fmem.GetParameter(0);
  ftimeSlope=ftime.GetParameter(0);
  picArray->AddLast(gr->Clone());
  pcstream->GetFile()->cd();
  gr->Write();
  (*pcstream)<<"summaryInfo"<<
    "ftimeSlope="<<ftimeSlope<<     // time slope per event
    "fmemSlope="<<fmemSlope<<       // mem slope  per event
    "fmemOffset="<<fmemOffset<<     // mem offset at the initialization
    "fmeandVM="<<fmeandVM<<         // mean delta VM per task
    "fmeandT="<<fmeandT<<           // mean delta T  per task
    "\n";
}

void GetSummaryTask(){
  //
  //
  //
  TH1 * hisVM =0;
  TH1 * hisT =0;

  chainWatch->Draw("1000*deltaVM:sname",cutStep1+"isTask&&isEvent","profile");
  hisVM = (TH1*)chainWatch->GetHistogram()->Clone();
  hisVM->SetMarkerStyle(20);
  hisVM->SetName("deltaVMperTask");
  hisVM->SetTitle("delta VM per Task");
  hisVM->GetYaxis()->SetTitle("#DeltaVM/Event (kBy)");
  hisVM->Draw();
  picArray->AddLast(hisVM);
  pcstream->GetFile()->cd();
  hisVM->Write();
  //
  chainWatch->Draw("1000*deltaT:sname",cutStep1+"isTask&&isEvent","profile");
  hisT = (TH1*)chainWatch->GetHistogram()->Clone();
  hisT->SetMarkerStyle(20);
  hisT->SetName("deltaTperTask");
  hisT->SetTitle("deltaT per Task");
  hisT->GetYaxis()->SetTitle("#Delta_{t}/Event (ms)");
  hisT->Draw();
  picArray->AddLast(hisT);
  pcstream->GetFile()->cd();
  hisT->Write();
  //
  //
  Int_t nbins = hisT->GetXaxis()->GetNbins();
  for (Int_t iname=0; iname<nbins; iname++){
    TObjString *hname= new TObjString(hisVM->GetXaxis()->GetLabels()->At(iname)->GetName());
    printf("%s\t%f\t%f\n",hname->GetName(), hisVM->GetBinContent(iname+1),hisT->GetBinContent(iname+1)); 
    Double_t vmev = hisVM->GetBinContent(iname+1);
    Double_t tev  = hisT->GetBinContent(iname+1);
    Double_t vmevErr = hisVM->GetBinError(iname+1);
    Double_t tevErr  = hisT->GetBinError(iname+1);
    char hstring[1000];
    sprintf(hstring,"%s",hname->GetName());
    (*pcstream)<<"taskInfo"<<
      "taskName.="<<hname<<   // host name
      "vmev="<< vmev<<          // memory per task per even
      "tev="<<  tev<<           // time per event per task
      "vmevErr="<< vmevErr<<          // memory per task per even
      "teverr="<<  tevErr<<           // time per event per task
      // MEAN summary (ALL)
      "ftimeSlope="<<ftimeSlope<<     // time slope per event
      "fmemSlope="<<fmemSlope<<       // mem slope  per event
      "fmemOffset="<<fmemOffset<<     // mem offset at the initialization
      "fmeandVM="<<fmeandVM<<         // mean delta VM per task
      "fmeandT="<<fmeandT<<           // mean delta T  per task
      "\n";
  }
}


void GetSummaryHost(){
  //
  // Get summary information per hosts
  // 2 histograms
  // 
  Int_t nbins=0;
  TH1 * hisIOT=0;
  TH1 * hisT=0;
  chainWatch->Draw("0.0000001*fileBytesRead/T:hostName","EventNumber>0&&isReading","prof");
  hisIOT = (TH1*)chainWatch->GetHistogram()->Clone();
  hisIOT->SetDirectory(0);
  chainWatch->Draw("1000*T/EventNumber:hostName","EventNumber>0&&isReading","prof");
  hisT = (TH1*)chainWatch->GetHistogram()->Clone();
  hisT->SetDirectory(0);
  //
  hisIOT->SetTitle("input MBy/s per host");
  hisIOT->SetName("input-MBysPerHost");
  hisIOT->SetTitle("MBy/s");
  hisIOT->GetYaxis()->SetTitle("MBy/s");
  hisIOT->SetMarkerStyle(22);
  hisIOT->Draw();
  picArray->AddLast(hisIOT);
  pcstream->GetFile()->cd();
  hisIOT->Write();
  //
  hisT->SetTitle("Event/ms per host");
  hisT->SetName("Event/ms per Host");
  hisT->GetYaxis()->SetTitle("Events/ms");
  hisT->SetMarkerStyle(22);
  hisT->Draw();
  picArray->AddLast(hisT);
  pcstream->GetFile()->cd();
  hisT->Write();
  //
  nbins = hisIOT->GetXaxis()->GetNbins();
  for (Int_t iname=0; iname<nbins; iname++){
    TObjString *hname= new TObjString(hisIOT->GetXaxis()->GetLabels()->At(iname)->GetName());
    printf("%s\t%f\t%f\n",hname->GetName(), hisIOT->GetBinContent(iname+1),hisT->GetBinContent(iname+1)); 
    Double_t iot= hisIOT->GetBinContent(iname+1);
    Double_t tev= hisT->GetBinContent(iname+1);
    char hstring[1000];
    sprintf(hstring,"%s",hname->GetName());
    (*pcstream)<<"hostInfo"<<
      "hostName.="<<hname<<   // host name
      "iot="<< iot<<          // reading per second
      "tev="<< tev<<          // events per milisecond 
      // mean summary (all)
      "ftimeSlope="<<ftimeSlope<<     // time slope per event
      "fmemSlope="<<fmemSlope<<       // mem slope  per event
      "fmemOffset="<<fmemOffset<<     // mem offset at the initialization
      "fmeandVM="<<fmeandVM<<         // mean delta VM per task
      "fmeandT="<<fmeandT<<           // mean delta T  per task
      "\n";
  }

}
