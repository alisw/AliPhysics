//--------------------------------------------------------------------------
// Macro for QA monitoring.
//
// In case it is not run with full aliroot, it needs the following libraries to compile:
//  - libSTEERBase.so
//  - libESD.so
//  - libAOD.so
//  - libANALYSIS.so
//  - libANALYSISalice.so
//  - libCORRFW.so
//  - libPWG3muon.so
//
// TString includePath = "-I${ALICE_ROOT}/PWG3/base/ ";  gSystem->SetIncludePath(includePath.Data());
//
// The macro reads results of the QA task and produce monitoring plots.
//
// Author: Philippe Pillot - SUBATECH Nantes
// Modified by Christophe Suire, Cynthia Hadjidakis - IPN Orsay
//--------------------------------------------------------------------------

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>

// ROOT includes
#include "TEnv.h"
#include "TMath.h"
#include "TGrid.h"
#include "TGridResult.h"
#include "THashList.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TF1.h"

// ALIROOT includes
#include "AliCounterCollection.h"

#endif


// .x PlotMuonQApp.C("alien:///alice/cern.ch/user/c/cynthia/muon/QA/pp/LHC10e/pass2_test/results",0,kFALSE)
// .x PlotMuonQApp.C("/Users/cynthia/Documents/alice/data/MuonQA/results",0,kFALSE)
// .x PlotMuonQApp.C("/Users/cynthia/Documents/alice/data/MuonQA/results","/Users/cynthia/Documents/alice/data/MuonQA/LHC10e/pass2/runlist_period3_test3_3runs.txt",kFALSE)
//--------------------------------------------------------------------------
void PlotMuonQApp(const char* baseDir, const char* runList = 0x0, const char * triggerList = 0x0, Bool_t selectPhysics = kFALSE, const char *LHCPeriod = "LHC11c", const char *QAFileName = "AnalysisResults.root")
{
  /// Macro for QA monitoring.
  /// Example: baseDir = "alien:///alice/cern.ch/user/p/ppillot/pp7TeV/LHC10d/MuonQA/pass1/results/".
  /// If runList != 0x0: only the given runs will be used. Otherwise use all runs found in baseDir.
  
#if defined(__CINT__) && !defined(__MAKECINT__)
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWG3base");
  gSystem->Load("libPWG3muon");
#endif
  
  // Cosmetics and configuration
  gStyle->SetFillColor(10);
  gStyle->SetPadGridX(kTRUE);
  gStyle->SetPadGridY(kTRUE);
  gStyle->SetPadRightMargin(0.01);
  
  Int_t colorInd[10]={1,4,2,3,6};
	
  TString OutFileName = "QA_";  OutFileName += LHCPeriod;
  TString OutFileNamePDF=  OutFileName.Data();  OutFileNamePDF+= ".pdf";
  TString OutFileNamePDF_open = OutFileNamePDF.Data(); OutFileNamePDF_open += "[";  
  TString OutFileNamePDF_close= OutFileNamePDF.Data(); OutFileNamePDF_close += "]";  
  TString OutFileNameROOT=  OutFileName.Data();  OutFileNameROOT+= ".root";
  
  Int_t PRINTSTAT = 1;
  //Trigger in the trigger list for matched track asymmetry and run per run statistics. 
  Int_t kCMUS = 2;//2;
	
  if (0){ // Equivalent to the fast read option
    gEnv->SetValue("XNet.ConnectTimeout",10);
    gEnv->SetValue("XNet.RequestTimeout",10);
    gEnv->SetValue("XNet.MaxRedirectCount",2);
    gEnv->SetValue("XNet.ReconnectTimeout",10);
    gEnv->SetValue("XNet.FirstConnectMaxCnt",1);
  }
  
  TH1::AddDirectory(kFALSE);
  
  TString alienBaseDir = baseDir;
  
  
  if (alienBaseDir.Contains("alien:") && !TGrid::Connect("alien://")) {
    Error("MergeQA","cannot connect to grid");
    return;
  }
  
  
  
  //---------------------------------- //
  //          Run selection            //
  //---------------------------------- //
  
  
  // list runs to be analyzed
  TString selectRuns = "run:";
  TObjArray runs, runs2;
  runs.SetOwner();

  if (runList) {
    // only the ones in the runList
    ifstream inFile(runList);
    if (!inFile.is_open()) {
      Error("PlotMuonQApp","unable to open file %s", runList);
      return;
    }
    
    TString currRun;
    while (!inFile.eof()) {
      currRun.ReadLine(inFile, kTRUE);
      if (currRun.IsNull()) continue;
      if (!currRun.IsDigit()) {
	Error("PlotMuonQApp","invalid run number: %s", currRun.Data());
	return;
      }
      runs.AddLast(new TObjString(Form("%09d", currRun.Atoi())));
			runs2.AddLast(new TObjString(Form("%d", currRun.Atoi())));
      selectRuns += Form("%s,",currRun.Data());
    }
    selectRuns.Remove(TString::kTrailing, ',');
    inFile.close();
    
  } else {
    // all runs
    runs.AddLast(new TObjString("*"));
		runs2.AddLast(new TObjString("*"));
  }
  
  // physics selection
  TString select = selectPhysics ? "selected:yes" : "";
  
  cout<<"//---------------------------------- //"<<endl;
  cout<<"//        Trigger selection          //"<<endl;
  cout<<"//---------------------------------- //"<<endl;

  TObjArray triggers, triggersB, triggersAC, triggersE;
  TString selectAllTriggers = "", selectAllTriggersAC = "", selectAllTriggersE = "";
  triggers.SetOwner();

  TString sTrigList = triggerList;
  if ( ! sTrigList.IsNull() ) {
    
    // only the ones in the triggerList
    ifstream inFile(triggerList);
    if (!inFile.is_open()) {
      Error("PlotMuonQApp","unable to open file %s", triggerList);
      return;
    }
    
    TString currTrig;
    while (!inFile.eof()) {
      currTrig.ReadLine(inFile, kTRUE);
      if (currTrig.IsNull()) continue;
      if (!currTrig.IsAlnum()) {
        Error("PlotMuonQApp","invalid trigger name: %s", currTrig.Data());
        return;
      }
      triggers.AddLast(new TObjString(currTrig));
    }
    inFile.close();
  } else {
    // by default all triggers from new period in LHC11c
    triggers.AddLast(new TObjString("CINT7"));
    triggers.AddLast(new TObjString("CMUSH7"));
    //triggers.AddLast(new TObjString("CMUS7"));
    triggers.AddLast(new TObjString("CMUL7"));
    triggers.AddLast(new TObjString("CMUU7"));
  }
  
  cout<<" Nr of triggers read "<<triggers.GetEntriesFast()<<endl;
  for(Int_t i = 0; i < triggers.GetEntriesFast(); i++){
    TString triggerName = ( (TObjString*) triggers.At(i) )->GetString();
    if(triggerName == "ANY" || triggerName =="OTHER" || triggerName == "CINT7I"){
      selectAllTriggers += Form("%s,", triggerName.Data());
      selectAllTriggersAC += Form("%s,", triggerName.Data());//fake
      selectAllTriggersE += Form("%s,", triggerName.Data());//fake
      continue;
    }
    selectAllTriggers += Form("%sB,", triggerName.Data());
    selectAllTriggersAC += Form("%sAC,", triggerName.Data());
    selectAllTriggersE += Form("%sE,", triggerName.Data());
  }
  selectAllTriggers.Remove(TString::kTrailing, ',')	;
  selectAllTriggersAC.Remove(TString::kTrailing, ',')	;
  selectAllTriggersE.Remove(TString::kTrailing, ',')	;
  
  for(Int_t i = 0; i <= triggers.GetEntriesFast(); i++){
    TString triggerName = "";
    if(i==0){
      triggersB.AddLast(new TObjString(Form("%s", selectAllTriggers.Data())));
      triggersAC.AddLast(new TObjString(Form("%s", selectAllTriggersAC.Data())));
      triggersE.AddLast(new TObjString(Form("%s", selectAllTriggersE.Data())));		
    }
    else{
      triggerName = ( (TObjString*) triggers.At(i-1) )->GetString();
      if(triggerName == "ANY" || triggerName =="OTHER" || triggerName == "CINT7I"){
	triggersB.AddLast(new TObjString(Form("%s,", triggerName.Data())));
	triggersAC.AddLast(new TObjString(Form("%s", triggerName.Data())));//fake
	triggersE.AddLast(new TObjString(Form("%s", triggerName.Data())));//fake
	continue;
      }
      triggersB.AddLast(new TObjString(Form("%sB", triggerName.Data())));
      triggersAC.AddLast(new TObjString(Form("%sAC", triggerName.Data())));
      triggersE.AddLast(new TObjString(Form("%sE", triggerName.Data())));
    }
    cout<<i<<" "<< ( (TObjString*) triggersB.At(i) )->GetString()<<endl;
    
  }
  
  cout<<"//---------------------------------- //"<<endl;
  cout<<"//        plot global counter        //"<<endl;
  cout<<"//---------------------------------- //"<<endl;
  
  TFile *globalFile = TFile::Open(Form("%s/%s", baseDir,QAFileName));
  if (!globalFile || ! globalFile->IsOpen()) {
    Error("PlotQA", "failed to open file: %s/%s", baseDir, QAFileName);
    return;
  }
  globalFile->Cd("MUON_QA");
  
  TString selection;
	
  // get counters
  AliCounterCollection* eventCounters = static_cast<AliCounterCollection*>(globalFile->FindObjectAny("eventCounters"));
  AliCounterCollection* trackCounters = static_cast<AliCounterCollection*>(globalFile->FindObjectAny("trackCounters"));
  if (!runList) selectRuns += trackCounters->GetKeyWords("run");
  
  //declare a default canvas c1 
  TString CanvasName = "c1";
  TCanvas *c1 = new TCanvas(CanvasName.Data(),CanvasName.Data());
  c1->cd();
  
  TH1* hBNoPS[10]={}; 
  TH1* hBWithPS[10]={};
  TH1* hB[10]={};
  TH1* hTriggerB[10], *hTrackerB[10], *hMatchedB[10], *hAllTracksB[10], *hMatchedLowPtB[10], *hMatchedHighPtB[10];
  TH1* hMatchedLowPtBNoPS[10], *hMatchedHighPtBNoPS[10];
  TH1* hPosMatchedB[10], *hNegMatchedB[10], *hAllMatchedB[10];
  TH1 *hACWithPS[10]={}; 
  TH1 *hACNoPS[10]={};
  //TH1 *hEWithPS[10]={};
  TH1 *hENoPS[10]={};
  
  if(triggers.GetEntriesFast()>=10){
    cout<<"Too many triggers"<<endl;
    return;
  }
  //Loop on trigger (first is all triggers, then for each defined trigger)
  for(Int_t i = 0; i <= triggers.GetEntriesFast(); i++){
    
    TString triggerName = ( (TObjString*) triggersB.At(i) )->GetString();
    
    // Histo trigger without Phys. Sel. 
    selection = Form("trigger:%s/%s", triggerName.Data(), selectRuns.Data());
    cout<<selection<<endl;
    hBNoPS[i] = (TH1*) eventCounters->Draw("run",selection);
    if(!hBNoPS[i]) return;
    hBNoPS[i]->Sumw2();
    hBNoPS[i]->LabelsOption("a");
    
    TString triggerNameAC = ( (TObjString*) triggersAC.At(i) )->GetString();
    // Histo trigger without Phys. Sel. AC
    selection = Form("trigger:%s/%s", triggerNameAC.Data(), selectRuns.Data());
    cout<<selection<<endl;
    hACNoPS[i] = (TH1*) eventCounters->Draw("run",selection);
    if(!hACNoPS[i]) hACNoPS[i] = new TH1D("hACNoPS","",10,0,10);
    else{
      hACNoPS[i]->Sumw2();
      hACNoPS[i]->LabelsOption("a");
    }
    // Histo trigger with Phys. Sel. AC
    selection = Form("trigger:%s/%s/selected:yes", triggerNameAC.Data(), selectRuns.Data());
    cout<<selection<<endl;
    hACWithPS[i] = (TH1*) eventCounters->Draw("run",selection);
    if(!hACWithPS[i]) hACWithPS[i] = new TH1D("hACWithPS","",10,0,10);
    else{
      hACWithPS[i]->Sumw2();
      hACWithPS[i]->LabelsOption("a");
    }
    /*TString triggerNameE = ( (TObjString*) triggersE.At(i) )->GetString();
    // Histo trigger without Phys. Sel. E
    selection = Form("trigger:%s/%s", triggerNameE.Data(), selectRuns.Data());
    cout<<selection<<endl;
    hENoPS[i] = (TH1*) eventCounters->Draw("run",selection);
    if(!hENoPS[i]) hENoPS[i] = new TH1D("hENoPS","",10,0,10);
    else{
    hENoPS[i]->Sumw2();
    hENoPS[i]->LabelsOption("a");
    }
    // Histo trigger with Phys. Sel. E
    selection = Form("trigger:%s/%s/selected:yes", triggerNameE.Data(), selectRuns.Data());
    cout<<selection<<endl;
    hEWithPS[i] = (TH1*) eventCounters->Draw("run",selection);
    if(!hEWithPS[i]) hEWithPS[i]=new TH1D("hEWithPS","",10,0,10);
    else{
    hEWithPS[i]->Sumw2();
    hEWithPS[i]->LabelsOption("a");
    }*/
    
    // Histo trigger with Phys. Sel. 
    selection = Form("trigger:%s/%s/selected:yes", triggerName.Data(), selectRuns.Data());
    cout<<selection<<endl;
    hBWithPS[i] = (TH1*) eventCounters->Draw("run",selection);
    if(!hBWithPS[i]) return;
    hBWithPS[i]->Sumw2();
    hBWithPS[i]->LabelsOption("a");
    
    // Histo trigger : Phys. Sel.  is selected or not depending on the macro arguments
    selection = Form("trigger:%s/%s/%s", triggerName.Data(), selectRuns.Data(), select.Data());
    hB[i] = (TH1*) eventCounters->Draw("run", selection.Data());
    hB[i]->Sumw2();
    hB[i]->LabelsOption("a");	
    
    // Histo tracking : Phys. Sel.  is selected or not depending on the macro arguments
    selection = Form("track:triggeronly/trigger:%s/%s/%s", triggerName.Data(), selectRuns.Data(), select.Data());
    hTriggerB[i] = (TH1*) trackCounters->Draw("run",selection);
    hTriggerB[i]->Sumw2();
    hTriggerB[i]->LabelsOption("a");
    selection = Form("track:trackeronly/trigger:%s/%s/%s", triggerName.Data(), selectRuns.Data(), select.Data());
    hTrackerB[i] = (TH1*) trackCounters->Draw("run",selection);
    hTrackerB[i]->Sumw2();
    hTrackerB[i]->LabelsOption("a");
    selection = Form("track:matched/trigger:%s/%s/%s", triggerName.Data(), selectRuns.Data(), select.Data());
    hMatchedB[i] = (TH1*) trackCounters->Draw("run",selection);
    hMatchedB[i]->Sumw2();
    hMatchedB[i]->LabelsOption("a");
    selection = Form("trigger:%s/%s/%s", triggerName.Data(), selectRuns.Data(), select.Data());
    hAllTracksB[i] = (TH1*) trackCounters->Draw("run",selection);
    hAllTracksB[i]->Sumw2();  
    hAllTracksB[i]->LabelsOption("a");
    selection = Form("track:matched/trigger:%s/%s/%s/pt:low/acc:in", triggerName.Data() ,selectRuns.Data(), select.Data());
    hMatchedLowPtB[i] = (TH1*) trackCounters->Draw("run", selection);
    hMatchedLowPtB[i]->Sumw2();  
    hMatchedLowPtB[i]->LabelsOption("a");
    selection = Form("track:matched/trigger:%s/%s/%s/pt:high/acc:in", triggerName.Data() ,selectRuns.Data(), select.Data());
    hMatchedHighPtB[i] = (TH1*) trackCounters->Draw("run", selection);
    hMatchedHighPtB[i]->Sumw2();  
    hMatchedHighPtB[i]->LabelsOption("a");
    selection = Form("track:matched/trigger:%s/%s/pt:low/acc:in", triggerName.Data() ,selectRuns.Data());
    hMatchedLowPtBNoPS[i] = (TH1*) trackCounters->Draw("run", selection);
    hMatchedLowPtBNoPS[i]->Sumw2();  
    hMatchedLowPtBNoPS[i]->LabelsOption("a");
    selection = Form("track:matched/trigger:%s/%s/pt:high/acc:in", triggerName.Data() ,selectRuns.Data());
    hMatchedHighPtBNoPS[i] = (TH1*) trackCounters->Draw("run", selection);
    hMatchedHighPtBNoPS[i]->Sumw2();  
    hMatchedHighPtBNoPS[i]->LabelsOption("a");
    selection = Form("track:matched/trigger:%s/%s/charge:pos/%s/acc:in",triggerName.Data(), select.Data(),selectRuns.Data());
    hPosMatchedB[i] =  (TH1*) trackCounters->Draw("run",selection);
    hPosMatchedB[i]->Sumw2();  
    hPosMatchedB[i]->LabelsOption("a");
    selection = Form("track:matched/trigger:%s/%s/charge:neg/%s/acc:in",triggerName.Data(), select.Data(),selectRuns.Data());
    hNegMatchedB[i] =  (TH1*) trackCounters->Draw("run",selection);
    hNegMatchedB[i]->Sumw2();  
    hNegMatchedB[i]->LabelsOption("a");
    selection = Form("track:matched/trigger:%s/%s/%s/acc:in",triggerName.Data(), select.Data(),selectRuns.Data());
    hAllMatchedB[i] =  (TH1*) trackCounters->Draw("run",selection);
    hAllMatchedB[i]->Sumw2();  
    hAllMatchedB[i]->LabelsOption("a"); 
  }
	
  TH1* hAll = (TH1*) eventCounters->Draw("trigger","run",Form("run:any/%s",select.Data()));
  hAll->LabelsOption("a");
  hAll->Draw("TEXT");
  
  Int_t NumOfBNoPS[10];
  Int_t NumOfBWithPS[10];
  Int_t NumOfACNoPS[10];
  Int_t NumOfENoPS[10];
  Int_t NumOfACWithPS[10];
  Int_t NumOfEWithPS[10];
  for(Int_t i = 0; i <= triggers.GetEntriesFast(); i++){
    NumOfBNoPS[i] = hBNoPS[i]->Integral();
    NumOfBWithPS[i] = hBWithPS[i]->Integral();
    NumOfACNoPS[i] = 0;//hACNoPS[i]->Integral();
    NumOfENoPS[i] = 0;//hENoPS[i]->Integral();
    NumOfACWithPS[i] = 0;//hACWithPS[i]->Integral();
    NumOfEWithPS[i] = 0;//hEWithPS[i]->Integral();
  }
  
  cout<<"//==================================================================================="<<endl;
  cout<<"// Put all plots in a ps file, easier to publish (Twiki)"<<endl;
  cout<<"//==================================================================================="<<endl;
  
  c1->Print(OutFileNamePDF_open.Data());
  
  TFile *rootFileOut = TFile::Open(OutFileNameROOT.Data(),"RECREATE");
  
  cout<<"//==================================================================================="<<endl;
  cout<<"// new canvas with the total number of trigger"<<endl;
  cout<<"//==================================================================================="<<endl;
  CanvasName =  LHCPeriod ; 
  CanvasName += "_TriggerContent"; 
  TCanvas *cTriggerContent = new TCanvas(CanvasName.Data(),"cTriggerContent",1200,900);
  cTriggerContent->SetTopMargin(0.05);
  cTriggerContent->SetRightMargin(0.01);
  cTriggerContent->SetGridy(1);
  cTriggerContent->cd();
  
  TLegend* legcTC = new TLegend(0.2,0.15,0.50,0.40);
  legcTC->SetHeader("Physics Selection");
  legcTC->AddEntry(".","applied :","");
  
  for(Int_t i = 0; i <= triggers.GetEntriesFast(); i++){
    if(i==0) continue;
    if(i==1){
      hBWithPS[i]->SetMinimum(0.001);
      hBWithPS[i]->Draw();
      hBWithPS[i]->GetYaxis()->SetTitle("Trigger content w/ Phys. Sel."); 
    }
    else hBWithPS[i]->Draw("same");
    hBWithPS[i]->SetLineColor(colorInd[i]);
    legcTC->AddEntry(hBWithPS[i],(( (TObjString*) triggersB.At(i) )->GetString()).Data(),"l");
  }
  legcTC->Draw("same");
  
  cTriggerContent->Print(OutFileNamePDF.Data());
  cTriggerContent->Write();
  
  cout<<"//==================================================================================="<<endl;
  cout<<"// new canvas with the relative content of each trigger w/ and w/o  physics selection"<<endl;
  cout<<"//==================================================================================="<<endl;
  CanvasName =  LHCPeriod ; 
  CanvasName += "_RelativeTriggerContent"; 
  TCanvas *cRelativeTriggerContent = new TCanvas(CanvasName.Data(),"cRelativeTriggerContent",1200,900);
  cRelativeTriggerContent->SetTopMargin(0.05);
  cRelativeTriggerContent->SetRightMargin(0.01);
  cRelativeTriggerContent->SetGridy(1);
  cRelativeTriggerContent->SetLogy(1);
  cRelativeTriggerContent->cd();
  
  
  TH1* ratioB[10], *ratioBNoPS[10];
  TH1*ratioACNoPS[10];
  TH1*ratioENoPS[10];
  TLegend* legcRTC = new TLegend(0.2,0.15,0.50,0.40);
  legcRTC->SetHeader("Physics Selection");
  
  TString hName, hTriggerName;
  for(Int_t i = 0; i <= triggers.GetEntriesFast(); i++){
    hName = "ratio";
    hName += ( (TObjString*) triggersB.At(i) )->GetString();
    ratioB[i] = static_cast<TH1*> (hBWithPS[i]->Clone(hName));
    ratioB[i]->Divide(hBNoPS[0]);
    ratioB[i]->SetLineWidth(2);
    ratioB[i]->SetLineColor(colorInd[i]);
    hName = "ratioNoPS";
    hName += ( (TObjString*) triggersB.At(i) )->GetString();
    ratioBNoPS[i] = static_cast<TH1*> (hBNoPS[i]->Clone(hName));
    ratioBNoPS[i]->Divide(hBNoPS[0]);
    ratioBNoPS[i]->SetLineWidth(0);
    ratioBNoPS[i]->SetLineStyle(2);
    ratioBNoPS[i]->SetMarkerStyle(24+i);
    ratioBNoPS[i]->SetMarkerSize(1);
    ratioBNoPS[i]->SetLineColor(colorInd[i]);
    ratioBNoPS[i]->SetMarkerColor(colorInd[i]);
    if(NumOfACNoPS[0]>0){
      hName = "ratioACNoPS";
      hName += ( (TObjString*) triggersAC.At(i) )->GetString();
      ratioACNoPS[i] = static_cast<TH1*> (hACNoPS[i]->Clone(hName));
      ratioACNoPS[i]->Divide(hACNoPS[0]);
      ratioACNoPS[i]->SetLineWidth(0);
      ratioACNoPS[i]->SetLineStyle(3);
      ratioACNoPS[i]->SetMarkerStyle(24+i);
      ratioACNoPS[i]->SetMarkerSize(1);
      ratioACNoPS[i]->SetLineColor(colorInd[i]);
      ratioACNoPS[i]->SetMarkerColor(colorInd[i]);
    }
    if(NumOfENoPS[0]>0){
      hName = "ratioENoPS";
      hName += ( (TObjString*) triggersE.At(i) )->GetString();
      ratioENoPS[i] = static_cast<TH1*> (hENoPS[i]->Clone(hName));
      ratioENoPS[i]->Divide(hENoPS[0]);
      ratioENoPS[i]->SetLineWidth(0);
      ratioENoPS[i]->SetLineStyle(3);
      ratioENoPS[i]->SetMarkerStyle(24+i);
      ratioENoPS[i]->SetMarkerSize(1);
      ratioENoPS[i]->SetLineColor(colorInd[i]);
      ratioENoPS[i]->SetMarkerColor(colorInd[i]);
    }
    
    if(i==0) continue;
    else if(i==1){
      ratioB[i]->SetMaximum(1.5);
      ratioB[i]->SetMinimum(0.001);
      ratioB[i]->SetLabelSize(0.02);
      ratioB[i]->GetYaxis()->SetTitle("Relative trigger content w/ and w/o Phys. Sel."); 
      ratioB[i]->Draw("E");
      ratioBNoPS[i]->Draw("EPSAME");
      if(NumOfACNoPS[0]>0) ratioACNoPS[i]->Draw("EPSAME");
      if(NumOfENoPS[0]>0)ratioENoPS[i]->Draw("EPSAME");
    }
    else{
      ratioB[i]->Draw("ESAME");
      ratioBNoPS[i]->Draw("EPSAME");
      if(NumOfACNoPS[0]>0) ratioACNoPS[i]->Draw("EPSAME");
      if(NumOfENoPS[0]>0) ratioENoPS[i]->Draw("EPSAME");
    }
  }
  
  legcRTC->AddEntry(".","applied :","");
  for(Int_t i = 1; i <= triggers.GetEntriesFast(); i++){
    legcRTC->AddEntry(ratioB[i],(( (TObjString*) triggersB.At(i) )->GetString()).Data(),"l");
  }
  legcRTC->AddEntry(".","not applied :","");
  for(Int_t i = 1; i <= triggers.GetEntriesFast(); i++){
    legcRTC->AddEntry(ratioBNoPS[i],(( (TObjString*) triggersB.At(i) )->GetString()).Data(),"p");
    if(NumOfACNoPS[i]) legcRTC->AddEntry(ratioACNoPS[i],(( (TObjString*) triggersAC.At(i) )->GetString()).Data(),"p");
    if(NumOfENoPS[i]) legcRTC->AddEntry(ratioENoPS[i],(( (TObjString*) triggersE.At(i) )->GetString()).Data(),"p");
  }
  legcRTC->Draw("same");
  
  cRelativeTriggerContent->Print(OutFileNamePDF.Data());
  cRelativeTriggerContent->Write();

  cout<<"//==========================================="<<endl;
  cout<<" // Draw ratio of tracks over CINT1B, CMUS1B, ... type versus run"<<endl;
  cout<<"//==========================================="<<endl;
  
  TH1 *hTrackerPerB[10], *hTriggerPerB[10], *hMatchedPerB[10], *hAllTracksPerB[10];
  
  TH1* hSumTriggerOverB[10], *hSumTrackerOverB[10]; 
  
  TH1* hTrackerOverTriggerB[10], *hMatchedOverTriggerB[10], *hMatchedOverTrackerB[10];	
  
  for(Int_t i = 0; i <= triggers.GetEntriesFast(); i++){
    
    hName = "hTrackerPer";
    hName += ( (TObjString*) triggersB.At(i) )->GetString();
    hTrackerPerB[i] = static_cast<TH1*>(hTrackerB[i]->Clone(hName));
    hTrackerPerB[i]->Divide(hB[i]);
    hTrackerPerB[i]->SetLineWidth(2);
    hTrackerPerB[i]->SetLineColor(kRed);
    
    hName = "hTriggerPer";
    hName += ( (TObjString*) triggersB.At(i) )->GetString();
    hTriggerPerB[i] = static_cast<TH1*>(hTriggerB[i]->Clone(hName));
    hTriggerPerB[i]->Divide(hB[i]);
    hTriggerPerB[i]->SetLineWidth(2);
    hTriggerPerB[i]->SetLineColor(kBlue);
    
    hName = "hMatchedPer";
    hName += ( (TObjString*) triggersB.At(i) )->GetString();
    hMatchedPerB[i] = static_cast<TH1*>(hMatchedB[i]->Clone(hName));
    hMatchedPerB[i]->Divide(hB[i]);
    hMatchedPerB[i]->SetLineWidth(2);
    hMatchedPerB[i]->SetLineColor(kViolet);
    
    hName = "hAllTracksPer";
    hName += ( (TObjString*) triggersB.At(i) )->GetString();
    hAllTracksPerB[i] = static_cast<TH1*>(hAllTracksB[i]->Clone(hName));
    hAllTracksPerB[i]->Divide(hB[i]);
    hAllTracksPerB[i]->SetLineWidth(3);
    hAllTracksPerB[i]->SetLineColor(kBlack);
    hAllTracksPerB[i]->SetTitle(Form("Ratio (Number of Tracks)/%s",(( (TObjString*) triggersB.At(i) )->GetString()).Data()));
    hAllTracksPerB[i]->SetMinimum(0.0001);
    hAllTracksPerB[i]->SetLabelSize(0.02);
    
    hName = Form("hSumTriggerOver%s",(( (TObjString*) triggersB.At(i) )->GetString()).Data());
    hSumTriggerOverB[i] = static_cast<TH1*>(hTriggerB[i]->Clone(hName));
    hSumTriggerOverB[i]->Add(hMatchedB[i]);
    hSumTriggerOverB[i]->Divide(hB[i]);
    
    hTriggerName = ( (TObjString*) triggersB.At(i) )->GetString();
    hName = Form("Sum of trigger tracks (matched+trigger-only) in %s events / # %s events",hTriggerName.Data(),hTriggerName.Data());
    hSumTriggerOverB[i]->SetTitle(hName);
    hSumTriggerOverB[i]->SetLabelSize(0.02);
    hSumTriggerOverB[i]->SetLineWidth(2);
    hSumTriggerOverB[i]->SetLineColor(kBlue);
    hName = Form("hSumTrackerOver%s",(( (TObjString*) triggersB.At(i) )->GetString()).Data());
    hSumTrackerOverB[i] = static_cast<TH1*>(hTrackerB[i]->Clone(hName));
    hSumTrackerOverB[i]->Add(hMatchedB[i]);
    hSumTrackerOverB[i]->Divide(hB[i]);
    hName = Form("Sum tracker tracks (matched+tracker-only) in %s events / # %s events",hTriggerName.Data(),hTriggerName.Data());
    hSumTrackerOverB[i]->SetTitle(hName);
    //hSumTrackerOverCINT1B->LabelsOption("u");
    hSumTrackerOverB[i]->SetLabelSize(0.02);
    hSumTrackerOverB[i]->SetLineWidth(2);
    hSumTrackerOverB[i]->SetLineColor(kBlue);
    
    hName = Form("hTrackerOverTrigger%s",(( (TObjString*) triggersB.At(i) )->GetString()).Data());
    hTrackerOverTriggerB[i] = static_cast<TH1*>(hTrackerB[i]->Clone(hName));
    hTrackerOverTriggerB[i]->Divide(hTriggerB[i]);
    hTriggerName = ( (TObjString*) triggersB.At(i) )->GetString();
    hName = Form("# tracker tracks / # trigger tracks in %s",hTriggerName.Data());
    hTrackerOverTriggerB[i]->SetTitle(hName);
    //hTrackerOverTriggerCINT1B->LabelsOption("u");
    hTrackerOverTriggerB[i]->SetLabelSize(0.02);
    hTrackerOverTriggerB[i]->SetLineWidth(2);
    hTrackerOverTriggerB[i]->SetLineColor(kBlue);
    
    hName = Form("hMatchedOverTrigger%s",(( (TObjString*) triggersB.At(i) )->GetString()).Data());	
    hMatchedOverTriggerB[i] = static_cast<TH1*>(hMatchedB[i]->Clone(hName));
    hMatchedOverTriggerB[i]->Divide(hTriggerB[i]);
    hTriggerName = ( (TObjString*) triggersB.At(i) )->GetString();
    hName = Form("# matched tracks / # trigger tracks in %s",hTriggerName.Data());
    hMatchedOverTriggerB[i]->SetTitle(hName);
    //hMatchedOverTriggerCINT1B->LabelsOption("u");
    hMatchedOverTriggerB[i]->SetLabelSize(0.02);
    hMatchedOverTriggerB[i]->SetLineWidth(2);
    hMatchedOverTriggerB[i]->SetLineColor(kBlue);
    
    hName = Form("hMatchedOverTracker%s",(( (TObjString*) triggersB.At(i) )->GetString()).Data());
    hMatchedOverTrackerB[i] = static_cast<TH1*>(hMatchedB[i]->Clone(hName));
    hMatchedOverTrackerB[i]->Divide(hTrackerB[i]);
    hTriggerName = ( (TObjString*) triggersB.At(i) )->GetString();
    hName = Form("# matched tracks / # tracker tracks in %s",hTriggerName.Data());
    hMatchedOverTrackerB[i]->SetTitle(hName);
    //hMatchedOverTrackerCINT1B->LabelsOption("u");
    hMatchedOverTrackerB[i]->SetLabelSize(0.02);
    hMatchedOverTrackerB[i]->SetLineWidth(2);
    hMatchedOverTrackerB[i]->SetLineColor(kBlue);
  }
  
  TCanvas *cRatioTrackTypesB[10];
  TLegend* legcTTCINT1B[10]; 
  
  TCanvas* cTrackMultB[10];
  TCanvas* cRatioTrackB[10];
  
  for(Int_t i = 1; i <= triggers.GetEntriesFast(); i++){
    CanvasName =  LHCPeriod ; 
    CanvasName += Form("_RatioTrackTypes%s",(( (TObjString*) triggersB.At(i) )->GetString()).Data()); 
    cout<<CanvasName<<endl;
    
    cRatioTrackTypesB[i] = new TCanvas(CanvasName.Data(),Form("cTrackTypes%s",(( (TObjString*) triggersB.At(i) )->GetString()).Data()),1200,900);
    cRatioTrackTypesB[i]->SetRightMargin(0.01);
    cRatioTrackTypesB[i]->SetGridy(1);
    cRatioTrackTypesB[i]->cd();	
    
    hAllTracksPerB[i]->Draw("E");
    hTrackerPerB[i]->Draw("Esame");
    hMatchedPerB[i]->Draw("Esame");
    hTriggerPerB[i]->Draw("Esame");
    
    legcTTCINT1B[i] = new TLegend(0.70,0.5,0.90,0.70);
    legcTTCINT1B[i]->AddEntry(hAllTracksPerB[i],"All tracks","l");
    legcTTCINT1B[i]->AddEntry(hTrackerPerB[i],"Tracking (only) tracks","l");
    legcTTCINT1B[i]->AddEntry(hMatchedPerB[i],"Matched tracks","l");
    legcTTCINT1B[i]->AddEntry(hTriggerPerB[i],"Trigger (only) tracks","l");
    legcTTCINT1B[i]->Draw("same");
    
    cRatioTrackTypesB[i]->Print(OutFileNamePDF.Data());
    cRatioTrackTypesB[i]->Write();
    
    if(i!=1&&i!=2) continue;
    
    CanvasName = Form("cTrackMult%s",(( (TObjString*) triggersB.At(i) )->GetString()).Data());
    cout<<CanvasName<<endl;
    cTrackMultB[i] = new TCanvas(CanvasName,CanvasName,1200,900);
    cTrackMultB[i]->Divide(1,2);
    cTrackMultB[i]->cd(1);
    hSumTriggerOverB[i]->Draw("e");
    cTrackMultB[i]->cd(2);
    hSumTrackerOverB[i]->Draw("e");
    
    cTrackMultB[i]->Print(OutFileNamePDF.Data());
    cTrackMultB[i]->Write();
    
    CanvasName = Form("cRatioTrack%s",(( (TObjString*) triggersB.At(i) )->GetString()).Data());
    cout<<CanvasName<<endl;	
    cRatioTrackB[i] = new TCanvas(CanvasName,CanvasName,1200,900);
    cRatioTrackB[i]->Divide(1,3);
    cRatioTrackB[i]->cd(1);
    hTrackerOverTriggerB[i]->Draw("e");	
    cRatioTrackB[i]->cd(2);
    hMatchedOverTriggerB[i]->Draw("e");	
    cRatioTrackB[i]->cd(3);
    hMatchedOverTrackerB[i]->Draw("e");	
    
    cRatioTrackB[i]->Print(OutFileNamePDF.Data());
    cRatioTrackB[i]->Write();
  }
  
  cout<<"  //=================================================="<<endl;
  cout<<"// Draw matched tracks asymmetry for cmus type  "<<endl;
  cout<<"  //=================================================="<<endl;
  
  TH1 *hDiffMatchedCMUS1B= static_cast<TH1*>(hPosMatchedB[kCMUS]->Clone("hDiffMatchedCMUS1B"));
  hDiffMatchedCMUS1B->Add(hNegMatchedB[kCMUS],-1);
  hDiffMatchedCMUS1B->Sumw2();
  
  TH1 *hAsymMatchedCMUS1B= static_cast<TH1*>(hDiffMatchedCMUS1B->Clone("hAsymMatchedCMUS1B"));
  hAsymMatchedCMUS1B->Divide(hAllMatchedB[kCMUS]);
  hAsymMatchedCMUS1B->SetLineColor(kRed);
  hAsymMatchedCMUS1B->SetLineWidth(2);
  hAsymMatchedCMUS1B->SetMinimum(-0.1);
  hAsymMatchedCMUS1B->SetMaximum(0.1);
  hAsymMatchedCMUS1B->SetLabelSize(0.02);
  hName = Form("Matched tracks asymmetry for %s with acc. cuts",(( (TObjString*) triggersB.At(kCMUS) )->GetString()).Data());
  hAsymMatchedCMUS1B->SetTitle(hName);
  
  CanvasName =  LHCPeriod ; 
  CanvasName += "_AsymMatched"; 
  TCanvas *cAsymMatched = new TCanvas(CanvasName.Data(),"cAsymMatched",1200,900);
  cAsymMatched->SetRightMargin(0.01);
  cAsymMatched->SetGridy(1);
  cAsymMatched->cd();
  hAsymMatchedCMUS1B->GetYaxis()->SetTitle("Asymmetry");  
  hAsymMatchedCMUS1B->Draw("EH");
  
  cAsymMatched->Print(OutFileNamePDF.Data());
  cAsymMatched->Write();
  
  
  cout<<"//=================================================="<<endl;
  cout<<"// Draw high pt tracks per trigger"<<endl;
  cout<<"  //=================================================="<<endl;
  
  TH1* hMatchedLowPtPerB[10], *hMatchedHighPtPerB[10];
  TH1* hMatchedLowPtPerBNoPS[10], *hMatchedHighPtPerBNoPS[10];
  
  for(Int_t i = 1; i <= triggers.GetEntriesFast(); i++){
    
    hName = Form("hMatchedLowPtPer%s ",(( (TObjString*) triggersB.At(i) )->GetString()).Data());
    hMatchedLowPtPerB[i] = static_cast<TH1*> (hMatchedLowPtB[i]->Clone(hName));
    hMatchedLowPtPerB[i]->Sumw2();
    hMatchedLowPtPerB[i]->Divide(hB[i]);
    hMatchedLowPtPerB[i]->SetLineWidth(2);
    hMatchedLowPtPerB[i]->SetLineColor(kBlue);
    hMatchedLowPtPerB[i]->SetTitle("");
    hName = Form("Ratio per %s ",(( (TObjString*) triggersB.At(i) )->GetString()).Data());
    hMatchedLowPtPerB[i]->GetYaxis()->SetTitle(hName);
    //hMatchedLowPtPerB[i]->SetMaximum(0.15);
    hMatchedLowPtPerB[i]->SetMinimum(0.0001);
    hMatchedLowPtPerB[i]->SetLabelSize(0.02);
    
    hName = Form("hMatchedHighPtPer%s ",(( (TObjString*) triggersB.At(i) )->GetString()).Data());
    hMatchedHighPtPerB[i] = static_cast<TH1*> (hMatchedHighPtB[i]->Clone(hName));
    hMatchedHighPtPerB[i]->Sumw2();
    hMatchedHighPtPerB[i]->Divide(hB[i]);
    hMatchedHighPtPerB[i]->SetLineWidth(2);
    hMatchedHighPtPerB[i]->SetLineColor(kRed);
    
    hName = Form("hMatchedLowPtNoPSPer%s ",(( (TObjString*) triggersB.At(i) )->GetString()).Data());
    hMatchedLowPtPerBNoPS[i] = static_cast<TH1*> (hMatchedLowPtBNoPS[i]->Clone(hName));
    hMatchedLowPtPerBNoPS[i]->Sumw2();
    hMatchedLowPtPerBNoPS[i]->Divide(hB[i]);
    hMatchedLowPtPerBNoPS[i]->SetLineStyle(3);
    hMatchedLowPtPerBNoPS[i]->SetLineColor(kBlue);
    
    hName = Form("hMatchedHighPtNoPSPer%s ",(( (TObjString*) triggersB.At(i) )->GetString()).Data());
    hMatchedHighPtPerBNoPS[i] = static_cast<TH1*> (hMatchedHighPtBNoPS[i]->Clone(hName));
    hMatchedHighPtPerBNoPS[i]->Sumw2();
    hMatchedHighPtPerBNoPS[i]->Divide(hB[i]);
    hMatchedHighPtPerBNoPS[i]->SetLineStyle(3);
    hMatchedHighPtPerBNoPS[i]->SetLineColor(kRed);
    
  }
  
  TCanvas *cHighPtMuons[10];
  TLegend* legcHPM[10];
  
  for(Int_t i = 1; i <= triggers.GetEntriesFast(); i++){
    CanvasName =  LHCPeriod ; 
    CanvasName += "_HighPtMuons"; 
    CanvasName+=i+1;
    cHighPtMuons[i]	 = new TCanvas(CanvasName.Data(),CanvasName.Data(),1200,900);
    cHighPtMuons[i]->SetTopMargin(0.05);
    cHighPtMuons[i]->SetRightMargin(0.01);
    cHighPtMuons[i]->SetGridy(1);
    cHighPtMuons[i]->cd();
    
    hMatchedLowPtPerB[i]->Draw("E");
    hMatchedHighPtPerB[i]->Draw("Esame");
    
    legcHPM[i] = new TLegend(0.60,0.45,0.98,0.65);
    hName = Form("Number of matched track per %s (include Vtx and R_{Abs} cuts)",(( (TObjString*) triggersB.At(i) )->GetString()).Data());
    legcHPM[i]->SetHeader(hName);
    legcHPM[i]->AddEntry(".","Physics selection applied :","");	
    legcHPM[i]->AddEntry(hMatchedLowPtPerB[i]," p_{T} > 1 GeV/c ","l");
    legcHPM[i]->AddEntry(hMatchedHighPtPerB[i]," p_{T} >  2 GeV/c ","l");
    legcHPM[i]->Draw("same");
    
    cHighPtMuons[i]->Print(OutFileNamePDF.Data());
    cHighPtMuons[i]->Write(); 
  }
  
  
  // close merged file	
  globalFile->Close();
  
  
  // close the PDF file
  //	c1->Print(OutFileNamePDF_close.Data());
  //		rootFileOut->Close();
  
  //===================================================================================
  //Print out the number of trigger without and with Phys. Sel.
  //===================================================================================
  
  cout << endl << endl;
  //====================================================
  if (PRINTSTAT){
    // set the format to print labels
    THashList* labels = hBWithPS[kCMUS]->GetXaxis()->GetLabels();
    TString format(Form("\n%%%ds %%9d",0));
    Int_t nRuns=0;
    
    // print value for each label
    TObjString* label = 0x0;
    TIter nextLabel(labels);
    cout << "-------------------------------------------------" << endl;
    cout << "Run Number" << "\t Number of "<< ( (TObjString*) triggersB.At(kCMUS) )->GetString()<<" after Phys. Sel. " << endl ;  
    while ((label = static_cast<TObjString*>(nextLabel()))) {
      nRuns++;
      Int_t bin = (Int_t) label->GetUniqueID();
      printf(format.Data(), label->String().Data(), (Int_t) hBWithPS[kCMUS]->GetBinContent(bin));
    }
    printf("\n========== Total #runs = %d ==============\n",nRuns);
    printf("\n\n");
    
    
    cout << "-------------------------------------------------" << endl;
    cout << "Total statistic" << endl; 
    cout << " " << endl ; 
    
    for(Int_t i = 0; i <= triggers.GetEntriesFast(); i++){
      cout << "Number of "<< ( (TObjString*) triggersB.At(i) )->GetString() << endl ;
      cout << "\t before selection " << NumOfBNoPS[i] << "\t  after selection " <<   NumOfBWithPS[i] << " --> rejection = ";
      if(NumOfBNoPS[i]>0) cout <<  (Double_t) (NumOfBNoPS[i]-NumOfBWithPS[i])/(NumOfBNoPS[i])*100. << "%" << endl ;
      else cout << "-1" << endl; 
      
      if(NumOfACNoPS[i]>0){
	cout << "Number of "<< ( (TObjString*) triggersAC.At(i) )->GetString() << endl ;
	cout << "\t before selection " << NumOfACNoPS[i] << "\t  after selection " <<   NumOfACWithPS[i] << " --> rejection = ";
	cout <<  (Double_t) (NumOfACNoPS[i]-NumOfACWithPS[i])/(NumOfACNoPS[i])*100. << "%" << endl ;
      }
      //else cout<<NumOfACNoPS[i]<<endl;
      if(NumOfENoPS[i]){	
	cout << "Number of "<< ( (TObjString*) triggersE.At(i) )->GetString() << endl ;
	cout << "\t before selection " << NumOfENoPS[i] << "\t  after selection " <<   NumOfEWithPS[i] << " --> rejection = ";
	if(NumOfENoPS[i]>0) cout <<  (Double_t) (NumOfENoPS[i]-NumOfEWithPS[i])/(NumOfENoPS[i])*100. << "%" << endl ;
      }
      //else cout<<NumOfENoPS[i]<<endl;
      cout << " " << endl ; 
    }
    
  }
  
  
  //return;
  
  //--------------------------------------------- //
  //        monitor quantities run per run        //
  //--------------------------------------------- //
  
  TH1F* hTriggerCutVsRun[2], *hTriggerCutWidthVsRun[2];
  for ( Int_t ihisto=0; ihisto<2; ++ihisto ) {
    TString cutName = ( ihisto == 0 ) ? "Lpt" : "Hpt";
    hTriggerCutVsRun[ihisto] = new TH1F(Form("hTriggerCutVsRun%s", cutName.Data()), Form("Trigger %s cut per run", cutName.Data()), 10000,1,10000);
    hTriggerCutWidthVsRun[ihisto] = (TH1F*)hTriggerCutVsRun[ihisto]->Clone(Form("hTriggerCutWidthVsRun%s", cutName.Data()));
  }
  TF1* fitMatchTrig = new TF1("fitMatchTrig","[3] + [0] * ( 1. + TMath::Erf((x - [1]) / [2]))", 0.1, 6.);
  TH1F* hNClustersPerTrackVsRun_Mean = new TH1F("hNClustersPerTrackVsRun_Mean", "averaged number of associated clusters per track;run;<n_{clusters}>",10000,1,10000);
  TH1F* hNClustersPerTrackVsRun_Sigma = new TH1F("hNClustersPerTrackVsRun_Sigma", "dispersion of the number of associated clusters per track;run;#sigma_{n_{clusters}}",10000,1,10000);
  TH1F* hNChamberHitPerTrack_Mean = new TH1F("hNChamberHitPerTrack_Mean", "averaged number of chambers hit per track;run;<n_{chamber hit}>",10000,1,10000);
  TH1F* hNChamberHitPerTrack_Sigma = new TH1F("hNChamberHitPerTrack_Sigma", "dispersion of the number of chambers hit per track;run;#sigma_{n_{chamber hit}}",10000,1,10000);
  TH1F* hChi2_Mean = new TH1F("hChi2_Mean", "averaged normalized #chi^{2} distribution;run;<#chi^{2} / ndf>",10000,1,10000);
  TH1F* hChi2_Sigma = new TH1F("hChi2_Sigma", "dispersion of normalized #chi^{2} distribution;run;#sigma_{#chi^{2} / ndf}",10000,1,10000);
  TH1F* hNClustersInCh[10];
  for (Int_t ich=0; ich<10; ich++){
    hNClustersInCh[ich] = new TH1F(Form("hNClustersInCh%d",ich+1), Form("averaged number of clusters in chamber %d per track;run;<n_{clusters}>",ich+1),10000,1,10000);
  }
  TH1F* hClusterHitMapXInCh[10];
  for (Int_t ich=0; ich<10; ich++){
    hClusterHitMapXInCh[ich] = new TH1F(Form("hClusterHitMapXInCh%d",ich+1), Form("averaged cluster position distribution in chamber %d;X (cm)",ich+1),10000,1,10000);
  }
  TH1F* hClusterHitMapYInCh[10];
  for (Int_t ich=0; ich<10; ich++){
    hClusterHitMapYInCh[ich] = new TH1F(Form("hClusterHitMapYInCh%d",ich+1), Form("averaged cluster position distribution in chamber %d;Y (cm)",ich+1),10000,1,10000);
  }
  Int_t ibin = 1;
  
  // Are the runs stored locally or in alien?
  Int_t isAlienFile = 0;
  if(alienBaseDir.Contains("alien:")){
    isAlienFile = 1;
    alienBaseDir.ReplaceAll("alien://","");
  }
  cout<<"============================================================"<<endl;
  cout<< "Monitoring quantities run per run: loop over "<<runs.GetEntriesFast()<<" runs."<<endl;
  cout<<"  directory = "<<alienBaseDir.Data()<<endl;
  cout<<"============================================================"<<endl;
  // Loop over runs
  for ( Int_t irun=0; irun<runs.GetEntriesFast(); irun++ ) {
    
    TString run = ((TObjString*)runs.UncheckedAt(irun))->GetString();
    TString run2 = ((TObjString*)runs2.UncheckedAt(irun))->GetString();
    // get the file (or list of files) to be analyzed
    TString command;
    TGridResult *res = 0;
    TObjString *objs = 0;	
    
    if(isAlienFile){
      command = Form("find %s/ %s/%s", alienBaseDir.Data(), run.Data(), QAFileName);
      res = gGrid->Command(command);
      if (!res) {
	Error("PlotMuonQApp","no result for the command: %s",command.Data());
	return;
      }
    }
    else{
      res = new TGridResult();	
      
      if(runList){
	objs = new TObjString(Form("%s/%s/%s", alienBaseDir.Data(), run2.Data(), QAFileName));
	res->Add(objs);
      }
      else {
	//loop over the directory to find the root files
	void *dir = gSystem->OpenDirectory(alienBaseDir.Data());
	TString sDirFilename;
	Int_t iEntry=0,	iFile=0;
	while(kTRUE){
	  iEntry++;
	  const	char* dirFilename = gSystem->GetDirEntry(dir);
	  if(!dirFilename) break;
	  sDirFilename = dirFilename;
	  if(!sDirFilename.IsDigit()) continue;
	  iFile++;
	  objs = new TObjString(Form("%s/%s/%s", alienBaseDir.Data(), sDirFilename.Data(), QAFileName));
	  res->Add(objs);
	}
      }
    }
    
    // Loop over 'find' results and get next LFN
    TIter nextmap(res);
    TMap *map = 0;
    
    //some checks
    Int_t iLoop=0, iLoopMax=200;
    while (kTRUE){
      
      // get the current file url
      if(isAlienFile){
	map=(TMap*)nextmap();
	if(!map) break;
	objs = dynamic_cast<TObjString*>(map->GetValue("turl"));
      }
      else{
	objs=(TObjString*)nextmap();
	if(!objs) break;
      } 
      //in case of infinite loop
      iLoop++;
      if(iLoop>iLoopMax) break;
      
      if (!objs || !objs->GetString().Length()) {
	Error("PlotMuonQApp","turl/obj not found for the run %s... SKIPPING", run.Data());
	continue;
      }
      
      // open the outfile for this run
      TFile *runFile = TFile::Open(objs->GetString());
      if (!runFile || ! runFile->IsOpen()) {
	Error("PlotMuonQApp","failed to open file: %s", objs->GetName());
	continue;//return;
      }
      runFile->Cd("MUON_QA");
      
      // get interesting histos
      TObjArray* general1 = static_cast<TObjArray*>(runFile->FindObjectAny("general1"));
      TObjArray* general2 = static_cast<TObjArray*>(runFile->FindObjectAny("general2"));
      TObjArray* expert = static_cast<TObjArray*>(runFile->FindObjectAny("expert"));
      
      if (!general1 || !general2 || !expert){
	Error("PlotMuonQApp","All objects not here !!! ===> Skipping...for %s",objs->GetName());		
	continue;
      }
      
      TH1* hNClustersPerTrack = static_cast<TH1*>(general1->FindObject("hNClustersPerTrack"));
      TH1* hNChamberHitPerTrack = static_cast<TH1*>(general1->FindObject("hNChamberHitPerTrack"));
      TH1* hChi2 = static_cast<TH1*>(general1->FindObject("hChi2"));
      TH1* hNClustersPerCh = static_cast<TH1*>(general2->FindObject("hNClustersPerCh"));
      TH1* hPtDistrib = static_cast<TH1*>(general1->FindObject("hPt"));
      TH1* hPtDistribLpt = static_cast<TH1*>(general1->FindObject("hPtMatchLpt"));
      TH1* hPtDistribHpt = static_cast<TH1*>(general1->FindObject("hPtMatchHpt"));
      if ( hPtDistrib && hPtDistribLpt && hPtDistribHpt ) {
        if ( hPtDistrib->GetSumw2N() == 0 ) hPtDistrib->Sumw2();
        TH1* histoMatch[2] = {hPtDistribLpt, hPtDistribHpt};
        for ( Int_t ihisto=0; ihisto<2; ++ihisto ) {
          if ( histoMatch[ihisto]->GetSumw2N() == 0 ) histoMatch[ihisto]->Sumw2();
          if ( histoMatch[ihisto]->GetEntries() == 0 ) continue;
          histoMatch[ihisto]->Divide(hPtDistrib);
          Double_t minEff = 99999., maxEff = -1.;
          Double_t ptMinFit = 0.1;
          Double_t ptMaxFit = 6.;
          Int_t ptBinLow = histoMatch[ihisto]->GetXaxis()->FindBin(ptMinFit);
          Int_t ptBinHigh = histoMatch[ihisto]->GetXaxis()->FindBin(ptMaxFit);
          for ( Int_t currBin=ptBinLow; currBin<=ptBinHigh; currBin++ ) {
            Double_t currEff = histoMatch[ihisto]->GetBinContent(currBin);
            Double_t currPt = histoMatch[ihisto]->GetXaxis()->GetBinCenter(currBin);
            if ( currPt < 1.5 && minEff > currEff ) {
              ptMinFit = currPt;
              minEff = currEff;
            }
            if ( currPt > 0.5 && maxEff < currEff ) {
              ptMaxFit = currPt;
              maxEff = currEff;
            }
          } // loop on histo bins
          fitMatchTrig->SetParameters(0.5, 0.5, 0.8, 0.2);
          fitMatchTrig->SetParLimits(0,0.,1.);
          fitMatchTrig->SetParLimits(1,0.,5.);
          fitMatchTrig->SetParLimits(2,0.,5.);
          fitMatchTrig->SetParLimits(3,0.,0.5);
          histoMatch[ihisto]->Fit(fitMatchTrig,"RQ0","",ptMinFit,ptMaxFit);          
          Double_t ptCut = fitMatchTrig->GetParameter(1);
          Double_t ptCutErr = fitMatchTrig->GetParError(1);
          Double_t ptCutWidth = fitMatchTrig->GetParameter(2);
          if ( ptCut < 0 || ptCut > 10. ) {
            ptCut = ptCutErr = ptCutWidth = 0.;
          }
          hTriggerCutVsRun[ihisto]->SetBinContent(ibin, ptCut);
          hTriggerCutVsRun[ihisto]->SetBinError(ibin, ptCutErr);
          hTriggerCutWidthVsRun[ihisto]->SetBinContent(ibin, ptCut);
          hTriggerCutWidthVsRun[ihisto]->SetBinError(ibin, ptCutWidth);
        } // loop on match histos
      }
      TH2* hClusterHitMapInCh[10];
      for(Int_t ich=0; ich<10; ich++) hClusterHitMapInCh[ich] = static_cast<TH2*>(expert->FindObject(Form("hClusterHitMapInCh%d",ich+1)));
      
      // skip empty runs... not anymore ! cs !
      if (!hNClustersPerCh) {
	Warning("PlotMuonQApp","File: %s has empty histograms !", objs->GetName());
	hNClustersPerTrackVsRun_Mean->SetBinContent(ibin, 0.);
	hNClustersPerTrackVsRun_Mean->SetBinError(ibin, 1.);
	hNClustersPerTrackVsRun_Sigma->SetBinContent(ibin, 0.);
	hNClustersPerTrackVsRun_Sigma->SetBinError(ibin, 1.);
	hNChamberHitPerTrack_Mean->SetBinContent(ibin, 0.);
	hNChamberHitPerTrack_Mean->SetBinError(ibin, 1.);
	hNChamberHitPerTrack_Sigma->SetBinContent(ibin, 0.);
	hNChamberHitPerTrack_Sigma->SetBinError(ibin, 1.);
	hChi2_Mean->SetBinContent(ibin, 0.);
	hChi2_Mean->SetBinError(ibin, 1.);
	hChi2_Sigma->SetBinContent(ibin, 0.);
	hChi2_Sigma->SetBinError(ibin, 1.);
	for (Int_t ich=0; ich<10; ich++) {
	  hNClustersInCh[ich]->SetBinContent(ibin,0.);
	  hNClustersInCh[ich]->SetBinError(ibin,1.);
	  hClusterHitMapXInCh[ich]->SetBinContent(ibin,0.);
	  hClusterHitMapXInCh[ich]->SetBinError(ibin,1.);
	  hClusterHitMapYInCh[ich]->SetBinContent(ibin,0.);
	  hClusterHitMapYInCh[ich]->SetBinError(ibin,1.);	
	}
	//runFile->Close();
	//continue;
      }
      else {
	// fill monitoring plots
	hNClustersPerTrackVsRun_Mean->SetBinContent(ibin, hNClustersPerTrack->GetMean());
	hNClustersPerTrackVsRun_Mean->SetBinError(ibin, hNClustersPerTrack->GetMeanError());
	hNClustersPerTrackVsRun_Sigma->SetBinContent(ibin, hNClustersPerTrack->GetRMS());
	hNClustersPerTrackVsRun_Sigma->SetBinError(ibin, hNClustersPerTrack->GetRMSError());
	hNChamberHitPerTrack_Mean->SetBinContent(ibin, hNChamberHitPerTrack->GetMean());
	hNChamberHitPerTrack_Mean->SetBinError(ibin, hNChamberHitPerTrack->GetMeanError());
	hNChamberHitPerTrack_Sigma->SetBinContent(ibin, hNChamberHitPerTrack->GetRMS());
	hNChamberHitPerTrack_Sigma->SetBinError(ibin, hNChamberHitPerTrack->GetRMSError());
	hChi2_Mean->SetBinContent(ibin, hChi2->GetMean());
	hChi2_Mean->SetBinError(ibin, hChi2->GetMeanError());
	hChi2_Sigma->SetBinContent(ibin, hChi2->GetRMS());
	hChi2_Sigma->SetBinError(ibin, hChi2->GetRMSError());
	for (Int_t ich=0; ich<10; ich++) {
	  hNClustersInCh[ich]->SetBinContent(ibin,hNClustersPerCh->GetBinContent(ich+1));
	  hNClustersInCh[ich]->SetBinError(ibin,hNClustersPerCh->GetBinError(ich+1));
	  hClusterHitMapXInCh[ich]->SetBinContent(ibin,hClusterHitMapInCh[ich]->GetMean(1));
	  hClusterHitMapXInCh[ich]->SetBinError(ibin,hClusterHitMapInCh[ich]->GetMeanError(1));
	  hClusterHitMapYInCh[ich]->SetBinContent(ibin,hClusterHitMapInCh[ich]->GetMean(2));
	  hClusterHitMapYInCh[ich]->SetBinError(ibin,hClusterHitMapInCh[ich]->GetMeanError(2));	  
	}
      }
      
      // set labels
      run = objs->GetString();
      run.ReplaceAll(Form("/%s",QAFileName), "");
      run.ReplaceAll(alienBaseDir, "");
      run.Remove(TString::kLeading, '/');
      run.Remove(TString::kLeading, '0');
      hNClustersPerTrackVsRun_Mean->GetXaxis()->SetBinLabel(ibin, run.Data());
      hNClustersPerTrackVsRun_Sigma->GetXaxis()->SetBinLabel(ibin, run.Data());
      hNChamberHitPerTrack_Mean->GetXaxis()->SetBinLabel(ibin, run.Data());
      hNChamberHitPerTrack_Sigma->GetXaxis()->SetBinLabel(ibin, run.Data());
      hChi2_Mean->GetXaxis()->SetBinLabel(ibin, run.Data());
      hChi2_Sigma->GetXaxis()->SetBinLabel(ibin, run.Data());
      for (Int_t ich=0; ich<10; ich++){
	hNClustersInCh[ich]->GetXaxis()->SetBinLabel(ibin, run.Data());
	hClusterHitMapXInCh[ich]->GetXaxis()->SetBinLabel(ibin, run.Data());
	hClusterHitMapYInCh[ich]->GetXaxis()->SetBinLabel(ibin, run.Data());
      }
      for ( Int_t ihisto=0; ihisto<2; ++ihisto) {
        hTriggerCutVsRun[ihisto]->GetXaxis()->SetBinLabel(ibin, run.Data());
        hTriggerCutWidthVsRun[ihisto]->GetXaxis()->SetBinLabel(ibin, run.Data());
      }
      // close outfile for this run
      runFile->Close();
      ibin++;      
    }
    
    delete res;
  }
  
  //sort label
  hNClustersPerTrackVsRun_Mean->LabelsOption("a");
  hNClustersPerTrackVsRun_Sigma->LabelsOption("a");
  hNChamberHitPerTrack_Mean->LabelsOption("a");
  hNChamberHitPerTrack_Sigma->LabelsOption("a");
  hChi2_Mean->LabelsOption("a");
  hChi2_Sigma->LabelsOption("a");
  
  for(Int_t ich=0; ich<10; ich++){
    hNClustersInCh[ich]->LabelsOption("a");
    hClusterHitMapXInCh[ich]->LabelsOption("a");
    hClusterHitMapYInCh[ich]->LabelsOption("a");
  }  

  TString dirToGo =  OutFileNameROOT.Data(); dirToGo+=":/";
  gDirectory->Cd(dirToGo.Data());
  //==================================================
  //Display Mean and Sigma of the number of associated clusters to a track 
  TLegend *lNClusters = new TLegend(0.75,0.85,0.99,0.99);
  lNClusters->AddEntry(hNClustersPerTrackVsRun_Mean,"clusters","PL");
  lNClusters->AddEntry(hNChamberHitPerTrack_Mean,"chamber hit","PL");
  
  TCanvas* cNClusters = new TCanvas("cNClusters","cNClusters",1200,900);
  cNClusters->Divide(1,2);
  cNClusters->cd(1);
  //hNClustersPerTrackVsRun_Mean->SetMaximum(11);
  hNClustersPerTrackVsRun_Mean->SetMinimum(7);
  hNClustersPerTrackVsRun_Mean->SetStats(kFALSE);
  hNClustersPerTrackVsRun_Mean->GetXaxis()->SetRange(1,ibin-1);
  hNClustersPerTrackVsRun_Mean->GetXaxis()->SetNdivisions(1,kFALSE);
  //hNClustersPerTrackVsRun_Mean->LabelsOption("u");
  hNClustersPerTrackVsRun_Mean->SetLabelSize(0.04);
  hNClustersPerTrackVsRun_Mean->SetTitle("averaged number of associated clusters or of the number of chamber hit per track");
  hNClustersPerTrackVsRun_Mean->SetLineWidth(2);
  hNClustersPerTrackVsRun_Mean->Draw("e");
  hNChamberHitPerTrack_Mean->SetLineColor(kRed);
  hNChamberHitPerTrack_Mean->SetLineWidth(2);
  hNChamberHitPerTrack_Mean->Draw("esame");
  lNClusters->Draw("same");
  
  cNClusters->cd(2);
  //hNClustersPerTrackVsRun_Sigma->SetMaximum(1.1);
  hNClustersPerTrackVsRun_Sigma->SetMinimum(0.4);
  hNClustersPerTrackVsRun_Sigma->SetStats(kFALSE);
  hNClustersPerTrackVsRun_Sigma->GetXaxis()->SetRange(1,ibin-1);
  hNClustersPerTrackVsRun_Sigma->GetXaxis()->SetNdivisions(1,kFALSE);
  //hNClustersPerTrackVsRun_Sigma->LabelsOption("u");
  hNClustersPerTrackVsRun_Sigma->SetLabelSize(0.04);
  hNClustersPerTrackVsRun_Sigma->SetTitle("dispersion of the number of associated clusters or of the number of chamber hit per track");
  hNClustersPerTrackVsRun_Sigma->SetLineWidth(2);
  hNClustersPerTrackVsRun_Sigma->Draw("e");
  hNChamberHitPerTrack_Sigma->SetLineWidth(2);
  hNChamberHitPerTrack_Sigma->SetLineColor(kRed);
  hNChamberHitPerTrack_Sigma->Draw("esame");
  lNClusters->Draw("same");

  cNClusters->Print(OutFileNamePDF.Data());
  cNClusters->Write();


  //==================================================
  // Display average number of cluster per chamber
  TLegend *lNClustersPerCh = new TLegend(0.92,0.45,0.99,0.99);
  TCanvas* cNClustersPerCh = new TCanvas("cNClustersPerCh","cNClustersPerCh",1200,900);
  cNClustersPerCh->cd();
  cNClustersPerCh->SetRightMargin(0.1);
  hNClustersInCh[0]->SetStats(kFALSE);
  hNClustersInCh[0]->GetXaxis()->SetRange(1,ibin-1);
  hNClustersInCh[0]->GetXaxis()->SetNdivisions(1,kFALSE);
  //hNClustersInCh[0]->LabelsOption("u");
  hNClustersInCh[0]->SetLabelSize(0.02);
  hNClustersInCh[0]->SetTitle("averaged number of clusters in chamber i per track");
  hNClustersInCh[0]->SetMaximum(1.2);
  hNClustersInCh[0]->SetMinimum(0.01);
  for (Int_t ich=0; ich<10; ich++) {
    hNClustersInCh[ich]->SetLineColor(ich+1+ich/9);
    hNClustersInCh[ich]->SetLineWidth(2);
    if (ich == 0) hNClustersInCh[ich]->Draw("e");
    else hNClustersInCh[ich]->Draw("esame");
    lNClustersPerCh->AddEntry(hNClustersInCh[ich],Form("ch%d",ich+1),"PL");
  }
  lNClustersPerCh->Draw("same");
  
  cNClustersPerCh->Print(OutFileNamePDF.Data());
  cNClustersPerCh->Write();

  //==================================================
  // Display average X and Y position of clusters per chamber
  TLegend *lClusterHitMapPerCh = new TLegend(0.92,0.45,0.99,0.99);
  TCanvas* cClusterHitMapPerCh = new TCanvas("cClusterHitMapPerCh","cClusterHitMapPerCh",1200,900);
  cClusterHitMapPerCh->Divide(1,2);
  cClusterHitMapPerCh->GetPad(1)->SetRightMargin(0.1);
  cClusterHitMapPerCh->GetPad(2)->SetRightMargin(0.1);
	
  cClusterHitMapPerCh->cd(1);
  hClusterHitMapXInCh[0]->SetStats(kFALSE);
  hClusterHitMapXInCh[0]->GetXaxis()->SetRange(1,ibin-1);
  hClusterHitMapXInCh[0]->GetXaxis()->SetNdivisions(1,kFALSE);
  //hNClustersInCh[0]->LabelsOption("u");
  hClusterHitMapXInCh[0]->SetLabelSize(0.04);
  hClusterHitMapXInCh[0]->SetTitle("<X> of clusters - associated to a track - in chamber i");
  hClusterHitMapXInCh[0]->SetMaximum(30);
  hClusterHitMapXInCh[0]->SetMinimum(-30);
  for (Int_t ich=0; ich<10; ich++) {
    hClusterHitMapXInCh[ich]->SetLineColor(ich+1+ich/9);
    hClusterHitMapXInCh[ich]->SetLineWidth(2);
    if (ich == 0) hClusterHitMapXInCh[ich]->Draw("e");
    else hClusterHitMapXInCh[ich]->Draw("esame");
    
    lClusterHitMapPerCh->AddEntry(hClusterHitMapXInCh[ich],Form("ch%d",ich+1),"PL");
  }
  lClusterHitMapPerCh->Draw("same");
  
  cClusterHitMapPerCh->cd(2);
  hClusterHitMapYInCh[0]->SetStats(kFALSE);
  hClusterHitMapYInCh[0]->GetXaxis()->SetRange(1,ibin-1);
  hClusterHitMapYInCh[0]->GetXaxis()->SetNdivisions(1,kFALSE);
  //hNClustersInCh[0]->LabelsOption("u");
  hClusterHitMapYInCh[0]->SetLabelSize(0.04);
  hClusterHitMapYInCh[0]->SetTitle("<Y> of clusters - associated to a track - in chamber i");
  hClusterHitMapYInCh[0]->SetMaximum(30);
  hClusterHitMapYInCh[0]->SetMinimum(-30);
  for (Int_t ich=0; ich<10; ich++) {
    hClusterHitMapYInCh[ich]->SetLineColor(ich+1+ich/9);
    hClusterHitMapYInCh[ich]->SetLineWidth(2);
    if (ich == 0) hClusterHitMapYInCh[ich]->Draw("e");
    else hClusterHitMapYInCh[ich]->Draw("esame");
  }
  lClusterHitMapPerCh->Draw("same");

  cClusterHitMapPerCh->Print(OutFileNamePDF.Data());
  cClusterHitMapPerCh->Write();


  //==================================================
  // Display tracks ChiSquare 
  TCanvas* cChi2 = new TCanvas("cChi2","cChi2",1200,900);
  cChi2->Divide(1,2);
  cChi2->cd(1);
  hChi2_Mean->SetStats(kFALSE);
  hChi2_Mean->GetXaxis()->SetRange(1,ibin-1);
  hChi2_Mean->GetXaxis()->SetNdivisions(1,kFALSE);
  //hChi2_Mean->LabelsOption("u");
  hChi2_Mean->SetLabelSize(0.04);
  hChi2_Mean->SetLineWidth(2);
  hChi2_Mean->Draw("e");

  cChi2->cd(2);
  hChi2_Sigma->SetStats(kFALSE);
  hChi2_Sigma->GetXaxis()->SetRange(1,ibin-1);
  hChi2_Sigma->GetXaxis()->SetNdivisions(1,kFALSE);
  //hChi2_Sigma->LabelsOption("u");
  hChi2_Sigma->SetLabelSize(0.04);
  hChi2_Sigma->SetLineWidth(2);
  hChi2_Sigma->Draw("e");

  cChi2->Print(OutFileNamePDF.Data());
  cChi2->Write();
  
  //==================================================
  // Display track Lpt/Hpt 
  if ( hTriggerCutVsRun[0] && hTriggerCutVsRun[1] ) {
    TCanvas* cLptHpt = new TCanvas("cLptHpt","cLptHpt",1200,900);
    cLptHpt->Divide(1,2);
    TLegend* legLptHpt = new TLegend(0.72,0.7,0.9,0.85);
    legLptHpt->SetBorderSize(1);
    for ( Int_t ihisto=0; ihisto<2; ++ihisto) {
      cLptHpt->cd(ihisto+1);
      TH1* currHistos[2] = {hTriggerCutVsRun[ihisto], hTriggerCutWidthVsRun[ihisto]};
      for ( Int_t jhisto=0; jhisto<2; jhisto++ ) {
        currHistos[jhisto]->GetXaxis()->SetRange(1,ibin-1);
        currHistos[jhisto]->GetYaxis()->SetRangeUser(0.,5.);
        currHistos[jhisto]->LabelsOption("a");
        currHistos[jhisto]->SetStats(kFALSE);
        currHistos[jhisto]->GetXaxis()->SetLabelSize(0.04);
        currHistos[jhisto]->SetLineWidth(2);
      }
      hTriggerCutWidthVsRun[ihisto]->SetLineColor(2);
      hTriggerCutWidthVsRun[ihisto]->SetMarkerColor(2);
      hTriggerCutWidthVsRun[ihisto]->SetFillColor(2);
      hTriggerCutWidthVsRun[ihisto]->SetFillStyle(3001);
      hTriggerCutWidthVsRun[ihisto]->Draw("e2");
      hTriggerCutVsRun[ihisto]->Draw("esame");
      if ( ihisto == 0 ) {
        legLptHpt->AddEntry(hTriggerCutWidthVsRun[ihisto],"Fit width","f");
        legLptHpt->AddEntry(hTriggerCutVsRun[ihisto],"pt cut from fit (stat error)","lp");
        legLptHpt->Draw("same");
      }
    }
    cLptHpt->Print(OutFileNamePDF.Data());
    cLptHpt->Write();
  }

  // close the PDF file
  c1->Print(OutFileNamePDF_close.Data());
  rootFileOut->Close();
	    			
  return;
	
}
	
 
