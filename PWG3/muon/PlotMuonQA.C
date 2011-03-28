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

// ALIROOT includes
#include "AliCounterCollection.h"

#endif


// .x PlotMuonQA.C("alien:///alice/cern.ch/user/c/cynthia/muon/QA/pp/LHC10e/pass2_test/results",0,kFALSE)
// .x PlotMuonQA.C("/Users/cynthia/Documents/alice/data/MuonQA/results",0,kFALSE)
// .x PlotMuonQA.C("/Users/cynthia/Documents/alice/data/MuonQA/results","/Users/cynthia/Documents/alice/data/MuonQA/LHC10e/pass2/runlist_period3_test3_3runs.txt",kFALSE)
//--------------------------------------------------------------------------
void PlotMuonQA(const char* baseDir, const char* runList = 0x0, Bool_t selectPhysics = kFALSE)
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
  
  TString LHCPeriod = "LHC11a"; 
  
  
  TString OutFileName = "QA_";  OutFileName += LHCPeriod;
  TString OutFileNamePDF=  OutFileName.Data();  OutFileNamePDF+= ".pdf";
  TString OutFileNamePDF_open = OutFileNamePDF.Data(); OutFileNamePDF_open += "[";  
  TString OutFileNamePDF_close= OutFileNamePDF.Data(); OutFileNamePDF_close += "]";  
  TString OutFileNameROOT=  OutFileName.Data();  OutFileNameROOT+= ".root";
  
  Int_t PRINTSTAT = 1; 
  
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
  TObjArray runs;
  runs.SetOwner();

  if (runList) {
    // only the ones in the runList
    ifstream inFile(runList);
    if (!inFile.is_open()) {
      Error("PlotQA",Form("unable to open file %s", runList));
      return;
    }
    
    TString currRun;
    while (!inFile.eof()) {
      currRun.ReadLine(inFile, kTRUE);
      if (currRun.IsNull()) continue;
      if (!currRun.IsDigit()) {
	Error("PlotQA","invalid run number: %s", currRun.Data());
	return;
      }
      runs.AddLast(new TObjString(Form("%09d", currRun.Atoi())));
      selectRuns += Form("%s,",currRun.Data());
    }
    selectRuns.Remove(TString::kTrailing, ',');
    inFile.close();
    
  } else {
    // all runs
    runs.AddLast(new TObjString("*"));
  }
  
  // physics selection
  TString select = selectPhysics ? "selected:yes" : "";
  
  
  //---------------------------------- //
  //        plot global counter        //
  //---------------------------------- //
  
  TFile *globalFile = TFile::Open(Form("%s/MergedAnalysisResults.root", baseDir));
  if (!globalFile || ! globalFile->IsOpen()) {
    Error("PlotQA", Form("failed to open file: %s/MergedAnalysisResults.root", baseDir));
    return;
  }
  globalFile->Cd("MUON_QA");
  
  // get counters
  AliCounterCollection* eventCounters = static_cast<AliCounterCollection*>(globalFile->FindObjectAny("eventCounters"));
  AliCounterCollection* trackCounters = static_cast<AliCounterCollection*>(globalFile->FindObjectAny("trackCounters"));
  if (!runList) selectRuns += trackCounters->GetKeyWords("run");
  
  //declare a default canvas c1 
  TString CanvasName = "c1";
  TCanvas *c1 = new TCanvas(CanvasName.Data());
  c1->cd();
  
  // Histo trigger without Phys. Sel. 
  TH1* hAllTriggersNoPS = eventCounters->Draw("run",Form("trigger:CINT1B,CMUS1B,CSH1B/%s", selectRuns.Data()));
	if(!hAllTriggersNoPS) return;
  hAllTriggersNoPS->Sumw2();
  TH1* hCINT1BNoPS = eventCounters->Draw("run",Form("trigger:CINT1B/%s", selectRuns.Data()));
	if(!hCINT1BNoPS) return;
  hCINT1BNoPS->Sumw2();
  Int_t NumOfCINT1BNoPS = hCINT1BNoPS->Integral();
  TH1* hCMUS1BNoPS = eventCounters->Draw("run",Form("trigger:CMUS1B/%s", selectRuns.Data()));
	if(!hCMUS1BNoPS) return;
  hCMUS1BNoPS->Sumw2();
  Int_t NumOfCMUS1BNoPS = hCMUS1BNoPS->Integral();
  TH1* hCSH1BNoPS = eventCounters->Draw("run",Form("trigger:CSH1B/%s", selectRuns.Data()));
	Int_t NumOfCSH1BNoPS = 0;
	if(!hCMUS1BNoPS) return;
	hCSH1BNoPS->Sumw2();
	NumOfCSH1BNoPS = hCSH1BNoPS->Integral();
	

  // Histo trigger with Phys. Sel. 
  TH1* hAllTriggersWithPS = eventCounters->Draw("run",Form("trigger:CINT1B,CMUS1B,CSH1B/%s/selected:yes", selectRuns.Data()));
  hAllTriggersWithPS->Sumw2();
  TH1* hCINT1BWithPS = eventCounters->Draw("run",Form("trigger:CINT1B/%s/selected:yes", selectRuns.Data()));
  hCINT1BWithPS->Sumw2();
  Int_t NumOfCINT1BWithPS = hCINT1BWithPS->Integral();
  TH1* hCMUS1BWithPS = eventCounters->Draw("run",Form("trigger:CMUS1B/%s/selected:yes", selectRuns.Data()));
  hCMUS1BWithPS->Sumw2();
  Int_t NumOfCMUS1BWithPS = hCMUS1BWithPS->Integral();
  TH1* hCSH1BWithPS = eventCounters->Draw("run",Form("trigger:CSH1B/%s/selected:yes", selectRuns.Data()));
  hCSH1BWithPS->Sumw2();
  Int_t NumOfCSH1BWithPS = hCSH1BWithPS->Integral();

  //Background estimator in CMUS1
  TH1* hCMUS1ACNoPS = eventCounters->Draw("run",Form("trigger:CMUS1AC/%s", selectRuns.Data()));
  hCMUS1ACNoPS->Sumw2();
  TH1* hCMUS1ENoPS = eventCounters->Draw("run",Form("trigger:CMUS1E/%s", selectRuns.Data()));
	Int_t NumOfCMUS1ENoPS = 0;
	if(hCMUS1ENoPS){
		hCMUS1ENoPS->Sumw2();
		NumOfCMUS1ENoPS = hCMUS1ENoPS->Integral();
	}

  // Histo trigger : Phys. Sel.  is selected or not depending on the macro arguments
  TH1* hAllTriggers = eventCounters->Draw("run",Form("trigger:CINT1B,CMUS1B,CSH1B/%s/%s", selectRuns.Data(), select.Data()));
  hAllTriggers->Sumw2();
  TH1* hCINT1B = eventCounters->Draw("run",Form("trigger:CINT1B/%s/%s", selectRuns.Data(), select.Data()));
  hCINT1B->Sumw2();
  TH1* hCMUS1B = eventCounters->Draw("run",Form("trigger:CMUS1B/%s/%s", selectRuns.Data(), select.Data()));
  hCMUS1B->Sumw2();
  TH1* hCSH1B = eventCounters->Draw("run",Form("trigger:CSH1B/%s/%s", selectRuns.Data(), select.Data()));
  hCSH1B->Sumw2();

  // Histo tracking : Phys. Sel.  is selected or not depending on the macro arguments
  TH1* hTriggerCINT1B = trackCounters->Draw("run",Form("track:triggeronly/trigger:CINT1B/%s/%s", selectRuns.Data(), select.Data()));
  hTriggerCINT1B->Sumw2();
  TH1* hTrackerCINT1B = trackCounters->Draw("run",Form("track:trackeronly/trigger:CINT1B/%s/%s", selectRuns.Data(), select.Data()));
  hTrackerCINT1B->Sumw2();
  TH1* hMatchedCINT1B = trackCounters->Draw("run",Form("track:matched/trigger:CINT1B/%s/%s", selectRuns.Data(), select.Data()));
  hMatchedCINT1B->Sumw2();
  TH1* hTriggerCMUS1B = trackCounters->Draw("run",Form("track:triggeronly/trigger:CMUS1B/%s/%s", selectRuns.Data(), select.Data()));
  hTriggerCMUS1B->Sumw2();
  TH1* hTrackerCMUS1B = trackCounters->Draw("run",Form("track:trackeronly/trigger:CMUS1B/%s/%s", selectRuns.Data(), select.Data()));
  hTrackerCMUS1B->Sumw2();
  TH1* hMatchedCMUS1B = trackCounters->Draw("run",Form("track:matched/trigger:CMUS1B/%s/%s", selectRuns.Data(), select.Data()));
  hMatchedCMUS1B->Sumw2();
  TH1* hAllTracksCINT1B = trackCounters->Draw("run",Form("trigger:CINT1B/%s/%s", selectRuns.Data(), select.Data()));
  hAllTracksCINT1B->Sumw2();
  TH1* hAllTracksCMUS1B = trackCounters->Draw("run",Form("trigger:CMUS1B/%s/%s", selectRuns.Data(), select.Data()));
  hAllTracksCMUS1B->Sumw2();
  
  TH1* hMatchedLowPtCMUS1B = trackCounters->Draw("run",Form("track:matched/trigger:CMUS1B/%s/%s/pt:low", selectRuns.Data(), select.Data()));
  hMatchedLowPtCMUS1B->Sumw2();
  TH1* hMatchedHighPtCMUS1B = trackCounters->Draw("run",Form("track:matched/trigger:CMUS1B/%s/%s/pt:high", selectRuns.Data(), select.Data()));
  hMatchedHighPtCMUS1B->Sumw2();

  TH1 *hPosMatched =  trackCounters->Draw("run",Form("track:matched/%s/charge:pos/%s",select.Data(),selectRuns.Data()));
  hPosMatched->Sumw2();
  TH1 *hNegMatched =  trackCounters->Draw("run",Form("track:matched/%s/charge:neg/%s",select.Data(),selectRuns.Data())); 
  hNegMatched->Sumw2();
  TH1 *hAllMatched=  trackCounters->Draw("run",Form("track:matched/%s/%s",select.Data(),selectRuns.Data())); 
  hAllMatched->Sumw2();
   // for CMUS1B only 
  TH1 *hPosMatchedCMUS1B =  trackCounters->Draw("run",Form("track:matched/%s/charge:pos/trigger:CMUS1B/%s",select.Data(),selectRuns.Data()));
  hPosMatchedCMUS1B->Sumw2();
  TH1 *hNegMatchedCMUS1B=  trackCounters->Draw("run",Form("track:matched/%s/charge:neg/trigger:CMUS1B/%s",select.Data(),selectRuns.Data())); 
  hNegMatchedCMUS1B->Sumw2();
  TH1 *hAllMatchedCMUS1B=  trackCounters->Draw("run",Form("track:matched/%s/trigger:CMUS1B/%s",select.Data(),selectRuns.Data())); 
  hAllMatchedCMUS1B->Sumw2();

  //TH1* hAll = eventCounters->Draw("trigger","run","run:any");
  //hAll->Draw();
  
  //===================================================================================
  // Put all plots in a ps file, easier to publish (Twiki)
  //TPostScript QAps(OutFileName.Data(),4112); 
  //QAps.Range(26,19);
  //c1->Print(PSopen.Data());
  c1->Print(OutFileNamePDF_open.Data());

  TFile *rootFileOut = TFile::Open(OutFileNameROOT.Data(),"RECREATE");
  
  //===================================================================================
  // new canvas with the relative content of each trigger w/ and w/o  physics selection
  TH1 *ratioCMUS1B = static_cast<TH1*>(hCMUS1BWithPS->Clone("ratioCMUS1B"));
  ratioCMUS1B->Divide(hAllTriggersNoPS);
  ratioCMUS1B->SetLineColor(kRed);
  ratioCMUS1B->SetLineWidth(2);
  TH1 *ratioCINT1B = static_cast<TH1*>(hCINT1BWithPS->Clone("ratioCINT1B"));
  ratioCINT1B->Divide(hAllTriggersNoPS);
  ratioCINT1B->SetLineColor(kBlue);
  ratioCINT1B->SetLineWidth(2);
  TH1 *ratioCSH1B = static_cast<TH1*>(hCSH1BWithPS->Clone("ratioCSH1B"));
  ratioCSH1B->Divide(hAllTriggersNoPS);
  ratioCSH1B->SetLineColor(kGreen);
  ratioCSH1B->SetLineWidth(2);

  TH1 *ratioCMUS1BNoPS = static_cast<TH1*>(hCMUS1BNoPS->Clone("ratioCMUS1BNoPS"));
  ratioCMUS1BNoPS->Divide(hAllTriggersNoPS);
  ratioCMUS1BNoPS->SetLineColor(kRed);
  ratioCMUS1BNoPS->SetLineWidth(0);
  ratioCMUS1BNoPS->SetLineStyle(2);
  ratioCMUS1BNoPS->SetMarkerColor(kRed);
  ratioCMUS1BNoPS->SetMarkerStyle(24);
  ratioCMUS1BNoPS->SetMarkerSize(1);
  TH1 *ratioCINT1BNoPS = static_cast<TH1*>(hCINT1BNoPS->Clone("ratioCINT1BNoPS"));
  ratioCINT1BNoPS->Divide(hAllTriggersNoPS);
  ratioCINT1BNoPS->SetLineColor(kBlue);
  ratioCINT1BNoPS->SetLineWidth(0);
  ratioCINT1BNoPS->SetLineStyle(2);
  ratioCINT1BNoPS->SetMarkerColor(kBlue);
  ratioCINT1BNoPS->SetMarkerStyle(25);
  ratioCINT1BNoPS->SetMarkerSize(1);
  TH1 *ratioCSH1BNoPS = static_cast<TH1*>(hCSH1BNoPS->Clone("ratioCSH1BNoPS"));
  ratioCSH1BNoPS->Divide(hAllTriggersNoPS);
  ratioCSH1BNoPS->SetLineColor(kGreen);
  ratioCSH1BNoPS->SetLineWidth(0);
  ratioCSH1BNoPS->SetLineStyle(2);
  ratioCSH1BNoPS->SetMarkerColor(kGreen);
  ratioCSH1BNoPS->SetMarkerStyle(26);
  ratioCSH1BNoPS->SetMarkerSize(1);

  CanvasName =  LHCPeriod.Data() ; 
  CanvasName += "_RelativeTriggerContent"; 
  TCanvas *cRelativeTriggerContent = new TCanvas(CanvasName.Data(),"cRelativeTriggerContent",1200,900);
  cRelativeTriggerContent->SetTopMargin(0.05);
  cRelativeTriggerContent->SetRightMargin(0.01);
  cRelativeTriggerContent->SetGridy(1);
  cRelativeTriggerContent->SetLogy(1);
  cRelativeTriggerContent->cd();

  ratioCMUS1B->SetMaximum(1.5);
  ratioCMUS1B->SetMinimum(0.001);
  ratioCMUS1B->SetTitle("");
  ratioCMUS1B->SetLabelSize(0.02);
  ratioCMUS1B->GetYaxis()->SetTitle("Relative trigger content w/ and w/o Phys. Sel."); 

  ratioCMUS1B->Draw();
  ratioCINT1B->Draw("ESAME");
  ratioCSH1B->Draw("ESAME");
  ratioCINT1BNoPS->Draw("EPSAME");
  ratioCMUS1BNoPS->Draw("EPSAME");
  ratioCSH1BNoPS->Draw("EPSAME"); 

  TLegend* legcRTC = new TLegend(0.2,0.15,0.50,0.40);
  legcRTC->SetHeader("Physics Selection");
  legcRTC->AddEntry(".","applied :","");
  legcRTC->AddEntry(ratioCMUS1B,"CMUS1B","l");
  legcRTC->AddEntry(ratioCINT1B,"CINT1B","l");
  legcRTC->AddEntry(ratioCSH1B,"CSH1B","l");
  legcRTC->AddEntry(".","not applied :","");
  legcRTC->AddEntry(ratioCMUS1BNoPS,"CMUS1B","p");
  legcRTC->AddEntry(ratioCINT1BNoPS,"CINT1B","p");
  legcRTC->AddEntry(ratioCSH1BNoPS,"CSH1B","p");
  legcRTC->Draw("same");
  
  cRelativeTriggerContent->Print(OutFileNamePDF.Data());
  cRelativeTriggerContent->Write();

  //===================================================================================
  // new canvas to estimate CMUS1B background percentage w/o  physics selection
  TH1 *ratioBckCMUS1BNoPS = static_cast<TH1*>(hCMUS1BNoPS->Clone("ratioBckCMUS1BNoPS"));
  ratioBckCMUS1BNoPS->Add(hCMUS1ACNoPS,-1);
  if(hCMUS1ENoPS) ratioBckCMUS1BNoPS->Add(hCMUS1ENoPS,2);
  ratioBckCMUS1BNoPS->Divide(hCMUS1BNoPS);
  ratioBckCMUS1BNoPS->SetLineColor(kRed);
  ratioBckCMUS1BNoPS->SetLineWidth(2);

  TH1 *ratioBckCMUS1BWithPS = static_cast<TH1*>(hCMUS1BWithPS->Clone("ratioBckCMUS1BWithPS"));
  ratioBckCMUS1BWithPS->Divide(hCMUS1BNoPS);
  ratioBckCMUS1BWithPS->SetLineColor(kBlue);
  ratioBckCMUS1BWithPS->SetLineWidth(2);

  CanvasName =  LHCPeriod.Data() ; 
  CanvasName += "_BackgroundInCMUS1B"; 
  TCanvas *cBackgroundInCMUS1B = new TCanvas(CanvasName.Data(),"cBackgroundInCMUS1B",1200,900);
  cBackgroundInCMUS1B->SetTopMargin(0.05);
  cBackgroundInCMUS1B->SetRightMargin(0.01);
  cBackgroundInCMUS1B->SetGridy(1);
  //cBackgroundInCMUS1B->SetLogy(1);
  cBackgroundInCMUS1B->cd();

  ratioBckCMUS1BNoPS->SetMaximum(1.2);
  ratioBckCMUS1BNoPS->SetMinimum(0.1);
  ratioBckCMUS1BNoPS->SetTitle("");
  ratioBckCMUS1BNoPS->SetLabelSize(0.02);
  ratioBckCMUS1BNoPS->GetYaxis()->SetTitle("Ratio"); 

  ratioBckCMUS1BNoPS->Draw();
  ratioBckCMUS1BWithPS->Draw("ESAME");

  TLegend* legcBTS = new TLegend(0.2,0.25,0.50,0.40);
  legcBTS->SetHeader("Background subtraction in CMUS1B: B-A-C+2E");
  legcBTS->AddEntry(ratioBckCMUS1BNoPS,"Background substracted/All ","l");
  legcBTS->AddEntry(ratioBckCMUS1BWithPS,"Physics Selected/All","l");
  legcBTS->Draw("same");
  
  cBackgroundInCMUS1B->Print(OutFileNamePDF.Data());
  cBackgroundInCMUS1B->Write();

  //===========================================
  // Draw ratio of tracks over CINT1B versus run
  TH1* hTrackerPerCINT1B= static_cast<TH1*>(hTrackerCINT1B->Clone("hTrackerPerCINT1B"));
  hTrackerPerCINT1B->Divide(hCINT1B);
  hTrackerPerCINT1B->SetLineWidth(2);
  hTrackerPerCINT1B->SetLineColor(kRed);

  TH1* hTriggerPerCINT1B= static_cast<TH1*>(hTriggerCINT1B->Clone("hTriggerPerCINT1B"));
  hTriggerPerCINT1B->Divide(hCINT1B);
  hTriggerPerCINT1B->SetLineWidth(2);
  hTriggerPerCINT1B->SetLineColor(kBlue);

  TH1* hMatchedPerCINT1B= static_cast<TH1*>(hMatchedCINT1B->Clone("hMatchedPerCINT1B"));
  hMatchedPerCINT1B->Divide(hCINT1B);
  hMatchedPerCINT1B->SetLineWidth(2);
  hMatchedPerCINT1B->SetLineColor(kViolet);


  TH1* hAllTracksPerCINT1B= static_cast<TH1*>(hAllTracksCINT1B->Clone("hAllTracksPerCINT1B"));
  hAllTracksPerCINT1B->Divide(hCINT1B);
  hAllTracksPerCINT1B->SetLineWidth(3);
  hAllTracksPerCINT1B->SetLineColor(kBlack);

  CanvasName =  LHCPeriod.Data() ; 
  CanvasName += "_RatioTrackTypesCINT1B"; 
  TCanvas *cRatioTrackTypesCINT1B = new TCanvas(CanvasName.Data(),"cRatioTrackTypesCINT1B",1200,900);
  cRatioTrackTypesCINT1B->SetRightMargin(0.01);
  cRatioTrackTypesCINT1B->SetGridy(1);
  cRatioTrackTypesCINT1B->cd();

  hAllTracksPerCINT1B->SetTitle("Ratio (Number of Tracks)/CINT1B");
  hAllTracksPerCINT1B->SetMinimum(0.0001);
  hAllTracksPerCINT1B->SetLabelSize(0.02);
  hAllTracksPerCINT1B->Draw("E");

  hTrackerPerCINT1B->Draw("Esame");
  hMatchedPerCINT1B->Draw("Esame");
  hTriggerPerCINT1B->Draw("Esame");

  TLegend* legcTTCINT1B = new TLegend(0.70,0.5,0.90,0.70);
  legcTTCINT1B->AddEntry(hAllTracksPerCINT1B,"All tracks","l");
  legcTTCINT1B->AddEntry(hTrackerPerCINT1B,"Tracking (only) tracks","l");
  legcTTCINT1B->AddEntry(hMatchedPerCINT1B,"Matched tracks","l");
  legcTTCINT1B->AddEntry(hTriggerPerCINT1B,"Trigger (only) tracks","l");
  legcTTCINT1B->Draw("same");

  cRatioTrackTypesCINT1B->Print(OutFileNamePDF.Data());
  cRatioTrackTypesCINT1B->Write();

  //===========================================
  // draw ratio of track (summed, i.e all trigger=(matched+trigger-only)) over CINT1B versus run
  TCanvas* cTrackMultCINT1B = new TCanvas("cTrackMultCINT1B","cTrackMultCINT1B",1200,900);
  // must be applied on the pads created by  divide
  cTrackMultCINT1B->Divide(1,2);
  cTrackMultCINT1B->cd(1);
  TH1* hSumTriggerOverCINT1B = static_cast<TH1*>(hTriggerCINT1B->Clone("hSumTriggerOverCINT1B"));
  hSumTriggerOverCINT1B->Add(hMatchedCINT1B);
  hSumTriggerOverCINT1B->Divide(hCINT1B);
  hSumTriggerOverCINT1B->SetTitle("Sum of trigger tracks (matched+trigger-only) in CINT1B events / # CINT1B events");
  //hSumTriggerOverCINT1B->LabelsOption("u");
  hSumTriggerOverCINT1B->SetLabelSize(0.02);
  hSumTriggerOverCINT1B->SetLineWidth(2);
  hSumTriggerOverCINT1B->SetLineColor(kBlue);
  hSumTriggerOverCINT1B->Draw("e");
  cTrackMultCINT1B->cd(2);
  TH1* hSumTrackerOverCINT1B = static_cast<TH1*>(hTrackerCINT1B->Clone("hSumTrackerOverCINT1B"));
  hSumTrackerOverCINT1B->Add(hMatchedCINT1B);
  hSumTrackerOverCINT1B->Divide(hCINT1B);
  hSumTrackerOverCINT1B->SetTitle("Sum tracker tracks (matched+tracker-only) in CINT1B events / # CINT1B events");
  //hSumTrackerOverCINT1B->LabelsOption("u");
  hSumTrackerOverCINT1B->SetLabelSize(0.02);
  hSumTrackerOverCINT1B->SetLineWidth(2);
  hSumTrackerOverCINT1B->SetLineColor(kBlue);
  hSumTrackerOverCINT1B->Draw("e");

  cTrackMultCINT1B->Print(OutFileNamePDF.Data());
  cTrackMultCINT1B->Write();
  
  //===========================================
  // draw mixed ratio of track over track versus run for CINT1B
  TCanvas* cRatioTrackCINT1B = new TCanvas("cRatioTrackCINT1B","cRatioTrackCINT1B",1200,900);
  cRatioTrackCINT1B->Divide(1,3);

  cRatioTrackCINT1B->cd(1);
  TH1* hTrackerOverTriggerCINT1B = static_cast<TH1*>(hTrackerCINT1B->Clone("hTrackerOverTriggerCINT1B"));
  hTrackerOverTriggerCINT1B->Divide(hTriggerCINT1B);
  hTrackerOverTriggerCINT1B->SetTitle("# tracker tracks / # trigger tracks in CINT1B");
  //hTrackerOverTriggerCINT1B->LabelsOption("u");
  hTrackerOverTriggerCINT1B->SetLabelSize(0.02);
  hTrackerOverTriggerCINT1B->SetLineWidth(2);
  hTrackerOverTriggerCINT1B->SetLineColor(kBlue);
  hTrackerOverTriggerCINT1B->Draw("e");

  cRatioTrackCINT1B->cd(2);
  TH1* hMatchedOverTriggerCINT1B = static_cast<TH1*>(hMatchedCINT1B->Clone("hMatchedOverTriggerCINT1B"));
  hMatchedOverTriggerCINT1B->Divide(hTriggerCINT1B);
  hMatchedOverTriggerCINT1B->SetTitle("# matched tracks / # trigger tracks in CINT1B");
  //hMatchedOverTriggerCINT1B->LabelsOption("u");
  hMatchedOverTriggerCINT1B->SetLabelSize(0.02);
  hMatchedOverTriggerCINT1B->SetLineWidth(2);
  hMatchedOverTriggerCINT1B->SetLineColor(kBlue);
  hMatchedOverTriggerCINT1B->Draw("e");

  cRatioTrackCINT1B->cd(3);
  TH1* hMatchedOverTrackerCINT1B = static_cast<TH1*>(hMatchedCINT1B->Clone("hMatchedOverTrackerCINT1B"));
  hMatchedOverTrackerCINT1B->Divide(hTrackerCINT1B);
  hMatchedOverTrackerCINT1B->SetTitle("# matched tracks / # tracker tracks in CINT1B");
  //hMatchedOverTrackerCINT1B->LabelsOption("u");
  hMatchedOverTrackerCINT1B->SetLabelSize(0.02);
  hMatchedOverTrackerCINT1B->SetLineWidth(2);
  hMatchedOverTrackerCINT1B->SetLineColor(kBlue);
  hMatchedOverTrackerCINT1B->Draw("e");
  
  cRatioTrackCINT1B->Print(OutFileNamePDF.Data());
  cRatioTrackCINT1B->Write();

  //===========================================
  // Draw ratio of tracks over CMUS1B versus run
  TH1* hTrackerPerCMUS1B= static_cast<TH1*>(hTrackerCMUS1B->Clone("hTrackerPerCMUS1B"));
  hTrackerPerCMUS1B->Divide(hCMUS1B);
  hTrackerPerCMUS1B->SetLineWidth(2);
  hTrackerPerCMUS1B->SetLineColor(kRed);

  TH1* hTriggerPerCMUS1B= static_cast<TH1*>(hTriggerCMUS1B->Clone("hTriggerPerCMUS1B"));
  hTriggerPerCMUS1B->Divide(hCMUS1B);
  hTriggerPerCMUS1B->SetLineWidth(2);
  hTriggerPerCMUS1B->SetLineColor(kBlue);

  TH1* hMatchedPerCMUS1B= static_cast<TH1*>(hMatchedCMUS1B->Clone("hMatchedPerCMUS1B"));
  hMatchedPerCMUS1B->Divide(hCMUS1B);
  hMatchedPerCMUS1B->SetLineWidth(2);
  hMatchedPerCMUS1B->SetLineColor(kViolet);

  TH1* hAllTracksPerCMUS1B= static_cast<TH1*>(hAllTracksCMUS1B->Clone("hAllTracksPerCMUS1B"));
  hAllTracksPerCMUS1B->Divide(hCMUS1B);
  hAllTracksPerCMUS1B->SetLineWidth(3);
  hAllTracksPerCMUS1B->SetLineColor(kBlack);


  CanvasName =  LHCPeriod.Data() ; 
  CanvasName += "_RatioTrackTypesCMUS1B"; 
  TCanvas *cRatioTrackTypesCMUS1B = new TCanvas(CanvasName.Data(),"cRatioTrackTypesCMUS1B",1200,900);
  cRatioTrackTypesCMUS1B->SetRightMargin(0.01);
  cRatioTrackTypesCMUS1B->SetGridy(1);
  cRatioTrackTypesCMUS1B->cd();

  hAllTracksPerCMUS1B->SetTitle("Ratio (Number of Tracks)/CMUS1B");
  hAllTracksPerCMUS1B->SetMinimum(0.01);
  hAllTracksPerCMUS1B->SetLabelSize(0.02);
  hAllTracksPerCMUS1B->Draw("E");

  hTrackerPerCMUS1B->Draw("Esame");
  hMatchedPerCMUS1B->Draw("Esame");
  hTriggerPerCMUS1B->Draw("Esame");

  TLegend* legcTTCMUS1B = new TLegend(0.75,0.55,0.90,0.75);
  legcTTCMUS1B->AddEntry(hAllTracksPerCMUS1B,"All tracks","l");
  legcTTCMUS1B->AddEntry(hTrackerPerCMUS1B,"Tracking (only) tracks","l");
  legcTTCMUS1B->AddEntry(hMatchedPerCMUS1B,"Matched tracks","l");
  legcTTCMUS1B->AddEntry(hTriggerPerCMUS1B,"Trigger (only) tracks","l");
  legcTTCMUS1B->Draw("same");

  cRatioTrackTypesCMUS1B->Print(OutFileNamePDF.Data());
  cRatioTrackTypesCMUS1B->Write();

  //===========================================
  // draw ratio of track (summed, i.e all trigger=(matched+trigger-only)) over CMUS1B versus run
  TCanvas* cTrackMultCMUS1B = new TCanvas("cTrackMultCMUS1B","cTrackMultCMUS1B",1200,900);
  cTrackMultCMUS1B->Divide(1,2);
  cTrackMultCMUS1B->cd(1);
  TH1* hSumTriggerOverCMUS1B = static_cast<TH1*>(hTriggerCMUS1B->Clone("hSumTriggerOverCMUS1B"));
  hSumTriggerOverCMUS1B->Add(hMatchedCMUS1B);
  hSumTriggerOverCMUS1B->Divide(hCMUS1B);
  hSumTriggerOverCMUS1B->SetTitle("Sum of trigger tracks (matched+trigger-only) in CMUS1B events / # CMUS1B events");
  //hSumTriggerOverCMUS1B->LabelsOption("u");
  hSumTriggerOverCMUS1B->SetLabelSize(0.02);
  hSumTriggerOverCMUS1B->SetLineWidth(2);
  hSumTriggerOverCMUS1B->SetLineColor(kRed);
  hSumTriggerOverCMUS1B->Draw("e");

  cTrackMultCMUS1B->cd(2);
  TH1* hSumTrackerOverCMUS1B = static_cast<TH1*>(hTrackerCMUS1B->Clone("hSumTrackerOverCMUS1B"));
  hSumTrackerOverCMUS1B->Add(hMatchedCMUS1B);
  hSumTrackerOverCMUS1B->Divide(hCMUS1B);
  hSumTrackerOverCMUS1B->SetTitle("Sum of tracker tracks (matched+tracker-only) in CMUS1B events / # CMUS1B events");
  //hSumTrackerOverCMUS1B->LabelsOption("u");
  hSumTrackerOverCMUS1B->SetLabelSize(0.02);
  hSumTrackerOverCMUS1B->SetLineWidth(2);
  hSumTrackerOverCMUS1B->SetLineColor(kRed);
  hSumTrackerOverCMUS1B->Draw("e");

  cTrackMultCMUS1B->Print(OutFileNamePDF.Data());
  cTrackMultCMUS1B->Write();
  
  //===========================================
  // draw mixed ratio of track over track versus run for CMUS1B
  TCanvas* cRatioTrackCMUS1B = new TCanvas("cRatioTrackCMUS1B","cRatioTrackCMUS1B",1200,900);
  cRatioTrackCMUS1B->Divide(1,3);

  cRatioTrackCMUS1B->cd(1);
  TH1* hTrackerOverTriggerCMUS1B = static_cast<TH1*>(hTrackerCMUS1B->Clone("hTrackerOverTriggerCMUS1B"));
  hTrackerOverTriggerCMUS1B->Divide(hTriggerCMUS1B);
  hTrackerOverTriggerCMUS1B->SetTitle("# tracker tracks / # trigger tracks in CMUS1B");
  //hTrackerOverTriggerCMUS1B->LabelsOption("u");
  hTrackerOverTriggerCMUS1B->SetLabelSize(0.02);
  hTrackerOverTriggerCMUS1B->SetLineWidth(2);
  hTrackerOverTriggerCMUS1B->SetLineColor(kRed);
  hTrackerOverTriggerCMUS1B->Draw("e");

  cRatioTrackCMUS1B->cd(2);
  TH1* hMatchedOverTriggerCMUS1B = static_cast<TH1*>(hMatchedCMUS1B->Clone("hMatchedOverTriggerCMUS1B"));
  hMatchedOverTriggerCMUS1B->Divide(hTriggerCMUS1B);
  hMatchedOverTriggerCMUS1B->SetTitle("# matched tracks / # trigger tracks in CMUS1B");
  //hMatchedOverTriggerCMUS1B->LabelsOption("u");
  hMatchedOverTriggerCMUS1B->SetLabelSize(0.02);
  hMatchedOverTriggerCMUS1B->SetLineWidth(2);
  hMatchedOverTriggerCMUS1B->SetLineColor(kRed);
  hMatchedOverTriggerCMUS1B->Draw("e");

  cRatioTrackCMUS1B->cd(3);
  TH1* hMatchedOverTrackerCMUS1B = static_cast<TH1*>(hMatchedCMUS1B->Clone("hMatchedOverTrackerCMUS1B"));
  hMatchedOverTrackerCMUS1B->Divide(hTrackerCMUS1B);
  hMatchedOverTrackerCMUS1B->SetTitle("# matched tracks / # tracker tracks in CMUS1B");
  //hMatchedOverTrackerCMUS1B->LabelsOption("u");
  hMatchedOverTrackerCMUS1B->SetLabelSize(0.02);
  hMatchedOverTrackerCMUS1B->SetLineWidth(2);
  hMatchedOverTrackerCMUS1B->SetLineColor(kRed);
  hMatchedOverTrackerCMUS1B->Draw("e");
  
  cRatioTrackCMUS1B->Print(OutFileNamePDF.Data());
  cRatioTrackCMUS1B->Write();

  //==================================================
  // Draw matched tracks asymetry :  
  // for all triggers  
  TH1 *hDiffMatched= static_cast<TH1*>(hPosMatched->Clone("hDiffMatched"));
  hDiffMatched->Add(hNegMatched,-1);
  hDiffMatched->Sumw2();
  
  TH1 *hAsymMatched= static_cast<TH1*>(hDiffMatched->Clone("hAsymMatched"));
  hAsymMatched->Divide(hAllMatched);
  hAsymMatched->SetTitle("Asymetry Matched tracks for Physics Selected events");
  hAsymMatched->SetMarkerStyle(20);
  

  TH1 *hDiffMatchedCMUS1B= static_cast<TH1*>(hPosMatchedCMUS1B->Clone("hDiffMatchedCMUS1B"));
  hDiffMatchedCMUS1B->Add(hNegMatchedCMUS1B,-1);
  hDiffMatchedCMUS1B->Sumw2();
  
  TH1 *hAsymMatchedCMUS1B= static_cast<TH1*>(hDiffMatchedCMUS1B->Clone("hAsymMatchedCMUS1B"));
  hAsymMatchedCMUS1B->Divide(hAllMatchedCMUS1B);
  hAsymMatchedCMUS1B->SetLineColor(kRed);
  hAsymMatchedCMUS1B->SetLineWidth(2);
  hAsymMatchedCMUS1B->SetMinimum(-0.1);
  hAsymMatchedCMUS1B->SetMaximum(0.1);
  hAsymMatchedCMUS1B->SetLabelSize(0.02);
  hAsymMatchedCMUS1B->SetTitle("Matched tracks asymetry in Physics Selected events for CMUS1B");
  
  CanvasName =  LHCPeriod.Data() ; 
  CanvasName += "_AsymMatched"; 
  TCanvas *cAsymMatched = new TCanvas(CanvasName.Data(),"cAsymMatched",1200,900);
  cAsymMatched->SetRightMargin(0.01);
  cAsymMatched->SetGridy(1);
  cAsymMatched->cd();
  hAsymMatchedCMUS1B->GetYaxis()->SetTitle("Asymetry");  
  hAsymMatchedCMUS1B->Draw("EH");

  cAsymMatched->Print(OutFileNamePDF.Data());
  cAsymMatched->Write();


  //==================================================
  // Draw high pt tracks per CMUS1B
  TH1* hMatchedLowPtPerCMUS1B = static_cast<TH1*> (hMatchedLowPtCMUS1B->Clone("hMatchedLowPtPerCMUS1B"));
  hMatchedLowPtPerCMUS1B->Sumw2();
  hMatchedLowPtPerCMUS1B->Divide(hCMUS1B);
  hMatchedLowPtPerCMUS1B->SetLineWidth(2);
  hMatchedLowPtPerCMUS1B->SetLineColor(kBlue);

  TH1* hMatchedHighPtPerCMUS1B = static_cast<TH1*> (hMatchedHighPtCMUS1B->Clone("hMatchedHighPtPerCMUS1B"));
  hMatchedHighPtPerCMUS1B->Sumw2();
  hMatchedHighPtPerCMUS1B->Divide(hCMUS1B);
  hMatchedHighPtPerCMUS1B->SetLineWidth(2);
  hMatchedHighPtPerCMUS1B->SetLineColor(kRed);

  CanvasName =  LHCPeriod.Data() ; 
  CanvasName += "_HighPtMuons"; 
  TCanvas *cHighPtMuons = new TCanvas(CanvasName.Data(),"cHighPtMuons",1200,900);
  cHighPtMuons->SetTopMargin(0.05);
  cHighPtMuons->SetRightMargin(0.01);
  cHighPtMuons->SetGridy(1);
  cHighPtMuons->cd();

  hMatchedLowPtPerCMUS1B->SetTitle("");
  hMatchedLowPtPerCMUS1B->GetYaxis()->SetTitle("Ratio per CMUS1B");
  //hMatchedLowPtPerCMUS1B->SetMaximum(0.15);
  hMatchedLowPtPerCMUS1B->SetMinimum(0.0001);
  hMatchedLowPtPerCMUS1B->SetLabelSize(0.02);

  hMatchedLowPtPerCMUS1B->Draw("E");
  hMatchedHighPtPerCMUS1B->Draw("Esame");
  
  TLegend* legcHPM = new TLegend(0.60,0.45,0.98,0.65);
  legcHPM->SetHeader("Number of matched track per CMUS1B (include Vtx and R_{Abs} cuts) ");
  legcHPM->AddEntry(hMatchedLowPtPerCMUS1B," p_{T} > 1 GeV/c ","l");
  legcHPM->AddEntry(hMatchedHighPtPerCMUS1B," p_{T} >  2 GeV/c ","l");
  legcHPM->Draw("same");

  cHighPtMuons->Print(OutFileNamePDF.Data());
  cHighPtMuons->Write(); 

  // close merged file
  globalFile->Close();

  //--------------------------------------------- //
  //        monitor quantities run per run        //
  //--------------------------------------------- //

  TH1F* hNClustersPerTrackVsRun_Mean = new TH1F("hNClustersPerTrackVsRun_Mean", "averaged number of associated clusters per track;run;<n_{clusters}>",10000,1,10000);
  TH1F* hNClustersPerTrackVsRun_Sigma = new TH1F("hNClustersPerTrackVsRun_Sigma", "dispersion of the number of associated clusters per track;run;#sigma_{n_{clusters}}",10000,1,10000);
  TH1F* hNChamberHitPerTrack_Mean = new TH1F("hNChamberHitPerTrack_Mean", "averaged number of chambers hit per track;run;<n_{chamber hit}>",10000,1,10000);
  TH1F* hNChamberHitPerTrack_Sigma = new TH1F("hNChamberHitPerTrack_Sigma", "dispersion of the number of chambers hit per track;run;#sigma_{n_{chamber hit}}",10000,1,10000);
  TH1F* hChi2_Mean = new TH1F("hChi2_Mean", "averaged normalized #chi^{2} distribution;run;<#chi^{2} / ndf>",10000,1,10000);
  TH1F* hChi2_Sigma = new TH1F("hChi2_Sigma", "dispersion of normalized #chi^{2} distribution;run;#sigma_{#chi^{2} / ndf}",10000,1,10000);
  TH1F* hNClustersInCh[10];
  for (Int_t ich=0; ich<10; ich++) hNClustersInCh[ich] = new TH1F(Form("hNClustersInCh%d",ich+1), Form("averaged number of clusters in chamber %d per track;run;<n_{clusters}>",ich+1),10000,1,10000);
  TH1F* hClusterHitMapXInCh[10];
  for (Int_t ich=0; ich<10; ich++) hClusterHitMapXInCh[ich] = new TH1F(Form("hClusterHitMapXInCh%d",ich+1), Form("averaged cluster position distribution in chamber %d;X (cm)",ich+1),10000,1,10000);
  TH1F* hClusterHitMapYInCh[10];
  for (Int_t ich=0; ich<10; ich++) hClusterHitMapYInCh[ich] = new TH1F(Form("hClusterHitMapYInCh%d",ich+1), Form("averaged cluster position distribution in chamber %d;Y (cm)",ich+1),10000,1,10000);
  
  Int_t ibin = 1;
  
  // Are the runs stored locally or in alien?
  Int_t isAlienFile = 0;
  if(alienBaseDir.Contains("alien:")){
    isAlienFile = 1;
    alienBaseDir.ReplaceAll("alien://","");
  }
  
  // Loop over runs
  for ( Int_t irun=0; irun<runs.GetEntriesFast(); irun++ ) {
    
    TString run = ((TObjString*)runs.UncheckedAt(irun))->GetString();
    cout << "processing run nr="<<run <<endl;
    // get the file (or list of files) to be analyzed
    TString command;
    TGridResult *res = 0;
    TObjString *objs = 0;	
    
    if(isAlienFile){
      command = Form("find %s/ %s/AnalysisResults.root", alienBaseDir.Data(), run.Data());
      res = gGrid->Command(command);
      if (!res) {
	Error("PlotQA",Form("no result for the command: %s",command.Data()));
	return;
      }
    }
    else{
      res = new TGridResult();	
      
      if(runList){
	objs = new TObjString(Form("%s/%s/AnalysisResults.root", alienBaseDir.Data(), run.Data()));
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
	  objs = new TObjString(Form("%s/%s/AnalysisResults.root", alienBaseDir.Data(), sDirFilename.Data()));
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
	Error("PlotQA","turl/obj not found for the run %s... SKIPPING", run.Data());
	continue;
      }
      
      // open the outfile for this run
      TFile *runFile = TFile::Open(objs->GetString());
      if (!runFile || ! runFile->IsOpen()) {
	Error("PlotQA", Form("failed to open file: %s", objs->GetName()));
	continue;//return;
      }
      runFile->Cd("MUON_QA");
      
      // get interesting histos
      TObjArray* general1 = static_cast<TObjArray*>(runFile->FindObjectAny("general1"));
			TObjArray* general2 = static_cast<TObjArray*>(runFile->FindObjectAny("general2"));
      TObjArray* expert = static_cast<TObjArray*>(runFile->FindObjectAny("expert"));
	
			if (!general1 || !general2 || !expert){
			 Error("PlotMUONQA", Form("All objects not here !!! ===> Skipping..."));		
			 continue;
			}
			TH1* hNClustersPerTrack = static_cast<TH1*>(general1->FindObject("hNClustersPerTrack"));
      TH1* hNChamberHitPerTrack = static_cast<TH1*>(general1->FindObject("hNChamberHitPerTrack"));
      TH1* hChi2 = static_cast<TH1*>(general1->FindObject("hChi2"));
			TH1* hNClustersPerCh = static_cast<TH1*>(general2->FindObject("hNClustersPerCh"));
      
				
      TH2* hClusterHitMapInCh[10];
      for(Int_t ich=0; ich<10; ich++) hClusterHitMapInCh[ich] = static_cast<TH2*>(expert->FindObject(Form("hClusterHitMapInCh%d",ich+1)));
      
      // skip empty runs... not anymore ! cs !
      if (!hNClustersPerCh) {
	Warning("PlotQA", Form("File: %s has empty histograms !", objs->GetName()));
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
      run.ReplaceAll(baseDir, "");
      run.Remove(TString::kLeading, '/');
      run.Remove(TString::kLeading, '0');
      run.ReplaceAll("/AnalysisResults.root", "");
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
      
      // close outfile for this run
      runFile->Close();
      ibin++;      
    }
    
    delete res;
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
  hNClustersPerTrackVsRun_Mean->SetLabelSize(0.02);
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
  hNClustersPerTrackVsRun_Sigma->SetLabelSize(0.02);
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
  hClusterHitMapXInCh[0]->SetLabelSize(0.02);
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
  hClusterHitMapYInCh[0]->SetLabelSize(0.02);
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
  hChi2_Mean->SetLabelSize(0.02);
  hChi2_Mean->SetLineWidth(2);
  hChi2_Mean->Draw("e");

  cChi2->cd(2);
  hChi2_Sigma->SetStats(kFALSE);
  hChi2_Sigma->GetXaxis()->SetRange(1,ibin-1);
  hChi2_Sigma->GetXaxis()->SetNdivisions(1,kFALSE);
  //hChi2_Sigma->LabelsOption("u");
  hChi2_Sigma->SetLabelSize(0.02);
  hChi2_Sigma->SetLineWidth(2);
  hChi2_Sigma->Draw("e");

  cChi2->Print(OutFileNamePDF.Data());
  cChi2->Write();


  // close the PDF file
  c1->Print(OutFileNamePDF_close.Data());
  rootFileOut->Close();
  
  
  //====================================================
  if (PRINTSTAT){
    // set the format to print labels
    THashList* labels = hCMUS1B->GetXaxis()->GetLabels();
    TString format(Form("\n%%%ds %%9d",0));
		Int_t nRuns=0;
	
    // print value for each label
    TObjString* label = 0x0;
    TIter nextLabel(labels);
    cout << "-------------------------------------------------" << endl;
    cout << "Run Number" << "\t Number of CMUS1B after Phys. Sel. " << endl ;  
    while ((label = static_cast<TObjString*>(nextLabel()))) {
			nRuns++;
      Int_t bin = (Int_t) label->GetUniqueID();
      printf(format.Data(), label->String().Data(), (Int_t) hCMUS1B->GetBinContent(bin));
    }
		printf("\n========== Total #runs = %d ==============\n",nRuns);
    printf("\n\n");


    cout << "-------------------------------------------------" << endl;
    cout << "Total statistic" << endl; 
    cout << " " << endl ; 
    cout << "Number of  NumOfCINT1B " << endl ;
    cout << "\t before selection " << NumOfCINT1BNoPS << "\t  after selection " <<   NumOfCINT1BWithPS << " --> rejection = " <<  (Double_t) (NumOfCINT1BNoPS-NumOfCINT1BWithPS)/(NumOfCINT1BNoPS)*100. << "%" << endl ; 
    cout << " " << endl ; 
    cout << "Number of  NumOfCMUS1B " << endl ;
    cout << "\t before selection " << NumOfCMUS1BNoPS << "\t  after selection " <<   NumOfCMUS1BWithPS << " --> rejection = " <<  (Double_t) (NumOfCMUS1BNoPS-NumOfCMUS1BWithPS)/(NumOfCMUS1BNoPS)*100. << "%" << endl ; 
    cout << " " << endl ; 
    if (NumOfCSH1BNoPS>0){
      cout << "Number of  NumOfCSH1B " << endl ;
      cout << "\t before selection " << NumOfCSH1BNoPS << "\t  after selection " <<   NumOfCSH1BWithPS << " --> rejection = " <<  (Double_t) (NumOfCSH1BNoPS-NumOfCSH1BWithPS)/(NumOfCSH1BNoPS)*100. << "%" << endl ; 
      cout << " " << endl ;  
    }
  }
}
