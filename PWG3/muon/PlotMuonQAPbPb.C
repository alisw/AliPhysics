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
// The macro reads results of the QA task and produce monitoring plots.
//
// Author: Philippe Pillot - SUBATECH Nantes, Christophe Suire - IPN Orsay 
// and Cynthia Hadjidakis - IPN Orsay
//--------------------------------------------------------------------------

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>

// ROOT includes
#include "TMath.h"
#include "TGrid.h"
#include "TFile.h"
#include "TH1.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPad.h"

// ALIROOT includes
#include "../PWG3/base/AliCounterCollection.h"

#endif

//--------------------------------------------------------------------------
void PlotMuonQAPbPb(const char* baseDir, const char* runList = 0x0, Bool_t selectPhysics = kFALSE)
{
  /// Macro for QA monitoring.
  /// Examples: 
  /// baseDir = "alien:///alice/cern.ch/user/p/ppillot/pp7TeV/LHC10d/MuonQA/pass1/results"
  /// baseDir = "results"
  /// If runList != 0x0: only the given runs will be used. Otherwise use all runs found in baseDir.
  /// Usage (local/grid)
  /// .x PlotMuonQAPbPb.C("alien:///alice/cern.ch/user/s/suire/QAP2LHC10h/output",0,kTRUE)
  

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

  TString LHCPeriod = "LHC10hpass2"; 

  TString OutFileName = "muonQA_";  OutFileName += LHCPeriod;
  TString OutFileNamePDF=  OutFileName.Data();  OutFileNamePDF+= ".pdf";
  TString OutFileNamePDF_open = OutFileNamePDF.Data(); OutFileNamePDF_open += "[";  
  TString OutFileNamePDF_close= OutFileNamePDF.Data(); OutFileNamePDF_close += "]";  
  TString OutFileNameROOT=  OutFileName.Data();  OutFileNameROOT+= ".root";

  //TString QAFileName = "QAresults.root";
  TString QAFileName = "AnalysisResults.root";
  TString MergedQAFileName = "Merged"; MergedQAFileName+= QAFileName.Data();

  Int_t PRINTSTAT = 1; 

  Int_t ALIENFASTREAD = 0 ; 

  if (ALIENFASTREAD){ // Equivalent to the fast read option
    gEnv->SetValue("XNet.ConnectTimeout",10);
    gEnv->SetValue("XNet.RequestTimeout",10);
    gEnv->SetValue("XNet.MaxRedirectCount",2);
    gEnv->SetValue("XNet.ReconnectTimeout",10);
    gEnv->SetValue("XNet.FirstConnectMaxCnt",1);
  }

  TH1::AddDirectory(kFALSE);
  
  TString alienBaseDir = baseDir;
  
  if (alienBaseDir.Contains("alien:") && !TGrid::Connect("alien://")) {
      Error("PlotMuonQA","cannot connect to grid");
      return;
  }
  
  Float_t LabelSize = 0.03; 
  
    
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
      Error("PlotMuonQAPbPb",Form("unable to open file %s", runList));
      return;
    }
    
    TString currRun;
    while (!inFile.eof()) {
      currRun.ReadLine(inFile, kTRUE);
      if (currRun.IsNull()) continue;
      if (!currRun.IsDigit()) {
	Error("PlotMuonQAPbPb","invalid run number: %s", currRun.Data());
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
  
  TFile *globalFile = TFile::Open(Form("%s/%s", baseDir,MergedQAFileName.Data()));
  if (!globalFile || ! globalFile->IsOpen()) {
    Error("PlotQA", Form("failed to open file: %s/%s", baseDir,MergedQAFileName.Data()));
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

  //TString allTriggers = "trigger:CMBAC-B,CMBS2A-B,CMBS2C-B,CMBACS2-B,C0SMH-B";
  //TString minBiasTrigers = "trigger:CMBAC-B,CMBS2A-B,CMBS2C-B,CMBACS2-B";
  TString allTriggers = "trigger:any"; //"trigger:CMBAC-B,CMBS2A-B,CMBS2C-B,C0SMH-B";
  TString minBiasTrigers = "trigger:any"; //CMBAC-B,CMBS2A-B,CMBS2C-B";
  TString bckTriggers = "trigger:CMBACS2-A,CMBACS2-C,CMBS2A-A,CMBS2A-C,CMBS2C-A,CMBS2C-C,CMBAC-A,CMBAC-A,C0SMH-A,C0SMH-C";

  // Histo trigger without Phys. Sel. 
  TH1* hAllTriggersNoPS = eventCounters->Draw("run",Form("%s/%s",allTriggers.Data(), selectRuns.Data()));
  hAllTriggersNoPS->Sumw2();
  Int_t NumOfAllTriggersNoPS = hAllTriggersNoPS->Integral();
  TH1* hAllTriggerslowNoPS = eventCounters->Draw("run",Form("v0mult:low/%s/%s",allTriggers.Data(), selectRuns.Data()));
  hAllTriggerslowNoPS->Sumw2();
  TH1* hAllTriggershighNoPS = eventCounters->Draw("run",Form("v0mult:high/%s/%s",allTriggers.Data(), selectRuns.Data()));
  hAllTriggershighNoPS->Sumw2();
  TH1* hMBBNoPS = eventCounters->Draw("run",Form("v0mult:low,int,high/%s/%s",minBiasTrigers.Data(), selectRuns.Data()));
  hMBBNoPS->Sumw2();
  Int_t NumOfMBBNoPS = hMBBNoPS->Integral();
  TH1* hMBBlowNoPS = eventCounters->Draw("run",Form("v0mult:low/%s/%s",minBiasTrigers.Data(), selectRuns.Data()));
  hMBBlowNoPS->Sumw2();
  Int_t NumOfMBBlowNoPS = hMBBlowNoPS->Integral();
  TH1* hMBBhighNoPS = eventCounters->Draw("run",Form("v0mult:high/%s/%s",minBiasTrigers.Data(), selectRuns.Data()));
  hMBBhighNoPS->Sumw2();
  Int_t NumOfMBBhighNoPS = hMBBhighNoPS->Integral();
  //TH1* hCMUS1BNoPS = eventCounters->Draw("run",Form("trigger:CMUS1B/%s", selectRuns.Data()));
  //hCMUS1BNoPS->Sumw2();
  //Int_t NumOfCMUS1BNoPS = hCMUS1BNoPS->Integral();
  //TH1* hCSH1BNoPS = eventCounters->Draw("run",Form("trigger:CSH1B/%s", selectRuns.Data()));
  //hCSH1BNoPS->Sumw2();
  //Int_t NumOfCSH1BNoPS = hCSH1BNoPS->Integral();

  TH1* hBckTriggersNoPS = eventCounters->Draw("run",Form("v0mult:low,int,high/%s/%s",bckTriggers.Data(), selectRuns.Data()));
  hBckTriggersNoPS->Sumw2();				      
					      

  // Histo trigger with Phys. Sel. 
  TH1* hAllTriggersWithPS = eventCounters->Draw("run",Form("%s/%s/selected:yes",allTriggers.Data(), selectRuns.Data()));
  hAllTriggersWithPS->Sumw2();
  Int_t  NumOfAllTriggersWithPS = hAllTriggersWithPS->Integral();
  TH1* hAllTriggersLowWithPS = eventCounters->Draw("run",Form("v0mult:low/%s/%s/selected:yes",allTriggers.Data(), selectRuns.Data()));
  hAllTriggersWithPS->Sumw2();
  TH1* hAllTriggersHighWithPS = eventCounters->Draw("run",Form("v0mult:high/%s/%s/selected:yes",allTriggers.Data(), selectRuns.Data()));
  hAllTriggersWithPS->Sumw2();
  TH1* hMBBWithPS = eventCounters->Draw("run",Form("v0mult:low,int,high/%s/%s/selected:yes",minBiasTrigers.Data(), selectRuns.Data()));
  hMBBWithPS->Sumw2();
  Int_t NumOfMBBWithPS = hMBBWithPS->Integral();
  TH1* hMBBlowWithPS = eventCounters->Draw("run",Form("v0mult:low/%s/%s/selected:yes",minBiasTrigers.Data(), selectRuns.Data()));
  hMBBlowWithPS->Sumw2();
  Int_t NumOfMBBlowWithPS = hMBBlowWithPS->Integral();
  TH1* hMBBhighWithPS = eventCounters->Draw("run",Form("v0mult:high/%s/%s/selected:yes",minBiasTrigers.Data(), selectRuns.Data()));
  hMBBhighWithPS->Sumw2();
  Int_t NumOfMBBhighWithPS = hMBBhighWithPS->Integral();
  
  // Histo trigger : Phys. Sel.  is selected or not depending on the macro arguments
  TH1* hAllTriggers = eventCounters->Draw("run",Form("%s/%s/%s",allTriggers.Data(), selectRuns.Data(), select.Data()));
  hAllTriggers->Sumw2();
  TH1* hMBB = eventCounters->Draw("run",Form("v0mult:low,int,high/%s/%s/%s",minBiasTrigers.Data(), selectRuns.Data(), select.Data()));
  hMBB->Sumw2();
  TH1* hMBBlow = eventCounters->Draw("run",Form("v0mult:low/%s/%s/%s",minBiasTrigers.Data(), selectRuns.Data(), select.Data()));
  hMBBlow->Sumw2();
  TH1* hMBBhigh = eventCounters->Draw("run",Form("v0mult:high/%s/%s/%s",minBiasTrigers.Data(), selectRuns.Data(), select.Data()));
  hMBBhigh->Sumw2();

  // Histo tracking : Phys. Sel.  is selected or not depending on the macro arguments
  TH1* hAllTracksMBB = trackCounters->Draw("run",Form("v0mult:low,int,high/%s/%s/%s",minBiasTrigers.Data(), selectRuns.Data(), select.Data()));
  hAllTracksMBB->Sumw2();
  TH1* hTriggerMBB = trackCounters->Draw("run",Form("v0mult:low,int,high/track:triggeronly/%s/%s/%s",minBiasTrigers.Data(), selectRuns.Data(), select.Data()));
  hTriggerMBB->Sumw2();
  TH1* hTrackerMBB = trackCounters->Draw("run",Form("v0mult:low,int,high/track:trackeronly/%s/%s/%s",minBiasTrigers.Data(), selectRuns.Data(), select.Data()));
  hTrackerMBB->Sumw2();
  TH1* hMatchedMBB = trackCounters->Draw("run",Form("v0mult:low,int,high/track:matched/%s/%s/%s",minBiasTrigers.Data(), selectRuns.Data(), select.Data()));
  hMatchedMBB->Sumw2();
  TH1* hAllTracksMBBlow = trackCounters->Draw("run",Form("v0mult:low/%s/%s/%s",minBiasTrigers.Data(), selectRuns.Data(), select.Data()));
  hAllTracksMBBlow->Sumw2();
  TH1* hTriggerMBBlow = trackCounters->Draw("run",Form("v0mult:low/track:triggeronly/%s/%s/%s",minBiasTrigers.Data(), selectRuns.Data(), select.Data()));
  hTriggerMBBlow->Sumw2();
  TH1* hTrackerMBBlow = trackCounters->Draw("run",Form("v0mult:low/track:trackeronly/%s/%s/%s",minBiasTrigers.Data(), selectRuns.Data(), select.Data()));
  hTrackerMBBlow->Sumw2();
  TH1* hMatchedMBBlow = trackCounters->Draw("run",Form("v0mult:low/track:matched/%s/%s/%s",minBiasTrigers.Data(), selectRuns.Data(), select.Data()));
  hMatchedMBBlow->Sumw2();
  TH1* hAllTracksMBBhigh = trackCounters->Draw("run",Form("v0mult:high/%s/%s/%s",minBiasTrigers.Data(), selectRuns.Data(), select.Data()));
  hAllTracksMBBhigh->Sumw2();
  TH1* hTriggerMBBhigh = trackCounters->Draw("run",Form("v0mult:high/track:triggeronly/%s/%s/%s",minBiasTrigers.Data(), selectRuns.Data(), select.Data()));
  hTriggerMBBhigh->Sumw2();
  TH1* hTrackerMBBhigh = trackCounters->Draw("run",Form("v0mult:high/track:trackeronly/%s/%s/%s",minBiasTrigers.Data(), selectRuns.Data(), select.Data()));
  hTrackerMBBhigh->Sumw2();
  TH1* hMatchedMBBhigh = trackCounters->Draw("run",Form("v0mult:high/track:matched/%s/%s/%s",minBiasTrigers.Data(), selectRuns.Data(), select.Data()));
  hMatchedMBBhigh->Sumw2();
  
  TH1* hMatchedLowPtMBB = trackCounters->Draw("run",Form("v0mult:low,int,high/track:matched/acc:in/%s/%s/%s/pt:low",minBiasTrigers.Data(), selectRuns.Data(), select.Data()));
  hMatchedLowPtMBB->Sumw2();
  TH1* hMatchedHighPtMBB = trackCounters->Draw("run",Form("v0mult:low,int,high/track:matched/acc:in/%s/%s/%s/pt:high",minBiasTrigers.Data(), selectRuns.Data(), select.Data()));
  hMatchedHighPtMBB->Sumw2();
  TH1* hMatchedLowPtMBBlow = trackCounters->Draw("run",Form("v0mult:low/track:matched/acc:in/%s/%s/%s/pt:low",minBiasTrigers.Data(), selectRuns.Data(), select.Data()));
  hMatchedLowPtMBBlow->Sumw2();
  TH1* hMatchedHighPtMBBlow = trackCounters->Draw("run",Form("v0mult:low/track:matched/acc:in/%s/%s/%s/pt:high",minBiasTrigers.Data(), selectRuns.Data(), select.Data()));
  hMatchedHighPtMBBlow->Sumw2();
  TH1* hMatchedLowPtMBBhigh = trackCounters->Draw("run",Form("v0mult:high/track:matched/acc:in/%s/%s/%s/pt:low",minBiasTrigers.Data(), selectRuns.Data(), select.Data()));
  hMatchedLowPtMBBhigh->Sumw2();
  TH1* hMatchedHighPtMBBhigh = trackCounters->Draw("run",Form("v0mult:high/track:matched/acc:in/%s/%s/%s/pt:high",minBiasTrigers.Data(), selectRuns.Data(), select.Data()));
  hMatchedHighPtMBBhigh->Sumw2();
  
  
  // for MBB only, low and high mult
  TH1 *hPosMatchedMBBlow =  trackCounters->Draw("run",Form("v0mult:low/track:matched/charge:pos/acc:in/%s/%s/%s",minBiasTrigers.Data(),select.Data(),selectRuns.Data()));
  hPosMatchedMBBlow->Sumw2();
  TH1 *hNegMatchedMBBlow=  trackCounters->Draw("run",Form("v0mult:low/track:matched/charge:neg/acc:in/%s/%s/%s",minBiasTrigers.Data(),select.Data(),selectRuns.Data())); 
  hNegMatchedMBBlow->Sumw2();
  TH1 *hAllMatchedMBBlow=  trackCounters->Draw("run",Form("v0mult:low/track:matched/acc:in/%s/%s/%s",minBiasTrigers.Data(),select.Data(),selectRuns.Data())); 
  hAllMatchedMBBlow->Sumw2();
  TH1 *hPosMatchedMBBhigh =  trackCounters->Draw("run",Form("v0mult:high/track:matched/charge:pos/acc:in/%s/%s/%s",minBiasTrigers.Data(),select.Data(),selectRuns.Data()));
  hPosMatchedMBBhigh->Sumw2();
  TH1 *hNegMatchedMBBhigh=  trackCounters->Draw("run",Form("v0mult:high/track:matched/charge:neg/acc:in/%s/%s/%s",minBiasTrigers.Data(),select.Data(),selectRuns.Data())); 
  hNegMatchedMBBhigh->Sumw2();
  TH1 *hAllMatchedMBBhigh=  trackCounters->Draw("run",Form("v0mult:high/track:matched/acc:in/%s/%s/%s",minBiasTrigers.Data(),select.Data(),selectRuns.Data())); 
  hAllMatchedMBBhigh->Sumw2();
  

  TH1* hAll = eventCounters->Draw("trigger","run","run:any/selected:yes");
  hAll->Draw("TEXT");
  

  //===================================================================================
  // Put all plots in a pdf file, easier to publish (Twiki)
  c1->Print(OutFileNamePDF_open.Data());

  TFile *rootFileOut = TFile::Open(OutFileNameROOT.Data(),"RECREATE");
  
  //===================================================================================
  // new canvas of the statistic wrt trigger w/  physics selection
  CanvasName =  LHCPeriod.Data() ; 
  CanvasName += "_StatByTrigger"; 
  TCanvas *cStatByTrigger = new TCanvas(CanvasName.Data(),"cStatByTrigger",1200,900);
  cStatByTrigger->SetTopMargin(0.05);
  cStatByTrigger->SetRightMargin(0.01);
  cStatByTrigger->SetGridy(1);
  cStatByTrigger->SetLogy(1);
  cStatByTrigger->cd();

  //hAllTriggersWithPS->SetLabelSize(LabelSize);
  //hAllTriggersWithPS->SetTitle("Statistic per trigger  w/  Phys. Sel."); 
  hAllTriggersWithPS->SetTitle(""); 
  hAllTriggersWithPS->SetLineWidth(2);
  hAllTriggersWithPS->SetMinimum(10);
  hAllTriggersWithPS->SetLabelSize(LabelSize);

  hAllTriggersWithPS->Draw("E");
  hMBBWithPS->SetLineWidth(2);
  hMBBWithPS->SetLineColor(kBlue);
  hMBBWithPS->Draw("ESAME");
  hMBBlowWithPS->SetMarkerStyle(27);
  hMBBlowWithPS->SetMarkerSize(2);
  hMBBlowWithPS->SetMarkerColor(kBlue);
  hMBBlowWithPS->SetLineColor(kBlue);
  hMBBlowWithPS->Draw("ESAME");
  hMBBhighWithPS->SetMarkerStyle(24);
  hMBBhighWithPS->SetMarkerSize(2);
  hMBBhighWithPS->SetMarkerColor(kBlue);
  hMBBhighWithPS->SetLineColor(kBlue);
  hMBBhighWithPS->Draw("ESAME");


  TLegend* legcSBT = new TLegend(0.12,0.15,0.42,0.3);
  legcSBT->SetHeader("Trigger Statistic w/ Phys. Sel.");
  legcSBT->AddEntry(hAllTriggersWithPS,"All collisions triggers","l");
  legcSBT->AddEntry(hMBBWithPS,"Min Bias [0-80]% from V0 amplitude","l");
  legcSBT->AddEntry(hMBBlowWithPS,"MB low multiplicity events [60-80]%","p"); 
  legcSBT->AddEntry(hMBBhighWithPS,"MB high multiplicity events [0-10]%","p"); 
  legcSBT->Draw("same");

  cStatByTrigger->Print(OutFileNamePDF.Data());
  cStatByTrigger->Write();


  //===================================================================================
  // new canvas with the relative content of each trigger w/ and w/o  physics selection
  TH1 *ratioMBB = static_cast<TH1*>(hMBBWithPS->Clone("ratioMBB"));
  ratioMBB->Divide(hMBBNoPS);
  ratioMBB->SetLineColor(kBlack);
  ratioMBB->SetLineWidth(2);
  TH1 *ratioMBBlow = static_cast<TH1*>(hMBBlowWithPS->Clone("ratioMBBlow"));
  ratioMBBlow->Divide(hMBBlowNoPS);
  ratioMBBlow->SetLineColor(kBlue);
  ratioMBBlow->SetLineWidth(2);
  ratioMBBlow->SetMarkerSize(0);
  ratioMBBlow->SetMarkerColor(kBlue);
  TH1 *ratioMBBhigh = static_cast<TH1*>(hMBBhighWithPS->Clone("ratioMBBhigh"));
  ratioMBBhigh->Divide(hMBBhighNoPS);
  ratioMBBhigh->SetLineColor(kRed);
  ratioMBBhigh->SetLineWidth(2);
  ratioMBBhigh->SetMarkerSize(0);
  ratioMBBhigh->SetMarkerColor(kRed);


  TH1 *ratioMBBNoPS = static_cast<TH1*>(hMBBNoPS->Clone("ratioMBBNoPS"));
  ratioMBBNoPS->Divide(hMBBNoPS);
  ratioMBBNoPS->SetLineColor(kBlack);
  ratioMBBNoPS->SetLineWidth(0);
  ratioMBBNoPS->SetLineStyle(2);
  ratioMBBNoPS->SetMarkerColor(kBlack);
  ratioMBBNoPS->SetMarkerStyle(24);
  ratioMBBNoPS->SetMarkerSize(1);
  TH1 *ratioMBBlowNoPS = static_cast<TH1*>(hMBBlowNoPS->Clone("ratioMBBlowNoPS"));
  ratioMBBlowNoPS->Divide(hMBBlowNoPS);
  ratioMBBlowNoPS->SetLineColor(kBlue);
  ratioMBBlowNoPS->SetLineWidth(0);
  ratioMBBlowNoPS->SetLineStyle(2);
  ratioMBBlowNoPS->SetMarkerColor(kBlue);
  ratioMBBlowNoPS->SetMarkerStyle(24);
  ratioMBBlowNoPS->SetMarkerSize(1);
  TH1 *ratioMBBhighNoPS = static_cast<TH1*>(hMBBhighNoPS->Clone("ratioMBBhighNoPS"));
  ratioMBBhighNoPS->Divide(hMBBhighNoPS);
  ratioMBBhighNoPS->SetLineColor(kRed);
  ratioMBBhighNoPS->SetLineWidth(0);
  ratioMBBhighNoPS->SetLineStyle(2);
  ratioMBBhighNoPS->SetMarkerColor(kRed);
  ratioMBBhighNoPS->SetMarkerStyle(24);
  ratioMBBhighNoPS->SetMarkerSize(1);
  

  TH1 *ratioBck = static_cast<TH1*>(hBckTriggersNoPS->Clone("ratioBck"));
  ratioBck->Divide(hMBBNoPS);
  ratioBck->SetLineColor(kBlack);
  ratioBck->SetLineWidth(3);
  ratioBck->SetLineStyle(3);
  ratioBck->Scale(10);

  CanvasName =  LHCPeriod.Data() ; 
  CanvasName += "_RelativeTriggerContent"; 
  TCanvas *cRelativeTriggerContent = new TCanvas(CanvasName.Data(),"cRelativeTriggerContent",1200,900);
  cRelativeTriggerContent->SetTopMargin(0.05);
  cRelativeTriggerContent->SetRightMargin(0.01);
  cRelativeTriggerContent->SetGridy(1);
  //cRelativeTriggerContent->SetLogy(1);
  cRelativeTriggerContent->cd();

  ratioMBB->SetMaximum(1.0);
  ratioMBB->SetMinimum(0.05);
  ratioMBB->SetTitle("");
  ratioMBB->SetLabelSize(LabelSize);
  ratioMBB->LabelsOption("v");
  ratioMBB->GetYaxis()->SetTitle("Relative trigger content w/ and w/o Phys. Sel."); 

  ratioMBB->Draw("E");
  ratioMBBlow->Draw("ESAME");
  ratioMBBhigh->Draw("ESAME");
  ratioMBBNoPS->Draw("EPSAME");
  ratioMBBlowNoPS->Draw("EPSAME");
  ratioMBBhighNoPS->Draw("EPSAME");
  ratioBck->Draw("EPSAME");

  TLegend* legcRTC = new TLegend(0.4,0.25,0.70,0.45);
  legcRTC->SetHeader("Physics Selection");
  legcRTC->AddEntry(".","applied :","");
  legcRTC->AddEntry(ratioMBB,"MBB","l");
  legcRTC->AddEntry(ratioMBBlow,"MBB low mult.","l");
  legcRTC->AddEntry(ratioMBBhigh,"MBB high mult.","l");
  legcRTC->AddEntry(".","not applied :","");
  legcRTC->AddEntry(ratioMBBNoPS,"MBB","p");
  legcRTC->AddEntry(ratioMBBlowNoPS,"MBB low mult.","p");
  legcRTC->AddEntry(ratioMBBhighNoPS,"MBB high mult.","p");
  legcRTC->AddEntry(ratioBck,"Background (x10)","l");
  legcRTC->Draw("same");
  
  cRelativeTriggerContent->Print(OutFileNamePDF.Data());
  cRelativeTriggerContent->Write();


  //===================================================================================
  // new canvas with the relative ratio of multiplicity bin after physic selection
  TH1 *relratioMBBlow = static_cast<TH1*>(hMBBlowWithPS->Clone("relratioMBBlow"));
  relratioMBBlow->Divide(hMBBWithPS);
  relratioMBBlow->SetLineColor(kBlue);
  relratioMBBlow->SetLineWidth(2);
  relratioMBBlow->SetMarkerSize(0);
  TH1 *relratioMBBhigh = static_cast<TH1*>(hMBBhighWithPS->Clone("relratioMBBhigh"));
  relratioMBBhigh->Divide(hMBBWithPS);
  relratioMBBhigh->SetLineColor(kRed);
  relratioMBBhigh->SetLineWidth(2);
  relratioMBBhigh->SetMarkerSize(0);


  TH1 *relratioMBBlowNoPS = static_cast<TH1*>(hMBBlowNoPS->Clone("relratioMBBlowNoPS"));
  relratioMBBlowNoPS->Divide(hMBBNoPS);
  relratioMBBlowNoPS->SetLineColor(kBlue);
  relratioMBBlowNoPS->SetLineWidth(0);
  relratioMBBlowNoPS->SetLineStyle(2);
  relratioMBBlowNoPS->SetMarkerColor(kBlue);
  relratioMBBlowNoPS->SetMarkerStyle(24);
  relratioMBBlowNoPS->SetMarkerSize(1);
  TH1 *relratioMBBhighNoPS = static_cast<TH1*>(hMBBhighNoPS->Clone("relratioMBBhighNoPS"));
  relratioMBBhighNoPS->Divide(hMBBNoPS);
  relratioMBBhighNoPS->SetLineColor(kRed);
  relratioMBBhighNoPS->SetLineWidth(0);
  relratioMBBhighNoPS->SetLineStyle(2);
  relratioMBBhighNoPS->SetMarkerColor(kRed);
  relratioMBBhighNoPS->SetMarkerStyle(24);
  relratioMBBhighNoPS->SetMarkerSize(1);
  
  CanvasName =  LHCPeriod.Data() ; 
  CanvasName += "_CentralityPercentileCheck"; 
  TCanvas *cCentralityPercentileCheck = new TCanvas(CanvasName.Data(),"CentralityPercentileCheck",1200,900);
  cCentralityPercentileCheck->SetTopMargin(0.05);
  cCentralityPercentileCheck->SetRightMargin(0.01);
  cCentralityPercentileCheck->SetGridy(1);
  //cCentralityPercentileCheck->SetLogy(1);
  cCentralityPercentileCheck->cd();

  relratioMBBlow->Scale(0.8);
  relratioMBBhigh->Scale(0.8);
  relratioMBBlowNoPS->Scale(0.8);
  relratioMBBhighNoPS->Scale(0.8);

  relratioMBBlow->SetMaximum(0.3);
  relratioMBBlow->SetMinimum(0.01);
  relratioMBBlow->SetTitle("");
  relratioMBBlow->SetLabelSize(LabelSize);
  relratioMBBlow->GetYaxis()->SetTitle("Centrality percentile check"); 

  relratioMBBlow->Draw();
  relratioMBBhigh->Draw("ESAME");
  relratioMBBlowNoPS->Draw("EPSAME");
  relratioMBBhighNoPS->Draw("EPSAME");

  TLegend* legcCPC = new TLegend(0.12,0.15,0.42,0.30);
  legcCPC->SetHeader("Physics Selection");
  legcCPC->AddEntry(".","applied :","");
  legcCPC->AddEntry(relratioMBBlow,"MBB low mult.","l");
  legcCPC->AddEntry(relratioMBBhigh,"MBB high mult.","l");
  legcCPC->AddEntry(".","not applied :","");
  legcCPC->AddEntry(relratioMBBlowNoPS,"MBB low mult.","p");
  legcCPC->AddEntry(relratioMBBhighNoPS,"MBB high mult.","p");
  legcCPC->Draw("same");
  
  cCentralityPercentileCheck->Print(OutFileNamePDF.Data());
  cCentralityPercentileCheck->Write();


  //====================================================
  // Draw ratio of tracks over MBB versus run 
  TH1* hTrackerPerMBB= static_cast<TH1*>(hTrackerMBB->Clone("hTrackerPerMBB"));
  hTrackerPerMBB->Divide(hMBB);
  hTrackerPerMBB->SetLineWidth(2);
  hTrackerPerMBB->SetLineColor(kRed);

  TH1* hTriggerPerMBB= static_cast<TH1*>(hTriggerMBB->Clone("hTriggerPerMBB"));
  hTriggerPerMBB->Divide(hMBB);
  hTriggerPerMBB->SetLineWidth(2);
  hTriggerPerMBB->SetLineColor(kBlue);


  TH1* hMatchedPerMBB= static_cast<TH1*>(hMatchedMBB->Clone("hMatchedPerMBB"));
  hMatchedPerMBB->Divide(hMBB);
  hMatchedPerMBB->SetLineWidth(2);
  hMatchedPerMBB->SetLineColor(kViolet);


  TH1* hAllTracksPerMBB= static_cast<TH1*>(hAllTracksMBB->Clone("hAllTracksPerMBB"));
  hAllTracksPerMBB->Divide(hMBB);
  hAllTracksPerMBB->SetLineWidth(3);
  hAllTracksPerMBB->SetLineColor(kBlack);


  CanvasName =  LHCPeriod.Data() ; 
  CanvasName += "_RatioTrackTypesMBB"; 
  TCanvas *cRatioTrackTypesMBB = new TCanvas(CanvasName.Data(),"cRatioTrackTypesMBB",1200,900);
  cRatioTrackTypesMBB->SetRightMargin(0.01);
  cRatioTrackTypesMBB->SetGridy(1);
  //cRatioTrackTypesMBB->SetLogy(1);
  cRatioTrackTypesMBB->cd();

  hAllTracksPerMBB->SetTitle("Ratio (Number of Tracks)/MBB for [0-80]% centrality");
  hAllTracksPerMBB->SetMaximum(4);
  hAllTracksPerMBB->SetMinimum(0.1);
  hAllTracksPerMBB->SetLabelSize(LabelSize);
  hAllTracksPerMBB->Draw("E");

  hTrackerPerMBB->Draw("Esame");
  hMatchedPerMBB->Draw("Esame");
  hTriggerPerMBB->Draw("Esame");

  TLegend* legcTTMBB = new TLegend(0.15,0.50,0.35,0.70);
  legcTTMBB->AddEntry(hAllTracksPerMBB,"All tracks","l");
  legcTTMBB->AddEntry(hTrackerPerMBB,"Tracking (only) tracks","l");
  legcTTMBB->AddEntry(hMatchedPerMBB,"Matched tracks","l");
  legcTTMBB->AddEntry(hTriggerPerMBB,"Trigger (only) tracks","l");
  legcTTMBB->Draw("same");

  cRatioTrackTypesMBB->Print(OutFileNamePDF.Data());
  cRatioTrackTypesMBB->Write();

  //==========================================================================
  // Draw ratio of tracks over MBB versus runs and low and high multiplicities 

  TH1* hTrackerPerMBBlow= static_cast<TH1*>(hTrackerMBBlow->Clone("hTrackerPerMBBlow"));
  hTrackerPerMBBlow->Divide(hMBBlow);
  hTrackerPerMBBlow->SetLineWidth(2);
  hTrackerPerMBBlow->SetLineStyle(2);
  hTrackerPerMBBlow->SetLineColor(kRed);
  hTrackerPerMBBlow->SetMarkerColor(kRed);
  hTrackerPerMBBlow->SetMarkerSize(2);
  hTrackerPerMBBlow->SetMarkerStyle(27);
  TH1* hTrackerPerMBBhigh= static_cast<TH1*>(hTrackerMBBhigh->Clone("hTrackerPerMBBhigh"));
  hTrackerPerMBBhigh->Divide(hMBBhigh);
  hTrackerPerMBBhigh->SetLineWidth(2);
  hTrackerPerMBBhigh->SetLineStyle(2);
  hTrackerPerMBBhigh->SetLineColor(kRed);
  hTrackerPerMBBhigh->SetMarkerColor(kRed);
  hTrackerPerMBBhigh->SetMarkerSize(2);
  hTrackerPerMBBhigh->SetMarkerStyle(24);


  TH1* hTriggerPerMBBlow= static_cast<TH1*>(hTriggerMBBlow->Clone("hTriggerPerMBBlow"));
  hTriggerPerMBBlow->Divide(hMBBlow);
  hTriggerPerMBBlow->SetLineWidth(2);
  hTriggerPerMBBlow->SetLineStyle(2);
  hTriggerPerMBBlow->SetLineColor(kBlue);
  hTriggerPerMBBlow->SetMarkerColor(kBlue);
  hTriggerPerMBBlow->SetMarkerSize(2);
  hTriggerPerMBBlow->SetMarkerStyle(27);
  TH1* hTriggerPerMBBhigh= static_cast<TH1*>(hTriggerMBBhigh->Clone("hTriggerPerMBBhigh"));
  hTriggerPerMBBhigh->Divide(hMBBhigh);
  hTriggerPerMBBhigh->SetLineWidth(2);
  hTriggerPerMBBhigh->SetLineStyle(2);
  hTriggerPerMBBhigh->SetLineColor(kBlue);
  hTriggerPerMBBhigh->SetMarkerColor(kBlue);
  hTriggerPerMBBhigh->SetMarkerSize(2);
  hTriggerPerMBBhigh->SetMarkerStyle(24);


  TH1* hMatchedPerMBBlow= static_cast<TH1*>(hMatchedMBBlow->Clone("hMatchedPerMBBlow"));
  hMatchedPerMBBlow->Divide(hMBBlow);
  hMatchedPerMBBlow->SetLineWidth(2);
  hMatchedPerMBBlow->SetLineStyle(2);
  hMatchedPerMBBlow->SetLineColor(kViolet);
  hMatchedPerMBBlow->SetMarkerColor(kViolet);
  hMatchedPerMBBlow->SetMarkerSize(2);
  hMatchedPerMBBlow->SetMarkerStyle(27);
  TH1* hMatchedPerMBBhigh= static_cast<TH1*>(hMatchedMBBhigh->Clone("hMatchedPerMBBhigh"));
  hMatchedPerMBBhigh->Divide(hMBBhigh);
  hMatchedPerMBBhigh->SetLineWidth(2);
  hMatchedPerMBBhigh->SetLineStyle(2);
  hMatchedPerMBBhigh->SetLineColor(kViolet);
  hMatchedPerMBBhigh->SetMarkerColor(kViolet);
  hMatchedPerMBBhigh->SetMarkerSize(2);
  hMatchedPerMBBhigh->SetMarkerStyle(24);


  TH1* hAllTracksPerMBBlow= static_cast<TH1*>(hAllTracksMBBlow->Clone("hAllTracksPerMBBlow"));
  hAllTracksPerMBBlow->Divide(hMBBlow);
  hAllTracksPerMBBlow->SetLineWidth(2);
  hAllTracksPerMBBlow->SetLineStyle(2);
  hAllTracksPerMBBlow->SetLineColor(kBlack);
  hAllTracksPerMBBlow->SetMarkerColor(kBlack);
  hAllTracksPerMBBlow->SetMarkerSize(2);
  hAllTracksPerMBBlow->SetMarkerStyle(27);
  TH1* hAllTracksPerMBBhigh= static_cast<TH1*>(hAllTracksMBBhigh->Clone("hAllTracksPerMBBhigh"));
  hAllTracksPerMBBhigh->Divide(hMBBhigh);
  hAllTracksPerMBBhigh->SetLineWidth(2);
  hAllTracksPerMBBhigh->SetLineStyle(2);
  hAllTracksPerMBBhigh->SetLineColor(kBlack);
  hAllTracksPerMBBhigh->SetMarkerColor(kBlack);
  hAllTracksPerMBBhigh->SetMarkerSize(2);
  hAllTracksPerMBBhigh->SetMarkerStyle(24);


  CanvasName =  LHCPeriod.Data() ; 
  CanvasName += "_RatioTrackTypesMBBVsMult"; 
  TCanvas *cRatioTrackTypesMBBVsMult = new TCanvas(CanvasName.Data(),"cRatioTrackTypesMBBVsMult",1200,900);
  cRatioTrackTypesMBBVsMult->SetRightMargin(0.01);
  cRatioTrackTypesMBBVsMult->SetGridy(1);
  //cRatioTrackTypesMBBVsMult->SetLogy(1);


  cRatioTrackTypesMBBVsMult->Divide(1,2);


  cRatioTrackTypesMBBVsMult->cd(1);

  hAllTracksPerMBBlow->SetTitle("Ratio (Number of Tracks)/MBB low mult");
  hAllTracksPerMBBlow->SetMaximum(0.30);
  hAllTracksPerMBBlow->SetMinimum(0.001);
  hAllTracksPerMBBlow->SetLabelSize(LabelSize);
  hAllTracksPerMBBlow->Draw("E");
  hTrackerPerMBBlow->Draw("Esame");
  hTriggerPerMBBlow->Draw("Esame");
  hMatchedPerMBBlow->Draw("Esame");

  TLegend* legcTTMBBlow = new TLegend(0.40,0.8,0.60,1.);
  legcTTMBBlow->AddEntry(hAllTracksPerMBBlow,"All tracks","p");
  legcTTMBBlow->AddEntry(hTrackerPerMBBlow,"Tracking (only) tracks","p");
  legcTTMBBlow->AddEntry(hMatchedPerMBBlow,"Matched tracks","p");
  legcTTMBBlow->AddEntry(hTriggerPerMBBlow,"Trigger (only) tracks","p");
  legcTTMBBlow->Draw("same");
 

  cRatioTrackTypesMBBVsMult->cd(2);
  hAllTracksPerMBBhigh->SetTitle("Ratio (Number of Tracks)/MBB high mult");
  hAllTracksPerMBBhigh->SetMaximum(15);
  hAllTracksPerMBBhigh->SetMinimum(0.01);
  hAllTracksPerMBBhigh->SetLabelSize(LabelSize);
  hAllTracksPerMBBhigh->Draw("E");
  hTrackerPerMBBhigh->Draw("Esame");
  hTriggerPerMBBhigh->Draw("Esame");
  hMatchedPerMBBhigh->Draw("Esame");



  TLegend* legcTTMBBhigh = new TLegend(0.40,0.8,0.60,1.);
  legcTTMBBhigh->AddEntry(hAllTracksPerMBBhigh,"All tracks","p");
  legcTTMBBhigh->AddEntry(hTrackerPerMBBhigh,"Tracking (only) tracks","p");
  legcTTMBBhigh->AddEntry(hMatchedPerMBBhigh,"Matched tracks","p");
  legcTTMBBhigh->AddEntry(hTriggerPerMBBhigh,"Trigger (only) tracks","p");
  legcTTMBBhigh->Draw("same");
    

  //   cRatioTrackTypesMBBVsMult->cd();
  //   TPaveText pt1(0.2,0.2,0.7,0.9);
  //   TText* t1=pt1.AddText("All tracks (tracking-only + trigger-only + matched)");
  //   t1->SetTextColor(kBlack);
  //   pt1.Draw();
  //   break;

  cRatioTrackTypesMBBVsMult->Print(OutFileNamePDF.Data());
  cRatioTrackTypesMBBVsMult->Write();


  //===========================================
  // draw ratio of track (summed, i.e all trigger=(matched+trigger-only)) for MB  versus run for low mult
  CanvasName =  LHCPeriod.Data() ; 
  CanvasName += "_TrackMultBBlow";   
  TCanvas* cTrackMultBBlow = new TCanvas("cTrackMultBBlow","cTrackMultBBlow",1200,900);
  // must be applied on the pads created by  divide
  cTrackMultBBlow->Divide(1,2);
  cTrackMultBBlow->cd(1);

  TH1* hSumTriggerOverMBBlow = static_cast<TH1*>(hTriggerMBBlow->Clone("hSumTriggerOverMBBlow"));
  hSumTriggerOverMBBlow->Add(hMatchedMBBlow);
  hSumTriggerOverMBBlow->Divide(hMBBlow);
  hSumTriggerOverMBBlow->SetTitle("Sum of trigger tracks (matched+trigger-only) in MBB events / # MBB events in low mult.");
  //hSumTriggerOverMBBlow->LabelsOption("u");
  hSumTriggerOverMBBlow->SetLabelSize(LabelSize);
  hSumTriggerOverMBBlow->SetLineWidth(2);
  hSumTriggerOverMBBlow->SetLineColor(kBlue);
  hSumTriggerOverMBBlow->Draw("e");
  cTrackMultBBlow->cd(2);
  TH1* hSumTrackerOverMBBlow = static_cast<TH1*>(hTrackerMBBlow->Clone("hSumTrackerOverMBBlow"));
  hSumTrackerOverMBBlow->Add(hMatchedMBBlow);
  hSumTrackerOverMBBlow->Divide(hMBBlow);
  hSumTrackerOverMBBlow->SetTitle("Sum tracker tracks (matched+tracker-only) in MBB events / # MBB events in low mult.");
  //hSumTrackerOverMBBlow->LabelsOption("u");
  hSumTrackerOverMBBlow->SetLabelSize(LabelSize);
  hSumTrackerOverMBBlow->SetLineWidth(2);
  hSumTrackerOverMBBlow->SetLineColor(kBlue);
  hSumTrackerOverMBBlow->Draw("e");


  cTrackMultBBlow->Print(OutFileNamePDF.Data());
  cTrackMultBBlow->Write();



  //===========================================
  // draw ratio of track (summed, i.e all trigger=(matched+trigger-only)) for MB  versus run for high mult
  CanvasName =  LHCPeriod.Data() ; 
  CanvasName += "_TrackMultBBhigh";   
  TCanvas* cTrackMultBBhigh = new TCanvas("cTrackMultBBhigh","cTrackMultBBhigh",1200,900);
  // must be applied on the pads created by  divide
  cTrackMultBBhigh->Divide(1,2);
  cTrackMultBBhigh->cd(1);

  TH1* hSumTriggerOverMBBhigh = static_cast<TH1*>(hTriggerMBBhigh->Clone("hSumTriggerOverMBBhigh"));
  hSumTriggerOverMBBhigh->Add(hMatchedMBBhigh);
  hSumTriggerOverMBBhigh->Divide(hMBBhigh);
  hSumTriggerOverMBBhigh->SetTitle("Sum of trigger tracks (matched+trigger-only) in MBB events / # MBB events in high mult.");
  //hSumTriggerOverMBBhigh->LabelsOption("u");
  hSumTriggerOverMBBhigh->SetLabelSize(LabelSize);
  hSumTriggerOverMBBhigh->SetLineWidth(2);
  hSumTriggerOverMBBhigh->SetLineColor(kBlue);
  hSumTriggerOverMBBhigh->Draw("e");
  cTrackMultBBhigh->cd(2);
  TH1* hSumTrackerOverMBBhigh = static_cast<TH1*>(hTrackerMBBhigh->Clone("hSumTrackerOverMBBhigh"));
  hSumTrackerOverMBBhigh->Add(hMatchedMBBhigh);
  hSumTrackerOverMBBhigh->Divide(hMBBhigh);
  hSumTrackerOverMBBhigh->SetTitle("Sum tracker tracks (matched+tracker-only) in MBB events / # MBB events in high mult.");
  //hSumTrackerOverMBBhigh->LabelsOption("u");
  hSumTrackerOverMBBhigh->SetLabelSize(LabelSize);
  hSumTrackerOverMBBhigh->SetLineWidth(2);
  hSumTrackerOverMBBhigh->SetLineColor(kBlue);
  hSumTrackerOverMBBhigh->Draw("e");


  cTrackMultBBhigh->Print(OutFileNamePDF.Data());
  cTrackMultBBhigh->Write();

  
  
  //===========================================
  // draw mixed ratio of track over track versus run for MBB in low mult 
  CanvasName =  LHCPeriod.Data() ; 
  CanvasName += "_RatioTrackMBBlow";   
  TCanvas* cRatioTrackMBBlow = new TCanvas("cRatioTrackMBBlow","cRatioTrackMBBlow",1200,900);
  cRatioTrackMBBlow->Divide(1,3);

  cRatioTrackMBBlow->cd(1);
  TH1* hTrackerOverTriggerMBBlow = static_cast<TH1*>(hTrackerMBBlow->Clone("hTrackerOverTriggerMBBlow"));
  hTrackerOverTriggerMBBlow->Divide(hTriggerMBBlow);
  hTrackerOverTriggerMBBlow->SetTitle("# tracker tracks / # trigger tracks in MBB low mult.");
  //hTrackerOverTriggerMBBlow->LabelsOption("u");
  hTrackerOverTriggerMBBlow->SetLabelSize(LabelSize);
  hTrackerOverTriggerMBBlow->SetLineWidth(2);
  hTrackerOverTriggerMBBlow->SetLineColor(kBlue);
  hTrackerOverTriggerMBBlow->Draw("e");

  cRatioTrackMBBlow->cd(2);
  TH1* hMatchedOverTriggerMBBlow = static_cast<TH1*>(hMatchedMBBlow->Clone("hMatchedOverTriggerMBBlow"));
  hMatchedOverTriggerMBBlow->Divide(hTriggerMBBlow);
  hMatchedOverTriggerMBBlow->SetTitle("# matched tracks / # trigger tracks in MBB low mult.");
  //hMatchedOverTriggerMBBlow->LabelsOption("u");
  hMatchedOverTriggerMBBlow->SetLabelSize(LabelSize);
  hMatchedOverTriggerMBBlow->SetLineWidth(2);
  hMatchedOverTriggerMBBlow->SetLineColor(kBlue);
  hMatchedOverTriggerMBBlow->Draw("e");

  cRatioTrackMBBlow->cd(3);
  TH1* hMatchedOverTrackerMBBlow = static_cast<TH1*>(hMatchedMBBlow->Clone("hMatchedOverTrackerMBBlow"));
  hMatchedOverTrackerMBBlow->Divide(hTrackerMBBlow);
  hMatchedOverTrackerMBBlow->SetTitle("# matched tracks / # tracker tracks in MBB low mult.");
  //hMatchedOverTrackerMBBlow->LabelsOption("u");
  hMatchedOverTrackerMBBlow->SetLabelSize(LabelSize);
  hMatchedOverTrackerMBBlow->SetLineWidth(2);
  hMatchedOverTrackerMBBlow->SetLineColor(kBlue);
  hMatchedOverTrackerMBBlow->Draw("e");
  
  cRatioTrackMBBlow->Print(OutFileNamePDF.Data());
  cRatioTrackMBBlow->Write();

  
  //===========================================
  // draw mixed ratio of track over track versus run for MBB in high mult 
  CanvasName =  LHCPeriod.Data() ; 
  CanvasName += "_RatioTrackMBBhigh";   
  TCanvas* cRatioTrackMBBhigh = new TCanvas("cRatioTrackMBBhigh","cRatioTrackMBBhigh",1200,900);
  cRatioTrackMBBhigh->Divide(1,3);

  cRatioTrackMBBhigh->cd(1);
  TH1* hTrackerOverTriggerMBBhigh = static_cast<TH1*>(hTrackerMBBhigh->Clone("hTrackerOverTriggerMBBhigh"));
  hTrackerOverTriggerMBBhigh->Divide(hTriggerMBBhigh);
  hTrackerOverTriggerMBBhigh->SetTitle("# tracker tracks / # trigger tracks in MBB high mult.");
  //hTrackerOverTriggerMBBhigh->LabelsOption("u");
  hTrackerOverTriggerMBBhigh->SetLabelSize(LabelSize);
  hTrackerOverTriggerMBBhigh->SetLineWidth(2);
  hTrackerOverTriggerMBBhigh->SetLineColor(kBlue);
  hTrackerOverTriggerMBBhigh->Draw("e");

  cRatioTrackMBBhigh->cd(2);
  TH1* hMatchedOverTriggerMBBhigh = static_cast<TH1*>(hMatchedMBBhigh->Clone("hMatchedOverTriggerMBBhigh"));
  hMatchedOverTriggerMBBhigh->Divide(hTriggerMBBhigh);
  hMatchedOverTriggerMBBhigh->SetTitle("# matched tracks / # trigger tracks in MBB high mult.");
  //hMatchedOverTriggerMBBhigh->LabelsOption("u");
  hMatchedOverTriggerMBBhigh->SetLabelSize(LabelSize);
  hMatchedOverTriggerMBBhigh->SetLineWidth(2);
  hMatchedOverTriggerMBBhigh->SetLineColor(kBlue);
  hMatchedOverTriggerMBBhigh->Draw("e");

  cRatioTrackMBBhigh->cd(3);
  TH1* hMatchedOverTrackerMBBhigh = static_cast<TH1*>(hMatchedMBBhigh->Clone("hMatchedOverTrackerMBBhigh"));
  hMatchedOverTrackerMBBhigh->Divide(hTrackerMBBhigh);
  hMatchedOverTrackerMBBhigh->SetTitle("# matched tracks / # tracker tracks in MBB high mult.");
  //hMatchedOverTrackerMBBhigh->LabelsOption("u");
  hMatchedOverTrackerMBBhigh->SetLabelSize(LabelSize);
  hMatchedOverTrackerMBBhigh->SetLineWidth(2);
  hMatchedOverTrackerMBBhigh->SetLineColor(kBlue);
  hMatchedOverTrackerMBBhigh->Draw("e");
  
  cRatioTrackMBBhigh->Print(OutFileNamePDF.Data());
  cRatioTrackMBBhigh->Write();





//   //===========================================
//   // Draw ratio of tracks over CMUS1B versus run
//   TH1* hTrackerPerCMUS1B= static_cast<TH1*>(hTrackerCMUS1B->Clone("hTrackerPerCMUS1B"));
//   hTrackerPerCMUS1B->Divide(hCMUS1B);
//   hTrackerPerCMUS1B->SetLineWidth(2);
//   hTrackerPerCMUS1B->SetLineColor(kRed);

//   TH1* hTriggerPerCMUS1B= static_cast<TH1*>(hTriggerCMUS1B->Clone("hTriggerPerCMUS1B"));
//   hTriggerPerCMUS1B->Divide(hCMUS1B);
//   hTriggerPerCMUS1B->SetLineWidth(2);
//   hTriggerPerCMUS1B->SetLineColor(kBlue);

//   TH1* hMatchedPerCMUS1B= static_cast<TH1*>(hMatchedCMUS1B->Clone("hMatchedPerCMUS1B"));
//   hMatchedPerCMUS1B->Divide(hCMUS1B);
//   hMatchedPerCMUS1B->SetLineWidth(2);
//   hMatchedPerCMUS1B->SetLineColor(kViolet);


//   TH1* hAllTracksPerCMUS1B= static_cast<TH1*>(hAllTracksCMUS1B->Clone("hAllTracksPerCMUS1B"));
//   hAllTracksPerCMUS1B->Divide(hCMUS1B);
//   hAllTracksPerCMUS1B->SetLineWidth(3);
//   hAllTracksPerCMUS1B->SetLineColor(kBlack);


//   CanvasName =  LHCPeriod.Data() ; 
//   CanvasName += "_RatioTrackTypesCMUS1B"; 
//   TCanvas *cRatioTrackTypesCMUS1B = new TCanvas(CanvasName.Data(),"cRatioTrackTypesCMUS1B",1200,900);
//   cRatioTrackTypesCMUS1B->SetRightMargin(0.01);
//   cRatioTrackTypesCMUS1B->SetGridy(1);
//   cRatioTrackTypesCMUS1B->cd();

//   hAllTracksPerCMUS1B->SetTitle("Ratio (Number of Tracks)/CMUS1B");
//   hAllTracksPerCMUS1B->SetMinimum(0.01);
//   hAllTracksPerCMUS1B->SetLabelSize(LabelSize);
//   hAllTracksPerCMUS1B->Draw("E");

//   hTrackerPerCMUS1B->Draw("Esame");
//   hMatchedPerCMUS1B->Draw("Esame");
//   hTriggerPerCMUS1B->Draw("Esame");

//   TLegend* legcTTCMUS1B = new TLegend(0.75,0.55,0.90,0.75);
//   legcTTCMUS1B->AddEntry(hAllTracksPerCMUS1B,"All tracks","l");
//   legcTTCMUS1B->AddEntry(hTrackerPerCMUS1B,"Tracking (only) tracks","l");
//   legcTTCMUS1B->AddEntry(hMatchedPerCMUS1B,"Matched tracks","l");
//   legcTTCMUS1B->AddEntry(hTriggerPerCMUS1B,"Trigger (only) tracks","l");
//   legcTTCMUS1B->Draw("same");

//   cRatioTrackTypesCMUS1B->Print(OutFileNamePDF.Data());
//   cRatioTrackTypesCMUS1B->Write();


//   //===========================================
//   // draw ratio of track (summed, i.e all trigger=(matched+trigger-only)) over CMUS1B versus run
//   TCanvas* cTrackMultCMUS1B = new TCanvas("cTrackMultCMUS1B","cTrackMultCMUS1B",1200,900);
//   cTrackMultCMUS1B->Divide(1,2);
//   cTrackMultCMUS1B->cd(1);
//   TH1* hSumTriggerOverCMUS1B = static_cast<TH1*>(hTriggerCMUS1B->Clone("hSumTriggerOverCMUS1B"));
//   hSumTriggerOverCMUS1B->Add(hMatchedCMUS1B);
//   hSumTriggerOverCMUS1B->Divide(hCMUS1B);
//   hSumTriggerOverCMUS1B->SetTitle("Sum of trigger tracks (matched+trigger-only) in CMUS1B events / # CMUS1B events");
//   //hSumTriggerOverCMUS1B->LabelsOption("u");
//   hSumTriggerOverCMUS1B->SetLabelSize(LabelSize);
//   hSumTriggerOverCMUS1B->SetLineWidth(2);
//   hSumTriggerOverCMUS1B->SetLineColor(kRed);
//   hSumTriggerOverCMUS1B->Draw("e");

//   cTrackMultCMUS1B->cd(2);
//   TH1* hSumTrackerOverCMUS1B = static_cast<TH1*>(hTrackerCMUS1B->Clone("hSumTrackerOverCMUS1B"));
//   hSumTrackerOverCMUS1B->Add(hMatchedCMUS1B);
//   hSumTrackerOverCMUS1B->Divide(hCMUS1B);
//   hSumTrackerOverCMUS1B->SetTitle("Sum of tracker tracks (matched+tracker-only) in CMUS1B events / # CMUS1B events");
//   //hSumTrackerOverCMUS1B->LabelsOption("u");
//   hSumTrackerOverCMUS1B->SetLabelSize(LabelSize);
//   hSumTrackerOverCMUS1B->SetLineWidth(2);
//   hSumTrackerOverCMUS1B->SetLineColor(kRed);
//   hSumTrackerOverCMUS1B->Draw("e");

//   cTrackMultCMUS1B->Print(OutFileNamePDF.Data());
//   cTrackMultCMUS1B->Write();
  
//   //===========================================
//   // draw mixed ratio of track over track versus run for CMUS1B
//   TCanvas* cRatioTrackCMUS1B = new TCanvas("cRatioTrackCMUS1B","cRatioTrackCMUS1B",1200,900);
//   cRatioTrackCMUS1B->Divide(1,3);

//   cRatioTrackCMUS1B->cd(1);
//   TH1* hTrackerOverTriggerCMUS1B = static_cast<TH1*>(hTrackerCMUS1B->Clone("hTrackerOverTriggerCMUS1B"));
//   hTrackerOverTriggerCMUS1B->Divide(hTriggerCMUS1B);
//   hTrackerOverTriggerCMUS1B->SetTitle("# tracker tracks / # trigger tracks in CMUS1B");
//   //hTrackerOverTriggerCMUS1B->LabelsOption("u");
//   hTrackerOverTriggerCMUS1B->SetLabelSize(LabelSize);
//   hTrackerOverTriggerCMUS1B->SetLineWidth(2);
//   hTrackerOverTriggerCMUS1B->SetLineColor(kRed);
//   hTrackerOverTriggerCMUS1B->Draw("e");

//   cRatioTrackCMUS1B->cd(2);
//   TH1* hMatchedOverTriggerCMUS1B = static_cast<TH1*>(hMatchedCMUS1B->Clone("hMatchedOverTriggerCMUS1B"));
//   hMatchedOverTriggerCMUS1B->Divide(hTriggerCMUS1B);
//   hMatchedOverTriggerCMUS1B->SetTitle("# matched tracks / # trigger tracks in CMUS1B");
//   //hMatchedOverTriggerCMUS1B->LabelsOption("u");
//   hMatchedOverTriggerCMUS1B->SetLabelSize(LabelSize);
//   hMatchedOverTriggerCMUS1B->SetLineWidth(2);
//   hMatchedOverTriggerCMUS1B->SetLineColor(kRed);
//   hMatchedOverTriggerCMUS1B->Draw("e");

//   cRatioTrackCMUS1B->cd(3);
//   TH1* hMatchedOverTrackerCMUS1B = static_cast<TH1*>(hMatchedCMUS1B->Clone("hMatchedOverTrackerCMUS1B"));
//   hMatchedOverTrackerCMUS1B->Divide(hTrackerCMUS1B);
//   hMatchedOverTrackerCMUS1B->SetTitle("# matched tracks / # tracker tracks in CMUS1B");
//   //hMatchedOverTrackerCMUS1B->LabelsOption("u");
//   hMatchedOverTrackerCMUS1B->SetLabelSize(LabelSize);
//   hMatchedOverTrackerCMUS1B->SetLineWidth(2);
//   hMatchedOverTrackerCMUS1B->SetLineColor(kRed);
//   hMatchedOverTrackerCMUS1B->Draw("e");
  
//   cRatioTrackCMUS1B->Print(OutFileNamePDF.Data());
//   cRatioTrackCMUS1B->Write();

   //==================================================
   // Draw matched tracks charge asymetry in low mult for MBB  triggers  
   TH1 *hDiffMatchedMBBlow= static_cast<TH1*>(hPosMatchedMBBlow->Clone("hDiffMatchedMBBlow"));
   hDiffMatchedMBBlow->Add(hNegMatchedMBBlow,-1);
   hDiffMatchedMBBlow->Sumw2();
  
   TH1 *hAsymMatchedMBBlow= static_cast<TH1*>(hDiffMatchedMBBlow->Clone("hAsymMatchedMBBlow"));
   hAsymMatchedMBBlow->Divide(hAllMatchedMBBlow);
   hAsymMatchedMBBlow->SetLineColor(kBlue);
   hAsymMatchedMBBlow->SetLineWidth(2);
   hAsymMatchedMBBlow->SetMinimum(-0.2);
   hAsymMatchedMBBlow->SetMaximum(0.2);
   hAsymMatchedMBBlow->SetLabelSize(LabelSize);
   hAsymMatchedMBBlow->SetTitle("Matched tracks charge asymmetry in Physics Selected events for MBB low (blue) and high (red) mult.");

   TH1 *hDiffMatchedMBBhigh= static_cast<TH1*>(hPosMatchedMBBhigh->Clone("hDiffMatchedMBBhigh"));
   hDiffMatchedMBBhigh->Add(hNegMatchedMBBhigh,-1);
   hDiffMatchedMBBhigh->Sumw2();
  
   TH1 *hAsymMatchedMBBhigh= static_cast<TH1*>(hDiffMatchedMBBhigh->Clone("hAsymMatchedMBBhigh"));
   hAsymMatchedMBBhigh->Divide(hAllMatchedMBBhigh);
   hAsymMatchedMBBhigh->SetLineColor(kRed);
   hAsymMatchedMBBhigh->SetLineWidth(2);
   hAsymMatchedMBBhigh->SetMinimum(-0.2);
   hAsymMatchedMBBhigh->SetMaximum(0.2);
   hAsymMatchedMBBhigh->SetLabelSize(LabelSize);
   hAsymMatchedMBBhigh->SetTitle("Matched tracks asymetry in Physics Selected events");
  

   CanvasName =  LHCPeriod.Data() ; 
   CanvasName += "_AsymMatched"; 
   TCanvas *cAsymMatched = new TCanvas(CanvasName.Data(),"cAsymMatched",1200,900);
   cAsymMatched->SetRightMargin(0.01);
   cAsymMatched->SetGridy(1);
   cAsymMatched->cd();
   hAsymMatchedMBBlow->GetYaxis()->SetTitle("Charge asymmetry");  
   hAsymMatchedMBBlow->Draw("EH");
   hAsymMatchedMBBhigh->Draw("EHsame");
   
   TLegend* legcAMT = new TLegend(0.60,0.25,0.98,0.45);
   legcAMT->SetHeader("Charge asymmetry  of matched track per MBB (include Vtx, #eta  and R_{Abs} cuts) ");
   legcAMT->AddEntry(hAsymMatchedMBBlow," Low mult. events","l");
   legcAMT->AddEntry(hAsymMatchedMBBhigh," High mult. events ","l");
   legcAMT->Draw("same");

   cAsymMatched->Print(OutFileNamePDF.Data());
   cAsymMatched->Write();


  //=========================================================
  // Draw low/high pt tracks in acceptance  per MBB  all mult 
  TH1* hMatchedLowPtPerMBB = static_cast<TH1*> (hMatchedLowPtMBB->Clone("hMatchedLowPtPerMBB"));
  hMatchedLowPtPerMBB->Sumw2();
  hMatchedLowPtPerMBB->Divide(hMBB);
  hMatchedLowPtPerMBB->SetLineWidth(2);
  hMatchedLowPtPerMBB->SetLineColor(kBlue);

  TH1* hMatchedHighPtPerMBB = static_cast<TH1*> (hMatchedHighPtMBB->Clone("hMatchedHighPtPerMBB"));
  hMatchedHighPtPerMBB->Sumw2();
  hMatchedHighPtPerMBB->Divide(hMBB);
  hMatchedHighPtPerMBB->SetLineWidth(2);
  hMatchedHighPtPerMBB->SetLineColor(kRed);

  TH1* hMatchedLowPtPerMBBlow = static_cast<TH1*> (hMatchedLowPtMBBlow->Clone("hMatchedLowPtPerMBBlow"));
  hMatchedLowPtPerMBBlow->Sumw2();
  hMatchedLowPtPerMBBlow->Divide(hMBBlow);
  hMatchedLowPtPerMBBlow->SetLineWidth(2);
  hMatchedLowPtPerMBBlow->SetMarkerSize(2);
  hMatchedLowPtPerMBBlow->SetMarkerStyle(27);
  hMatchedLowPtPerMBBlow->SetMarkerColor(kBlue);
  hMatchedLowPtPerMBBlow->SetLineWidth(2);
  hMatchedLowPtPerMBBlow->SetLineColor(kBlue);

  TH1* hMatchedHighPtPerMBBlow = static_cast<TH1*> (hMatchedHighPtMBBlow->Clone("hMatchedHighPtPerMBBlow"));
  hMatchedHighPtPerMBBlow->Sumw2();
  hMatchedHighPtPerMBBlow->Divide(hMBB);
  hMatchedHighPtPerMBBlow->SetLineWidth(2);
  hMatchedHighPtPerMBBlow->SetLineColor(kRed);
  hMatchedHighPtPerMBBlow->SetMarkerSize(2);
  hMatchedHighPtPerMBBlow->SetMarkerStyle(27);
  hMatchedHighPtPerMBBlow->SetMarkerColor(kRed);

  TH1* hMatchedLowPtPerMBBhigh = static_cast<TH1*> (hMatchedLowPtMBBhigh->Clone("hMatchedLowPtPerMBBhigh"));
  hMatchedLowPtPerMBBhigh->Sumw2();
  hMatchedLowPtPerMBBhigh->Divide(hMBBhigh);
  hMatchedLowPtPerMBBhigh->SetLineWidth(2);
  hMatchedLowPtPerMBBhigh->SetLineColor(kBlue);
  hMatchedLowPtPerMBBhigh->SetMarkerSize(2);
  hMatchedLowPtPerMBBhigh->SetMarkerStyle(24);
  hMatchedLowPtPerMBBhigh->SetMarkerColor(kBlue);


  TH1* hMatchedHighPtPerMBBhigh = static_cast<TH1*> (hMatchedHighPtMBBhigh->Clone("hMatchedHighPtPerMBBhigh"));
  hMatchedHighPtPerMBBhigh->Sumw2();
  hMatchedHighPtPerMBBhigh->Divide(hMBB);
  hMatchedHighPtPerMBBhigh->SetLineWidth(2);
  hMatchedHighPtPerMBBhigh->SetLineColor(kRed);
  hMatchedHighPtPerMBBhigh->SetMarkerSize(2);
  hMatchedHighPtPerMBBhigh->SetMarkerStyle(24);
  hMatchedHighPtPerMBBhigh->SetMarkerColor(kRed);


  CanvasName =  LHCPeriod.Data() ; 
  CanvasName += "_HighPtMuons"; 
  TCanvas *cHighPtMuons = new TCanvas(CanvasName.Data(),"cHighPtMuons",1200,900);
  cHighPtMuons->SetTopMargin(0.05);
  cHighPtMuons->SetRightMargin(0.01);
  //cHighPtMuons->SetLogy(1);
  cHighPtMuons->SetGridy(1);
  cHighPtMuons->cd();

  hMatchedLowPtPerMBB->SetTitle("");
  hMatchedLowPtPerMBB->GetYaxis()->SetTitle("Ratio per MBB");
  hMatchedLowPtPerMBB->SetMaximum(0.13);
  hMatchedLowPtPerMBB->SetMinimum(0.02);
  hMatchedLowPtPerMBB->SetLabelSize(LabelSize);

  hMatchedLowPtPerMBB->Draw("E");
  hMatchedHighPtPerMBB->Scale(5);
  hMatchedHighPtPerMBB->Draw("Esame");
    
  TLegend* legcHPM = new TLegend(0.60,0.50,0.98,0.65);
  legcHPM->SetHeader("Number of matched track per MBB (include Vtx, #eta  and R_{Abs} cuts) ");
  legcHPM->AddEntry(hMatchedLowPtPerMBB," p_{T} > 1 GeV/c ","l");
  legcHPM->AddEntry(hMatchedHighPtPerMBB," (x5) p_{T} >  2 GeV/c ","l");
  legcHPM->Draw("same");

  cHighPtMuons->Print(OutFileNamePDF.Data());
  cHighPtMuons->Write(); 


  //=====================================================================
  // Draw low/high pt tracks in acceptance per MBB  in low and high  mult 
  CanvasName =  LHCPeriod.Data() ; 
  CanvasName += "_HighPtMuonsVsMult"; 
  TCanvas *cHighPtMuonsVsMult = new TCanvas(CanvasName.Data(),"cHighPtMuonsVsMult",1200,900);
  cHighPtMuonsVsMult->SetTopMargin(0.05);
  cHighPtMuonsVsMult->SetRightMargin(0.01);
  cHighPtMuonsVsMult->SetLogy(1);
  cHighPtMuonsVsMult->SetGridy(1);
  cHighPtMuonsVsMult->Divide(1,2);


  cHighPtMuonsVsMult->cd(1);

  hMatchedLowPtPerMBBlow->SetTitle("");
  hMatchedLowPtPerMBBlow->GetYaxis()->SetTitle("Ratio per MBB in low mult");
  hMatchedLowPtPerMBBlow->SetMaximum(0.015);
  hMatchedLowPtPerMBBlow->SetMinimum(0.000001);
  hMatchedLowPtPerMBBlow->SetLabelSize(LabelSize);

  hMatchedLowPtPerMBBlow->Draw("E");
  hMatchedHighPtPerMBBlow->Scale(10);
  hMatchedHighPtPerMBBlow->Draw("Esame");
  TLegend* legcHPMlow = new TLegend(0.60,0.35,0.98,0.5);
  legcHPMlow->SetHeader("Number of matched track per MBB (include Vtx, #eta  and R_{Abs} cuts) ");
  legcHPMlow->AddEntry(hMatchedLowPtPerMBBlow," p_{T} > 1 GeV/c in low mult","p");
  legcHPMlow->AddEntry(hMatchedHighPtPerMBBlow,"(x10) p_{T} >  2 GeV/c in low mult","p");
  legcHPMlow->Draw("same");


  cHighPtMuonsVsMult->cd(2);
  hMatchedLowPtPerMBBhigh->SetTitle("");
  hMatchedLowPtPerMBBhigh->GetYaxis()->SetTitle("Ratio per MBB in high mult");
  hMatchedLowPtPerMBBhigh->SetMinimum(0.005);
  hMatchedLowPtPerMBBhigh->Draw("E");
  hMatchedHighPtPerMBBhigh->Scale(10);
  hMatchedHighPtPerMBBhigh->Draw("Esame");
    
  TLegend* legcHPMhigh = new TLegend(0.60,0.4,0.98,0.55);
  legcHPMhigh->SetHeader("Number of matched track per MBB (include Vtx, #eta  and R_{Abs} cuts) ");
  legcHPMhigh->AddEntry(hMatchedLowPtPerMBBhigh," p_{T} > 1 GeV/c in high mult","p");
  legcHPMhigh->AddEntry(hMatchedHighPtPerMBBhigh,"(x10)  p_{T} >  2 GeV/c in high mult","p");
  legcHPMhigh->Draw("same");
  
  cHighPtMuonsVsMult->Print(OutFileNamePDF.Data());
  cHighPtMuonsVsMult->Write(); 

  // close merged file
  globalFile->Close();


  //--------------------------------------------- //
  //        monitor quantities run per run        //
  //--------------------------------------------- //
  

  // Are the runs stored locally or in alien?
  Int_t isAlienFile = 0;
  if(alienBaseDir.Contains("alien:")){
    isAlienFile = 1;
    alienBaseDir.ReplaceAll("alien://","");
  }

  /*   if ( ! alienBaseDir.Contains("alien:")){
       Info("PlotMuonQA","Working locally, stopping here");
       rootFileOut->Close();
       return;
       }
  */ 
  //   else (!TGrid::Connect("alien://")) {
  //     Error("PlotMuonQA","cannot connect to grid. It's really needed here");
  //     c1->Print(OutFileNamePDF_close.Data());
  //     rootFileOut->Close();
  //     return;
  //   }

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
  
  // loop over runs
  for ( Int_t irun=0; irun<runs.GetEntriesFast(); irun++ ) { //loop over runs 
    
    TString run = ((TObjString*)runs.UncheckedAt(irun))->GetString();
    // get the file (or list of files) to be analyzed
    TString command;
    TGridResult *res = 0;
    TObjString *objs = 0;
  
    if(isAlienFile){
      // get the file (or list of files) to be analyzed
      command = Form("find %s/ %s/%s", alienBaseDir.Data(), run.Data(), QAFileName.Data());
      res = gGrid->Command(command);
      if (!res) {
	Error("PlotMuonQAPbPb",Form("no result for the command: %s",command.Data()));
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
        Int_t iEntry=0, iFile=0;

        while(kTRUE){
          iEntry++;
          const char* dirFilename = gSystem->GetDirEntry(dir);
          if(!dirFilename) break;
          sDirFilename = dirFilename;
          if(!sDirFilename.IsDigit()) continue;
          iFile++;
          objs = new TObjString(Form("%s/%s/AnalysisResults.root", alienBaseDir.Data(), sDirFilename.Data()));
          res->Add(objs);
        }
      }
    }

    // Loop over the 'find' results and get next LFN
    TIter nextmap(res);
    TMap *map = 0;
    TObjString *objs;

    //some checks
    Int_t iLoop=0, iLoopMax=1000;

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
	Error("PlotMuonQA","turl not found for the run %s... SKIPPING", run.Data());
	continue;
      }
      
      // open the outfile for this run
      TFile *runFile = TFile::Open(objs->GetString());
      if (!runFile || ! runFile->IsOpen()) {
	Error("PlotMuonQA", Form("failed to open file: %s", objs->GetName()));
	return;
      }

      // change here
      runFile->Cd("MUON_QA");
      
      // get interesting histos
      TObjArray* general1 = static_cast<TObjArray*>(runFile->FindObjectAny("general1"));
      TObjArray* general2 = static_cast<TObjArray*>(runFile->FindObjectAny("general2"));
      TObjArray* expert = static_cast<TObjArray*>(runFile->FindObjectAny("expert"));

      if (!general1 || !general2 || !expert){
        Error("PlotMUONQAPbPb", Form("All objects not here !!! ===> Skipping...for %s",objs->GetName()));
        continue;
      }

      TH1* hNClustersPerTrack = static_cast<TH1*>(general1->FindObject("hNClustersPerTrack"));
      TH1* hNChamberHitPerTrack = static_cast<TH1*>(general1->FindObject("hNChamberHitPerTrack"));
      TH1* hChi2 = static_cast<TH1*>(general1->FindObject("hChi2"));
      TH1* hNClustersPerCh = static_cast<TH1*>(general2->FindObject("hNClustersPerCh"));

      TH2* hClusterHitMapInCh[10];
      for(Int_t ich=0; ich<10; ich++) hClusterHitMapInCh[ich] = static_cast<TH2*>(expert->FindObject(Form("hClusterHitMapInCh%d",ich+1)));

      // do not skip empty runs butrather set values to 0
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
      run.ReplaceAll(QAFileName.Data(),"");
      run.Remove(TString::kTrailing, '/');
      hNClustersPerTrackVsRun_Mean->GetXaxis()->SetBinLabel(ibin, run.Data());
      hNClustersPerTrackVsRun_Sigma->GetXaxis()->SetBinLabel(ibin, run.Data());
      hNChamberHitPerTrack_Mean->GetXaxis()->SetBinLabel(ibin, run.Data());
      hNChamberHitPerTrack_Sigma->GetXaxis()->SetBinLabel(ibin, run.Data());
      hChi2_Mean->GetXaxis()->SetBinLabel(ibin, run.Data());
      hChi2_Sigma->GetXaxis()->SetBinLabel(ibin, run.Data());
      for (Int_t ich=0; ich<10; ich++) {
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
  cNClusters->SetLogy(0);
  cNClusters->Divide(1,2);
  cNClusters->cd(1);
  //hNClustersPerTrackVsRun_Mean->SetMaximum(11);
  hNClustersPerTrackVsRun_Mean->SetMinimum(7);
  hNClustersPerTrackVsRun_Mean->SetStats(kFALSE);
  hNClustersPerTrackVsRun_Mean->GetXaxis()->SetRange(1,ibin-1);
  hNClustersPerTrackVsRun_Mean->GetXaxis()->SetNdivisions(1,kFALSE);
  //hNClustersPerTrackVsRun_Mean->LabelsOption("u");
  hNClustersPerTrackVsRun_Mean->SetLabelSize(LabelSize);
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
  hNClustersPerTrackVsRun_Sigma->SetLabelSize(LabelSize);
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
  hNClustersInCh[0]->SetLabelSize(LabelSize);
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
  hClusterHitMapXInCh[0]->SetLabelSize(LabelSize);
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
  hClusterHitMapYInCh[0]->SetLabelSize(LabelSize);
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
  hChi2_Mean->SetLabelSize(LabelSize);
  hChi2_Mean->SetLineWidth(2);
  hChi2_Mean->Draw("e");

  cChi2->cd(2);
  hChi2_Sigma->SetStats(kFALSE);
  hChi2_Sigma->GetXaxis()->SetRange(1,ibin-1);
  hChi2_Sigma->GetXaxis()->SetNdivisions(1,kFALSE);
  //hChi2_Sigma->LabelsOption("u");
  hChi2_Sigma->SetLabelSize(LabelSize);
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
    THashList* labels = hMBB->GetXaxis()->GetLabels();
    TString format(Form("\n%%%ds %%9d",0));
  
    // print value for each label
    TObjString* label = 0x0;
    TIter nextLabel(labels);
    cout << "-------------------------------------------------" << endl;
    cout << "Run Number" << "\t Number of CMBB after Phys. Sel. " << endl ;  
    while ((label = static_cast<TObjString*>(nextLabel()))) {
      Int_t bin = (Int_t) label->GetUniqueID();
      printf(format.Data(), label->String().Data(), (Int_t) hMBB->GetBinContent(bin));
    }
    printf("\n\n");


    cout << "-------------------------------------------------" << endl;
    cout << "Total statistic" << endl; 
    cout << " " << endl ; 
    cout << "Number of  NumOfAllTriggers " << endl ;
    cout << "\t before selection " << NumOfAllTriggersNoPS << "\t  after selection " <<   NumOfAllTriggersWithPS << " --> rejection = " <<  (Double_t) (NumOfAllTriggersNoPS-NumOfAllTriggersWithPS)/(NumOfAllTriggersNoPS)*100. << "%" << endl ; 
    cout << "Number of  NumOfCMBB " << endl ;
    cout << "\t before selection " << NumOfMBBNoPS << "\t  after selection " <<   NumOfMBBWithPS << " --> rejection = " <<  (Double_t) (NumOfMBBNoPS-NumOfMBBWithPS)/(NumOfMBBNoPS)*100. << "%" << endl ; 
    cout << " " << endl ; 
    
  }
}
