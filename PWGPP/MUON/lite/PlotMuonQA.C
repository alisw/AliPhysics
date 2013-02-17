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
//  - libPWGmuon.so
//
// TString includePath = "-I${ALICE_ROOT}/PWGmuon ";  gSystem->SetIncludePath(includePath.Data());
//
// The macro reads results of the QA task and produce monitoring plots.
//
// Authors: Philippe Pillot - SUBATECH Nantes, Christophe Suire, Cynthia Hadjidakis - IPN Orsay
// To be done:
// - reorganize last part (reading and extracting info from run per run histos)

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
#include "TF1.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TObjArray.h"
#include "TObjString.h"

// ALIROOT includes
#include "AliCounterCollection.h"

#endif

TString GetRunList(const char *runList, TObjArray *runs, TObjArray *runs2);
Bool_t GetTriggerLists(const char *triggerList, TString listFromContainer, TObjArray *triggersB=0, TObjArray *triggersAC=0, TObjArray *triggersE=0, TObjArray *triggersShortName=0);
void SetCanvas(TCanvas *canvas, Int_t logy=1);

TH1* ProcessHisto( AliCounterCollection* counter, TString variable, TString selection, TString hName="", TString xName="", TString yName="", Int_t color=1);
TH2* ProcessHisto2D( AliCounterCollection* counter, TString hVariable, TString hVariable2, TString hSelection, TString hName);

TCanvas *ProcessCanvasTriggerContent(TObjArray *array, TH1 **histo, TH1 **histo2, TString canvasName);
TCanvas *ProcessCanvasRelativeTriggerContent(TObjArray *array, TH1 **histo, TString canvasName, TArrayI *colorInd);
TCanvas *ProcessCanvasPhysSelCut(TObjArray *triggersB, TObjArray *triggersAC, TObjArray *triggersE, TH1 **hBNoPS, TH1 **hACNoPS,TH1 **hENoPS, TH1 **hBWithPS, TString canvasName, TArrayI *colorInd, Bool_t isHeavyIon = kFALSE);
TCanvas *ProcessCanvasPhysSelCutCentrality(TObjArray *triggersB, TObjArray *triggersAC, TObjArray *triggersE, TH1 **hBNoPS, TH1 **hACNoPS, TH1 **hENoPS, TH1 **hBWithPS, Int_t k, TString canvasName, TArrayI *colorInd, TString *legendHeader, Bool_t isHeavyIon = kFALSE);
TCanvas *ProcessCanvasCentralityPercentile(TObjArray *triggersB, TH1 **hBNoPSCent, TH1 **hBWithPSCent, Int_t k, TString canvasName, TArrayI *colorInd, TString *legendHeader);
TCanvas *ProcessCanvasTracksoverTrigger(TObjArray *triggersB, TH1 **hB, TH1 **hTrackerB, TH1 **hTriggerB, TH1 **hMatchedB, TH1 **hAllTracksB, Int_t indTrigger, TString canvasName,TString legendHeader="");
TCanvas *ProcessCanvasTrackMultB(TObjArray *triggersB, TH1 **hB, TH1 **hTrackerB, TH1 **hTriggerB, TH1 **hMatchedB, TH1 **hAllTracksB, Int_t indTrigger, TString canvasName,TString legendHeader="");
TCanvas *ProcessCanvasRatioTrackB(TObjArray *triggersB, TH1 **hB, TH1 **hTrackerB, TH1 **hTriggerB, TH1 **hMatchedB, TH1 **hAllTracksB, Int_t indTrigger, TString canvasName,TString legendHeader="");
TCanvas *ProcessCanvasAsymMatched(TObjArray *triggersB, TH1 **hPosMatchedB, TH1 **hNegMatchedB, TH1 **hAllMatchedB, Int_t indTrigger, TString canvasName,TString legendHeader="");
TCanvas *ProcessCanvasHighPtMuons(TObjArray *triggersB, TH1 **hB, TH1 **hMatchedLowPtB, TH1 **hAllMatchedHightPtB, Int_t indTrigger, TString canvasName,TString legendHeader="");
TCanvas *ProcessCanvasBeamGasMatched(TObjArray *triggersB, TH1 **hBeamGasMatchedB, TH1** hBeamGasMatchedHighPtB, TH1 **hAllMatchedB, TH1** hMatchedHighPtB, Int_t indTrigger, TString canvasName,TString legendHeader="");

Bool_t IsTrigger(TObjArray *array, Int_t index, TString name);
Bool_t IsTriggerSelectedForMuonPhysics(TObjArray *array, Int_t index);
Bool_t IsHeavyIonCollision(AliCounterCollection *eventCounters);

const Int_t kNMaxTriggers = 25;

//--------------------------------------------------------------------------
void PlotMuonQA(const char* baseDir, const char* runList = 0x0, const char * triggerList = 0x0, Bool_t selectPhysics = kFALSE, const char *LHCPeriod = "LHC11c", const char *QAFileName = "QAresults.root") {
	
  /// Macro for QA monitoring.
  /// Example: baseDir = "alien:///alice/cern.ch/user/p/ppillot/pp7TeV/LHC10d/MuonQA/pass1/results/".
  /// If runList != 0x0: only the given runs will be used. Otherwise use all runs found in baseDir.
  /// If triggerList !=0x0: only the given triggers are displayed. Otherwise use the default list of triggers (see GetTriggerLists)
	
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
  gSystem->Load("libPWGHFbase");
  gSystem->Load("libPWGmuon");
#endif
  
  // Cosmetics and configuration
  gStyle->SetFillColor(10);
  gStyle->SetPadGridX(kTRUE);
  gStyle->SetPadGridY(kTRUE);
  gStyle->SetPadRightMargin(0.01);
  
  TString OutFileName = "QA_";  OutFileName += LHCPeriod;
  TString OutFileNamePDF=  OutFileName.Data();  OutFileNamePDF+= ".pdf";
  TString OutFileNamePDF_open = OutFileNamePDF.Data(); OutFileNamePDF_open += "[";  
  TString OutFileNamePDF_close= OutFileNamePDF.Data(); OutFileNamePDF_close += "]";  
  TString OutFileNameROOT=  OutFileName.Data();  OutFileNameROOT+= ".root";
  
  Int_t PRINTSTAT = 1;
  Int_t kCMUS = 2;
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
  
  cout<<"//---------------------------------- //"<<endl;
  cout<<"//          Run selection            //"<<endl;
  cout<<"//---------------------------------- //"<<endl;
  
  TObjArray *runs = new TObjArray();
  runs->SetOwner(kTRUE);
  TObjArray *runs2 = new TObjArray();
  runs2->SetOwner(kTRUE);
  TString selectRuns = GetRunList(runList,runs,runs2);
	
  // physics selection
  TString select = selectPhysics ? "selected:yes" : "";
  	
  cout<<"//---------------------------------- //"<<endl;
  cout<<"//        Get global counter        //"<<endl;
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
	
  cout<<"//---------------------------------- //"<<endl;
  cout<<"//        Trigger selection          //"<<endl;
  cout<<"//---------------------------------- //"<<endl;

  TObjArray *triggersB, *triggersAC, *triggersE, *triggersShortName;
  triggersB = new TObjArray();
  triggersB->SetOwner();
  triggersAC = new TObjArray();
  triggersAC->SetOwner();
  triggersE = new TObjArray();
  triggersE->SetOwner();
  triggersShortName = new TObjArray();
  triggersShortName->SetOwner();

  TString listFromContainer = eventCounters->GetKeyWords("trigger");
  Bool_t success = GetTriggerLists(triggerList, listFromContainer, triggersB, triggersAC, triggersE, triggersShortName);
  if(!success) return;


	
  cout<<"//---------------------------------- //"<<endl;
  cout<<"//      Set collision type ?          //"<<endl;
  cout<<"//---------------------------------- //"<<endl;

  Bool_t isHeavyIon = kTRUE;
  isHeavyIon = IsHeavyIonCollision(eventCounters);
	
  cout<<"//---------------------------------- //"<<endl;
  cout<<"//        Trigger plots              //"<<endl;
  cout<<"//---------------------------------- //"<<endl;
	
  //plot all trigger from event counters without any selection
  TString CanvasName = "cAll";
  TCanvas *cAll = new TCanvas(CanvasName.Data(),CanvasName.Data());
  cAll->SetLeftMargin(0.18);
  cAll->SetRightMargin(0.18);
  cAll->SetLogz(1);
  cAll->cd();
  //TH2* hAll = (TH2*) ProcessHisto2D(eventCounters, "trigger", "run", Form("run:any/%s",select.Data()) , "");
  TH2* hAll = (TH2*) ProcessHisto2D(eventCounters, "trigger", "run", "run:any" , "");
  for ( Int_t ibin=1; ibin<=hAll->GetYaxis()->GetNbins(); ++ibin ) {
    TString currLabel = hAll->GetYaxis()->GetBinLabel(ibin);
    TObjArray* labelArray = currLabel.Tokenize("-");
    labelArray->SetOwner();
    //cout<<currLabel<<endl;
    //labelArray->Print();
    TString newLabel = labelArray->At(0)->GetName();
    if ( labelArray->GetEntries() >= 2 ) newLabel = Form("%s-%s", newLabel.Data(), labelArray->At(1)->GetName());
    hAll->GetYaxis()->SetBinLabel(ibin, newLabel.Data());
    delete labelArray;
  }
  hAll->Draw("COLZ");

  //declare a default canvas c1 
  CanvasName = "c1";
  TCanvas *c1 = new TCanvas(CanvasName.Data(),CanvasName.Data());
  c1->cd();
	
  //loop on centrality
  Int_t centBin = 3;
  Int_t const centBinMax = 4;
  Int_t centBinMaxLoop = centBinMax;
  if(!isHeavyIon) centBinMaxLoop = 1;
  TString centBinName[centBinMax] = {"v0mult:any/","v0mult:low,int,high/","v0mult:low/","v0mult:high/"};
  TString centLegendName[centBinMax] ={"All collisions","[0-80%] from V0 amplitude","low mult. [60-80%] from V0 amplitude","high mult. [0-10%] from V0 amplitude"};
  TString centLegendNameShort[centBinMax] ={"All","[0-80%]","[60-80%]","[0-10%]"};
  TString selectionCent;
	       
  Int_t *colorTab = new Int_t[triggersB->GetEntriesFast()];
  Int_t const colorNrFirst = 8;
  Int_t colorIndex = 0;
  Int_t colorTabFirst[colorNrFirst] = {kGray,kRed,kBlue,kGreen,kOrange,kCyan,kMagenta,kYellow};
  for (Int_t i = 0; i< triggersB->GetEntriesFast(); i++ ) {
    colorTab[i] = colorTabFirst[i] + colorIndex;
    if ( i%colorNrFirst == 0 ) colorIndex++;
  } 
  TArrayI *colorInd = new TArrayI( triggersB->GetEntriesFast() );
  for(Int_t i=0; i< triggersB->GetEntriesFast(); i++) colorInd->AddAt(colorTab[i],i); 

  TH1* hBNoPS[centBinMax][kNMaxTriggers];
  TH1* hBWithPS[centBinMax][kNMaxTriggers];
  TH1* hBWithPST0Flag[centBinMax][kNMaxTriggers]; 
  TH1* hBWithPST0SPDFlag[centBinMax][kNMaxTriggers]; 
  TH1* hB[centBinMax][kNMaxTriggers]={};
  TH1* hTriggerB[centBinMax][kNMaxTriggers], *hTrackerB[centBinMax][kNMaxTriggers], *hMatchedB[centBinMax][kNMaxTriggers], *hAllTracksB[centBinMax][kNMaxTriggers], *hMatchedLowPtB[centBinMax][kNMaxTriggers], *hMatchedHighPtB[centBinMax][kNMaxTriggers];
  TH1* hMatchedLowPtBNoPS[centBinMax][kNMaxTriggers], *hMatchedHighPtBNoPS[centBinMax][kNMaxTriggers];
  TH1* hPosMatchedB[centBinMax][kNMaxTriggers], *hNegMatchedB[centBinMax][kNMaxTriggers], *hAllMatchedB[centBinMax][kNMaxTriggers];
  TH1* hBeamGasMatchedB[centBinMax][kNMaxTriggers], *hBeamGasMatchedHighPtB[centBinMax][kNMaxTriggers];
  TH1 *hACWithPS[centBinMax][kNMaxTriggers]={}; 
  TH1 *hACNoPS[centBinMax][kNMaxTriggers]={};
  TH1 *hEWithPS[centBinMax][kNMaxTriggers]={};
  TH1 *hENoPS[centBinMax][kNMaxTriggers]={};
  
  if(triggersB->GetEntriesFast()>=kNMaxTriggers){
    cout<<"Too many triggers = "<<triggersB->GetEntriesFast()<<endl;
    return;
  }
  
  //loop on centrality
  for(centBin = 0; centBin < centBinMaxLoop; centBin++){
    selectionCent = centBinName[centBin];
		
    //Loop on trigger (last is all triggers, then for each defined trigger)
    for(Int_t i = 0; i < triggersB->GetEntriesFast(); i++){
    
      TString histoNameBase = "h_trig", histoName;
      histoNameBase+= i+1;
 		
      TString triggerName = ( (TObjString*) triggersB->At(i) )->GetString();
      if(triggerName.EqualTo(" ")) continue;
      // Histo trigger without Phys. Sel. 
      selection = selectionCent; selection += Form("trigger:%s/%s", triggerName.Data(), selectRuns.Data());		
      // cout<<selection<<endl;
      histoName = histoNameBase;
      histoName += "BNoPS";
      hBNoPS[centBin][i] = (TH1*) ProcessHisto(eventCounters, "run", selection, histoName, "", "Trigger content w/o Phys. Sel.", colorInd->At(i));
      // Histo trigger with Phys. Sel. 
      selection = selectionCent; selection += Form("trigger:%s/%s/selected:yes", triggerName.Data(), selectRuns.Data());
      histoName = histoNameBase;
      histoName += "BWithPS";
      hBWithPS[centBin][i] = (TH1*) ProcessHisto(eventCounters, "run", selection, histoName, "", "Trigger content w/ Phys. Sel.", colorInd->At(i));
    // Histo trigger with Phys. Sel. and T0 pile up not flagged
      selection = selectionCent; selection += Form("trigger:%s/%s/selected:yes/t0pileup:no", triggerName.Data(), selectRuns.Data());
      histoName = histoNameBase;
      histoName += "BWithPST0Flag";
      hBWithPST0Flag[centBin][i] = (TH1*) ProcessHisto(eventCounters, "run", selection, histoName, "", "Trigger content w/ Phys. Sel. and no pile up from T0 flag", colorInd->At(i));
    // Histo trigger with Phys. Sel. and T0 + SPD pile up not flagged
      selection = selectionCent; selection += Form("trigger:%s/%s/selected:yes/t0pileup:no/spdpileup:no", triggerName.Data(), selectRuns.Data());
      histoName = histoNameBase;
      histoName += "BWithPST0SPDFlag";
      hBWithPST0SPDFlag[centBin][i] = (TH1*) ProcessHisto(eventCounters, "run", selection, histoName, "", "Trigger content w/ Phys. Sel. and no pile up from T0 and SPD flag", colorInd->At(i));
      // Histo trigger : Phys. Sel.  is selected or not depending on the macro arguments
      selection = selectionCent; selection += Form("trigger:%s/%s/%s", triggerName.Data(), selectRuns.Data(), select.Data());
      histoName = histoNameBase;
      histoName += "B";
      hB[centBin][i] = (TH1*) ProcessHisto(eventCounters, "run", selection, histoName);
		
      TString triggerNameAC = ( (TObjString*) triggersAC->At(i) )->GetString();
      // Histo trigger without Phys. Sel. AC
      histoName = histoNameBase;
      histoName += "ACNoPS";
      selection = selectionCent; selection += Form("trigger:%s/%s", triggerNameAC.Data(), selectRuns.Data());
      hACNoPS[centBin][i] =  (TH1*) ProcessHisto(eventCounters, "run", selection, histoName);
      // Histo trigger with Phys. Sel. AC
      selection = selectionCent; selection += Form("trigger:%s/%s/selected:yes", triggerNameAC.Data(), selectRuns.Data());
      histoName = histoNameBase;
      histoName += "ACWithPS";
      hACWithPS[centBin][i] =  (TH1*) ProcessHisto(eventCounters, "run", selection, histoName);
    
      TString triggerNameE = ( (TObjString*) triggersE->At(i) )->GetString();
      // Histo trigger without Phys. Sel. E
      selection = selectionCent; selection += Form("trigger:%s/%s", triggerNameE.Data(), selectRuns.Data());
      histoName = histoNameBase;
      histoName += "ENoPS";
      hENoPS[centBin][i] =  (TH1*) ProcessHisto(eventCounters, "run", selection, histoName);
      // Histo trigger with Phys. Sel. E
      selection = selectionCent; selection += Form("trigger:%s/%s/selected:yes", triggerNameE.Data(), selectRuns.Data());
      histoName = histoNameBase;
      histoName += "EWithPS";
      hEWithPS[centBin][i] =  (TH1*) ProcessHisto(eventCounters, "run", selection, histoName);

      // Histo tracking : Phys. Sel.  is selected or not depending on the macro arguments
      selection = selectionCent; selection += Form("track:triggeronly/trigger:%s/%s/%s", triggerName.Data(), selectRuns.Data(), select.Data());
      hTriggerB[centBin][i] = (TH1*) ProcessHisto(trackCounters, "run", selection, "");
		
      selection = selectionCent; selection += Form("track:trackeronly/trigger:%s/%s/%s", triggerName.Data(), selectRuns.Data(), select.Data());
      hTrackerB[centBin][i] = (TH1*) ProcessHisto(trackCounters, "run", selection, "");
		
      selection = selectionCent; selection += Form("track:matched/trigger:%s/%s/%s", triggerName.Data(), selectRuns.Data(), select.Data());
      hMatchedB[centBin][i] = (TH1*) ProcessHisto(trackCounters, "run", selection, "");
		
      selection = selectionCent; selection += Form("trigger:%s/%s/%s", triggerName.Data(), selectRuns.Data(), select.Data());
      hAllTracksB[centBin][i] = (TH1*) ProcessHisto(trackCounters, "run", selection, "");
		
      selection = selectionCent; selection += Form("track:matched/trigger:%s/%s/%s/pt:low/acc:in", triggerName.Data() ,selectRuns.Data(), select.Data());
      hMatchedLowPtB[centBin][i] = (TH1*) ProcessHisto(trackCounters, "run", selection, "");
		
      selection = selectionCent; selection += Form("track:matched/trigger:%s/%s/%s/pt:high/acc:in", triggerName.Data() ,selectRuns.Data(), select.Data());
      hMatchedHighPtB[centBin][i] = (TH1*) ProcessHisto(trackCounters, "run", selection, "");
		
      selection = selectionCent; selection += Form("track:matched/trigger:%s/%s/pt:low/acc:in", triggerName.Data() ,selectRuns.Data());
      hMatchedLowPtBNoPS[centBin][i] = (TH1*) ProcessHisto(trackCounters, "run", selection, "");
		
      selection = selectionCent; selection += Form("track:matched/trigger:%s/%s/pt:high/acc:in", triggerName.Data() ,selectRuns.Data());
      hMatchedHighPtBNoPS[centBin][i] = (TH1*) ProcessHisto(trackCounters, "run", selection, "");
		
      selection = selectionCent; selection += Form("track:matched/trigger:%s/%s/charge:pos/%s/acc:in",triggerName.Data(), select.Data(),selectRuns.Data());
      hPosMatchedB[centBin][i] = (TH1*) ProcessHisto(trackCounters, "run", selection, "");
    
      selection = selectionCent; selection += Form("track:matched/trigger:%s/%s/charge:neg/%s/acc:in",triggerName.Data(), select.Data(),selectRuns.Data());
      hNegMatchedB[centBin][i] =  (TH1*) ProcessHisto(trackCounters, "run", selection, "");
		
      selection = selectionCent; selection += Form("track:matched/trigger:%s/%s/%s/acc:in",triggerName.Data(), select.Data(),selectRuns.Data());
      hAllMatchedB[centBin][i] =  (TH1*) ProcessHisto(trackCounters, "run", selection, "");

      selection = selectionCent; selection += Form("track:matched/trigger:%s/%s/%s/acc:in/tagTrack:beamGas",triggerName.Data(), select.Data(),selectRuns.Data());
      hBeamGasMatchedB[centBin][i] =  (TH1*) ProcessHisto(trackCounters, "run", selection, "");

      selection = selectionCent; selection += Form("track:matched/trigger:%s/%s/%s/acc:in/tagTrack:beamGas/pt:high",triggerName.Data(), select.Data(),selectRuns.Data());
      hBeamGasMatchedHighPtB[centBin][i] =  (TH1*) ProcessHisto(trackCounters, "run", selection, "");

    }
  }
	
  //if there is not B triggers just stop now
  Int_t count_trigger=0;
  centBin = 0;
  for(Int_t i = 0; i < triggersB->GetEntriesFast(); i++){
    count_trigger += hBNoPS[centBin][i]->GetEntries();
  }
  if(count_trigger<=0) return;
	
	
  Int_t NumOfBNoPS[centBinMax][kNMaxTriggers];
  Int_t NumOfBWithPS[centBinMax][kNMaxTriggers];
  Int_t NumOfBWithPST0Flag[centBinMax][kNMaxTriggers];
  Int_t NumOfBWithPST0SPDFlag[centBinMax][kNMaxTriggers];
  Int_t NumOfACNoPS[centBinMax][kNMaxTriggers];
  Int_t NumOfENoPS[centBinMax][kNMaxTriggers];
  Int_t NumOfACWithPS[centBinMax][kNMaxTriggers];
  Int_t NumOfEWithPS[centBinMax][kNMaxTriggers];
	
  for(centBin = 0; centBin < centBinMaxLoop; centBin++){
    for(Int_t i = 0; i < triggersB->GetEntriesFast(); i++){
      NumOfBNoPS[centBin][i] = hBNoPS[centBin][i]->Integral();
      NumOfBWithPS[centBin][i] = hBWithPS[centBin][i]->Integral();
      NumOfBWithPST0Flag[centBin][i] = hBWithPST0Flag[centBin][i]->Integral();
      NumOfBWithPST0SPDFlag[centBin][i] = hBWithPST0Flag[centBin][i]->Integral();
      NumOfACNoPS[centBin][i] = hACNoPS[centBin][i]->Integral();
      NumOfENoPS[centBin][i] = hENoPS[centBin][i]->Integral();
      NumOfACWithPS[centBin][i] = hACWithPS[centBin][i]->Integral();
      NumOfEWithPS[centBin][i] = hEWithPS[centBin][i]->Integral();
    }
  }
  centBin = 0;
	
  cout<<"//==================================================================================="<<endl;
  cout<<"// Put all plots in a ps file, easier to publish (Twiki)"<<endl;
  cout<<"//==================================================================================="<<endl;
  
  c1->Print(OutFileNamePDF_open.Data());
  TFile *rootFileOut = TFile::Open(OutFileNameROOT.Data(),"RECREATE");
  rootFileOut->cd();
  TDirectoryFile *dirGlobal = new TDirectoryFile("Global","Global","",(TDirectory*)rootFileOut->GetMotherDir());
  cout<<"dirGlobal mother "<<(dirGlobal->GetMotherDir())->GetName()<<endl;
  //An array of TDirectoryFile
  TObjArray *dirTrigger = new TObjArray;
  dirTrigger->SetOwner();
  TObjArray *dirCent = new TObjArray;
  dirCent->SetOwner();
  for ( Int_t i = 0; i < triggersB->GetEntriesFast()-1 ; i++) {
    TString currTrigName = ( (TObjString*) triggersShortName->At(i) )->GetString();
    TDirectoryFile *dirFile = new TDirectoryFile( currTrigName.Data(),currTrigName.Data(),"",(TDirectory*)rootFileOut->GetMotherDir() );
    dirTrigger->AddLast( dirFile );
    for( Int_t j = 0; j < centBinMaxLoop; j++) {
      TString centName = centLegendNameShort[j];
      TDirectoryFile *dirFileCent = new TDirectoryFile( centName.Data(),centName.Data(),"",dirFile );
      dirCent->AddLast( dirFileCent );
    }
  }

  cAll->Print(OutFileNamePDF.Data());
  dirGlobal->cd();
  cAll->Write();
  	
  cout<<"//==================================================================================="<<endl;
  cout<<"// new canvas with the total number of trigger with and without Phys. Sel."<<endl;
  cout<<"//==================================================================================="<<endl;
	
  TCanvas *cTriggerContent = ProcessCanvasTriggerContent(triggersB, hBNoPS[centBin], hBWithPS[centBin], "TriggerContent");
  cTriggerContent->Draw(); 
  cTriggerContent->Print(OutFileNamePDF.Data());
  dirGlobal->cd();
  cTriggerContent->Write();
  
  cout<<"//==================================================================================="<<endl;
  cout<<"// new canvas with the relative content of each trigger w/o physics selection"<<endl;
  cout<<"//==================================================================================="<<endl;
	
  TCanvas *cRelativeTriggerContent = ProcessCanvasRelativeTriggerContent(triggersB, hBNoPS[centBin], "RelativeTriggerContent", colorInd);
  cRelativeTriggerContent->Draw();
  cRelativeTriggerContent->Print(OutFileNamePDF.Data());
  dirGlobal->cd();
  cRelativeTriggerContent->Write();
	
  cout<<"//==================================================================================="<<endl;
  cout<<"// new canvas with effect from physics selection for each trigger and background trigger "<<endl;
  cout<<"//==================================================================================="<<endl;
	
  TCanvas *cPhysSelCut = 0x0;
  cPhysSelCut = ProcessCanvasPhysSelCut(triggersB, triggersAC, triggersE, hBNoPS[centBin], hACNoPS[centBin], hENoPS[centBin], hBWithPS[centBin], "PhysSelCutOnCollTrigger", colorInd,isHeavyIon);
  cPhysSelCut->Draw();
  cPhysSelCut->Print(OutFileNamePDF.Data());
  dirGlobal->cd();
  cPhysSelCut->Write();
  rootFileOut->cd();
	
  cout<<"//==================================================================================="<<endl;
  cout<<"// new canvas with effect from T0 pile-up flag for each trigger and background trigger "<<endl;
  cout<<"//==================================================================================="<<endl;
	
  TCanvas *cPhysSelCutT0Flag = 0x0;
  cPhysSelCutT0Flag = ProcessCanvasPhysSelCut(triggersB, triggersAC, triggersE, hBWithPS[centBin], hACNoPS[centBin], hENoPS[centBin], hBWithPST0Flag[centBin], "PhysSelCutOnCollTriggerT0Flag", colorInd,isHeavyIon);
  cPhysSelCutT0Flag->Draw();
  cPhysSelCutT0Flag->Print(OutFileNamePDF.Data());
  dirGlobal->cd();
  cPhysSelCutT0Flag->Write();
  rootFileOut->cd();
	

 cout<<"//==================================================================================="<<endl;
  cout<<"// new canvas with effect from T0 + SPD pile-up flag for each trigger and background trigger "<<endl;
  cout<<"//==================================================================================="<<endl;
	
  TCanvas *cPhysSelCutT0SPDFlag = 0x0;
  cPhysSelCutT0SPDFlag = ProcessCanvasPhysSelCut(triggersB, triggersAC, triggersE, hBWithPS[centBin], hACNoPS[centBin], hENoPS[centBin], hBWithPST0SPDFlag[centBin], "PhysSelCutOnCollTriggerT0Flag", colorInd,isHeavyIon);
  cPhysSelCutT0SPDFlag->Draw();
  cPhysSelCutT0SPDFlag->Print(OutFileNamePDF.Data());
  dirGlobal->cd();
  cPhysSelCutT0SPDFlag->Write();
  rootFileOut->cd();



  cout<<"//==================================================================================="<<endl;
  cout<<"// new canvas with effect from physics selection for each trigger and centrality bin (only in PbPb) "<<endl;
  cout<<"//==================================================================================="<<endl;
		
  Int_t k=0;
  TString canvasName;
  TString triggerName;	
  TString legendHeader="";
	
	
  if ( isHeavyIon ){
    TCanvas *cPhysSelCutCentrality;
	
    //loop over trigger
    for(k = 0; k < triggersB->GetEntriesFast(); k++){
      //skip sum of all triggers
      if(k == (triggersB->GetEntriesFast()-1)) continue;

      canvasName = "PhysSel_trigger";
      canvasName += ( (TObjString*) triggersShortName->At(k) )->GetString();
	
      TH1* hBNoPSCent[centBinMax-1]={hBNoPS[1][k],hBNoPS[2][k],hBNoPS[3][k]}; 
      TH1* hACNoPSCent[centBinMax-1]={hACNoPS[1][k],hACNoPS[2][k],hACNoPS[3][k]}; 
      TH1* hENoPSCent[centBinMax-1]={hENoPS[1][k],hENoPS[2][k],hENoPS[3][k]}; 
      TH1* hBWithPSCent[centBinMax-1]={hBWithPS[1][k],hBWithPS[2][k],hBWithPS[3][k]}; 
		
      cPhysSelCutCentrality = ProcessCanvasPhysSelCutCentrality(triggersB, triggersAC, triggersE, hBNoPSCent, hACNoPSCent, hENoPSCent, hBWithPSCent, k, canvasName, colorInd, centLegendNameShort+1, isHeavyIon);
      cPhysSelCutCentrality->Draw();
      cPhysSelCutCentrality->Print(OutFileNamePDF.Data());
      ( (TDirectoryFile*) dirTrigger->At(k) )->cd();
      cPhysSelCutCentrality->Write();		
    }
  }
  rootFileOut->cd();
	
  cout<<"//==================================================================================="<<endl;
  cout<<"// new canvas for centrality percentile check (only in PbPb) "<<endl;
  cout<<"//==================================================================================="<<endl;
	
  if ( isHeavyIon ){
    TCanvas *cCentralityCheck;
		
    //loop over trigger
    for(k = 0; k < triggersB->GetEntriesFast(); k++){
      //skip sum of all triggers
      if(k == (triggersB->GetEntriesFast()-1)) continue;
 			
      canvasName = "CentralityCheck_trigger";
      canvasName +=( (TObjString*) triggersShortName->At(k) )->GetString();
			
      TH1* hBNoPSCent[centBinMax-1]={hBNoPS[1][k],hBNoPS[2][k],hBNoPS[3][k]}; 
      TH1* hBWithPSCent[centBinMax-1]={hBWithPS[1][k],hBWithPS[2][k],hBWithPS[3][k]}; 
			
      cCentralityCheck = ProcessCanvasCentralityPercentile(triggersB,hBNoPSCent,hBWithPSCent,k,canvasName,colorInd,centLegendNameShort+1 );

      cCentralityCheck->Draw();
      cCentralityCheck->Print(OutFileNamePDF.Data());
      ( (TDirectoryFile*) dirTrigger->At(k) )->cd();
      cCentralityCheck->Write();
			
    }
  }
  rootFileOut->cd();
	

  cout<<"//==================================================================================="<<endl;
  cout<<"// Ratio of tracks over trigger type (3 canvases) "<<endl;
  cout<<"//==================================================================================="<<endl;

  //Print a canvas per trigger type
  TCanvas *cTracksoverTrigger;
  TCanvas* cTrackMultB;
  TCanvas* cRatioTrackB;
 	
  //loop on centrality bin
  for ( centBin = 0; centBin < centBinMaxLoop; centBin++){
    if ( isHeavyIon ){
      legendHeader = "for ";
      legendHeader += centLegendName[centBin];
    }
    else legendHeader ="";
    //loop over trigger
    for(k = 0; k < triggersB->GetEntriesFast(); k++){
      //skip sum of all triggers
      if(k == (triggersB->GetEntriesFast()-1)) continue;

      ( (TDirectoryFile*) dirCent->At( k*centBinMaxLoop+centBin ) )->cd();

      canvasName = "RatioTrackTypes_cent";
      canvasName += centBin;
      canvasName +="trigger";
      canvasName += ( (TObjString*) triggersShortName->At(k) )->GetString();
      cTracksoverTrigger = ProcessCanvasTracksoverTrigger(triggersB, hB[centBin], hTrackerB[centBin], hTriggerB[centBin], hMatchedB[centBin], hAllTracksB[centBin], k, canvasName,legendHeader);
      cTracksoverTrigger->Draw();
      cTracksoverTrigger->Print(OutFileNamePDF.Data());
      cTracksoverTrigger->Write();

      canvasName = "TrackMult_cent";
      canvasName += centBin;
      canvasName +="trigger";
      canvasName +=( (TObjString*) triggersShortName->At(k) )->GetString();		
      cTrackMultB= ProcessCanvasTrackMultB(triggersB, hB[centBin], hTrackerB[centBin], hTriggerB[centBin], hMatchedB[centBin], hAllTracksB[centBin], k, canvasName, legendHeader);
      cTrackMultB->Draw();
      cTrackMultB->Print(OutFileNamePDF.Data());
      cTrackMultB->Write();
	
      canvasName = "RatioTrackB_cent";
      canvasName += centBin;
      canvasName +="trigger";
      canvasName +=( (TObjString*) triggersShortName->At(k) )->GetString();		
      cRatioTrackB = ProcessCanvasRatioTrackB(triggersB, hB[centBin], hTrackerB[centBin], hTriggerB[centBin], hMatchedB[centBin], hAllTracksB[centBin], k, canvasName, legendHeader);
      cRatioTrackB->Draw();
      cRatioTrackB->Print(OutFileNamePDF.Data());
      cRatioTrackB->Write();
    }
  }
  rootFileOut->cd();

  cout<<"//===================================================="<<endl;
  cout<<"// Draw matched tracks asymmetry for mus type trigger "<<endl;
  cout<<"//===================================================="<<endl;
	
  //Print a canvas per trigger type
  TCanvas *cAsymMatched;

  //Loop on centrality
  for ( centBin = 0; centBin < centBinMaxLoop; centBin++){
    if ( isHeavyIon ){
      legendHeader = "for ";
      legendHeader += centLegendName[centBin];
    }
    else legendHeader ="";
    //loop over trigger
    for(k = 0; k < triggersB->GetEntriesFast(); k++){
      //skip sum of all triggers
      if(k == (triggersB->GetEntriesFast()-1)) continue;

      ( (TDirectoryFile*) dirCent->At( k*centBinMaxLoop+centBin ) )->cd();
 
      canvasName = "AsymMatched";
      canvasName += centBin;
      canvasName +="trigger";
      canvasName +=( (TObjString*) triggersShortName->At(k) )->GetString();
      cAsymMatched = ProcessCanvasAsymMatched(triggersB, hPosMatchedB[centBin], hNegMatchedB[centBin], hAllMatchedB[centBin], k, canvasName,legendHeader);
      cAsymMatched->Draw();
      cAsymMatched->Print(OutFileNamePDF.Data());
      cAsymMatched->Write();
    }
  }
  rootFileOut->cd();
  legendHeader = ""; 

  cout<<"//===================================================================="<<endl;
  cout<<"// Draw beam gas contribution to matched tracks  for mus type trigger "<<endl;
  cout<<"//===================================================================="<<endl;
	
  //Print a canvas per trigger type
  TCanvas *cBeamGasMatched;

  //loop over trigger
  for(k = 0; k < triggersB->GetEntriesFast(); k++){
    //skip sum of all triggers
    if(k == (triggersB->GetEntriesFast()-1)) continue;

    ( (TDirectoryFile*) dirTrigger->At(k) )->cd();
 
    canvasName = "BeamGasMatched";
    canvasName +="trigger";
    canvasName +=( (TObjString*) triggersShortName->At(k) )->GetString();
    centBin = 0;
    cBeamGasMatched= ProcessCanvasBeamGasMatched(triggersB, hBeamGasMatchedB[centBin], hBeamGasMatchedHighPtB[centBin], hAllMatchedB[centBin], hMatchedHighPtB[centBin], k, canvasName,legendHeader);
    cBeamGasMatched->Draw();
    cBeamGasMatched->Print(OutFileNamePDF.Data());
    cBeamGasMatched->Write();
  }
  rootFileOut->cd();

  cout<<"//=================================================="<<endl;
  cout<<"// Draw high pt tracks per trigger"<<endl;
  cout<<"//=================================================="<<endl;

  //Print a canvas per trigger type
  TCanvas *cHighPtMuons;
	
  //Loop on centrality
  for ( centBin = 0; centBin < centBinMaxLoop; centBin++){
    if ( isHeavyIon ){
      legendHeader = "for ";
      legendHeader += centLegendName[centBin];
    }
    else legendHeader ="";
    //loop over trigger
    for(k = 0; k < triggersB->GetEntriesFast(); k++){
      //skip sum of all triggers
      if(k == (triggersB->GetEntriesFast()-1)) continue;
	
      ( (TDirectoryFile*) dirCent->At( k*centBinMaxLoop+centBin ) )->cd();

      canvasName = "HighPtMuons";
      canvasName += centBin;
      canvasName +="trigger";
      canvasName +=( (TObjString*) triggersShortName->At(k) )->GetString();
			
      cHighPtMuons = ProcessCanvasHighPtMuons(triggersB, hB[centBin], hMatchedLowPtB[centBin], hMatchedHighPtB[centBin], k, canvasName,legendHeader);
      cHighPtMuons->Draw();
      cHighPtMuons->Print(OutFileNamePDF.Data());
      cHighPtMuons->Write();
    }
  }
  rootFileOut->cd();
	
  // close merged file	
  globalFile->Close();
  
  //===================================================================================
  //Print out the number of trigger without and with Phys. Sel.
  //===================================================================================
  
  centBin = 0;
  cout << endl << endl;
  //====================================================
  if (PRINTSTAT){
    if ( triggersB->At(kCMUS) ) { 
    
      // set the format to print labels
      THashList* labels = hBWithPS[centBin][kCMUS]->GetXaxis()->GetLabels();
      TString format(Form("\n%%%ds %%9d",0));
      Int_t nRuns=0;
    
      // print value for each label
      TObjString* label = 0x0;
      TIter nextLabel(labels);
      cout << "-------------------------------------------------" << endl;
      cout << "Run Number" << "\t Number of "<< ( (TObjString*) triggersB->At(kCMUS) )->GetString()<<" after Phys. Sel. " << endl ;  
      while ((label = static_cast<TObjString*>(nextLabel()))) {
	nRuns++;
	Int_t bin = (Int_t) label->GetUniqueID();
	printf(format.Data(), label->String().Data(), (Int_t) hBWithPS[centBin][kCMUS]->GetBinContent(bin));
      }
      printf("\n========== Total #runs = %d ==============\n",nRuns);
      printf("\n\n");
    
    
      cout << "-------------------------------------------------" << endl;
      cout << "Total statistic" << endl; 
      cout << " " << endl ; 
    
      cout << "-------------------------------------------------------------------" << endl;
      cout<<"Number of triggers w/o Phys. Sel./ w/ Phys. Sel (Phys. Sel. cut in %)"<<endl;
      for(Int_t i = 0; i < triggersB->GetEntriesFast()-1; i++){
	TString triggerNameB = ( (TObjString*) triggersB->At(i) )->GetString();
	TString triggerNameAC = ( (TObjString*) triggersAC->At(i) )->GetString();
	TString triggerNameE = ( (TObjString*) triggersE->At(i) )->GetString();
			
	Int_t	cutinpercent	=	0;
	printf("%10s %30s %30s\n",triggerNameB.Data(),triggerNameAC.Data(),triggerNameE.Data());
	if(NumOfBNoPS[centBin][i]) cutinpercent = (Int_t) ((Double_t)(NumOfBNoPS[centBin][i]-NumOfBWithPS[centBin][i])/(NumOfBNoPS[centBin][i])*100.);
	printf("%5.2e / %.2e (%d%%)", (Double_t) NumOfBNoPS[centBin][i],(Double_t) NumOfBWithPS[centBin][i],cutinpercent);
	cutinpercent = 0;
	if(NumOfACNoPS[centBin][i]) cutinpercent = (Int_t) ((Double_t)(NumOfACNoPS[centBin][i]-NumOfACWithPS[centBin][i])/(NumOfACNoPS[centBin][i])*100.);
	printf("%15.2e / %.2e (%d%%)", (Double_t)NumOfACNoPS[centBin][i],(Double_t)NumOfACWithPS[centBin][i],cutinpercent);
	cutinpercent = 0;
	if(NumOfENoPS[centBin][i]) cutinpercent = (Int_t) ((Double_t)(NumOfENoPS[centBin][i]-NumOfEWithPS[centBin][i])/(NumOfENoPS[centBin][i])*100.);
	printf("%15.2e  / %.2e (%d%%)\n", (Double_t)NumOfENoPS[centBin][i],(Double_t)NumOfEWithPS[centBin][i],cutinpercent);

      }
    } 
  }
  
  //temporary
  //    return;

  
  //--------------------------------------------- //
  //        monitor quantities run per run        //
  //--------------------------------------------- //
  Int_t const nMuonTriggerCase = 7;
  TH1F *hMuonTrigger[nMuonTriggerCase];
  TString muonTriggerName[nMuonTriggerCase] = {"Unlike-only","Like-only","Hpt-only","Unlike&Like","Unlike&Hpt","Like&Hpt","Unlike&Like&Hpt"};
  for ( Int_t i = 0; i < nMuonTriggerCase; i++ ) {
    hMuonTrigger[i] = new TH1F(Form("hMuonTrigger_%s", muonTriggerName[i].Data()),Form("Trigger %s per run",muonTriggerName[i].Data()),10000,1,10000);
  }

  TH1F* hTriggerCutVsRun[2], *hTriggerCutWidthVsRun[2];
  for ( Int_t ihisto=0; ihisto<2; ++ihisto ) {
    TString cutName = ( ihisto == 0 ) ? "Lpt" : "Hpt";
    hTriggerCutVsRun[ihisto] = new TH1F(Form("hTriggerCutVsRun%s", cutName.Data()), Form("Trigger %s cut per run", cutName.Data()), 10000,1,10000);
    hTriggerCutWidthVsRun[ihisto] = (TH1F*)hTriggerCutVsRun[ihisto]->Clone(Form("hTriggerCutWidthVsRun%s", cutName.Data()));
  }
  TF1* fitMatchTrig = new TF1("fitMatchTrig","[3] + [0] * ( 1. + TMath::Erf((x - [1]) / [2] / TMath::Sqrt(2.)))", 0.1, 6.);	
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
  TH1F* hClusterChargeMeanInCh[10], *hClusterChargeSigmaInCh[10];
  for (Int_t ich=0; ich<10; ich++){
    hClusterChargeMeanInCh[ich] = new TH1F(Form("hClusterChargeMeanInCh%d",ich+1), Form("averaged cluster charge -mean- per track in chamber %d",ich+1),10000,1,10000);
    hClusterChargeSigmaInCh[ich] = new TH1F(Form("hClusterChargeSigmaInCh%d",ich+1), Form("averaged cluster charge -sigma- per track in chamber %d",ich+1),10000,1,10000);
  }
  TH1F* hClusterSizeMeanInCh[10], *hClusterSizeSigmaInCh[10];
  for (Int_t ich=0; ich<10; ich++){
    hClusterSizeMeanInCh[ich] = new TH1F(Form("hClusterSizeMeanInCh%d",ich+1), Form("averaged cluster size -mean- per track in chamber %d",ich+1),10000,1,10000);
    hClusterSizeSigmaInCh[ich] = new TH1F(Form("hClusterSizeSigmaInCh%d",ich+1), Form("averaged cluster size -sigma- per track in chamber %d",ich+1),10000,1,10000);
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
  cout<< "Monitoring quantities run per run: loop over "<<runs->GetEntriesFast()<<" runs."<<endl;
  cout<<"  directory = "<<alienBaseDir.Data()<<endl;
  cout<<"============================================================"<<endl;
  // Loop over runs
  for ( Int_t irun=0; irun<runs->GetEntriesFast(); irun++ ) {
    
    TString run = ((TObjString*)runs->UncheckedAt(irun))->GetString();
    TString run2 = ((TObjString*)runs2->UncheckedAt(irun))->GetString();
    // get the file (or list of files) to be analyzed
    TString command;
    TGridResult *res = 0;
    TObjString *objs = 0;	
    
    //cout<<irun<<" "<<run<<endl;

    if(isAlienFile){
      command = Form("find %s/ %s/%s", alienBaseDir.Data(), run.Data(), QAFileName);
      res = gGrid->Command(command);
      if (!res) {
	Error("PlotMuonQA","no result for the command: %s",command.Data());
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
	Error("PlotMuonQA","turl/obj not found for the run %s... SKIPPING", run.Data());
	continue;
      }
      
      run = objs->GetString();
      run.ReplaceAll(Form("/%s",QAFileName), "");
      run.ReplaceAll(alienBaseDir, "");
      run.Remove(TString::kLeading, '/');
      run.Remove(TString::kLeading, '0');
      
      if ( ! selectRuns.Contains(run.Data()) ) continue;      
      
      // open the outfile for this run
      TFile *runFile = TFile::Open(objs->GetString());
      if (!runFile || ! runFile->IsOpen()) {
	Error("PlotMuonQA","failed to open file: %s", objs->GetName());
	continue;//return;
      }
      runFile->Cd("MUON_QA");
      
      // get interesting histos
      TObjArray* general1 = static_cast<TObjArray*>(runFile->FindObjectAny("general1"));
      TObjArray* general2 = static_cast<TObjArray*>(runFile->FindObjectAny("general2"));
      TObjArray* expert = static_cast<TObjArray*>(runFile->FindObjectAny("expert"));
      
      if (!general1 || !general2 || !expert){
	Error("PlotMuonQA","All objects not here !!! ===> Skipping...for %s",objs->GetName());		
	continue;
      }
      
      TH1* hNClustersPerTrack = static_cast<TH1*>(general1->FindObject("hNClustersPerTrack"));
      TH1* hNChamberHitPerTrack = static_cast<TH1*>(general1->FindObject("hNChamberHitPerTrack"));
      TH1* hChi2 = static_cast<TH1*>(general1->FindObject("hChi2"));
      TH1* hNClustersPerCh = static_cast<TH1*>(general2->FindObject("hNClustersPerCh"));
      TH1 *hClusterChargePerChMean =static_cast<TH1*>(general2->FindObject("hClusterChargePerChMean")); 
      TH1 *hClusterChargePerChSigma =static_cast<TH1*>(general2->FindObject("hClusterChargePerChSigma")); 
      TH1 *hClusterSizePerChMean =static_cast<TH1*>(general2->FindObject("hClusterSizePerChMean")); 
      TH1 *hClusterSizePerChSigma =static_cast<TH1*>(general2->FindObject("hClusterSizePerChSigma")); 
      TH1* hMuonTriggers = static_cast<TH1*>(general1->FindObject("hMuonTriggers"));
      for (Int_t ihisto=0; ihisto<nMuonTriggerCase; ihisto++){
	Int_t val = 0;
	if ( hMuonTriggers ) val = hMuonTriggers->GetBinContent(ihisto+2);
	if ( hMuonTrigger[ihisto]->GetSumw2N() == 0 ) hMuonTrigger[ihisto]->Sumw2();
	hMuonTrigger[ihisto]->SetBinContent(ibin,val);
	hMuonTrigger[ihisto]->SetBinError(ibin,TMath::Sqrt(val));
	hMuonTrigger[ihisto]->GetXaxis()->SetBinLabel(ibin,run.Data());
      }
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
	Warning("PlotMUONQA","File: %s has empty histograms !", objs->GetName());
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
	  hClusterChargeMeanInCh[ich]->SetBinContent(ibin,hClusterChargePerChMean->GetBinContent(ich+1));
          hClusterChargeMeanInCh[ich]->SetBinError(ibin,hClusterChargePerChMean->GetBinError(ich+1));
	  hClusterChargeSigmaInCh[ich]->SetBinContent(ibin,hClusterChargePerChSigma->GetBinContent(ich+1));
          hClusterChargeSigmaInCh[ich]->SetBinError(ibin,hClusterChargePerChSigma->GetBinError(ich+1));
	  hClusterSizeMeanInCh[ich]->SetBinContent(ibin,hClusterSizePerChMean->GetBinContent(ich+1));
          hClusterSizeMeanInCh[ich]->SetBinError(ibin,hClusterSizePerChMean->GetBinError(ich+1));
	  hClusterSizeSigmaInCh[ich]->SetBinContent(ibin,hClusterSizePerChSigma->GetBinContent(ich+1));
          hClusterSizeSigmaInCh[ich]->SetBinError(ibin,hClusterSizePerChSigma->GetBinError(ich+1));
	  hNClustersInCh[ich]->SetBinContent(ibin,hNClustersPerCh->GetBinContent(ich+1));
          hNClustersInCh[ich]->SetBinError(ibin,hNClustersPerCh->GetBinError(ich+1));
	  hClusterHitMapXInCh[ich]->SetBinContent(ibin,hClusterHitMapInCh[ich]->GetMean(1));
	  hClusterHitMapXInCh[ich]->SetBinError(ibin,hClusterHitMapInCh[ich]->GetMeanError(1));
	  hClusterHitMapYInCh[ich]->SetBinContent(ibin,hClusterHitMapInCh[ich]->GetMean(2));
	  hClusterHitMapYInCh[ich]->SetBinError(ibin,hClusterHitMapInCh[ich]->GetMeanError(2));	  
	}
      }
      
      // set labels
      for ( Int_t ihisto=0; ihisto<2; ++ihisto ) {
        hTriggerCutVsRun[ihisto]->GetXaxis()->SetBinLabel(ibin,run.Data());
        hTriggerCutWidthVsRun[ihisto]->GetXaxis()->SetBinLabel(ibin,run.Data());
      }
      hNClustersPerTrackVsRun_Mean->GetXaxis()->SetBinLabel(ibin, run.Data());
      hNClustersPerTrackVsRun_Sigma->GetXaxis()->SetBinLabel(ibin, run.Data());
      hNChamberHitPerTrack_Mean->GetXaxis()->SetBinLabel(ibin, run.Data());
      hNChamberHitPerTrack_Sigma->GetXaxis()->SetBinLabel(ibin, run.Data());
      hChi2_Mean->GetXaxis()->SetBinLabel(ibin, run.Data());
      hChi2_Sigma->GetXaxis()->SetBinLabel(ibin, run.Data());
      for (Int_t ich=0; ich<10; ich++){
	hNClustersInCh[ich]->GetXaxis()->SetBinLabel(ibin, run.Data());
	hClusterChargeMeanInCh[ich]->GetXaxis()->SetBinLabel(ibin, run.Data());
	hClusterChargeSigmaInCh[ich]->GetXaxis()->SetBinLabel(ibin, run.Data());
	hClusterSizeMeanInCh[ich]->GetXaxis()->SetBinLabel(ibin, run.Data());
	hClusterSizeSigmaInCh[ich]->GetXaxis()->SetBinLabel(ibin, run.Data());
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
    hClusterChargeMeanInCh[ich]->LabelsOption("a");
    hClusterChargeSigmaInCh[ich]->LabelsOption("a");
    hClusterSizeMeanInCh[ich]->LabelsOption("a");
    hClusterSizeSigmaInCh[ich]->LabelsOption("a");
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
  dirGlobal->cd();
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
  dirGlobal->cd();
  cNClustersPerCh->Write();

  //==================================================
  // Display average cluster charge per chamber
                                                                          
  TLegend *lClusterChargePerCh = new TLegend(0.92,0.45,0.99,0.99);
  TCanvas* cClusterChargePerCh = new TCanvas("cClustersChargePerCh","cClustersChargePerCh",1200,900);
  cClusterChargePerCh->SetRightMargin(0.1);
  cClusterChargePerCh->Divide(1,2);

  cClusterChargePerCh->cd(1);
  hClusterChargeMeanInCh[0]->SetStats(kFALSE);
  hClusterChargeMeanInCh[0]->GetXaxis()->SetRange(1,ibin-1);
  hClusterChargeMeanInCh[0]->GetXaxis()->SetNdivisions(1,kFALSE);
  //hClusterChargeInCh[0]->LabelsOption("u"); 
  hClusterChargeMeanInCh[0]->SetLabelSize(0.04);
  hClusterChargeMeanInCh[0]->SetTitle("Cluster charge mean (fC) per track in chamber i");
  hClusterChargeMeanInCh[0]->SetMaximum(150);
  hClusterChargeMeanInCh[0]->SetMinimum(30);
  for (Int_t ich=0; ich<10; ich++) {
    hClusterChargeMeanInCh[ich]->SetLineColor(ich+1+ich/9);
    hClusterChargeMeanInCh[ich]->SetLineWidth(2);
    if (ich == 0) hClusterChargeMeanInCh[ich]->Draw("e");
    else hClusterChargeMeanInCh[ich]->Draw("esame");
    lClusterChargePerCh->AddEntry(hClusterChargeMeanInCh[ich],Form("ch%d",ich+1),"PL");
  }
  lClusterChargePerCh->Draw("same");

  cClusterChargePerCh->cd(2);
  hClusterChargeSigmaInCh[0]->SetStats(kFALSE);
  hClusterChargeSigmaInCh[0]->GetXaxis()->SetRange(1,ibin-1);
  hClusterChargeSigmaInCh[0]->GetXaxis()->SetNdivisions(1,kFALSE);
  //hClusterChargeInCh[0]->LabelsOption("u");                                                                               
  hClusterChargeSigmaInCh[0]->SetLabelSize(0.04);
  hClusterChargeSigmaInCh[0]->SetTitle("Cluster charge sigma per track in chamber i");
  hClusterChargeSigmaInCh[0]->SetMaximum(250);
  hClusterChargeSigmaInCh[0]->SetMinimum(50);
  for (Int_t ich=0; ich<10; ich++) {
    hClusterChargeSigmaInCh[ich]->SetLineColor(ich+1+ich/9);
    hClusterChargeSigmaInCh[ich]->SetLineWidth(2);
    if (ich == 0) hClusterChargeSigmaInCh[ich]->Draw("e");
    else hClusterChargeSigmaInCh[ich]->Draw("esame");
  }
  lClusterChargePerCh->Draw("same");

  cClusterChargePerCh->Print(OutFileNamePDF.Data());
  dirGlobal->cd();
  cClusterChargePerCh->Write();

  //==================================================
  // Display average cluster size per chamber     

  TLegend *lClusterSizePerCh = new TLegend(0.92,0.45,0.99,0.99);
  TCanvas* cClusterSizePerCh = new TCanvas("cClustersSizePerCh","cClustersSizePerCh",1200,900);
  cClusterSizePerCh->SetRightMargin(0.1);
  cClusterSizePerCh->Divide(1,2);

  cClusterSizePerCh->cd(1);
  hClusterSizeMeanInCh[0]->SetStats(kFALSE);
  hClusterSizeMeanInCh[0]->GetXaxis()->SetRange(1,ibin-1);
  hClusterSizeMeanInCh[0]->GetXaxis()->SetNdivisions(1,kFALSE);

  hClusterSizeMeanInCh[0]->SetLabelSize(0.04);
  hClusterSizeMeanInCh[0]->SetTitle("Cluster size mean (npads) per track in chamber i");
  hClusterSizeMeanInCh[0]->SetMaximum(18);
  hClusterSizeMeanInCh[0]->SetMinimum(0);
  for (Int_t ich=0; ich<10; ich++) {
    hClusterSizeMeanInCh[ich]->SetLineColor(ich+1+ich/9);
    hClusterSizeMeanInCh[ich]->SetLineWidth(2);
    if (ich == 0) hClusterSizeMeanInCh[ich]->Draw("e");
    else hClusterSizeMeanInCh[ich]->Draw("esame");
    lClusterSizePerCh->AddEntry(hClusterSizeMeanInCh[ich],Form("ch%d",ich+1),"PL");
  }
  lClusterSizePerCh->Draw("same");

  cClusterSizePerCh->cd(2);
  hClusterSizeSigmaInCh[0]->SetStats(kFALSE);
  hClusterSizeSigmaInCh[0]->GetXaxis()->SetRange(1,ibin-1);
  hClusterSizeSigmaInCh[0]->GetXaxis()->SetNdivisions(1,kFALSE);

  hClusterSizeSigmaInCh[0]->SetLabelSize(0.04);
  hClusterSizeSigmaInCh[0]->SetTitle("Cluster size sigma (npads) per track in chamber i");
  hClusterSizeSigmaInCh[0]->SetMaximum(7);
  hClusterSizeSigmaInCh[0]->SetMinimum(0);
  for (Int_t ich=0; ich<10; ich++) {
    hClusterSizeSigmaInCh[ich]->SetLineColor(ich+1+ich/9);
    hClusterSizeSigmaInCh[ich]->SetLineWidth(2);
    if (ich == 0) hClusterSizeSigmaInCh[ich]->Draw("e");
    else hClusterSizeSigmaInCh[ich]->Draw("esame");
  }
  lClusterSizePerCh->Draw("same");

  cClusterSizePerCh->Print(OutFileNamePDF.Data());
  dirGlobal->cd();
  cClusterSizePerCh->Write();


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
  dirGlobal->cd();
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
  dirGlobal->cd();
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
    dirGlobal->cd();
    cLptHpt->Write();
 }


  //==================================================
  // Display muon trigger
  //==================================================
  
  TCanvas* cMuonTriggerUnlikeLike = new TCanvas("cMuonTriggerUnlikeLike","cMuonTriggerUnlikeLike",1200,900);
  TCanvas* cMuonTriggerHpt = new TCanvas("cMuonTriggerHpt","cMuonTriggerHpt",1200,900);
  cMuonTriggerHpt->Divide(1,2);
  
  Int_t const ntrig = 3;
  Int_t  const nhist = 3;
  TH1F *hRelMuonTrigger[ntrig], *hTot;
  
  //3 configurations (3x3 histos for 3 canvases)
  Int_t conf[ntrig][nhist]={{0,1,3},{0,2,4},{1,2,5}};
  TString sconf[ntrig][nhist]={{"hTriggerUnlikeOnly","hTriggerLikeOnly","hTriggerUnlikeAndLike"},{"hTriggerUnlikeOnly","hTriggerHptOnly","hTriggerUnlikeAndHpt"},{"hTriggerLikeOnly","hTriggerHptOnly","hTriggerLikeAndHpt"}};
  TString sname;
  
  for ( Int_t i=0; i < ntrig; i++ ) {
    for ( Int_t j=0; j < nhist; j++ ) {
      sname = sconf[i][j]; sname +=(i+1);
      hMuonTrigger[conf[i][j]]->GetXaxis()->SetRange(1,ibin-1);
      hRelMuonTrigger[j] = (TH1F*)hMuonTrigger[conf[i][j]]->Clone(sname);
      hRelMuonTrigger[j]->SetLineColor(j+1);
      hRelMuonTrigger[j]->SetStats(kFALSE);
      hRelMuonTrigger[j]->GetXaxis()->SetNdivisions(1,kFALSE);
      hRelMuonTrigger[j]->LabelsOption("a");                 
      hRelMuonTrigger[j]->SetLabelSize(0.02);
      //hRelMuonTrigger[j]->GetXaxis()->SetLabelSize(0.04);
      hRelMuonTrigger[j]->SetLineWidth(2);
      hRelMuonTrigger[j]->SetTitle("");
      
    }
    sname = "hTot";
    sname += (i+1);
    hTot = (TH1F*) hRelMuonTrigger[0]->Clone(sname);
    hTot->Add(hRelMuonTrigger[1]);
    hTot->Add(hRelMuonTrigger[2],-1);
    hRelMuonTrigger[0]->Add(hRelMuonTrigger[2],-1);
    hRelMuonTrigger[1]->Add(hRelMuonTrigger[2],-1);
    
    for(Int_t j=0; j < nhist; j++) hRelMuonTrigger[j]->Divide(hTot);
    
    if(i==0) cMuonTriggerUnlikeLike->cd();
    else cMuonTriggerHpt->cd(i);
    
    if (i==0) hRelMuonTrigger[i]->SetTitle("Relative muon triggers content");
    TLegend *leg = new TLegend(0.72,0.7,0.9,0.85);
    leg->SetBorderSize(1);
    
    for(Int_t j=0; j<nhist; j++){
     if(j==0){
       hRelMuonTrigger[j]->SetMaximum(1);
       hRelMuonTrigger[j]->SetMinimum(0);
       hRelMuonTrigger[j]->Draw("e");
     }
     else hRelMuonTrigger[j]->Draw("esame");
     sname = sconf[i][j];
     leg->AddEntry(hRelMuonTrigger[j],sname,"l");
    }
    leg->Draw("same");
  }
  cMuonTriggerUnlikeLike->Print(OutFileNamePDF.Data());
  dirGlobal->cd();
  cMuonTriggerUnlikeLike->Write();
  cMuonTriggerHpt->Print(OutFileNamePDF.Data());
  cMuonTriggerHpt->Write(); 
  
  // close the PDF file
 c1->Print(OutFileNamePDF_close.Data());
 rootFileOut->Close();
 
 delete runs;
 delete runs2;
 delete triggersB;
 delete triggersAC;
 delete triggersE;
 delete dirTrigger;
 delete dirCent;
 delete colorInd;

 return;
 
}

void SetCanvas(TCanvas *canvas, Int_t logy){
  
  if(!canvas) return;
  canvas->SetTopMargin(0.05);
  canvas->SetRightMargin(0.01);
  canvas->SetGridy(1);
  canvas->SetLogy(logy);
	
  return;
}

Bool_t IsTrigger(TObjArray *array, Int_t index, TString name){
	
  Bool_t process = kFALSE;
	
  if( !array) return process;
	
  TString triggerName = (( (TObjString*) array->At(index) )->GetString());
	
  if( triggerName.Contains(name) ) process = kTRUE;
	
  return process;
}

Bool_t IsTriggerSelectedForMuonPhysics(TObjArray *triggersB, Int_t k){
	
  if ( !triggersB ) return kFALSE;
  Bool_t selected;
  selected =  (IsTrigger(triggersB,k,"CMB") || IsTrigger(triggersB, k, "CPBI") || IsTrigger(triggersB, k, "CVHN") || IsTrigger(triggersB, k, "CVLN") || IsTrigger(triggersB, k, "CCENT") || IsTrigger(triggersB, k, "CSEMI") );
  if (!selected) selected =  (IsTrigger(triggersB,k,"CINT7") || IsTrigger(triggersB, k, "CMUSH7") ||IsTrigger(triggersB, k, "CMUS7") || IsTrigger(triggersB, k, "CMUL7") || IsTrigger(triggersB, k, "CMUU7") );

  return selected;
}


Bool_t IsHeavyIonCollision(AliCounterCollection *eventCounters){
	
  if(!eventCounters) return kFALSE;
	
  Double_t sum = eventCounters->GetSum("v0mult:low,int,high");
  Bool_t result =		kTRUE;
  if(sum<=0) result = kFALSE;
	
  cout<<" Collision type is set to ";
  if( result == kFALSE) cout<<"p-p"<<endl;
  else cout<<"heavy-ion"<<endl;
	
  return result;
}

TCanvas *ProcessCanvasTriggerContent(TObjArray *array, TH1 **hBNoPS, TH1 **hBWithPS, TString canvasName){
 
  if(!array || !hBNoPS || !hBWithPS) return 0x0;
	
  TString cName =  "c"; 
  cName += canvasName; 
  TCanvas *cTriggerContent = new TCanvas(canvasName,cName,1200,900);
  SetCanvas(cTriggerContent);
  cTriggerContent->cd();
 
  TLegend* legcTC = new TLegend(0.2,0.15,0.50,0.40);
  legcTC->SetHeader("Physics Selection");
  legcTC->AddEntry(".","applied :","");
 
  for(Int_t i = 0; i < array->GetEntriesFast(); i++){
    //skip the sum of all triggers
    if( i == (array->GetEntriesFast()-1) ) continue;
    hBNoPS[i]->SetLineStyle(2);
    if(i==0){
      hBNoPS[i]->SetMinimum(1e-3);
      hBNoPS[i]->Draw();
      hBWithPS[i]->Draw("same");
    }
    else{
      hBNoPS[i]->Draw("same");
      hBWithPS[i]->Draw("same");
    }
    legcTC->AddEntry(hBWithPS[i],(( (TObjString*) array->At(i) )->GetString()).Data(),"l");
  }
  legcTC->AddEntry(".","not applied :","");
	
  for(Int_t i = 0; i < array->GetEntriesFast(); i++){
    legcTC->AddEntry(hBNoPS[i],(( (TObjString*) array->At(i) )->GetString()).Data(),"l");	 
  }
	
  legcTC->Draw("same");
 
  return cTriggerContent;
}

TCanvas *ProcessCanvasRelativeTriggerContent(TObjArray *triggersB, TH1 **histo, TString canvasName, TArrayI *colorInd){
	
  if(!triggersB || !histo ) return 0x0;
	
  TString cName =  "c" ; 
  cName += canvasName;
  TCanvas *cRelativeTriggerContent = new TCanvas(canvasName,cName,1200,900);
  SetCanvas(cRelativeTriggerContent);
  cRelativeTriggerContent->cd();
	
  TH1** ratio = new TH1*[triggersB->GetEntriesFast()];
  TLegend* legcRTC = new TLegend(0.2,0.15,0.50,0.40);
  legcRTC->SetHeader("Physics Selection");
	
  TString hName, hTriggerName;
  Int_t indAllTrig = triggersB->GetEntriesFast()-1;

  for(Int_t i = 0; i < triggersB->GetEntriesFast()-1; i++){
    hName = "ratio";
    hName += ( (TObjString*) triggersB->At(i) )->GetString();
    ratio[i] = static_cast<TH1*> (histo[i]->Clone(hName));
    ratio[i]->Divide(histo[indAllTrig]);
    ratio[i]->SetLineWidth(2);
    ratio[i]->SetLineColor(colorInd->At(i));
    if(i==0){
      ratio[i]->SetMaximum(1.5);
      ratio[i]->SetMinimum(0.001);
      ratio[i]->SetLabelSize(0.04);
      ratio[i]->GetYaxis()->SetTitle("Relative trigger content"); 
      ratio[i]->Draw("E");
    }
    else{
      ratio[i]->Draw("ESAME");
    }
  }
	
  legcRTC->AddEntry(".","not applied :","");
  for(Int_t i = 0; i < triggersB->GetEntriesFast()-1; i++){
    legcRTC->AddEntry(ratio[i],(( (TObjString*) triggersB->At(i) )->GetString()).Data(),"l");
  }
  legcRTC->Draw("same");
	
  return cRelativeTriggerContent;
}

TCanvas *ProcessCanvasPhysSelCut(TObjArray *triggersB, TObjArray *triggersAC, TObjArray *triggersE, TH1 **hBNoPS, TH1 **hACNoPS, TH1 **hENoPS, TH1 **hBWithPS, TString canvasName, TArrayI *colorInd, Bool_t isHeavyIon){
	
  if(!triggersB || !triggersE || !triggersAC || !hBNoPS || !hACNoPS || !hENoPS || !hBWithPS) return 0x0;
	
  //Multiplicative factor for A,C and E triggers
  Float_t scale = 3; //3 for pp - 10 for PbPb
  if(isHeavyIon) scale = 10;
	
  TString cName = "c";
  cName += canvasName;
  TCanvas *c1 = new TCanvas(canvasName,cName,1200,900);
  SetCanvas(c1);
  c1->cd();
	 
  TH1** ratioB = new TH1*[triggersB->GetEntriesFast()], **ratioBNoPS = new TH1*[triggersB->GetEntriesFast()];
  TH1** ratioACNoPS = new TH1*[triggersB->GetEntriesFast()];
  TH1** ratioENoPS = new TH1*[triggersB->GetEntriesFast()];
  TLegend* legcRTC = new TLegend(0.2,0.15,0.50,0.40);
  legcRTC->SetHeader("Physics Selection");

  TString hName;
  for(Int_t i = 0; i < triggersB->GetEntriesFast()-1; i++){
		
    //scale A, C and E triggers
    hACNoPS[i]->Scale(scale);
    hENoPS[i]->Scale(scale);
		
    hName = "ratio";
    hName += ( (TObjString*) triggersB->At(i) )->GetString();
    ratioB[i] = static_cast<TH1*> (hBWithPS[i]->Clone(hName));
    ratioB[i]->Divide(hBNoPS[i]);
    ratioB[i]->SetLineWidth(2);
    ratioB[i]->SetLineColor(colorInd->At(i));
    hName = "ratioNoPS";
    hName += ( (TObjString*) triggersB->At(i) )->GetString();
    ratioBNoPS[i] = static_cast<TH1*> (hBNoPS[i]->Clone(hName));
    ratioBNoPS[i]->Divide(hBNoPS[i]);
    ratioBNoPS[i]->SetLineWidth(0);
    ratioBNoPS[i]->SetLineStyle(1);
    ratioBNoPS[i]->SetMarkerStyle(24+i);
    ratioBNoPS[i]->SetMarkerSize(1);
    ratioBNoPS[i]->SetLineColor(colorInd->At(i));
    ratioBNoPS[i]->SetMarkerColor(colorInd->At(i));
		 
    hName = "ratioACNoPS";
    hName += ( (TObjString*) triggersAC->At(i) )->GetString();
    ratioACNoPS[i] = static_cast<TH1*> (hACNoPS[i]->Clone(hName));
    if ( ratioACNoPS[i]->GetEntries() > 0 ) ratioACNoPS[i]->Divide(hBNoPS[i]);
    ratioACNoPS[i]->SetLineWidth(0);
    ratioACNoPS[i]->SetLineStyle(2);
    ratioACNoPS[i]->SetMarkerStyle(24+i);
    ratioACNoPS[i]->SetMarkerSize(1);
    ratioACNoPS[i]->SetLineColor(colorInd->At(i));
    ratioACNoPS[i]->SetMarkerColor(colorInd->At(i));
		
    hName = "ratioENoPS";
    hName += ( (TObjString*) triggersE->At(i) )->GetString();
    ratioENoPS[i] = static_cast<TH1*> (hENoPS[i]->Clone(hName));
    if ( ratioENoPS[i]->GetEntries() > 0 ) ratioENoPS[i]->Divide(hBNoPS[i]);
    ratioENoPS[i]->SetLineWidth(0);
    ratioENoPS[i]->SetLineStyle(3);
    ratioENoPS[i]->SetMarkerStyle(24+i);
    ratioENoPS[i]->SetMarkerSize(1);
    ratioENoPS[i]->SetLineColor(colorInd->At(i));
    ratioENoPS[i]->SetMarkerColor(colorInd->At(i));
		 
	 
    if(i==0){
      ratioB[i]->SetMaximum(1.5);
      ratioB[i]->SetMinimum(0.05);
      ratioB[i]->SetLabelSize(0.02);
      ratioB[i]->GetYaxis()->SetTitle("Accepted / All from Phys. Sel."); 
      ratioB[i]->SetTitle("Phys. Sel. for all selected triggers"); 
      ratioB[i]->Draw("E");
      //ratioBNoPS[i]->Draw("EPSAME");
      ratioACNoPS[i]->Draw("EPSAME");
      ratioENoPS[i]->Draw("EPSAME");
    }
    else{
      ratioB[i]->Draw("ESAME");
      //ratioBNoPS[i]->Draw("EPSAME");
      ratioACNoPS[i]->Draw("EPSAME");
      ratioENoPS[i]->Draw("EPSAME");
    }
  }
	 
  legcRTC->AddEntry(".","applied :","");
  for(Int_t i = 0; i < triggersB->GetEntriesFast()-1; i++){
    legcRTC->AddEntry(ratioB[i],(( (TObjString*) triggersB->At(i) )->GetString()).Data(),"l");
  }
  legcRTC->AddEntry(".","not applied :","");
  for(Int_t i = 0; i < triggersB->GetEntriesFast()-1; i++){
    //legcRTC->AddEntry(ratioBNoPS[i],(( (TObjString*) triggersB->At(i) )->GetString()).Data(),"pl");
    TString textLegend = ( (TObjString*) triggersAC->At(i) )->GetString();
    if( textLegend.CompareTo(" ") ){
      textLegend += " x";
      textLegend += scale;
      legcRTC->AddEntry(ratioACNoPS[i],textLegend.Data(),"pl");
    }
    textLegend = ( (TObjString*) triggersE->At(i) )->GetString();
    if( textLegend.CompareTo(" ") ){
      //cout<<"trigger="<<textLegend.Data()<<"-"<<endl;
      textLegend += " x";
      textLegend += scale;
      legcRTC->AddEntry(ratioENoPS[i],textLegend.Data(),"pl");
    }
  }
  legcRTC->Draw("same");
	
	
  return c1;
}	

TCanvas *ProcessCanvasPhysSelCutCentrality(TObjArray *triggersB, TObjArray *triggersAC, TObjArray *triggersE, TH1 **hBNoPSCent, TH1 **hACNoPSCent, TH1 **hENoPSCent, TH1 **hBWithPSCent, Int_t k, TString canvasName, TArrayI *colorInd, TString *legendHeader, Bool_t isHeavyIon){
	
  if(!triggersB || !triggersE || !triggersAC || !hBNoPSCent || !hACNoPSCent || !hENoPSCent || !hBWithPSCent || !legendHeader) return 0x0;
	
  //Multiplicative factor for A,C and E triggers
  Float_t scale = 3; //3 for pp - 10 for PbPb
  if(isHeavyIon) scale = 10;
	
  TString cName = "c";
  cName += canvasName;
  TCanvas *c1 = new TCanvas(canvasName,cName,1200,900);
  SetCanvas(c1);
  c1->cd();
	
  Int_t const centBinMax =3;
  TH1* ratioB[centBinMax], *ratioBNoPS[centBinMax];
  TH1* ratioACNoPS[centBinMax];
  TH1* ratioENoPS[centBinMax];
  TLegend* legcRTC = new TLegend(0.2,0.15,0.50,0.40);
  legcRTC->SetHeader("Physics Selection");
	
  TString hName;
	
  Float_t yMin = 0.05, yMax = 2;
	
  for(Int_t centBin = 0; centBin < centBinMax; centBin++){
    //scale A, C and E triggers
    hACNoPSCent[centBin]->Scale(scale);
    hENoPSCent[centBin]->Scale(scale);
		
    hName = "ratio";
    hName += ( (TObjString*) triggersB->At(k) )->GetString();
    ratioB[centBin] = static_cast<TH1*> (hBWithPSCent[centBin]->Clone(hName));
    ratioB[centBin]->Divide(hBNoPSCent[centBin]);
    ratioB[centBin]->SetLineWidth(2);
    ratioB[centBin]->SetLineColor(colorInd->At(centBin));
    hName = "ratioNoPS";
    hName += ( (TObjString*) triggersB->At(k) )->GetString();
    ratioBNoPS[centBin] = static_cast<TH1*> (hBNoPSCent[centBin]->Clone(hName));
    ratioBNoPS[centBin]->Divide(hBNoPSCent[centBin]);
    ratioBNoPS[centBin]->SetLineWidth(0);
    ratioBNoPS[centBin]->SetLineStyle(1);
    ratioBNoPS[centBin]->SetMarkerStyle(24+centBin);
    ratioBNoPS[centBin]->SetMarkerSize(1);
    ratioBNoPS[centBin]->SetLineColor(colorInd->At(centBin));
    ratioBNoPS[centBin]->SetMarkerColor(colorInd->At(centBin));
		
    hName = "ratioACNoPS";
    hName += ( (TObjString*) triggersAC->At(k) )->GetString();
    ratioACNoPS[centBin] = static_cast<TH1*> (hACNoPSCent[centBin]->Clone(hName));
    if ( ratioACNoPS[centBin]->GetEntries() > 0 )  ratioACNoPS[centBin]->Divide(hBNoPSCent[centBin]);
    ratioACNoPS[centBin]->SetLineWidth(0);
    ratioACNoPS[centBin]->SetLineStyle(2);
    ratioACNoPS[centBin]->SetMarkerStyle(24+centBin);
    ratioACNoPS[centBin]->SetMarkerSize(1);
    ratioACNoPS[centBin]->SetLineColor(colorInd->At(centBin));
    ratioACNoPS[centBin]->SetMarkerColor(colorInd->At(centBin));
		
		
    hName = "ratioENoPS";
    hName += ( (TObjString*) triggersE->At(k) )->GetString();
    ratioENoPS[centBin] = static_cast<TH1*> (hENoPSCent[centBin]->Clone(hName));
    if ( ratioENoPS[centBin]->GetEntries() > 0 ) ratioENoPS[centBin]->Divide(hBNoPSCent[centBin]);
    ratioENoPS[centBin]->SetLineWidth(0);
    ratioENoPS[centBin]->SetLineStyle(3);
    ratioENoPS[centBin]->SetMarkerStyle(24+centBin);
    ratioENoPS[centBin]->SetMarkerSize(1);
    ratioENoPS[centBin]->SetLineColor(colorInd->At(centBin));
    ratioENoPS[centBin]->SetMarkerColor(colorInd->At(centBin));
		
		
    if(centBin==0){
      ratioB[centBin]->SetMaximum(yMax);
      ratioB[centBin]->SetMinimum(yMin);
      ratioB[centBin]->SetLabelSize(0.02);
      ratioB[centBin]->GetYaxis()->SetTitle("Accepted / All from Phys. Sel.");
      TString sTitle = "for ", sTitle2 = (( (TObjString*) triggersB->At(k) )->GetString()).Data();
      if ( !sTitle2.IsNull() ) sTitle += sTitle2;
      else sTitle = "";
      ratioB[centBin]->SetTitle(Form("Phys. Sel. %s - Multiplicity from V0 amplitude",sTitle.Data()));
      ratioB[centBin]->Draw("E");
      //ratioBNoPS[centBin]->Draw("EPSAME");
      ratioACNoPS[centBin]->Draw("EPSAME");
      ratioENoPS[centBin]->Draw("EPSAME");
    }
    else{
      ratioB[centBin]->Draw("ESAME");
      //ratioBNoPS[centBin]->Draw("EPSAME");
      ratioACNoPS[centBin]->Draw("EPSAME");
      ratioENoPS[centBin]->Draw("EPSAME");
    }
  }
	
  legcRTC->AddEntry(".","applied :","");
  for(Int_t centBin = 0; centBin < centBinMax; centBin++){
    legcRTC->AddEntry(ratioB[centBin],(legendHeader[centBin]).Data(),"l");
  }
  legcRTC->AddEntry(".","not applied :","");
  for(Int_t centBin = 0; centBin < centBinMax; centBin++){
    TString textLegend = ( (TObjString*) triggersAC->At(k) )->GetString();
    if( textLegend.CompareTo(" ") && ratioACNoPS[centBin]->GetMaximum() > yMin ){
      textLegend += " x";
      textLegend += scale;
      legcRTC->AddEntry(ratioACNoPS[centBin],textLegend.Data(),"pl");
    }
    textLegend = ( (TObjString*) triggersE->At(k) )->GetString();
    if( textLegend.CompareTo(" ") && ratioENoPS[centBin]->GetMaximum() > yMin){
      //cout<<"trigger="<<textLegend.Data()<<"-"<<endl;
      textLegend += " x";
      textLegend += scale;
      legcRTC->AddEntry(ratioENoPS[centBin],textLegend.Data(),"pl");
    }
  }
  legcRTC->Draw("same");
	
	
  return c1;
}	

TCanvas *ProcessCanvasCentralityPercentile(TObjArray *triggersB, TH1 **hBNoPSCent, TH1 **hBWithPSCent, Int_t k, TString canvasName, TArrayI *colorInd, TString *legendHeader){
	
  if(!triggersB || !hBNoPSCent || !hBWithPSCent || !legendHeader) return 0x0;
		
  TString cName = "c";
  cName += canvasName;
  TCanvas *c1 = new TCanvas(canvasName,cName,1200,900);
  SetCanvas(c1,0);
  c1->cd();
	
  Int_t const centBinMax =2;
  TH1* ratioB[centBinMax], *ratioBNoPS[centBinMax];
  TLegend* legcRTC = new TLegend(0.2,0.15,0.50,0.40);
  legcRTC->SetHeader("Physics Selection");
	
  TString hName;
	
  Float_t yMin = 0., yMax = 0.3;
	
  for(Int_t centBin = 0; centBin < centBinMax; centBin++){
			
    hName = "ratio";
    hName += ( (TObjString*) triggersB->At(k) )->GetString();
    ratioB[centBin] = static_cast<TH1*> (hBWithPSCent[centBin+1]->Clone(hName));
    ratioB[centBin]->Divide(hBWithPSCent[0]);
    ratioB[centBin]->Scale(0.8);
    ratioB[centBin]->SetLineWidth(2);
    ratioB[centBin]->SetLineColor(colorInd->At(centBin+1));
    hName = "ratioNoPS";
    hName += ( (TObjString*) triggersB->At(k) )->GetString();
    ratioBNoPS[centBin] = static_cast<TH1*> (hBNoPSCent[centBin+1]->Clone(hName));
    ratioBNoPS[centBin]->Divide(hBNoPSCent[0]);
    ratioBNoPS[centBin]->Scale(0.8);
    ratioBNoPS[centBin]->SetLineWidth(0);
    ratioBNoPS[centBin]->SetLineStyle(1);
    ratioBNoPS[centBin]->SetMarkerStyle(24+centBin+1);
    ratioBNoPS[centBin]->SetMarkerSize(1);
    ratioBNoPS[centBin]->SetLineColor(colorInd->At(centBin+1));
    ratioBNoPS[centBin]->SetMarkerColor(colorInd->At(centBin+1));
		
    if(centBin==0){
      ratioB[centBin]->SetMaximum(yMax);
      ratioB[centBin]->SetMinimum(yMin);
      ratioB[centBin]->SetLabelSize(0.02);
      ratioB[centBin]->GetYaxis()->SetTitle("Centrality percentile check"); 
      TString sTitle = "for ", sTitle2 = (( (TObjString*) triggersB->At(k) )->GetString()).Data();
      if ( !sTitle2.IsNull() ) sTitle += sTitle2;
      else sTitle = "";
      ratioB[centBin]->SetTitle(Form("Centrality percentile check %s - Multiplicity from V0 amplitude",sTitle.Data()));
      ratioB[centBin]->Draw("E");
      ratioBNoPS[centBin]->Draw("EPSAME");
    }
    else{
      ratioB[centBin]->Draw("ESAME");
      ratioBNoPS[centBin]->Draw("EPSAME");
    }
  }
	
  legcRTC->AddEntry(".","applied :","");
  for(Int_t centBin = 0; centBin < centBinMax; centBin++){
    legcRTC->AddEntry(ratioB[centBin],(legendHeader[centBin+1]).Data(),"l");
  }
  legcRTC->AddEntry(".","not applied :","");
  for(Int_t centBin = 0; centBin < centBinMax; centBin++){
    legcRTC->AddEntry(ratioBNoPS[centBin],(legendHeader[centBin+1]).Data(),"l");
  }
  legcRTC->Draw("same");
	
	
  return c1;
}	

TCanvas *ProcessCanvasTracksoverTrigger(TObjArray *triggersB, TH1 **hB, TH1 **hTrackerB, TH1 **hTriggerB, TH1 **hMatchedB, TH1 **hAllTracksB, Int_t indTrigger, TString canvasName, TString legendHeader){
	 
  if(!triggersB || !hB || !hTrackerB || !hTriggerB || !hMatchedB || !hAllTracksB || indTrigger<0 ) return 0x0;
	
  TH1 *hTrackerPerB, *hTriggerPerB, *hMatchedPerB, *hAllTracksPerB;
		 
  TString hName, hNameBase;
  hNameBase =( (TObjString*) triggersB->At(indTrigger) )->GetString();
		
  hName = "hTrackerPer";
  hName += hNameBase;
  hTrackerPerB = static_cast<TH1*>(hTrackerB[indTrigger]->Clone(hName));
  if ( hTrackerPerB->GetEntries() > 0 )  hTrackerPerB->Divide(hB[indTrigger]);
  hTrackerPerB->SetLineWidth(2);
  hTrackerPerB->SetLineColor(kRed);
	 
  hName = "hTriggerPer";
  hName += hNameBase;
  hTriggerPerB = static_cast<TH1*>(hTriggerB[indTrigger]->Clone(hName));
  if ( hTriggerPerB->GetEntries() > 0 ) hTriggerPerB->Divide(hB[indTrigger]);
  hTriggerPerB->SetLineWidth(2);
  hTriggerPerB->SetLineColor(kBlue);
	 
  hName = "hMatchedPer";
  hName += hNameBase;
  hMatchedPerB = static_cast<TH1*>(hMatchedB[indTrigger]->Clone(hName));
  if ( hMatchedPerB->GetEntries() > 0 ) hMatchedPerB->Divide(hB[indTrigger]);
  hMatchedPerB->SetLineWidth(2);
  hMatchedPerB->SetLineColor(kViolet);
	 
  hName = "hAllTracksPer";
  hName += hNameBase;
  hAllTracksPerB = static_cast<TH1*>(hAllTracksB[indTrigger]->Clone(hName));
  if ( hAllTracksPerB->GetEntries() > 0 ) hAllTracksPerB->Divide(hB[indTrigger]);
  hAllTracksPerB->SetLineWidth(3);
  hAllTracksPerB->SetLineColor(kBlack);
  hAllTracksPerB->SetTitle(Form("Number of Tracks /%s %s",hNameBase.Data(),legendHeader.Data()));
  hAllTracksPerB->SetMinimum(0.0001);
  hAllTracksPerB->SetLabelSize(0.02);
	

  TString cName = "c";
  cName += canvasName;
  hNameBase =( (TObjString*) triggersB->At(indTrigger) )->GetString();
  cName += hNameBase;	
  TCanvas *cRatioTrackTypesB = new TCanvas(canvasName,cName,1200,900);
  SetCanvas(cRatioTrackTypesB,0);
  cRatioTrackTypesB->cd();
	
  TLegend* legcTTCINT1B; 
	 	 
  hAllTracksPerB->Draw("E");
  hTrackerPerB->Draw("Esame");
  hMatchedPerB->Draw("Esame");
  hTriggerPerB->Draw("Esame");
	 
  legcTTCINT1B = new TLegend(0.70,0.5,0.90,0.70);
  legcTTCINT1B->AddEntry(hAllTracksPerB,"All tracks","l");
  legcTTCINT1B->AddEntry(hTrackerPerB,"Tracking (only) tracks","l");
  legcTTCINT1B->AddEntry(hMatchedPerB,"Matched tracks","l");
  legcTTCINT1B->AddEntry(hTriggerPerB,"Trigger (only) tracks","l");
  legcTTCINT1B->Draw("same");
	
  return cRatioTrackTypesB;
	
}


TCanvas *ProcessCanvasTrackMultB(TObjArray *triggersB, TH1 **hB, TH1 **hTrackerB, TH1 **hTriggerB, TH1 **hMatchedB, TH1 **hAllTracksB, Int_t indTrigger, TString canvasName,TString legendHeader){
	
  if(!triggersB || !hB || !hTrackerB || !hTriggerB || !hMatchedB || !hAllTracksB || indTrigger<0 ) return 0x0;
	
  TString cName = "c";
  cName += canvasName;
  TString hNameBase =( (TObjString*) triggersB->At(indTrigger) )->GetString();
  cName += hNameBase;	
  TCanvas *cTrackMultB = new TCanvas(canvasName,cName,1200,900);
  SetCanvas(cTrackMultB,0);

  cTrackMultB->Divide(1,2);
  cTrackMultB->cd(1);
	
	
  TH1* hSumTriggerOverB, *hSumTrackerOverB; 

  TString hName; 

  hName = Form("hSumTriggerOver%s",hNameBase.Data());
  hSumTriggerOverB = static_cast<TH1*>(hTriggerB[indTrigger]->Clone(hName));
  hSumTriggerOverB->Add(hMatchedB[indTrigger]);
  hSumTriggerOverB->Divide(hB[indTrigger]);
    
  hName = Form("Sum of trigger tracks (matched+trigger-only) / # events in %s %s",hNameBase.Data(),legendHeader.Data());
  hSumTriggerOverB->SetTitle(hName);
  hSumTriggerOverB->SetLabelSize(0.04);
  hSumTriggerOverB->SetLineWidth(2);
  hSumTriggerOverB->SetLineColor(kBlue);
  hName = Form("hSumTrackerOver%s",hNameBase.Data());
  hSumTrackerOverB = static_cast<TH1*>(hTrackerB[indTrigger]->Clone(hName));
  hSumTrackerOverB->Add(hMatchedB[indTrigger]);
  hSumTrackerOverB->Divide(hB[indTrigger]);
  hName = Form("Sum of tracker tracks (matched+tracker-only) / # events in %s %s",hNameBase.Data(),legendHeader.Data());
  hSumTrackerOverB->SetTitle(hName);
  //hSumTrackerOverCINT1B->LabelsOption("u");
  hSumTrackerOverB->SetLabelSize(0.04);
  hSumTrackerOverB->SetLineWidth(2);
  hSumTrackerOverB->SetLineColor(kBlue);
		
	
	
  hSumTriggerOverB->Draw("e");
  cTrackMultB->cd(2);
  hSumTrackerOverB->Draw("e");
	
  return cTrackMultB;
	
}

TCanvas *ProcessCanvasRatioTrackB(TObjArray *triggersB, TH1 **hB, TH1 **hTrackerB, TH1 **hTriggerB, TH1 **hMatchedB, TH1 **hAllTracksB, Int_t indTrigger, TString canvasName,TString legendHeader){
	
  if(!triggersB || !hB || !hTrackerB || !hTriggerB || !hMatchedB || !hAllTracksB || indTrigger<0 ) return 0x0;
	
  TString cName = "c";
  cName += canvasName;
  TString hNameBase =( (TObjString*) triggersB->At(indTrigger) )->GetString();
  cName += hNameBase;	
  TCanvas *cRatioTrackB = new TCanvas(canvasName,cName,1200,900);
  SetCanvas(cRatioTrackB,0);
	
  TH1* hTrackerOverTriggerB, *hMatchedOverTriggerB, *hMatchedOverTrackerB;	
  
  TString hName = Form("hTrackerOverTrigger%s",hNameBase.Data());
  hTrackerOverTriggerB = static_cast<TH1*>(hTrackerB[indTrigger]->Clone(hName));
  hTrackerOverTriggerB->Divide(hTriggerB[indTrigger]);
  hName = Form("# tracker tracks / # trigger tracks in %s %s",hNameBase.Data(),legendHeader.Data());
  hTrackerOverTriggerB->SetTitle(hName);
  //hTrackerOverTriggerCINT1B->LabelsOption("u");
  hTrackerOverTriggerB->SetLabelSize(0.02);
  hTrackerOverTriggerB->SetLineWidth(2);
  hTrackerOverTriggerB->SetLineColor(kBlue);
    
  hName = Form("hMatchedOverTrigger%s",hNameBase.Data());	
  hMatchedOverTriggerB = static_cast<TH1*>(hMatchedB[indTrigger]->Clone(hName));
  hMatchedOverTriggerB->Divide(hTriggerB[indTrigger]);
   
  hName = Form("# matched tracks / # trigger tracks in %s %s",hNameBase.Data(),legendHeader.Data());
  hMatchedOverTriggerB->SetTitle(hName);
  //hMatchedOverTriggerCINT1B->LabelsOption("u");
  hMatchedOverTriggerB->SetLabelSize(0.02);
  hMatchedOverTriggerB->SetLineWidth(2);
  hMatchedOverTriggerB->SetLineColor(kBlue);
    
  hName = Form("hMatchedOverTracker%s",hNameBase.Data());
  hMatchedOverTrackerB = static_cast<TH1*>(hMatchedB[indTrigger]->Clone(hName));
  hMatchedOverTrackerB->Divide(hTrackerB[indTrigger]);
  hName = Form("# matched tracks / # tracker tracks in %s %s",hNameBase.Data(),legendHeader.Data());
  hMatchedOverTrackerB->SetTitle(hName);
  //hMatchedOverTrackerCINT1B->LabelsOption("u");
  hMatchedOverTrackerB->SetLabelSize(0.02);
  hMatchedOverTrackerB->SetLineWidth(2);
  hMatchedOverTrackerB->SetLineColor(kBlue);
  
	
  cRatioTrackB->Divide(1,3);
  cRatioTrackB->cd(1);
  hTrackerOverTriggerB->Draw("e");	
  cRatioTrackB->cd(2);
  hMatchedOverTriggerB->Draw("e");	
  cRatioTrackB->cd(3);
  hMatchedOverTrackerB->Draw("e");	
    
  return cRatioTrackB;
	
}

TCanvas *ProcessCanvasAsymMatched(TObjArray *triggersB, TH1 **hPosMatchedB, TH1 **hNegMatchedB, TH1 **hAllMatchedB, Int_t indTrigger, TString canvasName,TString legendHeader){
	
  if(!triggersB || !hPosMatchedB || !hNegMatchedB || !hAllMatchedB || indTrigger<0 ) return 0x0;

  TString hName, hNameBase = (( (TObjString*) triggersB->At(indTrigger) )->GetString());
	
  TString cName =	"c";	
  cName += canvasName;
  cName += hNameBase;	
  TCanvas *cAsymMatched = new TCanvas(canvasName.Data(),cName,1200,900);
  SetCanvas(cAsymMatched,0);
  cAsymMatched->cd();
	
	
  TH1 *hDiffMatchedCMUS1B= static_cast<TH1*>(hPosMatchedB[indTrigger]->Clone("hDiffMatchedCMUS1B"));
  hDiffMatchedCMUS1B->Add(hNegMatchedB[indTrigger],-1);
  if ( hDiffMatchedCMUS1B->GetSumw2N() == 0 ) hDiffMatchedCMUS1B->Sumw2();
	 
  TH1 *hAsymMatchedCMUS1B= static_cast<TH1*>(hDiffMatchedCMUS1B->Clone("hAsymMatchedCMUS1B"));
  hAsymMatchedCMUS1B->Divide(hAllMatchedB[indTrigger]);
  hAsymMatchedCMUS1B->SetLineColor(kRed);
  hAsymMatchedCMUS1B->SetLineWidth(2);
  hAsymMatchedCMUS1B->SetMinimum(-0.3);
  hAsymMatchedCMUS1B->SetMaximum(0.3);
  hAsymMatchedCMUS1B->SetLabelSize(0.02);
  hName = Form("Matched tracks charge asymmetry for %s with acc. cuts %s",hNameBase.Data(),legendHeader.Data());
  hAsymMatchedCMUS1B->SetTitle(hName);
	 
  hAsymMatchedCMUS1B->GetYaxis()->SetTitle("Charged tracks asymmetry");  
  hAsymMatchedCMUS1B->Draw("EH");
	
  return cAsymMatched;
	
}

TCanvas *ProcessCanvasHighPtMuons(TObjArray *triggersB, TH1 **hB, TH1 **hMatchedLowPtB, TH1 **hMatchedHighPtB, Int_t indTrigger, TString canvasName, TString legendHeader){
	
  if(!triggersB || !hB || !hMatchedLowPtB || !hMatchedHighPtB || indTrigger<0 ) return 0x0;
	
  TString hName, hNameBase = (( (TObjString*) triggersB->At(indTrigger) )->GetString());
	
  TString cName =	"c";	
  cName += canvasName;
  cName += hNameBase;	
  TCanvas *cHighPtMuons = new TCanvas(canvasName.Data(),cName,1200,900);
  SetCanvas(cHighPtMuons,0);
  cHighPtMuons->cd();
	
  TLegend* legcHPM;
	
  TH1* hMatchedLowPtPerB, *hMatchedHighPtPerB;
  hName = Form("hMatchedLowPtPer%s ",hNameBase.Data());
  hMatchedLowPtPerB = static_cast<TH1*> (hMatchedLowPtB[indTrigger]->Clone(hName));
  if ( hMatchedLowPtPerB->GetSumw2N() == 0 ) hMatchedLowPtPerB->Sumw2();
  hMatchedLowPtPerB->Divide(hB[indTrigger]);
  hMatchedLowPtPerB->SetLineWidth(2);
  hMatchedLowPtPerB->SetLineColor(kBlue);
  hMatchedLowPtPerB->SetTitle("");
  hName = Form("Ratio per %s ",hNameBase.Data());
  hMatchedLowPtPerB->GetYaxis()->SetTitle(hName);
  //hMatchedLowPtPerB->SetMaximum(0.15);
  hMatchedLowPtPerB->SetMinimum(0.0001);
  hMatchedLowPtPerB->SetLabelSize(0.02);
	 
  hName = Form("hMatchedHighPtPer%s ",hNameBase.Data());
  hMatchedHighPtPerB = static_cast<TH1*> (hMatchedHighPtB[indTrigger]->Clone(hName));
  if ( hMatchedHighPtPerB->GetSumw2N() == 0 ) hMatchedHighPtPerB->Sumw2();
  hMatchedHighPtPerB->Divide(hB[indTrigger]);
  hMatchedHighPtPerB->SetLineWidth(2);
  hMatchedHighPtPerB->SetLineColor(kRed);
	 	 
  hName = Form("Number of matched track per %s (include Vtx and R_{Abs} cuts) %s",hNameBase.Data(),legendHeader.Data());
  hMatchedLowPtPerB->SetTitle(hName);
  hMatchedLowPtPerB->Draw("E");
  hMatchedHighPtPerB->Draw("Esame");
	 
  legcHPM = new TLegend(0.60,0.45,0.98,0.65);
  //legcHPM->SetHeader(hName);
  legcHPM->AddEntry(".","Physics selection applied :","");	
  legcHPM->AddEntry(hMatchedLowPtPerB," p_{T} > 1 GeV/c ","l");
  legcHPM->AddEntry(hMatchedHighPtPerB," p_{T} >  2 GeV/c ","l");
  legcHPM->Draw("same");
	
  return cHighPtMuons;
}

TCanvas *ProcessCanvasBeamGasMatched(TObjArray *triggersB, TH1 **hBeamGasMatchedB, TH1 **hBeamGasMatchedHighPtB, TH1 **hAllMatchedB, TH1 **hMatchedHighPtB, Int_t indTrigger, TString canvasName,TString legendHeader){
	
  if(!triggersB || !hBeamGasMatchedB || !hBeamGasMatchedHighPtB || !hAllMatchedB || indTrigger<0 || !hMatchedHighPtB )
    return 0x0;

  TString hName, hNameBase = (( (TObjString*) triggersB->At(indTrigger) )->GetString());
	
  TString cName = "c";	
  cName += canvasName;
  cName += hNameBase;	
  TCanvas *cBeamGasMatched = new TCanvas(canvasName.Data(),cName,1200,900);
  SetCanvas(cBeamGasMatched,0);
  cBeamGasMatched->cd();
		
  hName = Form("hBeamGasMatchedPer%s ",hNameBase.Data());
  TH1 *hBeamGasMatchedCMUS1B= static_cast<TH1*>(hBeamGasMatchedB[indTrigger]->Clone(hName));
  hBeamGasMatchedCMUS1B->Divide(hAllMatchedB[indTrigger]);
  hBeamGasMatchedCMUS1B->SetLineColor(kBlack);
  hBeamGasMatchedCMUS1B->SetLineWidth(2);
  hBeamGasMatchedCMUS1B->SetMinimum(0.0);
  hBeamGasMatchedCMUS1B->SetMaximum(1.1);
  hBeamGasMatchedCMUS1B->SetLabelSize(0.02);
  
  hName = Form("hBeamGasMatchedHightPtPer%s ",hNameBase.Data());
  TH1 *hBeamGasMatchedHighPtCMUS1B= static_cast<TH1*>(hBeamGasMatchedHighPtB[indTrigger]->Clone(hName));
  hBeamGasMatchedHighPtCMUS1B->Divide(hMatchedHighPtB[indTrigger]);
  hBeamGasMatchedHighPtCMUS1B->SetLineColor(kRed);
  hBeamGasMatchedHighPtCMUS1B->SetLineWidth(2);

  hName = Form("Identified beam-gas tracks (pxDCA cuts) in matched tracks for %s",hNameBase.Data());
  if(!legendHeader.IsNull()) hName += Form(" %s",legendHeader.Data());

  hBeamGasMatchedCMUS1B->SetTitle(hName);
	 
  hBeamGasMatchedCMUS1B->GetYaxis()->SetTitle("Relative beam-gas tracks");  
  hBeamGasMatchedCMUS1B->Draw("EH");
  hBeamGasMatchedHighPtCMUS1B->Draw("EHsame");

  TLegend *leg = new TLegend(0.60,0.45,0.98,0.65);
  leg->AddEntry(".","Physics selection applied :","");	
  leg->AddEntry(hBeamGasMatchedCMUS1B," All p_{T}","l");
  leg->AddEntry(hBeamGasMatchedHighPtCMUS1B," p_{T} >  2 GeV/c ","l");
  leg->Draw("same");
	
  return cBeamGasMatched;
	
}

TH1* ProcessHisto( AliCounterCollection* counter, TString hVariable, TString hSelection, TString hName, TString xName, TString yName, Int_t color){
  
  
  TH1* h1 = 0x0;
  if( !counter ) return h1;

  //cout<<"ProcessHisto selection "<<hSelection<<endl;
	
  if ( !hSelection.Contains("trigger: /") && !hSelection.Contains("trigger:/") ) h1 = (TH1*) counter->Draw(hVariable,hSelection);
  //cout<<"ProcessHisto selection2 "<<h1<<endl;
  if ( !h1 ) h1 = new TH1D(hName,"",10,0,10);
  else {
    if ( h1->GetSumw2N() == 0 ) h1->Sumw2();
    h1->LabelsOption("a");
    if(hName.Sizeof()>1) h1->SetName(hName);
    if(xName.Sizeof()>1) h1->GetXaxis()->SetTitle(xName);
    if(yName.Sizeof()>1) h1->GetYaxis()->SetTitle(yName);
    if(color>0) h1->SetLineColor(color);
    
  }
  
  return h1;
}

TH2* ProcessHisto2D( AliCounterCollection* counter, TString hVariable, TString hVariable2, TString hSelection, TString hName){
	
  
  TH2* h1 = 0x0;
  if( !counter ) return h1;
  Bool_t setName = kTRUE;
  
  if(hName.Sizeof()==1) setName = kFALSE;
  
  h1 = (TH2*) counter->Draw(hVariable,hVariable2,hSelection);
  if ( !h1 ) h1 = new TH2D(hName,"",10,0,10,10,0,10);
  else {
    if ( h1->GetSumw2N() == 0 ) h1->Sumw2();
    h1->LabelsOption("a");
    if(setName) h1->SetName(hName);
  }
  
  return h1;
}

Bool_t GetTriggerLists(const char* triggerList, TString listFromContainer, TObjArray *triggersB, TObjArray *triggersAC, TObjArray *triggersE, TObjArray *triggersShortName){
	
  //Get the trigger list from a file
  //The file should consist of a line for each trigger with the following layout:
  //        MB triggernameB triggerNameAC triggerNameE
  //     or MUONUNLIKE triggernameB triggerNameA,triggerNameC triggerNameE
  //     or NOSHOW triggernameB triggerNameACE notrigger
  //if filename is 0, then default trigger names (pp 2011) are used
	
  if( !triggersB || !triggersAC || !triggersE || !triggersShortName) return kFALSE;
  Int_t const nColumn = 4;
  TObjArray* triggers[nColumn] = {triggersShortName, triggersB, triggersAC, triggersE};
  
  TString trigSuffix[nColumn] = {"","B", "AC", "E"};
  TString currTrigName = "";
  TObjArray* fullTriggerList[nColumn];
	
  for ( Int_t ibeam=0; ibeam<nColumn; ++ibeam ) {
    fullTriggerList[ibeam] = new TObjArray;
    fullTriggerList[ibeam]->SetOwner();
  }
  
  // Build trigger list (from file or default)
  if ( triggerList ) {
    // only the ones in the triggerList
    ifstream inFile(triggerList);
    if (!inFile.is_open()) {
      Error("PlotMuonQA","unable to open file %s", triggerList);
      return kFALSE;
    }
    
    while ( !inFile.eof() ) {
      Bool_t isGoodB = kTRUE;
      for ( Int_t ibeam=0; ibeam<nColumn; ++ibeam ) {
        currTrigName.ReadToken( inFile );
        if ( ! isGoodB ) continue;
        if ( currTrigName.IsNull() || ! currTrigName.IsAscii() ) {
          if ( ibeam==0 || ibeam == 1 ) {
            isGoodB = kFALSE;
            continue;
          }
          currTrigName = "notrigger";
        }
        fullTriggerList[ibeam]->AddLast(new TObjString(currTrigName));
      }
    }
    inFile.close();
  }
  else {
    TString baseTrigName[4] = {"CINT7", "CMUSH7", "CMUL7", "CMUU7"};
    for ( Int_t ibase=0; ibase<4; ++ibase ) {
      for ( Int_t ibeam=0; ibeam<nColumn; ++ibeam ) {
        // by default all triggers from new period in LHC11c
        currTrigName = baseTrigName[ibase] + trigSuffix[ibeam];
        fullTriggerList[ibeam]->AddLast(new TObjString(currTrigName));
      }
    }
  }
  //
  // Select only existing triggers in container
  //
  TObjArray *triggersFromContainer = listFromContainer.Tokenize(",");
  TObjString* trigName = 0x0;
	
  TString selectAllTriggers[nColumn] = {"", "", "", ""};
  for ( Int_t itrig=0; itrig<fullTriggerList[0]->GetEntries(); ++itrig ) {
    Bool_t isBadTrig = kFALSE;
    for ( Int_t ibeam=0; ibeam<nColumn; ++ibeam ) {
      currTrigName = fullTriggerList[ibeam]->At(itrig)->GetName();
      //condition on trigger name from trigger list
      if ( ibeam == 0 && currTrigName.Contains("NOSHOW") ) {
	break;
      }
      if ( currTrigName.Contains("notrigger") ){
        isBadTrig = kTRUE;
        if ( ibeam == 0 || ibeam == 1 ){
	  if ( ibeam == 1) if ( triggers[0]->GetLast() ) triggers[0]->RemoveAt(triggers[0]->GetLast());
	  break;
	}
        currTrigName = " ";
      }
      //select only the existing triggers in the container 
      //note that the trigger in the trigger file can be a list of different trigger
      if ( ibeam > 0 && triggersFromContainer ) {
	TIter nextTrigger( triggersFromContainer );
	isBadTrig = kTRUE;
	while ( ( trigName = static_cast<TObjString*>(nextTrigger()) ) ) {
	  if ( currTrigName.Contains(trigName->GetString()) ){
	    isBadTrig = kFALSE;
	  }
	}
	if ( isBadTrig == kTRUE ){ 
	  if ( ibeam == 1){
	    if ( triggers[0]->GetLast() != (triggers[0]->LowerBound()-1) ) triggers[0]->RemoveAt(triggers[0]->GetLast());
	    break;
	  }
	  currTrigName = " ";
	}
      }
      triggers[ibeam]->AddLast(new TObjString(currTrigName));
      if ( isBadTrig ) continue;
      if ( ! selectAllTriggers[ibeam].IsNull() ) selectAllTriggers[ibeam] += ",";
      selectAllTriggers[ibeam] += currTrigName;
    }
  }
  if(triggersFromContainer) delete triggersFromContainer;
  if(trigName) delete trigName;
	
	
  // Complete trigger list and print values
  cout<<" Nr of triggers read "<<triggers[0]->GetEntriesFast()<<endl;
  for ( Int_t ibeam=0; ibeam<nColumn; ++ibeam ) {
    triggers[ibeam]->AddLast(new TObjString(selectAllTriggers[ibeam]));
    printf(" %s triggers:\n", trigSuffix[ibeam].Data());
    triggers[ibeam]->Print();
    delete fullTriggerList[ibeam];
  }
	
  return kTRUE;
}

TString GetRunList(const char *runList, TObjArray *runs, TObjArray *runs2){

  // list of runs to be analyzed
  TString selectRuns = "run:";
  
  if(runList) {
    // only the ones in the runList
    ifstream inFile(runList);
    if (!inFile.is_open()) {
      Error("PlotMuonQA","unable to open file %s", runList);
      return selectRuns;
    }
    
    TString currRun;
    while (!inFile.eof()) {
      currRun.ReadLine(inFile, kTRUE);
      if (currRun.IsNull()) continue;
      if (!currRun.IsDigit()) {
	Error("PlotMuonQA","invalid run number: %s", currRun.Data());
	return selectRuns;
      }
      if(runs) runs->AddLast(new TObjString(Form("%09d", currRun.Atoi())));
      if(runs2) runs2->AddLast(new TObjString(Form("%d", currRun.Atoi())));
      selectRuns += Form("%s,",currRun.Data());
    }
    selectRuns.Remove(TString::kTrailing, ',');
    inFile.close();
    
  } else {
    // all runs
    cout<<"runList is not set"<<endl;
    if(runs) runs->AddLast(new TObjString("*"));
    if(runs2) runs2->AddLast(new TObjString("*"));
  }
  
  printf("selected runs from runlist %s: %s\n",runList, selectRuns.Data());
	
  return selectRuns;
}
