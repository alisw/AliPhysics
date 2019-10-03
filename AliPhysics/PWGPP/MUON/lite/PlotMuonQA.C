//--------------------------------------------------------------------------
// Macro for QA monitoring.
//
// In case it is not run with full aliroot, it needs the following libraries to compile:
//  - libSTEERBase
//  - libESD
//  - libAOD
//  - libANALYSIS
//  - libANALYSISalice
//  - libCORRFW
//  - libPWGmuon
//
// TString includePath = "-I${ALICE_ROOT}/PWGmuon ";  gSystem->SetIncludePath(includePath.Data());
//
// The macro reads results of the QA task and produce monitoring plots.
//
// Authors: Cynthia Hadjidakis - IPN Orsay
// QA histos from counters (event, track) and run per run histos
// To be done:
// - reorganize last part (reading and extracting info from run per run histos)
// - remove trigger selection when muon QA task modified (now a selection is done one triggers' name)
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
#include "TFileMerger.h"
// ALIROOT includes
#include "AliCounterCollection.h"
#endif

Int_t GetRunNumber(TString);
TString GetRunList(const char *runList, TObjArray *runs);
Bool_t MergeOutputs(const char*,const char*);
Bool_t GetTriggerLists(const char *triggerList, TString listFromContainer, TObjArray *triggersB=0, TObjArray *triggersShortName=0);
void SetCanvas(TCanvas *canvas, Int_t logy=1);
TH1* ProcessHisto( AliCounterCollection* counter, TString variable, TString selection, TString hName="", TString xName="", TString yName="", Int_t color=1);
TH2* ProcessHisto2D( AliCounterCollection* counter, TString hVariable, TString hVariable2, TString hSelection, TString hName);
Int_t GetIndex(TObjArray *triggersB, Int_t trigNr, Int_t centNr);
TCanvas *ProcessCanvasAllTrigger(AliCounterCollection *counter, TString canvasName);
TCanvas *ProcessCanvasTriggerContent(TObjArray *trigName, TObjArray trigNoPS, TObjArray trigWithPS, TString canvasName);
TCanvas *ProcessCanvasRelativeTriggerContent(TObjArray *array, TObjArray trigNoPS, TString canvasName);
TCanvas *ProcessCanvasPhysSelCut(TObjArray *triggersB, TObjArray trigNoPS, TObjArray trigWithPS, TString canvasName);
TCanvas *ProcessCanvasPhysSelCutCentrality(TObjArray *triggersB, TObjArray trigNoPS, TObjArray trigBWithPS, Int_t trigNr, TString canvasName, TString *legendHeader);
TCanvas *ProcessCanvasCentralityPercentile(TObjArray *triggersB, TObjArray trigNoPS, TObjArray trigWithPS, Int_t trigNr, TString canvasName, TString *legendHeader);
TCanvas *ProcessCanvasTracksoverTrigger(TObjArray *triggersB, TObjArray trigSel, TObjArray trackTracker, TObjArray trackTrigger, TObjArray trackMatched, TObjArray trackAll, Int_t trigNr, Int_t centNr, TString canvasName,TString legendHeader=""); 
TCanvas *ProcessCanvasTrackMultB(TObjArray *triggersB, TObjArray trigSel, TObjArray trackTracker, TObjArray trackTrigger, TObjArray trackMatched, Int_t trigNr, Int_t centNr, TString canvasName,TString legendHeader="");
TCanvas *ProcessCanvasRatioTrackB(TObjArray *triggersB, TObjArray trigSel, TObjArray trackTracker, TObjArray trackTrigger, TObjArray trackMatched, Int_t trigNr, Int_t centNr, TString canvasName,TString legendHeader="");
TCanvas *ProcessCanvasAsymMatched(TObjArray *triggersB, TObjArray trackPosMatched, TObjArray trackNegMatched, TObjArray trackhAllMatched, Int_t trigNr, Int_t centNr, TString canvasName,TString legendHeader="");
TCanvas *ProcessCanvasHighPtMuons(TObjArray *triggersB, TObjArray trigSel, TObjArray trackMatchedLowPt, TObjArray trackAllMatchedHightPt, Int_t trigNr, Int_t centNr, TString canvasName,TString legendHeader="");
TCanvas *ProcessCanvasBeamGasMatched(TObjArray *triggersB, TObjArray trackBeamGasMatched, TObjArray trackBeamGasMatchedHighPt, TObjArray trackAllMatched, TObjArray trackMatchedHighPt, Int_t trigNr, Int_t centNr, TString canvasName,TString legendHeader="");
Bool_t IsHeavyIonCollision(AliCounterCollection *eventCounters);

//loop on centrality
const Int_t kCentBinMax = 4;
TString kCentBinName[kCentBinMax] = {"v0mult:any/","v0mult:low,int,high/","v0mult:low/","v0mult:high/"};
TString kCentLegendName[kCentBinMax] ={"all collisions","[0-80%] from V0 amplitude","low mult. [60-80%] from V0 amplitude","high mult. [0-10%] from V0 amplitude"};
TString kCentLegendNameShort[kCentBinMax] ={"all","[0-80%]","[60-80%]","[0-10%]"};

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
  Int_t kCMUS = 1;

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
  TString selectRuns = GetRunList(runList,runs);
	
  // physics selection
  TString select = selectPhysics ? "selected:yes" : "";
  	
  cout<<"//---------------------------------- //"<<endl;
  cout<<"//        Get global counter        //"<<endl;
  cout<<"//---------------------------------- //"<<endl;
  
  TString mergedFilename = Form("%s/%s", baseDir,QAFileName);
  if (runList) MergeOutputs(runList, mergedFilename);

  TFile *globalFile = TFile::Open(mergedFilename.Data());
  if (!globalFile || ! globalFile->IsOpen()) {
    Error("PlotQA", "failed to open file: %s", mergedFilename.Data());
    return;
  }
  globalFile->Cd("MUON_QA");
  
  TString selection;
	
  // get counters
  AliCounterCollection* eventCounters = static_cast<AliCounterCollection*>(globalFile->FindObjectAny("eventCounters"));
  AliCounterCollection* trackCounters = static_cast<AliCounterCollection*>(globalFile->FindObjectAny("trackCounters"));

  // run list from counters
  if (!runList) selectRuns += trackCounters->GetKeyWords("run");
	
  cout<<"//---------------------------------- //"<<endl;
  cout<<"//        Trigger selection          //"<<endl;
  cout<<"//---------------------------------- //"<<endl;

  TObjArray *triggersB, *triggersShortName;
  triggersB = new TObjArray();
  triggersB->SetOwner();
  triggersShortName = new TObjArray();
  triggersShortName->SetOwner();

  TString listFromContainer = eventCounters->GetKeyWords("trigger");
  Bool_t success = GetTriggerLists(triggerList, listFromContainer, triggersB, triggersShortName);
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
  TString CanvasName = "AllTriggers";
  TCanvas *cAll = ProcessCanvasAllTrigger(eventCounters, CanvasName);

  cout<<"//---------------------------------- //"<<endl;
  cout<<"//   Processing event counters       //"<<endl;
  cout<<"//---------------------------------- //"<<endl;

  //declare a default canvas c1 
  CanvasName = "c1";
  TCanvas *c1 = new TCanvas(CanvasName.Data(),CanvasName.Data());
  c1->cd();
  
  Int_t nCentBin = ( isHeavyIon ) ? kCentBinMax : 1;
  TString selectionCent;


  TArrayI colorInd( triggersB->GetEntriesFast() );
  Int_t const colorNrFirst = 8;
  Int_t colorIndex = 0;
  Int_t colorTabFirst[colorNrFirst] = {kGray+2,kRed,kBlue,kGreen,kOrange,kCyan,kMagenta,kYellow};
  for ( Int_t i = 0; i < triggersB->GetEntriesFast(); i++ ) {
    colorInd[i] = colorTabFirst[i%colorNrFirst] + colorIndex;
    if ( i%colorNrFirst == 0 ) colorIndex++;
  } 


  Int_t nTrig = triggersB->GetEntriesFast();

  TObjArray trigNoPS(nTrig*nCentBin);
  TObjArray trigWithPS(nTrig*nCentBin);
  TObjArray trigWithPST0Flag(nTrig*nCentBin);
  TObjArray trigWithPST0SPDFlag(nTrig*nCentBin);

  TObjArray trigSel;
  TObjArray trackTrigger(nTrig*nCentBin); 
  TObjArray trackTracker(nTrig*nCentBin); 
  TObjArray trackMatched(nTrig*nCentBin);
  TObjArray trackAll(nTrig*nCentBin); 

  TObjArray trackAllMatched(nTrig*nCentBin); //matched tracks + additional geometrical cut
  TObjArray trackMatchedLowPt(nTrig*nCentBin);
  TObjArray trackMatchedHighPt(nTrig*nCentBin);

  TObjArray trackPosMatched(nTrig*nCentBin);
  TObjArray trackNegMatched(nTrig*nCentBin);

  TObjArray trackBeamGasMatched(nTrig*nCentBin);
  TObjArray trackBeamGasMatchedHighPt(nTrig*nCentBin);

  cout<<Form("Processing for %d triggers...",triggersB->GetEntriesFast())<<endl;
    
  //loop on centrality
  for ( Int_t iCentBin = 0; iCentBin < nCentBin; iCentBin++){
    selectionCent = kCentBinName[iCentBin];
		
    //Loop on trigger (last is all triggers)
    for(Int_t iTrig = 0; iTrig < triggersB->GetEntriesFast(); iTrig++){
    
      TString histoNameBase = "h_trig", histoName;
      histoNameBase+= iTrig+1;
 		
      Int_t index = GetIndex(triggersB,iTrig,iCentBin);

      TH1* histo = 0;

      TString triggerName = ( (TObjString*) triggersB->At(iTrig) )->GetString();
      if(triggerName.EqualTo(" ")) continue;
      // Histo trigger without Phys. Sel. 
      selection = selectionCent; selection += Form("trigger:%s/%s", triggerName.Data(), selectRuns.Data());		
      // cout<<selection<<endl;
      histoName = histoNameBase;
      histoName += "BNoPS";
      histo = (TH1*) ProcessHisto(eventCounters, "run", selection, histoName, "", "Trigger content w/o Phys. Sel.", colorInd[iTrig]);
      trigNoPS.AddAt(histo,index);
      // Histo trigger with Phys. Sel. 
      selection = selectionCent; selection += Form("trigger:%s/%s/selected:yes", triggerName.Data(), selectRuns.Data());
      histoName = histoNameBase;
      histoName += "BWithPS";
      histo = (TH1*) ProcessHisto(eventCounters, "run", selection, histoName, "", "Trigger content w/ Phys. Sel.", colorInd[iTrig]);
      trigWithPS.AddAt(histo,index);
    // Histo trigger with Phys. Sel. and T0 pile up not flagged
      selection = selectionCent; selection += Form("trigger:%s/%s/selected:yes/t0pileup:no", triggerName.Data(), selectRuns.Data());
      histoName = histoNameBase;
      histoName += "BWithPST0Flag";
      histo = (TH1*) ProcessHisto(eventCounters, "run", selection, histoName, "", "Trigger content w/ Phys. Sel. and no pile up from T0 flag", colorInd[iTrig]);
      trigWithPST0Flag.AddAt(histo,index);
    // Histo trigger with Phys. Sel. and T0 + SPD pile up not flagged
      selection = selectionCent; selection += Form("trigger:%s/%s/selected:yes/t0pileup:no/spdpileup:no", triggerName.Data(), selectRuns.Data());
      histoName = histoNameBase;
      histoName += "BWithPST0SPDFlag";
      histo = (TH1*) ProcessHisto(eventCounters, "run", selection, histoName, "", "Trigger content w/ Phys. Sel. and no pile up from T0 and SPD flag", colorInd[iTrig]);
      trigWithPST0SPDFlag.AddAt(histo,index);
		
      // Histo tracking : Phys. Sel.  is selected or not depending on the macro arguments
      selection = selectionCent; selection += Form("track:triggeronly/trigger:%s/%s/%s", triggerName.Data(), selectRuns.Data(), select.Data());
      histo = (TH1*) ProcessHisto(trackCounters, "run", selection, "");
      trackTrigger.AddAt(histo,index);

      selection = selectionCent; selection += Form("track:trackeronly/trigger:%s/%s/%s", triggerName.Data(), selectRuns.Data(), select.Data());
      histo = (TH1*) ProcessHisto(trackCounters, "run", selection, "");
      trackTracker.AddAt(histo,index);

      selection = selectionCent; selection += Form("track:matched/trigger:%s/%s/%s", triggerName.Data(), selectRuns.Data(), select.Data());
      histo = (TH1*) ProcessHisto(trackCounters, "run", selection, "");
      trackMatched.AddAt(histo,index);
		
      histo = (TH1*) ( (TH1*) trackTrigger.At(index))->Clone("");
      histo->Add( (TH1*) trackTracker.At(index) );
      histo->Add( (TH1*) trackMatched.At(index) );
      trackAll.AddAt(histo,index);

      //for the following, only integrated over centrality
      if ( iCentBin > 0 ) continue;
      
      selection = selectionCent; selection += Form("track:matched/trigger:%s/%s/%s/pt:low/acc:in", triggerName.Data() ,selectRuns.Data(), select.Data());
      histo = (TH1*) ProcessHisto(trackCounters, "run", selection, "");
      trackMatchedLowPt.AddAt(histo,index);

      selection = selectionCent; selection += Form("track:matched/trigger:%s/%s/%s/pt:high/acc:in", triggerName.Data() ,selectRuns.Data(), select.Data());
      histo = (TH1*) ProcessHisto(trackCounters, "run", selection, "");
      trackMatchedHighPt.AddAt(histo,index);

      selection = selectionCent; selection += Form("track:matched/trigger:%s/%s/charge:pos/%s/acc:in",triggerName.Data(), select.Data(),selectRuns.Data());
      histo = (TH1*) ProcessHisto(trackCounters, "run", selection, "");
      trackPosMatched.AddAt(histo,index);
    
      selection = selectionCent; selection += Form("track:matched/trigger:%s/%s/charge:neg/%s/acc:in",triggerName.Data(), select.Data(),selectRuns.Data());
      histo =  (TH1*) ProcessHisto(trackCounters, "run", selection, "");
      trackNegMatched.AddAt(histo,index);

      selection = selectionCent; selection += Form("track:matched/trigger:%s/%s/%s/acc:in",triggerName.Data(), select.Data(),selectRuns.Data());
      histo =  (TH1*) ProcessHisto(trackCounters, "run", selection, "");
      trackAllMatched.AddAt(histo,index);

      selection = selectionCent; selection += Form("track:matched/trigger:%s/%s/%s/acc:in/tagTrack:beamGas",triggerName.Data(), select.Data(),selectRuns.Data());
      histo =  (TH1*) ProcessHisto(trackCounters, "run", selection, "");
      trackBeamGasMatched.AddAt(histo,index);

      selection = selectionCent; selection += Form("track:matched/trigger:%s/%s/%s/acc:in/tagTrack:beamGas/pt:high",triggerName.Data(), select.Data(),selectRuns.Data());
      histo =  (TH1*) ProcessHisto(trackCounters, "run", selection, "");
      trackBeamGasMatchedHighPt.AddAt(histo,index);
    }
  }

  if ( selectPhysics) trigSel = trigWithPS;
  else trigSel = trigNoPS;
	
  //if there is not Beam triggers just stop now
  Int_t count_trigger=0;
  for(Int_t iTrig = 0; iTrig < triggersB->GetEntriesFast(); iTrig++){
    Int_t centBinNr = 0;
    Int_t index = GetIndex(triggersB,iTrig,centBinNr);
    count_trigger += (static_cast<TH1*>(trigNoPS.At(index)))->GetEntries();
  }
  if(count_trigger<=0) return;
	
  Double_t NumOfBNoPS[nCentBin*nTrig];
  Double_t NumOfBWithPS[nCentBin*nTrig];
	
  for ( Int_t iCentBin = 0; iCentBin < nCentBin; iCentBin++){
    for ( Int_t iTrig = 0; iTrig < nTrig; iTrig++){
      Int_t index = GetIndex(triggersB,iTrig,iCentBin);
      NumOfBNoPS[index] = (static_cast<TH1*>(trigNoPS.At(index)))->Integral();
      NumOfBWithPS[index] = (static_cast<TH1*>(trigWithPS.At(index)))->Integral();
    }
  }
	
  cout<<"//==================================================================================="<<endl;
  cout<<"// Put all plots in a ps file, easier to publish (Twiki)"<<endl;
  cout<<"//==================================================================================="<<endl;
  
  c1->Print(OutFileNamePDF_open.Data());
  TFile *rootFileOut = TFile::Open(OutFileNameROOT.Data(),"RECREATE");
  rootFileOut->cd();
  TDirectoryFile *dirGlobal = new TDirectoryFile("Global","Global","",(TDirectory*)rootFileOut->GetMotherDir());
  cout<<"dirGlobal mother "<<(dirGlobal->GetMotherDir())->GetName()<<endl;
  //An array of TDirectoryFile
  TObjArray *dirTrigger = new TObjArray();
  dirTrigger->SetOwner();
  TObjArray *dirCent = new TObjArray();
  dirCent->SetOwner();
  for ( Int_t i = 0; i < triggersB->GetEntriesFast()-1 ; i++) {
    TString currTrigName = ( (TObjString*) triggersShortName->At(i) )->GetString();
    TDirectoryFile *dirFile = new TDirectoryFile( currTrigName.Data(),currTrigName.Data(),"",(TDirectory*)rootFileOut->GetMotherDir() );
    dirTrigger->AddLast( dirFile );
    for( Int_t j = 0; j < nCentBin; j++) {
      TString centName = kCentLegendNameShort[j];
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
	
  TCanvas *cTriggerContent = ProcessCanvasTriggerContent(triggersB, trigNoPS, trigWithPS, "TriggerContent");
  cTriggerContent->Print(OutFileNamePDF.Data());
  dirGlobal->cd();
  cTriggerContent->Write();
  cTriggerContent->Close();

  cout<<"//==================================================================================="<<endl;
  cout<<"// new canvas with the relative content of each trigger w/o physics selection"<<endl;
  cout<<"//==================================================================================="<<endl;
	
  TCanvas *cRelativeTriggerContent = ProcessCanvasRelativeTriggerContent(triggersB, trigNoPS, "RelativeTriggerContent");
  cRelativeTriggerContent->Print(OutFileNamePDF.Data());
  dirGlobal->cd();
  cRelativeTriggerContent->Write();
  cRelativeTriggerContent->Close();
  

  cout<<"//==================================================================================="<<endl;
  cout<<"// new canvas with effect from physics selection for each trigger "<<endl;
  cout<<"//==================================================================================="<<endl;
	
  TCanvas *cPhysSelCut = 0x0;
  cPhysSelCut = ProcessCanvasPhysSelCut(triggersB, trigNoPS, trigWithPS, "PhysSelCutOnCollTrigger");
  cPhysSelCut->Print(OutFileNamePDF.Data());
  dirGlobal->cd();
  cPhysSelCut->Write();
  cPhysSelCut->Close();
  rootFileOut->cd();
	
  cout<<"//==================================================================================="<<endl;
  cout<<"// new canvas with effect from T0 pile-up flag for each trigger "<<endl;
  cout<<"//==================================================================================="<<endl;
	
  TCanvas *cPhysSelCutT0Flag = 0x0;
  cPhysSelCutT0Flag = ProcessCanvasPhysSelCut(triggersB, trigNoPS, trigWithPST0Flag, "PhysSelCutOnCollTriggerT0Flag");
  cPhysSelCutT0Flag->Print(OutFileNamePDF.Data());
  dirGlobal->cd();
  cPhysSelCutT0Flag->Write();
  cPhysSelCutT0Flag->Close();
  rootFileOut->cd();
	

 cout<<"//==================================================================================="<<endl;
  cout<<"// new canvas with effect from T0 + SPD pile-up flag for each trigger "<<endl;
  cout<<"//==================================================================================="<<endl;
	
  TCanvas *cPhysSelCutT0SPDFlag = 0x0;
  cPhysSelCutT0SPDFlag = ProcessCanvasPhysSelCut(triggersB, trigNoPS, trigWithPST0SPDFlag, "PhysSelCutOnCollTriggerT0SPDFlag");
  cPhysSelCutT0SPDFlag->Print(OutFileNamePDF.Data());
  dirGlobal->cd();
  cPhysSelCutT0SPDFlag->Write();
  cPhysSelCutT0SPDFlag->Close();
  rootFileOut->cd();

  TString canvasName;
  TString triggerName;	
  TString legendHeader="";
	
  if ( isHeavyIon ){

  cout<<"//==================================================================================="<<endl;
  cout<<"// new canvas with effect from physics selection for each trigger and centrality bin (only in heavy-ion) "<<endl;
  cout<<"//==================================================================================="<<endl;

    TCanvas *cPhysSelCutCentrality;
	
    //loop over trigger
    for( Int_t iTrig = 0; iTrig < triggersB->GetEntriesFast(); iTrig++){
      //skip sum of all triggers
      if(iTrig == (triggersB->GetEntriesFast()-1)) continue;

      canvasName = "PhysSel_trigger";
      canvasName += ( (TObjString*) triggersShortName->At(iTrig) )->GetString();
	
      cPhysSelCutCentrality = ProcessCanvasPhysSelCutCentrality(triggersB, trigNoPS, trigWithPS, iTrig, canvasName, kCentLegendNameShort+1);
      cPhysSelCutCentrality->Print(OutFileNamePDF.Data());
      ( (TDirectoryFile*) dirTrigger->At(iTrig) )->cd();
      cPhysSelCutCentrality->Write();		
      cPhysSelCutCentrality->Close();		
    }
  }
  rootFileOut->cd();
	
  if ( isHeavyIon ){
    cout<<"//==================================================================================="<<endl;
    cout<<"// new canvas for centrality percentile check (only in PbPb) "<<endl;
    cout<<"//==================================================================================="<<endl;
   
    TCanvas *cCentralityCheck;
		
    //loop over trigger
    for ( Int_t iTrig = 0; iTrig < triggersB->GetEntriesFast(); iTrig++){
      //skip sum of all triggers
      if ( iTrig == (triggersB->GetEntriesFast()-1) ) continue;
 			
      canvasName = "CentralityCheck_trigger";
      canvasName += ( (TObjString*) triggersShortName->At(iTrig) )->GetString();
			
      cCentralityCheck = ProcessCanvasCentralityPercentile(triggersB,trigNoPS,trigWithPS,iTrig,canvasName,kCentLegendNameShort);

      cCentralityCheck->Print(OutFileNamePDF.Data());
      ( (TDirectoryFile*) dirTrigger->At(iTrig) )->cd();
      cCentralityCheck->Write();
      cCentralityCheck->Close();
    }
  }
  rootFileOut->cd();
	

  cout<<"//==================================================================================="<<endl;
  cout<<"// Ratio of tracks over trigger type (2 canvases) "<<endl;
  cout<<"//==================================================================================="<<endl;

  //Print a canvas per trigger type
  TCanvas *cTracksoverTrigger;
  TCanvas* cTrackMultB;
 	
  //loop over centrality bins
  for ( Int_t iCentBin = 0; iCentBin < nCentBin; iCentBin++){
    if ( isHeavyIon ){
      legendHeader = "for ";
      legendHeader += kCentLegendName[iCentBin];
    }
    else legendHeader ="";
    //loop over triggers
    for(Int_t iTrig = 0; iTrig < triggersB->GetEntriesFast(); iTrig++){
      //skip sum of all triggers
      if( iTrig == (triggersB->GetEntriesFast()-1) ) continue;

      ( (TDirectoryFile*) dirCent->At( iTrig*nCentBin+iCentBin ) )->cd();

      canvasName = "RatioTrackTypes_cent";
      canvasName += iCentBin;
      canvasName +="trigger";
      canvasName += ( (TObjString*) triggersShortName->At(iTrig) )->GetString();
      cTracksoverTrigger = ProcessCanvasTracksoverTrigger(triggersB, trigSel, trackTracker, trackTrigger, trackMatched, trackAll, iTrig, iCentBin, canvasName,legendHeader);
      cTracksoverTrigger->Print(OutFileNamePDF.Data());
      cTracksoverTrigger->Write();
      cTracksoverTrigger->Close();

      canvasName = "TrackMult_cent";
      canvasName += iCentBin;
      canvasName +="trigger";
      canvasName +=( (TObjString*) triggersShortName->At(iTrig) )->GetString();		
      cTrackMultB= ProcessCanvasTrackMultB(triggersB, trigSel, trackTracker, trackTrigger, trackMatched, iTrig, iCentBin, canvasName, legendHeader);
      cTrackMultB->Print(OutFileNamePDF.Data());
      cTrackMultB->Write();
      cTrackMultB->Close();
    }
  }
  rootFileOut->cd();

  cout<<"//===================================================="<<endl;
  cout<<"// Draw matched tracks asymmetry for each trigger "<<endl;
  cout<<"//===================================================="<<endl;
	
  //Print a canvas per trigger type
  TCanvas *cAsymMatched;

  legendHeader = "for all collisions";
  for ( Int_t iTrig = 0; iTrig < triggersB->GetEntriesFast(); iTrig++){
    //skip sum of all triggers
    if ( iTrig == (triggersB->GetEntriesFast()-1) ) continue;
    
    ( (TDirectoryFile*) dirTrigger->At( iTrig ) )->cd();
    
    canvasName = "AsymMatched";
    canvasName +="trigger";
    canvasName +=( (TObjString*) triggersShortName->At(iTrig) )->GetString();
    Int_t iCentBin = 0;
    cAsymMatched = ProcessCanvasAsymMatched(triggersB, trackPosMatched, trackNegMatched, trackAllMatched, iTrig, iCentBin, canvasName,legendHeader);
    cAsymMatched->Print(OutFileNamePDF.Data());
    cAsymMatched->Write();
    cAsymMatched->Close();
  }
  rootFileOut->cd();
  legendHeader = ""; 

  cout<<"//===================================================================="<<endl;
  cout<<"// Draw beam gas contribution to matched tracks for each trigger "<<endl;
  cout<<"//===================================================================="<<endl;
	
  //Print a canvas per trigger type
  TCanvas *cBeamGasMatched;

  //loop over trigger
  for ( Int_t iTrig = 0; iTrig < triggersB->GetEntriesFast(); iTrig++){
    //skip sum of all triggers
    if ( iTrig == (triggersB->GetEntriesFast()-1)) continue;

    ( (TDirectoryFile*) dirTrigger->At(iTrig) )->cd();
    canvasName = "BeamGasMatched";
    canvasName +="trigger";
    canvasName +=( (TObjString*) triggersShortName->At(iTrig) )->GetString();
    Int_t iCentBin = 0;
    cBeamGasMatched= ProcessCanvasBeamGasMatched(triggersB, trackBeamGasMatched, trackBeamGasMatchedHighPt, trackAllMatched, trackMatchedHighPt, iTrig, iCentBin, canvasName,legendHeader);
    cBeamGasMatched->Print(OutFileNamePDF.Data());
    cBeamGasMatched->Write();
    cBeamGasMatched->Close();
  }
  rootFileOut->cd();

  cout<<"//=================================================="<<endl;
  cout<<"// Draw high pt tracks per trigger"<<endl;
  cout<<"//=================================================="<<endl;

  //Print a canvas per trigger type
  TCanvas *cHighPtMuons;
  legendHeader = "for all collisions";
  //loop over trigger
  for ( Int_t iTrig = 0; iTrig < triggersB->GetEntriesFast(); iTrig++){
    //skip sum of all triggers
    if( iTrig == (triggersB->GetEntriesFast()-1)) continue;
    
    ( (TDirectoryFile*) dirTrigger->At( iTrig ) )->cd();
	
    canvasName = "HighPtMuons";
    canvasName +="trigger";
    canvasName +=( (TObjString*) triggersShortName->At(iTrig) )->GetString();

    Int_t iCentBin = 0;
    cHighPtMuons = ProcessCanvasHighPtMuons(triggersB, trigSel, trackMatchedLowPt, trackMatchedHighPt, iTrig, iCentBin, canvasName,legendHeader);
    cHighPtMuons->Print(OutFileNamePDF.Data());
    cHighPtMuons->Write();
    cHighPtMuons->Close();
  }
  rootFileOut->cd();
	
  // close merged file	
  globalFile->Close();
  
  //===================================================================================
  //Print out the number of trigger without and with Phys. Sel.
  //===================================================================================
  
  cout << endl << endl;
  //====================================================
  if (PRINTSTAT){
    if ( triggersB->At(kCMUS) ) { 
    
      // set the format to print labels
      Int_t centBinNr = 0;
      THashList* labels = (static_cast<TH1*>(trigWithPS.At(triggersB->GetEntriesFast()*centBinNr+kCMUS)))->GetXaxis()->GetLabels();
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
	printf(format.Data(), label->String().Data(), (Int_t) (static_cast<TH1*>(trigWithPS.At(triggersB->GetEntriesFast()*centBinNr+kCMUS)))->GetBinContent(bin));
      }
      printf("\n========== Total #runs = %d ==============\n",nRuns);
      printf("\n\n");
    
    
      cout << "-------------------------------------------------" << endl;
      cout << "Total statistic" << endl; 
      cout << " " << endl ; 
    
      cout << "-------------------------------------------------------------------" << endl;
      cout<<"Number of triggers w/o Phys. Sel./ w/ Phys. Sel (Phys. Sel. cut in %)"<<endl;
      for(Int_t iTrig = 0; iTrig < triggersB->GetEntriesFast()-1; iTrig++){
	TString triggerNameB = ( (TObjString*) triggersB->At(iTrig) )->GetString();
			
	Double_t	cutinpercent	=	0;
	printf("%10s\t",triggerNameB.Data());
	Int_t index = GetIndex(triggersB,iTrig,centBinNr);
	if(NumOfBNoPS[index] > 0.) cutinpercent = (NumOfBNoPS[index]-NumOfBWithPS[index])/NumOfBNoPS[index]*100.;
	printf("%g / %g (%g%%)\n", NumOfBNoPS[index],NumOfBWithPS[index],cutinpercent);
      }
    } 
  }
  
  //temporary
  //    return;

  
  //--------------------------------------------- //
  //        monitor quantities run per run        //
  //--------------------------------------------- //

  //output histos
  //const Int_t kNHistos =  100;
  //TObjArray outHist(kNHistos);

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

    // get the file (or list of files) to be analyzed
    TString command;
    TGridResult *res = 0;
    
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
      
      command = Form("find %s/*%s/ -name %s | xargs", alienBaseDir.Data(), run.Data(), QAFileName);
      TString foundFiles = gSystem->GetFromPipe(command.Data());
      TObjArray* arr = foundFiles.Tokenize(" ");

      for ( Int_t iarr=0; iarr<arr->GetEntries(); iarr++ ) {
        res->Add(new TObjString(arr->At(iarr)->GetName()));
      }
      delete arr;
    }
    
    // Loop over 'find' results and get next LFN
    TIter nextmap(res);
    TObjString *objs = 0;
    TObject* currObj = 0x0;
    Bool_t searchRunNr = kFALSE;
    if ( !isAlienFile && run.Contains("*") ) searchRunNr = kTRUE;

    //some checks
    while ( ( currObj = nextmap() ) ){
      
      // get the current file url
      if(isAlienFile){
        objs = static_cast<TObjString*>(static_cast<TMap*>(currObj)->GetValue("turl"));
      }
      else{
        objs=static_cast<TObjString*>(currObj);
      }
      
      if (objs->GetString().IsNull()) {
	Error("PlotMuonQA","turl/obj not found for the run %s... SKIPPING", run.Data());
	continue;
      }
    
      if ( searchRunNr ) {
	Int_t runNr = GetRunNumber(objs->GetString());
	if (runNr > 0) run = Form("%i",runNr);
      }
  
      if ( run.IsDigit() && ! selectRuns.Contains(Form("%i",run.Atoi())) ) continue;
      
      // open the outfile for this run
      TFile *runFile = TFile::Open(objs->GetString());
      if (!runFile || ! runFile->IsOpen()) {
	Error("PlotMuonQA","failed to open file: %s", objs->GetName());
	continue;
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
  }//end loop over runs
  
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
  cout<<"=================================================="<<endl;
  cout<<"Display Mean and Sigma of the number of associated clusters to a track "<<endl;
  cout<<"=================================================="<<endl;
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
  
  cout<<"=================================================="<<endl;
  cout<<" Display average number of cluster per chamber "<<endl;
  cout<<"=================================================="<<endl;

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

  cout<<"=================================================="<<endl;
  cout<<" Display average cluster charge per chamber "<<endl;
  cout<<"=================================================="<<endl;                                                                   
       
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

  cout<<"=================================================="<<endl;
  cout<<" Display average cluster size per chamber     "<<endl;
  cout<<"=================================================="<<endl;

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


  cout<<"=================================================="<<endl;
  cout<<"Display average X and Y position of clusters per chamber"<<endl;
  cout<<"=================================================="<<endl;
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


  cout<<"=================================================="<<endl;
  cout<<" Display tracks ChiSquare "<<endl;
  cout<<"=================================================="<<endl;
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

  cout<<"=================================================="<<endl;
  cout<<" Display track Lpt/Hpt "<<endl;
  cout<<"=================================================="<<endl;

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


  cout<<"=================================================="<<endl;
  cout<<" Display muon trigger "<<endl;
  cout<<"=================================================="<<endl;
  
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
 c1->Close();
 //Note: closing the file delete all related TDirectoryFile (dirCent and dirTrigger)
 rootFileOut->Close(); 
 

 delete runs;
 delete triggersB;
  delete triggersShortName;

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

Bool_t IsHeavyIonCollision(AliCounterCollection *eventCounters){
	
  if(!eventCounters) return kFALSE;
	
  Double_t sum = eventCounters->GetSum("v0mult:low,int,high");
  Bool_t result = kTRUE;
  if(sum<=0) result = kFALSE;
	
  cout<<" Collision type is set to ";
  if( result == kFALSE) cout<<"p-p"<<endl;
  else cout<<"heavy-ion"<<endl;
	
  return result;
}

TCanvas *ProcessCanvasAllTrigger( AliCounterCollection *eventCounters, TString canvasName) {

  if ( !eventCounters ) return 0;

  TString cName = Form("c%s",canvasName.Data());
  TCanvas *cAll = new TCanvas(canvasName.Data(),canvasName.Data());
  cAll->SetLeftMargin(0.18);
  cAll->SetRightMargin(0.18);
  cAll->SetLogz(1);
  cAll->cd();

  TH2* hAll = (TH2*) ProcessHisto2D(eventCounters, "trigger", "run", "run:any" , "");
  for ( Int_t ibin=1; ibin <= hAll->GetYaxis()->GetNbins(); ++ibin ) {
    TString currLabel = hAll->GetYaxis()->GetBinLabel(ibin);
    TObjArray* labelArray = currLabel.Tokenize("-");
    labelArray->SetOwner();
    //cout<<currLabel<<endl;
    TString newLabel = labelArray->At(0)->GetName();
    if ( labelArray->GetEntries() >= 2 ) newLabel = Form("%s-%s", newLabel.Data(), labelArray->At(1)->GetName());
    hAll->GetYaxis()->SetBinLabel(ibin, newLabel.Data());
    delete labelArray;
  }
  hAll->Draw("COLZ");

  return cAll;
}

TCanvas *ProcessCanvasTriggerContent(TObjArray *array, TObjArray trigNoPS, TObjArray trigWithPS, TString canvasName){
 
  if(!array) return 0x0;
	
  TString cName =  "c"; 
  cName += canvasName; 
  TCanvas *cTriggerContent = new TCanvas(canvasName,cName,1200,900);
  SetCanvas(cTriggerContent);
  cTriggerContent->cd();
 
  TLegend* legcTC = new TLegend(0.2,0.15,0.50,0.40);
  legcTC->SetHeader("Physics Selection");
  legcTC->AddEntry(".","applied :","");
 
  for(Int_t iTrig = 0; iTrig < array->GetEntriesFast(); iTrig++){
    //skip the sum of all triggers
    if( iTrig == (array->GetEntriesFast()-1) ) continue;
    TH1* hNoPS = static_cast<TH1*>(trigNoPS.At(iTrig));
    TH1* hWithPS = static_cast<TH1*>(trigWithPS.At(iTrig));
    if (!hNoPS ||!hWithPS) continue;
    hNoPS->SetLineStyle(2);
    if(iTrig==0){
      hNoPS->SetMinimum(1e-3);
      hNoPS->Draw();
      hWithPS->Draw("same");
    }
    else{
      hNoPS->Draw("same");
      hWithPS->Draw("same");
    }
    legcTC->AddEntry(hWithPS,(( (TObjString*) array->At(iTrig) )->GetString()).Data(),"l");
  }
  legcTC->AddEntry(".","not applied :","");
	
  for(Int_t iTrig = 0; iTrig < array->GetEntriesFast(); iTrig++){
   if( iTrig == (array->GetEntriesFast()-1) ) continue;
    TH1* hNoPS = static_cast<TH1*>(trigNoPS.At(iTrig));
    if(hNoPS) legcTC->AddEntry(hNoPS,(( (TObjString*) array->At(iTrig) )->GetString()).Data(),"l");	 
  }
	
  legcTC->Draw("same");
 
  return cTriggerContent;
}

TCanvas *ProcessCanvasRelativeTriggerContent(TObjArray *triggersB, TObjArray trigNoPS, TString canvasName){
	
  if(!triggersB) return 0x0;
	
  TString cName =  "c" ; 
  cName += canvasName;
  TCanvas *cRelativeTriggerContent = new TCanvas(canvasName,cName,1200,900);
  SetCanvas(cRelativeTriggerContent);
  cRelativeTriggerContent->cd();
	
  TObjArray relTrigNoPS(triggersB->GetEntriesFast());
  TLegend* legcRTC = new TLegend(0.2,0.15,0.50,0.40);
	
  TString hName, hTriggerName;
  Int_t indAllTrig = triggersB->GetEntriesFast()-1;
  TH1* hAllTrig = static_cast<TH1*> (trigNoPS.At(indAllTrig));
  if(!hAllTrig) return 0;

  for(Int_t iTrig = 0; iTrig < triggersB->GetEntriesFast()-1; iTrig++){
    hName = "ratio";
    hName += ( (TObjString*) triggersB->At(iTrig) )->GetString();
    TH1* histo = static_cast<TH1*> (trigNoPS.At(iTrig));
    if(!histo) continue;
    TH1* hRatio = (TH1*) histo->Clone(hName);
    hRatio->Divide(hAllTrig);
    hRatio->SetLineStyle(1);
    if(iTrig==0){
      hRatio->SetMaximum(1.5);
      hRatio->SetMinimum(0.001);
      hRatio->SetLabelSize(0.04);
      hRatio->GetYaxis()->SetTitle("Relative trigger content"); 
      hRatio->Draw("E");
    }
    else{
      hRatio->Draw("ESAME");
    }
    relTrigNoPS.AddAt(hRatio,iTrig);
  }

  legcRTC->SetHeader("Physics Selection not applied");	
  for(Int_t iTrig = 0; iTrig < triggersB->GetEntriesFast()-1; iTrig++){
    TH1* hRatio = static_cast<TH1*>(relTrigNoPS.At(iTrig));
    if (!hRatio) continue;
    legcRTC->AddEntry(hRatio,(( (TObjString*) triggersB->At(iTrig) )->GetString()).Data(),"l");
  }
  legcRTC->Draw("same");
	
  return cRelativeTriggerContent;
}

TCanvas *ProcessCanvasPhysSelCut(TObjArray *triggersB, TObjArray trigNoPS, TObjArray trigWithPS, TString canvasName){
	
  if(!triggersB) return 0x0;

  TString cName = "c";
  cName += canvasName;
  TCanvas *c1 = new TCanvas(canvasName,cName,1200,900);
  SetCanvas(c1);
  c1->cd();
	 
  TObjArray trigRatio(triggersB->GetEntriesFast());
  TLegend* legcRTC = new TLegend(0.2,0.15,0.50,0.40);
  TString header = "Physics Selection Cut on selected triggers:";
  if (canvasName.Contains("T0Flag")) header += " and T0 pile-up event selection"; 
  if (canvasName.Contains("T0SPDFlag")) header += " and T0, SPD pile-up event selection"; 
  legcRTC->SetHeader(header.Data());

  TString hName;
  for(Int_t iTrig = 0; iTrig < triggersB->GetEntriesFast()-1; iTrig++){
		
    hName = "ratio";
    hName += ( (TObjString*) triggersB->At(iTrig) )->GetString();
    TH1 * hWithPS = (TH1*) (static_cast<TH1*> (trigWithPS.At(iTrig)));
    TH1 * hNoPS = (TH1*) (static_cast<TH1*> (trigNoPS.At(iTrig)));
    if (!hNoPS || !hWithPS) continue;
    TH1 *hRatio = (TH1*) hWithPS->Clone(hName);
    hRatio->Divide(hNoPS);
    hName = "ratioNoPS";
    hName += ( (TObjString*) triggersB->At(iTrig) )->GetString();

    if(iTrig==0){
      hRatio->SetMaximum(1.5);
      hRatio->SetMinimum(0.05);
      hRatio->SetLabelSize(0.02);
      hRatio->GetYaxis()->SetTitle("Accepted / All from Phys. Sel."); 
      hRatio->SetTitle("Physics Selection Cut"); 
      hRatio->Draw("E");
    }
    else{
      hRatio->Draw("ESAME");
    }
    trigRatio.AddAt(hRatio,iTrig);
  }
	 
  for(Int_t iTrig = 0; iTrig < triggersB->GetEntriesFast()-1; iTrig++){
    TH1 * histo = static_cast<TH1*> (trigRatio.At(iTrig));
    if (!histo) continue;
    legcRTC->AddEntry(histo,(( (TObjString*) triggersB->At(iTrig) )->GetString()).Data(),"l");
  }
  legcRTC->Draw("same");
		
  return c1;
}	

TCanvas *ProcessCanvasPhysSelCutCentrality(TObjArray *triggersB, TObjArray trigNoPS, TObjArray trigWithPS, Int_t trigNr, TString canvasName, TString *legendHeader){
	
  if(!triggersB || !legendHeader) return 0x0;

  TString cName = "c";
  cName += canvasName;
  TCanvas *c1 = new TCanvas(canvasName,cName,1200,900);
  SetCanvas(c1);
  c1->cd();
	
  Int_t const centBinMaxi = kCentBinMax - 1;
  TObjArray trigWithPSRatio(centBinMaxi);
  TLegend* legcRTC = new TLegend(0.2,0.15,0.50,0.40);
  legcRTC->SetHeader("V0 amplitude bins:");
  Int_t colorTab[3] = {kBlack,kRed,kBlue};
  TArrayI colorInd(centBinMaxi);
  for (Int_t i = 0; i < centBinMaxi; i++ ) {
    if(i<3) colorInd.AddAt(colorTab[i],i);
    else colorInd.AddAt(colorTab[2]+i-2,i);
  }

  TString hName;
	
  Float_t yMin = 0.05, yMax = 2;
	
  for(Int_t iCentBin = 0; iCentBin < centBinMaxi; iCentBin++){

    Int_t index = GetIndex(triggersB,trigNr,iCentBin);
    TH1 *hWithPS = static_cast<TH1*> (trigWithPS.At(index));
    TH1 *hNoPS = static_cast<TH1*> (trigNoPS.At(index));
    if (!hNoPS || !hWithPS) continue;
		
    hName = "ratio";
    hName += ( (TObjString*) triggersB->At(trigNr) )->GetString();
    TH1 *hRatio  = (TH1*) hWithPS->Clone(hName);
    hRatio->Divide(hNoPS);
    hRatio->SetLineColor(colorInd.At(iCentBin));
    if ( iCentBin == 0 ) {
      hRatio->SetMaximum(yMax);
      hRatio->SetMinimum(yMin);
      hRatio->SetLabelSize(0.02);
      hRatio->GetYaxis()->SetTitle("Accepted / All from Phys. Sel.");
      TString sTitle = "for ", sTitle2 = (( (TObjString*) triggersB->At(trigNr) )->GetString()).Data();
      if ( !sTitle2.IsNull() ) sTitle += sTitle2;
      else sTitle = "";
      hRatio->SetTitle(Form("Phys. Sel. %s - Multiplicity from V0 amplitude",sTitle.Data()));
      hRatio->Draw("E");
    }
    else{
      hRatio->Draw("ESAME");
    }
    trigWithPSRatio.AddAt(hRatio,iCentBin);
  }
	
  for ( Int_t centBin = 0; centBin < centBinMaxi; centBin++ ){
    TH1 *hRatio = static_cast<TH1*> (trigWithPSRatio.At(centBin));
    if ( !hRatio ) continue;
    legcRTC->AddEntry(hRatio,(legendHeader[centBin]).Data(),"l");
  }
  legcRTC->Draw("same");
	
	
  return c1;
}	

TCanvas *ProcessCanvasCentralityPercentile(TObjArray *triggersB, TObjArray trigNoPS, TObjArray trigWithPS, Int_t trigNr, TString canvasName, TString *legendHeader){
	
  if(!triggersB || !legendHeader) return 0x0;
		
  TString cName = "c";
  cName += canvasName;
  TCanvas *c1 = new TCanvas(canvasName,cName,1200,900);
  SetCanvas(c1,0);
  c1->cd();
	
  Int_t const centBinMaxi = 2;
  TObjArray trigRatio(centBinMaxi), trigRatioNoPS(centBinMaxi);
  TLegend* legcRTC = new TLegend(0.2,0.15,0.50,0.40);
  legcRTC->SetHeader("Physics Selection");
	
 Int_t colorTab[2] = {kRed,kBlue};
 TArrayI colorInd(centBinMaxi);
 for (Int_t i = 0; i < centBinMaxi; i++ ) {
   if(i<2) colorInd.AddAt(colorTab[i],i);
   else colorInd.AddAt(colorTab[1]+i-1,i);
 }
 
 TString hName;
 
 Float_t yMin = 0., yMax = 0.3;
 
 //process centrality bin 0-10% (centBin=3) and 60-80% (centBin=3) and compare it to 0-80% (centBin=1)	
 Int_t centBinMin = 2;
 TH1 *hWithPSAll = static_cast<TH1*> (trigWithPS.At(1*triggersB->GetEntriesFast()+trigNr));
 TH1 *hNoPSAll = static_cast<TH1*> (trigNoPS.At(1*triggersB->GetEntriesFast()+trigNr));
 if (!hNoPSAll || !hWithPSAll) return 0;
 
 for ( Int_t centBin = centBinMin; centBin < centBinMaxi+centBinMin; centBin++){
   
   if ( centBin > kCentBinMax ) continue;
   
   Int_t index = GetIndex(triggersB,trigNr,centBin);
   TH1 *hWithPS = static_cast<TH1*> (trigWithPS.At(index));
   TH1 *hNoPS = static_cast<TH1*> (trigNoPS.At(index));
   if (!hNoPS || !hWithPS) continue;
   
   hName = "ratio";
   hName += ( (TObjString*) triggersB->At(trigNr) )->GetString();
   TH1 *hRatio = (TH1*) hWithPS->Clone(hName);
   hRatio->Divide(hWithPSAll);
   hRatio->Scale(0.8);
   hRatio->SetLineColor(colorInd.At(centBin-centBinMin));
   hName = "ratioNoPS";
   hName += ( (TObjString*) triggersB->At(trigNr) )->GetString();
   TH1 *hRatioNoPS = (TH1*) (hNoPS->Clone(hName));
   hRatioNoPS->Divide(hNoPSAll);
   hRatioNoPS->Scale(0.8);
   hRatioNoPS->SetLineStyle(2);
   hRatioNoPS->SetLineColor(colorInd.At(centBin-centBinMin));

   if ( centBin == centBinMin ){
     hRatio->SetMaximum(yMax);
     hRatio->SetMinimum(yMin);
     hRatio->SetLabelSize(0.02);
     hRatio->GetYaxis()->SetTitle("Centrality percentile check"); 
     TString sTitle = "for ", sTitle2 = (( (TObjString*) triggersB->At(trigNr) )->GetString()).Data();
     if ( !sTitle2.IsNull() ) sTitle += sTitle2;
     else sTitle = "";
     hRatio->SetTitle(Form("Centrality percentile check %s - Multiplicity from V0 amplitude",sTitle.Data()));
     hRatio->Draw("E");
     hRatioNoPS->Draw("EPSAME");
   }
   else{
     hRatio->Draw("ESAME");
     hRatioNoPS->Draw("EPSAME");
   }
   trigRatio.AddAt(hRatio,centBin-centBinMin);
   trigRatioNoPS.AddAt(hRatioNoPS,centBin-centBinMin);
 }
	
  legcRTC->AddEntry(".","applied :","");
  for(Int_t centBin = centBinMin; centBin < centBinMaxi+centBinMin; centBin++){
    TH1 *hRatio = static_cast<TH1*> (trigRatio.At(centBin-centBinMin));
    if (!hRatio) continue;
    legcRTC->AddEntry(hRatio,(legendHeader[centBin]).Data(),"l");
  }
  legcRTC->AddEntry(".","not applied :","");
  for(Int_t centBin = centBinMin; centBin < centBinMaxi+centBinMin; centBin++){
    TH1 *hRatioNoPS = static_cast<TH1*> (trigRatioNoPS.At(centBin-centBinMin));
    if (!hRatioNoPS) continue;
    legcRTC->AddEntry(hRatioNoPS,(legendHeader[centBin]).Data(),"l");
  }
  legcRTC->Draw("same");
	
	
  return c1;
}	

TCanvas *ProcessCanvasTracksoverTrigger(TObjArray *triggersB, TObjArray trigSel, TObjArray trackTracker, TObjArray trackTrigger, TObjArray trackMatched, TObjArray trackAll, Int_t trigNr, Int_t centNr, TString canvasName,TString legendHeader){
	 
  if(!triggersB || trigNr<0 || centNr<0 ) return 0x0;
	
  Int_t index = GetIndex(triggersB,trigNr,centNr);

  TH1 *hTrackerPerB, *hTriggerPerB, *hMatchedPerB, *hAllTracksPerB, *hTrigSel, *histo;
		 
  TString hName, hNameBase;
  hNameBase =( (TObjString*) triggersB->At(trigNr) )->GetString();
		
  hTrigSel = static_cast<TH1*> (trigSel.At(index));
  if (!hTrigSel) return 0;

  hName = Form("hTrackerPer%s",hNameBase.Data());
  histo = static_cast<TH1*> (trackTracker.At(index));
  if (!histo) return 0;
  hTrackerPerB = (TH1*) histo->Clone(hName);
  if ( hTrackerPerB->GetEntries() ) hTrackerPerB->Divide(hTrigSel);
  hTrackerPerB->SetLineColor(kRed);
	 
  hName = Form("hTriggerPer%s",hNameBase.Data());
  histo = static_cast<TH1*> (trackTrigger.At(index));
  if (!histo) return 0;
  hTriggerPerB = (TH1*) histo->Clone(hName);
  if ( hTriggerPerB->GetEntries() > 0 ) hTriggerPerB->Divide(hTrigSel);
  hTriggerPerB->SetLineColor(kBlue);
	 
  hName = Form("hMatchedPer%s",hNameBase.Data());
  histo = static_cast<TH1*> (trackMatched.At(index));
  if (!histo) return 0;
  hMatchedPerB = (TH1*) histo->Clone(hName);
  if ( hMatchedPerB->GetEntries() > 0 ) hMatchedPerB->Divide(hTrigSel);
  hMatchedPerB->SetLineColor(kViolet);
	 
  hName = Form("hAllTracksPer%s",hNameBase.Data());
  histo = static_cast<TH1*> (trackAll.At(index));
  if (!histo) return 0;
  hAllTracksPerB = (TH1*) histo->Clone(hName);
  if ( hAllTracksPerB->GetEntries() > 0 ) hAllTracksPerB->Divide(hTrigSel);
 
  hAllTracksPerB->SetLineWidth(3);
  hAllTracksPerB->SetLineColor(kBlack);
  hAllTracksPerB->SetTitle(Form("Number of Tracks /%s %s",hNameBase.Data(),legendHeader.Data()));
  hAllTracksPerB->SetMinimum(0.0001);
  if ( hAllTracksPerB->GetEntries() == 0 ) hAllTracksPerB->SetMaximum(0.1);
  hAllTracksPerB->SetLabelSize(0.02);
	
  TString cName = "c";
  cName += canvasName;
  hNameBase = ( (TObjString*) triggersB->At(trigNr) )->GetString();
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

TCanvas *ProcessCanvasTrackMultB(TObjArray *triggersB, TObjArray trigSel, TObjArray trackTracker, TObjArray trackTrigger, TObjArray trackMatched, Int_t trigNr, Int_t centNr, TString canvasName,TString legendHeader){

  if( !triggersB || centNr < 0 || trigNr < 0 ) return 0x0;
	
  Int_t index = GetIndex(triggersB,trigNr,centNr);
  TString hNameBase = ( (TObjString*) triggersB->At(trigNr) )->GetString();

  TString cName = Form("c%s%s",canvasName.Data(),hNameBase.Data());
  TCanvas *cTrackMultB = new TCanvas(canvasName,cName,1200,900);
  SetCanvas(cTrackMultB,0);
  cTrackMultB->Divide(1,2);
  cTrackMultB->cd(1);
	
  TH1* hSumTriggerOverB, *hSumTrackerOverB, *hTrigSel, *hTracker, *hTrigger, *hMatched; 
  TString hName; 

  hTrigSel = static_cast<TH1*> (trigSel.At(index));
  if (!hTrigSel) return 0;
  hTracker = static_cast<TH1*> (trackTracker.At(index));
  if (!hTracker) return 0;
  hTrigger = static_cast<TH1*> (trackTrigger.At(index));
  if (!hTrigger) return 0;
  hMatched = static_cast<TH1*> (trackMatched.At(index));
  if (!hMatched) return 0;
  
  hName = Form("hSumTriggerOver%s",hNameBase.Data());
  hSumTriggerOverB = (TH1*) hTrigger->Clone(hName);
  hSumTriggerOverB->Add(hMatched);
  hSumTriggerOverB->Divide(hTrigSel);
  hName = Form("Sum of trigger tracks (matched + trigger-only) / # events in %s %s",hNameBase.Data(),legendHeader.Data());
  hSumTriggerOverB->SetTitle(hName);
  hSumTriggerOverB->SetLabelSize(0.04);
  hSumTriggerOverB->SetLineColor(kBlue);

  hName = Form("hSumTrackerOver%s",hNameBase.Data());
  hSumTrackerOverB = (TH1*) hTracker->Clone(hName);
  hSumTrackerOverB->Add(hMatched);
  hSumTrackerOverB->Divide(hTrigSel);
  hName = Form("Sum of tracker tracks (matched + tracker-only) / # events in %s %s",hNameBase.Data(),legendHeader.Data());
  hSumTrackerOverB->SetTitle(hName);
  //hSumTrackerOverCINT1B->LabelsOption("u");
  hSumTrackerOverB->SetLabelSize(0.04);
  hSumTrackerOverB->SetLineColor(kBlue);
      	
  hSumTriggerOverB->Draw("e");
  cTrackMultB->cd(2);
  hSumTrackerOverB->Draw("e");
	
  return cTrackMultB;
	
}

TCanvas *ProcessCanvasRatioTrackB(TObjArray *triggersB, TObjArray trigSel, TObjArray trackTracker, TObjArray trackTrigger, TObjArray trackMatched, Int_t trigNr, Int_t centNr, TString canvasName,TString legendHeader) {

  if(!triggersB || trigNr < 0 || centNr < 0 ) return 0x0;
  
  Int_t index = GetIndex(triggersB,trigNr,centNr);
  TString hNameBase =( (TObjString*) triggersB->At(trigNr) )->GetString();
  
  TString cName = Form("c%s%s",canvasName.Data(),hNameBase.Data());
  TCanvas *cRatioTrackB = new TCanvas(canvasName,cName,1200,900);
  SetCanvas(cRatioTrackB,0);
  
  TH1* hTrackerOverTriggerB, *hMatchedOverTriggerB, *hMatchedOverTrackerB, *hTrigSel, *hTracker, *hTrigger, *hMatched; 	
  hTrigSel = static_cast<TH1*> (trigSel.At(index));
  if (!hTrigSel) return 0;
  hTracker = static_cast<TH1*> (trackTracker.At(index));
  if (!hTracker) return 0;
  hTrigger = static_cast<TH1*> (trackTrigger.At(index));
  if (!hTrigger) return 0;
  hMatched = static_cast<TH1*> (trackMatched.At(index));
  if (!hMatched) return 0;
  
  TString hName = Form("hTrackerOverTrigger%s",hNameBase.Data());
  hTrackerOverTriggerB = (TH1*) hTracker->Clone(hName);
  hTrackerOverTriggerB->Divide(hTrigger);
  hName = Form("# tracker tracks / # trigger tracks in %s %s",hNameBase.Data(),legendHeader.Data());
  hTrackerOverTriggerB->SetTitle(hName);
  //hTrackerOverTriggerCINT1B->LabelsOption("u");
  hTrackerOverTriggerB->SetLabelSize(0.02);
  hTrackerOverTriggerB->SetLineColor(kBlue);
    
  hName = Form("hMatchedOverTrigger%s",hNameBase.Data());	
  hMatchedOverTriggerB = (TH1*) hMatched->Clone(hName);
  hMatchedOverTriggerB->Divide(hTrigger);
  hName = Form("# matched tracks / # trigger tracks in %s %s",hNameBase.Data(),legendHeader.Data());
  hMatchedOverTriggerB->SetTitle(hName);
  //hMatchedOverTriggerCINT1B->LabelsOption("u");
  hMatchedOverTriggerB->SetLabelSize(0.02);
  hMatchedOverTriggerB->SetLineColor(kBlue);
    
  hName = Form("hMatchedOverTracker%s",hNameBase.Data());
  hMatchedOverTrackerB = (TH1*) hMatched->Clone(hName);
  hMatchedOverTrackerB->Divide(hTracker);
  hName = Form("# matched tracks / # tracker tracks in %s %s",hNameBase.Data(),legendHeader.Data());
  hMatchedOverTrackerB->SetTitle(hName);
  //hMatchedOverTrackerCINT1B->LabelsOption("u");
  hMatchedOverTrackerB->SetLabelSize(0.02);
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

TCanvas *ProcessCanvasAsymMatched(TObjArray *triggersB, TObjArray trackPosMatched, TObjArray trackNegMatched, TObjArray trackAllMatched, Int_t trigNr, Int_t centNr, TString canvasName,TString legendHeader) {

  if(!triggersB || trigNr < 0 || centNr < 0 ) return 0x0;
  
  Int_t index = GetIndex(triggersB,trigNr,centNr);
  TString hName, hNameBase = (( (TObjString*) triggersB->At(trigNr) )->GetString());
  
  TString cName = Form("c%s%s",canvasName.Data(),hNameBase.Data());	
  TCanvas *cAsymMatched = new TCanvas(canvasName.Data(),cName,1200,900);
  SetCanvas(cAsymMatched,0);
  cAsymMatched->cd();
	
  TH1* hPosMatched, *hNegMatched, *hAllMatched;
 
  hPosMatched = static_cast<TH1*> (trackPosMatched.At(index));
  if (!hPosMatched) return 0;
  hNegMatched = static_cast<TH1*> (trackNegMatched.At(index));
  if (!hNegMatched) return 0;
  hAllMatched = static_cast<TH1*> (trackAllMatched.At(index));
  if (!hAllMatched) return 0;

  hName = Form("hAsyMatchedFor%s",hNameBase.Data());
  TH1 *hAsymMatched = (TH1*) hPosMatched->Clone(hName); 
  hAsymMatched->Add(hNegMatched,-1);
  hAsymMatched->Divide(hAllMatched);
  hAsymMatched->SetLineColor(kRed);
  hAsymMatched->SetMinimum(-0.3);
  hAsymMatched->SetMaximum(0.3);
  hAsymMatched->SetLabelSize(0.02);
  hName = Form("Matched tracks charge asymmetry for %s with acc. cuts %s",hNameBase.Data(),legendHeader.Data());
  hAsymMatched->SetTitle(hName);
  hAsymMatched->GetYaxis()->SetTitle("Charged tracks asymmetry");  
  hAsymMatched->Draw("EH");
	
  return cAsymMatched;

}

TCanvas *ProcessCanvasHighPtMuons(TObjArray *triggersB, TObjArray trigSel, TObjArray trackMatchedLowPt, TObjArray trackAllMatchedHighPt, Int_t trigNr, Int_t centNr, TString canvasName,TString legendHeader) {

  if ( !triggersB || trigNr < 0 || centNr < 0 ) return 0x0;
  
  Int_t index = GetIndex(triggersB,trigNr,centNr);
  TString hName, hNameBase = (( (TObjString*) triggersB->At(trigNr) )->GetString());
  
  TString cName = Form("c%s%s",canvasName.Data(),hNameBase.Data());	
  TCanvas *cHighPtMuons = new TCanvas(canvasName.Data(),cName,1200,900);
  SetCanvas(cHighPtMuons,0);
  cHighPtMuons->cd();
  TLegend* legcHPM;
  
  TH1* hRelMatchedLowPt, *hRelMatchedHighPt, *hTrigSel, *hMatchedLowPt, *hAllMatchedHighPt;  
  hTrigSel = static_cast<TH1*> (trigSel.At(index));
  if (!hTrigSel) return 0;
  hMatchedLowPt = static_cast<TH1*> (trackMatchedLowPt.At(index));
  if (!hMatchedLowPt) return 0;
  hAllMatchedHighPt = static_cast<TH1*> (trackAllMatchedHighPt.At(index));
  if (!hAllMatchedHighPt) return 0;

  hName = Form("hMatchedLowPtPer%s ",hNameBase.Data());
  hRelMatchedLowPt = (TH1*) hMatchedLowPt->Clone(hName);
  hRelMatchedLowPt->Divide(hTrigSel);
  hRelMatchedLowPt->SetLineColor(kBlue);
  hRelMatchedLowPt->SetTitle("");
  hName = Form("Ratio per %s ",hNameBase.Data());
  hRelMatchedLowPt->GetYaxis()->SetTitle(hName);
  hRelMatchedLowPt->SetMinimum(0.0001);
  hRelMatchedLowPt->SetLabelSize(0.02);
	 
  hName = Form("hMatchedHighPtPer%s ",hNameBase.Data());
  hRelMatchedHighPt = (TH1*) hAllMatchedHighPt->Clone(hName);
  hRelMatchedHighPt->Divide(hTrigSel);
  hRelMatchedHighPt->SetLineColor(kRed);
	 	 
  hName = Form("Number of matched track per %s (include Vtx and R_{Abs} cuts) %s",hNameBase.Data(),legendHeader.Data());
  hRelMatchedLowPt->SetTitle(hName);
  hRelMatchedLowPt->Draw("E");
  hRelMatchedHighPt->Draw("Esame");
	 
  legcHPM = new TLegend(0.60,0.45,0.98,0.65);
  //legcHPM->SetHeader(hName);
  legcHPM->AddEntry(".","Physics selection applied :","");	
  legcHPM->AddEntry(hRelMatchedLowPt," p_{T} > 1 GeV/c ","l");
  legcHPM->AddEntry(hRelMatchedHighPt," p_{T} >  2 GeV/c ","l");
  legcHPM->Draw("same");
	
  return cHighPtMuons;
}

TCanvas *ProcessCanvasBeamGasMatched ( TObjArray *triggersB, TObjArray trackBeamGasMatched, TObjArray trackBeamGasMatchedHighPt, TObjArray trackAllMatched, TObjArray trackMatchedHighPt, Int_t trigNr, Int_t centNr, TString canvasName,TString legendHeader) {

  if(!triggersB || trigNr < 0 || centNr < 0 )    return 0x0;

  Int_t index = GetIndex(triggersB,trigNr,centNr);
  TString hName, hNameBase = (( (TObjString*) triggersB->At(trigNr) )->GetString());
	
  TString cName = Form("c%s%s",canvasName.Data(),hNameBase.Data());	
  TCanvas *cBeamGasMatched = new TCanvas(canvasName.Data(),cName,1200,900);
  SetCanvas(cBeamGasMatched,0);
  cBeamGasMatched->cd();
		
  TH1* hBeamGasMatched, *hBeamGasMatchedHighPt, *hAllMatched, *hMatchedHighPt;

  hBeamGasMatched = static_cast<TH1*> (trackBeamGasMatched.At(index));
  if (!hBeamGasMatched) return 0;
  hBeamGasMatchedHighPt = static_cast<TH1*> (trackBeamGasMatchedHighPt.At(index));
  if (!hBeamGasMatchedHighPt) return 0;
  hAllMatched = static_cast<TH1*> (trackAllMatched.At(index));
  if (!hAllMatched) return 0;
  hMatchedHighPt = static_cast<TH1*> (trackMatchedHighPt.At(index));
  if (!hMatchedHighPt) return 0;

  hName = Form("hBeamGasMatchedPer%s ",hNameBase.Data());
  TH1 *hRelBeamGasMatched = (TH1*) hBeamGasMatched->Clone(hName);
  hRelBeamGasMatched->Divide(hAllMatched);
  hRelBeamGasMatched->SetLineColor(kBlack);
  hRelBeamGasMatched->SetMinimum(0.0);
  hRelBeamGasMatched->SetMaximum(1.1);
  hRelBeamGasMatched->SetLabelSize(0.02);
  
  hName = Form("hBeamGasMatchedHightPtPer%s ",hNameBase.Data());
  TH1 *hRelBeamGasMatchedHighPt = (TH1*) hBeamGasMatchedHighPt->Clone(hName);
  hRelBeamGasMatchedHighPt->Divide(hMatchedHighPt);
  hRelBeamGasMatchedHighPt->SetLineColor(kRed);

  hName = Form("Identified beam-gas tracks (pxDCA cuts) in matched tracks for %s",hNameBase.Data());
  if(!legendHeader.IsNull()) hName += Form(" %s",legendHeader.Data());
  hRelBeamGasMatched->SetTitle(hName);	 
  hRelBeamGasMatched->GetYaxis()->SetTitle("Relative beam-gas tracks");  
  hRelBeamGasMatched->Draw("EH");
  hRelBeamGasMatchedHighPt->Draw("EHsame");

  TLegend *leg = new TLegend(0.60,0.45,0.98,0.65);
  leg->AddEntry(".","Physics selection applied :","");	
  leg->AddEntry(hRelBeamGasMatched," All p_{T}","l");
  leg->AddEntry(hRelBeamGasMatchedHighPt," p_{T} >  2 GeV/c ","l");
  leg->Draw("same");
	
  return cBeamGasMatched;
}

Int_t GetIndex( TObjArray *triggers, Int_t trigNr, Int_t centNr ) {

  return ( centNr * triggers->GetEntriesFast() + trigNr );
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
    h1->SetLineWidth(2);
    
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

Bool_t GetTriggerLists(const char* triggerList, TString listFromContainer, TObjArray *triggersB, TObjArray *triggersShortName){
	
  //Get the trigger list from a file
  //The file should consist of a line for each trigger with the following layout:
  //        MB triggernameB 
  //     or MUONUNLIKE triggernameB 
  //     or NOSHOW triggernameB 
  //if filename is 0, then all the triggers stored are used	
  if( !triggersB || !triggersShortName) return kFALSE;
  Int_t const nColumn = 2;
  TObjArray* triggers[nColumn] = {triggersShortName, triggersB};
  
  TString trigSuffix[nColumn] = {"","B"};
  TString currTrigName = "";
  TObjArray* fullTriggerList[nColumn];
	
  for ( Int_t ibeam=0; ibeam<nColumn; ++ibeam ) {
    fullTriggerList[ibeam] = new TObjArray();
    fullTriggerList[ibeam]->SetOwner();
  }
  
  // Build trigger list (from file or use all trigger stored)
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
    TObjArray *triggersInContainer = listFromContainer.Tokenize(",");
    Int_t nTrig = 0;
    for ( Int_t iTrig = 0; iTrig < triggersInContainer->GetEntriesFast(); iTrig++ ) {
      currTrigName = triggersInContainer->At(iTrig)->GetName();
      Bool_t keep = kTRUE;
      TString rejectPatterns = "-E- -A- -C- WU UP SPI PHI EMC ZED TRUE SHM TPC BEAM 1A 1C";
      TObjArray* rejArray = rejectPatterns.Tokenize(" ");
      for ( Int_t ipat=0; ipat<rejArray->GetEntriesFast(); ipat++ ) {
        if ( currTrigName.Contains(rejArray->At(ipat)->GetName()) ) {
          keep = kFALSE;
          break;
        }
      }
      delete rejArray;
      if (!keep) continue;
      nTrig++;
      for (Int_t ibeam = 0; ibeam < nColumn; ibeam++) {
	fullTriggerList[ibeam]->AddLast( new TObjString(currTrigName) );
      }
    }
    //if no triggers are kept, then keep all of them
    if (nTrig == 0) {
      printf("INFO: no trigger selected over %d triggers: all triggers kept!!\n",
	     triggersInContainer->GetEntriesFast());
      for ( Int_t iTrig = 0; iTrig < triggersInContainer->GetEntriesFast(); iTrig++ ) {
	currTrigName = triggersInContainer->At(iTrig)->GetName();
	for (Int_t ibeam = 0; ibeam < nColumn; ibeam++) {
	  fullTriggerList[ibeam]->AddLast( new TObjString(currTrigName) );
	}
      }
    }
    if ( triggersInContainer ) delete triggersInContainer;
  }

  //
  // Select only existing triggers in container
  //
  TObjArray *triggersFromContainer = listFromContainer.Tokenize(",");
  TObjString* trigName = 0x0;
	
  TString selectAllTriggers[nColumn];
  for ( Int_t ibeam=0; ibeam<nColumn; ++ibeam ) selectAllTriggers[ibeam]= "";
  
  for ( Int_t itrig=0; itrig<fullTriggerList[0]->GetEntries(); ++itrig ) {
    Bool_t isBadTrig = kFALSE;
    for ( Int_t ibeam=0; ibeam<nColumn; ++ibeam ) {
      currTrigName = fullTriggerList[ibeam]->At(itrig)->GetName();
      //condition on trigger name from trigger list
      if ( ibeam == 0 && currTrigName.Contains("NOSHOW") ) {
	break;
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

Int_t GetRunNumber(TString filePath)
{
  /// Get run number from file path
  TObjArray* array = filePath.Tokenize("/");
  array->SetOwner();
  TString auxString = "";
  Int_t runNum = -1;
  for ( Int_t ientry=0; ientry<array->GetEntries(); ientry++ ) {
    auxString = array->At(ientry)->GetName();
    if ( auxString.IsDigit() && auxString.Length()>=6 && auxString.Length()<=9 ) {
      runNum = auxString.Atoi();
      break;
    }
  }
  delete array;

  if ( runNum < 0 ) {
    array = auxString.Tokenize("_");
    array->SetOwner();
    auxString = array->Last()->GetName();
    auxString.ReplaceAll(".root","");
    if ( auxString.IsDigit() ) runNum = auxString.Atoi();
    delete array;
  }

  return runNum;
}

TString GetRunList(const char *runList, TObjArray *runs){

  // list of runs to be analyzed
  TString selectRuns = "run:";
  
  if(runList) {
    // only the ones in the runList
    ifstream inFile(runList);
    if (!inFile.is_open()) {
      Error("PlotMuonQA","unable to open file %s", runList);
      return selectRuns;
    }
    
    TString currLine;
    while (!inFile.eof()) {
      currLine.ReadLine(inFile);
      if ( currLine.IsNull() ) continue;
      Int_t currRun = GetRunNumber(currLine);
      if (currRun<0) {
        Warning("PlotMuonQA","invalid run number: %s", currLine.Data());
        continue;
      }
      if(runs) runs->AddLast(new TObjString(Form("%d", currRun)));
      selectRuns += Form("%i,",currRun);
    }
    selectRuns.Remove(TString::kTrailing, ',');
    inFile.close();
    
  } else {
    // all runs
    cout<<"runList is not set"<<endl;
    if(runs) runs->AddLast(new TObjString("*"));
  }
  
  printf("selected runs from runlist %s: %s\n",runList, selectRuns.Data());
	
  return selectRuns;
}

Bool_t MergeOutputs(const char* inputList,const char* outFilename)
{
  ifstream inFile(inputList);
  if ( ! inFile.is_open()) {
    printf("Error: cannot find inputList %s\n", inputList);
    return kFALSE;
  }
  TString currLine = "";
  TFileMerger fileMerger;
  Int_t mergeType = ( TFileMerger::kRegular | TFileMerger::kAll | TFileMerger::kOnlyListed );
  fileMerger.AddObjectNames("MUON_QA");
  while ( ! inFile.eof() ) {
    currLine.ReadLine(inFile);
    if ( ! currLine.EndsWith(".root") ) continue;
    fileMerger.AddFile(currLine.Data());
  }
  inFile.close();
  if ( fileMerger.GetMergeList()->GetEntries() == 0 ) return kFALSE;
  fileMerger.OutputFile(outFilename,kTRUE,1); // needed when merging single files for specific directories
  fileMerger.PartialMerge(mergeType);

  return kTRUE;
}
