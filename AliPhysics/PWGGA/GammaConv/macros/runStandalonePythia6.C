// runmcgenThermalModel.C
// =====================
// This macro runs jet analysis on pythia events with AliMCGen framework
//
// Author: M. Verweij
 
#include <ctime>
#include "TGrid.h"

//Analysis tasks
const Int_t iEmcalSetup  = 0;
const Int_t iJetTagger   = 0;
const Int_t iJetMass     = 0;
const Int_t iJetMassBkg  = 0;

const Bool_t   saveManager         = kFALSE;//kFALSE;//
UInt_t         pSel                = AliVEvent::kAny;
const char    *kCentEst            = "V0M";

void runStandalonePythia6(Long64_t nEvents = 50000)
{

	AliAnalysisManager *mgr = new AliAnalysisManager("MCGenThermalModel");

	// Handlers
	AliDummyHandler *dumH = new AliDummyHandler();
	//mgr->SetInputEventHandler(dumH);

	// AliDummyHandler *dumHM = static_cast<AliDummyHandler*>mgr->GetInputEventHandler();
	// if(dumHM) Printf("Found dummy handler");
	// if(!dumHM->GetEvent()) Printf("handler has no event");
	AliESDEvent *esdE = new AliESDEvent();
	esdE->CreateStdContent();
	AliESDVertex *vtx = new AliESDVertex(0.,0.,100);
	vtx->SetName("VertexTracks");
	vtx->SetTitle("VertexTracks");
	esdE->SetPrimaryVertexTracks(vtx);
	if(esdE->GetPrimaryVertex()) Printf("vtx set");
	dumH->SetEvent(esdE);
	// if(dumHM->GetEvent()) Printf("handler has event");
	mgr->SetInputEventHandler(dumH);

	AliMCGenHandler* mcInputHandler = new AliMCGenHandler();  
	mgr->SetMCtruthEventHandler(mcInputHandler);

	// Generator
	// gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddMCGenHijing.C");
	// AliGenerator* gener = AddMCGenHijing();
	// mcInputHandler->SetGenerator(gener);
	// mcInputHandler->SetSeedMode(2);
	gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddMCGenPythia.C");
	
	AliGenerator* gener = AddMCGenPythia(2760., 16., 21.);
	((AliGenPythiaPlus*)gener)->SetEventListRange(-1, -1);
// 	((AliGenPythiaPlus*)gener)->SetProcess(kPyJets);
// 	((AliGenPythiaPlus*)gener)->SetPtHard(5.,7.);

	mcInputHandler->SetGenerator(gener);
	mcInputHandler->SetSeedMode(2);

  // Analysis Tasks
	gROOT->LoadMacro("$ALICE_PHYSICS/PWGGA/GammaConv/macros/AddTask_GammaPureMC.C");
	AliAnalysisTask *taskA = AddTask_GammaPureMC();		

	// Run analysis
	if(saveManager) {
		TFile *fileMgr = new TFile("AnalysisManager.root","RECREATE");
		mgr->Write();

		fileMgr->Write();
		fileMgr->Close();

		runUseSavedManager("AnalysisManager.root",nEvents);
	}
	else {
		mgr->InitAnalysis();
		mgr->PrintStatus();
		mgr->EventLoop(nEvents);
	}
	(AliPythia6::Instance())->PrintStatistics();

}

void runUseSavedManager(TString fileMgr = "AnalysisManager.root", Int_t nEvents = 1000, Bool_t useGrid = kFALSE) {

  Printf("run saved");

  LoadLibs();

  TFile *f = new TFile(fileMgr.Data());
  AliAnalysisManager *mgr = dynamic_cast<AliAnalysisManager*>f->Get("MCGenThermalModel");
  if(!mgr) {
    Printf("Did not find manager");
    return;
  }

  // AliDummyHandler *dumH = static_cast<AliDummyHandler*>mgr->GetInputEventHandler();
  // AliESDEvent *esdE = new AliESDEvent();
  // esdE->CreateStdContent();
  // AliESDVertex *vtx = new AliESDVertex(0.,0.,100);
  // vtx->SetName("VertexTracks");
  // vtx->SetTitle("VertexTracks");
  // esdE->SetPrimaryVertexTracks(vtx);
  // if(esdE->GetPrimaryVertex()) Printf("vtx set");
  // dumH->SetEvent(esdE);

  mgr->SetDebugLevel(11);
  mgr->InitAnalysis();
  mgr->PrintStatus();
  mgr->EventLoop(nEvents);
}


