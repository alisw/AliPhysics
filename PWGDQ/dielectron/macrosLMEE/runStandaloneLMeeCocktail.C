// Adapted by Oton Vazquez Doce, from original:
//  runmcgenThermalModel.C
//  This macro runs jet analysis on pythia events with AliMCGen framework
//  Author: M. Verweij
 
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

// Define here number of events
void runStandaloneLMeeCocktail(Long64_t nEvents = 10)
{

	AliAnalysisManager *mgr = new AliAnalysisManager("MCGenEMCocktail");

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

  //=========================== Generator settings ==============================================
  Int_t selectedMothers   = 63; // pi0, eta, rho0, omega, eta', phi (111111)_2
    // Full Example (
    // (62591)_10 = (1111010001111111)_2
    // includes: pi0, eta, rho0, omega, etaprime, phi, J/psi, Delta+, Delta0, Rho+, Rho-, K0*
  Int_t numberOfParticles = 1000;
  Int_t collisionSystem  = 200;
  Int_t centrality        = 0;
  Int_t decayMode         = 3; // ee
  Int_t selectedMothers   = 63; // pi0, eta, rho0, omega, eta', phi
  Double_t minPt          = 0.;
  Double_t maxPt          = 20;
  Int_t pythiaErrorTolerance  = 2000;
  Bool_t ExternalDecayer      = 1; //0 pythia, 1 exodus
  Bool_t decayLongLived       = 0;

  // cocktail generator
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/Cocktail/macros/AddMCEMCocktailV2.C");
  //AliGenerator* gener = AddMCEMCocktailV2(collisionSystem,centrality,decayMode,selectedMothers,"$ALICE_PHYSICS/PWG/Cocktail/parametrisations/pp_7TeV_default.root","7TeV_default",numberOfParticles,minPt,maxPt,pythiaErrorTolerance,ExternalDecayer); // 7TeV
  AliGenerator* gener = AddMCEMCocktailV2(collisionSystem,centrality,decayMode,selectedMothers,"$ALICE_PHYSICS/PWGDQ/dielectron/files/pp_13TeV_defalt_parametrizations.root","13TeV_default",numberOfParticles,minPt,maxPt,pythiaErrorTolerance,ExternalDecayer); //13TeV
  

  mcInputHandler->SetGenerator(gener);
  mcInputHandler->SetSeedMode(2);

  // Analysis Tasks
  //---------------
  // LMee:
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/AddTask_LMeeCocktailMC.C");

  // Task parameters:
  Float_t MaxY=0.8;
  Float_t MinPt=0.2;
  if(collisionSystem>=400) MinPt=0.4;
  Bool_t WriteTTree = kFALSE;
  // MinPt=0.2;

  AliAnalysisTask *taskA = AddTask_LMeeCocktailMC(collisionSystem,MaxY,MinPt,WriteTTree);   
  
	// Run analysis
	//-------------
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
// 	(AliPythia6::Instance())->PrintStatistics();
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

  mgr->SetDebugLevel(11);
  mgr->InitAnalysis();
  mgr->PrintStatus();
  mgr->EventLoop(nEvents);
}


