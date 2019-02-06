// runStandaloneCocktail.C
// =====================
// This macro runs LMee cocktail generation with AliMCGen framework
// orginal author: M. Verweij
// modified by: F.Bock, L.Altenkaemper, N.Schmidt, O. Vazquez Doce

#include <ctime>
#include "TGrid.h"

const Bool_t saveManager = kFALSE;

void runStandaloneLMeeCocktail( Long64_t nEvents        = 400,
                                Int_t CollisionSystem  = 200,      // 200: pp 7TeV
                                Int_t motherSelect      = 63,       // pi0, eta, etaprime, rho, omega, phi. Explanation in AliPhysics/PWG/Cocktail/AliGenEMCocktailV2.h
                                Int_t decayMode         = 3,        // 0: kAll, 1: kGammaEM, 2: kElectronEM, 3: kDiElectronEM
                                Int_t numberOfParticles = 1000,
                                Int_t centrality        = 0,
                                Double_t minPt          = 0.,
                                Double_t maxPt          = 30,
                                Int_t pythiaErrorTolerance  = 2000,
                                Bool_t ExternalDecayer      = 1,    //0 pythia, 1 exodus
                                Bool_t decayLongLived       = 1
){
    AliAnalysisManager* mgr         = new AliAnalysisManager("MCGenEMCocktail");
    AliDummyHandler*    dumH        = new AliDummyHandler();

    // create blank ESD event
    AliESDEvent *esdE               = new AliESDEvent();
    esdE->CreateStdContent();
    AliESDVertex *vtx               = new AliESDVertex(0.,0.,100);
    vtx->SetName("VertexTracks");
    vtx->SetTitle("VertexTracks");
    esdE->SetPrimaryVertexTracks(vtx);
    if(esdE->GetPrimaryVertex()) Printf("vtx set");
    dumH->SetEvent(esdE);
    mgr->SetInputEventHandler(dumH);

    // greate MC input Handler
    AliMCGenHandler* mcInputHandler = new AliMCGenHandler();
    mgr->SetMCtruthEventHandler(mcInputHandler);

    // cocktail generator
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/Cocktail/macros/AddMCEMCocktailV2.C");
    AliGenerator* gener             = 0x0;
    gener = AddMCEMCocktailV2(CollisionSystem,centrality,decayMode,motherSelect,"$ALICE_PHYSICS/PWGDQ/dielectron/files/pp13TeV.root","13TeV_default",numberOfParticles,minPt,maxPt,pythiaErrorTolerance,ExternalDecayer,decayLongLived);

    mcInputHandler->SetGenerator(gener);
    mcInputHandler->SetSeedMode(2);


    // Analysis Tasks
    //---------------
    // LMee:
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/AddTask_LMeeCocktailMC.C");
    // Task parameters:
    Float_t MaxAcceptEta=0.8;
    Float_t MinAcceptPt=0.2;
    Bool_t WriteTTree = kFALSE;
    Int_t ResolType=2; // 1 = RunI, 2 = RunII
    Int_t ALTweightType = 1; // Weight type for alternative histos. 0=pt efficiency weight, 1=HM simulation.
    AliAnalysisTask *taskA = AddTask_LMeeCocktailMC(CollisionSystem,MaxAcceptEta,MinAcceptPt,WriteTTree,ResolType,ALTweightType);


    // give some indication everything is running fine
    mgr->SetUseProgressBar(1, 10);

    // Run analysis
    if(saveManager) {
        TFile *fileMgr = new TFile("AnalysisManager.root","RECREATE");
        mgr->Write();
        fileMgr->Write();
        fileMgr->Close();
        runUseSavedManager("AnalysisManager.root",nEvents);
    } else {
        mgr->InitAnalysis();
        mgr->PrintStatus();
        mgr->EventLoop(nEvents);
    }
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



