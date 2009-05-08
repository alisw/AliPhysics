#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TChain.h>
#include <TSystem.h>
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliAODHandler.h"
#include "AliAnalysisTaskESDfilter.h"
#include "AliAnalysisDataContainer.h"
#endif

void CreateAODfromESD(const char *inFileName = "AliESDs.root",
		      const char *outFileName = "AliAODs.root",
		      Bool_t bKineFilter = kTRUE) 
{
  
    gSystem->Load("libTree");
    gSystem->Load("libGeom");
    gSystem->Load("libPhysics");
    gSystem->Load("libVMC");
    gSystem->Load("libSTEERBase");
    gSystem->Load("libESD");
    gSystem->Load("libAOD");
    
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
    gSystem->Load("libPWG3muon");

    TChain *chain = new TChain("esdTree");
    // Steering input chain
    chain->Add(inFileName);
    AliAnalysisManager *mgr  = new AliAnalysisManager("ESD to AOD", "Analysis Manager");

    // Input
    AliESDInputHandler* inpHandler = new AliESDInputHandler();
    inpHandler->SetReadTags();
    mgr->SetInputEventHandler  (inpHandler);
    // Output
    AliAODHandler* aodHandler   = new AliAODHandler();
    aodHandler->SetOutputFileName(outFileName);
    mgr->SetOutputEventHandler(aodHandler);

    // MC Truth
    if(bKineFilter){
	AliMCEventHandler* mcHandler = new AliMCEventHandler();
	mgr->SetMCtruthEventHandler(mcHandler);
    }


    // Tasks
    // Filtering of MC particles (decays conversions etc)
    // this task is also needed to set the MCEventHandler
    // to the AODHandler, this will not be needed when
    // AODHandler goes to ANALYSISalice
    AliAnalysisTaskMCParticleFilter *kinefilter = new AliAnalysisTaskMCParticleFilter("Particle Filter");
    if (bKineFilter) mgr->AddTask(kinefilter);
    
    // Barrel Tracks
    AliAnalysisTaskESDfilter *filter = new AliAnalysisTaskESDfilter("Filter");
    mgr->AddTask(filter);
    // Muons
    AliAnalysisTaskESDMuonFilter *esdmuonfilter = new AliAnalysisTaskESDMuonFilter("ESD Muon Filter");
    mgr->AddTask(esdmuonfilter);

    // Cuts on primary tracks
    AliESDtrackCuts* esdTrackCutsL = new AliESDtrackCuts("AliESDtrackCuts", "Standard");
    esdTrackCutsL->SetMinNClustersTPC(50);
    esdTrackCutsL->SetMaxChi2PerClusterTPC(3.5);
    esdTrackCutsL->SetMaxCovDiagonalElements(2, 2, 0.5, 0.5, 2);
    esdTrackCutsL->SetRequireTPCRefit(kTRUE);
    esdTrackCutsL->SetMaxDCAToVertexXY(3.0);
    esdTrackCutsL->SetMaxDCAToVertexZ(3.0);
    esdTrackCutsL->SetDCAToVertex2D(kTRUE);
    esdTrackCutsL->SetRequireSigmaToVertex(kFALSE);
    esdTrackCutsL->SetAcceptKinkDaughters(kFALSE);
    // ITS stand-alone tracks
    AliESDtrackCuts* esdTrackCutsITSsa = new AliESDtrackCuts("AliESDtrackCuts", "ITS stand-alone");
    esdTrackCutsITSsa->SetRequireITSStandAlone(kTRUE);

    AliAnalysisFilter* trackFilter = new AliAnalysisFilter("trackFilter");
    trackFilter->AddCuts(esdTrackCutsL);
    trackFilter->AddCuts(esdTrackCutsITSsa);

    // Cuts on V0s
    AliESDv0Cuts*   esdV0Cuts = new AliESDv0Cuts("AliESDv0Cuts", "Standard pp");
    esdV0Cuts->SetMinRadius(0.2);
    esdV0Cuts->SetMaxRadius(100);
    esdV0Cuts->SetMinDcaPosToVertex(0.05);
    esdV0Cuts->SetMinDcaNegToVertex(0.05);
    esdV0Cuts->SetMaxDcaV0Daughters(0.5);
    esdV0Cuts->SetMinCosinePointingAngle(0.99);
    AliAnalysisFilter* v0Filter = new AliAnalysisFilter("v0Filter");
    v0Filter->AddCuts(esdV0Cuts);


//
    filter->SetTrackFilter(trackFilter);
    filter->SetV0Filter(v0Filter);


//  Create AOD Tags
    AliAnalysisTaskTagCreator* tagTask = new AliAnalysisTaskTagCreator("AOD Tag Creator");
    mgr->AddTask(tagTask);

    // Pipelining
    AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();    
    AliAnalysisDataContainer *coutput1 = mgr->GetCommonOutputContainer();
    
    
    AliAnalysisDataContainer *coutputT
	= mgr->CreateContainer("cTag",  TTree::Class(), AliAnalysisManager::kOutputContainer, "AOD.tag.root");

    coutput1->SetSpecialOutput();
    coutputT->SetSpecialOutput();
    
    if(bKineFilter) {
	mgr->ConnectInput  (kinefilter,     0, cinput1  );
	mgr->ConnectOutput (kinefilter,     0, coutput1 );
    }

    mgr->ConnectInput (filter, 0, cinput1 );
    mgr->ConnectOutput(filter, 0, coutput1);

    mgr->ConnectInput (esdmuonfilter, 0, cinput1 );
//    mgr->ConnectOutput(esdmuonfilter, 0, coutput1);

    mgr->ConnectInput (tagTask, 0, cinput1);
    mgr->ConnectOutput(tagTask, 1, coutputT);

    //
    // Run the analysis
    //
    mgr->InitAnalysis();
    mgr->PrintStatus();
    mgr->StartAnalysis("local", chain);
}

