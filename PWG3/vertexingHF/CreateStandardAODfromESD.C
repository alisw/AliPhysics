#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TChain.h>
#include <TSystem.h>
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliAODHandler.h"
#include "AliAnalysisTaskESDfilter.h"
#include "AliAnalysisDataContainer.h"
#endif

void CreateStandardAODfromESD() 
{
 
  const char *inFileName = "AliESDs.root";
  const char *outFileName = "AliAOD.root";
  Bool_t writeKineToAOD = kTRUE;
  TString mode="local"; // "grid" 

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
  
  TChain *chain=0;
  if(mode=="local") { // local
    chain = new TChain("esdTree");
    // Steering input chain
    chain->Add(inFileName);
  } else if(mode=="grid") { // grid
    const char *collectionfile = "esd_coll1.xml";
    TGrid::Connect("alien:",0,0,"t") ;
    //Create an AliRunTagCuts and an AliEventTagCuts Object and impose some selection criteria
    AliRunTagCuts      *runCuts   = new AliRunTagCuts();
    AliEventTagCuts    *eventCuts = new AliEventTagCuts();
    AliLHCTagCuts      *lhcCuts   = new AliLHCTagCuts();
    AliDetectorTagCuts *detCuts   = new AliDetectorTagCuts();
    // eventCuts->SetMultiplicityRange(0,20000);
    //Create an AliTagAnalysis Object and chain the tags
    AliTagAnalysis   *tagAna = new AliTagAnalysis();
    tagAna->SetType("ESD");
    TAlienCollection *coll   = TAlienCollection::Open(collectionfile);
    TGridResult      *tagResult = coll->GetGridResult("",0,0);
    tagResult->Print();
    tagAna->ChainGridTags(tagResult);
    //Create a new esd chain and assign the chain that is returned by querying the tags
    chain = tagAna->QueryTags(runCuts,lhcCuts,detCuts,eventCuts);
  } else {
    printf("ERROR: mode has to be \"local\" or \"grid\"\n");
    return;
  }

  AliAnalysisManager *mgr  = new AliAnalysisManager("ESD to AOD", "Analysis Manager");
  
  // Input
  AliESDInputHandler* inpHandler = new AliESDInputHandler();
  mgr->SetInputEventHandler (inpHandler);
  
  // Output
  AliAODHandler* aodHandler   = new AliAODHandler();
  aodHandler->SetOutputFileName(outFileName);
  mgr->SetOutputEventHandler(aodHandler);
 
  // MC Truth
  AliMCEventHandler* mcHandler = new AliMCEventHandler();
  if(writeKineToAOD) mgr->SetMCtruthEventHandler(mcHandler);

  // Tasks
  
  // Filtering of MC particles (decays conversions etc)
  // this task is also needed to set the MCEventHandler
  // to the AODHandler, this will not be needed when
  // AODHandler goes to ANALYSISalice
  AliAnalysisTaskMCParticleFilter *kinefilter = new AliAnalysisTaskMCParticleFilter("Particle Filter");
  if(writeKineToAOD) mgr->AddTask(kinefilter);
  
  // Barrel Tracks
  AliAnalysisTaskESDfilter *filter = new AliAnalysisTaskESDfilter("Filter");
  mgr->AddTask(filter);

  // Muons
  AliAnalysisTaskESDMuonFilter *esdmuonfilter = new AliAnalysisTaskESDMuonFilter("ESD Muon Filter");
  mgr->AddTask(esdmuonfilter);

  AliESDtrackCuts* esdTrackCutsHF = new AliESDtrackCuts("AliESDtrackCuts", "Heavy flavour");
  esdTrackCutsHF->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);

  AliAnalysisFilter* trackFilter = new AliAnalysisFilter("trackFilter");
  trackFilter->AddCuts(esdTrackCutsHF);

  filter->SetTrackFilter(trackFilter);

  // Pipelining
  mgr->ConnectInput(filter,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(filter,0,mgr->GetCommonOutputContainer());
  mgr->ConnectInput(esdmuonfilter,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(esdmuonfilter,0,mgr->GetCommonOutputContainer());
  if(writeKineToAOD) {
    mgr->ConnectInput(kinefilter,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(kinefilter,0,mgr->GetCommonOutputContainer());
  }
  /*
  // before v4-17-Release  
  AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("cchain", TChain::Class(),
							   AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("tree", TTree::Class(),
							    AliAnalysisManager::kOutputContainer,
							    "default");
  mgr->ConnectInput(filter,0,cinput1);
  mgr->ConnectOutput(filter,0,coutput1);
  mgr->ConnectInput(esdmuonfilter,0,cinput1);
  mgr->ConnectOutput(esdmuonfilter,0,coutput1);
  if(writeKineToAOD) {
    mgr->ConnectInput(kinefilter,0,cinput1);
    mgr->ConnectOutput(kinefilter,0,coutput1);
  }
  */
  //
  // Run the analysis
  //
  mgr->InitAnalysis();
  mgr->PrintStatus();
  mgr->StartAnalysis(mode.Data(),chain);

  return;
}
