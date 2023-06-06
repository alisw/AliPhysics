#if !defined(__CINT__) || defined(__CLING__)
#include "AliMCEventHandler.h"
#include "AliESDInputHandler.h"
#include "AliAODInputHandler.h"
#include "AliDummyHandler.h"
#include "AliAnalysisAlien.h"
#include "AliAnalysisManager.h"

R__ADD_INCLUDE_PATH($ALICE_ROOT)

R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#include "AliAnalysisTaskGammaPythia.h"

#endif
#include "TGrid.h"
 
/* #include <ctime> */
/* #include "TGrid.h" */

//Analysis tasks
const Bool_t   saveManager         = kFALSE;//kFALSE;//
UInt_t         pSel                = AliVEvent::kAny;
const char    *kCentEst            = "V0M";

void runStandalonePythia8(Long64_t nEvents = 1)
{
  // since we will compile a class, tell root where to look for headers
#if !defined (__CINT__) || defined (__CLING__)
  gInterpreter->ProcessLine(".include $ROOTSYS/include");
  gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
#else
  gROOT->ProcessLine(".include $ROOTSYS/include");
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
#endif
  
  AliAnalysisManager* mgr         = new AliAnalysisManager("MCgen");
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
  
  cout << __LINE__ << endl;
  AliMCGenHandler* mcInputHandler = new AliMCGenHandler();
  mgr->SetMCtruthEventHandler(mcInputHandler);
  
#if !defined (__CINT__) || defined (__CLING__)
  AliGenerator* gener=reinterpret_cast<AliGenerator*>(gInterpreter->ExecuteMacro("AddMCGen_Pythia8_13TeV_Monash.C()"));
#else

#endif
  
  mcInputHandler->SetGenerator(gener);
  mcInputHandler->SetSeedMode(2); // check what this does

  // since we will compile a class, tell root where to look for headers
#if !defined (__CINT__) || defined (__CLING__)
  gInterpreter->LoadMacro("AliAnalysisTaskGammaPythia.cxx++g");
  AliAnalysisTaskGammaPythia *task = reinterpret_cast<AliAnalysisTaskGammaPythia*>(gInterpreter->ExecuteMacro("AddTask_GammaPythia.C()"));
#else
  gROOT->LoadMacro("AddTask_GammaPythia.C");
  AliAnalysisTask *taskA = AddTask_GammaPythia();		 
#endif

  (AliPythia8::Instance())->PrintStatistics();

  mgr->InitAnalysis();
  mgr->PrintStatus();
  mgr->EventLoop(nEvents);

}
