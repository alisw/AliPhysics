#include "TObject.h"
#include "TFile.h"
#include "AliPhysicsSelectionTask.h"
#include "AliPhysicsSelection.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliMultiInputEventHandler.h"
#include "AliAnalysisDataContainer.h"

ClassImp(AliPhysicsSelectionTask)

AliPhysicsSelectionTask::AliPhysicsSelectionTask() :
  AliAnalysisTaskSE("AliPhysicsSelectionTask"),
  fOutput(0),
  fOption(""),
  fUseSpecialOutput(kFALSE),
  fPhysicsSelection(0)
{
  //
  // Default event handler
  //
}

AliPhysicsSelectionTask::AliPhysicsSelectionTask(const char* opt) :
  AliAnalysisTaskSE("AliPhysicsSelectionTask"),
  fOutput(0),
  fOption(opt),
  fUseSpecialOutput(kFALSE),
  fPhysicsSelection(new AliPhysicsSelection("blablabla"))
{
  //
  // Constructor. Initialization of pointers
  //
  AliInputEventHandler* handler = dynamic_cast<AliInputEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (handler) {
    handler->SetEventSelection(fPhysicsSelection);
    AliInfo("Physics Event Selection enabled.");
  } else {
    AliError("No input event handler connected to analysis manager. No Physics Event Selection.");
  }

  //
  TString opts = opt;
  opts.ToLower();
  if (opts.Contains("specialoutput")) fUseSpecialOutput = kTRUE;

  // Define input and output slots here
  DefineOutput(1, TList::Class());
  
  AliLog::SetClassDebugLevel("AliPhysicsSelectionTask", AliLog::kWarning);
}

AliPhysicsSelectionTask::~AliPhysicsSelectionTask(){
  //
  // Destructor
  //

  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor

  if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fOutput;
    fOutput = 0;
  }
}

void AliPhysicsSelectionTask::UserCreateOutputObjects(){
  // create result objects and add to output list

  Printf("AliPhysicsSelectionTask::CreateOutputObjects");

  if (fUseSpecialOutput) OpenFile(1);

  fOutput = new TList;
  fOutput->SetOwner();
  
  if (!fPhysicsSelection)
    fPhysicsSelection = new AliPhysicsSelection;
  
  fOutput->Add(fPhysicsSelection);
  // All tasks must post data once for all outputs (AG)
  PostData(1, fOutput);
}

void AliPhysicsSelectionTask::UserExec(Option_t*){
  // process the event

  // AliPhysicsSelection::IsCollisionCandidate is called from the event handler
  // post the data here anyway!
  PostData(1, fOutput);
}

void AliPhysicsSelectionTask::FinishTaskOutput(){
// This gets called at the end of the processing on the worker. It allows dumping
// statistics printed by the physics selection object to the statistics message
// handled by the analysis manager.
   if (!fPhysicsSelection) return;
   fPhysicsSelection->FillStatistics();
//   fPhysicsSelection->Print("STAT");
}

void AliPhysicsSelectionTask::Terminate(Option_t *){
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput)
    Printf("ERROR: fOutput not available");
    
  if (fOutput)
  {
    fPhysicsSelection = dynamic_cast<AliPhysicsSelection*> (fOutput->FindObject("AliPhysicsSelection"));
  }

  TFile* fout = new TFile("event_stat.root", "RECREATE");

  if (fPhysicsSelection)
  {
    fPhysicsSelection->Print();
    fPhysicsSelection->SaveHistograms();
  }
    
  fout->Write();
  fout->Close();
  
  Printf("Writing result to event_stat.root");
}

AliPhysicsSelectionTask* AliPhysicsSelectionTask::AddTaskPhysicsSelection( const Bool_t mCAnalysisFlag, const Bool_t applyPileupCuts, const UInt_t deprecatedFlag2, const Bool_t useSpecialOutput) {
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPhysicsSelection", "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPhysicsSelection", "This task requires an input event handler");
    return NULL;
  }

  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();

  TString inputDataType = inputHandler->GetDataType(); // can be "ESD" or "AOD"

  // Configure analysis
  //===========================================================================
  AliPhysicsSelectionTask *task = new AliPhysicsSelectionTask("");
  task->SetUseSpecialOutput(useSpecialOutput); // RS: optionally use special output
  // this makes physics selection to work using AliMultiInputEventHandler
  if (inputHandler && (inputHandler->IsA() == AliMultiInputEventHandler::Class())) {
    AliMultiInputEventHandler *multiInputHandler=(AliMultiInputEventHandler*)inputHandler;
    AliInputEventHandler *ih = multiInputHandler->GetFirstInputEventHandler();
    if (!ih) {
      ::Error("AddTaskPhysicsSelection","ESD or AOD input handler is missing");
      return NULL;
    }
    ih->SetEventSelection(multiInputHandler->GetEventSelection());
    inputDataType = ih->GetDataType(); // can be "ESD" or "AOD"
  }

  mgr->AddTask(task);

  AliPhysicsSelection* physSel = task->GetPhysicsSelection();
  if (mCAnalysisFlag)  { physSel->SetAnalyzeMC(); }
  if (deprecatedFlag2) { ::Error("AddTaskPhysicsSelection", "ComputeBG functionality is deprecated"); }

  physSel->ApplyPileupCuts(applyPileupCuts);

  AliAnalysisDataContainer *cinput0 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("cstatsout",
                TList::Class(),
                AliAnalysisManager::kOutputContainer,
                "EventStat_temp.root");
  //
  if (useSpecialOutput) { coutput1->SetSpecialOutput(); } // RS: optionally use special output
  //
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,coutput1);

  return task;
}

