/* $Id$ */

#include "AliPhysicsSelectionTask.h"

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>

#include <AliLog.h>
#include <AliESDEvent.h>
#include <AliHeader.h>

#include "AliPhysicsSelection.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"

//#include "AliBackgroundSelection.h"

ClassImp(AliPhysicsSelectionTask)

AliPhysicsSelectionTask::AliPhysicsSelectionTask() :
  AliAnalysisTaskSE("AliPhysicsSelectionTask"),
  fOutput(0),
  fOption(""),
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
  fPhysicsSelection(new AliPhysicsSelection())
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
  // Define input and output slots here
  DefineOutput(1, TList::Class());
  fBranchNames = "ESD:AliESDRun.,AliESDHeader.,AliMultiplicity.,AliESDVZERO.,"
                 "AliESDZDC.,SPDVertex.,PrimaryVertex.,TPCVertex.,Tracks,SPDPileupVertices";
  
  AliLog::SetClassDebugLevel("AliPhysicsSelectionTask", AliLog::kWarning);
}

AliPhysicsSelectionTask::~AliPhysicsSelectionTask()
{
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

void AliPhysicsSelectionTask::UserCreateOutputObjects()
{
  // create result objects and add to output list

  Printf("AliPhysicsSelectionTask::CreateOutputObjects");

  fOutput = new TList;
  fOutput->SetOwner();
  
  if (!fPhysicsSelection)
    fPhysicsSelection = new AliPhysicsSelection;
  
  fOutput->Add(fPhysicsSelection);
  // All tasks must post data once for all outputs (AG)
  PostData(1, fOutput);
}

void AliPhysicsSelectionTask::UserExec(Option_t*)
{
  // process the event

  // AliPhysicsSelection::IsCollisionCandidate is called from the event handler
  // post the data here anyway!
  PostData(1, fOutput);
}

void AliPhysicsSelectionTask::FinishTaskOutput()
{
// This gets called at the end of the processing on the worker. It allows dumping
// statistics printed by the physics selection object to the statistics message
// handled by the analysis manager.
   if (fPhysicsSelection) fPhysicsSelection->Print("STAT");
}

void AliPhysicsSelectionTask::Terminate(Option_t *)
{
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
  
  Printf("Writting result to event_stat.root");
}
