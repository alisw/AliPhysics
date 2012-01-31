/* $Id: AliEventStatsTask.cxx 35782 2009-10-22 11:54:31Z jgrosseo $ */

#include "AliEventStatsTask.h"

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>

#include <AliLog.h>
#include <AliESDEvent.h>
#include <AliHeader.h>

#include "AliPhysicsSelection.h"
//#include "AliBackgroundSelection.h"

ClassImp(AliEventStatsTask)

AliEventStatsTask::AliEventStatsTask(const char* opt) :
  AliAnalysisTaskSE("AliEventStatsTask"),
  fOutput(0),
  fOption(opt),
  fPhysicsSelection(0)
{
  //
  // Constructor. Initialization of pointers
  //

  // Define input and output slots here
  DefineOutput(1, TList::Class());
  
  AliLog::SetClassDebugLevel("AliEventStatsTask", AliLog::kWarning);
}

AliEventStatsTask::~AliEventStatsTask()
{
  //
  // Destructor
  //

  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor

  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }
}

void AliEventStatsTask::UserCreateOutputObjects()
{
  // create result objects and add to output list

  Printf("AliEventStatsTask::CreateOutputObjects");

  fOutput = new TList;
  fOutput->SetOwner();
  
  if (!fPhysicsSelection)
    fPhysicsSelection = new AliPhysicsSelection;
  
  fOutput->Add(fPhysicsSelection);
}

void AliEventStatsTask::UserExec(Option_t*)
{
  // process the event

  // post the data already here
  PostData(1, fOutput);

  AliESDEvent* esd = dynamic_cast<AliESDEvent*> (InputEvent());

  if (!esd)
  {
    AliError("ESD branch not available");
    return;
  }
  
  fPhysicsSelection->IsCollisionCandidate(esd);
}

void AliEventStatsTask::Terminate(Option_t *)
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
    fPhysicsSelection->SaveHistograms("physics_selection");
  }
    
  fout->Write();
  fout->Close();
  
  Printf("Writting result to event_stat.root");
}
