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
  //
  TString opts = opt;
  opts.ToLower();
  if (opts.Contains("specialoutput")) fUseSpecialOutput = kTRUE;

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

  if (fUseSpecialOutput) OpenFile(1);

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
  
  Printf("Writing result to event_stat.root");
}

void AliPhysicsSelectionTask::NotifyRun(){
  if (fPhysicsSelection->IsMC()) return;
  TObject* prodInfoData = fInputHandler->GetUserInfo()->FindObject("alirootVersion");
  TString filePath;
  if (prodInfoData) {
    // take filePath from UserInfo - available only from ~LHC12d period
    TString str(prodInfoData->GetTitle());
    TObjArray* tokens = str.Tokenize(";");
    for (Int_t i=0;i<=tokens->GetLast();i++) {
      TObjString* stObj = (TObjString*) tokens->At(i);
      TString s = stObj->GetString();
      if (s.Contains("OutputDir")) {
        filePath = s;
        break;
      }
    }
    delete tokens;
  } else {
    // guess name from the input filename
    // may be a problem for local analysis
    filePath = fInputHandler->GetTree()->GetCurrentFile()->GetName();
  }

  TString passName="";

  TObjArray* tokens = filePath.Tokenize("/");
  for (Int_t i=0;i<=tokens->GetLast();i++) {
    TObjString* stObj = (TObjString*) tokens->At(i);
    TString s = stObj->GetString();
    if (s.Contains("pass")) {
      passName = s;
      break;
    }
  }
  delete tokens;

  if (!passName.Contains("pass")){
    AliError(" ******** Failed to find reconstruction pass name *********");
    AliError(" ******** Default parameters loaded: parameters unreliable ******");
    AliError("      --> If these are MC data: please set kTRUE first argument of AddTaskPhysicsSelection");
    AliError("      --> If these are real data: ");
    AliError("          (a) please insert pass number inside the path of your local file OR");
    AliError("          (b) specify reconstruction pass number when adding physics selection task");
    AliError(" Using default pass parameters for physics selection (PS will probably fail for LHC10h pass1 and pass2 data)!");
  }
  fPhysicsSelection->SetPassName(passName);
}
