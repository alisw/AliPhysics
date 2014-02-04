// $Id$
//
// Task to debug problems
//
// Author: C.Loizides

#include "AliEmcalDebugTask.h"
#include <TClonesArray.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TSystem.h>
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"

ClassImp(AliEmcalDebugTask)

//________________________________________________________________________
AliEmcalDebugTask::AliEmcalDebugTask() : 
  AliAnalysisTaskSE(),
  fId(0),
  fFileTest(),
  fPrintEnv(0),
  fOutput(0),
  fFileName(),
  fRand(0)
{
  // Constructor.
}

//________________________________________________________________________
AliEmcalDebugTask::AliEmcalDebugTask(const char *name) : 
  AliAnalysisTaskSE(name),
  fId(0),
  fFileTest(),
  fPrintEnv(0),
  fOutput(0),
  fFileName(),
  fRand(0)
{
  // Constructor.

  DefineOutput(1, TList::Class()); 
  fBranchNames = "ESD:AliESDHeader.,AliESDRun.,Tracks";
}

//________________________________________________________________________
AliEmcalDebugTask::~AliEmcalDebugTask()
{
  // Destructor.
}

//________________________________________________________________________
void AliEmcalDebugTask::UserCreateOutputObjects()
{
  // Create user objects

  fOutput = new TList();
  fOutput->SetOwner();

  TRandom3 r(0);
  fRand = r.Integer(kMaxUInt);
  fOutput->Add(new TNamed(Form("%u",fId),Form("%u",fRand)));

  AliInfo(Form("AliEmcalDebug: %u %u",fId,fRand));
  if (fPrintEnv)
    gSystem->Exec("env");

  PostData(1, fOutput); 
}

//________________________________________________________________________
void AliEmcalDebugTask::UserExec(Option_t *) 
{
  // Main loop, called for each event.

  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
  if (!am) {
    AliError("Manager zero, returning");
    return;
  }

  TString filename;

  TTree *t = am->GetTree();
  if (t) {
    TFile *f = t->GetCurrentFile();
    if (f) {
      filename = f->GetName();
    }
  }

  if (filename==fFileName)
    return;

  if (fFileTest.Length()>0) {
    if (!fFileName.Contains(fFileTest)) {
      AliError(Form("Filename %s does not contain %s", fFileName.Data(), fFileTest.Data()));
      return;
    }
  }
  fFileName = filename;

  AliInfo(Form("New file: %s", fFileName.Data()));
  fOutput->Add(new TNamed(Form("%u:%u",fId,fRand),fFileName.Data()));
}
