//
// Class AliRsnReaderTask
//
// An AnalysisTask object to convert any kind of source event type (ESD/AOD/MC)
// into the RSN internal format (AliRsnEvent).
// The output of this task is a TTree with converted events, which is saved in a file
// and can then be processed as many times as desired, to build invariant mass spectra.
// ---
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#include "AliLog.h"

#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODInputHandler.h"
#include "AliAODHandler.h"

#include "AliRsnEvent.h"
#include "AliRsnReader.h"
#include "AliRsnPID.h"
#include "AliRsnReaderTask.h"

ClassImp(AliRsnReaderTask)

//_____________________________________________________________________________
AliRsnReaderTask::AliRsnReaderTask() :
    AliRsnAnalysisTaskBase(),
    fRsnTree(0x0)
{
//=========================================================
// Default constructor (not recommended)
//=========================================================
}

//_____________________________________________________________________________
AliRsnReaderTask::AliRsnReaderTask(const char *name) :
    AliRsnAnalysisTaskBase(name),
    fRsnTree(0x0)
{
//=========================================================
// Working constructor (recommended)
//=========================================================
  InitIOVars();
  DefineOutput(0, TTree::Class());
}

//_____________________________________________________________________________
AliRsnReaderTask::AliRsnReaderTask(const AliRsnReaderTask& obj) :
    AliRsnAnalysisTaskBase(obj),
    fRsnTree(0x0)
{
//=========================================================
// Copy constructor (not recommended)
//=========================================================
}

//_____________________________________________________________________________
AliRsnReaderTask& AliRsnReaderTask::operator= (const AliRsnReaderTask& /*obj*/)
{
//=========================================================
// Assignment operator (not recommended)
//=========================================================

  AliInfo("Not implemented. Avoid using the assignment operator");
  return *this;
}

//_____________________________________________________________________________
void AliRsnReaderTask::InitIOVars()
{
//=========================================================
// Init values for input and output data
//=========================================================
  AliDebug(AliLog::kDebug, "<-");
  AliRsnAnalysisTaskBase::InitIOVars();
  AliDebug(AliLog::kDebug, "->");
}

//_____________________________________________________________________________
void AliRsnReaderTask::CreateOutputObjects()
{
//=========================================================
// Create the output container
//=========================================================

  AliDebug(1, "Creating output objects");
  OpenFile(0);
  fRsnTree = new TTree("aodTree", "AliRsnEvents");

  fRSN[0] = new AliRsnEvent();
  fRSN[0]->SetName("rsnEvents");
  fRSN[0]->Init();
  fRsnTree->Branch("rsnEvents","AliRsnEvent",&fRSN[0]);
//     AddAODBranch("AliRsnEvent", &fRsnEvent);
}

//_____________________________________________________________________________
void AliRsnReaderTask::LocalInit()
{
//=========================================================
// Initialization
//=========================================================

  AliDebug(1, "Initializing");
}

Bool_t AliRsnReaderTask::Notify()
{
  return kTRUE;
}

//_____________________________________________________________________________
void AliRsnReaderTask::Exec(Option_t */*option*/)
{
//=========================================================
// Loops on input container to store data of all tracks.
// Uses the AliRsnReader methods to save them in the output.
//=========================================================

  fRSN[0] = GetRsnEventFromInputType(0);
  if (!fRSN[0]) return;

  AliDebug(1, Form("Collected %d tracks", fRSN[0]->GetMultiplicity()));

  fRsnTree->Fill();
  PostData(0, fRsnTree);
}

//_____________________________________________________________________________
void AliRsnReaderTask::Terminate(Option_t */*option*/)
{
//=========================================================
// Terminate analysis
//=========================================================

  AliDebug(1, "Terminating");
}

