//
// Class AliRsnReaderTaskSE
//
// An AnalysisTask object to convert any kind of source event type (ESD/AOD/MC)
// into the RSN internal format (AliRsnEvent).
// The output of this task is a TTree with converted events, which is saved in a file
// and can then be processed as many times as desired, to build invariant mass spectra.
// ---
// original author: A. Pulvirenti (alberto.pulvirenti@ct.infn.it)
// adapted for Analysis Framework by: R. Vernet (renaud.vernet@cern.ch)
//

#include "AliLog.h"

#include "AliAnalysisManager.h"

#include "AliMCEvent.h"
#include "AliMCEventHandler.h"

#include "AliESDInputHandler.h"
#include "AliAODInputHandler.h"
#include "AliAODHandler.h"

#include "AliRsnPID.h"
#include "AliRsnEvent.h"
#include "AliRsnReader.h"
#include "AliRsnReaderTaskSE.h"

ClassImp(AliRsnReaderTaskSE)

//_____________________________________________________________________________
AliRsnReaderTaskSE::AliRsnReaderTaskSE() :
    AliRsnAnalysisTaskSEBase(),
    fRsnEvent(0x0)
{
//
// Default constructor (not recommended)
//
}

//_____________________________________________________________________________
AliRsnReaderTaskSE::AliRsnReaderTaskSE(const char *name) :
    AliRsnAnalysisTaskSEBase(name),
    fRsnEvent(0x0)
{
//
// Working constructor (recommended)
//
}

//_____________________________________________________________________________
void AliRsnReaderTaskSE::UserCreateOutputObjects()
{
//
// Instantiates the output object (AliRsnEvent) and adds a branch
// to the non-standard AOD output TTree to include it.
// Checks that the necessary data member objects for
// conversion and PID are allocated. If this is not the case,
// raises a fatal error which breaks the AliRoot session.
//

  fRsnEvent = new AliRsnEvent();
  fRsnEvent->SetName("rsnEvents");
  fRsnEvent->Init();
  AddAODBranch("AliRsnEvent", &fRsnEvent);
}

//_____________________________________________________________________________
void AliRsnReaderTaskSE::Init()
{
//
// Inherited function.
// Here it does not need to do anything, so it is left dummy.
//
}

//_____________________________________________________________________________
void AliRsnReaderTaskSE::UserExec(Option_t */*option*/)
{
//
// Execution core of the class.
// Uses the AliRsnReader and AliRsnPID methods to convert input data
// and store them in the output AOD event, with all required computations.
//

  AliDebug(1,Form("Reading event %d", ++fEntry));

  // before adding new data, the ones from previous event
  // must be cleared explicitly
  fRsnEvent->Clear();


  // step 1: conversion
  Bool_t ok = kFALSE;
  switch (fInputType[0])
  {
    case kAOD:
      AliDebug(5, "Reading AOD event...");
      ok = fReader.FillFromAOD(fRsnEvent, (AliAODEvent*)fInputEvent, fMCEvent);
      AliDebug(5, "...done");
      break;
    case kESD:
      AliDebug(5, "Reading ESD event...");
      ok = fReader.FillFromESD(fRsnEvent, (AliESDEvent*)fInputEvent, fMCEvent);
      AliDebug(5, "...done");
      break;
    case kESDMC:
      AliDebug(5, "Reading ESD event with MC...");
      ok = fReader.FillFromESD(fRsnEvent, (AliESDEvent*)fInputEvent, fMCEvent);
      AliDebug(5, "...done");
      break;
    case kMC:
      AliDebug(5, "Reading MC only event...");
      ok = fReader.FillFromMC(fRsnEvent, fMCEvent);
      AliDebug(5, "...done");
      break;
    default:
      AliError("Type not supported ...");
      return;
  }
  if (!ok) AliWarning("Failed reading");

  // step 2: PID probability computation
  //if (!fPID.Process(fRsnEvent)) AliWarning("Failed PID");

  AliDebug(1,Form("Collected %d tracks", fRsnEvent->GetMultiplicity()));
}

//_____________________________________________________________________________
void AliRsnReaderTaskSE::Terminate(Option_t */*option*/)
{
//
// Inherited function.
// Here it does not need to do anything, so it is left dummy.
//
}
