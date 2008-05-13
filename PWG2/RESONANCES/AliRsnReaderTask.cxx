/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//----------------------------------------------------------------------------------
//  Class AliRsnReaderTask
// ------------------------
// Reader for conversion of ESD output into the internal format
// used for resonance study.
// ---
// original author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
// ---
// adapted for Analysis Framework
// by    : R. Vernet                          (email: renaud.vernet@cern.ch)
//----------------------------------------------------------------------------------

#include <Riostream.h>

#include "AliLog.h"

#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODInputHandler.h"
#include "AliAODHandler.h"

#include "AliRsnEvent.h"
#include "AliRsnReaderTask.h"

ClassImp(AliRsnReaderTask)

//_____________________________________________________________________________
AliRsnReaderTask::AliRsnReaderTask(ESource source) :
  AliAnalysisTaskSE(),
  fSource(source),
  fReader(0x0),
  fPID(0x0),
  fRsnEvents(0x0)
{
//=========================================================
// Default constructor (not recommended)
//=========================================================
}

//_____________________________________________________________________________
AliRsnReaderTask::AliRsnReaderTask(const char *name, ESource source) : 
  AliAnalysisTaskSE(name),
  fSource(source),
  fReader(0x0),
  fPID(0x0),
  fRsnEvents(0x0)
{
//=========================================================
// Working constructor (recommended)
//=========================================================
}

//_____________________________________________________________________________
AliRsnReaderTask::AliRsnReaderTask(const AliRsnReaderTask& obj) :
  AliAnalysisTaskSE(obj),
  fSource(obj.fSource),
  fReader(obj.fReader),
  fPID(obj.fPID),
  fRsnEvents(0x0)
{
//=========================================================
// Copy constructor (not recommended)
//=========================================================
}

//_____________________________________________________________________________
AliRsnReaderTask& AliRsnReaderTask::operator=(const AliRsnReaderTask& /*obj*/) 
{
//=========================================================
// Assignment operator (not recommended)
//=========================================================

    AliInfo("Not implemented. Avoid using the assignment operator");
	return *this;
}

//_____________________________________________________________________________
void AliRsnReaderTask::UserCreateOutputObjects()
{
//=========================================================
// Create the output container
//=========================================================

    AliDebug(1, "Creating USER output objects");
    fRsnEvents = new TClonesArray("AliRsnEvent", 0);
    fRsnEvents->SetName("AliRsnEvents");
    AddAODBranch("TClonesArray", fRsnEvents);
}

//_____________________________________________________________________________
void AliRsnReaderTask::Init()
{
//=========================================================
// Initialization
//=========================================================

    AliDebug(1, "Initializing");
}

//_____________________________________________________________________________
void AliRsnReaderTask::UserExec(Option_t */*option*/)
{
//=========================================================
// Loops on input container to store data of all tracks.
// Uses the AliRsnReader methods to save them in the output.
//=========================================================

    // static counter
    static Int_t ientry = 0;
    AliInfo(Form("Reading event %d", ++ientry));
    
    // clear previous sample
    fRsnEvents->Clear();
    
    // check for existence of reader, otherwise abort
    if (!fReader) {
        AliError("Event reader not initialized. Impossible to continue");
        return;
    }
    if (!fPID) {
        AliError("PID manager not initialized. Impossible to continue");
        return;
    }
    

	// get MC reference event
	AliMCEventHandler* mcHandler = (AliMCEventHandler*)((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler());
	AliMCEvent *mcEvent = 0x0;
	if (mcHandler) mcEvent = mcHandler->MCEvent();
	
	// work-flow variables
	Bool_t successRead;
 	Int_t  nextIndex = fRsnEvents->GetEntries();
	TClonesArray &array = *fRsnEvents;
    AliRsnEvent *event = new(array[nextIndex]) AliRsnEvent;
    if (!event) {
        AliError("Problem occurred while creating new AliRsnEvent object. Aborting");
        return;
    }
    event->Init();
    AliESDInputHandler *handlerESD;
    AliAODInputHandler *handlerAOD;
    AliESDEvent *esd;
    AliAODEvent *aod;
	
	// read source event according to reader settings
	switch (fSource) {
        case kESD:  // read from ESD event
            handlerESD = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
            if (!handlerESD) {
                AliError("Source set to 'ESD' but ESD handler not initialized");
                return;
            }
            esd = handlerESD->GetEvent();
            if (!esd) {
                AliError("Source set to 'ESD' but ESD event not found in handler");
                return;
            }
            if (!mcEvent) AliWarning("MC info not present");
            successRead = fReader->FillFromESD(event, esd, mcEvent);
            break;
        case kAOD:  // read from AOD event
            handlerAOD = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
            if (!handlerAOD) {
                AliError("Source set to 'AOD' but AOD handler not initialized");
                return;
            }
            aod = handlerAOD->GetEvent();
            if (!aod) {
                AliError("Source set to 'AOD' but AOD event not found in handler");
                return;
            }
            if (!mcEvent) AliWarning("MC info not present");
            successRead = fReader->FillFromAOD(event, aod, mcEvent);
            break;
        case kMC: // read from MC truth only
            if (!mcEvent) {
                AliError("Required usage of MC event for reading. Not possible without MC info. Impossible to continue");
                return;
            }
            successRead = fReader->FillFromMC(event, mcEvent);
            break;
        default: // unrecognized option
            AliError("Source flag not properly set");
            return;
    }
    if (!successRead) {
        array.RemoveAt(nextIndex);
        return;
    }
    
    // if the event reading is successful, perform particle identification
    if (!fPID->Identify(event)) AliWarning(Form("Failed PID for event %d", ientry));
    AliInfo(Form("Event %d: collected %d tracks", ientry, event->GetMultiplicity()));
}

//_____________________________________________________________________________
void AliRsnReaderTask::Terminate(Option_t */*option*/)
{
//=========================================================
// Terminate analysis
//=========================================================

    AliDebug(1, "Terminating");
}
