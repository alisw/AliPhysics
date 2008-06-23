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
//  Class AliRsnReaderTaskSE
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
#include "AliRsnReader.h"
#include "AliRsnPID.h"
#include "AliRsnReaderTaskSE.h"

ClassImp(AliRsnReaderTaskSE)

//_____________________________________________________________________________
AliRsnReaderTaskSE::AliRsnReaderTaskSE() :
  AliAnalysisTaskSE(),
  fReader(0x0),
  fPID(0x0),
  fRsnEvent(0x0)
{
//=========================================================
// Default constructor (not recommended)
//=========================================================
}

//_____________________________________________________________________________
AliRsnReaderTaskSE::AliRsnReaderTaskSE(const char *name) :
  AliAnalysisTaskSE(name),
  fReader(0x0),
  fPID(0x0),
  fRsnEvent(0x0)
{
//=========================================================
// Working constructor (recommended)
// Initializes the reader and PID objects.
//=========================================================

    fReader = new AliRsnReader;
    fPID = new AliRsnPID;
}

//_____________________________________________________________________________
AliRsnReaderTaskSE::AliRsnReaderTaskSE(const AliRsnReaderTaskSE& obj) :
  AliAnalysisTaskSE(obj),
  fReader(obj.fReader),
  fPID(obj.fPID),
  fRsnEvent(0x0)
{
//=========================================================
// Copy constructor (not recommended)
//=========================================================
}

//_____________________________________________________________________________
AliRsnReaderTaskSE& AliRsnReaderTaskSE::operator=(const AliRsnReaderTaskSE& /*obj*/)
{
//=========================================================
// Assignment operator (not recommended)
//=========================================================

    AliInfo("Not implemented. Avoid using the assignment operator");
	return *this;
}

//_____________________________________________________________________________
void AliRsnReaderTaskSE::UserCreateOutputObjects()
{
//=========================================================
// Create the output container
//=========================================================

    AliDebug(1, "Creating USER output objects");

    // check for existence of reader, otherwise abort
    if (!fReader) {
        AliFatal("Event reader not initialized. Impossible to continue");
        return;
    }
    if (!fPID) {
        AliFatal("PID manager not initialized. Impossible to continue");
        return;
    }
    else {
        // the PID object is always used in realistic mode here
        // to fill the "realistic" index array in AliRsnEvent
        fPID->SetMethod(AliRsnPID::kRealistic);
    }

    fRsnEvent = new AliRsnEvent();
    fRsnEvent->SetName("rsnEvents");
    fRsnEvent->Init();
    AddAODBranch("AliRsnEvent", &fRsnEvent);
}

//_____________________________________________________________________________
void AliRsnReaderTaskSE::Init()
{
//=========================================================
// Initialization
//=========================================================

    AliDebug(1, "Initializing");
}

//_____________________________________________________________________________
void AliRsnReaderTaskSE::UserExec(Option_t */*option*/)
{
//=========================================================
// Loops on input container to store data of all tracks.
// Uses the AliRsnReader methods to save them in the output.
//=========================================================

    // static counter
    AliInfo(Form("Reading event %d", ++fEntry));

    // clear previous sample
    fRsnEvent->Clear();

    // read event, identify
    if (!fReader->Fill(fRsnEvent, fInputEvent, fMCEvent)) AliWarning("Failed reading");
    if (!fPID->Identify(fRsnEvent)) AliWarning("Failed PID");
    AliInfo(Form("Collected %d tracks", fRsnEvent->GetMultiplicity()));
}

//_____________________________________________________________________________
void AliRsnReaderTaskSE::Terminate(Option_t */*option*/)
{
//=========================================================
// Terminate analysis
//=========================================================

    AliDebug(1, "Terminating");
}
