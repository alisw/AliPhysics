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
//  Class AliRsnAnalysisSimpleTask
// ------------------------
// Reader for conversion of ESD output into the internal format
// used for resonance study.
// ---
// original author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
// ---
// adapted for Analysis Framework
// by    : R. Vernet                          (email: renaud.vernet@cern.ch)
//----------------------------------------------------------------------------------

#include <TChain.h>

#include "AliLog.h"
#include "AliVEvent.h"
#include "AliMCEventHandler.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"

#include "AliMCEvent.h"

#include "AliRsnEvent.h"
#include "AliRsnReader.h"
#include "AliRsnPID.h"
#include "AliRsnPairSimple.h"
#include "AliRsnAnalyzerSimple.h"
#include "AliRsnAnalysisSimpleTask.h"

ClassImp(AliRsnAnalysisSimpleTask)

//_____________________________________________________________________________
AliRsnAnalysisSimpleTask::AliRsnAnalysisSimpleTask() :
  AliAnalysisTask(),
  fEvent(0x0),
  fMC(0x0),
  fReader(0x0),
  fPID(0x0),
  fAnalyzer(0x0),
  fCurrEvent(0x0),
  fHistograms(0x0)
{
//
// Default constructor (not recommended)
//
}

//_____________________________________________________________________________
AliRsnAnalysisSimpleTask::AliRsnAnalysisSimpleTask(const char *name) :
  AliAnalysisTask(name, ""),
  fEvent(0x0),
  fMC(0x0),
  fReader(0x0),
  fPID(0x0),
  fAnalyzer(0x0),
  fCurrEvent(0x0),
  fHistograms(0x0)
{
//
// Working constructor (recommended)
//

    DefineInput (0, TChain::Class());
    DefineOutput (0, TList::Class());
}

//_____________________________________________________________________________
AliRsnAnalysisSimpleTask::AliRsnAnalysisSimpleTask(const AliRsnAnalysisSimpleTask& obj) :
  AliAnalysisTask(obj),
  fEvent(0x0),
  fMC(0x0),
  fReader(obj.fReader),
  fPID(obj.fPID),
  fAnalyzer(obj.fAnalyzer),
  fCurrEvent(0x0),
  fHistograms(0x0)
{
//
// Copy constructor (not recommended)
//
}

//_____________________________________________________________________________
AliRsnAnalysisSimpleTask& AliRsnAnalysisSimpleTask::operator=(const AliRsnAnalysisSimpleTask& /*obj*/)
{
//
// Assignment operator (not recommended)
//

    AliInfo("Not implemented. Avoid using the assignment operator");
    return *this;
}

//_____________________________________________________________________________
void AliRsnAnalysisSimpleTask::CreateOutputObjects()
{
//
// Create the output container
//

    // check for presence of NECESSARY data-members
    if (!fReader) {
        AliFatal("Event reader not initialized. Impossible to continue. Aborting with fatal error.");
        return;
    }
    if (!fPID) {
        AliFatal("PID manager not initialized. Impossible to continue. Aborting with fatal error.");
        return;
    }
    if (!fAnalyzer) {
        AliFatal("Analysis manager not initialized. Impossible to continue. Aborting with fatal error.");
        return;
    }

    // OpenFile (0);

    // output histogram list
    fHistograms = new TList(0);

    // initialize analyzer
    fAnalyzer->Init();

    // store all histograms in the pairs into the list
    TObjArray *array = fAnalyzer->GetPairs();
    AliRsnPairSimple *pair;
    if (array) {
        TObjArrayIter iter(array);
        while ( (pair = (AliRsnPairSimple*)iter.Next()) ) {
            fHistograms->AddLast((TObject*)pair->GetHistogram());
            fHistograms->AddLast((TObject*)pair->GetHistogramMC());
        }
    }
    array = fAnalyzer->GetMixPairs();
    if (array) {
        TObjArrayIter iter(array);
        while ( (pair = (AliRsnPairSimple*)iter.Next()) ) {
            fHistograms->AddLast((TObject*)pair->GetHistogram());
            fHistograms->AddLast((TObject*)pair->GetHistogramMC());
        }
    }
}

//_____________________________________________________________________________
void AliRsnAnalysisSimpleTask::ConnectInputData(Option_t *)
{
//
// Connect the input data
//

    // connect ESD
    AliInputEventHandler *inH = dynamic_cast<AliInputEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    if (!inH) {
		AliError("Could not get InputHandler");
	}
	else {
		fEvent = inH->GetEvent();
	}

	// connect MC
    AliMCEventHandler* mcH = dynamic_cast<AliMCEventHandler*>((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler());
    if (!mcH) {
		AliError("Could not get MCInputHandler");
	}
    else {
		fMC = mcH->MCEvent();
	}
}

//_____________________________________________________________________________
void AliRsnAnalysisSimpleTask::Exec(Option_t */*option*/)
{
//
// Loops on input container to store data of all tracks.
// Uses the AliRsnReader methods to save them in the output.
//

    // check ESD and MC events
    if (!fEvent) {
        AliWarning("Event not available here");
        return;
    }
    if (!fMC) {
        AliWarning("MC event not available here");
        return;
    }

    // clear previous event
    if (!fCurrEvent) {
        fCurrEvent = new AliRsnEvent;
        fCurrEvent->Init();
    }
    else fCurrEvent->Clear();

    // read event, identify
    if (!fReader->Fill(fCurrEvent, fEvent, fMC)) AliWarning("Failed reading");
    if (!fPID->Identify(fCurrEvent)) AliWarning("Failed PID");
    AliInfo(Form("Collected %d tracks", fCurrEvent->GetMultiplicity()));

    // process event with analyzer
    fPID->Identify(fCurrEvent);
    fAnalyzer->Process(fCurrEvent);

    // post collected data
    PostData (0, fHistograms);
}

//_____________________________________________________________________________
void AliRsnAnalysisSimpleTask::Terminate(Option_t */*option*/)
{
//
// Terminate analysis
//

    AliDebug(1, "Terminating");
    AliAnalysisTask::Terminate();
}
