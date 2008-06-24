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

//-------------------------------------------------------------------------
//                     Class AliRsnAnalysisSimple
//-------------------------------------------------------------------------
// This class is a manager to process one or more pair analysis with a
// given event sample, specified as a TTree passed to this object.
// Each single pair analysis must be defined with its specifications, like
// histogram binning, cuts and everything else.
// This class utiliy consists in defining a unique event sample which is
// used for all pairs analysis, and all kinds of event mixing (if required).
// All histograms computed in a single execution, are then saved into a file.
// This object contains a two TObjArray's:
//  - one to contain single event pair analysis objects (signal, like-sign)
//  - one to contain all event-mixing pair analysis objects
// When a new pair is added, it must be then specified in what container it
// must be placed, in order to avoid meaningless results.
//
// author: A. Pulvirenti
// email : alberto.pulvirenti@ct.infn.it
//-------------------------------------------------------------------------

#include <TH1.h>
#include <TTree.h>
#include <TFile.h>
#include <TArray.h>
#include <TClonesArray.h>

#include "AliLog.h"
#include "AliRsnPID.h"
#include "AliRsnDaughter.h"
#include "AliRsnEvent.h"
#include "AliRsnEventBuffer.h"
#include "AliRsnPairSimple.h"
#include "AliRsnAnalyzerSimple.h"
#include "AliRsnAnalysisSimple.h"

ClassImp(AliRsnAnalysisSimple)

//_____________________________________________________________________________
AliRsnAnalysisSimple::AliRsnAnalysisSimple(AliRsnAnalyzerSimple *ana, AliRsnPID *pid) :
  TObject(),
  fInitialized(kFALSE),
  fStep(1000),
  fTree(0x0),
  fPID(pid),
  fAnalyzer(ana)
{
//
// Constructor
// Initializes all pointers and collections to NULL.
//

    strcpy(fFileName, "default.root");
}


//_____________________________________________________________________________
void AliRsnAnalysisSimple::Clear(Option_t* /*option*/)
{
//
// Clear heap
//
}

//_____________________________________________________________________________
void AliRsnAnalysisSimple::SetEventsTree(TTree *tree)
{
//
// Set the tree containing the events to be processed.
// Counts also the number of events and stores it in a private datamember.
// This can avoid the time-wasting entries count in a long TChain.
//
    fTree = tree;
    AliInfo(Form("Total number of events to be processed: %d", fTree->GetEntries()));
}

//_____________________________________________________________________________
Bool_t AliRsnAnalysisSimple::Initialize()
{
//
// Various initialization processes
//
    // check process objects
    if (!fPID) {
        AliError("PID not initialized");
        return kFALSE;
    }
    if (!fAnalyzer) {
        AliError("Analyzer not initialized");
        return kFALSE;
    }

    // set reference to PID method to all pairs
    AliRsnPairSimple *pair = 0;
    TObjArrayIter pairIterator(fAnalyzer->GetPairs());
    while ( (pair = (AliRsnPairSimple*)pairIterator.Next()) ) {
        pair->SetPIDMethod(fPID->GetMethod());
    }
    if (fAnalyzer->GetMixPairs() && !fAnalyzer->GetMixPairs()->IsEmpty()) {
        TObjArrayIter mixIterator(fAnalyzer->GetMixPairs());
        while ( (pair = (AliRsnPairSimple*)mixIterator.Next()) ) {
            pair->SetPIDMethod(fPID->GetMethod());
        }
    }

    // initialize analyzer
    fAnalyzer->Init();

    // at the end, update flag for initialization
    fInitialized = kTRUE;

    return kTRUE;
}

//_____________________________________________________________________________
Stat_t AliRsnAnalysisSimple::Process()
{
//
// Computes all invariant mass distributions defined in the AliRsnPair objects collected.
// Depending on the kind of background evaluation method, computes also this one.
//

    // check initialization
    if (!fInitialized) {
        AliError("Analysis not initialized. Use method 'Initialize()'");
        return 0.0;
    }

    // set cursor object
    AliRsnEvent *event = 0x0;
    fTree->SetBranchAddress("rsnEvents", &event);

    // output counter
    Stat_t nPairs = 0.0;

	// loop on events
	Int_t i, nEvents = (Int_t)fTree->GetEntries();
    for (i = 0; i < nEvents; i++) {
        // message
        if ((i % fStep) == 0) AliInfo(Form("Processing event %d", i));
        // get entry
        fTree->GetEntry(i);
        if (!event) continue;
        fPID->Identify(event);
        nPairs += fAnalyzer->Process(event);
    }

    return nPairs;
}

//_____________________________________________________________________________
void AliRsnAnalysisSimple::SaveOutput() const
{
//
// Writes histograms in current directory
//
    TFile *file = TFile::Open(fFileName, "RECREATE");
    AliRsnPairSimple *pair = 0;
    TH1D *hist = 0;
    TObjArrayIter pairIterator(fAnalyzer->GetPairs());
    while ( (pair = (AliRsnPairSimple*)pairIterator.Next()) ) {
        hist = pair->GetHistogram();
        if (hist) hist->Write();
        hist = pair->GetHistogramMC();
        if (hist) hist->Write();
    }
    if (fAnalyzer->GetMixPairs() && !fAnalyzer->GetMixPairs()->IsEmpty()) {
    	TObjArrayIter mixIterator(fAnalyzer->GetMixPairs());
        while ( (pair = (AliRsnPairSimple*)mixIterator.Next()) ) {
            hist = pair->GetHistogram();
            if (hist) hist->Write();
            hist = pair->GetHistogramMC();
            if (hist) hist->Write();
        }
	}
	file->Close();
}
