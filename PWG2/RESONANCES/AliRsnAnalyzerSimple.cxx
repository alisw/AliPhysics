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

//=========================================================================
// Class AliRsnAnalyzerSimple
//
// Implementation of the event processing which returns histograms of
// invariant mass for resonances and backgrounds evaluation.
//
// author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
//=========================================================================

#include <TTree.h>
#include <TFile.h>
#include <TArray.h>
#include <TClonesArray.h>

#include "AliLog.h"
#include "AliRsnDaughter.h"
#include "AliRsnEvent.h"
#include "AliRsnEventBuffer.h"
#include "AliRsnPairSimple.h"
#include "AliRsnAnalyzerSimple.h"

ClassImp(AliRsnAnalyzerSimple)

//_____________________________________________________________________________
AliRsnAnalyzerSimple::AliRsnAnalyzerSimple(Int_t bufferSize) :
  TObject(),
  fBufferSize(bufferSize),
  fMixMultCut(10),
  fMixVzCut(0.5),
  fNMix(10),
  fPairs(0x0),
  fMixPairs(0x0),
  fBuffer(0x0)
{
//
// Constructor
// Initializes all pointers and collections to NULL.
//
}


//_____________________________________________________________________________
void AliRsnAnalyzerSimple::Clear(Option_t *option)
{
//
// Clear heap
//
    fPairs->Clear(option);
    fMixPairs->Clear(option);
    delete fBuffer;
    fBuffer = 0;
}

//_____________________________________________________________________________
void AliRsnAnalyzerSimple::AddPair(AliRsnPairSimple *pair)
{
//
// Add a pair of particle types to be analyzed.
// Second argument tells if the added object is for event mixing.
//
    Bool_t mixing = pair->IsForMixing();
    TObjArray* &target = mixing ? fMixPairs : fPairs;
    if (!target) target = new TObjArray(0);
    target->AddLast(pair);
}

//_____________________________________________________________________________
void AliRsnAnalyzerSimple::Init()
{
//
// Initialize what needs to.
//

    // buffer
    fBuffer = new AliRsnEventBuffer(fBufferSize, kFALSE);

    // histograms
    AliRsnPairSimple *pair = 0;
    if (fPairs) {
        TObjArrayIter pairIterator(fPairs);
        while ( (pair = (AliRsnPairSimple*)pairIterator.Next()) ) pair->InitHistograms();
    }
    if (fMixPairs) {
        TObjArrayIter mixIterator(fMixPairs);
        while ( (pair = (AliRsnPairSimple*)mixIterator.Next()) ) pair->InitHistograms();
    }
}

//_____________________________________________________________________________
Stat_t AliRsnAnalyzerSimple::Process(AliRsnEvent *event)
{
//
// Computes all invariant mass distributions defined in the AliRsnPair objects collected.
// Depending on the kind of background evaluation method, computes also this one.
//

    // skip empty events
    if (event->GetMultiplicity() < 1) return 0.0;

    // create buffer if not already present
    if (!fBuffer) Init();

    // initialize output values and utility variables
    Stat_t nPairs = 0;
    AliRsnPairSimple *pair = 0x0;

    // break here if NULL argument is passed
    if (!event) {
        AliError("NULL event passed");
        return 0.0;
    }

    // loop over the collection of pair defs
    if (fPairs) {
        TObjArrayIter pairIterator(fPairs);
        while ( (pair = (AliRsnPairSimple*)pairIterator.Next()) ) {
            nPairs += pair->Process(event, event);
        }
    }
    if (fMixPairs) {
        // variables for mixing
        Int_t    mult1, mult2, i, j, count = 0;
        Double_t vz1, vz2;
        // add this event to buffer
        fBuffer->AddEvent(event);
        // event mixing
        vz1 = event->GetPrimaryVertexZ();
        mult1 = event->GetMultiplicity();
        if (!fBuffer->NEmptySlots()) {
            i = fBuffer->GetEventsBufferIndex();
            j = i+1;
            for (;;j++) {
                if (count >= fNMix) break;
                if (j == fBufferSize) j = 0;
                if (j == i) break;
                AliRsnEvent *evmix = fBuffer->GetEvent(j);
                if (!evmix) continue;
                vz2 = evmix->GetPrimaryVertexZ();
                mult2 = evmix->GetMultiplicity();
                if (TMath::Abs(vz1 - vz2) <= fMixVzCut && TMath::Abs(mult1- mult2) <= fMixMultCut) {
                    // loop over the collection of pair defs
                    TObjArrayIter mixIterator(fMixPairs);
                    while ( (pair = (AliRsnPairSimple*)mixIterator.Next()) ) {
                        nPairs += pair->Process(event, evmix);
                        nPairs += pair->Process(evmix, event);
                    }
                    count++;
                }
            }
        }
    }

    return nPairs;
}
