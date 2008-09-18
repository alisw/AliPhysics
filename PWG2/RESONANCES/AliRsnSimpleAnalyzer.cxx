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
// Class AliRsnSimpleAnalyzer
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
#include <TDirectory.h>

#include "AliLog.h"

#include "AliRsnDaughter.h"
#include "AliRsnEvent.h"
#include "AliRsnEventBuffer.h"
#include "AliRsnSimpleFunction.h"

#include "AliRsnSimpleAnalyzer.h"

ClassImp(AliRsnSimpleAnalyzer)

//_____________________________________________________________________________
AliRsnSimpleAnalyzer::AliRsnSimpleAnalyzer(Int_t bufferSize) :
    TNamed("RsnSimpleAnalyzer", ""),
    fBufferSize(bufferSize),
    fMixMultCut(10),
    fMixVzCut(0.5),
    fNMix(10),
    fSingle(0x0),
    fMix(0x0),
    fBuffer(0x0)
{
//
// Constructor
// Initializes all pointers and collections to NULL.
// Adds this object to the global Directory.
//
  gDirectory->Append(this, kTRUE);
}


//_____________________________________________________________________________
void AliRsnSimpleAnalyzer::Clear(Option_t *option)
{
//
// Clear heap
//

  fSingle->Clear(option);
  fMix->Clear(option);
  delete fBuffer;
  fBuffer = 0;
}

//_____________________________________________________________________________
void AliRsnSimpleAnalyzer::Add(AliRsnSimpleFunction *fcn)
{
//
// Add a pair of particle types to be analyzed.
// Second argument tells if the added object is for event mixing.
//

  Bool_t mixing = fcn->MixFlag();
  TObjArray* &target = mixing ? fMix : fSingle;
  if (!target) target = new TObjArray(0);
  target->AddLast(fcn);
}

//_____________________________________________________________________________
void AliRsnSimpleAnalyzer::Init()
{
//
// Initialize what needs to.
//

  // buffer
  fBuffer = new AliRsnEventBuffer(fBufferSize, kFALSE);

  // histograms
  AliRsnSimpleFunction *fcn = 0;
  if (fSingle)
  {
    TObjArrayIter iter(fSingle);
    while ((fcn = (AliRsnSimpleFunction*)iter.Next()))
    {
      fcn->Init();
    }
  }
  if (fMix)
  {
    TObjArrayIter iter(fMix);
    while ((fcn = (AliRsnSimpleFunction*)iter.Next()))
    {
      fcn->Init();
    }
  }
}

//_____________________________________________________________________________
void AliRsnSimpleAnalyzer::Process(AliRsnEvent *event)
{
//
// Computes all invariant mass distributions defined in the AliRsnPair objects collected.
// Depending on the kind of background evaluation method, computes also this one.
//

  // loop over the collection of pair defs
  ProcessEvents(fSingle, event, 0x0);

  if (fMix)
  {
    // variables for mixing
    Int_t    mult1, mult2, i, j, count = 0;
    Double_t vz1, vz2;
    // add this event to buffer
    fBuffer->AddEvent(event);
    // event mixing
    vz1 = event->GetPrimaryVertexZ();
    mult1 = event->GetMultiplicity();
    if (!fBuffer->NEmptySlots())
    {
      i = fBuffer->GetEventsBufferIndex();
      j = i+1;
      for (;;j++)
      {
        if (count >= fNMix) break;
        if (j == fBufferSize) j = 0;
        if (j == i) break;
        AliRsnEvent *evmix = fBuffer->GetEvent(j);
        if (!evmix) continue;
        vz2 = evmix->GetPrimaryVertexZ();
        mult2 = evmix->GetMultiplicity();
        if (TMath::Abs(vz1 - vz2) <= fMixVzCut && TMath::Abs(mult1- mult2) <= fMixMultCut)
        {
          // loop over the collection of pair defs
          ProcessEvents(fMix, event, evmix);
          ProcessEvents(fMix, evmix, event);
          count++;
        }
      }
    }
  }
}

//_____________________________________________________________________________
void AliRsnSimpleAnalyzer::ProcessEvents
(TObjArray *array, AliRsnEvent *event1, AliRsnEvent *event2)
{
//
// Takes all AliRsnSimpleFunction objects in the passed array
// and processes the given pair of AliRsnEvents with all of them
// If the array is NULL nothing is done
//

  if (!array) return;

  AliRsnSimpleFunction *fcn = 0;
  TObjArrayIter iter(array);
  while ((fcn = (AliRsnSimpleFunction*)iter.Next()))
  {
    if (event2) fcn->ProcessTwo(event1, event2);
    else fcn->ProcessOne(event1);
  }
}
