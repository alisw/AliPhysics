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
//                     Class AliRsnManager
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

#include "AliRsnEvent.h"
#include "AliRsnPair.h"
#include "AliRsnManager.h"

ClassImp(AliRsnManager)

//_____________________________________________________________________________
AliRsnManager::AliRsnManager() :
  TObject(),
  fUsePID(kTRUE),
  fStep(100),
  fMixEvents(10),
  fMixMultCut(10),
  fMixVzCut(0.5),
  fQueuePos(-1),
  fPairs(0x0),
  fMixPairs(0x0),
  fBuffer(0x0),
  fOutputList(0x0)

{
//
// Constructor
//
}

//_____________________________________________________________________________
void AliRsnManager::Clear(Option_t* /*option*/)
{
//
// Clear heap
//
	
	if (fPairs) {
	   fPairs->Delete();
	   delete fPairs;
    }
    if (fMixPairs) {
	   fMixPairs->Delete();
	   delete fMixPairs;
    }
	if (fBuffer) delete fBuffer;
	if (fOutputList) {
	   fOutputList->Delete();
	   delete fOutputList;
    }
}

//_____________________________________________________________________________
void AliRsnManager::SetQueueSize(Int_t size)
{
    if (fBuffer) {
        fBuffer->Delete();
        delete fBuffer;
    }
    fBuffer = new TObjArray(size);
    fBuffer->SetOwner();
    fQueuePos = -1;
}

//_____________________________________________________________________________
void AliRsnManager::AddPair(AliRsnPair *pair)
{
//
// Add a pair of particle types to be analyzed.
// Second argument tells if the added object is for event mixing.
//

    Bool_t mixing = pair->IsForMixing();
	TObjArray* &target = mixing ? fMixPairs : fPairs;
	if (!target) target = new TObjArray(0);
	target->AddLast(pair);
	
	if (!fOutputList) fOutputList = new TList;
	
	fOutputList->Add(pair->GetHistogram());
}

//_____________________________________________________________________________
Stat_t AliRsnManager::Process(AliRsnEvent *event)
{
//
// For each AliRsnPair definition, fills its histogram 
// taking particles from the passed event.
//
    Stat_t nPairs = 0;
    if (!event) return nPairs;
    
    // check pair list
    if (!fPairs) {
    	AliError("Uninitialized array");
        return 0.0;
    }
    if (fPairs->IsEmpty()) {
    	AliError("No pairs defined");
        return 0.0;
    }
    if (!fBuffer) {
        AliError("Buffer uninitialized");
        return 0.0;
    }

	// create an iterator which run over all pairs
	TObjArrayIter pairIterator(fPairs);
	AliRsnPair *pair = 0;
    while ( (pair = (AliRsnPair*)pairIterator.Next()) ) {
        nPairs += pair->Process(event, event);
    }
    
    // if there are mix pairs, do event mixing
    if (fMixPairs) Mix(event);
    
    // add this event to the queue
    fQueuePos++;
    if (fQueuePos >= fBuffer->GetSize()) fQueuePos = 0;
    AliRsnEvent *prev = (AliRsnEvent*) fBuffer->At(fQueuePos);
    if (prev) fBuffer->RemoveAt(fQueuePos);
    fBuffer->AddAt(event, fQueuePos);
    
    return nPairs;
}


//_____________________________________________________________________________
Stat_t AliRsnManager::Mix(AliRsnEvent *event)
{
//
// Performs event mixing.
// It takes the array of fMixPairDefs and stores results in fMixHistograms.
// First parameter defines how many events must be mixed with each one.
// Events to be mixed are chosen using the other parameters:
//
// - multDiffMax defines the maximum allowed difference in particle multiplicity,
//   a) when last argument is kFALSE (default) the multiplicity comparison is made with respect
//      to the particles of 'type 2' in the pair
//   b) when last argument is kTRUE, the multiplicity comparison is made with respect of total
//      particle multiplicity in the events
//
// - vzDiffMax defines maximum allowed difference in Z coordinate of prim. vertex.
//
// If one wants to exchange the particle types, one has to add both combinations of particles
// as different pair defs.
//
    Stat_t nPairs = 0;

	if (!fBuffer) {
        AliError("Uninitialized queue");
        return 0.0;
    }
	
	// iterator for pairs
	TObjArrayIter pairIterator(fMixPairs);
	AliRsnPair *pair = 0;
	
	// iterator for events
	TObjArrayIter eventIterator(fBuffer);
	AliRsnEvent *mixed = 0;
	    
    // loop on array and mix with passed event
    while ( (mixed = (AliRsnEvent*)eventIterator.Next()) ) {
        if (!CanBeMixed(event, mixed)) continue;
        pairIterator.Reset();
        while ( (pair = (AliRsnPair*)pairIterator.Next()) ) {
            nPairs += pair->Process(event, mixed);
		}
	}
	
	return nPairs;
}

//_____________________________________________________________________________
Bool_t AliRsnManager::CanBeMixed(AliRsnEvent *ev1, AliRsnEvent *ev2)
{
//
// Checks if two events are within the mixing cuts
//
    
    Int_t diffMult = TMath::Abs(ev1->GetMultiplicity() - ev2->GetMultiplicity());
    Bool_t diffVz = TMath::Abs(ev1->GetPrimaryVertexZ() - ev2->GetPrimaryVertexZ());
    if (diffMult <= fMixMultCut && diffVz <= fMixVzCut) {
        return kTRUE;
    }
    else {
        return kFALSE;
    }
}
