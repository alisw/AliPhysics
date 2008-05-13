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
//                     Class AliRsnAnalysis
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
#include "AliRsnDaughterCut.h"
#include "AliRsnDaughterCutPair.h"
#include "AliRsnEvent.h"
#include "AliRsnPair.h"
#include "AliRsnAnalysis.h"

ClassImp(AliRsnAnalysis)

//_____________________________________________________________________________
AliRsnAnalysis::AliRsnAnalysis(const char *branchName) :
  TObject(),
  fSkipUnbalanced(kFALSE),
  fStep(1000),
  fBranchName(branchName),
  fMixMultCut(10),
  fMixVzCut(0.5),
  fNEvents(0),
  fPID(0x0),
  fPairs(0x0),
  fMixPairs(0x0),
  fTree(0x0),
  fMatches(0x0),
  fEvents(0x0)
{
//
// Constructor
// Initializes all pointers and collections to NULL.
//
}


//_____________________________________________________________________________
void AliRsnAnalysis::Clear(Option_t *option)
{
//
// Clear heap
//
	fPairs->Clear(option);
    delete [] fMatches;
    fMatches = 0;
    delete fPID;
    fPID = 0;
}

//_____________________________________________________________________________
void AliRsnAnalysis::SetEventsTree(TTree *tree)
{
//
// Set the tree containing the events to be processed.
// Counts also the number of events and stores it in a private datamember.
// This can avoid the time-wasting entries count in a long TChain.
//
	fTree = tree;
    fNEvents = tree->GetEntries();
    AliInfo(Form("Total number of events to be processed: %d", fNEvents));
    
    // link branch
    if (fEvents) delete fEvents;
    fEvents = new TClonesArray("AliRsnEvent", 0);
    fTree->SetBranchAddress(fBranchName.Data(), &fEvents);
}


//_____________________________________________________________________________
void AliRsnAnalysis::AddPair(AliRsnPair *pair)
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
Stat_t AliRsnAnalysis::Process()
{
//
// Computes all invariant mass distributions defined in the AliRsnPair objects collected.
// Depending on the kind of background evaluation method, computes also this one.
//
    AliInfo("Processing");

	// create an iterator which run over all pairs
	TObjArrayIter pairIterator(fPairs);
	AliRsnPair *pair = 0;
		
    Bool_t usePID;
	Stat_t nPairs = 0;
    Int_t  i;
	
	// loop on events
    for (i = 0; i < fNEvents; i++) {
        // get entry
        AliRsnEvent *event = Evt(i);
        if (!event) continue;
        // if requested, skip unbalanced events
        if (fSkipUnbalanced && (event->GetNPos() == 0 || event->GetNNeg() == 0)) continue;
        // adjust PID
        usePID = AdjustPID(event);
        // loop over the collection of pair defs
        pairIterator.Reset();
        while ( (pair = (AliRsnPair*)pairIterator.Next()) ) {
            nPairs += pair->Process(event, event, usePID);
        }
        // message
        if ((i+1) % fStep == 0) AliInfo(Form("%d events processed: %d pairs collected", i+1, (Int_t)nPairs));
    }
	
	return nPairs;
}


//_____________________________________________________________________________
Stat_t AliRsnAnalysis::EventMixing(Int_t nEventsToMix)
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

    AliInfo("Processing");
	// create an iterator which run over all pairs
    if (!fMixPairs) {
    	AliError("Uninitialized MIX array");
        return 0.0;
    }
    if (fMixPairs->IsEmpty()) {
    	AliError("No mixing pairs defined");
        return 0.0;
    }
	TObjArrayIter pairIterator(fMixPairs);
	AliRsnPair *pair = 0;
    
    // find matches for all events
    if (fMatches) delete [] fMatches;
    FindMatches(nEventsToMix);
    		
	// define variables to store data about particles
	Bool_t usePID;
	Stat_t nPairs = 0;
	Int_t iev, imatch, imatched, nEvents = fNEvents;
	for (iev = 0; iev < nEvents; iev++) {
		// get event
		AliRsnEvent *evtry = Evt(iev);
		if (!evtry) continue;
		AdjustPID(evtry);
		AliRsnEvent *event1 = new AliRsnEvent(*evtry);
        // loop on matches
        for (imatch = 0; imatch < nEventsToMix; imatch++) {
        	imatched = (fMatches[iev])[imatch];
            if (imatched < 0 || imatched > nEvents) continue;
            // get other event
            AliRsnEvent *event2 = Evt(imatched);
            if (!event2) continue;
            usePID = AdjustPID(event2);
            if (fPID) event2->FillPIDArrays();
            // loop on pair definitions
            pairIterator.Reset();
            while ( (pair = (AliRsnPair*)pairIterator.Next()) ) {
            	// processing is done twice, exchanging the event arguments
                // mix particles of type #1 from event #1 with particles of type #2 in event #2
                nPairs += pair->Process(event1, event2, usePID);
                // mix particles of type #1 from event #2 with particles of type #2 in event #1
                nPairs += pair->Process(event2, event1, usePID);
			}
		}
		delete event1;
		event1 = 0;
		// message
        if ((iev+1) % fStep == 0) AliInfo(Form("%d events processed: %d pairs collected", iev+1, (Int_t)nPairs));
	} // end loop on events
	
	return nPairs;
}


//_____________________________________________________________________________
void AliRsnAnalysis::SaveOutput(const char *fileName, const char *fileOpt) const
{
//
// Writes histograms in current directory
//
	TFile *file = TFile::Open(fileName, fileOpt);
	AliRsnPair *pair = 0;
    TObjArrayIter pairIterator(fPairs);
    while ( (pair = (AliRsnPair*)pairIterator.Next()) ) pair->GetHistogram()->Write();
    if (fMixPairs && !fMixPairs->IsEmpty()) {
    	TObjArrayIter mixIterator(fMixPairs);
        while ( (pair = (AliRsnPair*)mixIterator.Next()) ) pair->GetHistogram()->Write();
	}       
	file->Close();
}

//_____________________________________________________________________________
AliRsnEvent* AliRsnAnalysis::Evt(Int_t i)
{
//
// Return the event stored in the TClonesArray corresponding
// to the entry number passed as argument
//

    if (!fEvents) return 0x0;
    if (i < 0 || i >= fNEvents) return 0x0;
    
    fTree->GetEntry(i);
    TObjArrayIter iter(fEvents);
    AliRsnEvent *evt = (AliRsnEvent*)iter.Next();
    return evt;
}

//_____________________________________________________________________________
void AliRsnAnalysis::FindMatches(Int_t nEventsToMatch)
{
//
// Initializes the "fMatch" array and computes, for each event, the good matches 
// for event mixing computations (if requested).
//
	Int_t    *mult    = new Int_t[fNEvents];
    Double_t *vz      = new Double_t[fNEvents];
    
    // loop on events
    Int_t iev;
	for (iev = 0; iev < fNEvents; iev++) {
		AliRsnEvent *event = Evt(iev);
		if (!event) continue;
        if (fSkipUnbalanced && (event->GetNPos() == 0 || event->GetNNeg() == 0)) continue;
        vz[iev] = event->GetPrimaryVertexZ();
        mult[iev] = event->GetMultiplicity();
	}
    
    // initialize arrays of matches and compute
    Int_t imatch, nmatch, diffMult;
    Double_t diffVz;
    fMatches = new TArrayI[fNEvents];
    for (iev = 0; iev < fNEvents; iev++) {
    	fMatches[iev].Set(nEventsToMatch);
        imatch = iev;
        nmatch = 0;
        while (nmatch < nEventsToMatch) {
        	imatch++;
            if (imatch >= fNEvents) imatch -= fNEvents;
            if (imatch == iev) break;
            diffMult = TMath::Abs(mult[iev] - mult[imatch]);
            diffVz = TMath::Abs(vz[iev] - vz[imatch]);
            if (diffMult <= fMixMultCut && diffVz <= fMixVzCut) {
                (fMatches[iev])[nmatch] = imatch;
                nmatch++;
            }
        }
    }
    
    delete [] mult;
    delete [] vz;
}


//_____________________________________________________________________________
Bool_t AliRsnAnalysis::AdjustPID(AliRsnEvent* &event)
{
//
// if a PID different from "native" has been chosen, particles are identified again
// the AliRsnPair::Process method wants an argument (3rd) which is TRUE when PID must
// be taken into account, and FALSE otherwise.
// Returns de corresponding value of this flag for that method, according to the choice done.
//
//
    
    // if the fPID object is NULL, the native PID is used
    if (!fPID) {
        return kTRUE;
    }
    else {
        fPID->Identify(event);
        AliRsnPID::EMethod method = fPID->GetMethod();
        if (method == AliRsnPID::kNone) return kFALSE;
        else return kTRUE;
    }
}
