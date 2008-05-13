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
//                     Class AliRsnPair
//-------------------------------------------------------------------------
// This class computes the invariant mass spectrum of a specified pair of
// particles, throughout a list of AliRsnEvents, and returns it as a TH1D.
// This object is not supposed to be used directly: an AliRsnAnalysis
// should be initialized in a macro and filled with one or more AliRsnPair's
// which are then processed with a given sample of events.
//   
// author: A. Pulvirenti
// email : alberto.pulvirenti@ct.infn.it
//-------------------------------------------------------------------------

#include <Riostream.h>

#include <TH1.h>
#include <TString.h>
#include <TRefArray.h>
#include <TClonesArray.h>

#include "AliLog.h"
#include "AliRsnParticle.h"
#include "AliRsnDaughter.h"
#include "AliRsnDaughterCut.h"
#include "AliRsnDaughterCutPair.h"
#include "AliRsnEvent.h"

#include "AliRsnPair.h"

ClassImp(AliRsnPair)
//--------------------------------------------------------------------------------------------------------
AliRsnPair::AliRsnPair() :
  TNamed(),
  fForMixing(kFALSE),
  fStoreOnlyTrue(kFALSE),
  fTrueMotherPDG(0),
  fPtMin(0.0),
  fPtMax(0.0),
  fVtMax(0.0),
  fCutsPair(0x0),
  fHistogram(0x0)
{
//
// Empty constructor.
// Initializes the data members to default values:
//  - default (empty string) name and title for the object
//  - switch off the 'fStoreOnlyTrue' flag;
//  - assume working with real data (no PDG code of the mother);
//  - no cuts of any kind;
//  - no definition of particles in the pair;
//  - histogram undefined.
// When using this constructor, all analysis elements (particles, histogram)
// must be defined before starting event processing.
//

	Int_t i;
	for (i = 0; i < 2; i++) {
		fMass[i] = 0.0;
		fCharge[i] = '0';
		fType[i] = AliRsnPID::kUnknown;
		fCutsSingle[i] = 0x0;
	}
}
//--------------------------------------------------------------------------------------------------------
AliRsnPair::AliRsnPair
(const char *name, const char *title, Int_t nbins, Double_t min, Double_t max, 
 Double_t ptmin, Double_t ptmax, Double_t dmax) :
  TNamed(name, title),
  fForMixing(kFALSE),
  fStoreOnlyTrue(kFALSE),
  fTrueMotherPDG(0),
  fPtMin(ptmin),
  fPtMax(ptmax),
  fVtMax(dmax),
  fCutsPair(0x0),
  fHistogram(0x0)
{
//
// Constructor with arguments.
// This constructor allows to define some of the initialization values:
//  - name and title of the object
//  - histogram binning and edges
// The other parameters are initialized as in the default constructor.
//

	Int_t i;
	for (i = 0; i < 2; i++) {
		fMass[i] = 0.0;
		fCharge[i] = '0';
		fType[i] = AliRsnPID::kUnknown;
		fCutsSingle[i] = 0x0;
	}
	fHistogram = new TH1D(Form("histo_%s", name), "AliRsnPair::fHistogram", nbins, min, max);
}
//--------------------------------------------------------------------------------------------------------
AliRsnPair::AliRsnPair(const AliRsnPair &copy) :
  TNamed(copy),
  fForMixing(copy.fForMixing),
  fStoreOnlyTrue(copy.fStoreOnlyTrue),
  fTrueMotherPDG(copy.fTrueMotherPDG),
  fPtMin(copy.fPtMin),
  fPtMax(copy.fPtMax),
  fVtMax(copy.fVtMax),
  fCutsPair(0x0),
  fHistogram(0x0)
{
//
// Copy constructor.
// Default behavior as a copy constructor for what concerns non-array data-members.
// The arrays are cloned if they are not NULL.
//

	if (copy.fHistogram) fHistogram = (TH1D*)(copy.fHistogram->Clone());
	if (copy.fCutsPair) fCutsPair = (TObjArray*)(copy.fCutsPair->Clone());
	
	Int_t i;
	for (i = 0; i < 2; i++) {
		fMass[i] = copy.fMass[i];
		fCharge[i] = copy.fCharge[i];
		fType[i] = copy.fType[i];
		if (copy.fCutsSingle[i]) fCutsSingle[i] = (TObjArray*)copy.fCutsSingle[i]->Clone();
	}
}
//--------------------------------------------------------------------------------------------------------
const AliRsnPair& AliRsnPair::operator=(const AliRsnPair &copy)
{
//
// Assignment operator.
// Default behavior like copy constructor.
//

	fHistogram = 0x0;
	fCutsPair = 0x0;
	fCutsSingle[0] = fCutsSingle[1] = 0x0;

    fForMixing = copy.fForMixing;
	fStoreOnlyTrue = copy.fStoreOnlyTrue;
	fTrueMotherPDG = copy.fTrueMotherPDG;
	if (copy.fHistogram) fHistogram = (TH1D*)(copy.fHistogram->Clone());
	if (copy.fCutsPair) fCutsPair = (TObjArray*)(copy.fCutsPair->Clone());
    
    fPtMin = copy.fPtMin;
    fPtMax = copy.fPtMax;
    fVtMax = copy.fVtMax;
	
	Int_t i;
	for (i = 0; i < 2; i++) {
		fMass[i] = copy.fMass[i];
		fCharge[i] = copy.fCharge[i];
		fType[i] = copy.fType[i];
		if (copy.fCutsSingle[i]) fCutsSingle[i] = (TObjArray*)copy.fCutsSingle[i]->Clone();
	}
	
	return (*this);
}
//--------------------------------------------------------------------------------------------------------
void AliRsnPair::Clear(Option_t* /*option*/)
{
//
// Clear arrays and histogram.
// For the sake of security, all pointers are also set explicitly to NULL.
//
	fCutsSingle[0]->Delete();
	fCutsSingle[1]->Delete();
	fCutsPair->Delete();
	
	delete fCutsSingle[0];
	delete fCutsSingle[1];
	delete fCutsPair;
	delete fHistogram;
	
	fCutsSingle[0] = fCutsSingle[1] = fCutsPair = 0;
	fHistogram = 0;
}
//--------------------------------------------------------------------------------------------------------
void AliRsnPair::SetPair
(Char_t charge1, AliRsnPID::EType type1, Char_t charge2, AliRsnPID::EType type2)
{
//
// This method allows to set at once all the parameters of the particles in the pair.
// The mass must not be specified, and it is retrieved from TDatabasePDG,
// using a static method defined in AliRsnDaughter class.
//

	fCharge[0] = charge1;
	fType[0] = type1;
	SetMass(0, AliRsnPID::ParticleMass(type1));

	fCharge[1] = charge2;
	fType[1] = type2;
	SetMass(1, AliRsnPID::ParticleMass(type2));
}
//--------------------------------------------------------------------------------------------------------
void AliRsnPair::AddCutPair(AliRsnDaughterCutPair *cut)
{
//
// Add a pair cut.
// If the cut array is NULL, it is initialized here.
//

	if (!fCutsPair) fCutsPair = new TObjArray(0);
	fCutsPair->AddLast(cut);
}
//--------------------------------------------------------------------------------------------------------
void AliRsnPair::AddCutSingle(Int_t i, AliRsnDaughterCut *cut)
{
//
// Add a single particle cut.
// If the cut array is NULL, it is initialized here.
//

	if (i < 0 || i > 1) return;
	if (!fCutsSingle[i]) fCutsSingle[i] = new TObjArray(0);
	fCutsSingle[i]->AddLast(cut);
}
//--------------------------------------------------------------------------------------------------------
Stat_t AliRsnPair::Process(AliRsnEvent *event1, AliRsnEvent *event2, Bool_t usePID)
{
//
// Scans the two events specified in argument to fill the histogram.
// This method essentially calls the AliRsnPair::Fill() method one or many times.
// When the "noPID" argument is kFALSE, the analysis is done with identified particles
// and this causes the Fill() method to be called only once, for the two lists of
// identified particles of the two kinds specified in AliRsnPair datamembers.
// When the "noPID" argument is kTRUE, the analysis is done with all collections
// of particles of the same sign as specified in the two arguments of the pair.
// ---
// Particles of type #1 are taken in 'event1', and particles of type #2 are taken in 'event2'.
// When doing single-event analysis (for resonance signal or like-sign background),
// the second argument can be simply skipped.
// When doing event mixing, the two arguments must be not null and different.
// If argument #1 is NULL, an error is raised, while if argument #2 is NULL, no error is raised,
// and 'event2' argument is set equal to 'event1' (= single event processing).
// ---
// Return value is the total number of pairs processed.
//

	// preliminary checks
	if (!event1) {
		// argument #1 cannot be NULL
		AliError("Argument #1 cannot be NULL.");
		return 0.0;
	}
	if (!event2) {
		// if argument #2 is NULL, it is put equal to argument #1
		event2 = event1;
	}
    
    // define if same indexes must be summed or not depending if
    // the two events are the same or not
    Bool_t skipSameIndex = (event1 == event2);
    
    TRefArray *list1, *list2;
    if (usePID) {
        // analysis with PID: only two collections are processed
    	list1 = event1->GetTracks(fCharge[0], fType[0]);
    	list2 = event2->GetTracks(fCharge[1], fType[1]);
    }
    else {
        // analysis without PID: directly take the two arrays with all particles of a given sign
        list1 = event1->GetCharged(fCharge[0]);
        list2 = event2->GetCharged(fCharge[1]);
    }
    
    /*
    TRefArray *list[2];
    Short_t i;
    for (i = 0; i < 2; i++) {
        if (usePID) {
            // analysis with PID: only two collections are processed
    	   list[i] = event1->GetTracks(fCharge[i], fType[i]);
        }
        else {
            // analysis without PID: directly take the two arrays with all particles of a given sign
    	   list[i] = event1->GetCharged(fCharge[i]);
        }
    }
    Stat_t nPairs = Fill(list[0], list[1], skipSameIndex);
    */
    
    Stat_t nPairs = Fill(list1, list2, skipSameIndex);
    return nPairs;
}
//--------------------------------------------------------------------------------------------------------
Stat_t AliRsnPair::Fill(TRefArray *list1, TRefArray *list2, Bool_t skipSameIndex)
{
//
// [PRIVATE METHOD]
// This is the core of the pair work flow.
// It loops on all particles contained in the two lists passed as arguments,
// and for each pair it computes the invariant mass and fills its histogram.
// Third argument (skipSameIndex) is a security flag: 
// when the two lists come from the same event, two tracks with the same index 
// must not be summed, and setting to kTRUE this flag this is ensured.
// When it is kFALSE, this check is not done (event mixing).
// ---
// Before starting the operation, it checks that the two arguments are meaningful.
// ---
// Return value is the number of pairs processed.
//

	// define two 'cursor' objects
	AliRsnDaughter *track1 = 0, *track2 = 0;
	
	// preliminary checks
	if (!list1) {
		AliError("List #1 cannot be NULL.");
		return 0.0;
	}
    if (!list2) {
		AliError("List #2 cannot be NULL.");
		return 0.0;
	}
    
    // create histogram if it is not present
	if (!fHistogram) {
		AliError("Histogram undefined.");
		return 0.0;
	}
	
	// define iterators for the two collections
	TRefArrayIter iter1(list1);
	TRefArrayIter iter2(list2);
	
	// define temporary variables for better code readability
	Stat_t nPairs = 0;
	Int_t  pdgRef = TMath::Abs(fTrueMotherPDG);
	
	// loop on particle of type 1 (in 'event1')
	while ( (track1 = (AliRsnDaughter*)iter1.Next()) ) {
		// check against impact parameter cut
        if (fVtMax > 0.0 && track1->Vt() > fVtMax) continue;
        // check against 1-particle cuts (particle type #1)
        if (!SingleCutCheck(0, track1)) continue;
		// loop on particles of type #2 (in 'event2')
		iter2.Reset();
		while ( (track2 = (AliRsnDaughter*)iter2.Next()) ) {
			// skip the case when particle 2 is the same as particle 1
			if (skipSameIndex && (track1->Index() == track2->Index())) continue;
            // check against impact parameter cut
            if (fVtMax > 0.0 && track2->Vt() > fVtMax) continue;
			// check against 1-particle cuts (particle type #2)
			if (!SingleCutCheck(1, track2)) continue;
			// check against 2-particle cuts
			if (!PairCutCheck(track1, track2)) continue;
			// compute total 4-momentum
			track1->SetM(fMass[0]);
			track2->SetM(fMass[1]);
			AliRsnDaughter sum = AliRsnDaughter::Sum(*track1, *track2);
            // reject wrong pairs with mass = 0
            if (sum.M() <= 1E-4) continue;
            // check transverse momentum bin
            if (fPtMin > 0.0 && sum.Pt() < fPtMin) continue;
            if (fPtMax > 0.0 && sum.Pt() > fPtMax) continue;
			// if the choice to get ONLY true pairs is selected, a check is made on the mothers
			 if (fStoreOnlyTrue) {
			     AliRsnParticle *part = sum.GetParticle();
			     if (!part) continue;
			     if (TMath::Abs(part->PDG()) != pdgRef) continue;
            }
			// fill histogram
			fHistogram->Fill(sum.M());
			nPairs++;
		} // end loop 2
	} // end loop 1
	
	return nPairs;
}
//--------------------------------------------------------------------------------------------------------
Bool_t AliRsnPair::SingleCutCheck(Int_t ipart, AliRsnDaughter *track) const
{
//
// [PRIVATE METHOD]
// Checks a track against single particle cuts (if defined)
//

	if (ipart < 0 || ipart > 1) {
		AliError(Form("Required cuts collection for particle index %d [allowed 0 or 1]", ipart));
		return kFALSE;
	}
	if (!fCutsSingle[ipart]) return kTRUE;
    if (fCutsSingle[ipart]->IsEmpty()) return kTRUE;
		
	TObjArrayIter iter(fCutsSingle[ipart]);
	AliRsnDaughterCut *cut = 0;
	while ( (cut = (AliRsnDaughterCut*)iter.Next()) ) {
		if (!cut->Pass(track)) return kFALSE;
	}
	
	return kTRUE;
}
//--------------------------------------------------------------------------------------------------------
Bool_t AliRsnPair::PairCutCheck(AliRsnDaughter *track1, AliRsnDaughter *track2) const
{
//
// [PRIVATE METHOD]
// Checks a pair against pair cuts (if defined)
//
	if (fCutsPair == 0x0) return kTRUE;
    if (fCutsPair->IsEmpty()) return kTRUE;
	
	TObjArrayIter iter(fCutsPair);
	AliRsnDaughterCutPair *cut = 0;
	while ( (cut = (AliRsnDaughterCutPair*)iter.Next()) ) {
		if (!cut->Pass(track1, track2)) return kFALSE;
	}
	
	return kTRUE;
}
