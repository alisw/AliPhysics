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
//            Reconstruction and analysis of a binary resonance
// ........................................
// ........................................
// ........................................
// ........................................
// ........................................
// 
// author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
//-------------------------------------------------------------------------

#include <Riostream.h>

#include <TH1.h>
#include <TTree.h>
#include <TDatabasePDG.h>

#include "AliRsnDaughter.h"
#include "AliRsnDaughterCut.h"
#include "AliRsnEvent.h"
#include "AliRsnAnalysis.h"

ClassImp(AliRsnAnalysis)

//--------------------------------------------------------------------------------------------------------
AliRsnAnalysis::AliRsnAnalysis() 
{
//
// Constructor
// Initializes all pointers and collections to NULL.
//
	fMixHistograms = 0;
	fHistograms = 0;
	fEventsTree = 0;
	fMixPairDefs = 0;
	fPairDefs = 0;
	fPairCuts = 0;
	
	Int_t i;
	for (i = 0; i < AliPID::kSPECIES; i++) fCuts[i] = 0;
}
//--------------------------------------------------------------------------------------------------------
void AliRsnAnalysis::AddCutPair(AliRsnDaughterCut *cut)
{
//
// Add a cut on pairs.
// This cut is global for all pairs.
//
	if (!cut->IsPairCut()) {
		Warning("AddCutPair", "This is a single cut, cannot be added");
		return;
	}
	
	if (!fPairCuts) fPairCuts = new TObjArray(0);

	fPairCuts->AddLast(cut);
}
//--------------------------------------------------------------------------------------------------------
void AliRsnAnalysis::AddCutSingle(AliPID::EParticleType type, AliRsnDaughterCut *cut)
{
//
// Add a cut on single particles.
// This cut must be specified for each particle type.
//
	if (cut->IsPairCut()) {
		Warning("AddCutSingle", "This is a pair cut, cannot be added");
		return;
	}
	
	if (type >= AliPID::kElectron && type <= AliPID::kProton) {
		if (!fCuts[type]) fCuts[type] = new TObjArray(0);
		fCuts[type]->AddLast(cut);
	}
}
//--------------------------------------------------------------------------------------------------------
void AliRsnAnalysis::AddMixPairDef
(AliPID::EParticleType p1, Char_t sign1, AliPID::EParticleType p2, Char_t sign2)
{
//
// Adds a new pair definition to create a new histogram,
// for the event mixing step.
// If the pair defs array is NULL, it is initialized here.
//
	if (!fMixPairDefs) fMixPairDefs = new TObjArray(0);
	fMixPairDefs->AddLast( new AliPairDef(p1, sign1, p2, sign2, 0, kFALSE) );
}
//--------------------------------------------------------------------------------------------------------
void AliRsnAnalysis::AddPairDef
(AliPID::EParticleType p1, Char_t sign1, AliPID::EParticleType p2, Char_t sign2, Bool_t onlyTrue)
{
//
// Adds a new pair definition to create a new histogram,
// for the signal evaluation (same event) step.
// If the pair defs array is NULL, it is initialized here.
//
	if (!fPairDefs) fPairDefs = new TObjArray(0);
	fPairDefs->AddLast( new AliPairDef(p1, sign1, p2, sign2, fTrueMotherPDG, onlyTrue) );
}
//--------------------------------------------------------------------------------------------------------
void AliRsnAnalysis::Clear(Option_t* /* option */)
{
//
// Clear heap
//
	fHistograms->Clear("C");
	fPairDefs->Clear("C");
		
	fMixHistograms->Clear("C");
	fMixPairDefs->Clear("C");
	
	Int_t i;
	for (i = 0; i < AliPID::kSPECIES; i++) fCuts[i]->Clear("C");
	fPairCuts->Clear("C");
}
//--------------------------------------------------------------------------------------------------------
Stat_t AliRsnAnalysis::EventMix(Int_t nmix, Int_t multDiffMax, Double_t vzDiffMax, Bool_t compareTotMult)
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
// EXAMPLE:
// analysis->AddMixPairDef(AliRsnEvent::kPion, '+', AliRsnEvent::kKaon, '-');
// analysis->AddMixPairDef(AliRsnEvent::kKaon, '-', AliRsnEvent::kPion, '+');
//
	// allocate the histograms array
	Int_t i, npairdefs = (Int_t)fMixPairDefs->GetEntries();
	fMixHistograms = new TObjArray(npairdefs);
	TObjArray &histos = *fMixHistograms;
	AliPairDef *pd = 0;
	for (i = 0; i < npairdefs; i++) {
		pd = (AliPairDef*)fMixPairDefs->At(i);
		histos[i] = new TH1D(Form("hmix_%s", pd->GetName()), pd->GetTitle(), fNBins, fHistoMin, fHistoMax);
	}
	
	// Link events branch
	AliRsnEvent *event2 = new AliRsnEvent;
	fEventsTree->SetBranchAddress("events", &event2);
	
	// define variables to store data about particles
	Stat_t nPairs = 0;
	Int_t iev, ievmix, nEvents = (Int_t)fEventsTree->GetEntries();
	
	// loop on events
	Int_t nmixed, mult1, mult2, diffMult;
	Double_t vz1, vz2, diffVz;
	for (iev = 0; iev < nEvents; iev++) {
	
		// get event
		event2->Clear("DELETE");
		fEventsTree->GetEntry(iev);
		
		// copy this event into 'event 1'
		AliRsnEvent *event1 = new AliRsnEvent(*event2);
		
		// get Z position of primary vertex
		vz1 = event1->GetPrimaryVertexZ();
		
		// if total multiplicities must be used
		// it is computed here
		if (compareTotMult) {
			mult1 = event1->GetMultiplicity();
		}
		else {
			mult1 = 0;
		}
		
		// message
		if (iev % 10 == 0) cout << "\rMixing event " << iev << flush;
		
		// loop on pair definitions
		for (i = 0; i < npairdefs; i++) {
			
			// get pair definition
			pd = (AliPairDef*)fMixPairDefs->At(i);
			
			// get histogram reference
			TH1D *histogram = (TH1D*)histos[i];
			
			// get multiplicity of particle type '2' in event '1'
			if (!mult1) {
				mult1 = (Int_t)event1->GetTracks(pd->GetSign2(), pd->GetParticle2())->GetEntries();
			}
			
			// loop on events for mixing
			nmixed = 0;
			ievmix = iev;
			while (nmixed < nmix) {
				
				// get other event (it starts from iev + 1, and 
				// loops again to 0 if reachs end of tree)
				ievmix++;
				if (ievmix >= nEvents) ievmix -= nEvents;
				
				// if event-2-mix is equal to 'iev', then
				// a complete loop has been done, and no more
				// mixing can be done, even if they are too few
				// then the loop breaks anyway
				if (ievmix == iev) break;
				
				// get other event
				event2->Clear("DELETE");
				fEventsTree->GetEntry(ievmix);
				
				// get event stats related to event 2
				vz2 = event2->GetPrimaryVertexZ();
				if (compareTotMult) {
					mult2 = event2->GetMultiplicity();
				}
				else {
					mult2 = event2->GetTracks(pd->GetSign2(), pd->GetParticle2())->GetEntries();
				}
				
				// fill histogram if checks are passed
				diffMult = TMath::Abs(mult1 - mult2);
				diffVz = TMath::Abs(vz1 - vz2);
				if (diffVz <= vzDiffMax && diffMult <= multDiffMax) {
					//cout << ievmix << " " << flush;
					nPairs += Compute(pd, histogram, event1, event2);
					nmixed++;
				}
			}
		}
		
		event1->Clear("DELETE");
		delete event1;
		event1 = 0;
		
	} // end loop on events
	cout << endl;
	
	return nPairs;
}
//--------------------------------------------------------------------------------------------------------
Stat_t AliRsnAnalysis::Process()
{
//
// Reads the list 'fPairDefs', and builds an inv-mass histogram for each definition.
// For each event, particle of 'type 1' are combined with particles of 'type 2' as 
// defined in each pair definition specified in the list, taking the list of both types
// from the same event.
// This method is supposed to evaluate the 'signal' of the resonance, containing the peak.
// It can also be used when one wants to evaluate the 'background' with the 'wrong sign'
// particle pairs.
//
	// allocate the histograms array in order to contain 
	// as many objects as the number of pair definitionss
	Int_t i, npairdefs = (Int_t)fPairDefs->GetEntries();
	fHistograms = new TObjArray(npairdefs);
	TObjArray &histos = *fHistograms;
	
	// loop over pair defs in order to retrieve the particle species chosen
	// which are used to set an unique name for the output histogram
	// There will be a direct correspondence between the inder of pairdef in its array
	// and the corresponding histogram in the 'fHistograms' list.
	AliPairDef *pd = 0;
	for (i = 0; i < npairdefs; i++) {
		pd = (AliPairDef*)fPairDefs->At(i);
		histos[i] = new TH1D(Form("h_%s", pd->GetName()), pd->GetTitle(), fNBins, fHistoMin, fHistoMax);
	}
	
	// Link events branch
	AliRsnEvent *event = new AliRsnEvent;
	fEventsTree->SetBranchAddress("events", &event);
	
	// define variables to store data about particles
	Stat_t nPairs = 0;
	Int_t iev, nEvents = (Int_t)fEventsTree->GetEntries();
	
	// loop on events
	for (iev = 0; iev < nEvents; iev++) {
	
		// get event
		event->Clear("DELETE");
		fEventsTree->GetEntry(iev);
		if (iev % 1000 == 0) cout << "\rProcessing event " << iev << flush;
			
		// loop over the collection of pair defs
		for (i = 0; i < npairdefs; i++) {
			pd = (AliPairDef*)fPairDefs->At(i);
			TH1D *histogram = (TH1D*)histos[i];
			// invoke AliPairDef method to fill histogram
			nPairs += Compute(pd, histogram, event, event);
		}
		
	} // end loop on events
	cout << endl;
	
	return nPairs;
}
//--------------------------------------------------------------------------------------------------------
void AliRsnAnalysis::WriteHistograms() const
{
//
// Writes histograms in current directory
//
	TH1D *histogram;
	TObjArrayIter iter(fHistograms);
	while ( (histogram = (TH1D*)iter.Next()) ) {
		histogram->Write();
	}
	
	if (fMixHistograms) {
		TObjArrayIter iterMix(fMixHistograms);
		while ( (histogram = (TH1D*)iterMix.Next()) ) {
			histogram->Write();
		}
	}
}
//--------------------------------------------------------------------------------------------------------
Stat_t AliRsnAnalysis::Compute
(AliPairDef *pd, TH1D* &histogram, AliRsnEvent *event1, AliRsnEvent *event2)
{
//
// Adds to the specified histogram the invariant mass spectrum calculated taking
// particles of type 1 from event 1, and particles of type 2 from event 2.
// Events can be equal (signal) or different (background with event mixing).
//
	// define two 'cursor' objects
	AliRsnDaughter *track1 = 0, *track2 = 0;
	
	// define iterators for the two collections
	if (!event1 || !event1->GetTracks(pd->GetSign1(), pd->GetParticle1())) {
		Error("Compute", "Null pointer for particle 1 collection");
		return 0.0;
	}
	if (!event2 || !event2->GetTracks(pd->GetSign2(), pd->GetParticle2())) {
		Error("Compute", "Null pointer for particle 2 collection");
		return 0.0;
	}
	TObjArrayIter iter1(event1->GetTracks(pd->GetSign1(), pd->GetParticle1()));
	TObjArrayIter iter2(event2->GetTracks(pd->GetSign2(), pd->GetParticle2()));
	
	// define temporary variables for better code readability
	Stat_t nPairs = 0;
	
	// loop on particle of type 1 (in event 1)
	while ( (track1 = (AliRsnDaughter*)iter1.Next()) ) {

		if (fRejectFakes && (track1->GetLabel() < 0)) continue;
		if (!SingleCutCheck(pd->GetParticle1(), track1)) continue;
		
		iter2.Reset();
		
		// loop on particles of type 2 (in event 2)
		while ( (track2 = (AliRsnDaughter*)iter2.Next()) ) {
			
			if (fRejectFakes && (track2->GetLabel() < 0)) continue;
			if (event1 == event2 && track1->GetIndex() == track2->GetIndex()) continue;
			if (!SingleCutCheck(pd->GetParticle2(), track2)) continue;
			if (!PairCutCheck(track1, track2)) continue;
					
			/*
			// check
			Char_t sign1 = (track1->GetSign() > 0) ? '+' : '-';
			Char_t sign2 = (track2->GetSign() > 0) ? '+' : '-';
			Info("Compute", "Particle 1: PDG = %d, sign = %c --- Particle 2: PDG = %d, sign = %c", track1->GetTruePDG(), sign1, track2->GetTruePDG(), sign2);
			*/
				
			// total 4-momentum
			track1->SetMass(pd->GetMass1());
			track2->SetMass(pd->GetMass2());
			AliRsnDaughter sum = AliRsnDaughter::Sum(*track1, *track2);
				
			// if the choice to get ONLY true pairs is selected, a check is made on the mothers
			if (pd->GetOnlyTrue()) {
				// the 'sum' AliRsnDaughter object is initialized with
				// the PDG code of the common mother (if it is the case)
				if (TMath::Abs(sum.GetMotherPDG()) == TMath::Abs(fTrueMotherPDG)) {
					histogram->Fill(sum.GetMass());
					nPairs++;
				}
			}
			else {
				histogram->Fill(sum.GetMass());
				nPairs++;
			}	
				
		} // end loop 2
				
	} // end loop 1
	
	return nPairs;
}
//--------------------------------------------------------------------------------------------------------
Bool_t AliRsnAnalysis::SingleCutCheck(Int_t itype, AliRsnDaughter *track) const
{
//
// Checks a track against single particle cuts (if defined)
//
	if (!fCuts[itype]) return kTRUE;
	
	TObjArrayIter iter(fCuts[itype]);
	AliRsnDaughterCut *cut = 0;
	while ( (cut = (AliRsnDaughterCut*)iter.Next()) ) {
		if (!cut->Pass(track)) return kFALSE;
	}
	
	return kTRUE;
}
//--------------------------------------------------------------------------------------------------------
Bool_t AliRsnAnalysis::PairCutCheck(AliRsnDaughter *track1, AliRsnDaughter *track2) const
{
//
// Checks a pair against pair cuts (if defined)
//
	if (!fPairCuts) return kTRUE;
	
	TObjArrayIter iter(fPairCuts);
	AliRsnDaughterCut *cut = 0;
	while ( (cut = (AliRsnDaughterCut*)iter.Next()) ) {
		if (!cut->Pass(track1, track2)) return kFALSE;
	}
	
	return kTRUE;
}
//--------------------------------------------------------------------------------------------------------
AliRsnAnalysis::AliPairDef::AliPairDef
(AliPID::EParticleType p1, Char_t sign1, 
 AliPID::EParticleType p2, Char_t sign2, Int_t pdgMother, Bool_t onlyTrue)
{
//
// Constructor for nested class
//
	fOnlyTrue = onlyTrue;
	fTrueMotherPDG = 0;
	if (fOnlyTrue) fTrueMotherPDG = pdgMother;
	
	fSign1 = sign1;
	fParticle1 = p1;
	
	fSign2 = sign2;
	fParticle2 = p2;
	
	// instance a PDG database to recovery true masses of particles
	TDatabasePDG *db = TDatabasePDG::Instance();
	Int_t pdg1 = AliPID::ParticleCode((Int_t)p1);
	Int_t pdg2 = AliPID::ParticleCode((Int_t)p2);
	fMass1 = db->GetParticle(pdg1)->Mass();
	fMass2 = db->GetParticle(pdg2)->Mass();
	
	// set name according to the chosen particles
	fName.Append(Form("%s(%c)_%s(%c)", ParticleName(p1), sign1, ParticleName(p2), sign2));
	fTitle.Append(Form("Inv. mass for %s(%c) and %s(%c)", ParticleName(p1), sign1, ParticleName(p2), sign2));
	if (onlyTrue) {
		fName.Append("_true");
		fTitle.Append(" (true pairs)");
	}
}
//--------------------------------------------------------------------------------------------------------
Text_t* AliRsnAnalysis::AliPairDef::ParticleName(AliPID::EParticleType part) const
{
//
// [PRIVATE]
// Returns the name of the particle in text format
//
	if (part == AliPID::kElectron) return ("E");
	else if (part == AliPID::kMuon) return ("Mu");
	else if (part == AliPID::kPion) return ("Pi");
	else if (part == AliPID::kKaon) return ("K");
	else if (part == AliPID::kProton) return ("P");
	else {
		Warning("ParticleName", "Unrecognized value of EParticle argument");
		return ("???");
	}
}
//--------------------------------------------------------------------------------------------------------
