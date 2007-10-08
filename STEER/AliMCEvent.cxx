/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

//-------------------------------------------------------------------------
//     Class for Kinematic Events
//     Author: Andreas Morsch, CERN
//-------------------------------------------------------------------------
#include <TArrow.h>
#include <TMarker.h>
#include <TH2F.h>
#include <TTree.h>
#include <TFile.h>
#include <TParticle.h>
#include <TClonesArray.h>
#include <TObjArray.h>

#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliTrackReference.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"


AliMCEvent::AliMCEvent():
    AliVEvent(),
    fStack(0),
    fMCParticles(new TClonesArray("AliMCParticle",1000)),
    fMCParticleMap(0),
    fHeader(0),
    fTrackReferences(0),
    fTreeTR(0),
    fTmpTreeTR(0),
    fTmpFileTR(0),
    fNprimaries(-1),
    fNparticles(-1)
{
    // Default constructor
}

AliMCEvent::AliMCEvent(const AliMCEvent& mcEvnt) :
    AliVEvent(mcEvnt),
    fStack(0),
    fMCParticles(0),
    fMCParticleMap(0),
    fHeader(0),
    fTrackReferences(0),
    fTreeTR(0),
    fTmpTreeTR(0),
    fTmpFileTR(0),
    fNprimaries(-1),
    fNparticles(-1) 
{ 
// Copy constructor
}


AliMCEvent& AliMCEvent::operator=(const AliMCEvent& mcEvnt)
{
    // assignment operator
    if (this!=&mcEvnt) { 
	AliVEvent::operator=(mcEvnt); 
    }
  
    return *this; 
}

void AliMCEvent::ConnectTreeE (TTree* tree)
{
    // Connect the event header tree
    tree->SetBranchAddress("Header", &fHeader);
}

void AliMCEvent::ConnectTreeK (TTree* tree)
{
    // Connect the kinematics tree to the stack
    fStack = fHeader->Stack();
    fStack->ConnectTree(tree);
    //
    // Load the event
    fStack->GetEvent();
    fNparticles = fStack->GetNtrack();
    fNprimaries = fStack->GetNprimary();
    AliInfo(Form("AliMCEvent: Number of particles: %5d (all) %5d (primaries)\n", 
		 fNparticles, fNprimaries));
 
    // This is a cache for the TParticles converted to MCParticles on user request
    if (fMCParticleMap) {
	fMCParticleMap->Clear();
	if (fNparticles>0) fMCParticleMap->Expand(fNparticles);}
    else
	fMCParticleMap = new TObjArray(fNparticles);
    fMCParticles->Clear();
}

void AliMCEvent::ConnectTreeTR (TTree* tree)
{
    // Connect the track reference tree
    fTreeTR = tree;
    
    if (fTreeTR->GetBranch("AliRun")) {
	if (fTmpFileTR) {
	    fTmpFileTR->Close();
	    delete fTmpFileTR;
	}
	// This is an old format with one branch per detector not in synch with TreeK
	ReorderAndExpandTreeTR();
    } else {
	// New format 
	fTreeTR->SetBranchAddress("TrackReferences", &fTrackReferences);
    }
}

Int_t AliMCEvent::GetParticleAndTR(Int_t i, TParticle*& particle, TClonesArray*& trefs)
{
    // Retrieve entry i
    if (i < 0 || i >= fNparticles) {
	AliWarning(Form("AliMCEventHandler::GetEntry: Index out of range"));
	particle = 0;
	trefs    = 0;
	return (-1);
    }
    particle = fStack->Particle(i);
    if (fTreeTR) {
	fTreeTR->GetEntry(fStack->TreeKEntry(i));
	trefs    = fTrackReferences;
	return trefs->GetEntries();
    } else {
	trefs = 0;
	return -1;
    }
}


void AliMCEvent::Clean()
{
    // Clean-up before new trees are connected
    
    if (fHeader) {
	delete fHeader;
	fHeader = 0;
    }
    
    delete fStack;

    // Clear TR
    if (fTrackReferences) {
	fTrackReferences->Clear();
	delete fTrackReferences;
	fTrackReferences = 0;
    }
}

void AliMCEvent::FinishEvent()
{
    // Clean-up after event
     fStack->Reset(0);
}


void AliMCEvent::DrawCheck(Int_t i, Int_t search)
{
    //
    // Simple event display for debugging
    if (!fTreeTR) {
	AliWarning("No Track Reference information available");
	return;
    } 
    
    if (i > -1 && i < fNparticles) {
	fTreeTR->GetEntry(fStack->TreeKEntry(i));
    } else {
	AliWarning("AliMCEvent::GetEntry: Index out of range");
    }
    
    Int_t nh = fTrackReferences->GetEntries();
    
    
    if (search) {
	while(nh <= search && i < fNparticles - 1) {
	    i++;
	    fTreeTR->GetEntry(fStack->TreeKEntry(i));
	    nh =  fTrackReferences->GetEntries();
	}
	printf("Found Hits at %5d\n", i);
    }
    TParticle* particle = fStack->Particle(i);
    
    TH2F*    h = new TH2F("", "", 100, -500, 500, 100, -500, 500);
    Float_t x0 = particle->Vx();
    Float_t y0 = particle->Vy();

    Float_t x1 = particle->Vx() + particle->Px() * 50.;
    Float_t y1 = particle->Vy() + particle->Py() * 50.;
    
    TArrow*  a = new TArrow(x0, y0, x1, y1, 0.01);
    h->Draw();
    a->SetLineColor(2);
    
    a->Draw();
    
    for (Int_t ih = 0; ih < nh; ih++) {
	AliTrackReference* ref = (AliTrackReference*) fTrackReferences->At(ih);
	TMarker* m = new TMarker(ref->X(), ref->Y(), 20);
	m->Draw();
	m->SetMarkerSize(0.4);
	
    }
}


void AliMCEvent::ReorderAndExpandTreeTR()
{
//
//  Reorder and expand the track reference tree in order to match the kinematics tree.
//  Copy the information from different branches into one
//
//  TreeTR

    fTmpFileTR = new TFile("TrackRefsTmp.root", "recreate");
    fTmpTreeTR = new TTree("TreeTR", "TrackReferences");
    if (!fTrackReferences)  fTrackReferences = new TClonesArray("AliTrackReference", 100);
    fTmpTreeTR->Branch("TrackReferences", "TClonesArray", &fTrackReferences, 32000, 0);
    

//
//  Activate the used branches only. Otherwisw we get a bad memory leak.
    fTreeTR->SetBranchStatus("*",        0);
    fTreeTR->SetBranchStatus("AliRun.*", 1);
    fTreeTR->SetBranchStatus("ITS.*",    1);
    fTreeTR->SetBranchStatus("TPC.*",    1);
    fTreeTR->SetBranchStatus("TRD.*",    1);
    fTreeTR->SetBranchStatus("TOF.*",    1);
    fTreeTR->SetBranchStatus("FRAME.*",  1);
    fTreeTR->SetBranchStatus("MUON.*",   1);
//
//  Connect the active branches
    TClonesArray* trefs[7];
    for (Int_t i = 0; i < 7; i++) trefs[i] = 0;
    if (fTreeTR){
	// make branch for central track references
	if (fTreeTR->GetBranch("AliRun")) fTreeTR->SetBranchAddress("AliRun", &trefs[0]);
	if (fTreeTR->GetBranch("ITS"))    fTreeTR->SetBranchAddress("ITS",    &trefs[1]);
	if (fTreeTR->GetBranch("TPC"))    fTreeTR->SetBranchAddress("TPC",    &trefs[2]);
	if (fTreeTR->GetBranch("TRD"))    fTreeTR->SetBranchAddress("TRD",    &trefs[3]);
	if (fTreeTR->GetBranch("TOF"))    fTreeTR->SetBranchAddress("TOF",    &trefs[4]);
	if (fTreeTR->GetBranch("FRAME"))  fTreeTR->SetBranchAddress("FRAME",  &trefs[5]);
	if (fTreeTR->GetBranch("MUON"))   fTreeTR->SetBranchAddress("MUON",   &trefs[6]);
    }

    Int_t np = fStack->GetNprimary();
    Int_t nt = fTreeTR->GetEntries();
    
    //
    // Loop over tracks and find the secondaries with the help of the kine tree
    Int_t ifills = 0;
    Int_t it     = 0;
    Int_t itlast = 0;
    TParticle* part;

    for (Int_t ip = np - 1; ip > -1; ip--) {
	part = fStack->Particle(ip);
//	printf("Particle %5d %5d %5d %5d %5d %5d \n", 
//	       ip, part->GetPdgCode(), part->GetFirstMother(), part->GetFirstDaughter(), 
//	       part->GetLastDaughter(), part->TestBit(kTransportBit));

	// Determine range of secondaries produced by this primary during transport	
	Int_t dau1  = part->GetFirstDaughter();
	if (dau1 < np) continue;  // This particle has no secondaries produced during transport
	Int_t dau2  = -1;
	if (dau1 > -1) {
	    Int_t inext = ip - 1;
	    while (dau2 < 0) {
		if (inext >= 0) {
		    part = fStack->Particle(inext);
		    dau2 =  part->GetFirstDaughter();
		    if (dau2 == -1 || dau2 < np) {
			dau2 = -1;
		    } else {
			dau2--;
		    }
		} else {
		    dau2 = fStack->GetNtrack() - 1;
		}
		inext--;
	    } // find upper bound
	}  // dau2 < 0
	

//	printf("Check (1) %5d %5d %5d %5d %5d \n", ip, np, it, dau1, dau2);
//
// Loop over reference hits and find secondary label
// First the tricky part: find the entry in treeTR than contains the hits or
// make sure that no hits exist.
//
	Bool_t hasHits   = kFALSE;
	Bool_t isOutside = kFALSE;

	it = itlast;
	while (!hasHits && !isOutside && it < nt) {
	    fTreeTR->GetEntry(it++);
	    for (Int_t ib = 0; ib < 7; ib++) {
		if (!trefs[ib]) continue;
		Int_t nh = trefs[ib]->GetEntries();
		for (Int_t ih = 0; ih < nh; ih++) {
		    AliTrackReference* tr = (AliTrackReference*) trefs[ib]->At(ih);
		    Int_t label = tr->Label();
		    if (label >= dau1 && label <= dau2) {
			hasHits = kTRUE;
			itlast = it - 1;
			break;
		    }
		    if (label > dau2 || label < ip) {
			isOutside = kTRUE;
			itlast = it - 1;
			break;
		    }
		} // hits
		if (hasHits || isOutside) break;
	    } // branches
	} // entries

	if (!hasHits) {
	    // Write empty entries
	    for (Int_t id = dau1; (id <= dau2); id++) {
		fTmpTreeTR->Fill();
		ifills++;
	    } 
	} else {
	    // Collect all hits
	    fTreeTR->GetEntry(itlast);
	    for (Int_t id = dau1; (id <= dau2) && (dau1 > -1); id++) {
		for (Int_t ib = 0; ib < 7; ib++) {
		    if (!trefs[ib]) continue;
		    Int_t nh = trefs[ib]->GetEntries();
		    for (Int_t ih = 0; ih < nh; ih++) {
			AliTrackReference* tr = (AliTrackReference*) trefs[ib]->At(ih);
			Int_t label = tr->Label();
			// Skip primaries
			if (label == ip) continue;
			if (label > dau2 || label < dau1) 
			    printf("AliMCEventHandler::Track Reference Label out of range !: %5d %5d %5d %5d \n", 
				   itlast, label, dau1, dau2);
			if (label == id) {
			    // secondary found
			    tr->SetDetectorId(ib-1);
			    Int_t nref =  fTrackReferences->GetEntriesFast();
			    TClonesArray &lref = *fTrackReferences;
			    new(lref[nref]) AliTrackReference(*tr);
			}
		    } // hits
		} // branches
		fTmpTreeTR->Fill();
		fTrackReferences->Clear();
		ifills++;
	    } // daughters
	} // has hits
    } // tracks

    //
    // Now loop again and write the primaries
    //
    it = nt - 1;
    for (Int_t ip = 0; ip < np; ip++) {
	Int_t labmax = -1;
	while (labmax < ip && it > -1) {
	    fTreeTR->GetEntry(it--);
	    for (Int_t ib = 0; ib < 7; ib++) {
		if (!trefs[ib]) continue;
		Int_t nh = trefs[ib]->GetEntries();
		// 
		// Loop over reference hits and find primary labels
		for (Int_t ih = 0; ih < nh; ih++) {
		    AliTrackReference* tr = (AliTrackReference*)  trefs[ib]->At(ih);
		    Int_t label = tr->Label();
		    if (label < np && label > labmax) {
			labmax = label;
		    }
		    
		    if (label == ip) {
			tr->SetDetectorId(ib-1);
			Int_t nref = fTrackReferences->GetEntriesFast();
			TClonesArray &lref = *fTrackReferences;
			new(lref[nref]) AliTrackReference(*tr);
		    }
		} // hits
	    } // branches
	} // entries
	it++;
	fTmpTreeTR->Fill();
	fTrackReferences->Clear();
	ifills++;
    } // tracks
    // Check


    // Clean-up
    delete fTreeTR; fTreeTR = 0;
    
    for (Int_t ib = 0; ib < 7; ib++) {
	if (trefs[ib]) {
	    trefs[ib]->Clear();
	    delete trefs[ib];
	    trefs[ib] = 0;
	}
    }

    if (ifills != fStack->GetNtrack()) 
	printf("AliMCEvent:Number of entries in TreeTR (%5d) unequal to TreeK (%5d) \n", 
	       ifills, fStack->GetNtrack());

    fTmpTreeTR->Write();
    fTreeTR = fTmpTreeTR;
}

AliMCParticle* AliMCEvent::GetTrack(Int_t i) const
{
    // Get MC Particle i
    AliMCParticle *mcParticle;
    TParticle     *particle;

    // Out of range check
    if (i < 0 || i >= fNparticles) {
	AliWarning(Form("AliMCEvent::GetEntry: Index out of range"));
	mcParticle = 0;
	return (mcParticle);
    }
    

    particle   = fStack->Particle(i);
    //
    // First check of the MC Particle has been already cached
    if(!fMCParticleMap->At(i)) {
	Int_t nentries = fMCParticles->GetEntriesFast();
	new ((*fMCParticles)[nentries]) AliMCParticle(particle);
	mcParticle = dynamic_cast<AliMCParticle*>((*fMCParticles)[nentries]);
	fMCParticleMap->AddAt(mcParticle, i);
    } else {
	mcParticle = dynamic_cast<AliMCParticle*>(fMCParticleMap->At(i));
    }

    return mcParticle;
}

 AliGenEventHeader* AliMCEvent::GenEventHeader() {return (fHeader->GenEventHeader());}

ClassImp(AliMCEvent)
