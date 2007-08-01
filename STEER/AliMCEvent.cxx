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

/* $Id$ */
//---------------------------------------------------------------------------------
//                          Class AliMCEvent
// This class gives access to MC truth during the analysis.
// Monte Carlo truth is containe in the kinematics tree (produced particles) and 
// the tree of reference hits.
//      
// Origin: Andreas Morsch, CERN, andreas.morsch@cern.ch 
//---------------------------------------------------------------------------------



#include "AliMCEvent.h"
#include "AliTrackReference.h"
#include "AliHeader.h"
#include "AliStack.h"

#include <TTree.h>
#include <TFile.h>
#include <TParticle.h>
#include <TClonesArray.h>
#include <TDirectoryFile.h>
#include <TArrow.h>
#include <TMarker.h>
#include <TH2F.h>


ClassImp(AliMCEvent)

AliMCEvent::AliMCEvent() :
    AliVirtualEventHandler(),
    fFileE(0),
    fFileK(0),
    fFileTR(0),
    fTreeE(0),
    fTreeK(0),
    fTreeTR(0),
    fStack(0),
    fHeader(0),
    fTrackReferences(0),
    fNEvent(-1),
    fEvent(-1),
    fNprimaries(-1),
    fNparticles(-1)
{
    // Default constructor
}

AliMCEvent::AliMCEvent(const char* name, const char* title) :
    AliVirtualEventHandler(name, title),
    fFileE(0),
    fFileK(0),
    fFileTR(0),
    fTreeE(0),
    fTreeK(0),
    fTreeTR(0),
    fStack(0),
    fHeader(new AliHeader()),
    fTrackReferences(new TClonesArray("AliTrackReference", 200)),
    fNEvent(-1),
    fEvent(-1),
    fNprimaries(-1),
    fNparticles(-1)
{
    // Constructor
}
AliMCEvent::~AliMCEvent()
{ 
    // Destructor
    delete fFileE;
    delete fFileK;
    delete fFileTR;

    delete fTreeE;
    delete fTreeK;
    delete fTreeTR;
    
}

Bool_t AliMCEvent::InitIO(Option_t* /*opt*/) 
{ 
    // Initialize input
    fFileE = new TFile("galice.root");
    fFileE->GetObject("TE", fTreeE);
    fTreeE->SetBranchAddress("Header", &fHeader);
    fNEvent = fTreeE->GetEntries();
    //
    // Tree K
    fFileK = new TFile("Kinematics.root");
    //
    // Tree TR
    fFileTR = new TFile("TrackRefs.root");
    //
    // Reset the event number
    fEvent = -1;
    printf("Number of events in this directory %5d \n", fNEvent);
    return kTRUE;
    
}

Bool_t AliMCEvent::BeginEvent()
{ 
    // Read the next event
    fEvent++;
    char folder[20];
    sprintf(folder, "Event%d", fEvent);
    // TreeE
    fTreeE->GetEntry(fEvent);
    fStack = fHeader->Stack();
    // Tree K
    TDirectoryFile* dirK  = 0;
    fFileK->GetObject(folder, dirK);
    dirK->GetObject("TreeK", fTreeK);
    fStack->ConnectTree(fTreeK);
    fStack->GetEvent();
    //Tree TR 
    TDirectoryFile* dirTR = 0;
    fFileTR->GetObject(folder, dirTR);
    dirTR->GetObject("TreeTR", fTreeTR);
    fTreeTR->SetBranchAddress("TrackReferences", &fTrackReferences);
    //
    fNparticles = fStack->GetNtrack();
    fNprimaries = fStack->GetNprimary();
    printf("Event#: %5d Number of particles %5d \n", fEvent, fNparticles);
    
    return kTRUE;
}

Int_t AliMCEvent::GetParticleAndTR(Int_t i, TParticle*& particle, TClonesArray*& trefs)
{
    // Retrieve entry i
    if (i > -1 && i < fNparticles) {
	fTreeTR->GetEntry(fStack->TreeKEntry(i));
    } else {
	printf("AliMCEvent::GetEntry: Index out of range");
    }
    particle = fStack->Particle(i);
    trefs    = fTrackReferences;
    printf("Returning %5d \n", particle->GetPdgCode());
    
    return trefs->GetEntries();
    
}

void AliMCEvent::DrawCheck(Int_t i)
{
    // Retrieve entry i and draw momentum vector and hits
    if (i > -1 && i < fNparticles) {
	fTreeTR->GetEntry(fStack->TreeKEntry(i));
    } else {
	printf("AliMCEvent::GetEntry: Index out of range");
    }

    TParticle* particle = fStack->Particle(i);
    Int_t nh =  fTrackReferences->GetEntries();
    
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

    

			    
Bool_t AliMCEvent::FinishEvent()
{
    // Dummy 
    return kTRUE;
}

Bool_t AliMCEvent::Terminate()
{ 
    // Dummy 
    return kTRUE;
}

Bool_t AliMCEvent::TerminateIO()
{ 
    // Dummy
    return kTRUE;
}
    
