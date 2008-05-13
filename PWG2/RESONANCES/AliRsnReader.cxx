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

//
// ==== Class AliRsnReader ========
//
// This object reads a 'standard' event and converts it into the internal
// format used for resonance analysis (AliRsnEvent).
// 'Standard' event means ESD, standard AOD and MC event.
//
// The input-2-AliRsnEvent conversion is done through a class which reads
// from AliAnalysisTaskSE, which is the standard analysis object. 
// This class creates the AliRsnEvent's before the input event is read, 
// so this class has not to 'create' a new outpu event, but instead it has 
// to 'fill' one which has already been created elsewhere.
// Then, the methods provided here accept an AliRsnEvent as argument passed
// by reference, and they 'fill' this object using the data from the inputs
// passed to them.
// 
// author: A. Pulvirenti
// email : alberto.pulvirenti@ct.infn.it
//

#include "AliLog.h"

#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"

#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODVertex.h"

#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliGenEventHeader.h"

#include "AliRsnParticle.h"
#include "AliRsnDaughter.h"
#include "AliRsnEvent.h"
#include "AliRsnReader.h"

ClassImp(AliRsnReader)
		
//_____________________________________________________________________________
AliRsnReader::AliRsnReader(Bool_t checkSplit, Bool_t rejectFakes) :
  TObject(),
  fCheckSplit(checkSplit),
  fRejectFakes(rejectFakes)
{
//=========================================================
// Constructor.
// Initializes the base-type data members:
//   - management of fake tracks
//=========================================================
}

//_____________________________________________________________________________
AliRsnReader::AliRsnReader(const AliRsnReader &copy) : 
  TObject(copy),
  fCheckSplit(copy.fCheckSplit),
  fRejectFakes(copy.fRejectFakes)
{
//=========================================================
// Copy constructor.
//=========================================================
}

//_____________________________________________________________________________
AliRsnReader& AliRsnReader::operator=(const AliRsnReader &copy)
{
//=========================================================
// Assignment operator.
//=========================================================
	
	fCheckSplit = copy.fCheckSplit;
	fRejectFakes = copy.fRejectFakes;
	return (*this);
}

//_____________________________________________________________________________
Bool_t AliRsnReader::FillFromESD(AliRsnEvent *rsn, AliESDEvent *esd, AliMCEvent *mc)
{
//=========================================================
// Filler from an ESD event.
// Stores all tracks (if a filter is defined, it will store
// only the ones which survive the cuts).
// If a reference MC event is provided, it is used to store
// the MC informations for each track (true PDG code, 
// GEANT label of mother, PDG code of mother, if any).
// When this is used, the 'source' flag of the output
// AliRsnEvent object will be set to 'kESD'.
//=========================================================
	
	// set source flag
	rsn->SetSource(AliRsnEvent::kESD);
	
	// retrieve stack (if possible)
	AliStack *stack = 0x0;
	if (mc) stack = mc->Stack();
	
	// get number of tracks
	Int_t ntracks = esd->GetNumberOfTracks();
	if (!ntracks) {
	   AliWarning("No tracks in this event");
	   return kFALSE;
    }
    
    // if required with the flag, scans the event
    // and searches all split tracks (= 2 tracks with the same label);
    // for each pair of split tracks, only the better (best chi2) is kept
    // and the other is rejected: this info is stored into a Boolean array
    Int_t i1, i2, lab1, lab2;
    Bool_t *accept = new Bool_t[ntracks];
    for (i1 = 0; i1 < ntracks; i1++) accept[i1] = kTRUE;
    if (fCheckSplit) {
        for (i1 = 0; i1 < ntracks; i1++) {
            AliESDtrack *trk1 = esd->GetTrack(i1);
            lab1 = TMath::Abs(trk1->GetLabel());
            for (i2 = i1+1; i2 < ntracks; i2++) {
                AliESDtrack *trk2 = esd->GetTrack(i2);
                lab2 = TMath::Abs(trk2->GetLabel());
                // check if labels are equal
                if (lab1 == lab2) {
                    if (trk1->GetConstrainedChi2() > trk2->GetConstrainedChi2()) {
                        accept[i1] = kTRUE;
                        accept[i2] = kFALSE;
                    }
                    else {
                        accept[i1] = kFALSE;
                        accept[i2] = kTRUE;
                    }
                }
            }
        }
    }
	
	// get primary vertex
	Double_t vertex[3];
	vertex[0] = esd->GetVertex()->GetXv();
	vertex[1] = esd->GetVertex()->GetYv();
	vertex[2] = esd->GetVertex()->GetZv();
	rsn->SetPrimaryVertex(vertex[0], vertex[1], vertex[2]);
	
	// store tracks from ESD
    Int_t  index, label, labmum;
    Bool_t check;
    AliRsnDaughter temp;
    for (index = 0; index < ntracks; index++) {
        // skip track recognized as the worse one in a splitted pair
        if (!accept[index]) {
            AliInfo(Form("Rejecting split track #%d in this event", index));
            continue;
        }
        // get and check ESD track
        AliESDtrack *esdTrack = esd->GetTrack(index);
        label = esdTrack->GetLabel();
        if (fRejectFakes && (label < 0)) continue;
        // copy ESD track data into RsnDaughter
        // if unsuccessful, this track is skipped
        check = temp.Adopt(esdTrack);
        if (!check) continue;
        // if stack is present, copy MC info
        if (stack) {
            TParticle *part = stack->Particle(TMath::Abs(label));
            if (part) {
                temp.InitParticle(part);
                labmum = part->GetFirstMother();
                if (labmum >= 0) {
                    TParticle *mum = stack->Particle(labmum);
                    temp.GetParticle()->SetMotherPDG(mum->GetPdgCode());
                }
            }
        }
        // set index and label and add this object to the output container
        temp.SetIndex(index);
        temp.SetLabel(label);
        AliRsnDaughter *ptr = rsn->AddTrack(temp);
        // if problems occurred while storing, that pointer is NULL
        if (!ptr) AliWarning(Form("Failed storing track#%d", index));
    }
    
    // compute total multiplicity
    if (rsn->GetMultiplicity() <= 0) {
        AliWarning("Zero Multiplicity in this event");
        return kFALSE;
    }
    
    return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliRsnReader::FillFromAOD(AliRsnEvent *rsn, AliAODEvent *aod, AliMCEvent *mc)
{
//=========================================================
// Filler from an AOD event.
// Stores all tracks (if a filter is defined, it will store
// only the ones which survive the cuts).
// If a reference MC event is provided, it is used to store
// the MC informations for each track (true PDG code, 
// GEANT label of mother, PDG code of mother, if any).
// When this is used, the 'source' flag of the output
// AliRsnEvent object will be set to 'kAOD'.
//=========================================================
	
    // set source flag
	rsn->SetSource(AliRsnEvent::kAOD);
    
    // retrieve stack (if possible)
	AliStack *stack = 0x0;
	if (mc) stack = mc->Stack();
    
    // get number of tracks
	Int_t ntracks = aod->GetNTracks();
	if (!ntracks) {
	   AliWarning("No tracks in this event");
	   return kFALSE;
    }
	
	// get primary vertex
	Double_t vertex[3];
	vertex[0] = aod->GetPrimaryVertex()->GetX();
	vertex[1] = aod->GetPrimaryVertex()->GetY();
	vertex[2] = aod->GetPrimaryVertex()->GetZ();
	rsn->SetPrimaryVertex(vertex[0], vertex[1], vertex[2]);
	
	// store tracks from ESD
    Int_t  index, label, labmum;
    Bool_t check;
    AliAODTrack *aodTrack = 0;
    AliRsnDaughter temp;
    TObjArrayIter iter(aod->GetTracks());
    while ( (aodTrack = (AliAODTrack*)iter.Next()) ) {
        // retrieve index
        index = aod->GetTracks()->IndexOf(aodTrack);
        label = aodTrack->GetLabel();
        if (fRejectFakes && (label < 0)) continue;
        // copy ESD track data into RsnDaughter
        // if unsuccessful, this track is skipped
        check = temp.Adopt(aodTrack);
        if (!check) continue;
        // if stack is present, copy MC info
        if (stack) {
            TParticle *part = stack->Particle(TMath::Abs(label));
            if (part) {
                temp.InitParticle(part);
                labmum = part->GetFirstMother();
                if (labmum >= 0) {
                    TParticle *mum = stack->Particle(labmum);
                    temp.GetParticle()->SetMotherPDG(mum->GetPdgCode());
                }
            }
        }
        // set index and label and add this object to the output container
        temp.SetIndex(index);
        temp.SetLabel(label);
        AliRsnDaughter *ptr = rsn->AddTrack(temp);
        // if problems occurred while storin, that pointer is NULL
        if (!ptr) AliWarning(Form("Failed storing track#%d"));
    }
    
    // compute total multiplicity
    if (rsn->GetMultiplicity() <= 0) {
        AliWarning("Zero multiplicity in this event");
        return kFALSE;
    }
    
    return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliRsnReader::FillFromMC(AliRsnEvent *rsn, AliMCEvent *mc)
{
//=========================================================
// Filler from an ESD event.
// Stores all tracks which generate at least one 
// TrackReference (a point in a sensitive volume).
// In this case, the MC info is stored by default and 
// perfect particle identification is the unique available.
// When this is used, the 'source' flag of the output
// AliRsnEvent object will be set to 'kMC'.
//=========================================================
	
	// set source flag
	rsn->SetSource(AliRsnEvent::kMC);
	
	// get number of tracks
	Int_t ntracks = mc->GetNumberOfTracks();
	if (!ntracks) {
	   AliWarning("No tracks in this event");
	   return kFALSE;
    }
    
    AliStack *stack = mc->Stack();
	
	// get primary vertex
	TArrayF fvertex(3);
	Double_t vertex[3];
	mc->GenEventHeader()->PrimaryVertex(fvertex);
	vertex[0] = (Double_t)fvertex[0];
	vertex[1] = (Double_t)fvertex[1];
	vertex[2] = (Double_t)fvertex[2];
	rsn->SetPrimaryVertex(vertex[0], vertex[1], vertex[2]);
	
	// store tracks from MC
    Int_t  index, labmum;
    Bool_t check;
    AliRsnDaughter temp;
    for (index = 0; index < ntracks; index++) {
        // get and check MC track
        AliMCParticle *mcTrack = mc->GetTrack(index);
        // if particle has no track references, it is rejected
        if (mcTrack->GetNumberOfTrackReferences() <= 0) continue;
        // try to insert in the RsnDaughter its data
        check = temp.Adopt(mcTrack);
        if (!check) continue;
        labmum = temp.GetParticle()->Mother();
        if (labmum >= 0) {
            TParticle *mum = stack->Particle(labmum);
            temp.GetParticle()->SetMotherPDG(mum->GetPdgCode());
        }
        // if successful, set other data and stores it
        temp.SetIndex(index);
        temp.SetLabel(mcTrack->Label());
        AliRsnDaughter *ptr = rsn->AddTrack(temp);
        // if problems occurred while storin, that pointer is NULL
        if (!ptr) AliWarning(Form("Failed storing track#%d", index));
    }
    
    // compute total multiplicity
    if (rsn->GetMultiplicity() <= 0) {
        AliWarning("Zero multiplicity in this event");
        return kFALSE;
    }
    
    return kTRUE;
}
