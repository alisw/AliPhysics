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
//                      Class AliRsnEvent
//                     -------------------
//           Simple collection of reconstructed tracks
//           selected from an ESD event
//           to be used for analysis.
//           .........................................
// 
// author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
//-------------------------------------------------------------------------

#include <Riostream.h>

#include <TString.h>
#include <TObjArray.h>
#include <TClonesArray.h>

#include "AliRsnDaughter.h"
#include "AliRsnEvent.h"

ClassImp(AliRsnEvent)

//--------------------------------------------------------------------------------------------------------
AliRsnEvent::AliRsnEvent() :
 fIsESD(kTRUE),
 fPath(""),
 fPVx(0.0),
 fPVy(0.0),
 fPVz(0.0),
 fMultiplicity(-1)
{
//
// Default constructor
//
	Int_t i;
	for (i = 0; i < AliPID::kSPECIES; i++) {
		fPos[i] = NULL;
		fNeg[i] = NULL;
	}
}
//--------------------------------------------------------------------------------------------------------
AliRsnEvent::AliRsnEvent(const AliRsnEvent &event) :
 TObject((TObject)event),
 fIsESD(event.fIsESD),
 fPath(event.fPath),
 fPVx(event.fPVx), 
 fPVy(event.fPVy), 
 fPVz(event.fPVz),
 fMultiplicity(event.fMultiplicity)
{
//
// Copy constructor.
// Creates new instances of all collections to store a copy of all objects.
//
	// clone tracks collections
	Int_t i;
	for (i = 0; i < AliPID::kSPECIES; i++) {
		fPos[i] = 0;
		fNeg[i] = 0;
		if (event.fPos[i]) fPos[i] = (TClonesArray*)event.fPos[i]->Clone();
		if (event.fNeg[i]) fNeg[i] = (TClonesArray*)event.fNeg[i]->Clone();
	}
}
AliRsnEvent& AliRsnEvent::operator=(const AliRsnEvent &event)
{
//
// Assignment operator.
// Creates new instances of all collections to store a copy of all objects.
//
	fIsESD = event.fIsESD;
	fPath = event.fPath;
	fPVx = event.fPVx; 
	fPVy = event.fPVy; 
	fPVz = event.fPVz;
	fMultiplicity = event.fMultiplicity;

	// clone tracks collections
	Int_t i;
	for (i = 0; i < AliPID::kSPECIES; i++) {
		fPos[i] = 0;
		fNeg[i] = 0;
		if (event.fPos[i]) fPos[i] = (TClonesArray*)event.fPos[i]->Clone();
		if (event.fNeg[i]) fNeg[i] = (TClonesArray*)event.fNeg[i]->Clone();
	}
	
	return (*this);
}
//--------------------------------------------------------------------------------------------------------//--------------------------------------------------------------------------------------------------------
void AliRsnEvent::AddTrack(AliRsnDaughter track)
{
//
// Stores a track into the correct array
//
	Int_t pdg, sign, ifirst, ilast;
	
	// if sign is zero, track is not stored
	sign = (Int_t)track.GetSign();
	if (!sign) return;
	
	// if PDG code is assigned, track is stored in the corresponding collection
	// otherwise, a copy of track is is stored in each collection (undefined track)
	pdg = track.GetPDG();
	if (pdg != 0) {
		ifirst = PDG2Enum(pdg);
		ilast = ifirst;
		if (ifirst < AliPID::kElectron || ifirst > AliPID::kProton) return;
	}
	else {
		ifirst = AliPID::kElectron;
		ilast = AliPID::kProton;
	}
	
	// track is stored
	Int_t i, index;
	for (i = ifirst; i <= ilast; i++) {
		if (sign > 0) {
			index = fPos[i]->GetEntries();
			TClonesArray &array = *fPos[i];
			new(array[index]) AliRsnDaughter(track);
		}
		else {
			index = fNeg[i]->GetEntries();
			TClonesArray &array = *fNeg[i];
			new(array[index]) AliRsnDaughter(track);
		}
	}
}
//--------------------------------------------------------------------------------------------------------
void AliRsnEvent::Clear(Option_t *option)
{
//
// Clears list of tracks and references.
// If the string "DELETE" is specified, the collection classes
// are also cleared from heap.
//
	// evaluate option
	TString opt(option);
	Bool_t deleteCollections = opt.Contains("DELETE", TString::kIgnoreCase);
	
	Int_t i;
	for (i = 0; i < AliPID::kSPECIES; i++) {
		if (fPos[i]) fPos[i]->Delete();
		if (fNeg[i]) fNeg[i]->Delete();
		if (deleteCollections) {
			delete fPos[i];
			delete fNeg[i];
			fPos[i] = 0;
			fNeg[i] = 0;
		}
	}
}
//--------------------------------------------------------------------------------------------------------
Int_t AliRsnEvent::GetMultiplicity(Bool_t recalc)
{
//
// Computes multiplicity.
// If it is already computed (fMultiplicity > -1), it returns that value,
// unless one sets the argument to kTRUE.
//
	if (fMultiplicity < 0) recalc = kTRUE;
	if (recalc) {
		fMultiplicity = 0;
		Int_t i;
		for (i = 0; i < AliPID::kSPECIES; i++) {
			fMultiplicity += (Int_t)fPos[i]->GetEntries();
			fMultiplicity += (Int_t)fNeg[i]->GetEntries();
		}
	}
	
	return fMultiplicity;
}
//--------------------------------------------------------------------------------------------------------
const char * AliRsnEvent::GetOriginFileName() const
{
//
// Returns the path where input file was stored
//
	TString str(fPath);
	if (fIsESD) {
		str.Append("/AliESDs.root");
	}
	else {
		str.Append("/galice.root");
	}
	
	return str.Data();
}
//--------------------------------------------------------------------------------------------------------
TClonesArray* AliRsnEvent::GetTracks(Char_t sign, AliPID::EParticleType type)
{
//
// Returns the particle collection specified in argument
//
	Int_t itype = (Int_t)type;
	if (itype >= 0 && itype < AliPID::kSPECIES) {
		if (sign == '+') return fPos[type]; else return fNeg[type];
	}
	else {
		return NULL;
	}
}
//--------------------------------------------------------------------------------------------------------
void AliRsnEvent::Init()
{
//
// Action 1: define default values for some data members (including pointers).
// Action 2: if 'ntracks' > 0 allocates memory to store tracks.
//
	Int_t i;
	for (i = 0; i < AliPID::kSPECIES; i++) {
		fPos[i] = new TClonesArray("AliRsnDaughter", 0);
		fNeg[i] = new TClonesArray("AliRsnDaughter", 0);
		fPos[i]->BypassStreamer(kFALSE);
		fNeg[i]->BypassStreamer(kFALSE);
	}
}
//--------------------------------------------------------------------------------------------------------
Int_t AliRsnEvent::PDG2Enum(Int_t pdgcode)
{
//
// Converts a PDG code into the correct slot in the EParticleType enumeration in AliPID
//
	Int_t i;
	for (i = 0; i < AliPID::kSPECIES; i++) {
		if (AliPID::ParticleCode((AliPID::EParticleType)i) == TMath::Abs(pdgcode)) {
			return i;
		}
	}
	
	return -1;
}
//--------------------------------------------------------------------------------------------------------
void AliRsnEvent::PrintTracks()
{
//
// Print data for particles stored in this event
//
	cout << endl;
	
	AliRsnDaughter *track = 0;
	
	Int_t type;
	for (type = 0; type < AliPID::kSPECIES; type++) {
		TObjArrayIter iterPos(fPos[type]);
		cout << "Positive " << AliPID::ParticleName(type) << "s" << endl;
		if (!fPos[type]) cout << "NOT INITIALIZED" << endl;
		while ( (track = (AliRsnDaughter*)iterPos.Next()) ) {
			cout << "Index, Sign, PDG = ";
			cout << track->GetIndex() << " ";
			cout << (Int_t)track->GetSign() << " ";
			cout << track->GetPDG() << endl;
		}
		TObjArrayIter iterNeg(fNeg[type]);
		cout << "Negative " << AliPID::ParticleName(type) << "s" << endl;
		if (!fNeg[type]) cout << "NOT INITIALIZED" << endl;
		while ( (track = (AliRsnDaughter*)iterNeg.Next()) ) {
			cout << "Index, Sign, PDG = ";
			cout << track->GetIndex() << " ";
			cout << (Int_t)track->GetSign() << " ";
			cout << track->GetPDG() << endl;
		}
	}
}
//--------------------------------------------------------------------------------------------------------
