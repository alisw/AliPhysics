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
#include <TRefArray.h>
#include <TClonesArray.h>

#include "AliLog.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliRsnDaughter.h"

#include "AliRsnEvent.h"

ClassImp(AliRsnEvent)

//_____________________________________________________________________________
AliRsnEvent::AliRsnEvent() :
  fSource(kUnknown),
  fPVx(0.0),
  fPVy(0.0),
  fPVz(0.0),
  fTracks(0x0),
  fPos(0x0),
  fNeg(0x0)
{
//=========================================================
// Default constructor 
// (implemented but not recommended for direct use)
//=========================================================

    Int_t i;
    for (i = 0; i <= AliRsnPID::kSpecies; i++) {
        fPosID[i] = 0x0;
        fNegID[i] = 0x0;
    }
}

//_____________________________________________________________________________
AliRsnEvent::AliRsnEvent(const AliRsnEvent &event) :
  TObject((TObject)event),
  fSource(event.fSource),
  fPVx(event.fPVx),
  fPVy(event.fPVy),
  fPVz(event.fPVz),
  fTracks(0x0),
  fPos(0x0),
  fNeg(0x0)
{
//=========================================================
// Copy constructor.
// Creates new instances of all collections 
// to store a copy of all objects.
//=========================================================
    
    // initialize arrays
    Init();
    
    // duplcate entries
    AliRsnDaughter *track = 0;
    TObjArrayIter iter(event.fTracks);
    while ( (track = (AliRsnDaughter*)iter.Next()) ) {
        AliRsnDaughter *ref = AddTrack(*track);
        if (!ref) AliWarning(Form("Problem occurred when copying track #%d", fTracks->IndexOf(ref)));
    }
    
    // fill PID arrays
    FillPIDArrays();
}

//_____________________________________________________________________________
AliRsnEvent& AliRsnEvent::operator=(const AliRsnEvent &event)
{
//=========================================================
// Assignment operator.
// Creates new instances of all collections 
// to store a copy of all objects.
//=========================================================
    
    // copy source info
    fSource = event.fSource;
    
    // copy primary vertex and initialize track counter to 0
    fPVx = event.fPVx;
    fPVy = event.fPVy;
    fPVz = event.fPVz;
        
    // initialize with size of argument
    Init();

    // loop on collection of argument's tracks and store a copy here
    AliRsnDaughter *track = 0;
    TObjArrayIter iter(event.fTracks);
    while ( (track = (AliRsnDaughter*)iter.Next()) ) {
        AliRsnDaughter *ref = AddTrack(*track);
        if (!ref) AliWarning(Form("Problem occurred when copying track #%d", fTracks->IndexOf(ref)));
    }
    
    // fill PID arrays
    FillPIDArrays();
    
    // return the newly created object
    return (*this);
}

//_____________________________________________________________________________
AliRsnEvent::~AliRsnEvent()
{
//=========================================================
// Destructor.
// Deletes the collection objects.
// If statements are present because if the event is 
// destroyed before that any track is added to it, then
// its collection classes will not be initialized.
//=========================================================
    
    Clear();
    
    if (fTracks) delete fTracks;
    if (fPos) delete fPos;
    if (fNeg) delete fNeg;
    
    Int_t i;
    for (i = 0; i <= AliRsnPID::kSpecies; i++) {
        if (fPosID[i]) delete fPosID[i];
        if (fNegID[i]) delete fNegID[i];
    }
}

//_____________________________________________________________________________
void AliRsnEvent::Init()
{
//=========================================================
// Initialize arrays
//=========================================================
    
    fTracks = new TClonesArray("AliRsnDaughter", 0);
    fTracks->BypassStreamer(kFALSE);
    fPos = new TRefArray;
    fNeg = new TRefArray;
    
    Int_t i;
    for (i = 0; i <= AliRsnPID::kSpecies; i++) {
        fPosID[i] = new TRefArray;
        fNegID[i] = new TRefArray;
    }
}

//_____________________________________________________________________________
AliRsnDaughter* AliRsnEvent::AddTrack(AliRsnDaughter track)
{
//=========================================================
// Stores a track into the array and proper references.
//=========================================================
    
    Int_t nextIndex = fTracks->GetEntriesFast();
    TClonesArray &tracks = (*fTracks);
    AliRsnDaughter *copy = new (tracks[nextIndex]) AliRsnDaughter(track);
    if (!copy) return 0x0;
    if (copy->Charge() > 0) {
        fPos->Add(copy);
        return copy;
    }
    else if (copy->Charge() < 0) {
        fNeg->Add(copy);
        return copy;
    }
    else {
        return 0x0;
    }
}

//_____________________________________________________________________________
void AliRsnEvent::Clear(Option_t* /*option*/)
{
//=========================================================
// Empties the collections (does not delete the objects).
// The track collection is emptied only at the end.
// Again, since the objects could be uninitialized, some
// if statement are used.
//=========================================================
    
    if (fPos) fPos->Delete();
    if (fNeg) fNeg->Delete();
    
    Int_t i;
    for (i = 0; i <= AliRsnPID::kSpecies; i++) {
        if (fPosID[i]) fPosID[i]->Delete();
        if (fNegID[i]) fNegID[i]->Delete();
    }
    
    if (fTracks) fTracks->Delete();
}

//_____________________________________________________________________________
void AliRsnEvent::Print(Option_t *option) const
{
//=========================================================
// Lists the details of the event, and the ones of each
// contained track, usind the Dump method of the track.
// The options are passed to AliRsnDaughter::Print().
// Look at that method to understand option values.
//=========================================================

    cout << "...Multiplicity     : " << fTracks->GetEntries() << endl;
    cout << "...Primary vertex   : " << fPVx << ' ' << fPVy << ' ' << fPVz << endl;
    
    TObjArrayIter iter(fTracks);
    AliRsnDaughter *d = 0;
    while ( (d = (AliRsnDaughter*)iter.Next()) ) {
        cout << "....Track #" << fTracks->IndexOf(d) << endl;
        d->Print(option);
    }
}

//_____________________________________________________________________________
Int_t AliRsnEvent::GetMultiplicity() const 
{
//=========================================================
// Get number of all tracks
//=========================================================
    
    if (!fTracks) return 0;
    return fTracks->GetEntries();
}

//_____________________________________________________________________________
Int_t AliRsnEvent::GetNPos() const 
{
//=========================================================
// Get number of positive tracks
//=========================================================
    
    if (!fPos) return 0;
    return fPos->GetEntries();
}

//_____________________________________________________________________________
Int_t AliRsnEvent::GetNNeg() const 
{
//=========================================================
// Get number of negative tracks
//=========================================================
    
    if (!fNeg) return 0;
    return fNeg->GetEntries();
}

//_____________________________________________________________________________
TRefArray* AliRsnEvent::GetTracks(Char_t sign, AliRsnPID::EType refType)
{
//=========================================================
// Returns the particle collection specified in argument.
// Arguments :
//   1) sign of particle ('+' or '-')
//   2) PID of particle (from AliRsnPID::EType)
//=========================================================
    
    if (sign == '+') {
        if (refType >= AliRsnPID::kElectron && refType <= AliRsnPID::kSpecies) {
            return fPosID[refType];
        }
        else {
            AliError(Form("Index %d out of range", refType));
            return 0x0;
        }
    }
    else if (sign == '-') {
        if (refType >= AliRsnPID::kElectron && refType <= AliRsnPID::kSpecies) {
            return fNegID[refType];
        }
        else {
            AliError(Form("Index %d out of range", refType));
            return 0x0;
        }
    }
    else {
        AliError(Form("Character '%c' not recognized as charge sign", sign));
        return 0x0;
    }
}

//_____________________________________________________________________________
void AliRsnEvent::FillPIDArrays()
{
//=========================================================
// Initializes and fills all the TRefArrays containing 
// references to particles identified as each available
// PID type (from AliRsnPID).
// This method is the unique way to do this, because it 
// cannot be assumed that the track PID is known when
// it is added to this event by AddTrack.
// Of course, if tracks have not been identified, this 
// method will do nothing and all tracks will be placed
// in the TRefArray of 'Unknown' particles and the 
// GetTracks() method will not work properly.
//=========================================================

    // clear arrays if they are present
    // initialize them if they are still 0x0
    Int_t i;
    for (i = 0; i <= AliRsnPID::kSpecies; i++) {
        if (fPosID[i]) fPosID[i]->Delete(); else fPosID[i] = new TRefArray;
        if (fNegID[i]) fNegID[i]->Delete(); else fNegID[i] = new TRefArray;
    }
    
    // loop on tracks and create references
    Short_t           charge;
    AliRsnPID::EType  type;
    AliRsnDaughter   *track = 0;
    TObjArrayIter     iter(fTracks);
    while ( (track = (AliRsnDaughter*)iter.Next()) ) {
        charge = track->Charge();
        type = track->PIDType();
        i = (Int_t)type;
        if (charge > 0) {
            fPosID[i]->Add(track);
        }
        else if (charge < 0) {
            fNegID[i]->Add(track);
        }
        else {
            AliError("Found particle with ZERO charge!!!");
        }
    }
}

