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
// *** Class AliRsnEvent ***
//
// A container for a collection of AliRsnDaughter objects from an event.
// Contains also the primary vertex, useful for some cuts.
// In order to retrieve easily the tracks which have been identified
// as a specific type and charge, there is an array of indexes which
// allows to avoid to loop on all tracks and have only the neede ones.
//
// authors: A. Pulvirenti (email: alberto.pulvirenti@ct.infn.it)
//          M. Vala (email: martin.vala@cern.ch)
//

#include <Riostream.h>

#include "AliLog.h"

#include "AliRsnDaughter.h"
#include "AliRsnEvent.h"
#include "AliRsnMCInfo.h"

ClassImp(AliRsnEvent)

//_____________________________________________________________________________
AliRsnEvent::AliRsnEvent() :
    TNamed("rsnEvent", ""),
    fPVx(0.0),
    fPVy(0.0),
    fPVz(0.0),
    fTracks(0x0),
    fNoPID(0x0),
    fPerfectPID(0x0),
    fRealisticPID(0x0)
{
//
// Default constructor
// (implemented but not recommended for direct use)
//
}

//_____________________________________________________________________________
AliRsnEvent::AliRsnEvent(const AliRsnEvent &event) :
    TNamed(event),
    fPVx(event.fPVx),
    fPVy(event.fPVy),
    fPVz(event.fPVz),
    fTracks(0x0),
    fNoPID(0x0),
    fPerfectPID(0x0),
    fRealisticPID(0x0)
{
//
// Copy constructor.
// Copies all the tracks from the argument's collection
// to this' one, and then recreates the PID index arrays,
// trusting on the PID informations in the copied tracks.
//

  // during track copy, counts how many faults happen
  Int_t errors = Fill(event.fTracks);
  if (errors) AliWarning(Form("%d errors occurred in copy", errors));

  // fill PID index arrays
  // FillPIDArrays();

  if (event.fNoPID) fNoPID = new AliRsnPIDIndex(* (event.fNoPID));
  if (event.fPerfectPID) fPerfectPID = new AliRsnPIDIndex(* (event.fPerfectPID));
  if (event.fRealisticPID) fRealisticPID = new AliRsnPIDIndex(* (event.fRealisticPID));
}

//_____________________________________________________________________________
AliRsnEvent& AliRsnEvent::operator= (const AliRsnEvent &event)
{
//
// Works in the same way as the copy constructor.
//
  // copy name and title
  SetName(event.GetName());
  SetTitle(event.GetTitle());

  // copy primary vertex and initialize track counter to 0
  fPVx = event.fPVx;
  fPVy = event.fPVy;
  fPVz = event.fPVz;

  // add tracks from array of argument
  Int_t errors = Fill(event.fTracks);
  if (errors) AliWarning(Form("%d errors occurred in copy", errors));

  // fill PID arrays
  // FillPIDArrays();
  if (event.fNoPID)
  {
    if (!fNoPID) fNoPID = new AliRsnPIDIndex(* (event.fNoPID));
    else (*fNoPID) = * (event.fNoPID);
  }
  if (event.fPerfectPID)
  {
    if (!fPerfectPID) fPerfectPID = new AliRsnPIDIndex(* (event.fPerfectPID));
    else (*fPerfectPID) = * (event.fPerfectPID);
  }
  if (event.fRealisticPID)
  {
    if (!fRealisticPID) fRealisticPID = new AliRsnPIDIndex(* (event.fRealisticPID));
    else (*fRealisticPID) = * (event.fRealisticPID);
  }

  // return this object
  return (*this);
}

//_____________________________________________________________________________
AliRsnEvent::~AliRsnEvent()
{
//
// Destructor.
// Deletes the TClonesArray, after clearing its content.
// Other memory-allocating arrays are cleared by their
// destructor, which is automatically called from here.
//

  Clear();
  if (fTracks) delete fTracks;
}

//_____________________________________________________________________________
void AliRsnEvent::Init()
{
//
// Initialize TClonesArray data-member.
//

  fTracks = new TClonesArray("AliRsnDaughter", 1);
  //fTracks->BypassStreamer (kFALSE);
}

//_____________________________________________________________________________
void AliRsnEvent::Clear(Option_t* /*option*/)
{
//
// Empties the collections (does not delete the objects).
// The track collection is emptied only at the end.
// Since some objects could be uninitialized, some
// "if" statement are used.
//

  if (fTracks) fTracks->Delete();
  delete fNoPID;
  fNoPID = 0x0;
  delete fPerfectPID;
  fPerfectPID = 0x0;
  delete fRealisticPID;
  fRealisticPID = 0x0;
}

//_____________________________________________________________________________
AliRsnDaughter* AliRsnEvent::AddTrack(AliRsnDaughter track)
{
//
// Stores a new track into the array and returns
// a reference pointer to it (which is NULL in case of errors).
//

  Int_t nextIndex = fTracks->GetEntriesFast();
  TClonesArray &tracks = (*fTracks);
  AliRsnDaughter *copy = new(tracks[nextIndex]) AliRsnDaughter(track);
  return copy;
}

//_____________________________________________________________________________
AliRsnDaughter* AliRsnEvent::GetTrack(Int_t index)
{
//
// Returns one track in the collection
// given the absolute index in the global TClonesArray
//
  return (AliRsnDaughter*) fTracks->UncheckedAt(index);
}

//_____________________________________________________________________________
AliRsnDaughter* AliRsnEvent::GetLeadingParticle
(Double_t ptMin, AliRsnPID::EType type, Bool_t realistic)
{
//
// Returns the particle in this event with largest transverse momentum,
// provided that this momentum is larger than the first argument
// and that the PID type correspond to the second argument.
// If one specifies "AliRsnPID::kUnknown" as second arguments, the PID check is not done.
//

  Double_t prob;
  AliRsnPID::EType trackType;
  AliRsnDaughter *track, *leading = 0x0;
  TObjArrayIter iter(fTracks);
  while ((track = (AliRsnDaughter*) iter.Next()))
  {
    if (track->Pt() < ptMin) continue;
    if (realistic)
    {
      AliRsnDaughter::SetPIDMethod(AliRsnDaughter::kRealistic);
    }
    else
    {
      AliRsnDaughter::SetPIDMethod(AliRsnDaughter::kPerfect);
    }
    trackType = track->PIDType(prob);
    if (type != AliRsnPID::kUnknown && trackType != type) continue;
    ptMin = track->Pt();
    leading = track;
  }

  return leading;
}

//_____________________________________________________________________________
Int_t AliRsnEvent::GetLastFastTrack
(Double_t ptMin, AliRsnPID::EType type, Bool_t realistic)
{
//
// Loops on the list of tracks (eventually skipping the ones which do not match
// the given PID type with the specified PID method) and returns the index of the last
// one whose transverse momentum is still larger than a specified value.
// When no tracks are found this way, the value "-1" is returned.
//

  if (realistic)
  {
    AliRsnDaughter::SetPIDMethod(AliRsnDaughter::kRealistic);
  }
  else
  {
    AliRsnDaughter::SetPIDMethod(AliRsnDaughter::kPerfect);
  }

  Double_t prob;
  Int_t i, nTracks = fTracks->GetEntries(), lastIndex = -1;
  for (i = 0; i < nTracks; i++)
  {
    AliRsnDaughter *d = (AliRsnDaughter*) fTracks->At(i);
    AliRsnPID::EType trackType = d->PIDType(prob);
    if (type != AliRsnPID::kUnknown && trackType != type) continue;
    if (d->Pt() >= ptMin) lastIndex = i;
  }

  return lastIndex;
}

//_____________________________________________________________________________
TArrayI* AliRsnEvent::GetCharged(Char_t sign)
{
//
// Returns an array with the indexes of all tracks with a given charge
// (arg can be '+' or '-'), irrespective of its PID.
// When the argument is wrong, a NULL pointer is returned.
//
  if (fNoPID) return fNoPID->GetTracksArray(sign, AliRsnPID::kUnknown);
  return 0x0;
}

//_____________________________________________________________________________
TArrayI * AliRsnEvent::GetTracksArray
(AliRsnDaughter::EPIDMethod pidtype, Char_t sign, AliRsnPID::EType type)
{
//
// Returns an array of indexes of all tracks in this event
// which match the charge sign and PID type in the arguments,
// according to one of the allowed PID methods (perfect or realistic).
// It retrieves this array from the AliRsnPIDIndex data members.
// If the arguments are wrong a NULL pointer is returned.
//

  switch (pidtype)
  {
    case AliRsnDaughter::kRealistic:
      if (fRealisticPID)
      {
        return fRealisticPID->GetTracksArray(sign, type);
      }
      break;
    case AliRsnDaughter::kPerfect:
      if (fPerfectPID)
      {
        return fPerfectPID->GetTracksArray(sign, type);
      }
      break;
    case AliRsnDaughter::kNoPID:
      if (fNoPID)
      {
        return fNoPID->GetTracksArray(sign, AliRsnPID::kUnknown);
      }
      break;
    default:
      AliError("Handled PID methods here are only fNoPID,kPerfect and kRealistic. Nothing done.");
      return 0x0;
  }
  return 0x0;
}

//_____________________________________________________________________________
void AliRsnEvent::FillPIDArrays(Int_t arraySizeInit)
{
//
// Initializes and fills the AliRsnPIDIndex objects containing
// arrays of indexes for each possible charge and PID type.
// This method is the unique way to do this, for safety reasons.
//

  if (fNoPID) delete fNoPID;
  if (fPerfectPID) delete fPerfectPID;
  if (fRealisticPID) delete fRealisticPID;
  fNoPID = new AliRsnPIDIndex(arraySizeInit);
  fPerfectPID = new AliRsnPIDIndex(arraySizeInit);
  fRealisticPID = new AliRsnPIDIndex(arraySizeInit);

  // set the default type to Realistic
  AliRsnDaughter::SetPIDMethod(AliRsnDaughter::kRealistic);

  // loop on tracks and create references
  Double_t prob;
  Int_t i, icharge, type;
  Short_t charge;
  AliRsnMCInfo *mcinfo = 0;
  AliRsnDaughter *track = 0;
  TObjArrayIter iter(fTracks);
  while ((track = (AliRsnDaughter*) iter.Next()))
  {
    charge = track->Charge();
    type = (Int_t) track->PIDType(prob);
    i = fTracks->IndexOf(track);
    mcinfo = track->GetMCInfo();
    if (charge > 0) icharge = 0;
    else if (charge < 0) icharge = 1;
    else
    {
      AliError("Found particle with ZERO charge!!!");
      continue;
    }
    // add to charged array
    fNoPID->AddIndex(i, icharge, (Int_t) AliRsnPID::kUnknown);
    // add to realistic PID array
    fRealisticPID->AddIndex(i, icharge, (Int_t) type);
    // add to perfect PID array (needs MCInfo present)
    if (mcinfo)
    {
      fPerfectPID->AddIndex(i, icharge, (Int_t) AliRsnPID::InternalType(mcinfo->PDG()));
    }
  }

  // adjusts the size of arrays
  if (fNoPID) fNoPID->SetCorrectIndexSize();
  if (fPerfectPID) fPerfectPID->SetCorrectIndexSize();
  if (fRealisticPID) fRealisticPID->SetCorrectIndexSize();
}

//_____________________________________________________________________________
void AliRsnEvent::Print(Option_t *option) const
{
//
// Lists the details of the event, and the ones of each
// contained track.
// The options are passed to AliRsnDaughter::Print().
// Look at that method to understand option values.
//

  cout << "...Multiplicity     : " << fTracks->GetEntries() << endl;
  cout << "...Primary vertex   : " << fPVx << ' ' << fPVy << ' ' << fPVz << endl;

  TObjArrayIter iter(fTracks);
  AliRsnDaughter *d = 0;
  while ((d = (AliRsnDaughter*) iter.Next()))
  {
    cout << "....Track #" << fTracks->IndexOf(d) << endl;
    d->Print(option);
  }
}

//_____________________________________________________________________________
void AliRsnEvent::MakeComputations()
{
//
// Computes all required overall variables:
// - multiplicity
// - mean phi of tracks
//

  if (!fTracks) 
  {
    fMult = 0; 
    fPhiMean = 1000.0;
  }
  else 
  {
    fMult = fTracks->GetEntries();
    if (fMult < 1) {
      fPhiMean = 1000.0;
    }
    else
    {
      fPhiMean = 0.0;
      AliRsnDaughter *d = 0;
      TObjArrayIter next(fTracks);
      while ( (d = (AliRsnDaughter*)next()) ) fPhiMean += d->Phi();
      fPhiMean /= (Double_t)fMult;
    }
  }
}

//_____________________________________________________________________________
void AliRsnEvent::CorrectByPrimaryVertex()
{
//
// Corrects the X,Y,Z position of DCA vertex of all tracks
// by the amount of stored primary vertex
//

    TObjArrayIter iter(fTracks);
    AliRsnDaughter *d = 0;
    while ((d = (AliRsnDaughter*) iter.Next())) {
        d->ShiftZero(fPVx, fPVy, fPVz);
    }
}

/*
//_____________________________________________________________________________
Int_t AliRsnEvent::GetMultiplicity() const
{
//
// Get number of all tracks
//

  if (!fTracks) return 0;
  return fTracks->GetEntries();
}
*/

//_____________________________________________________________________________
Int_t AliRsnEvent::GetNCharged(Char_t sign)
{
//
// Get number of charged tracks
//

  Int_t icharge;
  icharge = ChargeIndex(sign);
  if (icharge < 0) return 0;
  TArrayI *charged = GetCharged(sign);
  if (!charged) return 0;
  return charged->GetSize();
}

//_____________________________________________________________________________
Int_t AliRsnEvent::Fill(TObjArray *array)
{
//
// Fills the data-member TClonesArray of tracks with
// the ones stored in the array passed as argument.
// If this data-member is already present, it is cleared.
// Returns the number of tracks which raised problems
// while attempting to add them. Zero is the best.
//

  // clear the array if it is already instantiated,
  // create if otherwise
  if (fTracks) fTracks->Delete();
  else Init();

  // copy argument entries into data-member
  Int_t errors = 0;
  AliRsnDaughter *track = 0;
  TObjArrayIter iter(array);
  while ((track = (AliRsnDaughter*) iter.Next()))
  {
    AliRsnDaughter *ref = AddTrack(*track);
    if (!ref)
    {
      AliWarning(Form("Problem occurred when copying track #%d from passed array", array->IndexOf(track)));
      errors++;
    }
  }

  return errors;
}

//_____________________________________________________________________________
Int_t AliRsnEvent::ChargeIndex(Char_t sign) const
//
// Returns the array index corresponding to charge
// 0 for positive, 1 for negative
//
{
  if (sign == '+') return 0;
  else if (sign == '-') return 1;
  else
  {
    AliError(Form("Character '%c' not recognized as charge sign", sign));
    return -1;
  }
}
