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

#include "AliMUONDigit.h"

//-----------------------------------------------------------------------------
/// \class AliMUONDigit
/// A class representing a digit (with MC information if possible)
/// in the MUON spectrometer either in tracking or trigger chambers.
///
/// A digit holds the signal (proportional to a charge) on a pad
/// (or strip).
/// 
/// This class is used to represent either sdigits (purely simulated digit, 
/// with no electronic noise whatsoever) or digits (simulated ones but 
/// including electronic noise and de-calibration, to closely ressemble real ones).
//-----------------------------------------------------------------------------

/// \cond CLASSIMP
ClassImp(AliMUONDigit)
/// \endcond

//_____________________________________________________________________________
AliMUONDigit::AliMUONDigit()
: 
AliMUONVDigit(),
fDetElemId(0),
fManuId(0),
fManuChannel(0),
fSignal(0.0),
fPadX(-1),
fPadY(-1),
fCathode(0),
fADC(0),
fFlags(0),
fNtracks(0),
fTcharges(0x0),
fTracks(0x0),
fHit(0),
fStatusMap(0)
{
  /// Default constructor
}

//_____________________________________________________________________________
AliMUONDigit::AliMUONDigit(Int_t detElemId, Int_t manuId,
                           Int_t manuChannel, Int_t cathode)
: 
AliMUONVDigit(detElemId,manuId,manuChannel,cathode),
fDetElemId(detElemId),
fManuId(manuId),
fManuChannel(manuChannel),
fSignal(0.0),
fPadX(-1),
fPadY(-1),
fCathode(cathode),
fADC(0),
fFlags(0),
fNtracks(0),
fTcharges(0x0),
fTracks(0x0),
fHit(0),
fStatusMap(0)
{
  /// Normal constructor
}


//_____________________________________________________________________________
AliMUONDigit::AliMUONDigit(const AliMUONDigit& digit)
: AliMUONVDigit(),
fDetElemId(0),
fManuId(0),
fManuChannel(0),
fSignal(0.0),
fPadX(-1),
fPadY(-1),
fCathode(0),
fADC(0),
fFlags(0),
fNtracks(0),
fTcharges(0x0),
fTracks(0x0),
fHit(0),
fStatusMap(0)
{
  /// Copy constructor

   (static_cast<const AliMUONDigit&>(digit)).Copy(*this);
}

//_____________________________________________________________________________
AliMUONDigit::~AliMUONDigit()
{
  /// Destructor 

  delete[] fTcharges;
  delete[] fTracks;
}

//_____________________________________________________________________________
void
AliMUONDigit::AddTrack(Int_t trackNumber, Float_t trackCharge)
{
  /// Add 1 track information to the track list we keep.
  /// The implementation below is dumb, you've been warned !
  
  // First check if track is already there, in which
  // case we simply increment its charge.
  for ( Int_t i = 0; i < Ntracks(); ++i )
  {
      if ( Track(i) == trackNumber ) 
      {
        fTcharges[i] += trackCharge;
        return;
      }
  }
  
  // Nope. It's a brand new track. Make a new array to get space
  // for it, copy the old array into new one, and add the track.
  Int_t* newTracks = new Int_t[fNtracks+1];
  Float_t* newTcharges = new Float_t[fNtracks+1];
  
  for ( Int_t i = 0; i < fNtracks; ++i )
  {
    newTracks[i] = fTracks[i];
    newTcharges[i] = fTcharges[i];
  }
  
  newTracks[fNtracks] = trackNumber;
  newTcharges[fNtracks] = trackCharge;
  
  delete[] fTracks;
  delete[] fTcharges;
  
  fTracks = newTracks;
  fTcharges = newTcharges;
  
  ++fNtracks;
}

//_____________________________________________________________________________
void 
AliMUONDigit::Clear(Option_t*)
{
  /// Reset this digit, in particular the internal arrays are deleted.

  delete[] fTracks;
  delete[] fTcharges;
  fTracks=0x0;
  fTcharges=0x0;
  fNtracks=0;
}

//______________________________________________________________________________
void 
AliMUONDigit::Copy(TObject& obj) const
{
  /// Copy this line to line.

  TObject::Copy(obj);
  AliMUONDigit& digit = static_cast<AliMUONDigit&>(obj);
  
  digit.fDetElemId = fDetElemId;
  digit.fManuId = fManuId;
  digit.fManuChannel = fManuChannel;
  digit.fSignal = fSignal;
  
  digit.fPadX = fPadX;
  digit.fPadY = fPadY;
  digit.fCathode = fCathode;
  digit.fADC = fADC;
  digit.fFlags = fFlags;
  
  digit.fNtracks = fNtracks;
  
  delete[] digit.fTcharges;
  delete[] digit.fTracks;
  
  if ( fNtracks )
  {
    digit.fTcharges = new Float_t[fNtracks];
    digit.fTracks = new Int_t[fNtracks];
  }
  
  for ( Int_t i=0; i<fNtracks; ++i ) 
  {
    digit.fTcharges[i] = fTcharges[i];
    digit.fTracks[i] = fTracks[i];
  }
  
  digit.fHit = fHit;
  digit.fStatusMap = fStatusMap;
}


//_____________________________________________________________________________
Bool_t
AliMUONDigit::IsNoiseOnly() const
{
  /// Whether this (simulated only) digit is only due to noise.

  return (fFlags & fgkNoiseOnlyMask );
}

//_____________________________________________________________________________
Bool_t
AliMUONDigit::IsSaturated() const
{
  /// Whether this digit is saturated or not.

  return (fFlags & fgkSaturatedMask );
}

//_____________________________________________________________________________
Bool_t
AliMUONDigit::IsCalibrated() const
{
  /// Whether this digit is calibrated or not
  
  return (fFlags & fgkCalibratedMask );
}


//_____________________________________________________________________________
Bool_t
AliMUONDigit::IsUsed() const
{
  /// Whether this digit is used or not (in a cluster, for instance)
  
  return (fFlags & fgkUsedMask );
}

//_____________________________________________________________________________
Bool_t
AliMUONDigit::IsEfficiencyApplied() const
{
  /// Whether this digit had efficiency applied or not
  
  return (fFlags & fgkEfficiencyMask );
}

//_____________________________________________________________________________
void
AliMUONDigit::Used(Bool_t value)
{
  /// Set the Used status of this digit.
  
  if ( value )
  {
    fFlags |= fgkUsedMask;
  }
  else
  {
    fFlags &= ~fgkUsedMask;
  }
}

//_____________________________________________________________________________
void
AliMUONDigit::Calibrated(Bool_t value)
{
  /// Set the Calibrated status of this digit.
  
  if ( value )
  {
    fFlags |= fgkCalibratedMask;
  }
  else
  {
    fFlags &= ~fgkCalibratedMask;
  }
}

//_____________________________________________________________________________
void
AliMUONDigit::EfficiencyApplied(Bool_t value)
{
  /// Set the EfficiencyApplied status of this digit.
  
  if ( value )
  {
    fFlags |= fgkEfficiencyMask;
  }
  else
  {
    fFlags &= ~fgkEfficiencyMask;
  }
}

//_____________________________________________________________________________
Bool_t
AliMUONDigit::MergeWith(const AliMUONVDigit& src)
{
  /// Merge with src.
  
  Bool_t check = ( src.DetElemId() == DetElemId() &&
                   src.PadX() == PadX() &&
                   src.PadY() == PadY() &&
                   src.Cathode() == Cathode() );
  if (!check)
  {
    return kFALSE;
  }
  
  AddCharge(src.Charge());
  for ( Int_t i = 0; i < src.Ntracks(); ++i )
  {
    AddTrack(src.Track(i),src.TrackCharge(i));
  }
  return kTRUE;
}

//_____________________________________________________________________________
void
AliMUONDigit::NoiseOnly(Bool_t value)
{
  /// Set the NoiseOnly status of this digit.

  if ( value )
  {
    fFlags |= fgkNoiseOnlyMask;
  }
  else
  {
    fFlags &= ~fgkNoiseOnlyMask;
  }
}

//_____________________________________________________________________________
AliMUONDigit& 
AliMUONDigit::operator=(const AliMUONDigit& digit)
{
  /// Assignement operator.

  if ( this != &digit ) 
  {
    digit.Copy(*this);
  }
  return *this;
}

//_____________________________________________________________________________
void
AliMUONDigit::PatchTracks(Int_t mask)
{
  /// Add mask to each track number.

  for ( Int_t i = 0; i < Ntracks(); ++i )
  {
    fTracks[i] += mask;
  }
}

//_____________________________________________________________________________
void
AliMUONDigit::Saturated(Bool_t value)
{
  /// Set the saturation status of this digit.

  if ( value )
  {
    fFlags |= fgkSaturatedMask;
  }
  else
  {
    fFlags &= ~fgkSaturatedMask;
  }
}

//_____________________________________________________________________________
Int_t
AliMUONDigit::Track(Int_t i) const
{
  /// Return the i-th track number (if i is >=0 and < Ntracks()) or -1.

  if ( i >= 0 && i < fNtracks ) 
  {
    return fTracks[i];
  }

  return -1;
}

//_____________________________________________________________________________
Float_t
AliMUONDigit::TrackCharge(Int_t i) const
{
  /// Return the i-th track charge (if i is >=0 and < Ntracjs()) or -1.

  if ( i >= 0 && i < fNtracks ) 
  {
    return fTcharges[i];
  }

  return -1;
}

//_____________________________________________________________________________
UInt_t
AliMUONDigit::GetUniqueID() const
{
  /// Return a single integer with id information

  return BuildUniqueID(DetElemId(),ManuId(),ManuChannel(),Cathode());
}
