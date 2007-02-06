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

#include "Riostream.h"
#include "TString.h"

/// \class AliMUONDigit
/// A class representing a digit in the MUON spectrometer
/// either in tracking or trigger chambers.
///
/// A digit holds the signal (proportional to a charge) on a pad
/// (or strip).
/// 
/// This class is used to represent either sdigits (purely simulated digit, 
/// with no electronic noise whatsoever) or digits (either real digits or
/// simulated ones but including electronic noise and de-calibration, to 
/// closely ressemble real ones).

/// \cond CLASSIMP
ClassImp(AliMUONDigit)
/// \endcond

//_____________________________________________________________________________
AliMUONDigit::AliMUONDigit()
: 
TObject(),
fDetElemId(-1),
fManuId(-1),
fManuChannel(-1),
fSignal(0.0),
fPadX(-1),
fPadY(-1),
fCathode(-1),
fADC(0),
fFlags(0),
fNtracks(0),
fTcharges(0x0),
fTracks(0x0),
fPhysics(0.0),
fHit(0),
fStatusMap(0)
{
  /// Default constructor
}

//_____________________________________________________________________________
AliMUONDigit::AliMUONDigit(const AliMUONDigit& digit)
: TObject(digit),
fDetElemId(-1),
fManuId(-1),
fManuChannel(-1),
fSignal(0.0),
fPadX(-1),
fPadY(-1),
fCathode(-1),
fADC(0),
fFlags(0),
fNtracks(0),
fTcharges(0x0),
fTracks(0x0),
fPhysics(0.0),
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

//_____________________________________________________________________________
Int_t AliMUONDigit::Compare(const TObject *obj) const
{
  /// The order defined below is first by DE, then Signal, then 
  /// manuId, and then manuChannel, i.e. it should be a total ordering...

  const AliMUONDigit* d = static_cast<const AliMUONDigit*>(obj);
  
  if ( DetElemId() > d->DetElemId() ) 
  {
    return 1;
  }
  else if ( DetElemId() < d->DetElemId() )
  {
    return -1;
  }
  else
  {
    if ( Signal() > d->Signal() )
    {
      return 1;
    }
    else if ( Signal() < d->Signal() )
    {
      return -1;
    }
    else
    {
      if ( ManuId() < d->ManuId() )
      {
        return 1;
      }
      else if ( ManuId() > d->ManuId() )
      {
        return -1;
      }
      else
      {
        return ( ManuChannel() < d->ManuChannel() ) ? 1 : -1;
      }
    }
  }
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
  
  digit.fPhysics = fPhysics;
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
    fFlags ^= fgkNoiseOnlyMask;
  }
}

//_____________________________________________________________________________
AliMUONDigit& 
AliMUONDigit::operator=(const AliMUONDigit& digit)
{
  /// Assignement operator.

  AliMUONDigit a(digit);
  a.Copy(*this);
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
AliMUONDigit::Print(Option_t* opt) const
{
  /// Dump to screen.
  /// If opt=="tracks", info on tracks are printed too.

  cout << Form("<AliMUONDigit>: DE %4d Cath %d (Ix,Iy)=(%3d,%3d) (Manu,Channel)=(%4d,%2d)"
               ", Signal=%7.2f Physics=%7.2f",
               DetElemId(),Cathode(),PadX(),PadY(),ManuId(),ManuChannel(),Signal(),
               Physics());
  
  if ( IsSaturated() ) 
  {
    cout << "(S)";
  }
  else
  {
    cout << "   ";
  }
  cout << " ADC=" << setw(4) << ADC();
  cout << " Flags=0x" << setw(4) << hex << setfill('0') << fFlags << dec
    << setfill(' ');
  cout << " StatusMap=0x" << setw(4) << hex << setfill('0') << StatusMap() << dec
    << setfill(' ');

  TString options(opt);
  options.ToLower();
  if ( options.Contains("tracks") )
  {
    cout << " Hit " << setw(3) << Hit();
    Int_t ntracks = Ntracks();
    if (ntracks) 
    {
      cout << " Tracks : " << setw(2) << ntracks;
      for ( Int_t i = 0; i < ntracks; ++i )
      {
        cout << " Track(" << i << ")=" << setw(3) << Track(i)
        << " Charge(" << i << ")=" << setw(5) << TrackCharge(i);
      }
    }
    else
    {
      cout << " no track info.";
    }
  }
  cout << endl;  
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
    fFlags ^= fgkSaturatedMask;
  }
}

//_____________________________________________________________________________
void
AliMUONDigit::SetElectronics(Int_t manuId, Int_t manuChannel)
{
  //
  //FIXME: should we check that the values are ok here ??
  //
  fManuId=manuId;
  fManuChannel=manuChannel;
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
