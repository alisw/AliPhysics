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

ClassImp(AliMUONDigit)

//_____________________________________________________________________________
AliMUONDigit::AliMUONDigit()
: TObject(),
fPadX(0),
fPadY(0),
fCathode(0),
fSignal(0),
fPhysics(0),
fHit(0),
fDetElemId(0),
fManuId(-1),
fManuChannel(-1),
fADC(0),
fIsSaturated(kFALSE)
{
  // Default constructor
  
  for ( Int_t i=0; i<kMAXTRACKS; ++i ) 
  {
    fTcharges[i] = 0;
    fTracks[i] = 0;
  }
}

//_____________________________________________________________________________
AliMUONDigit::AliMUONDigit(const AliMUONDigit& digit)
  : TObject(digit)
{
  // copy constructor
 
   (static_cast<const AliMUONDigit&>(digit)).Copy(*this);
}


//_____________________________________________________________________________
AliMUONDigit::~AliMUONDigit()
{
  // Destructor 
}

//_____________________________________________________________________________
AliMUONDigit& 
AliMUONDigit::operator=(const AliMUONDigit& digit)
{
  AliMUONDigit a(digit);
  a.Copy(*this);
  return *this;
}

//______________________________________________________________________________
void 
AliMUONDigit::Copy(TObject& obj) const
{
  // Copy this line to line.
  
  TObject::Copy(obj);
  AliMUONDigit& digit = static_cast<AliMUONDigit&>(obj);
  digit.fPadX = fPadX;
  digit.fPadY = fPadY;
  digit.fCathode = fCathode;
  digit.fSignal = fSignal;
  digit.fHit = fHit;
  digit.fDetElemId = fDetElemId;
  digit.fManuId = fManuId;
  digit.fManuChannel = fManuChannel;
  digit.fADC = fADC;
  digit.fIsSaturated = fIsSaturated;
}


//_____________________________________________________________________________
AliMUONDigit::AliMUONDigit(Int_t *digits)
{
  //
  // Creates a MUON digit object to be updated
  //
    fPadX        = digits[0];
    fPadY        = digits[1];
    fCathode     = digits[2];
    fSignal      = digits[3];
    fPhysics     = digits[4];
    fHit         = digits[5];
    fDetElemId   = digits[6];
    fManuId = -1;
    fManuChannel = -1;
    fADC=0;
    fIsSaturated = kFALSE;
}
//_____________________________________________________________________________
AliMUONDigit::AliMUONDigit(Int_t *tracks, Int_t *charges, Int_t *digits)
{
    //
    // Creates a MUON digit object
    //
    fPadX        = digits[0];
    fPadY        = digits[1];
    fCathode     = digits[2];
    fSignal      = digits[3];
    fPhysics     = digits[4];
    fHit         = digits[5];
    fDetElemId   = digits[6];
    fManuId = -1;
    fManuChannel = -1;
    fADC=0;
    for(Int_t i=0; i<kMAXTRACKS; i++) {
      fTcharges[i]  = charges[i];
      fTracks[i]    = tracks[i];
    }
    fIsSaturated=kFALSE;
}

//_____________________________________________________________________________
Int_t AliMUONDigit::Compare(const TObject *obj) const
{
// sort by idDE

 AliMUONDigit* d = (AliMUONDigit*) obj;

 return ( fDetElemId > d->DetElemId()) ? 1 : -1;

}

//_____________________________________________________________________________
void
AliMUONDigit::Print(Option_t* opt) const
{
  cout << "DetEle " << setw(5) << DetElemId()
  << " Cath " << setw(2) << Cathode()
  << " (Ix,Iy)=(" << setw(3) << PadX() << "," << setw(3) << PadY()
  << ") "
  << " (Manu,Channel)=(" << setw(4) << ManuId() 
  << "," << setw(3) << ManuChannel() << ")"
  << " Signal=" << setw(6) << Signal()
  << " Physics=" << setw(4) << Physics();
  if ( fIsSaturated ) 
  {
    cout << "(S)";
  }
  else
  {
    cout << "   ";
  }
  cout << " ADC=" << setw(4) << ADC();
  TString options(opt);
  options.ToLower();
  if ( options.Contains("tracks") )
  {
    cout << " Track0=" << setw(3) << Track(0) 
    << " Charge0=" << setw(4) << TrackCharge(0)
    << " Track1=" << setw(3) << Track(1) 
    << " Charge1=" << setw(4) << TrackCharge(1);
  }
  cout << endl;  
}

//_____________________________________________________________________________
Bool_t
AliMUONDigit::IsSaturated() const
{
  return fIsSaturated;
}

//_____________________________________________________________________________
Int_t
AliMUONDigit::ManuId() const
{
  return fManuId;
}

//_____________________________________________________________________________
Int_t
AliMUONDigit::ManuChannel() const
{
  return fManuChannel;
}

//_____________________________________________________________________________
Int_t 
AliMUONDigit::ADC() const
{
  return fADC;
}

//_____________________________________________________________________________
void
AliMUONDigit::SetADC(Int_t adc)
{
  fADC = adc;
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
void
AliMUONDigit::Saturated(Bool_t saturated)
{
  fIsSaturated=saturated;
}


