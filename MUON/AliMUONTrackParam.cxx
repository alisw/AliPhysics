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

///////////////////////////////////////////////////
//
// Track parameters
// in
// ALICE
// dimuon
// spectrometer
//
///////////////////////////////////////////////////

#include <Riostream.h>

#include "AliMUONTrackParam.h"
#include "AliESDMuonTrack.h"
#include "AliLog.h"
#include "AliMUONHitForRec.h"

/// \cond CLASSIMP
ClassImp(AliMUONTrackParam) // Class implementation in ROOT context
/// \endcond

  //_________________________________________________________________________
AliMUONTrackParam::AliMUONTrackParam()
  : TObject(),
    fInverseBendingMomentum(0.),
    fBendingSlope(0.),
    fNonBendingSlope(0.),
    fZ(0.),
    fBendingCoor(0.),
    fNonBendingCoor(0.),
    fHitForRecPtr(0x0)
{
  /// Constructor
}

  //_________________________________________________________________________
AliMUONTrackParam::AliMUONTrackParam(const AliMUONTrackParam& theMUONTrackParam)
  : TObject(theMUONTrackParam),
    fInverseBendingMomentum(theMUONTrackParam.fInverseBendingMomentum), 
    fBendingSlope(theMUONTrackParam.fBendingSlope),
    fNonBendingSlope(theMUONTrackParam.fNonBendingSlope),
    fZ(theMUONTrackParam.fZ),
    fBendingCoor(theMUONTrackParam.fBendingCoor),
    fNonBendingCoor(theMUONTrackParam.fNonBendingCoor),
    fHitForRecPtr(theMUONTrackParam.fHitForRecPtr)
{
  /// Copy constructor
}

  //_________________________________________________________________________
AliMUONTrackParam& AliMUONTrackParam::operator=(const AliMUONTrackParam& theMUONTrackParam)
{
  /// Asignment operator
  if (this == &theMUONTrackParam)
    return *this;

  // base class assignement
  TObject::operator=(theMUONTrackParam);

  fInverseBendingMomentum =  theMUONTrackParam.fInverseBendingMomentum; 
  fBendingSlope           =  theMUONTrackParam.fBendingSlope; 
  fNonBendingSlope        =  theMUONTrackParam.fNonBendingSlope; 
  fZ                      =  theMUONTrackParam.fZ; 
  fBendingCoor            =  theMUONTrackParam.fBendingCoor; 
  fNonBendingCoor         =  theMUONTrackParam.fNonBendingCoor;
  fHitForRecPtr           =  theMUONTrackParam.fHitForRecPtr;

  return *this;
}

  //__________________________________________________________________________
AliMUONTrackParam::~AliMUONTrackParam()
{
/// Destructor
/// Update the number of TrackHit's connected to the attached HitForRec if any
  if (fHitForRecPtr) fHitForRecPtr->SetNTrackHits(fHitForRecPtr->GetNTrackHits() - 1); // decrement NTrackHits of hit
}

  //__________________________________________________________________________
void AliMUONTrackParam::SetTrackParam(AliMUONTrackParam& theMUONTrackParam)
{
  /// Set track parameters from "TrackParam" leaving pointer to fHitForRecPtr unchanged
  fInverseBendingMomentum =  theMUONTrackParam.fInverseBendingMomentum; 
  fBendingSlope           =  theMUONTrackParam.fBendingSlope; 
  fNonBendingSlope        =  theMUONTrackParam.fNonBendingSlope; 
  fZ                      =  theMUONTrackParam.fZ; 
  fBendingCoor            =  theMUONTrackParam.fBendingCoor; 
  fNonBendingCoor         =  theMUONTrackParam.fNonBendingCoor;
  
}

  //__________________________________________________________________________
AliMUONHitForRec* AliMUONTrackParam::GetHitForRecPtr(void) const
{
/// return pointer to HitForRec attached to the current TrackParam
/// this method should not be called when fHitForRecPtr == NULL
  if (!fHitForRecPtr) AliWarning("AliMUONTrackParam::GetHitForRecPtr: fHitForRecPtr == NULL");
  return fHitForRecPtr;
}

  //__________________________________________________________________________
Int_t AliMUONTrackParam::Compare(const TObject* TrackParam) const
{
/// "Compare" function to sort with decreasing Z (spectro. muon Z <0).
/// Returns 1 (0, -1) if Z of current TrackHit
/// is smaller than (equal to, larger than) Z of TrackHit
  if (fHitForRecPtr->GetZ() < ((AliMUONTrackParam*)TrackParam)->fHitForRecPtr->GetZ()) return(1);
  else if (fHitForRecPtr->GetZ() == ((AliMUONTrackParam*)TrackParam)->fHitForRecPtr->GetZ()) return(0);
  else return(-1);
}

  //_________________________________________________________________________
void AliMUONTrackParam::GetParamFrom(const AliESDMuonTrack& esdMuonTrack)
{
  /// assigned value form ESD track.
  fInverseBendingMomentum =  esdMuonTrack.GetInverseBendingMomentum();
  fBendingSlope           =  TMath::Tan(esdMuonTrack.GetThetaY());
  fNonBendingSlope        =  TMath::Tan(esdMuonTrack.GetThetaX());
  fZ                      =  esdMuonTrack.GetZ(); 
  fBendingCoor            =  esdMuonTrack.GetBendingCoor(); 
  fNonBendingCoor         =  esdMuonTrack.GetNonBendingCoor();
}

  //_________________________________________________________________________
void AliMUONTrackParam::SetParamFor(AliESDMuonTrack& esdMuonTrack)
{
  /// assigned value form ESD track.
  esdMuonTrack.SetInverseBendingMomentum(fInverseBendingMomentum);
  esdMuonTrack.SetThetaX(TMath::ATan(fNonBendingSlope));
  esdMuonTrack.SetThetaY(TMath::ATan(fBendingSlope));
  esdMuonTrack.SetZ(fZ); 
  esdMuonTrack.SetBendingCoor(fBendingCoor); 
  esdMuonTrack.SetNonBendingCoor(fNonBendingCoor);
}

  //__________________________________________________________________________
Double_t AliMUONTrackParam::Px() const
{
  /// return px from track paramaters
  Double_t pYZ, pZ, pX;
  pYZ = 0;
  if (  TMath::Abs(fInverseBendingMomentum) > 0 )
    pYZ = TMath::Abs(1.0 / fInverseBendingMomentum);
  pZ = -pYZ / (TMath::Sqrt(1.0 + fBendingSlope * fBendingSlope));  // spectro. (z<0)
  pX = pZ * fNonBendingSlope; 
  return pX;
}

  //__________________________________________________________________________
Double_t AliMUONTrackParam::Py() const
{
  /// return px from track paramaters
  Double_t pYZ, pZ, pY;
  pYZ = 0;
  if (  TMath::Abs(fInverseBendingMomentum) > 0 )
    pYZ = TMath::Abs(1.0 / fInverseBendingMomentum);
  pZ = -pYZ / (TMath::Sqrt(1.0 + fBendingSlope * fBendingSlope));  // spectro. (z<0)
  pY = pZ * fBendingSlope; 
  return pY;
}

  //__________________________________________________________________________
Double_t AliMUONTrackParam::Pz() const
{
  /// return px from track paramaters
  Double_t pYZ, pZ;
  pYZ = 0;
  if (  TMath::Abs(fInverseBendingMomentum) > 0 )
    pYZ = TMath::Abs(1.0 / fInverseBendingMomentum);
  pZ = -pYZ / (TMath::Sqrt(1.0 + fBendingSlope * fBendingSlope));  // spectro. (z<0)
  return pZ;
}

  //__________________________________________________________________________
Double_t AliMUONTrackParam::P() const
{
  /// return p from track paramaters
  Double_t  pYZ, pZ, p;
  pYZ = 0;
  if (  TMath::Abs(fInverseBendingMomentum) > 0 )
    pYZ = TMath::Abs(1.0 / fInverseBendingMomentum);
  pZ = -pYZ / (TMath::Sqrt(1.0 + fBendingSlope * fBendingSlope));  // spectro. (z<0)
  p = TMath::Abs(pZ) * 
    TMath::Sqrt(1.0 + fBendingSlope * fBendingSlope + fNonBendingSlope * fNonBendingSlope);
  return p;
  
}

//_____________________________________________-
void AliMUONTrackParam::Print(Option_t* opt) const
{
  /// Printing TrackParam information 
  /// "full" option for printing all the information about the TrackParam
  TString sopt(opt);
  sopt.ToUpper();
 
  if ( sopt.Contains("FULL") ) { 
    cout << "<AliMUONTrackParam> Bending P=" << setw(5) << setprecision(3)  << 1./GetInverseBendingMomentum() << 
      ", NonBendSlope=" << setw(5) << setprecision(3)  << GetNonBendingSlope()*180./TMath::Pi() <<
      ", BendSlope=" << setw(5) << setprecision(3)     << GetBendingSlope()*180./TMath::Pi()  << 
      ", (x,y,z)_IP=(" <<  setw(5) << setprecision(3) << GetNonBendingCoor() <<
      "," <<  setw(5) << setprecision(3) << GetBendingCoor() <<
      "," <<  setw(5) << setprecision(3) << GetZ() <<
      ") cm, (px,py,pz)=(" << setw(5) << setprecision(3) << Px() <<
      "," << setw(5) << setprecision(3) << Py() <<
      "," << setw(5) << setprecision(3) << Pz() << ") GeV/c" << endl;
  }
  else {
    cout << "<AliMUONTrackParam>"  << endl;
  }
    
}
