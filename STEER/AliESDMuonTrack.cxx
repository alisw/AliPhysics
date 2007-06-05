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

///////////////////////////////////////////////////////////////////////////////
//
//  Class to describe the MUON tracks
//  in the Event Summary Data class
//  This is where the results of reconstruction
//  are stored for the muons
//  Author: G.Martinez
//
///////////////////////////////////////////////////////////////////////////////

#include "AliESDMuonTrack.h"

#include <TLorentzVector.h>
#include <TMath.h>

ClassImp(AliESDMuonTrack)

//_____________________________________________________________________________
AliESDMuonTrack::AliESDMuonTrack ():
  TObject(),
  fInverseBendingMomentum(0),
  fThetaX(0),
  fThetaY(0),
  fZ(0),
  fBendingCoor(0),
  fNonBendingCoor(0),
  fInverseBendingMomentumUncorrected(0),
  fThetaXUncorrected(0),
  fThetaYUncorrected(0),
  fZUncorrected(0),
  fBendingCoorUncorrected(0),
  fNonBendingCoorUncorrected(0),
  fChi2(0),
  fNHit(0),
  fLocalTrigger(234),
  fChi2MatchTrigger(0),
  fHitsPatternInTrigCh(0)
{
  // Default constructor
}


//_____________________________________________________________________________
AliESDMuonTrack::AliESDMuonTrack (const AliESDMuonTrack& MUONTrack):
  TObject(MUONTrack),
  fInverseBendingMomentum(MUONTrack.fInverseBendingMomentum),
  fThetaX(MUONTrack.fThetaX),
  fThetaY(MUONTrack.fThetaY),
  fZ(MUONTrack.fZ),
  fBendingCoor(MUONTrack.fBendingCoor),
  fNonBendingCoor(MUONTrack.fNonBendingCoor),
  fInverseBendingMomentumUncorrected(MUONTrack.fInverseBendingMomentumUncorrected),
  fThetaXUncorrected(MUONTrack.fThetaXUncorrected),
  fThetaYUncorrected(MUONTrack.fThetaYUncorrected),
  fZUncorrected(MUONTrack.fZUncorrected),
  fBendingCoorUncorrected(MUONTrack.fBendingCoorUncorrected),
  fNonBendingCoorUncorrected(MUONTrack.fNonBendingCoorUncorrected),
  fChi2(MUONTrack.fChi2),
  fNHit(MUONTrack.fNHit),
  fLocalTrigger(MUONTrack.fLocalTrigger),
  fChi2MatchTrigger(MUONTrack.fChi2MatchTrigger),
  fHitsPatternInTrigCh(MUONTrack.fHitsPatternInTrigCh)
{
  //
  // Copy constructor
  // Deep copy implemented
  //
}

//_____________________________________________________________________________
AliESDMuonTrack& AliESDMuonTrack::operator=(const AliESDMuonTrack& MUONTrack)
{
  // 
  // Equal operator for a deep copy
  //
  if (this == &MUONTrack)
    return *this;

  fInverseBendingMomentum = MUONTrack.fInverseBendingMomentum; 
  fThetaX                 = MUONTrack.fThetaX;           
  fThetaY                 = MUONTrack.fThetaY;           
  fZ                      = MUONTrack.fZ;                
  fBendingCoor            = MUONTrack.fBendingCoor;      
  fNonBendingCoor         = MUONTrack.fNonBendingCoor;   
  
  fInverseBendingMomentumUncorrected = MUONTrack.fInverseBendingMomentumUncorrected; 
  fThetaXUncorrected                 = MUONTrack.fThetaXUncorrected;           
  fThetaYUncorrected                 = MUONTrack.fThetaYUncorrected;           
  fZUncorrected                      = MUONTrack.fZUncorrected;                
  fBendingCoorUncorrected            = MUONTrack.fBendingCoorUncorrected;      
  fNonBendingCoorUncorrected         = MUONTrack.fNonBendingCoorUncorrected;   
  
  fChi2                   = MUONTrack.fChi2;             
  fNHit                   = MUONTrack.fNHit; 

  fLocalTrigger           = MUONTrack.fLocalTrigger;  
  fChi2MatchTrigger       = MUONTrack.fChi2MatchTrigger; 

  fHitsPatternInTrigCh    = MUONTrack.fHitsPatternInTrigCh;
 
  return *this;
}

//_____________________________________________________________________________
Double_t AliESDMuonTrack::Px() const
{
  // return p_x from track parameters
  Double_t nonBendingSlope = TMath::Tan(fThetaX);
  Double_t bendingSlope    = TMath::Tan(fThetaY);
  Double_t pYZ = (fInverseBendingMomentum != 0.) ? TMath::Abs(1. / fInverseBendingMomentum) : 0.;
  Double_t pZ  = -pYZ / TMath::Sqrt(1.0 + bendingSlope*bendingSlope);  // spectro. (z<0)
  return pZ * nonBendingSlope;
}

//_____________________________________________________________________________
Double_t AliESDMuonTrack::Py() const
{
  // return p_y from track parameters
  Double_t bendingSlope = TMath::Tan(fThetaY);
  Double_t pYZ = (fInverseBendingMomentum != 0.) ? TMath::Abs(1. / fInverseBendingMomentum) : 0.;
  Double_t pZ  = -pYZ / TMath::Sqrt(1.0 + bendingSlope*bendingSlope);  // spectro. (z<0)
  return pZ * bendingSlope;
}

//_____________________________________________________________________________
Double_t AliESDMuonTrack::Pz() const
{
  // return p_z from track parameters
  Double_t bendingSlope = TMath::Tan(fThetaY);
  Double_t pYZ = (fInverseBendingMomentum != 0.) ? TMath::Abs(1. / fInverseBendingMomentum) : 0.;
  return -pYZ / TMath::Sqrt(1.0 + bendingSlope*bendingSlope);  // spectro. (z<0)
}

//_____________________________________________________________________________
Double_t AliESDMuonTrack::P() const
{
  // return p from track parameters
  Double_t nonBendingSlope = TMath::Tan(fThetaX);
  Double_t bendingSlope    = TMath::Tan(fThetaY);
  Double_t pYZ = (fInverseBendingMomentum != 0.) ? TMath::Abs(1. / fInverseBendingMomentum) : 0.;
  Double_t pZ  = -pYZ / TMath::Sqrt(1.0 + bendingSlope*bendingSlope);  // spectro. (z<0)
  return -pZ * TMath::Sqrt(1.0 + bendingSlope*bendingSlope + nonBendingSlope*nonBendingSlope);
}

//_____________________________________________________________________________
void AliESDMuonTrack::LorentzP(TLorentzVector& vP) const
{
  // return Lorentz momentum vector from track parameters
  Double_t muonMass = 0.105658369;
  Double_t nonBendingSlope = TMath::Tan(fThetaX);
  Double_t bendingSlope    = TMath::Tan(fThetaY);
  Double_t pYZ = (fInverseBendingMomentum != 0.) ? TMath::Abs(1. / fInverseBendingMomentum) : 0.;
  Double_t pZ  = -pYZ / TMath::Sqrt(1.0 + bendingSlope*bendingSlope);  // spectro. (z<0)
  Double_t pX  = pZ * nonBendingSlope;
  Double_t pY  = pZ * bendingSlope;
  Double_t e   = TMath::Sqrt(muonMass*muonMass + pX*pX + pY*pY + pZ*pZ);
  vP.SetPxPyPzE(pX, pY, pZ, e);
}

//_____________________________________________________________________________
Double_t AliESDMuonTrack::PxUncorrected() const
{
  // return p_x from track parameters
  Double_t nonBendingSlope = TMath::Tan(fThetaXUncorrected);
  Double_t bendingSlope    = TMath::Tan(fThetaYUncorrected);
  Double_t pYZ = (fInverseBendingMomentumUncorrected != 0.) ? TMath::Abs(1. / fInverseBendingMomentumUncorrected) : 0.;
  Double_t pZ  = -pYZ / TMath::Sqrt(1.0 + bendingSlope*bendingSlope);  // spectro. (z<0)
  return pZ * nonBendingSlope;
}

//_____________________________________________________________________________
Double_t AliESDMuonTrack::PyUncorrected() const
{
  // return p_y from track parameters
  Double_t bendingSlope = TMath::Tan(fThetaYUncorrected);
  Double_t pYZ = (fInverseBendingMomentumUncorrected != 0.) ? TMath::Abs(1. / fInverseBendingMomentumUncorrected) : 0.;
  Double_t pZ  = -pYZ / TMath::Sqrt(1.0 + bendingSlope*bendingSlope);  // spectro. (z<0)
  return pZ * bendingSlope;
}

//_____________________________________________________________________________
Double_t AliESDMuonTrack::PzUncorrected() const
{
  // return p_z from track parameters
  Double_t bendingSlope = TMath::Tan(fThetaYUncorrected);
  Double_t pYZ = (fInverseBendingMomentumUncorrected != 0.) ? TMath::Abs(1. / fInverseBendingMomentumUncorrected) : 0.;
  return -pYZ / TMath::Sqrt(1.0 + bendingSlope*bendingSlope);  // spectro. (z<0)
}

//_____________________________________________________________________________
Double_t AliESDMuonTrack::PUncorrected() const
{
  // return p from track parameters
  Double_t nonBendingSlope = TMath::Tan(fThetaXUncorrected);
  Double_t bendingSlope    = TMath::Tan(fThetaYUncorrected);
  Double_t pYZ = (fInverseBendingMomentumUncorrected != 0.) ? TMath::Abs(1. / fInverseBendingMomentumUncorrected) : 0.;
  Double_t pZ  = -pYZ / TMath::Sqrt(1.0 + bendingSlope*bendingSlope);  // spectro. (z<0)
  return -pZ * TMath::Sqrt(1.0 + bendingSlope*bendingSlope + nonBendingSlope*nonBendingSlope);
}

//_____________________________________________________________________________
void AliESDMuonTrack::LorentzPUncorrected(TLorentzVector& vP) const
{
  // return Lorentz momentum vector from track parameters
  Double_t muonMass = 0.105658369;
  Double_t nonBendingSlope = TMath::Tan(fThetaXUncorrected);
  Double_t bendingSlope    = TMath::Tan(fThetaYUncorrected);
  Double_t pYZ = (fInverseBendingMomentumUncorrected != 0.) ? TMath::Abs(1. / fInverseBendingMomentumUncorrected) : 0.;
  Double_t pZ  = -pYZ / TMath::Sqrt(1.0 + bendingSlope*bendingSlope);  // spectro. (z<0)
  Double_t pX  = pZ * nonBendingSlope;
  Double_t pY  = pZ * bendingSlope;
  Double_t e   = TMath::Sqrt(muonMass*muonMass + pX*pX + pY*pY + pZ*pZ);
  vP.SetPxPyPzE(pX, pY, pZ, e);
}

//_____________________________________________________________________________
Int_t AliESDMuonTrack::GetMatchTrigger() const
{
  //  backward compatibility after replacing fMatchTrigger by fLocalTrigger
  //  0 track does not match trigger
  //  1 track match but does not pass pt cut
  //  2 track match Low pt cut
  //  3 track match High pt cut

  if (LoCircuit() == -1) {
    return 0;
  } else if (LoLpt() == 0 && LoHpt() == 0) {
    return 1;
  } else if (LoLpt() >  0 && LoHpt() == 0) {
    return 2;
  } else {
    return 3;
  }

}

