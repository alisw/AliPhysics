#ifndef ALIESDMUONTRACK_H
#define ALIESDMUONTRACK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//  Class to describe the MUON tracks
//  in the Event Summary Data class
//  Author: G.Martinez


#include "TObject.h"

class AliESDMuonTrack : public TObject {
public:
  AliESDMuonTrack(){} //Constructor
  virtual ~AliESDMuonTrack(){} // Destructor
  AliESDMuonTrack(const AliESDMuonTrack& esdm);
  AliESDMuonTrack& operator=(const AliESDMuonTrack& esdm);


 // Get and Set methods for data
  Double_t GetInverseBendingMomentum(void) const {return fInverseBendingMomentum;}
  void SetInverseBendingMomentum(Double_t InverseBendingMomentum) 
    {fInverseBendingMomentum = InverseBendingMomentum;}
  Double_t GetThetaX(void) const {return fThetaX;}
  void SetThetaX(Double_t ThetaX) {fThetaX = ThetaX;}
  Double_t GetThetaY(void) const {return fThetaY;}
  void SetThetaY(Double_t ThetaY) {fThetaY = ThetaY;}
  Double_t GetZ(void) const {return fZ;}
  void SetZ(Double_t Z) {fZ = Z;}
  Double_t GetBendingCoor(void) const {return fBendingCoor;}
  void SetBendingCoor(Double_t BendingCoor) {fBendingCoor = BendingCoor;}
  Double_t GetNonBendingCoor(void) const {return fNonBendingCoor;}
  void SetNonBendingCoor(Double_t NonBendingCoor) {fNonBendingCoor = NonBendingCoor;}
  Double_t GetChi2(void) const {return fChi2;}
  void SetChi2(Double_t Chi2) {fChi2 = Chi2;}
  UInt_t GetNHit(void) const {return fNHit;}
  void SetNHit(UInt_t NHit) {fNHit = NHit;}

  Bool_t GetMatchTrigger() const {return fMatchTrigger;}
  void SetMatchTrigger(Bool_t MatchTrigger) {fMatchTrigger = MatchTrigger;}
  Double_t GetChi2MatchTrigger() const {return fChi2MatchTrigger;}
  void SetChi2MatchTrigger(Double_t Chi2MatchTrigger) {fChi2MatchTrigger = Chi2MatchTrigger;}

protected:
  // tracking chamber
  Double_t fInverseBendingMomentum; // Inverse bending momentum (GeV/c ** -1) times the charge 
  Double_t fThetaX;           // Angle of track at vertex in X direction (rad)
  Double_t fThetaY;           // Angle of track at vertex in Y direction (rad)
  Double_t fZ;                // Z coordinate (cm)
  Double_t fBendingCoor;      // bending coordinate (cm)
  Double_t fNonBendingCoor;   // non bending coordinate (cm)
  Double_t fChi2;             // chi2 in the MUON track fit
  UInt_t   fNHit;              // number of hit in the track

  // trigger matching
  Bool_t   fMatchTrigger; // 1 if track matches with trigger track, 0 if not
  Double_t fChi2MatchTrigger; // chi2 of trigger/track matching 


  ClassDef(AliESDMuonTrack,2)  //MUON ESD track class 
};

#endif 
