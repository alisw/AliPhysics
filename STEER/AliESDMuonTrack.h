#ifndef ALIESDMUONTRACK_H
#define ALIESDMUONTRACK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//  Class to describe the MUON tracks
//  in the Event Summary Data class
//  Author: G.Martinez


#include "TObject.h"

class TLorentzVector;

class AliESDMuonTrack : public TObject {
public:
  AliESDMuonTrack(); //Constructor
  virtual ~AliESDMuonTrack(){} // Destructor
  AliESDMuonTrack(const AliESDMuonTrack& esdm);
  AliESDMuonTrack& operator=(const AliESDMuonTrack& esdm);


 // Get and Set methods for data at vertex
  Double_t GetInverseBendingMomentum(void) const {return fInverseBendingMomentum;}
  void     SetInverseBendingMomentum(Double_t InverseBendingMomentum) 
		{fInverseBendingMomentum = InverseBendingMomentum;}
  Double_t GetThetaX(void) const {return fThetaX;}
  void     SetThetaX(Double_t ThetaX) {fThetaX = ThetaX;}
  Double_t GetThetaY(void) const {return fThetaY;}
  void     SetThetaY(Double_t ThetaY) {fThetaY = ThetaY;}
  Double_t GetZ(void) const {return fZ;}
  void     SetZ(Double_t Z) {fZ = Z;}
  Double_t GetBendingCoor(void) const {return fBendingCoor;}
  void     SetBendingCoor(Double_t BendingCoor) {fBendingCoor = BendingCoor;}
  Double_t GetNonBendingCoor(void) const {return fNonBendingCoor;}
  void     SetNonBendingCoor(Double_t NonBendingCoor) {fNonBendingCoor = NonBendingCoor;}
  
 // Get and Set methods for data at first station
  Double_t GetInverseBendingMomentumUncorrected(void) const {return fInverseBendingMomentumUncorrected;}
  void     SetInverseBendingMomentumUncorrected(Double_t InverseBendingMomentum) 
		{fInverseBendingMomentumUncorrected = InverseBendingMomentum;}
  Double_t GetThetaXUncorrected(void) const {return fThetaXUncorrected;}
  void     SetThetaXUncorrected(Double_t ThetaX) {fThetaXUncorrected = ThetaX;}
  Double_t GetThetaYUncorrected(void) const {return fThetaYUncorrected;}
  void     SetThetaYUncorrected(Double_t ThetaY) {fThetaYUncorrected = ThetaY;}
  Double_t GetZUncorrected(void) const {return fZUncorrected;}
  void     SetZUncorrected(Double_t Z) {fZUncorrected = Z;}
  Double_t GetBendingCoorUncorrected(void) const {return fBendingCoorUncorrected;}
  void     SetBendingCoorUncorrected(Double_t BendingCoor) {fBendingCoorUncorrected = BendingCoor;}
  Double_t GetNonBendingCoorUncorrected(void) const {return fNonBendingCoorUncorrected;}
  void     SetNonBendingCoorUncorrected(Double_t NonBendingCoor) {fNonBendingCoorUncorrected = NonBendingCoor;}
  
 // Get and Set methods for global tracking info
  Double_t GetChi2(void) const {return fChi2;}
  void     SetChi2(Double_t Chi2) {fChi2 = Chi2;}
  UInt_t   GetNHit(void) const {return fNHit;}
  void     SetNHit(UInt_t NHit) {fNHit = NHit;}

 // Get and Set methods for trigger matching
  Int_t    GetMatchTrigger() const;
  Double_t GetChi2MatchTrigger() const {return fChi2MatchTrigger;}
  void     SetChi2MatchTrigger(Double_t Chi2MatchTrigger) {fChi2MatchTrigger = Chi2MatchTrigger;}
  UShort_t GetHitsPatternInTrigCh() const {return fHitsPatternInTrigCh;}
  void     SetHitsPatternInTrigCh(UShort_t hitsPatternInTrigCh) {fHitsPatternInTrigCh = hitsPatternInTrigCh;}
  void     SetLocalTrigger(Int_t locTrig) { fLocalTrigger = locTrig; }
  Int_t    LoCircuit(void) const
  { Int_t circ = fLocalTrigger & 0xFF; return (circ == 234) ? -1 : circ; }
  Int_t    LoStripX(void) const  { return fLocalTrigger >>  8 & 0x1F; }
  Int_t    LoStripY(void) const  { return fLocalTrigger >> 13 & 0x0F; }
  Int_t    LoDev(void)    const  { return fLocalTrigger >> 17 & 0x1F; }
  Int_t    LoLpt(void)    const  { return fLocalTrigger >> 22 & 0x03; }
  Int_t    LoHpt(void)    const  { return fLocalTrigger >> 24 & 0x03; }
  
 // Methods to compute track momentum
  Double_t Px() const;
  Double_t Py() const;
  Double_t Pz() const;
  Double_t P() const;
  void     LorentzP(TLorentzVector& vP) const;
  Double_t PxUncorrected() const;
  Double_t PyUncorrected() const;
  Double_t PzUncorrected() const;
  Double_t PUncorrected() const;
  void     LorentzPUncorrected(TLorentzVector& vP) const;
  
  
protected:
 // parameters at vertex
  Double32_t fInverseBendingMomentum; // Inverse bending momentum (GeV/c ** -1) times the charge 
  Double32_t fThetaX;		    // Angle of track at vertex in X direction (rad)
  Double32_t fThetaY;		    // Angle of track at vertex in Y direction (rad)
  Double32_t fZ;			    // Z coordinate (cm)
  Double32_t fBendingCoor;	    // bending coordinate (cm)
  Double32_t fNonBendingCoor;	    // non bending coordinate (cm)
  
 // parameters at first tracking station
  Double32_t fInverseBendingMomentumUncorrected; // Inverse bending momentum (GeV/c ** -1) times the charge 
  Double32_t fThetaXUncorrected;		       // Angle of track at vertex in X direction (rad)
  Double32_t fThetaYUncorrected;		       // Angle of track at vertex in Y direction (rad)
  Double32_t fZUncorrected;		       // Z coordinate (cm)
  Double32_t fBendingCoorUncorrected;	       // bending coordinate (cm)
  Double32_t fNonBendingCoorUncorrected;	       // non bending coordinate (cm)
  
 // global tracking info
  Double32_t fChi2; // chi2 in the MUON track fit
  UInt_t   fNHit; // number of hit in the track

  Int_t fLocalTrigger;    ///< packed local trigger information
  
  Double32_t fChi2MatchTrigger; // chi2 of trigger/track matching
  
  UShort_t fHitsPatternInTrigCh; ///< Word containing info on the hits left in trigger chambers


  ClassDef(AliESDMuonTrack,5)  //MUON ESD track class 
};

#endif 
