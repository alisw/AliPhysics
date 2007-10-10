#ifndef ALIESDMUONTRACK_H
#define ALIESDMUONTRACK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//  Class to describe the MUON tracks
//  in the Event Summary Data class
//  Author: G.Martinez


#include <TMath.h>
#include <TMatrixD.h>
#include <TDatabasePDG.h>

#include "AliVParticle.h"

class TLorentzVector;

class AliESDMuonTrack : public AliVParticle {
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

 // Get and Set methods for covariance matrix of data at first station
  void     GetCovariances(TMatrixD& cov) const;
  void     SetCovariances(const TMatrixD& cov);
  void     GetCovarianceXYZPxPyPz(Double_t cov[21]) const;

 // Get and Set methods for global tracking info
  Double_t GetChi2(void) const {return fChi2;}
  void     SetChi2(Double_t Chi2) {fChi2 = Chi2;}
  UChar_t  GetNHit(void) const {return fNHit;}
  void     SetNHit(UInt_t NHit) {fNHit = NHit;}

 // Get and Set methods for trigger matching
  Int_t    GetMatchTrigger() const;
  Double_t GetChi2MatchTrigger() const {return fChi2MatchTrigger;}
  void     SetChi2MatchTrigger(Double_t Chi2MatchTrigger) {fChi2MatchTrigger = Chi2MatchTrigger;}
  UShort_t GetHitsPatternInTrigCh() const {return fHitsPatternInTrigCh;}
  void     SetHitsPatternInTrigCh(UShort_t hitsPatternInTrigCh) {fHitsPatternInTrigCh = hitsPatternInTrigCh;}
  void     SetLocalTrigger(Int_t locTrig) { fLocalTrigger = locTrig; }
  Int_t    LoCircuit(void) const { return fLocalTrigger & 0xFF;       }
  Int_t    LoStripX(void) const  { return fLocalTrigger >>  8 & 0x1F; }
  Int_t    LoStripY(void) const  { return fLocalTrigger >> 13 & 0x0F; }
  Int_t    LoDev(void)    const  { return fLocalTrigger >> 17 & 0x1F; }
  Int_t    LoLpt(void)    const  { return fLocalTrigger >> 22 & 0x03; }
  Int_t    LoHpt(void)    const  { return fLocalTrigger >> 24 & 0x03; }
  
 // Get and Set methods for muon cluster map
  UInt_t   GetMuonClusterMap() const {return fMuonClusterMap;}
  void     SetMuonClusterMap(UInt_t muonClusterMap) {fMuonClusterMap = muonClusterMap;}
  void     AddInMuonClusterMap(Int_t chamber);
  Bool_t   IsInMuonClusterMap(Int_t chamber) const;
  
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
  
  // additional methods to comply with AliVParticle
  Double_t Pt() const { return TMath::Sqrt(Px()*Px() + Py()*Py()); }
  Double_t OneOverPt() const { return 1./Pt(); }
  Double_t Phi() const { return TMath::ATan2(Py(), Px()); }
  Double_t Theta() const { return TMath::ATan2(Pt(), Pz()); }
  Double_t E() const { return TMath::Sqrt(M()*M() + P()*P()); }
  Double_t M() const { return TDatabasePDG::Instance()->GetParticle("mu-")->Mass(); }
  Double_t Eta() const { return -TMath::Log(TMath::Tan(0.5 * Theta()));}
  Double_t Y() const { return TMath::ATanH(Pz()/E()); }
  Short_t  Charge() const { return (Short_t)TMath::Sign(1., GetInverseBendingMomentum()); }
  const Double_t *PID() const { return (Double_t*)0x0; }
    
  
protected:
 // parameters at vertex
  Double32_t fInverseBendingMomentum; ///< Inverse bending momentum (GeV/c ** -1) times the charge 
  Double32_t fThetaX;                 ///< Angle of track at vertex in X direction (rad)
  Double32_t fThetaY;                 ///< Angle of track at vertex in Y direction (rad)
  Double32_t fZ;                      ///< Z coordinate (cm)
  Double32_t fBendingCoor;            ///< bending coordinate (cm)
  Double32_t fNonBendingCoor;         ///< non bending coordinate (cm)
  
 // parameters at first tracking station
  Double32_t fInverseBendingMomentumUncorrected; ///< Inverse bending momentum (GeV/c ** -1) times the charge 
  Double32_t fThetaXUncorrected;                 ///< Angle of track at vertex in X direction (rad)
  Double32_t fThetaYUncorrected;                 ///< Angle of track at vertex in Y direction (rad)
  Double32_t fZUncorrected;                      ///< Z coordinate (cm)
  Double32_t fBendingCoorUncorrected;            ///< bending coordinate (cm)
  Double32_t fNonBendingCoorUncorrected;         ///< non bending coordinate (cm)
  
 /// reduced covariance matrix of UNCORRECTED track parameters, ordered as follow:      <pre>
 /// [0] =  <X,X>
 /// [1] =<X,ThetaX>  [2] =<ThetaX,ThetaX>
 /// [3] =  <X,Y>     [4] =  <Y,ThetaX>     [5] =  <Y,Y>
 /// [6] =<X,ThetaY>  [7] =<ThetaX,ThetaY>  [8] =<Y,ThetaY>  [9] =<ThetaY,ThetaY>
 /// [10]=<X,InvP_yz> [11]=<ThetaX,InvP_yz> [12]=<Y,InvP_yz> [13]=<ThetaY,InvP_yz> [14]=<InvP_yz,InvP_yz>  </pre>
  Double32_t fCovariances[15]; ///< \brief reduced covariance matrix of parameters AT FIRST CHAMBER
  
 // global tracking info
  Double32_t fChi2;                ///< chi2 in the MUON track fit
  Double32_t fChi2MatchTrigger;    ///< chi2 of trigger/track matching
  Int_t      fLocalTrigger;        ///< packed local trigger information

  UInt_t     fMuonClusterMap;      ///< Map of clusters in tracking chambers
  UShort_t   fHitsPatternInTrigCh; ///< Word containing info on the hits left in trigger chambers
  UChar_t    fNHit;                ///< number of hit in the track


  ClassDef(AliESDMuonTrack,8) // MUON ESD track class 
};

#endif 
