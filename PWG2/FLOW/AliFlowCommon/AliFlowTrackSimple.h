/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef ALIFLOWTRACKSIMPLE_H
#define ALIFLOWTRACKSIMPLE_H

#include "TObject.h"
#include "TBits.h"
class TParticle;

// AliFlowTrackSimple:
// A simple track class to the the AliFlowEventSimple for flow analysis
// author: N. van der Kolk (kolk@nikhef.nl)
// mods: Mikolaj Krzewicki (mikolaj.krzewicki@cern.ch)

class AliFlowTrackSimple: public TObject {

public:
  AliFlowTrackSimple();
  AliFlowTrackSimple(TParticle* p);
  AliFlowTrackSimple(const AliFlowTrackSimple& aTrack);
  AliFlowTrackSimple& operator=(const AliFlowTrackSimple& aTrack);
  virtual  ~AliFlowTrackSimple();
  virtual AliFlowTrackSimple* Clone(const char* option="") const;
  
  Bool_t  IsFolder() const {return kTRUE;};
  //  void Browse(TBrowser *b); 
  virtual void Print(Option_t* option = "") const;

  Double_t Eta() const; 
  Double_t Pt()  const; 
  Double_t Phi() const;
  Double_t Weight() const; 
  Int_t Charge() const;
  Int_t PID() const {return 0;}
  

  Bool_t InRPSelection() const; 
  Bool_t InPOISelection() const; 
  Bool_t InSubevent(Int_t i) const;
  Bool_t IsDead() const {return (fFlowBits.CountBits()==0);}
      
  void SetEta(Double_t eta);
  void SetPt(Double_t pt); 
  void SetPhi(Double_t phi);
  void SetWeight(Double_t weight);
  void SetCharge(Int_t charge);
  void SetForRPSelection(Bool_t b=kTRUE); 
  void SetForPOISelection(Bool_t b=kTRUE); 
  void TagRP(Bool_t b=kTRUE) {SetForRPSelection(b);} 
  void TagPOI(Bool_t b=kTRUE) {SetForPOISelection(b);} 
  void SetForSubevent(Int_t i); 
  void ResetFlowTags() {fFlowBits.ResetAllBits();}
  void ResetSubEventTags() {fSubEventBits.ResetAllBits();}
  
  void ResolutionPt(Double_t resolution);

  void AddV1( Double_t v1,
              Double_t reactionPlaneAngle,
              Double_t precision,
              Int_t maxNumberOfIterations=100 );
  void AddV2( Double_t v2,
              Double_t reactionPlaneAngle,
              Double_t precision,
              Int_t maxNumberOfIterations=100 );
  void AddV3( Double_t v3,
              Double_t reactionPlaneAngle,
              Double_t precision,
              Int_t maxNumberOfIterations=100 );
  void AddV4( Double_t v4,
              Double_t reactionPlaneAngle,
              Double_t precision,
              Int_t maxNumberOfIterations=100 );
  void AddV5( Double_t v5,
              Double_t reactionPlaneAngle,
              Double_t precision,
              Int_t maxNumberOfIterations=100 );
  void AddFlow( Double_t v1,
                Double_t v2,
                Double_t v3,
                Double_t v4,
                Double_t reactionPlaneAngle,
                Double_t precision,
                Int_t maxNumberOfIterations=100 );
  void AddFlow( Double_t v1,
                Double_t v2,
                Double_t v3,
                Double_t v4,
                Double_t v5,
                Double_t reactionPlaneAngle,
                Double_t precision,
                Int_t maxNumberOfIterations=100 );
  void AddFlow( Double_t v1,
                Double_t v2,
                Double_t v3,
                Double_t v4,
                Double_t v5,
                Double_t rp1,
                Double_t rp2,
                Double_t rp3,
                Double_t rp4,
                Double_t rp5,
                Double_t precision,
                Int_t maxNumberOfIterations=100 );

  const TBits* GetFlowBits() const {return &fFlowBits;}

  void  SetID(Int_t i) {fID=i;}
  Int_t GetID() {return fID;}

 private:
  AliFlowTrackSimple(Double_t phi, Double_t eta, Double_t pt, Double_t weight, Int_t charge);
  Double_t fEta;         // eta
  Double_t fPt;          // pt
  Double_t fPhi;         // phi
  Double_t fTrackWeight; // weight
  Int_t fCharge;         //charge
  TBits    fFlowBits;    // bits to set if track is selected
  TBits    fSubEventBits;// bits to set if track is selected for a subevent
  Int_t    fID;          // Unique track ID, point back to the ESD track

  ClassDef(AliFlowTrackSimple,1)                 // macro for rootcint

};

//Getters
inline Double_t AliFlowTrackSimple::Eta() const { 
  return this->fEta; }
inline Double_t AliFlowTrackSimple::Pt() const {  
  return this->fPt;}
inline Double_t AliFlowTrackSimple::Phi() const { 
  return this->fPhi; }
inline Double_t AliFlowTrackSimple::Weight() const { 
  return this->fTrackWeight; }
inline Int_t AliFlowTrackSimple::Charge() const { 
  return this->fCharge; }
//TBits
inline Bool_t AliFlowTrackSimple::InRPSelection() const { 
  return this->fFlowBits.TestBitNumber(0); }
inline Bool_t AliFlowTrackSimple::InPOISelection() const { 
  return this->fFlowBits.TestBitNumber(1); }
inline Bool_t AliFlowTrackSimple::InSubevent(Int_t i) const { 
  return this->fSubEventBits.TestBitNumber(i); }

//Setters
inline void AliFlowTrackSimple::SetEta(Double_t val) {
  fEta = val; }
inline void AliFlowTrackSimple::SetPt(Double_t val) {
  fPt = val; }
inline void AliFlowTrackSimple::SetPhi(Double_t val) {
  fPhi = val; }
inline void AliFlowTrackSimple::SetWeight(Double_t val) {
  fTrackWeight = val; }
inline void AliFlowTrackSimple::SetCharge(Int_t val) {
  fCharge = val; }
//TBits
inline void AliFlowTrackSimple::SetForRPSelection(Bool_t val) {
  fFlowBits.SetBitNumber(0,val); }
inline void AliFlowTrackSimple::SetForPOISelection(Bool_t val) {
  fFlowBits.SetBitNumber(1,val); }
inline void AliFlowTrackSimple::SetForSubevent(Int_t i) {
  fSubEventBits.SetBitNumber(i,kTRUE); }

#endif

