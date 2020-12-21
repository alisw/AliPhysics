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
  enum poiTypes { kInvalid=-1,
                 kRP=0,
                 kPOI=1,
                 kPOI1=2,
                 kPOI2=3,
               };

  AliFlowTrackSimple();
  AliFlowTrackSimple(TParticle* p);
  AliFlowTrackSimple(const AliFlowTrackSimple& aTrack);
  AliFlowTrackSimple& operator=(const AliFlowTrackSimple& aTrack);
  virtual  ~AliFlowTrackSimple();
  virtual AliFlowTrackSimple* Clone(const char* option="") const;
  
  Bool_t  IsFolder() const {return kTRUE;};
  //  void Browse(TBrowser *b); 
  virtual void Print(Option_t* option = "") const;

  void Set(TParticle* p);

  Double_t Eta() const; 
  Double_t Pt()  const; 
  Double_t Phi() const;
  Double_t Weight() const; 
  Int_t Charge() const;
  Double_t Mass() const;
  Int_t ITStype() const;
  Int_t PID() const {return 0;}
  
  Bool_t InRPSelection() const; 
  Bool_t InPOISelection(Int_t poiType=1) const; 
  Bool_t IsPOItype(Int_t poiType) const;
  void SetPOItype(Int_t poiType, Bool_t b=kTRUE);
  Bool_t InSubevent(Int_t i) const;
  void TagRP(Bool_t b=kTRUE) {SetForRPSelection(b);} 
  void TagPOI(Bool_t b=kTRUE) {SetForPOISelection(b);} 
  void Tag(Int_t n, Bool_t b=kTRUE) {fPOItype.SetBitNumber(n,b);}
  Bool_t CheckTag(Int_t n) {return fPOItype.TestBitNumber(n);}
  void SetForSubevent(Int_t i); 
  void ResetPOItype() {fPOItype.ResetAllBits();}
  void ResetSubEventTags() {fSubEventBits.ResetAllBits();}
  Bool_t IsDead() const {return (fPOItype.CountBits()==0);}
      
  void SetEta(Double_t eta);
  void SetPt(Double_t pt); 
  void SetPhi(Double_t phi);
  void SetWeight(Double_t weight);
  void SetCharge(Int_t charge);
  void SetMass(Double_t mass);
  void SetITStype(Int_t val);
  void SetForRPSelection(Bool_t b=kTRUE); 
  void SetForPOISelection(Bool_t b=kTRUE); 
  virtual void Clear(Option_t* o="");
  
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

  const TBits* GetPOItype() const {return &fPOItype;}
  const TBits* GetFlowBits() const {return GetPOItype();}

  void  SetID(Int_t i) {fID=i;}
  Int_t GetID() const {return fID;}

  virtual Int_t GetNDaughters() const {return 0;}
  virtual void  AddDaughter(Int_t /*value*/) {}
  virtual Int_t GetIDDaughter(Int_t /*value*/) const {return 0;}
  virtual void SetDaughter(Int_t /*value*/, AliFlowTrackSimple* /*track*/) {}
  virtual AliFlowTrackSimple *GetDaughter(Int_t /*value*/) const {return NULL;}

 private:
  AliFlowTrackSimple(Double_t phi, Double_t eta, Double_t pt, Double_t weight, Int_t charge, Double_t mass=-1);
  Double_t fEta;         // eta
  Double_t fPt;          // pt
  Double_t fPhi;         // phi
  Double_t fTrackWeight; // weight
  Int_t    fCharge;      //charge
  Double_t fMass;        // mass
  TBits    fPOItype;     // bits to set if track is selected
  TBits    fSubEventBits;// bits to set if track is selected for a subevent
  Int_t    fID;          // Unique track ID, point back to the ESD track
  Int_t    fITStype;     // ITS hits identifier (test purpose only)

  ClassDef(AliFlowTrackSimple,2)                 // macro for rootcint

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
inline Double_t AliFlowTrackSimple::Mass() const { 
  return this->fMass; }
inline Int_t AliFlowTrackSimple::ITStype() const {
  return this->fITStype; }
//TBits
inline Bool_t AliFlowTrackSimple::InRPSelection() const { 
  return fPOItype.TestBitNumber(kRP); }
inline Bool_t AliFlowTrackSimple::InPOISelection(Int_t poiType) const { 
  return fPOItype.TestBitNumber(poiType); }
inline Bool_t AliFlowTrackSimple::IsPOItype(Int_t poiType) const {
  return fPOItype.TestBitNumber(poiType); }
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
inline void AliFlowTrackSimple::SetMass(Double_t val) {
  fMass = val; }
inline void AliFlowTrackSimple::SetITStype(Int_t val) {
  fITStype = val; }

  //TBits
inline void AliFlowTrackSimple::SetForRPSelection(Bool_t val) {
  fPOItype.SetBitNumber(kRP,val); }
inline void AliFlowTrackSimple::SetForPOISelection(Bool_t val) {
  fPOItype.SetBitNumber(kPOI,val); }
inline void AliFlowTrackSimple::SetForSubevent(Int_t i) {
  fSubEventBits.SetBitNumber(i,kTRUE); }

inline void AliFlowTrackSimple::SetPOItype(Int_t poiType, Bool_t b) {
  fPOItype.SetBitNumber(poiType,b); }

#endif

