/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef ALIFLOWTRACKSIMPLE_H
#define ALIFLOWTRACKSIMPLE_H

#include "TNamed.h"
#include "TBits.h"
class TParticle;

// AliFlowTrackSimple:
// A simple track class to the the AliFlowEventSimple for flow analysis
// author: N. van der Kolk (kolk@nikhef.nl)
// mods: Mikolaj Krzewicki (mikolaj.krzewicki@cern.ch)

class AliFlowTrackSimple: public TNamed {

 public:
  AliFlowTrackSimple();
  AliFlowTrackSimple(const TParticle* p);
  AliFlowTrackSimple(const AliFlowTrackSimple& aTrack);
  AliFlowTrackSimple& operator=(const AliFlowTrackSimple& aTrack);
  virtual  ~AliFlowTrackSimple();
  
  Bool_t  IsFolder() const {return kTRUE;};
  //  void Browse(TBrowser *b); 
  //  void Print(Option_t* option = "") const;      //method to print stats
 

  Double_t Eta() const; 
  Double_t Pt()  const; 
  Double_t Phi() const; 
  TBits    GetFlowBits() const;
  Bool_t   InRPSelection() const; 
  Bool_t   InPOISelection() const; 
  Bool_t   InSubevent(Int_t i) const; 
      
  void SetEta(Double_t eta);
  void SetPt(Double_t pt); 
  void SetPhi(Double_t phi);
  void SetForRPSelection(Bool_t b=kTRUE); 
  void SetForPOISelection(Bool_t b=kTRUE); 
  void SetForSubevent(Int_t i); 
  
  void ResolutionPt(Double_t resolution);

  void AddFlow( Double_t flow, Double_t planeAngle);
    
 private:
  Double_t fEta;    // eta
  Double_t fPt;     // pt
  Double_t fPhi;    // phi
  TBits fFlowBits;  // bits to set if track is selected
  TBits fSubEventBits;  // bits to set if track is selected for a subevent
  

  ClassDef(AliFlowTrackSimple,1)                 // macro for rootcint

};

inline Double_t AliFlowTrackSimple::Eta() const { 
  return this->fEta; }
inline Double_t AliFlowTrackSimple::Pt() const {  
  //  cout << "Returned pt:" << fPt << endl; 
  return this->fPt;}
inline Double_t AliFlowTrackSimple::Phi() const { 
  return this->fPhi; }
//TBits
inline TBits AliFlowTrackSimple::GetFlowBits() const { 
  return this->fFlowBits; }
inline Bool_t AliFlowTrackSimple::InRPSelection() const { 
  return this->fFlowBits.TestBitNumber(0); }
inline Bool_t AliFlowTrackSimple::InPOISelection() const { 
  return this->fFlowBits.TestBitNumber(1); }
inline Bool_t AliFlowTrackSimple::InSubevent(Int_t i) const { 
  return this->fSubEventBits.TestBitNumber(i); }

inline void AliFlowTrackSimple::SetEta(Double_t val) {
  fEta = val; }
inline void AliFlowTrackSimple::SetPt(Double_t val) {
  fPt = val; }
  //  cout << "pt set to:" << fPt << endl;}
inline void AliFlowTrackSimple::SetPhi(Double_t val) {
  fPhi = val; }
//TBits
inline void AliFlowTrackSimple::SetForRPSelection(Bool_t val) {
  fFlowBits.SetBitNumber(0,val); }
inline void AliFlowTrackSimple::SetForPOISelection(Bool_t val) {
  fFlowBits.SetBitNumber(1,val); }
inline void AliFlowTrackSimple::SetForSubevent(Int_t i) {
  fSubEventBits.SetBitNumber(i,kTRUE); }

#endif

