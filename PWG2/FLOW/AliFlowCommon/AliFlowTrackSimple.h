/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliFlowTrackSimple_H
#define AliFlowTrackSimple_H

#include "TBits.h"

// AliFlowTrackSimple:
// A simple track class to the the AliFlowEventSimple for flow analysis
// author: N. van der Kolk (kolk@nikhef.nl)

class AliFlowTrackSimple: public TObject {

 public:
  AliFlowTrackSimple();
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
  //Bool_t   UseForIntegratedFlow() const;
  //Bool_t   UseForDifferentialFlow() const;
  Bool_t   InRPSelection() const; //will replace UseForIntegratedFlow()
  Bool_t   InPOISelection() const; //will replace UseForDifferentialFlow()
  
  void SetEta(Double_t eta);
  void SetPt(Double_t pt); 
  void SetPhi(Double_t phi);
  void SetForIntegratedFlow(Bool_t b);
  void SetForDifferentialFlow(Bool_t b);
  void SetForRPSelection(Bool_t b); //will replace SetForIntegratedFlow
  void SetForPOISelection(Bool_t b); //will replace SetForDifferentialFlow
  
 private:
  Double_t fEta;    // eta
  Double_t fPt;     // pt
  Double_t fPhi;    // phi
  TBits fFlowBits;  // bits to set if track is selected
  

  ClassDef(AliFlowTrackSimple,1)                 // macro for rootcint

};

inline Double_t AliFlowTrackSimple::Eta() const { 
  return this->fEta; }
inline Double_t AliFlowTrackSimple::Pt() const {  
  //  cout << "Returned pt:" << fPt << endl; 
  return this->fPt;}
inline Double_t AliFlowTrackSimple::Phi() const { 
  return this->fPhi; }
inline TBits AliFlowTrackSimple::GetFlowBits() const { 
  return this->fFlowBits; }
//inline Bool_t AliFlowTrackSimple::UseForIntegratedFlow() const { 
//  return this->fFlowBits.TestBitNumber(0); }
//inline Bool_t AliFlowTrackSimple::UseForDifferentialFlow() const { 
//  return this->fFlowBits.TestBitNumber(1); }
inline Bool_t AliFlowTrackSimple::InRPSelection() const { 
  return this->fFlowBits.TestBitNumber(0); }
inline Bool_t AliFlowTrackSimple::InPOISelection() const { 
  return this->fFlowBits.TestBitNumber(1); }

inline void AliFlowTrackSimple::SetEta(Double_t val) {
  fEta = val; }
inline void AliFlowTrackSimple::SetPt(Double_t val) {
  fPt = val; }
  //  cout << "pt set to:" << fPt << endl;}
inline void AliFlowTrackSimple::SetPhi(Double_t val) {
  fPhi = val; }
//inline void AliFlowTrackSimple::SetForIntegratedFlow(Bool_t val) { //replace
//  fFlowBits.SetBitNumber(0,val); }
//inline void AliFlowTrackSimple::SetForDifferentialFlow(Bool_t val) { //replace
//  fFlowBits.SetBitNumber(1,val); }
inline void AliFlowTrackSimple::SetForRPSelection(Bool_t val) {
  fFlowBits.SetBitNumber(0,val); }
inline void AliFlowTrackSimple::SetForPOISelection(Bool_t val) {
  fFlowBits.SetBitNumber(1,val); }
#endif

