/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

// AliFlowTrackSimpleCuts:
// A simple track cut class to the the AliFlowTrackSimple for basic
// kinematic cuts
// author: N. van der Kolk (kolk@nikhef.nl)
// mods: Mikolaj Krzewicki (mikolaj.krzewicki@cern.ch), Redmer Bertens (rbertens@cern.ch)

#ifndef ALIFLOWTRACKSIMPLECUTS_H
#define ALIFLOWTRACKSIMPLECUTS_H

#include "TNamed.h"

class AliFlowTrackSimple;
class TParticle;
class TH2;

class AliFlowTrackSimpleCuts : public TNamed {

 public:
  AliFlowTrackSimpleCuts(const char* name="");
  //AliFlowTrackSimpleCuts(const AliFlowTrackSimpleCuts& someCuts);
  //AliFlowTrackSimpleCuts& operator=(const AliFlowTrackSimpleCuts& someCuts);
  virtual  ~AliFlowTrackSimpleCuts() {}
  
  //setters
  void SetPtMax(Double_t max)   {this->fPtMax = max; fCutPt=kTRUE; }
  void SetPtMin(Double_t min)   {this->fPtMin = min; fCutPt=kTRUE;  }
  void SetEtaMax(Double_t max)  {this->fEtaMax = max; fCutEta=kTRUE; }
  void SetEtaMin(Double_t min)  {this->fEtaMin = min; fCutEta=kTRUE; }
  void SetEtaGap(Double_t min, Double_t max)
        {fEtaGapMin = min, fEtaGapMax = max, fCutEtaGap = kTRUE; }
  void SetPhiMax(Double_t max)  {this->fPhiMax = max; fCutPhi=kTRUE; }
  void SetPhiMin(Double_t min)  {this->fPhiMin = min; fCutPhi=kTRUE; }
  void SetPID(Int_t pid)        {this->fPID = pid; fCutPID=kTRUE; }
  void SetCharge(Int_t c)       {this->fCharge = c; fCutCharge=kTRUE; }
  void SetMassMax(Double_t max) {this->fMassMax = max; fCutMass=kTRUE; }
  void SetMassMin(Double_t min) {this->fMassMin = min; fCutMass=kTRUE; }
  void SetEtaPhiEff(TH2* eff)   {this->fEtaPhiEff = eff; fCutEtaPhiEff=kTRUE; }
  
  //getters
  Double_t GetPtMax() const     {return this->fPtMax; }
  Double_t GetPtMin() const     {return this->fPtMin; }
  Double_t GetEtaMax() const    {return this->fEtaMax; }
  Double_t GetEtaMin() const    {return this->fEtaMin; }
  Double_t GetPhiMax() const    {return this->fPhiMax; }
  Double_t GetEtaGapMin() const {return this->fEtaGapMin; }
  Double_t GetEtaGapMax() const {return this->fEtaGapMax; }
  Double_t GetPhiMin() const    {return this->fPhiMin; }
  Int_t    GetPID() const       {return this->fPID; }
  Int_t    GetCharge() const    {return this->fCharge; }
  Double_t GetMassMax() const   {return this->fMassMax; }
  Double_t GetMassMin() const   {return this->fMassMin; }
  
  //simple method to check if the simple track passes the simple cuts:
  Bool_t PassesCuts(const AliFlowTrackSimple *track) const;
  Bool_t PassesCuts(TParticle* p) const;

  Int_t GetPOItype() const {return fPOItype;}
  void SetPOItype(Int_t t) {fPOItype=t;}

  virtual Bool_t IsSelected(TObject* obj, Int_t id=-1);

 protected:
  Bool_t   fCutPt; //cut on pt?
  Double_t fPtMax; //max pt
  Double_t fPtMin; //min pt
  Bool_t   fCutEta; //cut on eta?
  Double_t fEtaMax; //max eta
  Double_t fEtaMin; //min eta
  Bool_t   fCutEtaGap; //have eta gap (i.e. discontinuous rp selection in eta?)
  Double_t fEtaGapMax; // max eta gap
  Double_t fEtaGapMin; // min eta gap
  Bool_t   fCutPhi; //cut on phi?
  Double_t fPhiMax; //max phi
  Double_t fPhiMin; //min phi
  Bool_t   fCutPID; //cut on pid?
  Int_t    fPID;    //pid
  Bool_t   fCutCharge; //cut on charge?
  Int_t    fCharge;    //charge
  Bool_t   fCutMass; // cut on mass?
  Double_t fMassMax; //max mass
  Double_t fMassMin; //min mass
  TH2*     fEtaPhiEff; //eta, phi efficiency map
  Bool_t   fCutEtaPhiEff; //apply eta, phi efficiency map?
  Int_t    fPOItype; //which poi type do we produce? (RP=0, POI=1,2,3,4,5...)
  
  // not implemented: copy and assignment operator
  AliFlowTrackSimpleCuts(const AliFlowTrackSimpleCuts&);        
  AliFlowTrackSimpleCuts& operator=(const AliFlowTrackSimpleCuts&);
  
  ClassDef(AliFlowTrackSimpleCuts,4)
};

#endif


