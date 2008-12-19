/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

// AliFlowTrackSimpleCuts:
// A simple track cut class to the the AliFlowTrackSimple for basic
// kinematic cuts
// author: N. van der Kolk (kolk@nikhef.nl)

#ifndef AliFlowTrackSimpleCuts_H
#define AliFlowTrackSimpleCuts_H

class TObject;

class AliFlowTrackSimpleCuts : public TObject {

 public:
  AliFlowTrackSimpleCuts();
  AliFlowTrackSimpleCuts(const AliFlowTrackSimpleCuts& someCuts);
  AliFlowTrackSimpleCuts& operator=(const AliFlowTrackSimpleCuts& someCuts);
  virtual  ~AliFlowTrackSimpleCuts();
  
  //setters
  void SetPtMax(Double_t max)   {this->fPtMax = max; }
  void SetPtMin(Double_t min)   {this->fPtMin = min; }
  void SetEtaMax(Double_t max)  {this->fEtaMax = max; }
  void SetEtaMin(Double_t min)  {this->fEtaMin = min; }
  void SetPhiMax(Double_t max)  {this->fPhiMax = max; }
  void SetPhiMin(Double_t min)  {this->fPhiMin = min; }
  void SetPID(Int_t pid)        {this->fPID = pid; }
  
  //getters
  Double_t GetPtMax() const     {return this->fPtMax; }
  Double_t GetPtMin() const     {return this->fPtMin; }
  Double_t GetEtaMax() const    {return this->fEtaMax; }
  Double_t GetEtaMin() const    {return this->fEtaMin; }
  Double_t GetPhiMax() const    {return this->fPhiMax; }
  Double_t GetPhiMin() const    {return this->fPhiMin; }
  Int_t    GetPID() const       {return this->fPID; }

 private:
  Double_t fPtMax;
  Double_t fPtMin;
  Double_t fEtaMax;
  Double_t fEtaMin;
  Double_t fPhiMax;
  Double_t fPhiMin;
  Int_t    fPID;

  ClassDef(AliFlowTrackSimpleCuts,1)
};

#endif


