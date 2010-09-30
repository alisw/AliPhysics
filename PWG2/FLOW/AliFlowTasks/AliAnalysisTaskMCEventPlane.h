/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef ALIANALYSISTASKMCEVENTPLANE_H
#define ALIANALYSISTASKMCEVENTPLANE_H

// AliAnalysisTaskMCEventPlane:
// analysis task for 
// Monte Carlo Event Plane
// Author: 
//        Naomi van der Kolk (kolk@nikhef.nl)

class AliFlowEventSimple;
class AliFlowAnalysisWithMCEventPlane;
class TList;

#include "TString.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskMCEventPlane : public AliAnalysisTaskSE {
 public:

  AliAnalysisTaskMCEventPlane();
  AliAnalysisTaskMCEventPlane(const char *name);
  virtual ~AliAnalysisTaskMCEventPlane();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  void SetHarmonic(Int_t const harmonic) {this->fHarmonic = harmonic;};
  Int_t GetHarmonic() const {return this->fHarmonic;};
  
  // Objects needed for mixed harmonics study:
  void SetEvaluateMixedHarmonics(Bool_t const emh) {this->fEvaluateMixedHarmonics = emh;};
  Bool_t GetEvalauteMixedHarmonics() const {return this->fEvaluateMixedHarmonics;};
  void SetNinCorrelator(Int_t const n) {this->fNinCorrelator = n;};
  Int_t GetNinCorrelator() const {return this->fNinCorrelator;};     
  void SetMinCorrelator(Int_t const m) {this->fMinCorrelator = m;};
  Int_t GetMinCorrelator() const {return this->fMinCorrelator;};     
  void SetXinPairAngle(Double_t const xipa) {this->fXinPairAngle = xipa;};
  Double_t GetXinPairAngle() const {return this->fXinPairAngle;};       
  void SetnBinsMult(Int_t const nbm) {this->fnBinsMult = nbm;};
  Int_t GetnBinsMult() const {return this->fnBinsMult;};  
  void SetMinMult(Double_t const minm) {this->fMinMult = minm;};
  Double_t GetMinMult() const {return this->fMinMult;};
  void SetMaxMult(Double_t const maxm) {this->fMaxMult = maxm;};
  Double_t GetMaxMult() const {return this->fMaxMult;};     

 private:
 
  AliAnalysisTaskMCEventPlane(const AliAnalysisTaskMCEventPlane& aAnalysis);
  AliAnalysisTaskMCEventPlane& operator=(const AliAnalysisTaskMCEventPlane& aAnalysis);
  
  AliFlowEventSimple*               fEvent;      //input event
  AliFlowAnalysisWithMCEventPlane*  fMc;         // MC EP analysis object
  TList*                            fListHistos; // collection of output
  Int_t                             fHarmonic;   // harmonic 
  // Objects needed for mixed harmonics study:
  Bool_t fEvaluateMixedHarmonics; // evaluate and store objects relevant for mixed harmonics
  Int_t fnBinsMult; // number of multiplicity bins for mixed harmonics analysis versus multiplicity  
  Double_t fMinMult; // minimal multiplicity for mixed harmonics analysis versus multiplicity  
  Double_t fMaxMult; // maximal multiplicity for mixed harmonics analysis versus multiplicity    
  Int_t fNinCorrelator; // n in <cos[m*phi_{pair}-n*RP]> and <sin[m*phi_{pair}-n*RP]>, where phi_{pair} = x*phi1+(1-x)*phi2
  Int_t fMinCorrelator; // m in <cos[m*phi_{pair}-n*RP]> and <sin[m*phi_{pair}-n*RP]>, where phi_{pair} = x*phi1+(1-x)*phi2   
  Double_t fXinPairAngle; // x in definition phi_{pair} = x*phi1+(1-x)*phi2

  ClassDef(AliAnalysisTaskMCEventPlane, 1); // AliAnalysisTaskMCEventPlane class object
};

#endif

