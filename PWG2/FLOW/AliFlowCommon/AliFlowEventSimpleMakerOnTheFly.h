/* 
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 
 * See cxx source for full Copyright notice 
 * $Id$ 
 */

/************************************ 
 * Create an event and perform full *
 * flow analysis 'on the fly'.      * 
 *                                  * 
 * author: Ante Bilandzic           * 
 *         (abilandzic@gmail.com)   *
 ************************************/ 

#ifndef ALIFLOWEVENTSIMPLEMAKERONTHEFLY_H
#define ALIFLOWEVENTSIMPLEMAKERONTHEFLY_H

class TF1;
class TRandom3;

#include "AliFlowEventSimple.h" // needed as include
#include "AliFlowTrackSimpleCuts.h"
    
class AliFlowEventSimpleMakerOnTheFly{
 public:
  AliFlowEventSimpleMakerOnTheFly(UInt_t uiSeed = 0); // constructor
  virtual ~AliFlowEventSimpleMakerOnTheFly(); // destructor
  virtual void Init();   
  Bool_t AcceptOrNot(AliFlowTrackSimple *pTrack);  
  AliFlowEventSimple* CreateEventOnTheFly(AliFlowTrackSimpleCuts *cutsRP, AliFlowTrackSimpleCuts *cutsPOI); 
  // Setters and getters:
  void SetMinMult(Int_t iMinMult) {this->fMinMult = iMinMult;}
  Int_t GetMinMult() const {return this->fMinMult;} 
  void SetMaxMult(Int_t iMaxMult) {this->fMaxMult = iMaxMult;}
  Int_t GetMaxMult() const {return this->fMaxMult;} 
  void SetMass(Double_t dMass) {this->fMass = dMass;}
  Double_t GetMass() const {return this->fMass;} 
  void SetTemperature(Double_t dT) {this->fTemperature = dT;}
  Double_t GetTemperature() const {return this->fTemperature;} 
  void SetV1(Double_t dV1) {this->fV1 = dV1;}
  Double_t GetV1() const {return this->fV1;} 
  void SetV2(Double_t dV2) {this->fV2 = dV2;}
  Double_t GetV2() const {return this->fV2;} 
  void SetV3(Double_t dV3) {this->fV3 = dV3;}
  Double_t GetV3() const {return this->fV3;} 
  void SetV4(Double_t dV4) {this->fV4 = dV4;}
  Double_t GetV4() const {return this->fV4;} 
  void SetUniformFluctuationsV2(Bool_t b) {this->fUniformFluctuationsV2 = b;}
  Bool_t GetUniformFluctuationsV2() const {return this->fUniformFluctuationsV2;} 
  void SetMinV2(Double_t dMinV2) {this->fMinV2 = dMinV2;}
  Double_t GetMinV2() const {return this->fMinV2;} 
  void SetMaxV2(Double_t dMaxV2) {this->fMaxV2 = dMaxV2;}
  Double_t GetMaxV2() const {return this->fMaxV2;} 
  void SetPtDependentV2(Bool_t b) {this->fPtDependentV2 = b;}
  Bool_t GetPtDependentV2() const {return this->fPtDependentV2;} 
  void SetV2vsPtCutOff(Double_t dV2vsPtCutOff) {this->fV2vsPtCutOff = dV2vsPtCutOff;}
  Double_t GetV2vsPtCutOff() const {return this->fV2vsPtCutOff;} 
  void SetV2vsPtMax(Double_t dV2vsPtMax) {this->fV2vsPtMax = dV2vsPtMax;}
  Double_t GetV2vsPtMax() const {return this->fV2vsPtMax;} 
  void SetSubeventEtaRange(Double_t minA, Double_t maxA, Double_t minB, Double_t maxB) 
  {this->fEtaMinA = minA;this->fEtaMaxA = maxA;this->fEtaMinB = minB;this->fEtaMaxB = maxB;};
  void SetNTimes(Int_t nt) {this->fNTimes = nt;}
  Int_t GetNTimes() const {return this->fNTimes;} 
  void SetUniformAcceptance(Bool_t ua) {this->fUniformAcceptance = ua;}
  Bool_t GetUniformAcceptance() const {return this->fUniformAcceptance;} 
  void SetFirstSectorPhiMin(Double_t dPhiMin1) {this->fPhiMin1 = dPhiMin1;}
  Double_t GetFirstSectorPhiMin() const {return this->fPhiMin1;} 
  void SetFirstSectorPhiMax(Double_t dPhiMax1) {this->fPhiMax1 = dPhiMax1;}
  Double_t GetFirstSectorPhiMax() const {return this->fPhiMax1;}  
  void SetFirstSectorProbability(Double_t dProbability1) {this->fProbability1 = dProbability1;}
  Double_t GetFirstSectorProbability() const {return this->fProbability1;}    
  void SetSecondSectorPhiMin(Double_t dPhiMin2) {this->fPhiMin2 = dPhiMin2;}
  Double_t GetSecondSectorPhiMin() const {return this->fPhiMin2;} 
  void SetSecondSectorPhiMax(Double_t dPhiMax2) {this->fPhiMax2 = dPhiMax2;}
  Double_t GetSecondSectorPhiMax() const {return this->fPhiMax2;}
  void SetSecondSectorProbability(Double_t dProbability2) {this->fProbability2 = dProbability2;}
  Double_t GetSecondSectorProbability() const {return this->fProbability2;}       
 private:
  AliFlowEventSimpleMakerOnTheFly(const AliFlowEventSimpleMakerOnTheFly& anAnalysis); // copy constructor
  AliFlowEventSimpleMakerOnTheFly& operator=(const AliFlowEventSimpleMakerOnTheFly& anAnalysis); // assignment operator
  Int_t fCount; // count number of events 
  Int_t fMinMult; // uniformly sampled multiplicity is >= iMinMult
  Int_t fMaxMult; // uniformly sampled multiplicity is < iMaxMult
  TF1 *fPtSpectra; // transverse momentum distribution (pt is sampled from hardwired Boltzmann distribution)
  Double_t fMass; // mass in pt distribution (hardwired is Boltzmann pt distribution)
  Double_t fTemperature; // "temperature" in pt distribution (hardwired is Boltzmann pt distribution)   
  TF1 *fPhiDistribution; // azimuthal distribution (phi is sampled from hardwired Fourier-like distribution)
  Double_t fV1; // harmonic v1
  Double_t fV2; // harmonic v2
  Double_t fV3; // harmonic v3
  Double_t fV4; // harmonic v4
  Bool_t fUniformFluctuationsV2; // v2 is sampled uniformly for each event and for all particles from [fMinV2,fMaxV2] 
  Double_t fMinV2; // if v2 is sampled uniformly for each event, this is lower boundary on its value  
  Double_t fMaxV2; // if v2 is sampled uniformly for each event, this is upper boundary on its value
  Bool_t fPtDependentV2; // v2 is pt-dependent
  Double_t fV2vsPtCutOff; // if v2 is pt-dependent: for v2 < fV2vsPtCutOff v2 is growing linearly, otherwise v2 = fV2vsPtMax
  Double_t fV2vsPtMax; // if v2 is pt-dependent: v2 = fV2vsPtMax for v2 >= fV2vsPtCutOff 
  Double_t fEtaMinA; // minimum eta of subevent A
  Double_t fEtaMaxA; // maximum eta of subevent A
  Double_t fEtaMinB; // minimum eta of subevent B
  Double_t fEtaMaxB; // maximum eta of subevent B 
  Int_t fNTimes; // number of times to use the same particle in the analysis (simulating nonflow)
  Bool_t fUniformAcceptance; // detector has uniform azimuthal acceptance or not
  Double_t fPhiMin1; // first sector with non-uniform acceptance starts at azimuth fPhiMin1
  Double_t fPhiMax1; // first sector with non-uniform acceptance ends at azimuth fPhiMax1
  Double_t fProbability1; // particles emitted in fPhiMin1 < phi < fPhiMax1 are taken with probability fProbability1 
  Double_t fPhiMin2; // second sector with non-uniform acceptance starts at azimuth fPhiMin2
  Double_t fPhiMax2; // second sector with non-uniform acceptance ends at azimuth fPhiMax2
  Double_t fProbability2; // particles emitted in fPhiMin2 < phi < fPhiMax2 are taken with probability fProbability2
  Double_t fPi; // pi

  ClassDef(AliFlowEventSimpleMakerOnTheFly,1) // macro for rootcint
};
 
#endif



