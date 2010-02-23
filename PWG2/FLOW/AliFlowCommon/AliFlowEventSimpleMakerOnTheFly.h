/* 
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 
 * See cxx source for full Copyright notice 
 * $Id$ 
 */

/********************************** 
 * create an event and perform    *
 * flow analysis 'on the fly'     * 
 *                                * 
 * authors: Raimond Snellings     *
 *           (snelling@nikhef.nl) * 
 *          Ante Bilandzic        * 
 *           (anteb@nikhef.nl)    *
 *********************************/ 

#ifndef ALIFLOWEVENTSIMPLEMAKERONTHEFLY_H
#define ALIFLOWEVENTSIMPLEMAKERONTHEFLY_H

class TF1;
class TRandom3;

#include "AliFlowEventSimple.h"  //needed as include
    
class AliFlowEventSimpleMakerOnTheFly {

 public:

  AliFlowEventSimpleMakerOnTheFly(UInt_t);    // constructor
  virtual ~AliFlowEventSimpleMakerOnTheFly(); // destructor

  virtual void Init(); 
  
  AliFlowEventSimple* CreateEventOnTheFly();  // create an event on the fly
 
    
  //                        *****************************
  //                        **** SETTERS AND GETTERS ****
  //                        *****************************
  //................................................................................................
  // setters and getters for global parameters:
  void SetMultDistrOfRPsIsGauss(Bool_t const mdorig) {this->fMultDistrOfRPsIsGauss = mdorig;};
  Bool_t GetMultDistrOfRPsIsGauss() const {return this->fMultDistrOfRPsIsGauss;};
  
  void SetMultiplicityOfRP(Int_t multRP) {this->fMultiplicityOfRP = multRP;}
  Int_t GetMultiplicityOfRP() const {return this->fMultiplicityOfRP;} 
  
  void SetMultiplicitySpreadOfRP(Double_t multSpreadRP) {this->fMultiplicitySpreadOfRP = multSpreadRP;}
  Double_t GetMultiplicitySpreadOfRP() const {return this->fMultiplicitySpreadOfRP;} 
  
  void SetMinMultOfRP(Int_t minmr) {this->fMinMultOfRP = minmr;}
  Int_t GetMinMultOfRP() const {return this->fMinMultOfRP;} 
  
  void SetMaxMultOfRP(Int_t maxmr) {this->fMaxMultOfRP = maxmr;}
  Int_t GetMaxMultOfRP() const {return this->fMaxMultOfRP;} 
  
  void SetTemperatureOfRP(Double_t temperatureRP) {this->fTemperatureOfRP = temperatureRP;}
  Double_t GetTemperatureOfRP() const {return this->fTemperatureOfRP;} 
  
  void SetPtDependentHarmonics(Bool_t const pdh) {this->fPtDependentHarmonics = pdh;};
  Bool_t GetPtDependentHarmonics() const {return this->fPtDependentHarmonics;};
  
  void SetEtaDependentHarmonics(Bool_t const edh) {this->fEtaDependentHarmonics = edh;};
  Bool_t GetEtaDependentHarmonics() const {return this->fEtaDependentHarmonics;};
  
  // constant harmonics:  
  void SetV1RP(Double_t dV1RP) {this->fV1RP = dV1RP;}
  Double_t GetV1RP() const {return this->fV1RP;} 
  
  void SetV1SpreadRP(Double_t dV1SpreadRP) {this->fV1SpreadRP = dV1SpreadRP;}
  Double_t GetV1SpreadRP() const {return this->fV1SpreadRP;} 
  
  void SetConstantV2IsSampledFromGauss(Bool_t const cV2isfg) {this->fConstantV2IsSampledFromGauss = cV2isfg;};
  Bool_t GetConstantV2IsSampledFromGauss() const {return this->fConstantV2IsSampledFromGauss;};
  
  void SetV2RP(Double_t dV2RP) {this->fV2RP = dV2RP;}
  Double_t GetV2RP() const {return this->fV2RP;} 
  
  void SetV2SpreadRP(Double_t dV2SpreadRP) {this->fV2SpreadRP = dV2SpreadRP;}
  Double_t GetV2SpreadRP() const {return this->fV2SpreadRP;} 
  
  void SetMinV2RP(Double_t dMinV2RP) {this->fMinV2RP = dMinV2RP;}
  Double_t GetMinV2RP() const {return this->fMinV2RP;} 
  
  void SetMaxV2RP(Double_t dMaxV2RP) {this->fMaxV2RP = dMaxV2RP;}
  Double_t GetMaxV2RP() const {return this->fMaxV2RP;} 
  
  void SetV4RP(Double_t dV4RP) {this->fV4RP = dV4RP;}
  Double_t GetV4RP() const {return this->fV4RP;} 
  
  void SetV4SpreadRP(Double_t dV4SpreadRP) {this->fV4SpreadRP = dV4SpreadRP;}
  Double_t GetV4SpreadRP() const {return this->fV4SpreadRP;} 
  
  // (pt,eta) dependent harmonics:
  void SetV2RPMax(Double_t dV2RPMax) {this->fV2vsPtEtaMax = dV2RPMax;}
  Double_t GetV2RPMax() const {return this->fV2vsPtEtaMax;} 
  
  void SetV2vsPtEtaMax(Double_t dV2vsPtEtaMax) {this->fV2vsPtEtaMax = dV2vsPtEtaMax;}
  Double_t GetV2vsPtEtaMax() const {return this->fV2vsPtEtaMax;} 
  
  void SetPtCutOff(Double_t dPtCutOff) {this->fPtCutOff = dPtCutOff;}
  Double_t GetPtCutOff() const {return this->fPtCutOff;} 
  
  void SetV2vsEtaSpread(Double_t dV2vsEtaSpread) {this->fV2vsEtaSpread = dV2vsEtaSpread;}
  Double_t GetV2vsEtaSpread() const {return this->fV2vsEtaSpread;} 
  
  void SetFirstSectorPhiMin(Double_t dPhiMin1) {this->fPhiMin1 = dPhiMin1;}
  Double_t GetFirstSectorPhiMin() const {return this->fPhiMin1;} 
  
  void SetFirstSectorPhiMax(Double_t dPhiMax1) {this->fPhiMax1 = dPhiMax1;}
  Double_t GetFirstSectorPhiMax() const {return this->fPhiMax1;}
  
  void SetFirstSectorProbability(Double_t dProbability1) {this->fProbability1 = dProbability1;}
  Double_t GetFirstProbability() const {return this->fProbability1;}  
  
  void SetSecondSectorPhiMin(Double_t dPhiMin2) {this->fPhiMin2 = dPhiMin2;}
  Double_t GetSecondSectorPhiMin() const {return this->fPhiMin2;} 
  
  void SetSecondSectorPhiMax(Double_t dPhiMax2) {this->fPhiMax2 = dPhiMax2;}
  Double_t GetSecondSectorPhiMax() const {return this->fPhiMax2;}
  
  void SetSecondSectorProbability(Double_t dProbability2) {this->fProbability2 = dProbability2;}
  Double_t GetSecondProbability() const {return this->fProbability2;}  
  //................................................................................................
  
  void SetNoOfLoops(Int_t noofl) {this->fNoOfLoops = noofl;}
  Int_t GetNoOfLoops() const {return this->fNoOfLoops;} 
  void SetPhiRange(Double_t phr) {this->fPhiRange = phr;}
  Double_t GetPhiRange() const {return this->fPhiRange;}   
  void SetPtRange(Double_t pr) {this->fPtRange = pr;}
  Double_t GetPtRange() const {return this->fPtRange;}   
  void SetEtaRange(Double_t er) {this->fEtaRange = er;}
  Double_t GetEtaRange() const {return this->fEtaRange;} 
  void SetSubeventEtaRange(Double_t minA,Double_t maxA,Double_t minB,Double_t maxB) 
  {this->fEtaMinA = minA; this->fEtaMaxA = maxA;this->fEtaMinB = minB; this->fEtaMaxB = maxB;};

 private:
 
  AliFlowEventSimpleMakerOnTheFly(const AliFlowEventSimpleMakerOnTheFly& anAnalysis);            // copy constructor
  AliFlowEventSimpleMakerOnTheFly& operator=(const AliFlowEventSimpleMakerOnTheFly& anAnalysis); // assignment operator
  
  //................................................................................................
  // global parameters:
  Bool_t    fMultDistrOfRPsIsGauss;  // 1.) if kTRUE  = multiplicitiy of RPs is sampled e-b-e from Gaussian distribution with
                                     //                 mean = fMultiplicityOfRP and spread = fMultiplicitySpreadOfRP
                                     // 2.) if kFALSE = multiplicitiy of RPs is sampled e-b-e uniformly from 
                                     //                 interval [fMinMultOfRP,fMaxMultOfRP]
  Int_t     fMultiplicityOfRP;       // mean multiplicity of RPs (if sampled from Gaussian)
  Double_t  fMultiplicitySpreadOfRP; // multiplicity spread of RPs (if sampled from Gaussian)
  Int_t     fMinMultOfRP;            // minimal multiplicity of RPs (if sampled uniformly)
  Int_t     fMaxMultOfRP;            // maximum multiplicity of RPs (if sampled uniformly)
  Double_t  fTemperatureOfRP;        // "temperature" of RPs in GeV/c (increase this parameter to get more high pt RPs) 
  Bool_t    fPtDependentHarmonics;   // harmonic V2 is a function of pt     
  Bool_t    fEtaDependentHarmonics;  // harmonic V2 is a function of eta     
  // constant harmonics: 
  Double_t  fV1RP;                   // directed flow of RPs
  Double_t  fV1SpreadRP;             // directed flow spread of RPs
  
  Bool_t    fConstantV2IsSampledFromGauss; // 1.) if kTRUE  = elliptic flow of RPs is sampled e-b-e from Gaussian distribution with
                                           //                 mean = fV2RP and spread = fV2SpreadRP
                                           // 2.) if kFALSE = elliptic flow of RPs is sampled e-b-e uniformly from 
                                           //                 interval [fMinV2RP,fMaxV2RP]
  Double_t  fV2RP;                   // mean elliptic flow of RPs (if sampled from Gaussian)
  Double_t  fV2SpreadRP;             // elliptic flow spread of RPs (if sampled from Gaussian)
  Double_t  fMinV2RP;                // minimal elliptic flow of RPs (if sampled uniformly)
  Double_t  fMaxV2RP;                // minimal elliptic flow of RPs (if sampled uniformly)
  
  Double_t  fV4RP;                   // harmonic V4 of RPs
  Double_t  fV4SpreadRP;             // harmonic V4's spread of RPs
  // (pt,eta) dependent harmonics:
  Double_t  fV2vsPtEtaMax;           // max value of (pt,eta) dependent V2
  Double_t  fPtCutOff;               // V2(pt) is linear up to pt = dPtCutOff and for pt > dPtCutOff it is constant, V2(pt) = dV2vsPtEtaMax:
  Double_t  fV2vsEtaSpread;          // V2(eta) is Gaussian centered at midrapidity (eta=0) with spread = dV2vsEtaSpread
  // non-uniform acceptance:
  Double_t  fPhiMin1;                // first non-uniform sector starts at azimuth fPhiMin1
  Double_t  fPhiMax1;                // first non-uniform sector ends at azimuth fPhiMax1
  Double_t  fProbability1;           // particles emitted in fPhiMin1 < phi < fPhiMax1 are taken with probability fProbability1 
  Double_t  fPhiMin2;                // second non-uniform sector starts at azimuth fPhiMin2
  Double_t  fPhiMax2;                // second non-uniform sector starts at azimuth fPhiMax2
  Double_t  fProbability2;           // particles emitted in fPhiMin2 < phi < fPhiMax2 are taken with probability fProbability2
  //................................................................................................
  
  //................................................................................................
  // equations for distributions: 
  TF1*      fPtSpectra;  // transverse momentum distribution
  TF1*      fPhiDistribution; // azimuthal distribution
  //................................................................................................
  
  TRandom3* fMyTRandom3; // our TRandom3 generator
  Int_t     fCount;      // count number of events 
  Int_t     fNoOfLoops;  // number of times to use the same particle (nonflow)
  Double_t  fPhiRange;   // splitted track phi range (+/- from original track's phi) for uniform sampling
  Double_t  fPtRange;    // splitted track pt range (+/- from original track's pt) for uniform sampling
  Double_t  fEtaRange;   // splitted track eta range (+/- from original track's eta) for uniform sampling
  Double_t  fEtaMinA;    // minimum eta of subevent A eta range
  Double_t  fEtaMaxA;    // maximum eta of subevent A eta range
  Double_t  fEtaMinB;    // minimum eta of subevent B eta range
  Double_t  fEtaMaxB;    // maximum eta of subevent B eta range  

  ClassDef(AliFlowEventSimpleMakerOnTheFly,0) // macro for rootcint
};
 
#endif



