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
//class TUnuran;

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
  void SetMultiplicityOfRP(Int_t multRP) {this->fMultiplicityOfRP = multRP;}
  Int_t GetMultiplicityOfRP() const {return this->fMultiplicityOfRP;} 
  
  void SetMultiplicitySpreadOfRP(Double_t multSpreadRP) {this->fMultiplicitySpreadOfRP = multSpreadRP;}
  Double_t GetMultiplicitySpreadOfRP() const {return this->fMultiplicitySpreadOfRP;} 
  
  void SetUseConstantHarmonics(Bool_t const uch) {this->fUseConstantHarmonics = uch;};
  Bool_t GetUseConstantHarmonics() const {return this->fUseConstantHarmonics;};
  
  // constant harmonics:  
  void SetV1RP(Double_t dV1RP) {this->fV1RP = dV1RP;}
  Double_t GetV1RP() const {return this->fV1RP;} 
  
  void SetV1SpreadRP(Double_t dV1SpreadRP) {this->fV1SpreadRP = dV1SpreadRP;}
  Double_t GetV1SpreadRP() const {return this->fV1SpreadRP;} 
  
  void SetV2RP(Double_t dV2RP) {this->fV2RP = dV2RP;}
  Double_t GetV2RP() const {return this->fV2RP;} 
  
  void SetV2SpreadRP(Double_t dV2SpreadRP) {this->fV2SpreadRP = dV2SpreadRP;}
  Double_t GetV2SpreadRP() const {return this->fV2SpreadRP;} 
  
  void SetV4RP(Double_t dV4RP) {this->fV4RP = dV4RP;}
  Double_t GetV4RP() const {return this->fV4RP;} 
  
  void SetV4SpreadRP(Double_t dV4SpreadRP) {this->fV4SpreadRP = dV4SpreadRP;}
  Double_t GetV4SpreadRP() const {return this->fV4SpreadRP;} 
  
  // (pt,eta) dependent harmonics:
  void SetV2RPMax(Double_t dV2RPMax) {this->fV2RPMax = dV2RPMax;}
  Double_t GetV2RPMax() const {return this->fV2RPMax;} 
  
  void SetPtCutOff(Double_t dPtCutOff) {this->fPtCutOff = dPtCutOff;}
  Double_t GetPtCutOff() const {return this->fPtCutOff;} 
  //................................................................................................
  
  void SetNoOfLoops(Int_t noofl) {this->fNoOfLoops = noofl;}
  Int_t GetNoOfLoops() const {return this->fNoOfLoops;} 

 private:
 
  AliFlowEventSimpleMakerOnTheFly(const AliFlowEventSimpleMakerOnTheFly& anAnalysis);            // copy constructor
  AliFlowEventSimpleMakerOnTheFly& operator=(const AliFlowEventSimpleMakerOnTheFly& anAnalysis); // assignment operator
  
  //................................................................................................
  // global parameters:
  Int_t     fMultiplicityOfRP;       // multiplicity of RPs
  Double_t  fMultiplicitySpreadOfRP; // multiplicity spread of RPs 
  Bool_t    fUseConstantHarmonics;      // harmonics V1, V2, V4... are constant (kTRUE) or functions of pt and eta (kFALSE)     
  // constant harmonics: 
  Double_t  fV1RP;                   // directed flow of RPs
  Double_t  fV1SpreadRP;             // directed flow spread of RPs
  Double_t  fV2RP;                   // elliptic flow of RPs
  Double_t  fV2SpreadRP;             // elliptic flow spread of RPs
  Double_t  fV4RP;                   // harmonic V4 of RPs
  Double_t  fV4SpreadRP;             // harmonic V4's spread of RPs
  // (pt,eta) dependent harmonics:
  Double_t  fV2RPMax;                // elliptic flow of RPs
  Double_t  fPtCutOff;               // elliptic flow spread of RPs
  //................................................................................................
  
  //................................................................................................
  // equations for distributions: 
  TF1*      fPtSpectra;  // transverse momentum distribution
  TF1*      fPhiDistribution; // azimuthal distribution
  //................................................................................................
  
  TRandom3* fMyTRandom3; // our TRandom3 generator
  //TUnuran*  fMyUnuran;   // our TUnuran generator
  Int_t     fCount;      // count number of events 
  Int_t     fNoOfLoops;  // number of times to use the same particle (nonflow)

  ClassDef(AliFlowEventSimpleMakerOnTheFly,0) // macro for rootcint
};
 
#endif



