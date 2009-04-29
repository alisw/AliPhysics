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

  void SetV1RP(Double_t dV1RP) {this->fV1RP = dV1RP;}
  Double_t GetV1RP() const {return this->fV1RP;} 
  
  void SetV1SpreadRP(Double_t dV1SpreadRP) {this->fV1SpreadRP = dV1SpreadRP;}
  Double_t GetV1SpreadRP() const {return this->fV1SpreadRP;} 
  
  void SetV2RP(Double_t dV2RP) {this->fV2RP = dV2RP;}
  Double_t GetV2RP() const {return this->fV2RP;} 
  
  void SetV2SpreadRP(Double_t dV2SpreadRP) {this->fV2SpreadRP = dV2SpreadRP;}
  Double_t GetV2SpreadRP() const {return this->fV2SpreadRP;} 
  //................................................................................................
  
 private:
 
  AliFlowEventSimpleMakerOnTheFly(const AliFlowEventSimpleMakerOnTheFly& anAnalysis);            // copy constructor
  AliFlowEventSimpleMakerOnTheFly& operator=(const AliFlowEventSimpleMakerOnTheFly& anAnalysis); // assignment operator
  
  //................................................................................................
  // global parameters:
  Int_t fMultiplicityOfRP;          // multiplicity of RPs
  Double_t fMultiplicitySpreadOfRP; // multiplicity spread of RPs  
  Double_t fV1RP;                   // directed flow of RPs
  Double_t fV1SpreadRP;             // directed flow spread of RPs
  Double_t fV2RP;                   // elliptic flow of RPs
  Double_t fV2SpreadRP;             // elliptic flow spread of RPs
  //................................................................................................
  
  //................................................................................................
  // equations for distributions: 
  TF1 *fPtSpectra;  // transverse momentum distribution
  TF1 *fPhiDistribution; // azimuthal distribution
  //................................................................................................
  
  TRandom3* fMyTRandom3; // our random generator
  Int_t fCount;

  ClassDef(AliFlowEventSimpleMakerOnTheFly,0) // macro for rootcint
};
 
#endif



