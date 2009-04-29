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

#include "AliFlowEventSimple.h"  //needed as include
    
class AliFlowEventSimpleMakerOnTheFly {

 public:

  AliFlowEventSimpleMakerOnTheFly();          // constructor
  virtual ~AliFlowEventSimpleMakerOnTheFly(); // destructor

  AliFlowEventSimple* CreateEventOnTheFly();  // create an event on the fly
 
    
  //                        *****************************
  //                        **** SETTERS AND GETTERS ****
  //                        *****************************
  //................................................................................................
  // setters and getters for global parameters:
  void SetMultiplicityOfRP(Int_t multRP) {this->fMultiplicityOfRP = multRP;}
  Int_t GetMultiplicityOfRP() const {return this->fMultiplicityOfRP;} 
  //................................................................................................
  
  //................................................................................................
  // setters and getters for event-by-event parameters:
  void SetMultiplicitySpreadOfRP(Double_t multSpreadRP) {this->fMultiplicitySpreadOfRP = multSpreadRP;}
  Double_t GetMultiplicitySpreadOfRP() const {return this->fMultiplicitySpreadOfRP;} 
  //................................................................................................
   
 private:
 
  AliFlowEventSimpleMakerOnTheFly(const AliFlowEventSimpleMakerOnTheFly& anAnalysis);            // copy constructor
  AliFlowEventSimpleMakerOnTheFly& operator=(const AliFlowEventSimpleMakerOnTheFly& anAnalysis); // assignment operator
  
  //................................................................................................
  // global parameters:
  Int_t fMultiplicityOfRP; // multiplicity of RPs
  //................................................................................................
  
  //................................................................................................
  // event-by-event parameters:
  Double_t fMultiplicitySpreadOfRP; // multiplicity spread of RPs
  //................................................................................................
  
  //................................................................................................
  // equations for distributions: 
  TF1 *fPtFormula;  // transverse momentum distribution
  TF1 *fPhiFormula; // azimuthal distribution
  //................................................................................................
  
  ClassDef(AliFlowEventSimpleMakerOnTheFly,0) // macro for rootcint
};
 
#endif



