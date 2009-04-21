/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#ifndef AliFlowEventSimpleMaker_H
#define AliFlowEventSimpleMaker_H

#include "AliFlowEventSimple.h"  //needed as include
#include "AliFlowTrackSimpleCuts.h"

class TTree;
class AliCFManager;
class AliMCEvent;
class AliESDEvent;
class AliAODEvent;

// AliFlowEventSimpleMaker:
// Class to fill the AliFlowEventSimple with AliFlowTrackSimple objects
// author: N. van der Kolk (kolk@nikhef.nl),  Ante Bilandzic (anteb@nikhef.nl),  Raimond Snellings (Raimond.Snellings@nikhef.nl) 
          
class AliFlowEventSimpleMaker {

 public:

  AliFlowEventSimpleMaker();             //constructor
  virtual ~AliFlowEventSimpleMaker();    //destructor

  void SetMCReactionPlaneAngle(Double_t fPhiRP)  { this->fMCReactionPlaneAngle = fPhiRP; } 
  //TTree
  AliFlowEventSimple* FillTracks(TTree* anInput, AliFlowTrackSimpleCuts* rpCuts, AliFlowTrackSimpleCuts* poiCuts);   //use own cut class
  //AliMCEvent
  AliFlowEventSimple* FillTracks(AliMCEvent* anInput);   //use own cuts
  AliFlowEventSimple* FillTracks(AliMCEvent* anInput, AliCFManager* rpCFManager, AliCFManager* poiCFManager ); //use CF(2x)
  //AliESDEvent
  AliFlowEventSimple* FillTracks(AliESDEvent* anInput);   //use own cuts
  AliFlowEventSimple* FillTracks(AliESDEvent* anInput,  AliCFManager* rpCFManager, AliCFManager* poiCFManager); //use CF(2x)
  //AliESDEvent & AliMCEvent
  AliFlowEventSimple* FillTracks(AliESDEvent* anInput, AliMCEvent* anInputMc, Int_t anOption);  //use own cuts
  AliFlowEventSimple* FillTracks(AliESDEvent* anInput, AliMCEvent* anInputMc, AliCFManager* rpCFManager, AliCFManager* poiCFManager, Int_t anOption);  //use CF(2x)
  // anOption = 0 : kine from ESD
  // anOption = 1 : kine from MC
  //AliAODEvent
  AliFlowEventSimple* FillTracks(AliAODEvent* anInput); //use own cuts
  AliFlowEventSimple* FillTracks(AliAODEvent* anInput, AliCFManager* rpCFManager, AliCFManager* poiCFManager);  //use CF(2x)
  
  void SetNoOfLoops(Int_t noofl) {this->fNoOfLoops = noofl;}
  Int_t GetNoOfLoops() const {return this->fNoOfLoops;} 
  
  void SetEllipticFlowValue(Double_t elfv) {this->fEllipticFlowValue = elfv;}
  Double_t GetEllipticFlowValue() const {return this->fEllipticFlowValue;} 
  
  void SetMultiplicityOfEvent(Int_t multevnt) {this->fMultiplicityOfEvent = multevnt;}
  Int_t GetMultiplicityOfEvent() const {return this->fMultiplicityOfEvent;} 
  
  
 private:
  AliFlowEventSimpleMaker(const AliFlowEventSimpleMaker& anAnalysis);            //copy constructor
  AliFlowEventSimpleMaker& operator=(const AliFlowEventSimpleMaker& anAnalysis); //assignment operator
  Double_t  fMCReactionPlaneAngle;   // the angle of the reaction plane from the MC truth
  Int_t     fCount;   // counter for the number of events processed

  Int_t     fNoOfLoops; // number of times to use the same particle (nonflow) 
  Double_t  fEllipticFlowValue; // Add Flow. Must be in range [0,1].
  Int_t     fMultiplicityOfEvent; // Set maximal multiplicity.
      
  ClassDef(AliFlowEventSimpleMaker,0)    // macro for rootcint
};
 
     
#endif

