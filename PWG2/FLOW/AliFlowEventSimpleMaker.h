/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#ifndef AliFlowEventSimpleMaker_H
#define AliFlowEventSimpleMaker_H

#include "AliFlowEventSimple.h"  //needed as include
//class AliFlowEventSimple; //does not compile

#include "../../CORRFW/AliCFManager.h"
//class AliCFManager;
#include "AliFlowTrackSimpleCuts.h"

class TTree;
class AliMCEvent;
class AliESDEvent;
class AliAODEvent;

// AliFlowEventSimpleMaker:
// Class to fill the AliFlowEventSimple with AliFlowTrackSimple objects
// author: N. van der Kolk (kolk@nikhef.nl)
          
class AliFlowEventSimpleMaker {

 public:

  AliFlowEventSimpleMaker();             //constructor
  virtual ~AliFlowEventSimpleMaker();    //destructor
  
  //TTree
  AliFlowEventSimple* FillTracks(TTree* anInput, AliFlowTrackSimpleCuts* intCuts, AliFlowTrackSimpleCuts* diffCuts);   //use own cut class
  //AliMCEvent
  AliFlowEventSimple* FillTracks(AliMCEvent* anInput);   //use own cuts
  AliFlowEventSimple* FillTracks(AliMCEvent* anInput, AliCFManager* intCFManager, AliCFManager* diffCFManager ); //use CF(2x)
  //AliESDEvent
  AliFlowEventSimple* FillTracks(AliESDEvent* anInput);   //use own cuts
  AliFlowEventSimple* FillTracks(AliESDEvent* anInput,  AliCFManager* intCFManager, AliCFManager* diffCFManager); //use CF(2x)
  //AliESDEvent & AliMCEvent
  AliFlowEventSimple* FillTracks(AliESDEvent* anInput, AliMCEvent* anInputMc, Int_t anOption);  //use own cuts
  AliFlowEventSimple* FillTracks(AliESDEvent* anInput, AliMCEvent* anInputMc, AliCFManager* intCFManager, AliCFManager* diffCFManager, Int_t anOption);  //use CF(2x)
  // anOption = 0 : kine from ESD
  // anOption = 1 : kine from MC
  //AliAODEvent
  AliFlowEventSimple* FillTracks(AliAODEvent* anInput); //use own cuts
  AliFlowEventSimple* FillTracks(AliAODEvent* anInput, AliCFManager* intCFManager, AliCFManager* diffCFManager);  //use CF(2x)
    
 private:
  AliFlowEventSimpleMaker(const AliFlowEventSimpleMaker& anAnalysis);            //copy constructor
  AliFlowEventSimpleMaker& operator=(const AliFlowEventSimpleMaker& anAnalysis); //assignment operator
       
  ClassDef(AliFlowEventSimpleMaker,0)    // macro for rootcint
};
 
     
#endif

