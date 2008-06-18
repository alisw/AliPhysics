/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#ifndef AliFlowEventSimpleMaker_H
#define AliFlowEventSimpleMaker_H

#include "AliFlowEventSimple.h"  //needed as include
//class AliFlowEventSimple; //does not compile

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
  virtual  ~AliFlowEventSimpleMaker();   //destructor
  
  AliFlowEventSimple* FillTracks(TTree* anInput);
  AliFlowEventSimple* FillTracks(AliMCEvent* anInput);
  AliFlowEventSimple* FillTracks(AliESDEvent* anInput);
  AliFlowEventSimple* FillTracks(AliESDEvent* anInput, AliMCEvent* anInputMc, Int_t anOption);
  // anOption = 0 : kine from ESD
  // anOption = 1 : kine from MC
  AliFlowEventSimple* FillTracks(AliAODEvent* anInput);
    
 private:
  AliFlowEventSimpleMaker(const AliFlowEventSimpleMaker& anAnalysis);            //copy constructor
  AliFlowEventSimpleMaker& operator=(const AliFlowEventSimpleMaker& anAnalysis); //assignment operator
          
  ClassDef(AliFlowEventSimpleMaker,0)    // macro for rootcint
};
 
     
#endif

