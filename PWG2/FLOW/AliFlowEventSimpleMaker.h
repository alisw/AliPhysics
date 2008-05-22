/* copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#ifndef AliFlowEventSimpleMaker_H
#define AliFlowEventSimpleMaker_H

#include "AliFlowTrackSimple.h" //needed as include
#include "AliFlowEventSimple.h"  //needed as include

//class AliFlowTrackSimple; //does not compile
//class AliFlowEventSimple; //does not compile

class TTree;
class TParticle;
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
  
  AliFlowEventSimple* FillTracks(TTree* fInput);
  AliFlowEventSimple* FillTracks(AliMCEvent* fInput);
  AliFlowEventSimple* FillTracks(AliESDEvent* fInput);
  AliFlowEventSimple* FillTracks(AliESDEvent* fInput, AliMCEvent* fInputMc, Int_t fOption);
  // fOption = 0 : kine from ESD
  // fOption = 1 : kine from MC
  AliFlowEventSimple* FillTracks(AliAODEvent* fInput);
    
 private:
  
  AliFlowEventSimple*   fEvent;      //!
  AliFlowTrackSimple*   fTrack;      //!
  TParticle*            fParticle;   //!
  
  ClassDef(AliFlowEventSimpleMaker,0)    // macro for rootcint
};
 
     
#endif

