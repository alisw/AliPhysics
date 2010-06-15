/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

/*****************************************************************
  AliFlowEvent: Event container for flow analysis                  
                                     
  origin:   Mikolaj Krzewicki  (mikolaj.krzewicki@cern.ch)
*****************************************************************/

#ifndef ALIFLOWEVENT_H
#define ALIFLOWEVENT_H

class TTree;
class AliFlowTrackSimpleCuts;
class AliFlowTrack;
class AliCFManager;
class AliMCEvent;
class AliESDEvent;
class AliMCEvent;
class AliAODEvent;
class AliMultiplicity;

#include "AliFlowEventSimple.h"

class AliFlowEvent: public AliFlowEventSimple {
public:

  enum KineSource { kNoKine, kESDkine, kMCkine };

  AliFlowEvent();
  AliFlowEvent(const AliFlowEvent& event);
  AliFlowEvent& operator=(const AliFlowEvent& event);
  virtual  ~AliFlowEvent() {}

  AliFlowEvent( const AliMCEvent* anInput,
                const AliCFManager* rpCFManager=NULL,
                const AliCFManager* poiCFManager=NULL ); //use CF(2x)
  AliFlowEvent( const AliESDEvent* anInput,  
                const AliCFManager* rpCFManager=NULL, 
                const AliCFManager* poiCFManager=NULL ); //use CF(2x)
  AliFlowEvent( const AliAODEvent* anInput, 
                const AliCFManager* rpCFManager=NULL, 
                const AliCFManager* poiCFManager=NULL );  //use CF(2x)
  AliFlowEvent( const AliESDEvent* anInput, 
                const AliMCEvent* anInputMc, 
                KineSource anOption=kNoKine,
                const AliCFManager* rpCFManager=NULL, 
                const AliCFManager* poiCFManager=NULL );  //use CF(2x)
  AliFlowEvent( const AliESDEvent* anInput,
                const AliMultiplicity* anInputTracklets,
                const AliCFManager* poiCFManager );

  void SetMCReactionPlaneAngle(const AliMCEvent* mcEvent);

  AliFlowTrack* GetTrack( Int_t i );

  ClassDef(AliFlowEvent,1)
};

#endif


