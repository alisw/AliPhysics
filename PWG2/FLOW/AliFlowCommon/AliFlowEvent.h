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
class AliCFManager;
class AliMCEvent;
class AliESDEvent;
class AliMCEvent;
class AliAODEvent;

#include "AliFlowEventSimple.h"

class AliFlowEvent: public AliFlowEventSimple {
public:
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
                Int_t anOption=-1,
                const AliCFManager* rpCFManager=NULL, 
                const AliCFManager* poiCFManager=NULL );  //use CF(2x)

  void SetMCReactionPlaneAngle(const AliMCEvent* mcEvent);

  ClassDef(AliFlowEvent,1)
};

#endif


