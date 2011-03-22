/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

/*****************************************************************
  AliFlowEvent: Event container for flow analysis                  
                                     
  origin:   Mikolaj Krzewicki  (mikolaj.krzewicki@cern.ch)
*****************************************************************/

#ifndef ALIFLOWEVENT_H
#define ALIFLOWEVENT_H

class AliFlowTrackCuts;
class AliFlowTrack;
class AliCFManager;
class AliVEvent;
class AliMCEvent;
class AliESDEvent;
class AliAODEvent;
class AliMultiplicity;
class AliESDPmdTrack;
class TH2F;

#include "AliFlowEventSimple.h"

class AliFlowEvent: public AliFlowEventSimple {
public:

  enum KineSource { kNoKine, kESDkine, kMCkine };

  AliFlowEvent();
  AliFlowEvent(Int_t n);
  AliFlowEvent(const AliFlowEvent& event);
  AliFlowEvent& operator=(const AliFlowEvent& event);
  virtual  ~AliFlowEvent() {}

  //deprecated
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
                const AliCFManager* poiCFManager,
                Bool_t hybrid);
  AliFlowEvent( const AliESDEvent* anInput, 
                const AliMCEvent* anInputMc, 
                KineSource anOption=kNoKine,
                const AliCFManager* rpCFManager=NULL, 
                const AliCFManager* poiCFManager=NULL );  //use CF(2x)
  AliFlowEvent( const AliESDEvent* anInput,
                const AliMultiplicity* anInputTracklets,
                const AliCFManager* poiCFManager );
  AliFlowEvent( const AliESDEvent* anInput,
                const TH2F* anInputFMDhist,
                const AliCFManager* poiCFManager );
  //pmd
  AliFlowEvent( const AliESDEvent* anInput,
                const AliESDPmdTrack *pmdtracks,
                const AliCFManager* poiCFManager );
  //pmd
  //end of deprecated

  AliFlowEvent( AliFlowTrackCuts* rpCuts,
                AliFlowTrackCuts* poiCuts );
  
  void Fill( AliFlowTrackCuts* rpCuts,
             AliFlowTrackCuts* poiCuts );

  void SetMCReactionPlaneAngle(const AliMCEvent* mcEvent);
  using AliFlowEventSimple::SetMCReactionPlaneAngle;

  AliFlowTrack* GetTrack( Int_t i );

protected:
  AliFlowTrack* ReuseTrack( Int_t i);

  ClassDef(AliFlowEvent,1)
};

#endif


