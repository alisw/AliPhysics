/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

/*****************************************************************
  AliFlowEventStar: Event container for flow analysis of star data
                                     
  origin:   Mikolaj Krzewicki  (mikolaj.krzewicki@cern.ch)
*****************************************************************/

#ifndef ALIFLOWEVENTSTAR_H
#define ALIFLOWEVENTSTAR_H

#include "AliFlowEventSimple.h"

class AliStarEvent;
class AliStarTrackCuts;

class AliFlowEventStar: public AliFlowEventSimple {

 public:
  AliFlowEventStar();
  AliFlowEventStar( const AliFlowEventStar& event );
  AliFlowEventStar( const AliStarEvent* event,
                    const AliStarTrackCuts* rpCuts=NULL,
                    const AliStarTrackCuts* poiCuts=NULL);
  AliFlowEventStar& operator=( const AliFlowEventStar& event );
  virtual  ~AliFlowEventStar() {}

  ClassDef(AliFlowEventStar,1)
};

#endif


