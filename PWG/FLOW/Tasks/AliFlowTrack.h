/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef ALIFLOWTRACK_H
#define ALIFLOWTRACK_H

#include "AliFlowTrackSimple.h"
class AliVParticle;

// AliFlowTrack:
// A track class to the the AliFlowEvent for flow analysis
// origin: Mikolaj Krzewicki (mikolaj.krzewicki@cern.ch)

class AliFlowTrack: public AliFlowTrackSimple {

public:
  enum trackSource { kFromESD=0,
                     kFromMC=1,
                     kFromAOD=2,
                     kFromTracklet=3,
                     kFromFMD=4,
                     kFromPMD=5,
                     kFromV0=6 };
  AliFlowTrack();
  AliFlowTrack(const AliVParticle* p);
  AliFlowTrack& operator=(const AliFlowTrack& aTrack);
  //virtual AliFlowTrackSimple& operator=(const AliFlowTrackSimple& aTrack);
  AliFlowTrack(const AliFlowTrack& aTrack);
  virtual  ~AliFlowTrack();
  virtual AliFlowTrack* Clone(const char* option="") const;

  void Set(const AliVParticle* p);
 
  void SetSource( trackSource s )
                  { fTrackSourceBits.SetBitNumber(UInt_t(s),kTRUE); }
  Bool_t IsSource( trackSource s ) const
                 { return fTrackSourceBits.TestBitNumber(s); }

private:
  TBits fTrackSourceBits; //where do i come from?
  
  ClassDef(AliFlowTrack,1);
};

#endif

