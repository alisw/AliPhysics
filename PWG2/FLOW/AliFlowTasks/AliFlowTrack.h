/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef ALIFLOWTRACK_H
#define ALIFLOWTRACK_H

#include "AliFlowTrackSimple.h"

// AliFlowTrack:
// A track class to the the AliFlowEvent for flow analysis
// origin: Mikolaj Krzewicki (mikolaj.krzewicki@cern.ch)

class AliFlowTrack: public AliFlowTrackSimple {

public:
  enum trackSource {kFromESD=0,
                    kFromMC=1,
                    kFromAOD=2,
                    kFromTracklet=3};
  AliFlowTrack();
  AliFlowTrack& operator=(const AliFlowTrack& aTrack);
  AliFlowTrack(const AliFlowTrack& aTrack);
  virtual  ~AliFlowTrack();
 
  void SetFMDMultiplicity( const Float_t m ) {fFMDmultiplicity=m;} 
  Float_t GetFMDMultiplicity() const {return fFMDmultiplicity;}
  void SetSource( trackSource s )
                  { fTrackSourceBits.SetBitNumber(UInt_t(s),kTRUE); }
  Bool_t IsSource( trackSource s ) const
                 { return fTrackSourceBits.TestBitNumber(s); }



private:
  TBits fTrackSourceBits; //where do i come from?
  Float_t fFMDmultiplicity; //FMD multiplicity
  

  ClassDef(AliFlowTrack,1);
};

#endif

