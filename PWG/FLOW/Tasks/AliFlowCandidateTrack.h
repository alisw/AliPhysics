/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id: $ */

#ifndef ALIFLOWCANDIDATETRACK_H
#define ALIFLOWCANDIDATETRACK_H

#include "AliFlowTrack.h"

////////////////////////////////////////////////////
// AliFlowCandidateTrack:
// Class for reconstructed particles to be used in flow analysis
// Author: Carlos Perez (cperez@cern.ch)
////////////////////////////////////////////////////

class AliFlowCandidateTrack : public AliFlowTrack {
  public:
    AliFlowCandidateTrack();
    AliFlowCandidateTrack(const AliFlowCandidateTrack& );
    AliFlowCandidateTrack& operator=(const AliFlowCandidateTrack& );
    ~AliFlowCandidateTrack();

    Int_t GetNDaughters(void)        const { return fNDaughters; }
    void  AddDaughter(Int_t value)  { if(fNDaughters<3) fDaughter[fNDaughters++]=value; }
    Int_t GetIDDaughter(Int_t value) const { return fDaughter[value]; }

    void SetDaughter(Int_t value, AliFlowTrack *track) { fTrack[value]=track; }
    AliFlowTrack *GetDaughter(Int_t value) const { return fTrack[value]; }

  protected:
    Int_t fNDaughters;        // number of daughters (5 max)
    Int_t fDaughter[5];       // fID of daughter, points back to ESD track
    AliFlowTrack *fTrack[5];  // pointer to daughter in FlowEvent
    

    ClassDef(AliFlowCandidateTrack, 2);
};

#endif
