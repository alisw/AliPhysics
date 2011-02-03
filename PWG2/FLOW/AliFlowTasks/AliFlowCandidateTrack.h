/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id: $ */

#ifndef AliFlowCandidateTrack_H
#define AliFlowCandidateTrack_H

#include "AliFlowTrack.h"

////////////////////////////////////////////////////
// AliFlowCandidateTrack:
// Class for reconstructed particles to be used in flow analysis
// Author: Carlos Perez (cperez@cern.ch)
////////////////////////////////////////////////////

class AliFlowCandidateTrack : public AliFlowTrack {
  protected:
    Double_t fMass;           // mass
    Int_t fNDaughters;        // number of daughters (5 max)
    Int_t fDaughter[5];       // fID of daughter, points back to ESD track
    AliFlowTrack *fTrack[5];  // pointer to daughter in FlowEvent
    
  public:
    AliFlowCandidateTrack();
    AliFlowCandidateTrack(const AliFlowCandidateTrack& );
    AliFlowCandidateTrack& operator=(const AliFlowCandidateTrack& );
    ~AliFlowCandidateTrack();

    Double_t Mass(void)            { return fMass; }
    void SetMass(Double_t value)   { fMass=value; }

    Int_t GetNDaughters(void)           { return fNDaughters; }
    void  AddDaughter(Int_t value)      { if(fNDaughters<3) fDaughter[fNDaughters++]=value; }
    Int_t GetIDDaughter(Int_t value) { return fDaughter[value]; }

    void SetDaughter(Int_t value, AliFlowTrack *track) { fTrack[value]=track; }
    AliFlowTrack *GetDaughter(Int_t value) { return fTrack[value]; }


    ClassDef(AliFlowCandidateTrack, 1);
};

#endif
