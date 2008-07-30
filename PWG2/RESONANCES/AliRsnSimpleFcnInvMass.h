//
// Class AliRsnSimpleFcnInvMass
//
// This is the most fundamental AliRsnSimpleFunction,
// which computes the invariant mass spectrum of a resonance,
// by correlating pairs of tracks from an event (or mixing, for BG)
//

#ifndef AliRsnSimpleFcnInvMass_H
#define AliRsnSimpleFcnInvMass_H

#include "AliRsnDaughter.h"
#include "AliRsnSimpleFunction.h"

class AliRsnSimpleFcnInvMass : public AliRsnSimpleFunction
{
public:

    AliRsnSimpleFcnInvMass();
    AliRsnSimpleFcnInvMass
      (const char* name, AliRsnDaughter::EPIDMethod method, AliRsnPairDef *pd, 
       AliRsnHistoDef *hd, AliRsnCutMgr *cuts=0, Option_t *option="");

    void               SetUseMCValues(Bool_t value) {fUseMCValues = value;}
    Bool_t             UseMCValues() {return fUseMCValues;}

    // virtual working routines
    virtual Bool_t ProcessOne(AliRsnEvent *event);
    virtual Bool_t ProcessTwo(AliRsnEvent *event1, AliRsnEvent *event2);

private:

    // dummy required methods
    AliRsnSimpleFcnInvMass(const AliRsnSimpleFcnInvMass &copy) :
      AliRsnSimpleFunction(copy),fUseMCValues(kFALSE),fPIDMethod(AliRsnDaughter::kNoPID)
      { /*nothing */ }
    AliRsnSimpleFcnInvMass& operator=(const AliRsnSimpleFcnInvMass& /*copy*/)
      { /* nothing */ return (*this); }

    // utility methods
    Bool_t Add(AliRsnDaughter *t1, AliRsnDaughter *t2);

    Bool_t                     fUseMCValues;  // flag to switch between MC/rec values
    AliRsnDaughter::EPIDMethod fPIDMethod;    // flag to PID method used

    ClassDef(AliRsnSimpleFcnInvMass, 1); // dictionary
};

#endif
