//
// Class AliRsnSimpleFcnResolution
//
// This is the most fundamental AliRsnSimpleFunction,
// which computes the invariant mass spectrum of a resonance,
// by correlating pairs of tracks from an event (or mixing, for BG)
//

#ifndef AliRsnSimpleFcnResolution_H
#define AliRsnSimpleFcnResolution_H

#include "AliRsnSimpleFunction.h"

class AliRsnSimpleFcnResolution : public AliRsnSimpleFunction
{
  public:

    AliRsnSimpleFcnResolution();
    AliRsnSimpleFcnResolution
    (const char* name, AliRsnPairDef *pd,
     AliRsnHistoDef *hd, AliRsnCutMgr *cuts=0, Option_t *option="");

    // virtual working routines
    virtual Bool_t ProcessOne(AliRsnEvent *event);

  private:

    // dummy required methods
    AliRsnSimpleFcnResolution(const AliRsnSimpleFcnResolution &copy) :
        AliRsnSimpleFunction(copy) { /*nothing */ }
    AliRsnSimpleFcnResolution& operator=(const AliRsnSimpleFcnResolution& /*copy*/)
    { /* nothing */ return (*this); }

    // utility methods
    Bool_t Add(AliRsnDaughter *t1, AliRsnDaughter *t2);

    ClassDef(AliRsnSimpleFcnResolution, 1); // dictionary
};

#endif
