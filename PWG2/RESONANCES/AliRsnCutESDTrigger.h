//
// Class AliRsnCutRange
//
// General implementation of cuts which check a value inside a range.
// This range can be defined by two integers or two doubles.
// A user-friendly enumeration allows to define what is checked.
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#ifndef ALIRSNCUTESDTrigger_H
#define ALIRSNCUTESDTrigger_H

#include <TString.h>
#include "AliRsnCut.h"

class AliRsnCutESDTrigger : public AliRsnCut
{
  public:

    AliRsnCutESDTrigger();
    AliRsnCutESDTrigger(const char *name, const char *mask);
    virtual ~AliRsnCutESDTrigger() {;};
      
    void    SetTriggerMask(const char *mask) {fTrigger = mask;}

    virtual Bool_t IsSelected(AliRsnCut::ETarget tgt, AliRsnDaughter* daughter);
    virtual Bool_t IsSelected(ETarget tgt, AliRsnPairParticle *pair);
    virtual Bool_t IsSelected(ETarget tgt, AliRsnEvent * const event);
    virtual Bool_t IsSelected(ETarget tgt, AliRsnEvent *ev1, AliRsnEvent *ev2);

  protected:

    TString fTrigger;  // trigger mask

    ClassDef(AliRsnCutESDTrigger, 1)
};

#endif
