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

#ifndef ALIRSNCUTESDPRIMARY_H
#define ALIRSNCUTESDPRIMARY_H

#include "AliESDtrackCuts.h"
#include "AliRsnCut.h"

class AliRsnCutESDPrimary : public AliRsnCut
{
  public:

    AliRsnCutESDPrimary();
    AliRsnCutESDPrimary(const char *name);
    virtual ~AliRsnCutESDPrimary() {;};

    AliESDtrackCuts* GetCuts() {return &fCuts;}
    virtual Bool_t   IsSelected(TObject *object);

  protected:

    AliESDtrackCuts fCuts;  // set of ESD track cuts

    ClassDef(AliRsnCutESDPrimary, 1)
};

#endif
