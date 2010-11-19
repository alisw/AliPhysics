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

#ifndef ALIRSNCUTESDCUTMULTIPLICITY_H
#define ALIRSNCUTESDCUTMULTIPLICITY_H

#include "AliESDtrackCuts.h"
#include "AliRsnCut.h"

class AliRsnCutESDCutMultiplicity : public AliRsnCut
{
  public:

    AliRsnCutESDCutMultiplicity();
    AliRsnCutESDCutMultiplicity(const char *name, Int_t min, Int_t max);
    virtual ~AliRsnCutESDCutMultiplicity() {;};

    const AliESDtrackCuts* GetCuts() const {return &fCuts;}
    virtual Bool_t IsSelected(TObject *object);

  protected:

    AliESDtrackCuts fCuts;  // set of ESD track cuts

    ClassDef(AliRsnCutESDCutMultiplicity, 1)
};

#endif
