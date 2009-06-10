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

  AliESDtrackCuts GetCuts() {return fCuts;}

  virtual Bool_t IsSelected(ETarget tgt, AliRsnDaughter *daughter);
  virtual Bool_t IsSelected(ETarget tgt, AliRsnPairParticle *pair);
  virtual Bool_t IsSelected(ETarget tgt, AliRsnEvent *event);
  virtual Bool_t IsSelected(ETarget tgt, AliRsnEvent *ev1, AliRsnEvent *ev2);

protected:

  AliESDtrackCuts fCuts;  // set of ESD track cuts

  ClassDef(AliRsnCutESDPrimary, 1)
};

#endif
