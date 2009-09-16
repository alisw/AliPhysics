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

#ifndef ALIRSNCUTPRIMARYVERTEX_H
#define ALIRSNCUTPRIMARYVERTEX_H

#include "AliRsnCut.h"
#include "AliPID.h"

#include "Riostream.h"
#include "TMath.h"

#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"

#include "AliRsnEvent.h"

class AliRsnEvent;
class AliRsnDaughter;
class AliRsnPairParticle;
class AliRsnCutPrimaryVertex : public AliRsnCut
{
  public:

    AliRsnCutPrimaryVertex();
    AliRsnCutPrimaryVertex(const char *name, Int_t minContributors);
    virtual ~AliRsnCutPrimaryVertex() {;};

    virtual Bool_t IsSelected(AliRsnCut::ETarget tgt, AliRsnDaughter *daughter);
    virtual Bool_t IsSelected(AliRsnCut::ETarget tgt, AliRsnPairParticle *pair);
    virtual Bool_t IsSelected(ETarget tgt, AliRsnEvent *event);
    virtual Bool_t IsSelected(ETarget tgt, AliRsnEvent *ev1, AliRsnEvent *ev2);

  protected:

    ClassDef(AliRsnCutPrimaryVertex, 1)
};

#endif
