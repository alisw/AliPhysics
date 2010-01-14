//
// Class AliRsnCutPrimaryVertex
//
// This cut implementation checks the quality of event primary vertex.
// It currently works only with ESD events (not AOD).
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
    AliRsnCutPrimaryVertex(const char *name, Int_t minContributors, Bool_t acceptTPC = kFALSE);
    virtual ~AliRsnCutPrimaryVertex() {;};

    virtual Bool_t IsSelected(AliRsnCut::ETarget tgt, AliRsnDaughter *daughter);
    virtual Bool_t IsSelected(AliRsnCut::ETarget tgt, AliRsnPairParticle *pair);
    virtual Bool_t IsSelected(ETarget tgt, AliRsnEvent *event);
    virtual Bool_t IsSelected(ETarget tgt, AliRsnEvent *ev1, AliRsnEvent *ev2);

  protected:

    Bool_t fAcceptTPC;  // if kTRUE, the TPC primary vertexes are accepted

    ClassDef(AliRsnCutPrimaryVertex, 1)
};

#endif
