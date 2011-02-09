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
class AliRsnCutPrimaryVertex : public AliRsnCut {
public:

   AliRsnCutPrimaryVertex(const char *name = "cutPrimVert", Double_t maxVz = 10.0, Int_t minContributors = 1, Bool_t acceptTPC = kFALSE);
   virtual ~AliRsnCutPrimaryVertex() {;};

   void           SetCheckPileUp(Bool_t doit = kTRUE) {fCheckPileUp = doit;}
   virtual Bool_t IsSelected(TObject *object);
   virtual void   Print(const Option_t *option = "") const;

protected:

   Bool_t fAcceptTPC;   // if kTRUE, the TPC primary vertexes are accepted
   Bool_t fCheckPileUp; // check and reject pileupped events (pp)

   ClassDef(AliRsnCutPrimaryVertex, 1)
};

#endif
