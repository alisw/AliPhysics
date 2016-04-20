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
class AliVVertex;

class AliRsnCutPrimaryVertex : public AliRsnCut {
public:

   AliRsnCutPrimaryVertex(const char *name = "cutPrimVert", Double_t maxVz = 10.0, Int_t minContributors = 1, Bool_t acceptTPC = kFALSE, Bool_t acceptSPD = kTRUE);
   virtual ~AliRsnCutPrimaryVertex() {;};

   void           SetCheckPileUp(Bool_t doit = kTRUE) {fCheckPileUp = doit;}
   virtual Bool_t IsSelected(TObject *object);
   virtual void   Print(const Option_t *option = "") const;

   Bool_t         GetCheckZResolutionSPD(){return fCheckZResolutionSPD;}
   void           SetCheckZResolutionSPD(Bool_t doit=kTRUE){fCheckZResolutionSPD=doit;}
   Float_t        GetMaxZResolutionSPD(){return fMaxZResolutionSPD;}
   void           SetMaxZResolutionSPD(Float_t val=0.25){fMaxZResolutionSPD=val;}
   Bool_t         GoodZResolutionSPD(const AliESDVertex* v);
   Bool_t         GoodZResolutionSPD(AliAODVertex* v);
   Bool_t         GoodZResolutionSPD();

   Bool_t         GetCheckDispersionSPD(){return fCheckDispersionSPD;}
   void           SetCheckDispersionSPD(Bool_t doit=kTRUE){fCheckDispersionSPD=doit;}
   Float_t        GetMaxDispersionSPD(){return fMaxDispersionSPD;}
   void           SetMaxDispersionSPD(Float_t val=0.04){fMaxDispersionSPD=val;}
   Bool_t         GoodDispersionSPD(const AliESDVertex* v);
   Bool_t         GoodDispersionSPD();

   Bool_t         GetCheckZDifferenceSPDTrack(){return fCheckZDifferenceSPDTrack;}
   void           SetCheckZDifferenceSPDTrack(Bool_t doit=kTRUE){fCheckZDifferenceSPDTrack=doit;}
   Float_t        GetMaxZDifferenceSPDTrack(){return fMaxZDifferenceSPDTrack;}
   void           SetMaxZDifferenceSPDTrack(Float_t val=0.5){fMaxZDifferenceSPDTrack=val;}
   Bool_t         GoodZDifferenceSPDTrack(const AliESDVertex* vTrk, const AliESDVertex* vSPD);
   Bool_t         GoodZDifferenceSPDTrack(AliAODVertex* vTrk, AliAODVertex* vSPD);
   Bool_t         GoodZDifferenceSPDTrack();

protected:

   Bool_t         CheckVertex(AliVVertex *vert);

   Bool_t         fAcceptTPC;   // if kTRUE, the TPC primary vertexes are accepted
   Bool_t         fAcceptSPD;   // if kTRUE, the SPD primary vertexes are accepted
   Bool_t         fCheckPileUp; // check and reject pileupped events (pp)
   Bool_t         fCheckZResolutionSPD;// check the z resolution of SPD vertex
   Float_t        fMaxZResolutionSPD;// maximum z resolution of SPD vertex (ESD only)
   Bool_t         fCheckDispersionSPD;// check the dispersion of SPD vertex
   Float_t        fMaxDispersionSPD;// maximum dispersion of SPD vertex (ESD only)
   Bool_t         fCheckZDifferenceSPDTrack;// check the z position difference between track and SPD vertices (ESD only)
   Float_t        fMaxZDifferenceSPDTrack;// maximum z position difference between track and SPD vertices
   ClassDef(AliRsnCutPrimaryVertex, 2)
};

#endif
