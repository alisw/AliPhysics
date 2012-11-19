#ifndef ALIRSNCUTTOFMATCH_H
#define ALIRSNCUTTOFMATCH_H

//
// Class for TOF-matching cut.
// Author: Francesca Bellini (fbellini@cern.ch)

#include <TMath.h>
#include <TClonesArray.h>

#include "AliESDtrack.h"
#include "AliRsnCut.h"

class AliVTrack;

class AliRsnCutTOFMatch : public AliRsnCut {
public:

   AliRsnCutTOFMatch();
   AliRsnCutTOFMatch(const char *name);
   virtual ~AliRsnCutTOFMatch() { }

   Bool_t   MatchTOF(const AliVTrack *vtrack) const;
   Bool_t   IsSelected(TObject *object);

   ClassDef(AliRsnCutTOFMatch, 1)

};

#endif
