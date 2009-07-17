//
// Class AliRsnCutMgr
//
// The cut manager: contains a complete set of cut definitions
// to be applied to all possible targets (one for each target),
// in order to ease the set-up procedure of cuts and allow to
// pass them at once to each object which must use them
//
// author: Martin Vala (martin.vala@cern.ch)
//

#ifndef ALIRSNCUTMGR_H
#define ALIRSNCUTMGR_H

#include <TNamed.h>

#include "AliRsnCut.h"

class AliRsnCutSet;

class AliRsnCutMgr : public TNamed
{
  public:

    AliRsnCutMgr();
    AliRsnCutMgr(const char *name, const char* title);
    ~AliRsnCutMgr();

    void          SetCutSet(AliRsnCut::ETarget type, AliRsnCutSet*const cutset);
    AliRsnCutSet* GetCutSet(AliRsnCut::ETarget type) {return fCutSets[type];}
    Bool_t        IsSelected(AliRsnCut::ETarget type, TObject *const obj);

  private:

    // dummy constructors
    AliRsnCutMgr(const AliRsnCutMgr &cut) : TNamed(cut) {}
    AliRsnCutMgr& operator=(const AliRsnCutMgr& /*cut*/) {return *this;}

    AliRsnCutSet *fCutSets[AliRsnCut::kLastCutTarget];  // cut definitions for all targets

    ClassDef(AliRsnCutMgr, 1)  // dictionary
};

#endif
