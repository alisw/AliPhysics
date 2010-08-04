//
// Class AliRsnCutManager
//
// The cut manager: contains a complete set of cut definitions
// to be applied to all possible targets (one for each target),
// in order to ease the set-up procedure of cuts and allow to
// pass them at once to each object which must use them
//
// author: Martin Vala (martin.vala@cern.ch)
//

#ifndef ALIRSNCUTMANAGER_H
#define ALIRSNCUTMANAGER_H

#include <TNamed.h>

#include "AliRsnCut.h"
#include "AliRsnDaughter.h"
#include "AliRsnMother.h"
#include "AliRsnCutSet.h"

class AliRsnCutManager : public TNamed
{
  public:

    AliRsnCutManager();
    AliRsnCutManager(const char *name, const char* title = "");
    AliRsnCutManager(const AliRsnCutManager &cut);
    AliRsnCutManager& operator=(const AliRsnCutManager& cut);
    ~AliRsnCutManager();
    
    AliRsnCutSet*    GetCommonDaughterCuts() {return fDaughterCutsCommon;}
    AliRsnCutSet*    GetDaughter1Cuts() {return fDaughterCuts1;}
    AliRsnCutSet*    GetDaughter2Cuts() {return fDaughterCuts2;}
    AliRsnCutSet*    GetMotherCuts() {return fMotherCuts;}

    void    SetCommonDaughterCuts(AliRsnCutSet *cuts) {fDaughterCutsCommon = cuts;}
    void    SetDaughter1Cuts(AliRsnCutSet *cuts) {fDaughterCuts1 = cuts;}
    void    SetDaughter2Cuts(AliRsnCutSet *cuts) {fDaughterCuts2 = cuts;}
    void    SetMotherCuts(AliRsnCutSet *cuts) {fMotherCuts = cuts;}
    void    SetEvent(AliRsnEvent *event);
    
    Bool_t  PassCommonDaughterCuts(AliRsnDaughter *daughter) const {if (fDaughterCutsCommon) return fDaughterCutsCommon->IsSelected(daughter); return kTRUE;}
    Bool_t  PassDaughter1Cuts(AliRsnDaughter *daughter) const {if (fDaughterCuts1) return fDaughterCuts2->IsSelected(daughter); return kTRUE;}
    Bool_t  PassDaughter2Cuts(AliRsnDaughter *daughter) const {if (fDaughterCuts2) return fDaughterCuts2->IsSelected(daughter); return kTRUE;}
    Bool_t  PassMotherCuts(AliRsnMother *mother) {if (fMotherCuts) return fMotherCuts->IsSelected(mother); return kTRUE;}

  private:

    AliRsnCutSet  *fDaughterCutsCommon; //  single-track cuts common to both daughters
    AliRsnCutSet  *fDaughterCuts1;      //  single-track cuts for only first daughter
    AliRsnCutSet  *fDaughterCuts2;      //  single-track cuts for only second daughter
    AliRsnCutSet  *fMotherCuts;         //  mother cuts (on relations between daughters)

    ClassDef(AliRsnCutManager, 1)  // dictionary
};

#endif
