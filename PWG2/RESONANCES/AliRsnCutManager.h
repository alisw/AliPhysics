//
// *** Class AliRsnCutManager ***
//
// This class is used both in normal analysis and efficiency computation
// as a collection of all cuts which could be needed in a single job.
// It allocates an AliRsnCutSet for each possible target:
//  - one with all cuts common to all tracks
//  - one with all cuts for first candidate daughter (definition #1 in pairDef)
//  - one with all cuts for second candidate daughter (definition #2 in pairDef)
//  - one with all cuts on the pair
// -----
// This object is used to define a step in efficiency CORRFW container
// and also is contained in all AliRsnPair objects to decide if two candidates
// can be accepted or not.
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@cern.ch)
//

#ifndef ALIRSNCUTMANAGER_H
#define ALIRSNCUTMANAGER_H

#include <TNamed.h>

#include "AliRsnDaughter.h"
#include "AliRsnMother.h"
#include "AliRsnCutSet.h"

class AliRsnCut;

class AliRsnCutManager : public TNamed {
public:

   AliRsnCutManager();
   AliRsnCutManager(const char *name, const char* title = "");
   AliRsnCutManager(const AliRsnCutManager &cut);
   AliRsnCutManager& operator=(const AliRsnCutManager& cut);
   ~AliRsnCutManager();

   AliRsnCutSet*  GetCommonDaughterCuts() {return &fDaughterCutsCommon;}
   AliRsnCutSet*  GetDaughter1Cuts()      {return &fDaughterCuts1;}
   AliRsnCutSet*  GetDaughter2Cuts()      {return &fDaughterCuts2;}
   AliRsnCutSet*  GetMotherCuts()         {return &fMotherCuts;}

   Bool_t         IsSelected(TObject *object);
   Bool_t         PassCommonDaughterCuts(AliRsnDaughter *daughter) {return fDaughterCutsCommon.IsSelected(daughter);}
   Bool_t         PassDaughter1Cuts(AliRsnDaughter *daughter)      {return fDaughterCuts1.IsSelected(daughter);}
   Bool_t         PassDaughter2Cuts(AliRsnDaughter *daughter)      {return fDaughterCuts2.IsSelected(daughter);}
   Bool_t         PassMotherCuts(AliRsnMother *mother)             {return fMotherCuts.IsSelected(mother);}
   Bool_t         PassSpecificDaughterCuts(Bool_t first, AliRsnDaughter *daughter)
   {if (first) return PassDaughter1Cuts(daughter); else return PassDaughter2Cuts(daughter);}

private:

   AliRsnCutSet  fDaughterCutsCommon; // single-track cuts common to both daughters
   AliRsnCutSet  fDaughterCuts1;      // single-track cuts for only first daughter
   AliRsnCutSet  fDaughterCuts2;      // single-track cuts for only second daughter
   AliRsnCutSet  fMotherCuts;         // mother cuts (on relations between daughters)

   ClassDef(AliRsnCutManager, 2)      // dictionary
};

inline Bool_t AliRsnCutManager::IsSelected(TObject *object)
{
//
// Check all selection cuts
//

   if (object->InheritsFrom(AliRsnDaughter::Class())) {
      return PassCommonDaughterCuts((AliRsnDaughter*)object);
   } else if (object->InheritsFrom(AliRsnMother::Class())) {
      AliRsnMother *mother = (AliRsnMother*)object;
      if (!PassCommonDaughterCuts(mother->GetDaughter(0))) return kFALSE;
      if (!PassCommonDaughterCuts(mother->GetDaughter(1))) return kFALSE;
      if (!PassDaughter1Cuts(mother->GetDaughter(0))) return kFALSE;
      if (!PassDaughter2Cuts(mother->GetDaughter(1))) return kFALSE;
      if (!PassMotherCuts(mother)) return kFALSE;
      return kTRUE;
   } else {
      AliError("AliRsnCutManager can check only daughters and mothers");
      return kFALSE;
   }
}

#endif
