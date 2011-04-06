#ifndef ALIRSNPIDMANAGER_H
#define ALIRSNPIDMANAGER_H

#include "TNamed.h"
#include <TObjArray.h>

class AliRsnDaughter;
class AliRsnCutSet;
class TEntryList;
class AliRsnEvent;
class AliRsnPIDManager : public TNamed {

public:
   AliRsnPIDManager(const char *name = "name", const char *title = "title");
   virtual ~AliRsnPIDManager();

   virtual void Print(Option_t* option = "") const;
   void Init();
   void Reset();

   void SetCut(AliRsnCutSet* cut, const Int_t& pidId);
   TEntryList *GetParticles(const Int_t &charge, const Int_t &pidId);

   void ApplyCuts(AliRsnEvent*ev);

private:

   TObjArray fIdParticles[2];
   TObjArray fCuts;


   void CheckTrack(AliRsnDaughter* d, const Int_t& id);

   AliRsnPIDManager(const AliRsnPIDManager& copy);
   AliRsnPIDManager &operator=(const AliRsnPIDManager &copy);

   ClassDef(AliRsnPIDManager, 1)
};

#endif
