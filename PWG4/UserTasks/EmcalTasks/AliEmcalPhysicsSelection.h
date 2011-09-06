#ifndef ALIEMCALPHYSICSSELECTION_H
#define ALIEMCALPHYSICSSELECTION_H

// $Id$

#include "AliPhysicsSelection.h"

class AliEmcalPhysicsSelection: public AliPhysicsSelection
{
 public:
  AliEmcalPhysicsSelection() : AliPhysicsSelection(), fExcludeFastOnly(0) {;}
  virtual ~AliEmcalPhysicsSelection() {;}

  virtual UInt_t GetSelectionMask(const TObject* obj);

  void           SetExcludeFastOnly(Bool_t b) { fExcludeFastOnly = b; }

 protected:
  Bool_t         fExcludeFastOnly; //=true then exclude FastOnly events (only for LHC11a)

  ClassDef(AliEmcalPhysicsSelection, 1); // My physics selection
};
#endif
