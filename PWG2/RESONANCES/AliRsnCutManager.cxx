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

#include "AliLog.h"

#include "AliRsnCutManager.h"

ClassImp(AliRsnCutManager)

//_____________________________________________________________________________
AliRsnCutManager::AliRsnCutManager() :
  TNamed("defaultName", "defaultTitle"),
  fMotherCuts(0x0)
{
//
// Constructor without arguments.
//

  Int_t i;
  for (i = 0; i < 3; i++) fDaughterCuts[i] = 0x0;
}

//_____________________________________________________________________________
AliRsnCutManager::AliRsnCutManager(const char *name, const char *title) :
  TNamed(name, title),
  fMotherCuts(0x0)
{
//
// Constructor with name and title.
//

  Int_t i;
  for (i = 0; i < 3; i++) fDaughterCuts[i] = 0x0;
}

//_____________________________________________________________________________
AliRsnCutManager::AliRsnCutManager(const AliRsnCutManager &cut) :
  TNamed(cut),
  fMotherCuts(cut.fMotherCuts)
{
//
// Constructor with name and title.
//

  Int_t i;
  for (i = 0; i < 3; i++) fDaughterCuts[i] = cut.fDaughterCuts[i];
}

AliRsnCutManager& AliRsnCutManager::operator=(const AliRsnCutManager &cut)
{
//
// Assign operator
//

  SetName(cut.GetName());
  SetTitle(cut.GetTitle());
  
  Int_t i;
  for (i = 0; i < 3; i++) fDaughterCuts[i] = cut.fDaughterCuts[i];
  
  fMotherCuts = cut.fMotherCuts;
  
  return (*this);
}

//_____________________________________________________________________________
AliRsnCutManager::~AliRsnCutManager()
{
//
// Destructor.
// Deletes all cut definitions.
//

  Int_t i;
  for (i = 0; i < 3; i++) delete fDaughterCuts[i];
  
  delete fMotherCuts;
}

//_____________________________________________________________________________
void AliRsnCutManager::SetEvent(AliRsnEvent *event)
{
//
// Sets reference event in all cut sets
//

  Int_t i;
  for (i = 0; i < 3; i++) if (fDaughterCuts[i]) fDaughterCuts[i]->SetEvent(event);
  
  if (fMotherCuts) fMotherCuts->SetEvent(event);
}
