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
  fDaughterCutsCommon(0x0),
  fDaughterCuts1(0x0),
  fDaughterCuts2(0x0),
  fMotherCuts(0x0)
{
//
// Constructor without arguments.
//
}

//_____________________________________________________________________________
AliRsnCutManager::AliRsnCutManager(const char *name, const char *title) :
  TNamed(name, title),
  fDaughterCutsCommon(0x0),
  fDaughterCuts1(0x0),
  fDaughterCuts2(0x0),
  fMotherCuts(0x0)
{
//
// Constructor with name and title.
//
}

//_____________________________________________________________________________
AliRsnCutManager::AliRsnCutManager(const AliRsnCutManager &cut) :
  TNamed(cut),
  fDaughterCutsCommon(cut.fDaughterCutsCommon),
  fDaughterCuts1(cut.fDaughterCuts1),
  fDaughterCuts2(cut.fDaughterCuts2),
  fMotherCuts(cut.fMotherCuts)
{
//
// Constructor with name and title.
//
}

AliRsnCutManager& AliRsnCutManager::operator=(const AliRsnCutManager &cut)
{
//
// Assign operator
//

  SetName(cut.GetName());
  SetTitle(cut.GetTitle());
  
  fDaughterCuts2      = cut.fDaughterCuts2;
  fDaughterCuts1      = cut.fDaughterCuts1;
  fDaughterCutsCommon = cut.fDaughterCutsCommon;
  fMotherCuts         = cut.fMotherCuts;
  
  return (*this);
}

//_____________________________________________________________________________
AliRsnCutManager::~AliRsnCutManager()
{
//
// Destructor.
// Deletes all cut definitions.
//

  delete fDaughterCuts2;
  delete fDaughterCuts1;
  delete fDaughterCutsCommon;
  delete fMotherCuts;
}

//_____________________________________________________________________________
void AliRsnCutManager::SetEvent(AliRsnEvent *event)
{
//
// Sets reference event in all cut sets
//

  if (fDaughterCuts2     ) fDaughterCuts2      ->SetEvent(event);
  if (fDaughterCuts1     ) fDaughterCuts1      ->SetEvent(event);
  if (fDaughterCutsCommon) fDaughterCutsCommon ->SetEvent(event);
  if (fMotherCuts        ) fMotherCuts         ->SetEvent(event);
}
