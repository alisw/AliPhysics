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

#include "AliLog.h"

#include "AliRsnCut.h"
#include "AliRsnCutManager.h"

ClassImp(AliRsnCutManager)

//_____________________________________________________________________________
AliRsnCutManager::AliRsnCutManager() :
  TNamed("defaultName", "defaultTitle"),
  fDaughterCutsCommon("defaultCommon", AliRsnTarget::kDaughter),
  fDaughterCuts1("defaultD1", AliRsnTarget::kDaughter),
  fDaughterCuts2("defaultD2", AliRsnTarget::kDaughter),
  fMotherCuts("defaultPair", AliRsnTarget::kMother)
{
//
// Constructor without arguments.
//
}

//_____________________________________________________________________________
AliRsnCutManager::AliRsnCutManager(const char *name, const char *title) :
  TNamed(name, title),
  fDaughterCutsCommon(Form("common_%s", name), AliRsnTarget::kDaughter),
  fDaughterCuts1(Form("d1_%s", name), AliRsnTarget::kDaughter),
  fDaughterCuts2(Form("d2_%s", name), AliRsnTarget::kDaughter),
  fMotherCuts(Form("pair_%s", name), AliRsnTarget::kMother)
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

//_____________________________________________________________________________
AliRsnCutManager& AliRsnCutManager::operator=(const AliRsnCutManager &cut)
{
//
// Assign operator
//

  SetName(cut.GetName());
  SetTitle(cut.GetTitle());
  
  fDaughterCuts2 = cut.fDaughterCuts2;
  fDaughterCuts1 = cut.fDaughterCuts1;
  fDaughterCutsCommon = cut.fDaughterCutsCommon;
  fMotherCuts = cut.fMotherCuts;
  
  return (*this);
}

//_____________________________________________________________________________
AliRsnCutManager::~AliRsnCutManager()
{
//
// Destructor.
// Does nothing.
//
}
