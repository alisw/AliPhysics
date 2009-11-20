// $Id$
// Author: Matevz Tadel 2009

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveConfigManager.h"

#include <TEveManager.h>
#include <TEveBrowser.h>
#include <TGMenu.h>

#include <RVersion.h>

//______________________________________________________________________________
// Full description of AliEveConfigManager
//

ClassImp(AliEveConfigManager)

AliEveConfigManager* AliEveConfigManager::fgMaster = 0;

namespace
{
 enum EAliEveMenu_e
 {
   kAEMTest
 };
}
 
//______________________________________________________________________________
AliEveConfigManager* AliEveConfigManager::InitializeMaster()
{
  // Get main instance.

  static const TEveException kEH("AliEveConfigManager::InitializeMaster ");

  if (fgMaster)
    throw kEH + "Master already initialized.";

  fgMaster = new AliEveConfigManager;
  return fgMaster;
}

//______________________________________________________________________________
AliEveConfigManager* AliEveConfigManager::GetMaster()
{
  // Get main instance.

  static const TEveException kEH("AliEveConfigManager::GetMaster ");

  if (!fgMaster)
    throw kEH + "Master not initialized.";

  return fgMaster;
}

//______________________________________________________________________________
AliEveConfigManager::AliEveConfigManager() :
  TObject(),
  fAliEvePopup(0)
{
  // Constructor.
  // Expected TEveManager is already initialized.

  fAliEvePopup = new TGPopupMenu(gClient->GetRoot());
  fAliEvePopup->AddEntry("&Test", kAEMTest);

  fAliEvePopup->AddSeparator();

  fAliEvePopup->Connect("Activated(Int_t)", "AliEveConfigManager",
                        this, "AliEvePopupHandler(Int_t)");

  TGMenuBar *mBar = 0;
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,25,4) || defined XXX_LATEST_ROOT
  mBar = gEve->GetBrowser()->GetMenuBar();
#else
  // Uber hack as TRootBrowser does not provede manu-bar getter.
  TGFrameElement   *xxFE = (TGFrameElement*)   gEve->GetBrowser()->GetList()->First();
  TGCompositeFrame *xxCF = (TGCompositeFrame*) xxFE->fFrame;
  xxFE = (TGFrameElement*)   xxCF->GetList()->First();
  xxCF = (TGCompositeFrame*) xxFE->fFrame;
  xxFE = (TGFrameElement*)   xxCF->GetList()->First();
  xxCF = (TGCompositeFrame*) xxFE->fFrame;
  xxFE = (TGFrameElement*)   xxCF->GetList()->First();
  mBar = (TGMenuBar*) xxFE->fFrame;
#endif
  mBar->AddPopup("&AliEve", fAliEvePopup, new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0));
  ((TGCompositeFrame*)mBar->GetParent()->GetParent())->Layout();
}

//==============================================================================

void AliEveConfigManager::AliEvePopupHandler(Int_t id)
{
  // Handle user selections from AliEve popup.

  static const TEveException kEH("AliEveConfigManager::AliEvePopupHandler ");

  switch (id)
  {
    case kAEMTest:
    {
      printf("Test!\n");
      break;
    }

    default:
    {
      Warning(kEH, "Unknown menu entry.");
      break;
    }
  }
}
