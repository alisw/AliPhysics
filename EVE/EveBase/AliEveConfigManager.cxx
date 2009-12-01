// $Id$
// Author: Matevz Tadel 2009

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveConfigManager.h"

#include <AliEveMultiView.h>
#include <TEveManager.h>
#include <TEveBrowser.h>
#include <TGFileDialog.h>
#include <TGMenu.h>

class AliEveMacroExecutor;
class TEveProjectionManager;
class TEveGeoShape;
class TEveUtil;

#include <TSystem.h>
#include <TPRegexp.h>
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
   kAEMSave, kAEMLoad, kAEMDefault, kAEMScreen, kAEMProjector, kAEMTransparentDark, kAEMTransparentLight, kAEMTransparentMonoDark, kAEMTransparentMonoLight, kAEMTpc
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
  fAliEvePopup->AddEntry("&Save", kAEMSave);
  fAliEvePopup->AddEntry("&Load", kAEMLoad);

  fAliEvePopup->AddSeparator();

  fAliEvePopup->AddEntry("&Default", kAEMDefault);
  fAliEvePopup->AddEntry("&Screen", kAEMScreen);
  fAliEvePopup->AddEntry("&Projector", kAEMProjector);

  fAliEvePopup->AddSeparator();

  fAliEvePopup->AddEntry("&Transparent screen", kAEMTransparentDark);
  fAliEvePopup->AddEntry("&Transparent projector", kAEMTransparentLight);
  fAliEvePopup->AddEntry("&Transparent mono dark", kAEMTransparentMonoDark);
  fAliEvePopup->AddEntry("&Transparent mono light", kAEMTransparentMonoLight);

  fAliEvePopup->AddSeparator();

  fAliEvePopup->AddEntry("&TPC", kAEMTpc);

  fAliEvePopup->AddSeparator();

  fAliEvePopup->Connect("Activated(Int_t)", "AliEveConfigManager",
                        this, "AliEvePopupHandler(Int_t)");

#if ROOT_VERSION_CODE >= ROOT_VERSION(5,25,4)
  TGMenuBar *mBar = gEve->GetBrowser()->GetMenuBar();
  mBar->AddPopup("&AliEve", fAliEvePopup, new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0));
  gEve->GetBrowser()->GetTopMenuFrame()->Layout();
#else
  // Uber hack as TRootBrowser does not provede manu-bar getter.
  TGFrameElement   *xxFE = (TGFrameElement*)   gEve->GetBrowser()->GetList()->First();
  TGCompositeFrame *xxCF = (TGCompositeFrame*) xxFE->fFrame;
  xxFE = (TGFrameElement*)   xxCF->GetList()->First();
  xxCF = (TGCompositeFrame*) xxFE->fFrame;
  xxFE = (TGFrameElement*)   xxCF->GetList()->First();
  xxCF = (TGCompositeFrame*) xxFE->fFrame;
  xxFE = (TGFrameElement*)   xxCF->GetList()->First();
  TGMenuBar *mBar = (TGMenuBar*) xxFE->fFrame;
  mBar->AddPopup("&AliEve", fAliEvePopup, new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0));
  ((TGCompositeFrame*)mBar->GetParent()->GetParent())->Layout();
#endif
}

//==============================================================================

namespace
{
const char *gMacroSaveAsTypes[] = {"CINT Macro", "*.C",
                                   0, 0};
}

void AliEveConfigManager::AliEvePopupHandler(Int_t id)
{
  // Handle user selections from AliEve popup.

  static const TEveException kEH("AliEveConfigManager::AliEvePopupHandler ");

  switch (id)
  {
    case kAEMSave:
    {
      TGFileInfo fi;
      fi.fFileTypes   = gMacroSaveAsTypes;
      fi.fIniDir      = StrDup(""); // current directory
      fi.fFileTypeIdx = 0;
      fi.fOverwrite   = kTRUE;
      new TGFileDialog(gClient->GetDefaultRoot(), gEve->GetMainWindow(), kFDSave, &fi);
      if (!fi.fFilename) return;

      TPMERegexp filere(".*/([^/]+$");
      if (filere.Match(fi.fFilename) != 2)
      {
        Warning("AliEvePopupHandler", "file '%s' bad.", file.Data());
        return;
      }
      printf("Saving...\n");

      TString file(filere[1]);
      if (!file.EndsWith(".C"))
        file += ".C";
      gSystem->ChangeDirectory(fi.fIniDir);
      gEve->SaveVizDB(file);

      printf("Saved\n");
      break;

    }

    case kAEMLoad:
    {

      printf("Loading...\n");

      gEve->LoadVizDB("VizDB_custom.C", kTRUE, kTRUE);

      printf("Loaded\n");
      break;

    }

    case kAEMDefault:
    {

      gEve->GetScenes()->FirstChild()->DestroyElements();

      gEve->GetScenes()->FindChild("RPhi Geometry")->LastChild()->DestroyElements();

      gEve->GetScenes()->FindChild("RhoZ Geometry")->LastChild()->DestroyElements();

      gEve->LoadVizDB("geom_gentle_default.C", kTRUE, kTRUE);

      gEve->LoadVizDB("VizDB_scan.C", kTRUE, kTRUE);

      if(gEve->GetViewers()->UseLightColorSet())
        gEve->GetViewers()->SwitchColorSet();

      gEve->FullRedraw3D();   
      
      break;
    }

    case kAEMScreen:
    {

      gEve->GetScenes()->FirstChild()->DestroyElements();

      gEve->GetScenes()->FindChild("RPhi Geometry")->LastChild()->DestroyElements();

      gEve->GetScenes()->FindChild("RhoZ Geometry")->LastChild()->DestroyElements();

      gEve->LoadVizDB("geom_gentle.C", kTRUE, kTRUE);

      gEve->LoadVizDB("VizDB_scan_screen.C", kTRUE, kTRUE);

      if(gEve->GetViewers()->UseLightColorSet())
        gEve->GetViewers()->SwitchColorSet();

      gEve->FullRedraw3D();    

      break;
    }

    case kAEMProjector:
    {

      gEve->GetScenes()->FirstChild()->DestroyElements();

      gEve->GetScenes()->FindChild("RPhi Geometry")->LastChild()->DestroyElements();

      gEve->GetScenes()->FindChild("RhoZ Geometry")->LastChild()->DestroyElements();

      gEve->LoadVizDB("geom_gentle.C", kTRUE, kTRUE);

      gEve->LoadVizDB("VizDB_scan_projector.C", kTRUE, kTRUE);

      if(!gEve->GetViewers()->UseLightColorSet())
        gEve->GetViewers()->SwitchColorSet();

      gEve->FullRedraw3D();    

      break;
    }

    case kAEMTransparentDark:
    {

      gEve->GetScenes()->FirstChild()->DestroyElements();

      gEve->GetScenes()->FindChild("RPhi Geometry")->LastChild()->DestroyElements();

      gEve->GetScenes()->FindChild("RhoZ Geometry")->LastChild()->DestroyElements();

      gEve->LoadVizDB("geom_gentle_transparent.C", kTRUE, kTRUE);

      gEve->LoadVizDB("VizDB_scan_transparentdark.C", kTRUE, kTRUE);

      if(gEve->GetViewers()->UseLightColorSet())
        gEve->GetViewers()->SwitchColorSet();

      gEve->FullRedraw3D();   

      break;
    }

    case kAEMTransparentLight:
    {

      gEve->GetScenes()->FirstChild()->DestroyElements();

      gEve->GetScenes()->FindChild("RPhi Geometry")->LastChild()->DestroyElements();

      gEve->GetScenes()->FindChild("RhoZ Geometry")->LastChild()->DestroyElements();

      gEve->LoadVizDB("geom_gentle_transparent.C", kTRUE, kTRUE);

      gEve->LoadVizDB("VizDB_scan_transparentlight.C", kTRUE, kTRUE);

      if(!gEve->GetViewers()->UseLightColorSet())
        gEve->GetViewers()->SwitchColorSet();

      gEve->FullRedraw3D();    

      break;
    }

    case kAEMTransparentMonoDark:
    {

      gEve->GetScenes()->FirstChild()->DestroyElements();

      gEve->GetScenes()->FindChild("RPhi Geometry")->LastChild()->DestroyElements();

      gEve->GetScenes()->FindChild("RhoZ Geometry")->LastChild()->DestroyElements();

      gEve->LoadVizDB("geom_gentle_transparentdark.C", kTRUE, kTRUE);

      gEve->LoadVizDB("VizDB_scan_transparentdark.C", kTRUE, kTRUE);

      if(gEve->GetViewers()->UseLightColorSet())
        gEve->GetViewers()->SwitchColorSet();

      gEve->FullRedraw3D();   

      break;
    }

    case kAEMTransparentMonoLight:
    {

      gEve->GetScenes()->FirstChild()->DestroyElements();

      gEve->GetScenes()->FindChild("RPhi Geometry")->LastChild()->DestroyElements();

      gEve->GetScenes()->FindChild("RhoZ Geometry")->LastChild()->DestroyElements();

      gEve->LoadVizDB("geom_gentle_transparentlight.C", kTRUE, kTRUE);

      gEve->LoadVizDB("VizDB_scan_transparentlight.C", kTRUE, kTRUE);

      if(!gEve->GetViewers()->UseLightColorSet())
        gEve->GetViewers()->SwitchColorSet();

      gEve->FullRedraw3D();    

      break;
    }

    case kAEMTpc:
    {

      gEve->GetScenes()->FirstChild()->DestroyElements();

      gEve->GetScenes()->FindChild("RPhi Geometry")->LastChild()->DestroyElements();

      gEve->GetScenes()->FindChild("RhoZ Geometry")->LastChild()->DestroyElements();

      gEve->LoadVizDB("geom_gentle_tpc.C", kTRUE, kTRUE);

      gEve->LoadVizDB("VizDB_scan_tpc.C", kTRUE, kTRUE);

      if(!gEve->GetViewers()->UseLightColorSet())
        gEve->GetViewers()->SwitchColorSet();

      gEve->FullRedraw3D();

      break;
    }

    default:
    {
      Warning(kEH, "Unknown menu entry.");
      break;
    }
  }
}
