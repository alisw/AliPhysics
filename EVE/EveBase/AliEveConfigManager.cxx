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
#include <TEveWindow.h>
#include <TGFileDialog.h>
#include <TGMenu.h>

#include "AliEveEventManager.h"
#include "AliEveMacroExecutor.h"
#include "AliEveMacro.h"
#include "AliEvePreferencesWindow.h"
#include "AliEveCutsWindow.h"
#include "AliEveMomentumVectors.h"

//Storage Manager:
#ifdef ZMQ
#include "AliStorageAdministratorPanelListEvents.h"
#include "AliStorageAdministratorPanelMarkEvent.h"
#endif

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
	 kAEMDefault, kAEMScreen, kAEMProjector, kAEMNotransparency, kAEMTransparentDark, kAEMTransparentLight, kAEMTransparentMonoDark, kAEMTransparentMonoLight, kAEMGreen, kAEMBright, kAEMYellow, kAEMTpc, kAEMAll, kAEM3d, kAEMRphi, kAEMRhoz, kAEMAllhr, kAEM3dhr, kAEMRphihr, kAEMRhozhr, kAEMSave, kAEMOpen, kAEMSetDefault, kAEMResiduals,  kAEMCuts, kAEMVectors, kAEMGui, kStorageListEvents, kStorageMarkEvent, kPreferences
 };
}
 
//______________________________________________________________________________
AliEveConfigManager* AliEveConfigManager::InitializeMaster(bool storageManager)
{
  // Get main instance.

  static const TEveException kEH("AliEveConfigManager::InitializeMaster ");

  if (fgMaster)
    throw kEH + "Master already initialized.";

  fgMaster = new AliEveConfigManager(storageManager);
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
AliEveConfigManager::AliEveConfigManager(bool storageManager) :
  TObject(),
  fAnalysisPopup(0),
  fAliEvePopup(0),
  fAliEveGeometries(0),
  fAliEvePictures(0),
  fAliEvePicturesHR(0),
  fAliEveVizDBs(0),
  fLoadCheck(kFALSE)
{
  // Constructor.
  // Expected TEveManager is already initialized.


  fAliEveGeometries = new TGPopupMenu(gClient->GetRoot());
  fAliEveGeometries->AddEntry("&Default", kAEMDefault);
  fAliEveGeometries->AddEntry("&Screen", kAEMScreen);
  fAliEveGeometries->AddEntry("&Projector", kAEMProjector);

  fAliEveGeometries->AddSeparator();

  fAliEveGeometries->AddEntry("&Low transparency", kAEMNotransparency);

  fAliEveGeometries->AddSeparator();

  fAliEveGeometries->AddEntry("&Transparent screen", kAEMTransparentDark);
  fAliEveGeometries->AddEntry("&Transparent projector", kAEMTransparentLight);
  fAliEveGeometries->AddEntry("&Transparent mono dark", kAEMTransparentMonoDark);
  fAliEveGeometries->AddEntry("&Transparent mono light", kAEMTransparentMonoLight);

  fAliEveGeometries->AddSeparator();

  fAliEveGeometries->AddEntry("&First collision setup", kAEMGreen);
  fAliEveGeometries->AddEntry("&Bright", kAEMBright);

  fAliEveGeometries->AddSeparator();

  fAliEveGeometries->AddEntry("&TPC Yellow", kAEMYellow);
  fAliEveGeometries->AddEntry("&TPC Blue", kAEMTpc);

  fAliEveGeometries->AddSeparator();

  fAliEvePictures = new TGPopupMenu(gClient->GetRoot());

  fAliEvePictures->AddEntry("&Save all views", kAEMAll);
  fAliEvePictures->AddEntry("&Save 3D View",   kAEM3d);
  fAliEvePictures->AddEntry("&Save RPhi View", kAEMRphi);
  fAliEvePictures->AddEntry("&Save RhoZ View", kAEMRhoz);

  fAliEvePictures->AddSeparator();

  fAliEvePicturesHR = new TGPopupMenu(gClient->GetRoot());

  fAliEvePicturesHR->AddEntry("&Save all views HR", kAEMAllhr);
  fAliEvePicturesHR->AddEntry("&Save 3D View HR",   kAEM3dhr);
  fAliEvePicturesHR->AddEntry("&Save RPhi View HR", kAEMRphihr);
  fAliEvePicturesHR->AddEntry("&Save RhoZ View HR", kAEMRhozhr);

  fAliEvePicturesHR->AddSeparator();

  fAliEveVizDBs = new TGPopupMenu(gClient->GetRoot());

  fAliEveVizDBs->AddEntry("&Save VizDB", kAEMSave);
  fAliEveVizDBs->AddEntry("&Load VizDB", kAEMOpen);

  fAliEveVizDBs->AddSeparator();

  fAliEvePopup = new TGPopupMenu(gClient->GetRoot());
    fAliEvePopup->AddEntry("&Preferences", kPreferences);
    fAliEvePopup->AddSeparator();
  fAliEvePopup->AddEntry("&Set Default Settings", kAEMSetDefault);
  fAliEvePopup->AddSeparator();
  fAliEvePopup->AddPopup("&Geometries/VizDBs", fAliEveGeometries);
  fAliEvePopup->AddSeparator();
  fAliEvePopup->AddPopup("&Pictures", fAliEvePictures);
  fAliEvePopup->AddSeparator();
  fAliEvePopup->AddPopup("&PicturesHR", fAliEvePicturesHR);
  fAliEvePopup->AddSeparator();
  fAliEvePopup->AddPopup("&VizDBs",  fAliEveVizDBs);
  fAliEvePopup->AddSeparator();

  fAnalysisPopup = new TGPopupMenu(gClient->GetRoot());  
//  fAnalysisPopup->AddEntry("&Residuals", kAEMResiduals);
//  fAnalysisPopup->AddSeparator();
  fAnalysisPopup->AddEntry("&Cuts", kAEMCuts);
  fAnalysisPopup->AddSeparator();
  fAnalysisPopup->AddEntry("&Momentum vectors", kAEMVectors);
  fAnalysisPopup->AddSeparator();
//  fAnalysisPopup->AddEntry("&Gui Mode", kAEMGui);
//  fAnalysisPopup->AddSeparator();

  fAliEvePopup->Connect("Activated(Int_t)", "AliEveConfigManager",
                        this, "AliEvePopupHandler(Int_t)");

  fAnalysisPopup->Connect("Activated(Int_t)", "AliEveConfigManager",
                        this, "AliEvePopupHandler(Int_t)");

  //Storage Manager:
#ifdef ZMQ
    if(storageManager)
    {
        gEve->GetBrowser()->StartEmbedding(0);
        AliStorageAdministratorPanelListEvents* fListEventsTab = AliStorageAdministratorPanelListEvents::GetInstance();
        gEve->GetBrowser()->StopEmbedding("List");

        fListEventsTab->Connect("SelectedEvent()","AliEveConfigManager",this,"SetEventInEventManager()");
    }
#endif
  
  fLoadCheck = kFALSE;

#if ROOT_VERSION_CODE >= ROOT_VERSION(5,25,4)
  TGMenuBar *mBar = gEve->GetBrowser()->GetMenuBar();
  mBar->AddPopup("&AliEve", fAliEvePopup, new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0));
  ((TGCompositeFrame*)mBar->GetParent()->GetParent())->Layout();
  mBar->AddPopup("&Tools", fAnalysisPopup, new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0));
    ((TGCompositeFrame*)mBar->GetParent()->GetParent())->Layout();
  
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
  mBar->AddPopup("&Tools", fAnalysisPopup, new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0));
  ((TGCompositeFrame*)mBar->GetParent()->GetParent())->Layout();
  mBar->AddPopup("&Storage Manager",fStoragePopup,new TGLayoutHints(kLHintsTop | kLHintsLeft,0,4,0,0));
    ((TGCompositeFrame*)mBar->GetParent()->GetParent())->Layout();
#endif
}

//==============================================================================

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

namespace
{
const char *gMacroSaveAsTypes[] = {"CINT Macro", "*.C",
                                   0, 0}; //for saving/loading macros

const char *gPictureSaveAsTypes[] = {"PNG Image", "*.png",
                                   0, 0}; //for saving pictures

}

void AliEveConfigManager::ConnectEventManagerSignals()
{
    AliEveEventManager *manager = AliEveEventManager::GetMaster();
    manager->Connect("StorageManagerOk()","AliEveConfigManager",this,"StorageManagerChangedState(=1)");
    manager->Connect("StorageManagerDown()","AliEveConfigManager",this,"StorageManagerChangedState(=0)");
}

void AliEveConfigManager::AliEvePopupHandler(Int_t id)
{
  // Handle user selections from AliEve popup.

  static const TEveException kEH("AliEveConfigManager::AliEvePopupHandler ");
    
  switch (id)
  {

    case kAEMDefault: //default geometry and VizDB
    {
      AliEveMultiView *mv = AliEveMultiView::Instance();

      mv->DestroyAllGeometries(); //destroy RPhi and Rhoz geometries before putting new

      gEve->LoadVizDB("geom_gentle_default.C", kTRUE, kTRUE); //loading geometry

      gEve->LoadVizDB("VizDB_scan.C", kTRUE, kTRUE); //loading VizDB

      if(!gEve->GetViewers()->UseLightColorSet())
        gEve->GetViewers()->SwitchColorSet(); //white background

      gEve->FullRedraw3D();   
      
      break;
    }

    case kAEMScreen: //default geometry with black background
    {
      AliEveMultiView *mv = AliEveMultiView::Instance();

      mv->DestroyAllGeometries();

      gEve->LoadVizDB("geom_gentle_default.C", kTRUE, kTRUE);

      gEve->LoadVizDB("VizDB_scan_screen.C", kTRUE, kTRUE);

      if(gEve->GetViewers()->UseLightColorSet())
        gEve->GetViewers()->SwitchColorSet();

      gEve->FullRedraw3D();    

      break;
    }

    case kAEMProjector: //default geometry with white background
    {
      AliEveMultiView *mv = AliEveMultiView::Instance();

      mv->DestroyAllGeometries();

      gEve->LoadVizDB("geom_gentle_projector.C", kTRUE, kTRUE);

      gEve->LoadVizDB("VizDB_scan_projector.C", kTRUE, kTRUE);

      if(!gEve->GetViewers()->UseLightColorSet())
        gEve->GetViewers()->SwitchColorSet();

      gEve->FullRedraw3D();    

      break;
    }

    case kAEMNotransparency: //default geometry with low transparency (5%)
    {
      AliEveMultiView *mv = AliEveMultiView::Instance();

      mv->DestroyAllGeometries();

      gEve->LoadVizDB("geom_gentle_notransparency.C", kTRUE, kTRUE);

      gEve->LoadVizDB("VizDB_scan.C", kTRUE, kTRUE);

      if(!gEve->GetViewers()->UseLightColorSet())
        gEve->GetViewers()->SwitchColorSet();

      gEve->FullRedraw3D();    

      break;
    }


    case kAEMTransparentDark: //default geometry with black background, high transparency (80%)
    {
      AliEveMultiView *mv = AliEveMultiView::Instance();

      mv->DestroyAllGeometries();

      gEve->LoadVizDB("geom_gentle_transparent.C", kTRUE, kTRUE);

      gEve->LoadVizDB("VizDB_scan_transparentdark.C", kTRUE, kTRUE);

      if(gEve->GetViewers()->UseLightColorSet())
        gEve->GetViewers()->SwitchColorSet();

      gEve->FullRedraw3D();   

      break;
    }

    case kAEMTransparentLight: //default geometry with white background, high transparency (80%)
    {
      AliEveMultiView *mv = AliEveMultiView::Instance();

      mv->DestroyAllGeometries();

      gEve->LoadVizDB("geom_gentle_transparent.C", kTRUE, kTRUE);

      gEve->LoadVizDB("VizDB_scan_transparentlight.C", kTRUE, kTRUE);

      if(!gEve->GetViewers()->UseLightColorSet())
        gEve->GetViewers()->SwitchColorSet();

      gEve->FullRedraw3D();    

      break;
    }

    case kAEMTransparentMonoDark:
    {
      AliEveMultiView *mv = AliEveMultiView::Instance();

      mv->DestroyAllGeometries();

      gEve->LoadVizDB("geom_gentle_transparentdark.C", kTRUE, kTRUE);

      gEve->LoadVizDB("VizDB_scan_transparentdark.C", kTRUE, kTRUE);

      if(gEve->GetViewers()->UseLightColorSet())
        gEve->GetViewers()->SwitchColorSet();

      gEve->FullRedraw3D();   

      break;
    }

    case kAEMTransparentMonoLight:
    {
      AliEveMultiView *mv = AliEveMultiView::Instance();

      mv->DestroyAllGeometries();

      gEve->LoadVizDB("geom_gentle_transparentlight.C", kTRUE, kTRUE);

      gEve->LoadVizDB("VizDB_scan_transparentlight.C", kTRUE, kTRUE);

      if(!gEve->GetViewers()->UseLightColorSet())
        gEve->GetViewers()->SwitchColorSet();

      gEve->FullRedraw3D();    

      break;
    }

    case kAEMGreen:
    {
      AliEveMultiView *mv = AliEveMultiView::Instance();

      mv->DestroyAllGeometries();

      gEve->LoadVizDB("geom_gentle_green.C", kTRUE, kTRUE);

      gEve->LoadVizDB("VizDB_scan.C", kTRUE, kTRUE);

      if(!gEve->GetViewers()->UseLightColorSet())
        gEve->GetViewers()->SwitchColorSet();

      gEve->FullRedraw3D();    

      break;
    }

    case kAEMBright:
    {
      AliEveMultiView *mv = AliEveMultiView::Instance();

      mv->DestroyAllGeometries();

      gEve->LoadVizDB("geom_gentle_bright.C", kTRUE, kTRUE);

      gEve->LoadVizDB("VizDB_scan.C", kTRUE, kTRUE);

      if(gEve->GetViewers()->UseLightColorSet())
        gEve->GetViewers()->SwitchColorSet();

      gEve->FullRedraw3D();

      break;
    }

    case kAEMYellow:
    {
      AliEveMultiView *mv = AliEveMultiView::Instance();

      mv->DestroyAllGeometries();

      gEve->LoadVizDB("geom_gentle_yellow.C", kTRUE, kTRUE);

      gEve->LoadVizDB("VizDB_scan_yellow.C", kTRUE, kTRUE);

      if(!gEve->GetViewers()->UseLightColorSet())
        gEve->GetViewers()->SwitchColorSet();

      gEve->FullRedraw3D();    

      break;
    }

    case kAEMTpc:
    {
      AliEveMultiView *mv = AliEveMultiView::Instance();

      mv->DestroyAllGeometries();

      gEve->LoadVizDB("geom_gentle_tpc.C", kTRUE, kTRUE);

      gEve->LoadVizDB("VizDB_scan_tpc.C", kTRUE, kTRUE);

      if(!gEve->GetViewers()->UseLightColorSet())
        gEve->GetViewers()->SwitchColorSet();

      gEve->FullRedraw3D();

      break;
    }

    case kAEMAll: //saving pictures from all three viewers
    {
      
     TGFileInfo fi;
     fi.fFileTypes   = gPictureSaveAsTypes;
     fi.fIniDir      = StrDup(""); // current directory
     fi.fFileTypeIdx = 0;
     fi.fOverwrite   = kTRUE;
     new TGFileDialog(gClient->GetDefaultRoot(),
     gEve->GetMainWindow(), kFDSave, &fi); // dialog 
     if (!fi.fFilename) return;

     TPMERegexp filere(".*/([^/]+$)");
     if (filere.Match(fi.fFilename) != 2)
     {
       Warning("AliEvePopupHandler", "file '%s' bad.", fi.fFilename);
       return;
     }

     TString file1(filere[1]);
     TString file2(filere[1]);
     TString file3(filere[1]);
     
     if (!file1.EndsWith(".png"))
       file1 += "_3D.png"; // adding extensions
       
     if (!file2.EndsWith(".png"))
       file2 += "_RPhi.png"; // adding extensions
       
     if (!file3.EndsWith(".png"))
       file3 += "_RhoZ.png"; // adding extensions

     gSystem->ChangeDirectory(fi.fIniDir);
      
     printf("Saving...\n");

     TEveViewerList *viewers = gEve->GetViewers();
     TEveElement::List_i i = viewers->BeginChildren();
     i++;
     TEveViewer* view3d = ((TEveViewer*)*i);
     view3d->GetGLViewer()->SavePicture(file1); // saving pictures
     i++;
     TEveViewer* viewrphi = ((TEveViewer*)*i);
     viewrphi->GetGLViewer()->SavePicture(file2); // saving pictures
     i++;
     TEveViewer* viewrhoz = ((TEveViewer*)*i);
     viewrhoz->GetGLViewer()->SavePicture(file3); // saving pictures
     
     printf("Done.\n"); 
      
      break;
    }
    
    case kAEM3d: // saving only 3d view
    {
      
     TGFileInfo fi;
     fi.fFileTypes   = gPictureSaveAsTypes;
     fi.fIniDir      = StrDup(""); // current directory
     fi.fFileTypeIdx = 0;
     fi.fOverwrite   = kTRUE;
     new TGFileDialog(gClient->GetDefaultRoot(),
     gEve->GetMainWindow(), kFDSave, &fi);
     if (!fi.fFilename) return;

     TPMERegexp filere(".*/([^/]+$)");
     if (filere.Match(fi.fFilename) != 2)
     {
       Warning("AliEvePopupHandler", "file '%s' bad.", fi.fFilename);
       return;
     }

     TString file1(filere[1]);
     
     if (!file1.EndsWith(".png"))
       file1 += ".png";

     gSystem->ChangeDirectory(fi.fIniDir);
      
     printf("Saving...\n");

     TEveViewerList *viewers = gEve->GetViewers();
     TEveElement::List_i i = viewers->BeginChildren();
     i++;
     TEveViewer* view3d = ((TEveViewer*)*i);
     view3d->GetGLViewer()->SavePicture(file1);
     
     printf("Done.\n"); 
      
      break;
    }
    
     case kAEMRphi: // saving only RPhi view
    {
      
     TGFileInfo fi;
     fi.fFileTypes   = gPictureSaveAsTypes;
     fi.fIniDir      = StrDup(""); // current directory
     fi.fFileTypeIdx = 0;
     fi.fOverwrite   = kTRUE;
     new TGFileDialog(gClient->GetDefaultRoot(),
     gEve->GetMainWindow(), kFDSave, &fi);
     if (!fi.fFilename) return;

     TPMERegexp filere(".*/([^/]+$)");
     if (filere.Match(fi.fFilename) != 2)
     {
       Warning("AliEvePopupHandler", "file '%s' bad.", fi.fFilename);
       return;
     }

     TString file1(filere[1]);
     
     if (!file1.EndsWith(".png"))
       file1 += ".png";
     
     gSystem->ChangeDirectory(fi.fIniDir);
      
     printf("Saving...\n");

     TEveViewerList *viewers = gEve->GetViewers();
     TEveElement::List_i i = viewers->BeginChildren();
     i++;
     i++;
     TEveViewer* viewrphi = ((TEveViewer*)*i);
     viewrphi->GetGLViewer()->SavePicture(file1);
     
     printf("Done.\n"); 
      
      break;
    }
    
     case kAEMRhoz: // saving only RhoZ view
    {
    
     TGFileInfo fi;
     fi.fFileTypes   = gPictureSaveAsTypes;
     fi.fIniDir      = StrDup(""); // current directory
     fi.fFileTypeIdx = 0;
     fi.fOverwrite   = kTRUE;
     new TGFileDialog(gClient->GetDefaultRoot(),
     gEve->GetMainWindow(), kFDSave, &fi);
     if (!fi.fFilename) return;

     TPMERegexp filere(".*/([^/]+$)");
     if (filere.Match(fi.fFilename) != 2)
     {
       Warning("AliEvePopupHandler", "file '%s' bad.", fi.fFilename);
       return;
     }

     TString file1(filere[1]);
       
     if (!file1.EndsWith(".png"))
       file1 += ".png";

     gSystem->ChangeDirectory(fi.fIniDir);
     
     printf("Saving...\n");

     TEveViewerList *viewers = gEve->GetViewers();
     TEveElement::List_i i = viewers->BeginChildren();
     i++;
     i++;
     i++;
     TEveViewer* viewrhoz = ((TEveViewer*)*i);
     viewrhoz->GetGLViewer()->SavePicture(file1);
     
     printf("Done.\n"); 
      
      break;
    }

    case kAEMAllhr: // saving all three views in high resolution
    {
      
     TGFileInfo fi;
     fi.fFileTypes   = gPictureSaveAsTypes;
     fi.fIniDir      = StrDup(""); // current directory
     fi.fFileTypeIdx = 0;
     fi.fOverwrite   = kTRUE;
     new TGFileDialog(gClient->GetDefaultRoot(),
     gEve->GetMainWindow(), kFDSave, &fi);
     if (!fi.fFilename) return;

     TPMERegexp filere(".*/([^/]+$)");
     if (filere.Match(fi.fFilename) != 2)
     {
       Warning("AliEvePopupHandler", "file '%s' bad.", fi.fFilename);
       return;
     }

     TString file1(filere[1]);
     TString file2(filere[1]);
     TString file3(filere[1]);
     
     if (!file1.EndsWith(".png"))
       file1 += "_3D.png";
       
     if (!file2.EndsWith(".png"))
       file2 += "_RPhi.png";
       
     if (!file3.EndsWith(".png"))
       file3 += "_RhoZ.png";

     gSystem->ChangeDirectory(fi.fIniDir);
      
     printf("Saving...\n");

     TEveViewerList *viewers = gEve->GetViewers();
     TEveElement::List_i i = viewers->BeginChildren();
     i++;
     TEveViewer* view3d = ((TEveViewer*)*i);
     view3d->GetGLViewer()->SavePictureScale(file1,6.0); // getting high resolution
     i++;
     TEveViewer* viewrphi = ((TEveViewer*)*i);
     viewrphi->GetGLViewer()->SavePictureScale(file2,6.0);
     i++;
     TEveViewer* viewrhoz = ((TEveViewer*)*i);
     viewrhoz->GetGLViewer()->SavePictureScale(file3,6.0);
     
     printf("Done.\n"); 
      
      break;
    }
    
    case kAEM3dhr: // saving only 3d view in high resolution
    {
      
     TGFileInfo fi;
     fi.fFileTypes   = gPictureSaveAsTypes;
     fi.fIniDir      = StrDup(""); // current directory
     fi.fFileTypeIdx = 0;
     fi.fOverwrite   = kTRUE;
     new TGFileDialog(gClient->GetDefaultRoot(),
     gEve->GetMainWindow(), kFDSave, &fi);
     if (!fi.fFilename) return;

     TPMERegexp filere(".*/([^/]+$)");
     if (filere.Match(fi.fFilename) != 2)
     {
       Warning("AliEvePopupHandler", "file '%s' bad.", fi.fFilename);
       return;
     }

     TString file1(filere[1]);
     
     if (!file1.EndsWith(".png"))
       file1 += ".png";

     gSystem->ChangeDirectory(fi.fIniDir);
      
     printf("Saving...\n");

     TEveViewerList *viewers = gEve->GetViewers();
     TEveElement::List_i i = viewers->BeginChildren();
     i++;
     TEveViewer* view3d = ((TEveViewer*)*i);
     view3d->GetGLViewer()->SavePictureScale(file1,4.0);
     
     printf("Done.\n"); 
      
      break;
    }
    
     case kAEMRphihr: // saving only RPhi view in high resolution
    {
      
     TGFileInfo fi;
     fi.fFileTypes   = gPictureSaveAsTypes;
     fi.fIniDir      = StrDup(""); // current directory
     fi.fFileTypeIdx = 0;
     fi.fOverwrite   = kTRUE;
     new TGFileDialog(gClient->GetDefaultRoot(),
     gEve->GetMainWindow(), kFDSave, &fi);
     if (!fi.fFilename) return;

     TPMERegexp filere(".*/([^/]+$)");
     if (filere.Match(fi.fFilename) != 2)
     {
       Warning("AliEvePopupHandler", "file '%s' bad.", fi.fFilename);
       return;
     }

     TString file1(filere[1]);
     
     if (!file1.EndsWith(".png"))
       file1 += ".png";
     
     gSystem->ChangeDirectory(fi.fIniDir);
      
     printf("Saving...\n");

     TEveViewerList *viewers = gEve->GetViewers();
     TEveElement::List_i i = viewers->BeginChildren();
     i++;
     i++;
     TEveViewer* viewrphi = ((TEveViewer*)*i);
     viewrphi->GetGLViewer()->SavePictureScale(file1,4.0);
     
     printf("Done.\n"); 
      
      break;
    }
    
     case kAEMRhozhr: // saving only RhoZ view in high resolution
    {
    
     TGFileInfo fi;
     fi.fFileTypes   = gPictureSaveAsTypes;
     fi.fIniDir      = StrDup(""); // current directory
     fi.fFileTypeIdx = 0;
     fi.fOverwrite   = kTRUE;
     new TGFileDialog(gClient->GetDefaultRoot(),
     gEve->GetMainWindow(), kFDSave, &fi);
     if (!fi.fFilename) return;

     TPMERegexp filere(".*/([^/]+$)");
     if (filere.Match(fi.fFilename) != 2)
     {
       Warning("AliEvePopupHandler", "file '%s' bad.", fi.fFilename);
       return;
     }

     TString file1(filere[1]);
       
     if (!file1.EndsWith(".png"))
       file1 += ".png";

     gSystem->ChangeDirectory(fi.fIniDir);
     
     printf("Saving...\n");

     TEveViewerList *viewers = gEve->GetViewers();
     TEveElement::List_i i = viewers->BeginChildren();
     i++;
     i++;
     i++;
     TEveViewer* viewrhoz = ((TEveViewer*)*i);
     viewrhoz->GetGLViewer()->SavePictureScale(file1,4.0);
     
     printf("Done.\n"); 
      
      break;
    }
    case kAEMSave://saving VizDB
    {
      TGFileInfo fi;
      fi.fFileTypes   = gMacroSaveAsTypes;
      fi.fIniDir      = StrDup(""); // current directory
      fi.fFileTypeIdx = 0;
      fi.fOverwrite   = kTRUE;
      new TGFileDialog(gClient->GetDefaultRoot(), gEve->GetMainWindow(), kFDSave, &fi);
      if (!fi.fFilename) return;

      TPMERegexp filere(".*/([^/]+$)");
      if (filere.Match(fi.fFilename) != 2)
      {
        Warning("AliEvePopupHandler", "file '%s' bad.", fi.fFilename);
        return;
      }
      printf("Saving...\n");

      TString file(filere[1]);
      if (!file.EndsWith(".C"))
        file += ".C";
      gSystem->ChangeDirectory(fi.fIniDir);
      gEve->SaveVizDB(file);

//Last line "gEve->SaveVizDB(file);" gives macro with many unnecessary
//lines like "x038->SetMinPt(0);" tahat are not interpreted properly later

      string text;
      string all;

      ifstream myfile1(file);
      if(myfile1.is_open())
        {
        while(!myfile1.eof())
          {
            getline(myfile1,text);
            TString check(text);
            if(!(check.EndsWith("MinPt(0);")||check.EndsWith("MaxPt(0);")
               ||check.EndsWith("LimPt(0);")||check.EndsWith("MinP(0);")
               ||check.EndsWith("MaxP(0);")||check.EndsWith("LimP(0);")))
              {
              all += text; //Cut off unnecessary lines and bring everything together
              all += "\n";
              }
          }
        myfile1.close();
        }

      ofstream myfile2(file); //Replacing old file with the one without "bad" lines
      myfile2 << all;
      myfile2.close();

      printf("Done.\n");
      break;

    }

    case kAEMOpen://Opening VizDB
    {
      TGFileInfo fi;
      fi.fFileTypes   = gMacroSaveAsTypes;
      fi.fIniDir      = StrDup(""); // current directory
      fi.fFileTypeIdx = 0;
      fi.fOverwrite   = kTRUE;
      new TGFileDialog(gClient->GetDefaultRoot(), gEve->GetMainWindow(), kFDOpen, &fi);
      if (!fi.fFilename) return;

      TPMERegexp filere(".*/([^/]+$)");
      if (filere.Match(fi.fFilename) != 2)
      {
        Warning("AliEvePopupHandler", "file '%s' bad.", fi.fFilename);
        return;
      }
      printf("Opening...\n");

      TString file(filere[1]);

      gSystem->ChangeDirectory(fi.fIniDir);

      gEve->LoadVizDB(file, kTRUE, kTRUE);

      gEve->Redraw3D(kTRUE);

      printf("Done.\n");
      break;

    }

    case kAEMSetDefault://Restore default settings
    {

      printf("Setting...\n");

      TEveBrowser *browser = gEve->GetBrowser();
      browser->ShowCloseTab(kFALSE);

      if(fLoadCheck)
        browser->RemoveTab(TRootBrowser::kRight, 5);//remove the tab with previous DataSelection window
      else
        browser->RemoveTab(TRootBrowser::kRight, 2);

      TEveUtil::Macro("geom_gentle_default.C");
      gEve->LoadVizDB("VizDB_scan.C", kTRUE, kTRUE);
      TEveUtil::Macro("DataSelection_init.C");

      if(!gEve->GetViewers()->UseLightColorSet())
        gEve->GetViewers()->SwitchColorSet(); //white background

      AliEveEventManager *eman = AliEveEventManager::GetMaster();//reload event (gEve->Refresh() crashes)
      Int_t ev = eman->GetEventId();
      eman->Close();
      eman->Open();
      eman->GotoEvent(ev);

      printf("Done.\n");

      fLoadCheck = kTRUE;

      gEve->Redraw3D(kTRUE);


    }

/*
    case kAEMResiduals:
    {

      TEveUtil::Macro("make_residuals.C");

      break;

    }
*/

    case kAEMCuts:
    {

        TEveBrowser *browser = gEve->GetBrowser();
        browser->StartEmbedding(TRootBrowser::kLeft);
        new AliEveCutsWindow();
        browser->StopEmbedding("Cuts");
        
      break;

    }

    case kAEMVectors:
      {
        new AliEveMomentumVectors();

      break;

    }

/*
    case kAEMGui:
    {

      TEveUtil::Macro("alieve_gui_mode.C");

      break;

    }
*/
          /*
      case kStorageListEvents:
      {
#ifdef ZMQ
          fListEventsWindow =
          AliStorageAdministratorPanelListEvents::GetInstance();
          
          fListEventsWindow->Connect("SelectedEvent()","AliEveConfigManager",this,"SetEventInEventManager()");
#endif
          break;
          
      }
           */
      case kStorageMarkEvent:
      {
#ifdef ZMQ
          AliStorageAdministratorPanelMarkEvent *markEventWindow =
          AliStorageAdministratorPanelMarkEvent::GetInstance();
#endif
          break;
      }
      case kPreferences:
      {
          AliEvePreferencesWindow::Instance();
          
          break;
      }
          
      default:
      {
          Warning(kEH, "Unknown menu entry.");
          break;
      }
  }
}

void AliEveConfigManager::SetEventInEventManager()
{
#ifdef ZMQ

    AliEveEventManager *manager = AliEveEventManager::GetMaster();
    AliStorageAdministratorPanelListEvents* fListEventsTab = AliStorageAdministratorPanelListEvents::GetInstance();
    AliESDEvent *event = fListEventsTab->GetSelectedEvent();
    
    if(event)
    {
	    cout<<"SETTING EVENT IN ED"<<endl;

        manager->SetAutoLoad(kFALSE);
        AliEveDataSource *dataSource = manager->GetCurrentDataSource();
        dataSource->SetEventFromStorageManager(event);
        dataSource->NextEvent();
    }
#endif
}

void AliEveConfigManager::StorageManagerChangedState(int state)
{
#ifdef ZMQ
    AliEveEventManager *manager = AliEveEventManager::GetMaster();
    AliStorageAdministratorPanelListEvents* listEventsTab = AliStorageAdministratorPanelListEvents::GetInstance();
    
//    if (manager->IsOnlineMode()) {
        if (state == 0)// storage manager is down
        {
            listEventsTab->SetOfflineMode(kTRUE);
        }
        else if(state == 1)// storage manager is alive
        {
            listEventsTab->SetOfflineMode(kFALSE);
        }
//    }
#endif
}


