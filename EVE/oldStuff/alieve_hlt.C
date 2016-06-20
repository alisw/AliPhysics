/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveEventManager.h"
#include "AliEveEventManagerEditor.h"
#include "AliEveMultiView.h"
#include "AliEveMacroExecutor.h"
#include "AliEveMacro.h"
#include "AliSysInfo.h"
#include "AliCDBManager.h"

#include <TH2D.h>
#include <TTimeStamp.h>
#include <TROOT.h>
#include <TEveManager.h>
#include <TEveBrowser.h>
#include <TEveMacro.h>
#include <TGTab.h>
#include <TGFileBrowser.h>
#include <TInterpreter.h>
#include <TEnv.h>
#include <TPRegexp.h>
#include <TFolder.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TList.h>
#include "TString.h"

#include <iostream>

using namespace std;

class TEveProjectionManager;
class TEveGeoShape;
class TEveUtil;
class AliTriggerAnalysis;
class AliSysInfo;

TH2D *V0StateHistogram;
Bool_t gCenterProjectionsAtPrimaryVertex = kFALSE;

Int_t      g_pic_id  = 0;
Int_t      g_pic_max = 100;
TTimeStamp g_pic_prev(0, 0);

void alieve_init_import_macros();
void alieve_online_on_new_event();

void alieve_hlt(TString ocdbStorage="local://OCDB")
{
    if (gSystem->Getenv("ocdbStorage"))
      ocdbStorage=gSystem->Getenv("ocdbStorage");
    AliEveEventManager::SetCdbUri(ocdbStorage);         // current OCDB snapshot
    
    cout<<"Creating multiview...";
    AliEveMultiView *multiView = new AliEveMultiView(kTRUE);
    cout<<"created"<<endl;
    
    Info("alieve_init", "Adding standard macros.");
    alieve_init_import_macros();
    gROOT->ProcessLine(".L geom_gentle.C+");
    gROOT->ProcessLine(".L geom_gentle_trd.C+");
    gROOT->ProcessLine(".L geom_gentle_muon.C+");

    AliEveEventManager* eventManager = new AliEveEventManager(AliEveEventManager::kSourceHLT);
    gEve->AddEvent(eventManager);
    cout<<"Event manager created"<<endl;
    
    TEveUtil::AssertMacro("VizDB_scan.C");
    gSystem->ProcessEvents();
    cout<<"VizDB_scan loaded"<<endl;
    TEveBrowser *browser = gEve->GetBrowser();
    browser->ShowCloseTab(kFALSE);
    cout<<"browser created"<<endl;
    
    TEveUtil::LoadMacro("geom_gentle.C");
    cout<<"geom gentle loaded"<<endl;
    
    multiView->InitGeomGentle(geom_gentle(),
                              geom_gentle_rphi(),
                              geom_gentle_rhoz(),
                              geom_gentle_rhoz());
    cout<<"geom gentl inited"<<endl;
    
    //TEveUtil::LoadMacro("geom_gentle_trd.C");
    multiView->InitGeomGentleTrd(geom_gentle_trd());

    //TEveUtil::LoadMacro("geom_gentle_muon.C");
    multiView->InitGeomGentleMuon(geom_gentle_muon(), kFALSE, kFALSE, kTRUE);
 
    printf("============ Setting macro executor ============\n");
    
    AliEveMacroExecutor *exec = AliEveEventManager::Instance()->GetExecutor();
    printf("exec created\n");
    /*
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX",         "primary_vertex.C", "primary_vertex",             "",                kTRUE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX Ellipse", "primary_vertex.C", "primary_vertex_ellipse",     "",                kTRUE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX Box",     "primary_vertex.C", "primary_vertex_box",         "kFALSE, 3, 3, 3", kFALSE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX",         "primary_vertex.C", "primary_vertex_spd",         "",                kTRUE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX Ellipse", "primary_vertex.C", "primary_vertex_ellipse_spd", "",                kTRUE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX Box",     "primary_vertex.C", "primary_vertex_box_spd",     "kFALSE, 3, 3, 3", kFALSE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX",         "primary_vertex.C", "primary_vertex_tpc",         "",                kFALSE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX Ellipse", "primary_vertex.C", "primary_vertex_ellipse_tpc", "",                kFALSE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX Box",     "primary_vertex.C", "primary_vertex_box_tpc",     "kFALSE, 3, 3, 3", kFALSE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Track",      "esd_tracks.C",        "esd_tracks",             "", kFALSE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Track",      "esd_tracks.C",        "esd_tracks_MI",          "", kFALSE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Track MUON", "esd_muon_tracks.C", "esd_muon_tracks",        "kTRUE,kFALSE", kTRUE));
    */
    
    
    // default appearance:
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Tracks by category",  "esd_tracks.C", "esd_tracks_by_category",  "", kTRUE));
    
    // preset for cosmics:
    //exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Tracks by category",  "esd_tracks.C", "esd_tracks_by_category",  "kGreen,kGreen,kGreen,kGreen,kGreen,kGreen,kGreen,kGreen,kGreen,kFALSE", kTRUE));
    
    //exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC FMD",        "fmd_esd.C",           "fmd_esd",                "", kTRUE));//huge leak
    //
    
    // ???
    // exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC TRD", "trd_detectors.C", "trd_detectors",         "", kFALSE));
    // trd_tracks disabled due to memory leaks
    cout<<"macros added to exec"<<endl;

    //----------------------------------------------------------------------------
    /* something is wrong here:
    slot = TEveWindow::CreateWindowInTab(browser->GetTabRight());
    slot->StartEmbedding();
    AliEveMacroExecutorWindow* exewin = new AliEveMacroExecutorWindow(exec);
    slot->StopEmbedding("DataSelection");
    exewin->PopulateMacros();
    */
    //============================================================================
    // Final GUI setup
    //============================================================================
    
    browser->GetTabRight()->SetTab(1);
    browser->StartEmbedding(TRootBrowser::kBottom);
    new AliEveEventManagerWindow(AliEveEventManager::Instance());
    browser->StopEmbedding("EventCtrl");
    
    gEve->FullRedraw3D(kTRUE);
    gSystem->ProcessEvents();
    
    // move and rotate sub-views
    TGLViewer *glv1 = multiView->Get3DView()->GetGLViewer();
    TGLViewer *glv2 = multiView->GetRPhiView()->GetGLViewer();
    TGLViewer *glv3 = multiView->GetRhoZView()->GetGLViewer();
    TGLViewer *glv4 = multiView->GetMuonView()->GetGLViewer();
    
    glv1->CurrentCamera().RotateRad(-0.4, 0.6);
    glv2->CurrentCamera().Dolly(90, kFALSE, kFALSE);
    glv3->CurrentCamera().Dolly(2300, kFALSE, kFALSE);
    glv4->CurrentCamera().Dolly(1, kFALSE, kFALSE);
    
    AliEveEventManager::Instance()->AddNewEventCommand("alieve_online_on_new_event();");
    
    gEve->FullRedraw3D();
    gSystem->ProcessEvents();
    gEve->Redraw3D(kTRUE);
    
    // set autoload by default
    AliEveEventManager::Instance()->SetAutoLoad(true);
}

void alieve_online_on_new_event()
{
  if (!AliEveEventManager::HasESD()) return;

  AliESDEvent* esd = AliEveEventManager::AssertESD();
	
  Double_t x[3] = { 0, 0, 0 };
  esd->GetPrimaryVertex()->GetXYZ(x);
  
  TTimeStamp ts(esd->GetTimeStamp());
  TString win_title("Eve Main Window -- Timestamp: ");
  win_title += ts.AsString("s");
  win_title += "; Event # in ESD file: ";
  win_title += esd->GetEventNumberInFile();
  win_title += " | run: ";
  win_title += esd->GetRunNumber();
  gEve->GetBrowser()->SetWindowName(win_title);
    
	TEveElement* top = gEve->GetCurrentEvent();
    
	AliEveMultiView *mv = AliEveMultiView::Instance();
    
	mv->DestroyEventRPhi();
	if (gCenterProjectionsAtPrimaryVertex){
	  mv->SetCenterRPhi(x[0], x[1], x[2]);
	}
	mv->ImportEventRPhi(top);
    
	mv->DestroyEventRhoZ();
	if (gCenterProjectionsAtPrimaryVertex){
	  mv->SetCenterRhoZ(x[0], x[1], x[2]);
	}
	mv->ImportEventRhoZ(top);
    
	if (gCenterProjectionsAtPrimaryVertex)
	  mv->SetCenterMuon(x[0], x[1], x[2]);
	mv->ImportEventMuon(top);
}

void alieve_init_import_macros()
{
  // Put macros in the list of browsables, add a macro browser to
  // top-level GUI.

  TString macdir("$(ALICE_ROOT)/EVE/macros");

  if (gSystem->Getenv("ALICE_ROOT") != 0)
  {
    gInterpreter->AddIncludePath(Form("%s/MUON", gSystem->Getenv("ALICE_ROOT")));
    gInterpreter->AddIncludePath(Form("%s/MUON/mapping", gSystem->Getenv("ALICE_ROOT")));
    gSystem->ExpandPathName(macdir);
  }

  TFolder* f = gEve->GetMacroFolder();
  void* dirhandle = gSystem->OpenDirectory(macdir.Data());
  if (dirhandle != 0)
  {
    char* filename;
    TPMERegexp re("\\.C$");
    TObjArray names;
    while ((filename = (char*)(gSystem->GetDirEntry(dirhandle))) != 0)
    {
      if (re.Match(filename))
        names.AddLast(new TObjString(filename));
    }
    names.Sort();

    for (Int_t ii=0; ii<names.GetEntries(); ++ii)
    {
      TObjString * si = (TObjString*) names.At(ii);
      f->Add(new TEveMacro(Form("%s/%s", macdir.Data(), (si->GetString()).Data())));
    }
  }
  gSystem->FreeDirectory(dirhandle);

  gROOT->GetListOfBrowsables()->Add(new TSystemDirectory(macdir.Data(), macdir.Data()));

  {
    TEveBrowser   *br = gEve->GetBrowser();
    TGFileBrowser *fb = 0;
    fb = br->GetFileBrowser();
    fb->GotoDir(macdir);
    {
      br->StartEmbedding(0);
      fb = br->MakeFileBrowser();
      fb->BrowseObj(f);
      fb->Show();
      br->StopEmbedding();
      br->SetTabTitle("Macros", 0);
      br->SetTab(0, 0);
    }
  }
}

