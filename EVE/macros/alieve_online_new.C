/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include <../../MONITOR/AliSocket.h>

class TEveProjectionManager;
class TEveGeoShape;
class TEveUtil;
class AliTriggerAnalysis;
class AliSysInfo;

TH2D* V0StateHistogram;

Bool_t gCenterProjectionsAtPrimaryVertex = kFALSE;

void alieve_online_new()
{
    printf("alieve_online_init() ...\n");
 printf("================================ Correct macro r ...\n");
 
	if (gSystem->Getenv("ALICE_ROOT") != 0)
  {
    gInterpreter->AddIncludePath(Form("%s/MUON", gSystem->Getenv("ALICE_ROOT")));
    gInterpreter->AddIncludePath(Form("%s/MUON/mapping", gSystem->Getenv("ALICE_ROOT")));
  }
  
   AliEveEventManager::SetCdbUri("local://$ALICE_ROOT/OCDB");
  
  Info("alieve_init", "Adding standard macros.");
  TString  hack = gSystem->pwd(); // Problem with TGFileBrowser cding
  alieve_init_import_macros();
  gSystem->cd(hack);
  
  new AliEveEventManager("online", -1);
  //AliEveEventManager::GetMaster()->AddNewEventCommand("alieve_online_on_new_event()");
  gEve->AddEvent(AliEveEventManager::GetMaster());
  
  
  TEveUtil::AssertMacro("VizDB_scan.C");
  
  gSystem->ProcessEvents();
   
  
  
 
  
  AliEveMacroExecutor *exec  = AliEveEventManager::GetMaster()->GetExecutor();
  TEveBrowser         *browser = gEve->GetBrowser();
  browser->ShowCloseTab(kFALSE);
  

  // Gentle-geom loading changes gGeoManager.
  //TEveGeoManagerHolder mgrRestore;

  AliEveMultiView *multiView = new AliEveMultiView(kTRUE);

  TEveUtil::LoadMacro("geom_gentle.C");
  multiView->InitGeomGentle(geom_gentle(),
                             geom_gentle_rphi(), 
                             geom_gentle_rhoz(),
                             geom_gentle_rhoz());

  TEveUtil::LoadMacro("geom_gentle_trd.C");
  multiView->InitGeomGentleTrd(geom_gentle_trd());

  TEveUtil::LoadMacro("geom_gentle_muon.C");
  multiView->InitGeomGentleMuon(geom_gentle_muon(), kFALSE, kFALSE, kTRUE);

  //============================================================================
  // Standard macros to execute -- not all are enabled by default.
  //============================================================================

  
   printf("============ Setting macro executor\n");

  AliEveMacroExecutor *exec = AliEveEventManager::GetMaster()->GetExecutor();

  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX",         "primary_vertex.C", "primary_vertex",             "",                kTRUE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX Ellipse", "primary_vertex.C", "primary_vertex_ellipse",     "",                kTRUE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX Box",     "primary_vertex.C", "primary_vertex_box",         "kFALSE, 3, 3, 3", kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX",         "primary_vertex.C", "primary_vertex_spd",         "",                kTRUE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX Ellipse", "primary_vertex.C", "primary_vertex_ellipse_spd", "",                kTRUE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX Box",     "primary_vertex.C", "primary_vertex_box_spd",     "kFALSE, 3, 3, 3", kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX",         "primary_vertex.C", "primary_vertex_tpc",         "",                kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX Ellipse", "primary_vertex.C", "primary_vertex_ellipse_tpc", "",                kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX Box",     "primary_vertex.C", "primary_vertex_box_tpc",     "kFALSE, 3, 3, 3", kFALSE));

  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clus ITS",   "its_clusters.C",   "its_clusters"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clus TPC",   "tpc_clusters.C",   "tpc_clusters"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clus TRD",   "trd_clusters.C",   "trd_clusters"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clus TOF",   "tof_clusters.C",   "tof_clusters"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clus HMPID", "hmpid_clusters.C", "hmpid_clusters"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clus MUON",  "muon_clusters.C",  "muon_clusters"));

  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "DIG EMCAL",   "emcal_digits.C",   "emcal_digits"));

  exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW ITS",     "its_raw.C",     "its_raw"));
  //  exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW TPC",     "tpc_raw.C",     "tpc_raw"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW TOF",     "tof_raw.C",     "tof_raw"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW VZERO",   "vzero_raw.C",   "vzero_raw", "", kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW ACORDE",  "acorde_raw.C",  "acorde_raw", "", kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW MUON",    "muon_raw.C",  "muon_raw"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW FMD",     "fmd_raw.C",     "fmd_raw"));

  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Track",      "esd_tracks.C",        "esd_tracks",             "", kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Track",      "esd_tracks.C",        "esd_tracks_MI",          "", kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Track",      "esd_tracks.C",        "esd_tracks_by_category", "", kTRUE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Track MUON", "esd_muon_tracks.C", "esd_muon_tracks",        "kTRUE,kFALSE", kTRUE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC FMD",        "fmd_esd.C",           "fmd_esd",                "", kTRUE));

  // ???
  // exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC TRD", "trd_detectors.C", "trd_detectors",         "", kFALSE));
  // trd_tracks disabled due to memory leaks

  //----------------------------------------------------------------------------

  slot = TEveWindow::CreateWindowInTab(browser->GetTabRight());
  slot->StartEmbedding();
  AliEveMacroExecutorWindow* exewin = new AliEveMacroExecutorWindow(exec);
  slot->StopEmbedding("DataSelection");
  exewin->PopulateMacros();

  //============================================================================
  // Final GUI setup
  //============================================================================

  browser->GetTabRight()->SetTab(1);

  browser->StartEmbedding(TRootBrowser::kBottom);
  new AliEveEventManagerWindow(AliEveEventManager::GetMaster());
  browser->StopEmbedding("EventCtrl");

  browser->MoveResize(0, 0, gClient->GetDisplayWidth(),
		      gClient->GetDisplayHeight() - 32);

  gEve->GetViewers()->SwitchColorSet();
  gEve->FullRedraw3D(kTRUE);
	gSystem->ProcessEvents();

  TGLViewer *glv1 = multiView->Get3DView()->GetGLViewer();
  TGLViewer *glv2 = multiView->GetRPhiView()->GetGLViewer();
  TGLViewer *glv3 = multiView->GetRhoZView()->GetGLViewer();

  glv1->CurrentCamera().RotateRad(-0.4, -1.8);
  glv2->CurrentCamera().Dolly(450, kFALSE, kFALSE);
  glv3->CurrentCamera().Dolly(1500, kFALSE, kFALSE);

  gEve->FullRedraw3D();
  gSystem->ProcessEvents();
  
  // Register command to call on each event.
  // AliEveEventManager::GetMaster()->AddNewEventCommand("alieve_online_on_new_event();");
  //AliEveEventManager::GetMaster()->GotoEvent(-1);
 
  printf("================================ Connecting to Server ...\n");
   
    //AliEveEventManager::ConnectToServer("tcp://137.138.55.173", 5024);
  if(AliEveEventManager::ConnectToServer("tcp://137.138.93.150", 5024))
  {
	  printf("\nconnected\n");
  }	else printf("not connected\n");
    
    AliSocket* subscriber = AliEveEventManager::AssertSubscriber();
  
  	if(subscriber ==0) {
  		printf("===================== Not connected! ====================\n");
    }

  

  gEve->Redraw3D(kTRUE);
    
}


Int_t      g_pic_id  = 0;
Int_t      g_pic_max = 100;
TTimeStamp g_pic_prev(0, 0);

void alieve_online_on_new_event()
{
  AliSysInfo::AddStamp("on_new_event_start");
 
	AliSysInfo::AddStamp("on_new_event_end");
}

void alieve_init_import_macros()
{
  // Put macros in the list of browsables, add a macro browser to
  // top-level GUI.

  TString macdir("$(ALICE_ROOT)/EVE/alice-macros");
  gSystem->ExpandPathName(macdir);

  TFolder* f = gEve->GetMacroFolder();
  void* dirhandle = gSystem->OpenDirectory(macdir.Data());
  if (dirhandle != 0)
  {
    char* filename;
    TPMERegexp re("\\.C$");
    TObjArray names;
    while ((filename = gSystem->GetDirEntry(dirhandle)) != 0)
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
