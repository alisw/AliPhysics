// $Id$

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

/// \ingroup evemacros
/// \file muon_init.C
///
/// \author P. Pillot, L. Aphecetche; Subatech

class AliEveMacroExecutor;
class TEveProjectionManager;
class TEveGeoShape;
class TEveUtil;
class TSystem;
class TInterpreter;

Bool_t gShowMuonRPhi = kFALSE;
Bool_t gShowMuonRhoZ = kTRUE;
Bool_t gShowMuon = kTRUE;

Bool_t gCenterProjectionsAtPrimaryVertex = kFALSE;


void muon_init(const TString& cdburi = "",
               const TString& path   = ".",
	       Bool_t showBarrel = kFALSE)
{
  if (gSystem->Getenv("ALICE_ROOT") != 0)
  {
    gInterpreter->AddIncludePath(Form("%s/MUON", gSystem->Getenv("ALICE_ROOT")));
    gInterpreter->AddIncludePath(Form("%s/MUON/mapping", gSystem->Getenv("ALICE_ROOT")));
  }
  
  if (cdburi.IsNull() && ! AliCDBManager::Instance()->IsDefaultStorageSet())
  {
    gEnv->SetValue("Root.Stacktrace", "no");
    Fatal("muon_init.C", "OCDB path MUST be specified as the first argument.");
  }
  
  TEveUtil::LoadMacro("alieve_init.C");
  path.Remove(TString::kTrailing, '/');
  if (path.BeginsWith("alien:")) AliEveEventManager::SearchRawForCentralReconstruction();
  alieve_init(cdburi, path, -1);
  
  TEveUtil::AssertMacro("VizDB_scan.C");
  
  AliEveMacroExecutor *exec    = AliEveEventManager::Instance()->GetExecutor();
  TEveBrowser         *browser = gEve->GetBrowser();
  browser->ShowCloseTab(kFALSE);
  
  
  //==============================================================================
  // Geometry, scenes, projections and viewers
  //==============================================================================
  
  AliEveMultiView *mv = new AliEveMultiView;
  
  mv->SetDepth(-10);
  
  TEveUtil::LoadMacro("geom_gentle.C");
  TEveUtil::LoadMacro("geom_gentle_muon.C+");

  mv->InitGeomGentle(geom_gentle(), geom_gentle_rphi(), geom_gentle_rhoz(), geom_gentle_muon(kFALSE));
  
  mv->InitGeomGentleMuon(geom_gentle_muon(kFALSE), gShowMuonRPhi, gShowMuonRhoZ, gShowMuon);
  
  mv->SetDepth(0);
  
  //==============================================================================
  // Registration of per-event macros
  //==============================================================================
  
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "SIM Track","kine_tracks.C+",   "kine_tracks",  "", kFALSE));

  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "SIM TrackRef","muon_trackRefs.C+","muon_trackRefs","kTRUE", kTRUE));
  
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW MUON", "muon_raw.C+",     "muon_raw",     "", kTRUE));

  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "DIG MUON", "muon_digits.C+",  "muon_digits",  "", kFALSE));
  
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "CLU MUON", "muon_clusters.C+","muon_clusters","", kTRUE));

  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Track", "esd_muon_tracks.C+", "esd_muon_tracks","kTRUE,kTRUE", kTRUE));

  if (showBarrel) {
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX", "primary_vertex.C+", "primary_vertex", "", kTRUE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX SPD", "primary_vertex.C+", "primary_vertex_spd", "", kTRUE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Tracks by category", "esd_tracks.C+", "esd_tracks_by_category", "", kTRUE));
  }
  
  //==============================================================================
  // Additional GUI components
  //==============================================================================
  
  // Macro / data selection
  slot = TEveWindow::CreateWindowInTab(browser->GetTabRight());
  slot->StartEmbedding();
  AliEveMacroExecutorWindow* exewin = new AliEveMacroExecutorWindow(exec);
  slot->StopEmbedding("DataSelection");
  exewin->PopulateMacros();
  
  // Event selection tab
  slot = TEveWindow::CreateWindowInTab(browser->GetTabRight());
  slot->StartEmbedding();
  new AliEveEventSelectorWindow(gClient->GetRoot(), 600, 400, AliEveEventManager::Instance()->GetEventSelector());
  slot->StopEmbedding("Selections");
  
  // QA viewer
  slot = TEveWindow::CreateWindowInTab(browser->GetTabRight());
  slot->StartEmbedding();
  new AliQAHistViewer(gClient->GetRoot(), 600, 400, kTRUE);
  slot->StopEmbedding("QA histograms");
  
  browser->GetTabRight()->SetTab(1);
  
  browser->StartEmbedding(TRootBrowser::kBottom);
  new AliEveEventManagerWindow(AliEveEventManager::Instance());
  browser->StopEmbedding("EventCtrl");
  
  slot = TEveWindow::CreateWindowInTab(browser->GetTabRight());
  TEveWindowTab *store_tab = slot->MakeTab();
  store_tab->SetElementNameTitle("WindowStore",
				 "Undocked windows whose previous container is not known\n"
				 "are placed here when the main-frame is closed.");
  gEve->GetWindowManager()->SetDefaultContainer(store_tab);
  
  
  //==============================================================================
  // AliEve objects - global tools
  //==============================================================================
  
  AliEveTrackCounter* g_trkcnt = new AliEveTrackCounter("Primary Counter");
  gEve->AddToListTree(g_trkcnt, kFALSE);
  
  
  //==============================================================================
  // Final stuff
  //==============================================================================
  
  // A refresh to show proper window.
  //gEve->GetViewers()->SwitchColorSet();
  gEve->Redraw3D(kTRUE);
  gSystem->ProcessEvents();
  
  // Register command to call on each event.
  AliEveEventManager::Instance()->AddNewEventCommand("on_new_event();");
  AliEveEventManager::Instance()->GotoEvent(0);
  
  gEve->EditElement(g_trkcnt);
  
  gEve->Redraw3D(kTRUE);
  
  // Assure 3D view rotates around the origin.
  gSystem->ProcessEvents();
  AliEveMultiView::Instance()->Get3DView()->GetGLViewer()->CurrentCamera().SetCenterVec(0,0,0);
  AliEveMultiView::Instance()->Get3DView()->GetGLViewer()->RequestDraw();
}

/******************************************************************************/

void on_new_event()
{
  Double_t x[3] = { 0, 0, 0 };
  
  if (AliEveEventManager::HasESD())
  {
    AliESDEvent* esd = AliEveEventManager::Instance()->AssertESD();
    esd->GetPrimaryVertex()->GetXYZ(x);
    
    TTimeStamp ts(esd->GetTimeStamp());
    TString win_title("Eve Main Window -- Timestamp: ");
    win_title += ts.AsString("s");
    win_title += "; Event # in ESD file: ";
    win_title += esd->GetEventNumberInFile();
    gEve->GetBrowser()->SetWindowName(win_title);
  }
  
  TEveElement* top = gEve->GetCurrentEvent();
  
  AliEveMultiView *mv = AliEveMultiView::Instance();
  
  mv->DestroyEventRPhi();
  if (gCenterProjectionsAtPrimaryVertex)
    mv->SetCenterRPhi(x[0], x[1], x[2]);
  mv->ImportEventRPhi(top);
  
  mv->DestroyEventRhoZ();
  if (gCenterProjectionsAtPrimaryVertex)
    mv->SetCenterRhoZ(x[0], x[1], x[2]);
  mv->ImportEventRhoZ(top);
}
