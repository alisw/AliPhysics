// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

class AliEveMacroExecutor;
class MultiView;
class TEveProjectionManager;
class TEveGeoShape;
class TEveUtil;


Bool_t gShowTrd      = kTRUE;
Bool_t gShowMuonRPhi = kFALSE;
Bool_t gShowMuonRhoZ = kTRUE;

Bool_t gCenterProjectionsAtPrimaryVertex = kFALSE;


void visscan_init(const TString& cdburi = "",
		  const TString& path   = ".", Bool_t show_extra_geo = kFALSE)
{
  if (cdburi.IsNull() && ! AliCDBManager::Instance()->IsDefaultStorageSet())
  {
    gEnv->SetValue("Root.Stacktrace", "no");
    Fatal("visscan_init.C", "OCDB path MUST be specified as the first argument.");
  }

  if (gROOT->LoadMacro("MultiView.C+") != 0)
  {
    gEnv->SetValue("Root.Stacktrace", "no");
    Fatal("visscan_init.C", "Failed loading MultiView.C in compiled mode.");
  }

  if (!show_extra_geo)
  {
    gShowTrd = gShowMuonRPhi = gShowMuonRhoZ = kFALSE;
  }

  AliEveEventManager::AddAODfriend("AliAOD.VertexingHF.root");

  TEveUtil::LoadMacro("alieve_init.C");
  alieve_init(cdburi, path, -1);

  // TEveLine::SetDefaultSmooth(1);

  TEveUtil::AssertMacro("VizDB_scan.C");

  AliEveMacroExecutor *exec    = AliEveEventManager::GetMaster()->GetExecutor();
  TEveBrowser         *browser = gEve->GetBrowser();
  browser->ShowCloseTab(kFALSE);


  //==============================================================================
  // Geometry, scenes, projections and viewers
  //==============================================================================

  gMultiView = new MultiView;

  TEveUtil::LoadMacro("geom_gentle.C");
  gMultiView->InitGeomGentle(geom_gentle(),
                             geom_gentle_rphi(), 
                             geom_gentle_rhoz());

  if (gShowTrd) {
    TEveUtil::LoadMacro("geom_gentle_trd.C");
    gMultiView->InitGeomGentleTrd(geom_gentle_trd());
  }

  if (gShowMuonRPhi || gShowMuonRhoZ) {
    TEveUtil::LoadMacro("geom_gentle_muon.C");
    gMultiView->InitGeomGentleMuon(geom_gentle_muon(kFALSE), gShowMuonRPhi, gShowMuonRhoZ);
  }


  //==============================================================================
  // Registration of per-event macros
  //==============================================================================

  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "SIM Track",   "kine_tracks.C", "kine_tracks", "", kFALSE));

  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "SIM Hit ITS", "its_hits.C",    "its_hits",    "", kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "SIM Hit TPC", "tpc_hits.C",    "tpc_hits",    "", kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "SIM Hit T0",  "t0_hits.C",     "t0_hits",     "", kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "SIM Hit FMD", "fmd_hits.C",    "fmd_hits",    "", kFALSE));

  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "DIG ITS",     "its_digits.C",  "its_digits",  "", kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "DIG FMD",     "fmd_digits.C",  "fmd_digits",  "", kFALSE));

  exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW ITS",     "its_raw.C",     "its_raw",     "", kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW TPC",     "tpc_raw.C",     "tpc_raw",     "", kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW T0",      "t0_raw.C",      "t0_raw",      "", kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW FMD",     "fmd_raw.C",     "fmd_raw",     "", kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW VZERO",   "vzero_raw.C",   "vzero_raw",   "", kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW ACORDE",  "acorde_raw.C",  "acorde_raw",  "", kFALSE));

  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX",         "primary_vertex.C", "primary_vertex",             "",                kTRUE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX Ellipse", "primary_vertex.C", "primary_vertex_ellipse",     "",                kTRUE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX Box",     "primary_vertex.C", "primary_vertex_box",         "kFALSE, 3, 3, 3", kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX",         "primary_vertex.C", "primary_vertex_spd",         "",                kTRUE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX Ellipse", "primary_vertex.C", "primary_vertex_ellipse_spd", "",                kTRUE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX Box",     "primary_vertex.C", "primary_vertex_box_spd",     "kFALSE, 3, 3, 3", kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX",         "primary_vertex.C", "primary_vertex_tpc",         "",                kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX Ellipse", "primary_vertex.C", "primary_vertex_ellipse_tpc", "",                kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX Box",     "primary_vertex.C", "primary_vertex_box_tpc",     "kFALSE, 3, 3, 3", kFALSE));

  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC V0",   "esd_V0_points.C",       "esd_V0_points_onfly"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC V0",   "esd_V0_points.C",       "esd_V0_points_offline"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC V0",   "esd_V0.C",              "esd_V0"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC CSCD", "esd_cascade_points.C",  "esd_cascade_points"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC CSCD", "esd_cascade.C",         "esd_cascade"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC KINK", "esd_kink_points.C",     "esd_kink_points"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC KINK", "esd_kink.C",            "esd_kink"));

  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Track", "esd_tracks.C", "esd_tracks",             "", kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Track", "esd_tracks.C", "esd_tracks_MI",          "", kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Track", "esd_tracks.C", "esd_tracks_by_category", "", kTRUE));

  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Tracklet", "esd_spd_tracklets.C", "esd_spd_tracklets", "", kFALSE));

  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC ZDC",      "esd_zdc.C", "esd_zdc", "", kFALSE));

  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clus",     "clusters.C+",     "clusters", "", kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clus ITS", "its_clusters.C+", "its_clusters"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clus TPC", "tpc_clusters.C+", "tpc_clusters"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clus TRD", "trd_clusters.C+", "trd_clusters"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clus TOF", "tof_clusters.C+", "tof_clusters"));

  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clus TPC", "vplot_tpc.C+",    "vplot_tpc", "", kFALSE));

  exec->AddMacro(new AliEveMacro(AliEveMacro::kAOD, "ANA HF",   "aod_HF.C",   "aod_HF",   "", kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kAOD, "ANA Jets", "jetplane.C", "jetplane", "", kFALSE));


  //==============================================================================
  // Additional GUI components
  //==============================================================================

  slot = TEveWindow::CreateWindowInTab(browser->GetTabRight());
  slot->StartEmbedding();
  AliEveMacroExecutorWindow* exewin = new AliEveMacroExecutorWindow(exec);
  slot->StopEmbedding("DataSelection");
  exewin->PopulateMacros();

  slot = TEveWindow::CreateWindowInTab(browser->GetTabRight());
  slot->StartEmbedding();
  new AliQAHistViewer(gClient->GetRoot(), 600, 400, kTRUE);
  slot->StopEmbedding("QA histograms");

  //event selection tab
  slot = TEveWindow::CreateWindowInTab(browser->GetTabRight());
  slot->StartEmbedding();
  new AliEveEventSelectorWindow(gClient->GetRoot(), 600, 400, AliEveEventManager::GetMaster()->GetEventSelector());
  slot->StopEmbedding("Selections");

  browser->GetTabRight()->SetTab(1);

  browser->StartEmbedding(TRootBrowser::kBottom);
  new AliEveEventManagerWindow(AliEveEventManager::GetMaster());
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

  AliEveTrackFitter* fitter = new AliEveTrackFitter();
  gEve->AddToListTree(fitter, 1);
  gEve->AddElement(fitter, gEve->GetEventScene());

  AliEveTrackCounter* g_trkcnt = new AliEveTrackCounter("Primary Counter");
  gEve->AddToListTree(g_trkcnt, kFALSE);


  //==============================================================================
  // Final stuff
  //==============================================================================

  // A refresh to show proper window.
  gEve->Redraw3D(kTRUE);
  gSystem->ProcessEvents();

  // Register command to call on each event.
  AliEveEventManager::GetMaster()->AddNewEventCommand("on_new_event();");
  AliEveEventManager::GetMaster()->GotoEvent(0);

  gEve->EditElement(g_trkcnt);

  gEve->Redraw3D(kTRUE);
}

/******************************************************************************/

void on_new_event()
{
  AliEveTrackCounter* g_trkcnt = AliEveTrackCounter::fgInstance;
  g_trkcnt->Reset();
  g_trkcnt->SetEventId(AliEveEventManager::GetMaster()->GetEventId());

  if (g_esd_tracks_by_category_container != 0)
  {
    TEveElementList* cont = g_esd_tracks_by_category_container;

    // Here we expect several TEveTrackList containers.
    // First two have reasonable primaries (sigma-to-prim-vertex < 5).
    // Others are almost certainly secondaries.
    Int_t count = 1;
    TEveElement::List_i i = cont->BeginChildren();
    while (i != cont->EndChildren())
    {
      TEveTrackList* l = dynamic_cast<TEveTrackList*>(*i);
      if (l != 0)
      {
	g_trkcnt->RegisterTracks(l, (count <= 2));
	++count;
      }
      ++i;
    }

    // Set it to zero, so that we do not reuse an old one.
    g_esd_tracks_by_category_container = 0;
  }
  else
  {
    Warning("on_new_event", "g_esd_tracks_by_category_container not initialized.");
  }

  Double_t x[3] = { 0, 0, 0 };

  if (AliEveEventManager::HasESD())
  {
    AliESDEvent* esd = AliEveEventManager::AssertESD();
    esd->GetPrimaryVertex()->GetXYZ(x);

    TTimeStamp ts(esd->GetTimeStamp());
    TString win_title("Eve Main Window -- Timestamp: ");
    win_title += ts.AsString("s");
    win_title += "; Event # in ESD file: ";
    win_title += esd->GetEventNumberInFile();
    gEve->GetBrowser()->SetWindowName(win_title);
  }

  TEveElement* top = gEve->GetCurrentEvent();

  gMultiView->DestroyEventRPhi();
  if (gCenterProjectionsAtPrimaryVertex)
    gMultiView->SetCenterRPhi(x[0], x[1], x[2]);
  gMultiView->ImportEventRPhi(top);

  gMultiView->DestroyEventRhoZ();
  if (gCenterProjectionsAtPrimaryVertex)
    gMultiView->SetCenterRhoZ(x[0], x[1], x[2]);
  gMultiView->ImportEventRhoZ(top);
}
