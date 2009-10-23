/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

class TEveProjectionManager;
class TEveGeoShape;
class TEveUtil;

Bool_t gCenterProjectionsAtPrimaryVertex = kFALSE;


void alieve_online_init()
{
  if (gROOT->LoadMacro("MultiView.C++") != 0)
  {
    gEnv->SetValue("Root.Stacktrace", "no");
    Fatal("alieve_online.C", "Failed loading MultiView.C in compiled mode.");
  }

  TEveUtil::AssertMacro("VizDB_scan.C");

  TEveBrowser *browser = gEve->GetBrowser();
  browser->ShowCloseTab(kFALSE);

  // Gentle-geom loading changes gGeoManager.
  TEveGeoManagerHolder mgrRestore;

  gMultiView = new MultiView;

  TEveUtil::LoadMacro("geom_gentle.C");
  gMultiView->InitGeomGentle(geom_gentle(),
                             geom_gentle_rphi(), 
                             geom_gentle_rhoz());

  // See visscan_init.C for how to add TRD / MUON geometry.

  //============================================================================
  // Standard macros to execute -- not all are enabled by default.
  //============================================================================

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

  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clus ITS",   "its_clusters.C++",   "its_clusters"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clus TPC",   "tpc_clusters.C++",   "tpc_clusters"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clus TOF",   "tof_clusters.C++",   "tof_clusters"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clus HMPID", "hmpid_clusters.C++", "hmpid_clusters"));

  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clus TOF",   "emcal_digits.C++",   "emcal_digits"));

  exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW ITS",     "its_raw.C",     "its_raw"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW TPC",     "tpc_raw.C",     "tpc_raw"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW TOF",     "tof_raw.C",     "tof_raw"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW VZERO",   "vzero_raw.C",   "vzero_raw"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW ACORDE",  "acorde_raw.C",  "acorde_raw"));

  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Track", "esd_tracks.C", "esd_tracks",             "", kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Track", "esd_tracks.C", "esd_tracks_MI",          "", kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Track", "esd_tracks.C", "esd_tracks_by_category", "", kTRUE));

  // ???
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC TRD", "trd_detectors.C++", "trd_detectors",         "", kFALSE));
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

  TString autoRun(gSystem->Getenv("ONLINERECO_AUTORUN"));
  if (autoRun == "1" || autoRun.CompareTo("true", TString::kIgnoreCase) == 0)
  {
    AliEveEventManager::GetMaster()->SetAutoLoad(kTRUE);
  }

  {
    TGTab *tab = gEve->GetBrowser()->GetTab(2);

    TGHorizontalFrame *hf = (TGHorizontalFrame*) tab->GetParent();
    TGVerticalFrame   *vf = (TGVerticalFrame*)   hf ->GetParent();

    hf->Resize(hf->GetWidth(), hf->GetHeight() + 80);
    vf->Layout();
  }

  gEve->FullRedraw3D(kTRUE);

  TGLViewer *glv = gMultiView->f3DView->GetGLViewer();
  glv->CurrentCamera().RotateRad(-0.4, 1);
  glv->DoDraw();
}

void alieve_online_on_new_event()
{
  AliESDEvent* esd = AliEveEventManager::AssertESD();
  Double_t x[3];
  esd->GetPrimaryVertex()->GetXYZ(x);

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
