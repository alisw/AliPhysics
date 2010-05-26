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
  
  if (gSystem->Getenv("ALICE_ROOT") != 0)
  {
    gInterpreter->AddIncludePath(Form("%s/MUON", gSystem->Getenv("ALICE_ROOT")));
    gInterpreter->AddIncludePath(Form("%s/MUON/mapping", gSystem->Getenv("ALICE_ROOT")));
  }
  
  TEveUtil::AssertMacro("VizDB_scan.C");

  TEveBrowser *browser = gEve->GetBrowser();
  browser->ShowCloseTab(kFALSE);

  // Gentle-geom loading changes gGeoManager.
  TEveGeoManagerHolder mgrRestore;

  AliEveMultiView *multiView = new AliEveMultiView;

  TEveUtil::LoadMacro("geom_gentle.C");
  multiView->InitGeomGentle(geom_gentle(),
                             geom_gentle_rphi(), 
                             geom_gentle_rhoz());

  TEveUtil::LoadMacro("geom_gentle_trd.C");
  multiView->InitGeomGentleTrd(geom_gentle_trd());

  TEveUtil::LoadMacro("geom_gentle_muon.C");
  multiView->InitGeomGentleMuon(geom_gentle_muon(), kFALSE, kTRUE);

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
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clus TRD",   "trd_clusters.C++",   "trd_clusters"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clus TOF",   "tof_clusters.C++",   "tof_clusters"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clus HMPID", "hmpid_clusters.C++", "hmpid_clusters"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clus MUON",  "muon_clusters.C++",  "muon_clusters"));

  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clus TOF",   "emcal_digits.C++",   "emcal_digits"));

  exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW ITS",     "its_raw.C",     "its_raw"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW TPC",     "tpc_raw.C",     "tpc_raw"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW TOF",     "tof_raw.C",     "tof_raw"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW VZERO",   "vzero_raw.C",   "vzero_raw"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW ACORDE",  "acorde_raw.C",  "acorde_raw"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW MUON",    "muon_raw.C++",  "muon_raw"));

  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Track", "esd_tracks.C", "esd_tracks",             "", kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Track", "esd_tracks.C", "esd_tracks_MI",          "", kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Track", "esd_tracks.C", "esd_tracks_by_category", "", kTRUE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Track MUON", "esd_muon_tracks.C++", "esd_muon_tracks", "kTRUE,kFALSE", kTRUE));

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

  TGLViewer *glv = multiView->Get3DView()->GetGLViewer();
  glv->CurrentCamera().RotateRad(-0.4, 1);
  glv->DoDraw();
}

//   multiView->Get3DView()->GetGLViewer()->CurrentCamera().RotateRad(-1, 1)

Int_t      g_pic_id  = 0;
Int_t      g_pic_max = 100;
TTimeStamp g_pic_prev(0, 0);

void alieve_online_on_new_event()
{
  AliESDEvent* esd = AliEveEventManager::AssertESD();
  Double_t x[3];
  esd->GetPrimaryVertex()->GetXYZ(x);

  TEveElement* top = gEve->GetCurrentEvent();

  AliEveMultiView *multiView = AliEveMultiView::Instance();

  TGLViewer *glv = (dynamic_cast<TEveViewer*>(gEve->GetViewers()->FindChild("3D View")))->GetGLViewer();
  
  if(gEve->GetScenes()->FirstChild()->FindChild("Gentle MUON"))
  {
    if (esd->GetNumberOfMuonTracks() == 0 && !gEve->GetKeepEmptyCont())
    {
      gEve->GetScenes()->FirstChild()->FindChild("Gentle MUON")->SetRnrChildren(kFALSE);
      
      if(gEve->GetEventScene()->FirstChild()->FindChild("MUON Clusters"))
        gEve->GetEventScene()->FirstChild()->FindChild("MUON Clusters")->SetRnrSelf(kFALSE);
      if(gEve->GetEventScene()->FirstChild()->FindChild("MUON Raw digits"))
        gEve->GetEventScene()->FirstChild()->FindChild("MUON Raw digits")->SetRnrChildren(kFALSE);

      gEve->FullRedraw3D(kTRUE);
      glv->CurrentCamera().RotateRad(-0.4, -1.8);
    }
    else
    {
      gEve->GetScenes()->FirstChild()->FindChild("Gentle MUON")->SetRnrChildren(kTRUE);

      if(gEve->GetEventScene()->FirstChild()->FindChild("MUON Clusters"))
        gEve->GetEventScene()->FirstChild()->FindChild("MUON Clusters")->SetRnrSelf(kTRUE);
      if(gEve->GetEventScene()->FirstChild()->FindChild("MUON Raw digits"))
        gEve->GetEventScene()->FirstChild()->FindChild("MUON Raw digits")->SetRnrChildren(kTRUE);

      gEve->FullRedraw3D(kTRUE);
      glv->CurrentCamera().RotateRad(-0.4, 1);
    }
  }

  glv->DoDraw();

  multiView->DestroyEventRPhi();
  if (gCenterProjectionsAtPrimaryVertex)
    multiView->SetCenterRPhi(x[0], x[1], x[2]);
  multiView->ImportEventRPhi(top);

  multiView->DestroyEventRhoZ();
  if (gCenterProjectionsAtPrimaryVertex)
    multiView->SetCenterRhoZ(x[0], x[1], x[2]);
  multiView->ImportEventRhoZ(top);

  // Register image to amore.
  const TString pichost("aldaqacrs3");
  TTimeStamp now;
  Double_t delta = now.AsDouble() - g_pic_prev.AsDouble();

  printf("Pre image dump: host='%s', delta=%f.\n",
	 gSystem->HostName(), delta);

  if (pichost == gSystem->HostName() && delta >= 30)
  {
    TString id;      id.Form("online-viz-%03d", g_pic_id);
    TString pic(id); pic += ".png";

    printf("In image dump: file='%s'.\n", pic.Data());

    gEve->GetBrowser()->RaiseWindow();
    gEve->FullRedraw3D();
    gSystem->ProcessEvents();

    Int_t status;

    status = gSystem->Exec(TString::Format("xwd -id %u | convert - %s",
			   gEve->GetBrowser()->GetId(), pic.Data()));

    printf("Post capture -- status=%d.\n", status);

    status = gSystem->Exec(TString::Format("SendImageToAmore %s %s %d",
		          id.Data(), pic.Data(),
		          AliEveEventManager::AssertRawReader()->GetRunNumber()));

    printf("Post AMORE reg -- status=%d, run=%d.\n", status,
	   AliEveEventManager::AssertRawReader()->GetRunNumber());

    if (++g_pic_id >= g_pic_max)
      g_pic_id = 0;
    g_pic_prev.Set();
  }
}
