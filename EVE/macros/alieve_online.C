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

  gROOT->LoadMacro("primary_vertex.C");
  gROOT->LoadMacro("esd_tracks.C");
  //  Disabled due to memory leaks
  //  gROOT->LoadMacro("trd_tracks.C++");
  gROOT->LoadMacro("trd_detectors.C++");

  gROOT->LoadMacro("its_clusters.C++");
  gROOT->LoadMacro("tpc_clusters.C++");
  gROOT->LoadMacro("tof_clusters.C++");
  gROOT->LoadMacro("hmpid_clusters.C++");
  gROOT->LoadMacro("emcal_digits.C++");

  gROOT->LoadMacro("acorde_raw.C");
  gROOT->LoadMacro("its_raw.C");
  gROOT->LoadMacro("tpc_raw.C");
  gROOT->LoadMacro("tof_raw.C");
  gROOT->LoadMacro("vzero_raw.C");

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
  if (AliEveEventManager::HasRawReader())
    its_raw();
  its_clusters();

  if (AliEveEventManager::HasRawReader())
    tpc_raw();
  tpc_clusters();

  if (AliEveEventManager::HasRawReader())
    tof_raw();
  tof_clusters();

  hmpid_clusters();

  if (AliEveEventManager::HasRawReader())
    acorde_raw();

  if (AliEveEventManager::HasRawReader())
    vzero_raw();

  emcal_digits();

  primary_vertex();
  esd_tracks();

  //  Disabled due to memory leaks
  //  if (AliEveEventManager::HasESDfriend()) trd_tracks();
  //  AliSysInfo::AddStamp("EveTRDTr");
  trd_detectors();

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
