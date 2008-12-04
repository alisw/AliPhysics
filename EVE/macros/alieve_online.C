/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

TEveGeoShape *gGeomGentle = 0;

void alieve_online_init()
{
  gROOT->LoadMacro("geom_gentle.C");

  gROOT->LoadMacro("primary_vertex.C");
  gROOT->LoadMacro("esd_tracks.C");
  gROOT->LoadMacro("trd_tracks.C++");
  gROOT->LoadMacro("trd_detectors.C++");

  gROOT->LoadMacro("its_clusters.C++");
  gROOT->LoadMacro("tpc_clusters.C++");
  gROOT->LoadMacro("hmpid_clusters.C++");

  gROOT->LoadMacro("acorde_raw.C");
  gROOT->LoadMacro("its_raw.C");
  gROOT->LoadMacro("tpc_raw.C");

  TEveUtil::AssertMacro("VizDB_scan.C");

  // Temp fix !!!
  TGeoManager *man = gGeoManager;
  gGeomGentle = geom_gentle();
  // Temp fix !!!
  gGeoManager = man;

  gROOT->ProcessLine(".L SplitGLView.C++g"); // !!!! debug-mode
  TEveBrowser* browser = gEve->GetBrowser();
  browser->ExecPlugin("SplitGLView", 0, "new SplitGLView(gClient->GetRoot(), 600, 450, kTRUE)");

  if (gRPhiMgr) {
    TEveProjectionAxes* a = new TEveProjectionAxes(gRPhiMgr);
    a->SetNumTickMarks(3);
    a->SetText("R-Phi");
    a->SetFontFile("comicbd");
    a->SetFontSize(10);
    gEve->GetScenes()->FindChild("R-Phi Projection")->AddElement(a);
  }
  if (gRhoZMgr) {
    TEveProjectionAxes* a = new TEveProjectionAxes(gRhoZMgr);
    a->SetNumTickMarks(3);
    a->SetText("Rho-Z");
    a->SetFontFile("comicbd");
    a->SetFontSize(10);
    gEve->GetScenes()->FindChild("Rho-Z Projection")->AddElement(a);
  }

  TEveBrowser* browser = gEve->GetBrowser();

  browser->StartEmbedding(TRootBrowser::kBottom);
  new AliEveEventManagerWindow(AliEveEventManager::GetMaster());
  browser->StopEmbedding("EventCtrl");

  gEve->Redraw3D(kTRUE);
}

void alieve_online_on_new_event()
{
  if (AliEveEventManager::HasRawReader())
    its_raw();
  its_clusters();

  if (AliEveEventManager::HasRawReader())
    tpc_raw();
  tpc_clusters();

  hmpid_clusters();

  if (AliEveEventManager::HasRawReader())
    acorde_raw();

  primary_vertex();
  esd_tracks();

  if (AliEveEventManager::HasESDfriend()) trd_tracks();
  trd_detectors();

  AliESDEvent* esd = AliEveEventManager::AssertESD();
  Double_t x[3];
  esd->GetPrimaryVertex()->GetXYZ(x);

  TEveElement* top = gEve->GetCurrentEvent();

  if (gRPhiMgr && top) {
    gRPhiMgr->DestroyElements();
    gRPhiMgr->SetCenter(x[0], x[1], x[2]);
    gRPhiMgr->ImportElements(gGeomGentle);
    gRPhiMgr->ImportElements(top);
  }
  if (gRhoZMgr && top) {
    gRhoZMgr->DestroyElements();
    gRhoZMgr->SetCenter(x[0], x[1], x[2]);
    gRhoZMgr->ImportElements(gGeomGentle);
    gRhoZMgr->ImportElements(top);
  }

  gROOT->ProcessLine("SplitGLView::UpdateSummary()");
}
