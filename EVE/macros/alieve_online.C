/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

class TEveProjectionManager;
class TEveGeoShape;
class TEveUtil;

TEveGeoShape *gGeomGentle = 0;

TEveScene *gRPhiGeomScene  = 0;
TEveScene *gRhoZGeomScene  = 0;
TEveScene *gRPhiEventScene = 0;
TEveScene *gRhoZEventScene = 0;

TEveProjectionManager *gRPhiMgr = 0;
TEveProjectionManager *gRhoZMgr = 0;

TEveViewer *g3DView   = 0;
TEveViewer *gRPhiView = 0;
TEveViewer *gRhoZView = 0;

Bool_t gCenterProjectionsAtPrimaryVertex = kFALSE;

void alieve_online_init()
{
  gROOT->LoadMacro("geom_gentle.C");

  gROOT->LoadMacro("primary_vertex.C");
  gROOT->LoadMacro("esd_tracks.C");
  gROOT->LoadMacro("trd_tracks.C++");
  gROOT->LoadMacro("trd_detectors.C++");

  gROOT->LoadMacro("its_clusters.C++");
  gROOT->LoadMacro("tpc_clusters.C++");
  gROOT->LoadMacro("tof_clusters.C++");
  gROOT->LoadMacro("hmpid_clusters.C++");

  gROOT->LoadMacro("acorde_raw.C");
  gROOT->LoadMacro("its_raw.C");
  gROOT->LoadMacro("tpc_raw.C");
  gROOT->LoadMacro("tof_raw.C");

  TEveUtil::AssertMacro("VizDB_scan.C");

  TEveBrowser         *browser = gEve->GetBrowser();
  browser->ShowCloseTab(kFALSE);

  // Temp fix !!!
  TGeoManager *man = gGeoManager;
  gGeomGentle = geom_gentle();
  gGeomGentleRPhi = geom_gentle_rphi(); gGeomGentleRPhi->IncDenyDestroy();
  gGeomGentleRhoZ = geom_gentle_rhoz(); gGeomGentleRhoZ->IncDenyDestroy();
  // Temp fix !!!
  gGeoManager = man;

  // Scenes
  gRPhiGeomScene  = gEve->SpawnNewScene("RPhi Geometry",
                    "Scene holding projected geometry for the RPhi view.");
  gRhoZGeomScene  = gEve->SpawnNewScene("RhoZ Geometry",
		    "Scene holding projected geometry for the RhoZ view.");
  gRPhiEventScene = gEve->SpawnNewScene("RPhi Event Data",
		    "Scene holding projected geometry for the RPhi view.");
  gRhoZEventScene = gEve->SpawnNewScene("RhoZ Event Data",
		    "Scene holding projected geometry for the RhoZ view.");

  gRPhiMgr = new TEveProjectionManager();
  gRPhiMgr->SetProjection(TEveProjection::kPT_RPhi);
  gEve->AddToListTree(gRPhiMgr, kFALSE);
  {
    TEveProjectionAxes* a = new TEveProjectionAxes(gRPhiMgr);
    a->SetMainColor(kWhite);
    a->SetTitle("R-Phi");
    a->SetTitleSize(0.05);
    a->SetTitleFont(102);
    a->SetLabelSize(0.025);
    a->SetLabelFont(102);
    gRPhiGeomScene->AddElement(a);
  }
  gRPhiMgr->SetCurrentDepth(-10);
  gRPhiMgr->ImportElements(gGeomGentleRPhi, gRPhiGeomScene);
  gRPhiMgr->SetCurrentDepth(0);

  gRhoZMgr = new TEveProjectionManager();
  gRhoZMgr->SetProjection(TEveProjection::kPT_RhoZ);
  gEve->AddToListTree(gRhoZMgr, kFALSE);
  {
    TEveProjectionAxes* a = new TEveProjectionAxes(gRhoZMgr);
    a->SetMainColor(kWhite);
    a->SetTitle("Rho-Z");
    a->SetTitleSize(0.05);
    a->SetTitleFont(102);
    a->SetLabelSize(0.025);
    a->SetLabelFont(102);
    gRhoZGeomScene->AddElement(a);
  }
  gRhoZMgr->SetCurrentDepth(-10);
  gRhoZMgr->ImportElements(gGeomGentleRhoZ, gRhoZGeomScene);
  gRhoZMgr->SetCurrentDepth(0);

  // Viewers
  TEveWindowSlot *slot = 0;
  TEveWindowPack *pack = 0;

  slot = TEveWindow::CreateWindowInTab(browser->GetTabRight());
  pack = slot->MakePack();
  pack->SetElementName("Multi View");
  pack->SetHorizontal();
  pack->SetShowTitleBar(kFALSE);
  pack->NewSlot()->MakeCurrent();
  g3DView = gEve->SpawnNewViewer("3D View", "");
  g3DView->AddScene(gEve->GetGlobalScene());
  g3DView->AddScene(gEve->GetEventScene());

  pack = pack->NewSlot()->MakePack();
  pack->SetShowTitleBar(kFALSE);
  pack->NewSlot()->MakeCurrent();
  gRPhiView = gEve->SpawnNewViewer("RPhi View", "");
  gRPhiView->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
  gRPhiView->AddScene(gRPhiGeomScene);
  gRPhiView->AddScene(gRPhiEventScene);

  pack->NewSlot()->MakeCurrent();
  gRhoZView = gEve->SpawnNewViewer("RhoZ View", "");
  gRhoZView->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
  gRhoZView->AddScene(gRhoZGeomScene);
  gRhoZView->AddScene(gRhoZEventScene);


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

  if (AliEveEventManager::HasRawReader())
    tof_raw();
  tof_clusters();

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

  if (gRPhiMgr && top)
  {
    gRPhiEventScene->DestroyElements();
    if (gCenterProjectionsAtPrimaryVertex)
      gRPhiMgr->SetCenter(x[0], x[1], x[2]);
    gRPhiMgr->ImportElements(top, gRPhiEventScene);
  }
  if (gRhoZMgr && top)
  {
    gRhoZEventScene->DestroyElements();
    if (gCenterProjectionsAtPrimaryVertex)
      gRhoZMgr->SetCenter(x[0], x[1], x[2]);
    gRhoZMgr->ImportElements(top, gRhoZEventScene);
  }
}
