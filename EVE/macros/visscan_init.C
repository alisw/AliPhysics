// $Id: NLT_trackcount_init.C 24927 2008-04-04 13:46:04Z mtadel $
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

class TEveProjectionManager;
class TEveGeoShape;
class TEveUtil;

R__EXTERN TEveProjectionManager *gRPhiMgr;
R__EXTERN TEveProjectionManager *gRhoZMgr;

TEveGeoShape *gGeomGentle     = 0;
TEveGeoShape *gGeomGentleRPhi = 0;
TEveGeoShape *gGeomGentleRhoZ = 0;
TEveGeoShape *gGeomGentleTRD  = 0;

Bool_t gShowTRD = kFALSE;

void visscan_init()
{
  TEveUtil::LoadMacro("alieve_init.C");
  alieve_init(".", -1, 0, 0, 0, 0, kFALSE, kTRUE, kFALSE, kFALSE);

  AliEveTrackFitter* fitter = new AliEveTrackFitter();
  gEve->AddToListTree(fitter, 1);
  gEve->AddElement(fitter, gEve->GetEventScene());

  AliEveTrackCounter* g_trkcnt = new AliEveTrackCounter("Primary Counter");
  gEve->AddToListTree(g_trkcnt, kFALSE);

  // Geometry
  TEveUtil::LoadMacro("geom_gentle.C");
  gGeomGentle     = geom_gentle();
  gGeomGentleRPhi = geom_gentle_rphi(); gGeomGentleRPhi->IncDenyDestroy();
  gGeomGentleRhoZ = geom_gentle_rhoz(); gGeomGentleRhoZ->IncDenyDestroy();
  if (gShowTRD) {
    TEveUtil::LoadMacro("geom_gentle_trd.C");
    gGeomGentleTRD = geom_gentle_trd();
  }

  // Per event data
  TEveUtil::LoadMacro("primary_vertex.C");
  TEveUtil::LoadMacro("esd_V0_points.C");
  TEveUtil::LoadMacro("esd_V0.C");
  TEveUtil::LoadMacro("esd_cascade_points.C");
  TEveUtil::LoadMacro("esd_cascade.C");
  TEveUtil::LoadMacro("esd_tracks.C");
  TEveUtil::LoadMacro("its_clusters.C+");
  TEveUtil::LoadMacro("tpc_clusters.C+");
  TEveUtil::LoadMacro("trd_clusters.C+");
  TEveUtil::LoadMacro("tof_clusters.C+");

  // TEveLine::SetDefaultSmooth(1);

  TEveBrowser* browser = gEve->GetBrowser();
  browser->ShowCloseTab(kFALSE);

  gROOT->ProcessLine(".L SplitGLView.C+");
  browser->ExecPlugin("SplitGLView", 0, "new SplitGLView(gClient->GetRoot(), 600, 450, kTRUE)");

  browser->ShowCloseTab(kTRUE);

  browser->StartEmbedding(TRootBrowser::kBottom);
  new AliEveEventManagerWindow(AliEveEventManager::GetMaster());
  browser->StopEmbedding("EventCtrl");

  // Projections
  if (gRPhiMgr) {
    TEveProjectionAxes* a = new TEveProjectionAxes(gRPhiMgr);
    a->SetMainColor(kWhite);
    a->SetTitle("R-Phi");
    a->SetTitleSize(0.05);
    a->SetTitleFontName("comicbd");
    a->SetLabelSize(0.025);
    a->SetLabelFontName("comicbd");
    gEve->GetScenes()->FindChild("R-Phi Projection")->AddElement(a);
  }
  if (gRhoZMgr) {
    TEveProjectionAxes* a = new TEveProjectionAxes(gRhoZMgr);
    a->SetMainColor(kWhite);
    a->SetTitle("Rho-Z");
    a->SetTitleSize(0.05);
    a->SetTitleFontName("comicbd");
    a->SetLabelSize(0.025);
    a->SetLabelFontName("comicbd");
    gEve->GetScenes()->FindChild("Rho-Z Projection")->AddElement(a);
  }

  // Event
  AliEveEventManager::GetMaster()->AddNewEventCommand("on_new_event();");
  AliEveEventManager::GetMaster()->GotoEvent(0);

  gEve->EditElement(g_trkcnt);

  gEve->Redraw3D(kTRUE);
}

/******************************************************************************/

void on_new_event()
{
  try {
    TEvePointSet* itsc = its_clusters();
    if (itsc) {
      itsc->SetMarkerColor(5);
    }

    TEvePointSet* tpcc = tpc_clusters();
    if (tpcc) {
      tpcc->SetMarkerColor(4);
    }

    TEvePointSet* trdc = trd_clusters();
    if (trdc) {
      trdc->SetMarkerColor(7);
      trdc->SetMarkerStyle(4);
      trdc->SetMarkerSize(0.5);
    }

    TEvePointSet* tofc = tof_clusters();
    if (tofc) {
      tofc->SetMarkerColor(kOrange);
      tofc->SetMarkerStyle(4);
      tofc->SetMarkerSize(0.5);
    }
  }
  catch(TEveException& exc) {
    printf("Exception loading ITS/TPC clusters: %s\n", exc.Data());
  }

  primary_vertex();
  primary_vertex_ellipse();
  primary_vertex_spd();
  primary_vertex_ellipse_spd();

  esd_V0_points();
  esd_V0();

  esd_cascade_points();
  esd_cascade();

  AliEveTrackCounter* g_trkcnt = AliEveTrackCounter::fgInstance;
  g_trkcnt->Reset();
  g_trkcnt->SetEventId(AliEveEventManager::GetMaster()->GetEventId());

  TEveElementList* cont = esd_tracks_by_category();

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

  AliESDEvent* esd = AliEveEventManager::AssertESD();
  {
    TTimeStamp ts(esd->GetTimeStamp());
    TString win_title("Eve Main Window -- Timestamp: ");
    win_title += ts.AsString("s");
    win_title += "; Event # in ESD file: ";
    win_title += esd->GetEventNumberInFile();
    gEve->GetBrowser()->SetWindowName(win_title);
  }
  Double_t x[3];
  esd->GetPrimaryVertex()->GetXYZ(x);

  TEveElement* top = gEve->GetCurrentEvent();

  if (gRPhiMgr && top) {
    gRPhiMgr->DestroyElements();
    gRPhiMgr->SetCenter(x[0], x[1], x[2]);
    gRPhiMgr->ImportElements(gGeomGentleRPhi);
    if (gShowTRD) gRPhiMgr->ImportElements(gGeomGentleTRD);
    gRPhiMgr->ImportElements(top);
  }
  if (gRhoZMgr && top) {
    gRhoZMgr->DestroyElements();
    gRhoZMgr->SetCenter(x[0], x[1], x[2]);
    gRhoZMgr->ImportElements(gGeomGentleRhoZ);
    if (gShowTRD) gRhoZMgr->ImportElements(gGeomGentleTRD);
    gRhoZMgr->ImportElements(top);
  }

  gROOT->ProcessLine("SplitGLView::UpdateSummary()");
}
