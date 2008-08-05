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

TEveGeoShape *gGeomGentle    = 0;
TEveGeoShape *gGeomGentleTRD = 0;

Bool_t gShowTRD = kFALSE;

void visscan_init()
{
  TEveUtil::LoadMacro("alieve_init.C");
  alieve_init(".", -1);

  TEveUtil::LoadMacro("geom_gentle.C");
  if (gShowTRD) TEveUtil::LoadMacro("geom_gentle_trd.C");

  TEveUtil::LoadMacro("primary_vertex.C");
  TEveUtil::LoadMacro("esd_V0_points.C");
  TEveUtil::LoadMacro("esd_V0.C");
  TEveUtil::LoadMacro("esd_tracks.C");
  TEveUtil::LoadMacro("its_clusters.C+");
  TEveUtil::LoadMacro("tpc_clusters.C+");
  TEveUtil::LoadMacro("trd_clusters.C+");
  TEveUtil::LoadMacro("tof_clusters.C+");

  TEveLine::SetDefaultSmooth(1);

  /*
  TEveViewer* nv = gEve->SpawnNewViewer("NLT Projected");
  TEveScene*  ns = gEve->SpawnNewScene("NLT");
  nv->AddScene(ns);
  TGLViewer* v = nv->GetGLViewer();
  v->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
  TGLCameraMarkupStyle* mup = v->GetCameraMarkup();
  if(mup) mup->SetShow(kFALSE);
  */

  AliEveTrackFitter* fitter = new AliEveTrackFitter();
  gEve->AddToListTree(fitter, 1);
  gEve->AddElement(fitter, gEve->GetEventScene());

  TEveTrackCounter* g_trkcnt = new TEveTrackCounter("Primary Counter");
  gEve->AddToListTree(g_trkcnt, kFALSE);


  // geometry
  gGeomGentle = geom_gentle();
  if (gShowTRD) gGeomGentleTRD = geom_gentle_trd();


  gROOT->ProcessLine(".L SplitGLView.C+");
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

  // event
  gAliEveEvent->AddNewEventCommand("on_new_event();");
  gAliEveEvent->GotoEvent(0);

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

  primary_vertex(1, 1);
  esd_V0_points();
  esd_V0();

  TEveElementList* cont = esd_tracks_vertex_cut();

  // Here we expect five TEveTrackList containers.
  // First two have reasonable primaries (sigma-to-prim-vertex < 5).
  // Other three are almost certainly secondaries.
  Int_t count = 1;
  TEveTrackCounter* g_trkcnt = TEveTrackCounter::fgInstance;
  g_trkcnt->Reset();
  g_trkcnt->SetEventId(gAliEveEvent->GetEventId());
  TEveElement::List_i i = cont->BeginChildren();
  while (i != cont->EndChildren()) {
    TEveTrackList* l = dynamic_cast<TEveTrackList*>(*i);
    if (l != 0) {
      // l->SetLineWidth(2);
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
    gRPhiMgr->ImportElements(gGeomGentle);
    if (gShowTRD) gRPhiMgr->ImportElements(gGeomGentleTRD);
    gRPhiMgr->ImportElements(top);
  }
  if (gRhoZMgr && top) {
    gRhoZMgr->DestroyElements();
    gRhoZMgr->SetCenter(x[0], x[1], x[2]);
    gRhoZMgr->ImportElements(gGeomGentle);
    if (gShowTRD) gRhoZMgr->ImportElements(gGeomGentleTRD);
    gRhoZMgr->ImportElements(top);
  }

  gROOT->ProcessLine("SplitGLView::UpdateSummary()");
}

/******************************************************************************/

TParticle* id(Int_t label=0, Bool_t showParents=kTRUE)
{
  AliRunLoader* rl = AliEveEventManager::AssertRunLoader();
  rl->LoadKinematics();
  AliStack* stack = rl->Stack();

  printf("Number primaries %d, all particles %d, label %d\n",
	 stack->GetNprimary(), stack->GetNtrack(), label);
  if (label < 0 || label >= stack->GetNtrack()) {
    printf("  Label exceeds available range.\n");
    return 0;
  }

  TParticle* part = stack->Particle(label);
  if (part != 0) {
    part->Print();
    if (showParents) {
      while (part->GetMother(0) >= 0) {
	part = stack->Particle(part->GetMother(0));
	part->Print();
      }
    }
  }
  return stack->Particle(label);
}
