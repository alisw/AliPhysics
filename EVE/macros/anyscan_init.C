// $Id: NLT_trackcount_init.C 24927 2008-04-04 13:46:04Z mtadel $
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

class AliEveMacroExecutor;

class TEveProjectionManager;
class TEveGeoShape;
class TEveUtil;

R__EXTERN TEveProjectionManager *gRPhiMgr;
R__EXTERN TEveProjectionManager *gRhoZMgr;

TEveGeoShape *gGeomGentle    = 0;
TEveGeoShape *gGeomGentleTRD = 0;

Bool_t gShowTRD = kFALSE;

void anyscan_init()
{
  TEveUtil::LoadMacro("alieve_init.C");
  alieve_init(".", -1);

  TEveLine::SetDefaultSmooth(1);

  TEveUtil::AssertMacro("VizDB_scan.C");


  AliEveTrackFitter* fitter = new AliEveTrackFitter();
  gEve->AddToListTree(fitter, 1);
  gEve->AddElement(fitter, gEve->GetEventScene());

  AliEveTrackCounter* g_trkcnt = new AliEveTrackCounter("Primary Counter");
  gEve->AddToListTree(g_trkcnt, kFALSE);


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

  // geometry
  TEveUtil::LoadMacro("geom_gentle.C");
  gGeomGentle = geom_gentle();
  if (gShowTRD) {
    TEveUtil::LoadMacro("geom_gentle_trd.C");
    gGeomGentleTRD = geom_gentle_trd();
  }


  AliEveMacroExecutor *exec = gAliEveEvent->GetExecutor();

  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "SIM Track",   "kine_tracks.C", "kine_tracks", "", kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "SIM Hit ITS", "its_hits.C",    "its_hits",    "", kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "SIM Hit TPC", "tpc_hits.C",    "tpc_hits",    "", kFALSE));

  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX", "primary_vertex.C", "primary_vertex"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC V0",   "esd_V0_points.C",  "esd_V0_points"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC V0",   "esd_V0.C",         "esd_V0"));
  // exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "esd_tracks.C",     ""));
  TEveUtil::LoadMacro("esd_tracks.C");

  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clus ITS", "its_clusters.C+", "its_clusters"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clus TPC", "tpc_clusters.C+", "tpc_clusters"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clus TRD", "trd_clusters.C+", "trd_clusters"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clus TOF", "tof_clusters.C+", "tof_clusters"));

  TEveBrowser* browser = gEve->GetBrowser();

  browser->StartEmbedding(TRootBrowser::kRight);
  AliEveMacroExecutorWindow* exewin = new AliEveMacroExecutorWindow(exec);
  browser->StopEmbedding("DataSelection");
  exewin->PopulateMacros();

  browser->StartEmbedding();
  new AliQAHistViewer(gClient->GetRoot(), 600, 400, kTRUE);
  browser->StopEmbedding("QA histograms");

  browser->GetTabRight()->SetTab(1);

  browser->StartEmbedding(TRootBrowser::kBottom);
  new AliEveEventManagerWindow;
  browser->StopEmbedding("EventCtrl");

  // event
  gAliEveEvent->AddNewEventCommand("on_new_event();");
  gAliEveEvent->GotoEvent(0);

  gEve->EditElement(g_trkcnt);

  gEve->Redraw3D(kTRUE);
}

/******************************************************************************/

void on_new_event()
{
  printf("on_new_event() entered ...\n");

  TEveElementList* cont = esd_tracks_vertex_cut();

  // Here we expect five TEveTrackList containers.
  // First two have reasonable primaries (sigma-to-prim-vertex < 5).
  // Other three are almost certainly secondaries.
  Int_t count = 1;
  AliEveTrackCounter* g_trkcnt = AliEveTrackCounter::fgInstance;
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
