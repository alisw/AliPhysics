// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

class TEveProjectionManager;
class TEveGeoShape;
class TEveUtil;

TEveProjectionManager  *proj = 0;
TEveGeoShape *geom = 0;

void trackcount_init()
{
  TEveUtil::LoadMacro("alieve_init.C");
  alieve_init(".", -1);

  TEveUtil::LoadMacro("geom_gentle.C");

  TEveUtil::LoadMacro("primary_vertex.C");
  TEveUtil::LoadMacro("esd_tracks.C");
  TEveUtil::LoadMacro("its_clusters.C+");
  TEveUtil::LoadMacro("tpc_clusters.C+");

  TEveViewer* nv = gEve->SpawnNewViewer("Projected View");
  TEveScene*  ns = gEve->SpawnNewScene("Projected Scene");
  nv->AddScene(ns);
  TGLViewer* v = nv->GetGLViewer();
  v->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
  TGLCameraMarkupStyle* mup = v->GetCameraMarkup();
  if(mup) mup->SetShow(kFALSE);

  AliEveTrackFitter* fitter = new AliEveTrackFitter();
  gEve->AddToListTree(fitter, 1);
  gEve->AddElement(fitter, gEve->GetEventScene());

  AliEveTrackCounter* g_trkcnt = new AliEveTrackCounter("Primary Counter");
  gEve->AddToListTree(g_trkcnt, kFALSE);

  TEveProjectionManager* pm = new TEveProjectionManager(); 
  proj = pm;
  gEve->AddToListTree(proj, kTRUE);
  gEve->AddElement(proj, ns);
  
  TEveProjectionAxes* axes = new TEveProjectionAxes(proj);
  axes->SetText("Projected Track-Count");
  axes->SetFontFile("comicbd");
  axes->SetFontSize(20);
  gEve->AddElement(axes, ns);
  gEve->AddToListTree(axes, kTRUE);

  // geometry
  TEveGeoShape* gg = geom_gentle();
  geom = gg;

  // event
  gAliEveEvent->AddNewEventCommand("on_new_event();");
  gAliEveEvent->GotoEvent(0);

  gEve->Redraw3D(kTRUE);
}

/******************************************************************************/

void on_new_event()
{
  try {
    TEvePointSet* itsc = its_clusters();
    itsc->SetMarkerColor(5);

    TEvePointSet* tpcc = tpc_clusters();
    tpcc->SetMarkerColor(4);
  }
  catch(TEveException& exc) {
    printf("Exception loading ITS/TPC clusters: %s\n", exc.Data());
  }

  primary_vertex(1, 1);

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
      l->SetLineWidth(2);
      g_trkcnt->RegisterTracks(l, (count <= 2));
      ++count;
    }
    ++i;
  }
  TEveElement* top = gEve->GetCurrentEvent();
  proj->DestroyElements();
  AliESDEvent* esd = AliEveEventManager::AssertESD();
  Double_t x[3];
  esd->GetPrimaryVertex()->GetXYZ(x);
  proj->SetCenter(x[0], x[1], x[2]);

  // geom
  proj->ImportElements(geom);
  // event
  proj->ImportElements(top);
  // top->SetRnrState(kFALSE);
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
