namespace Reve
{
class NLTProjector;
class GeoShapeRnrEl;
class RnrElement*;
}

Reve::NLTProjector  * proj = 0;
Reve::GeoShapeRnrEl * geom = 0;

void NLT_trackcount_init()
{
  Reve::LoadMacro("alieve_init.C");
  alieve_init(".", -1);

  Reve::LoadMacro("geom_gentle.C");

  Reve::LoadMacro("primary_vertex.C");
  Reve::LoadMacro("esd_tracks.C");
  Reve::LoadMacro("its_clusters.C+");
  Reve::LoadMacro("tpc_clusters.C+");

  Reve::Viewer* nv = gReve->SpawnNewViewer("NLT Projected");
  Reve::Scene*  ns = gReve->SpawnNewScene("NLT"); 
  nv->AddScene(ns);
  TGLViewer* v = nv->GetGLViewer();
  v->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
  TGLCameraMarkupStyle* mup = v->GetCameraMarkup();
  if(mup) mup->SetShow(kFALSE);

  Reve::TrackCounter* g_trkcnt = new Reve::TrackCounter("Primary Counter");
  gReve->AddToListTree(g_trkcnt, kFALSE);

  Reve::NLTProjector* p = new Reve::NLTProjector; proj = p;
  gReve->AddToListTree(p, kTRUE);
  gReve->AddRenderElement(proj, ns);

  // geometry
  Reve::GeoShapeRnrEl* gg = geom_gentle();
  geom = gg;

  // event
  Alieve::gEvent->AddNewEventCommand("on_new_event();");
  Alieve::gEvent->GotoEvent(0);

  gReve->Redraw3D(kTRUE);
}

/**************************************************************************/

void on_new_event()
{
  Reve::PointSet* itsc = its_clusters();
  itsc->SetMarkerColor(5);

  Reve::PointSet* tpcc = tpc_clusters();
  tpcc->SetMarkerColor(4);

  primary_vertex(1, 1);

  Reve::RenderElementList* cont = esd_tracks_vertex_cut();

  // Here we expect five TrackList containers.
  // First two have reasonable primaries (sigma-to-prim-vertex < 5).
  // Other three are almost certainly secondaries.
  Int_t count = 1;
  Reve::TrackCounter* g_trkcnt = Reve::TrackCounter::fgInstance;
  g_trkcnt->Reset();
  g_trkcnt->SetEventId(Alieve::gEvent->GetEventId());
  Reve::RenderElement::List_i i = cont->BeginChildren();
  while (i != cont->EndChildren()) {
    Reve::TrackList* l = dynamic_cast<Reve::TrackList*>(*i);
    if (l != 0) {
      l->SetLineWidth(2);
      g_trkcnt->RegisterTracks(l, (count <= 2));
      ++count;
    }
    ++i;
  }
  Reve::RenderElement* top = gReve->GetCurrentEvent();
  proj->DestroyElements();
  AliESDEvent* esd = Alieve::Event::AssertESD();
  Double_t x[3];
  esd->GetPrimaryVertex()->GetXYZ(x);
  proj->SetCenter(x[0], x[1], x[2]);

  // geom
  proj->ImportElements(geom);
  // event
  proj->ImportElements(top);
  // top->SetRnrState(kFALSE);
}

/**************************************************************************/

TParticle* id(Int_t label=0, Bool_t showParents=kTRUE)
{
  AliRunLoader* rl = Alieve::Event::AssertRunLoader();
  rl->LoadKinematics();
  AliStack* stack = rl->Stack();

  printf("Number primaries %d, all particles %d, label %d\n",
	 stack->GetNprimary(), stack->GetNtrack(), label);
  if (label < 0 || label >= stack->GetNtrack()) {
    printf("  Label exceeds available range.\n");
    return 0;
  }

  TParticle* part = stack->Particle(label);
  if(part != 0) {
    part->Print();
    if(showParents) {
      while(part->GetMother(0) >= 0) {
	part = stack->Particle(part->GetMother(0));
	part->Print();
      }
    }
  }
  return stack->Particle(label);
}
