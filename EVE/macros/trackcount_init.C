void trackcount_init()
{
  Reve::LoadMacro("alieve_init.C");
  alieve_init(".", -1);

  Reve::LoadMacro("primary_vertex.C");
  Reve::LoadMacro("esd_tracks.C");
  Reve::LoadMacro("its_clusters.C+");
  Reve::LoadMacro("tpc_clusters.C+");

  {
    TGLViewer* glv = (TGLViewer *)gReve->GetGLViewer();
    glv->SetIgnoreSizesOnUpdate(kTRUE);
    // The size of ortho cameras can not be set in advance.
    glv->SetOrthoCamera(TGLViewer::kCameraOrthoXOY, -0.1, 0.1, 0.1, -0.1);
    glv->SetOrthoCamera(TGLViewer::kCameraOrthoZOY, -22, 22, 22, -22);
  }

  Reve::TrackCounter* g_trkcnt = new Reve::TrackCounter("Primary Counter");
  gReve->AddGlobalRenderElement(g_trkcnt);

  Alieve::gEvent->AddNewEventCommand("on_new_event();");
  Alieve::gEvent->GotoEvent(0);

  gReve->Redraw3D(kTRUE);
}

void on_new_event()
{
  Reve::PointSet* itsc = its_clusters();
  itsc->SetMarkerColor(5);

  Reve::PointSet* tpcc = tpc_clusters();
  tpcc->SetMarkerColor(4);

  primary_vertex(1, 1);

  Reve::RenderElementList* cont = esd_tracks_vertex_cut();
  TGListTree* lt = gReve->GetListTree();
  TGListTreeItem* ti = cont->FindListTreeItem(lt);
  lt->OpenItem(ti);

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
}

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
