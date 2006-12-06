Reve::TrackCounter* g_trkcnt = 0;

void trackcount_init()
{
  Reve::LoadMacro("alieve_init.C");
  alieve_init(".", -1);

  {
    TGLViewer* glv = dynamic_cast<TGLViewer*>(gReve->GetGLCanvas()->GetViewer3D());
    glv->SetIgnoreSizesOnUpdate(kTRUE);
    // The size of ortho cameras can not be set in advance.
    // glv->SetOrthoCamera(TGLViewer::kCameraOrthoXOY, -0.1, 0.1, 0.1, -0.1);
    // glv->SetOrthoCamera(TGLViewer::kCameraOrthoZOY, -22, 22, 22, -22);
  }

  g_trkcnt = new Reve::TrackCounter("Primary Counter");
  gReve->AddGlobalRenderElement(g_trkcnt);

  Alieve::gEvent->AddNewEventCommand("on_new_event();");
  Alieve::gEvent->GotoEvent(0);
}

void on_new_event()
{
  Reve::LoadMacro("primary_vertex.C");
  Reve::LoadMacro("esd_tracks.C");
  Reve::LoadMacro("its_hits.C");
  Reve::LoadMacro("tpc_hits.C");

  Reve::PointSet* p = its_hits();
  p->SetMarkerStyle(4);

  primary_vertex(1, 1);

  Reve::RenderElementList* cont = esd_tracks_vertex_cut();
  TGListTree* lt = gReve->GetListTree();
  TGListTreeItem* ti = cont->FindListTreeItem(lt);
  lt->OpenItem(ti);

  // Here we expect five TrackList containers.
  // First two have reasonable primaries (sigma-to-prim-vertex < 5).
  // Other three are almost certainly secondaries.
  Int_t count = 1;
  g_trkcnt->Reset();
  Reve::RenderElement::List_i i = cont->BeginChildren();
  while (i != cont->EndChildren()) {
    Reve::TrackList* l = dynamic_cast<Reve::TrackList*>(*i);
    if (l != 0) {
      l->SetWidth(2);
      g_trkcnt->RegisterTracks(l, (count <= 2));
      ++count;
    }
    ++i;
  }
}

void id(Int_t label=0, Bool_t showParents=kFALSE)
{
  AliRunLoader* rl = Alieve::Event::AssertRunLoader();
  rl->LoadKinematics();
  AliStack* stack = rl->Stack();

  printf("Number primaries %d %d\n", stack->GetNprimary(), stack->GetNtrack());

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
}
