// Import tracks from kinematics-tree / particle-stack.
// Preliminary/minimal solution.

Reve::TrackList* kine_tracks(Double_t min_pt=0.5, Double_t max_pt=100)
{
  AliRunLoader* rl =  Alieve::Event::AssertRunLoader();
  rl->LoadKinematics();
  AliStack* stack = rl->Stack();
  if (!stack) {
    Error("kine_tracks.C", "can not get kinematics.");
    return 0;
  }

  Reve::TrackList* cont = new Reve::TrackList("Kine Tracks"); 
  cont->SetMainColor(Color_t(6));
  Reve::TrackRnrStyle* rnrStyle = cont->GetRnrStyle();
  rnrStyle->fColor = 8;
  // !!! Watch the '-', apparently different sign convention then for ESD.
  rnrStyle->SetMagField( - gAlice->Field()->SolenoidField() );

  gReve->AddRenderElement(cont);

  Int_t count = 0;
  Int_t     N = stack->GetNtrack();
  for (Int_t i=0; i<N; ++i) {
    TParticle* p = stack->Particle(i);
    Double_t  pT = p->Pt();
    if (pT<min_pt || pT>max_pt) continue;

    ++count;
    Reve::Track* track = new Reve::Track(p, i, rnrStyle);
    track->SetName(Form("%s [%d]", p->GetName(), i));
    gReve->AddRenderElement(cont, track);
  }
  
  const Text_t* tooltip = Form("pT ~ (%.2lf, %.2lf), N=%d", min_pt, max_pt, count);
  cont->SetTitle(tooltip); // Not broadcasted automatically ...
  cont->UpdateItems();

  cont->MakeTracks();
  cont->MakeMarkers();
  gReve->Redraw3D();

  return cont;
}
