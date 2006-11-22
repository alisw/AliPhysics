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
    //PH The line below is replaced waiting for a fix in Root
    //PH which permits to use variable siza arguments in CINT
    //PH on some platforms (alphalinuxgcc, solariscc5, etc.)
    //PH    track->SetName(Form("%s [%d]", p->GetName(), i));
    char form[1000];
    sprintf(form,"%s [%d]", p->GetName(), i);
    track->SetName(form);
    gReve->AddRenderElement(cont, track);
  }
  
  //PH  const Text_t* tooltip = Form("pT ~ (%.2lf, %.2lf), N=%d", min_pt, max_pt, count);
  char tooltip[1000];
  sprintf(tooltip,"pT ~ (%.2lf, %.2lf), N=%d", min_pt, max_pt, count);
  cont->SetTitle(tooltip); // Not broadcasted automatically ...
  cont->UpdateItems();

  cont->MakeTracks();
  cont->MakeMarkers();
  gReve->Redraw3D();

  return cont;
}
