// Import daughters of track with given label.

Reve::TrackList*
kine_daughter_tracks(Int_t label, Reve::TrackList* cont = 0)
{
  using namespace Reve;

  gSystem->IgnoreSignal(kSigSegmentationViolation, true);
  if (label < 0) {
    Warning("kine_daughter_tracks", "label not set.");
    return 0;
  }
 
  //Bool_t create_cont;
  //create_cont = cont ? 0:1;
  AliRunLoader* rl =  Alieve::Event::AssertRunLoader();
  rl->LoadKinematics();
  AliStack* stack = rl->Stack();

  TParticle* p = stack->Particle(label);
  if (p->GetNDaughters()) 
  {
    if (cont == 0) 
    {
      printf("create container \n");
      cont =  
	new TrackList(Form("%d kine_daughter_tracks %d", p->GetNDaughters(), label));

      Reve::TrackRnrStyle* rnrStyle = cont->GetRnrStyle();
      // !!! Watch the '-', apparently different sign convention then for ESD.
      rnrStyle->SetMagField( - gAlice->Field()->SolenoidField() );
      char tooltip[1000];
      sprintf(tooltip,"%d, mother %d", p->GetNDaughters(), label);
      cont->SetTitle(tooltip);  
      cont->SelectByPt(0, 100);
      rnrStyle->fColor = 8;  
      rnrStyle->fMaxOrbs = 8;  
      gReve->AddRenderElement(cont);
    }

    for (int d=p->GetFirstDaughter(); d>0 && d<=p->GetLastDaughter(); ++d) 
    {	
      TParticle* dp = stack->Particle(d);
      Track* track = new Reve::Track(dp, d, cont->GetRnrStyle());  
      char form[1000];
      sprintf(form,"%s [%d]", dp->GetName(), d);
      track->SetName(form);
      gReve->AddRenderElement(cont, track);
    }

    cont->UpdateItems(); // update list tree
    cont->MakeTracks();
    cont->MakeMarkers();
    gReve->Redraw3D();
    return cont;
  }
  else 
  {
    printf("Particle %d does not have daughters.", label);
    return 0;
  }

}
