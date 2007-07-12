Reve::TrackList* primary_vertex_tracks()
{
  Reve::LoadMacro("esd_tracks.C");
  AliESDEvent   *esd = Alieve::Event::AssertESD();
  AliESDVertex *pv  = esd->GetPrimaryVertex();

  Reve::TrackList* cont = new Reve::TrackList("Tracks for Primary Vertex"); 
  cont->SetMainColor(Color_t(7));
  Reve::TrackRnrStyle* rnrStyle = cont->GetRnrStyle();
  rnrStyle->SetMagField( esd->GetMagneticField() );

  gReve->AddRenderElement(cont);

  for (Int_t n=0; n<pv->GetNIndices(); n++)
  {
    AliESDtrack* at = esd->GetTrack(pv->GetIndices()[n]);
    Reve::Track* track = esd_make_track(rnrStyle, n, at, at);
    track->SetLineWidth(4);
    track->SetLineStyle(7);
    gReve->AddRenderElement(cont, track);
  }

  cont->MakeTracks();
  cont->MakeMarkers();
  gReve->Redraw3D();

  return cont;
}
