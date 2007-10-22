Reve::TrackList* primary_vertex_tracks()
{
  Reve::LoadMacro("esd_tracks.C");
  AliESDEvent   *esd = Alieve::Event::AssertESD();
  AliESDVertex *pv  = esd->GetPrimaryVertex();

  Reve::TrackList* cont = new Reve::TrackList("Tracks for Primary Vertex"); 
  cont->SetMainColor(Color_t(7));
  Reve::TrackRnrStyle* rnrStyle = cont->GetRnrStyle();
  rnrStyle->SetMagField( esd->GetMagneticField() );
  rnrStyle->fRnrFV = kTRUE;
  rnrStyle->fFVAtt->SetMarkerColor(2);
  gReve->AddRenderElement(cont);

  for (Int_t n=0; n<pv->GetNIndices(); n++)
  {
    AliESDtrack* at = esd->GetTrack(pv->GetIndices()[n]);
    Reve::Track* track = esd_make_track(rnrStyle, n, at, at);
    track->SetLineWidth(4);
    track->SetLineColor(cont->GetMainColor());
    track->SetLineStyle(7);
    gReve->AddRenderElement(track, cont);
  }

  cont->MakeTracks();
  gReve->Redraw3D();

  return cont;
}
