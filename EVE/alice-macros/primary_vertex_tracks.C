TEveTrackList* primary_vertex_tracks()
{
  TEveUtil::LoadMacro("esd_tracks.C");
  AliESDEvent   *esd = Alieve::Event::AssertESD();
  AliESDVertex *pv  = esd->GetPrimaryVertex();

  TEveTrackList* cont = new TEveTrackList("Tracks for Primary Vertex"); 
  cont->SetMainColor(Color_t(7));
  TEveTrackPropagator* rnrStyle = cont->GetPropagator();
  rnrStyle->SetMagField( esd->GetMagneticField() );
  rnrStyle->fRnrFV = kTRUE;
  rnrStyle->fFVAtt->SetMarkerColor(2);
  gEve->AddElement(cont);

  for (Int_t n=0; n<pv->GetNIndices(); n++)
  {
    AliESDtrack* at = esd->GetTrack(pv->GetIndices()[n]);
    TEveTrack* track = esd_make_track(rnrStyle, n, at, at);
    track->SetLineWidth(4);
    track->SetLineColor(cont->GetMainColor());
    track->SetLineStyle(7);
    gEve->AddElement(track, cont);
  }

  cont->MakeTracks();
  gEve->Redraw3D();

  return cont;
}
