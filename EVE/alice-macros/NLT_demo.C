// $Id$

void NLT_demo(Int_t type = Reve::NLTProjection::RhoZ, Float_t distortion = 0.)
{
  Reve::LoadMacro("alieve_init.C");
  Reve::LoadMacro("esd_tracks.C");
  Reve::LoadMacro("its_hits.C");
  Reve::LoadMacro("NLT_geo_demo.C");

  Reve::NLTProjector* nlt = new  Reve::NLTProjector();
  nlt->SetProjection(type, distortion);

  // geometry
  make_geo(nlt);

  // event
  alieve_init();
  gReve->DisableRedraw();
  {
    // tacks
    Reve::TrackList* tl  =  esd_tracks();
    nlt->ProjectTrackList(tl);
    tl->SetStyle(2);
    // hits
    Reve::PointSet* ih = its_hits();
    nlt->ProjectPointSet(ih);
  }
  gReve->EnableRedraw();

  TGLViewer* glv = dynamic_cast<TGLViewer*>(gReve->GetGLCanvas()->GetViewer3D());
  glv->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
  gReve->Redraw3D();
}
