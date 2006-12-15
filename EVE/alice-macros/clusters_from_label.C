// $Id$

Reve::PointSet* clusters_from_label(Int_t label=0)
{
  AliESD* esd = Alieve::Event::AssertESD();
  Reve::PointSet* clusters = new Reve::PointSet(64);
  clusters->SetOwnIds(kTRUE);

  for (Int_t n=0; n<esd->GetNumberOfTracks(); n++)
  {
    AliESDtrack* at = esd->GetTrack(n);
    if (at->GetLabel() == label) {
      const AliTrackPointArray* pArr = at->GetTrackPointArray();
      if (pArr == 0) {
	Warning("clusters_from_label", "TrackPointArray not stored with ESD track.");
	continue;
      }
      Int_t np =  pArr->GetNPoints();
      const Float_t* x = pArr->GetX();
      const Float_t* y = pArr->GetY();
      const Float_t* z = pArr->GetZ();
      for (Int_t i=0; i<np; ++i) {
	clusters->SetNextPoint(x[i], y[i], z[i]);
	AliTrackPoint *atp = new AliTrackPoint;
	pArr->GetPoint(*atp, i);
	clusters->SetPointId(atp);
      }
    }
  }
  clusters->SetMarkerStyle(2);
  clusters->SetMarkerSize(0.5);
  clusters->SetMarkerColor(4);
  //PH The line below is replaced waiting for a fix in Root
  //PH which permits to use variable siza arguments in CINT
  //PH on some platforms (alphalinuxgcc, solariscc5, etc.)
  //PH  clusters->SetName(Form("Clusters lab=%d", label));
  char form[1000];
  sprintf(form,"Clusters lab=%d", label);
  clusters->SetName(form);

  using namespace Reve;
  gReve->AddRenderElement(clusters);
  gReve->Redraw3D();

  return clusters;
}
