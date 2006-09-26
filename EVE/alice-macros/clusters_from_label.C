// $Id$

void clusters_from_label(Int_t label=0)
{
  AliESD* esd = Alieve::Event::AssertESD();
  Reve::PointSet* clusters = new Reve::PointSet(64);
  clusters->SetOwnIds(kTRUE);

  for (Int_t n=0; n<esd->GetNumberOfTracks(); n++) {
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
  clusters->SetMarkerSize(1);
  clusters->SetMarkerColor(4);
  clusters->SetName(Form("Clusters lab=%d", label));

  using namespace Reve;
  Color_t* colp = FindColorVar(clusters, "fMarkerColor");
  RenderElementObjPtr* rnrEl = new RenderElementObjPtr(clusters, *colp);
  gReve->AddRenderElement(rnrEl);
  gReve->Redraw3D();
}
