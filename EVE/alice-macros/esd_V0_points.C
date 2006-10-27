// $Id$

Reve::PointSet*
esd_V0_points()
{
  AliESD* esd = Alieve::Event::AssertESD();

  Int_t NV0s = esd->GetNumberOfV0s();
  Reve::PointSet* points = new Reve::PointSet("V0 CA points", NV0s);

  for (Int_t n =0; n<NV0s; n++) {
    AliESDv0* av = esd->GetV0(n);
    points->SetNextPoint(av->GetXr(0), av->GetXr(1), av->GetXr(2));
    points->SetPointId(av);
  }

  points->SetTitle(Form("N=%d", points->Size()));
  points->SetMarkerStyle(4);
  points->SetMarkerSize(1);
  points->SetMarkerColor((Color_t)30);

  gReve->AddRenderElement(points);
  gReve->Redraw3D();

  return points;
}
