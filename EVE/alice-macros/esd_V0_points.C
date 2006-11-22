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

  //PH The line below is replaced waiting for a fix in Root
  //PH which permits to use variable siza arguments in CINT
  //PH on some platforms (alphalinuxgcc, solariscc5, etc.)
  //PH  points->SetTitle(Form("N=%d", points->Size()));
  char form[1000];
  sprintf(form,"N=%d", points->Size());
  points->SetTitle(form);
  points->SetMarkerStyle(4);
  points->SetMarkerSize(1);
  points->SetMarkerColor((Color_t)30);

  gReve->AddRenderElement(points);
  gReve->Redraw3D();

  return points;
}
