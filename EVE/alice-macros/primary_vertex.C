// $Header$

TPolyMarker3D* make_vertex_marker(AliESDVertex* v, const Text_t* name)
{
  Double_t x[3], e[3];
  v->GetTruePos(x);
  v->GetSigmaXYZ(e);

  printf("%s: %f %f %f   -   %f %f %f\n", name,
	 x[0], x[1], x[2], e[0], e[1], e[2]);

  TPolyMarker3D* m = new TPolyMarker3D(8);
  m->SetName(name);
  m->SetPoint(0, x[0], x[1], x[2]);

  return m;
}

void register_vertex_marker(TPolyMarker3D* m)
{
  using namespace Reve;
  Color_t* colp = FindColorVar(m, "fMarkerColor");
  RenderElementObjPtr* rnrEl = new RenderElementObjPtr(m, *colp);
  gReve->AddRenderElement(rnrEl);
  gReve->Redraw3D();
}

void primary_vertex(Bool_t showSPD=kTRUE)
{
  AliESD* esd = Alieve::Event::AssertESD();

  AliESDVertex*  pv  = esd->GetPrimaryVertex();
  TPolyMarker3D* pvm = make_vertex_marker(pv, "Primary Vertex");
  pvm->SetMarkerStyle(5);
  pvm->SetMarkerColor(3);
  register_vertex_marker(pvm);

  if(showSPD) {
    AliESDVertex*  spdv  = esd->GetVertex();
    TPolyMarker3D* spdvm = make_vertex_marker(spdv, "SPD Vertex");
    spdvm->SetMarkerStyle(2);
    spdvm->SetMarkerColor(7);
    register_vertex_marker(spdvm);
  }
}
