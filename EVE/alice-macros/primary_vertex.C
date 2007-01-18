// $Header$

TPolyMarker3D* make_vertex_marker(AliESDVertex* v, const Text_t* name)
{
  Double_t x[3], e[3];
  v->GetXYZ(x);
  v->GetSigmaXYZ(e);

  printf("%16s: %f %f %f   -   %f %f %f\n", name,
	 x[0], x[1], x[2], e[0], e[1], e[2]);

  TPolyMarker3D* m = new TPolyMarker3D(1);
  m->SetName(name);
  m->SetPoint(0, x[0], x[1], x[2]);

  return m;
}

Reve::BoxSet* make_vertex_boxes(AliESDVertex* v)
{
  Double_t x[3], e[3];
  v->GetTruePos(x);
  v->GetSigmaXYZ(e);

  Reve::BoxSet* bs;

  bs = new BoxSet("+- 10 x 10 x 20mm");
  bs->SetRenderMode(Reve::BoxSet::RM_Line);
  bs->AddBox(Reve::Box(-1, x[0], x[1], x[2], 1, 1, 2));
  bs->SetMainColor((Color_t) 2);
  gReve->AddRenderElement(bs);

  bs = new BoxSet("+- 30 sigma_r x 10 sigma_z");
  bs->SetRenderMode(Reve::BoxSet::RM_Line);
  bs->AddBox(Reve::Box(-1, x[0], x[1], x[2], 30*e[0], 30*e[1], 10*e[2]));
  bs->SetMainColor((Color_t) 3);
  gReve->AddRenderElement(bs);

  gReve->Redraw3D();
}

void register_vertex_marker(TPolyMarker3D* m)
{
  using namespace Reve;
  Color_t* colp = FindColorVar(m, "fMarkerColor");
  RenderElementObjPtr* rnrEl = new RenderElementObjPtr(m, *colp);
  gReve->AddRenderElement(rnrEl);
  gReve->Redraw3D();
}

void primary_vertex(Bool_t showSPD=kTRUE, Bool_t showBoxes=kFALSE)
{
  AliESD* esd = Alieve::Event::AssertESD();

  AliESDVertex*  pv  = esd->GetPrimaryVertex();
  TPolyMarker3D* pvm = make_vertex_marker(pv, "Primary Vertex");
  pvm->SetMarkerStyle(5);
  pvm->SetMarkerColor(5);
  register_vertex_marker(pvm);

  if(showSPD) {
    AliESDVertex*  spdv  = esd->GetVertex();
    TPolyMarker3D* spdvm = make_vertex_marker(spdv, "SPD Vertex");
    spdvm->SetMarkerStyle(2);
    spdvm->SetMarkerColor(6);
    register_vertex_marker(spdvm);
  }

  if(showBoxes)
    make_vertex_boxes(pv);
}
