// $Header$

TPolyMarker3D* make_vertex_marker(AliESDVertex* v, const Text_t* name )
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
  v->GetXYZ(x);
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

void primary_vertex_primitive(Bool_t showSPD=kTRUE, Bool_t showBoxes=kFALSE)
{
  AliESDEvent* esd = Alieve::Event::AssertESD();

  AliESDVertex*  pv  = esd->GetPrimaryVertex();
  TPolyMarker3D* pvm = make_vertex_marker(pv, "Primary Vertex");
  pvm->SetMarkerStyle(5);
  pvm->SetMarkerColor(5);
  pvm->SetMarkerSize(1.4);
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

/**************************************************************************/

Reve::StraightLineSet* ESDvertex_lineset(AliESDVertex* v, const Text_t* name)
{ 
  using namespace Reve;

  Double_t x[3], e[3];
  v->GetXYZ(x); v->GetSigmaXYZ(e);
  printf("%16s: %f %f %f   -   %f %f %f\n", name,
	 x[0], x[1], x[2], e[0], e[1], e[2]);

  // dimensions
  Reve::StraightLineSet* ls = new Reve::StraightLineSet();
  ls->SetName(name);
  ls->AddLine(e[0], 0, 0, -e[0], 0, 0); 
  ls->AddLine(0, e[1], 0, 0, -e[1], 0);
  ls->AddLine(0, 0, e[2], 0, 0, -e[2]);
  for(Int_t i =0; i < 3; i++)
  {
    ls->AddMarker(i, 0);
    ls->AddMarker(i, 1);
  }

  // centre marker
  ls->AddMarker(0, 0.5);
  ls->RefHMTrans().SetPos(x);
  return ls;  
}

void make_vertex_ellipses(Reve::StraightLineSet* ls, AliESDVertex* v, Bool_t ellipseUseSigma)
{
  using namespace Reve;

  Double_t x[3], e[3];
  v->GetXYZ(x); v->GetSigmaXYZ(e);

  if(ellipseUseSigma)
  {
    e[0] *= 30; e[1] *= 30; e[2] *= 10;
    ls->SetMarkerStyle(5);
    ls->SetMarkerColor(5);
    ls->SetMarkerSize(1.4);
    ls->SetLineColor(7);
    ls->SetTitle("+- 30 sigma_r x 10 sigma_z");
  }
  else
  {
    e[0] = 1; e[1] = 1; e[2] = 2;
    ls->SetMarkerStyle(2);
    ls->SetMarkerColor(6);
    ls->SetLineColor(6);
    ls->SetTitle("+- 10 x 10 x 20mm");
  }
  Int_t N = 32;
  Float_t S = 2*TMath::Pi()/N;
  Float_t b, a, phi;

  a = e[0]; b = e[1];
  for(Int_t i = 0; i<N; i++)
    ls->AddLine(a*TMath::Cos(i*S)  , b*TMath::Sin(i*S)  , 0, 
		a*TMath::Cos(i*S+S), b*TMath::Sin(i*S+S), 0);

  a = e[0]; b = e[2];
  for(Int_t i = 0; i<N; i++)
    ls->AddLine(a*TMath::Cos(i*S)  , 0, b*TMath::Sin(i*S), 
		a*TMath::Cos(i*S+S), 0, b*TMath::Sin(i*S+S));

  a = e[1]; b = e[2];
  for(Int_t i = 0; i<N; i++)
    ls->AddLine(0, a*TMath::Cos(i*S)  ,  b*TMath::Sin(i*S), 
		0, a*TMath::Cos(i*S+S),  b*TMath::Sin(i*S+S));
}

void primary_vertex(Bool_t showSPD=kTRUE, Bool_t rnrEllipse=kTRUE)
{ 
  AliESDEvent* esd = Alieve::Event::AssertESD();
  Reve::StraightLineSet* ls;

  AliESDVertex* PV  =  esd->GetPrimaryVertex();
  ls = ESDvertex_lineset(PV, "Primary Vertex");
  if(rnrEllipse) make_vertex_ellipses(ls, PV, kTRUE);
  gReve->AddRenderElement(ls);

  if(showSPD) 
  {
    AliESDVertex*  SPDV  = esd->GetVertex();
    ls = ESDvertex_lineset(SPDV, "SPD Vertex");
    if(rnrEllipse) make_vertex_ellipses(ls, SPDV, kFALSE);
    gReve->AddRenderElement(ls);
  }
}
