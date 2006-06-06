// $Header$

void region_marker(Float_t a=10, Float_t b=10, Float_t c=20,
		   Float_t x=0, Float_t y=0, Float_t z=0)
{
  TPolyMarker3D* mark = new TPolyMarker3D(8);
  mark->SetName("Origin marker");
  mark->SetMarkerStyle(3);
  mark->SetMarkerColor(6);
  mark->SetPoint(0, x+a, y+b, z+c);
  mark->SetPoint(1, x+a, y-b, z+c);
  mark->SetPoint(2, x-a, y-b, z+c);
  mark->SetPoint(3, x-a, y+b, z+c);

  mark->SetPoint(4, x+a, y+b, z-c);
  mark->SetPoint(5, x+a, y-b, z-c);
  mark->SetPoint(6, x-a, y+b, z-c);
  mark->SetPoint(7, x-a, y-b, z-c);

  using namespace Reve;
  Color_t* colp = FindColorVar(mark, "fMarkerColor");
  RenderElementObjPtr* rnrEl = new RenderElementObjPtr(mark, *colp);
  gReve->AddGlobalRenderElement(rnrEl);
  gReve->DrawRenderElement(rnrEl);
}
