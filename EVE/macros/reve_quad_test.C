// $Header$

void reve_quad_test()
{
  Alieve::ITSModule* qs = new Alieve::ITSModule("QuadSet Test");
  qs->Test(200);
  gReve->AddRenderElement(qs);
  gReve->DrawRenderElement(qs);
}
