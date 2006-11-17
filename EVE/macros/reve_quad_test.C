// $Header$

void reve_quad_test()
{
  TRandom r(0);

  gStyle->SetPalette(1, 0);

  Reve::QuadSet* q = new Reve::QuadSet("Pepe");
  q->Reset(Reve::QuadSet::QT_AxisAligned, kFALSE, 32);
  for (Int_t i=0; i<128; ++i) {
    q->AddQuad(r.Uniform(-10, 10), r.Uniform(-10, 10), r.Uniform(-10, 10),
	       r.Uniform(-1, 1), r.Uniform(-1, 1));
    q->QuadValue(r.Uniform(0, 130));
  }
  q->RefitPlex();

  gReve->AddRenderElement(q);
  gReve->Redraw3D();
}
