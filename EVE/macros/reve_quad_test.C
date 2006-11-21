// $Header$

Reve::QuadSet* reve_quad_test(Float_t x=0, Float_t y=0, Float_t z=0)
{
  TRandom r(0);

  gStyle->SetPalette(1, 0);

  Reve::QuadSet* q = new Reve::QuadSet("Pepe");
  q->Reset(Reve::QuadSet::QT_AxisAligned, kFALSE, 32);
  for (Int_t i=0; i<10; ++i) {
    q->AddQuad(r.Uniform(-10, 10), r.Uniform(-10, 10), r.Uniform(-10, 10),
	       r.Uniform(-1, 1), r.Uniform(-1, 1));
    q->QuadValue(r.Uniform(0, 130));
  }
  q->RefitPlex();

  Reve::ZTrans& t = q->RefHMTrans();
  t.SetPos(x, y, z);

  gReve->AddRenderElement(q);
  gReve->Redraw3D();

  return q;
}

Reve::QuadSet* reve_quad_test_circ()
{
  TRandom r(0);

  gStyle->SetPalette(1, 0);

  Reve::QuadSet* q = new Reve::QuadSet("Pepe");
  q->Reset(Reve::QuadSet::QT_AxisAligned, kFALSE, 32);

  Float_t R = 10, dW = 1, dH = .5;
  for (Int_t i=0; i<12; ++i) {
    Float_t x = R * TMath::Cos(TMath::TwoPi()*i/12);
    Float_t y = R * TMath::Sin(TMath::TwoPi()*i/12);
    q->AddQuad(x-dW, y-dH, r.Uniform(-1, 1), 2*dW, 2*dH);
    q->QuadValue(r.Uniform(0, 130));
  }
  q->RefitPlex();

  Reve::ZTrans& t = q->RefHMTrans();
  t.SetPos(0, 0, 300);

  gReve->AddRenderElement(q);
  gReve->Redraw3D();

  return q;
}
