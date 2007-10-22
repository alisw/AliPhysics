// $Header$

Reve::QuadSet* reve_quad_test(Float_t x=0, Float_t y=0, Float_t z=0,
			      Int_t num=100, Bool_t register=kTRUE)
{
  TRandom r(0);

  gStyle->SetPalette(1, 0);

  Reve::RGBAPalette* pal = new Reve::RGBAPalette(0, 130);

  Reve::QuadSet* q = new Reve::QuadSet("RectangleXY");
  q->SetPalette(pal);
  q->Reset(Reve::QuadSet::QT_RectangleXY, kFALSE, 32);
  for (Int_t i=0; i<num; ++i) {
    q->AddQuad(r.Uniform(-10, 10), r.Uniform(-10, 10), r.Uniform(-10, 10),
	       r.Uniform(0.2, 1), r.Uniform(0.2, 1));
    q->QuadValue(r.Uniform(0, 130));
  }
  q->RefitPlex();

  Reve::ZTrans& t = q->RefHMTrans();
  t.SetPos(x, y, z);

  if (register)
  {
    gReve->AddRenderElement(q);
    gReve->Redraw3D();
  }

  return q;
}

Reve::QuadSet* reve_quad_test_emc(Float_t x=0, Float_t y=0, Float_t z=0,
				  Int_t num=100)
{
  TRandom r(0);

  gStyle->SetPalette(1, 0);

  Reve::QuadSet* q = new Reve::QuadSet("EMC Supermodule");
  q->SetOwnIds(kTRUE);
  q->Reset(Reve::QuadSet::QT_RectangleXZFixedDimY, kFALSE, 32);
  q->SetDefWidth(8);
  q->SetDefHeight(8);

  for (Int_t i=0; i<num; ++i) {
    q->AddQuad(r.Uniform(-100, 100), r.Uniform(-100, 100));
    q->QuadValue(r.Uniform(0, 130));
    q->AddId(new TNamed(Form("Cell %d", i)));
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
  q->Reset(Reve::QuadSet::QT_RectangleXY, kFALSE, 32);

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

Reve::QuadSet* reve_quad_test_hex(Float_t x=0, Float_t y=0, Float_t z=0,
				  Int_t num=100, Bool_t register=kTRUE)
{
  TRandom r(0);

  gStyle->SetPalette(1, 0);

  {
    Reve::QuadSet* q = new Reve::QuadSet("HexagonXY");
    q->Reset(Reve::QuadSet::QT_HexagonXY, kFALSE, 32);
    for (Int_t i=0; i<num; ++i) {
      q->AddHexagon(r.Uniform(-10, 10), r.Uniform(-10, 10), r.Uniform(-10, 10),
		    r.Uniform(0.2, 1));
      q->QuadValue(r.Uniform(0, 120));
    }
    q->RefitPlex();

    Reve::ZTrans& t = q->RefHMTrans();
    t.SetPos(x, y, z);

    if (register)
    {
      gReve->AddRenderElement(q);
      gReve->Redraw3D();
    }
  }

  {
    Reve::QuadSet* q = new Reve::QuadSet("HexagonYX");
    q->Reset(Reve::QuadSet::QT_HexagonYX, kFALSE, 32);
    for (Int_t i=0; i<num; ++i) {
      q->AddHexagon(r.Uniform(-10, 10), r.Uniform(-10, 10), r.Uniform(-10, 10),
		    r.Uniform(0.2, 1));
      q->QuadValue(r.Uniform(0, 120));
    }
    q->RefitPlex();

    Reve::ZTrans& t = q->RefHMTrans();
    t.SetPos(x, y, z);

    if (register)
    {
      gReve->AddRenderElement(q);
      gReve->Redraw3D();
    }
  }

  return q;
}

Reve::QuadSet* reve_quad_test_hexid(Float_t x=0, Float_t y=0, Float_t z=0,
				    Int_t num=100, Bool_t register=kTRUE)
{
  TRandom r(0);

  gStyle->SetPalette(1, 0);

  {
    Reve::QuadSet* q = new Reve::QuadSet("HexagonXY");
    q->SetOwnIds(kTRUE);
    q->Reset(Reve::QuadSet::QT_HexagonXY, kFALSE, 32);
    for (Int_t i=0; i<num; ++i) {
      q->AddHexagon(r.Uniform(-10, 10), r.Uniform(-10, 10), r.Uniform(-10, 10),
		    r.Uniform(0.2, 1));
      q->QuadValue(r.Uniform(0, 120));
      q->QuadId(new TNamed(Form("Quad with idx=%d", i), "This title is not confusing."));
    }
    q->RefitPlex();

    Reve::ZTrans& t = q->RefHMTrans();
    t.SetPos(x, y, z);

    if (register)
    {
      gReve->AddRenderElement(q);
      gReve->Redraw3D();
    }
  }

  return q;
}

void reve_quad_test_hierarchy(Int_t n=4)
{
  gStyle->SetPalette(1, 0);

  Reve::RGBAPalette* pal = new Reve::RGBAPalette(20, 100);
  pal->SetLimits(0, 120);

  Reve::FrameBox*    box = new Reve::FrameBox();
  box->SetAABox(-10, -10, -10, 20, 20, 20);
  box->SetFrameColor((Color_t) 33);

  Reve::RenderElementList* l = new Reve::RenderElementList("Parent/Dir");
  l->SetTitle("Tooltip");
  //  l->SetMainColor((Color_t)3);
  gReve->AddRenderElement(l);

  // PMD: obtain digit-tree from run-loader, loop over entries.

  for (Int_t i=0; i<n; ++i)
  {
    Reve::QuadSet* qs = reve_quad_test_hexid(0, 0, 50*i, 50, kFALSE);
    // PMD: loop over clones-array, create hexagons above threshold.
    qs->SetPalette(pal);
    qs->SetFrame(box);
    gReve->AddRenderElement(qs, l);
  }

  gReve->Redraw3D();
}
