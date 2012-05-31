void
TestESDPhi()
{
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliCDBManager::Instance()->SetRun(0);
  
  AliGeomManager::LoadGeometry();
  
  AliFMDGeometry* geom = AliFMDGeometry::Instance();
  geom->Init();
  geom->InitTransformations();

  AliESDFMD esd;
  
  for (UShort_t d = 1; d <= 3; d++) { 
    UShort_t nRng = (d == 1 ? 1 : 2);
    for (UShort_t q = 0; q < nRng; q++) { 
      Char_t   r    = (q == 0 ? 'I' : 'O');
      UShort_t nSec = (q == 0 ?  20 :  40);
      for (UShort_t s = 0; s < nSec; s++) {
	Double_t x, y, z;
	geom->Detector2XYZ(d, r, s, 0, x, y, z);
	Double_t a = TMath::ATan2(y, x);
	Double_t p = esd.Phi(d, r, s, 0);
	
	if (a < 0) a+= 2 * TMath::Pi();
	a *= 180 / TMath::Pi();

	Printf("FMD%d%c[%2d]: Geom: %5.1f, ESD: %5.1f", d, r, s, a, p);
      }
    }
  }
  for (UShort_t q = 0; q < nRng; q++) { 
    Char_t   r    = (q == 0 ? 'I' : 'O');
    UShort_t nStr = (q == 0 ? 512 : 256);
    for (UShort_t t = 0; t < nStr; t++) {
      Double_t x, y, z;
      geom->Detector2XYZ(2, r, 0, t, x, y, z);
      Double_t l = TMath::Sqrt(x * x + y * y);
      Double_t v = esd.R(2, r, 0, t);

      Printf("FMD%c[%3d]: Geom: %6.3f, ESD: %6.3f", r, t, l, v);
    }
  }
}
