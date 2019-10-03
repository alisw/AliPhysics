void
Compare(UShort_t d, Char_t r, UShort_t s, UShort_t t,
	const char* what, 
	Double_t fromGeom, Double_t fromUtil, UShort_t type=0,
	Double_t eps=1e-4)
{
  if (TMath::Abs(fromGeom-fromUtil) < eps) return;
  Double_t c = 1;
  if (type==1) c = TMath::RadToDeg();
  Double_t g = c * fromGeom;
  Double_t u = c * fromUtil;
  Printf("  FMD%d%c[%2d,%3d]: %20s G:%10.4f U:%10.4f -> D:%10.4f > %g %s",
	 d, r, s, t, what, g, u, (g-u), eps,
	 type==0 ? "" : type==1 ? "[degrees]" : "[cm]" );
}
void
TestGeomUtils(Bool_t init=false, ULong_t run=138190)
{
  if (!gROOT->GetGlobal("gGeoManager")) init = true;
  AliCDBManager* cdb  = AliCDBManager::Instance();
  if (init) {
    cdb->SetDefaultStorageFromRun(run);
    cdb->SetRun(run);
    AliGeomManager::LoadGeometry();
    AliGeomManager::ApplyAlignObjsFromCDB("FMD");
  }
  AliFMDGeometry* fmd = AliFMDGeometry::Instance();
  if (init) {    
    fmd->Init();
    fmd->InitTransformations();
  }

  TVector3 ip(0,0,0);
  for (UShort_t d = 1; d <= 3; d++) {
    Printf("FMD%d", d);
    Int_t nQ = (d == 1  ? 1 : 2);
    for (UShort_t q=0; q<nQ; q++) {
      char r = (q == 0 ? 'I' : 'O');
      UShort_t nSec = (q == 0 ?  20 :  40);
      UShort_t nStr = (q == 0 ? 512 : 256);
      printf(" FMD%d%c ", d, r);
      for (UShort_t s = 0; s < nSec; s++) {
	printf(".");
	for (UShort_t t = 0; t < nStr; t++) {
	  UShort_t m = (t % 4);
	  Char_t   c = (m == 0 ? '/' : m == 1 ? '-' : m == 2 ? '\\' : '|');
	  printf("%c\b", c);
	  Double_t x, y, z;
	  fmd->Detector2XYZ(d, r, s, t, x, y, z);

	  TVector3 pos;
	  AliForwardUtil::GetXYZ(d, r, s, t, ip, pos);
	  // tolerance set to 1/2 mm for X,Y
	  Double_t tolXY = 5e-2;
	  Compare(d, r, s, t, "X", x, pos.X(), 2, tolXY);
	  Compare(d, r, s, t, "Y", y, pos.Y(), 2, tolXY);
	  Compare(d, r, s, t, "Z", z, pos.Z(), 2);
	  
	  Double_t gRadius, gEta, gPhi, gTheta;
	  fmd->XYZ2REtaPhiTheta(x,y,z,gRadius,gEta,gPhi,gTheta);
	  // if (gPhi < 0) gPhi += TMath::TwoPi();

	  Double_t uEta, uPhi;
	  AliForwardUtil::GetEtaPhi(d, r, s, t, ip, uEta, uPhi);
	  // Tolerance for eta set to 0.01
	  Double_t tolEta = 0.01;
	  Compare(d,r,s,t,"eta",gEta, uEta, 0, tolEta);
	  Compare(d,r,s,t,"phi",gPhi, uPhi, 1);
	} // for t
      } // for s
      Printf("");
    }
  }
}
