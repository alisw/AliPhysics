class TCanvas;

TH1* Make(UShort_t d, Char_t r)
{
  TH1* h = new TH1D(Form("fmd%d%c", d, r),
		    Form("FMD%d%c", d, r),
		    200, -4, 6);
  Color_t col = ((d == 1 ? kRed+2 :
		  d == 2 ? kGreen+(r == 'I'?2:-2) :
		  kBlue+(r == 'I'?2:-2)));
  // (r == 'I' || r == 'i') ? -2 : +2);
  h->SetDirectory(0);
  h->SetXTitle("#eta");
  h->SetYTitle("");
  h->SetFillColor(col);
  h->SetLineColor(col);
  h->SetMarkerColor(col);
  h->SetMarkerStyle(19+d);
  h->Sumw2();
  return h;
}

void Event(const TVector3& ip, TList* list, Bool_t fake)
{
  for (UShort_t d = 1; d <= 3+(fake?1:0); d++) {
    // Printf("FMD%d", d);
    Int_t nQ = (d == 1  ? 1 : 2);
    for (UShort_t q=0; q<nQ; q++) {
      char r = (q == 0 ? 'i' : 'o');
      UShort_t nSec = (q == 0 ?  20 :  40);
      UShort_t nStr = (q == 0 ? 512 : 256);
      UShort_t idx  = (d-1)*2-(1-q)+(d==1?1:0);
      TH1* h = static_cast<TH1*>(list->At(idx));
      if (!h) {
	h = Make(d,r);
	Info("", "Adding %s @ %d", h->GetName(), idx);
	list->AddAt(h, idx);
      }
      // Info("", "Got histogram %p", h);
      // printf(" FMD%d%c ", d, r);
      for (UShort_t s = 0; s < nSec; s++) {
	// printf(".");
	for (UShort_t t = 0; t < nStr; t++) {
	  // UShort_t m = (t % 4);
	  // Char_t   c = (m == 0 ? '/' : m == 1 ? '-' : m == 2 ? '\\' : '|');
	  // printf("%c\b", c);
	  // Double_t x, y, z;
	  // fmd->Detector2XYZ(d, r, s, t, x, y, z);

	  // TVector3 pos;
	  // AliForwardUtil::GetXYZ(d, r, s, t, ip, pos);
	  // tolerance set to 1/2 mm for X,Y
	  // Double_t tolXY = 5e-2;
	  // Compare(d, r, s, t, "X", x, pos.X(), 2, tolXY);
	  // Compare(d, r, s, t, "Y", y, pos.Y(), 2, tolXY);
	  // Compare(d, r, s, t, "Z", z, pos.Z(), 2);
	  
	  // Double_t gRadius, gEta, gPhi, gTheta;
	  // fmd->XYZ2REtaPhiTheta(x,y,z,gRadius,gEta,gPhi,gTheta);
	  // if (gPhi < 0) gPhi += TMath::TwoPi();

	  Double_t uEta, uPhi;
	  AliForwardUtil::GetEtaPhi(d, r, s, t, ip, uEta, uPhi);
	  h->Fill(uEta);
	  // Tolerance for eta set to 0.01
	  // Double_t tolEta = 0.01;
	  // Compare(d,r,s,t,"eta",gEta, uEta, 0, tolEta);
	  // Compare(d,r,s,t,"phi",gPhi, uPhi, 1);
	} // for t
      } // for s
      // Printf("");
    }
  }
}
		  
void TestEtaPhi()
{

  TList list;
  Int_t nEv = 500;
  Double_t mean = .5;
  Double_t var  = 4;
  Int_t iEv = 0;
  
  do {
    Double_t ipZ = gRandom->Gaus(mean, var);
    if (TMath::Abs(ipZ) > 10) continue;
    if (iEv != 0 && (iEv % 100 == 0)) Printf("Event # %4d", iEv);
    // Info("", "Event %d ipZ=%f", iEv, ipZ);
    Double_t ipX = gRandom->Gaus(.1, .003);
    Double_t ipY = gRandom->Gaus(.33, .003);
    TVector3 ip(ipX, ipY, ipZ);
    Event(ip, &list, false);
    iEv++;
  } while (iEv < nEv);

  THStack* stack = new THStack("stack","");
  TIter next(&list);
  TH1*  hist = 0;
  while ((hist = static_cast<TH1*>(next())))
    stack->Add(hist);
  stack->Draw("nostack p");
    
}
