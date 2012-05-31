void
res(Char_t i)
{
  Double_t pitch  = (i == 'I' ?  0.025 :  0.048);
  Double_t open   = (i == 'I' ? 18.0   :  9.0  ) / 180 * TMath::Pi();
  Double_t rmin   = (i == 'I' ?  4.3   : 15.6  );
  Double_t rmax   = (i == 'I' ? 17.2   : 28.0  );
  Double_t phimin = 0;
  Double_t phimax = 2 * TMath::Pi();
  TF2* xres = new TF2("xres", "sqrt(pow(cos(x),2)*[0]+y*y*pow(sin(x),2)*[1])",
		      phimin,phimax,rmin,rmax);
  TF2* yres = new TF2("yres", "sqrt(pow(sin(x),2)*[0]+y*y*pow(cos(x),2)*[1])",
		      phimin,phimax,rmin,rmax);
  xres->SetParameters(pitch*pitch, open*open); 
  yres->SetParameters(pitch*pitch, open*open);
  xres->GetHistogram()->SetXTitle("#phi [radians]");
  xres->GetHistogram()->SetYTitle("r [cm]");
  xres->GetHistogram()->SetZTitle("#delta x [cm]");
  yres->GetHistogram()->SetXTitle("#phi [radians]");
  yres->GetHistogram()->SetYTitle("r [cm]");
  yres->GetHistogram()->SetZTitle("#delta y [cm]");
  xres->SetLineColor(i == 'I' ? 2 : 6);
  yres->SetLineColor(i == 'I' ? 3 : 7);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  xres->Draw("surf");
  yres->Draw("same surf");

  TLegend* l = new TLegend(.7,.8,.95,.95, i == 'I' ? "Inner" : "Outer");
  l->SetBorderSize(0);
  l->AddEntry(xres, "#delta x", "l");
  l->AddEntry(yres, "#delta y", "l");
  l->Draw();
}

