#include <TF1.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TMath.h>
#include <TH2.h>
#include <FMD/flow/AliFMDFlowBessel.h>

//____________________________________________________________________
double myI(double* xp, double* pp)
{
  double n = pp[0];
  double x = xp[0];
  return AliFMDFlowBessel::I(n, x);
}

//____________________________________________________________________
double myDI(double* xp, double* pp)
{
  double n = pp[0];
  double x = xp[0];
  return AliFMDFlowBessel::DiffI(n, x);
}

//____________________________________________________________________
double _rootI(double n, double x)
{
  if (fabs(x) < 0.00001) return 0;
  if (n == 0.5) 
    return TMath::Sqrt(2 * x / TMath::Pi()) * TMath::SinH(x)/x;
  if (n == 1.5) 
    return TMath::Sqrt(2 * x / TMath::Pi()) * 
      (TMath::CosH(x) / x - TMath::SinH(x) / (x * x));
  if (n == 2.5) return _rootI(0.5, x) - 3 * _rootI(1.5, x) / x;
  return TMath::BesselI(int(fabs(n)), x);
}

//____________________________________________________________________
double rootI(double* xp, double* pp)
{
  return _rootI(pp[0], xp[0]);
}

//____________________________________________________________________
const char* ntonm(double n)
{
  if (int(2 * fabs(n)) % 2 == 1) 
    return Form("I%d_2", int(2 * n));
  return Form("I%d", int(n));
}

//____________________________________________________________________
const char* ntostr(double n)
{
  if (int(2 * fabs(n)) % 2 == 1) 
    return Form("I_{%d/2}", int(2 * n));
  return Form("I_{%d}", int(n));
}

//____________________________________________________________________
TCanvas* c = 0;
TCanvas* d = 0;

//____________________________________________________________________
void compare(double n, int col)
{
  Double_t min = 0.00001;
  Double_t max = 4;
  if (!c) {
    c = new TCanvas("c_I", "I_{#nu}");
    c->SetFillColor(0);
    c->cd();
    TH2* f = new TH2F("f_I", "I_{#nu}", 100, min, max, 100, -1, 13);
    f->SetXTitle("x");
    f->SetYTitle("I_{#nu}"); 
    f->SetStats(0);
    f->Draw();
  }
  TString  str = ntostr(n);
  TString  nm  = ntonm(n);

  c->cd();
  TF1* r = new TF1(Form("root_%s", nm.Data()), rootI, min, max, 1);
  r->SetTitle(Form("%s (ROOT)", str.Data()));
  r->SetParameter(0, n);
  r->SetLineColor(col);
  r->Draw("same");

  c->cd();
  TF1* m = new TF1(Form("my_%s", nm.Data()), myI, min, max, 1);
  m->SetTitle(Form("%s (My)", str.Data()));
  m->SetParameter(0, n);
  m->SetLineColor(col);
  m->SetLineStyle(2);
  m->Draw("same");

  c->Modified();
  c->Update();
  c->cd();

}
//____________________________________________________________________
void dcompare(double n, int col)
{
  Double_t min = 0;
  Double_t max = 4;
  if (!d) { 
    d = new TCanvas("c_dI", "dI_{#nu}/dx");
    d->SetFillColor(0);
    d->cd();
    TH2* f = new TH2F("f_dI", "dI_{#nu}/dx", 100, min, max, 100, -1, 12);
    f->SetXTitle("x");
    f->SetYTitle("dI_{#nu}/dx");
    f->SetStats(0);
    f->Draw();
  }

  TString  str = ntostr(n);
  TString  nm  = ntonm(n);

  d->cd();
  TF1* r = new TF1(Form("root_d%s", nm.Data()), rootI, min, max, 1);
  r->SetTitle(Form("d%s/dx (ROOT)", str.Data()));
  r->SetParameter(0, n);
  r->SetLineColor(col);
  TGraph* g = new TGraph(r, "d");
  g->SetTitle(r->GetTitle());
  g->Draw("same c");

  d->cd();
  TF1* m = new TF1(Form("my_d%s", nm.Data()), myDI, min, max, 1);
  m->SetTitle(Form("d%s/dx (My)", str.Data()));
  m->SetParameter(0, n);
  m->SetLineColor(col);
  m->SetLineStyle(2);
  m->Draw("same");

  d->Modified();
  d->Update();
  d->cd();
}
  
//____________________________________________________________________
void
TestBessel()
{
  gStyle->SetPalette(1);
  // gStyle->SetOptTitle(0);
  // double n[] = { 0, 0.5, 1, 1.5, 2, 2.5, 3 };
  double n[] = { -1, 0, 0.5, 1, 1.5, 2, 2.5, 3 };
  // double n[] = { 0.5, 1.5, 2.5 };
  // double n[] = { 0.5, 1.5, 2.5 };
  size_t  m  = sizeof(n)/sizeof(double);
  Int_t   nc = gStyle->GetNumberOfColors();
  for (size_t i = 0; i < m; i++) { 
    int col = gStyle->GetColorPalette(int(float(i) / m * nc));
    compare(n[i], col);
    dcompare(n[i], col);
  }
  c->cd();
  c->SetTopMargin(0.05);
  c->SetRightMargin(0.05);
  TLegend* l = c->BuildLegend(.1, .35, .4, .95);
  l->SetFillColor(0);
  l->SetBorderSize(1);
  // c->Print("I.ps");
  d->cd();
  d->SetTopMargin(0.05);
  d->SetRightMargin(0.05);
  l = d->BuildLegend(.1, .35, .4, .95);
  l->SetFillColor(0);
  l->SetBorderSize(1);
  // d->Print("dI.ps");
}

#ifndef __CINT__
int 
main()
{
  TApplication app("app", 0, 0);
  testBessel();
  app.Run();
  return 0;
}
#endif
