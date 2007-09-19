#include <FMD/flow/AliFMDFlowBinned1D.h>
#include <FMD/flow/AliFMDFlowTrue.h>
#include <FMD/flow/AliFMDFlowUtil.h>
#include <TRandom.h>
#include <TCanvas.h>
#include <TVirtualPad.h>
#include <TArrayD.h>
#include <TBrowser.h>
#include <iostream>
#include <TMath.h>
#include <TH2.h>
#include <TStyle.h>

//____________________________________________________________________
struct Generator 
{
  Generator(Double_t psi=-1, 
	    Double_t v1=.05,  Double_t v2=.05, 
	    UInt_t   min=100, UInt_t   max=1000) 
    : fPsi(psi), fV1(v1), fV2(v2), fMin(min), fMax(max) 
  {}
  UInt_t Prepare() 
  {
    // Generate a uniform random direction 
    fN = 0;
    if (fPsi >= 0) fPsiR = fPsi;
    else           fPsiR = gRandom->Uniform(0, 2 * TMath::Pi());
    return unsigned(gRandom->Uniform(fMin, fMax));
  }
  Double_t operator()() 
  {
    // Generate a uniform random direction 
    Double_t phi  =  gRandom->Uniform(0, 2 * TMath::Pi());
    Double_t rel  =  phi - fPsiR;
    Double_t dphi =  -2 * TMath::Sin(rel) * fV1;
    dphi          -= TMath::Sin(2 * rel) * fV2;
    phi           += dphi;
    return phi;
  }
  Double_t Psi() const { return fPsiR; }
  Double_t fPsi;
  Double_t fPsiR;
  Double_t fV1;
  Double_t fV2;
  Double_t fMin;
  Double_t fMax;
  UInt_t   fN;
};

void
TestFlow(UInt_t n_events=100, Int_t seg=-1, UInt_t n_max=20000)
{
  Generator           generator(-1, 0.00, 0.05, n_max, n_max);
  AliFMDFlowAxis      a(10, -5, 5);
  AliFMDFlowBinned1D* flow = new AliFMDFlowBinned1D(2, a);
  AliFMDFlowTrue1D*   real = new AliFMDFlowTrue1D(2, a);
  TH2D*               hist = new TH2D("hist","hist",a.N(),a.Bins(),
				      (seg>0?seg:90),0, 2* TMath::Pi());
  Double_t            dphi = (seg <= 0 ? 0 : 2 * TMath::Pi() / seg);
  TArrayD             phis;
  TArrayD             tphis;
  TArrayD             xs;
  
  std::cout << std::endl;
  for (UInt_t i = 0; i < n_events; i++) {
    std::cout << "\rEvent # " << i << std::flush;
    // Generate an event, get the number of objects, get the true
    // event plane, shuffle the phis around randomly. 
    UInt_t   n_obs  = generator.Prepare();
    Double_t rpsi_r = NormalizeAngle(generator.Psi());
    phis.Set(n_obs);
    phis.Reset(0);
    tphis.Set(n_obs);
    tphis.Reset(0);
    xs.Set(n_obs);
    xs.Reset(0);
    
    real->SetPsi(rpsi_r);
    std::cout << " (" << n_obs << " observations)" << std::flush;
    for (UInt_t j = 0; j < n_obs; j++) {
      if (j % 2000 == 0) std::cout << "." << std::flush;
      Double_t x   = gRandom->Gaus(0, 3);
      Double_t phi = generator();
      tphis[j]     = phi;
      if (seg >= 0) { 
	Int_t iseg = Int_t(phi / dphi);
	phi        = dphi * (iseg + .5);
      }
      hist->Fill(x, phi);
      phis[j]      = phi;
      xs[j]        = x;
    }
    flow->Event(phis.fArray,  xs.fArray, 0, phis.fN);
    real->Event(tphis.fArray, xs.fArray, 0, tphis.fN);
  }
  std::cout << std::endl;
  flow->Print("s");
  real->Print("s");

  gStyle->SetPalette(1);
  
  TCanvas* c = new TCanvas("C");
  c->SetFillColor(0);
  c->Divide(2,1);
  TVirtualPad* p = c->cd(1);
  p->Divide(1,2);
  p = p->cd(1);
  p->SetFillColor(0);
  flow->Draw("bnst colz");
  p = c->cd(1);
  p = p->cd(2);
  p->SetFillColor(0);
  flow->Draw("bnstr colz");
  c->cd(2);
  hist->Draw("colz");
  new TBrowser("b", flow);
}

#ifndef __CINT__
#include <TApplication.h>

int
main()
{
  TApplication app("app", 0, 0);
  TestFlow();
  app.Run();
}

#endif
  

