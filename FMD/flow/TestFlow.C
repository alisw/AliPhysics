/* Copyright (C) 2007 Christian Holm Christensen <cholm@nbi.dk>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1 of
 * the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 * USA
 */
//___________________________________________________________________
/** @file 
    @brief Example 
    
    Compile and run like 
    @verbatim 
    Root> gSystem->Load("libFMDflow.so");
    Root> gSystem->AddIncludePath("-I$ALICE_ROOT/FMD")
    Root> .x FMD/flow/TestFlow.C 
    @endverbatim 

    or 
    @verbatim 
    $ g++ `root-config --cflags --libs` -I$ALICE_ROOT/include \
       -I$ALICE_ROOT/FMD $ALICE_ROOT/FMD/flow/TestFlow.C -o testflow
    $ ./testflow 
    @endverbatim 
*/
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
#include <TFile.h>

//____________________________________________________________________
/** Generate events */
struct Generator 
{
  /** Constructor 
      @param psi Possibly fixed event plane (if <0 -> random)
      @param v1  value of v1 
      @param v2  value of v2 
      @param min Minimum number of observations
      @param max Maximum number of observations */
  Generator(Double_t psi=-1, 
	    Double_t v1=.05,  Double_t v2=.05, 
	    UInt_t   min=100, UInt_t   max=1000) 
    : fPsi(psi), fV1(v1), fV2(v2), fMin(min), fMax(max) 
  {}
  /** Prepare for the next event. 
      @return Number of observations */
  UInt_t Prepare() 
  {
    // Generate a uniform random direction 
    fN = 0;
    if (fPsi >= 0) fPsiR = fPsi;
    else           fPsiR = gRandom->Uniform(0, 2 * TMath::Pi());
    return unsigned(gRandom->Uniform(fMin, fMax));
  }
  /** Create an observation */
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
  /** Get the event plane */
  Double_t Psi() const { return fPsiR; }
  Double_t fPsi;
  Double_t fPsiR;
  Double_t fV1;
  Double_t fV2;
  Double_t fMin;
  Double_t fMax;
  UInt_t   fN;
};

/** Run test program */
void
TestFlow(UInt_t n_events=100, Int_t seg=-1, UInt_t n_max=20000,
	 Float_t v2=0.05)
{
  Generator           generator(-1, 0.00, v2, n_max, n_max);
  AliFMDFlowAxis      a(10, -5, 5);
  AliFMDFlowBinned1D  flow("test","Test",2, 1, a);
  AliFMDFlowTrue1D    real("true","true", 2, a);
  TH2D                hist("hist","hist",a.N(),a.Bins(),
			   (seg>0?seg:90),0, 2* TMath::Pi());
  TH1D                rela("rela","real",
			   (seg>0?seg:90),0, 2* TMath::Pi());
  Double_t            dphi = (seg <= 0 ? 0 : 2 * TMath::Pi() / seg);
  TArrayD             phis(n_max);
  TArrayD             tphis(n_max);
  TArrayD             xs(n_max);
  flow.SetLineColor(kRed+2);
  real.SetLineColor(kBlue+2);
  rela.SetXTitle("#Psi_{R}");
  rela.SetYTitle("#Psi_{2}");

  std::cout << std::endl;
  for (UInt_t i = 0; i < n_events; i++) {
    std::cout << "\rEvent # " << i << std::flush;
    // Generate an event, get the number of objects, get the true
    // event plane, shuffle the phis around randomly. 
    UInt_t   n_obs  = generator.Prepare();
    Double_t rpsi_r = NormalizeAngle(generator.Psi());
    
    real.SetPsi(rpsi_r);
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
      phis[j]      = phi;
      xs[j]        = x;
      hist.Fill(x, NormalizeAngle(phi-rpsi_r));
      rela.Fill(NormalizeAngle(phi-rpsi_r));
    }
    flow.Event(n_obs, phis.fArray,  xs.fArray);
    real.Event(n_obs, tphis.fArray, xs.fArray);
  }
  std::cout << std::endl;
  flow.Print("s");
  real.Print("s");

  gStyle->SetPalette(1);
  
  TCanvas* c = new TCanvas("C");
  c->SetFillColor(0);
  c->Divide(2,2);
  TVirtualPad* p = c->cd(1);
  // p->Divide(1,2);
  // p = p->cd(1);
  p->SetFillColor(0);
  flow.Draw("th:");
  real.Draw("th:same");
  // p = c->cd(1);
  p = c->cd(2);
  rela.Scale(1. / n_events);
  rela.Draw("");
  p = c->cd(3);
  p->SetFillColor(0);
  flow.Draw("tr:");
  real.Draw("tr:same");
  c->cd(4);
  hist.Scale(1. / n_events);
  hist.Draw("lego2");
  TBrowser* b = new TBrowser("b");
  b->Add(&flow);
  b->Add(&real);

  TFile* file = TFile::Open("flow_test.root", "RECREATE");
  flow.Write();
  real.Write();
  file->Close();
  delete file;
}

#ifndef __CINT__
#include <TApplication.h>
/** Entry point for test program */
int
main()
{
  TApplication app("app", 0, 0);
  TestFlow();
  app.Run();
}

#endif
  

