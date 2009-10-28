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
       -I$ALICE_ROOT/FMD -L$ALICE_ROOT/lib/tgt_${ALICE_TARGET} \
       -lFMDflow $ALICE_ROOT/FMD/flow/TestFlow.C -o testflow
    $ ./testflow 
    @endverbatim 
*/
#include <flow/AliFMDFlowBinned1D.h>
#include <flow/AliFMDFlowTrue.h>
#include <flow/AliFMDFlowUtil.h>
#include <TRandom.h>
#include <TCanvas.h>
#include <TVirtualPad.h>
#include <TLegend.h>
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

struct Tester 
{
  Tester(Float_t v2=0.05, Int_t seg=-1, UInt_t n_max=20000)
    : fGenerator(-1, 0.00, v2, n_max, n_max),
      fAxis(10, -5, 5), 
      fFlow("flow", "Analysed", 2, 1, fAxis), 
      fReal("true", "Truth", 2, fAxis), 
      fHist("hist", "Histogram", fAxis.N(), fAxis.Bins(),
	    (seg > 0 ? seg : 90), 0, 2*TMath::Pi()),
      fRela("rela", "Truth", (seg > 0 ? seg : 90), 0, 2*TMath::Pi()),
      fDPhi(seg <= 0 ? 0 : 2 * TMath::Pi() / seg), 
      fPhis(n_max), 
      fTruePhis(n_max), 
      fXs(n_max)
  {
    fFlow.SetLineColor(kRed+2);
    fReal.SetLineColor(kBlue+2);
    fRela.SetXTitle("#Psi_{R}");
    fRela.SetYTitle("#Psi_{2}");
  }
  void Run(UInt_t n_events) 
  {
    for (UInt_t i = 0; i < n_events; i++) Event(i);
    std::cout << std::endl;
    PrintResult();
    DrawResult(n_events);
    StoreResult();
  }
  void Event(Int_t i)
  {
    std::cout << "\rEvent # " << i << std::flush;
    // Generate an event, get the number of objects, get the true
    // event plane, shuffle the phis around randomly. 
    UInt_t   n_obs  = fGenerator.Prepare();
    Double_t rpsi_r = NormalizeAngle(fGenerator.Psi());
    
    fReal.SetPsi(rpsi_r);
    std::cout << " (" << n_obs << " observations)" << std::flush;
    for (UInt_t j = 0; j < n_obs; j++) {
      Observation(j, rpsi_r);
    }
    fFlow.Event(n_obs, fPhis.fArray,     fXs.fArray);
    fReal.Event(n_obs, fTruePhis.fArray, fXs.fArray);
  }
  void Observation(Int_t j, Double_t rpsi_r)
  {
    if (j % 2000 == 0) std::cout << "." << std::flush;
    Double_t x   = gRandom->Gaus(0, 3);
    Double_t phi = fGenerator();
    fTruePhis[j] = phi;
    if (fDPhi != 0) { 
      Int_t iseg = Int_t(phi / fDPhi);
      phi        = fDPhi * (iseg + .5);
    }
    fPhis[j]     = phi;
    fXs[j]       = x;
    fHist.Fill(x, NormalizeAngle(phi-rpsi_r));
    fRela.Fill(NormalizeAngle(phi-rpsi_r));
  }
  void PrintResult()
  {
    fFlow.Print("s");
    fReal.Print("s");
  }
  void DrawResult(Int_t n_events, Bool_t showHists=kFALSE)
  {
    gStyle->SetPalette(1);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    
    TCanvas* c = new TCanvas("C", "TestFlow", 600, (showHists ? 800 : 600));
    c->SetFillColor(0);
    c->SetBorderSize(0);
    c->SetBorderMode(0);
    Double_t x2 = (showHists ? 0.5 : 1);
    TPad* pflow = new TPad("pflow", "PFlow", 0.0, 0.5, x2, 1.0, 0, 0, 0);
    TPad* pres  = new TPad("pres",  "PRes",  0.0, 0.0, x2, 0.5, 0, 0, 0);
    TPad* phist = 0;
    TPad* p2d   = 0;
    if (showHists) {
      phist = new TPad("phist", "PHist", x2, 0.5, 1.0, 1.0, 0, 0, 0);
      p2d   = new TPad("p2d",   "P2d",   x2, 0.0, 1.0, 0.5, 0, 0, 0);
    }
    pflow->SetTopMargin(0.05);
    pflow->SetBottomMargin(0.00);
    pflow->SetRightMargin(0.05);
    pres->SetTopMargin(0.00);
    pres->SetRightMargin(0.05);

    c->cd();
    pflow->Draw();
    pflow->cd();
    TH1* fth =fFlow.MakeHistogram(AliFMDFlowBin::kTdr,AliFMDFlowBin::kHarmonic);
    TH1* tth =fReal.MakeHistogram(AliFMDFlowBin::kTdr,AliFMDFlowBin::kHarmonic);
    Double_t min = TMath::Min(fth->GetMinimum(), tth->GetMinimum());
    Double_t max = TMath::Max(fth->GetMaximum(), tth->GetMaximum());
    fth->SetMinimum(0.9 * min); 
    tth->SetMinimum(0.9 * min); 
    fth->SetMaximum(1.2 * max); 
    tth->SetMaximum(1.2 * max); 
    fth->Draw();
    tth->Draw("same");
    TLegend* l = pflow->BuildLegend(0.5, 0.7, 0.94, 0.94);
    l->SetFillColor(0);
    l->SetBorderSize(0);

    c->cd();
    pres->Draw();
    pres->cd();
    TH1* ftr = fFlow.MakeHistogram(AliFMDFlowBin::kTdr,
				   AliFMDFlowBin::kResolution);
    TH1* ttr = fReal.MakeHistogram(AliFMDFlowBin::kTdr,
				   AliFMDFlowBin::kResolution);
    min      = TMath::Min(ftr->GetMinimum(), ttr->GetMinimum());
    max      = TMath::Max(ftr->GetMaximum(), ttr->GetMaximum());
    ftr->SetMinimum(0.8 * min); 
    ttr->SetMinimum(0.8 * min); 
    ftr->SetMaximum(1.2 * max); 
    ttr->SetMaximum(1.2 * max); 
    ftr->Draw();
    ttr->Draw("same");
    l = pres->BuildLegend(0.5, 0.7, 0.94, 0.94);
    l->SetFillColor(0);
    l->SetBorderSize(0);

    if (p2d) {
      p2d->Draw();
      p2d->cd();
      p2d->SetTopMargin(0.05);
      p2d->SetRightMargin(0.05);
      p2d->SetFillColor(0);
      fRela.Scale(1. / n_events);
      fRela.DrawCopy("");
    }

    if (phist) { 
      phist->Draw();
      phist->cd();
      phist->SetTopMargin(0.05);
      phist->SetRightMargin(0.05);
      fHist.Scale(1. / n_events);
      fHist.DrawCopy("lego2");
    }
    c->Modified();
    c->Update();
    c->cd();
  }
  void BrowseResult()
  {
    TBrowser* b = new TBrowser("b");
    b->Add(&fFlow);
    b->Add(&fReal);
  }
  void StoreResult() 
  {
    TFile* file = TFile::Open("flow_test.root", "RECREATE");
    fFlow.Write();
    fReal.Write();
    file->Close();
    delete file;
  }
  Generator          fGenerator;   // Event plane generator
  AliFMDFlowAxis     fAxis;        // Bins
  AliFMDFlowBinned1D fFlow;        // Our histogram
  AliFMDFlowTrue1D   fReal;        // Histogram of true values
  TH2D               fHist;        // Histogram of data
  TH1D               fRela;        // Histogram of resolution
  Double_t           fDPhi;        // Phi segmentation
  TArrayD            fPhis;        // Cache of phis
  TArrayD            fTruePhis;    // Cache of true phis
  TArrayD            fXs;          // Cache of X's
};
//____________________________________________________________________
void
TestFlow(UInt_t n_events=100, Int_t seg=-1, UInt_t n_max=20000,
	 Float_t v2=0.05)
{
  std::cout << "Flow:                   " << v2 << "\n"
	    << "# of events:            " << n_events << "\n" 
	    << "# of phi segments:      ";
  if (seg < 0) std::cout << "none\n"; 
  else         std::cout << seg << "\n";
  std::cout << "Maximum # observations: " << n_max << std::endl;
		 
  Tester t(v2, seg, n_max);
  t.Run(n_events);
}

#ifndef __CINT__
#include <sstream>
#include <TApplication.h>
/** Entry point for test program */
void 
usage(std::ostream& o, const char* argv0)
{
  o << "Usage: " << argv0 << " [OPTIONS]\n\n"
    << "Options:\n"
    << "   -n N_EVENT          Set the number of events\n"
    << "   -s N_SEGMENTS       Set number of phi segments\n" 
    << "   -o MAX_OBSERVATIONS Set maximum number of observations per event\n"
    << "   -v FLOW             Set the flow\n"
    << "   -h                  Show this help\n"
    << std::endl;
}

template <typename T>
void str2val(const char* str, T& val)
{
  std::stringstream s(str);
  s >> val;
}

int
main(int argc, char** argv)
{
  UInt_t  n_events = 100;
  Int_t   seg      = -1;
  UInt_t  n_max    = 20000;
  Float_t v2      = 0.05;
  for (int i = 1; i < argc; i++) { 
    if (argv[i][0] == '-') { 
      switch (argv[i][1]) { 
      case 'h': usage(std::cout, argv[0]); return 0;
      case 'n': str2val(argv[++i], n_events); break;
      case 'o': str2val(argv[++i], n_max); break;
      case 's': str2val(argv[++i], seg); break;
      case 'v': str2val(argv[++i], v2); break;
      default:
	std::cerr << "Unknown option: " << argv[i] << std::endl;
	return 1;
      }
    }
  }
  TApplication app("app", 0, 0);
  TestFlow(n_events, seg, n_max, v2);
  app.Run();
}

#endif
  

