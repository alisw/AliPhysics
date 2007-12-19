#include <TH1D.h>
#include <TH2D.h>
#include <TBrowser.h>
#include <TArrayF.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TROOT.h>
#include "AliFMDAnaESD.h"
#include "AliFMDAnaFlowRing.h"
#include <TFile.h>
#if 1
# include "background.C"
// # include "background2.C"
#else
TH2* bgFMD1i() { return 0; }
TH2* bgFMD2i() { return 0; }
TH2* bgFMD2o() { return 0; }
TH2* bgFMD3o() { return 0; }
TH2* bgFMD3o() { return 0; }
#endif
#if 1
# include "reactionplane.C"
#else
Float_t reactionplane(int ev) { return 0; }
#endif
#include <iostream>


//====================================================================
class AliFMDAnaFlow : public AliFMDAnaESD
{
public:
  //__________________________________________________________________
  /** Constructor */
  AliFMDAnaFlow(Int_t n=2, Bool_t bg=false, Float_t cut0=.68) 
    : fAxis(n, -3.7, 5.3),
      fV2("all", "Full FMD", 2, 1, fAxis)
  {
    for (int i = 0; i < 5; i++) { 
      UShort_t det = (i == 0 ? 1 : (i <= 2 ? 2 : 3));
      Char_t   rng = (i == 0 ? 'I' :  (i % 2 == 1 ? 'I' : 'O'));
      AddRing(new AliFMDAnaFlowRing(det,rng,n,EtaMin(i)+.1,EtaMax(i)-.1,
			   Bg(bg,i),cut0,Cut1(i)));
    }
  }
  //__________________________________________________________________
  /** Constructor */
  AliFMDAnaFlow(Float_t deta, Bool_t bg=false, Float_t cut0=.68) 
    : fAxis(Int_t((5.3- -3.7)/deta), -3.7, 5.3),
      fV2("all", "Full FMD", 2, 1, fAxis)
  {
    for (int i = 0; i < 5; i++) { 
      UShort_t det = (i == 0 ? 1 : i <= 2 ? 2 : 3);
      Char_t   rng = (i % 2 == 0 ? 'I' : 'O');
      AddRing(new AliFMDAnaFlowRing(det,rng,fAxis,Bg(bg,i),cut0,Cut1(i)));
    }
  }
  //__________________________________________________________________
  /** Initialize */
  Bool_t Init()
  {
    fPhis.Set(51200);
    fEtas.Set(51200);
    fWeights.Set(51200);
    fNObserv = 0;

    Float_t pi2 = 2 * TMath::Pi();

    fPsi  = new TH2D("all_psi", "Event planes", 
		     fAxis.N(), fAxis.Bins(), 40, 0, pi2);
    fPsi->SetXTitle("#eta");
    fPsi->SetYTitle("#Psi");
    fPsi->SetZTitle("Frequency");
    fPsi->SetMarkerStyle(20);

    fDPsi = new TProfile("all_dpsi", "Distance to real #Psi_{R}", 
			 fAxis.N(), fAxis.Bins(), -pi2, pi2);
    fDPsi->SetXTitle("#eta");
    fDPsi->SetYTitle("#Delta#Psi=|#Psi-#Psi_{R}|");

    fDPhi = new TProfile2D("all_dphi","Distribution of phi signals", 
			  fAxis.N(), -3.7, 5.3, 40, 0, pi2, 0, 20);
    fDPhi->SetXTitle("#eta");
    fDPhi->SetYTitle("#varphi");
    fDPhi->SetZTitle("signal");
    fDPhi->SetMarkerStyle(20);

    return AliFMDAnaESD::Init();
  }
  //__________________________________________________________________
  /** Called at the beginning of an event */
  Bool_t Begin(Int_t ev) 
  {
    fNObserv = 0;
    Float_t psir = reactionplane(fNEvents);
    AliFMDAnaFlowRing* fring = 0;
    for (size_t i = 0; i < 5; i++) {
      if (!fRing[i]) continue;
      fring = static_cast<AliFMDAnaFlowRing*>(fRing[i]);
      fring->fTrue.SetPsi(psir);
    }
    return AliFMDAnaESD::Begin(ev);
  }
  //__________________________________________________________________
  /** Fill in one strip */
  void Fill(Float_t phi, Float_t eta, Float_t mult)
  {
    Int_t imult = Int_t(mult+.33);
    if (fNObserv + imult >= fEtas.fN) { 
      ULong_t n = Int_t(fEtas.fN * 1.5)+imult;
      fEtas.Set(n);
      fPhis.Set(n);
      fWeights.Set(n);
    }
    for (Int_t i = 0; i < imult; i++) { 
      fEtas[fNObserv]    = eta;
      fPhis[fNObserv]    = phi;
      fWeights[fNObserv] = 1;
      fNObserv++;
    }
    fDPhi->Fill(eta, phi);
  }
  //__________________________________________________________________
  /** Called at the end of an event */
  Bool_t End() 
  {
    float       psir    = reactionplane(fNEvents-1);
    // std::cout << "PsiR=" << psir << std::endl;
    fV2.Event(fNObserv, fPhis.fArray, fEtas.fArray, 
	      fWeights.fArray, fWeights.fArray);
    for (UShort_t i = 0; i < fAxis.N(); i++) { 
      AliFMDFlowBin* bin = fV2.GetBin(i);
      if (bin && bin->Counts() > 0) {
	fPsi->Fill(fAxis.BinCenter(i),  bin->Psi() - psir);
	fDPsi->Fill(fAxis.BinCenter(i), bin->Psi() - psir);
      }
    }
    return AliFMDAnaESD::End();
  }
  //__________________________________________________________________
  /** Make a frame */
  TH1* MakeFrame(const char* xname, const char* yname, 
		 Float_t ymin, Float_t ymax)
  {
    static int id = 0;
    id++;
    TH1* frame = new TH1D(Form("frame%d", id), "Frame", 
			  fAxis.N(), fAxis.Bins()[0], fAxis.Bins()[fAxis.N()]);
    frame->SetMinimum(ymin);
    frame->SetMaximum(ymax);
    frame->SetXTitle(xname);
    frame->SetYTitle(yname);
    frame->SetLabelFont(132, "XYZ");
    frame->SetTitleFont(132, "XYZ");
    frame->SetLabelSize(.08, "XYZ");
    frame->SetTitleSize(.08, "XYZ");
    frame->SetTitleOffset(.5, "XYZ");
    frame->SetNdivisions(8, "XYZ");
    frame->Draw();
    return frame;
  }
  //__________________________________________________________________
  /** called at the end of a run */
  Bool_t Finish()
  {
    gStyle->SetPalette(1);
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    AliFMDAnaFlowRing* fring[5];
    for (size_t i = 0; i < 5; i++) 
      fring[i] = static_cast<AliFMDAnaFlowRing*>(fRing[i]);
    
    TCanvas* c = new TCanvas("c1", "Result", 600, 750);
    c->cd();
    c->SetGridy();
    c->SetTopMargin(0.02);
    c->SetRightMargin(0.02);
    c->SetFillColor(0);
    c->SetBorderSize(0);
    c->SetBorderMode(0);
    c->SetHighLightColor(0);
    c->Divide(1,3,0,0);

    TVirtualPad* p = 0;
    TH1*         f = 0;

    p = c->cd(1);
    p->SetGridx();
    p->SetGridy();
    p->SetFillColor(0);
    f = MakeFrame("#eta", "|#Psi_{2}-#Psi_{R}|", -TMath::Pi(), TMath::Pi());
    for (int i = 0; i < 5; i++) fring[i]->fDPsi.Draw("same");
    fDPsi->Draw("same");

    p = c->cd(2);
    p->SetGridx();
    p->SetGridy();
    p->SetFillColor(0);
    f = MakeFrame("#eta", "R_{k}", .7, 1);
    for (int i = 0; i < 5; i++) fring[i]->fV2.Draw("tr:same");
    fV2.Draw("tr:same");

    p = c->cd(3);
    p->SetGridx();
    p->SetGridy();
    p->SetFillColor(0);
    f = MakeFrame("#eta", "v_{2}", 0, .1);
    for (int i = 0; i < 5; i++) { 
      fring[i]->fV2.Draw("th:same");
      fring[i]->fTrue.Draw("th:same");
    }
    fV2.Draw("th:same");


    for (int i = 0; i < 5; i++) { 
      fring[i]->fV2.Print("s");
      fring[i]->fTrue.Print("s");
    }
    fV2.Print("s");
    gROOT->GetListOfSpecials()->Add(this);
    return kTRUE;
  }
  void ToFile(const char* file) { 
    TFile* out = TFile::Open(file, "RECREATE");
    for (int i = 0; i < 5; i++) fRing[i]->Write(fRing[i]->Name());
    fV2.Write();
    if (fPsi) fPsi->Write();
    if (fDPsi) fDPsi->Write();
    if (fDPhi) fDPhi->Write();
    out->Close();
  }
    
    

  //__________________________________________________________________
  /** Get min pseudo rapidity of ring idx */
  Float_t EtaMin(Int_t idx) const { 
    switch (idx) { 
    case 0: return  3.61728; break; // FMD1I
    case 1: return  2.18122; break; // FMD2I
    case 2: return  1.79821; break; // FMD2O
    case 3: return -3.39913; break; // FMD3I
    case 4: return -2.28925; break; // FMD3O
    }
    return 0;
  }
  //__________________________________________________________________
  /** Get max pseudo rapidity of ring idx */
  Float_t EtaMax(Int_t idx) const { 
    switch (idx) { 
    case 0: return  5.02643; break; // FMD1I
    case 1: return  3.57899; break; // FMD2I
    case 2: return  2.39085; break; // FMD2O
    case 3: return -2.00644; break; // FMD3I
    case 4: return -1.70080; break; // FMD3O
    }
    return 0;
  }
  //__________________________________________________________________
  /** Get cut 1 of ring idx */
  Float_t Cut1(Int_t idx) const { 
    switch (idx) { 
    case 0: return 1.33857; break; // FMD1I
    case 1: return 1.30464; break; // FMD2I
    case 2: return 1.27263; break; // FMD2O
    case 3: return 1.28052; break; // FMD3I
    case 4: return 1.23041; break; // FMD3O
    }
    return 0;
  }
  //__________________________________________________________________
  /** Get background correction of ring idx */
  TH2* Bg(Bool_t use, Int_t idx) const { 
    if (!use) return 0;
    switch (idx) { 
    case 0: return bgFMD1i(); break;
    case 1: return bgFMD2i(); break;
    case 2: return bgFMD2o(); break;
    case 3: return bgFMD3i(); break;
    case 4: return bgFMD3o(); break;
    }
    return 0;
  }
  /** Flow axis */
  AliFMDFlowAxis fAxis;
  /** Flow object */
  AliFMDFlowBinned1D fV2;
  /** Cache of phi */
  TArrayD            fPhis;
  /** Cache of eta */
  TArrayD            fEtas;
  /** Cache of weights */
  TArrayD            fWeights;
  /** Number of observables */
  Int_t              fNObserv;
  /** Histogram of Psi */
  TH2D*              fPsi;       // Event planes
  /** Histogram of dPsi */
  TProfile*          fDPsi;   // Mean distance to Psi
  /** Histogram of dPhi */
  TProfile2D*        fDPhi;
  ClassDef(AliFMDAnaFlow,1) 
};

//____________________________________________________________________
//
// EOF
//
