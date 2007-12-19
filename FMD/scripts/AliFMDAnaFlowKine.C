#include <AliFMDHit.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TParticle.h>
#include <TBrowser.h>
#include <TArrayF.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TROOT.h>
#include <AliFMDAnaESD.h>
#include "AliFMDAnaFlowRing.h"
#if 1
# include "background.C"
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
class AliFMDAnaFlowKine : public AliFMDInput
{
public:
  //__________________________________________________________________
  /** Constructor */
  AliFMDAnaFlowKine(Int_t n=1, Bool_t segments=true, 
		    Bool_t primary_only=true, Float_t cut0=.68) 
    : fAxis(n, -3.7, 5.3),
      fV2("all", "Full FMD", 2, fAxis), 
      fPrimaryOnly(primary_only),
      fSegmented(segments)
  {
    AddLoad(kHits);
    AddLoad(kKinematics);
    for (int i = 0; i < 5; i++) { 
      UShort_t det = (i == 0 ? 1 : (i <= 2 ? 2 : 3));
      Char_t   rng = (i == 0 ? 'I' :  (i % 2 == 1 ? 'I' : 'O'));
      fRing[i] = new AliFMDAnaFlowRing(det,rng,n,EtaMin(i)-.1,EtaMax(i)+.1,
				       0,cut0,Cut1(i));
    }
  }
  //__________________________________________________________________
  /** Constructor */
  AliFMDAnaFlowKine(Float_t deta, Bool_t segments=true, 
		    Bool_t primary_only=true, Float_t cut0=.68) 
    : fAxis(Int_t((5.3- -3.7)/deta), -3.7, 5.3),
      fV2("all", "Full FMD", 2, fAxis),
      fPrimaryOnly(primary_only),
      fSegmented(segments)
  {
    for (int i = 0; i < 5; i++) { 
      UShort_t det = (i == 0 ? 1 : i <= 2 ? 2 : 3);
      Char_t   rng = (i % 2 == 0 ? 'I' : 'O');
      fRing[i] = new AliFMDAnaFlowRing(det,rng,fAxis,0,cut0,Cut1(i));
    }
  }
  //__________________________________________________________________
  Int_t FindRing(UShort_t det, Char_t ring) const
  { 
    Int_t idx = -1;
    switch (det) { 
    case 1: idx = 0; break;
    case 2: idx = 1 + (ring == 'O' || ring == 'o' ? 1 : 0); break;
    case 3: idx = 3 + (ring == 'O' || ring == 'o' ? 1 : 0); break;
    }
    return idx;
  }  
  //__________________________________________________________________
  /** Initialize */
  Bool_t Init()
  {
    fPhis.Set(51200);
    fEtas.Set(51200);
    fWeights.Set(51200);
    fNObserv = 0;
    fNEvents = 0;

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

    for (size_t i = 0; i < 5; i++) fRing[i]->Init();
    return AliFMDInput::Init();
  }
  //__________________________________________________________________
  /** Called at the beginning of an event */
  Bool_t Begin(Int_t ev) 
  {
    fNObserv = 0;
    fNEvents++;
    Float_t psir  = reactionplane(fNEvents-1);
    for (size_t i = 0; i < 5; i++) { 
      fRing[i]->Begin();
      fRing[i]->fTrue.SetPsi(psir);
    }
    return AliFMDInput::Begin(ev);
  }
  //__________________________________________________________________
  /** Called at the beginning of an event */
  Bool_t ProcessHit(AliFMDHit* hit, TParticle* p)
  {
    if (fPrimaryOnly && !p->IsPrimary()) return kTRUE;
    if (hit->IsStop())   return kTRUE;
    UShort_t           det  = hit->Detector();
    Char_t             rng  = hit->Ring();
    UShort_t           sec  = hit->Sector();
    AliFMDAnaFlowRing* frn  = fRing[FindRing(det, rng)];
    Double_t           the  = p->Theta();
    Double_t           phi  = p->Phi();
    Double_t           eta  = -TMath::Log(TMath::Tan(the/2));
    // Phi 
    if (fSegmented) { 
      Double_t         dPhi =  2*TMath::Pi() / frn->NSeq();
      if (det == 3)    phi  =  TMath::Pi() - (sec + .5) * dPhi;
      else     	       phi  =  (sec + .5) * dPhi;
      if (phi < 0)     phi  += 2 * TMath::Pi();
    }
    
    Double_t           mult = 1; // hit->Edep();
    frn->Fill(phi, eta, mult);
    Fill(phi, eta, mult);
    return kTRUE;
  }

  //__________________________________________________________________
  /** Fill in one strip */
  void Fill(Float_t phi, Float_t eta, Float_t mult)
  {
    if (fNObserv >= fEtas.fN) { 
      ULong_t n = Int_t(fEtas.fN * 1.5)+1;
      fEtas.Set(n);
      fPhis.Set(n);
      fWeights.Set(n);
    }
    fEtas[fNObserv]    = eta;
    fPhis[fNObserv]    = phi;
    fWeights[fNObserv] = mult;
    fNObserv++;
    fDPhi->Fill(eta, phi);
  }
  //__________________________________________________________________
  /** Called at the end of an event */
  Bool_t End() 
  {
    for (size_t i = 0; i < 5; i++) fRing[i]->End();
    float       psir    = reactionplane(fNEvents-1)/2;
    // std::cout << "PsiR=" << psir << std::endl;
    fV2.Event(fPhis.fArray, fEtas.fArray, fWeights.fArray, fNObserv);
    for (UShort_t i = 0; i < fAxis.N(); i++) { 
      AliFMDFlowBin* bin = fV2.GetBin(i);
      if (bin && bin->Counts() > 0) {
	fPsi->Fill(fAxis.BinCenter(i),  bin->Psi() - psir);
	fDPsi->Fill(fAxis.BinCenter(i), bin->Psi() - psir);
      }
    }
    return AliFMDInput::End();
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
    f = MakeFrame("#eta", "v_{2}", 0, .1);
    for (int i = 0; i < 5; i++) { 
      fring[i]->fV2.Draw("th:same");
      fring[i]->fTrue.Draw("th:same");
    }
    fV2.Draw("th:same");

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
    f = MakeFrame("#eta", "Counts", 0, 2e5);
    f->GetYaxis()->SetNoExponent();
    for (int i = 0; i < 5; i++) fring[i]->fV2.Draw("tc:same");
    // fV2.Draw("tc:same");

    c = new TCanvas("c2", "Event plane");
    c->cd();
    c->SetGridy();
    c->SetTopMargin(0.02);
    c->SetRightMargin(0.02);
    c->SetFillColor(0);
    c->SetBorderSize(0);
    c->SetBorderMode(0);
    c->SetHighLightColor(0);
    float pi2 = 2*TMath::Pi();
    TH2* h = new TH2D("h", "h", fAxis.N(), fAxis.Bins(), 80, -pi2, pi2);
    h->Draw();
    for (int i = 0; i < 5; i++) fring[i]->fDPsi.Draw("same col");
    // fDPhi->Draw("same");

    for (int i = 0; i < 5; i++) { 
      fring[i]->fV2.Print("s");
      fring[i]->fTrue.Print("s");
    }
    fV2.Print("s");
    // gROOT->GetListOfSpecials()->Add(this);
    return kTRUE;
  }
  //__________________________________________________________________
  Bool_t IsFolder() const { return kTRUE; }
  //__________________________________________________________________
  void Browse(TBrowser* b) 
  {
    for (size_t i = 0; i < 5; i++) 
      b->Add(fRing[i], fRing[i]->Name());
    b->Add(&fV2);
    if (fPsi)  b->Add(fPsi);
    if (fDPsi) b->Add(fDPsi);
    if (fDPhi) b->Add(fDPhi);
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
  /** Number of events */ 
  Int_t fNEvents;
  /** Flow rings */
  AliFMDAnaFlowRing* fRing[5];
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
  /** Wether to include primaries only */ 
  Bool_t             fPrimaryOnly;
  /** Whether to segment the phi space according to detector */
  Bool_t             fSegmented;
  ClassDef(AliFMDAnaFlowKine,1) 
};

//____________________________________________________________________
//
// EOF
//
