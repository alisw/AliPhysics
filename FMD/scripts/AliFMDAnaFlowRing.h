// -*- mode: C++ -*-
#ifndef ALIFMDANAFLOWRING_H
#define ALIFMDANAFLOWRING_H
#include <flow/AliFMDFlowBinned1D.h>
#include <flow/AliFMDFlowTrue.h>
#include <flow/AliFMDFlowBin.h>
#include <flow/AliFMDFlowAxis.h>
#include <flow/AliFMDFlowUtil.h>
#include <AliFMDAnaRing.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TMath.h>
#include <TBrowser.h>

//====================================================================
class AliFMDAnaFlowRing : public AliFMDAnaRing
{
public:
  /** Constructor 
      @param d      Detector number 
      @param r      Ring identifier 
      @param n      Number of bins
      @param etamin Minimum @f$\eta@f$
      @param etamax Maximum @f$\eta@f$
      @param bg     Background 
      @param c1     Lower cut 
      @param c2     higher cut */
  AliFMDAnaFlowRing(UShort_t d=0, Char_t r='\0', 
		    Int_t n=0, Float_t etamin=0, Float_t etamax=0,
		    TH2* bg=0, Float_t c1=0, Float_t c2=0)
    : AliFMDAnaRing(d, r, bg, c1, c2), 
      fAxis(n, etamin, etamax),
      fV2(Name(), Name(), 2, 1, fAxis),
      fTrue(Form("True%s", Name()), Name(), 2, fAxis),
      fNObserv(0), 
      fPsi(Form("%s_psi", fName), "Event planes", 
	   fAxis.N(), fAxis.Bins(), 2*fNSeq, -TMath::Pi(), TMath::Pi()),
      fDPsi(Form("%s_dpsi", fName), "Distance to real #Psi_{R}", 
	    fAxis.N(), fAxis.Bins(), -TMath::Pi(), TMath::Pi()),
      fMult(Form("%s_mult", fName), "Multiplicity per hit", 
	    fAxis.N(), fAxis.Bins(), 100, 0, 10)
  {}
  /** Constructor 
      @param d     Detector number 
      @param r     Ring identifier 
      @param a     Axis 
      @param bg    Background 
      @param c1    Lower cut 
      @param c2    higher cut */
  AliFMDAnaFlowRing(UShort_t d, Char_t r, const AliFMDFlowAxis& a,
		    TH2* bg, Float_t c1, Float_t c2)
    : AliFMDAnaRing(d, r, bg, c1, c2), 
      fAxis(a),
      fV2(Name(), Name(), 2, 1, fAxis),
      fTrue(Form("True%s", Name()), Name(), 2, fAxis),
      fNObserv(0), 
      fPsi(Form("%s_psi", fName), "Event planes", 
	   fAxis.N(), fAxis.Bins(), 2*fNSeq, -TMath::Pi(), TMath::Pi()),
      fDPsi(Form("%s_dpsi", fName), "Distance to real #Psi_{R}", 
	    fAxis.N(), fAxis.Bins(), -TMath::Pi(), TMath::Pi()),
      fMult(Form("%s_mult", fName), "Multiplicity per hit", 
	    fAxis.N(), fAxis.Bins(), 100, 0, 10)
  {}
  /** Initialize */
  void Init()
  {
    fPhis.Set(10240);
    fEtas.Set(10240);
    fWeights.Set(10240);
    fNObserv = 0;

    fV2.SetLineColor(Color());
    fV2.SetMarkerColor(Color());
    fV2.SetFillColor(Color());
    fV2.SetFillStyle(3003+(fRing == 'I' ? 1 : 2));

    fPsi.SetXTitle("#eta");
    fPsi.SetYTitle("#Psi");
    fPsi.SetZTitle("Frequency");

    fDPsi.SetXTitle("#eta");
    fDPsi.SetYTitle("#Delta#Psi=|#Psi-#Psi_{R}|");

    fMult.SetXTitle("#eta");
    fMult.SetYTitle("M");

    ModHist(fPsi);
    ModHist(fDPsi);
    ModHist(fMult);

    AliFMDAnaRing::Init();
  }
  /** Modify histogram 
      @param h Histogram */
  void ModHist(TH1& h)
  {
    h.GetXaxis()->SetTitleFont(132);
    h.GetXaxis()->SetLabelFont(132);
    h.GetXaxis()->SetNdivisions(8);
    h.GetYaxis()->SetTitleFont(132);
    h.GetYaxis()->SetLabelFont(132);
    h.GetYaxis()->SetNdivisions(8);
    h.GetZaxis()->SetTitleFont(132);
    h.GetZaxis()->SetLabelFont(132);
    h.GetZaxis()->SetNdivisions(8);
    h.SetLineColor(Color());
    h.SetMarkerColor(Color());
    h.SetFillColor(Color());
    h.SetMarkerStyle(20);
    h.SetDirectory(0);
  }
  /** Called at beginning of event */
  void Begin() { fNObserv = 0; }
  /** Fill in one strip 
      @param phi   Azimuthal angle @f$ \varphi@f$ 
      @param eta   Pseudo rapidity @f$ \eta@f$
      @param mult  Signal */
  void Fill(Float_t phi, Float_t eta, Float_t mult)
  {
    // if (mult < 0.2 || mult > 2) return;
    // if (fBg && mult < fCut0) return;
    if (mult <= 0.01) return;
    Int_t imult = Int_t(mult+.1);
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
      fMult.Fill(eta, mult);
      fNObserv++;
    }
  }
  /** Called at end of event */
  void End() 
  { 
    fTrue.Event(fNObserv, fPhis.fArray, fEtas.fArray, 
		fWeights.fArray, 0); // fWeights.fArray);
    fV2.Event(fNObserv, fPhis.fArray, fEtas.fArray, 
	      fWeights.fArray, 0); // fWeights.fArray);
    for (UShort_t i = 0; i < fAxis.N(); i++) { 
      AliFMDFlowBin* bin  = fV2.GetBin(i);
      AliFMDFlowBin* tbin = fTrue.GetBin(i);
      if (bin && bin->Counts() > 0) {
	Double_t psir = tbin->Psi() / fV2.PsiOrder();
	fPsi.Fill(fAxis.BinCenter(i), bin->Psi() - psir);
	fDPsi.Fill(fAxis.BinCenter(i), bin->Psi() - psir);
      }
    }
  }  
  /** Browse this object 
      @param b Browser */
  void Browse(TBrowser* b)
  {
    AliFMDAnaRing::Browse(b);
    b->Add(&fV2,   "V2");
    b->Add(&fTrue, "TrueV2");
    b->Add(&fPsi);
    b->Add(&fDPsi);
    b->Add(&fMult);
  }
protected:
  /** Flow axis */
  AliFMDFlowAxis     fAxis;
  /** Flow histogram */
  AliFMDFlowBinned1D fV2;
  /** Flow histogram */
  AliFMDFlowTrue1D   fTrue;
  /** Cache of phi */
  TArrayD            fPhis;
  /** Cache of eta */
  TArrayD            fEtas;
  /** Cache of weights */
  TArrayD            fWeights;
  /** Number of observables */
  Int_t              fNObserv;
  /** Histogram of Psi */
  TH2D               fPsi;       // Event planes
  /** Histogram of dPsi */
  TProfile           fDPsi;   // Mean distance to Psi
  /** Histogram of multiplicities */
  TH2D               fMult;

  ClassDef(AliFMDAnaFlowRing,1) // Flow ring
};

#endif
//
// EOF
//
