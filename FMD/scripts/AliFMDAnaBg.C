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
#include <AliFMDInput.h>
#include <AliFMDAnaRing.h>
#include <AliFMDGeometry.h>
#include <AliCDBManager.h>
#include <iostream>

//====================================================================
class AliFMDAnaBgRing : public AliFMDAnaRing
{
public:
  AliFMDAnaBgRing(UShort_t det, Char_t ring, 
		  Float_t deta, Float_t etaMin, Float_t etaMax)
    : AliFMDAnaRing(det, ring, 0, 0, 0),
      fTotal(Form("total_%s", fName), "Total number of hits", 
	     Int_t((etaMax-etaMin)/deta+.5), etaMin, etaMax, 
	     fNSeq, 0, 2*TMath::Pi()),
      fPrimary(Form("primary_%s", fName), "Number of hits from primary", 
	       Int_t((etaMax-etaMin)/deta+.5), etaMin, etaMax, 
	       fNSeq, 0, 2*TMath::Pi()),
      fRatio(Form("ratio_%s", fName), "Correction factors", 
	       Int_t((etaMax-etaMin)/deta+.5), etaMin, etaMax, 
	     fNSeq, 0, 2*TMath::Pi()),
      fInverse(Form("inverse_%s", fName), "Inverse Correction factors", 
	       Int_t((etaMax-etaMin)/deta+.5), etaMin, etaMax, 
	       fNSeq, 0, 2*TMath::Pi())
  {
    fTotal.SetDirectory(0);
    fTotal.Sumw2();
    fTotal.SetXTitle("#eta");
    fTotal.SetYTitle("#varphi");
    fTotal.SetZTitle("# of hits");
    fPrimary.SetDirectory(0);
    fPrimary.Sumw2();
    fPrimary.SetXTitle("#eta");
    fPrimary.SetYTitle("#varphi");
    fPrimary.SetZTitle("# of primary");
    fRatio.SetDirectory(0);
    fRatio.Sumw2();
    fRatio.SetXTitle("#eta");
    fRatio.SetYTitle("#varphi");
    fRatio.SetZTitle("# of primary/# of hits");
    fInverse.SetDirectory(0);
    fInverse.Sumw2();
    fInverse.SetXTitle("#eta");
    fInverse.SetYTitle("#varphi");
    fInverse.SetZTitle("# of hits/# of primary");
  }
  void Fill(Float_t phi, Float_t eta, Float_t mult)
  {
    fTotal.Fill(eta, phi);
    if (mult >= 1) fPrimary.Fill(eta,phi);
  }
  void Finish()
  {
    fRatio.Add(&fPrimary);
    fRatio.Divide(&fTotal);
    fInverse.Add(&fTotal);
    fInverse.Divide(&fPrimary);
  }
  void Browse(TBrowser* b) 
  {
    // AliFMDAnaRing::Browse(b);
    b->Add(&fTotal);
    b->Add(&fPrimary);
    b->Add(&fRatio);
    b->Add(&fInverse);
  }
  TH2D fTotal;
  TH2D fPrimary;
  TH2D fRatio;
  TH2D fInverse;
};

//====================================================================
class AliFMDAnaBg : public AliFMDInput
{
public:
  //__________________________________________________________________
  /** Constructor */
  AliFMDAnaBg(Float_t deta=.1) 
    : fEtaMin(-3.7), 
      fEtaMax(5.3),
      fTotal("total", "Total number of hits", 
	     Int_t((fEtaMax-fEtaMin)/deta+.5), fEtaMin, fEtaMax, 
	     40, 0, 2*TMath::Pi()),
      fPrimary("primary", "Number of hits from primary", 
	       Int_t((fEtaMax-fEtaMin)/deta+.5), fEtaMin, fEtaMax, 
	       40, 0, 2*TMath::Pi()),
      fRatio("ratio", "Correction factors", 
	     Int_t((fEtaMax-fEtaMin)/deta+.5), fEtaMin, fEtaMax, 
	     40, 0, 2*TMath::Pi())
  {
    fTotal.SetDirectory(0);
    fTotal.Sumw2();
    fTotal.SetXTitle("#eta");
    fTotal.SetYTitle("#varphi");
    fTotal.SetZTitle("# of hits");
    fPrimary.SetDirectory(0);
    fPrimary.Sumw2();
    fPrimary.SetXTitle("#eta");
    fPrimary.SetYTitle("#varphi");
    fPrimary.SetZTitle("# of primary");
    fRatio.SetDirectory(0);
    fRatio.Sumw2();
    fRatio.SetXTitle("#eta");
    fRatio.SetYTitle("#varphi");
    fRatio.SetZTitle("# of primary/# of hits");

    AliCDBManager::Instance()->SetRun(0);

    AddLoad(kHits);
    AddLoad(kKinematics);
    AddLoad(kGeometry);
    for (int i = 0; i < 5; i++) { 
      UShort_t det = (i == 0 ? 1 : (i <= 2 ? 2 : 3));
      Char_t   rng = (i == 0 ? 'I' :  (i % 2 == 1 ? 'I' : 'O'));
      fRing[i] = new AliFMDAnaBgRing(det,rng,deta,EtaMin(i)-.2,EtaMax(i)+.12);
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
    Bool_t ret = AliFMDInput::Init();
    for (size_t i = 0; i < 5; i++) fRing[i]->Init();
    AliFMDGeometry::Instance()->Init();
    AliFMDGeometry::Instance()->InitTransformations();
    return ret;
  }
  //__________________________________________________________________
  /** Called at the beginning of an event */
  Bool_t Begin(Int_t ev) 
  {
    for (size_t i = 0; i < 5; i++) fRing[i]->Begin();
    return AliFMDInput::Begin(ev);
  }
  //__________________________________________________________________
  /** Called at the beginning of an event */
  Bool_t ProcessHit(AliFMDHit* hit, TParticle* p)
  {
    Float_t          mult = p->IsPrimary() ? 1 : 0;
    UShort_t         det  = hit->Detector();
    Char_t           rng  = hit->Ring();
    UShort_t         sec  = hit->Sector();
    UShort_t         str  = hit->Strip();
    AliFMDAnaBgRing* frn  = fRing[FindRing(det, rng)];
    Double_t x, y, z;
    AliFMDGeometry::Instance()->Detector2XYZ(det,rng,sec,str,x,y,z);
    Double_t         phi  = TMath::ATan2(y,x);
    Double_t         r    = TMath::Sqrt(x * x + y * y);
    Double_t         the  = TMath::ATan2(r, z);
    Double_t         eta  = -TMath::Log(TMath::Tan(the/2));
    if (phi < 0) phi += 2*TMath::Pi();
    frn->Fill(phi, eta, mult);
    Fill(phi, eta, mult);
    return kTRUE;
  }

  //__________________________________________________________________
  /** Fill in one strip */
  void Fill(Float_t phi, Float_t eta, Float_t mult)
  {
    fTotal.Fill(eta, phi);
    if (mult >= 1) fPrimary.Fill(eta,phi);
  }
  //__________________________________________________________________
  /** Called at the end of an event */
  Bool_t End() 
  {
    for (UShort_t i = 0; i < 5; i++) fRing[i]->End();
    return AliFMDInput::End();
  }
  //__________________________________________________________________
  /** called at the end of a run */
  Bool_t Finish()
  {
    gStyle->SetPalette(1);
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

    fRatio.Add(&fPrimary);
    fRatio.Divide(&fTotal);

    for (UShort_t i = 0; i < 5; i++) fRing[i]->Finish();

    TCanvas* c = new TCanvas("c1", "Result", 800, 600);
    c->cd();
    c->SetGridy();
    c->SetTopMargin(0.02);
    c->SetRightMargin(0.02);
    c->SetFillColor(0);
    c->SetBorderSize(0);
    c->SetBorderMode(0);
    c->SetHighLightColor(0);
    c->Divide(5,1,0,0);
    for (UShort_t i = 0; i < 5; i++) { 
      c->cd(i+1);
      fRing[i]->fRatio.Draw("colz");
    }

    c = new TCanvas("c2", "All", 800, 600);
    c->cd();
    c->SetGridy();
    c->SetTopMargin(0.02);
    c->SetRightMargin(0.02);
    c->SetFillColor(0);
    c->SetBorderSize(0);
    c->SetBorderMode(0);
    c->SetHighLightColor(0);
    c->Divide(5,1,0,0);
    fRatio.Draw("colz");

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
  /** Flow rings */
  AliFMDAnaBgRing* fRing[5];
  Float_t fEtaMin;
  Float_t fEtaMax;
  TH2D fTotal;
  TH2D fPrimary;
  TH2D fRatio;
  ClassDef(AliFMDAnaBg,1) 
};

//____________________________________________________________________
//
// EOF
//
