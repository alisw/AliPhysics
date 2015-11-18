#ifndef __CINT__
# include <TH1.h>
# include <TH2.h>
# include <TMath.h>
# include <TParticle.h>
# include <TObjArray.h>
# include <TString.h>
#else
class TH1;
class TParticle;
#endif
#include "FastAnalysis.C"

struct spectraAnalysis : public FastAnalysis
{
  enum {
    kPiSty = 20,
    kKSty  = 21,
    kPSty  = 22,
    kOSty  = 23
  };
  enum {
    kPiColor = kRed+2,
    kKColor  = kGreen+2,
    kPColor  = kBlue+2,
    kOColor  = kMagenta+2
  };
  spectraAnalysis(Bool_t verbose=false, Int_t monitor=0)
    : FastAnalysis(verbose,monitor),
      fNpi(0),
      fNK(0),
      fNP(0),
      fNO(0),      
      fdNdptPi(0),
      fdNdptK(0),
      fdNdptP(0),
      fdNdptO(0),
      fK2piRatio(0),
      fP2piRatio(0),
      fO2piRatio(0)      
  {
  }
  TH1* CreateSpectra(UInt_t pdg, Color_t col, Style_t sty)
  {
    TH1* h = new TH1D(dNdptName(pdg), "",   100, 0, 10);
    h->SetXTitle("p_{T}");
    h->SetYTitle("dN/dp_{T}");
    h->SetDirectory(0);
    h->SetMarkerColor(col);
    h->SetLineColor(col);
    h->SetFillColor(col);
    h->SetMarkerStyle(sty);
    h->Sumw2();
    fOutput->Add(h);
    return h;
  }
  TH1* CreateRatio(UInt_t pdg, Color_t col, Style_t sty)
  {
    TH1* h = new TH1D(RatioName(pdg), "", 100, 0, 10);
    h->SetXTitle(Form("%d/211",pdg));
    h->SetYTitle("Events");
    h->SetDirectory(0);
    h->SetMarkerColor(col);
    h->SetLineColor(col);
    h->SetFillColor(col);
    h->SetMarkerStyle(sty);
    h->Sumw2();
    fOutput->Add(h);
    return h;
  }
  /** 
   * Called on each slave at start of processing
   */
  virtual void SlaveBegin(TTree*)
  {
    Info("SlaveBegin", "Making dN/deta histogram");
    fdNdptPi = CreateSpectra(211,   kPiColor,kPiSty);
    fdNdptK  = CreateSpectra(321,   kKColor, kKSty);
    fdNdptP  = CreateSpectra(2212,  kPColor, kPSty);
    fdNdptO  = CreateSpectra(0,     kOColor, kOSty);

    fK2piRatio = CreateRatio(321,  kKColor, kKSty);
    fP2piRatio = CreateRatio(2212, kPColor, kPSty);
    fO2piRatio = CreateRatio(0,    kOColor, kOSty);
  }
  Bool_t ProcessHeader() { return true; }
  void ProcessParticles()
  {
    fNpi = fNK = fNP = fNO = 0;
    FastAnalysis::ProcessParticles();
    if (fNpi <= 0) return;
    fK2piRatio->Fill(Double_t(fNK)/fNpi);
    fP2piRatio->Fill(Double_t(fNP)/fNpi);
    fO2piRatio->Fill(Double_t(fNO)/fNpi);
  }
  /** 
   * Process a particle.  
   * 
   * @param p Particle to process
   */    
  virtual Bool_t ProcessParticle(const TParticle* p)
  {
    Double_t pT    = p->Pt();
    Int_t    pdg = TMath::Abs(p->GetPdgCode());
    if      (pdg == 211)     { fdNdptPi->Fill(pT); fNpi++; }
    else if (pdg == 321)     { fdNdptK ->Fill(pT); fNK++; }
    else if (pdg == 2212)    { fdNdptP ->Fill(pT); fNP++; }
    else                     { fdNdptO ->Fill(pT); fNO++; }
    return true;
  }
  static const char* dNdptName(UInt_t pdg) 
  {
    switch (pdg) {
    case 211:  return "dNdptPi";
    case 321:  return "dNdptK";
    case 2212: return "dNdptP";
    }
    return "dNdptO";
  }
  static const char* RatioName(UInt_t pdg) 
  {
    switch (pdg) {
    case 321:  return "KtoPiRatio";
    case 2212: return "PtoPiRatio";
    }
    return "OtoPiRatio";
  }
  Bool_t Scale(UInt_t pdg)
  {
    TH1* h = static_cast<TH1*>(GetOutputObject(dNdptName(pdg), TH1::Class()));
    if (!h) {
      SetStatus(-1);
      Warning("Terminate", "No %s histogram found", dNdptName(pdg));
      return false;
    }
    h->Scale(1./fOK, "width");
    return true;
  }
  Bool_t Ratio(TH1* r, Int_t bin, UInt_t pdg)
  {
    TH1* h = static_cast<TH1*>(GetOutputObject(RatioName(pdg),
					       TH1::Class()));
    if (!h) {
      SetStatus(-1);
      Warning("Terminate", "No %s histogram found", RatioName(pdg));
      return false;
    }
    r->SetBinContent(bin, h->GetMean());
    r->SetBinError(bin, h->GetRMS());
    return true;
  }
  /** 
   * Final processing.  Scales the histogram to the nubmer of events
   * and the bin width.
   */
  virtual void Terminate()
  {
    fOK = GetEventCount();
    if (fOK <= 0) {
      SetStatus(-1);
      Warning("Terminate", "No events selected");
      return;
    }
    Printf("A total of %ld events", fOK);

    Scale(211);
    Scale(321);
    Scale(2212);
    Scale(0);

    TH1* ratios = new TH1D("ratios","Ratios to #pi",
			   3, .5, 3.5);
    ratios->GetXaxis()->SetBinLabel(1, "K/#pi");
    ratios->GetXaxis()->SetBinLabel(2, "p/#pi");
    ratios->GetXaxis()->SetBinLabel(3, "other/#pi");
    ratios->SetFillColor(kCyan+2);
    ratios->SetLineColor(kCyan+2);
    ratios->SetMarkerColor(kCyan+2);
    ratios->SetMarkerStyle(20);
    Ratio(ratios, 1, 321);
    Ratio(ratios, 2, 2212);
    Ratio(ratios, 3, 0);

    fOutput->Add(ratios);
  }    
  Long_t fNpi;
  Long_t fNK;
  Long_t fNP;
  Long_t fNO;

  TH1* fdNdptPi;
  TH1* fdNdptK;
  TH1* fdNdptP;
  TH1* fdNdptO;
  
  TH1* fK2piRatio;
  TH1* fP2piRatio;
  TH1* fO2piRatio;
  
  ClassDef(spectraAnalysis,0);
};

//====================================================================
/*
 * The function to make our analyser 
 */
struct spectraMaker : public FastAnalysis::Maker
{
  spectraMaker() : FastAnalysis::Maker("spectra") {}
  
  FastAnalysis*  Make(const TString& subtype,
		      Int_t          monitor,
		      Bool_t         verbose,
		      TString&       uopt)
  {
    return new spectraAnalysis(verbose,monitor);
  }
  void List() const
  {
  }
  const char* Script() const { return __FILE__; }
};



  
