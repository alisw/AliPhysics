#ifndef SPECTRAANALYSIS_C
#define SPECTRAANALYSIS_C
#ifndef __CINT__
# include <TH1.h>
# include <TH2.h>
# include <TMath.h>
# include <TParticle.h>
# include <TObjArray.h>
# include <TString.h>
# include <TParticlePDG.h>
# include <TGraphErrors.h>
#else
class TH1;
class TH2;
class TParticle;
class TMap;
#endif
#include "FastAnalysis.C"
#include "FastCentHelper.C"

namespace spectra {
  /** 
   * Analyser to analyse the pT spectra and particle ratios 
   */
  struct Base : public FastAnalysis
  {
    enum {
      kPiSty = 20,
      kKSty  = 21,
      kPSty  = 22,
      kOSty  = 23,
      kASty  = 24
    };
    enum {
      kPiColor = kRed+2,
      kKColor  = kGreen+2,
      kPColor  = kBlue+2,
      kOColor  = kMagenta+2,
      kAColor  = kBlack
    };
    /** 
     * Constructor 
     * 
     * @param verbose 
     * @param monitor 
     */
    Base(Bool_t verbose=false, Int_t monitor=0)
      : FastAnalysis(verbose,monitor),
	fNpi(0),
	fNK(0),
	fNP(0),
	fNO(0),
	fdNdpt(0),
	fdNdptPi(0),
	fdNdptK(0),
	fdNdptP(0),
	fdNdptO(0),
	fdNdm(0),
	fK2piRatio(0),
	fP2piRatio(0),
	fO2piRatio(0),
	fd2Ndptdy(0),
	fd2Ndmdy(0),
	fMaxY(10)
    {
    }
    /** 
     * Create a spectrum histogram 
     * 
     * @param pdg 
     * @param col 
     * @param sty 
     * 
     * @return 
     */
    static TH1* CreateSpectra(UInt_t pdg, Color_t col, Style_t sty)
    {
      TH1* h = new TH1D(dNdptName(pdg), "", 1000, 0, 10);
      h->SetXTitle("p_{T}");
      h->SetYTitle("dN/dp_{T}");
      h->SetDirectory(0);
      h->SetMarkerColor(col);
      h->SetLineColor(col);
      h->SetFillColor(col);
      h->SetMarkerStyle(sty);
      h->Sumw2();
      return h;
    }
    static TH1* CreateMt(UInt_t pdg, Color_t col, Style_t sty)
    {
      TH1* h = new TH1D(dNdmtName(pdg), "", 1000, 0, 10);
      h->SetXTitle("m_{T}");
      h->SetYTitle("dN/dm_{T}");
      h->SetDirectory(0);
      h->SetMarkerColor(col);
      h->SetLineColor(col);
      h->SetFillColor(col);
      h->SetMarkerStyle(sty);
      h->Sumw2();
      return h;
    }
    static TH1* CreateMass()
    {
      TH1* h = new TH1D("dNdm", "Mass distribution", 300, 0, 3);
      h->SetXTitle("m\\ \\mathrm{(GeV)}");
      h->SetYTitle("\\mathrm{d}N/\\mathrm{d}m");
      h->SetDirectory(0);
      h->SetMarkerColor(kBlack);
      h->SetLineColor(kBlack);
      h->SetFillColor(kBlack);
      h->SetMarkerStyle(20);
      h->Sumw2();
      return h;
    }      
    /** 
     * Create a ratio histogram 
     * 
     * @param pdg 
     * @param col 
     * @param sty 
     * 
     * @return 
     */
    static TH1* CreateRatio(UInt_t pdg, Color_t col, Style_t sty)
    {
      TH1* h = new TH1D(RatioName(pdg), "", 100, 0, 1);
      h->SetXTitle(Form("%d/211",pdg));
      h->SetYTitle("Events");
      h->SetDirectory(0);
      h->SetMarkerColor(col);
      h->SetLineColor(col);
      h->SetFillColor(col);
      h->SetMarkerStyle(sty);
      h->Sumw2();
      return h;
    }
    static TObject* CreateOutput()
    {
      TList* ret = new TList;
      ret->SetName("all");
      ret->SetOwner(true);
      ret->Add(CreateSpectra(-1,    kAColor, kASty));
      ret->Add(CreateSpectra(211,   kPiColor,kPiSty));
      ret->Add(CreateSpectra(321,   kKColor, kKSty));
      ret->Add(CreateSpectra(2212,  kPColor, kPSty));
      ret->Add(CreateSpectra(0,     kOColor, kOSty));
      ret->Add(CreateMass());
      ret->Add(CreateRatio(321,  kKColor, kKSty));
      ret->Add(CreateRatio(2212, kPColor, kPSty));
      ret->Add(CreateRatio(0,    kOColor, kOSty));

      TH2* d2Ndptdy = new TH2D("d2Ndptdy","Spectra", 20, -10, +10, 100, 0, 10);
      d2Ndptdy->SetDirectory(0);
      d2Ndptdy->SetXTitle("\\mathit{y}");
      d2Ndptdy->SetYTitle("p_{\\mathrm{T}}");
      d2Ndptdy->SetZTitle("\\mathrm{d}^{2}N_{\\mathrm{ch}}/"
			  "\\mathrm{d}p_{\\mathrm{T}}\\mathrm{d}y");
      ret->Add(d2Ndptdy);

      TH2* d2Ndmdy = new TH2D("d2Ndmdy", "Mass", 20, -10, +10, 300, 0, 3);
      d2Ndmdy->SetDirectory(0);
      d2Ndmdy->SetXTitle("\\mathit{y}");
      d2Ndmdy->SetYTitle("\\mathit{m}\\ \\hbox{GeV}");
      d2Ndmdy->SetZTitle("Events");
      ret->Add(d2Ndmdy);
      
      return ret;
    }      
    /** 
     * Called on each slave at start
     * of processing.  Histograms
     * should be added to the output
     * here.
     */
    virtual void SlaveBegin(TTree*)
    {
      Info("SlaveBegin", "Making dN/dpt histogram");
      TList*  l  = static_cast<TList*>(CreateOutput());
      if (!l) {
	Warning("SlaveBegin", "Failed to make histograms");
	SetStatus(-1);
	return;
      }
      l->SetOwner(false);
      l->ls();
      fdNdpt     = static_cast<TH1*>(l->At(0));
      fdNdptPi   = static_cast<TH1*>(l->At(1));
      fdNdptK    = static_cast<TH1*>(l->At(2));
      fdNdptP    = static_cast<TH1*>(l->At(3));
      fdNdptO    = static_cast<TH1*>(l->At(4));
      fdNdm      = static_cast<TH1*>(l->At(5));
      fK2piRatio = static_cast<TH1*>(l->At(6));
      fP2piRatio = static_cast<TH1*>(l->At(7));
      fO2piRatio = static_cast<TH1*>(l->At(8));
      fd2Ndptdy  = static_cast<TH2*>(l->At(9));
      fd2Ndmdy   = static_cast<TH2*>(l->At(10));
      fOutput->AddAll(l);
    }
    /** 
     * Process header.  Implement to make event selection 
     * 
     * 
     * @return true if event is accepted 
     */
    Bool_t ProcessHeader() { return true; }
    /** 
     * Process all particles 
     * 
     */
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
      Double_t y   = p->Y();
      Double_t pT  = p->Pt();
      Double_t m   = p->GetPDG()->Mass();
      Int_t    pdg = TMath::Abs(p->GetPdgCode());
      fd2Ndptdy->Fill(y, pT);
      fd2Ndmdy ->Fill(y, m);

      if (TMath::Abs(y) > fMaxY) return true;

      fdNdpt->Fill(pT);
      fdNdm ->Fill(m);
      if      (pdg == 211)     { fdNdptPi->Fill(pT); fNpi++; }
      else if (pdg == 321)     { fdNdptK ->Fill(pT); fNK++; }
      else if (pdg == 2212)    { fdNdptP ->Fill(pT); fNP++; }
      else                     { fdNdptO ->Fill(pT); fNO++; }
      return true;
    }
    /** 
     * Name of spectrum histogram 
     * 
     * @param pdg Particle code 
     * 
     * @return Name 
     */
    static const char* dNdptName(UInt_t pdg) 
    {
      switch (pdg) {
      case -1:   return "dNdpt";
      case 211:  return "dNdptPi";
      case 321:  return "dNdptK";
      case 2212: return "dNdptP";
      }
      return "dNdptO";
    }
    /** 
     * Name of spectrum histogram 
     * 
     * @param pdg Particle code 
     * 
     * @return Name 
     */
    static const char* dNdmtName(UInt_t pdg) 
    {
      switch (pdg) {
      case -1:   return "dNdmt";
      case 211:  return "dNdmtPi";
      case 321:  return "dNdmtK";
      case 2212: return "dNdmtP";
      }
      return "dNdmtO";
    }
    /** 
     * Name of ratio histograms
     * 
     * @param pdg Particle code 
     * 
     * @return Name 
     */
    static const char* RatioName(UInt_t pdg) 
    {
      switch (pdg) {
      case 321:  return "KtoPiRatio";
      case 2212: return "PtoPiRatio";
      }
      return "OtoPiRatio";
    }
    /** 
     * Scale a spectrum by the number of events and the bin width
     * 
     * @param pdg Particle number 
     * 
     * @return true on succes
     */
    static Bool_t Scale(TList* l, UInt_t pdg, Int_t nOK)
    {
      TH1* h = static_cast<TH1*>(l->FindObject(dNdptName(pdg)));
      if (!h) {
	::Warning("Scale", "No %s histogram found", dNdptName(pdg));
	return false;
      }
      Double_t intg = h->Integral();
      h->Scale(1./nOK/intg, "width");
      return true;
    }
    /** 
     * Fill the mean and rms of a ratio into a histogram 
     * 
     * @param r    Histogram to fill 
     * @param bin  Bin number 
     * @param pdg  Particle ID
     * 
     * @return 
     */
    static Bool_t Ratio(TList* l, TH1* r, Int_t bin, UInt_t pdg)
    {
      TH1* h = static_cast<TH1*>(l->FindObject(RatioName(pdg)));
      if (!h) {
	::Warning("Ratio", "No %s histogram found", RatioName(pdg));
	return false;
      }
      r->SetBinContent(bin, h->GetMean());
      r->SetBinError(bin, h->GetRMS());
      return true;
    }
    static TH1* Normalize(TObject* o, Int_t nEv)
    {
      // Printf("Normalising spectra to %d events");
      TList* l = static_cast<TList*>(o);
      Scale(l, -1,  nEv);
      Scale(l, 211, nEv);
      Scale(l, 321, nEv);
      Scale(l, 2212,nEv);
      Scale(l, 0,   nEv);

      TH1* ratios = new TH1D("ratios","Ratios to #pi", 3, .5, 3.5);
      ratios->GetXaxis()->SetBinLabel(1, "K/#pi");
      ratios->GetXaxis()->SetBinLabel(2, "p/#pi");
      ratios->GetXaxis()->SetBinLabel(3, "other/#pi");
      ratios->SetFillColor(kCyan+2);
      ratios->SetLineColor(kCyan+2);
      ratios->SetMarkerColor(kCyan+2);
      ratios->SetMarkerStyle(20);
      ratios->SetDirectory(0);
      Ratio(l, ratios, 1, 321);
      Ratio(l, ratios, 2, 2212);
      Ratio(l, ratios, 3, 0);
      l->Add(ratios);

      TH2*      d2Ndptdy = static_cast<TH2D*>(l->FindObject("d2Ndptdy"));
      TProfile* meanPt   = d2Ndptdy->ProfileX("meanPt");
      meanPt->SetTitle("\\langle p_{\\mathrm{T}}\\rangle");
      meanPt->SetDirectory(0);
      meanPt->SetXTitle("\\mathit{y}");
      meanPt->SetYTitle("\\langle p_{\\mathrm{T}}\\rangle");
      l->Add(meanPt);

      TH2*      d2Ndmdy = static_cast<TH2D*>(l->FindObject("d2Ndmdy"));
      TProfile* meanM   = d2Ndmdy->ProfileX("meanM");
      meanM->SetTitle("\\langle p_{\\mathrm{T}}\\rangle");
      meanM->SetDirectory(0);
      meanM->SetXTitle("\\mathit{y}");
      meanM->SetYTitle("\\langle p_{\\mathrm{T}}\\rangle");
      l->Add(meanM);

      TH1* mpvPt    = new TH1D("mpvPt","\\lambda(p_{\\mathrm{T}}",
			       meanPt->GetNbinsX(),
			       meanPt->GetXaxis()->GetXmin(),
			       meanPt->GetXaxis()->GetXmax());
      mpvPt->SetDirectory(0);
      mpvPt->SetStats(0);
      mpvPt->SetXTitle("\\mathit{y}");
      mpvPt->SetYTitle(mpvPt->GetTitle());
      mpvPt->SetLineColor(meanPt->GetLineColor());
      mpvPt->SetMarkerColor(meanPt->GetMarkerColor());
      mpvPt->SetMarkerStyle(meanPt->GetMarkerStyle()+4);
      mpvPt->SetMarkerSize(meanPt->GetMarkerSize());
      for (Int_t i = 1; i <= mpvPt->GetNbinsX(); i++) {
	d2Ndptdy->GetXaxis()->SetRange(i,i);
	Int_t    ix, iy, iz;
	Int_t    maxBin = d2Ndptdy->GetMaximumBin(ix, iy, iz);
	Double_t mpv    = d2Ndptdy->GetYaxis()->GetBinCenter(iy);
	Double_t empv   = d2Ndptdy->GetYaxis()->GetBinWidth(iy);
	mpvPt->SetBinContent(i, mpv);
	mpvPt->SetBinError  (i, empv);
      }
      l->Add(mpvPt);
      
      TH1* mpvM    = new TH1D("mpvM","\\lambda(m)",
			       meanM->GetNbinsX(),
			       meanM->GetXaxis()->GetXmin(),
			       meanM->GetXaxis()->GetXmax());
      mpvM->SetDirectory(0);
      mpvM->SetStats(0);
      mpvM->SetXTitle("\\mathit{y}");
      mpvM->SetYTitle(mpvM->GetTitle());
      mpvM->SetLineColor(meanM->GetLineColor());
      mpvM->SetMarkerColor(meanM->GetMarkerColor());
      mpvM->SetMarkerStyle(meanM->GetMarkerStyle()+4);
      mpvM->SetMarkerSize(meanM->GetMarkerSize());
      for (Int_t i = 1; i <= mpvM->GetNbinsX(); i++) {
	d2Ndmdy->GetXaxis()->SetRange(i,i);
	Int_t    ix, iy, iz;
	Int_t    maxBin = d2Ndmdy->GetMaximumBin(ix, iy, iz);
	Double_t mpv    = d2Ndmdy->GetYaxis()->GetBinCenter(iy);
	Double_t empv   = d2Ndmdy->GetYaxis()->GetBinWidth(iy);
	mpvM->SetBinContent(i, mpv);
	mpvM->SetBinError  (i, empv);
      }
      l->Add(mpvM);
      
      return ratios;
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

      Normalize(fOutput, fOK);
    }    
    /** 
     * Get the list of monitored objects 
     * 
     * @return The list of monitored objects 
     */
    virtual TList* GetMonitorObjects()
    {
      TObject* m1 = new TNamed(dNdptName(-1), "");
      m1->SetUniqueID(0x8); // Scale 

      TList* ret = new TList;
      ret->Add(m1);
    
      return ret;
    }
    void Print(Option_t* option="") const
    {
      FastAnalysis::Print(option);
      Printf(" Max rapidity:       +/-%f", fMaxY);
    }
    Long_t fNpi;
    Long_t fNK;
    Long_t fNP;
    Long_t fNO;

    TH1* fdNdpt;       //!
    TH1* fdNdptPi;     //!
    TH1* fdNdptK;      //!
    TH1* fdNdptP;      //!
    TH1* fdNdptO;      //!

    TH1* fdNdm;        //!
  
    TH1* fK2piRatio;   //!
    TH1* fP2piRatio;   //!
    TH1* fO2piRatio;   //!

    TH2* fd2Ndptdy;    //!
    TH2* fd2Ndmdy;     //! 
    
    Double_t fMaxY;
    
    ClassDef(Base,1);
  };

  //====================================================================
  /**
   * Processes V0AND events and build the @f$ 1/N dN_{ch}/d\eta@f$ 
   * 
   */
  struct V0AND : public Base
  {
    /** 
     * Constructor 
     * 
     * @param verbose  Whether to be verbose 
     * @param monitor  Frequency 
     */  
    V0AND(Bool_t verbose=false, Int_t monitor=0)
      : Base(verbose,monitor)
    {
      SetTrigger("V0AND");
    }
    /** 
     * Process the header.  
     * 
     * @return Always true - by definition INEL is all events. 
     */
    virtual Bool_t ProcessHeader()
    {
      return CheckTrigger();
    }
    ClassDef(V0AND,1);
  };
  //====================================================================
  struct Cent : public Base 
  {
    FastCentHelper fHelper;
    const Long_t   fMinEvents;
    Int_t          fCentBin;
    /** 
     * Constructor 
     * 
     * @param method   Centrality method 
     * @param verbose  Verbosity flag 
     * @param monitor  Frequency 
     */
    Cent(const char* method="V0M", Bool_t verbose=true, Int_t monitor=0)
      : Base(verbose,monitor),
	fHelper(method),
	fMinEvents(100),
	fCentBin(0)	
    {
      fTrigMask = 0x0;
      fCentMethod = fHelper.fCentMeth;
    }
    /** 
     * Called on each slave at start of processing
     * 
     * @param t Ignored
     */
    virtual void SlaveBegin(TTree* t)
    {
      Base::SlaveBegin(t);
      // We use the parent fdNdeta histogram as a cache 
      fHelper.CreateHistos(fOutput, Base::CreateOutput);
      fOutput->ls();
    }
    virtual Bool_t SetupEstimator()
    {
      if (!Base::SetupEstimator()) return false;
    
      fHelper.CreateDiagnostics(fOutput, fCentHist);
      return true;
    }    
		     
    /** 
     * Clear our internal caches
     */
    virtual void Clear(Option_t* option="") 
    {
      Base::Clear(option);
      fdNdpt    ->Reset();
      fdNdptPi  ->Reset();
      fdNdptK   ->Reset();
      fdNdptP   ->Reset();
      fdNdptO   ->Reset();
      fdNdm        ->Reset();
      fK2piRatio->Reset();
      fP2piRatio->Reset();
      fO2piRatio->Reset();
      fd2Ndptdy ->Reset();
      fd2Ndmdy  ->Reset();
    }
    /** 
     * Process the header.  Accepts events in range 
     * 
     * @return true 
     */
    virtual Bool_t ProcessHeader()
    {
      if (!CheckTrigger()) return false;
      Double_t cent = GetCentrality();
      fCentBin      = fHelper.CheckCentrality(cent,
					      fEventMult,
					      fHeader->fB,
					      fHeader->fNtgt+fHeader->fNproj,
					      fHeader->fNbin);
      return fCentBin >= 0;
    }
    virtual void FillOut(TList* out)
    {
      static_cast<TH1*>(out->At(0)) ->Add(fdNdpt);
      static_cast<TH1*>(out->At(1)) ->Add(fdNdptPi);
      static_cast<TH1*>(out->At(2)) ->Add(fdNdptK);
      static_cast<TH1*>(out->At(3)) ->Add(fdNdptP);
      static_cast<TH1*>(out->At(4)) ->Add(fdNdptO);
      static_cast<TH1*>(out->At(5)) ->Add(fdNdm);
      static_cast<TH1*>(out->At(6)) ->Add(fK2piRatio);
      static_cast<TH1*>(out->At(7)) ->Add(fP2piRatio);
      static_cast<TH1*>(out->At(8)) ->Add(fO2piRatio);
      static_cast<TH1*>(out->At(9)) ->Add(fd2Ndptdy);
      static_cast<TH1*>(out->At(10))->Add(fd2Ndmdy);
    }
    /** 
     * Process a single event.  
     *
     * First we fill the internal cache using the base class methods.
     * Then we find the histogram for this particular reference
     * multiplicity and add our event cache to that bin.  The number of
     * events in each bin is counted in the unique ID of each bin.
     */
    virtual void ProcessParticles()
    {
      // Check we got a bin 
      if (fCentBin < 0) return;

      // Find the histogram to update 
      TList* out  = static_cast<TList*>(fHelper.fCentList->At(fCentBin-1));
      // If we still have no histogram, return immediately 
      if (!out) return;

      // Use parent function to fill cache 
      Base::ProcessParticles();

      // Add our cache to the appropriate bin 
      FillOut(out);
    }
    virtual void SetAttr(TH1* h, const char* title, Color_t col, Style_t sty)
    {
      h->SetDirectory(0);
      h->SetTitle(title);
      h->SetYTitle(title);
      h->SetMarkerColor(col);
      h->SetMarkerStyle(sty);
      h->SetMarkerSize(2);
      h->SetLineColor(col);
      h->SetStats(0);
    }
	
    /** 
     * Final processing. 
     * 
     * Normalize each bin to the number of events in each bin (stored in
     * the unique ID of that bin).
     */
    virtual void Terminate()
    {
      fOK = GetEventCount();
      if (fOK <= 0) {
	SetStatus(-1);
	Warning("Terminate", "No events selected");
	return;
      }

      if (!fHelper.Finalize(fOutput, fMinEvents, Base::Normalize))
	SetStatus(-1);

      
      TH1* meanPt = static_cast<TH1*>(fHelper.fCentAll->Clone("cmeanPt"));
      meanPt->Reset();
      SetAttr(meanPt, "\\langle p_{\\mathrm{T}}\\rangle", kKColor, kKSty);
      fOutput->Add(meanPt);

      TH1* mpvPt = static_cast<TH1*>(meanPt->Clone("cmpvPt"));
      mpvPt->Reset();
      SetAttr(mpvPt,"\\lambda\\left(p_{\\mathrm{T}}\\right)",kKColor,kKSty+4);
      fOutput->Add(mpvPt);

      TH1* meanM = static_cast<TH1*>(fHelper.fCentAll->Clone("cmeanM"));
      meanM->Reset();
      SetAttr(meanM,"\\langle m\\rangle",kPColor,kPSty);
      fOutput->Add(meanM);

      TH1* mpvM = static_cast<TH1*>(meanM->Clone("cmpvM"));
      mpvM->Reset();
      SetAttr(mpvM,"\\lambda(m)",kPColor,kPSty+4);
      fOutput->Add(mpvM);
      
      TH1* pt2m = static_cast<TH1*>(fHelper.fCentAll->Clone("cpt2m"));
      pt2m->Reset();
      SetAttr(pt2m,"\\langle p_{\\mathrm{T}}\\rangle/\\langle m\\rangle",
	      kPiColor, kPiSty);
      fOutput->Add(pt2m);

      TH1* lpt2m = static_cast<TH1*>(pt2m->Clone("clpt2m"));
      lpt2m->Reset();
      SetAttr(lpt2m, "\\lambda\\left(p_{\\mathrm{T}}\\right)/\\lambda(m)",
	      kKColor, kKSty);
      fOutput->Add(lpt2m);

      TH1* mpt2m = static_cast<TH1*>(pt2m->Clone("cmpt2m"));
      mpt2m->Reset();
      SetAttr(mpt2m,"\\lambda\\left(p_{\\mathrm{T}}\\right)/\\langle m\\rangle",
	      kPColor, kPSty);
      fOutput->Add(mpt2m);

      TH1* k2pi = static_cast<TH1*>(fHelper.fCentAll->Clone("ck2pi"));
      k2pi->Reset();
      SetAttr(k2pi,"\\langle K/\\pi\\rangle",kKColor,kKSty);
      fOutput->Add(k2pi);

      TH1* p2pi = static_cast<TH1*>(fHelper.fCentAll->Clone("cp2pi"));
      p2pi->Reset();
      SetAttr(p2pi,"\\langle p/\\pi \\rangle",kPColor,kPSty);
      fOutput->Add(p2pi);

      TH1* o2pi = static_cast<TH1*>(fHelper.fCentAll->Clone("co2pi"));
      o2pi->Reset();
      SetAttr(o2pi,"\\langle X/\\pi \\rangle",kOColor,kOSty);
      fOutput->Add(o2pi);

      for (Int_t i = 1; i <= meanPt->GetNbinsX(); i++) {
	TList* out     = static_cast<TList*>(fHelper.CentCollection(i));
	if (!out) continue;
	TH1*  dNdpt     = static_cast<TH1*>(out->At(0));
	TH1*  dNdptPi   = static_cast<TH1*>(out->At(1));
	TH1*  dNdptK    = static_cast<TH1*>(out->At(2));
	TH1*  dNdptP    = static_cast<TH1*>(out->At(3));
	TH1*  dNdptO    = static_cast<TH1*>(out->At(4));
	TH1*  dNdm      = static_cast<TH1*>(out->At(5));
	TH1*  k2piRatio = static_cast<TH1*>(out->At(6));
	TH1*  p2piRatio = static_cast<TH1*>(out->At(7));
	TH1*  o2piRatio = static_cast<TH1*>(out->At(8));
	TH2*  d2Ndptdy  = static_cast<TH2*>(out->At(9));
	TH2*  d2Ndmdy   = static_cast<TH2*>(out->At(10));
	TH1*  rat       = static_cast<TH1*>(out->At(11));
	Int_t bPt       = dNdpt->GetMaximumBin();
	Int_t bM        = dNdm ->GetMaximumBin();
	meanPt->SetBinContent(i, dNdpt->GetMean());
	meanPt->SetBinError  (i, dNdpt->GetMeanError());
	meanM ->SetBinContent(i, dNdm ->GetMean());
	meanM ->SetBinError  (i, dNdm ->GetMeanError());
	mpvPt->SetBinContent (i, dNdpt->GetXaxis()->GetBinCenter(bPt));
	mpvPt->SetBinError   (i, dNdpt->GetXaxis()->GetBinWidth (bPt));
	mpvM ->SetBinContent (i, dNdm ->GetXaxis()->GetBinCenter(bM));
	mpvM ->SetBinError   (i, dNdm ->GetXaxis()->GetBinWidth (bM)); 
	k2pi  ->SetBinContent(i, rat  ->GetBinContent(1));
	k2pi  ->SetBinError  (i, rat  ->GetBinError  (1));
	p2pi  ->SetBinContent(i, rat  ->GetBinContent(2));
	p2pi  ->SetBinError  (i, rat  ->GetBinError  (2));
	o2pi  ->SetBinContent(i, rat  ->GetBinContent(3));
	o2pi  ->SetBinError  (i, rat  ->GetBinError  (3));
      }
      pt2m->Add(meanPt);
      pt2m->Divide(meanM);
      lpt2m->Add(mpvPt);
      lpt2m->Divide(mpvM);
      mpt2m->Add(mpvPt);
      mpt2m->Divide(meanM);

      Double_t fit[] = { 1.61,1.60,1.57,1.50,1.47,1.48,1.51,1.52,1.49,1.47,0};
      Double_t err[] = { 0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.03,0.03,0.04,0};
      TProfile* nPart =
	static_cast<TProfile*>(fOutput->FindObject("centNPartMean"));
      TGraphErrors* ptm = new TGraphErrors;
      ptm->SetName("cptm");
      ptm->SetTitle("\\langle p_{\\mathrm{T}}\\rangle/\\langle m\\rangle");
      ptm->SetMarkerColor(kPiColor);
      ptm->SetMarkerStyle(kPiSty);
      ptm->SetLineColor(kPiColor);
      fOutput->Add(ptm);
      TGraphErrors* lptm = new TGraphErrors;
      lptm->SetName("clptm");
      lptm->SetTitle("\\lambda\\left(p_{\\mathrm{T}}\\right)/\\lambda(m)");
      lptm->SetMarkerColor(kPiColor);
      lptm->SetMarkerStyle(kPiSty+4);
      lptm->SetLineColor(kPiColor);
      fOutput->Add(lptm);      
      for (Int_t i = 1; i <= nPart->GetNbinsX(); i++) {
	Double_t nP    = nPart->GetBinContent(i);
	Double_t eP    = nPart->GetBinError  (i);
	Double_t p2m   = pt2m ->GetBinContent(i);
	Double_t ep2m  = pt2m ->GetBinError  (i);
	Double_t lp2m  = lpt2m->GetBinContent(i);
	Double_t lep2m = lpt2m->GetBinError  (i);
	ptm->SetPoint      (i-1, nP, p2m);
	ptm->SetPointError (i-1, eP, ep2m);
	lptm->SetPoint     (i-1, nP, lp2m);
	lptm->SetPointError(i-1, eP, lep2m);
	Printf("%2d  Npart=%5.1f+/-%4.1f  <pT>/<m>=%5.3f+/-%5.3f  "
	       "l(pT)/l(m)=%5.3f+/-%5.3f  a=%5.3f+/-%5.3f",
	       i, nP, eP, p2m, ep2m, lp2m, lep2m, fit[i-1], err[i-1]);
      }
      TMultiGraph* mg = new TMultiGraph;
      mg->SetName("mg");
      mg->SetTitle("Estimates");
      mg->Add(ptm);
      mg->Add(lptm);
      fOutput->Add(mg);
      
      THStack* means = new THStack("meansCent","Particle means");
      means->Add(meanPt);
      means->Add(meanM);
      fOutput->Add(means);

      THStack* ptms = new THStack("ptmsCent","Particle ptms");
      ptms->Add(pt2m);
      ptms->Add(lpt2m);
      ptms->Add(mpt2m);
      fOutput->Add(ptms);
      
      THStack* ratios = new THStack("ratiosCent","Particle ratios");
      ratios->Add(k2pi);
      ratios->Add(p2pi);
      ratios->Add(o2pi);
      fOutput->Add(ratios);

      ptm->SaveAs("ptm.C");
    }
    /** 
     * Get the list of monitored objects 
     * 
     * @return The list of monitored objects 
     */
    virtual TList* GetMonitorObjects()
    {
      TObject* m1 = new TNamed("cent",     "hist text30");
      TObject* m2 = new TNamed("centAcc",  "hist text30");
      TObject* m3 = new TNamed("byCent",   "e");
      TObject* m0 = new TNamed(dNdptName(-1), "");
      m0->SetUniqueID(0x8); // Scale     
      m3->SetUniqueID(0x8); // Scale 
      TList* ret = new TList;
      ret->Add(m0);
      ret->Add(m1);
      ret->Add(m2);
      ret->Add(m3);
    
      return ret;
    }
    void Print(Option_t* option="") const
    {
      Base::Print(option);
      Printf(" Least # events:        %d", fMinEvents);
      fHelper.Print(option);
    }
    ClassDef(Cent,1);
  };

}
//====================================================================
/*
 * The function to make our analyser 
 */
struct spectraMaker : public FastAnalysis::Maker
{
  /**
   * CTOR 
   */
  spectraMaker() : FastAnalysis::Maker("spectra") {}
  /** 
   * Create an analyser or return 0. 
   * 
   * @param subtype  Sub-type to make  
   * @param monitor  Monitor period (in seconds)
   * @param verbose  Verbosity 
   * @param uopt     Possibly additiona options 
   *  
   * @return Newly allocated analyser or null
   */
  FastAnalysis*  Make(const TString& subtype,
		      Int_t          monitor,
		      Bool_t         verbose,
		      TMap&          uopt)
  {
    spectra::Base* ret = 0;
    TString t(subtype);
    if      (t.EqualTo("V0AND"))   ret = new spectra::V0AND(verbose,monitor);
    else if (t.BeginsWith("CENT")) {
      TString w(t(4, t.Length()-4));
      if (!(w.BeginsWith("RefMult") ||
	    w.BeginsWith("ZNA") || 
	    w.BeginsWith("ZNC") || 
	    w.BeginsWith("ZPA") || 
	    w.BeginsWith("ZPC") || 
	    w.BeginsWith("V0M") ||
	    w.BeginsWith("V0A") ||
	    w.BeginsWith("V0C") ||
	    w.BeginsWith("B")   ||
	    w.IsNull())) {
	Printf("Warning: spectraMaker::Make: Unknown estimator: %s",
	       w.Data());
	return 0;
      }
      ret = new spectra::Cent(w, verbose, monitor);
    }
    else
      ret = new spectra::Base(verbose, monitor);
    
    if (ret) {
      TPair* tp = static_cast<TPair*>(uopt.FindObject("trig"));
      if (tp) ret->SetTrigger(tp->Value()->GetName());
      TPair* my = static_cast<TPair*>(uopt.FindObject("maxy"));
      if (my) ret->fMaxY = TString(my->Value()->GetName()).Atof();
    }
    else 
      Printf("Error: spectraMaker::Run: Invalid spec: %s", t.Data());
    return ret;
  }
  /** 
   * Show list of possible sub-types
   * 
   */
  void List() const
  {
    Printf(" <default>       - Inelastic");
    Printf(" V0AND           - Visible X-section");
    Printf(" CENT<est>       - Centrality classes. <est> is one of ");
    Printf("  ZNA            - ZNA signal");
    Printf("  ZNC            - ZNC signal");
    Printf("  ZPA            - ZPA signal");
    Printf("  ZPC            - ZPC signal");
    Printf("  V0M            - V0-A + -C");
    Printf("  V0A            - V0-A");
    Printf("  V0C            - V0-C");
  }
  /** 
   * Get name of script to load 
   */
  const char* Script() const { return __FILE__; }
};


// ------------------------------------------------------------------
// Create instance of maker 
spectraMaker* _spectraMaker = new spectraMaker;

#endif
//
// EOF
// 
