#ifndef FASTCENTESTIMATOR_C
#define FASTCENTESTIMATOR_C
#ifndef __CINT__
# include <TObject.h>
# include <TString.h>
# include <TH1.h>
# include <TMath.h>
# include <TCollection.h>
# include <TTree.h>
# include <TParticle.h>
# include <TError.h>
# include "FastShortHeader.C"
# include <TH2.h>
#else
class TH1;
// class TH2;
class TCollection;
class TTree;
class TParticle;
class FastShortHeader;
#endif 

enum {
  kNeutron = 2112,
  kProton = 2212
};

//====================================================================
/** 
 * Base class for centrality estimators 
 */
struct FastCentEstimator : public TObject
{
  TString fName;
  Bool_t fVerbose;
  /** 
   * Constructor 
   * 
   * @param name Name of the estimator 
   */
  FastCentEstimator(const char* name="")
    : TObject(), fName(name), fVerbose(false)
  {}
  /** 
   * Destructor 
   */
  virtual ~FastCentEstimator() {}
  /**
   * Get the name 
   *
   * @return The name 
   */
  const char* GetName() const { return fName.Data(); }
  void SetVerbose(Bool_t verb) { fVerbose = verb; }
  /** 
   * Set-up this estimator.  Output objects should be stored in @a
   * out, and a branch can be registerd in the TTree.
   * 
   * @param out   Output list to add stuff to 
   * @param tree  Output tree
   * @param sNN   Collision energy in GeV
   * @param tgtA  True if target is a nucleus 
   * @param projA True if projectile is a nucleus 
   */
  virtual void Setup(TCollection* out, TTree* tree, UShort_t sNN,
		     Bool_t tgtA, Bool_t projA) = 0;
  /** 
   * Called before the start of an event 
   * 
   */
  virtual void PreEvent() {}
  virtual void ProcessHeader(FastShortHeader&) {}
  /** 
   * Called for each particle produced by the event generator.
   * Sub-classes should decide if they want to take the information
   * from the passed particle, and then process that information.
   * 
   * @param p Generated particle 
   */
  virtual void Process(const TParticle* p) = 0;
  /** 
   * Called at the end of an event 
   */
  virtual void PostEvent() {}
  /** 
   * Do the final calculations 
   * 
   * @param out Output list to add stuff to 
   */
  virtual void Terminate(TCollection* out) = 0;
  virtual void Print(Option_t* option="") const
  {
    Printf("%s: %s", ClassName(), GetName());
  }
  /** 
   * Get the particle polar angle 
   * 
   * @param p Particle 
   * 
   * @return Polar angle 
   */
  static Double_t Theta(const TParticle* p)
  {
    Double_t pT    = p->Pt();
    Double_t pZ    = p->Pz();
    Double_t theta = TMath::ATan2(pT, pZ);
    return theta;
  }
  /** 
   * Get the particle pseudo-rapidity @f$\eta@f$ 
   * 
   * @param p The particle 
   * 
   * @return Pseudo-rapidity @f$\eta@f$ 
   */
  static Double_t Eta(const TParticle* p)
  {
    Double_t theta = Theta(p);
    Double_t tanth = TMath::Tan(theta/2);
    if (tanth < 1e-100) return +1e10;
    if (tanth > 1e100)  return -1e10;
    Double_t eta   = -TMath::Log(tanth);
    return eta;
  }
  /** 
   * Get the particle azimuth angle @f$\varphi@f$ 
   * 
   * @param p The particle 
   * 
   * @return Azimuth angle @f$\varphi@f$ 
   */
  static Double_t Phi(const TParticle* p)
  {
    Double_t px  = p->Px();
    Double_t py  = p->Py();
    Double_t phi = TMath::ATan2(px, py);
    return phi;
  }
  /** 
   * Check if a particle is primary 
   * 
   * @param p Particle 
   * 
   * @return true if primary 
   */
  static Bool_t IsPrimary(const TParticle* p)
  {
    return p->TestBit(BIT(14));
  }
  /** 
   * Check if particle is from weak decay 
   * 
   * @param p Particle 
   * 
   * @return true if from weak decay 
   */
  static Bool_t IsWeakDecay(const TParticle* p)
  {
    return p->TestBit(BIT(15));
  }
  /** 
   * Check if particle is charged 
   * 
   * @param p Particle 
   * 
   * @return true if particle is charged 
   */
  static Bool_t IsCharged(const TParticle* p)
  {
    return p->TestBit(BIT(16));
  }
  ClassDef(FastCentEstimator,1);
};




//____________________________________________________________________
/** 
 * A 1-dimensional centrality estimator 
 */
struct Fast1DCentEstimator : public FastCentEstimator
{
  /** Sum of signals for a given event */ 
  ULong64_t fCache;
  /** 
   * Our histogram. The Setup member function must be overridden to
   * define this member. 
   */
  TH1* fHistogram; //!
  Bool_t fFromTop; 
  /** 
   * Constructor 
   * 
   * @param name Name of the estimator 
   */
  Fast1DCentEstimator(const char* name="")
    : FastCentEstimator(name), fHistogram(0), fFromTop(true)
  {}
  /** 
   * Destructor 
   */
  virtual ~Fast1DCentEstimator() {}
  virtual const char* MultSpec() const { return "l"; }
  /** 
   * Set-up this object.  Defines the internal histogram and add to
   * output
   * 
   * @param l Output list
   * @param tree Tree to add branch to 
   */
  virtual void Setup(TCollection* l, TTree* tree, UShort_t,
		     Bool_t, Bool_t)
  {
    if (fHistogram && l) {
      Info("Setup", "Adding histogram %s to output",
	   fHistogram->GetName());
      l->Add(fHistogram);
    }
    TString leaves; leaves.Form("value/%s", MultSpec());
    if (tree) tree->Branch(GetName(), &fCache, leaves.Data());
  }
  /** 
   * Called before each event.  Zeros the cache variable 
   */
  virtual void PreEvent()
  {
    fCache = 0;
  }
  /** 
   * Fills the summed observable into the histogram 
   */
  virtual void PostEvent()
  {
    if (fVerbose) Info("PostEvent", " Got %lld %s particles",
		       fCache, GetName());
    fHistogram->Fill(fCache);
  }
  virtual TH1* GetHistogram(TCollection* l) = 0;
  /** 
   * Called at the end of the processing.  The member function
   * produces a histogram like the defined observable histogram, but
   * the bin content corresponds to the centrality corresponding to
   * that bin.  In that way, we can do a simple look-up in the output
   * histogram to extract the centrality
   * 
   * @param out Output list to add information to. 
   */
  virtual void Terminate(TCollection* out)
  {
    TH1* h    = GetHistogram(out);
    TH1* cent = static_cast<TH1*>(h->Clone(GetName()));
    cent->SetDirectory(0);
    cent->SetYTitle("Centrality [%]");
    cent->SetTitle(Form("%s mapping", GetName()));
    cent->Reset();
    out->Add(cent);
    
    Int_t    nX         = h->GetNbinsX();
    Double_t total      = h->Integral(1,nX);
    Int_t    dBin       = (fFromTop ? -1 : 1);
    Int_t    end        = (fFromTop ? 0  : nX+1);
    Int_t    start      = (fFromTop ? nX : 1);
    if (fVerbose)
      Info("Teminate", "Integrating %s from bin %d to 1",
	   h->GetName(), nX);
    for (Int_t i = start; i != end; i += dBin) {
      Double_t curInt  = (fFromTop ?
			  h->Integral(i, start) :
			  h->Integral(start,i));
      if (curInt < 0) continue;
      Double_t curCent = curInt / total * 100;
      cent->SetBinContent(i, curCent);
      if (fVerbose)
	Info("Terminate", "Bin %3d -> %9f/%9f -> %5.1f%%",
	     i, curInt, total, curCent);
    }    
  }
  virtual void Print(Option_t* option="nh") const
  {
    TString opt(option); opt.ToLower();
    if (opt.Contains("n")) 
      Printf("1D Estimator: %s (%screasing)",
	     GetName(),(fFromTop ? "de" : "in"));
    if (opt.Contains("h") && fHistogram) {
      Int_t nBin = fHistogram->GetNbinsX();
      Printf("  %d bins between %f and %f",nBin, 
	     fHistogram->GetXaxis()->GetXmin(),
	     fHistogram->GetXaxis()->GetXmax());
    }
  }
  ClassDef(Fast1DCentEstimator,1);
};

//____________________________________________________________________
struct BCentEstimator : public Fast1DCentEstimator
{
  // TH2* fBvsC;
  const Double_t fkFactor;
  BCentEstimator()
    : Fast1DCentEstimator("B"),
      // fBvsC(0),
      fkFactor(1000)
  {}
  /** 
   * 
   * 
   * @param out 
   * @param tree 
   * @param sNN 
   * @param tgtA 
   * @param projA 
   */
  virtual void Setup(TCollection* out, TTree* tree, UShort_t sNN,
		     Bool_t tgtA, Bool_t projA)
  {
    fHistogram = MakeHistogram(sNN, tgtA, projA);
    Fast1DCentEstimator::Setup(out, tree, sNN, tgtA, projA);
    // if (tree) tree->Branch(GetName(), &fB, "value/D");    
  }
  /** 
   * Get the histogram to accumulate the observable in.  
   * 
   * @return Pointer to the histogram. 
   */
  virtual TH1* GetHistogram(TCollection* l)
  {
    return static_cast<TH1*>(l->FindObject(Form("raw%s",GetName())));
  }
  TH1* MakeHistogram(UShort_t sNN, Bool_t tgtA, Bool_t projA)
  {
    TArrayD cents;
    TArrayD bins; // In B
    if (tgtA && projA) { // Pb-Pb
      if (sNN == 2760) {
	// PbPb @ 2.76TeV only 
	// Updated 4th of November 2014 from 
	// cern.ch/twiki/bin/view/ALICE/CentStudies
	//        #Tables_with_centrality_bins_AN1
	Double_t bs[] = { 0,      1.57,  2.22,  2.71,  3.13,
			  3.50,   4.94,  6.05,  6.98,  7.81,
			  8.55,   9.23,  9.88, 10.47, 11.04,
			  11.58, 12.09, 12.58, 13.05, 13.52,
			  13.97, 14.43, 14.96, 15.67, 20.00 };
	Double_t cs[] = { 0.5,   1.5,   2.5,   3.5,   4.5,
			  7.5,   12.5,  17.5,  22.5,  27.5,
			  32.5,  37.5,  42.5,  47.5,  52.5,
			  57.5,  62.5,  67.5,  72.5,  77.5,
			  82.5,  87.5,  92.5,  97.5 };
	cents.Set(24,cs);
	bins.Set(25,bs);
      }
      else if (sNN == 5023) {
	// PbPb @ 5.02TeV only 
	// https://twiki.cern.ch/twiki/bin/viewauth/ALICE/CentralityCodeSnippets
	Double_t bs[] = { 0.00, 1.56, 2.22, 2.71, 3.13,
			  3.51, 3.84, 4.15, 4.43, 4.71,
			  4.96, 6.08, 7.01, 7.84, 8.59,
			  9.27, 9.92, 10.5, 11.1, 11.6,
			  12.1, 12.6, 13.1, 13.6, 14.0,
			  14.5, 15.0, 15.7, 19.6  }; // 29
	Double_t cs[] = { 0.5,   1.5,  2.5,  3.5,  4.5,
			  5.5,   6.5,  7.5,  8.5,  9.5,
			  12.5, 17.5, 22.5, 27.5, 32.5,
			  37.5, 42.5, 47.5, 52.5, 57.5,
			  62.5, 67.5, 72.5, 77.5, 82.5,
			  87.5, 92.5, 97.5  };	  
	cents.Set(28,cs);
	bins.Set(29,bs);
      }
      else if (sNN == 5440) {
	// Xe-Xe @ 5.44TeV
	// https://twiki.cern.ch/twiki/bin/viewauth/ALICE/XeXeCentStudies
	Double_t bs[] = {
	  0.00, 2.12, 3.00, 3.66, 4.23, 5.18, 5.98, 6.68, 7.32, 7.91, 8.46,
	  8.97, 9.45, 9.91, 10.4, 10.8, 11.2, 11.6, 12.0, 12.4, 12.9, 18.3
	};
	Double_t cs[] = {
	  1.25, 3.75, 6.25, 8.75, 12.50, 17.50, 22.50, 27.50, 32.50, 37.50,
	  42.50, 47.50, 52.50, 57.50, 62.50, 67.50, 72.50, 77.50, 82.50,
	  87.50, 95.50
	};
	cents.Set(21, cs);
	bins.Set(22, bs);
      }
    }
    else if (tgtA || projA) { // p-Pb or Pb-p
      if (sNN == 5023) {
	Double_t cs[] = { 2.5,     7.5,     15.,      30.,
			  50.,     70.,     90. };
	Double_t bs[] = { 0,       1.83675, 2.59375,  3.66875,
			  5.18625, 6.35475, 7.40225, 13.8577 };
	cents.Set(7, cs);
	bins.Set(8, bs);
      }
    }
    if (bins.GetSize() <= 0 || cents.GetSize() <= 0 ) {
      // Nothing defined
      Warning("MakeHistogram", "No bins defined for sNN=%d (%c-%c)",
	      sNN, tgtA ? 'A' : 'p', projA ? 'A' : 'p');
      return 0;
    }
    printf("b bins: ");
    for (Int_t i = 0; i < bins.GetSize(); i++) {
      bins[i] *= fkFactor; // Scale to 1/1000 fm
      printf("%s%7.1f", i!=0 ? "-" : "", bins[i]);
    }
    Printf("");
    TH1* h = new TH1D(Form("raw%s",GetName()), "B to Centrality",
		      bins.GetSize()-1, bins.GetArray());
    h->SetDirectory(0);
    h->SetXTitle("b\\hbox{ [10^{3}fm]}");
    h->SetYTitle("c\\hbox{ [\\%]}");
    h->SetBinContent(0,1);
    for (Int_t i = 1; i <= cents.GetSize(); i++) {
      h->SetBinContent(i, cents[i-1]);
    }
    return h;
  }
  /** 
   * Reset cache 
   * 
   */
  void PreEvent() { fCache = 100*fkFactor; };
  /** 
   * Process header 
   * 
   * @param h Header
   */
  void ProcessHeader(FastShortHeader& h)
  {
    // In Ap (EPOS) we have many spec in the target, meaning they will
    // be detected on the C side (p is the projectile, A is the target)
    // if (!fSpectators) return;
    fCache = h.fB * fkFactor;
    // Info("", "Cache=%ld", fCache);
    // if (fBvsC) fBvsC->Fill(fCache, h.fC);
  }
  /** 
   * Do nothing here 
   */
  virtual void Process(const TParticle*) {};
  /** 
   * Do Nothing 
   * 
   */
  void PostEvent() {};
  /** 
   * Called at the end of the processing.  Just copy mapping to output
   * 
   * @param out Output list to add information to. 
   */
  virtual void Terminate(TCollection* out)
  {
    TH1* h    = GetHistogram(out);
    if (!h) {
      Warning("Terminate", "No histogram on input");
      out->ls();
      return;
    }
    TH1* cent = static_cast<TH1*>(h->Clone(GetName()));
    cent->SetDirectory(0);
    cent->SetYTitle("Centrality [%]");
    cent->SetTitle(Form("%s mapping", GetName()));
    out->Add(cent);
    // Scale to number of workers
    Double_t scale = cent->GetBinContent(0);
    cent->Scale(1./scale);
  }
  /** 
   * Special function for returning centrality early 
   * 
   * 
   * @return Event centrality 
   */
  Double_t GetCentrality(Double_t b) const
  {
    Double_t ret = (fHistogram ?
		    fHistogram->GetBinContent(fHistogram->FindBin(b*fkFactor)) :
		    200);
    // Info("", "Look-up of %f (%f) -> %5.1f%%", b*fkFactor, b, ret);
    return ret;
  }
  ClassDef(BCentEstimator,2);
};
  
//____________________________________________________________________
/**
 * Centrality estimator using charged particles 
 */
struct FastNchCentEstimator : public Fast1DCentEstimator
{
  /** 
   * Constructor 
   * 
   * @param name Name 
   */
  FastNchCentEstimator(const char* name="")
    : Fast1DCentEstimator(name)
  {}
  /** 
   * Destructor 
   */
  virtual ~FastNchCentEstimator() {}
  /** 
   * Process a single particle. 
   * 
   * @param p Particle 
   */
  virtual void Process(const TParticle* p)
  {
    if (!IsCharged(p)) return;
    if (!Accept(p))    return;
    fCache++;
  }
  /** 
   * Must be overloaded. Should return true if we accept the particle 
   *
   * @param p Particle to investigate 
   *
   * @return true if we are to count this particle 
   */
  virtual Bool_t Accept(const TParticle* p) = 0;
  virtual void Print(Option_t* option="nh") const
  {
    TString opt(option); opt.ToLower();
    if (opt.Contains("n")) 
      Printf("1D Nch Estimator: %s (%screasing)",
	     GetName(),(fFromTop ? "de" : "in"));
    opt.ReplaceAll("n", "");
    Fast1DCentEstimator::Print(opt);
  }
  ClassDef(FastNchCentEstimator,1);
};

//____________________________________________________________________
/**
 * Centrality estimator using the V0 signal 
 */
struct V0CentEstimator : public FastNchCentEstimator
{
  /** Mode: Negative, use C side, positive use A side, otherwise sum */
  Short_t  fMode;
  Bool_t   fOnlyPrimary;
  Double_t fAEtaMin;
  Double_t fAEtaMax;
  Double_t fCEtaMin;
  Double_t fCEtaMax;
  /** 
   * Constructor 
   * 
   * @param mode Mode: Negative, use C side, positive use A side, otherwise sum 
   * @param onlyPrimary IF true, only investigate primaries  
   */
  V0CentEstimator(Short_t mode=0, Bool_t onlyPrimary=false) 
    : FastNchCentEstimator(Form("%s%s", 
				(mode < 0 ? "V0C" : 
				 mode > 0 ? "V0A" : "V0M"), 
				(onlyPrimary ? "P" : ""))),
      fMode(mode),
      fOnlyPrimary(onlyPrimary)			   
  {
    fAEtaMin = +2.8;
    fAEtaMax = +5.1;
    fCEtaMin = -3.7;
    fCEtaMax = -1.7;
  }
  void Flip(Bool_t onlySign=true)
  {
    if (!onlySign) {
      std::swap(fAEtaMin, fCEtaMax); // a=[-1.7,*] c=[*,+2.8]
      std::swap(fAEtaMax, fCEtaMin); // a=[*,-3.7] c=[+5.1,*]
    }
    else {
      std::swap(fAEtaMin, fAEtaMax);  // a=[+5.1,+2.8]
      std::swap(fCEtaMin, fCEtaMax);  // c=[-1.7,-3.7]
    }
    //                 onlySign     !onlySign
    fAEtaMin *= -1; // a=[-5.1,*]   a=[+1.7,*]
    fAEtaMax *= -1; // a=[*,-2.8]   a=[*,+3.7]
    fCEtaMin *= -1; // c=[+1.7,*]   c=[-5.1,*]
    fCEtaMax *= -1; // c=[*,+3.7]   c=[*,-2.8]
    Printf("Flipped %s", (onlySign ? "sign" : "acceptance"));
  }
  /** 
   * Set-up this object.  Defines the internal histogram and add to
   * output
   * 
   * @param l Output list
   * @param tree Tree to add branch to 
   * @param sNN   Collision energy in GeV
   * @param tgtA  True if target is a nucleus 
   * @param projA True if projectile is a nucleus 
   */
  void Setup(TCollection* l, TTree* tree, UShort_t sNN,
	     Bool_t tgtA, Bool_t projA)
  {
    Bool_t  isAA  = (tgtA && projA);
    Bool_t  isPA  = (tgtA ^ projA); // XOR
    UInt_t  max   = (isAA ? 13000 : isPA ? 800 : 300);
    UInt_t  dBin  = (isAA ? 10    : isPA ?   1 :   1);
    Color_t color = (fMode < 0 ? kRed : fMode > 0 ? kBlue : kGreen)+2;
    fHistogram = new TH1D(Form("raw%s",GetName()),
			  Form("%s #it{N}_{ch} distribution", GetName()),
			  max/dBin, 0, (fMode == 0 ? 2 : 1)*max);
    fHistogram->SetXTitle("#it{N}_{ch}");
    fHistogram->SetYTitle("Raw #it{P}(#it{N}_{ch})");
    fHistogram->SetDirectory(0);
    fHistogram->SetLineColor(color);
    fHistogram->SetFillColor(color);
    fHistogram->SetMarkerColor(color);
    fHistogram->SetMarkerStyle(20);
    fHistogram->SetFillStyle(3002);

    Fast1DCentEstimator::Setup(l, tree, sNN, tgtA, projA);
  }
  /** 
   * Whether we should accept a particle.  We accept a particle if it
   * falls within the acceptance of the V0.
   * 
   * @param p Particle to investigate 
   * 
   * @return true if to be used 
   */
  Bool_t Accept(const TParticle* p)
  {
    if (fOnlyPrimary && !IsPrimary(p)) return false;
    Double_t eta = Eta(p);
    Bool_t   v0A = ((eta >= fAEtaMin) && (eta <= fAEtaMax));
    Bool_t   v0C = ((eta >= fCEtaMin) && (eta <= fCEtaMax));
    if (fMode < 0) return v0C;
    if (fMode > 0) return v0A;
    return v0A || v0C;
  }
  /** 
   * Get the histogram to accumulate the observable in.  
   * 
   * @return Pointer to the histogram. 
   */
  virtual TH1* GetHistogram(TCollection* l)
  {
    return static_cast<TH1*>(l->FindObject(Form("raw%s",GetName())));
  }
  virtual void Print(Option_t* option="nah") const
  {
    TString opt(option); opt.ToLower();
    if (opt.Contains("n")) 
      Printf("V0 Estimator: %s (%screasing) mode=%d",
	     GetName(),(fFromTop ? "de" : "in"), fMode);
    if (opt.Contains("a"))
	Printf("  A: eta=[%f,%f], C: eta=[%f,%f]",
	       fAEtaMin, fAEtaMax, fCEtaMin, fCEtaMax);
    opt.ReplaceAll("n", "");
    FastNchCentEstimator::Print(opt);
  }
  ClassDef(V0CentEstimator,1);
};
//____________________________________________________________________
/**
 * Centrality estimator using the V0 signal 
 */
struct RefMultEstimator : public FastNchCentEstimator
{
  /** Mode: Negative, use C side, positive use A side, otherwise sum */
  Double_t fEtaCut;
  /** 
   * Constructor 
   * 
   * @param etaCut Cut on eta 
   */
  RefMultEstimator(Double_t etaCut=0.8) 
    : FastNchCentEstimator(Form("RefMult%02dd%02d",
				Int_t(etaCut), Int_t(100*etaCut)%100)),
      fEtaCut(etaCut)
  {}
  /** 
   * Set-up this object.  Defines the internal histogram and add to
   * output
   * 
   * @param l Output list
   * @param tree Tree to add branch to 
   * @param sNN   Collision energy in GeV
   * @param tgtA  Target atomic weight 
   * @param projA Projectile atomic weight 
   */
  void Setup(TCollection* l, TTree* tree, UShort_t sNN,
	     Bool_t tgtA, Bool_t projA)
  {
    Bool_t  isAA  = (tgtA && projA);
    Bool_t  isPA  = (tgtA ^ projA); // XOR
    UInt_t  max   = (isAA ? 15000 : isPA ? 900 : 200);
    UInt_t  dBin  = (isAA ? 10    : isPA ?   1 :   1);
    Color_t color = kMagenta;
    fHistogram = new TH1D(Form("raw%s",GetName()),
			  Form("#it{N}_{ch} |#it{#eta}|<%5.2f distribution",
			       fEtaCut), max/dBin, 0, max);
    fHistogram->SetXTitle("#it{N}_{ch}");
    fHistogram->SetYTitle("Raw #it{P}(#it{N}_{ch})");
    fHistogram->SetDirectory(0);
    fHistogram->SetLineColor(color);
    fHistogram->SetFillColor(color);
    fHistogram->SetMarkerColor(color);
    fHistogram->SetMarkerStyle(20);
    fHistogram->SetFillStyle(3002);

    Fast1DCentEstimator::Setup(l, tree, sNN, tgtA, projA);
  }
  /** 
   * Whether we should accept a particle.  We accept a particle if it
   * falls within the acceptance of the V0.
   * 
   * @param p Particle to investigate 
   * 
   * @return true if to be used 
   */
  Bool_t Accept(const TParticle* p)
  {
    Double_t eta = Eta(p);
    if (!IsPrimary(p)) return false;
    if (TMath::Abs(eta) > fEtaCut) return false;
    return true;
  }
  /** 
   * Get the histogram to accumulate the observable in.  
   * 
   * @return Pointer to the histogram. 
   */
  virtual TH1* GetHistogram(TCollection* l)
  {
    return static_cast<TH1*>(l->FindObject(Form("raw%s",GetName())));
  }
  ClassDef(RefMultEstimator,1);
};
//____________________________________________________________________
/**
 * Centrality estimator using the V0 signal 
 */
struct ZNCentEstimator : public Fast1DCentEstimator
{
  // What to put in 
  /** Mode: Negative, use C side, positive use A side, otherwise sum */
  Short_t   fMode;
  Bool_t    fNeutrons;
  Bool_t    fSpectators;
  Bool_t    fPrimary;
  Double_t  fMinEta;
  Double_t  fMaxEta;
  Double_t  fMaxPhi;
  ULong64_t fNSpec;
  /** 
   * Constructor 
   * 
   * @param mode Mode: Negative, use C side, positive use A side, otherwise sum 
   * @param neutrons 
   * @param spectators 
   * @param primary 
   */
  ZNCentEstimator(Short_t mode=0,
		  Bool_t neutrons=true,
		  Bool_t spectators=false,
		  Bool_t primary=false)
    : Fast1DCentEstimator(Form("Z%c%c%c%c",
			       (neutrons   ? 'N' : 'P'),
			       (mode < 0   ? 'C' : (mode > 0 ? 'A' : 'M')),
			       (spectators ? 'S' : 'E'),
			       (primary    ? 'P' : 'A'))),
      fMode(mode),
      fNeutrons(neutrons),
      fSpectators(spectators),
      fPrimary(primary),
      fMinEta(0),
      fMaxEta(0),
      fMaxPhi(-1)
  {
    fFromTop = (fSpectators ? false : true);
    if (fNeutrons) {
      const Double_t zN = 1161.3; // distance to ZN cm
      const Double_t dN = 7;      // Area of ZN
      const Double_t rN = TMath::Sqrt(2*dN*dN); // Radius of ZN
      const Double_t tN = TMath::ATan2(rN,zN); // Largest angle ZN
      fMinEta           = -TMath::Log(TMath::Tan(tN/2));
      fMaxEta           = TMath::Infinity();
      fMaxPhi           = -1;
    }
    else {
      const Double_t  zP  = 1156.3; // distance to ZC cm
      const Double_t  wP  = 20.8; // Width of cal
      const Double_t  hP  = 12; // Width of cal
      const Double_t  oP  = 19;
      const Double_t  rO  = TMath::Sqrt(TMath::Power(oP+wP/2,2)+hP*hP/4);
      const Double_t  rI  = TMath::Sqrt(TMath::Power(oP-wP/2,2)+hP*hP/4);
      const Double_t  tO  = TMath::ATan2(rO, zP);
      const Double_t  tI  = TMath::ATan2(rI, zP);
      fMinEta             = -TMath::Log(TMath::Tan(tO/2));
      fMaxEta             = -TMath::Log(TMath::Tan(tI/2));
      fMaxPhi             = TMath::ATan2(hP/2,oP);
    }
  }
  /** 
   * Set-up this object.  Defines the internal histogram and add to
   * output
   * 
   * @param l Output list
   * @param tree Tree to add branch to 
   * @param sNN   Collision energy in GeV
   * @param tgtA  True if target is a nucleus 
   * @param projA True if projectile is a nucleus 
   */
  void Setup(TCollection* l, TTree* tree, UShort_t sNN,
	     Bool_t tgtA, Bool_t projA)
  {
    Bool_t  isAA  = (tgtA && projA);
    Bool_t  isPA  = (tgtA ^ projA); // XOR
    UInt_t  max   = (isAA ? 2000 : isPA ? 300 : 30);
    UInt_t  dBin  = (isAA ? 10   : isPA ?   1 :  1);
    Color_t color = (fMode < 0 ? kRed : fMode > 0 ? kBlue : kGreen)+2;

    TString nTxt;
    nTxt.Form("#it{N}_{%s%c} (%s)",
	      fSpectators ? "spec," : "", fNeutrons ? 'n' : 'p',
	      fPrimary ? "primary" : "all");
    fHistogram = new TH1D(Form("raw%s",GetName()),
			  Form("%s %s  distribution", GetName(), nTxt.Data()),
			  max/dBin, 0, (fMode == 0 ? 2 : 1)*max);
    fHistogram->SetXTitle(nTxt);
    fHistogram->SetYTitle(Form("Raw #it{P}(%s)", nTxt.Data()));
    fHistogram->SetDirectory(0);
    fHistogram->SetLineColor(color);
    fHistogram->SetFillColor(color);
    fHistogram->SetMarkerColor(color);
    fHistogram->SetMarkerStyle(20);
    fHistogram->SetFillStyle(3002);

    Fast1DCentEstimator::Setup(l, tree, sNN, tgtA, projA);
  }
  /** 
   * Process a single particle. 
   * 
   * @param p Particle 
   */
  virtual void Process(const TParticle* p)
  {
    // Check if we should count emitted 
    if (fSpectators) return;

    // Check if we should count only primaries
    if      (fPrimary && !IsPrimary(p)) return; 
    else if (!fPrimary && p->GetStatusCode()==4) return; 
    // if ((fFlags & kPrimary) != 0 && !IsPrimary(p)) return;

    // Check if this is a spectator 
    // Int_t status = p->GetStatusCode();
    // if ((fFlags & kSpectators) == 0 && (status==13 || status==14)) return;

    // Check particle type 
    Int_t aPdg = TMath::Abs(p->GetPdgCode());
    if (fNeutrons && aPdg != kNeutron)
      // Looking for neutrons 
      return;
    if (!fNeutrons && aPdg != kProton)
      // Looking for protons 
      return;

    // Check side
    Double_t eta = Eta(p);
    if (fMode < 0 && eta > 0) return; // Wrong side
    if (fMode > 0 && eta < 0) return; // Wrong side

    // Check acceptance 
    Double_t aeta = TMath::Abs(eta);
    if (aeta > fMaxEta || aeta < fMinEta) return; // Not acceptance

    if (fMaxPhi > 0) { 
      Double_t phi = TMath::Abs(Phi(p) - (eta < 0 ? 0 : TMath::Pi()));
      if (phi > fMaxPhi) return; // Not acceptance 
    }

    // Decrement
    fCache++;
    // Printf("%s: increment from %d (%f) -> %lld",
    //        GetName(), aPdg, eta, fCache);
  }
  void ProcessHeader(FastShortHeader& h)
  {
    // In Ap (EPOS) we have many spec in the target, meaning they will
    // be detected on the C side (p is the projectile, A is the target)
    // if (!fSpectators) return;
    fCache   = 0;
    fNSpec   = 0;
    Int_t nC = 0;
    Int_t nA = 0;
    if (fNeutrons) {
      nC = h.fNSpecNtgt;
      nA = h.fNSpecNproj;
    }
    else {
      nC = h.fNSpecPtgt;
      nA = h.fNSpecPproj;
    }
    if      (fMode < 0) fNSpec = nC;
    else if (fMode > 0) fNSpec = nA;
    else                fNSpec = nC + nA;
    // Printf("%s: Initial cache value: %d", GetName(), fNSpec);
  }
  virtual void PostEvent()
  {
    // Info(GetName(), "Nspec=%lld Nneu=%lld", fNSpec, fCache);
    if (fSpectators) fCache = fNSpec;
    else {
      if (fNSpec > fCache) {
	// Warning("PostEvent", "Nspec (%d) > Nneu (%d)", fNSpec, fCache);
	fNSpec = fCache;
      }
      if (!fPrimary) fCache -= fNSpec;
    }
    Fast1DCentEstimator::PostEvent();
  }
  /** 
   * Get the histogram to accumulate the observable in.  
   * 
   * @return Pointer to the histogram. 
   */
  virtual TH1* GetHistogram(TCollection* l)
  {
    return static_cast<TH1*>(l->FindObject(Form("raw%s",GetName())));
  }
  virtual void Print(Option_t* option="nah") const
  {
    TString opt(option); opt.ToLower();
    if (opt.Contains("n")) 
      Printf("ZN Estimator: %s (%screasing) mode=%d - %s %s",
	     GetName(),(fFromTop ? "de" : "in"), fMode,
	     (fSpectators ? "spectator" : "emitted"),
	     (fPrimary    ? "primary"   : "all"),
	     (fNeutrons   ? "neutrons"  : "protons"));
    if (opt.Contains("a"))
      Printf("  eta=[%f,%f], phi=%f", fMinEta, fMaxEta, fMaxPhi);
    opt.ReplaceAll("n", "");
    Fast1DCentEstimator::Print(opt);
  }
  ClassDef(ZNCentEstimator,1);
};
#endif
//
// EOF
//
