#ifndef DNDETAANALYSIS_C
#define DNDETAANALYSIS_C
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
class TMap;
#endif
#include "FastAnalysis.C"
#include "FastCentHelper.C"
namespace dNdeta {
  //====================================================================
  /** 
   * Base class for making @f$ 1/N dN_{ch}/d\eta@f$ 
   */
  struct Base : public FastAnalysis
  {
    TH1* fdNdeta; //!
  
    Base(Bool_t verbose=false, Int_t monitor=0)
      : FastAnalysis(verbose,monitor), fdNdeta(0)
    {}
    /** 
     * Static member function to create a histogram 
     * 
     * @return Newly created histogram
     */  
    static TH1D* CreatedNdeta()
    {
      Double_t maxEta = 9; // 10;
      Int_t    nEta   = Int_t(200 * maxEta/5+.5);
      TH1D*    eta    = new TH1D("dNdeta",
				 "Charged particle pseudo-rapidity density",
				 nEta, -maxEta, +maxEta);
      eta->Sumw2();
      eta->SetXTitle("#it{#eta}");
      eta->SetYTitle("1/N d#it{N}_{ch}/d#it{#eta}");
      eta->SetMarkerColor(kRed+2);
      eta->SetMarkerStyle(20);
      eta->SetDirectory(0);
    
      return eta;
    }
    static TObject* CreateOutput()
    {
      return CreatedNdeta();
    }
    /** 
     * Called on each slave at start of processing
     */
    virtual void SlaveBegin(TTree*)
    {
      Info("SlaveBegin", "Making dN/deta histogram");
      fdNdeta = CreatedNdeta();
      fdNdeta->SetMarkerColor(kBlack);
      fdNdeta->SetMarkerStyle(21);
      fdNdeta->SetTitle(GetName());
      fOutput->Add(fdNdeta);
    }
    /** 
     * Process a particle.  
     * 
     * @param p Particle to process
     */    
    virtual Bool_t ProcessParticle(const TParticle* p)
    {
      Double_t   pT    = p->Pt();
      Double_t   pZ    = p->Pz();
      Double_t   theta = TMath::ATan2(pT, pZ);
      Double_t   eta   = (pT < 1e-10 ? 1024 : -TMath::Log(TMath::Tan(theta/2)));
      if (TMath::Abs(eta) > 1000) return false;

      Fill(eta);
      return true;
    }
    virtual void Fill(Double_t eta) { fdNdeta->Fill(eta); }
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
      fdNdeta = static_cast<TH1*>(GetOutputObject("dNdeta", TH1::Class()));
      if (!fdNdeta) {
	SetStatus(-1);
	Warning("Terminate", "No dN/deta histogram found");
	return;
      }
      Printf("A total of %ld events", fOK);
      TH1* copy = static_cast<TH1*>(fdNdeta->Clone("before"));
      fOutput->Add(copy);
						
      fdNdeta->Scale(1. / fOK, "width");
      fdNdeta->Draw();
      SetStatus(0);
    }
    /** 
     * Get the list of monitored objects 
     * 
     * @return The list of monitored objects 
     */
    virtual TList* GetMonitorObjects()
    {
      TObject* m1 = new TNamed("dNdeta", "");
      m1->SetUniqueID(0x8); // Scale 

      TList* ret = new TList;
      ret->Add(m1);
    
      return ret;
    }
    ClassDef(Base,1);
  };

  //====================================================================
  /** 
   * Select NSD events and builds the @f$ 1/N dN_{ch}/d\eta@f$ 
   * 
   */
  struct NSD : public Base
  {
    /** 
     * Constructor 
     * 
     * @param verbose  Whether to be verbose 
     * @param monitor  Frequency 
     */
    NSD(Bool_t verbose=false, Int_t monitor=0)
      : Base(verbose, monitor)
    {}
    /** 
     * Process the header.  Return 
     * 
     * @return True in case the event is flagged as NSD, false
     * otherwise.
     */
    virtual Bool_t ProcessHeader()
    {
      if (fHeader->fType & 0x1) return false;
      return true;
    }
    ClassDef(NSD,1);
  };

  //====================================================================
  /**
   * Processes INEL events and build the @f$ 1/N dN_{ch}/d\eta@f$ 
   * 
   */
  struct INEL : public Base
  {
    /** 
     * Constructor 
     * 
     * @param verbose  Whether to be verbose 
     * @param monitor  Frequency 
     */  
    INEL(Bool_t verbose=false, Int_t monitor=0)
      : Base(verbose,monitor)
    {}
    /** 
     * Process the header.  
     * 
     * @return Always true - by definition INEL is all events. 
     */
    virtual Bool_t ProcessHeader() { return true; }
    ClassDef(INEL,1);
  };

  //====================================================================
  /**
   * Processes INEL events and build the @f$ 1/N dN_{ch}/d\eta@f$ 
   * 
   */
  struct INELGt0 : public Base
  {
    /** 
     * Constructor 
     * 
     * @param verbose  Whether to be verbose 
     * @param monitor  Frequency 
     */  
    INELGt0(Bool_t verbose=false, Int_t monitor=0)
      : Base(verbose,monitor)
    {
      SetTrigger("INEL>0");
    }
    /** 
     * Process the header.  
     * 
     * @return Always true - by definition INEL is all events. 
     */
    virtual Bool_t ProcessHeader()
    {
      if (fHeader->fType & 0x4) return true;
      return false;
    }
    ClassDef(INELGt0,1);
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
     * Interface for helper 
     * 
     * @return Newly allocated histogram 
     */
    virtual TH1* CreateSingle() { return Base::CreatedNdeta(); }
    /** 
     * Called on each slave at start of processing
     * 
     * @param t Ignored
     */
    virtual void SlaveBegin(TTree* t)
    {
      Base::SlaveBegin(t);
      // We use the parent fdNdeta histogram as a cache 
      fdNdeta->SetName("cache");
      fdNdeta->SetTitle("Cache");

      fHelper.CreateHistos(fOutput, Base::CreateOutput);
      fOutput->ls();
    }
    /** 
     * Set-up centrality estimator
     *
     * @return true on success 
     */
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
      fdNdeta->Reset(); // Clear the cache
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
      TH1*     out  = static_cast<TH1*>(fHelper.fCentList->At(fCentBin-1));
      // If we still have no histogram, return immediately 
      if (!out) return;

      // Use parent function to fill cache 
      Base::ProcessParticles();
      if (fVerbose) {
	Double_t n0   = fdNdeta->GetBinContent(fdNdeta->GetNbinsX()/2);
	Double_t d0   = fdNdeta->GetXaxis()->GetBinWidth(fdNdeta->GetNbinsX()
							 /2);
	Double_t eta0 = (n0 / d0);
	Printf("Centrality %6.2f-%6.2f (bin %4d) "
	       "Nch=%8.1f/%4.2f=%8.1f out=%s (%d)",
	       fHelper.fCentAxis->GetBinLowEdge(fCentBin),
	       fHelper.fCentAxis->GetBinUpEdge(fCentBin),
	       fCentBin,
	       n0, d0, eta0,
	       out->GetName(),
	       Int_t(out->GetBinContent(0)));
      }

      // Make sure we have nothing in the underflow bin. 
      fdNdeta->SetBinContent(0,0);

      // Add our cache to the appropriate bin 
      out->Add(fdNdeta);
      // Increment underflow bin for the event count 
      out->SetBinContent(0, out->GetBinContent(0)+1);
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

      if (!fHelper.Finalize(fOutput, fMinEvents, Cent::Normalize))
	SetStatus(-1);
    }
    static TH1* Normalize(TObject* o, Int_t n)
    {
      if (!o->IsA()->InheritsFrom(TH1::Class())) return 0;
      TH1*  h = static_cast<TH1*>(o);
      Int_t m = h->GetBinContent(0);
      if (m != n) 
	::Warning("Normalize", "Different event count here:%d helper:%d",m,n);
      h->SetBinContent(0,0);
      h->Scale(1. / n, "width");
      return h;
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
    
      m3->SetUniqueID(0x8); // Scale 
      TList* ret = new TList;
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
  
  //====================================================================
  /**
   * Processes events and build the @f$ 1/N dN_{ch}/d\eta@f$ for each
   * bin in reference multiplicity.  The reference multiplicity is the
   * number of charged particles with @f$|\eta|\le0.8@f$
   */
  struct Mult : public Cent
  {
    /** 
     * Constructor. 
     * 
     * @param method   Method 
     * @param monitor  Frequency 
     * @param verbose  Whether to verbose
     */
    Mult(const char* method="RefMult00d80",
	 Bool_t verbose=false,
	 Int_t monitor=0)
      : Cent(method, verbose,monitor)
    {
      //              +1 +2 +3 +3  +5, +5, +5, +5,+10,+10,+10,+10,+10,+10,+10
      Double_t bins[]={ 0, 3, 6, 9, 14, 19, 24, 29, 39, 49, 59, 69, 79, 89, 99};
      fHelper.SetCentAxis(14, bins);
    }
    ClassDef(Mult,1);
  };
} // End of namespace 
//====================================================================
/*
 * The function to make our analyser 
 */
struct dNdetaMaker : public FastAnalysis::Maker
{
  /** 
   * CTOR 
   */
  dNdetaMaker() : FastAnalysis::Maker("dNdeta") {}
  
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
    FastAnalysis* ret = 0;
    TString t(subtype);
    if      (t.EqualTo("INEL"))    ret = new dNdeta::INEL   (verbose,monitor);
    else if (t.EqualTo("NSD"))     ret = new dNdeta::NSD    (verbose,monitor);
    else if (t.EqualTo("INELGt0")) ret = new dNdeta::INELGt0(verbose,monitor);
    else if (t.EqualTo("V0AND"))   ret = new dNdeta::V0AND  (verbose,monitor);
    else if (t.BeginsWith("MULT") || t.BeginsWith("CENT")) {
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
	Printf("Warning: dNdetaMaker::Make: Unknown estimator: %s",
	       w.Data());
	return 0;
      }
      if (t.BeginsWith("MULT"))
	ret = new dNdeta::Mult(w, verbose, monitor);
      else
	ret = new dNdeta::Cent(w, verbose, monitor);
    }
    
    if (ret) {
      TPair* tp = static_cast<TPair*>(uopt.FindObject("trig"));
      if( tp) ret->SetTrigger(tp->Value()->GetName());
    }
    else 
      Printf("Error: dNdetaMaker::Make: Invalid spec: %s", t.Data());
    return ret;
  }
  /** 
   * Show list of possible sub-types
   * 
   */
  void List() const
  {
    Printf(" INEL            - inelastic");
    Printf(" INELGt0         - inelastic with at least 1 particle in |eta|<1");
    Printf(" NSD             - Non-single diffractive");
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
  const char* Script() const {   return __FILE__;  }
};

// ------------------------------------------------------------------
// Create instance of maker
dNdetaMaker* _dNdetaMaker = new dNdetaMaker;

#ifdef __MAKECINT__
#pragma link C++ nestedclasses;
#pragma link C++ namespace dNdeta;
#endif
#endif
//
// EOF
//
