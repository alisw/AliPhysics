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
#include "FastCentHelper.C"
namespace dNdy {
  //====================================================================
  /** 
   * Base class for making @f$ 1/N dN_{ch}/dy@f$ 
   */
  struct Base : public FastAnalysis
  {
    TH1* fdNdy; //!
  
    Base(Bool_t verbose=false, Int_t monitor=0)
      : FastAnalysis(verbose,monitor), fdNdy(0)
    {}
    /** 
     * Static member function to create a histogram 
     * 
     * @return Newly created histogram
     */  
    static TH1D* CreatedNdy()
    {
      Double_t maxY   = 9; // 10;
      Double_t dY     = .25; 
      TH1D*    h      = new TH1D("dNdy", "Charged particle-rapidity density",
				 Int_t(2*maxY/dY+.5), -maxY, +maxY);
      h->Sumw2();
      h->SetXTitle("\\mathit{y}");
      h->SetYTitle("\\hbox{d}N_{\\hbox{ch}}/\\hbox{d}y");
      h->SetMarkerColor(kRed+2);
      h->SetMarkerStyle(20);
      h->SetDirectory(0);
    
      return h;
    }
    /** 
     * Called on each slave at start of processing
     */
    virtual void SlaveBegin(TTree*)
    {
      Info("SlaveBegin", "Making dN/dy histogram");
      fdNdy = CreatedNdy();
      fdNdy->SetMarkerColor(kBlack);
      fdNdy->SetMarkerStyle(21);
      fdNdy->SetTitle(GetName());
      fOutput->Add(fdNdy);
    }
    /** 
     * Process a particle.  
     * 
     * @param p Particle to process
     */    
    virtual Bool_t ProcessParticle(const TParticle* p)
    {
      Fill(p->Y());
      return true;
    }
    virtual void Fill(Double_t y) { fdNdy->Fill(y); }
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
      fdNdy = static_cast<TH1*>(GetOutputObject("dNdy", TH1::Class()));
      if (!fdNdy) {
	SetStatus(-1);
	Warning("Terminate", "No dN/deta histogram found");
	return;
      }
      Printf("A total of %ld events", fOK);
      TH1* copy = static_cast<TH1*>(fdNdy->Clone("before"));
      fOutput->Add(copy);
						
      fdNdy->Scale(1. / fOK, "width");
      fdNdy->Draw();
      SetStatus(0);
    }
    /** 
     * Get the list of monitored objects 
     * 
     * @return The list of monitored objects 
     */
    virtual TList* GetMonitorObjects()
    {
      TObject* m1 = new TNamed("dNdy", "");
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
     * @param filename File to open 
     * @param verbose  Whether to be verbose 
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
   * Processes INEL events and build the @f$ dN_{ch}/dy@f$ 
   * 
   */
  struct INEL : public Base
  {
    /** 
     * Constructor 
     * 
     * @param verbose  Whether to be verbose 
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
   * Processes INEL events and build the @f$ dN_{ch}/dy@f$ 
   * 
   */
  struct INELGt0 : public Base
  {
    /** 
     * Constructor 
     * 
     * @param verbose  Whether to be verbose 
     */  
    INELGt0(Bool_t verbose=false, Int_t monitor=0)
      : Base(verbose,monitor)
    {}
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
     */
    Cent(const char* method="V0M", Bool_t verbose=true, Int_t monitor=0)
      : Base(verbose,monitor),
	fHelper(method),
	fMinEvents(100),
	fCentBin(0)	
    {
      fCentMethod = fHelper.fCentMeth;
    }
    /** 
     * Interface for helper 
     * 
     * @return Newly allocated histogram 
     */
    virtual TH1* CreateSingle() { return Base::CreatedNdy(); }
    /** 
     * Called on each slave at start of processing
     * 
     * @param t Ignored
     */
    virtual void SlaveBegin(TTree* t)
    {
      Base::SlaveBegin(t);
      // We use the parent fdNdy histogram as a cache 
      fdNdy->SetName("cache");
      fdNdy->SetTitle("Cache");

      fHelper.CreateHistos(fOutput, Base::CreatedNdy);
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
      fdNdy->Reset(); // Clear the cache
    }
    /** 
     * Process the header.  Accepts events in range 
     * 
     * @return true 
     */
    virtual Bool_t ProcessHeader()
    {
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
	Double_t n0   = fdNdy->GetBinContent(fdNdy->GetNbinsX()/2);
	Double_t d0   = fdNdy->GetXaxis()->GetBinWidth(fdNdy->GetNbinsX()/2);
	Double_t y0   = (n0 / d0);
	Printf("Centrality %6.2f-%6.2f (bin %4d) "
	       "Nch=%8.1f/%4.2f=%8.1f out=%s (%d)",
	       fHelper.fCentAxis->GetBinLowEdge(fCentBin),
	       fHelper.fCentAxis->GetBinUpEdge(fCentBin),
	       fCentBin,
	       n0, d0, y0,
	       out->GetName(),
	       Int_t(out->GetBinContent(0)));
      }

      // Make sure we have nothing in the underflow bin. 
      fdNdy->SetBinContent(0,0);

      // Add our cache to the appropriate bin 
      out->Add(fdNdy);
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

      if (!fHelper.Finalize(fOutput, fMinEvents))
	SetStatus(-1);
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
    ClassDef(Cent,1);
  };
  
  //====================================================================
  /**
   * Processes events and build the @f$ dN_{ch}/dy@f$ for each
   * bin in reference multiplicity.  The reference multiplicity is the
   * number of charged particles with @f$|\eta|\le0.8@f$
   */
  struct Mult : public Cent
  {
    /** 
     * Constructor. 
     * 
     * @param filename File to open
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
struct dNdyMaker : public FastAnalysis::Maker
{
  /** 
   * CTOR 
   */
  dNdyMaker() : FastAnalysis::Maker("dNdy") {}
  
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
		      TString&       uopt)
  {
    TString t(subtype);
    if      (t.EqualTo("INEL"))    return new dNdy::INEL(verbose,monitor);
    else if (t.EqualTo("NSD"))     return new dNdy::NSD(verbose,monitor);
    else if (t.EqualTo("INELGt0")) return new dNdy::INELGt0(verbose,monitor);
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
	Printf("Warning: dNdyMaker::Make: Unknown estimator: %s",
	       w.Data());
	return 0;
      }
      if (t.BeginsWith("MULT"))
	return new dNdy::Mult(w, verbose, monitor);
      else
	return new dNdy::Cent(w, verbose, monitor);
    }
    Printf("Error: dNdyMaker::Run: Invalid spec: %s", t.Data());
    return 0;
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
dNdyMaker* _dNdyMaker = new dNdyMaker;

#ifdef __MAKECINT__
#pragma link C++ nestedclasses;
#pragma link C++ namespace dNdy;
#endif

//
// EOF
//
