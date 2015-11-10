#ifndef FASTSIM_H
#define FASTSIM_H
#include <TSelector.h>
#include <TQObject.h>
#ifndef __CINT__
# include "AliGenerator.h"
# include "AliRunLoader.h"
# include "AliStack.h"
# include "AliHeader.h"
# include "AliGenEventHeader.h"
# include "AliRun.h"
# include "AliCollisionGeometry.h"
# include "AliGenPythiaEventHeader.h"
# include "AliGenDPMjetEventHeader.h"
# include "AliGenGeVSimEventHeader.h"
# include "AliGenHerwigEventHeader.h"
# include "AliGenCocktailEventHeader.h"
# include "AliPDG.h"
# include <TROOT.h>
# include <TString.h>
# include <TMath.h>
# include <TParticle.h>
# include <TH1.h>
# include <TTree.h>
# include <TClonesArray.h>
# include <TList.h>
# include <TProof.h>
# include <TChain.h>
# include <TParticlePDG.h>
# include <TStopwatch.h>
# include <TFile.h>
# include <TProofOutputFile.h>
# include <TCanvas.h>
# include <TTimer.h>
# include <TRandom.h>
# include <TUrl.h>
# include <TMacro.h>
# include <TSystemDirectory.h>
# include <TFileCollection.h>
# include <TPRegexp.h>
# include <fstream>
# include <TLeaf.h>
# include <algorithm>
#else
class AliGenerator;
class AliRunLoader;
class AliStack;
class AliHeader;
class AliGenEventHeader;
class TH1;
class TTree;
class TClonesArray;
class TBrowser;
class TList;
class TFile;
class TProofOutputFile;
class TCanvas;
class TVirtualPad;
class TTimer;
class TUrl;
class TAxis;
class TParticle;
class TMacro;
class TLeaf;
#endif

enum {
  kNeutron = 2112,
  kProton = 2212
};

  
/** To get DPMJEt common block */
typedef struct {
   Double_t	rproj;
   Double_t	rtarg;
   Double_t	bimpac;
   Int_t	nwtsam;
   Int_t	nwasam;
   Int_t	nwbsam;
   Int_t	nwtacc;
   Int_t	nwtaac;
   Int_t	nwtbac;
   Int_t        ncp;
   Int_t        nct;
} DtglcpCommon;
DtglcpCommon* _dtglcp = 0;

//====================================================================
/** 
 * Header structure 
 */
struct ShortHeader {
  UInt_t   fRunNo;
  UInt_t   fEventId;
  UInt_t   fNtgt;
  UInt_t   fNproj;
  UInt_t   fNbin;
  UInt_t   fType;
  Double_t fIpX;
  Double_t fIpY;
  Double_t fIpZ;
  Double_t fB;
  Double_t fC;
  Double_t fPhiR;
  UInt_t   fNSpecNproj;  // # of spectator neutrons in projectile
  UInt_t   fNSpecNtgt;   // # of spectator neutrons in target 
  UInt_t   fNSpecPproj;  // # of spectator protons in projectile
  UInt_t   fNSpecPtgt;   // # of spectator protons in target
  
  void Print()
  {
    Printf(" Run #/Event:          %9d/%9d", fRunNo, fEventId);
    Printf(" Participants/binary:  %4d/%4d/%3d", fNtgt, fNproj, fNbin);
    Printf(" Event type:           %7s%12s",(fType==1?"Non":
					     fType==2?"Single":
					     "Double"), "-diffractive");
    Printf(" IP:                   (%-5.1f,%-5.1f,%-5.1f)",fIpX,fIpY,fIpZ);
    Printf(" Impact par./cent.:    (%13f/%-3d)", fB, Int_t(fC));
    Printf(" Reaction plane:       %19f", fPhiR);
    Printf(" Specs (Nt,Np,Pt,Pp):  %4d/%4d/%4d/%4d",
	   fNSpecNtgt, fNSpecNproj, fNSpecPtgt, fNSpecPproj);
  }
  void Reset(UInt_t runNo, UInt_t eventNo)
  {
    fRunNo      = runNo;
    fEventId    = eventNo;
    fIpX        = 1024;
    fIpY        = 1024;
    fIpZ        = 1024;
    fNtgt       = -1;
    fNproj      = -1;
    fNbin       = -1;
    fPhiR       = -1;
    fB          = -1;
    fC          = -1;
    fNSpecNtgt  = -1;
    fNSpecNproj = -1;
    fNSpecPtgt  = -1;
    fNSpecPproj = -1;
    fEG         = kUnknown;
  }

  enum {
    kUnknown = 0,
    kPythia,
    kHijing,
    kDPMJet,
    kEPOS,
  } fEG;
};

//====================================================================
/** 
 * Monitor output objects
 */
struct FastMonitor : public TObject, public TQObject 
{
  /** 
   * Constructor 
   * 
   * 
   * @return 
   */
  FastMonitor(TSelector* s=0)
    : fName("FastMonitor"),
      fCanvas(0),
      fSelector(s)
  {
    if (gROOT->IsBatch()) {
      Warning("FastMonitor", "Batch processing, no monitoring");
      return;
    }

    if (gProof) {
      fName = gProof->GetSessionTag();
      gDirectory->Add(this);
      Bool_t ret = gProof->Connect("Feedback(TList *objs)",
				   "FastMonitor", this, 
				   "Feedback(TList *objs)");
      if (!ret) {
	Warning("FastMonitor", "Failed to connect to Proof");
	return;
      }
    }
    else if (!s) return;
    
    fCanvas = new TCanvas(fName, Form("Monitor %s", fName.Data()), 1000, 800);
    fCanvas->SetFillColor(0);
    fCanvas->SetFillStyle(0);
    fCanvas->SetTopMargin(0.01);
    fCanvas->SetRightMargin(0.01);

    fCanvas->Divide(3,2);
    RegisterDraw(1, "histograms/type",            "", 0);
    RegisterDraw(2, "histograms/b",               "", 0);
    RegisterDraw(3, "histograms/cent",            "", 0);
    RegisterDraw(4, "histograms/dNdeta",          "", 0x8);
    RegisterDraw(5, "estimators/rawV0M",          "", 0x2);
    RegisterDraw(6, "estimators/rawRefMult00d80", "", 0x2);
  }
  /** 
   * Register a draw of a an object 
   * 
   * @param i      Pad number 
   * @param name   Name of object 
   * @param option Drawing option
   * @param flags  Flags 
   *
   *  - 0x1   Log(x)
   *  - 0x2   Log(y)
   *  - 0x4   Log(z)
   *  - 0x8   Scale to events and bin width 
   */
  void RegisterDraw(Int_t i,
		    const char* name,
		    const char* option,
		    UShort_t    flags=0)
  {
    TVirtualPad* p = fCanvas->GetPad(i);
    if (!p) {
      Warning("RegisterDraw", "Not enough sub-pads (%d)", i);
      return;
    }
    p->SetFillColor(0);
    p->SetFillStyle(0);
    p->SetTopMargin(0.01);
    p->SetRightMargin(0.01);
    p->SetName(Form("p_%s", name));
    p->SetTitle(option);
    if (flags & 0x1) p->SetLogx();
    if (flags & 0x2) p->SetLogy();
    if (flags & 0x4) p->SetLogz();
    if (flags & 0x8) p->SetBit(BIT(15));
  }
  /** 
   * Desctructor 
   */
  virtual ~FastMonitor() 
  {
    if (!gProof) return;
    gProof->Disconnect("Feedback(TList *objs)",this, 
		       "Feedback(TList* objs)");
  }
  /** 
   * Set name of this object 
   * 
   * @param name Name 
   */
  void SetName(const char* name) { fName = name; }
  /** 
   * Get the name of this object 
   * 
   * @return Name 
   */
  const char* GetName() const { return fName.Data(); }
  /** 
   * Find pad corresponding to an object
   * 
   * @param name Name of object 
   * 
   * @return Pointer to pad or null
   */
  TVirtualPad* FindPad(const TString& name)
  {
    TVirtualPad* p = 0;
    Int_t        i = 1;
    TString      t = Form("p_%s", name.Data());
    while ((p = fCanvas->GetPad(i))) {
      if (t.EqualTo(p->GetName())) return p;
      i++;
    }
    return 0;
  }
  /** 
   * Find an object in the list @a l which corresponds to a registered
   * pad.
   * 
   * @param padName Pad's name 
   * @param l       Input collection 
   * 
   * @return The found object, or null
   */
  TObject* FindPadObject(const Char_t* padName, TCollection* l)
  {
    TString path(padName);
    path.Remove(0,2);
    if (path.Index("/") == kNPOS)
      return l->FindObject(path);
    TObjArray*   tokens  = path.Tokenize("/");
    Int_t        nTokens = tokens->GetEntriesFast();
    TCollection* current = l;
    TObject*     ret     = 0;
    for (Int_t i = 0; i < nTokens; i++) {
      TObject* o = current->FindObject(tokens->At(i)->GetName());
      if (!o) break;
      if (i == nTokens - 1) {
	ret = o;
	break;
      }
      if (!o->IsA()->InheritsFrom(TCollection::Class())) {
	Warning("FindPadObject", "Path object %s of %s is not a collection "
		"but a %s", o->GetName(), path.Data(), o->ClassName());
	break;
      }
      current = static_cast<TCollection*>(o);
    }
    delete tokens;
    // if (!ret) l->ls();
    return ret;
  }
    
  /** 
   * Called when we get notified of 
   * 
   * @param objs List of monitored objects
   */
  void Feedback(TList* objs)
  {
    // Info("FeedBack", "List is %p", objs);
    // if (objs) objs->ls();
    if (!fCanvas) return;

    // objs->ls();
    TList* l = static_cast<TList*>(objs->FindObject("list"));
    if (!l) {
      Warning("Feedback", "No list");
      return;
    }
    TList* hs = static_cast<TList*>(l->FindObject("histograms"));
    Int_t nEvents = 1;
    TObject* oIpz = hs->FindObject("ipZ");
    if (oIpz && oIpz->IsA()->InheritsFrom(TH1::Class())) 
      nEvents = static_cast<TH1*>(oIpz)->GetEntries();
    else 
      Warning("Feedback", "Histogram ipZ not found");

    Int_t        iPad = 1;
    TVirtualPad* p    = 0;
    while ((p = fCanvas->GetPad(iPad))) {
      TObject* o = FindPadObject(p->GetName(), l);
      if (!o) {
	Warning("Feedback", "Object correspondig to pad %s (%d) not found",
		p->GetName(), iPad);
	iPad++;
      }
      p->cd();
      if (o->IsA()->InheritsFrom(TH1::Class())) {
	TH1* h = static_cast<TH1*>(o);
	TH1* c = h->DrawCopy(p->GetTitle());
	c->SetDirectory(0);
	c->SetBit(TObject::kCanDelete);
	if (p->TestBit(BIT(15))) {
	  // Info("Feedback", "Scaling %s by 1./%d and width",
	  //      c->GetName(), nEvents);
	  c->Scale(1./nEvents, "width");
	}
      }
      else {
	TObject* c = o->DrawClone(p->GetTitle());
	c->SetBit(TObject::kCanDelete);
      }
      p->Modified();
      iPad++;
    }
    fCanvas->Modified();
    fCanvas->Update();
    fCanvas->cd();
  }
  /** 
   * Function to handle connect signals 
   * 
   */
  void Handle()
  {
    HandleTimer(0);
  }
  /**
   * Function to handle timer events 
   */
  Bool_t HandleTimer(TTimer*)
  {
    // Info("HandleTimer", "Selector=%p", fSelector);
    if (!fSelector) return false;
    Feedback(fSelector->GetOutputList());
    return true;
  }
  /** Our name */
  TString fName;
  /** Our canvas */
  TCanvas* fCanvas;
  /** Possibly link to selector */
  TSelector* fSelector;
  ClassDef(FastMonitor,1);
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
   * @param tgtA  True if target is a nucleus 
   * @param projA True if projectile is a nucleus 
   */
  virtual void Setup(TCollection* out, TTree* tree,
		     Bool_t tgtA, Bool_t projA) = 0;
  /** 
   * Called before the start of an event 
   * 
   */
  virtual void PreEvent() {}
  virtual void ProcessHeader(ShortHeader&) {}
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
  /** 
   * Set-up this object.  Defines the internal histogram and add to
   * output
   * 
   * @param l Output list
   * @param tree Tree to add branch to 
   * @param tgtA  True if target is a nucleus 
   * @param projA True if projectile is a nucleus 
   */
  void Setup(TCollection* l, TTree* tree,
	     Bool_t, Bool_t)
  {
    if (fHistogram && l) l->Add(fHistogram);
    if (tree) tree->Branch(GetName(), &fCache, "value/l");
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
   * @param tgtA  True if target is a nucleus 
   * @param projA True if projectile is a nucleus 
   */
  void Setup(TCollection* l, TTree* tree,
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

    Fast1DCentEstimator::Setup(l, tree, tgtA, projA);
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
   * @param mode Mode: Negative, use C side, positive use A side, otherwise sum 
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
   */
  void Setup(TCollection* l, TTree* tree,
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

    Fast1DCentEstimator::Setup(l, tree, tgtA, projA);
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
   * @param tgtA  True if target is a nucleus 
   * @param projA True if projectile is a nucleus 
   */
  void Setup(TCollection* l, TTree* tree,
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

    Fast1DCentEstimator::Setup(l, tree, tgtA, projA);
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
  void ProcessHeader(ShortHeader& h)
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
    Info(GetName(), "Nspec=%lld Nneu=%lld", fNSpec, fCache);
    if (fSpectators) fCache = fNSpec;
    else {
      if (fNSpec > fCache) {
	Warning("PostEvent", "Nspec (%d) > Nneu (%d)",
		fNSpec, fCache);
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
  

//====================================================================
/** 
 * Run a event generator simulation 
 */
struct FastSim : public TSelector
{
  /** 
   * Constructor 
   * 
   * @param eg     Event generator 
   * @param runNo  Run number to simulate 
   * @param bMin   Lease impact parameter 
   * @param bMax   Largest impact parameter 
   */
  FastSim(const char* eg="",
	  ULong_t runNo=0,
	  Double_t bMin=0,
	  Double_t bMax=20,
	  Long64_t nEvents=0)
    : TSelector(),
      fEGName(eg),
      fRunNo(runNo),
      fBMin(bMin),
      fBMax(bMax),
      fGRP(0),
      fOverrides(0),
      fNEvents(nEvents),
      fIsTgtA(false),
      fIsProjA(false),
      fGenerator(0),
      fRunLoader(0),
      fStack(0),
      fHeader(0),
      fTree(0),
      fParticles(0),
      fList(0),
      fHEta(0),
      fHIpz(0),
      fHType(0),
      fHCent(0),
      fHB(0),
      fHPhiR(0),
      fHTime(0),
      fCentEstimators(0),
      fProofFile(0),
      fAliceFile(0),
      fKineFile(0),
      fFile(0),
      fFileName(""),
      fVerbose(true)
  {}
  const char* FileName() const
  {
    static TString fn;
    if (fn.IsNull()) {
      if (!fFileName.IsNull())  fn = fFileName;
      else {
	TString egName = (fGenerator ? fGenerator->GetName() : "");
	if (egName.IsNull()) egName = fEGName;
	
	fn = Form("%s_%09d", egName.Data(), fRunNo);
	if (fGenerator) {
	  TString tgt, proj;
	  Int_t   tgtA, tgtZ, projA, projZ;
	  fGenerator->GetTarget(tgt, tgtA, tgtZ);
	  fGenerator->GetProjectile(proj, projA, projZ);
	  fn.Append(Form("_%s%s", tgt.Data(), proj.Data()));
	  fn.Append(Form("_%05d", Int_t(fGenerator->GetEnergyCMS())));
	}
	else {
	  Long_t en = gROOT->ProcessLine("grp->energy");
	  // Bool_t a1 = gROOT->ProcessLine("(GRPData*)%p->beam1.IsA()",fGRP);
	  // Bool_t a2 = gROOT->ProcessLine("(GRPData*)%p->beam2.IsA()",fGRP);
	  fn.Append(Form("_%s%s", (fIsTgtA ? "A" : "p"),
			 (fIsProjA ? "A" : "p")));
	  fn.Append(Form("_%05ld", (fGRP ? en : 0)));
	}

	if (fNEvents > 0) {
	  if (fNEvents >= 1000000)
	    fn.Append(Form("_%lldM", fNEvents/1000000));
	  else if (fNEvents >= 1000)
	    fn.Append(Form("_%lldk", fNEvents/1000));
	  else
	    fn.Append(Form("_%lld", fNEvents));
	}
	fn.Append(".root");
	fFileName = fn;
      }
    }
    return fn.Data();
  }
  /** 
   * Get the name of the selector 
   * 
   * @return Name of selector 
   */
  const char* GetName() const { return "FastSim"; }
  /** 
   * Get the title of the selector 
   * 
   * @return The title 
   */
  const char* GetTitle() const { return "ALICE Event Generator simulation"; }
  /** 
   * Create our outputs 
   * 
   * @return true on success 
   */
  Bool_t SetupOutput()
  {
    Printf("=== Output =============================\n"
	   " File name: %s\n"
	   "========================================", FileName());

    if (fVerbose) Info("SetupOutput", "First the file");
    Bool_t isProof = false;
    if (fInput && fInput->FindObject("PROOF_Ordinal"))
      isProof = true;
    if (isProof) {
      Info("SetupOutput", "Making Proof File");
      fProofFile = new TProofOutputFile(FileName(), "M");
      // TProofOutputFile::kMerge,
      // TProofOutputFile::kRemote);
      fFile = fProofFile->OpenFile("RECREATE");
    }
    else
      fFile = TFile::Open(FileName(), "RECREATE");

    if (fVerbose) Info("SetupOutput", "Making our tree");
    fTree      = new TTree("T", (fGenerator ?
				 fGenerator->GetTitle() :
				 "T"));
    fParticles = new TClonesArray("TParticle");
    fTree->Branch("header", &fShortHead,
		  "run/i:event:ntgt:nproj:nbin:type:"
		  "ipx/D:ipy:ipz:b:c:phir:"
		  "nsnp/i:nsnt:nspp:nspt");
    fTree->Branch("particles", &fParticles);
    fTree->SetAutoSave(500); // Save every on every 100 events 
    fTree->SetDirectory(fFile);
    fTree->SetAlias("primary", "(particles.fBits&(1<<14))");
    fTree->SetAlias("weak",    "(particles.fBits&(1<<15))");
    fTree->SetAlias("charged", "(particles.fBits&(1<<16))");
    fTree->SetAlias("pt",      "(sqrt(pow(particles.fPx,2)+"
		    /*       */"pow(particles.fPy,2)))");
    fTree->SetAlias("eta",     "(particles.Pt()<1e-100?"
		    "(particles.fPz>0?100:-100):"
		    "-log(tan(atan2(particles.Pt(),particles.fPz)/2)))");
    fTree->SetAlias("good",    "(primary&&charged&&abs(eta)<1000)");
    fTree->SetAlias("sd",      "(header.fType & 0x1)");
    fTree->SetAlias("dd",      "(header.fType & 0x2)");
    fTree->SetAlias("pion",    "(abs(particles.fPdgCode)==211)");
    fTree->SetAlias("kaon",    "(abs(particles.fPdgCode)==321)");
    fTree->SetAlias("proton",  "(abs(particles.fPdgCode)==2212)");
    fTree->SetAlias("electron","(abs(particles.fPdgCode)==11)");
    fTree->SetAlias("other",   "(!pion&&!kaon&&!proton&&!electron)");
    fTree->SetAlias("beta",    "(particles.P()/particle.Energy())");
    fTree->SetAlias("gamma",   "(1./sqrt(1-beta*beta))");
    fTree->SetAlias("npart",   "(header.ntgt+header.nproj)");

    if (fVerbose) Info("SetupOutput", "Making histograms");
    Double_t maxEta = 10;
    Double_t dEta   = 10./200;
    fHEta = new TH1D("dNdeta", "Charged particle pseudo-rapidity density",
		     Int_t(2*maxEta/dEta+.5), -maxEta, +maxEta);
    fHEta->Sumw2();
    fHEta->SetXTitle("#it{#eta}");
    fHEta->SetYTitle("1/N d#it{N}_{ch}/d#it{#eta}");
    fHEta->SetMarkerColor(kRed+2);
    fHEta->SetMarkerStyle(20);
    fHEta->SetDirectory(0);
    
    fHIpz = new TH1D("ipZ", "Z-coordinate of interaction point",
		     10, -10, 10);
    fHIpz->SetMarkerColor(kGreen+2);
    fHIpz->SetFillColor(kGreen+2);
    fHIpz->SetFillStyle(3001);
    fHIpz->SetXTitle("IP_{#it{z}} [cm]");
    fHIpz->SetDirectory(0);

    fHType = new TH1D("type", "Diffractive", 3, .5, 3.5);
    fHType->SetMarkerColor(kOrange+2);
    fHType->SetFillColor(kOrange+2);
    fHType->SetFillStyle(3001);
    fHType->SetDirectory(0);
    fHType->GetXaxis()->SetBinLabel(1, "Non");
    fHType->GetXaxis()->SetBinLabel(2, "Single");
    fHType->GetXaxis()->SetBinLabel(3, "Double");

    fHCent = new TH1D("cent", "Centrality", 20, 0, 100);
    fHCent->SetMarkerColor(kPink+2);
    fHCent->SetFillColor(kPink+2);
    fHCent->SetFillStyle(3001);
    fHCent->SetDirectory(0);
    fHCent->SetXTitle("Centrality [%]");
    
    fHB = new TH1D("b", "Impact parameter", 20, 0, 20);
    fHB->SetMarkerColor(kYellow+2);
    fHB->SetFillColor(kYellow+2);
    fHB->SetFillStyle(3001);
    fHB->SetXTitle("#it{b} [fm]");
    fHB->SetDirectory(0);

    fHPhiR = new TH1D("phiR", "Event plane angle", 360, 0, 360);
    fHPhiR->SetMarkerColor(kMagenta+2);
    fHPhiR->SetFillColor(kMagenta+2);
    fHPhiR->SetFillStyle(3001);
    fHPhiR->SetXTitle("#it{#Phi}_{R} [degrees]");
    fHPhiR->SetDirectory(0);

    fHTime = new TH1D("timing", "Timing of processing", 5,0.5,5.5);
    fHTime->SetMarkerColor(kBlue+2);
    fHTime->SetFillColor(kBlue+2);
    fHTime->SetFillStyle(3001);
    fHTime->GetXaxis()->SetBinLabel(1,"Reset");
    fHTime->GetXaxis()->SetBinLabel(2,"Generate");
    fHTime->GetXaxis()->SetBinLabel(3,"Header");
    fHTime->GetXaxis()->SetBinLabel(4,"Particles");
    fHTime->GetXaxis()->SetBinLabel(5,"Filling");
    fHTime->SetDirectory(0);
				    
    fList = new TList;
    fList->SetName("list");

    TList* histos = new TList;
    histos->SetName("histograms");
    histos->SetOwner(true);
    histos->Add(fHEta);
    histos->Add(fHIpz);
    histos->Add(fHType);
    histos->Add(fHCent);
    histos->Add(fHB);
    histos->Add(fHPhiR);
    histos->Add(fHTime);
    fList->Add(histos);
    
    TList* estimators = new TList;
    estimators->SetName("estimators");
    estimators->SetOwner(true);
    fList->Add(estimators);
    
    TIter next(fCentEstimators);
    FastCentEstimator* estimator = 0;
    while ((estimator = static_cast<FastCentEstimator*>(next()))) {
      estimator->Setup(estimators, fTree,fIsTgtA,fIsProjA);
      estimator->SetVerbose(fVerbose);
      estimator->Print("nah");
    }
    
    if (fVerbose) Info("SetupOutput", "Adding list ot outputs");
    fOutput->Add(fList);
    // fOutput->ls();
    
    return true;
  }
  /** 
   * Set our seed.  This will read a UInt_t worth of data from
   * /dev/urandom
   * 
   */
  virtual void SetupSeed()
  {
    std::ifstream in("/dev/urandom");
    UInt_t seed = 0;
    in.read(reinterpret_cast<char*>(&seed), sizeof(seed));
    in.close();
    Printf("=== Random =============================\n"
	   " Random number seed: %d\n"
	   "========================================", seed);
    gRandom->SetSeed(seed);
  }
  virtual Bool_t SetupGRP()
  {
    Printf(" === Setup ==============================");
    Printf("  Run #:                          %6d", fRunNo);
    Printf("  EG:     %30s", fEGName.Data());
    Printf("  B range:             %5.1ffm - %5.1ffm", fBMin, fBMax);
    Printf(" ========================================");
    Printf("Macro path: %s", gROOT->GetMacroPath());

    // --- Check if we shoud get the GRP line ------------------------
    if (!fGRP && fInput) {
      fGRP = fInput->FindObject("GRP");
      std::ofstream* pout = new std::ofstream("grp.dat");
      if (pout) {
	if (fVerbose) 
	  Info("SetupGRP", "Writing GRP line '%s' to \"grp.dat\"",
	       fGRP->GetTitle());
	std::ostream& out = *pout;
	out << fGRP->GetTitle() << std::endl;
	pout->close();
      }
    }
    if (fVerbose) 
      Info("SetupGRP", "Overrides: %p Input: %p", fOverrides, fInput);
    if (!fOverrides && fInput) {
      fOverrides = static_cast<TList*>(fInput->FindObject("overrides"));
      if (!fOverrides && fVerbose)
	Info("SetupGRP", "No GRP overrides found in input:");
    }
    
    // --- Load our settings -----------------------------------------
    if (fVerbose) Info("SetupGRP", "Loading scripts");
    // Check if we have the global "grp" already 
    if (gROOT->ProcessLine("grp") == 0) 
      gROOT->Macro(Form("GRP.C(%d)", fRunNo));
    if (fVerbose) Info("SetupGRP", "Perhaps override");
    OverrideGRP();
    gROOT->ProcessLine("grp->Print()");
    return true;
  }    
  /** 
   * Set-up the generator. 
   * 
   * @return true on success 
   */
  virtual Bool_t SetupGen()
  {
    if (fVerbose) Info("SetupGen", "Load base config");
    gROOT->Macro("BaseConfig.C");
    if (fVerbose) Info("SetupGen", "Load EG config");
    gROOT->Macro("EGConfig.C");

    gROOT->ProcessLine(Form("VirtualEGCfg::LoadGen(\"%s\")",fEGName.Data()));

    // --- Make our generator ----------------------------------------
    // Info("SetupGen", "Creating generator");
    TString egMk = Form("egCfg->MakeGenerator(\"%s\",%f,%f)",
			fEGName.Data(), fBMin, fBMax);
    Long64_t egPtr = gROOT->ProcessLine(egMk);
    if (egPtr == 0) {
      Error("SetupGen", "Failed to make generator");
      return false;
    }
    fGenerator = reinterpret_cast<AliGenerator*>(egPtr);
    TString tgt, proj;
    Int_t   tgtA=0, tgtZ=0, projA=0, projZ=0;
    fGenerator->GetTarget(tgt, tgtA, tgtZ);
    fGenerator->GetProjectile(proj, projA, projZ);
    fIsTgtA  = !(tgtA  == tgtZ  && tgtA == 1);
    fIsProjA = !(projA == projZ && projZ == 1);
    if (fVerbose) 
      Info("SetupGen", "tgt=%s (%3d,%2d) proj=%s (%3d,%2d) CMS=%fGeV",
	   tgt.Data(), tgtA, tgtZ, proj.Data(), projA, projZ,
	   fGenerator->GetEnergyCMS());
    Info("SetupGen", "Generator: %s", fGenerator->GetTitle());
    if (fFileName.IsNull()) FileName();

    return true;
  }
  /** 
   * Setup the generator etc. of the job 
   * 
   * @param nev Maximum number of events per file 
   * 
   * @return true on success 
   */
  virtual Bool_t SetupRun()
  {
    // --- gAlice (bare ROOT) ----------------------------------------
    if (!gAlice)
      new AliRun("gAlice", "The ALICE Off-line framework");

    Long64_t nev = (fNEvents <= 0 ? 0xFFFFFFFF : fNEvents);
    Printf("=== Run ================================\n"
	   " Number of events: %lld\n"
	   "========================================", nev);
    TObject* ord      = (fInput ? fInput->FindObject("PROOF_Ordinal") : 0);
    UShort_t saveMode = 0;
    TString  post     = "";
    TString  dir      = "";
    if (ord) {
      TObject* save = fInput->FindObject("PROOF_SaveGALICE");
      if (save && fVerbose) {
	Info("SetupRun", "Got save option:");
	save->Print();
      }
      TString optSave(save ? save->GetTitle() : "split");
      optSave.ToLower();
      if       (optSave.EqualTo("none"))   saveMode = 0;
      else if  (optSave.EqualTo("merge"))  saveMode = 1;
      else if  (optSave.EqualTo("split"))  saveMode = 2;
      if (fProofFile && saveMode > 0) 
	dir  = fProofFile->GetDir(true);
      if (saveMode > 1)
	post = Form("_%s", ord->GetTitle());	
    }
    TString  galiceName(Form("%sgalice.root",dir.Data()));
    TString  kineName(Form("%sKinematics.root",dir.Data()));
    
    // --- Run-loader, stack, etc  -----------------------------------
    // Info("SetupRun", "Set-up run Loader");    
    fRunLoader = AliRunLoader::Open(galiceName, "FASTRUN", "RECREATE");
    fRunLoader->SetKineFileName(kineName);
    fRunLoader->SetCompressionLevel(2);
    fRunLoader->SetNumberOfEventsPerFile(nev);
    fRunLoader->LoadKinematics("RECREATE");
    fRunLoader->MakeTree("E");
    gAlice->SetRunLoader(fRunLoader);
    fRunLoader->MakeStack();
    fStack  = fRunLoader->Stack();
    fHeader = fRunLoader->GetHeader();

    // --- Initialize generator --------------------------------------
    // Info("SetupRun", "Initializing generator");
    fGenerator->Init();
    fGenerator->SetStack(fStack);

    if (saveMode < 1) {
      if (ord) 
	Info("SetupRun", "Not saving galice.root and Kinematics.root");
      return true;
    }
    
    TString aliceOut = Form("galice%s.root", post.Data());
    fAliceFile = new TProofOutputFile(aliceOut, "M");

    TString kineOut = Form("Kinematics%s.root", post.Data());
    fKineFile = new TProofOutputFile(kineOut, "M");
    
    return true;
  }
  /** 
   * Read the previously created grp.dat file 
   * 
   */
  Bool_t ReadGRPLine()
  {
    std::ifstream* pin = new std::ifstream("grp.dat");
    if (!pin) {
      Warning("ReadGRPLine", "Failed to open \"grp.dat\"");
      return false;
    }
    std::istream&  in  = *pin;
    TString line;
    TString env;
    do {
      line.ReadLine(in);
      if (line.IsNull()) continue;
      if (line.BeginsWith("#")) continue;
      env = line;
      break;
    } while (!in.eof());
    pin->close();

    if (env.IsNull()) {
      Warning("ReadGRPLine", "Got no line from \"grp.dat\"");
      return false;
    }
    
    fGRP = new TNamed("GRP",env.Data());
    if (fVerbose) Info("ReadGRPLine", "Read \"%s\"", env.Data());
    return true;
  }
  /** 
   * Possibly override settings from GRP. 
   * 
   */
  void OverrideGRP()
  {
    Long_t ret = gROOT->ProcessLine("grp");
    if (ret == 0) {
      Warning("OverrideGRP", "GRP not set yet, cannot override");
      return;
    }
    if (!fOverrides) {
      if (fVerbose) Info("OverrideGRP", "No overrides defined");
      return;
    }
    TIter next(fOverrides);
    TObject* o = 0;
    while ((o = next())) {
      if (fVerbose)
	Info("OverrideGRP", "Overriding GRP setting %s with %s",
	     o->GetName(), o->GetTitle());
      gROOT->ProcessLine(Form("grp->%s = %s;",
			      o->GetName(), o->GetTitle()));
    }
    Info("OverrideGRP", "After overriding:");
    gROOT->ProcessLine("grp->Print()");
  }
  /** 
   * Add an item to the list of things from GRP to override
   * 
   * @param field Field name of GRPData
   * @param value Field value of GRPData
   */
  void AddOverride(const TString& field, const TString& value)
  {
    if (!fOverrides) {
      fOverrides = new TList;
      fOverrides->SetName("overrides");
    }
    fOverrides->Add(new TNamed(field, value));
  }
  /** 
   * Set up job 
   * 
   */
  virtual void Init(TTree*)
  {
  }
  /** 
   * Set up job 
   * 
   */
  virtual void Begin(TTree*)
  {
    // Make a monitor
    // Info("Begin", "gProof=%p Nomonitor=%p",
    //      gProof, (gProof ? gProof->GetParameter("NOMONITOR") : 0));
    if (fVerbose) Info("Begin", "Called for FastSim");
       
    if (gProof && !gProof->GetParameter("NOMONITOR")) { 
      new FastMonitor;
      gProof->AddFeedback("list");
      // Info("Begin", "Adding monitoring");
    }
    gROOT->Macro(Form("GRP.C(%d)", fRunNo));
    if (ReadGRPLine()) {
      if(gProof) {
	gProof->AddInput(fGRP);
	if (fOverrides) gProof->AddInput(fOverrides);
      }
    }
    if (fVerbose) Info("Begin", "Perhaps override");
    OverrideGRP();
    if (fVerbose) Info("Begin", "Defining centrality estimators");
    fCentEstimators = new TList;
    // fCentEstimators->Add(new V0CentEstimator(-1));           //V0C
    // fCentEstimators->Add(new V0CentEstimator( 0));           //V0M
    // fCentEstimators->Add(new V0CentEstimator(+1));           //V0A
    fCentEstimators->Add(new V0CentEstimator(-1,true));      //V0CP
    fCentEstimators->Add(new V0CentEstimator( 0,true));      //V0MP
    fCentEstimators->Add(new V0CentEstimator(+1,true));      //V0AP
    // fCentEstimators->Add(new ZNCentEstimator(-1,zN|zE));     //ZNCE
    // fCentEstimators->Add(new ZNCentEstimator(+1,zN|zE));     //ZNAE 
    fCentEstimators->Add(new ZNCentEstimator(-1,true,false)); //ZNCE
    fCentEstimators->Add(new ZNCentEstimator(+1,true,false)); //ZNAE 
    fCentEstimators->Add(new ZNCentEstimator(-1,true,false,true)); //ZNCEP
    fCentEstimators->Add(new ZNCentEstimator(+1,true,false,true)); //ZNAE 
    fCentEstimators->Add(new ZNCentEstimator(-1,true,true));  //ZNCS
    fCentEstimators->Add(new ZNCentEstimator(+1,true,true));  //ZNAS
    // fCentEstimators->Add(new ZNCentEstimator(-1,zE|zP));     //ZPCE
    // fCentEstimators->Add(new ZNCentEstimator(+1,zE|zP));     //ZPAE
    // fCentEstimators->Add(new ZNCentEstimator(-1,zS));        //ZPCS
    // fCentEstimators->Add(new ZNCentEstimator(+1,zS));        //ZPAS
    // fCentEstimators->Add(new RefMultEstimator(0.8));
    // fCentEstimators->Add(new RefMultEstimator(0.5));
    TIter next(fCentEstimators);
    FastCentEstimator* estimator = 0;
    while ((estimator = static_cast<FastCentEstimator*>(next()))) {
      estimator->Print("nah");
    }
    if (fVerbose) Info("Begin", "End of begin");
  }
  /** 
   * Set-up this sub-job 
   * 
   */
  void SlaveBegin(TTree*)
  {
     if (fVerbose) Info("SlavesBegin", "Called for FastSim");
    SetupSeed();
    SetupGRP();
    SetupGen();
    SetupOutput();
    SetupRun();
  }
  /* Reset internal caches etc. 
   * 
   * @param iEv Event number 
   * 
   * @return true on success
   */
  virtual Bool_t PreEvent(Long64_t iEv)
  {
    // --- Reset header ----------------------------------------------
    fShortHead.Reset(fRunNo, iEv);
    fParticles->Clear();
    // --- Reset header, etc.  ---------------------------------------
    fHeader->Reset(fRunNo, iEv);
    fRunLoader->SetEventNumber(iEv);
    fStack->Reset();
    fRunLoader->MakeTree("K");

    TIter next(fCentEstimators);
    FastCentEstimator* estimator = 0;
    while ((estimator = static_cast<FastCentEstimator*>(next())))
      estimator->PreEvent();
    
    return true;
  }
  virtual void PbPbCent(Double_t b, Double_t& c, Float_t& np, Float_t& nc) const
  {
    // PbPb @ 2.76TeV only 
    // Updated 4th of November 2014 from 
    // cern.ch/twiki/bin/view/ALICE/CentStudies
    //        #Tables_with_centrality_bins_AN1
    if      (0.00 >= b  && b < 1.57)  { c=0.5;  np=403.8; nc=1861; } 
    else if (1.57 >= b  && b < 2.22)  { c=1.5;  np=393.6; nc=1766; } 
    else if (2.22 >= b  && b < 2.71)  { c=2.5;  np=382.9; nc=1678; } 
    else if (2.71 >= b  && b < 3.13)  { c=3.5;  np=372;   nc=1597; }  
    else if (3.13 >= b  && b < 3.50)  { c=4.5;  np=361.1; nc=1520; } 
    else if (3.50 >= b  && b < 4.94)  { c=7.5;  np=329.4; nc=1316; } 
    else if (4.94 >= b  && b < 6.05)  { c=12.5; np=281.2; nc=1032; } 
    else if (6.05 >= b  && b < 6.98)  { c=17.5; np=239;   nc=809.8; }
    else if (6.98 >= b  && b < 7.81)  { c=22.5; np=202.1; nc=629.6; }
    else if (7.81 >= b  && b < 8.55)  { c=27.5; np=169.5; nc=483.7; }
    else if (8.55 >= b  && b < 9.23)  { c=32.5; np=141;   nc=366.7; }
    else if (9.23 >= b  && b < 9.88)  { c=37.5; np=116;   nc=273.4; }
    else if (9.88 >= b  && b < 10.47) { c=42.5; np=94.11; nc=199.4; } 
    else if (10.47 >= b && b < 11.04) { c=47.5; np=75.3;  nc=143.1; } 
    else if (11.04 >= b && b < 11.58) { c=52.5; np=59.24; nc=100.1; }
    else if (11.58 >= b && b < 12.09) { c=57.5; np=45.58; nc=68.46; }
    else if (12.09 >= b && b < 12.58) { c=62.5; np=34.33; nc=45.79; }
    else if (12.58 >= b && b < 13.05) { c=67.5; np=25.21; nc=29.92; }
    else if (13.05 >= b && b < 13.52) { c=72.5; np=17.96; nc=19.08; }
    else if (13.52 >= b && b < 13.97) { c=77.5; np=12.58; nc=12.07; }
    else if (13.97 >= b && b < 14.43) { c=82.5; np=8.812; nc=7.682; }
    else if (14.43 >= b && b < 14.96) { c=87.5; np=6.158; nc=4.904; }
    else if (14.96 >= b && b < 15.67) { c=92.5; np=4.376; nc=3.181; }
    else if (15.67 >= b && b < 20.00) { c=97.5; np=3.064; nc=1.994; }
  }
  void pPbNpartCent(Int_t nPart, Double_t& c) const
  {
    Double_t cm[33];
    cm[0]       = -1;
    cm[1]       = 1;
    cm[2]	= 0.823203;
    cm[3]	= 0.723342;
    cm[4]	= 0.647615;
    cm[5]	= 0.583512;
    cm[6]	= 0.524967;
    cm[7]	= 0.469422;
    cm[8]	= 0.416093;
    cm[9]	= 0.362399;
    cm[10]	= 0.30997;
    cm[11]	= 0.258351;
    cm[12]	= 0.209152;
    cm[13]	= 0.164198;
    cm[14]	= 0.124642;
    cm[15]	= 0.0912396;
    cm[16]	= 0.0639762;
    cm[17]	= 0.0432744;
    cm[18]	= 0.0279033;
    cm[19]	= 0.0174596;
    cm[20]	= 0.010496;
    cm[21]	= 0.00588846;
    cm[22]	= 0.00319697;
    cm[23]	= 0.00165676;
    cm[24]	= 0.000807567;
    cm[25]	= 0.000393674;
    cm[26]	= 0.000177213;
    cm[27]	= 7.96863e-05;
    cm[28]	= 3.44911e-05;
    cm[29]	= 1.42722e-05;
    cm[30]	= 8.32543e-06;
    cm[31]	= 3.56804e-06;
    cm[32]	= 1.18935e-06;
    Int_t n = TMath::Max(TMath::Min(32, nPart),0);
    c = cm[n] * 100;
  }
  void pPbCent(Double_t b, Double_t& c, Float_t&, Float_t&) const
  {
    // pPb/Pbp @ 5.02TeV 
#if 1
    // From Alberica
    const Double_t cl[] = {   0,    5,   10,   20,   40,   60,   80, -1 };
    const Double_t ch[] = {   5,   10,   20,   40,   60,   80,  100, -1 };
    const Double_t bm[] = { 1.83675, 2.59375,  3.66875, 5.18625, 6.35475,
			    7.40225, 13.8577, -1};
    for (Int_t i = 0; i < 7; i++) {
      if (b <= bm[i]) {
	c = (ch[i]+cl[i]) / 2;
	return;
      }
    }
    c = 100;
#elif 0 
    // From Glauber
    const Double_t cl[] = {   0,    5,   10,   20,   40,   60,   80, -1 };
    const Double_t ch[] = {   5,   10,   20,   40,   60,   80,  100, -1 };
    const Double_t bm[] = {3.12, 3.50, 3.85, 4.54, 5.57, 6.63, 7.51, 15 };
    const Double_t bs[] = {1.39, 1.48, 1.57, 1.69, 1.69, 1.45, 1.11, -1 };
    static Double_t bl[] = { -1, -1, -1, -1, -1, -1, -1, -1, -1 };
    if (bl[0] == -1) {
      for (Int_t i = 0; i < 7; i++)
	bl[i] = bm[i] + (bm[i+1]-bm[i])/2;
    }
    
    Double_t* t   = std::lower_bound(bl, bl+7, b);
    Int_t     idx = t-bl;
    if (idx < 0 || idx > 6) { c = ch[6]; return; }
    c = (ch[idx]+cl[idx])/2;
#else 
    // cern.ch/twiki/bin/viewauth/ALICE/PACentStudies
    //          #Ncoll_in_centrality_bins_from_CL
    Double_t ac[] = { 2.5, 7.5, 15., 30., 50., 70., 90. };
    Double_t ab[] = { 100., 14.4, 13.8, 12.7, 10.2, 6.30, 3.10, 1.44, 0 };
    for (Int_t i = 0; i < 7; i++) {
      Double_t bC = ab[i+1];
      Double_t bH = bC + (ab[i]   - bC) / 2;
      Double_t bL = bC - (bC - ab[i+2]) / 2;
      if (b >= bL && b < bH)  {
	c = ac[i];
	break;
      }
    }
#endif 
  }    
    
  /** 
   * Process the event header 
   * 
   * @return true if the event should be diagnosed
   */
  virtual Bool_t ProcessHeader()
  {
    // --- Copy to short header --------------------------------------
    fShortHead.fRunNo   = fHeader->GetRun();
    fShortHead.fEventId = fHeader->GetEvent();
    TArrayF ip;
    fHeader->GenEventHeader()->PrimaryVertex(ip);
    fShortHead.fIpX     = ip[0];
    fShortHead.fIpY     = ip[1];
    fShortHead.fIpZ     = ip[2];

    // --- Check header type -----------------------------------------
    AliGenEventHeader* genHeader = fHeader->GenEventHeader();
    AliCollisionGeometry* geometry = 
      dynamic_cast<AliCollisionGeometry*>(genHeader);
    AliGenPythiaEventHeader* pythia    = 
      dynamic_cast<AliGenPythiaEventHeader*>(genHeader);
    AliGenDPMjetEventHeader* dpm       = 
      dynamic_cast<AliGenDPMjetEventHeader*>(genHeader);
    AliGenGeVSimEventHeader* gev       = 
      dynamic_cast<AliGenGeVSimEventHeader*>(genHeader);
    AliGenHerwigEventHeader* herwig    = 
      dynamic_cast<AliGenHerwigEventHeader*>(genHeader);
    AliGenCocktailEventHeader* cocktail = 
      dynamic_cast<AliGenCocktailEventHeader*>(genHeader);
    AliGenHijingEventHeader* hijing = 
      dynamic_cast<AliGenHijingEventHeader*>(genHeader);
    if (cocktail) {
      TList* headers = cocktail->GetHeaders();
      if (!headers) Warning("", "No headers in cocktail!");
      TIter next(headers);
      AliGenEventHeader* header = 0;
      AliCollisionGeometry* geom = 0;
      while ((header = static_cast<AliGenEventHeader*>(next()))) {
	AliCollisionGeometry* g = dynamic_cast<AliCollisionGeometry*>(header);
	if (g) geom = g;
	hijing = dynamic_cast<AliGenHijingEventHeader*>(header);
      }
      if (geom) geometry = geom;
    }
    if (hijing)   fShortHead.fEG = ShortHeader::kHijing;
    if (dpm)      fShortHead.fEG = ShortHeader::kDPMJet;
    if (pythia)   fShortHead.fEG = ShortHeader::kPythia;
    
    if (geometry) {
      fShortHead.fB          = geometry->ImpactParameter();
      fShortHead.fNtgt       = geometry->TargetParticipants();
      fShortHead.fNproj      = geometry->ProjectileParticipants();
      fShortHead.fNbin       = geometry->NN();
      fShortHead.fPhiR       = geometry->ReactionPlaneAngle();
      fShortHead.fNSpecNproj = geometry->ProjSpectatorsn(); 
      fShortHead.fNSpecPproj = geometry->ProjSpectatorsp();
      fShortHead.fNSpecNtgt  = geometry->TargSpectatorsn(); 
      fShortHead.fNSpecPtgt  = geometry->TargSpectatorsp();
   }
    // --- Determine diffraction flags -------------------------------
    Bool_t sd = false;
    Bool_t dd = false;
    if (pythia) {
      Int_t type = pythia->ProcessType();
      if (type < 100) { // pythia6
	switch (type) {
	case 92: case 93: sd = true; break;
	case 94:          dd = true; break;
	}
      }
      else {
	switch (type) { // Pythia8
	case 103: case 104: sd = true; break;
	case 105:           dd = true; break;
	}
      }
      fShortHead.fB     = pythia->GetImpactParameter();
      fShortHead.fNtgt  = 1;
      fShortHead.fNproj = 1;
      fShortHead.fNbin  = 1;
    }
    if (dpm) {
      // Try to figure out number of spectators in either side 
      UInt_t nSpecNproj = 0;
      UInt_t nSpecNtgt  = 0;
      UInt_t nSpecPproj = 0;
      UInt_t nSpecPtgt  = 0;
      Int_t  nPart      = fStack->GetNprimary();
      for (Int_t iPart = 0; iPart < nPart; iPart++) {
	TParticle*    particle  = fStack->Particle(iPart);
	Int_t         ks        = particle->GetStatusCode();
	Int_t         side      = 0;
	if (ks == 13) side = -1;
	if (ks == 14) side = +1;
	if (side == 0) continue;
	Int_t         kf        = particle->GetPdgCode();
	if (kf == kNeutron) {
	  if (side < 0) nSpecNproj++;
	  else          nSpecNtgt++;
	}
	else if (kf == kProton) { 
	  if (side < 0) nSpecPproj++;
	  else          nSpecPtgt++;
	}
      }
      fShortHead.fNSpecNtgt   = nSpecNtgt; 
      fShortHead.fNSpecNproj  = nSpecNproj;
      fShortHead.fNSpecPtgt   = nSpecPtgt;
      fShortHead.fNSpecPproj  = nSpecPproj;
      // fShortHead.Print();

      Int_t type = dpm->ProcessType();
#ifndef NO_DPMJET_TYPE
      switch (type) {
      case 5: case 6: sd = true;
      case 7:         dd = true;
      }      
#else
      static bool   first  = true;
      
      if (first) {
	Func_t add = gSystem->DynFindSymbol("*", "dtglcp_");
	if (!add) 
	  Warning("", "Didn't find dtglcp_");
	else 
	  _dtglcp = (DtglcpCommon*)add;
      }
      // The below - or rather a different implementation with some
      // errors - was proposed by Cvetan - I don't think it's right
      // though.  See also
      //
      //   https://cern.ch/twiki/pub/ALICE/PAPaperCentrality/normalization.pdf
      //   https://cern.ch/twiki/bin/view/ALICE/PAMCProductionStudies
      //
      Int_t nsd1=0, nsd2=0, ndd=0;
      Int_t npP = dpm->ProjectileParticipants();
      Int_t npT = dpm->TargetParticipants();
      // Get the numbeer of single and double diffractive participants
      dpm->GetNDiffractive(nsd1,nsd2,ndd);
      // Check if all partipants are single/double diffractive 
      if      ((ndd == 0) && ((npP == nsd1) || (npT == nsd2)))	sd = true;
      else if (ndd == (npP + npT))                       	dd = true;
      Int_t ncp = dpm->NN();
      Int_t nct = dpm->NNw();
      Int_t nwp = dpm->NwN();
      Int_t nwt = dpm->NwNw();
      Int_t nwtacc = _dtglcp->nwtacc;
      Int_t nwtsam = _dtglcp->nwtsam;
      if (first) {
	Printf("@ Npp sd1 Npt sd2  dd tpe Ncp Nct Nwp Nwt acc sam");
	first = false;
      }
      Printf("@ %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d",
	     npP, nsd1, npT, nsd2, ndd, type, ncp, nct, nwp, nwt,
	     nwtacc, nwtsam);
#endif 
    }
    if (gev) fShortHead.fPhiR = gev->GetEventPlane();
    if (herwig) {
      Int_t type = herwig->ProcessType();
      switch (type) {
      case 5: case 6: sd = true; break;
      }
      fShortHead.fNtgt  = 1;
      fShortHead.fNproj = 1;
      fShortHead.fNbin  = 1;
    }
    fShortHead.fType = (sd ? 0x1 : 0) | (dd ? 0x2 : 0);

    // --- Check centrality -----------------------------------------
    Float_t  np = -1;
    Float_t  nc = -1;
    Double_t c  = -1;
    Double_t b  = fShortHead.fB;
    Info("ProcessHeader", "b=%f isProjA=%d isTgtA=%d cms=%f",
	 b, fIsProjA, fIsTgtA, fGenerator ? fGenerator->GetEnergyCMS() : -1);
    if (b >= 0 && fIsProjA && fIsTgtA &&
	TMath::Abs(fGenerator->GetEnergyCMS()-2760) < 10) {
      PbPbCent(b, c, np, nc);
    }
    else if (b >= 0 && (fIsTgtA || fIsProjA) &&
	     (fGenerator
	      ? TMath::Abs(fGenerator->GetEnergyCMS()-5023) < 10
	      : true)) {
      pPbCent(b, c, np, nc);
      Info("", "Getting pPb centrality b=%f -> %f%%", b, c);
    }
    if (c >= 0) fShortHead.fC = c;
    Double_t nb   = nc/2;
    // Be careful to round off
    if ((fShortHead.fNtgt+fShortHead.fNproj) <= 0 && np >= 0) {
      fShortHead.fNtgt  = Int_t(np-nb+.5);
      fShortHead.fNproj = Int_t(nb+.5);
    }
    if (fShortHead.fNbin  <= 0 && nb >= 0)
      fShortHead.fNbin  = Int_t(nb+.5);
    
    // --- Check if within vertex cut -------------------------------
    Bool_t selected = (fShortHead.fIpZ <= fHIpz->GetXaxis()->GetXmax() &&
		       fShortHead.fIpZ >= fHIpz->GetXaxis()->GetXmin());

    // --- Only update histograms if within IPz cut ------------------
    if (selected) {
      fHPhiR->Fill(fShortHead.fPhiR*TMath::RadToDeg());
      fHB->Fill(fShortHead.fB);
      fHIpz->Fill(fShortHead.fIpZ);
      if (dd) fHType->Fill(3);
      if (sd) fHType->Fill(2);
      if (!dd && !sd) fHType->Fill(1);
      fHCent->Fill(c);
      // fShortHead.Print();
    }
    TIter next(fCentEstimators);
    FastCentEstimator* estimator = 0;
    while ((estimator = static_cast<FastCentEstimator*>(next())))
      estimator->ProcessHeader(fShortHead);
    
    return selected;
  }
  /** 
   * Process all particles 
   * 
   * @param selected True if particle information should be diagnosed
   * 
   * @return true on success
   */
  virtual Bool_t ProcessParticles(Bool_t selected)
  {
    Int_t nPart = fStack->GetNprimary();
    for (Int_t iPart = 0; iPart < nPart; iPart++) {
      TParticle*    particle  = fStack->Particle(iPart);
      TParticlePDG* pdg       = particle->GetPDG();
      Bool_t        primary   = fStack->IsPhysicalPrimary(iPart);
      Bool_t        weakDecay = fStack->IsSecondaryFromWeakDecay(iPart);
      Bool_t        charged   = (pdg &&  TMath::Abs(pdg->Charge()) > 0);
      if (primary)   particle->SetBit(BIT(14));
      if (weakDecay) particle->SetBit(BIT(15));
      if (charged)   particle->SetBit(BIT(16));

      new ((*fParticles)[iPart]) TParticle(*particle);

      TIter next(fCentEstimators);
      FastCentEstimator* estimator = 0;
      while ((estimator = static_cast<FastCentEstimator*>(next())))
	estimator->Process(particle);
      
      if (!selected || !charged || !primary) continue;
      Double_t pT    = particle->Pt();
      if (pT < 1e-10) continue; /// Along beam axis 
      Double_t pZ    = particle->Pz();
      Double_t theta = TMath::ATan2(pT, pZ);
      Double_t eta   = -TMath::Log(TMath::Tan(theta/2));
      fHEta->Fill(eta);
    }
    return true;
  }
  /** 
   * Do final event processing (fill output)
   * 
   */
  virtual void PostEvent()
  {
    fHeader->SetNprimary(fStack->GetNprimary());
    fHeader->SetNtrack(fStack->GetNtrack());

    TIter next(fCentEstimators);
    FastCentEstimator* estimator = 0;
    while ((estimator = static_cast<FastCentEstimator*>(next())))
      estimator->PostEvent();

    fTree->Fill();
    
    fStack->FinishEvent();
    fHeader->SetStack(fStack);
    
    fRunLoader->TreeE()->Fill();
    fRunLoader->WriteKinematics("OVERWRITE");
  }
  virtual void Generate() { fGenerator->Generate(); }
 
  /** 
   * Process one event 
   * 
   * @param iEv Event number 
   * 
   * @return true on success, false otherwize 
   */
  virtual Bool_t Process(Long64_t iEv)
  {
    // --- The stopwatch ---------------------------------------------
    TStopwatch timer;
    timer.Start();
    PreEvent(iEv);
    fHTime->Fill(1, timer.RealTime());
    
    // --- Generate event --------------------------------------------
    timer.Start();
    Generate();
    fHTime->Fill(2, timer.RealTime());

    // --- Process the header ----------------------------------------
    timer.Start();
    Bool_t selected = ProcessHeader();
    fHTime->Fill(3, timer.RealTime());
    
    // --- Loop over particles ---------------------------------------
    timer.Start();
    ProcessParticles(selected);
    fHTime->Fill(4, timer.RealTime());

    // --- Do final stuff --------------------------------------------    
    timer.Start();
    PostEvent();
    fHTime->Fill(5, timer.RealTime());

    return true;
  }
  virtual void FinishRun()
  {
    fGenerator->FinishRun();
    fRunLoader->WriteHeader("OVERWRITE");
    fGenerator->Write();
    fRunLoader->Write();

  }
  /** 
   * Finalize this sub-job  
   * 
   */
  void SlaveTerminate()
  {
    FinishRun();
    if (fFile) {
      if (fProofFile) {
	if (fVerbose) fProofFile->Print();
	fOutput->Add(fProofFile);
	fOutput->Add(new TH1F("filename", fFileName.Data(),1,0,1));
      }
      // Flush out tree 
      fFile->cd();
      fTree->Write(0, TObject::kOverwrite);
      fFile->Close();
      fFile->Delete();
      fFile = 0;
    }
    if (fAliceFile) {
      TFile* galice = GetGAlice();
      if (galice) {
	if (fVerbose) fAliceFile->Print();
	fAliceFile->AdoptFile(galice);
	fAliceFile->SetOutputFileName(fAliceFile->GetName());
	fOutput->Add(fAliceFile);
	galice->Write();
      }
    }
    if (fKineFile) {
      TFile* kine = GetKine();
      if (kine) {
	if (fVerbose) fKineFile->Print();
	fKineFile->AdoptFile(kine);
	fKineFile->SetOutputFileName(fKineFile->GetName());
	fOutput->Add(fKineFile);
	kine->Write();
      }
    }

    if (fVerbose) {
      Info("SlaveTerminate", "Content of output list");
      gROOT->IncreaseDirLevel();
      fOutput->ls();
      gROOT->DecreaseDirLevel();

      gSystem->Exec("echo \"Content of working directory\"");
      gSystem->Exec("ls -l1 | sed 's/^/  /'");
    }
  }
  /** 
   * Write a collection to disk, transforming sub-collections to
   * directories.
   * 
   * @param c 
   * @param dir 
   */
  void FlushList(TCollection* c, TDirectory* dir)
  {
    dir->cd();
    TIter next(c);
    TObject* o = 0;
    while ((o = next())) {
      if (o->IsA()->InheritsFrom(TCollection::Class())) {
	if (fVerbose) Info("FlushList", "Got collection: %s", c->GetName());
	TDirectory* cur = dir->mkdir(o->GetName());
	FlushList(static_cast<TCollection*>(o), cur);
	dir->cd();
	continue;
      }
      o->Write();
    }
    dir->cd();
  }
      
  /** 
   * Final processing of the data 
   * 
   */
  void Terminate()
  {
    if (gProof) gProof->ClearFeedback();

    if (!fList)
      fList = static_cast<TList*>(fOutput->FindObject("histograms"));
    if (!fList) {
      Error("Terminate", "No output list");
      return;
    }
    
    if (!fProofFile) {
      TObject* fn = fOutput->FindObject("filename");
      if (fn) fFileName  = fn->GetTitle();
      fProofFile =
	static_cast<TProofOutputFile*>(fOutput->FindObject(FileName()));
    }
    if (fProofFile) 
      fFile = fProofFile->OpenFile("UPDATE");
    if (!fFile)
      fFile = TFile::Open(FileName(),"UPDATE");

    TList* estimators = static_cast<TList*>(fList->FindObject("estimators"));
    TList* histos     = static_cast<TList*>(fList->FindObject("histograms"));
    if (!histos) {
      Warning("Terminate", "No histogram list found in output");
      fList->ls();
    }
    TIter next(fCentEstimators);
    FastCentEstimator* estimator = 0;
    while ((estimator = static_cast<FastCentEstimator*>(next())))
      estimator->Terminate(estimators);
	
    fHEta  = static_cast<TH1*>(histos->FindObject("dNdeta"));
    fHIpz  = static_cast<TH1*>(histos->FindObject("ipZ"));
    fHType = static_cast<TH1*>(histos->FindObject("type"));
    fHCent = static_cast<TH1*>(histos->FindObject("cent"));
    fHB    = static_cast<TH1*>(histos->FindObject("b"));
    fHPhiR = static_cast<TH1*>(histos->FindObject("phiR"));
    fHTime = static_cast<TH1*>(histos->FindObject("timing"));

    if (!(fHEta && fHIpz && fHType && fHB && fHPhiR && fHTime)) {
      Warning("Terminate", "Missing histograms (%p,%p,%p,%p,%p,%p)",
	      fHEta, fHIpz, fHType, fHB, fHPhiR, fHTime);
      return;
    }

    Int_t nTotal = fHIpz->GetEntries();
    fHEta ->Scale(1./nTotal, "width");
    fHB   ->Scale(1./nTotal, "width");
    fHPhiR->Scale(1./nTotal, "width");
    fHTime->Scale(1./nTotal, "width");

    if (!fFile){
      Warning("Terminate", "No file to write to");
      return;
    }

    FlushList(fList, fFile); // ->Write();
    
    fTree = static_cast<TTree*>(fFile->Get("T"));
    if (!fTree)  Warning("Terminate", "No tree");

    if (fVerbose) fFile->ls();
    fFile->Close();

    MoveAliceFiles();
  }
  /** 
   * Retrieve the galice.root file from ROOT 
   * 
   * @return Pointer to file or null
   */
  TFile* GetGAlice()
  {
    if (!fRunLoader) return 0;
    TString galiceName = fRunLoader->GetFileName();
    TFile*  file       = gROOT->GetFile(galiceName);
    if (!file) {
      Warning("GetGAlice", "Didn't find galice file \"%s\"", galiceName.Data());
      gROOT->GetListOfFiles()->ls();
      return 0;
    }
    return file;
  }
  /** 
   * Retrieve the Kinematics.root file from ROOT 
   * 
   * @return Pointer to file or null
   */
  TFile* GetKine()
  {
    if (!fRunLoader) return 0;
    TString kineName   = "Kinematics.root";
    TString galiceName = fRunLoader->GetFileName();
    TString dir        = gSystem->DirName(galiceName);
    if (dir.EqualTo(".")) dir = "";
    if (!dir.IsNull() && dir[dir.Length()-1] != '/') dir.Append("/");
    kineName.Prepend(dir);

    TFile*  file = gROOT->GetFile(kineName);
    if (!file) {
      Warning("GetKine", "Didn't find kinematics file \"%s\"", kineName.Data());
      gROOT->GetListOfFiles()->ls();
      return 0;
    }
    return file;
  }
  /** 
   * Move retrieved ALICE files (galice.root and Kinematics.root) to
   * separate su-directories, and create a collection of the TE tree
   * stored in the galice.root files.
   * 
   */
  void MoveAliceFiles()
  {
    if (!fInput) return;

    TObject* save  = fInput->FindObject("PROOF_SaveGALICE");
    if (!save) return;
    
    TString  sMode = save->GetTitle();
    if (!sMode.EqualTo("split", TString::kIgnoreCase)) return;

    TList*            lst   = new TList;
    TSystemDirectory* dir   = new TSystemDirectory(".",
						   gSystem->WorkingDirectory());
    TList*            files = dir->GetListOfFiles();
    TSystemFile*      file  = 0;
    TIter             next(files);
    while ((file = static_cast<TSystemFile*>(next()))) {
      if (file->IsDirectory()) continue;
      TString fn(file->GetName());
      if (!fn.BeginsWith("galice") && !fn.BeginsWith("Kinematics"))
	continue;

      TPRegexp regex("(.*)_([^_]+)\\.root");
      TObjArray* matches = regex.MatchS(fn);
      if (matches->GetEntriesFast() < 3) {
	delete matches;
	continue;
      }
      TString ord = matches->At(2)->GetName();
      TString bse = matches->At(1)->GetName();

      if (gSystem->AccessPathName(ord,kFileExists))
	gSystem->MakeDirectory(ord);

      if (fVerbose) 
	Info("MoveAliceFiles", "Moving %s to %s/%s.root",
	     fn.Data(), ord.Data(), bse.Data());
      file->Move(Form("%s/%s.root", ord.Data(), bse.Data()));

      if (!bse.EqualTo("galice")) continue;
      TObjString* url = new TObjString(Form("file://%s/%s/%s.root?#TE",
					    file->GetTitle(),
					    ord.Data(),
					    bse.Data()));
      if (fVerbose) 
	Info("MoveAliceFiles", "Adding \"%s\" to file list",
	     url->GetName());
      lst->Add(url);
    }
    if (lst->GetEntries() <= 0) return;
    if (fVerbose) lst->ls();
    
    TFile* out   = TFile::Open("index.root","RECREATE");
    lst->Write("TE",TObject::kSingleKey);
    out->Write();
    out->Close();
    
  }
  /** 
   * Interface version used 
   * 
   * @return 1y
   */
  Int_t Version() const { return 1; }

  /**
   * @{ 
   * @name Parameters 
   */
  TString  fEGName;               // Name of event generator
  Int_t    fRunNo;                // Run to simulate 
  Double_t fBMin;                 // Least impact parameter 
  Double_t fBMax;                 // Largest impact parameter
  TObject* fGRP;                  //! GRP in one line
  TList*   fOverrides;            //! GRP setting to override
  Long64_t fNEvents;              //  Number of requested events
  Bool_t   fIsTgtA;               //! True if target beam is nuclei
  Bool_t   fIsProjA;              //! True if projectile beam is nuclei
  /* @} */
  /** 
   * @{ 
   * @name ALICE EG interface 
   */
  AliGenerator* fGenerator;       //! Event generator
  AliRunLoader* fRunLoader;       //! Loader of trees
  AliStack*     fStack;           //! Stack of particles
  AliHeader*    fHeader;          //! Header handler
  /* @} */
  /** 
   * @{ 
   * @name Custom output 
   */
  TTree*        fTree;            //! Custom tree 
  TClonesArray* fParticles;       //! List of particles
  /**
   * @{ 
   * @name Diagnostics 
   */
  TList* fList;                   //! List of outputs
  TH1*   fHEta;                   //! dN/deta
  TH1*   fHIpz;                   //! IPz histogram
  TH1*   fHType;                  //! Event type histogram
  TH1*   fHCent;                  //! Event type histogram
  TH1*   fHB;                     //! B histogram
  TH1*   fHPhiR;                  //! Reaction plane
  TH1*   fHTime;                  //! Timing 
  /* @} */
  /** 
   * @{ 
   * @name Centrality 
   */
  TList* fCentEstimators;          // Centrality estimators
  /* @} */
  /**
   * @{ 
   * @name Output files 
   */
  TProofOutputFile* fProofFile;   //! Proof output file
  TProofOutputFile* fAliceFile;   //! 
  TProofOutputFile* fKineFile;    //! 
  TFile*            fFile;        //! Output file
  mutable TString   fFileName;    //! Output file name 
  /* @} */
  Bool_t fVerbose; // Verbosity 

#ifndef __CINT__
  ShortHeader fShortHead;
#endif 
  /** 
   * Run this selector as a normal process
   * 
   * @param nev        Number of events
   * @param run        Run number to anchor in
   * @param gen        Generator 
   * @param bMin       Least impact parameter [fm]
   * @param bMax       Largest impact parameter [fm]
   * @param monitor    Monitor frequency [s]
   * 
   * @return true on succes
   */
  static Bool_t  LocalRun(Long64_t       nev,
			  UInt_t         run,
			  const TString& gen,
			  Double_t       bMin,
			  Double_t       bMax,
			  Int_t          monitor,
			  Bool_t         verbose,
			  const TString& overrides="")
			  
  {
    FastSim* sim = new FastSim(gen,run,bMin,bMax,nev);
    SetOverrides(sim, overrides);
    sim->fVerbose = verbose;
    sim->Begin(0);
    sim->SlaveBegin(0);

    TTimer* timer = 0;
    if (monitor > 0) {
      // timer = new TTimer(new FastMonitor(sim), monitor*1000,true);
      timer = new TTimer(1000);
      timer->Connect("Timeout()","FastMonitor",
		     new FastMonitor(sim), "Handle()");
      // ::Info("Run", "Turning on monitoring");
      timer->Start(-1,false);
    }
      
    for (Long64_t i=0; i <nev; i++) {
      Printf("=== Event # %6lld/%6lld ==========================",
	     i+1, nev);
      sim->Process(i);
      if (timer && (i > 0) && (i % 500 == 0)) {
	if (timer->CheckTimer(gSystem->Now()))
	  Printf("Fired timer");
      }
    }
    if (timer) timer->TurnOff();
    sim->SlaveTerminate();
    sim->Terminate();

    return true;
  }
  /**
   * Load needed libraries in a proof serssion 
   */
  static void ProofLoadLibs()
  {
    if (!gProof) return;

    // Remember to copy changes to RunFast.C
    TList clsLib;
    clsLib.Add(new TNamed("TVirtualMC",              "libVMC"));
    clsLib.Add(new TNamed("TLorentzVector",          "libPhysics"));
    clsLib.Add(new TNamed("TLinearFitter",           "libMinuit"));
    clsLib.Add(new TNamed("TTree",                   "libTree"));
    clsLib.Add(new TNamed("TProof",                  "libProof"));
    clsLib.Add(new TNamed("TGFrame",                 "libGui"));
    clsLib.Add(new TNamed("TSAXParser",              "libXMLParser"));
    clsLib.Add(new TNamed("AliVEvent",               "libSTEERBase"));
    clsLib.Add(new TNamed("AliESDEvent",             "libESD"));
    clsLib.Add(new TNamed("AliAODEvent",             "libAOD"));
    clsLib.Add(new TNamed("AliAnalysisManager",      "libANALYSIS"));
    clsLib.Add(new TNamed("AliCDBManager",           "libCDB"));
    clsLib.Add(new TNamed("AliRawVEvent",            "libRAWDatabase"));
    clsLib.Add(new TNamed("AliHit",                  "libSTEER"));
    clsLib.Add(new TNamed("AliGenMC",                "libEVGEN"));
    clsLib.Add(new TNamed("AliFastEvent",            "libFASTSIM"));

    TIter next(&clsLib);
    TObject* obj = 0;
    while ((obj = next())) {
      gProof->Exec(Form("gROOT->LoadClass(\"%s\",\"%s\");",
			obj->GetName(), obj->GetTitle()));
    }
  }
  /** 
   * Run this selector in PROOF(Lite)
   * 
   * @param url        Proof URL
   * @param nev        Number of events
   * @param run        Run number to anchor in
   * @param gen        Generator 
   * @param bMin       Least impact parameter [fm]
   * @param bMax       Largest impact parameter [fm]
   * @param monitor    Monitor frequency [s]
   * @param opt        Compilation options
   * 
   * @return true on succes
   */
  static Bool_t ProofRun(const TUrl&    url,
			 Long64_t       nev,
			 UInt_t         run,
			 const TString& gen,
			 Double_t       bMin,
			 Double_t       bMax,
			 Int_t          monitor=-1,
			 Bool_t         verbose=false,
			 const TString& overrides="",
			 const TString& save="none",
			 const char*    opt="")
  {
    TProof::Reset(url.GetUrl());
    TProof::Open(url.GetUrl());
    gProof->ClearCache();

    TString phy = gSystem->ExpandPathName("$(ALICE_PHYSICS)");
    TString ali = gSystem->ExpandPathName("$(ALICE_ROOT)");
    // TString fwd = gSystem->ExpandPathName("$ANA_SRC");
    TString fwd = phy + "/PWGLF/FORWARD/analysis2";

    gProof->AddIncludePath(Form("%s/include", ali.Data()));
    gProof->AddIncludePath(Form("%s/include", phy.Data()));
    ProofLoadLibs();
    gProof->Load(Form("%s/sim/GRP.C",fwd.Data()), true);
    gProof->Load(Form("%s/sim/BaseConfig.C",fwd.Data()), true);
    gProof->Load(Form("%s/sim/EGConfig.C",fwd.Data()), true);

    // gROOT->ProcessLine("gProof->SetLogLevel(5);");
    gProof->Load(Form("%s/sim/FastSim.C+%s", fwd.Data(), opt),true);

    if (monitor <= 0) gProof->SetParameter("NOMONITOR", true/*ignored*/);
    else              gProof->SetParameter("PROOF_FeedbackPeriod",
					   monitor*1000/*ms*/);
    gProof->SetParameter("PROOF_SaveGALICE", save);

    FastSim* sim = new FastSim(gen,run,bMin,bMax,nev);
    SetOverrides(sim, overrides);
    sim->fVerbose = verbose;
    gProof->Process(sim, nev, "");

    return true; // status >= 0;
  }
  /** 
   * Extract key value pair from string 
   * 
   * @param in  Input string 
   * @param key On return, the key 
   * @param val On return, the value
   * @param sep Separator between key an value 
   * 
   * @return false of separator isn't found in input 
   */
  static Bool_t Str2KeyVal(const TString& in,
			   TString&       key,
			   TString&       val,
			   const char     sep='=')
  {
    Int_t idx = in.Index(sep);
    if (idx == kNPOS) return false;

    key = in(0,idx);
    val = in(idx+1, in.Length()-idx-1);
    return true;
  }
  static void SetOverrides(FastSim* sim, const TString& override)
  {
    if (override.IsNull()) return;

    const char* valid[] = { "beamEnergy", // UInt_t [GeV]
			    "energy",     // UInt_t [GeV]
			    "period",     // String			     
			    "run",        // UInt_t
			    "beam1.a",    // UInt_t
			    "beam1.z",    // UInt_t
			    "beam2.a",    // UInt_t
			    "beam2.z",    // UInt_t
			    0 };
    TObjArray*  tokens = override.Tokenize(",");
    TObjString* token  = 0;
    TIter       next(tokens);
    while ((token = static_cast<TObjString*>(next()))) {
      TString& str = token->String();
      if (str.IsNull()) continue;

      TString  key, val;
      if (!Str2KeyVal(str,key,val, ':')) {
	Printf("Warning: FastSim::Run: incomplete override '%s'",str.Data());
	continue;
      }
      const char** pvalid = valid;
      while (*pvalid) {
	if (key.EqualTo(*pvalid, TString::kIgnoreCase)) {
	  break;
	}
	pvalid++;
      }
      if (!*pvalid) {
	Printf("Warning: FastSim::Run: Invalid override '%s'", key.Data());
	continue;
      }
      // Special case for a string 
      if (key.EqualTo("period",TString::kIgnoreCase))
	val = Form("\"%s\"", val.Data());
      sim->AddOverride(*pvalid, val);
    }
    // delete tokens;
  }
  /** 
   * Run a simulation. 
   * 
   * @a url is the execution URL of the form 
   * 
   * @verbatim 
   PROTOCOL://[HOST[:PORT]]/[?OPTIONS]
   @endverbatim 
   *
   * Where PROTOCOL is one of 
   * 
   * - local for local (single thread) execution 
   * - lite for Proof-Lite execution 
   * - proof for Proof exection 
   * 
   * HOST and PORT is only relevant for Proof. 
   *
   * Options is a list of & separated options 
   * 
   * - events=NEV  Set the number of events to process 
   * - run=RUNNO   Set the run number to anchor in 
   * - eg=NAME     Set the event generator 
   * - b=RANGE     Set the impact parameter range in fermi 
   * - monitor=SEC Set the update rate in seconds of monitor histograms 
   *
   * 
   * @param url Exection URL 
   * @param opt Optimization used when compiling 
   * 
   * @return true on success 
   */
  static Bool_t Run(const char*  url, const char* opt="")
  {
    Printf("Will run fast simulation with:\n\n\t%s\n\n",url);
    Long64_t     nev     = 10000;
    UInt_t       run     = 0;
    TString      eg      = "default";
    TString      override= "";
    TString      save    = "none";
    Double_t     bMin    = 0;
    Double_t     bMax    = 20;
    Int_t        monitor = -1;
    Bool_t       verbose = false;
    TUrl         u(url);
    TString      out;
    TObjArray*   opts    = TString(u.GetOptions()).Tokenize("&");
    TObjString*  token   = 0;
    TIter        nextToken(opts);
    while ((token = static_cast<TObjString*>(nextToken()))) {
      TString& str = token->String();
      if (str.IsNull()) continue;

      if (str.EqualTo("verbose")) { verbose = true; str = ""; }

      if (str.IsNull()) continue;
      
      TString  key, val;
      if (!Str2KeyVal(str,key,val)) {
	if (!out.IsNull()) out.Append("&");
	out.Append(str);
	continue;
      }

      if      (key.EqualTo("events"))   nev      = val.Atoll();
      else if (key.EqualTo("run"))      run      = val.Atoi();
      else if (key.EqualTo("eg"))       eg       = val;
      else if (key.EqualTo("override")) override = val;
      else if (key.EqualTo("save"))     save     = val;
      else if (key.EqualTo("monitor"))  monitor  = val.Atoi();
      else if (key.EqualTo("b")) {
	TString min, max;
	if (Str2KeyVal(val, min, max, '-')) {
	  bMin = min.Atof();
	  bMax = max.Atof();
	}
      }
      else {
	if (!out.IsNull()) out.Append("&");
	out.Append(str);
      }
    }
    opts->Delete();
    u.SetOptions(out);
    if (!u.IsValid()) {
      Printf("Error: FastSim::Run: URL %s is invalid", u.GetUrl());
      return false;
    }
    
    Bool_t isLocal = TString(u.GetProtocol()).EqualTo("local");

    Printf("Run %s for %lld events anchored at %d\n"
	   "  Impact paramter range:  %5.1f-%5.1f fm\n"
	   "  Monitor frequency:      %d sec\n"
	   "  Execution url:          %s",
	   eg.Data(), nev, run, bMin, bMax, monitor, u.GetUrl());


    TStopwatch timer;
    timer.Start();

    Bool_t ret = false;
    if (isLocal)
      ret = LocalRun(nev, run, eg, bMin, bMax, monitor, verbose, override);
    else 
      ret = ProofRun(u, nev, run, eg, bMin, bMax,
		     monitor, verbose, override, save, opt);
    timer.Print();

    return ret;
  }
		    
  ClassDef(FastSim,3); 
};


struct EPosSim : public FastSim
{
  EPosSim(UInt_t run=0)
    : FastSim("epos", run, 0, 20, 100000),
      fInTree(0),
      fInNTot(0),
      fInB(0),
      fInPDG(0),
      fInStatus(0),
      fInPx(0),
      fInPy(0),
      fInPz(0),
      fInE(0),
      fInM(0),
      fInNcollH(0),
      fInNpartP(0),
      fInNpartT(0),
      fInNcoll(0),
      fInNSpcPN(0),
      fInNSpcTN(0),
      fInNSpcPP(0),
      fInNSpcTP(0),
      fInPhiR(0)
  {
  }
  Bool_t SetupBranches()
  {
    fInNTot   = fInTree->GetLeaf("nPart");
    fInB      = fInTree->GetLeaf("ImpactParameter");
    fInPDG    = fInTree->GetLeaf("pdgid");
    fInStatus = fInTree->GetLeaf("status");
    fInPx     = fInTree->GetLeaf("px");
    fInPy     = fInTree->GetLeaf("py");
    fInPz     = fInTree->GetLeaf("pz");
    fInE      = fInTree->GetLeaf("E");
    fInM      = fInTree->GetLeaf("m");
    // These are probably EPOS-LHC specific 
    fInNcollH = fInTree->GetLeaf("Ncoll_hard");
    fInNpartP = fInTree->GetLeaf("Npart_proj");
    fInNpartT = fInTree->GetLeaf("Npart_targ");
    fInNcoll  = fInTree->GetLeaf("Ncoll");
    fInNSpcPN = fInTree->GetLeaf("Nspec_proj_neut");
    fInNSpcTN = fInTree->GetLeaf("Nspec_targ_neut");
    fInNSpcPP = fInTree->GetLeaf("Nspec_proj_prot");
    fInNSpcTP = fInTree->GetLeaf("Nspec_targ_prot");
    fInPhiR   = fInTree->GetLeaf("phiR");
    
    

    return (fInNTot && fInB && fInPDG && fInPx && fInPy && fInPz && fInE);
  }
  void Begin(TTree* tree)
  {
    FastSim::Begin(tree);
    TIter next(fCentEstimators);
    FastCentEstimator* estimator = 0;
    while ((estimator = static_cast<FastCentEstimator*>(next()))) {
      if (!estimator->IsA()->InheritsFrom(V0CentEstimator::Class()))
	continue;
      V0CentEstimator* v = static_cast<V0CentEstimator*>(estimator);
      v->Flip(false); // Flip detector acceptance of V0A/C
      // v->Print("nah");
    }
  }
  void Init(TTree* tree)
  {
    Info("Init", "Initializing with tree %p (%s)",
	 tree, (tree ? tree->ClassName() : ""));
    if (!tree) return;

    TFile* file = tree->GetCurrentFile();
    Info("Init", "Current file: (%p) %s", file,
	 (file ? file->GetName() : ""));
    
    fInTree = tree;
    if (!SetupBranches())
      Fatal("Init", "Failed to set-up branches");
    // if (!SetupEstimator())
    //   Fatal("Init", "Failed to set-up estimator");
  }    
  /** 
   * Called when the file changes 
   * 
   * @return true on success, false otherwise 
   */
  Bool_t Notify()
  {
    if (!fInTree) {
      Warning("Notify", "No tree set yet!");
      return false;
    }
    TFile* file = fInTree->GetCurrentFile();
    Info("Notify", "processing file: (%p) %s", file,
	 (file ? file->GetName() : ""));
    if (!file) return true;
    if (!SetupBranches()) {
      Warning("Notify", "Failed to set-up branches");
      return false;
    }
    return true;
  }
  void SetupSeed() {}
  Bool_t SetupGen()
  {
    fIsTgtA  = gROOT->ProcessLine("grp->beam1.IsA()");
    fIsProjA = gROOT->ProcessLine("grp->beam2.IsA()");
    return true;
  }
  Bool_t SetupRun() { return true; }
  Bool_t PreEvent(Long64_t iEv)
  {
    Int_t read = fInTree->GetTree()->GetEntry(iEv);
    if (read <= 0) return false;
    
    // --- Reset header ----------------------------------------------
    fShortHead.Reset(fRunNo, iEv);
    fParticles->Clear();
    // Reset input

    TIter next(fCentEstimators);
    FastCentEstimator* estimator = 0;
    while ((estimator = static_cast<FastCentEstimator*>(next())))
      estimator->PreEvent();
    
    return true;    
  }
  Bool_t ProcessHeader()
  {
    fShortHead.fIpX        = 0;
    fShortHead.fIpY        = 0;
    fShortHead.fIpZ        = 0;
    fShortHead.fB          = fInB->GetValue();
    fShortHead.fNtgt       = (fInNpartT ? fInNpartT->GetValue() : 0);
    fShortHead.fNproj      = (fInNpartP ? fInNpartP->GetValue() : 0);
    fShortHead.fNbin       = (fInNcoll  ? fInNcoll ->GetValue() : 0);
    fShortHead.fPhiR       = (fInPhiR   ? fInPhiR  ->GetValue() : 0);
    fShortHead.fNSpecNproj = (fInNSpcPN ? fInNSpcPN->GetValue() : 0);
    fShortHead.fNSpecNtgt  = (fInNSpcTN ? fInNSpcTN->GetValue() : 0);
    fShortHead.fNSpecPproj = (fInNSpcPP ? fInNSpcPP->GetValue() : 0);
    fShortHead.fNSpecPtgt  = (fInNSpcTP ? fInNSpcTP->GetValue() : 0);
    fShortHead.fEG         = ShortHeader::kEPOS;
    
    Double_t c = -1;
    Float_t np = -1, nc = -1;
    if (fIsProjA && fIsTgtA) 
      PbPbCent(fShortHead.fB, c, np, nc);
    else if (fIsProjA && !fIsTgtA)
      pPbCent(fShortHead.fB, c, np, nc);
    
    if (c >= 0) fShortHead.fC = c;
    Double_t nb   = nc/2;
    // Be careful to round off
    if ((fShortHead.fNtgt+fShortHead.fNproj) <= 0 && np >= 0) {
      fShortHead.fNtgt  = Int_t(np-nb+.5);
      fShortHead.fNproj = Int_t(nb+.5);
    }
    if (fShortHead.fNbin  <= 0 && nb >= 0)
      fShortHead.fNbin  = Int_t(nb+.5);

    // --- Check if within vertex cut -------------------------------
    Bool_t selected = (fShortHead.fIpZ <= fHIpz->GetXaxis()->GetXmax() &&
		       fShortHead.fIpZ >= fHIpz->GetXaxis()->GetXmin());

    // --- Only update histograms if within IPz cut ------------------
    if (selected) {
      fHPhiR->Fill(fShortHead.fPhiR*TMath::RadToDeg());
      fHB->Fill(fShortHead.fB);
      fHIpz->Fill(fShortHead.fIpZ);
      fHCent->Fill(c);
      // fShortHead.Print();
    }
    TIter next(fCentEstimators);
    FastCentEstimator* estimator = 0;
    while ((estimator = static_cast<FastCentEstimator*>(next())))
      estimator->ProcessHeader(fShortHead);
    return selected;
  }
  virtual Bool_t ProcessParticles(Bool_t selected)
  {
    Int_t nTot = fInNTot->GetValue();
    for (Int_t iPart = 0; iPart < nTot; iPart++) {
      Int_t    status = fInStatus->GetValue(iPart);
      Int_t    pdg    = fInPDG->GetValue(iPart);
      Double_t px     = fInPx->GetValue(iPart);
      Double_t py     = fInPy->GetValue(iPart);
      Double_t pz     = fInPz->GetValue(iPart);
      // Double_t pz     = -fInPz->GetValue(iPart); // flip sign on pZ
      Double_t e      = fInE->GetValue(iPart);
      // Double_t m      = fInM->GetValue(iPart);

      TParticle* particle =
	new ((*fParticles)[iPart]) TParticle(pdg, status,-1,-1,-1,-1,
					     px, py, pz, e, 0, 0, 0, 0);
      TParticlePDG* pdgP      = particle->GetPDG();
      Bool_t        primary   = status == 1;
      Bool_t        weakDecay = false;
      Bool_t        charged   = (pdgP &&  TMath::Abs(pdgP->Charge()) > 0);
      if (primary)   particle->SetBit(BIT(14));
      if (weakDecay) particle->SetBit(BIT(15));
      if (charged)   particle->SetBit(BIT(16));

      TIter next(fCentEstimators);
      FastCentEstimator* estimator = 0;
      while ((estimator = static_cast<FastCentEstimator*>(next())))
	estimator->Process(particle);
      
      if (!selected || !charged || !primary) continue;
      Double_t pT    = particle->Pt();
      if (pT < 1e-10) continue; /// Along beam axis 
      Double_t pZ    = particle->Pz();
      Double_t theta = TMath::ATan2(pT, pZ);
      Double_t eta   = -TMath::Log(TMath::Tan(theta/2));
      fHEta->Fill(eta);
    }
    return true;
  }
  void PostEvent()
  {
    fTree->Fill();

    TIter next(fCentEstimators);
    FastCentEstimator* estimator = 0;
    while ((estimator = static_cast<FastCentEstimator*>(next())))
      estimator->PostEvent();
  }    
  void Generate() {}
  void FinishRun() {}
  TTree* fInTree;   //!
  TLeaf* fInNTot;   //!
  TLeaf* fInB;      //!
  TLeaf* fInPDG;    //!
  TLeaf* fInStatus; //! 
  TLeaf* fInPx;     //!
  TLeaf* fInPy;     //! 
  TLeaf* fInPz;     //! 
  TLeaf* fInE;      //! 
  TLeaf* fInM;      //!
  TLeaf* fInNcollH; //!
  TLeaf* fInNpartP; //!
  TLeaf* fInNpartT; //!
  TLeaf* fInNcoll;  //!
  TLeaf* fInNSpcPN; //!
  TLeaf* fInNSpcTN; //!
  TLeaf* fInNSpcPP; //!
  TLeaf* fInNSpcTP; //!
  TLeaf* fInPhiR;   //!
  
  /** 
   * Run this selector as a normal process
   * 
   * @param nev        Number of events
   * @param run        Run number to anchor in
   * @param gen        Generator 
   * @param bMin       Least impact parameter [fm]
   * @param bMax       Largest impact parameter [fm]
   * @param monitor    Monitor frequency [s]
   * 
   * @return true on succes
   */
  static Bool_t  LocalRun(Long64_t       nev,
			  UInt_t         run,
			  Int_t          monitor,
			  Bool_t         verbose)
			  
  {
    EPosSim* sim = new EPosSim(run);
    sim->fVerbose = verbose;
    sim->Begin(0);
    sim->SlaveBegin(0);

    TTimer* timer = 0;
    if (monitor > 0) {
      // timer = new TTimer(new FastMonitor(sim), monitor*1000,true);
      timer = new TTimer(1000);
      timer->Connect("Timeout()","FastMonitor",
		     new FastMonitor(sim), "Handle()");
      // ::Info("Run", "Turning on monitoring");
      timer->Start(-1,false);
    }
      
    for (Long64_t i=0; i <nev; i++) {
      Printf("=== Event # %6lld/%6lld ==========================",
	     i+1, nev);
      sim->Process(i);
      if (timer && (i > 0) && (i % 500 == 0)) {
	if (timer->CheckTimer(gSystem->Now()))
	  Printf("Fired timer");
      }
    }
    if (timer) timer->TurnOff();
    sim->SlaveTerminate();
    sim->Terminate();

    return true;
  }
  /** 
   * Run this selector in PROOF(Lite)
   * 
   * @param url        Proof URL
   * @param nev        Number of events
   * @param run        Run number to anchor in
   * @param gen        Generator 
   * @param bMin       Least impact parameter [fm]
   * @param bMax       Largest impact parameter [fm]
   * @param monitor    Monitor frequency [s]
   * @param opt        Compilation options
   * 
   * @return true on succes
   */
  static Bool_t SetupProof(const TUrl&    url,
			   Int_t          monitor=-1,
			   const char*    opt="")
  {
    TProof::Reset(url.GetUrl());
    TProof::Open(url.GetUrl());
    gProof->ClearCache();

    TString phy = gSystem->ExpandPathName("$(ALICE_PHYSICS)");
    TString ali = gSystem->ExpandPathName("$(ALICE_ROOT)");
    // TString fwd = gSystem->ExpandPathName("$ANA_SRC");
    TString fwd = phy + "/PWGLF/FORWARD/analysis2";

    gProof->AddIncludePath(Form("%s/include", ali.Data()));
    gProof->AddIncludePath(Form("%s/include", phy.Data()));
    ProofLoadLibs();
    gProof->Load(Form("%s/sim/GRP.C",fwd.Data()), true);
    gProof->Load(Form("%s/sim/BaseConfig.C",fwd.Data()), true);
    gProof->Load(Form("%s/sim/EGConfig.C",fwd.Data()), true);

    // gROOT->ProcessLine("gProof->SetLogLevel(5);");
    gProof->Load(Form("%s/sim/FastSim.C+%s", fwd.Data(), opt),true);

    if (monitor <= 0) gProof->SetParameter("NOMONITOR", true/*ignored*/);
    else              gProof->SetParameter("PROOF_FeedbackPeriod",
					   monitor*1000/*ms*/);

    return true; // status >= 0;
  }
  /** 
   * Run a simulation. 
   * 
   * @a url is the execution URL of the form 
   * 
   * @verbatim 
   PROTOCOL://[HOST[:PORT]]/[?OPTIONS]
   @endverbatim 
   *
   * Where PROTOCOL is one of 
   * 
   * - local for local (single thread) execution 
   * - lite for Proof-Lite execution 
   * - proof for Proof exection 
   * 
   * HOST and PORT is only relevant for Proof. 
   *
   * Options is a list of & separated options 
   * 
   * - events=NEV  Set the number of events to process 
   * - run=RUNNO   Set the run number to anchor in 
   * - eg=NAME     Set the event generator 
   * - b=RANGE     Set the impact parameter range in fermi 
   * - monitor=SEC Set the update rate in seconds of monitor histograms 
   *
   * 
   * @param url Exection URL 
   * @param opt Optimization used when compiling 
   * 
   * @return true on success 
   */
  static Bool_t Run(const char*  url, const char* opt="")
  {
    Printf("Will run fast simulation with:\n\n\t%s\n\n",url);
    UInt_t       run     = 0;
    Long64_t     nev     = -1;
    Int_t        monitor = -1;
    Bool_t       verbose = false;
    TUrl         u(url);
    TString      out;
    TObjArray*   opts    = TString(u.GetOptions()).Tokenize("&");
    TObjString*  token   = 0;
    TIter        nextToken(opts);
    while ((token = static_cast<TObjString*>(nextToken()))) {
      TString& str = token->String();
      if (str.IsNull()) continue;

      if (str.EqualTo("verbose")) { verbose = true; str = ""; }

      if (str.IsNull()) continue;
      
      TString  key, val;
      if (!Str2KeyVal(str,key,val)) {
	if (!out.IsNull()) out.Append("&");
	out.Append(str);
	continue;
      }

      if      (key.EqualTo("run"))      run      = val.Atoi();
      else if (key.EqualTo("monitor"))  monitor  = val.Atoi();
      else if (key.EqualTo("events"))   nev      = val.Atoi();
      else {
	if (!out.IsNull()) out.Append("&");
	out.Append(str);
      }
    }
    opts->Delete();
    u.SetOptions(out);
    if (!u.IsValid()) {
      Printf("Error: FastSim::Run: URL %s is invalid", u.GetUrl());
      return false;
    }
    
    Printf("Run EPos anchored at %d\n"
	   "  Monitor frequency:      %d sec\n"
	   "  Execution url:          %s",
	   run, monitor, u.GetUrl());

    TString treeName = u.GetAnchor();
    if (treeName.IsNull()) treeName = "Particle";
    TFile*  file = TFile::Open(u.GetFile(), "READ");
    if (!file) {
      Printf("Error: FastAnalysis::Run: Failed to open %s",
	     u.GetFile());
      return false;
    }

    TChain*   chain = new TChain(treeName, treeName);
    TObject*  o = file->Get(treeName);
    if (!o) {
      Printf("Error: FastAnalysis::Run: Couldn't get %s from %s",
	     treeName.Data(), u.GetFile());
      file->Close();
      return false;
    }
    Int_t cret = 0;
    if (o->IsA()->InheritsFrom(TChain::Class()))
      cret = chain->Add(static_cast<TChain*>(o));
    else if (o->IsA()->InheritsFrom(TTree::Class()))
      cret = chain->AddFile(u.GetFile());
    else if (o->IsA()->InheritsFrom(TCollection::Class()))
      cret = chain->AddFileInfoList(static_cast<TCollection*>(o));
    else if (o->IsA()->InheritsFrom(TFileCollection::Class()))
      cret = chain->AddFileInfoList(static_cast<TFileCollection*>(o)
				    ->GetList());
    file->Close();
    if (cret <= 0 || chain->GetListOfFiles()->GetEntries() <= 0) {
      Printf("Error: FastAnalysis::Run: Failed to create chain");
      return false;
    }

    TString       proto    = u.GetProtocol();
    Bool_t        isProof  = (proto.EqualTo("proof") || proto.EqualTo("lite"));
    if (isProof) {
      if (!SetupProof(u,monitor,opt)) return false;
      chain->SetProof();
    }

    EPosSim* sim = new EPosSim(run);
    sim->fVerbose = verbose;
    if (nev < 0) nev = TChain::kBigNumber;
    
    TStopwatch timer;
    timer.Start();

    Long64_t ret = chain->Process(sim, "", nev, 0);
    timer.Print();

    return ret > 0;
  }
  Int_t Version() const { return 2; }
  ClassDef(EPosSim,1);
};

#endif
//
// EOF
// 
