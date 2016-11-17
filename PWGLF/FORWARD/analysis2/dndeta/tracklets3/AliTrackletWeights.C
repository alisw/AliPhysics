/**
 * @file   AliTrackletWeights.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Wed Apr 27 16:50:58 2016
 * 
 * @brief  Encode simulation weights for 2nd pass 
 * 
 * 
 * @ingroup pwglf_forward_tracklets
 */
#ifndef ALITRACKLETWEIGHTS_C
#define ALITRACKLETWEIGHTS_C
#include <TNamed.h>
#include <map>
#ifndef __CINT__
# include <TH1.h>
# include <TH2.h>
# include <TH3.h>
# include <TAxis.h>
# include <TList.h>
# include <TBrowser.h>
# include <TCanvas.h>
# include <THStack.h>
# include <TROOT.h>
# include <TClass.h>
# include "AliAODTracklet.C"
# include <TParameter.h>
# include <TParticle.h>
#else
class TH1D;
class TH2D;
class TH1;
class TH2;
class TH3;
class TAxis;
class TList;
class TBrowser;
class THStack;
class AliAODTracklet;
class TCanvas; // Autoload
class TVirtualPad;
class TParticle;
#endif

//====================================================================
#ifndef __CINT__
namespace {
  template <typename T>
  T GetPar(const char* name, TCollection* l)
  {
    TObject* o = l->FindObject(name);
    if (!o) {
      Warning("GetPar", "Didn't find parameter %s in %s",
	      name, l->GetName());
      return T();
    }
    TClass* cls = TParameter<T>::Class();
    if (!o->IsA()->InheritsFrom(cls)) {
      Warning("GetPar", "Object %s is a %s, not a %s",
	      name, o->ClassName(), cls->GetName());
      return T();
    }
    TParameter<T>* p = static_cast<TParameter<T>*>(o);
    return p->GetVal();
  }
}
#endif

//====================================================================
/**
 * Base class for tracklets weights 
 */
class AliTrackletBaseWeights : public TNamed
{
public:
  /**						
   * Mode calculation 
   */
  enum ECalc {
    /** 
     * Tracklet weights calculated as 
     *
     * @f[
     w_{ij} = \left\{\begin{array}{cl}
       w_i & i=j \mbox{ good tracklets}\\
       w_i w_j & i\neq j \mbox{ fake tracklets}	
       \end{array}\right.
    @f] 
    */
    kProduct,
    /** 
     * Tracklet weights calculated as 
     *
     * @f[
     w_{ij} = \sqrt{w_i w_j} = \left\{\begin{array}{cl}
       w_i & i=j \mbox{ good tracklets}\\
       \sqrt{w_i w_j} & i\neq j \mbox{ fake tracklets}	
       \end{array}\right.
    @f] 
    */
    kSquare,
    /** 
     * Tracklet weights calculated as 
     *
     * @f[
     w_{ij} = 1+(w_i-1)+(w_j-1) 
    @f] 
    */
    kSum,
    /** 
     * Tracklet weights calculated as 
     *
     * @f[
     w_{ij} = 1+\frac{(w_i-1)+(w_j-1)}{2} = \frac{w_i + w_j}{2} 
     = \left\{\begin{array}{cl}
       w_i & i=j \mbox{ good tracklets}\\
       (w_i + w_j)/2 & i\neq j \mbox{ fake tracklets}	
       \end{array}\right.
    @f] 
    */
    kAverage
  };
  /**
   * Default constructor - ROOT I/O only
   * 
   */
  AliTrackletBaseWeights()
    : TNamed(),
      fCalc(kProduct),
      fMask(0xFF),
      fVeto(0x0),
      fInverse(false),
      fDebug(0)
  {}
  /**
   * Named constructor 
   * 
   * @param name  Name of object 
   * @param title (optional) free-form title 
   */
  AliTrackletBaseWeights(const char* name,
			 const char* title="Sim. tracklet weights")
    : TNamed(name,title),
      fCalc(kProduct),
      fMask(0xFF),
      fVeto(0x0),
      fInverse(false),
      fDebug(0)
  {}
  /**
   * Copy constructor 
   *
   * @param o Object to copy from  
   */
  AliTrackletBaseWeights(const AliTrackletBaseWeights& o)
    : TNamed(o),
      fCalc   (o.fCalc),
      fMask   (o.fMask),
      fVeto   (o.fVeto),
      fInverse(o.fInverse),
      fDebug  (o.fDebug)
  {}
  /**
   * Destructor 
   */
  virtual ~AliTrackletBaseWeights() {}

  /**
   * Assignment operator  
   *
   * @param o Object to assign from  
   *
   * @return Reference to this object 
   */
  AliTrackletBaseWeights& operator=(const AliTrackletBaseWeights& o)
  {
    if (&o == this) return *this;
    TNamed::operator=(o);
    fCalc    = o.fCalc;
    fMask    = o.fMask;
    fVeto    = o.fVeto;
    fInverse = o.fInverse;
    fDebug   = o.fDebug;
    return *this;
  }
  /** 
   * Set the mode for calculating the weight of a tracklet.  A
   * tracklet consist of two clusters with labels @f$ i@f$ and @f$
   * j@f$.  A weight (@f$ w_i@f$ and @f$ w_j@f$) is assigned to each
   * of the particles corresponding to the tracks with these labels.
   * The weight is calculated from the particle species and transverse
   * momentum of the primary mother particle of the tracks @f$ i@f$
   * and @f$ j@f$.  The weight of a tracklet can then be calculated in
   * four different ways
   *
   * @b Product 
   * @f[
   w = \left\{\begin{array}{cl}
   w_i & \mbox{for} i=j\\
   w_i w_j & \mbox{for} i\neq j\\
   \end{array}\right.\quad,
   @f] 
   * 
   * @b Square root of product of weights 
   *
   * @f[
   w = \left\{\begin{array}{cl}
   w_i & \mbox{for} i=j\\
   \sqrt(w_i w_j) & \mbox{for} i\neq j\\
   \end{array}\right.\quad,
   @f] 
   * 
   * @b Sum of weights 
   *
   * @f[ 
   w = 1 + \left\{\begin{array}{cl}
   2(w_i-1) & \mbox{for} i=j\\
   (w_i-1) + (w_j-1) & \mbox{for} i\neq j\\
   \end{array}\right.\quad,
   @f]    
   * 
   * @b Average of weights 
   *
   * @f[ 
   w = 1 + \left\{\begin{array}{cl}
   (w_i-1) & \mbox{for} i=j\\
   ((w_i-1) + (w_j-1))/2 & \mbox{for} i\neq j\\
   \end{array}\right.\quad.
   @f]    
   * 
   * @param mode Whether to take square root or not 
   */
  void SetCalc(UChar_t mode=kProduct) { fCalc = mode; }
  /** 
   * Set the tracklet mask 
   * 
   * @param mask Mask to use 
   */
  void SetMask(UChar_t mask) { fMask = mask; }
  /** 
   * Set the tracklet veto 
   * 
   * @param veto Veto to use 
   */
  void SetVeto(UChar_t veto) { fVeto = veto; }
  /** 
   * Inverse the weights calculated.  That is, if this option is
   * enabled, then the weight used is @f$1/w@f$ where @f$w@f$ is the
   * normal weight.
   * 
   * @param inv If true, inverse weights 
   */
  void SetInverse(Bool_t inv) { fInverse = inv; }
  /** 
   * Set the debug level 
   * 
   * @param lvl Debug level 
   */
  void SetDebug(Int_t lvl) { fDebug = lvl; }
  /** 
   * Check if tracklet is to be reweighed according to mask and veto 
   * 
   * @param tracklet Tracklet 
   * 
   * @return true if to be weighed, false otherwise 
   */
  Bool_t CheckTracklet(const AliAODTracklet* tracklet) const
  {
    UChar_t flags = tracklet->GetFlags();
    if (fMask != 0xFF && (fMask & flags) == 0) {
      if (fDebug > 3) 
	Info("LookupWeight", "Tracklet 0x%02x does not fullfill mask 0x%02x",
             flags, fMask);	  
      return false;
    }
    if ((fVeto & flags) != 0) {
      if (fDebug > 3) 
	Info("LookupWeight", "Tracklet 0x%02x vetoed by 0x%02x",flags, fVeto);
      return false;
    }
    return true;
  }    
  /** 
   * Find the tracklet weight 
   * 
   * @param tracklet Tracklet 
   * @param cent     Centrality 
   * @param ipz      Interaction point Z coordinate 
   * @param corr     Optional histogram to fill with correlation of weights
   * 
   * @return The weight
   */
  Double_t LookupWeight(AliAODTracklet* tracklet,
			Double_t        cent,
			Double_t        ipz,
			TH2*            corr=0) const
  {
    if (!CheckTracklet(tracklet)) return 1;
    Double_t w = CalcWeight(tracklet, cent, ipz, corr);
    if (fInverse) w = 1/w;
    return w;
  }
  /** 
   * Look-up weight of a particle 
   * 
   * @param particle Particle 
   * @param cent     Centrality 
   * @param ipz      Interaction point z-coordinate 
   * 
   * @return The weight
   */
  Double_t LookupWeight(TParticle* particle,
			Double_t   cent,
			Double_t   ipz) const
  {
    Double_t w = CalcWeight(particle, cent, ipz);
    if (fInverse) w = 1/w;
    return w;
  }
  /** 
   * Calculate the weight of a tracklet.  This member function must be
   * overloaded.
   * 
   * @param tracklet Tracklet 
   * @param cent     Centrality 
   * @param ipz      Interaction point Z coordinate 
   * @param corr     Optional histogram to fill with correlation of weights
   * 
   * @return The weight
   */
  virtual Double_t CalcWeight(AliAODTracklet* tracklet,
			      Double_t        cent,
			      Double_t        ipz,
			      TH2*            corr) const = 0;
  /** 
   * Calculate weight of a particle 
   * 
   * @param particle Particle 
   * @param cent     Centrality 
   * @param ipZ      Interaction point 
   * 
   * @return The weight
   */
  virtual Double_t CalcWeight(TParticle* particle,
			      Double_t   cent,
			      Double_t   ipZ) const = 0;
  /** 
   * Store values 
   * 
   * @param parent Parent container  
   */
  virtual TCollection* Store(TCollection* parent)
  {
    TList* top = new TList;
    top->SetName(GetName());
    top->SetOwner(true);
    parent->Add(top);
    top->Add(new TParameter<int> ("mask", fMask,    'f'));
    top->Add(new TParameter<int> ("veto", fVeto,    'f'));
    top->Add(new TParameter<int> ("calc", fCalc,    'f'));
    top->Add(new TParameter<bool>("inv",  fInverse, 'f'));
    return top;
  }
  /** 
   * Retrieve weights from a collection 
   * 
   * @param in Input collection 
   * 
   * @return Container read from or null
   */
  virtual TCollection* Retrieve(TCollection* in)
  {
    TCollection* top = static_cast<TCollection*>(in->FindObject(GetName()));
    if (!top) {
      Warning("Retrieve",
	      "Collection %s not found in %s", GetName(), in->GetName());
      in->ls();
      return 0;
    }
    fCalc   = GetPar<int>("calc", top);  
    fMask   = GetPar<int>("mask", top);
    fVeto   = GetPar<int>("veto", top);
    return top;
  }
  /** 
   * @return always true 
   */
  Bool_t IsFolder() const { return true; }  
  /** 
   * Print information on reweights 
   *
   * @param option Not used 
   */
  virtual void Print(Option_t* option="") const //*MENU*
  {
    gROOT->IndentLevel();
    Printf("%s : %s", ClassName(), GetName());
    Printf(" Weight calculation:  %s",
	   fCalc == kProduct ? "product" :
	   fCalc == kSquare  ? "square"  :
	   fCalc == kSum     ? "sum"     : "average");
    Printf(" Tracklet mask:       0x%02x", fMask);
    Printf(" Tracklet veto:       0x%02x", fVeto);
    Printf(" Take inverse:        %s", fInverse ? "yes" : "no");
  }
  UChar_t fCalc;        // Whether the square of the weight is calculated
  UChar_t fMask;        // Which partiles to take
  UChar_t fVeto;        // Which particles not to take
  Bool_t  fInverse;     // If true, do 1/w
  Int_t   fDebug;       // Debug level
  ClassDef(AliTrackletBaseWeights,1); // Base class of weights 
};

//====================================================================
/**
 * Simulation weights
 * 
 * @ingroup pwglf_forward_tracklets
 */
class AliTrackletPtPidStrWeights : public AliTrackletBaseWeights
{
public:
  /** 
   * Mode of weights 
   */
  enum EMode {
    /** If this bit is set, adjust weight up by error */    
    kUp       = (0x1) << 14,
    /** If this bit is set, adjust weight down by error */
    kDown     = (0x2) << 14,
    /** If this bit is set, the weight is disabled */
    kDisabled = (0x4) << 14
  };
  
  /** Map a particle species to a weight */
  typedef std::map<short,TH1D*> PdgMap;
  /**
   * Default constructor - ROOT I/O only
   * 
   */
  AliTrackletPtPidStrWeights()
    : AliTrackletBaseWeights(),
      fPt(0),
      fAbundance(),
      fStrangeness()
  {}
  /**
   * Named constructor 
   * 
   * @param name  Name of object 
   * @param title (optional) free-form title 
   */
  AliTrackletPtPidStrWeights(const char* name,
		     const char* title="Sim. tracklet weights");
  /**
   * Copy constructor 
   *
   * @param o Object to copy from  
   */
  AliTrackletPtPidStrWeights(const AliTrackletPtPidStrWeights& o);
  /**
   * Destructor 
   */
  virtual ~AliTrackletPtPidStrWeights() {}
  /**
   * Assignment operator  
   *
   * @param o Object to assign from  
   *
   * @return Reference to this object 
   */
  AliTrackletPtPidStrWeights& operator=(const AliTrackletPtPidStrWeights& o);

  /** 
   * Find the tracklet weight 
   * 
   * @param tracklet Tracklet 
   * @param cent     Centrality 
   * @param ipz      Interaction point Z coordinate 
   * @param corr     Optional histogram to fill with correlation of weights
   * 
   * @return The weight
   */
  virtual Double_t CalcWeight(AliAODTracklet* tracklet,
			      Double_t        cent,
			      Double_t        ipz,
			      TH2*            corr=0) const;
  virtual Double_t CalcWeight(TParticle* particle,
			      Double_t   cent,
			      Double_t   ipZ) const;
  /** 
   * Add a histogram to weight particle abundances 
   * 
   * @param pdg  Particle species 
   * @param h    Weight of particle specie as a function of centrality. 
   * @param mode Mode of this weight (+1: up, 0: normal, -1: down)
   * 
   * @return true on success. 
   */
  Bool_t AddAbundanceWeight(Short_t pdg, const TH1D* h, UShort_t mode=0)
  {
    AddPdgWeight(fAbundance, pdg, h, mode);
  }
  /** 
   * Add a histogram to weight strange particle abundances 
   * 
   * @param pdg  Particle species 
   * @param h    Weight of particle specie as a function of centrality. 
   * @param mode Mode of this weight (+1: up, 0: normal, -1: down)
   * 
   * @return true on success. 
   */
  Bool_t AddStrangenessWeight(Short_t pdg, const TH1D* h, UShort_t mode=0)
  {
    AddPdgWeight(fStrangeness, pdg, h, mode);
  }
  /** 
   * Set the weight histogram per transverse momentum and centrality 
   * 
   * @param h    histogram of (centrality vs pT) weights
   * @param mode Mode of this weight (+1: up, 0: normal, -1: down)
   *
   * @return true on success 
   */
  Bool_t SetPtWeight(const TH2D* h, UShort_t mode=0);
  /** 
   * Set the mode of a single particle abundance weight.  Mode is a
   * bit mask of 
   * 
   * - kDisabled:  The weight is disabled - i.e., will be 1 
   * - kUp:        The weight is adjusted up by the error of the weight 
   * - kDown:      The weight is adjusted down by the error of the weight 
   * 
   * @param pdg  Particle specie 
   * @param mode Mode of this weight (+1: up, 0: normal, -1: down)
   *
   */
  void SetAbundanceMode(Short_t pdg, UShort_t mode)
  {
    SetPdgMode(fAbundance, pdg, mode);
  }
  /** 
   * Set the mode of a single particle strangeness weight.  Mode is a
   * bit mask of 
   * 
   * - kDisabled:  The weight is disabled - i.e., will be 1 
   * - kUp:        The weight is adjusted up by the error of the weight 
   * - kDown:      The weight is adjusted down by the error of the weight 
   * 
   * @param pdg  Particle specie 
   * @param mode Mode 
   */
  void SetStrangenessMode(Short_t pdg, UShort_t mode)
  {
    SetPdgMode(fStrangeness, pdg, mode);
  }
  /** 
   * Set the mode of a transverse momentum weight.  Mode is a
   * bit mask of 
   * 
   * - kDisabled:  The weight is disabled - i.e., will be 1 
   * - kUp:        The weight for @f$ p_T<0.05@f$ is adjusted by +30% 
   * - kDown:      The weight for @f$ p_T<0.05@f$ is adjusted by -30% 
   * 
   * @param mode Mode 
   */
  void SetPtMode(UShort_t mode);
  /** 
   * Draw the weights 
   * 
   * @param option 
   */
  void Draw(Option_t* option=""); //*MENU*
  /** 
   * Print information on reweights 
   *
   * @param option Not used 
   */
  void Print(Option_t* option="") const; //*MENU*
  /** 
   * Store weights histograms in output of analysis 
   * 
   * @param out Collection to add histograms to 
   *
   * @return Output container 
   */
  TCollection* Store(TCollection* out);
  /** 
   * Retrieve weights from a collection 
   * 
   * @param in Input collection 
   * 
   * @return Container read from or null
   */
  TCollection* Retrieve(TCollection* in);
  /** 
   * Get the abundance weight of a given particle type for a given
   * centrality
   * 
   * @param apdg Absolute value of the particle PDG identifier 
   * @param cent Centrality 
   * 
   * @return The associated weight 
   */
  Double_t GetAbundanceWeight(UShort_t apdg, Double_t cent) const
  {
    return GetPdgWeight(fAbundance, apdg, cent);
  }
  /** 
   * Get the strangeness weight of a given particle type for a given
   * centrality
   * 
   * @param apdg Absolute value of the particle PDG identifier 
   * @param cent Centrality 
   * 
   * @return The associated weight 
   */
  Double_t GetStrangenessWeight(UShort_t apdg, Double_t cent) const
  {
    return GetPdgWeight(fStrangeness, apdg, cent);
  }

protected:  
  /** 
   * Look-up weight based on transverse momentum @a pT, particle
   * species @a pdf and event centrality @a cent.  This is used to
   * calculate the weight of a cluster.
   * 
   * @param pT    (mother) particle transverse momentum 
   * @param pdg   (mother) particle specie number 
   * @param cent  Event centrality 
   * 
   * @return The accumulated weight 
   */
  virtual Double_t CalcWeight(Double_t pT, Short_t pdg, Double_t cent) const;
  /** 
   * Get the histogram associated with a particle species in a given map. 
   * 
   * @param pdg Particle species 
   * @param m   Map 
   * 
   * @return Pointer to histogram or null 
   */
  TH1D* GetPdgHist(const PdgMap&  m, Short_t pdg) const;
  /** 
   * Get a PDG depedent weight. 
   *
   * @param m    Map to look up in 
   * @param apdg The absolute value of the PDG code 
   * @param cent The event centrality 
   *
   * @return the weight associated with the particle 
   */
  virtual Double_t GetPdgWeight(const PdgMap& m,
				UShort_t      apdg,
				Double_t      cent) const;
  /** 
   * Add a weight histogram 
   * 
   * @param map  Map of weights 
   * @param pdg  Particle species 
   * @param w    Weight histogram 
   * @param mode Mode of this weight (+1: up, 0: normal, -1: down)
   * 
   * @return true on success 
   */
  Bool_t AddPdgWeight(PdgMap& map, Short_t pdg, const TH1D* w, UShort_t mode=0);
  /** 
   * Set the mode of a single particle weight.  Mode is a bit mask of 
   * 
   * - kDisabled:  The weight is disabled - i.e., will be 1 
   * - kUp:        The weight is adjusted up by the error of the weight 
   * - kDown:      The weight is adjusted down by the error of the weight 
   * 
   * @param map  Map 
   * @param pdg  Particle specie 
   * @param mode Mode 
   */
  void SetPdgMode(PdgMap& map, Short_t pdg, UShort_t mode);
  /** 
   * Write a maps histograms to an output container 
   * 
   * @param parent Output container 
   * @param name   Name of output sub-container 
   * @param m      Map of weights 
   */
  void StoreMap(TCollection* parent, const char* name, PdgMap& m);
  /** 
   * Retrieve map contents from an input collection 
   * 
   * @param parent Input collection 
   * @param name   Name of collection 
   * @param m      Map to add to 
   * 
   * @return true on success 
   */
  Bool_t RetrieveMap(TCollection* parent, const char* name, PdgMap& m);
  /** 
   * Modify a drawn stack 
   * 
   * @param stack Stack to modify 
   */
  void ModStack(THStack* stack);
  /** 
   * Print a map 
   *
   * @param m       The map 
   * @param name    Name of map
   * @param options Options 
   */
  void PrintMap(const PdgMap& m, const char* name, Option_t* options="") const;
  /** 
   * Print info on a histogram
   * 
   * @param tag Tag line 
   * @param h   histogramT
   */
  void PrintHist(const char* tag, TH1* h) const;
  
  TH2D*   fPt;          // Weight by centrality and pT
  PdgMap  fAbundance;   // Map for abundance weight 
  PdgMap  fStrangeness; // Map for strangeness weight 
  ClassDef(AliTrackletPtPidStrWeights,3); // Weighs for tracklet analysis 
};

//____________________________________________________________________
AliTrackletPtPidStrWeights::AliTrackletPtPidStrWeights(const char* name,
						       const char* title)
  : AliTrackletBaseWeights(name, title),
    fPt(0),
    fAbundance(),
    fStrangeness()
{}

//____________________________________________________________________
AliTrackletPtPidStrWeights::
AliTrackletPtPidStrWeights(const AliTrackletPtPidStrWeights& o)
  : AliTrackletBaseWeights(o),
    fPt(0),
    fAbundance(),
    fStrangeness()
{
  UInt_t mask = kUp|kDown|kDisabled;
  SetPtWeight(o.fPt, o.fPt->TestBits(mask));
  for (PdgMap::const_iterator i = o.fAbundance.begin();
       i != o.fAbundance.end(); ++i) 
    AddAbundanceWeight(i->first, i->second, i->second->TestBits(mask));
  for (PdgMap::const_iterator i = o.fStrangeness.begin();
       i != o.fStrangeness.end(); ++i) 
    AddStrangenessWeight(i->first, i->second, i->second->TestBits(mask));
}
//____________________________________________________________________
AliTrackletPtPidStrWeights&
AliTrackletPtPidStrWeights::operator=(const AliTrackletPtPidStrWeights& o)
{
  if (&o == this) return *this;
  AliTrackletBaseWeights::operator=(o);
  fName    = o.fName;
  fTitle   = o.fTitle;
  UInt_t mask = kUp|kDown|kDisabled;
  SetPtWeight(o.fPt, o.fPt->TestBits(mask));
  for (PdgMap::const_iterator i = o.fAbundance.begin();
       i != o.fAbundance.end(); ++i) 
    AddAbundanceWeight(i->first, i->second, i->second->TestBits(mask));
  for (PdgMap::const_iterator i = o.fStrangeness.begin();
       i != o.fStrangeness.end(); ++i) 
    AddStrangenessWeight(i->first, i->second, i->second->TestBits(mask));
  return *this;
}
//____________________________________________________________________
Bool_t AliTrackletPtPidStrWeights::AddPdgWeight(PdgMap&     m,
						Short_t     pdg,
						const TH1D* w,
						UShort_t    mode)
{
  UShort_t         apdg = TMath::Abs(pdg);
  TH1D*            copy = static_cast<TH1D*>(w->Clone(Form("w%d", apdg)));
  copy->SetDirectory(0);
  copy->SetXTitle("Centrality [%]");
  copy->SetBit(mode);
  PdgMap::iterator i    = m.find(apdg);
  if (i != m.end()) {
    Warning("AddPdgWeight", "Replacing weight for %d", pdg);
    delete i->second;
    i->second = 0;
    m.erase(i);
  }
  m[apdg] = copy;
  return true;
}

//____________________________________________________________________
Bool_t AliTrackletPtPidStrWeights::SetPtWeight(const TH2D* h, UShort_t mode)
{
  if (fPt) {
    delete fPt;
    fPt = 0;
  }
  if (!h) return false;
  fPt = static_cast<TH2D*>(h->Clone("centPt"));
  fPt->SetDirectory(0);
  fPt->SetXTitle("Centrality [%]");
  fPt->SetYTitle("#it{p}_{T} [GeV]");
  fPt->SetBit(mode);
  return true;
}
//____________________________________________________________________
void AliTrackletPtPidStrWeights::SetPtMode(UShort_t mode)
{
  if (fPt) fPt->SetBit(mode);
}
//____________________________________________________________________
void AliTrackletPtPidStrWeights::SetPdgMode(PdgMap&  m,
					    Short_t  pdg,
					    UShort_t mode)
{
  TH1* h = GetPdgHist(m, pdg);
  if (h) h->SetBit(mode);
}
//____________________________________________________________________
TH1D* AliTrackletPtPidStrWeights::GetPdgHist(const PdgMap&  m, Short_t pdg) const
{
  UShort_t apdg = TMath::Abs(pdg);
  PdgMap::const_iterator i = m.find(apdg);
  if (i == m.end()) return 0;
  return i->second;
}

//____________________________________________________________________
Double_t AliTrackletPtPidStrWeights::GetPdgWeight(const PdgMap&  m,
						  UShort_t       apdg,
						  Double_t       cent) const
{
  if (m.size() < 1) return 1;
  TH1D* h = GetPdgHist(m, apdg);
  if (!h || h->TestBit(kDisabled)) return 1;

#if 1
  Int_t bC = h->GetXaxis()->FindBin(cent);
  if (bC < 1 || bC > h->GetNbinsX()) return 1;
#else
  Int_t bC = 1;
#endif
  Double_t add = (h->TestBit(kUp)   ? +1 :
		  h->TestBit(kDown) ? -1 : 0) * h->GetBinError(bC);  
  return h->GetBinContent(bC) + add;
}
//____________________________________________________________________
Double_t
AliTrackletPtPidStrWeights::CalcWeight(Double_t pT,
				       Short_t  pdg,
				       Double_t cent) const
{
  Double_t w = 1;
  if (fPt && !fPt->TestBit(kDisabled)) {
    Int_t    bC  =  fPt->GetXaxis()->FindBin(cent);
    Int_t    bpT =  fPt->GetYaxis()->FindBin(pT);
    if (bC  >= 1 && bC  <= fPt->GetXaxis()->GetNbins() &&
	bpT >= 1 && bpT <= fPt->GetYaxis()->GetNbins()) {
      Double_t fac = (pT >= 0.05 ?  1 :
		      (fPt->TestBit(kUp)   ? 1.3 :
		       fPt->TestBit(kDown) ? 0.7 : 1));
      Double_t ptW =  fPt->GetBinContent(bC, bpT);
      w            *= fac*ptW;
      if (fDebug > 3)
	Info("CalcWeight", "pT=%5.2f -> %4.2f * %6.4f = %6.4f",
	     pT, fac, ptW, fac*ptW);
    }
  }
  UShort_t apdg = TMath::Abs(pdg);
  Double_t aW   = GetPdgWeight(fAbundance,   apdg, cent);
  Double_t sW   = GetPdgWeight(fStrangeness, apdg, cent);
  w *= aW * sW;
  if (fDebug > 3) 
    Info("CalcWeight","pdg=%4d -> %6.4f * %6.4f = %6.4f -> %6.4f",
	 pdg, aW, sW, aW*sW, w);
  // Printf("Weight of pT=%6.3f pdg=%5d cent=%5.1f -> %f", pT, pdg, cent, w);
  return w;
}

//____________________________________________________________________
Double_t
AliTrackletPtPidStrWeights::CalcWeight(TParticle* particle,
				       Double_t   cent,
				       Double_t   ipz) const
{
  Int_t    pdg = particle->GetPdgCode();
  Double_t pT  = particle->Pt();
  return CalcWeight(pT, pdg, cent);
}

//____________________________________________________________________
Double_t
AliTrackletPtPidStrWeights::CalcWeight(AliAODTracklet* tracklet,
				       Double_t        cent,
				       Double_t        ipz,
				       TH2*            corr) const
{
#if 0
  if (!tracklet->IsSimulated()) {
    Warning("CalcWeight", "Not a simulated tracklet");
    return 1;
  }
#endif
  Double_t w1 = 1, w2 = 1;
  
  // AliAODMCTracklet* mc = static_cast<AliAODMCTracklet*>(tracklet);
  AliAODTracklet* mc = tracklet;
  Short_t pdg1 = mc->GetParentPdg();
  Short_t pdg2 = mc->GetParentPdg(true);
  if      (pdg1!=0)         w1 = CalcWeight(mc->GetParentPt(),    pdg1,cent);
  if      (pdg2!=0)         w2 = CalcWeight(mc->GetParentPt(true),pdg2,cent);
  else if (fCalc!=kProduct) w2 = w1;

  if (corr && mc->IsMeasured()) corr->Fill(w1, w2);

  Double_t    ret = 1;
  const char* cm  = "?";
  switch (fCalc) {
  case kProduct: ret = w1 * w2;              cm="product"; break;
  case kSquare:  ret = TMath::Sqrt(w1 * w2); cm="square";  break;
  case kSum:     ret = 1+(w1-1)+(w2-1);      cm="sum";     break;
  case kAverage: ret = 1+((w1-1)+(w2-1))/2;  cm="average"; break;
  }

  if (fDebug > 1)
    Printf("pdg1=%5d -> %6.4f  pdg2=%5d -> %6.4f  (%10s) -> %6.4f",
	   pdg1, w1, pdg2, w2, cm, ret);
  return ret;
}

//____________________________________________________________________
void AliTrackletPtPidStrWeights::ModStack(THStack* stack)
{
  if (!stack || !stack->GetHists()) return;
  TIter    next(stack->GetHists());
  TH1*     hist = 0;
  Color_t  colors[] = { kPink+2, kRed+2, kOrange+2, kYellow+2,
			kSpring+2, kGreen+2, kTeal+2, kCyan+2,
			kBlue+2, kViolet+2 };
  Int_t    idx      = 0;
  while ((hist = static_cast<TH1*>(next()))) {
    Color_t c = colors[idx % 10];
    idx++;
    hist->SetLineColor(c);
    hist->SetMarkerColor(c);
    hist->SetMarkerStyle(20+idx%20);
  }
}

//____________________________________________________________________
TCollection*
AliTrackletPtPidStrWeights::Store(TCollection* parent)
{
  TCollection* top = AliTrackletBaseWeights::Store(parent);
  if (!top) return 0;
	   
  if (fPt) {
    TH2* copy = static_cast<TH2*>(fPt->Clone("centPt"));
    copy->SetDirectory(0);
    copy->SetBinContent(0,0,1); // For counting merges
    top->Add(copy);
  }
  StoreMap(top, "abundance",   fAbundance);
  StoreMap(top, "strangeness", fStrangeness);
  return top;
}

//____________________________________________________________________
void
AliTrackletPtPidStrWeights::StoreMap(TCollection* parent,
				     const char*  name,
				     PdgMap&      m)
{
  TList* top = new TList;
  top->SetName(name);
  top->SetOwner(true);
  parent->Add(top);
  for (PdgMap::const_iterator i = m.begin(); i != m.end(); ++i) {
    TH1* copy = static_cast<TH1*>(i->second->Clone(Form("w%d",i->first)));
    copy->SetDirectory(0);
    copy->SetBinContent(0,1); // For counting merges
    top->Add(copy);
  }
}


//____________________________________________________________________
TCollection*
AliTrackletPtPidStrWeights::Retrieve(TCollection* parent)
{
  TCollection* top = AliTrackletBaseWeights::Retrieve(parent);
  if (!top) return 0;
  
  fPt     = static_cast<TH2D*>(top->FindObject("centPt"));
  if (!fPt)
    Warning("Retrieve","centPt histogram not found in %s", GetName());
  else {
    fPt->SetDirectory(0);
    Double_t scale = fPt->GetBinContent(0,0);
    fPt->Scale(1/scale); // Counting merges
    // fPt->SetBinContent(0,0,1); // Zero merger count 
  }
  if (!RetrieveMap(top, "abundance",   fAbundance))   return 0;
  if (!RetrieveMap(top, "strangeness", fStrangeness)) return 0;

  return top;
}
//____________________________________________________________________
Bool_t
AliTrackletPtPidStrWeights::RetrieveMap(TCollection* parent,
					const char* name,
					PdgMap& m)
{
  m.clear();
  TList* top = static_cast<TList*>(parent->FindObject(name));
  if (!top) {
    Warning("RetrieveMap",
	    "Collection %s not found in %s", name, parent->GetName());
    return true;
  }
  TIter    next(top);
  TObject* o = 0;
  while ((o = next())) {
    if (!o->IsA()->InheritsFrom(TH1D::Class())) continue;
    TH1D* copy = static_cast<TH1D*>(o);
    copy->SetDirectory(0);
    Double_t scale = copy->GetBinContent(0); 
    copy->Scale(1/scale); // Counting merges
    // copy->SetBinContent(0,1); // Zero merger count 
    TString nme(copy->GetName());
    nme.Remove(0,1);
    Int_t apdg = nme.Atoi();
    AddPdgWeight(m, apdg, copy, copy->TestBits(kUp|kDown|kDisabled));
  }
  return true;
}

//____________________________________________________________________
void
AliTrackletPtPidStrWeights::Draw(Option_t* option)
{
  TVirtualPad* master = TVirtualPad::Pad();
  if (!master) {
    Warning("Draw", "No current pad to draw in");
    return;
  }

  Int_t nPad = 3;
  Int_t iPad = 0;
  if (fAbundance.size()   <= 0) nPad--;
  if (fStrangeness.size() <= 0) nPad--;
  TString opt(option); opt.ToUpper();
  
  master->Divide(nPad,1);  
  if (fPt) {
    THStack* stack = new THStack(fPt, (opt.Contains("C") ? "x" : "y"));
    ModStack(stack);
    stack->SetTitle("#it{p}_{T} weights");
    master->cd(++iPad);
    stack->Draw("nostack");
    stack->GetHistogram()->SetXTitle((opt.Contains("C") ?
				      "Centrality [%]" : "#it{p}_{T}"));
    master->GetPad(iPad)->BuildLegend();
    master->GetPad(iPad)->Modified();
  }

  if (fAbundance.size() > 0) {
    THStack* ha = new THStack("abundance", "Abundance weights");
    for (PdgMap::const_iterator i = fAbundance.begin();
	 i != fAbundance.end(); ++i) {
      ha->Add(i->second);
    }
    ModStack(ha);
    master->cd(++iPad);
    ha->Draw("nostack");
    ha->GetHistogram()->SetXTitle("Centrality [%]");
    master->GetPad(iPad)->BuildLegend();
  }

  if (fStrangeness.size() > 0) {
    THStack* hs = new THStack("strangeness", "Strangeness weights");
    for (PdgMap::const_iterator i = fStrangeness.begin();
	 i != fStrangeness.end(); ++i) {
      hs->Add(i->second);
    }
    ModStack(hs);
    master->cd(++iPad);
    hs->Draw("nostack");
    hs->GetHistogram()->SetXTitle("Centrality [%]");
    master->GetPad(iPad)->BuildLegend();
  }
  master->Modified();
}

//____________________________________________________________________
void
AliTrackletPtPidStrWeights::Print(Option_t* option) const
{
  AliTrackletBaseWeights::Print(option);
  gROOT->IncreaseDirLevel();
  PrintHist("pT", fPt);
  PrintMap(fAbundance,   "Abundance", option);
  PrintMap(fStrangeness, "Strangeness", option);
  
  gROOT->DecreaseDirLevel();
}
//____________________________________________________________________
void
AliTrackletPtPidStrWeights::PrintHist(const char* tag, TH1* h) const
{
  if (!h) return;
  gROOT->IndentLevel();
  Printf("%10s (%c): %p %s/%s", tag, 
	 h->TestBit(kDisabled) ? '0' :
	 h->TestBit(kUp)       ? '+' :
	 h->TestBit(kDown)     ? '-' : '=',
	 h, h->GetName(), h->GetTitle());

}
//____________________________________________________________________
void
AliTrackletPtPidStrWeights::PrintMap(const PdgMap& m,
				     const char*   name,
				     Option_t* option) const  
{
  gROOT->IndentLevel();
  Printf("Map of PDG codes: %s", name);
  gROOT->IncreaseDirLevel();
  for (PdgMap::const_iterator i = m.begin(); i != m.end(); ++i) 
    PrintHist(Form("%10d", i->first), i->second);

  gROOT->DecreaseDirLevel();
}


//====================================================================
/**
 * Base class for tracklets weights 
 */
class AliTrackletDeltaWeights : public AliTrackletBaseWeights
{
public:
  /**
   * Default constructor - ROOT I/O only
   * 
   */
  AliTrackletDeltaWeights()
    : AliTrackletBaseWeights(),
      fHistos(),
      fCentAxis()
  {}
  /**
   * Named constructor 
   * 
   * @param name  Name of object 
   * @param title (optional) free-form title 
   */
  AliTrackletDeltaWeights(const char* name,
			  const char* title="Sim. tracklet weights");
  /**
   * Copy constructor 
   *
   * @param o Object to copy from  
   */
  AliTrackletDeltaWeights(const AliTrackletDeltaWeights& o);
  /**
   * Destructor 
   */
  virtual ~AliTrackletDeltaWeights() {}

  /**
   * Assignment operator  
   *
   * @param o Object to assign from  
   *
   * @return Reference to this object 
   */
  AliTrackletDeltaWeights& operator=(const AliTrackletDeltaWeights& o);
  /**
   * Set the centrality axis 
   * 
   * @param axis Axis to copy
   */
  void SetCentAxis(const TAxis& axis);
  /** 
   * Set the centrality axis 
   * 
   * @param n     number of bins 
   * @param bins  Bin limits (n+1 entries)
   */
  void SetCentAxis(Int_t n, const Double_t* bins);
  /** 
   * Set the centrality axis 
   * 
   * @param n    number of bins 
   * @param low  Least value 
   * @param high Largest value
   */
  void SetCentAxis(Int_t n, Double_t low, Double_t high);
  /** 
   * Set the histogram to use for the given centrality bin 
   * 
   * @param bin Bin number.  Start at 1
   * @param h   Histogram.  The histogram will be cloned
   * 
   * @return true on success 
   */
  Bool_t SetHisto(Int_t bin, TH1* h);
  /** 
   * Find the histogram corresponding to a centrality 
   * 
   * @param cent Centrality 
   * 
   * @return Pointer to histogram or null
   */
  TH1* FindHisto(Double_t cent) const;
  /** 
   * Find the bin on an axis that correspnd to the passed value.  If
   * the value is outside the defined range, the closest bin number is
   * returned.
   * 
   * @param value  Value 
   * @param axis   Axis
   * 
   * @return Bin number 
   */
  static Int_t FindBin(Double_t value, const TAxis* axis);
  /** 
   * Find the tracklet weight 
   * 
   * @param tracklet Tracklet 
   * @param cent     Centrality 
   * @param ipZ      Interaction point Z coordinate 
   * @param corr     Optional histogram to fill with correlation of weights
   * 
   * @return The weight
   */
  virtual Double_t CalcWeight(AliAODTracklet* tracklet,
			      Double_t        cent,
			      Double_t        ipZ,
			      TH2*            corr=0) const;
  virtual Double_t CalcWeight(TParticle* particle,
			      Double_t   cent,
			      Double_t   ipZ) const;
  /** 
   * Project 3D histogram on 1 axis 
   * 
   * @param h 
   * @param which 
   * @param name 
   * 
   * @return 
   */
  TH1* Project(TH1* h, char which, const char* name);
  /** 
   * Draw the weights 
   * 
   * @param option 
   */
  void Draw(Option_t* option=""); //*MENU*
  /** 
   * Print information on reweights 
   *
   * @param option Not used 
   */
  void Print(Option_t* option="") const; //*MENU*
  /** 
   * Draw one stack 
   * 
   * @param pad 
   * @param nPad
   * @param stack 
   * @param xtitle 
   * @param min 
   * @param max 
   */
  void DrawStack(TVirtualPad* pad,
		 Int_t        nPad, 
		 THStack*     stack,
		 const char*  xtitle,
		 Double_t     min=0,
		 Double_t     max=2.5);
  void DrawOne(Double_t cent=2.5,
	       Double_t eta=0,
	       Double_t ipz=0) const; //*MENU*
/** List of histograms */
  TObjArray fHistos;
  /** Centrality axis */
  TAxis fCentAxis;
  ClassDef(AliTrackletDeltaWeights,1); // Base class of weights 
};

//____________________________________________________________________
AliTrackletDeltaWeights::AliTrackletDeltaWeights(const char* name,
						 const char* title)
  : AliTrackletBaseWeights(name, title),
    fHistos(),
    fCentAxis()
{
  fHistos.SetOwner(true);
}
//____________________________________________________________________
AliTrackletDeltaWeights::AliTrackletDeltaWeights(const AliTrackletDeltaWeights& o)
  : AliTrackletBaseWeights(o),
    fHistos(),
    fCentAxis(o.fCentAxis)
{
  fHistos.SetOwner(true);
  for (Int_t i = 0; i < o.fHistos.GetEntriesFast(); i++) {
    TH1*  h = static_cast<TH1*>(o.fHistos.At(i));
    if (!h) continue;
    h = static_cast<TH1*>(h->Clone());
    h->SetDirectory(0);
    fHistos.Add(h);
  }    
}
//____________________________________________________________________
AliTrackletDeltaWeights&
AliTrackletDeltaWeights::operator=(const AliTrackletDeltaWeights& o)
{
  if (&o == this) return *this;
  AliTrackletBaseWeights::operator=(o);
  fCentAxis = o.fCentAxis;
  for (Int_t i = 0; i < o.fHistos.GetEntriesFast(); i++) {
    TH1*  h = static_cast<TH1*>(o.fHistos.At(i));
    if (!h) continue;
    h = static_cast<TH1*>(h->Clone());
    h->SetDirectory(0);
    fHistos.Add(h);
  }
  return *this;
}

//____________________________________________________________________
void
AliTrackletDeltaWeights::SetCentAxis(const TAxis& a)
{
  if (a.GetNbins() && a.GetXbins()->GetArray())
    SetCentAxis(a.GetNbins(), a.GetXbins()->GetArray());
  else
    SetCentAxis(a.GetNbins(), a.GetXmin(), a.GetXmax());
}

//____________________________________________________________________
void
AliTrackletDeltaWeights::SetCentAxis(Int_t n, const Double_t* bins)
{
  fCentAxis.Set(n, bins);
}
//____________________________________________________________________
void
AliTrackletDeltaWeights::SetCentAxis(Int_t n, Double_t l, Double_t h)
{
  fCentAxis.Set(n, l, h);
}

//____________________________________________________________________
Bool_t
AliTrackletDeltaWeights::SetHisto(Int_t bin, TH1* h)
{
  if (bin < 1 || bin > fCentAxis.GetNbins()) {
    Warning("SetHisto", "Centrality bin %d out of range [%d,%d]",
	    bin, 1, fCentAxis.GetNbins());
    return false;
  }
  Double_t c1 = fCentAxis.GetBinLowEdge(bin);
  Double_t c2 = fCentAxis.GetBinUpEdge(bin);  
  TString name;
  name.Form("cent%03dd%02d_%03dd%02d",
	    Int_t(c1), Int_t(c1 * 100) % 100,
	    Int_t(c2), Int_t(c2 * 100) % 100);
  TH1* cpy = static_cast<TH1*>(h->Clone(name));
  cpy->SetDirectory(0);
  cpy->SetTitle(Form("%6.2f%% - %6.2f%%", c1, c2));
  fHistos.AddAtAndExpand(cpy, bin-1);
  return true;
}
  
//____________________________________________________________________
TH1*
AliTrackletDeltaWeights::FindHisto(Double_t cent) const
{
  Int_t bin = FindBin(cent, &fCentAxis);
  if       (bin >= fHistos.GetEntriesFast()) return 0;
  return static_cast<TH1*>(fHistos.At(bin-1));
}

//____________________________________________________________________
Int_t
AliTrackletDeltaWeights::FindBin(Double_t value,
				 const TAxis*   axis) 
{
  if (!axis) return 1;
  Int_t bin = const_cast<TAxis*>(axis)->FindBin(value);
  if       (bin <  1)                bin = 1;
  else if  (bin >  axis->GetNbins()) bin = axis->GetNbins();
  return bin;
}
  
//____________________________________________________________________
Double_t
AliTrackletDeltaWeights::CalcWeight(TParticle* particle,
				    Double_t   cent,
				    Double_t   ipZ) const
{
  // Perhaps this method doesn't really make much sense 
  TH1* h = FindHisto(cent);
  if (!h) return 1;

  Double_t eta      = particle->Eta();
  Double_t delta    = 0;
  Int_t    etaBin   = FindBin(eta,   h->GetXaxis());
  Int_t    deltaBin = FindBin(delta, h->GetYaxis());
  Int_t    ipzBin   = FindBin(ipZ,   h->GetZaxis());

  return h->GetBinContent(etaBin, deltaBin, ipzBin);
}  
				   
//____________________________________________________________________
Double_t
AliTrackletDeltaWeights::CalcWeight(AliAODTracklet* tracklet,
				    Double_t        cent,
				    Double_t        ipZ,
				    TH2*            corr) const
{
  if (tracklet->IsGenerated()) return 1;
  TH1* h = FindHisto(cent);
  if (!h) return 1;

  Double_t eta      = tracklet->GetEta();
  Double_t delta    = tracklet->GetDelta();
  Int_t    etaBin   = FindBin(eta,   h->GetXaxis());
  Int_t    deltaBin = FindBin(delta, h->GetYaxis());
  Int_t    ipzBin   = FindBin(ipZ,   h->GetZaxis());

  return h->GetBinContent(etaBin, deltaBin, ipzBin);
}

//____________________________________________________________________
TH1*
AliTrackletDeltaWeights::Project(TH1* h, char which, const char* name)
{
  TString n; n.Form("%s_%s",h->GetName(),name);
  TString t(h->GetTitle());
  TH1*    ret = 0;
  TAxis*  a0  = 0;
  TAxis* a1  = 0;
  TAxis* a2  = 0;
  typedef
    TH1D* (TH3::*ProjFunc)(const char*,Int_t,Int_t,Int_t,Int_t,Option_t*) const;
  ProjFunc func = 0;
  switch (which) {
  case 'x':
    a0   = h->GetXaxis();
    a1   = h->GetYaxis();
    a2   = h->GetZaxis();
    func = &TH3::ProjectionX;
    break;
  case 'y':
    a0   = h->GetYaxis();
    a1   = h->GetXaxis();
    a2   = h->GetZaxis();
    func = &TH3::ProjectionY;
    break;
  case 'z':
    a0   = h->GetZaxis();
    a1   = h->GetXaxis();
    a2   = h->GetYaxis();
    func = &TH3::ProjectionZ;
    break;
  }
  Int_t   nx  = a0->GetNbins();
  if (nx <= 1) return 0;
  if (a0->GetXbins() && a0->GetXbins()->GetArray())
    ret = new TH1D(n,t,nx,a0->GetXbins()->GetArray());
  else
    ret = new TH1D(n,t,nx,a0->GetXmin(),a0->GetXmax());
  ret->Sumw2();
  ret->SetXTitle(a0->GetTitle());
  static_cast<const TAttAxis*>(a0)->Copy(*(ret->GetXaxis()));
  ret->SetDirectory(0);
  ret->SetLineColor(h->GetMarkerColor());
  ret->SetMarkerColor(h->GetMarkerColor());
  ret->SetFillColor(kWhite);// color);
  ret->SetFillStyle(0);
  ret->SetMarkerStyle(20);
  
  for (Int_t i0 = 1; i0 <= ret->GetNbinsX(); i0++) {
    Int_t    count = 0;
    Double_t sum   = 0;
    Double_t sumE2 = 0;
    for (Int_t i1 = 1; i1 <= a1->GetNbins(); i1++) {
      for (Int_t i2 = 1; i2 <= a2->GetNbins(); i2++) {
	Int_t bin = 0;
	switch (which) {
	case 'x': bin = h->GetBin(i0, i1, i1); break;
	case 'y': bin = h->GetBin(i1, i0, i2); break;
	case 'z': bin = h->GetBin(i1, i2, i0); break;
	}
	Double_t c = h->GetBinContent(bin);
	Double_t e = h->GetBinError  (bin);
	if (c < 1e-3 || e < 1e-6) continue;
	count++;
	sum   += c;
	sumE2 += e*e;
      }
    }
    // ret->SetBinContent(i0, ret->GetBinContent(i0)/count);
    // ret->SetBinError  (i0, ret->GetBinError  (i0)/count);
    ret->SetBinContent(i0, sum/count);
    ret->SetBinError  (i0, TMath::Sqrt(sumE2)/count);
  }
	  
      
  // ret->Scale(1./(a1->GetNbins()*a2->GetNbins()));
  return ret;
}
  
//____________________________________________________________________
void
AliTrackletDeltaWeights::Print(Option_t* option) const
{
  AliTrackletBaseWeights::Print(option);
  gROOT->IncreaseDirLevel();
  gROOT->IndentLevel();
  Printf("%d centrality bins", fCentAxis.GetNbins());
  for (Int_t i = 0; i < fHistos.GetEntriesFast(); i++) {
    TH1*  h = static_cast<TH1*>(fHistos.At(i));
    if (!h) continue;
    h->Print(option);
  }    
  gROOT->DecreaseDirLevel();
}
//____________________________________________________________________
void
AliTrackletDeltaWeights::DrawOne(Double_t cent,
				 Double_t eta,
				 Double_t ipz) const
{
  TH1* h = FindHisto(cent);
  if (!h) return;

  Int_t etaBin = FindBin(eta, h->GetXaxis());
  Int_t ipzBin = FindBin(ipz, h->GetZaxis());

  TH1* tmp = 0;
  if (h->IsA()->InheritsFrom(TH2::Class()))
    tmp = static_cast<TH2*>(h)->ProjectionY("tmp", etaBin, etaBin);
  if (h->IsA()->InheritsFrom(TH3::Class()))
    tmp = static_cast<TH3*>(h)->ProjectionY("tmp", etaBin, etaBin);
  else
    tmp = static_cast<TH1*>(h->Clone("tmp"));
  if (!tmp) return;
  tmp->SetDirectory(0);
  tmp->Draw();
}
				 
//____________________________________________________________________
void
AliTrackletDeltaWeights::DrawStack(TVirtualPad* pad,
				   Int_t        nPad,
				   THStack*     stack,
				   const char*  xtitle,
				   Double_t     min,
				   Double_t     max)
{
  if (!stack->GetHists()) return;
  stack->SetMinimum(min);
  stack->SetMaximum(max);
  pad->SetTicks();
  pad->SetGridy();
  pad->cd();
  stack->Draw("nostack");
  TH1* h = stack->GetHistogram();
  h->SetXTitle(xtitle);
  h->SetYTitle(stack->GetTitle());
  h->GetXaxis()->SetNdivisions(205);
  h->GetYaxis()->SetNdivisions(205);
  h->GetXaxis()->SetLabelOffset(0.005*pad->GetWNDC());
  h->GetXaxis()->SetTitleOffset(1.5*pad->GetWNDC());
  h->GetXaxis()->SetLabelSize(0.03/pad->GetWNDC()/*nPad*/);
  h->GetYaxis()->SetLabelSize(0.03/pad->GetWNDC()/*nPad*/);
  h->GetXaxis()->SetTitleSize(0.03/pad->GetWNDC()/*nPad*/);
  h->GetYaxis()->SetTitleSize(0.03/pad->GetWNDC()/*nPad*/);
  pad->Modified();
}
  
//____________________________________________________________________
void
AliTrackletDeltaWeights::Draw(Option_t* option)
{
  TVirtualPad* master = TVirtualPad::Pad();
  if (!master) {
    Warning("Draw", "No current pad to draw in");
    return;
  }

  TString  opt(option); opt.ToUpper();
  THStack* stackEta   = new THStack("eta",   "#LTw#GT_{#Delta,IP_{z}}");
  THStack* stackDelta = new THStack("delta", "#LTw#GT_{#eta,IP_{z}}");
  THStack* stackIPz   = new THStack("ipz",   "#LTw#GT_{#eta,#Delta}");
  
  for (Int_t i = 0; i < fHistos.GetEntriesFast(); i++) {
    TH3*  h = static_cast<TH3*>(fHistos.At(i));
    if (!h) continue;

    TH1* etaProj   = Project(h, 'x', "eta");
    TH1* deltaProj = Project(h, 'y', "delta");
    TH1* ipzProj   = Project(h, 'z', "ipz");
    stackEta  ->Add(etaProj);
    stackDelta->Add(deltaProj);
    stackIPz  ->Add(ipzProj);
  }

  Int_t   nPad = 0;
  if (stackEta  ->GetHists()) nPad++;
  if (stackDelta->GetHists()) nPad++;
  if (stackIPz  ->GetHists()) nPad++;
  
  master->SetTopMargin(0.01);
  master->SetRightMargin(0.01);
  master->SetLeftMargin(0.07*nPad);
  if (opt.Contains("V"))
    master->Divide(1,nPad);
  else
    master->Divide(nPad,1, 0, 0);

  if (nPad >= 1) DrawStack(master->cd(1), nPad, stackEta,   "#eta");
  if (nPad >= 2) DrawStack(master->cd(2), nPad, stackDelta, "#Delta");
  if (nPad >= 3) DrawStack(master->cd(3), nPad, stackIPz,   "IP_{z}");
  master->GetPad(nPad)->SetRightMargin(0.01);
  
  master->Modified();
}

#endif
//____________________________________________________________________
// 
// EOF
// 

