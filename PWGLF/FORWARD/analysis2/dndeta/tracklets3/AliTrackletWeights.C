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
# include <TList.h>
# include <TBrowser.h>
# include <TCanvas.h>
# include <THStack.h>
# include <TROOT.h>
# include <TClass.h>
# include "AliAODTracklet.C"
# include <TParameter.h>
#else
class TH1D;
class TH2D;
class TH2;
class TList;
class TBrowser;
class THStack;
class AliAODTracklet;
class TCanvas; // Autoload 
#endif

/**
 * Simulation weights
 * 
 * @ingroup pwglf_forward_tracklets
 */
class AliTrackletWeights : public TNamed
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
  
  /** Map a particle species to a weight */
  typedef std::map<short,TH1D*> PdgMap;
  /**
   * Default constructor - ROOT I/O only
   * 
   */
  AliTrackletWeights()
    : TNamed(),
      fPt(0),
      fAbundance(),
      fStrangeness(),
      fCalc(kProduct)
  {}
  /**
   * Named constructor 
   * 
   * @param name  Name of object 
   * @param title (optional) free-form title 
   */
  AliTrackletWeights(const char* name,
		     const char* title="Sim. tracklet weights");
  /**
   * Copy constructor 
   *
   * @param o Object to copy from  
   */
  AliTrackletWeights(const AliTrackletWeights& o);
  /**
   * Destructor 
   */
  virtual ~AliTrackletWeights() {}
  /**
   * Assignment operator  
   *
   * @param o Object to assign from  
   *
   * @return Reference to this object 
   */
  AliTrackletWeights& operator=(const AliTrackletWeights& o);

  /** 
   * Find the tracklet weight 
   * 
   * @param tracklet Tracklet 
   * @param cent     Centrality 
   * @param corr     Optional histogram to fill with correlation of weights
   * 
   * @return The weight
   */
  virtual Double_t LookupWeight(AliAODTracklet* tracklet,
				Double_t        cent,
				TH2*            corr=0) const;
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
   * @return always true 
   */
  Bool_t IsFolder() const { return true; }
  /** 
   * Store weights histograms in output of analysis 
   * 
   * @param out Collection to add histograms to 
   */
  void Store(TCollection* out);
  /** 
   * Retrieve weights from a collection 
   * 
   * @param in Input collection 
   * 
   * @return true on success
   */
  Bool_t Retrieve(TCollection* in);
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
   * species @a pdf and event centrality @a cent.
   * 
   * @param pT    (mother) particle transverse momentum 
   * @param pdg   (mother) particle specie number 
   * @param cent  Event centrality 
   * 
   * @return The accumulated weight 
   */
  virtual Double_t LookupWeight(Double_t pT, Short_t pdg, Double_t cent) const;
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
  UChar_t fCalc;        // Whether the square of the weight is calculated
  UChar_t fMask;        // Which particles to take
  UChar_t fVeto;        // Which particles not to take 
  ClassDef(AliTrackletWeights,3); // Weighs for tracklet analysis 
};

//____________________________________________________________________
AliTrackletWeights::AliTrackletWeights(const char* name,
				       const char* title)
  : TNamed(name, title),
    fPt(0),
    fAbundance(),
    fStrangeness(),
    fCalc(kProduct),
    fMask(0xFF),
    fVeto(0x0)
{}

//____________________________________________________________________
AliTrackletWeights::AliTrackletWeights(const AliTrackletWeights& o)
  : TNamed(o),
    fPt(0),
    fAbundance(),
    fStrangeness(),
    fCalc(o.fCalc),
    fMask(o.fMask),
    fVeto(o.fVeto)
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
AliTrackletWeights& AliTrackletWeights::operator=(const AliTrackletWeights& o)
{
  if (&o == this) return *this;
  fName    = o.fName;
  fTitle   = o.fTitle;
  fCalc    = o.fCalc;
  fMask    = o.fMask;
  fVeto    = o.fVeto;
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
Bool_t AliTrackletWeights::AddPdgWeight(PdgMap&     m,
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
Bool_t AliTrackletWeights::SetPtWeight(const TH2D* h, UShort_t mode)
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
void AliTrackletWeights::SetPtMode(UShort_t mode)
{
  if (fPt) fPt->SetBit(mode);
}
//____________________________________________________________________
void AliTrackletWeights::SetPdgMode(PdgMap&  m,
				    Short_t  pdg,
				    UShort_t mode)
{
  TH1* h = GetPdgHist(m, pdg);
  if (h) h->SetBit(mode);
}
//____________________________________________________________________
TH1D* AliTrackletWeights::GetPdgHist(const PdgMap&  m, Short_t pdg) const
{
  UShort_t apdg = TMath::Abs(pdg);
  PdgMap::const_iterator i = m.find(apdg);
  if (i == m.end()) return 0;
  return i->second;
}

//____________________________________________________________________
Double_t AliTrackletWeights::GetPdgWeight(const PdgMap&  m,
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
AliTrackletWeights::LookupWeight(Double_t pT, Short_t pdg, Double_t cent) const
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
      w            *= fac*fPt->GetBinContent(bC, bpT);
    }
  }
  UShort_t apdg = TMath::Abs(pdg);

  w *= GetPdgWeight(fAbundance,   apdg, cent);
  w *= GetPdgWeight(fStrangeness, apdg, cent);
  // Printf("Weight of pT=%6.3f pdg=%5d cent=%5.1f -> %f", pT, pdg, cent, w);
  return w;
}

//____________________________________________________________________
Double_t
AliTrackletWeights::LookupWeight(AliAODTracklet* tracklet,
				 Double_t        cent,
				 TH2*            corr) const
{
#if 0
  if (!tracklet->IsSimulated()) {
    Warning("LookupWeight", "Not a simulated tracklet");
    return 1;
  }
#endif
  UChar_t flags = tracklet->GetFlags();
  if (fMask != 0xFF && (fMask & flags) == 0) {
    // Info("LookupWeight", "Tracklet 0x%02x does not fullfill mask 0x%02x",
    //       flags, fMask);	  
    return 1;
  }
  if ((fVeto & flags) != 0) {
    // Info("LookupWeight", "Tracklet 0x%02x vetoed by 0x%02x",flags, fVeto);
    return 1;
  }
  
  Double_t w1 = 1, w2 = 1;
  
  // AliAODMCTracklet* mc = static_cast<AliAODMCTracklet*>(tracklet);
  AliAODTracklet* mc = tracklet;
  Short_t pdg1 = mc->GetParentPdg();
  Short_t pdg2 = mc->GetParentPdg(true);
  if      (pdg1>0)          w1 = LookupWeight(mc->GetParentPt(),    pdg1,cent);
  if      (pdg2>0)          w2 = LookupWeight(mc->GetParentPt(true),pdg2,cent);
  else if (fCalc!=kProduct) w2 = w1;

  if (corr && mc->IsMeasured()) corr->Fill(w1, w2);

  switch (fCalc) {
  case kProduct: return w1 * w2;
  case kSquare:  return TMath::Sqrt(w1 * w2);
  case kSum:     return 1+(w1-1)+(w2-1);
  case kAverage: return 1+((w1-1)+(w2-1))/2;
  }
  return 1;
}

//____________________________________________________________________
void AliTrackletWeights::ModStack(THStack* stack)
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
void
AliTrackletWeights::Store(TCollection* parent)
{
  TList* top = new TList;
  top->SetName(GetName());
  top->SetOwner(true);
  parent->Add(top);
  top->Add(new TParameter<int>("calc", fCalc, 'f'));
  top->Add(new TParameter<int>("mask", fMask, 'f'));
  top->Add(new TParameter<int>("veto", fVeto, 'f'));
	   
  if (fPt) {
    TH2* copy = static_cast<TH2*>(fPt->Clone("centPt"));
    copy->SetDirectory(0);
    copy->SetBinContent(0,0,1); // For counting merges
    top->Add(copy);
  }
  StoreMap(top, "abundance",   fAbundance);
  StoreMap(top, "strangeness", fStrangeness);
}

//____________________________________________________________________
void
AliTrackletWeights::StoreMap(TCollection* parent, const char* name, PdgMap& m)
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

#ifndef __CINT__
namespace {
  template <typename T>
  T GetPar(const char* name, TList* l)
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

//____________________________________________________________________
Bool_t
AliTrackletWeights::Retrieve(TCollection* parent)
{
  TList* top = static_cast<TList*>(parent->FindObject(GetName()));
  if (!top) {
    Warning("Retrieve",
	    "Collection %s not found in %s", GetName(), parent->GetName());
    parent->ls();
    return false;
  }
  fCalc   = GetPar<int>("calc", top);
  fMask   = GetPar<int>("mask", top);
  fVeto   = GetPar<int>("veto", top);
  
  fPt        = static_cast<TH2D*>(top->FindObject("centPt"));
  if (!fPt)
    Warning("Retrieve","centPt histogram not found in %s", GetName());
  else {
    fPt->SetDirectory(0);
    Double_t scale = fPt->GetBinContent(0,0);
    fPt->Scale(1/scale); // Counting merges
    // fPt->SetBinContent(0,0,1); // Zero merger count 
  }
  if (!RetrieveMap(top, "abundance",   fAbundance)) return false;
  if (!RetrieveMap(top, "strangeness", fStrangeness)) return false;

  return true;
}
//____________________________________________________________________
Bool_t
AliTrackletWeights::RetrieveMap(TCollection* parent,
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
AliTrackletWeights::Draw(Option_t* option)
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
AliTrackletWeights::Print(Option_t* option) const
{
  gROOT->IndentLevel();
  Printf("%s : %s", ClassName(), GetName());
  Printf(" Weight calculation:  %s",
	 fCalc == kProduct ? "product" :
	 fCalc == kSquare  ? "square"  :
	 fCalc == kSum     ? "sum"     : "average");
  Printf(" Tracklet mask:       0x%02x", fMask);
  Printf(" Tracklet veto:       0x%02x", fVeto);
  gROOT->IncreaseDirLevel();
  PrintHist("pT", fPt);
  PrintMap(fAbundance,   "Abundance", option);
  PrintMap(fStrangeness, "Strangeness", option);
  
  gROOT->DecreaseDirLevel();
}
//____________________________________________________________________
void
AliTrackletWeights::PrintHist(const char* tag, TH1* h) const
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
AliTrackletWeights::PrintMap(const PdgMap& m,
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
#endif
//____________________________________________________________________
//
// EOF
// 

