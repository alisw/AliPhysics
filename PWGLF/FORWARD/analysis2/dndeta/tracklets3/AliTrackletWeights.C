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
#else
class TH1D;
class TH2D;
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
  
  /** Map a particle species to a weight */
  typedef std::map<short,TH1D*> PdgMap;
  /**
   * Default constructor - ROOT I/O only
   * 
   */
  AliTrackletWeights() : TNamed(), fPt(0), fAbundance(), fStrangeness() {}
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
   * 
   * @return The weight
   */
  virtual Double_t LookupWeight(AliAODTracklet* tracklet, Double_t cent) const;
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
  void Draw(Option_t* option="");
  /** 
   * Print information on reweights 
   *
   * @param option Not used 
   */
  void Print(Option_t* option="") const;
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

  Double_t GetAbundanceWeight(UShort_t apdg, Double_t cent) const
  {
    return GetPdgWeight(fAbundance, apdg, cent);
  }
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
  
  TH2D*  fPt;     // Weight by centrality and pT
  PdgMap fAbundance;  // Map for abundance weight 
  PdgMap fStrangeness;    // Map for strangeness weight 
  
  
  ClassDef(AliTrackletWeights,1); // Weighs for tracklet analysis 
};

//____________________________________________________________________
AliTrackletWeights::AliTrackletWeights(const char* name,
				       const char* title)
  : TNamed(name, title),
    fPt(0),
    fAbundance(),
    fStrangeness()
{}

//____________________________________________________________________
AliTrackletWeights::AliTrackletWeights(const AliTrackletWeights& o)
  : TNamed(o),
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
AliTrackletWeights& AliTrackletWeights::operator=(const AliTrackletWeights& o)
{
  if (&o == this) return *this;
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
  if (!fPt) return;
  fPt->SetBit(mode);
}
//____________________________________________________________________
void AliTrackletWeights::SetPdgMode(PdgMap&  m,
				    Short_t  pdg,
				    UShort_t mode)
{
  TH1* h = GetPdgHist(m, pdg);
  if (!h) return;
  h->SetBit(mode);
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
AliTrackletWeights::LookupWeight(AliAODTracklet* tracklet, Double_t cent) const
{
  Double_t          w  = 1;
#if 0
  if (!tracklet->IsSimulated()) {
    Warning("LookupWeight", "Not a simulated tracklet");
    return w;
  }
#endif
  
  AliAODMCTracklet* mc = static_cast<AliAODMCTracklet*>(tracklet);
  if (mc->GetParentPdg()>0) w *= LookupWeight(mc->GetParentPt(),
					      mc->GetParentPdg(),
					      cent);
  // if (!mc->IsCombinatorics()) return w;
  if (mc->GetParentPdg(true)>0) w *= LookupWeight(mc->GetParentPt(true),
						  mc->GetParentPdg(true),
						  cent);
  // Printf("Got tracklet weight %f", w);
  return w;
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
  fPt        = static_cast<TH2D*>(top->FindObject("centPt"));
  if (!fPt) {
    Warning("Retrieve","centPt histogram not found in %s", GetName());
    return false;
  }
  fPt->SetDirectory(0);
  Double_t scale = fPt->GetBinContent(0,0);
  fPt->Scale(1/scale); // Counting merges
  
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
    parent->ls();
    return false;
  }
  TIter    next(top);
  TObject* o = 0;
  while ((o = next())) {
    if (!o->IsA()->InheritsFrom(TH1D::Class())) continue;
    TH1D* copy = static_cast<TH1D*>(o);
    copy->SetDirectory(0);
    Double_t scale = copy->GetBinContent(0); 
    copy->Scale(1/scale); // Counting merges
    TString nme(copy->GetName());
    nme.Remove(0,1);
    Int_t apdg = nme.Atoi();
    AddPdgWeight(m, apdg, copy, copy->TestBits(kUp|kDown|kDisabled));
  }
  return true;
}

//____________________________________________________________________
void
AliTrackletWeights::Draw(Option_t*)
{
  TVirtualPad* master = gPad;
  if (!master) return;
  master->Divide(3,1);
  if (fPt) {
    THStack* stack = new THStack(fPt, "y");
    ModStack(stack);
    stack->SetTitle("#it{p}_{T} weights");
    master->cd(1);
    stack->Draw("nostack");
    stack->GetHistogram()->SetXTitle("#it{p}_{T}");
    master->GetPad(1)->BuildLegend();
  }

  THStack* ha = new THStack("abundance", "Abundance weights");
  for (PdgMap::const_iterator i = fAbundance.begin();
       i != fAbundance.end(); ++i) {
    ha->Add(i->second);
  }
  ModStack(ha);
  master->cd(2);
  ha->Draw("nostack");
  ha->GetHistogram()->SetXTitle("Centrality [%]");
  master->GetPad(2)->BuildLegend();

  THStack* hs = new THStack("strangeness", "Strangeness weights");
  for (PdgMap::const_iterator i = fStrangeness.begin();
       i != fStrangeness.end(); ++i) {
    hs->Add(i->second);
  }
  ModStack(hs);
  master->cd(3);
  hs->Draw("nostack");
  hs->GetHistogram()->SetXTitle("Centrality [%]");
  master->GetPad(3)->BuildLegend();
}

//____________________________________________________________________
void
AliTrackletWeights::Print(Option_t* option) const
{
  gROOT->IndentLevel();
  Printf("%s : %s", ClassName(), GetName());
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

