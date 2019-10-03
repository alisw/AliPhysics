/**
 * @file   TupleSelector.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Thu Nov 13 14:12:11 2014
 * 
 * @brief A selector to draw stuff from the a track tuple as made by
 * AliFMDMCTrackELoss.
 * 
 * This class assumes that the files with the TTree's created by 
 * AliFMDMCTrackELoss lives in tuple/forward_tuple_N.root. 
 * 
 * The trees are made by AliFMDMCTrackELoss, which in turn is embedded
 * in a AliFMDMCTrackInspector object.  This code, which also does
 * fits of the MC 'true' energy loss spectra is contained in the task
 * AliFMDMCTrackInspectorTask.  This task can be executed via the
 * train MakeFMDMCTrackTrain.C
 *
 * To run this selector do 
 *
 * @code 
 void
 Run(Bool_t proof=true, Long64_t maxEvents=-1)
 {
   const char* fwd = "${ALICE_PHYSICS}/PWGLF/FORWARD/analysis2";
   gSystem->AddIncludePath("-I${ALICE_PHYSICS}/include");
   gROOT->Macro(Form("%s/scripts/LoadLibs.C"));
   gROOT->LoadMacro(Form("%s/TupleSelector.C++g",fwd));

   if (proof) TupleSelector::Proof(maxEvents);
   else       TupleSelector::Run(maxEvents);
 }
 * @endcode 
 * 
 * Here, $ANA_SRC is assumed to point to the source directory of 
 * PWGLF/FORWARD/analysis2 
 */

#ifndef SELECTOR_C
#define SELECTOR_C

#include <TSelector.h>
#include <TQObject.h>
#ifndef __CINT__
# include <TH1.h>
# include <TString.h>
# include <TCanvas.h>
# include <TLegend.h>
# include <TLegendEntry.h>
# include <TError.h>
# include <TMath.h>
# include <THStack.h>
# include <TTree.h>
# include <TClonesArray.h>
# include <TChain.h>
# include <TSystem.h>
# include <TOutputListSelectorDataMap.h>
# include <TLatex.h>
# include <TROOT.h>
# include <TFile.h>
# include <TDirectory.h>
# include <TSystemDirectory.h>
# include <TRegexp.h>
# include <TKey.h>
# include <TFileCollection.h>
# include <THashList.h>
# include <TDSet.h>
# include <TQConnection.h>
# include <iostream>
# include <TProof.h>
# include <TStopwatch.h>
# include "AliFMDMCTrackELoss.h"
#else
// Forward declarations of types in the interface, and to trigger
// autoloading of some libraries
class TTree;
class TChain;
class TCanvas;
class TH1;
class THStack;
class TLegend;
class TCanvas;
class TString;
class TDirectory;
class TSystemDirectory;
class TRegexp;
class TDSet;
class TProof;
class AliFMDMCTrackELoss;
class AliFMDMCTrackELoss::Hit;
#endif

//====================================================================
/** 
 * Container of primary/secondary histograms 
 */
struct Spectra : public TObject
{
  /** The name */
  TString fName;
  /** Primaries */
  TH1*    fPrimary;
  /** Secondaries */
  TH1*    fSecondary;
  /** 
   * I/O CTOR
   */
  Spectra() 
    : TObject(), fName(""), fPrimary(0), fSecondary(0)
  {}
  /** 
   * User CTOR 
   * 
   * @param name   Name 
   * @param title  Title 
   * @param color  Color 
   * @param bins   Bin definition 
   */
  Spectra(const char*    name, 
	  const char*    title, 
	  Color_t        color,
	  const TArrayD& bins)
    : TObject(), 
      fName(name),
      fPrimary(0), 
      fSecondary(0)
  {
    fPrimary = new TH1D(Form("primary%s", name), 
			title, bins.GetSize()-1, 
			bins.GetArray());
    fPrimary->SetXTitle("#beta#gamma");
    fPrimary->SetMarkerColor(color);
    fPrimary->SetLineColor(color);
    fPrimary->SetMarkerStyle(20);
    fPrimary->SetFillStyle(0);
    fPrimary->SetDirectory(0);
    fPrimary->Sumw2();
      
    fSecondary = static_cast<TH1*>(fPrimary->Clone(Form("secondary%s",name)));
    fSecondary->SetMarkerStyle(24);
    fSecondary->SetMarkerSize(1.2);
    fSecondary->SetDirectory(0);
  }
  /** 
   * Copy CTOR
   * 
   * @param o Object to copy from 
   */ 
  Spectra(const Spectra& o) 
    : TObject(o), 
      fName(o.fName),
      fPrimary(0), 
      fSecondary(0) 
  {
    if (o.fPrimary)
      fPrimary = static_cast<TH1*>(o.fPrimary->Clone());
    if (o.fSecondary)
      fSecondary = static_cast<TH1*>(o.fSecondary->Clone());
  }
  /** 
   * DTOR
   */
  virtual ~Spectra() 
  {
    if (fPrimary)   delete fPrimary;
    if (fSecondary) delete fSecondary;
  }
  /** 
   * Assigment operator 
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this 
   */
  Spectra& operator=(const Spectra& o) 
  {
    if (&o == this) return *this;
      
    TObject::operator=(o);
    fName = o.fName;
    if (fPrimary)   { delete fPrimary;   fPrimary   = 0; }
    if (fSecondary) { delete fSecondary; fSecondary = 0; }
    if (o.fPrimary)
      fPrimary = static_cast<TH1*>(o.fPrimary->Clone());
    if (o.fSecondary)
      fSecondary = static_cast<TH1*>(o.fSecondary->Clone());
      
    return *this;
  }
  /** 
   * Get the name 
   * 
   * @return name 
   */
  const char* GetName() const { return fName.Data(); }
  /** 
   * Fill into histogram 
   * 
   * @param primary true if for primary 
   * @param x       What to fill
   */
  void Fill(Bool_t primary, Double_t x) 
  {
    if (primary) fPrimary->Fill(x);
    else         fSecondary->Fill(x);
  }
  /** 
   * Get a histogram 
   * 
   * @param primary If true, get histogram for primaries, otherwise
   * for secondaries.
   * 
   * @return Pointer to histogram (possibly null)
   */
  TH1* Hist(Bool_t primary) const 
  {
    return (primary ? fPrimary : fSecondary);
  }
  /** 
   * Add the histograms to a stack and legend 
   * 
   * @param stack Stack to add to 
   * @param leg   Legend
   */ 
  virtual void Stack(THStack* stack, TLegend* leg)
  {
    TH1* hists[] = { fPrimary, fSecondary, 0 };
    for (Int_t i = 0; i < 2; i++) { 
      if (!hists[i]) continue;

      Int_t n = hists[i]->GetEntries();
      if (n <= 0) continue;
	
      if (stack) {
	hists[i]->Scale(1. / n, "width");
	stack->Add(hists[i]);

	if (i != 0) continue;
	if (!leg) continue;
	  
	leg->AddEntry(hists[i], hists[i]->GetTitle(), "p");
      }
      else if (leg) { 
	TLegendEntry* e = leg->AddEntry("", 
					(i == 0 ? "Primary" : "Secondary"), 
					"p");
	e->SetMarkerStyle(hists[i]->GetMarkerStyle());
	e->SetFillStyle(hists[i]->GetFillStyle());
	e->SetFillColor(kBlack);
      }
    }
  }
  /** 
   * Merge this object with similar objects in the list.  
   * 
   * @param list list of objects 
   * 
   * @return Number of merged objects 
   */
  Long64_t Merge(TCollection *list)
  {
    TIter    nxt(list);
    TObject* obj = 0;
    Long64_t cnt = 0;
    while ((obj = nxt())) { 
      if (!obj->IsA()->InheritsFrom(this->IsA())) {
	Warning("Merge", "Will not merge a Spectra with a %s", 
		obj->IsA()->GetName());
	continue;
      }
      Spectra* oth  = static_cast<Spectra*>(obj);
      if (Add(oth)) cnt++;
    }
    // Info("Merge", "%s merged %lld entries", fName.Data(), cnt);
    return cnt;
  }
  Bool_t Add(const Spectra* oth) 
  {
    if (!oth->fName.EqualTo(fName)) {
      Warning("Merge", "Will not merge %s with %s", 
	      oth->fName.Data(), fName.Data());
      return false;
    }
    // Int_t nPrmOld = fPrimary->GetEntries();
    // Int_t nSecOld = fPrimary->GetEntries();

    fPrimary  ->Add(oth->fPrimary);
    fSecondary->Add(oth->fSecondary);

    // Int_t nPrmNow = fPrimary->GetEntries();
    // Int_t nSecNow = fPrimary->GetEntries();

    // Info("Add", "%s: %d/%d merged to %d/%d", 
    //      fName.Data(), nPrmOld, nSecOld, nPrmNow, nSecNow);
    return true;
  }
  /** 
   * List this object 
   * 
   * @param option Not used
   */	
  void ls(Option_t* option="") const 
  {
    TObject::ls(option);
    gROOT->IncreaseDirLevel();
    TH1* hists[] = { fPrimary, 
		     fSecondary, 
		     0 };
    for (Int_t i = 0; i < 2; i++) if (hists[i]) hists[i]->ls();
    gROOT->DecreaseDirLevel();
  }
  /** 
   * Store this on file 
   *
   * @param dir Parent directory 
   */
  void Store(TDirectory* dir)
  {
    TDirectory* out = dir->mkdir(fName);
    out->cd();
    fPrimary->Clone("primary")->Write();
    fSecondary->Clone("secondary")->Write();
    dir->cd();
  }
  ClassDef(Spectra,1);
};
//====================================================================
/** 
 * Container to type histograms
 */
struct Type : public Spectra 
{
  /** 
   * I/O CTOR
   */
  Type() : Spectra() {}
  /** 
   * User CTOR 
   * 
   * @param name   Name 
   * @param title  Title 
   * @param color  Color
   */
  Type(const char*    name, 
       const char*    title, 
       Color_t        color)
    : Spectra() // Do not use the base class ctor 
  {
    fName = name;
    fPrimary = new TH1D(Form("primary%s", name),
			title, 6, .5, 6.5);
    fPrimary->SetXTitle("particle type");
    fPrimary->SetMarkerColor(color);
    fPrimary->SetMarkerStyle(20);
    fPrimary->SetFillColor(color);
    fPrimary->SetFillStyle(3001);
    fPrimary->SetLineColor(color);
    fPrimary->SetDirectory(0);
    fPrimary->Sumw2();
    fPrimary->GetXaxis()->SetBinLabel(1, "e^{#pm}");
    fPrimary->GetXaxis()->SetBinLabel(2, "#pi^{#pm}");
    fPrimary->GetXaxis()->SetBinLabel(3, "K^{#pm}");
    fPrimary->GetXaxis()->SetBinLabel(4, "p/#bar{p}");
    fPrimary->GetXaxis()->SetBinLabel(5, "strange");
    fPrimary->GetXaxis()->SetBinLabel(6, "other");
      
    fSecondary = 
      static_cast<TH1*>(fPrimary->Clone(Form("secondary%s",name)));
    fSecondary->SetMarkerStyle(24);
    fSecondary->SetMarkerSize(1.2);
    fSecondary->SetFillStyle(3002);
    fSecondary->SetDirectory(0);
  }
  ClassDef(Type,1); 
};

//====================================================================
struct Ring : public TObject
{

  /** Detector number */
  UShort_t fD;
  /** Ring identifier */
  Char_t   fR;
  /** Cached name */
  TString  fName;
  /** Collection of spectra */
  TObjArray* fSpectra; 

  /** Enumeration of our spectra */
  enum { 
    kBetaGamma = 0,
    kBetaGammaPi, 
    kBetaGammaElectron, 
    kTypes,
    kNSpectra
  };
    
  /** 
   * I/O constructor 
   */
  Ring() 
    : TObject(),
      fD(0), 
      fR('\0'),
      fName(),
      fSpectra(0)
  {}
  /** 
   * Constructor 
   * 
   * @param d Detector 
   * @param r Ring 
   */
  Ring(UShort_t d, Char_t r) 
    : fD(d), 
      fR(r),
      fName(Form("FMD%d%c", d, r)),
      fSpectra(0)
  {
    Color_t color = AliForwardUtil::RingColor(fD, fR);
    TArrayD bgArray(601);
    AliForwardUtil::MakeLogScale(300, -2, 5, bgArray);
    
    fSpectra   = new TObjArray(kNSpectra);
    fSpectra->SetOwner(true);
    fSpectra->AddAt(new Spectra("BetaGamma", fName, color, bgArray), 
		    kBetaGamma);
    fSpectra->AddAt(new Spectra("BetaGammaPi", fName, color, bgArray), 
		    kBetaGammaPi);
    fSpectra->AddAt(new Spectra("BetaGammaElectron", fName, color, bgArray), 
		    kBetaGammaElectron);
    fSpectra->AddAt(new Type("Type", fName, color), kTypes);
    // fSpectra->ls();
  }
  /** 
   * Copy constructor 
   * 
   * @param r Other to copy from 
   */
  Ring(const Ring& r)
    : TObject(r), 
      fD(r.fD), 
      fR(r.fR), 
      fName(r.fName),
      fSpectra(0)
  {
    if (r.fSpectra) fSpectra = new TObjArray(*r.fSpectra);
  }    
  /** 
   * Assignment operator 
   * 
   * @param r Other to assign from 
   * 
   * @return Reference to this obect 
   */
  Ring& operator=(const Ring& r) 
  {
    if (&r == this) return *this;
    
    fD                     = r.fD;
    fR                     = r.fR;
    fName                  = r.fName;
    if (fSpectra) { delete fSpectra; fSpectra = 0; }
    if (r.fSpectra) fSpectra = new TObjArray(*r.fSpectra);

    return *this;
  }
  const char* GetName() const { return fName.Data(); }
  /** 
   * Get a spectra 
   * 
   * @param which Identifier 
   * 
   * @return Pointer to Spectra object or null
   */
  Spectra* Get(UInt_t which) 
  {
    if (!fSpectra) return 0;
    Int_t i = which;
    if (i > fSpectra->GetEntriesFast()) return 0;
    return static_cast<Spectra*>(fSpectra->At(i));
  }
  /** 
   * Get a spectra 
   * 
   * @param which Identifier 
   * 
   * @return Pointer to Spectra object or null
   */
  const Spectra* Get(UInt_t which) const
  {
    if (!fSpectra) return 0;
    Int_t i = which;
    if (i > fSpectra->GetEntriesFast()) return 0;
    return static_cast<Spectra*>(fSpectra->At(i));
  }
  /** 
   * Fill histograms 
   * 
   * @param hit Hit structure 
   */
  void Fill(AliFMDMCTrackELoss::Hit* hit)
  {
    Bool_t prim = hit->IsPrimary();
    Int_t  apdg = hit->AbsPdg();
    Int_t  mfl  = Int_t(apdg/TMath::Power(10,Int_t(TMath::Log10(apdg))));
    Int_t  type = (hit->IsElectron() ? 1 : 
		   hit->IsPion()     ? 2 : 
		   hit->IsKaon()     ? 3 : 
		   hit->IsProton()   ? 4 : 
		   mfl == 3          ? 5 : 6);
    
    Get(kBetaGamma)        ->Fill(prim, hit->BetaGamma());
    Get(kTypes)            ->Fill(prim, type);
    if (hit->IsPion() || hit->IsProton() || hit->IsKaon() || apdg > 100)
      Get(kBetaGammaPi)      ->Fill(prim, hit->BetaGamma());
    if (hit->IsElectron())
      Get(kBetaGammaElectron)->Fill(prim, hit->BetaGamma());

  }
  /** 
   * Get histograms
   * 
   * @param primary If true, for primaries  
   * @param which   Which histogram to get 
   * 
   * @return Histogram or null
   */
  TH1* Hist(Bool_t primary, UShort_t which) const
  {
    const Spectra* spe = Get(which); 
    if (!spe) return 0;
    return spe->Hist(primary);
  }
  /** 
   * Merge this object with similar objects in the list.  
   * 
   * @param list list of objects 
   * 
   * @return Number of merged objects 
   */
  Long64_t Merge(TCollection *list)
  {
    TIter    nxt(list);
    TObject* obj = 0;
    Long64_t cnt = 0;
    while ((obj = nxt())) { 
      if (!obj->IsA()->InheritsFrom(this->IsA())) {
	Warning("Merge", "Will not merge a Ring with a %s", 
		obj->IsA()->GetName());
	continue;
      }
      Ring* oth = static_cast<Ring*>(obj);
      if (oth->fD != fD || oth->fR != fR) {
	Warning("Merge", "Will not merge FMD%d%c with FMD%d%c", 
		fD, fR, oth->fD, oth->fR);
	continue;
      }
      
      if (fSpectra) {
	// fSpectra->Merge(oth->fSpectra);
	TIter thsNext(fSpectra);
	TIter othNext(oth->fSpectra);
	Spectra* thsSpec = 0;
	Spectra* othSpec = 0;
	
	while ((thsSpec = static_cast<Spectra*>(thsNext())) && 
	       (othSpec = static_cast<Spectra*>(othNext()))) 
	  thsSpec->Add(othSpec);

      }
      cnt++;
    }
    // Info("Merge", "FMD%d%c merged %lld entries", fD, fR, cnt);
    return cnt;
  }
  void Draw(Option_t* opt="") 
  { 
    Info("Draw", "Should draw rings here");
    ls(opt); 
  }
  /** 
   * List this object 
   * 
   * @param option Not used
   */
  void ls(Option_t* option="") const 
  {
    TObject::ls(option);
    gROOT->IncreaseDirLevel();
    if (fSpectra) fSpectra->ls(option);
    gROOT->DecreaseDirLevel();
  }
  /** 
   * Store this on file 
   *
   * @param dir Parent directory 
   */
  void Store(TDirectory* dir)
  {
    TDirectory* out = dir->mkdir(fName);
    out->cd();
    TIter next(fSpectra);
    Spectra* spec = 0;
    while ((spec = static_cast<Spectra*>(next()))) 
      spec->Store(out);
    dir->cd();
  }
  
  ClassDef(Ring,2);
};

//====================================================================
struct Monitor : public TObject, public TQObject 
{
  Monitor()
  {
    if (!gProof) return;

    fName = gProof->GetSessionTag();
    gDirectory->Add(this);
    // _must_ specify signal _exactly_ like below, or we won't be called
    Bool_t ret = gProof->Connect("Feedback(TList *objs)", "Monitor", this, 
				 "Feedback(TList *objs)");
    if (!ret) 
      Warning("Monitor", "Failed to connect to Proof");
  }
  virtual ~Monitor() 
  {
    if (!gProof) return;
    gProof->Disconnect("Feedback(TList *objs)",this, 
		       "Feedback(TList* objs)");
  }
  void SetName(const char* name) { fName = name; }
  const char* GetName() const { return fName.Data(); }
  void Feedback(TList* objs)
  {
    Info("Feedback", "Got a list of objects (%p)", objs);
    if (!objs) return;
    objs->ls();
  }
  TString fName;
  ClassDef(Monitor,1);
};


//====================================================================
struct TupleSelector : public TSelector 
{
  TString       fTitle; //! Must not be persistent 
  TTree*        fTree;  //! Must not be persistent 
  TClonesArray* fHits;  //! Must not be persistent 
  Int_t         fI;     //! Not persistent 
  TObjArray*    fRings; //! Not persistent

  /** 
   * Constructor 
   * 
   * @param name Optional title 
   */
  TupleSelector(const char* name="") :  
    fTitle(name), fTree(0), fHits(0), fI(0), fRings(0) 
  {
    if (fTitle.IsNull()) 
      fTitle = gSystem->BaseName(gSystem->WorkingDirectory());
    // SetOption("fb=rings");
    
  }
  const char* GetName() const { return fTitle.Data(); }
  /** 
   * @{ 
   * @name Setup 
   */
  /** 
   * Get index for a particular ring 
   * 
   * @param d Detector 
   * @param r Ring 
   * 
   * @return Index, or 0xFFFF
   */
  UShort_t Index(UShort_t d, Char_t r) const
  {
    Bool_t inner = (r == 'i' || r == 'I');
    switch (d) { 
    case 1: return 0;
    case 2: return (inner ? 1 : 2);
    case 3: return (inner ? 4 : 3);
    }
    ::Warning("", "Unknown ring FMD%d%c", d, r);
    return 0xFFFF;
  }
  void Init(TTree* tree) 
  {
    Info("Init", "Got a tree: %p", tree);

    if (!tree) { 
      Warning("Init", "No tree passed");
      return;
    }
    
    fTree = tree;
    fHits = new TClonesArray("AliFMDMCTrackELoss::Hit");
    fTree->SetBranchAddress("hits", &fHits);
  }    
  void Begin(TTree* tree)
  {
    Info("Begin", "Called w/tree @ %p", tree);
    // if (tree) SlaveBegin(tree);
    new Monitor;
  }
  /** 
   * Begin on slave 
   * 
   * @param tree Tree to process 
   */  
  void SlaveBegin(TTree* tree)
  {
    Info("SlaveBegin", "Got a tree: %p", tree);
    fRings = new TObjArray(5);
    fRings->SetName("rings");
    fRings->SetOwner(true);
    fRings->AddAt(new Ring(1,'I'), Index(1,'I'));
    fRings->AddAt(new Ring(2,'I'), Index(2,'I'));
    fRings->AddAt(new Ring(2,'O'), Index(2,'O'));
    fRings->AddAt(new Ring(3,'O'), Index(3,'O'));
    fRings->AddAt(new Ring(3,'I'), Index(3,'I'));

    fOutput->Add(fRings);
    if (tree) {
      // In case of local mode, we get the tree here, 
      // and init isn't called 
      Init(tree);
    }
  }
  /* @} */

  /** 
   * @{ 
   * @name Event processing 
   */
  Bool_t Notify() 
  {
    TFile* file = (fTree ? fTree->GetCurrentFile() : 0);
    Info("Notify","processing file: %p (%s)", 
	 file, (file ? file->GetName() : "nil"));
    return true;
  }
  /** 
   * Process an entry 
   * 
   * @param entry Entry 
   * 
   * @return true on success 
   */
  Bool_t Process(Long64_t entry)
  {
    if (!fTree) {
      Warning("Process", "No tree");
      return false;
    }
    if (!fTree->GetTree()) {
      Warning("Process", "No real tree");
      return false;
    }
    fHits->Clear();
    fTree->GetTree()->GetEntry(entry);
    fI++;
    // if (fI % 100 == 0) {
    //   printf("Event # %6d (%6d hits)\r", fI, fHits->GetEntries());
    //   fflush(stdout);
    // }
    // Printf("Event # %7d, %6d hits", fI, fHits->GetEntries());
  
    TIter next(fHits);
    AliFMDMCTrackELoss::Hit* hit = 0;
    while ((hit = static_cast<AliFMDMCTrackELoss::Hit*>(next()))) 
      Fill(hit);
   return true;
  }
  /** 
   * Get a Ring 
   * 
   * @param d Detector 
   * @param r Ring 
   *
   * @return Pointer to ring object or null
   */
  Ring* Get(UShort_t d, Char_t r) 
  {
    return static_cast<Ring*>(fRings->At(Index(d,r)));
  }
  /** 
   * Get a Ring 
   * 
   * @param d Detector 
   * @param r Ring 
   *
   * @return Pointer to constant ring object or null
   */
  const Ring* Get(UShort_t d, Char_t r) const
  {
    return static_cast<Ring*>(fRings->At(Index(d,r)));
  }
  /** 
   * Fill in to a ring object 
   * 
   * @param hit Hit information
   */
  void Fill(AliFMDMCTrackELoss::Hit* hit)
  {
    Ring* r = Get(hit->D(), hit->R());
    r->Fill(hit);
  }
  /* @} */

  /** 
   * @{ 
   * @name Final processing 
   */
  void SlaveTerminate()
  {
  }
  /** 
   * Terminate the job 
   * 
   */
  void Terminate()
  {
    Printf("\nDone (rings=%p)", fRings);

    if (!fRings) 
      fRings = static_cast<TObjArray*>(fOutput->FindObject("rings"));
    if (!fRings) { 
      Error("Terminate", "No rings in output array");
      return;
    }
    // fRings->ls();
    
    TFile* out = TFile::Open("tuple_summary.root", "RECREATE");
    Store(out);

    PlotOne(Ring::kBetaGamma,        "#beta#gamma - all",     out);
    PlotOne(Ring::kBetaGammaPi,      "#beta#gamma - Hadrons", out);
    PlotOne(Ring::kBetaGammaElectron,"#beta#gamma - e^{#pm}", out);
    PlotOne(Ring::kTypes,            "Relative abundance",    out);
  }
  /** 
   * Plot one quantity 
   * 
   * @param which 
   * @param title Title 
   * @param dir   Directory 
   * @param opt   Options
   */
  void PlotOne(UShort_t which, const char* title, 
	       TDirectory* dir=0, Option_t* opt="")
  {
    static Int_t cnt  = 1;
    const char*  name = Form("c%02d", cnt++);
      
    TCanvas* c = new TCanvas(name, title, 1000, 1000);
    c->SetFillColor(0);
    c->SetFillStyle(0);
    c->SetTopMargin(0.01);
    c->SetRightMargin(0.03);

    TLegend* leg   = new TLegend(0.65, 0.65, 0.975, 0.975);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);

    TString xTitle = "";
    THStack* stack = new THStack(name, title);
    for (Int_t i = 0; i < 5; i++) { 
      Ring*    ring = static_cast<Ring*>(fRings->At(i));
      Spectra* spec = ring->Get(which);
      if (!spec) { 
	Warning("PlotOne", "No spectra for %d", which);
	continue;
      }

      // Add to our stack 
      spec->Stack(stack, leg);

      // Get X title 
      if (xTitle.IsNull() && spec->Hist(true))
	xTitle = spec->Hist(true)->GetXaxis()->GetTitle();
      
      // If we're not at the last one continue 
      if (i < 4) continue;

      // Add primary/secondary entry to legend 
      spec->Stack(0, leg);

    }
    c->cd();
    if (!stack->GetHists() || 
	stack->GetHists()->GetEntries() < 0) {
      Warning("PlotOne", "No histograms added");
      return;
    }
    stack->Draw(Form("nostack %s", opt));
    stack->GetHistogram()->SetXTitle(xTitle);
    stack->GetHistogram()->GetListOfFunctions()->Add(leg);
    
    leg->Draw();

    if (which != Ring::kTypes) {
      c->SetLogx();
      c->SetLogy();
    }

    TLatex* ltx = new TLatex(0.02, 0.02, fTitle);
    ltx->SetNDC();
    ltx->SetTextColor(kRed+2);
    ltx->SetTextAlign(11);
    ltx->SetTextSize(0.02);
    ltx->Draw();
    stack->GetHistogram()->GetListOfFunctions()->Add(ltx);

    if (dir) {
      dir->cd();
      stack->Clone()->Write();
    }
      

    c->Modified();
    c->Update();
    c->cd();

    c->Print(Form("%s.pdf", name));
  }
  /** 
   * Store this on file 
   *
   * @param dir Parent directory 
   */
  void Store(TDirectory* dir)
  {
    TDirectory* out = dir->mkdir(fTitle);
    out->cd();
    TIter next(fRings);
    Ring* ring = 0;
    while ((ring = static_cast<Ring*>(next()))) 
      ring->Store(out);
    dir->cd();
  }
  /* @} */

  /** 
   * @{
   * @name Other interface 
   */
  /** 
   * Get the version of the selector 
   * 
   * @return 1 
   */
  Int_t Version() const { return 1; }
  /** 
   * Get the status
   * 
   * @return Number of processedd events 
   */
  Long64_t GetStatus() const { return fI; }
  
  /* @} */
    
  //------------------------------------------------------------------
  /** 
   * @{ 
   * @name Service functions 
   */
  //------------------------------------------------------------------
  /** 
   * Get a chain or data set from the given file 
   * 
   * @param file File to look in 
   * 
   * @return Pointer to object or null
   */
  static TObject* GetChainOrDataSet(const char* file)
  {
    if (gSystem->AccessPathName(file)) return 0;

    TFile* in = TFile::Open(file, "READ");
    if (!in) return 0;

    TObject* ret = in->Get("tree");
    if (!ret) { 
      in->Close();
      return 0;
    }
    
    if (ret->IsA()->InheritsFrom(TChain::Class())) 
      static_cast<TChain*>(ret)->SetDirectory(0);
    else if (ret->IsA()->InheritsFrom(TDSet::Class())) 
      static_cast<TDSet*>(ret)->SetDirectory(0);
    else { 
      ::Warning("GetChainOrDataSet", "Found object is a %s", 
		ret->IsA()->GetName());
      ret = 0;
    }
    in->Close();
    return ret;
  }
  //------------------------------------------------------------------
  /** 
   * make our data set
   * 
   * @param src        Source directory 
   * @param recursive  Whether to scan recursively 
   * @param verbose    Be verbose
   * 
   * @return Data set or null
   */
  static TDSet* MakeDataSet(const TString& src=".", 
			    Bool_t recursive=false, 
			    Bool_t verbose=false)
  {
    TString dsFile(Form("%s/dataset.root", src.Data()));
    TDSet* dataset = static_cast<TDSet*>(GetChainOrDataSet(dsFile));
    if (dataset) {
      /// dataset->Print("a");
      return dataset;
    }

    TChain* c = DoMakeChain(src, recursive, verbose);
    if (!c) return 0;
    
    dataset = new TDSet(*c, false);
    dataset->SetName("tree");
    dataset->SetLookedUp();
    dataset->Validate();
    
    delete c;
    if (dataset) { 
      TFile* out = TFile::Open(dsFile, "RECREATE");
      dataset->Write();
      dataset->SetDirectory(0);
      out->Close();
    }
      
    return dataset;
    
  }
  //------------------------------------------------------------------
  /** 
   * Create our chain 
   * 
   * @param src        Source directory 
   * @param recursive  Whether to scan recursively 
   * @param verbose    Be verbose
   * 
   * @return Chain or null
   */
  static TChain* MakeChain(const TString& src=".", 
			   Bool_t recursive=false,
			   Bool_t verbose=false)
  {
    TString chFile(Form("%s/chain.root", src.Data()));
    if (!gSystem->AccessPathName(chFile)) { 
      TFile* in = TFile::Open(chFile, "READ");
      if (in) { 
	TChain* ret = static_cast<TChain*>(in->Get("tree"));
	if (ret) {
	  ret->SetDirectory(0);
	  // ret->GetListOfFiles()->ls();
	}
	in->Close();
	if (ret) return ret;
      }
    }
	  
    TChain* chain = DoMakeChain(src, recursive, verbose);
    if (chain) { 
      TFile* out = TFile::Open(chFile, "RECREATE");
      chain->Write();
      chain->SetDirectory(0);
      out->Close();
    }
      
    return chain;
  }
  //------------------------------------------------------------------
  /** 
   * Create our chain 
   * 
   * @param src        Source directory 
   * @param recursive  Whether to scan recursively 
   * @param verbose    Be verbose
   * 
   * @return Chain or null
   */
  static TChain* DoMakeChain(const TString& src=".", 
			     Bool_t recursive=false,
			     Bool_t verbose=false)
  {
    TChain* chain = new TChain("tree"); 
    TString savdir(gSystem->WorkingDirectory());
    TSystemDirectory d(gSystem->BaseName(src.Data()), src.Data());
    if (!ScanDirectory(&d, chain, "forward_tree_*", recursive, verbose)) {
      delete chain;
      chain = 0;
    }
    else if (!chain->GetListOfFiles() || 
	chain->GetListOfFiles()->GetEntries() <= 0) {
      delete chain;
      chain = 0;
    }
    return chain;
  }
  
  //------------------------------------------------------------------
  /** 
   * Scan directory @a dir (possibly recursive) for tree files to add
   * to the chain.    This does not follow sym-links
   * 
   * @param dir        Directory to scan
   * @param chain      Chain to add to
   * @param pattern    File name pattern 
   * @param recursive  Be reursive 
   * @param verbose    Verbosity flags
   *
   * @return true if any files where added 
   */
  static Bool_t ScanDirectory(TSystemDirectory* dir,
			      TChain*           chain, 
			      const TString&    pattern, 
			      Bool_t            recursive,
			      Bool_t            verbose) 
  {
    if (!dir) {
      ::Error("ScanDirectory", "No diretory passed");
      return false;
    }
    if (verbose) ::Info("ScanDirectory", "Scanning %s", dir->GetName());

    Bool_t  ret = false;
    TRegexp wild(pattern, true);
    TString oldDir(gSystem->WorkingDirectory());
    TList* files = dir->GetListOfFiles();

    if (!gSystem->ChangeDirectory(oldDir)) { 
      ::Error("ScanDirectory", "Failed to go back to %s", 
	      oldDir.Data());
      return false;
    }
    if (!files) {
      ::Warning("ScanDirectory", "No files");
      return false;
    }

    TList toAdd;
    toAdd.SetOwner();

    // Sort list of files and check if we should add it 
    files->Sort();
    TIter next(files);
    TSystemFile* file = 0;
    while ((file = static_cast<TSystemFile*>(next()))) {
      TString name(file->GetName());
      TString title(file->GetTitle());
      TString full(gSystem->ConcatFileName(file->GetTitle(), name.Data()));
      if (file->IsA()->InheritsFrom(TSystemDirectory::Class())) full = title;


      if (name.EqualTo(".")||name.EqualTo(".."))  {
	// Ignore special links 
	if (verbose) ::Info("ScanDirectory", "Ignoring special %s",
			    name.Data());
	continue;
      }

      if (verbose) ::Info("ScanDirectory", "Looking at %s", full.Data());

      FileStat_t fs;
      if (gSystem->GetPathInfo(full.Data(), fs)) {
	::Warning("ScanDirectory", "Cannot stat %s (%s)", 
		  full.Data(), gSystem->WorkingDirectory());
	continue;
      }
      // Check if this is a directory 
      if (file->IsDirectory(full)) { 
	if (verbose) ::Info("ScanDirectory", "Got a directory: %s", 
			    full.Data());
	if (recursive) {
	  // if (title[0] == '/') 
	  TSystemDirectory* d = new TSystemDirectory(file->GetName(),
						     full.Data());
	  if (ScanDirectory(d, chain, pattern, recursive, verbose))
	    ret = true;
	  delete d;
	}
	continue;
      }

      // If this is not a root file, ignore 
      if (!name.EndsWith(".root") && 
	  !name.EndsWith(".zip",TString::kIgnoreCase)) {
	if (verbose) ::Info("ScanDirectory", "Ignoring non-ROOT or ZIP %s",
			    name.Data());
	continue;
      }

      // If this file does not contain AliESDs, ignore 
      if (!name.Contains(wild)) {
	if (verbose) ::Info("ScanDirectory", "%s does not match %s",
			    name.Data(), pattern.Data());
	continue;
      }
      // Add 
      if (verbose) 
	::Info("::ScanDirectory", "Candidate %s", full.Data());
      toAdd.Add(new TObjString(full));
    }
    TIter nextAdd(&toAdd);
    TObjString* s = 0;
    Int_t added = 0;
    while ((s = static_cast<TObjString*>(nextAdd()))) {
      // Info("ChainBuilder::ScanDirectory", 
      //      "Adding %s", s->GetString().Data());
      TString fn = s->GetString();
      if (!CheckFile(fn, chain, true/*verbose*/)) continue;
       
      added++;
    }
    if (added > 0) ret = true;
    if (verbose) ::Info("ScanDirectory", "Added %d files", added);

    gSystem->ChangeDirectory(oldDir);
    return ret;  
  }
  //------------------------------------------------------------------
  /** 
   * Check if we can add a file to the chain 
   * 
   * @param path   Full path to file 
   * @param chain  Chain 
   * @param verbose Whether to be verbose 
   * 
   * @return true on success, false otherwise
   */
  static Bool_t CheckFile(const TString& path, 
			  TChain* chain,
			  Bool_t verbose)
  {
    // Info("", "Checking %s", path.Data());
    TString fn   = path;

    gSystem->RedirectOutput("/dev/null", "w");
    TFile*  test = TFile::Open(fn, "READ");
    gSystem->RedirectOutput(0);
    if (!test) { 
      ::Warning("CheckFile", "Failed to open %s", fn.Data());
      return false;
    }
    
    Bool_t           ok = false;
    TObject*         o  = test->Get(chain->GetName());
    TTree*           t  = dynamic_cast<TTree*>(o);
    TFileCollection* c  = dynamic_cast<TFileCollection*>(o);
    if (t) {
      Int_t n = t->GetEntries();
      test->Close();
      chain->Add(fn, n);
      ok = true;
      if (verbose) ::Info("CheckFile", "Added file %s (%d)", fn.Data(), n);
    } else if (c) {
      chain->AddFileInfoList(c->GetList());
      ok = true;
      if (verbose) ::Info("CheckFile", "Added file collection %s", fn.Data());
    } else {
      // Let's try to find a TFileCollection 
      TList* l = test->GetListOfKeys();
      TIter next(l);
      TKey* k = 0;
      while ((k = static_cast<TKey*>(next()))) {
	TString cl(k->GetClassName());
	if (!cl.EqualTo("TFileCollection")) continue;
	c = dynamic_cast<TFileCollection*>(k->ReadObj());
	if (!c) { 
	  ::Warning("CheckFile", "Returned collection invalid");
	  continue;
	}
	// Info("", "Adding file collection");
	chain->AddFileInfoList(c->GetList());
	ok = true;
	if (verbose) ::Info("CheckFile", "Added file collection %s", fn.Data());
      }
      test->Close();
    }

    if (!ok) 
      ::Warning("CheckFile", 
		"The file %s does not contain the tree %s "
		"or a file collection", 
		path.Data(), chain->GetName());
    return ok;
  }

  /** 
   * Run this selector on a chain locally 
   * 
   * @param maxEvents Maximum number of events 
   * @param title     Optional title 
   * 
   * @return true on sucess 
   */
  static Bool_t Run(Long64_t maxEvents=-1, 
		    const char* title="")
  {
    TStopwatch timer;
    timer.Start();
    TChain* chain = MakeChain("tuple");
    if (!chain) { 
      ::Error("Run", "No chain!");
      return false;
    }

    TupleSelector* s = new TupleSelector(title);
    Int_t status= chain->Process(s, "", maxEvents);
    timer.Print();
    return status >= 0;
  }
  /** 
   * Run this selector on a chain in Proof 
   * 
   * @param maxEvents Maximum number of events 
   * @param title     Optional title 
   * @param opt       Options
   * 
   * @return true on sucess 
   */
  static Bool_t Proof(Long64_t    maxEvents, 
		      const char* opt="",
		      const char* title="")
  {
    TStopwatch timer;
    timer.Start();
    TProof::Reset("lite:///?workers=8");
    TProof::Open("lite:///?workers=8");
    gProof->ClearCache();
    TString ali = gSystem->ExpandPathName("$(ALICE_PHYSICS)");
    TString fwd = ali + "/PWGLF/FORWARD/analysis2";
    gProof->AddIncludePath(Form("%s/include", ali.Data()));
    gProof->Load(Form("%s/scripts/LoadLibs.C",fwd.Data()), true);
    gProof->Exec("LoadLibs()");
    // gROOT->ProcessLine("gProof->SetLogLevel(5);");
    gProof->Load(Form("%s/scripts/TupleSelector.C+%s", fwd.Data(), opt),true);

    TDSet* dataset = MakeDataSet("tuple");
    if (!dataset) { 
      ::Error("Proof", "No dataset");
      return false;
    }
    
    TupleSelector* s = new TupleSelector(title);
    gProof->AddFeedback("rings");
    gProof->Process(dataset, s, "", maxEvents);
    timer.Print();
    return true; // status >= 0;
  }    

  /* @} */
  ClassDef(TupleSelector,2);

private:
  /** 
   * Copy constructor 
   */
  TupleSelector(const TupleSelector&); // {}
  /** 
   * Assignment operator
   * 
   * @return Reference to this object 
   */
  TupleSelector& operator=(const TupleSelector&); // { return *this; } 

};


  
#endif
// Local Variables:
//   mode: C++
// End:
// 
// EOF
// 
