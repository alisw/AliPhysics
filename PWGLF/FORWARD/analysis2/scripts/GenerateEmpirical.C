#ifndef __CINT__
#include <TFile.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TCollection.h>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TClass.h>
#include <TMultiGraph.h>
#include <THStack.h>
#include <TError.h>
#else
class TCollection;
class TDirectory;
class TFile;
class TH1;
class TH2;
class TGraph;
class THStack;
class TAxis;
class TMultiGraph; 
#endif

struct EmpiricalMaker {
  /** 
   * Open a file 
   * 
   * @param fname Filename of file to open 
   * @param rw    If true, update for read/write acess
   * 
   * @return Pointer to file handle or null
   */
  static TFile* OpenFile(const TString& fname, Bool_t rw=false)
  {
    TFile* file = TFile::Open(fname, rw ? "RECREATE" : "READ");  
    if (!file) { 
      ::Error("OpenFile", "Couldn't open %s for %s", 
	      fname.Data(), (rw ? "read/write" : "read"));
      return 0;
    }
    return file;
  }
  /** 
   * Check the type of an object
   * 
   * @param o     Object to check
   * @param cl    Class of object required
   * @param quiet If true, report no errors
   * 
   * @return true on success, false otherwise 
   */
  static Bool_t CheckType(TObject* o, TClass* cl, Bool_t quiet=false)
  { 
    if (!o || !cl) { 
      if (!quiet) ::Error("CheckType", "No object and/or class given");
      return false;
    }
    if (!o->IsA()->InheritsFrom(cl)) {				
      if (!quiet) ::Error("CheckType", "Object %s is not a %s, but a %s", 
	    o->GetName(), cl->GetName(), o->ClassName());
      return false;
    }
    return true;
  }
  /** 
   * Get an object from a directory 
   * 
   * @param d     Parent directory
   * @param name  Name of object
   * @param cl    Class pointer for testing 
   * @param quiet If true, report no errors
   * 
   * @return Pointer to object, or null
   */
  static TObject* GetObject(TDirectory* d, const TString& name, 
			    TClass* cl=0, Bool_t quiet=false) 
  {
    if (!d) {					
      if (!quiet) ::Error("GetObject", "No parent directory");
      return 0;
    }
    TObject* o = d->Get(name.Data());
    if (!o) {					
      if (!quiet) ::Error("GetObject", "didn't find the object %s in %s",
			  name.Data(), d->GetName());
      if (!quiet) d->ls();
      return 0;
    }
    if (cl && !CheckType(o, cl, quiet)) return 0;

    return o;
  }
  /** 
   * Get an object from a directory 
   * 
   * @param d     Parent directory
   * @param name  Name of object
   * @param quiet If true, report no errors
   * @param cl    Class pointer for testing 
   * 
   * @return Pointer to object, or null
   */
  static TObject* GetObject(TCollection* d, const TString& name, 
			    TClass* cl=0, Bool_t quiet=false) 
  {
    if (!d) {					
      if (!quiet) ::Error("GetObject", "No parent directory");
      return 0;
    }
    TObject* o = d->FindObject(name.Data());
    if (!o) {					
      if (!quiet) ::Error("GetObject", "didn't find the object %s in %s", 
			  name.Data(), d->GetName());
      if (!quiet) d->ls();
      return 0;
    }
    if (cl && !CheckType(o, cl, quiet)) return 0;

    return o;
  }
  /** 
   * Get a collection from a directory 
   * 
   * @param d      Parent 
   * @param name   Name of collection to get 
   * @param quiet  If true, report no errors
   * 
   * @return The found collection or null
   */
  static TCollection* GetCollection(TDirectory* d, const TString& name, 
				    Bool_t quiet=false) 
  {
    return static_cast<TCollection*>(GetObject(d, name, TCollection::Class(), 
					       quiet));
  }
  /** 
   * Get a collection from another collection
   * 
   * @param d      Parent 
   * @param name   Name of collection to get 
   * @param quiet  If true, report no errors
   * 
   * @return The found collection or null
   */
  static TCollection* GetCollection(TCollection* d, const TString& name, 
				    Bool_t quiet=false) 
  {
    return static_cast<TCollection*>(GetObject(d, name, TCollection::Class(), 
					       quiet));
  }
  /** 
   * Get a histogram from a directory 
   * 
   * @param d      Parent 
   * @param name   Name of collection to get 
   * @param quiet  If true, report no errors
   * 
   * @return The found histogram or null
   */
  static TH1* GetH1(TDirectory* d, const TString& name, Bool_t quiet=false) 
  {
    return static_cast<TH1*>(GetObject(d, name, TH1::Class(), quiet));
  }
  /** 
   * Get a histogram from a collection
   * 
   * @param d      Parent 
   * @param name   Name of collection to get 
   * @param quiet  If true, report no errors
   * 
   * @return The found histogram or null
   */
  static TH1* GetH1(TCollection* d, const TString& name, Bool_t quiet=false) 
  {
    return static_cast<TH1*>(GetObject(d, name, TH1::Class(), quiet));
  }
  /** 
   * Get an graph from a collection
   * 
   * @param d      Parent 
   * @param name   Name of collection to get 
   * @param quiet  If true, report no errors
   * 
   * @return The found graph or null
   */
  static TGraph* GetGraph(TCollection* d, const TString& name, Bool_t quiet=false) 
  {
    return static_cast<TGraph*>(GetObject(d, name, TGraph::Class(), quiet));
  }
  /** 
   * Get an axis from a collection
   * 
   * @param d      Parent 
   * @param name   Name of collection to get 
   * @param quiet  If true, report no errors
   * 
   * @return The found axis or null
   */
  static TAxis* GetAxis(TCollection* d, const TString& name, 
			Bool_t quiet=false) 
  {
    return static_cast<TAxis*>(GetObject(d, name, TAxis::Class(), quiet));
  }
  /**
   * Compute the ratio of @a h to @a g.  @a g is evaluated at the bin
   * centers of @a h 
   * 
   * @param h  Numerator 
   * @param g  Divisor 
   * 
   * @return h/g 
   */
  static TH1* RatioHG(const TH1* h, const TGraph* g) 
  {
    if (!h || !g) return 0;

    TH1* ret = static_cast<TH1*>(h->Clone("tmp"));
    ret->Reset();
    ret->SetDirectory(0);
    Double_t xlow  = g->GetX()[0];
    Double_t xhigh = g->GetX()[g->GetN()-1];
    if (xlow > xhigh) { Double_t t = xhigh; xhigh = xlow; xlow = t; }

    for (Int_t i = 1; i <= h->GetNbinsX(); i++) { 
      Double_t c = h->GetBinContent(i);
      if (c <= 0) continue;

      Double_t x = h->GetBinCenter(i);
      if (x < xlow || x > xhigh) continue; 

      Double_t f = g->Eval(x);
      if (f <= 0) continue; 

      ret->SetBinContent(i, c / f);
      ret->SetBinError(i, h->GetBinError(i) / f);
    }
    return ret;
  }
  //==================================================================
  /** 
   * Constructor
   */
  EmpiricalMaker() 
    : fRefLoaded(false)
  {}
  /** 
   * Load reference data as compiled script
   * 
   */
  void LoadReferences()
  {
#if 0
    if (fRefLoaded) return;
    // --- Set the macro pathand load other data script --------------
    TString savPath(gROOT->GetMacroPath());
    gROOT->SetMacroPath(Form("%s:$(ALICE_PHYSICS)/PWGLF/FORWARD/analysis2",
			     gROOT->GetMacroPath()));
    // Always recompile 
    if (!gROOT->GetClass("RefData"))
      gROOT->LoadMacro("OtherData.C+");
    gROOT->SetMacroPath(savPath);
    fRefLoaded = true;
#endif 
  }
  /** 
   * Get reference dN/deta in some centrality range
   * 
   * @param c1  Low value of centrality range 
   * @param c2  High value of centrallity range 
   * 
   * @return Reference data, or null
   */
  TMultiGraph* GetReference(UShort_t c1, UShort_t c2)
  {
    return 0;
#if 0
     
    if (!fRefLoaded) LoadReferences();

    TMultiGraph* other = 0;
    UShort_t     sys   = 2; // (fSysString  ? fSysString->GetUniqueID() : 0);
    UShort_t     trg   = 0; // (fTrigString ? fTrigString->GetUniqueID() : 0);
    UShort_t     snn   = 2760; // (fSNNString  ? fSNNString->GetUniqueID() : 0);
    Long_t       ret   = 
      gROOT->ProcessLine(Form("RefData::GetData(%d,%d,%d,%d,%d,%d);",
			      sys,snn,trg,c1,c2,0xF));
    if (!ret) {
      Warning("GetReference", 
	      "No other data for sys=%d trg=%d sNN=%d c=%3d%%-%3d%%",
	      sys, trg, snn, c1, c2);
      return 0;
    }

    other = reinterpret_cast<TMultiGraph*>(ret);    
    return other;
#endif 
  }
  /** 
   * Run the full thing 
   * 
   * @param fileName File to query results
   */
  void Run(const TString& fileName) 
  {
    TFile* file = OpenFile(fileName, false);
    if (!file) return;
    
    TFile* out = OpenFile("empirical.root", true);
    if (!out) return;

    ProcessComponent(file, out, "Forward");
    ProcessComponent(file, out, "Central");

    out->ls();
    out->Close();
    file->Close();
  }
  /** 
   * Process a single component
   * 
   * @param d     Parent
   * @param out   Output 
   * @param name  Name of component
   */
  void ProcessComponent(TDirectory* d, TDirectory* out, const TString& name)
  {
    TDirectory* store = out->mkdir(name);
    store->cd();

    TString resName(Form("%sdNdetaResults", name.Data()));
    TCollection* results = GetCollection(d, resName);
    if (!results) return;

    TAxis* centAxis = GetAxis(results, "centAxis");
    if (!centAxis) return;

    THStack* corrs = new THStack("empirical", "Empirical correction");
    // TMultiGraph* corrs = new TMultiGraph("empirical");
    corrs->SetTitle(Form("Empirical corrections for %s", name.Data()));
    for (Int_t i=1; i<=centAxis->GetNbins(); i++) { 
      UShort_t c1 = UShort_t(centAxis->GetBinLowEdge(i));
      UShort_t c2 = UShort_t(centAxis->GetBinUpEdge(i));
    
      ProcessCent(c1, c2, results, name, corrs, store);
    }
    
    if (!corrs->GetHists() ||
	corrs->GetHists()->GetEntries() <= 0) return;

    store->cd();

    TH1* sum = static_cast<TH1*>(corrs->GetHists()->At(0)->Clone("mean"));
    sum->SetTitle("mean");
    sum->Reset();
    sum->SetMarkerColor(kBlack);
    sum->SetDirectory(0);
    Int_t nHist = corrs->GetHists()->GetEntries();
    for (Int_t i = 0; i < nHist; i++) 	
      sum->Add(static_cast<TH1*>(corrs->GetHists()->At(i)));
    sum->Scale(1. / nHist);
    sum->Write("default");

    THStack* ratios = new THStack("ratios", "Ratio to mean");
    for (Int_t i = 0; i < nHist; i++) {
      TH1* h = static_cast<TH1*>(corrs->GetHists()->At(i));
      TH1* r = static_cast<TH1*>(h->Clone());
      r->SetDirectory(0);
      r->Divide(sum);
      ratios->Add(r);
    }
    

    TCanvas* c = new TCanvas(name);
    c->SetRightMargin(0.05);
    c->SetTopMargin(0.05);
    c->Divide(1,2,0,0);
    TVirtualPad* p = c->cd(1);
    corrs->Add(sum);
    corrs->DrawClone("nostack");
    p->BuildLegend(0.2, 0.05, 0.4, 0.4);
    corrs->Write();

    p = c->cd(2);
    THStack* cpy = static_cast<THStack*>(ratios->DrawClone("nostack"));
    cpy->SetMinimum(0.95);
    cpy->SetMaximum(1.05);
    p->BuildLegend(0.35,0.15,0.55,0.4);
    ratios->Write();

    c->Print(Form("%sEmpirical.pdf", name.Data()));
  }
  /** 
   * Process a single centrality bin
   * 
   * @param c1     Lower cut
   * @param c2     Higher cut
   * @param c      Parent
   * @param name   Name of component 
   * @param s      Stack 
   * @param d      Directory
   */
  void ProcessCent(UShort_t c1, UShort_t c2, TCollection* c, 
		   const TString& name, THStack* s, TDirectory* d) 
  {
    // @todo reimplement to get from GSE's 
    TMultiGraph* others = GetReference(c1, c2);
    if (!others) return;

    TGraph* alice = GetGraph(others->GetListOfGraphs(), 
			     Form("alice_pbpb2760"));
    if (!alice) return;

    TString folderName(Form("cent%03d_%03d", c1, c2));
    TCollection* centBin = GetCollection(c, folderName);
    if (!centBin) return;

    TH1* ana = GetH1(centBin, Form("dndeta%s", name.Data()));
    if (!ana) return;

    TH1* ratio = RatioHG(ana, alice);
    if (!ratio) return;

    ::Info("", "Adding %s/%s", name.Data(), folderName.Data());
    ratio->SetName(folderName);
    ratio->SetTitle(Form("%3d%% - %3d%%", c1, c2));

    // if (c1 == 0) ratio->Write("default");

    TDirectory* store = d->mkdir(folderName);
    store->cd();
    ana->Write();
    alice->Write();
    d->cd();

    s->Add(ratio);
  }

  Bool_t fRefLoaded;
};

void GenerateEmpirical(const Char_t*  fileName="PbPb_2760_dndeta_nosec_CENT_20140513_1349/forward_dndeta.root")
{
  EmpiricalMaker m;
  m.Run(fileName);
}
