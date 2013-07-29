#ifndef UNFOLDMULTDISTS_H
#define UNFOLDMULTDISTS_H
#include <TFile.h>
#include <TList.h>
#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TLatex.h>
#include <TParameter.h>
#include <TCanvas.h>
#include <THStack.h>
#include <TMultiGraph.h>
#include <TRegexp.h>
#include <TSystem.h>
#include "RooUnfold.h"
#include "RooUnfoldResponse.h"

/** 
 * Structure to do unfolding of multiplcity distributions 
 * 
 * It takes as input 2 files: the output file (forward_multdists.root)
 * of the analysis train MakeMultDistsTrain.C run on real and MC data.
 *
 * It generates the output file forward_unfolded.root with histograms
 * and stacks of the input data and the unfolding results.
 */
struct Unfolder 
{
  enum { 
    kMeasuredColor = kBlue+1,
    kUnfoldedColor = kRed+1,
    kTruthColor    = kGreen+2,
    kCMSColor      = kYellow+1,
    kALICEColor    = kMagenta+1,
    kSysColor      = kBlue-10
  };
  /** 
   * Constructor
   */
  Unfolder() {}
  virtual ~Unfolder() {}
  /** 
   * Run this code 
   * 
   * @param realFile Output file of MakeMultDistsTrain on real data
   * @param mcFile   Output file of MakeMultDistsTrain on MC data
   */
  void Run(const TString& method="Bayes", 
	   Double_t       regParam=-1e30, 
	   const TString& realFile="forward_multdists.root", 
	   const TString& mcFile="forward_mcmultdists.root")
  {
    if (!gROOT->GetClass("RooUnfold")) gSystem->Load("libRooUnfold.so");

    // --- Open files ------------------------------------------------
    TFile* realF = TFile::Open(realFile,"READ");
    if (!realF) { 
      Error("Run", "Couldn't open %s", realFile.Data());
      return;
    }
    TFile* mcF = TFile::Open(mcFile,"READ");
    if (!mcF) { 
      Error("Run", "Couldn't open %s", mcFile.Data());
      return;
    }
    TFile* outF = TFile::Open("forward_unfolded.root", "RECREATE");
    
    // --- Get top-level containers ----------------------------------
    TCollection* realTop = GetCollection(realF, "ForwardMultResults");
    if (!realTop) {
      Error("Run", "Failed to get real collection");
      return;
    }
    TCollection* mcTop   = GetCollection(mcF,   "ForwardMultResults");
    if (!mcTop) { 
      Error("Run", "Failed to get MC collection");
      return;
    }

    // --- Decode the method -----------------------------------------
    struct Method { 
      UInt_t  id;
      TString name;
    };
    const Method methods[] = { {RooUnfold::kNone,    "None"},
			       {RooUnfold::kBayes,   "Bayes"},
			       {RooUnfold::kSVD,     "SVD"},
			       {RooUnfold::kBinByBin,"BinByBin"},
			       {RooUnfold::kTUnfold, "TUnfold"},
			       {RooUnfold::kInvert,  "Invert"},
			       {RooUnfold::kDagostini,"Dagostini"}, 
			       {0xDeadBeef,           "unknown"} };
    const Method* pMethod = methods;
    while (pMethod->id != 0xDeadBeef) {
      if (method.BeginsWith(pMethod->name, TString::kIgnoreCase)) break;
      pMethod++;
    }
    if (pMethod->id == 0xDeadBeef) {
      Error("Run", "Unknown unfolding method: %s", method.Data());
      return;
    }
    TNamed* methText = new TNamed("method", pMethod->name.Data());
    methText->SetUniqueID(pMethod->id);
    outF->cd();
    methText->Write();
    
    // --- Loop over the kinds of bins we have -----------------------
    const char*  types[] = { "symmetric", "negative", "positive", "other", 0 };
    const char** ptype   = types;
    Double_t     regParm = regParam;
    while ((*ptype)) { 
      TCollection* realType = GetCollection(realTop, *ptype);
      TCollection* mcType   = GetCollection(mcTop, *ptype);
      TDirectory*  outType  = outF->mkdir(*ptype);
      if(realType && mcType) {
	// Restore default 
	regParm = regParam; 
	// regParm is possibly modified here.
	ScanType(pMethod->id, regParm, realType, mcType, outType);
      }
      // Info("run", "%s gave regularisation parameter %f", *ptype, regParm);
      ptype++;
      outF->cd();
    }
    // --- Write used regularisation parameter -----------------------
    // We get the value that's possibly modified
    // Info("Run", "Regularisation parameter was %f, now %f", regParam,regParm);
    TParameter<double>* regP = new TParameter<double>("regParam", regParm);
    regP->Write();

    // --- Close all files -------------------------------------------
    // outF->ls();
    outF->Write();
    outF->Close();
    realF->Close();
    mcF->Close();
  }
  /** 
   * Get an object from a collection 
   * 
   * @param c    Collection
   * @param name Name of object
   * @param cl   Possible class to check against
   * 
   * @return Pointer to object or null
   */
  TObject* GetObject(TCollection* c, const TString& name, TClass* cl) const
  {
    TObject* o = c->FindObject(name);
    if (!o) { 
      Warning("GetObject", "%s not found in %s", name.Data(), c->GetName());
      return 0;
    }
    if (cl && !o->IsA()->InheritsFrom(cl)) {
      Warning("GetCollection", "%s is not a %s but a %s", 
	      name.Data(), cl->GetName(), o->ClassName());
      return 0;
    }
    return o;
  }
  /** 
   * Get a collection
   * 
   * @param d 
   * @param name 
   * 
   * @return 
   */   
  TCollection* GetCollection(TDirectory* d, const TString& name) const
  {
    TObject* o = d->Get(name);
    if (!o) { 
      Warning("GetCollection", "%s not found in %s", name.Data(), d->GetName());
      return 0;
    }
    if (!o->IsA()->InheritsFrom(TCollection::Class())) { 
      Warning("GetCollection", "%s is not a collection", name.Data());
      return 0;
    }
    return static_cast<TCollection*>(o);
  }
  /** 
   * Get a collection 
   * 
   * @param c 
   * @param name 
   * 
   * @return 
   */
  TCollection* GetCollection(TCollection* c, const TString& name) const
  {
    return static_cast<TCollection*>(GetObject(c, name, TCollection::Class()));
  }
  /** 
   * Get a 1D histogram from a collection
   * 
   * @param c    Collection
   * @param name Nanme of histogram
   * 
   * @return Pointer to object or null
   */
  TH1* GetH1(TCollection* c, const TString& name) 
  {
    return static_cast<TH1*>(GetObject(c, name, TH1::Class()));
  }
  /** 
   * Get a 2D histogram from a collection
   * 
   * @param c    Collection
   * @param name Nanme of histogram
   * 
   * @return Pointer to object or null
   */
  TH2* GetH2(TCollection* c, const TString& name) 
  {
    return static_cast<TH2*>(GetObject(c, name, TH2::Class()));
  }
  TH1* Ratio(const TH1* num, const TGraph* denom,
	     Double_t etaMin, Double_t etaMax) 
  {
    TH1* ret = static_cast<TH1*>(num->Clone(Form("%s_%s",
						 num->GetName(), 
						 denom->GetName())));
    ret->SetTitle(Form("%+5.1f<#eta<%+5.1f to %s", etaMin, etaMax,  
		       denom->GetTitle()));
    ret->SetDirectory(0);
    ret->SetMarkerColor(denom->GetMarkerColor());
    for (Int_t i = 1; i <= ret->GetNbinsX(); i++) { 
      Double_t x    = ret->GetXaxis()->GetBinCenter(i);
      Double_t numY = ret->GetBinContent(i);
      Double_t numE = ret->GetBinError(i);
      Double_t denY = denom->Eval(x);
      
      if (denY <= 0) {
	ret->SetBinContent(i,0);
	ret->SetBinError(i,0);
	continue;
      }
      
      ret->SetBinContent(i,(numY-denY)/denY);
      ret->SetBinError(i,numE/denY);
    }
    return ret;
  }
  /** 
   * Scan type (symmetric, negative, ...) list for bins 
   * 
   * @param real Real list
   * @param mc   MC list
   * @param dir  Output directory 
   */
  void ScanType(UShort_t     method,
		Double_t&    regParam,
		TCollection* real, 
		TCollection* mc, 
		TDirectory*  dir) 
  {
    // --- Create container stack ------------------------------------
    TString tit(real->GetName());
    tit[0] = toupper(tit[0]);
    tit.Append(" bins");
    
    THStack* stack = new THStack("all", tit);
    dir->Add(stack);
    
    THStack* ratios = 0;

    TMultiGraph* mg = 0;

    // --- Create list of entries for legend -------------------------
    TList* l = new TList;
    l->SetName("legend");

    // --- Get the bins container ------------------------------------
    TList* le = static_cast<TList*>(real->FindObject("lowEdges"));
    TList* he = static_cast<TList*>(real->FindObject("highEdges"));
    if (!le || !he) 
      Warning("ScanType", "didn't get the bin low/high edge containers");
    else {
      // Copy to output 
      dir->cd();
      le->Clone()->Write("lowEdges", TObject::kSingleKey);
      he->Clone()->Write("highEdges", TObject::kSingleKey);
    }

    // --- Setup for markers -----------------------------------------
    const Int_t   nMarkers   = 7;
    const Int_t   cMarkers[] = { 20,  21,  22,  23,  29,  33,  34  };
    const Int_t   oMarkers[] = { 24,  25,  26,  32,  30,  27,  28  };
    const Float_t sMarkers[] = { 1.1, 1.0, 1.2, 1.2, 1.2, 1.2, 1.0 };
    Int_t iMarker            = 0;

    // --- Containers that allow us to stack the objects properly ----
    TList* mcTruth      = new TList;
    TList* mcMeasured   = new TList;
    TList* realMeasured = new TList;
    TList* mcUnfolded   = new TList;
    TList* realUnfolded = new TList;
    // --- Loop over the contained objects ---------------------------
    static TRegexp regex("[pm][0-9]d[0-9]*_[pm][0-9]d[0-9]*");
    TIter          next(real);
    TObject*       o = 0;
    Int_t          f = 1;
    Double_t       r = regParam;
    while ((o = next())) {
      // if not a collection, don't bother 
      if (!o->IsA()->InheritsFrom(TCollection::Class())) continue;
    
      // If it doesn't match our regular expression, don't bother 
      TString n(o->GetName());
      if (n.Index(regex) == kNPOS) { 
	// Warning("ScanType", "%s in %s doesn't match eta range regexp", 
	//         n.Data(), real->GetName());
	continue;
      }
      // Cast object and find corresponding MC object 
      TCollection* realBin = static_cast<TCollection*>(o);
      TCollection* mcBin   = GetCollection(mc, n.Data());
      if (!mcBin) { 
	Warning("ScanType", "No corresponding MC bin for %s found in %s", 
		n.Data(), mc->GetName());
	continue;
      }
      TDirectory* outBin = dir->mkdir(realBin->GetName());
      
      // Restore default
      r = regParam;
      // Now do the unfolding 
      THStack* bin = UnfoldEtaBin(method, r, realBin, mcBin, outBin);
      if (!bin) { dir->cd(); continue; }

      // Loop over histograms and set properties 
      TIter   nextH(bin->GetHists());
      TH1*    hist    = 0;
      Int_t   cMarker = cMarkers[iMarker % nMarkers];
      Int_t   oMarker = oMarkers[iMarker % nMarkers];
      Float_t sMarker = sMarkers[iMarker % nMarkers];
      TH1*    realH   = 0;
      while ((hist = static_cast<TH1*>(nextH()))) {
	TH1* out = static_cast<TH1*>(hist->Clone());
	if (out->GetMarkerColor() == kUnfoldedColor) {
          if (out->GetMarkerStyle() != 24)  {
	    realUnfolded->Add(out);
	    realH = out;
	  }
	  else 
	    mcUnfolded->Add(out);
	}
	else if (out->GetMarkerColor() == kMeasuredColor) {
	  if (out->GetMarkerStyle() != 24)  
	    realMeasured->Add(out);
	  else 
	    mcMeasured->Add(out);
	}
	else if (out->GetMarkerColor() == kTruthColor) 
	  mcTruth->Add(out);
	else 
	  Warning("", "Unknown color for %s", out->GetName());

	out->SetDirectory(0);
	out->Scale(f);
	out->SetMarkerStyle(out->GetMarkerStyle() == 24 ? oMarker : cMarker);
	out->SetMarkerSize(out->GetMarkerSize() * sMarker);
	out->SetOption(nextH.GetOption());
	// stack->Add(out, nextH.GetOption());
      }
      TString nn(bin->GetTitle());
      nn.Append(Form(" (#times%d)", f));
      TObjString* lee = new TObjString(nn);
      lee->SetUniqueID(cMarker);
      l->Add(lee);

      // Now try to get external data and make a multigraph 
      nn = o->GetName();
      nn.ReplaceAll("p", "+");
      nn.ReplaceAll("m", "-");
      nn.ReplaceAll("d", ".");
      nn.ReplaceAll("_", " ");
      TObjArray*  tokens = nn.Tokenize(" ");
      TObjString* sMin   = static_cast<TObjString*>(tokens->At(0));
      TObjString* sMax   = static_cast<TObjString*>(tokens->At(1));
      Double_t    etaMin = sMin->String().Atof();
      Double_t    etaMax = sMax->String().Atof();
      if (TMath::Abs(etaMax + etaMin) < 1e-6) { 
	// Symmetric bin
	Double_t aEta = TMath::Abs(etaMin);
	TGraphAsymmErrors* g1 = GetOther(0, aEta, 900, f, cMarker); 
	TGraphAsymmErrors* g2 = GetOther(1, aEta, 900, f, cMarker); 
	if (g1 || g2) {
	  if (!mg) mg = new TMultiGraph("other", "Other results");
	  if (!ratios) ratios = new THStack("ratios",tit);
	}
	if (g1) {
	  g1->SetTitle("CMS");
	  mg->Add(g1);
	  if (realH) ratios->Add(Ratio(realH, g1, etaMin, etaMax));
	}
	if (g2) {
	  g2->SetTitle("ALICE");
	  mg->Add(g2);
	  if (realH) ratios->Add(Ratio(realH, g2, etaMin, etaMax));
	}
      }

      // Increment scaling and marker 
      f *= 10;
      iMarker++;

      dir->cd();
    }
    regParam = r;
    // Info("ScanType", "Regularisation parameter was %f, now %f", regParam, r);
    TH1*  tmp = 0;
    TIter nextMCT(mcTruth);
    while ((tmp = static_cast<TH1*>(nextMCT()))) 
      stack->Add(tmp, tmp->GetOption());
    TIter nextMCM(mcMeasured);
    while ((tmp = static_cast<TH1*>(nextMCM()))) 
      stack->Add(tmp, tmp->GetOption());
    TIter nextRM(realMeasured);
    while ((tmp = static_cast<TH1*>(nextRM()))) 
      stack->Add(tmp, tmp->GetOption());
    TIter nextMCU(mcUnfolded);
    while ((tmp = static_cast<TH1*>(nextMCU()))) 
      stack->Add(tmp, tmp->GetOption());
    TIter nextRU(realUnfolded);
    while ((tmp = static_cast<TH1*>(nextRU()))) 
      stack->Add(tmp, tmp->GetOption());
    
    dir->cd();
    if (mg) mg->Write();
    if (ratios) dir->Add(ratios);
    l->Write(l->GetName(), TObject::kSingleKey);
  }
  /** 
   * Do unfolding in an @f$\eta@f$ bin 
   * 
   * @param real Real list
   * @param mc   MC list
   * @param dir  Output directory
   * 
   * @return Stack of histograms on success
   */  
  THStack* UnfoldEtaBin(UInt_t       method, 
			Double_t&    regParam,
			TCollection* real, 
			TCollection* mc, 
			TDirectory* dir)
  {
    TH1*  realRaw  = GetH1(real, "rawDist");
    TH1*  mcRaw    = GetH1(mc,   "rawDist");
    TH1*  mcTruth  = GetH1(mc,   "truth");
    TH2*  response = GetH2(mc,   "response");
    
    if (!realRaw) { 
      Warning("UnfoldEtaBin", "Real raw distribution not found in %s", 
	      real->GetName());
      return 0;
    }
    if (!mcRaw) { 
      Warning("UnfoldEtaBin", "MC raw distribution not found in %s", 
	      mc->GetName());
      return 0;
    }
    if (!mcTruth) { 
      Warning("UnfoldEtaBin", "MC true distribution not found in %s", 
	      mc->GetName());
      return 0;
    }
    if (!response) { 
      Warning("UnfoldEtaBin", "Response matrix not found in %s", 
	      mc->GetName());
      return 0;
    }
    
    Int_t    mN = realRaw->GetNbinsX();
    Double_t mL = realRaw->GetXaxis()->GetXmin();
    Double_t mH = realRaw->GetXaxis()->GetXmax();
    Int_t    tN = mcRaw->GetNbinsX();
    Double_t tL = mcRaw->GetXaxis()->GetXmin();
    Double_t tH = mcRaw->GetXaxis()->GetXmax();
    
    RooUnfoldResponse matrix(mN, mL, mH, tN, tL, tH);
    for (Int_t i = 1; i <= mN; i++) { 
      Double_t mX = response->GetYaxis()->GetBinCenter(i);
      for (Int_t j = 1; j <= tN; j++) { 
	Double_t tX = response->GetXaxis()->GetBinCenter(j);
	matrix.Fill(mX, tX, response->GetBinContent(j, i));
      }
    }

    Double_t regParm = regParam;
    RooUnfold::Algorithm alg = (RooUnfold::Algorithm)method;
    RooUnfold* realUnfold = RooUnfold::New(alg, &matrix, realRaw, regParm);
    realUnfold->SetVerbose(0);
    TH1* resReal = realUnfold->Hreco();
    TH1* outReal = static_cast<TH1*>(resReal->Clone("realUnfolded"));
    resReal->SetDirectory(0);
    regParam = realUnfold->GetRegParm();
    // Info("UnfoldEtaBin", "Used regularization parameter: %f", regParam);
    delete resReal;
    delete realUnfold;

    RooUnfold* mcUnfold = RooUnfold::New(alg, &matrix, mcRaw, regParm);
    mcUnfold->SetVerbose(0);
    TH1* resMC = mcUnfold->Hreco();
    TH1* outMC = static_cast<TH1*>(resMC->Clone("mcUnfolded"));
    resMC->SetDirectory(0);
    delete resMC;
    delete mcUnfold;

    TH1* outRealRaw  = static_cast<TH1*>(realRaw->Clone("realRaw"));
    TH1* outMcRaw    = static_cast<TH1*>(mcRaw->Clone("mcRaw"));
    TH1* outMcTruth  = static_cast<TH1*>(mcTruth->Clone("mcTruth"));
    TH2* outResponse = static_cast<TH2*>(matrix.Hresponse()->Clone("response"));

    outRealRaw ->SetTitle("Measured");
    outMcRaw   ->SetTitle("Measured (MC)");
    outMcTruth ->SetTitle("Truth (MC)");
    outResponse->SetTitle("Response matrix");
    outReal    ->SetTitle("Unfolded");
    outMC      ->SetTitle("Unfolded (MC)");

    outRealRaw ->SetDirectory(dir);
    outMcRaw   ->SetDirectory(dir);
    outMcTruth ->SetDirectory(dir);
    outResponse->SetDirectory(dir);
    outReal    ->SetDirectory(dir);
    outMC      ->SetDirectory(dir);
    
    outMcRaw  ->SetMarkerStyle(24); // Open circle
    outRealRaw->SetMarkerStyle(20); // Circle
    outMcTruth->SetMarkerStyle(24); // Open Circle
    outMC     ->SetMarkerStyle(24); // Open Circle
    outReal   ->SetMarkerStyle(20); // Circle
    
    outMcRaw  ->SetMarkerSize(1.2); // Open circle
    outMC     ->SetMarkerSize(1.2); // Open circle
    outMcTruth->SetMarkerSize(1.6); // Open square

    outMcRaw  ->SetMarkerColor(kMeasuredColor);
    outRealRaw->SetMarkerColor(kMeasuredColor);
    outMcTruth->SetMarkerColor(kTruthColor);
    outMC     ->SetMarkerColor(kUnfoldedColor);
    outReal   ->SetMarkerColor(kUnfoldedColor);

    outMcRaw  ->SetFillColor(kSysColor);
    outRealRaw->SetFillColor(kSysColor);
    outMcTruth->SetFillColor(kSysColor);
    outMC     ->SetFillColor(kSysColor);
    outReal   ->SetFillColor(kSysColor);

    outMcRaw  ->SetFillStyle(1001);
    outRealRaw->SetFillStyle(0);
    outMcTruth->SetFillStyle(1001);
    outMC     ->SetFillStyle(0);
    outReal   ->SetFillStyle(0);

    outMcRaw  ->SetLineColor(kBlack);
    outRealRaw->SetLineColor(kBlack);
    outMcTruth->SetLineColor(kBlack);
    outMC     ->SetLineColor(kBlack);
    outReal   ->SetLineColor(kBlack);

    TString tit(real->GetName());
    tit.ReplaceAll("m", "-");
    tit.ReplaceAll("p", "+");
    tit.ReplaceAll("d", ".");
    tit.ReplaceAll("_", "<#it{#eta}<");
    
    THStack* stack = new THStack("all", tit);
    stack->Add(outMcRaw,      "E2");
    stack->Add(outRealRaw,    "E1");
    stack->Add(outMcTruth,    "E2");
    stack->Add(outMC,         "E1");
    stack->Add(outReal,       "E1");

    dir->Add(stack);

    return stack;
  }
  /** 
   * Drwa output 
   * 
   * @param which 
   */
  void DrawType(const char* which="symmetric") 
  {
    const char* fname = "forward_unfolded.root";
    TFile*   outF = 0;
    TObject* outO = gROOT->GetListOfFiles()->FindObject(fname);
    if (outO) outF = static_cast<TFile*>(outO);
    else      outF = TFile::Open(fname, "READ");
    
    if (!outF) { 
      Warning("DrawType", "Failed to open file %s", fname);
      return;
    }
    
    TNamed* method = 0;
    outF->GetObject("method", method);
    TParameter<double>* regParam = 0;
    outF->GetObject("regParam", regParam);
    if (!method) Warning("DrawType", "Didn't find the method string");
    if (!regParam) 
      Warning("DrawType", "Didn't find the regularization parameter");
    

    THStack* stack = 0;
    outF->GetObject(Form("/%s/all", which), stack);
    if (!stack) { 
      Warning("DrawType", "Failed to get /%s/all from file %s", which, fname);
      return;
    }
    
    TList* leg = 0;
    outF->GetObject(Form("/%s/legend", which), leg);
    if (!leg) { 
      Warning("DrawType", "Failed to get /%s/legend from file %s", 
	      which, fname);
      return;
    }

    TMultiGraph* other = 0;
    outF->GetObject(Form("/%s/other", which), other);
    // if (!other) Warning("DrawType", "No other data found for %s", which);

    THStack* ratios = 0;
    outF->GetObject(Form("/%s/ratios", which), ratios);
    // if (!ratios) Warning("DrawType", "No ratios data found for %s", which);

    TCanvas* c = new TCanvas(which, Form("%s P(#it{N}_{ch})", which));
    c->SetLogy();
    c->SetTopMargin(0.01);
    c->SetRightMargin(0.02);
    c->cd();
    
    THStack* drawn = static_cast<THStack*>(stack->Clone(which));
    drawn->Draw("nostack");
    drawn->GetXaxis()->SetTitle("#it{#eta}");
    drawn->GetYaxis()->SetTitle("P(#it{N}_{ch})");
    drawn->GetHistogram()->SetDirectory(0);
    TIter nextH(drawn->GetHists());
    TH1*  h = 0;
    while ((h = static_cast<TH1*>(nextH()))) h->SetDirectory(0);

    Bool_t        hasCMS   = false;
    Bool_t        hasALICE = false;
    if (other) {
      TIter nextG(other->GetListOfGraphs());
      TGraph* g = 0;
      while ((g = static_cast<TGraph*>(nextG()))) {
	g->Draw("p same");
	if (g->GetMarkerColor() == kALICEColor) hasALICE = true;
	if (g->GetMarkerColor() == kCMSColor)   hasCMS   = true;
      }
    }
    TLegend* l = new TLegend(.65, .75, .975, .975);
    l->SetBorderSize(0);
    l->SetFillColor(0);
    l->SetFillStyle(0);
    
    TIter         next(leg);
    TObject*      o = 0;
    TLegendEntry* e = 0;
    Int_t         d = 0;
    while ((o = next())) { 
      e = l->AddEntry(Form("dummy%02d", d++), o->GetName(), "p");
      e->SetMarkerStyle(o->GetUniqueID());
      e->SetLineColor(kBlack);
      e->SetMarkerSize(1.3);
    }
    l->Draw();

    Double_t x1 = other ? .12 : .7;
    Double_t y1 = .1;
    Double_t w1 = other ? .20 : .15;
    Double_t h1 = other ? .3 : .25;
    TLegend* l2 = new TLegend(x1, y1, x1+w1, y1+h1);
    l2->SetBorderSize(0);
    l2->SetFillColor(0);
    l2->SetFillStyle(0);

    const char* oth[] = { "Measured", 
			  "Truth", 
			  "Unfolded", 
			  "CMS", 
			  "ALICE" };
    Int_t       col[] = { kMeasuredColor, 
			  kTruthColor, 
			  kUnfoldedColor, 
			  kCMSColor, 
			  kALICEColor };
    for (Int_t i = 0; i < 5; i++) {
      if (i == 3 && !hasCMS)   continue;
      if (i == 4 && !hasALICE) continue;
      e = l2->AddEntry(Form("dummy%02d", d++), oth[i], "f");
      e->SetFillColor(col[i]);
      e->SetFillStyle(1001);
    }
    e = l2->AddEntry(Form("dummy%02d", d++), "Real", "p");
    e->SetMarkerStyle(20);
    e->SetMarkerSize(1.3);
    e->SetMarkerColor(kBlack);
    e = l2->AddEntry(Form("dummy%02d", d++), "MC", "p");
    e->SetMarkerStyle(24);
    e->SetMarkerSize(1.3);
    e->SetMarkerColor(kBlack);
    l2->Draw();


    TString post;
    if (method && regParam) {
      TString mes = method->GetTitle();
      if      (method->GetUniqueID() == RooUnfold::kBayes || 
	       method->GetUniqueID() == RooUnfold::kDagostini)
	mes.Append(Form(" (N_{iter}=%d)", int(regParam->GetVal())));
      else if (method->GetUniqueID() == RooUnfold::kSVD) 
	mes.Append(Form(" (k=%d)", int(regParam->GetVal())));
      else if (method->GetUniqueID() == RooUnfold::kTUnfold) 
	mes.Append(Form(" (#tau=%f)", regParam->GetVal()));
      
      TLatex* meth = new TLatex(.97, .749, mes);
      meth->SetNDC();
      meth->SetTextAlign(33);
      meth->SetTextFont(42);
      meth->SetTextSize(.03);
      meth->Draw();

      post.Form("_%s", method->GetTitle());
      post.ToLower();
    }
    c->Modified();
    c->Update();
    c->cd();
    c->SaveAs(Form("pnch_%s%s.pdf", which, post.Data()));

    if (ratios) {
      c = new TCanvas(Form("r%s", which), 
		      Form("%s P(#it{N}_{ch})", which));
      c->SetTopMargin(0.01);
      c->SetRightMargin(0.02);
      c->cd();
      ratios->Draw("nostack");
      c->BuildLegend();
    }

    outF->cd();
  }
  void DrawAll()
  {
    DrawType("symmetric");
    DrawType("negative");
    DrawType("positive");
    DrawType("other");
  }
  TGraphAsymmErrors* GetOther(UShort_t type, Double_t eta, UShort_t sNN,
			      Int_t factor, Int_t marker)
  {
    TGraphAsymmErrors* g = 0;
    Int_t color = kBlack;
    switch (type) { 
    case 0: g = GetCMS(eta, sNN);   color = kCMSColor; break;
    case 1: g = GetALICE(eta, sNN); color = kALICEColor; break;
    }
    if (!g) {
      // Warning("GetOther", "No other data found for type=%d, eta=%f, sNN=%d", 
      //         type, eta, sNN);
      return g;
    }
    
    g->SetMarkerStyle(marker);
    g->SetMarkerColor(color);
    for (Int_t j = 0; j < g->GetN(); j++) { 
      g->SetPoint(j, g->GetX()[j], g->GetY()[j]*factor);
      g->SetPointError(j, g->GetEXlow()[j], g->GetEXhigh()[j], 
		       g->GetEYlow()[j]*factor, g->GetEYhigh()[j]*factor);
    }
    return g;
  }    

  TGraphAsymmErrors* GetCMS(Double_t eta, UShort_t sNN)
  {
    Double_t aEta = TMath::Abs(eta);
    TGraphAsymmErrors* g = 0;
    switch (sNN) { 
    case 900: 
      if      (aEta <= 0.51)  g = CMSsqrts900eta05();
      else if (aEta <= 1.01)  g = CMSsqrts900eta10();
      else if (aEta <= 1.51)  g = CMSsqrts900eta15();
      else if (aEta <= 2.01)  g = CMSsqrts900eta20();
      else if (aEta <= 2.41)  g = CMSsqrts900eta24();
      break;
    default: 
      return g;
    }
    if (g) {
      g->SetMarkerColor(kCMSColor);
      g->SetLineColor(kBlack);
    }
    return g;
  }
  /** 
   * CMS pp, @f$\sqrt{s} = 900GeV@f$, NSD, @f$|\eta|\lt0.5@f$ 
   * 
   * @return 
   */
  TGraphAsymmErrors* CMSsqrts900eta05()
  {
    // Plot: p8068_d2x1y1
    double x[] = { 0.0,	1.0,
		   2.0,	3.0,
		   4.0,	5.0,
		   6.0,	7.0,
		   8.0,	9.0,
		   10.0,	11.0,
		   12.0,	13.0,
		   14.0,	15.0,
		   16.0,	17.0,
		   18.0,	19.0,
		   20.0,	22.0,
		   26.0 };
    double xem[] = { 0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	1.5,
		     2.5 };
    double xep[] = { 0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	1.5,
		     2.5 };
    double y[] = { 0.193,	0.1504,
		   0.13507,	0.11647,
		   0.09566,	0.07567,
		   0.05817,	0.04427,
		   0.0337,	0.02572,
		   0.01948,	0.01457,
		   0.01069,	0.00769,
		   0.0055,	0.00396,
		   0.00289,	0.002112,
		   0.001544,	0.001118,
		   7.52E-4,	3.35E-4,
		   9.8E-5 };
    double yem[] = { 0.02189155088156159,	0.012180722474467595,
		     0.007909551188278638,	0.006072832946821442,
		     0.004838605170914444,	0.003828276374558138,
		     0.003016836091006603,	0.002365797962633327,
		     0.0018847015678881366,	0.0015610893632332517,
		     0.001308166656049603,	0.0011140017953306899,
		     9.219544457292888E-4,	7.203471385380802E-4,
		     5.508175741568165E-4,	4.3011626335213136E-4,
		     3.3837848631377264E-4,	2.746961958236772E-4,
		     2.2377220560203628E-4,	1.9070920271449933E-4,
		     1.6011246047700347E-4,	8.134494452638099E-5,
		     2.5059928172283337E-5 };
    double yep[] = { 0.031663227883461285,	0.015861904047118684,
		     0.010948538715280683,	0.008490583018850945,
		     0.007046453008429135,	0.005675746646917919,
		     0.004395201929377079,	0.003365070578754627,
		     0.0026049952015310893,	0.0020535335400231475,
		     0.0016522106403240478,	0.0013497407158413798,
		     0.001118033988749895,	9.060905032059435E-4,
		     7.061161377563892E-4,	5.457105459856901E-4,
		     4.148493702538308E-4,	3.257759966602819E-4,
		     2.6689698387205504E-4,	2.292705825002414E-4,
		     1.904205871222962E-4,	9.646242791885347E-5,
		     3.417601498127012E-5 };
    double ysm[] = { 0.002,	0.0014,
		     7.9E-4,	7.2E-4,
		     6.1E-4,	5.4E-4,
		     4.7E-4,	4.1E-4,
		     3.6E-4,	3.1E-4,
		     2.7E-4,	2.3E-4,
		     2.0E-4,	1.7E-4,
		     1.5E-4,	1.3E-4,
		     1.1E-4,	9.7E-5,
		     8.5E-5,	7.1E-5,
		     5.6E-5,	2.9E-5,
		     1.2E-5 };
    double ysp[] = { 0.002,	0.0014,
		     7.9E-4,	7.2E-4,
		     6.1E-4,	5.4E-4,
		     4.7E-4,	4.1E-4,
		     3.6E-4,	3.1E-4,
		     2.7E-4,	2.3E-4,
		     2.0E-4,	1.7E-4,
		     1.5E-4,	1.3E-4,
		     1.1E-4,	9.7E-5,
		     8.5E-5,	7.1E-5,
		     5.6E-5,	2.9E-5,
		     1.2E-5 };
    int np = 23;
    TGraphAsymmErrors* g = new TGraphAsymmErrors(np, 
						 x, y, 
						 xem, 
						 xep, 
						 yem, 
						 yep);
    g->SetName("/HepData/8068/d2x1y1");
    g->SetTitle("/HepData/8068/d2x1y1");

    return g;
  }
  TGraphAsymmErrors* CMSsqrts900eta10()
  {

    // Plot: p8068_d3x1y1
    double x[] = { 0.0,	1.0,
		   2.0,	3.0,
		   4.0,	5.0,
		   6.0,	7.0,
		   8.0,	9.0,
		   10.0,	11.0,
		   12.0,	13.0,
		   14.0,	15.0,
		   16.0,	17.0,
		   18.0,	19.0,
		   20.0,	21.0,
		   22.0,	23.0,
		   24.0,	25.0,
		   26.0,	27.0,
		   28.0,	29.0,
		   30.0,	31.0,
		   32.0,	33.0,
		   34.0,	36.0,
		   40.0,	45.0,
		   50.0,	55.0 };
    double xem[] = { 0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	1.5,
		     2.5,	2.5,
		     2.5,	2.5 };
    double xep[] = { 0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	1.5,
		     2.5,	2.5,
		     2.5,	2.5 };
    double y[] = { 0.1076,	0.0655,
		   0.07574,	0.08168,
		   0.0813,	0.07622,
		   0.06843,	0.06034,
		   0.05259,	0.0458,
		   0.03987,	0.03472,
		   0.03035,	0.02661,
		   0.02328,	0.02018,
		   0.01728,	0.0147,
		   0.0125,	0.01067,
		   0.00916,	0.00786,
		   0.00671,	0.0057,
		   0.00478,	0.00394,
		   0.00319,	0.00251,
		   0.00196,	0.00156,
		   0.001301,	0.001121,
		   9.84E-4,	8.25E-4,
		   6.39E-4,	3.78E-4,
		   1.96E-4,	6.0E-5,
		   1.25E-5,	1.6E-6 };
    double yem[] = { 0.015802847844613326,	0.009464142856064674,
		     0.008214195030555337,	0.00676426640516176,
		     0.0058916636020736966,	0.005200240379059414,
		     0.004441733895676326,	0.003671511950137164,
		     0.0029149442533262966,	0.0024150569351466646,
		     0.0020554318281081475,	0.001767059704707229,
		     0.0015500322577288513,	0.0013800362314084368,
		     0.0012684636376341263,	0.0011763077828527702,
		     0.001084158659975559,	9.82344135219425E-4,
		     8.805679985100526E-4,	7.692203845452875E-4,
		     6.705221845696084E-4,	5.913543776789007E-4,
		     5.375872022286245E-4,	5.154609587543949E-4,
		     5.124451190127582E-4,	4.904079934095691E-4,
		     4.4922154890432404E-4,	3.7E-4,
		     3.0083217912982643E-4,	2.4166091947189146E-4,
		     2.1108529081866412E-4,	1.8252944967867517E-4,
		     1.662077013859466E-4,	1.4663219291819924E-4,
		     1.2229881438509532E-4,	7.826237921249263E-5,
		     4.74236228055175E-5,	2.1095023109728988E-5,
		     6.161980201201558E-6,	1.3601470508735443E-6 };
    double yep[] = { 0.02725949375905576,	0.011154371340420759,
		     0.009740395269186974,	0.008248078564126314,
		     0.007285581651453781,	0.006454339625399333,
		     0.005595158621522718,	0.004804039966528172,
		     0.00405504623894722,	0.0034544174617437305,
		     0.0029644561052577585,	0.002545584412271571,
		     0.0021980445855350615,	0.0019087430418995638,
		     0.0016690416411821486,	0.001468911161370898,
		     0.0013179529581893276,	0.0011866338946785568,
		     0.001065129100156408,	9.339164844888434E-4,
		     8.052328855678958E-4,	6.868041933477111E-4,
		     5.94810894318522E-4,	5.440588203494177E-4,
		     5.220153254455275E-4,	5.192301994298867E-4,
		     4.972926703662542E-4,	4.368065933568311E-4,
		     3.5735136770411276E-4,	2.8792360097775935E-4,
		     2.4200206610688264E-4,	2.0310588371585893E-4,
		     1.670748335327616E-4,	1.364587849865299E-4,
		     1.1970797801316335E-4,	9.470480452437458E-5,
		     5.758472019555187E-5,	2.195449840010015E-5,
		     7.481310045707236E-6,	1.7029386365926402E-6 };
    double ysm[] = { 0.0018,	0.0011,
		     6.3E-4,	6.8E-4,
		     6.1E-4,	5.6E-4,
		     5.3E-4,	4.8E-4,
		     4.5E-4,	4.1E-4,
		     3.8E-4,	3.6E-4,
		     3.5E-4,	3.3E-4,
		     3.1E-4,	2.9E-4,
		     2.7E-4,	2.5E-4,
		     2.3E-4,	2.1E-4,
		     2.0E-4,	1.9E-4,
		     1.7E-4,	1.6E-4,
		     1.5E-4,	1.4E-4,
		     1.3E-4,	1.2E-4,
		     1.1E-4,	1.0E-4,
		     9.4E-5,	8.6E-5,
		     8.3E-5,	7.5E-5,
		     6.1E-5,	3.5E-5,
		     2.0E-5,	1.1E-5,
		     4.6E-6,	1.1E-6 };
    double ysp[] = { 0.0018,	0.0011,
		     6.3E-4,	6.8E-4,
		     6.1E-4,	5.6E-4,
		     5.3E-4,	4.8E-4,
		     4.5E-4,	4.1E-4,
		     3.8E-4,	3.6E-4,
		     3.5E-4,	3.3E-4,
		     3.1E-4,	2.9E-4,
		     2.7E-4,	2.5E-4,
		     2.3E-4,	2.1E-4,
		     2.0E-4,	1.9E-4,
		     1.7E-4,	1.6E-4,
		     1.5E-4,	1.4E-4,
		     1.3E-4,	1.2E-4,
		     1.1E-4,	1.0E-4,
		     9.4E-5,	8.6E-5,
		     8.3E-5,	7.5E-5,
		     6.1E-5,	3.5E-5,
		     2.0E-5,	1.1E-5,
		     4.6E-6,	1.1E-6 };
    int np = 40;
    TGraphAsymmErrors* g = new TGraphAsymmErrors(np,
						 x,
						 y,
						 xem,
						 xep,
						 yem,
						 yep);
    g->SetName("/HepData/8068/d3x1y1");
    g->SetTitle("/HepData/8068/d3x1y1");
    return g;
  }
  TGraphAsymmErrors* CMSsqrts900eta15()
  {
    // Plot: p8068_d4x1y1
    double x[] = { 0.0,	1.0,
		   2.0,	3.0,
		   4.0,	5.0,
		   6.0,	7.0,
		   8.0,	9.0,
		   10.0,	11.0,
		   12.0,	13.0,
		   14.0,	15.0,
		   16.0,	17.0,
		   18.0,	19.0,
		   20.0,	21.0,
		   22.0,	23.0,
		   24.0,	25.0,
		   26.0,	27.0,
		   28.0,	29.0,
		   30.0,	31.0,
		   32.0,	33.0,
		   34.0,	35.0,
		   36.0,	37.0,
		   38.0,	39.0,
		   40.0,	41.0,
		   42.0,	43.0,
		   44.0,	45.0,
		   46.5,	49.0,
		   53.0,	58.0,
		   63.0,	68.0 };
    double xem[] = { 0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     1.0,	1.5,
		     2.5,	2.5,
		     2.5,	2.5 };
    double xep[] = { 0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     1.0,	1.5,
		     2.5,	2.5,
		     2.5,	2.5 };
    double y[] = { 0.0753,	0.03557,
		   0.04347,	0.05068,
		   0.05566,	0.05786,
		   0.05835,	0.05688,
		   0.05348,	0.0489,
		   0.04422,	0.04002,
		   0.03647,	0.03337,
		   0.03055,	0.02797,
		   0.02556,	0.02329,
		   0.02122,	0.01936,
		   0.01765,	0.01604,
		   0.01451,	0.01304,
		   0.01168,	0.01045,
		   0.00937,	0.0084,
		   0.00751,	0.00667,
		   0.00592,	0.00526,
		   0.00466,	0.00413,
		   0.00362,	0.00316,
		   0.00272,	0.00233,
		   0.00203,	0.00178,
		   0.00156,	0.00136,
		   0.0012,	0.001053,
		   9.13E-4,	7.84E-4,
		   5.85E-4,	3.73E-4,
		   2.13E-4,	9.6E-5,
		   3.6E-5,	7.4E-6 };
    double yem[] = { 0.014524806367039803,	0.0061507154055443014,
		     0.006900384047283165,	0.0068944397886992964,
		     0.00601252858621063,	0.004686747272896204,
		     0.0042880881520789655,	0.0039679339712248235,
		     0.0035288099977187778,	0.0031197756329582422,
		     0.0027127108213003464,	0.0023229722340140013,
		     0.002033937068839643,	0.0017866449003649271,
		     0.001591257364476281,	0.0014504137340772805,
		     0.0013412307780542466,	0.0012490796611905903,
		     0.0011594826432508594,	0.0010893117092916976,
		     0.0010384603988597735,	9.972462083156796E-4,
		     9.588013350011566E-4,	8.984430978086481E-4,
		     8.160882305241266E-4,	7.178439941937245E-4,
		     6.390618123468182E-4,	5.88727441181401E-4,
		     5.571355310873647E-4,	5.255473337388365E-4,
		     4.846648326421054E-4,	4.6615448083226656E-4,
		     4.4384682042344296E-4,	4.3081318457076036E-4,
		     4.0853396431630995E-4,	3.992492955535426E-4,
		     3.676955262170047E-4,	3.2695565448543625E-4,
		     2.9068883707497264E-4,	2.505992817228334E-4,
		     2.1954498400100151E-4,	1.8601075237738274E-4,
		     1.8601075237738274E-4,	1.7042300314218148E-4,
		     1.7998055450520203E-4,	1.6058953888718904E-4,
		     1.4985659811966905E-4,	8.561541917201597E-5,
		     6.324555320336759E-5,	4.657252408878007E-5,
		     2.4166091947189146E-5,	5.590169943749474E-6 };
    double yep[] = { 0.03405304685340212,	0.00724989655098609,
		     0.007817985674072318,	0.007751728839426725,
		     0.006919566460407762,	0.00571192611997039,
		     0.005232991496266739,	0.004813003220443552,
		     0.00435332057170156,	0.0038838769290491172,
		     0.003445692963686695,	0.0030549959083442323,
		     0.0027252339349127445,	0.0024367396249907374,
		     0.002189794510907359,	0.001969568480657629,
		     0.0017906702655709678,	0.0016395731151735808,
		     0.0015002999700059986,	0.0013710215169719256,
		     0.0012712198865656563,	0.0011812705024675763,
		     0.0011035397591387453,	0.0010139033484509261,
		     9.121403400793104E-4,	8.132035415564789E-4,
		     7.242237223399962E-4,	6.545991139621257E-4,
		     6.135144660071188E-4,	5.724508712544685E-4,
		     5.692099788303082E-4,	5.126402247190518E-4,
		     4.531004303683677E-4,	4.4011362169330773E-4,
		     4.1785164831552356E-4,	4.0853396431630995E-4,
		     3.8626415831655934E-4,	3.54682957019364E-4,
		     3.1780497164141406E-4,	2.7730849247724094E-4,
		     2.459674775249769E-4,	2.3706539182259396E-4,
		     2.1954498400100151E-4,	2.0068881383873889E-4,
		     1.5864425612041553E-4,	1.6058953888718904E-4,
		     1.4338758663147937E-4,	1.0245974819410792E-4,
		     6.609841147864296E-5,	4.657252408878007E-5,
		     2.5079872407968906E-5,	5.941380311005179E-6 };
    double ysm[] = { 0.0019,	9.3E-4,
		     5.3E-4,	5.8E-4,
		     5.2E-4,	5.0E-4,
		     4.9E-4,	4.7E-4,
		     4.5E-4,	4.3E-4,
		     4.2E-4,	3.9E-4,
		     3.7E-4,	3.6E-4,
		     3.6E-4,	3.4E-4,
		     3.3E-4,	3.1E-4,
		     3.0E-4,	2.9E-4,
		     2.8E-4,	2.7E-4,
		     2.7E-4,	2.6E-4,
		     2.4E-4,	2.3E-4,
		     2.2E-4,	2.1E-4,
		     2.0E-4,	1.9E-4,
		     1.8E-4,	1.8E-4,
		     1.7E-4,	1.6E-4,
		     1.5E-4,	1.5E-4,
		     1.4E-4,	1.3E-4,
		     1.3E-4,	1.2E-4,
		     1.1E-4,	1.1E-4,
		     1.1E-4,	1.0E-4,
		     8.8E-5,	7.5E-5,
		     5.6E-5,	3.3E-5,
		     2.0E-5,	1.2E-5,
		     1.0E-5,	4.1E-6 };
    double ysp[] = { 0.0019,	9.3E-4,
		     5.3E-4,	5.8E-4,
		     5.2E-4,	5.0E-4,
		     4.9E-4,	4.7E-4,
		     4.5E-4,	4.3E-4,
		     4.2E-4,	3.9E-4,
		     3.7E-4,	3.6E-4,
		     3.6E-4,	3.4E-4,
		     3.3E-4,	3.1E-4,
		     3.0E-4,	2.9E-4,
		     2.8E-4,	2.7E-4,
		     2.7E-4,	2.6E-4,
		     2.4E-4,	2.3E-4,
		     2.2E-4,	2.1E-4,
		     2.0E-4,	1.9E-4,
		     1.8E-4,	1.8E-4,
		     1.7E-4,	1.6E-4,
		     1.5E-4,	1.5E-4,
		     1.4E-4,	1.3E-4,
		     1.3E-4,	1.2E-4,
		     1.1E-4,	1.1E-4,
		     1.1E-4,	1.0E-4,
		     8.8E-5,	7.5E-5,
		     5.6E-5,	3.3E-5,
		     2.0E-5,	1.2E-5,
		     1.0E-5,	4.1E-6 };
    int np = 52;
    TGraphAsymmErrors* g = new TGraphAsymmErrors(np,
						 x,
						 y,
						 xem,
						 xep,
						 yem,
						 yep);
    g->SetName("/HepData/8068/d4x1y1");
    g->SetTitle("/HepData/8068/d4x1y1");
    return g;
  }
  TGraphAsymmErrors* CMSsqrts900eta24()
  {
    // Plot: p8068_d5x1y1
    double x[] = { 0.0,	1.0,
		   2.0,	3.0,
		   4.0,	5.0,
		   6.0,	7.0,
		   8.0,	9.0,
		   10.0,	11.0,
		   12.0,	13.0,
		   14.0,	15.0,
		   16.0,	17.0,
		   18.0,	19.0,
		   20.0,	21.0,
		   22.0,	23.0,
		   24.0,	25.0,
		   26.0,	27.0,
		   28.0,	29.0,
		   30.0,	31.0,
		   32.0,	33.0,
		   34.0,	35.0,
		   36.0,	37.0,
		   38.0,	39.0,
		   40.0,	41.0,
		   42.0,	43.0,
		   44.0,	45.0,
		   46.0,	47.0,
		   48.0,	49.0,
		   50.0,	51.0,
		   52.0,	53.0,
		   54.0,	55.0,
		   56.0,	57.5,
		   60.0,	66.5,
		   76.5,	86.5 };
    double xem[] = { 0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	1.0,
		     1.5,	5.0,
		     5.0,	5.0 };
    double xep[] = { 0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	1.0,
		     1.5,	5.0,
		     5.0,	5.0 };
    double y[] = { 0.0561,	0.02321,
		   0.02645,	0.031,
		   0.0362,	0.04081,
		   0.04402,	0.04548,
		   0.04592,	0.04515,
		   0.04342,	0.04101,
		   0.03836,	0.03573,
		   0.03316,	0.03079,
		   0.02868,	0.02679,
		   0.02506,	0.02343,
		   0.02187,	0.02039,
		   0.01901,	0.0178,
		   0.01668,	0.01556,
		   0.01445,	0.01334,
		   0.01229,	0.01132,
		   0.01045,	0.00967,
		   0.00896,	0.00829,
		   0.00763,	0.00699,
		   0.00637,	0.00578,
		   0.00527,	0.00482,
		   0.00441,	0.00412,
		   0.00363,	0.00328,
		   0.00292,	0.00261,
		   0.00233,	0.00209,
		   0.00185,	0.00164,
		   0.00146,	0.00131,
		   0.00117,	0.00106,
		   9.7E-4,	8.9E-4,
		   8.2E-4,	6.63E-4,
		   4.53E-4,	2.12E-4,
		   7.7E-5,	1.12E-5 };
    double yem[] = { 0.01767766952966369,	0.003660737630587584,
		     0.005079094407470687,	0.006069810540700591,
		     0.006147235150862541,	0.005458580401532985,
		     0.004473712105176193,	0.0035970126494078384,
		     0.003456848275524976,	0.003228761372415125,
		     0.0029797483115189447,	0.002739361239413305,
		     0.0025401181074902798,	0.002313114783144148,
		     0.0020573040611440983,	0.001820109886792553,
		     0.0016324827717314506,	0.001484318025222358,
		     0.0013727709204379296,	0.001287866452703851,
		     0.001217538500417954,	0.0011279184367674819,
		     0.0010412012293500234,	9.741663102366044E-4,
		     9.360555539069249E-4,	9.139474820797965E-4,
		     9.013878188659973E-4,	8.697700845625814E-4,
		     8.190848551890091E-4,	7.433034373659253E-4,
		     6.741661516273269E-4,	6.239390995922599E-4,
		     5.961543424315552E-4,	5.8309518948453E-4,
		     5.738466694161429E-4,	5.515432893255072E-4,
		     5.292447448959697E-4,	4.976946855251721E-4,
		     4.609772228646444E-4,	4.295346318982906E-4,
		     4.0249223594996216E-4,	4.386342439892262E-4,
		     3.712142238654117E-4,	3.6674241641784496E-4,
		     3.6249137920783716E-4,	3.5341194094144583E-4,
		     3.3105890714493693E-4,	3.130495168499705E-4,
		     2.8178005607210744E-4,	2.5553864678361276E-4,
		     2.3323807579381201E-4,	2.1633307652783938E-4,
		     1.8601075237738274E-4,	1.6401219466856724E-4,
		     1.562049935181331E-4,	1.562049935181331E-4,
		     1.5480633061990717E-4,	1.395886814895821E-4,
		     1.0756393447619885E-4,	6.676076692189808E-5,
		     2.9966648127543395E-5,	6.551335741663679E-6 };
    double yep[] = { 0.09323352401362935,	0.004710169848317574,
		     0.005736898116578331,	0.006618164398078972,
		     0.006715772479767314,	0.006096630216767292,
		     0.00519042387479096,	0.004362247127341595,
		     0.004172217635742412,	0.003913693907295255,
		     0.003644283194264683,	0.0033638668225719043,
		     0.0031145144083789367,	0.0028765604460883488,
		     0.002629087294100369,	0.0023903974564912843,
		     0.002181604913819182,	0.0020026232796010335,
		     0.0018514858897652987,	0.0017280335644888382,
		     0.00160822883943797,	0.0014787156589419076,
		     0.0013612494260788505,	0.0012539936203984452,
		     0.0011763077828527702,	0.0011157060544785082,
		     0.001084158659975559,	0.001033247308247159,
		     9.630160954002794E-4,	9.244457799135652E-4,
		     8.35224520712844E-4,	6.989277502002622E-4,
		     6.519202405202648E-4,	6.203224967708329E-4,
		     6.10982814815605E-4,	5.88727441181401E-4,
		     5.664803615307418E-4,	5.348831648126533E-4,
		     4.976946855251721E-4,	4.477722635447622E-4,
		     4.204759208325728E-4,	3.9357337308308855E-4,
		     3.981205847478877E-4,	3.8483762809787716E-4,
		     3.807886552931954E-4,	3.716180835212409E-4,
		     3.584689665786984E-4,	3.4928498393145964E-4,
		     3.1780497164141406E-4,	2.9068883707497264E-4,
		     2.505992817228334E-4,	2.2472205054244234E-4,
		     2.0248456731316588E-4,	1.8027756377319944E-4,
		     1.6401219466856724E-4,	1.6401219466856724E-4,
		     1.5640012787718557E-4,	1.4223923509355636E-4,
		     1.0288342918079667E-4,	6.86804193347711E-5,
		     3.361547262794322E-5,	7.071067811865475E-6 };
    double ysm[] = { 0.0025,	8.1E-4,
		     4.4E-4,	4.9E-4,
		     4.6E-4,	4.5E-4,
		     4.6E-4,	4.4E-4,
		     4.3E-4,	4.3E-4,
		     4.2E-4,	4.0E-4,
		     3.9E-4,	3.9E-4,
		     3.9E-4,	3.8E-4,
		     3.7E-4,	3.6E-4,
		     3.4E-4,	3.1E-4,
		     3.0E-4,	2.9E-4,
		     2.9E-4,	2.9E-4,
		     2.9E-4,	2.8E-4,
		     2.7E-4,	2.6E-4,
		     2.5E-4,	2.5E-4,
		     2.4E-4,	2.3E-4,
		     2.3E-4,	2.2E-4,
		     2.2E-4,	2.1E-4,
		     2.0E-4,	1.9E-4,
		     1.9E-4,	1.8E-4,
		     1.8E-4,	1.8E-4,
		     1.7E-4,	1.6E-4,
		     1.5E-4,	1.5E-4,
		     1.4E-4,	1.4E-4,
		     1.3E-4,	1.3E-4,
		     1.2E-4,	1.2E-4,
		     1.1E-4,	1.0E-4,
		     1.0E-4,	1.0E-4,
		     9.4E-5,	6.6E-5,
		     3.7E-5,	1.9E-5,
		     1.3E-5,	3.4E-6 };
    double ysp[] = { 0.0025,	8.1E-4,
		     4.4E-4,	4.9E-4,
		     4.6E-4,	4.5E-4,
		     4.6E-4,	4.4E-4,
		     4.3E-4,	4.3E-4,
		     4.2E-4,	4.0E-4,
		     3.9E-4,	3.9E-4,
		     3.9E-4,	3.8E-4,
		     3.7E-4,	3.6E-4,
		     3.4E-4,	3.1E-4,
		     3.0E-4,	2.9E-4,
		     2.9E-4,	2.9E-4,
		     2.9E-4,	2.8E-4,
		     2.7E-4,	2.6E-4,
		     2.5E-4,	2.5E-4,
		     2.4E-4,	2.3E-4,
		     2.3E-4,	2.2E-4,
		     2.2E-4,	2.1E-4,
		     2.0E-4,	1.9E-4,
		     1.9E-4,	1.8E-4,
		     1.8E-4,	1.8E-4,
		     1.7E-4,	1.6E-4,
		     1.5E-4,	1.5E-4,
		     1.4E-4,	1.4E-4,
		     1.3E-4,	1.3E-4,
		     1.2E-4,	1.2E-4,
		     1.1E-4,	1.0E-4,
		     1.0E-4,	1.0E-4,
		     9.4E-5,	6.6E-5,
		     3.7E-5,	1.9E-5,
		     1.3E-5,	3.4E-6 };
    int np = 62;
    TGraphAsymmErrors* g = new TGraphAsymmErrors(np,
						 x,
						 y,
						 xem,
						 xep,
						 yem,
						 yep);
    g->SetName("/HepData/8068/d5x1y1");
    g->SetTitle("/HepData/8068/d5x1y1");
    return g;
  }
  TGraphAsymmErrors* CMSsqrts900eta20()
  {
    // Plot: p8068_d6x1y1
    double x[] = { 0.0,	1.0,
		   2.0,	3.0,
		   4.0,	5.0,
		   6.0,	7.0,
		   8.0,	9.0,
		   10.0,	11.0,
		   12.0,	13.0,
		   14.0,	15.0,
		   16.0,	17.0,
		   18.0,	19.0,
		   20.0,	21.0,
		   22.0,	23.0,
		   24.0,	25.0,
		   26.0,	27.0,
		   28.0,	29.0,
		   30.0,	31.0,
		   32.0,	33.0,
		   34.0,	35.0,
		   36.0,	37.0,
		   38.0,	39.0,
		   40.0,	41.0,
		   42.0,	43.0,
		   44.0,	45.0,
		   46.0,	47.0,
		   48.0,	49.0,
		   50.0,	51.0,
		   52.0,	53.0,
		   54.0,	55.0,
		   56.0,	57.0,
		   58.0,	59.0,
		   60.0,	61.0,
		   62.5,	64.5,
		   67.5,	74.5,
		   84.5,	94.5 };
    double xem[] = { 0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     1.0,	1.0,
		     2.0,	5.0,
		     5.0,	5.0 };
    double xep[] = { 0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     1.0,	1.0,
		     2.0,	5.0,
		     5.0,	5.0 };
    double y[] = { 0.0494,	0.01789,
		   0.01895,	0.02174,
		   0.0255,	0.02961,
		   0.0335,	0.03646,
		   0.03826,	0.03906,
		   0.0392,	0.03868,
		   0.03753,	0.03596,
		   0.0341,	0.03206,
		   0.03,	0.02802,
		   0.02616,	0.02448,
		   0.023,	0.0217,
		   0.02058,	0.01959,
		   0.01866,	0.01776,
		   0.01685,	0.01595,
		   0.01504,	0.01415,
		   0.01326,	0.01241,
		   0.01161,	0.01087,
		   0.01015,	0.00948,
		   0.00887,	0.00829,
		   0.00777,	0.00734,
		   0.00681,	0.00624,
		   0.00578,	0.00535,
		   0.00494,	0.00457,
		   0.00422,	0.00388,
		   0.00356,	0.00326,
		   0.00298,	0.00274,
		   0.00253,	0.00234,
		   0.00212,	0.00189,
		   0.00167,	0.0015,
		   0.00136,	0.00123,
		   0.00112,	0.001041,
		   9.17E-4,	7.64E-4,
		   4.92E-4,	2.45E-4,
		   9.9E-5,	1.64E-5 };
    double yem[] = { 0.031045933711196384,	0.0025678006153126453,
		     0.003808017857100988,	0.004928792955683978,
		     0.005475372133471842,	0.005444713031923721,
		     0.004928792955683978,	0.0041622710147226115,
		     0.003405994715204356,	0.0030985964564621835,
		     0.0029811742652854096,	0.002811547616527239,
		     0.0026404734423962684,	0.002492308167141455,
		     0.0023525518060183073,	0.0021654098919142305,
		     0.001967053634245899,	0.0017886866690396059,
		     0.0016107451691686056,	0.0014552319402761885,
		     0.001326989073052224,	0.001225275479229059,
		     0.001164817582284883,	0.0011168258592994702,
		     0.0010756393447619887,	0.001027861858422619,
		     9.7082439194738E-4,	9.202716990106781E-4,
		     8.825531145489205E-4,	8.509406559801923E-4,
		     8.228000972289684E-4,	8.099382692526635E-4,
		     7.912016177940992E-4,	7.725283166331187E-4,
		     7.316419889536139E-4,	6.815423684555495E-4,
		     6.407807737440318E-4,	5.950630218724736E-4,
		     5.636488268416782E-4,	5.636488268416782E-4,
		     5.322593352868506E-4,	4.919349550499538E-4,
		     4.7853944456021594E-4,	4.7423622805517506E-4,
		     4.651881339845203E-4,	4.428317965096906E-4,
		     4.248529157249601E-4,	3.9824615503479755E-4,
		     3.671511950137164E-4,	3.3615472627943223E-4,
		     3.138470965295043E-4,	3.138470965295043E-4,
		     3.176476034853718E-4,	3.3541019662496844E-4,
		     3.3105890714493693E-4,	2.996664812754339E-4,
		     2.641968962724581E-4,	2.3853720883753126E-4,
		     2.1633307652783938E-4,	1.9416487838947602E-4,
		     1.6401219466856724E-4,	1.5181897114655994E-4,
		     1.5966527487215246E-4,	1.3961733416735901E-4,
		     9.740636529508736E-5,	6.800735254367721E-5,
		     3.733630940518893E-5,	9.437160589923222E-6 };
    double yep[] = { 0.027627884464793896,	0.0038148263394288343,
		     0.004395600072800072,	0.005347317084295637,
		     0.005874325493194942,	0.0058935897380119695,
		     0.005446999173857106,	0.004729587719875803,
		     0.004031935019317648,	0.0037237615390892046,
		     0.0035759474269066098,	0.00339607125955861,
		     0.0032050585018061684,	0.003026549190084311,
		     0.0028468403537957655,	0.00266865134478073,
		     0.002469412885687608,	0.0022703523955544874,
		     0.0020813697413001857,	0.0019043371550227129,
		     0.00175524927004685,	0.001643471934655411,
		     0.0015533190271158081,	0.001475127113166862,
		     0.0014046351839534703,	0.0013267252918370102,
		     0.0012490796611905903,	0.0011691449867317568,
		     0.0011208925015361642,	0.0011472575996697516,
		     9.932774033471214E-4,	9.52102935611481E-4,
		     9.139474820797965E-4,	8.75956619930462E-4,
		     8.34865258589672E-4,	7.749193506423748E-4,
		     7.244998274671983E-4,	6.685057965343307E-4,
		     6.092618484691126E-4,	5.727128425310541E-4,
		     5.412947441089743E-4,	5.280151512977634E-4,
		     5.147815070493501E-4,	5.015974481593781E-4,
		     4.924428900898052E-4,	4.792702786528704E-4,
		     4.7010637094172637E-4,	4.428317965096906E-4,
		     4.204759208325728E-4,	3.8910152916687437E-4,
		     3.6674241641784496E-4,	3.488552708502481E-4,
		     3.2649655434629013E-4,	3.0886890422961E-4,
		     3.04138126514911E-4,	2.996664812754339E-4,
		     2.9068883707497264E-4,	2.641968962724581E-4,
		     2.418677324489565E-4,	2.1095023109728988E-4,
		     1.972308292331602E-4,	1.739482681718907E-4,
		     1.4475496537252186E-4,	1.2083045973594572E-4,
		     1.097679370308106E-4,	6.992138442565336E-5,
		     3.54682957019364E-5,	9.078546139112804E-6 };
    double ysm[] = { 0.0063,	8.0E-4,
		     3.7E-4,	4.3E-4,
		     4.1E-4,	4.0E-4,
		     4.3E-4,	4.3E-4,
		     4.2E-4,	4.2E-4,
		     4.3E-4,	4.2E-4,
		     4.0E-4,	4.0E-4,
		     3.9E-4,	3.9E-4,
		     3.8E-4,	3.7E-4,
		     3.6E-4,	3.6E-4,
		     3.5E-4,	3.3E-4,
		     3.2E-4,	3.2E-4,
		     3.1E-4,	3.1E-4,
		     3.1E-4,	3.0E-4,
		     3.0E-4,	2.9E-4,
		     2.9E-4,	2.8E-4,
		     2.8E-4,	2.8E-4,
		     2.7E-4,	2.6E-4,
		     2.5E-4,	2.5E-4,
		     2.4E-4,	2.4E-4,
		     2.3E-4,	2.2E-4,
		     2.1E-4,	2.0E-4,
		     2.0E-4,	1.9E-4,
		     1.9E-4,	1.9E-4,
		     1.8E-4,	1.7E-4,
		     1.6E-4,	1.6E-4,
		     1.5E-4,	1.5E-4,
		     1.4E-4,	1.3E-4,
		     1.3E-4,	1.3E-4,
		     1.2E-4,	1.1E-4,
		     1.0E-4,	9.3E-5,
		     7.3E-5,	5.8E-5,
		     3.2E-5,	2.0E-5,
		     1.3E-5,	4.1E-6 };
    double ysp[] = { 0.0063,	8.0E-4,
		     3.7E-4,	4.3E-4,
		     4.1E-4,	4.0E-4,
		     4.3E-4,	4.3E-4,
		     4.2E-4,	4.2E-4,
		     4.3E-4,	4.2E-4,
		     4.0E-4,	4.0E-4,
		     3.9E-4,	3.9E-4,
		     3.8E-4,	3.7E-4,
		     3.6E-4,	3.6E-4,
		     3.5E-4,	3.3E-4,
		     3.2E-4,	3.2E-4,
		     3.1E-4,	3.1E-4,
		     3.1E-4,	3.0E-4,
		     3.0E-4,	2.9E-4,
		     2.9E-4,	2.8E-4,
		     2.8E-4,	2.8E-4,
		     2.7E-4,	2.6E-4,
		     2.5E-4,	2.5E-4,
		     2.4E-4,	2.4E-4,
		     2.3E-4,	2.2E-4,
		     2.1E-4,	2.0E-4,
		     2.0E-4,	1.9E-4,
		     1.9E-4,	1.9E-4,
		     1.8E-4,	1.7E-4,
		     1.6E-4,	1.6E-4,
		     1.5E-4,	1.5E-4,
		     1.4E-4,	1.3E-4,
		     1.3E-4,	1.3E-4,
		     1.2E-4,	1.1E-4,
		     1.0E-4,	9.3E-5,
		     7.3E-5,	5.8E-5,
		     3.2E-5,	2.0E-5,
		     1.3E-5,	4.1E-6 };
    int np = 68;
    TGraphAsymmErrors* g = new TGraphAsymmErrors(np,
						 x,
						 y,
						 xem,
						 xep,
						 yem,
						 yep);
    g->SetName("/HepData/8068/d6x1y1");
    g->SetTitle("/HepData/8068/d6x1y1");
    return   g;
  }
  TGraphAsymmErrors* GetALICE(Double_t eta, UShort_t sNN)
  {
    Double_t aEta = TMath::Abs(eta);
    TGraphAsymmErrors* g = 0;
    switch (sNN) { 
    case 900:
      if      (aEta <= 0.51) g = ALICEsqrts900eta05();
      else if (aEta <= 1.01) g = ALICEsqrts900eta10();
      else if (aEta <= 1.31) g = ALICEsqrts900eta13();
      break;
    }
    if (g) { 
      g->SetMarkerColor(kALICEColor);
    }
    return g;
  }
  TGraphAsymmErrors* ALICEsqrts900eta05()
  {
    // Plot: p7742_d8x1y1
    double x[] = { 0.0,	1.0,
		   2.0,	3.0,
		   4.0,	5.0,
		   6.0,	7.0,
		   8.0,	9.0,
		   10.0,	11.0,
		   12.0,	13.0,
		   14.0,	15.0,
		   16.0,	17.0,
		   18.0,	19.0,
		   20.0,	21.5,
		   23.5,	25.5 };
    double xem[] = { 0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	1.0,
		     1.0,	1.0 };
    double xep[] = { 0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	1.0,
		     1.0,	1.0 };
    double y[] = { 0.179767,	0.155985,
		   0.140806,	0.120933,
		   0.097273,	0.075258,
		   0.057395,	0.043794,
		   0.033212,	0.024939,
		   0.018642,	0.013927,
		   0.010337,	0.007693,
		   0.005665,	0.004134,
		   0.002994,	0.002144,
		   0.001527,	0.001093,
		   8.03E-4,	4.82E-4,
		   2.56E-4,	1.02E-4 };
    double yem[] = { 0.028627623478032542,	0.006826681770816624,
		     0.004736912179891031,	0.005286866936097409,
		     0.004605827830043151,	0.003753897441326814,
		     0.003172958871463669,	0.0023132784095305087,
		     0.0017387768689512751,	0.0012885437516824954,
		     9.891637882575362E-4,	7.64476945368531E-4,
		     5.769662035162892E-4,	4.500277769204919E-4,
		     3.7107950630558943E-4,	3.1458703088334716E-4,
		     2.7004073766748605E-4,	2.210701246211256E-4,
		     1.7819090885900997E-4,	1.3576450198781714E-4,
		     1.2447088012864694E-4,	9.762171889492625E-5,
		     6.095900261651268E-5,	3.4828149534536E-5 };
    double yep[] = { 0.028627623478032542,	0.006826681770816624,
		     0.004736912179891031,	0.005286866936097409,
		     0.004605827830043151,	0.003753897441326814,
		     0.003172958871463669,	0.0023132784095305087,
		     0.0017387768689512751,	0.0012885437516824954,
		     9.891637882575362E-4,	7.64476945368531E-4,
		     5.769662035162892E-4,	4.500277769204919E-4,
		     3.7107950630558943E-4,	3.1458703088334716E-4,
		     2.7004073766748605E-4,	2.210701246211256E-4,
		     1.7819090885900997E-4,	1.3576450198781714E-4,
		     1.2447088012864694E-4,	9.762171889492625E-5,
		     6.095900261651268E-5,	3.4828149534536E-5 };
    double ysm[] = { 0.003595,	0.00312,
		     0.002816,	0.002419,
		     0.001945,	0.001505,
		     0.001148,	8.76E-4,
		     6.64E-4,	4.99E-4,
		     3.97E-4,	3.57E-4,
		     3.03E-4,	2.66E-4,
		     2.34E-4,	2.01E-4,
		     1.71E-4,	1.46E-4,
		     1.26E-4,	9.6E-5,
		     9.7E-5,	7.1E-5,
		     4.6E-5,	2.7E-5 };
    double ysp[] = { 0.003595,	0.00312,
		     0.002816,	0.002419,
		     0.001945,	0.001505,
		     0.001148,	8.76E-4,
		     6.64E-4,	4.99E-4,
		     3.97E-4,	3.57E-4,
		     3.03E-4,	2.66E-4,
		     2.34E-4,	2.01E-4,
		     1.71E-4,	1.46E-4,
		     1.26E-4,	9.6E-5,
		     9.7E-5,	7.1E-5,
		     4.6E-5,	2.7E-5 };
    int np = 24;
    TGraphAsymmErrors* g = new TGraphAsymmErrors(np,
						 x,
						 y,
						 xem,
						 xep,
						 yem,
						 yep);
    g->SetName("/HepData/7742/d8x1y1");
    g->SetTitle("/HepData/7742/d8x1y1");
    return g;
  }
  TGraphAsymmErrors* ALICEsqrts900eta10()
  {
    // Plot: p7742_d9x1y1
    double x[] = { 0.0,	1.0,
		   2.0,	3.0,
		   4.0,	5.0,
		   6.0,	7.0,
		   8.0,	9.0,
		   10.0,	11.0,
		   12.0,	13.0,
		   14.0,	15.0,
		   16.0,	17.0,
		   18.0,	19.0,
		   20.0,	21.0,
		   22.0,	23.0,
		   24.0,	25.0,
		   26.0,	27.0,
		   28.0,	29.0,
		   30.0,	31.0,
		   32.0,	33.0,
		   34.0,	35.0,
		   36.0,	37.0,
		   38.0,	39.0,
		   40.0,	41.5 };
    double xem[] = { 0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	1.0 };
    double xep[] = { 0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	1.0 };
    double y[] = { 0.092396,	0.076087,
		   0.073233,	0.078506,
		   0.080885,	0.077084,
		   0.070058,	0.062258,
		   0.054437,	0.047195,
		   0.040978,	0.035352,
		   0.03048,	0.026336,
		   0.022788,	0.019725,
		   0.017071,	0.014674,
		   0.012509,	0.010593,
		   0.008935,	0.007527,
		   0.006361,	0.005416,
		   0.004665,	0.004041,
		   0.003519,	0.00306,
		   0.002637,	0.002234,
		   0.001856,	0.001512,
		   0.001213,	9.62E-4,
		   7.61E-4,	6.02E-4,
		   4.84E-4,	3.93E-4,
		   3.27E-4,	2.76E-4,
		   2.37E-4,	1.7E-4 };
    double yem[] = { 0.021888152137629163,	0.0066076761421849355,
		     0.004451948000594796,	0.0034909631908686747,
		     0.002952765483407038,	0.0028889466938661224,
		     0.0029957812002881653,	0.0032637882897026274,
		     0.002922475833946279,	0.0026227742945209753,
		     0.0023487871338203465,	0.002091141554271255,
		     0.0017347841364273539,	0.0014762689456870655,
		     0.0012594490859101847,	0.0010506721658062519,
		     9.019456746389995E-4,	7.836938177630343E-4,
		     6.432239112470866E-4,	5.679656679765071E-4,
		     4.967293025381128E-4,	4.475187146924696E-4,
		     4.040519768544636E-4,	3.3646693745448454E-4,
		     2.9206163733020465E-4,	2.820017730440715E-4,
		     2.6982401672201087E-4,	2.5470178640912594E-4,
		     2.3519353732617736E-4,	2.1914607000811127E-4,
		     1.9727392123643713E-4,	1.7402586014727813E-4,
		     1.5060212481900778E-4,	1.3009611831257687E-4,
		     1.1269871339105873E-4,	9.693812459502196E-5,
		     8.48528137423857E-5,	7.427651041883969E-5,
		     6.726812023536856E-5,	6.030754513325841E-5,
		     5.544366510251645E-5,	4.601086828130936E-5 };
    double yep[] = { 0.021888152137629163,	0.0066076761421849355,
		     0.004451948000594796,	0.0034909631908686747,
		     0.002952765483407038,	0.0028889466938661224,
		     0.0029957812002881653,	0.0032637882897026274,
		     0.002922475833946279,	0.0026227742945209753,
		     0.0023487871338203465,	0.002091141554271255,
		     0.0017347841364273539,	0.0014762689456870655,
		     0.0012594490859101847,	0.0010506721658062519,
		     9.019456746389995E-4,	7.836938177630343E-4,
		     6.432239112470866E-4,	5.679656679765071E-4,
		     4.967293025381128E-4,	4.475187146924696E-4,
		     4.040519768544636E-4,	3.3646693745448454E-4,
		     2.9206163733020465E-4,	2.820017730440715E-4,
		     2.6982401672201087E-4,	2.5470178640912594E-4,
		     2.3519353732617736E-4,	2.1914607000811127E-4,
		     1.9727392123643713E-4,	1.7402586014727813E-4,
		     1.5060212481900778E-4,	1.3009611831257687E-4,
		     1.1269871339105873E-4,	9.693812459502196E-5,
		     8.48528137423857E-5,	7.427651041883969E-5,
		     6.726812023536856E-5,	6.030754513325841E-5,
		     5.544366510251645E-5,	4.601086828130936E-5 };
    double ysm[] = { 0.001848,	0.001522,
		     0.001465,	0.00157,
		     0.001618,	0.001542,
		     0.001401,	0.001245,
		     0.001089,	9.44E-4,
		     8.2E-4,	7.07E-4,
		     6.1E-4,	5.27E-4,
		     4.56E-4,	3.94E-4,
		     3.41E-4,	3.0E-4,
		     2.76E-4,	2.56E-4,
		     2.38E-4,	2.23E-4,
		     2.07E-4,	1.91E-4,
		     1.8E-4,	1.7E-4,
		     1.59E-4,	1.47E-4,
		     1.3E-4,	1.25E-4,
		     1.14E-4,	1.02E-4,
		     9.1E-5,	8.2E-5,
		     7.4E-5,	6.6E-5,
		     6.0E-5,	5.4E-5,
		     5.0E-5,	4.6E-5,
		     4.3E-5,	3.4E-5 };
    double ysp[] = { 0.001848,	0.001522,
		     0.001465,	0.00157,
		     0.001618,	0.001542,
		     0.001401,	0.001245,
		     0.001089,	9.44E-4,
		     8.2E-4,	7.07E-4,
		     6.1E-4,	5.27E-4,
		     4.56E-4,	3.94E-4,
		     3.41E-4,	3.0E-4,
		     2.76E-4,	2.56E-4,
		     2.38E-4,	2.23E-4,
		     2.07E-4,	1.91E-4,
		     1.8E-4,	1.7E-4,
		     1.59E-4,	1.47E-4,
		     1.3E-4,	1.25E-4,
		     1.14E-4,	1.02E-4,
		     9.1E-5,	8.2E-5,
		     7.4E-5,	6.6E-5,
		     6.0E-5,	5.4E-5,
		     5.0E-5,	4.6E-5,
		     4.3E-5,	3.4E-5 };
    int np = 42;
    TGraphAsymmErrors* g = new TGraphAsymmErrors(np,
						 x,
						 y,
						 xem,
						 xep,
						 yem,
						 yep);
    g->SetName("/HepData/7742/d9x1y1");
    g->SetTitle("/HepData/7742/d9x1y1");
    return g;
  }
  TGraphAsymmErrors* ALICEsqrts900eta13()
  {
    // Plot: p7742_d10x1y1
    double x[] = { 0.0,	1.0,
		   2.0,	3.0,
		   4.0,	5.0,
		   6.0,	7.0,
		   8.0,	9.0,
		   10.0,	11.0,
		   12.0,	13.0,
		   14.0,	15.0,
		   16.0,	17.0,
		   18.0,	19.0,
		   20.0,	21.0,
		   22.0,	23.0,
		   24.0,	25.0,
		   26.0,	27.0,
		   28.0,	29.0,
		   30.0,	31.0,
		   32.0,	33.0,
		   34.0,	35.0,
		   36.0,	37.0,
		   38.0,	39.0,
		   40.0,	41.5,
		   43.5,	45.5,
		   47.5,	49.5,
		   51.5,	53.5 };
    double xem[] = { 0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	1.0,
		     1.0,	1.0,
		     1.0,	1.0,
		     1.0,	1.0 };
    double xep[] = { 0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	0.5,
		     0.5,	1.0,
		     1.0,	1.0,
		     1.0,	1.0,
		     1.0,	1.0 };
    double y[] = { 0.068695,	0.056382,
		   0.049828,	0.056116,
		   0.063345,	0.066198,
		   0.065598,	0.061117,
		   0.055042,	0.04927,
		   0.044132,	0.039615,
		   0.035906,	0.03262,
		   0.029489,	0.026657,
		   0.023839,	0.021302,
		   0.018884,	0.016732,
		   0.014707,	0.012944,
		   0.01137,	0.009982,
		   0.008753,	0.007661,
		   0.006711,	0.005891,
		   0.005179,	0.004572,
		   0.004064,	0.003636,
		   0.003278,	0.002968,
		   0.002677,	0.002392,
		   0.002111,	0.001829,
		   0.001561,	0.001319,
		   0.001115,	7.93E-4,
		   5.24E-4,	3.51E-4,
		   2.2E-4,	1.56E-4,
		   1.27E-4,	8.5E-5 };
    double yem[] = { 0.021322725646595934,	0.007630831737104416,
		     0.004038016592338372,	0.0036354708085748677,
		     0.0035003511252444373,	0.0030936373413831173,
		     0.0026374775828431228,	0.002920906195001818,
		     0.00276560481631053,	0.002864668392676542,
		     0.0026843833556330957,	0.0022856071403458645,
		     0.002325611317481922,	0.002006865217198205,
		     0.0017203360136903488,	0.0014989629748596194,
		     0.0012394615766533467,	0.0011637718848640398,
		     9.151568171630477E-4,	9.078199160626517E-4,
		     7.98723982361867E-4,	6.967302203866285E-4,
		     6.308502199413107E-4,	5.846990678973245E-4,
		     5.4516236113657E-4,	5.460009157501477E-4,
		     4.930608481719068E-4,	4.3796118549478783E-4,
		     3.4394330928221294E-4,	3.2323366161339073E-4,
		     3.0869078379504626E-4,	2.9443165590676553E-4,
		     2.756320010448714E-4,	2.630304164920856E-4,
		     2.5628889948649746E-4,	2.4142286552851616E-4,
		     2.245551157288562E-4,	2.0678491240900531E-4,
		     1.8828170383762728E-4,	1.7255723688098393E-4,
		     1.5911316727411343E-4,	1.3023056476879765E-4,
		     1.0761505470890214E-4,	8.62670273047588E-5,
		     6.129437168288782E-5,	6.037383539249432E-5,
		     3.535533905932738E-5,      3.8078865529319545E-5 };
    double yep[] = { 0.021322725646595934,	0.007630831737104416,
		     0.004038016592338372,	0.0036354708085748677,
		     0.0035003511252444373,	0.0030936373413831173,
		     0.0026374775828431228,	0.002920906195001818,
		     0.00276560481631053,	0.002864668392676542,
		     0.0026843833556330957,	0.0022856071403458645,
		     0.002325611317481922,	0.002006865217198205,
		     0.0017203360136903488,	0.0014989629748596194,
		     0.0012394615766533467,	0.0011637718848640398,
		     9.151568171630477E-4,	9.078199160626517E-4,
		     7.98723982361867E-4,	6.967302203866285E-4,
		     6.308502199413107E-4,	5.846990678973245E-4,
		     5.4516236113657E-4,	5.460009157501477E-4,
		     4.930608481719068E-4,	4.3796118549478783E-4,
		     3.4394330928221294E-4,	3.2323366161339073E-4,
		     3.0869078379504626E-4,	2.9443165590676553E-4,
		     2.756320010448714E-4,	2.630304164920856E-4,
		     2.5628889948649746E-4,	2.4142286552851616E-4,
		     2.245551157288562E-4,	2.0678491240900531E-4,
		     1.8828170383762728E-4,	1.7255723688098393E-4,
		     1.5911316727411343E-4,	1.3023056476879765E-4,
		     1.0761505470890214E-4,	8.62670273047588E-5,
		     6.129437168288782E-5,	6.037383539249432E-5,
		     3.535533905932738E-5,      3.8078865529319545E-5 };
    double ysm[] = { 0.012402,	0.001128,
		     9.97E-4,	0.001122,
		     0.001267,	0.001324,
		     0.001312,	0.001222,
		     0.001101,	9.85E-4,
		     8.83E-4,	7.92E-4,
		     7.18E-4,	6.52E-4,
		     5.9E-4,	5.33E-4,
		     4.77E-4,	4.26E-4,
		     3.94E-4,	3.81E-4,
		     3.58E-4,	3.33E-4,
		     3.16E-4,	3.08E-4,
		     2.91E-4,	2.71E-4,
		     2.55E-4,	2.39E-4,
		     2.24E-4,	2.12E-4,
		     2.07E-4,	1.99E-4,
		     1.82E-4,	1.72E-4,
		     1.72E-4,	1.62E-4,
		     1.49E-4,	1.38E-4,
		     1.27E-4,	1.2E-4,
		     1.14E-4,	9.6E-5,
		     8.5E-5,	7.1E-5,
		     5.1E-5,	5.4E-5,
		     2.5E-5,    3.3E-5 };
    double ysp[] = { 0.012402,	0.001128,
		     9.97E-4,	0.001122,
		     0.001267,	0.001324,
		     0.001312,	0.001222,
		     0.001101,	9.85E-4,
		     8.83E-4,	7.92E-4,
		     7.18E-4,	6.52E-4,
		     5.9E-4,	5.33E-4,
		     4.77E-4,	4.26E-4,
		     3.94E-4,	3.81E-4,
		     3.58E-4,	3.33E-4,
		     3.16E-4,	3.08E-4,
		     2.91E-4,	2.71E-4,
		     2.55E-4,	2.39E-4,
		     2.24E-4,	2.12E-4,
		     2.07E-4,	1.99E-4,
		     1.82E-4,	1.72E-4,
		     1.72E-4,	1.62E-4,
		     1.49E-4,	1.38E-4,
		     1.27E-4,	1.2E-4,
		     1.14E-4,	9.6E-5,
		     8.5E-5,	7.1E-5,
		     5.1E-5,	5.4E-5,
		     2.5E-5,    3.3E-5 };
    int np = 48;
    TGraphAsymmErrors* g = new TGraphAsymmErrors(np,
						 x,
						 y,
						 xem,
						 xep,
						 yem,
						 yep);
    g->SetName("/HepData/7742/d10x1y1");
    g->SetTitle("/HepData/7742/d10x1y1");
    return g;
  }
  ClassDef(Unfolder,1);
};
void UnfoldMultDists(const char* method="Bayes", 
		     Double_t    regParam=4,
		     const char* realFile="forward_multdists.root",
                     const char* mcFile="forward_mcmultdists.root")
{
  Unfolder* u = new Unfolder;
  u->Run(method, regParam, realFile, mcFile);
  
  u->DrawAll();
}
#endif
