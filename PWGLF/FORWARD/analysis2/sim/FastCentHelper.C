#ifndef FASTCENTHELPER_H
#define FASTCENTHELPER_H
#ifndef __CINT__
# include <TAxis.h>
# include <TString.h>
# include <TList.h>
# include <TH1.h>
# include <TH2.h>
# include <TError.h>
# include <THStack.h>
# include <TClass.h>
# include <TStyle.h>
#else
class TAxis;
class TString;
class TList;
class TH1;
class TH2;
class TCollection;
class TClass;
class TH1D;
#endif

struct FastCentHelper
{
  TAxis*       fCentAxis;
  TString      fCentMeth;
  TList*       fCentList;
  TH1*         fCentAll; //!
  TH1*         fCentAcc; //!
  TH2*         fCentMult; //!
  TH2*         fCentNPart; //!
  TH2*         fCentNBin; //!
  TH2*         fCentB; //!
  FastCentHelper(const char* method)
    : fCentAxis(0),
      fCentMeth(""),
      fCentList(0),
      fCentAll(0),
      fCentAcc(0),
      fCentMult(0),
      fCentNPart(0),
      fCentNBin(0),
      fCentB(0)
  {
    TString    axis("default");
    TString    meth(method);
    TObjArray* tokens = meth.Tokenize(":");

    if (fCentMeth.IsNull()) {
      SetCentAxis(axis);
      return;
    }
    fCentMeth = tokens->At(0)->GetName();
    if (tokens->GetEntriesFast() > 1)
      axis = tokens->At(1)->GetName();
    
    if (fCentMeth.Contains("RefMult")) SetCentAxis("mult");
    else                               SetCentAxis(axis);
  }
  /** 
   * Set the centrality axis (equi-distant)
   * 
   * @param n     Number of bin
   * @param low   Lowest bound 
   * @param high  Highest bound 
   */
  void SetCentAxis(Int_t n, Double_t low, Double_t high)
  {
    if (fCentAxis) {
      delete fCentAxis;
      fCentAxis = 0;
    }
    fCentAxis = new TAxis(n, low, high);
  }
  /** 
   * Set the centrality axis.
   * 
   * @param n       Number of bins 
   * @param edges   Bin edges 
   */
  void SetCentAxis(Int_t n, Double_t* edges)
  {
    if (fCentAxis) {
      delete fCentAxis;
      fCentAxis = 0;
    }
    fCentAxis = new TAxis(n, edges);
  }
  /** 
   * Set the centrality axis.  Spec is either a pre-defined string, or
   * a list of colon separated bin edges. Pre-defined settings are 
   *
   * - pbpb, aa, default: For PbPb
   * - ppb, pbp, pa, ap: For pPb/Pbp 
   * - pp: for pp centrality 
   * 
   * Example of specs could be 
   *
   * - 0:5:10:20:30:40:50:60:80:100
   * - 0:0.1:1:5:10:20:40:70:100 
   * 
   * @param spec 
   */
  void SetCentAxis(const char* spec)
  {
    Info("SetCentAxis", "Setting centrality axis from %s", spec);
    TString    s(spec);
    s.ToLower();
    if (s.IsNull()) return;
    if (s.EqualTo("pbpb") || s.EqualTo("aa") || s.EqualTo("default")) {
      Printf("Setting centrality axis Pb-Pb");
      Double_t aa[] = { 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };
      SetCentAxis(11, aa);
      fCentAxis->SetUniqueID(2);
      return;
    }
    if (s.EqualTo("ppb") || s.EqualTo("pbp") ||
	s.EqualTo("pa") || s.EqualTo("ap")) {
      Printf("Setting centrality axis p-Pb/Pb-p");
      Double_t pa[] = { 0, 5, 10, 20, 40, 60, 80, 100 };
      SetCentAxis(7, pa);
      fCentAxis->SetUniqueID(3);
      return;
    }
    if (s.EqualTo("pp")) {
      Printf("Setting centrality axis pp");
      Double_t pp[] = { 0, 0.01, 0.1, 1, 5, 10, 15, 20, 30, 40, 50, 70, 100 };
      SetCentAxis(12, pp);
      fCentAxis->SetUniqueID(1);
      return;
    }
    
    TObjArray* tokens  = s.Tokenize(":");
    Int_t      nTokens = tokens->GetEntriesFast();
    TArrayD    edges(nTokens);
    for (Int_t i = 0; i < nTokens; i++) {
      TObjString* token = static_cast<TObjString*>(tokens->At(i));
      TString&    edge  = token->String();
      edges[i]          = edge.Atof();
    }
    SetCentAxis(edges.GetSize()-1, edges.GetArray());
    delete tokens;
  }
  /** 
   * Get the color associated with a centrality bin. 
   * 
   * @param low  Low edge 
   * @param high High edge 
   * 
   * @return Color identifier. 
   */
  Int_t GetCentralityColor(Double_t, Double_t high) const
  {
    Float_t  fc       = high / 100;
    Int_t    nCol     = gStyle->GetNumberOfColors();
    Int_t    icol     = TMath::Min(nCol-1,int(fc * nCol + .5));
    Int_t    col      = gStyle->GetColorPalette(icol);
    //Info("GetCentralityColor","%3d: %3d-%3d -> %3d",bin,centLow,centHigh,col);
    return col;
  }
  /** 
   * Get the histogram name 
   * 
   * @param low  Low edge 
   * @param high High edge
   * 
   * @return String
   */
  virtual const char* HistName(Double_t low, Double_t high) const
  {
    return Form("h%03dd%02d_%03dd%02d",
		Int_t(low),  Int_t(100*low) %100,
		Int_t(high), Int_t(100*high)%100);
  }
  /** 
   * Get the histogram title
   * 
   * @param low  Low edge 
   * @param high High edge
   * 
   * @return String 
   */
  virtual const char* HistTitle(Double_t low, Double_t high) const
  {
    if (fCentMeth.BeginsWith("RefMult")) {
      Int_t iLow  = low;
      Int_t iHigh = high;
      if (iLow == iHigh) {
	if (iLow == 0) return " 0";
	else           return Form("%2d+", iLow);
      }
      return Form("%2d-%2d", iLow+1, iHigh);
    }      
    return Form("%6.2f-%6.2f%%", low, high);
  }
  /** 
   * Modify a histogram.  Sets the line, marker, and fill styles and
   * colors, as well as the name and title. 
   * 
   * @param h    Histogram to modify 
   * @param low  Low edge of our bin
   * @param up   High edge of our bin 
   */
  void ModHist(TH1* h, Double_t low, Double_t high)
  {
    Color_t col = GetCentralityColor(low, high);
    h->SetLineColor(col);
    h->SetLineStyle(1);
    h->SetMarkerColor(col);
    h->SetMarkerStyle(24);
    h->SetFillColor(kWhite);
    h->SetFillStyle(0);
    h->SetName(HistName(low, high));
    h->SetTitle(HistTitle(low, high));
  }
  /** 
   * Create histograms 
   * 
   * @param output 
   */
  void CreateHistos(TCollection* output, TH1D* (*callback)())
  {
    // Create our stack 
    // fCentStack = new THStack("stack", "Stack of all dN_{ch}/d#eta");
    fCentList  = new TList;
    fCentList->SetName("byCent");
    for (Int_t i = 1; i <= fCentAxis->GetNbins(); i++) {
      Double_t low  = fCentAxis->GetBinLowEdge(i);
      Double_t high = fCentAxis->GetBinUpEdge(i);
      TH1*     hist = (*callback)();
      ModHist(hist, low, high);
      hist->SetMinimum(0);
      hist->Sumw2();
      fCentList->Add(hist);
    }
    if (fCentMeth.BeginsWith("RefMult")) {
      // Create null-bin 
      TH1* first = (*callback)();
      ModHist(first, 0, 0);
      fCentList->AddFirst(first);

      // Create overflow-bin
      TH1* last = (*callback)();
      ModHist(last, 100, 100);
      fCentList->AddLast(last);
    }      
    output->Add(fCentList);

    if (fCentAxis->GetXbins()->GetArray()) 
      fCentAll = new TH1D("cent", "All Centralities",
			   fCentAxis->GetNbins(),
			   fCentAxis->GetXbins()->GetArray());
    else
      fCentAll = new TH1D("cent", "All Centralities",
			   fCentAxis->GetNbins(),
			   fCentAxis->GetXmin(),
			   fCentAxis->GetXmax());
    fCentAll->SetXTitle("Centrality [%]");
    fCentAll->SetYTitle("Events");
    fCentAll->SetFillColor(kRed+2);
    fCentAll->SetFillStyle(3002);
    fCentAll->SetMarkerStyle(20);
    fCentAll->SetMarkerColor(kRed+2);
    fCentAll->SetLineColor(kRed+2);
    fCentAll->SetDirectory(0);
    fCentAll->SetMinimum(0);
    output->Add(fCentAll);

    fCentAcc = static_cast<TH1*>(fCentAll->Clone("centAcc"));
    fCentAcc->SetTitle("Accepted centralities");
    fCentAcc->SetFillColor(kGreen+2);
    fCentAcc->SetMarkerColor(kGreen+2);
    fCentAcc->SetLineColor(kGreen+2);
    fCentAcc->SetDirectory(0);
    output->Add(fCentAcc);
    // output->ls();
  }
  /** 
   * Create diagnostics histograms 
   * 
   * @param output 
   */  
  void CreateDiagnostics(TCollection* output, TH1* centHist)
  {
    if (fCentMult || fCentNPart) return;

    fCentMult  = 0;
    fCentNPart = 0;
    fCentNBin  = 0;
    Int_t maxNPart = 2*210;
    Int_t maxNBin  = 3*210;
    if (fCentAxis->GetXbins()->GetArray()) {
      fCentNPart = new TH2D("centNPart", "Centrality vs. N_{part}",
			    fCentAxis->GetNbins(), 
			    fCentAxis->GetXbins()->GetArray(),
			    maxNPart, 0, maxNPart);
      fCentNBin  = new TH2D("centNBin", "Centrality vs. N_{bin}", 
			    fCentAxis->GetNbins(), 
			    fCentAxis->GetXbins()->GetArray(),
			    maxNBin, 0, maxNBin);
      fCentB     = new TH2D("centB", "Centrality vs. b",
			    fCentAxis->GetNbins(), 
			    fCentAxis->GetXbins()->GetArray(),
			    200, 0, 20);
    } else {
      fCentNPart = new TH2D("centNPart", "Centrality vs. N_{part}",
			    fCentAxis->GetNbins(), fCentAxis->GetXmin(), 
			    fCentAxis->GetXmax(),
			    maxNPart, 0, maxNPart);
      fCentNBin  = new TH2D("centNBin", "Centrality vs. N_{bin}", 
			    fCentAxis->GetNbins(), fCentAxis->GetXmin(), 
			    fCentAxis->GetXmax(),
			    maxNBin, 0, maxNBin);
      fCentB     = new TH2D("centB", "Centrality vs. b",
			    fCentAxis->GetNbins(), fCentAxis->GetXmin(), 
			    fCentAxis->GetXmax(),
			    200, 0, 20);
    }

    fCentNPart->SetXTitle(Form("Centrality (%s) %%", fCentMeth.Data()));
    fCentNPart->SetYTitle("N_{part}");
    fCentNPart->SetDirectory(0);

    fCentNBin->SetXTitle(Form("Centrality (%s) %%", fCentMeth.Data()));
    fCentNBin->SetYTitle("N_{bin}");
    fCentNBin->SetDirectory(0);

    fCentB->SetXTitle(Form("Centrality (%s) %%", fCentMeth.Data()));
    fCentB->SetYTitle("b [fm]");
    fCentB->SetDirectory(0);
    
    output->Add(fCentNPart);
    output->Add(fCentNBin);
    output->Add(fCentB);

    if (!centHist) return;;
    if (fCentAxis->GetXbins()->GetArray()) {
      fCentMult  = new TH2D("centMult","Event multiplicity vs. centrality",
			    centHist->GetXaxis()->GetNbins(),
			    centHist->GetXaxis()->GetXmin(),
			    centHist->GetXaxis()->GetXmax(),
			    fCentAxis->GetNbins(),
			    fCentAxis->GetXbins()->GetArray());
    } else {
      fCentMult  = new TH2D("centMult","Event multiplicity vs. centrality",
			    centHist->GetXaxis()->GetNbins(),
			    centHist->GetXaxis()->GetXmin(),
			    centHist->GetXaxis()->GetXmax(),
			    fCentAxis->GetNbins(),
			    fCentAxis->GetXmin(), fCentAxis->GetXmax());
    }
    fCentMult->SetXTitle(Form("Event multiplicity (%s)", fCentMeth.Data()));
    fCentMult->SetYTitle(Form("Centrality (%s) %%", fCentMeth.Data()));
    fCentMult->SetDirectory(0);    
    output->Add(fCentMult);

    return;
  }    
  /** 
   * Check the centrality 
   * 
   * @param cent  Centrality 
   * @param mult  Estimator value 
   * @param b     Impact parameter 
   * @param nPart Number of participants 
   * @param nBin  Number of binary collisions 
   * 
   * @return Centrality bin, or -1 in case of problems
   */
  Int_t CheckCentrality(Double_t cent,
			Double_t mult,
			Double_t b,
			Double_t nPart,
			Double_t nBin)
  {
    Int_t ret = -1;
    fCentAll->Fill(cent);
    if (fCentMult) fCentMult->Fill(mult, cent);
    fCentB->Fill(cent, b);
    if (cent < 0 || cent > 999) {
      Warning("FillDiagnostics",
	      "Centrality is unreasonable: %f -> %f",mult, cent);
      return ret;
    }
    Int_t n = fCentAxis->GetNbins();
    ret = fCentAxis->FindBin(cent);
    if (ret-1 == n && cent == fCentAxis->GetXmax())
      ret = n;
    if (ret < 1 || ret > n) {
      Warning("ProcessHeader", "Centrality %f -> %f -> bin # %d",
	      mult, cent, ret);
      ret = -1;
      return ret;
    }
    fCentNPart->Fill(cent, nPart);
    fCentNBin->Fill(cent, nBin);
    if (ret == n) cent -= 0.001;
    fCentAcc->Fill(cent);
    return ret;
  }
  /** 
   * Get an object from the output list, possibly checking the type 
   * 
   * @param name Name of object 
   * @param cls  Possible class pointer 
   * 
   * @return Found object, or null
   */
  TObject* GetOutputObject(TCollection* output, const char* name, TClass* cls)
  {
    TObject* o = output->FindObject(name);
    if (!o) {
      Warning("GetOutputObject", "Object %s not found in output", name);
      output->ls();
      return 0;
    }
    if (cls && !o->IsA()->InheritsFrom(cls)) {
      Warning("GetOutputObject", "Object %s from output is not a %s, but a %s",
	      o->GetName(), cls->GetName(), o->ClassName());
      return o;
    }
    return o;
  }
  Bool_t Finalize(TCollection* output, Long_t minEvents)
  {
    
    fCentAll   = static_cast<TH1*>(GetOutputObject(output, "cent",
						   TH1::Class()));
    fCentAcc   = static_cast<TH1*>(GetOutputObject(output, "centAcc",
						   TH1::Class()));
    fCentMult  = static_cast<TH2*>(GetOutputObject(output, "centMult",
						   TH2::Class()));
    fCentNPart = static_cast<TH2*>(GetOutputObject(output, "centNPart",
						   TH2::Class()));
    fCentNBin  = static_cast<TH2*>(GetOutputObject(output, "centNBin",
						   TH2::Class()));
    fCentList = static_cast<TList*>(GetOutputObject(output, "byCent",
						    TList::Class()));
    if (!fCentList || !fCentAll || !fCentAcc) {
      Warning("Finalize", "Missing stack and histograms");
      return false;
    }
    Info("Terminate", "Accepted %d/%d=%6.2f%% events",
	 int(fCentAcc->GetEntries()), int(fCentAll->GetEntries()),
	 100*fCentAcc->GetEntries()/fCentAll->GetEntries());

    THStack*  stack = new THStack("all", "All");
    TList*    hists = fCentList; // ->GetHists();
    TObjLink* link  = hists->FirstLink();
    Int_t     bin   = 1;
    Long64_t  sum   = 0;
    Long64_t  total = 0;
    Long64_t  all   = fCentAll->GetEntries();
    while (link) {
      TObject* o = link->GetObject();
      if (!o) {
	link = link->Next();
	bin++;
	continue;
      }
      TH1*  h = static_cast<TH1*>(o);
      Int_t n = h->GetBinContent(0);
      Int_t m = fCentAcc->GetBinContent(bin);
      total += m;
      h->SetBinContent(0,0);
      printf("%9d (%9d) events in bin %s ...", n, m, o->GetTitle());
      if (n < minEvents) {
	// Too few event, remove this
	TObjLink* tmp = link->Next();
	hists->Remove(link);
	link = tmp;
	delete o;
	Printf(" removed");
	bin++;
	continue;
      }
      sum += m;
      // Scale
      h->Scale(1. / n, "width");
      stack->Add(h);
      Printf(" scaled");
      link = link->Next();
      bin++;
    }
    Printf("ana/acc/all: %9lld/%9lld/%9lld [%6.2f%%/%6.2f%%]",
	   sum, total, all, float(100*sum)/total, float(100*total)/all);
    output->Add(stack);
  }
  ClassDef(FastCentHelper,1);
};
#endif

//
// EOF
//
