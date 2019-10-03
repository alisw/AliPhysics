static Int_t PbPbBin(Double_t c1, Double_t c2)
{
  Double_t c = (c1+c2) / 2;
  if      (c <  5) return 0;
  else if (c < 10) return 1;
  else if (c < 20) return 2;
  else if (c < 30) return 3;
  else if (c < 40) return 4;
  else if (c < 50) return 5;
  else if (c < 60) return 6;
  else if (c < 70) return 7;
  else if (c < 80) return 8;
  else if (c < 90) return 9;
  return                  10;
}
/** 
 * Get the centrality color for PbPb 
 * 
 * @param c1 Lower edge
 * @param c2 Upper edge
 * 
 * @return Color 
 */
static Color_t PbPbColor(Double_t c1, Double_t c2)
{
  const Color_t cc[] = { kMagenta+2,
			 kBlue+2,
			 kAzure-1, // 10,
			 kCyan+2,
			 kGreen+1,
			 kSpring+5,//+10,
			 kYellow+1,
			 kOrange+5,//+10,
			 kRed+1,
			 kPink+5,//+10,
			 kBlack };
  Int_t bin = PbPbBin(c1,c2);
  return cc[bin];
}
//____________________________________________________________________
TObject*
GetObject(TDirectory* d, const char* name, TClass* cls=0, Bool_t verb=true)
{
  if (!d) {
    Error("GetOject", "No directory");
    return 0;
  }
  
  TObject* o = d->Get(name);
  if (!o) {
    if (verb) 
      Error("GetObject", "No object %s in directory %s", name, d->GetName());
    return 0;
  }

  if (!cls) return o;

  if (!o->IsA()->InheritsFrom(cls)) {
    if (verb) 
      Error("GetObject", "%s from %s is not a %s, but a %s",
	    o->GetName(), d->GetName(), cls->GetName(), o->ClassName());
    return 0;
  }

  return o;
}

//____________________________________________________________________
TObject*
GetObject(TCollection* d, const char* name, TClass* cls=0, Bool_t verb=true)
{
  if (!d) {
    Error("GetOject", "No collection");
    return 0;
  }
  
  TObject* o = d->FindObject(name);
  if (!o) {
    if (verb)
      Error("GetObject", "No object %s in collection %s", name, d->GetName());
    return 0;
  }

  if (!cls) return o;

  if (!o->IsA()->InheritsFrom(cls)) {
    if (verb)
      Error("GetObject", "%s from %s is not a %s, but a %s",
	    o->GetName(), d->GetName(), cls->GetName(), o->ClassName());
    return 0;
  }

  return o;
}
//____________________________________________________________________
TCollection*
GetCollection(TDirectory* d, const char* name, Bool_t verb=true)
{
  return static_cast<TCollection*>(GetObject(d, name, TCollection::Class(), verb));
}
//____________________________________________________________________
TCollection*
GetCollection(TCollection* d, const char* name, Bool_t verb=true)
{
  return static_cast<TCollection*>(GetObject(d, name, TCollection::Class(), verb));
}
//____________________________________________________________________
TH2*
GetH2(TCollection* c, const char* name, Bool_t verb=true)
{
  return static_cast<TH2*>(GetObject(c, name, TH2::Class(), verb));
}
//____________________________________________________________________
TH1*
GetH1(TCollection* c, const char* name, Bool_t verb=true)
{
  return static_cast<TH1*>(GetObject(c, name, TH1::Class(), verb));
}
//____________________________________________________________________
THStack*
GetStack(TCollection* c, const char* name, Bool_t verb=true)
{
  return static_cast<THStack*>(GetObject(c, name, THStack::Class(),verb));
}
//____________________________________________________________________
TAxis*
GetAxis(TCollection* c, const char* name, Bool_t verb=true)
{
  return static_cast<TAxis*>(GetObject(c, name, TAxis::Class(),verb));
}
//____________________________________________________________________
TParameter<double>*
GetParam(TCollection* c, const char* name, Bool_t verb=false)
{
  return static_cast<TParameter<double>*>(GetObject(c, name,
						    TParameter<double>::Class(),
						    verb));
}

//====================================================================
void
AddPath(const TString& dir, Bool_t before=true)
{
  const char* now = gROOT->GetMacroPath();
  const char* fst = (before ? dir.Data() : now);
  const char* lst = (before ? now        : dir.Data());
  gROOT->SetMacroPath(Form("%s:%s",fst,lst));

  gSystem->AddIncludePath(Form("-I%s",dir.Data()));
}

//____________________________________________________________________
void
AddRelPath(const char* rel, Bool_t first=true)
{
  TString alp(Form("$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/%s",rel));
  AddPath(alp, first);
  if (gSystem->Getenv("ANA_SRC")) {
    TString mine(Form("$ANA_SRC/%s", rel));
    AddPath(mine, first);
  }
}
  


//====================================================================
TH1* Rebin(TH1* h, Int_t rebin, Bool_t cutEdges) const
{ 
  if (rebin <= 1) return h;

  Int_t nBins = h->GetNbinsX();
  if(nBins % rebin != 0) {
    Warning("Rebin", "Rebin factor %d is not a devisor of current number "
	    "of bins %d in the histogram %s", rebin, nBins, h->GetName());
    return;
  }
    
  // Make a copy 
  TH1* tmp = static_cast<TH1*>(h->Clone("tmp"));
  tmp->Rebin(rebin);
  tmp->SetDirectory(0);
  tmp->Reset();
  // The new number of bins 
  Int_t nBinsNew = nBins / rebin;
  for(Int_t i = 1;i<= nBinsNew; i++) {
    Double_t content = 0;
    Double_t sumw    = 0;
    Double_t wsum    = 0;
    Int_t    nbins   = 0;
    for(Int_t j = 1; j<=rebin;j++) {
      Int_t    bin = (i-1)*rebin + j;
      Double_t c   =  h->GetBinContent(bin);

      if (c <= 0) continue;

      if ((cutEdges)) {
	if (h->GetBinContent(bin+1)<=0 || 
	    h->GetBinContent(bin-1)<=0) {
	  Warning("Rebin", "removing bin %d=%f of %s (%d=%f,%d=%f)", 
		  bin, c, h->GetName(), 
		  bin+1, h->GetBinContent(bin+1), 
		  bin-1, h->GetBinContent(bin-1));
	  continue;
	}	
      }
      Double_t e =  h->GetBinError(bin);
      Double_t w =  1 / (e*e); // 1/c/c
      content    += c;
      sumw       += w;
      wsum       += w * c;
      nbins++;
    }
      
    if(content > 0 && nbins > 1 ) {
      tmp->SetBinContent(i, wsum / sumw);
      tmp->SetBinError(i,1./TMath::Sqrt(sumw));
    }
  }

  // Finally, rebin the histogram, and set new content
  return tmp;
}



//____________________________________________________________________
TObject*
DoOne(TCollection*   c,
      const TString& dir,
      const TString& name,
      Int_t          rebin,
      Double_t       eff, 
      void*          oa,
      Bool_t         cutEdges=false,
      TLegend*       l=0)
{
  SysErrorAdder* adder = reinterpret_cast<SysErrorAdder*>(oa);

  TCollection* sub = GetCollection(c, dir);
  if (!sub) return 0;

  TH1* hist = GetH1(sub, name);
  if (!hist) return 0;

  hist = Rebin(hist, rebin, cutEdges);
  hist->SetMarkerColor(kBlack);
  hist->SetFillColor(kBlack);
  hist->SetLineColor(kBlack);

  if (dir.BeginsWith("cent"))
    hist->SetName(Form("%s_%s",hist->GetName(),dir.Data()));

  Info("DoOne","Make graph from %s with eff=%f", hist->GetName(), eff);
  GraphSysErr* gse = adder->Make(hist, l, eff, true);
  if (!gse) return 0;

  gse->SetSumOption(GraphSysErr::kBox);
  gse->SetSumFillColor(gse->GetMarkerColor());
  gse->SetSumFillStyle(3002);
  gse->SetCommonSumOption(GraphSysErr::kBox);
  gse->SetCommonSumFillColor(gse->GetMarkerColor());
  gse->SetCommonSumFillStyle(3001);

  return gse;
}

//____________________________________________________________________
void ModOne(TObject* o, TString& sys, UInt_t sNN, TString trg, Color_t col)
{
  GraphSysErr* g = static_cast<GraphSysErr*>(o);
  if (!sys.EqualTo("pp"))
    g->AddQualifier("SQRT(S)/NUCLEON IN GEV", Form("%d.0", sNN), true);
  else {
    g->AddQualifier("SQRT(S) IN GEV", Form("%d.0", sNN), true);
    g->AddQualifier("TRIGGER", trg);
  }

  g->SetMarkerColor(col);
  g->SetSumFillColor(kWhite);
  g->SetSumFillStyle(0);
  g->SetSumLineColor(col);
  g->SetCommonSumFillColor(col);
  g->SetCommonSumFillStyle(1001);
  g->SetCommonSumLineColor(col);
  g->Print("sys qual key");
}

//____________________________________________________________________
TH1*
DrawOne(TObject* o, Option_t* opt, Double_t& min, Double_t& max)
{
  GraphSysErr* g = static_cast<GraphSysErr*>(o);
  g->Draw(opt);
  Double_t mn, mx;
  g->GetMinMax(opt, mn, mx);
  min = TMath::Min(min, mn);
  max = TMath::Max(max, mx);
  TMultiGraph* mg = g->GetMulti();
  // Info("DrawOne", "Got  Multigraph %p", mg);
  TH1* hist = (mg ? mg->GetHistogram() : 0);
  // Info("DrawOne", "Got frame histogram %p", hist);
  return hist;
}
  

//____________________________________________________________________
TList*
ExtractGSEs(const char* filename="forward_dndeta.root",
	    Int_t       rebin=5,
	    Double_t    eff=1,
	    Bool_t      raw=false,
	    Bool_t      cutEdges=false,
	    const char* name="Forward")
{
  if (!gROOT->GetClass("GraphSysErr")) {
    AddPath("$HOME/GraphSysErr");
    AddRelPath("gse",false);

    gROOT->LoadMacro("GraphSysErr.C+g");
  }
  if (!gROOT->GetClass("SysErrorAdder")) {
    AddRelPath("dndeta", true);
    gROOT->LoadMacro("SysErrorAdder.C+g");
  }

  TFile* file = TFile::Open(filename,"READ");
  if (!file) {
    Warning("ExtractGSEs", "Failed to open %s", filename);
    return 0;
  }

  TCollection* results = GetCollection(file, Form("%sdNdetaResults",name));
  if (!results) return 0;


  // TH1*     empH   = GetH1    (results, "empirical");
  TObject* empO   = GetObject(results, "empirical");
  Bool_t   emp    = (empO != 0 && !raw);
  TObject* osNN   = GetObject(results, "sNN");
  TObject* osys   = GetObject(results, "sys");
  TObject* otrg   = GetObject(results, "trigger");
  TObject* ocentM = GetObject(results, "centEstimator");
  TAxis*   centA  = GetAxis  (results, "centAxis");
  TString  sys    = osys->GetTitle();
  UShort_t sNN    = osNN->GetUniqueID();
  TString  mth    = (ocentM && centA ? ocentM->GetTitle() : "");
  if (mth.EqualTo("none",TString::kIgnoreCase) ||
      (centA && centA->GetNbins() < 1)) mth = "";
  TString  trg    = (!mth.IsNull() && centA ? "CENT" : otrg->GetTitle());
  TString  hname  = Form("dndeta%s%s",name,emp?"Emp":"");

  Printf("Read settings\n"
	 "  Empirical:   %s (%p)\n"
	 "  System:      %s (%p)\n"
	 "  Energy:      %d (%s - %p)\n"
	 "  Trigger:     %s (%s - %p)\n"
	 "  Efficiency:  %6.4f\n"
	 "  Centrality:  %s (%p)\n"
	 "    Axis:      %p",
	 (emp  ? "yes" : "no"), empO,
	 sys.Data(), osys,
	 sNN, (osNN ? osNN->GetTitle() : ""), osNN,
	 trg.Data(), eff, otrg->GetTitle(), otrg, mth.Data(),
	 ocentM, centA);
  if (centA) Printf("    %d bins between %f and %f",
		    centA->GetNbins(), centA->GetXmin(), centA->GetXmax());

  // Possibly change trigger 
  if (trg.EqualTo("MBOR")  && eff > 0 && eff < 1) trg="INEL";
  if (trg.EqualTo("V0AND") && eff > 0 && eff < 1) trg="NSD";
  if (trg.EqualTo("VISX")  && eff > 0 && eff < 1) trg="NSD";
  TString strg(trg);
  if (strg.Contains("INEL>0")) trg = "INELGt0";
  
  
  SysErrorAdder* adder = SysErrorAdder::Create(trg,sys,sNN,mth);
  if (!adder) {
    Warning("ExtractGSEs", "Failed to make adder");
    return 0;
  }

  TCanvas* c = new TCanvas("c", "C");
  c->SetTopMargin(0.01);
  c->SetRightMargin(0.20);
  TLegend* l = new TLegend(0.8,0.1,0.99,0.99,
			   Form("%s @ %s (%s)",
				sys.Data(), osNN->GetTitle(), trg.Data()));
  l->SetFillColor(0);
  l->SetFillStyle(0);
  l->SetBorderSize(0);
  
  
  TList* ret   = new TList;
  Bool_t first = true;
  TH1*   frame = 0;
  Double_t min = 100000, max = -1000000;
  if (!centA || centA->GetNbins() < 1 || mth.IsNull()) {
    Info("ExtractGSEs", "Doing pp-like extraction"
	 " all bin, %s, rebin=%d, eff=%f, %p, %s",
	 hname.Data(), rebin, eff, adder,
	 (cutEdges ? "cut edges" : "edges stay"));
    TObject* all = DoOne(results,"all",hname,rebin,eff,adder,cutEdges,l);
    if (!all) {
      Warning("ExtractGSEs", "Nothing returned from DoOne(\"all\"...)");
      return 0;
    }
    GraphSysErr* gg = static_cast<GraphSysErr*>(all); 
    ModOne(gg, sys, sNN, trg, gg->GetMarkerColor());
    gg->SetName(strg);
    ret->Add(all);
    frame = DrawOne(all, "SUM QUAD AXIS", min, max);
  }
  else {
    for (Int_t i = 1; i <= centA->GetNbins(); i++) {
      Double_t low  = centA->GetBinLowEdge(i);
      Double_t high = centA->GetBinUpEdge(i);
      TString  dir;
      dir.Form("cent%03dd%02d_%03dd%02d",
	       Int_t(low),  Int_t(low *100)%100,
	       Int_t(high), Int_t(high*100)%100);
      TObject* g = DoOne(results, dir, hname, rebin, eff, adder, cutEdges,
			 (first ? l : 0));
      if (!g) continue;
      ret->Add(g);

      GraphSysErr* gse = static_cast<GraphSysErr*>(g);
      Color_t col = PbPbColor(low, high);
      ModOne(gse, sys, sNN, trg, col);
      if (first) frame = DrawOne(g, "SUM QUAD AXIS", min, max);
      else  	         DrawOne(g, "SUM QUAD", min, max);
      TLegendEntry* e = l->AddEntry("", Form("%3d-%3d%%",
					     int(low), int(high)),"F");
      e->SetFillColor(col);
      e->SetFillStyle(1001);
      e->SetLineColor(col);
      e->SetMarkerColor(col);
      first = false;
    }
  }
  if (!frame) {
    Warning("ExtractGSEs", "No frame given");
  }
  else {
    frame->SetMinimum(0.9*min);
    frame->SetMaximum(1.1*max);
  }
  l->Draw();
  
  TString outName;
  if (mth.EqualTo("V0M", TString::kIgnoreCase) &&
      sys.EqualTo("PbPb",TString::kIgnoreCase))
    mth = "";
  if (raw) mth = "RAW";
  outName.Form("%s_%05d_%s%s.root",sys.Data(),sNN,strg.Data(),mth.Data());
  TFile* out = TFile::Open(outName,"RECREATE");
  ret->Write("container",TObject::kSingleKey);
  out->Write();

  Info("ExtractGSEs", "Written to %s", outName.Data());
  return ret;
}
  
