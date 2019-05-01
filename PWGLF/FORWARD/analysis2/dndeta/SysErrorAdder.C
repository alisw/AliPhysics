#ifndef __CINT__
# include <TString.h>
# include <TH1.h>
# include "GraphSysErr.C"
# include <TLegend.h>
# include <TLegendEntry.h>
# include <TObjString.h>
# include <TColor.h>
#else
class GraphSysErr;
class TH1;
class TLegend;
#endif

struct SysErrorAdder
{
  /** 
   * Fill styles 
   */
  enum {
    kTriggerFill   = 1001, // 0, // 3144
    kMergingFill   = 3001, // 0, // 3001,
    kDensityFill   = 3002, // 0, // 3002,
    kEmpiricalFill = 3144, // 0, // 3244
    kHadronFill    = 3244, // 0,
    kEMFill        = 3344,
    kDiffRatFill   = 3001,
    kDiffMassFill  = 3002,
    kOtherFill     = 0
  };
  /** System */
  TString fSys;
  /** Energy */
  UShort_t fSNN;
  /** Trigger */
  TString fTrig;
  /** 
   * Constructor 
   * 
   * @param sys   System 
   * @param sNN   Collision energy 
   * @param trig  Trigger 
   */
  SysErrorAdder(const TString& sys,
		UShort_t       sNN,
		const TString& trig)
    : fSys(sys), fSNN(sNN), fTrig(trig)
  {}
  /** 
   * Destructor 
   */
  virtual ~SysErrorAdder() {}
  /** 
   * Trigger name
   * 
   * @return The trigger name 
   */
  virtual const char* GetTriggerName() const { return "Trigger"; }
  virtual const char* GetTriggerString() const { return fTrig; }
  Int_t ModColor(Color_t n, UShort_t off) const
  {
    ULong_t max   = 0xff;
    ULong_t pixel = TColor::Number2Pixel(n);
    ULong_t r     = TMath::Min((((pixel >> 16) & max) + off), max);
    ULong_t g     = TMath::Min((((pixel >>  8) & max) + off), max);
    ULong_t b     = TMath::Min((((pixel >>  0) & max) + off), max);
    ULong_t next  = (r << 16) | (g << 8) | (b << 0);
    // Printf("0x%08x -> 0x%08x", pixel, next);
    return TColor::GetColor(next);
  }
    
  /** 
   * Modify an error 
   * 
   * @param gse   Graph 
   * @param id    Identifier 
   * @param fill  Fill style 
   * @param l     Legend 
   * @param off   Off-set 
   */
  void ModError(GraphSysErr* gse, Int_t id, Style_t fill,
		TLegend* l, UShort_t off=64) const
  {
    // UShort_t off = 64; // (id-1) * 32;
    // switch (fill) {
    // case kMergingFill:    off = 16; break; 
    // case kDensityFill:    off = 32; break;
    // case kEmpiricalFill:  off = 48; break;
    // case kHadronFill:     off = 64; break;
    // }
    Color_t c = gse->GetMarkerColor(); // ModColor(gse->GetMarkerColor(), off);
	
    gse->SetSysFillColor(id, c);
    gse->SetSysLineColor(id, c);
    gse->SetSysLineWidth(id, 1);
    gse->SetSysLineStyle(id, id);
    gse->SetSysFillStyle(id, fill);
    // gse->SetSysOption(id, GraphSysErr::kBox);
    gse->SetSysOption(id, GraphSysErr::kRect);
    
    AddLegend(l, id, gse->GetSysTitle(id), fill);
  }
  /** 
   * Add legend entry 
   * 
   * @param l       Legend
   * @param id      Identifier 
   * @param title   Title 
   * @param fill    Fill style 
   */
  void AddLegend(TLegend* l, Int_t id, const char* title, Style_t fill) const
  {
    if (!l) return;
    TLegendEntry* e = l->AddEntry(Form("sys%02d", id), title, "f");
    e->SetFillStyle(fill);
    e->SetFillColor(kBlack);
    e->SetLineColor(kBlack);
    e->SetLineStyle(id);
    e->SetLineWidth(1);
  }
  /** 
   * Declare the systematic error from the trigger 
   * 
   * @param gse Graph 
   * @param l    Legend
   * 
   * @return Id of systmatic error or -1
   */
  virtual Int_t  MakeTrigger(GraphSysErr* gse, TLegend* l) const
  {
    Double_t low = 0,high = 0;
    GetTrigger(low,high);
    if (low == 0 && high == 0) return -1;
    
    Int_t id = gse->DefineCommon(GetTriggerName(), true, low, high);
    ModError(gse, id, kTriggerFill, l, 0);
    return id;
  }
  /** 
   * Declare the systematic error from the merging 
   * 
   * @param gse Graph 
   * @param l    Legend
   * 
   * @return Id of systmatic error or -1
   */
  Int_t MakeMerging(GraphSysErr* gse, TLegend* l) const
  {
    if (GetMerging(-1) <= 0 && GetMerging(+1) <= 0) return -1;
    Int_t id = gse->DeclarePoint2Point("Merging", true); 
    ModError(gse, id, kMergingFill, l);
    return id;
  }
  /** 
   * Declare the systematic error from the denisty calculation
   * 
   * @param gse Graph 
   * @param l    Legend
   * 
   * @return Id of systmatic error or -1
   */
  Int_t MakeDensity(GraphSysErr* gse, TLegend* l) const
  {
    if (GetDensity(-1) <= 0 && GetDensity(+1) <= 0) return -1;
    Int_t id = gse->DeclarePoint2Point("Density", true);
    ModError(gse, id, kDensityFill, l);
    return id;
  }
  /** 
   * Declare the systematic error from the empirical correction 
   * 
   * @param gse Graph 
   * @param l    Legend
   * 
   * @return Id of systmatic error or -1
   */
  virtual Int_t MakeEmpirical(GraphSysErr* gse, TLegend* l) const
  {
    if (GetEmpirical() <= 0) return -1;
    Int_t id = gse->DeclarePoint2Point("Empirical", true); 
    ModError(gse, id, kEmpiricalFill, l);
    return id;
  }
  /** 
   *  Declare the systematic error from the hadron chemistry
   * 
   * @param gse  Graph 
   * @param l    Legend
   * 
   * @return Identifier 
   */
  virtual Int_t MakeHadron(GraphSysErr* gse, TLegend* l) const
  {
    if (GetHadron(4) <= 0) return -1;
    Int_t id = gse->DeclarePoint2Point("Hadron chemistry", true); 
    ModError(gse, id, kHadronFill, l);
    return id;
  }
  virtual const char* MakeName(Double_t eff) const
  {
    return GetTriggerString();
  } 
  /** 
   * @return The systematic error from merging 
   */
  virtual Double_t GetMerging(Short_t ) const { return 0.01; }
  /** 
   * @return The systematic error from density calculator 
   */
  virtual Double_t GetDensity(Short_t w) const { return (w<0?0.02:0.01); }
  /** 
   * @return The systematic error from empirical correction 
   */
  virtual Double_t GetEmpirical() const { return 0.061; }
  /** 
   * @return The systematic error from empirical correction 
   */
  virtual Double_t GetHadron(Double_t eta) const
  {
    Double_t ret = (eta >= 3.6 ? 0.02 : 0);
    /// Printf("Eta=%f -> %f", eta, ret);
    return ret;
  }
  /** 
   * Get trigger systematic error 
   * 
   * @param low   On return, the low error 
   * @param high  On return, the high error
   */
  virtual void GetTrigger(Double_t& low, Double_t& high) const = 0;
  /** 
   * Create a graph 
   * 
   * @param h    histogram 
   * @param l    Option Legend to add to 
   * @param eff  Efficiency 
   * @param verb Be verbose 
   * 
   * @return Graph 
   */  
  virtual GraphSysErr* Make(TH1* h, TLegend* l, Double_t eff=1,
			    Bool_t verb=false)
  {
    // --- Reaction key ------------------------------------------------
    if (verb) Info("", "Making graph from %s %p eff=%f", h->GetName(), l, eff);
    TString reac = fSys;
    reac.ToUpper();
    reac.Insert(reac.Index("P", 1, 1, TString::kIgnoreCase), " ");
    reac.Append(" --> CHARGED X");
    
    GraphSysErr* gse = new GraphSysErr(10);
    gse->SetName(MakeName(eff));
    gse->SetSumLineColor(kRed+2);
    gse->SetSumLineWidth(2);
    gse->SetSumTitle("All errors");
    gse->SetSumOption(GraphSysErr::kHat);
    gse->SetLineColor(h->GetLineColor());
    gse->SetMarkerColor(h->GetMarkerColor());
    gse->SetFillColor(h->GetFillColor());
    gse->SetMarkerStyle(h->GetMarkerStyle());
    gse->SetFillStyle(h->GetFillStyle());
    gse->SetXTitle("\\mathit{\\eta}");
    gse->SetYTitle("\\mathrm{d}N_{ch}/\\mathrm{d}\\eta");
    gse->SetKey("laboratory", "CERN");
    gse->SetKey("accelerator", "LHC");
    gse->SetKey("detector", "FORWARD");
    gse->SetKey("reackey", reac);
    gse->SetKey("obskey", "DN/DETARAP");
    gse->SetKey("title", "Systematic study of dNch/deta over widest "
		"possible eta region at the LHC");
    gse->SetKey("author", "CHRISTENSEN");
    gse->SetKey("abstract", "We present dNch/deta over widest "
		"possible eta region at the LHC");
    gse->SetKey("dscomment", "The pseudo-rapidity density of charged particle");
    gse->AddQualifier(Form("SQRT(S)%s IN GEV",
			   fSys.EqualTo("pp",TString::kIgnoreCase) ?
			   "" : "/NUCLEON"), Form("%d", fSNN));
    
    MakeTrigger(gse, l);
    Int_t idMerge = MakeMerging(gse, l); 
    Int_t idDens  = MakeDensity(gse, l); 
    Int_t idEmp   = MakeEmpirical(gse, l);
    Int_t idHad   = MakeHadron(gse, l);

    Int_t cnt = 0;
    if (verb) 
      Info("", "Looping over histogram w/%d bins", h->GetNbinsX());
    for (Int_t i = 1; i <= h->GetNbinsX(); i++) {
      Double_t y  = h->GetBinContent(i);
      if (y < 1e-6) continue;
      Double_t ey = h->GetBinError(i);
      Double_t ex = h->GetXaxis()->GetBinWidth(i)/2;
      Double_t x  = h->GetXaxis()->GetBinCenter(i);

      gse->SetPoint(cnt, x, eff*y);
      gse->SetPointError(cnt, ex);
      gse->SetStatError(cnt, eff*ey);

      if (idMerge>=0)
	gse->SetSysError(idMerge, cnt, ex, ex, GetMerging(-1), GetMerging(1));
      if (idDens >=0)
	gse->SetSysError(idDens,  cnt, ex, ex, GetDensity(-1), GetDensity(1));
      if (idEmp  >=0)
	gse->SetSysError(idEmp,   cnt, ex, GetEmpirical());
      if (idHad  >=0)
	gse->SetSysError(idHad,   cnt, ex, GetHadron(x));
      cnt++;
    }
    return gse;
  }
  /** 
   * Create a systematic uncertainty adder 
   * 
   * @param t Trigger type 
   * @param s System type 
   * @param e Collision energy 
   * @param c Centrality method
   * 
   * @return Newly created object 
   */
  static SysErrorAdder* Create(const TString& t,
			       const TString& s,
			       UShort_t       e,
			       const TString& c);
};
/**
 * For pp INEL results
 * 
 */
struct OfflineAdder : public SysErrorAdder
{
  /** 
   * Constructor 
   * 
   * @param sys Collision system 
   * @param sNN Collision energy
   */
  OfflineAdder(const TString& sys, UShort_t sNN)
    : SysErrorAdder(sys, sNN, "OFFLINE")
  {
  }
  /** 
   * Get trigger systematic error 
   * 
   * @param low   On return, the low error 
   * @param high  On return, the high error
   */
  virtual void GetTrigger(Double_t& low, Double_t& high) const
  {
    low  = 0;
    high = 0;
  }
  Int_t  MakeTrigger(GraphSysErr* gse, TLegend* l) const
  {
    return -1;
  }
};

struct BareAdder : public SysErrorAdder
{
  BareAdder(const TString& sys, UShort_t sNN, TString& trig)
    : SysErrorAdder(sys,sNN,trig)
  {}
  /** 
   * Get trigger systematic error 
   * 
   * @param low   On return, the low error 
   * @param high  On return, the high error
   */
  virtual void GetTrigger(Double_t& low, Double_t& high) const
  {
    low  = high = 0;
  }
  Int_t  MakeTrigger(GraphSysErr* gse, TLegend* l) const
  {
    gse->AddQualifier("TRIGGER", fTrig);
    return SysErrorAdder::MakeTrigger(gse, l);
  }
};
  
/**
 * For pp INEL results
 * 
 */
struct INELAdder : public SysErrorAdder
{
  Double_t fLow;
  Double_t fHigh;
  /** 
   * Constructor 
   * 
   * @param sys Collision system 
   * @param sNN Collision energy
   * @param low  Low trigger uncertainty 
   * @param high High trigger uncertainty 
   */
  INELAdder(const TString& sys, UShort_t sNN, Double_t low=-1, Double_t high=-1)
    : SysErrorAdder(sys, sNN, "INEL"), fLow(low), fHigh(high)
  {
    if (low > 0 && high > 0) return;
    if (fSys.EqualTo("pp", TString::kIgnoreCase)) {
      switch (fSNN) {
      case  900:   fLow = 0.001;  fHigh = 0.003; break;
      case 2760:   fLow = 0.0035; fHigh = 0.006; break;
      case 5023:   fLow = fHigh = 0.028;         break; 
      case 7000: 
      case 8000:   fLow = 0.003;  fHigh = 0.006; break;
      case 13000:  fLow = fHigh = 0.028;         break;
      default: break;
      }
    }
  }
  /** 
   * Get trigger systematic error 
   * 
   * @param low   On return, the low error 
   * @param high  On return, the high error
   */
  virtual void GetTrigger(Double_t& low, Double_t& high) const
  {
    low  = fLow;
    high = fHigh;
  }
  Int_t  MakeTrigger(GraphSysErr* gse, TLegend* l) const
  {
    gse->AddQualifier("TRIGGER", "INEL");
    Int_t id = SysErrorAdder::MakeTrigger(gse, l);

    if (fSNN == 5023) {
      Int_t rId = gse->DefineCommon("Diff. ratio",true,
				    0.023,0.034,GraphSysErr::kBox);
      Int_t mId = gse->DefineCommon("Diff. mass",true,
				    0,0.065,GraphSysErr::kBox);
      ModError(gse, rId, kDiffRatFill, l);
      ModError(gse, mId, kDiffMassFill, l);
    }
    return id;
  }
};
/**
 * For pp INEL>0 results
 * 
 */
struct INELGt0Adder : public SysErrorAdder
{
  /** 
   * Constructor 
   * 
   * @param sys Collision system 
   * @param sNN Collision energy
   */
  INELGt0Adder(const TString& sys, UShort_t sNN)
    : SysErrorAdder(sys, sNN, "INEL>0"), fLow(0), fHigh(0)
  {
    switch (sNN) {
    case  5023:  fLow = fHigh = 0.023; break;
    case 13000:  fLow = fHigh = 0.023; break;
    }    
  }
  const char* MakeName(Double_t) const { return "INELGt0"; }
  /** 
   * Get trigger systematic error 
   * 
   * @param low   On return, the low error 
   * @param high  On return, the high error
   */
  virtual void GetTrigger(Double_t& low, Double_t& high) const
  {
    low = fLow; high = fHigh;
  }
  Int_t  MakeTrigger(GraphSysErr* gse, TLegend* l) const
  {
    gse->AddQualifier("TRIGGER", "INEL>0");
    gse->AddQualifier("N", ">0");
    Int_t id = SysErrorAdder::MakeTrigger(gse, l);
    if (fSNN == 5023) {
      Int_t rId = gse->DefineCommon("Diff. ratio",true,
				    0.005,0.005,GraphSysErr::kBox);
      Int_t mId = gse->DefineCommon("Diff. mass",true,
				    0.005,0,GraphSysErr::kBox);
      ModError(gse, rId, kDiffRatFill, l);
      ModError(gse, mId, kDiffMassFill, l);
    }
    return id;
  }
  Double_t fLow;
  Double_t fHigh;
};
/**
 * For pp NSD results
 * 
 */
struct NSDAdder : public SysErrorAdder
{
  Double_t fValue;
  /** 
   * Constructor 
   * 
   * @param sys   Collision system 
   * @param sNN   Collision energy
   * @param value Value of uncertainty 
   */
  NSDAdder(const TString& sys, UShort_t sNN, Double_t value=-1)
    : SysErrorAdder(sys, sNN, "NSD"), fValue(value)
  {
    if (value > 0) return;
    if (fSys.EqualTo("pp", TString::kIgnoreCase)) {
      switch (fSNN) {
      case  900:   fValue = 0.02;  break;
      case 2760:   fValue = 0.03;  break;
      case 5023:   fValue = 0.044; break;
      case 7000: 
      case 8000:
      case 13000:  fValue = 0.02; break;
      default: break;
      }
    }
  }
  const char* MakeName(Double_t e) const
  {
    // if (e < 1e-3 || TMath::Abs(e - 1) < 1e-3) return "V0AND";
    return "NSD";
  }
  /** 
   * Get trigger systematic error 
   * 
   * @param low   On return, the low error 
   * @param high  On return, the high error
   */
  virtual void GetTrigger(Double_t& low, Double_t& high) const
  {
    low = high = fValue;
  }
  Int_t  MakeTrigger(GraphSysErr* gse, TLegend* l) const
  {
    if (TString(gse->GetName()).EqualTo("V0AND")) {
      gse->AddQualifier("TRIGGER", "V0AND");
      return -1;
    }
    gse->AddQualifier("TRIGGER", "NSD");
    Int_t id = SysErrorAdder::MakeTrigger(gse, l);
    if (fSNN == 5023) {
      Int_t rId = gse->DefineCommon("Diff. ratio",true,
				    0.008,0.014,GraphSysErr::kBox);
      Int_t mId = gse->DefineCommon("Diff. mass",true,
				    0,0.043,GraphSysErr::kBox);
      ModError(gse, rId, kDiffRatFill, l);
      ModError(gse, mId, kDiffMassFill, l);
    }
    return id;
  }
  
};
/**
 * For centrality 
 * 
 */
struct CENTAdder : public SysErrorAdder
{
  Double_t fCent;
  Double_t fValue;
  // Double_t fMin;
  // Double_t fMax;
  TH1*     fLookup;
  Double_t fCMin;
  Double_t fCMax;
  Double_t fEM1;
  Double_t fEM2;
  /** 
   * Constructor 
   * 
   * @param sys Collision system 
   * @param sNN Collision energy
   * @param method MEthod to use  
   */
  CENTAdder(const TString& sys, UShort_t sNN, const TString& method)
    : SysErrorAdder(sys, sNN, method),
      fCent(0),
      fValue(0),
      fCMin(0),
      fCMax(0),
      fEM1(0),
      fEM2(0),
      fLookup(0)
  {
    Double_t off = .1;
    if (fSys.EqualTo("pPb", TString::kIgnoreCase) ||
	fSys.EqualTo("Pbp", TString::kIgnoreCase)) {
      TString   m;
      if (!method.BeginsWith("CENT")) m = "CENT";
      m.Append(method);
      Double_t  cent[] = { 0, 5, 10, 20, 40, 60, 80, 100+off };
      Double_t  zna[]  = { 1.5, 1.5, 1.5, 1.5, 2., 2., 3. };
      Double_t  v0a[]  = { 1.5, 1.5, 1.5, 1.5, 2., 2., 2. };
      Double_t  v0m[]  = { 1.5, 1.5, 1.5, 1.5, 2., 2., 4. };
      Double_t  cl1[]  = { 2.5, 2.5, 2.5, 4.0, 6., 10., 11. };
      Double_t* est    = 0;
      if      (m.Contains("CENTZNA") || m.Contains("CENTZNC")) est = zna;
      else if (m.Contains("CENTV0A") || m.Contains("CENTV0C")) est = v0a;
      else if (m.Contains("CENTV0M"))                          est = v0m;
      else if (m.Contains("CENTCL1"))                          est = cl1;
      if (!est) {
	Warning("","Couldn't find estimator from %s", method.Data());
	return;
      }
      fTrig   = m;
      fLookup = new TH1D("lookup", "Centrality error lookup", 7, cent);
      for (Int_t i = 1; i <= 7; i++) fLookup->SetBinContent(i,est[i-1]/100);
    }
    else {
      fTrig = "CENT";
      fLookup = new TH1D("lookup", "Centrality error lookup", 100, 0, 100);
      Double_t min = 0.004, max = 0.062, top = 100; // 0.02
      if (sNN == 5023) { 
	min = 0.005;
	max = 0.095;
	top = 90;
	fEM1 = 0.04;
	fEM2 = 0.04;
	// max = (7.5-min)/TMath::Power(80,2) * TMath::Power(100,2) + min;
      }
      else if (sNN == 5440 || sNN == 5100) { // Preliminary values for Xe-Xe
	min = 0.00086;
	max = 0.042;
	top = 90;
	fEM1 = 0.014;
	fEM2 = 0.040;
      }
      for (Int_t i = 1; i <= 100; i++) {
	Double_t c = fLookup->GetXaxis()->GetBinCenter(i);
	Double_t e = max * TMath::Power(c/top,2)+min;
	fLookup->SetBinContent(i, e);
      }
    }
    if (fLookup) fLookup->SetDirectory(0);
    
  }
  const char* GetTriggerName() const 
  {
    // static TString n;
    // n = fTrig;
    // n.ReplaceAll("CENT", "Centrality (");
    // n.Append(")");
    return "Centrality";
  }
  virtual const char* GetTriggerString() const
  {
    return Form("%s_%03dd%02d_%03dd%02d",
		fTrig.Data(),
		Int_t(fCMin), Int_t(fCMin*100)%100,
		Int_t(fCMax), Int_t(fCMax*100)%100);
  }
  /** 
   * Get trigger systematic error 
   * 
   * @param low   On return, the low error 
   * @param high  On return, the high error
   */
  virtual void GetTrigger(Double_t& low, Double_t& high) const
  {
    if (!fLookup) { low = high = 0; return; }
    Int_t bin = fLookup->FindBin(fCent);
    low = high = fLookup->GetBinContent(bin);
  }
  Int_t  MakeTrigger(GraphSysErr* gse, TLegend* l) const
  {
    gse->AddQualifier("TRIGGER", fTrig);
    gse->AddQualifier("CENTRALITY IN PCT", Form("%6.2f TO %6.2f",fCMin,fCMax));

    Int_t ret = SysErrorAdder::MakeTrigger(gse, l);

    Double_t em = fCMax > 80 ? fEM2 : fCMax > 70 ? fEM1 : 0;
    Int_t    ei = gse->DefineCommon("EM contamination", true, em, em);
    ModError(gse, ei, kEMFill, l, 0);
    
    return ret;
  }
  /** 
   * Get centrality 
   * 
   * @param h Histogram 
   * @param verb Be verbose 
   * 
   * @return Centrality 
   */
  Double_t GetCentrality(TH1* h, Bool_t verb=false) 
  {
    TString name(h->GetName());
    if (verb) 
      Info("", "Extracting centrality from %s", name.Data());
    Int_t idx = name.Index("_cent");
    if (idx == kNPOS) {
      Warning("GetCentrality", "Don't know how to parse %s",
	      name.Data());
      return -1;
    }
    name.Remove(0,idx+5);
    name.ReplaceAll("d", ".");
    name.ReplaceAll("_", " ");
    TObjArray* tokens  = name.Tokenize(" ");
    TObjString* first  = static_cast<TObjString*>(tokens->At(0));
    TObjString* second = static_cast<TObjString*>(tokens->At(1));
    fCMin              = first->String().Atof();
    fCMax              = second->String().Atof();
    fCent              = (fCMin+fCMax)/2;
    return fCent;
  }
  /** 
   * Declare the systematic error from the empirical correction 
   * 
   * @param gse Graph 
   * @param l   Optional legend 
   * 
   * @return Id of systmatic error or -1
   */
  virtual Int_t MakeEmpirical(GraphSysErr* gse, TLegend* l) const
  {
    if (GetEmpirical() <= 0) return -1;
    if (fCent >= 0)
      return SysErrorAdder::MakeEmpirical(gse, l);

    AddLegend(l, 0, "Empirical", kEmpiricalFill);
    return -1;
  }
  
  /** 
   * Create a graph 
   * 
   * @param h histogram 
   * @param l 
   * @param eff 
   * @param verb 
   * 
   * @return Graph 
   */  
  virtual GraphSysErr* Make(TH1*     h,
			    TLegend* l,
			    Double_t eff=1,
			    Bool_t   verb=false)
  {
    fCent = GetCentrality(h, verb);
    if (fCent < 0) return 0;
    return SysErrorAdder::Make(h, l, eff, verb);
  }
};


SysErrorAdder*
SysErrorAdder::Create(const TString& t,
		      const TString& s,
		      UShort_t       e,
		      const TString& c)
{
  TString tt(t);
  tt.ToUpper();

  SysErrorAdder* a = 0;
  if      (tt.EqualTo("OFFLINE"))  a=new OfflineAdder(s,e);
  else if (tt.EqualTo("UNKNOWN"))  a=new OfflineAdder(s,e);
  else if (tt.EqualTo("INEL"))     a=new INELAdder(s,e);
  else if (tt.EqualTo("INEL>0"))   a=new INELGt0Adder(s,e);
  else if (tt.EqualTo("INELGT0"))  a=new INELGt0Adder(s,e);
  else if (tt.EqualTo("NSD"))      a=new NSDAdder(s,e);
  else if (tt.EqualTo("MBOR"))     a=new BareAdder(s,e,tt);
  else if (tt.EqualTo("V0AND"))    a=new BareAdder(s,e,tt);
  else if (tt.EqualTo("VISX"))     a=new BareAdder(s,e,tt);
  else                             a=new CENTAdder(s,e,c);
  Info("Create", "Created %s adder for %s/%s/%hu/%s: %p",
       a->GetTriggerString(), t.Data(), s.Data(), e, c.Data(), a);
  return a;
}


//
// EOF
//
