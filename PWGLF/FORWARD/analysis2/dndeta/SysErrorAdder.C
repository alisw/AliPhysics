#ifndef __CINT__
# include <TString.h>
# include <TH1.h>
# include "GraphSysErr.C"
# include <TLegend.h>
# include <TLegendEntry.h>
# include <TObjString.h>
#else
class GraphSysErr;
class TH1;
class TLegend;
#endif

struct SysErrorAdder
{
  enum {
    kTriggerFill = 1001, // 3144
    kMergingFill = 3001, // 3001,
    kDensityFill = 3002, // 3002,
    kEmpiricalFill = 3144, // 3244
    kHadronFill = 3244
  };
  TString fSys;
  UShort_t fSNN;
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
  virtual ~SysErrorAdder() {}
  virtual const char* GetTriggerName() const { return "Trigger"; }
  void ModError(GraphSysErr* gse, Int_t id, Style_t fill, TLegend* l) const
  {
    gse->SetSysFillColor(id,gse->GetMarkerColor());
    gse->SetSysLineColor(id,gse->GetMarkerColor());
    gse->SetSysLineWidth(id, 1);
    gse->SetSysLineStyle(id, 1);
    gse->SetSysFillStyle(id, fill);
    gse->SetSysOption(id, GraphSysErr::kBox);
    
    AddLegend(l, id, gse->GetSysTitle(id), fill);
  }
  void AddLegend(TLegend* l, Int_t id, const char* title, Style_t fill) const
  {
    if (!l) return;
    TLegendEntry* e = l->AddEntry(Form("sys%02d", id), title, "f");
    e->SetFillStyle(fill);
    e->SetFillColor(kBlack);
    e->SetLineColor(kBlack);
    e->SetLineWidth(1);
  }
  /** 
   * Declare the systematic error from the trigger 
   * 
   * @param gse Graph 
   * 
   * @return Id of systmatic error or -1
   */
  Int_t  MakeTrigger(GraphSysErr* gse, TLegend* l) const
  {
    Double_t low = 0,high = 0;
    GetTrigger(low,high);
    if (low == 0 && high == 0) return -1;
    
    Int_t id = gse->DefineCommon(GetTriggerName(), true, low, high);
    ModError(gse, id, kTriggerFill, l);
    return id;
  }
  /** 
   * Declare the systematic error from the merging 
   * 
   * @param gse Graph 
   * 
   * @return Id of systmatic error or -1
   */
  Int_t MakeMerging(GraphSysErr* gse, TLegend* l) const
  {
    if (GetMerging() <= 0) return -1;
    Int_t id = gse->DeclarePoint2Point("Merging", true); 
    ModError(gse, id, kMergingFill, l);
    return id;
  }
  /** 
   * Declare the systematic error from the denisty calculation
   * 
   * @param gse Graph 
   * 
   * @return Id of systmatic error or -1
   */
  Int_t MakeDensity(GraphSysErr* gse, TLegend* l) const
  {
    if (GetDensity() <= 0) return -1;
    Int_t id = gse->DeclarePoint2Point("Density", true);
    ModError(gse, id, kDensityFill, l);
    return id;
  }
  /** 
   * Declare the systematic error from the empirical correction 
   * 
   * @param gse Graph 
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
  virtual Int_t MakeHadron(GraphSysErr* gse, TLegend* l) const
  {
    if (GetHadron(4) <= 0) return -1;
    Int_t id = gse->DeclarePoint2Point("Hadron chemistry", true); 
    ModError(gse, id, kHadronFill, l);
    return id;
  }
  /** 
   * @return The systematic error from merging 
   */
  virtual Double_t GetMerging() const { return 0.01; }
  /** 
   * @return The systematic error from density calculator 
   */
  virtual Double_t GetDensity() const { return 0.01; }
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
   * @param h histogram 
   * 
   * @return Graph 
   */  
  virtual GraphSysErr* Make(TH1* h, TLegend* l)
  {
    GraphSysErr* gse = new GraphSysErr(10);
    gse->SetSumLineColor(kRed+2);
    gse->SetSumLineWidth(2);
    gse->SetSumTitle("All errors");
    gse->SetSumOption(GraphSysErr::kHat);
    gse->SetLineColor(h->GetLineColor());
    gse->SetMarkerColor(h->GetMarkerColor());
    gse->SetFillColor(h->GetFillColor());
    gse->SetMarkerStyle(h->GetMarkerStyle());
    gse->SetFillStyle(h->GetFillStyle());
    gse->SetXTitle("#it{#eta}");
    gse->SetYTitle("1/#it{N} d#it{N}_{ch}/d#it{#eta}");
		      

    MakeTrigger(gse, l);
    Int_t idMerge = MakeMerging(gse, l); 
    Int_t idDens  = MakeDensity(gse, l); 
    Int_t idEmp   = MakeEmpirical(gse, l);
    Int_t idHad   = MakeHadron(gse, l);

    Int_t cnt = 0;
    for (Int_t i = 1; i <= h->GetNbinsX(); i++) {
      Double_t y  = h->GetBinContent(i);
      if (y < 1e-6) continue;
      Double_t ey = h->GetBinError(i);
      Double_t ex = h->GetXaxis()->GetBinWidth(i)/2;
      Double_t x  = h->GetXaxis()->GetBinCenter(i);

      gse->SetPoint(cnt, x, y);
      gse->SetPointError(cnt, ex);
      gse->SetStatError(cnt, ey);

      if (idMerge>=0) gse->SetSysError(idMerge, cnt, ex, GetMerging());
      if (idDens >=0) gse->SetSysError(idDens,  cnt, ex, GetDensity());
      if (idEmp  >=0) gse->SetSysError(idEmp,   cnt, ex, GetEmpirical());
      if (idHad  >=0) gse->SetSysError(idHad,   cnt, ex, GetHadron(x));
      cnt++;
    }
    return gse;
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
   */
  INELAdder(const TString& sys, UShort_t sNN)
    : SysErrorAdder(sys, sNN, "INEL"), fLow(0), fHigh(0)
  {
    if (fSys.EqualTo("pp", TString::kIgnoreCase)) {
      switch (fSNN) {
      case  900:   fLow = 0.001;  fHigh = 0.003; break;
      case 2760:   fLow = 0.0035; fHigh = 0.006; break;
      case 7000: 
      case 8000:   fLow = 0.003;  fHigh = 0.006; break;
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
};
/**
 * For pp INEL results
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
    : SysErrorAdder(sys, sNN, "INEL>0")
  {}
  /** 
   * Get trigger systematic error 
   * 
   * @param low   On return, the low error 
   * @param high  On return, the high error
   */
  virtual void GetTrigger(Double_t& low, Double_t& high) const
  {
    low = high = 0;
  }
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
   * @param sys Collision system 
   * @param sNN Collision energy
   */
  NSDAdder(const TString& sys, UShort_t sNN)
    : SysErrorAdder(sys, sNN, "NSD"), fValue(0)
  {
    if (fSys.EqualTo("pp", TString::kIgnoreCase)) {
      switch (fSNN) {
      case  900:   fValue = 0.02; break;
      case 2760:   fValue = 0.03; break;
      case 7000: 
      case 8000:   fValue = 0.02; break;
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
    low = high = fValue;
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
  Double_t fMin;
  Double_t fMax;
  /** 
   * Constructor 
   * 
   * @param sys Collision system 
   * @param sNN Collision energy
   * @param Method 
   */
  CENTAdder(const TString& sys, UShort_t sNN, const TString& method)
    : SysErrorAdder(sys, sNN, method), fCent(0), fValue(0)
  {
    if (fSys.EqualTo("pPb", TString::kIgnoreCase)) {
      fMin = 0.02; fMax = 0.04;
    }
    else if (fSys.EqualTo("Pbp", TString::kIgnoreCase)) {
      fMin = 0.04; fMax = 0.06;
    }
    else {
      fMin = 0.01; fMax = 0.02;
    }
  }
  const char* GetTriggerName() const 
  {
    static TString n;
    n = fTrig;
    n.ReplaceAll("CENT", "Centrality (");
    n.Append(")");
    return n.Data();
  }
  /** 
   * Get trigger systematic error 
   * 
   * @param low   On return, the low error 
   * @param high  On return, the high error
   */
  virtual void GetTrigger(Double_t& low, Double_t& high) const
  {
    low = high = ((fCent-2.5)/100) * (fMax-fMin) + fMin;
    // Printf("Trigger error for centrality = %f -> %f", fCent, low);
  }
  /** 
   * Get centrality 
   * 
   * @param h Histogram 
   * 
   * @return Centrality 
   */
  Double_t GetCentrality(TH1* h) 
  { 
    TString name(h->GetName());
    Int_t idx = name.Index("_cent");
    name.Remove(0,idx+5);
    name.ReplaceAll("d", ".");
    name.ReplaceAll("_", " ");
    TObjArray* tokens  = name.Tokenize(" ");
    TObjString* first  = static_cast<TObjString*>(tokens->At(0));
    TObjString* second = static_cast<TObjString*>(tokens->At(1));
    Double_t    low    = first->String().Atof();
    Double_t    high   = second->String().Atof();
    fCent              = (low+high)/2;
    return fCent;
  }
  /** 
   * Declare the systematic error from the empirical correction 
   * 
   * @param gse Graph 
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
   * 
   * @return Graph 
   */  
  virtual GraphSysErr* Make(TH1* h, TLegend* l)
  {
    fCent = GetCentrality(h);
    return SysErrorAdder::Make(h, l);
  }
};
