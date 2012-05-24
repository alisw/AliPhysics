#include "AliForwardFlowWeights.h"
#include <TGraph.h>
#include <TF1.h>
#include <TList.h>
#include <TMath.h>

//____________________________________________________________________
AliForwardFlowWeights::AliForwardFlowWeights()
  : fV22Pt(0), 
    fV24Pt(0), 
    fV24AltPt(0),
    fV2B(0), 
    fV2C(0)
{}

//____________________________________________________________________
AliForwardFlowWeights::AliForwardFlowWeights(const AliForwardFlowWeights& o)
  : TObject(o),
    fV22Pt(o.fV22Pt),
    fV24Pt(o.fV24Pt),
    fV24AltPt(o.fV24AltPt),
    fV2B(o.fV2B), 
    fV2C(o.fV2C)
{}

//____________________________________________________________________
AliForwardFlowWeights&
AliForwardFlowWeights::operator=(const AliForwardFlowWeights& o)
{
  if (&o == this) return *this;

  fV22Pt    = (o.fV22Pt ?    static_cast<TGraph*>(o.fV22Pt->Clone())    : 0);
  fV24Pt    = (o.fV24Pt ?    static_cast<TGraph*>(o.fV24Pt->Clone())    : 0);
  fV24AltPt = (o.fV24AltPt ? static_cast<TGraph*>(o.fV24AltPt->Clone()) : 0);
  fV2B      = (o.fV2B   ?    static_cast<TGraph*>(o.fV2B->Clone())      : 0);
  fV2C      = (o.fV2C   ?    static_cast<TGraph*>(o.fV2C->Clone())      : 0);

  return *this;
}
//____________________________________________________________________
AliForwardFlowWeights::~AliForwardFlowWeights()
{
}

namespace {
  const char* fgkPt2Name = "v22VsPt";
  const char* fgkPt4Name = "v24VsPt";
  const char* fgkPt4Alt  = "v24AltVsPt";
  const char* fgkBName   = "v2VsB";
  const char* fgkCName   = "v2VsC";
}

//____________________________________________________________________
void
AliForwardFlowWeights::Init(TList* l)
{
  Int_t          ptN     = 19;
  const Double_t ptX[]   = {0.00,     0.25,     0.350,    0.45, 
			    0.55,     0.650,    0.75,     0.85, 
			    0.950,    1.10,     1.30,     1.500,
			    1.70,     1.90,     2.250,    2.75, 
			    3.25,     3.750,    4.50};
  { 
    // v2{2} dependence on pt
    const Double_t y[] = {0.00000,  0.043400, 0.059911, 0.073516,
			  0.089756, 0.105486, 0.117391, 0.128199,
			  0.138013, 0.158271, 0.177726, 0.196383,
			  0.208277, 0.216648, 0.242954, 0.249961,
			  0.240131, 0.269006, 0.207796};
    
    fV22Pt = new TGraph(ptN, ptX, y);
    fV22Pt->SetName(fgkPt2Name);
    fV22Pt->SetMarkerStyle(20);
    fV22Pt->SetMarkerColor(kRed+1);
    l->Add(fV22Pt);
  }

  {
    const Double_t y[] = {0.000000, 0.038646, 0.049824, 0.066662,
			  0.075856, 0.081583, 0.099778, 0.104674,
			  0.118545, 0.131874, 0.152959, 0.155348,
			  0.169751, 0.179052, 0.178532, 0.198851,
			  0.185737, 0.239901, 0.186098};

    // v2{4} dependence on pt 
    fV24Pt = new TGraph(ptN, ptX, y);
    fV24Pt->SetName(fgkPt4Name);
    fV24Pt->SetMarkerStyle(20);
    fV24Pt->SetMarkerColor(kBlue+1);
    l->Add(fV24Pt);
  }

  {
    const Double_t y[] = {0.000000, 0.037071, 0.048566, 0.061083,
			  0.070910, 0.078831, 0.091396, 0.102026,
			  0.109691, 0.124449, 0.139819, 0.155561,
			  0.165701, 0.173678, 0.191149, 0.202015,
			  0.204540, 0.212560, 0.195885};
    // v2{4} dependence on pt (30-40%)
    fV24AltPt = new TGraph(ptN, ptX, y);
    fV24AltPt->SetName(fgkPt4Alt);
    fV24AltPt->SetMarkerStyle(20);
    fV24AltPt->SetMarkerColor(kBlue+1);
    l->Add(fV24AltPt);
  }
  Int_t nb            = 8;
  const Double_t by[] = {0.017855, 0.032440, 0.055818, 0.073137,
			 0.083898, 0.086690, 0.082040, 0.077777};
  {
    // V2 dependence on impact parameter
    const Double_t x[] = {1.75,     4.225,    5.965,    7.765,
			  9.215,    10.46,    11.565,   12.575};
    fV2B = new TGraph(nb, x, by);
    fV2B->SetName(fgkBName);
    fV2B->SetMarkerStyle(20);
    fV2B->SetMarkerColor(kGreen+1);
    l->Add(fV2B);
  }
  {
    // V2 dependence on impact parameter
    const Double_t x[] = { 2.5, 7.5, 15, 25, 35, 45, 55, 65};
    fV2C = new TGraph(nb, x, by);
    fV2C->SetName(fgkCName);
    fV2C->SetMarkerStyle(20);
    fV2C->SetMarkerColor(kGreen+1);
    l->Add(fV2C);
  }
}

//____________________________________________________________________
Double_t
AliForwardFlowWeights::CalcEtaWeight(Double_t eta, Int_t type) const
{
  if (type == 0) return 1;
  return 0.1 * TMath::Gaus(eta, 0, (type == 2 ? 3. : 
				    type == 3 ? 15 : 9));
}
//____________________________________________________________________
Double_t
AliForwardFlowWeights::CalcPidWeight(Int_t id, Int_t type) const
{
  if (type == 0) return 1;
  if (type == 2) return 1.207;
  switch (TMath::Abs(id)) {
  case 211:  return 1.3; break; // pions 
  case 2212: return 1.0; break; // protons 
  default:   return 0.7; break;
  }
  return 1;
}
//____________________________________________________________________
Double_t
AliForwardFlowWeights::CalcPtWeight(Double_t pt, Int_t type) const
{
  switch (type) { 
  case 0: return 1;
  case 2: return fV22Pt->Eval(pt);
  case 3: return fV24AltPt->Eval(pt); // From 30-40
  case 4: return fV24Pt->Eval(pt);
  }
  return 0.5 * (fV22Pt->Eval(pt) + fV24Pt->Eval(pt));
}

//____________________________________________________________________
Double_t
AliForwardFlowWeights::CalcBWeight(Double_t b) const
{
  return fV2B->Eval(b) / fV2B->Eval(10.46);
}

//____________________________________________________________________
Double_t
AliForwardFlowWeights::CalcCentWeight(Double_t c) const
{
  return fV2C->Eval(c) / fV2C->Eval(45);
}

//____________________________________________________________________
Double_t
AliForwardFlowWeights::CalcWeight(Double_t eta, 
				  Double_t pt, 
				  Double_t phi,  
				  Int_t    id, 
				  Double_t phiR, 
				  Double_t cOrB, 
				  Int_t    type, 
				  UShort_t order,
				  UShort_t what) const
{
  Double_t w = 1;
  if (what & kEta)  w *= CalcEtaWeight(eta, type);
  if (what & kPt)   w *= CalcPtWeight(pt, type);
  if (what & kPID)  w *= CalcPidWeight(id, type);
  if      (what & kCent) w *= CalcCentWeight(cOrB);
  else if (what & kB)    w *= CalcBWeight(cOrB);
  
  w *= 20 * 2. * TMath::Cos(order * (phi - phiR));

  return w;
}

//____________________________________________________________________
Double_t
AliForwardFlowWeights::CalcWeight(Double_t eta,  Double_t pt, 
				  Double_t phi,  Int_t id, 
				  Double_t phiR, Double_t b) const
{
  return CalcWeight(eta, pt, phi, id, phiR, b, 1, 2, kEta|kPt|kPID|kB);
  
}

namespace {
  TObject* GetListObject(TList* l, const char* name)
  {
    if (!name || name[0] == '\0') { 
      Error("GetListObject", "No object name");
      return 0;
    }
    if (!l) { 
      Error("GetListObject", "No list");
      return 0;
    }
    TObject* o = l->FindObject(name);
    if (!o) { 
      Error("GetListObject", "Object %s not found in list %s",  
	    name, l->GetName());
      return 0;
    }
    return o;
  }
}
//____________________________________________________________________
AliForwardFlowWeights*
AliForwardFlowWeights::FromList(TList* l)
{
  TObject* pt2 = GetListObject(l, fgkPt2Name);
  TObject* pt4 = GetListObject(l, fgkPt4Name);
  TObject* alt = GetListObject(l, fgkPt4Alt);
  TObject* b   = GetListObject(l, fgkBName);
  TObject* c   = GetListObject(l, fgkCName);

  if (!pt2 || !pt4 || !alt || !b || !c) {
    ::Error("FromList", "One or more histograms not found");
    return 0;
  }
  AliForwardFlowWeights* ret = new AliForwardFlowWeights;
  ret->fV22Pt    = static_cast<TGraph*>(pt2);
  ret->fV24Pt    = static_cast<TGraph*>(pt4);
  ret->fV24AltPt = static_cast<TGraph*>(alt);
  ret->fV2B      = static_cast<TGraph*>(b);
  ret->fV2C      = static_cast<TGraph*>(c);
  
  return ret;
}

//____________________________________________________________________
//
// EOF
//
