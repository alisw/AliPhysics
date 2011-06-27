#include <AliFMDCalibDrawer.h>
#include <TH1.h>
#include <TH2.h>
#include <TPad.h>
#include <TCanvas.h>
#include <AliFMDParameters.h>
#include <AliCDBManager.h>
#include <AliLog.h>
// #include <AliFMDCalibPedestal.h>
// #include <AliFMDCalibSampleRate.h>
// #include <AliFMDCalibStripRange.h>

//____________________________________________________________________
void
AliFMDCalibDrawer::Init(Int_t runNo, const char* ocdb)
{
  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetRun(runNo);
  if (ocdb && ocdb[0] != '\0') cdb->SetDefaultStorage(ocdb);

  AliFMDParameters* pars = AliFMDParameters::Instance();
  pars->Init(kTRUE);
}

//____________________________________________________________________
Double_t
AliFMDCalibDrawer::GetHistMax(EWhat what) const
{
  switch (what) {
  case kPedestal:  return 150;
  case kNoise:     return 5;
  case kGain:      return 4;
  case kDead:      return 1.5;
  case kRate:      return 5;
  case kRange:     return 128;
  case kZeroSuppression:  return 10;
  }
  return 1024;
}
//____________________________________________________________________
Double_t
AliFMDCalibDrawer::GetHistMin(EWhat what) const
{
  switch (what) {
  case kPedestal:  return 0;
  case kNoise:     return 0;
  case kGain:      return 0;
  case kDead:      return -.5;
  case kRate:      return 0;
  case kRange:     return -1;
  case kZeroSuppression:  return -1;
  }
  return -1;
}
  
//____________________________________________________________________
const char*
AliFMDCalibDrawer::GetHistName(EWhat what) const
{
  switch (what) {
  case kPedestal:  return "peds"; 
  case kNoise:     return "noise";
  case kGain:      return "gain";
  case kDead:      return "dead";
  case kRate:      return "rate";
  case kRange:     return "range";
  case kZeroSuppression:  return "thrs";
  }
  return "unknown";
}
//____________________________________________________________________
const char*
AliFMDCalibDrawer::GetHistTitle(EWhat what) const
{
  switch (what) {
  case kPedestal:  return "Pedestal & noise"; 
  case kNoise:     return "Noise";
  case kGain:      return "Gain";
  case kDead:      return "Dead map";
  case kRate:      return "Sample rate";
  case kRange:     return "Strip range";
  case kZeroSuppression:  return "ZS threshold";
  }
  return "Unknown";
}
  
//____________________________________________________________________
Int_t
AliFMDCalibDrawer::GetRingColor(UShort_t d, Char_t r) const
{
  return ((d == 1 ? kRed : (d == 2 ? kGreen : kBlue))
	  + ((r == 'I' || r == 'i') ? 2 : -3));
}
//____________________________________________________________________
void
AliFMDCalibDrawer::SetAttributes(TH1* ret, EWhat what, 
				 UShort_t d, Char_t r) const
{
  ret->SetFillColor(GetRingColor(d,r));
  ret->SetMarkerColor(GetRingColor(d,r));
  ret->SetLineColor(GetRingColor(d,r));
  ret->SetFillStyle(3001);
  ret->SetMarkerStyle(20);
  ret->SetLineStyle(1);
  ret->SetStats(0);
  ret->SetDirectory(0);
  ret->SetMinimum(GetHistMin(what));
  ret->SetMaximum(GetHistMax(what));
}
//____________________________________________________________________
TH2D*
AliFMDCalibDrawer::Make2D(EWhat what, UShort_t d, Char_t r) const
{
  UShort_t q  = (r == 'I' || r == 'i') ? 0 : 1;
  UShort_t nY = q == 0 ?  20 :  40;
  UShort_t nX = q == 0 ? 512 : 256;
  TString  n(Form("%s_FMD%d%c", GetHistName(what), d, r));
  TString  t(Form("%s for FMD%d%c", GetHistTitle(what),  d, r));
  TH2D*    ret = new TH2D(n, t, nX, 0, nX, nY, 0, nY);
  ret->SetXTitle("Strip #");
  ret->SetYTitle("Sector #");
  ret->SetZTitle(GetHistTitle(what));
  ret->GetXaxis()->SetNdivisions(1600 + (q == 0 ? 4 : 2),false);
  ret->GetYaxis()->SetNdivisions(nY, false);
  SetAttributes(ret, what, d, r);
  return ret;

}

//____________________________________________________________________
TH1D*
AliFMDCalibDrawer::Make1D(EWhat what, UShort_t d, Char_t r, UShort_t s) const  
{
  UShort_t q    = (r == 'I' || r == 'i') ? 0 : 1;
  UShort_t nStr = (r == 'I' || r == 'i') ? 512 : 256;
  TString  n(Form("%s_FMD%d%c_%02d", GetHistName(what), d, r, s));
  TString  t(Form("%s for FMD%d%c[%2d]", GetHistTitle(what), d, r, s));
  TH1D* ret = new TH1D(n, t, nStr, 0, nStr);
  ret->SetXTitle("Strip #");
  ret->SetYTitle(GetHistTitle(what));
  ret->GetXaxis()->SetNdivisions(1600 + (q == 0 ? 4 : 2),false);
  SetAttributes(ret, what, d, r);
  return ret;
}

//____________________________________________________________________
void
AliFMDCalibDrawer::GetNumber(EWhat what, UShort_t d, Char_t r, UShort_t s, 
			     UShort_t t, Double_t& val, Double_t& err) const
{
  AliFMDParameters* pars = AliFMDParameters::Instance();
  err = 0;
  switch (what) {
  case kPedestal:
    val = pars->GetPedestal(d, r, s, t);
    err = pars->GetPedestalWidth(d, r, s, t);
    break;
  case kNoise: val = pars->GetPedestalWidth(d, r, s, t); break;
  case kGain:  val = pars->GetPulseGain(d, r, s, t); break;
  case kDead:  val = pars->IsDead(d, r, s, t) ? 0 : 1; break;
  case kRate:  val = pars->GetSampleRate(d, r, s, t); break;
  case kRange:
    val = pars->GetMaxStrip(d,r,s,t);
    err = pars->GetMinStrip(d,r,s,t);
    break;
  case kZeroSuppression:  val = pars->GetZeroSuppression(d,r,s,t); break;
  default: 
    AliError(Form("Unknown quantity: %d - nothing returned", what));
    break;
  }
}

//____________________________________________________________________
TH1*
AliFMDCalibDrawer::FillRing(EWhat what, UShort_t d, Char_t r) const
{
  if (what == kRate || what == kRange) { 
    TH1D* ret = new TH1D(Form("%s_FMD%d%c", GetHistName(what), d, r), 
			 Form("%s for FMD%d%c", GetHistTitle(what), d, r), 
			 2, -.5, 1.5);
    UShort_t nSec = (r == 'I' || r == 'i') ? 20 : 40;
    ret->GetXaxis()->SetBinLabel(1, Form("Top [%02d-%02d]", 0,      nSec/2-1));
    ret->GetXaxis()->SetBinLabel(2, Form("Bottom [%02d-%02d]",nSec/2-1,nSec-1));
    ret->SetYTitle(GetHistTitle(what));
    SetAttributes(ret, what, d, r);
    if (what == kRate) { 
      ret->SetBarWidth(0.8);
      ret->SetBarOffset(0.1);
    }
    else 
      ret->SetMarkerSize(0);
    for (UShort_t s = 0; s < ret->GetNbinsX(); s++) { 
      Double_t val; 
      Double_t err; 
      GetNumber(what, d, r, s*(nSec/2), 0, val, err);
      if (what == kRange) { 
	Double_t min = err;
	err          = (val - err) / 2;
	val          = min + err;
      }
      ret->SetBinContent(s+1, val);
      ret->SetBinError(s+1, err);
    }
    return ret;
  }
  if (what == kZeroSuppression) { 
    UShort_t nY  = (r == 'I' || r == 'i') ?  20 :  40;
    UShort_t nX  = ((r == 'I' || r == 'i') ? 512 : 256) / 128;
    TH2D*    ret = new TH2D(Form("%s_FMD%d%c", GetHistName(what), d, r), 
			    Form("%s for FMD%d%c", GetHistTitle(what), d, r),
			    nX, 0, nX, nY, 0, nY);
    ret->SetXTitle("Channel #");
    ret->SetYTitle("Sector #");
    ret->SetZTitle(GetHistTitle(what));
    ret->GetXaxis()->SetNdivisions(nX, false);
    ret->GetYaxis()->SetNdivisions(nY, false);
    SetAttributes(ret, what, d, r);
    for (UShort_t s = 0; s < ret->GetNbinsY(); s++) { 
      for (UShort_t t = 0; t < ret->GetNbinsX(); t++) { 
	Double_t val; 
	Double_t err; 
	GetNumber(what, d, r, s, t*128, val, err);
	ret->SetBinContent(t+1, s+1, val);
	ret->SetBinError(t+1, s+1, err);
      }
    }
    return ret;
  }
    
  TH2D* ret = Make2D(what, d, r);
  for (UShort_t s = 0; s < ret->GetNbinsY(); s++) { 
    for (UShort_t t = 0; t < ret->GetNbinsX(); t++) { 
      Double_t val; 
      Double_t err; 
      GetNumber(what, d, r, s, t, val, err);
      ret->SetBinContent(t+1, s+1, val);
      ret->SetBinError(t+1, s+1, err);
    }
  }
  return ret;
}
//____________________________________________________________________
TH1*
AliFMDCalibDrawer::FillSector(EWhat what,UShort_t d, Char_t r, UShort_t s) const
{
  if (what == kZeroSuppression) { 
    UShort_t nX  = ((r == 'I' || r == 'i') ? 512 : 256) / 128;
    TH1D*    ret = new TH1D(Form("%s_FMD%d%c[%02d]", 
				 GetHistName(what), d, r, s), 
			    Form("%s for FMD%d%c[%02d]", 
				 GetHistTitle(what), d, r, s),
			    nX, 0, nX);
    ret->SetXTitle("Channel #");
    ret->SetYTitle(GetHistTitle(what));
    SetAttributes(ret, what, d, r);
    ret->GetXaxis()->SetNdivisions(nX, false);
    for (UShort_t t = 0; t < ret->GetNbinsX(); t++) { 
      Double_t val; 
      Double_t err; 
      GetNumber(what, d, r, s, t*128, val, err);
      ret->SetBinContent(t+1, s+1, val);
      if (err >= 0) ret->SetBinError(t+1, s+1, err);
    }
    return ret;
  }

  TH1D* ret = Make1D(what, d, r, s);
  for (UShort_t t = 0; t < ret->GetNbinsX(); t++) { 
    Double_t val; 
    Double_t err; 
    GetNumber(what, d, r, s, t, val, err);
    ret->SetBinContent(t+1, s+1, val);
    if (err >= 0) ret->SetBinError(t+1, s+1, err);
  }
  return ret;
}

//____________________________________________________________________
void
AliFMDCalibDrawer::DrawOne(EWhat what, Short_t d, Char_t r, 
			   Short_t s, Short_t /*t*/) const
{
  Info("DrawOne", "Will draw %s for D=%d, R=%c, S=%d",
       GetHistTitle(what), d, (r == '\0' ? '-' : r), s);
  TVirtualPad* tp = gPad;
  if (!tp) tp = TCanvas::MakeDefCanvas();
  if (d <= 0 || d > 3) { 
    // All detectors - need to split current pad
    tp->Divide(1,3,0,0);
    for (UShort_t dd = 1; dd <= 3; dd++) { 
      tp->cd(dd);
      // Call our selves 
      DrawOne(what, dd, '\0', -1, -1);
    }
  }
  else if ((r != 'I' && r != 'i') && (r != 'O' && r != 'o')) { 
    // All rings in a detector.  Split current pad in two 
    tp->SetName(Form("FMD%d", d));
    tp->Divide(2,1,0,0);
    for (UShort_t q = 0; q < 2; q++) { 
      tp->cd(q+1);
      Char_t       rr = q == 0 ? 'I' : 'O';
      // Call ourselves 
      DrawOne(what, d, rr, -1, -1);
    }
  }
  else if (what == kRate || what == kRange || s < 0 || s > 39) {
    // A single ring - and special case for sample rate and strip range
    tp->SetName(Form("FMD%d%c", d, r));
    tp->SetFillColor(0);
    tp->SetFillStyle(0);
    tp->SetRightMargin(0.10);
    if (d == 1 && (r == 'O' || r == 'o')) return;
    TH1* h = FillRing(what, d, r);
    switch (what) { 
    case kRate:  h->Draw("bar"); break;
    case kRange: h->Draw("e3"); break;
    default: h->Draw("colz"); break;
    }
  }
  else {
    // Finally, we're down to a single sector.  In this case, we just do
    // a single 1D histogram 
    TH1* h = FillSector(what, d, r, s);
    h->Draw("hist e");
  }
  tp->Update();
  tp->Modified();
  tp->cd();
}



//____________________________________________________________________
//
// EOF
// 
