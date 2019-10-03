#include "AliFMDMultCuts.h"
#include "AliForwardCorrectionManager.h"
#include "AliFMDCorrELossFit.h"
#include "AliForwardUtil.h"
#include <AliLog.h>
#include <iostream>
#include <TROOT.h>
#include <TParameter.h>
#include <TH2.h>
namespace { 
  Int_t fDebug = 1;
}
ClassImp(AliFMDMultCuts)
#if 0
; // This is for Emacs
#endif

//____________________________________________________________________
AliFMDMultCuts::AliFMDMultCuts()
  : TObject(),
    fMethod(kFixed)
{
  Reset();
}
//____________________________________________________________________
AliFMDMultCuts::AliFMDMultCuts(EMethod method, 
			       Double_t cut1,  
			       Double_t cut2,
			       Double_t cut3,
			       Double_t cut4,
			       Double_t cut5)
  : TObject(),
    fMethod(method)
{
  Set(method, cut1, cut2, cut3, cut4, cut5);
}

//____________________________________________________________________
AliFMDMultCuts::AliFMDMultCuts(const AliFMDMultCuts& o)
  : TObject(o),
    fMethod(o.fMethod)
{
  for (Int_t i = 0; i < 5; i++) fCuts[i] = o.fCuts[i];
}
//____________________________________________________________________
AliFMDMultCuts&
AliFMDMultCuts::operator=(const AliFMDMultCuts& o)
{
  if (&o == this) return *this; 
  fMethod  = o.fMethod;
  for (Int_t i = 0; i < 5; i++) fCuts[i] = o.fCuts[i];
  return *this;
}
//____________________________________________________________________
void
AliFMDMultCuts::Reset()
{
  fMethod = kFixed;
  for (Int_t i = 0; i < 5; i++) fCuts[i] = -1;
}
//____________________________________________________________________
void
AliFMDMultCuts::Set(EMethod method, 
		    Double_t fmd1i, 
		    Double_t fmd2i, 
		    Double_t fmd2o, 
		    Double_t fmd3i, 
		    Double_t fmd3o)
{
  // First, reset
  Reset();
  fMethod      = method; 
  Double_t oFac = 1;
  if      (fMethod == kFixed)       oFac = 1.2;
  else if (fMethod == kMPVFraction) oFac = 1.1;
  fCuts[0] = fmd1i;
  fCuts[1] = fmd2i >= 0 ? fmd2i : fmd1i;
  fCuts[2] = fmd2o >= 0 ? fmd2o : fmd1i * oFac;
  fCuts[3] = fmd3i >= 0 ? fmd3i : fmd1i;
  fCuts[4] = fmd3o >= 0 ? fmd3o : fmd1i * oFac;
}
//____________________________________________________________________
void
AliFMDMultCuts::DepSet(const char* what,
		       EMethod method, 
		       Double_t fmd1i, 
		       Double_t fmd2i, 
		       Double_t fmd2o, 
		       Double_t fmd3i, 
		       Double_t fmd3o)
{
  Warning(what, "*** DEPRECATED - use AliFMDMultCuts::Set instead ***");
  Set(method, fmd1i, fmd2i, fmd2o, fmd3i, fmd3o);
}
  
//____________________________________________________________________
void
AliFMDMultCuts::SetIncludeSigma(Bool_t in) 
{
  Warning("SetIncludeSigma",  
	  "*** DEPRECATED - use AliFMDMultCuts::Set instead ***");
  if (in) { 
    if (fMethod == kLandauWidth) fMethod = kLandauSigmaWidth;
  }
  else {
    if (fMethod == kLandauSigmaWidth) fMethod = kLandauWidth;
  }
}

//____________________________________________________________________
Double_t
AliFMDMultCuts::GetCutParam(UShort_t d, Char_t r) const
{
  // Int_t    idx = (d == 1 ? 0 : 2*(d - 2) + 1 + ((r=='I' || r=='i') ? 0 : 1));
  Int_t idx = -1;
  switch (d) { 
  case 1: idx = 0; break;
  case 2: idx = 1 + ((r == 'I' || r == 'i') ? 0 : 1); break;
  case 3: idx = 3 + ((r == 'I' || r == 'i') ? 0 : 1); break;
  }
  if (idx < 0) return -kBad;
  return fCuts[idx];
}
			    
//____________________________________________________________________
Double_t
AliFMDMultCuts::GetMultCut(UShort_t d, Char_t r, Int_t ieta, 
			   Bool_t errors) const
{
  // 
  // Get the multiplicity cut.  If the user has set fMultCut (via
  // SetMultCut) then that value is used.  If not, then the lower
  // value of the fit range for the enery loss fits is returned.
  // 
  // Return:
  //    Lower cut on multiplicity
  //
  DGUARD(fDebug, 5, "Get mult cut for FMD%d%c (method %d) @ etabin=%d", 
	 d, r, fMethod, ieta);
  Double_t param = GetCutParam(d, r);
  if (param < 0) {
    Warning("GetMultCut", "Got bad cut parameter for FMD%d%c ieta=%d",
	    d, r, ieta);
    return -kBad;
  }

  // If we're using a fixed cut, just return
  if (fMethod == kFixed) {
    DMSG(fDebug, 5, "-> %8.4f", param);
    return param;
  }

  // Bad value 
  Double_t rcut = -kBad;

  // Get the energy loss fits 
  AliForwardCorrectionManager&  fcm = 
    AliForwardCorrectionManager::Instance();
  const AliFMDCorrELossFit* fits = fcm.GetELossFit();
  if (fits) {
    switch (fMethod) { 
    case kMPVFraction:
      // Return fMPVFraction * mpv 
      // rcut = fits->GetLowerBound(d, r, ieta, param); break;
      rcut = fits->GetMpvCut(d,r,ieta,param); break;
    case kLandauWidth:
      // Return MPV - fNXi * xi
      // rcut = fits->GetLowerBound(d, r, ieta, param, errors,false);
      rcut = fits->GetXiCut(d, r, ieta, param);
      break;
    case kLandauSigmaWidth:
      // Return MPV - fNXi * (xi+sigma)
      // rcut = fits->GetLowerBound(d, r, ieta, param, errors,true);
      rcut = fits->GetXiSigmaCut(d, r, ieta, param);
      break;
    case kProbability:
      // Return probability cut 
      // rcut = fits->GetLowerBound(d, r, ieta, param, true); break;
      rcut = fits->GetProbabilityCut(d, r, ieta, param); break;
    case kLandauSigmaAvg:
      // Return MPV - fNxi*WeightedAverage(xi,sigma)
      rcut = fits->GetAvgXiSigmaCut(d,r,ieta,param); break;
    default:
      // Return lower fit boundary
      rcut = fits->GetLowCut(); break;
    }
  }
  else 
    Warning("GetMultCut", "No energy loss fits obtained from manager");

  DMSG(fDebug, 5, "-> %8.4f", rcut);
    
  return rcut;
}
    
//____________________________________________________________________
Double_t
AliFMDMultCuts::GetMultCut(UShort_t d, Char_t r, Double_t eta,
			   Bool_t errors) const
{
  // 
  // Get the multiplicity cut.  If the user has set fMultCut (via
  // SetMultCut) then that value is used.  If not, then the lower
  // value of the fit range for the enery loss fits is returned.
  // 
  // Return:
  //    Lower cut on multiplicity
  //
  DGUARD(fDebug, 5, "Get mult cut for FMD%d%c @ eta=%8.4f", 
	 d, r, eta);
  AliForwardCorrectionManager&  fcm  = AliForwardCorrectionManager::Instance();
  const AliFMDCorrELossFit*     fits = fcm.GetELossFit();
  Int_t                         iEta = fits ? fits->FindEtaBin(eta) : 1;
  DMSG(fDebug, 5, "bin=%4d", iEta);
  return GetMultCut(d, r, iEta, errors);
}
//____________________________________________________________________
const char*
AliFMDMultCuts::GetMethodString(Bool_t latex) const
{
  return Method2String(fMethod, latex);
}

//____________________________________________________________________
const char*
AliFMDMultCuts::Method2String(EMethod method, Bool_t latex)
{
  switch (method) {
  case kFixed:            
    return latex ? "c=X" : "fixed value";
  case kMPVFraction:      
    return latex ? "c=X#times#Delta_{p}":"fraction of MPV";
  case kFitRange:         
    return latex ? "range" : "c: lower fit bound";
  case kLandauWidth:      
    return latex ? "c=#Delta_{p}-X#times#xi" : "landau width";
  case kLandauSigmaWidth: 
    return latex ? "c=#Delta_{p}-X#times(#xi+#sigma)" : "landau+sigma width";
  case kProbability:      
    return latex ? "c:P(#Delta<c)<X" : "probability";
  case kLandauSigmaAvg:
    return latex ? "x=#Delta_{p}-X#times#bar{#xi+#sigma}" :
      "landau+sigma avg";
  }
  return latex ? "c:?" : "unknown";
} 
//____________________________________________________________________
AliFMDMultCuts::EMethod
AliFMDMultCuts::String2Method(const char* str)
{
  TString m(str);
  if (m.EqualTo("fixed value")     || m.Contains("fix")) 
    return kFixed;
  else if (m.EqualTo("fraction of mpv") || m.Contains("mpv")) 
    return kMPVFraction;
  else if (m.EqualTo("fit range")       || m.Contains("fit")) 
    return kFitRange;
  else if (m.EqualTo("landau width")    || m.Contains("xi") || 
	   m.Contains("width")) return kLandauWidth;
  else if (m.EqualTo("landau+sigma avg") || m.Contains("avg")) 
    return kLandauSigmaAvg;
  else if (m.EqualTo("landau+sigma width") || m.Contains("sig")) 
    return kLandauSigmaWidth;
  else if (m.EqualTo("probability")        || m.Contains("prob")) 
    return kProbability;
  return kFixed;
 } 
  
//____________________________________________________________________
void
AliFMDMultCuts::FillHistogram(TH2* h) const
{
  DGUARD(fDebug, 5, "Fill Histogram %s with cuts", h->GetName());
  // AliInfoF("Caching multiplicity cuts (%s)", h->GetName());
  TAxis* yAxis = h->GetYaxis(); 
  for (Int_t iy = 1; iy <= yAxis->GetNbins(); iy++) { 
    TString lab(yAxis->GetBinLabel(iy));
    lab.Remove(0,3);
    UShort_t det = lab.Atoi();
    Char_t   rng = lab[1];
    // Printf("Filling for FMD%d%c (bin # %d) %s", det, rng, iy, lab.Data());
    DMSG(fDebug, 5, "FMD%d%c", det, rng);
    // AliInfoF("FMD%d%c", det, rng);
    for (Int_t ix = 1; ix <= h->GetNbinsX(); ix++) {
      Double_t eta = h->GetXaxis()->GetBinCenter(ix);
      Double_t c   = GetMultCut(det, rng, eta,  false);
      DMSG(fDebug, 10, "FMD%s bin=%4d -> eta=%8.4f -> %8.4f", 
	   lab.Data(), ix, eta, c);
      // Double_t c = GetMultCut(det, rng, ix, false);
      if (c > 0) h->SetBinContent(ix, iy, c);
    }
  }
}

//____________________________________________________________________
void
AliFMDMultCuts::Output(TList* l, const char* name) const
{
  TList* ll = l;
  if (name && name[0] != '\0') { 
    ll = new TList;
    ll->SetName(name);
    ll->SetOwner();
    l->Add(ll);
  }
    

  ll->Add(AliForwardUtil::MakeParameter("method", UShort_t(fMethod)));
  ll->Add(AliForwardUtil::MakeParameter("fmd1i", fCuts[0]));
  ll->Add(AliForwardUtil::MakeParameter("fmd2i", fCuts[1]));
  ll->Add(AliForwardUtil::MakeParameter("fmd2o", fCuts[2]));
  ll->Add(AliForwardUtil::MakeParameter("fmd3i", fCuts[3]));
  ll->Add(AliForwardUtil::MakeParameter("fmd3o", fCuts[4]));
}
//____________________________________________________________________
Bool_t
AliFMDMultCuts::Input(TList* l, const char* name)
{
  if (!l) return false;
  TList* ll = l;
  if (name && name[0] != '\0') { 
    ll = static_cast<TList*>(l->FindObject(name));
  }
  if (!ll) return false;
  
  TObject* meth    = ll->FindObject("method");
  TObject* fmd1i   = ll->FindObject("fmd1i");
  TObject* fmd2i   = ll->FindObject("fmd2i");
  TObject* fmd2o   = ll->FindObject("fmd2o");
  TObject* fmd3i   = ll->FindObject("fmd3i");
  TObject* fmd3o   = ll->FindObject("fmd3o");
  
  UShort_t methNum = 0;
  
  AliForwardUtil::GetParameter(meth,  methNum);
  switch (methNum) { 
  case 0: fMethod = kFixed;            break;
  case 1: fMethod = kMPVFraction;      break;
  case 2: fMethod = kFitRange;         break;
  case 3: fMethod = kLandauWidth;      break;
  case 4: fMethod = kLandauSigmaWidth; break;
  case 5: fMethod = kProbability;      break;
  }

  AliForwardUtil::GetParameter(fmd1i, fCuts[0]);
  AliForwardUtil::GetParameter(fmd2i, fCuts[1]);
  AliForwardUtil::GetParameter(fmd2o, fCuts[2]);
  AliForwardUtil::GetParameter(fmd3i, fCuts[3]);
  AliForwardUtil::GetParameter(fmd3o, fCuts[4]);
  
  return true;
}
#define PF(N,V,...)					\
  AliForwardUtil::PrintField(N,V, ## __VA_ARGS__)
#define PFB(N,FLAG)				\
  do {									\
    AliForwardUtil::PrintName(N);					\
    std::cout << std::boolalpha << (FLAG) << std::noboolalpha << std::endl; \
  } while(false)
#define PFV(N,VALUE)					\
  do {							\
    AliForwardUtil::PrintName(N);			\
    std::cout << (VALUE) << std::endl; } while(false)
  
//____________________________________________________________________
void
AliFMDMultCuts::Print(Option_t*) const
{
  gROOT->IncreaseDirLevel();
  PFV("Method used", GetMethodString());
  gROOT->IncreaseDirLevel();
  PFV("FMD1i", GetCutParam(1,'I'));
  PFV("FMD2i", GetCutParam(2,'I'));
  PFV("FMD2o", GetCutParam(2,'O'));
  PFV("FMD3i", GetCutParam(3,'I'));
  PFV("FMD3o", GetCutParam(3,'O'));
  gROOT->DecreaseDirLevel();
  gROOT->DecreaseDirLevel();
}
//____________________________________________________________________
//
// EOF
//
