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
    fMPVFraction(0), 
    fNXi(0), 
    fIncludeSigma(false),
    fProbability(0)
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
    fMPVFraction(0), 
    fNXi(0), 
    fIncludeSigma(false),
    fProbability(0)
{
  Set(method, cut1, cut2, cut3, cut4, cut5);
}

//____________________________________________________________________
AliFMDMultCuts::AliFMDMultCuts(const AliFMDMultCuts& o)
  : TObject(o),
    fMPVFraction(o.fMPVFraction), 
    fNXi(o.fNXi), 
    fIncludeSigma(o.fIncludeSigma),
    fProbability(o.fProbability)
{
  for (Int_t i = 0; i < 5; i++) fMultCuts[i] = o.fMultCuts[i];
}
//____________________________________________________________________
AliFMDMultCuts&
AliFMDMultCuts::operator=(const AliFMDMultCuts& o)
{
  if (&o == this) return *this; 
  fMPVFraction  = o.fMPVFraction;
  fNXi          = o.fNXi;
  fIncludeSigma = o.fIncludeSigma;
  fProbability  = o.fProbability;
  for (Int_t i = 0; i < 5; i++) fMultCuts[i] = o.fMultCuts[i];
  return *this;
}
//____________________________________________________________________
void
AliFMDMultCuts::Reset()
{
  for (Int_t i = 0; i < 5; i++) fMultCuts[i] = -1;
  fMPVFraction  = -1;
  fNXi          = -1;
  fIncludeSigma = false;
  fProbability  = -1;
}
//____________________________________________________________________
void
AliFMDMultCuts::Set(EMethod method, 
		    Double_t cut1, 
		    Double_t cut2, 
		    Double_t cut3, 
		    Double_t cut4, 
		    Double_t cut5)
{
  // First, reset
  Reset();
  
  // Then switch on method 
  switch(method) { 
  case kFixed:  
    if (cut2 < 0) SetMultCuts(cut1, cut1, cut1*1.2, cut1*1.2, cut1);
    else          SetMultCuts(cut1, cut2, cut3, cut4, cut5);
    break;
  case kMPVFraction: 
    SetMPVFraction(cut1);
    break;
  case kFitRange:
    break;
  case kLandauWidth:
    SetNXi(cut1);
    SetIncludeSigma(cut2 > 0);
    break;
  case kProbability:
    SetProbability(cut1);
    break;
  }
}
  

//____________________________________________________________________
Double_t
AliFMDMultCuts::GetFixedCut(UShort_t d, Char_t r) const
{
  // Int_t    idx = (d == 1 ? 0 : 2*(d - 2) + 1 + ((r=='I' || r=='i') ? 0 : 1));
  Int_t idx = -1;
  switch (d) { 
  case 1: idx = 0; break;
  case 2: idx = 1 + ((r == 'I' || r == 'i') ? 0 : 1); break;
  case 3: idx = 3 + ((r == 'I' || r == 'i') ? 0 : 1); break;
  }
  if (idx < 0) return -1024;
  return fMultCuts[idx];
}

//____________________________________________________________________
void
AliFMDMultCuts::SetMultCuts(Double_t fmd1i, 
			    Double_t fmd2i, 
			    Double_t fmd2o, 
			    Double_t fmd3i, 
			    Double_t fmd3o)
{
  fMultCuts[0] = fmd1i;
  fMultCuts[1] = fmd2i >= 0 ? fmd2i : fmd1i;
  fMultCuts[2] = fmd2o >= 0 ? fmd2o : fmd1i;
  fMultCuts[3] = fmd3i >= 0 ? fmd3i : fmd1i;
  fMultCuts[4] = fmd3o >= 0 ? fmd3o : fmd1i;
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
  UShort_t meth = GetMethod();
  DGUARD(fDebug, 5, "Get mult cut for FMD%d%c (method %d) @ etabin=%d", 
	 d, r, meth, ieta);
  Double_t rcut = -1024;
  if (meth == kFixed) rcut = GetFixedCut(d, r);

  if (rcut < 0) {
    // Get the energy loss fits 
    AliForwardCorrectionManager&  fcm = 
      AliForwardCorrectionManager::Instance();
    const AliFMDCorrELossFit* fits = fcm.GetELossFit();
    if (fits) {
      switch (meth) { 
      case kMPVFraction:
	// Return fMPVFraction * mpv 
	rcut = fits->GetLowerBound(d, r, ieta, fMPVFraction); break;
      case kLandauWidth:
	// Return MPV - fNXi * xi
	rcut = fits->GetLowerBound(d, r, ieta, fNXi, errors,fIncludeSigma);
	break;
      case kProbability:
	// Return probability cut 
	rcut = fits->GetLowerBound(d, r, ieta, fProbability, true); break;
      default:
	// Return lower fit boundary
	rcut = fits->GetLowCut(); break;
      }
    }
    else 
      Warning("GetMultCut", "No energy loss fits obtained from manager");
  }
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
UShort_t
AliFMDMultCuts::GetMethod() const
{
  return (fMultCuts[0] >= 0 ? kFixed : // Fixed
	  fMPVFraction >  0 ? kMPVFraction : // Fraction MPV
	  fNXi         >  0 ? kLandauWidth :  // Width 
	  fProbability >  0 ? kProbability : 
	  kFitRange); // Fit range
}
//____________________________________________________________________
const char*
AliFMDMultCuts::GetMethodString() const
{
  switch (GetMethod()) {
  case kFixed:       return "fixed value";
  case kMPVFraction: return "fraction of MPV";
  case kFitRange:    return "fit range";
  case kLandauWidth: return "landau width";
  case kProbability: return "probability";
  }
  return "unknown";
 } 
  
//____________________________________________________________________
void
AliFMDMultCuts::FillHistogram(TH2* h) const
{
  DGUARD(fDebug, 5, "Fill Histogram %s with cuts", h->GetName());
  AliInfoF("Caching multiplicity cuts (%s)", h->GetName());
  TAxis* yAxis = h->GetYaxis(); 
  for (Int_t iy = 1; iy <= yAxis->GetNbins(); iy++) { 
    TString lab(yAxis->GetBinLabel(iy));
    lab.Remove(0,3);
    UShort_t det = lab.Atoi();
    Char_t   rng = lab[1];
    // Printf("Filling for FMD%d%c (bin # %d) %s", det, rng, iy, lab.Data());
    AliInfoF("FMD%d%c", det, rng);
    for (Int_t ix = 1; ix <= h->GetNbinsX(); ix++) {
      Double_t eta = h->GetXaxis()->GetBinCenter(ix);
      Double_t c   = GetMultCut(det, rng, eta,  false);
      DMSG(fDebug, 5, "FMD%s bin=%4d -> eta=%8.4f -> %8.4f", 
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
    

  ll->Add(AliForwardUtil::MakeParameter("nXi", fNXi));
  ll->Add(AliForwardUtil::MakeParameter("frac", fMPVFraction));
  ll->Add(AliForwardUtil::MakeParameter("sigma", fIncludeSigma));
  ll->Add(AliForwardUtil::MakeParameter("probability", fProbability));
  ll->Add(AliForwardUtil::MakeParameter("method", GetMethod()));
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
  
  TObject* nXi   = ll->FindObject("nXi");
  TObject* frac  = ll->FindObject("frac");
  TObject* sigma = ll->FindObject("sigma");
  TObject* prob  = ll->FindObject("probability");
  if (!nXi || !frac || !sigma) return false;
  AliForwardUtil::GetParameter(nXi, fNXi);
  AliForwardUtil::GetParameter(frac, fMPVFraction);
  AliForwardUtil::GetParameter(sigma, fIncludeSigma);
  AliForwardUtil::GetParameter(prob, fProbability);
  
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
  PF("Fixed cuts","");
  gROOT->IncreaseDirLevel();
  PFV("FMD1i", GetFixedCut(1,'I'));
  PFV("FMD2i", GetFixedCut(2,'I'));
  PFV("FMD2o", GetFixedCut(2,'O'));
  PFV("FMD3i", GetFixedCut(3,'I'));
  PFV("FMD3o", GetFixedCut(3,'O'));
  gROOT->DecreaseDirLevel();
  PFV("N xi factor",		fNXi);
  PFB("Include sigma in cut",	fIncludeSigma);
  PFV("MPV fraction",		fMPVFraction);
  PFV("Probability",		fProbability);
  gROOT->DecreaseDirLevel();
}
//____________________________________________________________________
//
// EOF
//
