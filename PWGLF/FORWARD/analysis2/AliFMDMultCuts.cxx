#include "AliFMDMultCuts.h"
#include "AliForwardCorrectionManager.h"
#include "AliFMDCorrELossFit.h"
#include <iostream>
#include <TROOT.h>
#include <TParameter.h>

ClassImp(AliFMDMultCuts)
#if 0
; // This is for Emacs
#endif

//____________________________________________________________________
AliFMDMultCuts::AliFMDMultCuts()
  : TObject(),
    fMPVFraction(0), 
    fNXi(0), 
    fIncludeSigma(true)
{
  for (Int_t i = 0; i < 5; i++) fMultCuts[i] = 0;
}
//____________________________________________________________________
AliFMDMultCuts::AliFMDMultCuts(const AliFMDMultCuts& o)
  : TObject(o),
    fMPVFraction(o.fMPVFraction), 
    fNXi(o.fNXi), 
    fIncludeSigma(o.fIncludeSigma)
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
  for (Int_t i = 0; i < 5; i++) fMultCuts[i] = o.fMultCuts[i];
  return *this;
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
  if (idx < 0) return 1024;
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
  Double_t rcut = GetFixedCut(d, r);
  if (rcut > 0) return rcut;

  AliForwardCorrectionManager&  fcm = AliForwardCorrectionManager::Instance();
  AliFMDCorrELossFit* fits = fcm.GetELossFit();
  if (fMPVFraction > 0) 
    return fits->GetLowerBound(d, r, ieta, fMPVFraction);

  if (fNXi < 0) return fits->GetLowCut();

  return fits->GetLowerBound(d, r, ieta, fNXi, errors, fIncludeSigma);
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
  AliForwardCorrectionManager&  fcm  = AliForwardCorrectionManager::Instance();
  AliFMDCorrELossFit*           fits = fcm.GetELossFit();
  Int_t                         iEta = fits->FindEtaBin(eta);
  
  return GetMultCut(d, r, iEta, errors);
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
    
  TParameter<double>* nXi    = new TParameter<double>("nXi", fNXi);
  TParameter<double>* frac   = new TParameter<double>("frac", fMPVFraction);
  TNamed*             sigma  = new TNamed("sigma", fIncludeSigma ? 
					  "included" : "excluded");
  sigma->SetUniqueID(fIncludeSigma);
  ll->Add(nXi);
  ll->Add(frac);
  ll->Add(sigma);
}
  
//____________________________________________________________________
void
AliFMDMultCuts::Print(Option_t*) const
{
  char ind[gROOT->GetDirLevel()+1];
  for (Int_t i = 0; i < gROOT->GetDirLevel(); i++) ind[i] = ' ';
  ind[gROOT->GetDirLevel()] = '\0';
  std::cout << std::boolalpha 
	    << ind << "  Method used:           " 
	    << (fMultCuts[0] >= 0 ? "fixed value" : 
		fMPVFraction > 0 ? "fraction of MPV" : 
		fNXi < 0 ? "fit range" : 
		"landau width") << '\n'
	    << ind << "  Fixed cuts:            "
	    << "FMD1i=" << GetFixedCut(1,'I') << " "
	    << "FMD2i=" << GetFixedCut(2,'I') << " "
	    << "FMD2o=" << GetFixedCut(2,'O') << " "
	    << "FMD3i=" << GetFixedCut(3,'I') << " "
	    << "FMD3o=" << GetFixedCut(3,'O') << "\n"
	    << ind << "  N xi factor:           " << fNXi    << '\n'
	    << ind << "  Include sigma in cut:  " << fIncludeSigma << '\n'
	    << ind << "  MPV fraction:          " << fMPVFraction 
	    << std::noboolalpha << std::endl;
}
//____________________________________________________________________
//
// EOF
//
