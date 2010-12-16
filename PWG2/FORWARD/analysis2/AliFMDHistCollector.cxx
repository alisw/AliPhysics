// Calculate the total energy loss 
#include "AliFMDHistCollector.h"
#include <AliESDFMD.h>
#include <TAxis.h>
#include <TList.h>
#include <TMath.h>
// #include "AliFMDAnaParameters.h"
#include "AliForwardCorrectionManager.h"
#include "AliLog.h"
#include <TH2D.h>
#include <TH1I.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TArrayI.h>
#include <TROOT.h>
#include <iostream>
#include <iomanip>

ClassImp(AliFMDHistCollector)
#if 0
; // For Emacs
#endif 


//____________________________________________________________________
AliFMDHistCollector&
AliFMDHistCollector::operator=(const AliFMDHistCollector& o)
{
  TNamed::operator=(o);

  fNCutBins       = o.fNCutBins;
  fCorrectionCut  = o.fCorrectionCut;
  fFirstBins      = o.fFirstBins;
  fLastBins       = o.fLastBins;
  fDebug          = o.fDebug;

  return *this;
}

//____________________________________________________________________
void
AliFMDHistCollector::Init(const TAxis& vtxAxis)
{
  // AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  AliForwardCorrectionManager& fcm = AliForwardCorrectionManager::Instance();

  UShort_t nVz = vtxAxis.GetNbins();
  fFirstBins.Set(5*nVz);
  fLastBins.Set(5*nVz);

  // Find the eta bin ranges 
  for (UShort_t iVz = 1; iVz <= nVz; iVz++) {
    // Find the first and last eta bin to use for each ring for 
    // each vertex bin.   This is instead of using the methods 
    // provided by AliFMDAnaParameters 
    for (Int_t iIdx = 0; iIdx < 5; iIdx++) {
      UShort_t d = 0;
      Char_t   r = 0;
      GetDetRing(iIdx, d, r);
      
      // Get the background object 
      // TH2F* bg    = pars->GetBackgroundCorrection(d,r,iVz);
      TH2D* bg    = fcm.GetSecondaryMap()->GetCorrection(d,r,iVz);
      Int_t nEta  = bg->GetNbinsX();
      Int_t first = nEta+1;
      Int_t last  = 0;

      // Loop over the eta bins 
      for (Int_t ie = 1; ie <= nEta; ie++) { 
	if (bg->GetBinContent(ie,1) < fCorrectionCut) continue;

	// Loop over the phi bins to make sure that we 
	// do not have holes in the coverage 
	bool ok = true;
	for (Int_t ip = 1; ip <= bg->GetNbinsY(); ip++) { 
	  if (bg->GetBinContent(ie,ip) < 0.001) {
	    ok = false;
	    continue;
	  }
	}
	if (!ok) continue;

	first = TMath::Min(ie, first);
	last  = TMath::Max(ie, last);
      }
      
      // Store the result for later use 
      fFirstBins[(iVz-1)*5+iIdx] = first;
      fLastBins[(iVz-1)*5+iIdx]  = last;
    } // for j 
  }
}

//____________________________________________________________________
Int_t
AliFMDHistCollector::GetIdx(UShort_t d, Char_t r) const
{
  Int_t idx = -1;
  switch (d) { 
  case 1: idx = 0; break;
  case 2: idx = 1 + (r == 'I' || r == 'i' ? 0 : 1); break;
  case 3: idx = 3 + (r == 'I' || r == 'i' ? 0 : 1); break;
  }
  return idx;
}
//____________________________________________________________________
void
AliFMDHistCollector::GetDetRing(Int_t idx, UShort_t& d, Char_t& r) const
{
  d = 0; 
  r = '\0';
  switch (idx) { 
  case 0: d = 1; r = 'I'; break;
  case 1: d = 2; r = 'I'; break;
  case 2: d = 2; r = 'O'; break;
  case 3: d = 3; r = 'I'; break;
  case 4: d = 3; r = 'O'; break;
  }
}

//____________________________________________________________________
void
AliFMDHistCollector::GetFirstAndLast(Int_t idx, UShort_t vtxbin, 
				     Int_t& first, Int_t& last) const
{
  first = 0; 
  last  = 0;
  
  if (idx    <  0) return;
  if (vtxbin <= 0) return;
  idx += (vtxbin-1) * 5;
      
  if (idx < 0 || idx >= fFirstBins.GetSize()) return;
  
  first = fFirstBins.At(idx)+fNCutBins;  
  last  = fLastBins.At(idx)-fNCutBins;
}

//____________________________________________________________________
Int_t
AliFMDHistCollector::GetFirst(Int_t idx, UShort_t v) const 
{
  Int_t f, l;
  GetFirstAndLast(idx,v,f,l);
  return f;
}


//____________________________________________________________________
Int_t
AliFMDHistCollector::GetLast(Int_t idx, UShort_t v) const 
{
  Int_t f, l;
  GetFirstAndLast(idx,v,f,l);
  return l;
}

//____________________________________________________________________
Int_t 
AliFMDHistCollector::GetOverlap(UShort_t d, Char_t r, 
				Int_t bin,  UShort_t vtxbin) const
{
  // Get the possibly overlapping histogram 
  Int_t other = -1;
  if (d == 1) {
    if (bin <= GetLast(2,'I',vtxbin)) other = GetIdx(2,'I');
  }
  else if (d == 2 && r == 'I') {
    if      (bin <= GetLast(2,  'O', vtxbin)) other = GetIdx(2, 'O');
    else if (bin >= GetFirst(1, 'I', vtxbin)) other = GetIdx(1, 'I');
  }
  else if (d == 2 && r == 'O') {
    if (bin >= GetFirst(2, 'I', vtxbin)) other = GetIdx(2,'I');
  }
  else if (d == 3 && r == 'O') {
    if (bin <= GetLast(3, 'I', vtxbin)) other = GetIdx(3, 'I');
  }
  else if (d == 3 && r == 'I') {
    if (bin >= GetFirst(3, 'O', vtxbin)) other = GetIdx(3, 'O');
  }
  // AliInfo(Form("FMD%d%c (%d) -> %d", d, r, GetIdx(d,r), other));
  return other;
}
//____________________________________________________________________
Int_t
AliFMDHistCollector::GetOverlap(Int_t idx, Int_t bin, UShort_t vtxbin) const
{
  UShort_t d = 0; 
  Char_t   r = '\0';
  GetDetRing(idx, d, r);
  if (d == 0 || r == '\0') return 0;

  return GetOverlap(d, r, bin, vtxbin);
}
  
  
//____________________________________________________________________
Bool_t
AliFMDHistCollector::Collect(AliForwardUtil::Histos& hists,
			     UShort_t                vtxbin, 
			     TH2D&                   out)
{
  for (UShort_t d=1; d<=3; d++) { 
    UShort_t nr = (d == 1 ? 1 : 2);
    for (UShort_t q=0; q<nr; q++) { 
      Char_t      r = (q == 0 ? 'I' : 'O');
      TH2D*       h = hists.Get(d,r);
      TH2D*       t = static_cast<TH2D*>(h->Clone(Form("FMD%d%c_tmp",d,r)));
      
      // Outer rings have better phi segmentation - rebin to same as inner. 
      if (q == 1) t->RebinY(2);

      Int_t first = 0;
      Int_t last  = 0;
      GetFirstAndLast(d, r, vtxbin, first, last);

      // Now update profile output 
      for (Int_t iEta = first; iEta <= last; iEta++) { 

	// Get the possibly overlapping histogram 
	Int_t overlap = GetOverlap(d,r,iEta,vtxbin);

	// Fill eta acceptance for this event into the phi underlow bin
	Float_t ooc      = out.GetBinContent(iEta,0);
	Float_t noc      = overlap >= 0? 0.5 : 1;
	out.SetBinContent(iEta, 0, ooc + noc);

	for (Int_t iPhi = 1; iPhi <= h->GetNbinsY(); iPhi++) { 
	  Double_t c = t->GetBinContent(iEta,iPhi);
	  Double_t e = t->GetBinError(iEta,iPhi);

	  // If there's no signal, continue 
	  // if (e <= 0) continue;
	  if (c <= 0 || e <= 0)     continue;
	  
	  // If there's no overlapping histogram (ring), then 
	  // fill in data and continue to the next phi bin 
	  if (overlap < 0) { 
	    out.SetBinContent(iEta,iPhi,c);
	    out.SetBinError(iEta,iPhi,e);
	    continue;
	  }

	  // Get the current bin content and error 
	  Double_t oc = out.GetBinContent(iEta,iPhi);
	  Double_t oe = out.GetBinError(iEta,iPhi);

#define USE_STRAIGHT_MEAN
// #define USE_STRAIGHT_MEAN_NONZERO
// #define USE_WEIGHTED_MEAN
// #define USE_MOST_CERTAIN
#if defined(USE_STRAIGHT_MEAN)
	  // calculate the average of old value (half the original), 
	  // and this value, as well as the summed squared errors 
	  // of the existing content (sqrt((e_1/2)^2=sqrt(e_1^2/4)=e_1/2) 
	  // and half the error of this.   
	  //
	  // So, on the first overlapping histogram we get 
	  // 
	  //    c = c_1 / 2
	  //    e = sqrt((e_1 / 2)^2) = e_1/2
	  // 
	  // On the second we get 
	  // 
	  //    c' = c_2 / 2 + c = c_2 / 2 + c_1 / 2 = (c_1+c_2)/2 
	  //    e' = sqrt(e^2 + (e_2/2)^2) 
	  //       = sqrt(e_1^2/4 + e_2^2/4) 
	  //       = sqrt(1/4 * (e_1^2+e_2^2)) 
	  //       = 1/2 * sqrt(e_1^2 + e_2^2)
	  out.SetBinContent(iEta,iPhi,oc + c/2);
	  out.SetBinError(iEta,iPhi,TMath::Sqrt(oe*oe+(e*e)/4));
#elif defined(USE_STRAIGHT_MEAN_NONZERO)
# define ZERO_OTHER
	  // If there's data in the overlapping histogram, 
	  // calculate the average and add the errors in 
	  // quadrature.  
	  // 
	  // If there's no data in the overlapping histogram, 
	  // then just fill in the data 
	  if (oe > 0) {
	    out.SetBinContent(iEta,iPhi,(oc + c)/2);
	    out.SetBinError(iEta,iPhi,TMath::Sqrt(oe*oe + e*e)/2);
	  }
	  else {
	    out.SetBinContent(iEta,iPhi,c);
	    out.SetBinError(iEta,iPhi,e);
	  }	    
#elif defined(USE_WEIGHTED_MEAN) 
	  // Calculate the weighted mean 
	  Double_t w  = 1/(e*e);
	  Double_t sc = w * c;
	  Double_t sw = w;
	  if (oe > 0) {
	    Double_t ow = 1/(oe*oe);
	    sc          += ow * oc;
	    sw          += ow;
	  }
	  Double_t nc = sc / sw;
	  Double_t ne = TMath::Sqrt(1 / sw);
	  out.SetBinContent(iEta,iPhi,nc,ne);
#elif defined(USE_MOST_CERTAIN)
# define ZERO_OTHER
	  if (e < oe) {
	    out.SetBinContent(iEta,iPhi,c);
	    out.SetBinError(iEta,iPhi,e);
	  }
	  else {
	    out.SetBinContent(iEta,iPhi,oc);
	    out.SetBinError(iEta,iPhi,oe);
	  }
#else 
#         error No method for defining content of overlapping bins defined
#endif
#if defined(ZERO_OTHER)
	  // Get the content of the overlapping histogram, 
	  // and zero the content so that we won't use it 
	  // again 
	  UShort_t od; Char_t oR; 
	  GetDetRing(overlap, od, oR);
	  TH2D* other = hists.Get(od,oR);
	  other->SetBinContent(iEta,iPhi,0);
	  other->SetBinError(iEta,iPhi,0);
#endif
	}
      }
      // Remove temporary histogram 
      delete t;
    } // for r
  } // for d 
  return kTRUE;
}

//____________________________________________________________________
void
AliFMDHistCollector::Print(Option_t* /* option */) const
{
  char ind[gROOT->GetDirLevel()+1];
  for (Int_t i = 0; i < gROOT->GetDirLevel(); i++) ind[i] = ' ';
  ind[gROOT->GetDirLevel()] = '\0';
  std::cout << ind << "AliFMDHistCollector: " << GetName() << '\n'
	    << ind << " # of cut bins:          " << fNCutBins << '\n'
	    << ind << " Correction cut:         " << fCorrectionCut << '\n'
	    << ind << " Bin ranges:\n" << ind << "  v_z bin";
  Int_t nVz = fFirstBins.fN / 5;
  for (UShort_t iVz = 1; iVz <= nVz; iVz++) 
    std::cout << " | " << std::setw(7) << iVz;
  std::cout << '\n' << ind << "  / ring ";
  for (UShort_t iVz = 1; iVz <= nVz; iVz++) 
    std::cout << "-+--------";
  std::cout << std::endl;
    
  for (Int_t iIdx = 0; iIdx < 5; iIdx++) {
    UShort_t d = 0;
    Char_t   r = 0;
    GetDetRing(iIdx, d, r);
    
    std::cout << ind << "  FMD" << d << r << "  ";
    for (UShort_t iVz = 1; iVz <= nVz; iVz++) {
      Int_t first, last;
      GetFirstAndLast(iIdx, iVz, first, last);
      std::cout << " | " << std::setw(3) << first << "-" 
		<< std::setw(3) << last;
    }
    std::cout << std::endl;
  }
}

//____________________________________________________________________
//
// EOF
//
	  


