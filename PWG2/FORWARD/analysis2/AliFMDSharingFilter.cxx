//
// Class to do the sharing correction.  That is, a filter that merges 
// adjacent strip signals presumably originating from a single particle 
// that impinges on the detector in such a way that it deposite energy 
// into two or more strips. 
//
// Input: 
//    - AliESDFMD object  - from reconstruction
//
// Output: 
//    - AliESDFMD object  - copy of input, but with signals merged 
//
// Corrections used: 
//    - AliFMDCorrELossFit
//
// Histograms: 
//    - For each ring (FMD1i, FMD2i, FMD2o, FMD3i, FMD3o) the distribution of 
//      signals before and after the filter.  
//    - For each ring (see above), an array of distributions of number of 
//      hit strips for each vertex bin (if enabled - see Init method)
// 
//
//
#include "AliFMDSharingFilter.h"
#include <AliESDFMD.h>
#include <TAxis.h>
#include <TList.h>
#include <TH1.h>
#include <TMath.h>
#include "AliForwardCorrectionManager.h"
#include "AliFMDCorrELossFit.h"
#include <AliLog.h>
#include <TROOT.h>
#include <THStack.h>
#include <iostream>
#include <iomanip>

ClassImp(AliFMDSharingFilter)
#if 0
; // This is for Emacs
#endif 

#define DBG(L,M) \
  do { if (L>fDebug)break; std::cout << (M) << std::flush;} while(false) 
#define DBGL(L,M) \
  do { if (L>fDebug)break; std::cout << (M) << std::endl;} while(false) 
			      


//____________________________________________________________________
AliFMDSharingFilter::AliFMDSharingFilter()
  : TNamed(), 
    fRingHistos(),
    fLowCut(0.),
    fCorrectAngles(kFALSE), 
    fNXi(1),
    fIncludeSigma(true),
    fSummed(0),
    fHighCuts(0),
    fOper(0),
    fDebug(0),
    fZeroSharedHitsBelowThreshold(0)
{
  // 
  // Default Constructor - do not use 
  //
}

//____________________________________________________________________
AliFMDSharingFilter::AliFMDSharingFilter(const char* title)
  : TNamed("fmdSharingFilter", title), 
    fRingHistos(), 
    fLowCut(0.),
    fCorrectAngles(kFALSE), 
    fNXi(1),
    fIncludeSigma(true),
    fSummed(0),
    fHighCuts(0),
    fOper(0),
    fDebug(0),
    fZeroSharedHitsBelowThreshold(0)
{
  // 
  // Constructor 
  // 
  // Parameters:
  //    title Title of object  - not significant 
  //
  fRingHistos.Add(new RingHistos(1, 'I'));
  fRingHistos.Add(new RingHistos(2, 'I'));
  fRingHistos.Add(new RingHistos(2, 'O'));
  fRingHistos.Add(new RingHistos(3, 'I'));
  fRingHistos.Add(new RingHistos(3, 'O'));
}

//____________________________________________________________________
AliFMDSharingFilter::AliFMDSharingFilter(const AliFMDSharingFilter& o)
  : TNamed(o), 
    fRingHistos(), 
    fLowCut(o.fLowCut),
    fCorrectAngles(o.fCorrectAngles), 
    fNXi(o.fNXi),
    fIncludeSigma(o.fIncludeSigma),
    fSummed(o.fSummed),
    fHighCuts(o.fHighCuts),
    fOper(o.fOper),
    fDebug(o.fDebug),
    fZeroSharedHitsBelowThreshold(o.fZeroSharedHitsBelowThreshold)
{
  // 
  // Copy constructor 
  // 
  // Parameters:
  //    o Object to copy from 
  //
  TIter    next(&o.fRingHistos);
  TObject* obj = 0;
  while ((obj = next())) fRingHistos.Add(obj);
}

//____________________________________________________________________
AliFMDSharingFilter::~AliFMDSharingFilter()
{
  // 
  // Destructor
  //
  fRingHistos.Delete();
}

//____________________________________________________________________
AliFMDSharingFilter&
AliFMDSharingFilter::operator=(const AliFMDSharingFilter& o)
{
  // 
  // Assignment operator 
  // 
  // Parameters:
  //    o Object to assign from 
  // 
  // Return:
  //    Reference to this 
  //
  TNamed::operator=(o);

  fLowCut                       = o.fLowCut;
  fCorrectAngles                = o.fCorrectAngles;
  fNXi                          = o.fNXi;
  fDebug                        = o.fDebug;
  fOper                         = o.fOper;
  fSummed                       = o.fSummed;
  fHighCuts                     = o.fHighCuts;
  fIncludeSigma                 = o.fIncludeSigma;
  fZeroSharedHitsBelowThreshold = o.fZeroSharedHitsBelowThreshold;

  fRingHistos.Delete();
  TIter    next(&o.fRingHistos);
  TObject* obj = 0;
  while ((obj = next())) fRingHistos.Add(obj);
  
  return *this;
}

//____________________________________________________________________
AliFMDSharingFilter::RingHistos*
AliFMDSharingFilter::GetRingHistos(UShort_t d, Char_t r) const
{
  // 
  // Get the ring histogram container 
  // 
  // Parameters:
  //    d Detector
  //    r Ring 
  // 
  // Return:
  //    Ring histogram container 
  //
  Int_t idx = -1;
  switch (d) { 
  case 1: idx = 0; break;
  case 2: idx = 1 + (r == 'I' || r == 'i' ? 0 : 1); break;
  case 3: idx = 3 + (r == 'I' || r == 'i' ? 0 : 1); break;
  }
  if (idx < 0 || idx >= fRingHistos.GetEntries()) return 0;
  
  return static_cast<RingHistos*>(fRingHistos.At(idx));
}

//____________________________________________________________________
void
AliFMDSharingFilter::Init()
{
  // Initialise 
  AliForwardCorrectionManager&  fcm = AliForwardCorrectionManager::Instance();
 
  // Get the high cut.  The high cut is defined as the 
  // most-probably-value peak found from the energy distributions, minus 
  // 2 times the width of the corresponding Landau.
  AliFMDCorrELossFit* fits = fcm.GetELossFit();
  const TAxis& eAxis = fits->GetEtaAxis();

  UShort_t nEta = eAxis.GetNbins();
  fHighCuts->SetBins(nEta, eAxis.GetXmin(), eAxis.GetXmax(), 5, .5, 5.5);
  fHighCuts->GetYaxis()->SetBinLabel(1, "FMD1i");
  fHighCuts->GetYaxis()->SetBinLabel(2, "FMD2i");
  fHighCuts->GetYaxis()->SetBinLabel(3, "FMD2o");
  fHighCuts->GetYaxis()->SetBinLabel(4, "FMD3i");
  fHighCuts->GetYaxis()->SetBinLabel(5, "FMD3o");

  UShort_t ybin = 0;
  for (UShort_t d = 1; d <= 3; d++) {
    UShort_t nr = (d == 1 ? 1 : 2);
    for (UShort_t q = 0; q < nr; q++) { 
      Char_t r = (q == 0 ? 'I' : 'O');
      ybin++;
      for (UShort_t e = 1; e <= nEta; e++) { 
	Double_t eta = eAxis.GetBinCenter(e);
	Double_t cut = GetHighCut(d, r, eta, false);
	if (cut <= 0) continue;
	fHighCuts->SetBinContent(e, ybin, cut);
      }
    }
  }
}

//____________________________________________________________________
Bool_t
AliFMDSharingFilter::Filter(const AliESDFMD& input, 
			    Bool_t           lowFlux,
			    AliESDFMD&       output)
{
  // 
  // Filter the input AliESDFMD object
  // 
  // Parameters:
  //    input     Input 
  //    lowFlux   If this is a low-flux event 
  //    output    Output AliESDFMD object 
  // 
  // Return:
  //    True on success, false otherwise 
  //
  output.Clear();
  TIter    next(&fRingHistos);
  RingHistos* o      = 0;
  while ((o = static_cast<RingHistos*>(next())))
    o->Clear();

  if (fOper) fOper->Reset(0);
  Int_t nNone      = 0;
  Int_t nCandidate = 0;
  Int_t nMerged    = 0;
  Int_t nSummed    = 0;

  Status status[512];
  
  for(UShort_t d = 1; d <= 3; d++) {
    Int_t nRings = (d == 1 ? 1 : 2);
    for (UShort_t q = 0; q < nRings; q++) {
      Char_t      r      = (q == 0 ? 'I' : 'O');
      UShort_t    nsec   = (q == 0 ?  20 :  40);
      UShort_t    nstr   = (q == 0 ? 512 : 256);
      RingHistos* histos = GetRingHistos(d, r);
      
      for(UShort_t s = 0; s < nsec;  s++) {
	for (UShort_t t = 0; t < nstr; t++) status[t] = kCandidate;
#ifdef USE_OLDER_MERGING
	Bool_t usedThis   = kFALSE;
	Bool_t usedPrev   = kFALSE;
#endif	
	for(UShort_t t = 0; t < nstr; t++) {
	  output.SetMultiplicity(d,r,s,t,0.);
	  Float_t mult = SignalInStrip(input,d,r,s,t);
	  // Get the pseudo-rapidity 
	  Double_t eta = input.Eta(d,r,s,t);
	  Double_t phi = input.Phi(d,r,s,t) * TMath::Pi() / 180.;
	  if (s == 0) output.SetEta(d,r,s,t,eta);
	  
	  // Keep dead-channel information. 
	  if(mult == AliESDFMD::kInvalidMult)
	    output.SetMultiplicity(d,r,s,t,AliESDFMD::kInvalidMult);

	  // If no signal or dead strip, go on. 
	  if (mult == AliESDFMD::kInvalidMult || mult == 0) {
	    if (mult == 0) histos->fSum->Fill(eta,phi,mult);
	    status[t] = kNone;
	    continue;
	  }

	  // Fill the diagnostics histogram 
	  histos->fBefore->Fill(mult);

	  // Get next and previous signal - if any 
	  Double_t prevE = 0;
	  Double_t nextE = 0;
	  Status   prevStatus = (t == 0      ? kNone : status[t-1]);
	  Status   thisStatus = status[t];
	  Status   nextStatus = (t == nstr-1 ? kNone : status[t+1]);
	  if (t != 0) {
	    prevE = SignalInStrip(input,d,r,s,t-1);
	    if (prevE == AliESDFMD::kInvalidMult) prevE = 0;
	  }
	  if (t != nstr - 1) {
	    nextE = SignalInStrip(input,d,r,s,t+1);
	    if (nextE == AliESDFMD::kInvalidMult) nextE = 0;
	  }
	  if (t != 0) histos->fNeighborsBefore->Fill(prevE, mult);
	  
#ifdef USE_OLDER_MERGING
	  Double_t mergedEnergy = MultiplicityOfStrip(mult,eta,prevE,nextE,
						      lowFlux,d,r,s,t, 
						      usedPrev,usedThis);
	  status[t] = (usedPrev ? kMergedWithOther : kNone);
	  if (t != nstr - 1) status[t] = (usedThis ? kMergedWithOther : kNone);
#else 
	  Double_t mergedEnergy = MultiplicityOfStrip(mult, prevE, nextE, 
						      eta, lowFlux, 
						      d, r, s, t, 
						      prevStatus, 
						      thisStatus, 
						      nextStatus);
	  if (t != 0)      status[t-1] = prevStatus;
	  if (t != nstr-1) status[t+1] = nextStatus;
	  status[t] = thisStatus;

#endif
	  // If we're processing on non-angle corrected data, we
	  // should do the angle correction here
	  if (!fCorrectAngles)
	    mergedEnergy = AngleCorrect(mergedEnergy, eta);
	  if (mergedEnergy > 0) histos->Incr();

	  if (t != 0) 
	    histos->fNeighborsAfter->Fill(output.Multiplicity(d,r,s,t-1), 
					  mergedEnergy);
	  histos->fBeforeAfter->Fill(mult, mergedEnergy);
	  histos->fAfter->Fill(mergedEnergy);
	  histos->fSum->Fill(eta,phi,mergedEnergy);

	  output.SetMultiplicity(d,r,s,t,mergedEnergy);
	} // for strip
	for (UShort_t t = 0; t < nstr; t++) {
	  if (fOper) fOper->operator()(d, r, s, t) = status[t];
	  switch (status[t]) { 
	  case kNone:            nNone++;      break;
	  case kCandidate:       nCandidate++; break;
	  case kMergedWithOther: nMerged++;    break;
	  case kMergedInto:      nSummed++;    break;
	  }
	}
      } // for sector
    } // for ring 
  } // for detector
  fSummed->Fill(kNone,            nNone);
  fSummed->Fill(kCandidate,       nCandidate);
  fSummed->Fill(kMergedWithOther, nMerged);
  fSummed->Fill(kMergedInto,      nSummed);

  DBGL(1, Form("none=%9d, candidate=%9d, merged=%9d, summed=%9d", 
	       nNone, nCandidate, nMerged, nSummed));
  next.Reset();
  while ((o = static_cast<RingHistos*>(next())))
    o->Finish();

  return kTRUE;
}

//_____________________________________________________________________
Double_t 
AliFMDSharingFilter::SignalInStrip(const AliESDFMD& input, 
				   UShort_t         d,
				   Char_t           r,
				   UShort_t         s,
				   UShort_t         t) const
{
  // 
  // Get the signal in a strip 
  // 
  // Parameters:
  //    fmd   ESD object
  //    d     Detector
  //    r     Ring 
  //    s     Sector 
  //    t     Strip
  // 
  // Return:
  //    The energy signal 
  //
  Double_t mult = input.Multiplicity(d,r,s,t);
  // In case of 
  //  - bad value (invalid or 0) 
  //  - we want angle corrected and data is 
  //  - we don't want angle corrected and data isn't 
  // just return read value  
  if (mult == AliESDFMD::kInvalidMult               || 
      mult == 0                                     ||
      (fCorrectAngles && input.IsAngleCorrected()) || 
      (!fCorrectAngles && !input.IsAngleCorrected()))
    return mult;

  // If we want angle corrected data, correct it, 
  // otherwise de-correct it 
  if (fCorrectAngles) mult = AngleCorrect(mult, input.Eta(d,r,s,t));
  else                mult = DeAngleCorrect(mult, input.Eta(d,r,s,t));
  return mult;
}
//_____________________________________________________________________
Double_t 
AliFMDSharingFilter::GetLowCut(UShort_t, Char_t, Double_t) const
{
  //
  // Get the low cut.  Normally, the low cut is taken to be the lower
  // value of the fit range used when generating the energy loss fits.
  // However, if fLowCut is set (using SetLowCit) to a value greater
  // than 0, then that value is used.
  //
  if (fLowCut > 0) return fLowCut;
  AliForwardCorrectionManager&  fcm = AliForwardCorrectionManager::Instance();
  AliFMDCorrELossFit* fits = fcm.GetELossFit();
  return fits->GetLowCut();
}
			
//_____________________________________________________________________
Double_t 
AliFMDSharingFilter::GetHighCut(UShort_t d, Char_t r, 
				Double_t eta, Bool_t errors) const
{
  //
  // Get the high cut.  The high cut is defined as the 
  // most-probably-value peak found from the energy distributions, minus 
  // 2 times the width of the corresponding Landau.
  //
  AliForwardCorrectionManager&  fcm = AliForwardCorrectionManager::Instance();

 
  // Get the high cut.  The high cut is defined as the 
  // most-probably-value peak found from the energy distributions, minus 
  // 2 times the width of the corresponding Landau.
  AliFMDCorrELossFit* fits = fcm.GetELossFit();
  
  return fits->GetLowerBound(d, r, eta, fNXi, errors, fIncludeSigma);
}

//_____________________________________________________________________
namespace { 
  const char* status2String(AliFMDSharingFilter::Status s)
  {
    switch(s) { 
    case AliFMDSharingFilter::kNone:            return "none     "; 
    case AliFMDSharingFilter::kCandidate:       return "candidate"; 
    case AliFMDSharingFilter::kMergedWithOther: return "merged   "; 
    case AliFMDSharingFilter::kMergedInto:      return "summed   "; 
    }
    return "unknown  ";
  }
}

//_____________________________________________________________________
Double_t 
AliFMDSharingFilter::MultiplicityOfStrip(Double_t thisE,
					 Double_t prevE,
					 Double_t nextE,
					 Double_t eta,
					 Bool_t   lowFlux,
					 UShort_t d,
					 Char_t   r,
					 UShort_t s,
					 UShort_t t,
					 Status&  prevStatus, 
					 Status&  thisStatus, 
					 Status&  nextStatus) const
{
  Double_t lowCut = GetLowCut(d, r, eta);

  DBG(3,Form("FMD%d%c[%2d,%3d]: this=%9f(%s), prev=%9f(%s), next=%9f(%s)", 
	     d, r, s, t,
	     thisE, status2String(thisStatus), 
	     prevE, status2String(prevStatus), 
	     nextE, status2String(nextStatus)));

  // If below cut, then modify zero signal and make sure the next
  // strip is considered a candidate.
  if (thisE < lowCut || thisE > 20) { 
    thisStatus = kNone;
    DBGL(3,Form(" %9f<%9f || %9f>20, 0'ed", thisE, lowCut, thisE));
    if (prevStatus == kCandidate) prevStatus = kNone;
    return 0;  
  }
  // It this strip was merged with the previous strip, then 
  // make the next strip a candidate and zero the value in this strip. 
  if (thisStatus == kMergedWithOther) { 
    DBGL(3,Form(" Merged with other, 0'ed"));
    return 0;
  }

  // Get the upper cut on the sharing 
  Double_t highCut = GetHighCut(d, r, eta ,false);
  if (highCut <= 0) {
    thisStatus = kNone;
    return 0;
  }

  // Only for low-flux  events 
  if (lowFlux) { 
    // If this is smaller than the next, and the next is larger 
    // then the cut, mark both as candidates for sharing. 
    // They should be be merged on the next iteration 
    if (thisE < nextE && nextE > highCut) { 
      thisStatus = kCandidate;
      nextStatus = kCandidate;
      return 0;
    }
  }
  
  // Variable to sum signal in 
  Double_t totalE = thisE;

  // If the previous signal was between the two cuts, and is still
  // marked as a candidate , then merge it in here.
  if (prevStatus == kCandidate && prevE > lowCut && prevE < highCut) {
    totalE     += prevE;
    prevStatus =  kMergedWithOther;
    DBG(3, Form(" Prev candidate %9f<%9f<%9f->%9f", 
		lowCut, prevE, highCut, totalE));
  }

  // If the next signal is between the two cuts, then merge it here
  if (nextE > lowCut && nextE < highCut) { 
    totalE     += nextE;
    nextStatus =  kMergedWithOther;
    DBG(3, Form(" Next %9f<%9f<%9f->%9f", lowCut, nextE, highCut,totalE));
  }
  
  // If the signal is bigger than the high-cut, then return the value 
  if (totalE > highCut) { 
    if (totalE > thisE) 
      thisStatus = kMergedInto;
    else 
      thisStatus = kNone;
    DBGL(3, Form(" %9f>%f9 (was %9f) -> %9f %s", 
		 totalE, highCut, thisE, totalE,status2String(thisStatus)));
    return totalE;
  }
  else {
    // This is below cut, so flag next as a candidate 
    DBG(3, Form(" %9f<=%9f, next candidate", totalE, highCut));
    nextStatus = kCandidate;
  }
  
  // If the total signal is smaller than low cut then zero this and kill this 
  if (totalE < lowCut)  {
    DBGL(3, Form(" %9f<%9f (was %9f), zeroed", totalE, lowCut, thisE));
    thisStatus = kNone;
    return 0;
  }

  // If total signal not above high cut or lower than low cut, 
  // mark this as a candidate for merging into the next, and zero signal 
  DBGL(3, Form(" %9f<%9f<%9f (was %9f), zeroed, candidate", 
		   lowCut, totalE, highCut, thisE));
  thisStatus = kCandidate;
  if(fZeroSharedHitsBelowThreshold)
    return 0;
  else return totalE; 
}
	   
//_____________________________________________________________________
Double_t 
AliFMDSharingFilter::MultiplicityOfStrip(Double_t mult,
					 Double_t eta,
					 Double_t prevE,
					 Double_t nextE,
					 Bool_t   lowFlux,
					 UShort_t d,
					 Char_t   r,
					 UShort_t /*s*/,
					 UShort_t /*t*/,
					 Bool_t&  usedPrev, 
					 Bool_t&  usedThis) const
{
  // 
  // The actual algorithm 
  // 
  // Parameters:
  //    mult      The unfiltered signal in the strip
  //    eta       Psuedo rapidity 
  //    prevE     Previous strip signal (or 0)
  //    nextE     Next strip signal (or 0) 
  //    lowFlux   Whether this is a low flux event 
  //    d         Detector
  //    r         Ring 
  //    s         Sector 
  //    t         Strip
  //    usedPrev  Whether the previous strip was used in sharing or not
  //    usedThis  Wether this strip was used in sharing or not. 
  // 
  // Return:
  //    The filtered signal in the strip
  //

  // IF the multiplicity is very large, or below the cut, reject it, and 
  // flag it as candidate for sharing 
  Double_t    lowCut = GetLowCut(d,r,eta);
  if(mult > 12 || mult < lowCut) {
    usedThis      = kFALSE;
    usedPrev      = kFALSE;
    return 0;
  }

  // If this was shared with the previous signal, return 0 and mark it as 
  // not a candiate for sharing 
  if(usedThis) {
    usedThis      = kFALSE;
    usedPrev      = kTRUE;
    return 0.;
  }

  //analyse and perform sharing on one strip
  
  // Get the high cut 
  Double_t highCut = GetHighCut(d, r, eta, false);
  if (highCut <= 0) {
    usedThis = kFALSE;
    usedPrev = kTRUE;
    return 0;
  }

  // If this signal is smaller than the next, and the next signal is smaller 
  // than then the high cut, and it's a low-flux event, then mark this and 
  // the next signal as candidates for sharing 
  // 
  // This is the test if the signal is the smaller of two possibly
  // shared signals.   
  if (mult  < nextE   && 
      nextE > highCut && 
      lowFlux) {
    usedThis      = kFALSE;
    usedPrev      = kFALSE;
    return 0;
  }
  
  Double_t totalE  = mult;
  
  // If the previous signal was larger than the low cut, and 
  // the previous was smaller than high cut, and the previous signal 
  // isn't marked as used, then add it's energy to this signal 
  if (prevE > lowCut && 
      prevE < highCut && 
      !usedPrev) 
    totalE += prevE;

  // If the next signal is larger than the low cut, and 
  // the next signal is smaller than the high cut, then add the next signal 
  // to this, and marked the next signal as used 
  if(nextE > lowCut && nextE < highCut ) {
    totalE      += nextE;
    usedThis =  kTRUE;
  }

  // If we're processing on non-angle corrected data, we should do the 
  // angle correction here 
  // if (!fCorrectAngles) 
  //   totalE = AngleCorrect(totalE, eta);

  // Fill post processing histogram
  // if(totalE > fLowCut)
  //   GetRingHistos(d,r)->fAfter->Fill(totalE);

  // Return value 
  Double_t mergedEnergy = 0;
  
  if(totalE > 0) {
    // If we have a signal, then this is used (sharedPrev=true) and
    // the signal is set to the result
    mergedEnergy = totalE;
    usedPrev     = kTRUE;
  }
  else {
    // If we do not have a signal (too low), then this is not shared, 
    // and the next strip is not shared either 
    usedThis  = kFALSE;
    usedPrev  = kFALSE;
  }
  
  return mergedEnergy; 
}
//____________________________________________________________________
Double_t
AliFMDSharingFilter::AngleCorrect(Double_t mult, Double_t eta) const
{
  // 
  // Angle correct the signal 
  // 
  // Parameters:
  //    mult Angle Un-corrected Signal 
  //    eta  Pseudo-rapidity 
  // 
  // Return:
  //    Angle corrected signal 
  //
  Double_t theta =  2 * TMath::ATan(TMath::Exp(-eta));
  if (eta < 0) theta -= TMath::Pi();
  return mult * TMath::Cos(theta);
}
//____________________________________________________________________
Double_t
AliFMDSharingFilter::DeAngleCorrect(Double_t mult, Double_t eta) const
{
  // 
  // Angle de-correct the signal 
  // 
  // Parameters:
  //    mult Angle corrected Signal 
  //    eta  Pseudo-rapidity 
  // 
  // Return:
  //    Angle un-corrected signal 
  //
  Double_t theta =  2 * TMath::ATan(TMath::Exp(-eta));
  if (eta < 0) theta -= TMath::Pi();
  return mult / TMath::Cos(theta);
}

//____________________________________________________________________
void
AliFMDSharingFilter::ScaleHistograms(const TList* dir, Int_t nEvents)
{
  // 
  // Scale the histograms to the total number of events 
  // 
  // Parameters:
  //    dir     Where the output is 
  //    nEvents Number of events 
  //
  if (nEvents <= 0) return;
  TList* d = static_cast<TList*>(dir->FindObject(GetName()));
  if (!d) return;

  TIter    next(&fRingHistos);
  RingHistos* o = 0;
  THStack* sums = new THStack("sums", "Sum of ring signals");
  while ((o = static_cast<RingHistos*>(next()))) {
    o->ScaleHistograms(d, nEvents);
    TH1D* sum = o->fSum->ProjectionX(o->GetName(), 1, o->fSum->GetNbinsY(),"e");
    sum->Scale(1., "width");
    sum->SetTitle(o->GetName());
    sum->SetDirectory(0);
    sum->SetYTitle("#sum #Delta/#Delta_{mip}");
    sums->Add(sum);
  }
  d->Add(sums);
}

//____________________________________________________________________
void
AliFMDSharingFilter::DefineOutput(TList* dir)
{
  // 
  // Define the output histograms.  These are put in a sub list of the
  // passed list.   The histograms are merged before the parent task calls 
  // AliAnalysisTaskSE::Terminate 
  // 
  // Parameters:
  //    dir Directory to add to 
  //
  TList* d = new TList;
  d->SetName(GetName());
  dir->Add(d);

  fSummed = new TH2I("operations", "Strip operations", 
		     kMergedInto, kNone-.5, kMergedInto+.5,
		     51201, -.5, 51200.5);
  fSummed->SetXTitle("Operation");
  fSummed->SetYTitle("# of strips");
  fSummed->SetZTitle("Events");
  fSummed->GetXaxis()->SetBinLabel(kNone,            "None");
  fSummed->GetXaxis()->SetBinLabel(kCandidate,       "Candidate");
  fSummed->GetXaxis()->SetBinLabel(kMergedWithOther, "Merged w/other");
  fSummed->GetXaxis()->SetBinLabel(kMergedInto,      "Merged into");
  fSummed->SetDirectory(0);
  d->Add(fSummed);

  fHighCuts = new TH2D("highCuts", "High cuts used", 1,0,1, 1,0,1);
  fHighCuts->SetXTitle("#eta");
  fHighCuts->SetDirectory(0);
  d->Add(fHighCuts);


  TIter    next(&fRingHistos);
  RingHistos* o = 0;
  while ((o = static_cast<RingHistos*>(next()))) {
    o->Output(d);
  }
}
//____________________________________________________________________
void
AliFMDSharingFilter::Print(Option_t* /*option*/) const
{
  // 
  // Print information
  // 
  // Parameters:
  //    option Not used 
  //
  char ind[gROOT->GetDirLevel()+1];
  for (Int_t i = 0; i < gROOT->GetDirLevel(); i++) ind[i] = ' ';
  ind[gROOT->GetDirLevel()] = '\0';
  std::cout << ind << ClassName() << ": " << GetName() << '\n'
	    << std::boolalpha 
	    << ind << " Debug:                  " << fDebug << "\n"
	    << ind << " Low cut:                " << fLowCut << '\n'
	    << ind << " N xi factor:            " << fNXi    << '\n'
	    << ind << " Include sigma in cut:   " << fIncludeSigma << '\n'
	    << ind << " Use corrected angles:   " << fCorrectAngles 
	    << std::noboolalpha << std::endl;
}
  
//====================================================================
AliFMDSharingFilter::RingHistos::RingHistos()
  : AliForwardUtil::RingHistos(), 
    fBefore(0), 
    fAfter(0), 
    fBeforeAfter(0),
    fNeighborsBefore(0),
    fNeighborsAfter(0),
    fSum(0),
    fHits(0),
    fNHits(0)
{
  // 
  // Default CTOR
  //

}

//____________________________________________________________________
AliFMDSharingFilter::RingHistos::RingHistos(UShort_t d, Char_t r)
  : AliForwardUtil::RingHistos(d,r), 
    fBefore(0), 
    fAfter(0),
    fBeforeAfter(0),
    fNeighborsBefore(0),
    fNeighborsAfter(0),
    fSum(0),
    fHits(0),
    fNHits(0)
{
  // 
  // Constructor
  // 
  // Parameters:
  //    d detector
  //    r ring 
  //
  fBefore = new TH1D("esdEloss", "Energy loss (reconstruction)", 
		     1000, 0, 25);
  fBefore->SetXTitle("#Delta E/#Delta E_{mip}");
  fBefore->SetYTitle("P(#Delta E/#Delta E_{mip})");
  fBefore->SetFillColor(Color());
  fBefore->SetFillStyle(3001);
  fBefore->SetDirectory(0);

  fAfter  = static_cast<TH1D*>(fBefore->Clone("anaEloss"));
  fAfter->SetTitle("Energy loss in %s (sharing corrected)");
  fAfter->SetFillColor(Color()+2);
  fAfter->SetDirectory(0);

  Double_t max = 15;
  Double_t min = -1;
  Int_t    n   = int((max-min) / (max / 300));
  fBeforeAfter = new TH2D("beforeAfter", "Before and after correlation", 
			  n, min, max, n, min, max);
  fBeforeAfter->SetXTitle("#Delta E/#Delta E_{mip} before");
  fBeforeAfter->SetYTitle("#Delta E/#Delta E_{mip} after");
  fBeforeAfter->SetZTitle("Correlation");
  fBeforeAfter->SetDirectory(0);

  fNeighborsBefore = static_cast<TH2D*>(fBeforeAfter->Clone("neighborsBefore"));
  fNeighborsBefore->SetTitle("Correlation of neighbors befire");
  fNeighborsBefore->SetXTitle("#Delta E_{i}/#Delta E_{mip}");
  fNeighborsBefore->SetYTitle("#Delta E_{i+1}/#Delta E_{mip}");
  fNeighborsBefore->SetDirectory(0);
  
  fNeighborsAfter = 
    static_cast<TH2D*>(fNeighborsBefore->Clone("neighborsAfter"));
  fNeighborsAfter->SetTitle("Correlation of neighbors after");
  fNeighborsAfter->SetDirectory(0);

  fSum = new TH2D("summed", "Summed signal", 200, -4, 6, 
		  (fRing == 'I' || fRing == 'i' ? 20 : 40), 0, 2*TMath::Pi());
  fSum->SetDirectory(0);
  fSum->Sumw2();
  fSum->SetMarkerColor(Color());
  // fSum->SetFillColor(Color());
  fSum->SetXTitle("#eta");
  fSum->SetYTitle("#varphi [radians]");
  fSum->SetZTitle("#sum #Delta/#Delta_{mip}(#eta,#varphi) ");
  
  fHits = new TH1D("hits", "Number of hits", 200, 0, 200000);
  fHits->SetDirectory(0);
  fHits->SetXTitle("# of hits");
  fHits->SetFillColor(kGreen+1);
}
//____________________________________________________________________
AliFMDSharingFilter::RingHistos::RingHistos(const RingHistos& o)
  : AliForwardUtil::RingHistos(o), 
    fBefore(o.fBefore), 
    fAfter(o.fAfter),
    fBeforeAfter(o.fBeforeAfter),
    fNeighborsBefore(o.fNeighborsBefore),
    fNeighborsAfter(o.fNeighborsAfter),
    fSum(o.fSum),
    fHits(o.fHits),
    fNHits(o.fNHits)
{
  // 
  // Copy constructor 
  // 
  // Parameters:
  //    o Object to copy from 
  //
}

//____________________________________________________________________
AliFMDSharingFilter::RingHistos&
AliFMDSharingFilter::RingHistos::operator=(const RingHistos& o)
{
  // 
  // Assignment operator 
  // 
  // Parameters:
  //    o Object to assign from 
  // 
  // Return:
  //    Reference to this 
  //
  AliForwardUtil::RingHistos::operator=(o);
  fDet = o.fDet;
  fRing = o.fRing;
  
  if (fBefore) delete  fBefore;
  if (fAfter)  delete  fAfter;
  if (fHits)   delete fHits;
  
  fBefore          = static_cast<TH1D*>(o.fBefore->Clone());
  fAfter           = static_cast<TH1D*>(o.fAfter->Clone());
  fBeforeAfter     = static_cast<TH2D*>(o.fBeforeAfter->Clone());
  fNeighborsBefore = static_cast<TH2D*>(o.fNeighborsBefore->Clone());
  fNeighborsAfter  = static_cast<TH2D*>(o.fNeighborsAfter->Clone());
  fHits            = static_cast<TH1D*>(o.fHits->Clone());
  fSum             = static_cast<TH2D*>(o.fSum->Clone());

  return *this;
}
//____________________________________________________________________
AliFMDSharingFilter::RingHistos::~RingHistos()
{
  // 
  // Destructor 
  //
}
//____________________________________________________________________
void
AliFMDSharingFilter::RingHistos::Finish()
{
  // 
  // Finish off 
  // 
  //
  fHits->Fill(fNHits);
}

//____________________________________________________________________
void
AliFMDSharingFilter::RingHistos::ScaleHistograms(const TList* dir, 
						 Int_t nEvents)
{
  // 
  // Scale the histograms to the total number of events 
  // 
  // Parameters:
  //    nEvents Number of events 
  //    dir     Where the output is 
  //
  TList* l = GetOutputList(dir);
  if (!l) return; 

  TH1D* before = static_cast<TH1D*>(l->FindObject("esdELoss"));
  TH1D* after  = static_cast<TH1D*>(l->FindObject("anaELoss"));
  TH2D* summed = static_cast<TH2D*>(l->FindObject("summed"));
  if (before) before->Scale(1./nEvents);
  if (after)  after->Scale(1./nEvents);
  if (summed) summed->Scale(1./nEvents);
  
}

//____________________________________________________________________
void
AliFMDSharingFilter::RingHistos::Output(TList* dir)
{
  // 
  // Make output 
  // 
  // Parameters:
  //    dir where to store 
  //
  TList* d = DefineOutputList(dir);

  d->Add(fBefore);
  d->Add(fAfter);
  d->Add(fBeforeAfter);
  d->Add(fNeighborsBefore);
  d->Add(fNeighborsAfter);
  d->Add(fHits);
  d->Add(fSum);

  dir->Add(d);
}

//____________________________________________________________________
//
// EOF
//
