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
#include "AliFMDStripIndex.h"
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
#include <TParameter.h>
#include <iostream>
#include <iomanip>

ClassImp(AliFMDSharingFilter)
#if 0
; // This is for Emacs
#endif 

#define DBG(L,M)							\
  do { if (L>fDebug)break; std::cout << (M) << std::flush;} while(false) 
#define DBGL(L,M)							\
  do { if (L>fDebug)break; std::cout << (M) << std::endl;} while(false) 
			      


//____________________________________________________________________
AliFMDSharingFilter::AliFMDSharingFilter()
  : TNamed(), 
    fRingHistos(),
    fCorrectAngles(kFALSE), 
    // fSummed(0),
    fHighCuts(0),
    fLowCuts(0),
    // fOper(0),
    fDebug(0),
    fZeroSharedHitsBelowThreshold(false),
    fLCuts(),
    fHCuts(),
    fUseSimpleMerging(false),
    fThreeStripSharing(true),
    fMergingDisabled(false),
    fIgnoreESDForAngleCorrection(false)
{
  // 
  // Default Constructor - do not use 
  //
  DGUARD(fDebug,1, "Default CTOR for AliFMDSharingFilter");
}

//____________________________________________________________________
AliFMDSharingFilter::AliFMDSharingFilter(const char* title)
  : TNamed("fmdSharingFilter", title), 
    fRingHistos(), 
    fCorrectAngles(kFALSE), 
    // fSummed(0),
    fHighCuts(0),
    fLowCuts(0),
    // fOper(0),
    fDebug(0),
    fZeroSharedHitsBelowThreshold(false),
    fLCuts(),
    fHCuts(),
    fUseSimpleMerging(false),
    fThreeStripSharing(true),
    fMergingDisabled(false),
    fIgnoreESDForAngleCorrection(false)
{
  // 
  // Constructor 
  // 
  // Parameters:
  //    title Title of object  - not significant 
  //
  DGUARD(fDebug,1, "Named CTOR for AliFMDSharingFilter: %s", title);
  fRingHistos.SetName(GetName());
  fRingHistos.SetOwner();
  
  fHCuts.Set(AliFMDMultCuts::kLandauSigmaWidth, 1);
  fLCuts.Set(AliFMDMultCuts::kFixed, .15);

  // fExtraDead.Reset(-1);
}

//____________________________________________________________________
AliFMDSharingFilter::~AliFMDSharingFilter()
{
  // 
  // Destructor
  //
  DGUARD(fDebug,3, "DTOR for AliFMDSharingFilter");
  // fRingHistos.Delete();
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
AliFMDSharingFilter::SetupForData(const TAxis& axis)
{
  // Initialise - called on first event
  DGUARD(fDebug,1, "Initialize for AliFMDSharingFilter");
  AliForwardCorrectionManager& fcm  = AliForwardCorrectionManager::Instance();
  const AliFMDCorrELossFit*    fits = fcm.GetELossFit();

  // Get the high cut.  The high cut is defined as the 
  // most-probably-value peak found from the energy distributions, minus 
  // 2 times the width of the corresponding Landau.
  
  TAxis eAxis(axis.GetNbins(),
	      axis.GetXmin(),
	      axis.GetXmax());
  if(fits) 
    eAxis.Set(fits->GetEtaAxis().GetNbins(), 
	      fits->GetEtaAxis().GetXmin(),
	      fits->GetEtaAxis().GetXmax());

  UShort_t nEta = eAxis.GetNbins();
	
  fHighCuts->SetBins(nEta, eAxis.GetXmin(), eAxis.GetXmax(), 5, .5, 5.5);
  fHighCuts->GetYaxis()->SetBinLabel(1, "FMD1i");
  fHighCuts->GetYaxis()->SetBinLabel(2, "FMD2i");
  fHighCuts->GetYaxis()->SetBinLabel(3, "FMD2o");
  fHighCuts->GetYaxis()->SetBinLabel(4, "FMD3i");
  fHighCuts->GetYaxis()->SetBinLabel(5, "FMD3o");

  fLowCuts->SetBins(nEta, eAxis.GetXmin(), eAxis.GetXmax(), 5, .5, 5.5);
  fLowCuts->GetYaxis()->SetBinLabel(1, "FMD1i");
  fLowCuts->GetYaxis()->SetBinLabel(2, "FMD2i");
  fLowCuts->GetYaxis()->SetBinLabel(3, "FMD2o");
  fLowCuts->GetYaxis()->SetBinLabel(4, "FMD3i");
  fLowCuts->GetYaxis()->SetBinLabel(5, "FMD3o");

  // Cache our cuts in histograms 
  fLCuts.FillHistogram(fLowCuts);
  fHCuts.FillHistogram(fHighCuts);
}

//____________________________________________________________________
#define ETA2COS(ETA)						\
  TMath::Cos(2*TMath::ATan(TMath::Exp(-TMath::Abs(ETA))))

Bool_t
AliFMDSharingFilter::Filter(const AliESDFMD& input, 
			    Bool_t           /*lowFlux*/,
			    AliESDFMD&       output, 
			    Double_t         /*zvtx*/)
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
  DGUARD(fDebug,1, "Filter event in AliFMDSharingFilter");
  output.Clear();
  TIter    next(&fRingHistos);
  RingHistos* o      = 0;
  while ((o = static_cast<RingHistos*>(next()))) o->Clear();

  Int_t nSingle    = 0;
  Int_t nDouble    = 0;
  Int_t nTriple    = 0;

  for(UShort_t d = 1; d <= 3; d++) {
    Int_t nRings = (d == 1 ? 1 : 2);
    for (UShort_t q = 0; q < nRings; q++) {
      Char_t      r      = (q == 0 ? 'I' : 'O');
      UShort_t    nsec   = (q == 0 ?  20 :  40);
      UShort_t    nstr   = (q == 0 ? 512 : 256);
      RingHistos* histos = GetRingHistos(d, r);
      
      for(UShort_t s = 0; s < nsec;  s++) {	
	// `used' flags if the _current_ strip was used by _previous_ 
	// iteration. 
	Bool_t   used            = kFALSE;
	// `eTotal' contains the current sum of merged signals so far 
	Double_t eTotal          = -1;
	// Int_t    nDistanceBefore = -1;
	// Int_t    nDistanceAfter  = -1;
	// `twoLow' flags if we saw two consequtive strips with a 
	// signal between the two cuts. 
	Bool_t   twoLow          = kFALSE;
        Int_t    nStripsAboveCut = 0;
	
	for(UShort_t t = 0; t < nstr; t++) {
	  // nDistanceBefore++;
	  // nDistanceAfter++;

	  output.SetMultiplicity(d,r,s,t,0.);
	  Float_t mult         = SignalInStrip(input,d,r,s,t);
	  Float_t multNext     = (t<nstr-1) ? SignalInStrip(input,d,r,s,t+1) :0;
	  Float_t multNextNext = (t<nstr-2) ? SignalInStrip(input,d,r,s,t+2) :0;
	  if (multNext     ==  AliESDFMD::kInvalidMult) multNext     = 0;
	  if (multNextNext ==  AliESDFMD::kInvalidMult) multNextNext = 0;
	  if(!fThreeStripSharing) multNextNext = 0;

	  // Get the pseudo-rapidity 
	  Double_t eta = input.Eta(d,r,s,t);
	  Double_t phi = input.Phi(d,r,s,t) * TMath::Pi() / 180.;
	  if (s == 0) output.SetEta(d,r,s,t,eta);
	  
	  // Keep dead-channel information - either from the ESD (but
	  // see above for older data) or from the settings in the
	  // ForwardAODConfig.C file.
	  if (mult == AliESDFMD::kInvalidMult) {
	    output.SetMultiplicity(d,r,s,t,AliESDFMD::kInvalidMult);
	    histos->fBefore->Fill(-1);
	    mult = AliESDFMD::kInvalidMult;
	  }
	  
	  Double_t lowCut  = GetLowCut(d, r, eta);
	  Double_t highCut = GetHighCut(d, r, eta, false);
	  if (mult != AliESDFMD::kInvalidMult && mult > lowCut) {
	    // Always fill the ESD sum histogram 
	    histos->fSumESD->Fill(eta, phi, mult);
	  }

	  // If no signal or dead strip, go on. 
	  if (mult == AliESDFMD::kInvalidMult || mult == 0) {
	    if (mult == 0) histos->fSum->Fill(eta,phi,mult);
	    // Flush a possible signal 
	    if (eTotal > 0 && t > 0) 
	      output.SetMultiplicity(d,r,s,t-1,eTotal);
	    // Reset states so we do not try to merge over a dead strip. 
	    eTotal = -1;
	    used   = false;
	    twoLow = false;
	    if (t > 0)	
	      histos->fNConsecutive->Fill(nStripsAboveCut);
	    if (mult == AliESDFMD::kInvalidMult)
	      // Why not fill immidiately here? 
	      nStripsAboveCut = -1;
	    else
	      // Why not fill immidiately here? 
	      nStripsAboveCut = 0;	
	    continue;
	  }

	  // Fill the diagnostics histogram 
	  histos->fBefore->Fill(mult);

	  Double_t mergedEnergy = mult;
	  // it seems to me that this logic could be condensed a bit
          if(mult > lowCut) {		  
	    if(nStripsAboveCut < 1) {
	      if(t > 0)
		histos->fNConsecutive->Fill(nStripsAboveCut);
	      nStripsAboveCut=0;
	    }
	    nStripsAboveCut++;
	  }	
	  else {
	    if (t > 0)
	      histos->fNConsecutive->Fill(nStripsAboveCut);
	    nStripsAboveCut=0;
	  }		

	  if (!fMergingDisabled) {
	    mergedEnergy = 0;

	    // The current sum
	    Float_t etot = 0;
	  
	    // Fill in neighbor information
	    if (t < nstr-1) histos->fNeighborsBefore->Fill(mult,multNext);

	    Bool_t thisValid = mult     > lowCut;
	    Bool_t nextValid = multNext > lowCut;
	    Bool_t thisSmall = mult     < highCut;
	    Bool_t nextSmall = multNext < highCut;
	  
	    // If this strips signal is above the high cut, reset distance
	    // if (!thisSmall) {
	    //    histos->fDistanceBefore->Fill(nDistanceBefore);
	    //    nDistanceBefore = -1;
	    // }
	  
	    // If the total signal in the past 1 or 2 strips are non-zero
	    // we need to check 
	    if (eTotal > 0) {
	      // Here, we have already flagged one strip as a candidate 
	    
	      // If 3-strip merging is enabled, then check the next 
	      // strip to see that it falls within cut, or if we have 
	      // two low signals 
	      if (fThreeStripSharing && nextValid && (nextSmall || twoLow)) {
		eTotal = eTotal + multNext;
		used = kTRUE;
		histos->fTriple->Fill(eTotal);
		nTriple++;
		twoLow = kFALSE;
	      }
	      // Otherwise, we got a double hit before, and that 
	      // should be stored. 
	      else {
		used = kFALSE;
		histos->fDouble->Fill(eTotal);
		nDouble++;
	      }
	      // Store energy loss and reset sum 
	      etot   = eTotal;
	      eTotal = -1;
	    } // if (eTotal>0)
	    else {
	      // If we have no current sum 
	    
	      // Check if this is marked as used, and if so, continue
	      if (used) {used = kFALSE; continue; }
	    
	      // If the signal is abvoe the cut, set current
	      if (thisValid) etot = mult;
	    
	      // If the signal is abiove the cut, and so is the next 
	      // signal and either of them are below the high cut, 
	      if (thisValid  && nextValid  && (thisSmall || nextSmall)) {
	      
		// If this is below the high cut, and the next is too, then 
		// we have two low signals 
		if (thisSmall && nextSmall) twoLow = kTRUE;
	      
		// If this signal is bigger than the next, and the 
		// one after that is below the low-cut, then update 
		// the sum
		if (mult>multNext && multNextNext < lowCut) {
		  etot = mult + multNext;
		  used = kTRUE;
		  histos->fDouble->Fill(etot);
		  nDouble++;
		}
		// Otherwise, we may need to merge with a third strip
		else {
		  etot   = 0;
		  eTotal = mult + multNext;
		}
	      }
	      // This is a signle hit 
	      else if(etot > 0) {
		histos->fSingle->Fill(etot);
		histos->fSinglePerStrip->Fill(etot,t);
		nSingle++;
	      }
	    } // else if (etotal >= 0)
	  
	    mergedEnergy = etot;
	    // if (mergedEnergy > GetHighCut(d, r, eta ,false)) {
	    //   histos->fDistanceAfter->Fill(nDistanceAfter);
	    //   nDistanceAfter    = -1;
	    // }
	    //if(mult>0 && multNext >0)
	    //  std::cout<<mult<<"  "<<multNext<<"  "<<mergedEnergy<<std::endl;
	  } // if (!fMergingDisabled)

	  if (!fCorrectAngles)
	    mergedEnergy = AngleCorrect(mergedEnergy, eta);
	  // if (mergedEnergy > 0) histos->Incr();
	  
	  if (t != 0) 
	    histos->fNeighborsAfter->Fill(output.Multiplicity(d,r,s,t-1), 
					  mergedEnergy);
	  histos->fBeforeAfter->Fill(mult, mergedEnergy);
	  if(mergedEnergy > 0)
	    histos->fAfter->Fill(mergedEnergy);
	  histos->fSum->Fill(eta,phi,mergedEnergy);
	  
	  output.SetMultiplicity(d,r,s,t,mergedEnergy);
	} // for strip
	histos->fNConsecutive->Fill(nStripsAboveCut); // fill the last sector 
      } // for sector
    } // for ring 
  } // for detector
  DMSG(fDebug, 3,"single=%9d, double=%9d, triple=%9d", 
       nSingle, nDouble, nTriple);
  next.Reset();
  // while ((o = static_cast<RingHistos*>(next()))) o->Finish();

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
      (fCorrectAngles && (fIgnoreESDForAngleCorrection || input.IsAngleCorrected())) || 
      (!fCorrectAngles && !fIgnoreESDForAngleCorrection && !input.IsAngleCorrected()))
    return mult;

  // If we want angle corrected data, correct it, 
  // otherwise de-correct it 
  if (fCorrectAngles) mult = AngleCorrect(mult, input.Eta(d,r,s,t));
  else                mult = DeAngleCorrect(mult, input.Eta(d,r,s,t));
  return mult;
}

namespace {
  Double_t Rng2Cut(UShort_t d, Char_t r, Double_t eta, TH2* h) {
    Double_t ret = 1024;
    Int_t ybin = 0;							
    switch(d) {								
    case 1: ybin = 1; break;						
    case 2: ybin = (r=='i' || r=='I') ? 2 : 3; break;			
    case 3: ybin = (r=='i' || r=='I') ? 4 : 5; break;			
    default: return ret;
    }									
    Int_t xbin = h->GetXaxis()->FindBin(eta);				
    if (xbin < 1 && xbin > h->GetXaxis()->GetNbins()) return ret;
    ret = h->GetBinContent(xbin,ybin);					
    return ret;
  }
}
    
  //_____________________________________________________________________
Double_t 
AliFMDSharingFilter::GetLowCut(UShort_t d, Char_t r, Double_t eta) const
{
  //
  // Get the low cut.  Normally, the low cut is taken to be the lower
  // value of the fit range used when generating the energy loss fits.
  // However, if fLowCut is set (using SetLowCit) to a value greater
  // than 0, then that value is used.
  //
  return Rng2Cut(d, r, eta, fLowCuts);
  // return fLCuts.GetMultCut(d,r,eta,false);
}
			
//_____________________________________________________________________
Double_t 
AliFMDSharingFilter::GetHighCut(UShort_t d, Char_t r, 
				Double_t eta, Bool_t /*errors*/) const
{
  //
  // Get the high cut.  The high cut is defined as the 
  // most-probably-value peak found from the energy distributions, minus 
  // 2 times the width of the corresponding Landau.
  //
  return Rng2Cut(d, r, eta, fHighCuts);
  // return fHCuts.GetMultCut(d,r,eta,errors); 
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
AliFMDSharingFilter::Terminate(const TList* dir, TList* output, Int_t nEvents)
{
  // 
  // Scale the histograms to the total number of events 
  // 
  // Parameters:
  //    dir     Where the output is 
  //    nEvents Number of events 
  //
  DGUARD(fDebug,1, "Scale histograms in AliFMDSharingFilter");
  if (nEvents <= 0) return;
  TList* d = static_cast<TList*>(dir->FindObject(GetName()));
  if (!d) return;

  TList* out = new TList;
  out->SetName(d->GetName());
  out->SetOwner();

  TParameter<int>* nFiles = 
    static_cast<TParameter<int>*>(d->FindObject("nFiles"));

  TH2* lowCuts  = static_cast<TH2*>(d->FindObject("lowCuts"));
  TH2* highCuts = static_cast<TH2*>(d->FindObject("highCuts"));
  if (lowCuts && nFiles) {
    TH1* oh = static_cast<TH1*>(lowCuts->Clone());
    oh->Scale(1. / nFiles->GetVal());	
    oh->SetBit(BIT(20));
    out->Add(oh);
  }
  else 
    AliWarning("low cuts histogram not found in input list");
  if (highCuts && nFiles) {
    TH1* oh = static_cast<TH1*>(highCuts->Clone());
    oh->Scale(1. / nFiles->GetVal());
    oh->SetBit(BIT(20));
    out->Add(oh);
  }
  else 
    AliWarning("high cuts histogram not found in input list");
  
  fRingHistos.Add(new RingHistos(1, 'I'));
  fRingHistos.Add(new RingHistos(2, 'I'));
  fRingHistos.Add(new RingHistos(2, 'O'));
  fRingHistos.Add(new RingHistos(3, 'I'));
  fRingHistos.Add(new RingHistos(3, 'O'));

  TIter    next(&fRingHistos);
  RingHistos* o = 0;
  THStack* sums = new THStack("sums", "Sum of ring signals");
  THStack* sumsESD = new THStack("sumsESD", "Sum of ring ESD signals");
  while ((o = static_cast<RingHistos*>(next()))) {
    o->Terminate(d, nEvents);
    if (!o->fSum) { 
      Warning("Terminate", "No sum histogram found for ring %s", o->GetName());
      continue;
    }
    TH1D* sum = o->fSum->ProjectionX(o->GetName(), 1, o->fSum->GetNbinsY(),"e");
    sum->Scale(1., "width");
    sum->SetTitle(o->GetName());
    sum->SetDirectory(0);
    sum->SetYTitle("#sum_{c} #Delta/#Delta_{mip}");
    sums->Add(sum);


    if (o->fSumESD) { 
      sum = o->fSumESD->ProjectionX(o->GetName(), 1, o->fSumESD->GetNbinsY(),"e");
      sum->Scale(1., "width");
      sum->SetTitle(o->GetName());
      sum->SetDirectory(0);
      sum->SetYTitle("#sum_{s} #Delta/#Delta_{mip}");
      sumsESD->Add(sum);
    }
  }
  out->Add(sums);
  out->Add(sumsESD);
  output->Add(out);
}

//____________________________________________________________________
void
AliFMDSharingFilter::CreateOutputObjects(TList* dir)
{
  // 
  // Define the output histograms.  These are put in a sub list of the
  // passed list.   The histograms are merged before the parent task calls 
  // AliAnalysisTaskSE::Terminate 
  // 
  // Parameters:
  //    dir Directory to add to 
  //
  DGUARD(fDebug,1, "Define output in AliFMDSharingFilter");
  TList* d = new TList;
  d->SetOwner();
  d->SetName(GetName());
  dir->Add(d);

#if 0
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
#endif

  fHighCuts = new TH2D("highCuts", "High cuts used", 1,0,1, 1,0,1);
  fHighCuts->SetXTitle("#eta");
  fHighCuts->SetDirectory(0);
  d->Add(fHighCuts);

  fLowCuts = new TH2D("lowCuts", "Low cuts used", 1,0,1, 1,0,1);
  fLowCuts->SetXTitle("#eta");
  fLowCuts->SetDirectory(0);
  d->Add(fLowCuts);

  // d->Add(lowCut);
  // d->Add(nXi);
  // d->Add(sigma);
  d->Add(AliForwardUtil::MakeParameter("angle", fCorrectAngles));
  d->Add(AliForwardUtil::MakeParameter("lowSignal", 
				       fZeroSharedHitsBelowThreshold));
  d->Add(AliForwardUtil::MakeParameter("simple", fUseSimpleMerging));
  d->Add(AliForwardUtil::MakeParameter("sumThree", fThreeStripSharing));
  d->Add(AliForwardUtil::MakeParameter("disabled", fMergingDisabled));
  TParameter<int>* nFiles = new TParameter<int>("nFiles", 1);
  nFiles->SetMergeMode('+');
  d->Add(nFiles);

  fLCuts.Output(d,"lCuts");
  fHCuts.Output(d,"hCuts");

  fRingHistos.Add(new RingHistos(1, 'I'));
  fRingHistos.Add(new RingHistos(2, 'I'));
  fRingHistos.Add(new RingHistos(2, 'O'));
  fRingHistos.Add(new RingHistos(3, 'I'));
  fRingHistos.Add(new RingHistos(3, 'O'));

  TIter    next(&fRingHistos);
  RingHistos* o = 0;
  while ((o = static_cast<RingHistos*>(next()))) {
    o->CreateOutputObjects(d);
  }
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
AliFMDSharingFilter::Print(Option_t* /*option*/) const
{
  // 
  // Print information
  // 
  // Parameters:
  //    option Not used 
  //
  AliForwardUtil::PrintTask(*this);
  gROOT->IncreaseDirLevel();

  PFB("Use corrected angles",  fCorrectAngles);
  PFB("Zero below threshold",  fZeroSharedHitsBelowThreshold);
  PFB("Use simple sharing",    fUseSimpleMerging);
  PFB("Allow 3 strip merging", fThreeStripSharing);
  PF("Low cuts",	"");
  fLCuts.Print();
  PF("High cuts",	"");
  fHCuts.Print();
  gROOT->DecreaseDirLevel();
}
  
//====================================================================
AliFMDSharingFilter::RingHistos::RingHistos()
  : AliForwardUtil::RingHistos(), 
    fBefore(0), 
    fAfter(0), 
    fSingle(0),
    fDouble(0),
    fTriple(0),
    fSinglePerStrip(0),
    // fDistanceBefore(0),
    // fDistanceAfter(0),
    fBeforeAfter(0),
    fNeighborsBefore(0),
    fNeighborsAfter(0),
    fSumESD(0),
    fSum(0),
    fNConsecutive(0)	
     // ,
    // fHits(0),
    // fNHits(0)
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
    fSingle(0),
    fDouble(0),
    fTriple(0),    
    fSinglePerStrip(0),
    // fDistanceBefore(0),
    // fDistanceAfter(0),
    fBeforeAfter(0),
    fNeighborsBefore(0),
    fNeighborsAfter(0),
    fSumESD(0),
    fSum(0),
    fNConsecutive(0)	
     //,
    // fHits(0),
    // fNHits(0)
{
  // 
  // Constructor
  // 
  // Parameters:
  //    d detector
  //    r ring 
  //
  fBefore = new TH1D("esdEloss", Form("Energy loss in %s (reconstruction)", 
				      GetName()), 640, -1, 15);
  fBefore->SetXTitle("#Delta E/#Delta E_{mip}");
  fBefore->SetYTitle("P(#Delta E/#Delta E_{mip})");
  fBefore->SetFillColor(Color());
  fBefore->SetFillStyle(3001);
  fBefore->SetLineColor(kBlack);
  fBefore->SetLineStyle(2);
  fBefore->SetDirectory(0);

  fAfter  = static_cast<TH1D*>(fBefore->Clone("anaEloss"));
  fAfter->SetTitle(Form("Energy loss in %s (sharing corrected)", GetName()));
  fAfter->SetFillColor(Color()+2);
  fAfter->SetLineStyle(1);
  fAfter->SetDirectory(0);
  
  fSingle = new TH1D("singleEloss", "Energy loss (single strips)", 
		     600, 0, 15);
  fSingle->SetXTitle("#Delta/#Delta_{mip}");
  fSingle->SetYTitle("P(#Delta/#Delta_{mip})");
  fSingle->SetFillColor(Color());
  fSingle->SetFillStyle(3001);
  fSingle->SetLineColor(kBlack);
  fSingle->SetLineStyle(2);
  fSingle->SetDirectory(0);

  fDouble = static_cast<TH1D*>(fSingle->Clone("doubleEloss"));
  fDouble->SetTitle("Energy loss (two strips)");
  fDouble->SetFillColor(Color()+1);
  fDouble->SetDirectory(0);
  
  fTriple = static_cast<TH1D*>(fSingle->Clone("tripleEloss"));
  fTriple->SetTitle("Energy loss (three strips)"); 
  fTriple->SetFillColor(Color()+2);
  fTriple->SetDirectory(0);
  
  //Int_t nBinsForInner = (r == 'I' ? 32 : 16);
  Int_t nBinsForInner = (r == 'I' ? 512 : 256);
  Int_t nStrips       = (r == 'I' ? 512 : 256);
  
  fSinglePerStrip = new TH2D("singlePerStrip", "SinglePerStrip", 
			     600,0,15, nBinsForInner,0,nStrips);
  fSinglePerStrip->SetXTitle("#Delta/#Delta_{mip}");
  fSinglePerStrip->SetYTitle("Strip #");
  fSinglePerStrip->SetZTitle("Counts");
  fSinglePerStrip->SetDirectory(0);

#if 0
  fDistanceBefore = new TH1D("distanceBefore", "Distance before sharing", 
			     nStrips , 0,nStrips );
  fDistanceBefore->SetXTitle("Distance");
  fDistanceBefore->SetYTitle("Counts");
  fDistanceBefore->SetFillColor(kGreen+2);
  fDistanceBefore->SetFillStyle(3001);
  fDistanceBefore->SetLineColor(kBlack);
  fDistanceBefore->SetLineStyle(2);
  fDistanceBefore->SetDirectory(0);

  fDistanceAfter = static_cast<TH1D*>(fDistanceBefore->Clone("distanceAfter"));
  fDistanceAfter->SetTitle("Distance after sharing"); 
  fDistanceAfter->SetFillColor(kGreen+1);
  fDistanceAfter->SetDirectory(0);
#endif
  
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
  fNeighborsBefore->SetTitle("Correlation of neighbors before");
  fNeighborsBefore->SetXTitle("#Delta E_{i}/#Delta E_{mip}");
  fNeighborsBefore->SetYTitle("#Delta E_{i+1}/#Delta E_{mip}");
  fNeighborsBefore->SetDirectory(0);
  
  fNeighborsAfter = 
    static_cast<TH2D*>(fNeighborsBefore->Clone("neighborsAfter"));
  fNeighborsAfter->SetTitle("Correlation of neighbors after");
  fNeighborsAfter->SetDirectory(0);

  fSumESD = new TH2D("summedESD", "Summed ESD signal", 200, -4, 6, 
		  NSector(), 0, 2*TMath::Pi());
  fSumESD->SetDirectory(0);
  fSumESD->Sumw2();
  fSumESD->SetMarkerColor(Color());
  // fSum->SetFillColor(Color());
  fSumESD->SetXTitle("#eta");
  fSumESD->SetYTitle("#varphi [radians]");
  fSumESD->SetZTitle("#sum_{strip} #Delta/#Delta_{mip}(#eta,#varphi) ");

  fSum = static_cast<TH2D*>(fSumESD->Clone("summed"));
  fSum->SetTitle("Summed cluster signal");
  fSum->SetZTitle("#sum_{cluster} #Delta/#Delta_{mip}(#eta,#varphi) ");
  fSum->SetDirectory(0);
 
  // Perhaps we need to ensure that this histogram has enough range to
  // accommondate all possible ranges - that is, from -1 to the number
  // of strips in this ring(-type) - i.e., NStrips().  Perhaps the
  // axis should be defined with increasin bin size - e.g.,
  //
  //   -1.5,-.5,.5,1.5,...,100.5,128.5,192.5,...,NStrips()
  // 
  fNConsecutive = new TH1D("nConsecutive","# consecutive strips above low cut",
			   201,-1.5,199.5);
  fNConsecutive->SetXTitle("N_{strips}");
  fNConsecutive->SetYTitle("N_{entries}");
  fNConsecutive->SetFillColor(kYellow+2);
  fNConsecutive->SetFillStyle(3001); 
  fNConsecutive->SetDirectory(0);
  
 
#if 0  
  fHits = new TH1D("hits", "Number of hits", 200, 0, 200000);
  fHits->SetDirectory(0);
  fHits->SetXTitle("# of hits");
  fHits->SetFillColor(kGreen+1);
#endif
}
//____________________________________________________________________
AliFMDSharingFilter::RingHistos::RingHistos(const RingHistos& o)
  : AliForwardUtil::RingHistos(o), 
    fBefore(o.fBefore), 
    fAfter(o.fAfter),
    fSingle(o.fSingle),
    fDouble(o.fDouble),
    fTriple(o.fTriple),
    fSinglePerStrip(o.fSinglePerStrip),
    // fDistanceBefore(o.fDistanceBefore),
    // fDistanceAfter(o.fDistanceAfter),    
    fBeforeAfter(o.fBeforeAfter),
    fNeighborsBefore(o.fNeighborsBefore),
    fNeighborsAfter(o.fNeighborsAfter),
    fSumESD(o.fSumESD), //,
    fSum(o.fSum),
    fNConsecutive(o.fNConsecutive)
     //,
    // fHits(o.fHits),
    // fNHits(o.fNHits)
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
  if (&o == this) return *this;
  AliForwardUtil::RingHistos::operator=(o);
  fDet = o.fDet;
  fRing = o.fRing;
  
  if (fBefore) 	       delete  fBefore;
  if (fAfter)  	       delete  fAfter;
  if (fSingle) 	       delete  fSingle;
  if (fDouble) 	       delete  fDouble;
  if (fTriple)         delete  fTriple;
  if (fSinglePerStrip) delete  fSinglePerStrip;
  if (fNConsecutive)   delete  fNConsecutive;
  // if (fDistanceBefore) delete fDistanceBefore;
  // if (fDistanceAfter)  delete fDistanceAfter;
  // if (fHits)   	       delete fHits;
  
  
  fBefore          = static_cast<TH1D*>(o.fBefore->Clone());
  fAfter           = static_cast<TH1D*>(o.fAfter->Clone());
  fSingle          = static_cast<TH1D*>(o.fSingle->Clone());
  fDouble          = static_cast<TH1D*>(o.fDouble->Clone());
  fTriple          = static_cast<TH1D*>(o.fTriple->Clone());
  fSinglePerStrip  = static_cast<TH2D*>(o.fSinglePerStrip->Clone());
  // fDistanceBefore  = static_cast<TH1D*>(o.fDistanceBefore->Clone());
  // fDistanceAfter   = static_cast<TH1D*>(o.fDistanceAfter->Clone());
  fBeforeAfter     = static_cast<TH2D*>(o.fBeforeAfter->Clone());
  fNeighborsBefore = static_cast<TH2D*>(o.fNeighborsBefore->Clone());
  fNeighborsAfter  = static_cast<TH2D*>(o.fNeighborsAfter->Clone());
  // fHits            = static_cast<TH1D*>(o.fHits->Clone());
  fSumESD          = static_cast<TH2D*>(o.fSumESD->Clone());
  fSum             = static_cast<TH2D*>(o.fSum->Clone());
  fNConsecutive    = static_cast<TH1D*>(o.fNConsecutive->Clone());

  return *this;
}
//____________________________________________________________________
AliFMDSharingFilter::RingHistos::~RingHistos()
{
  // 
  // Destructor 
  //
}
#if 0
//____________________________________________________________________
void
AliFMDSharingFilter::RingHistos::Finish()
{
  // 
  // Finish off 
  // 
  //
  // fHits->Fill(fNHits);
}
#endif
//____________________________________________________________________
void
AliFMDSharingFilter::RingHistos::Terminate(const TList* dir, Int_t nEvents)
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

#if 0
  TH1D* before = static_cast<TH1D*>(l->FindObject("esdEloss"));
  TH1D* after  = static_cast<TH1D*>(l->FindObject("anaEloss"));
  if (before) before->Scale(1./nEvents);
  if (after)  after->Scale(1./nEvents);
#endif

  TH2D* summed = static_cast<TH2D*>(l->FindObject("summed"));
  if (summed) summed->Scale(1./nEvents);
  fSum = summed;

  TH2D* summedESD = static_cast<TH2D*>(l->FindObject("summedESD"));
  if (summedESD) summedESD->Scale(1./nEvents);
  fSumESD = summedESD;

  TH1D* consecutive = static_cast<TH1D*>(l->FindObject("nConsecutive"));
  if (consecutive) consecutive->Scale(1./nEvents);
  fNConsecutive= consecutive;
}

//____________________________________________________________________
void
AliFMDSharingFilter::RingHistos::CreateOutputObjects(TList* dir)
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
  d->Add(fSingle);
  d->Add(fDouble);
  d->Add(fTriple);
  d->Add(fSinglePerStrip);
  // d->Add(fDistanceBefore);
  // d->Add(fDistanceAfter);
  d->Add(fBeforeAfter);
  d->Add(fNeighborsBefore);
  d->Add(fNeighborsAfter);
  // d->Add(fHits);
  d->Add(fSumESD);
  d->Add(fSum);
  d->Add(fNConsecutive);

  // Removed to avoid doubly adding the list which destroys 
  // the merging
  //dir->Add(d);
}

//____________________________________________________________________
//
// EOF
//
