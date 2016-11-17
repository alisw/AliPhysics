// This class calculates the inclusive charged particle density
// in each for the 5 FMD rings. 
//
#include "AliFMDDensityCalculator.h"
#include <AliESDFMD.h>
#include <TAxis.h>
#include <TList.h>
#include <TMath.h>
#include "AliForwardCorrectionManager.h"
#include "AliFMDCorrDoubleHit.h"
#include "AliFMDCorrELossFit.h"
#include "AliLog.h"
#include "AliForwardUtil.h"
#include <TH2D.h>
#include <TProfile.h>
#include <THStack.h>
#include <TROOT.h>
#include <TVector3.h>
#include <TStopwatch.h>
#include <TParameter.h>
#include <iostream>
#include <iomanip>
#include <cstring>

ClassImp(AliFMDDensityCalculator)
#if 0
; // For Emacs
#endif 

//____________________________________________________________________
const char* AliFMDDensityCalculator::fgkFolderName = "fmdDensityCalculator";

//____________________________________________________________________
AliFMDDensityCalculator::AliFMDDensityCalculator()
  : TNamed(), 
    fRingHistos(),
    fSumOfWeights(0),
    fWeightedSum(0),
    fCorrections(0),
    fMaxParticles(5),
    fUsePoisson(false),
    fUsePhiAcceptance(kPhiCorrectNch),
    fAccI(0),
    fAccO(0),
    fFMD1iMax(0),
    fFMD2iMax(0),
    fFMD2oMax(0),
    fFMD3iMax(0),
    fFMD3oMax(0),
    fMaxWeights(0),
    fLowCuts(0),
    fEtaLumping(32), 
    fPhiLumping(4),    
    fDebug(0),
    fCuts(),
    fRecalculatePhi(false),
    fMinQuality(AliFMDCorrELossFit::kDefaultQuality),
    fHitThreshold(0.9),
    fCache(),
    fDoTiming(false),
    fHTiming(0), 
    fMaxOutliers(0.05),
    fOutlierCut(0.50)
{
  // 
  // Constructor 
  //
  DGUARD(fDebug, 3, "Default CTOR of FMD density calculator");
}

//____________________________________________________________________
AliFMDDensityCalculator::AliFMDDensityCalculator(const char* title)
  : TNamed(fgkFolderName, title), 
    fRingHistos(), 
    fSumOfWeights(0),
    fWeightedSum(0),
    fCorrections(0),
    fMaxParticles(5),
    fUsePoisson(false),
    fUsePhiAcceptance(kPhiCorrectNch),
    fAccI(0),
    fAccO(0),
    fFMD1iMax(0),
    fFMD2iMax(0),
    fFMD2oMax(0),
    fFMD3iMax(0),
    fFMD3oMax(0),
    fMaxWeights(0),
    fLowCuts(0),
    fEtaLumping(32), 
    fPhiLumping(4),
    fDebug(0),
    fCuts(),
    fRecalculatePhi(false),
    fMinQuality(AliFMDCorrELossFit::kDefaultQuality),
    fHitThreshold(0.9),
    fCache(),
    fDoTiming(false),
    fHTiming(0), 
    fMaxOutliers(0.05),
    fOutlierCut(0.50)
{
  // 
  // Constructor 
  // 
  // Parameters:
  //    name Name of object
  //
  DGUARD(fDebug, 3, "Named CTOR of FMD density calculator: %s", title);
  fRingHistos.SetName(GetName());
  fRingHistos.SetOwner();
  fSumOfWeights = new TH1D("sumOfWeights", "Sum of Landau weights",
			   200, 0, 20);
  fSumOfWeights->SetFillColor(kRed+1);
  fSumOfWeights->SetXTitle("#sum_{i} a_{i} f_{i}(#Delta)");
  fWeightedSum  = new TH1D("weightedSum", "Weighted sum of Landau propability",
			   200, 0, 20);
  fWeightedSum->SetFillColor(kBlue+1);
  fWeightedSum->SetXTitle("#sum_{i} i a_{i} f_{i}(#Delta)");
  fCorrections  = new TH1D("corrections", "Distribution of corrections", 
			   60, 0, 1.2);
  fCorrections->SetFillColor(kBlue+1);
  fCorrections->SetXTitle("correction");

  fAccI = GenerateAcceptanceCorrection('I');
  fAccO = GenerateAcceptanceCorrection('O');

  fMaxWeights = new TH2D("maxWeights", "Maximum i of a_{i}'s to use", 
			 1, 0, 1, 1, 0, 1);
  fMaxWeights->SetXTitle("#eta");
  fMaxWeights->SetDirectory(0);

  fLowCuts = new TH2D("lowCuts", "Low cuts used", 1, 0, 1, 1, 0, 1);
  fLowCuts->SetXTitle("#eta");
  fLowCuts->SetDirectory(0);

}

//____________________________________________________________________
AliFMDDensityCalculator::AliFMDDensityCalculator(const 
						 AliFMDDensityCalculator& o)
  : TNamed(o), 
    fRingHistos(), 
    fSumOfWeights(o.fSumOfWeights),
    fWeightedSum(o.fWeightedSum),
    fCorrections(o.fCorrections),
    fMaxParticles(o.fMaxParticles),
    fUsePoisson(o.fUsePoisson),
    fUsePhiAcceptance(o.fUsePhiAcceptance),
    fAccI(o.fAccI),
    fAccO(o.fAccO),
    fFMD1iMax(o.fFMD1iMax),
    fFMD2iMax(o.fFMD2iMax),
    fFMD2oMax(o.fFMD2oMax),
    fFMD3iMax(o.fFMD3iMax),
    fFMD3oMax(o.fFMD3oMax),
    fMaxWeights(o.fMaxWeights),
    fLowCuts(o.fLowCuts),
    fEtaLumping(o.fEtaLumping), 
    fPhiLumping(o.fPhiLumping),
    fDebug(o.fDebug),
    fCuts(o.fCuts),
    fRecalculatePhi(o.fRecalculatePhi),
    fMinQuality(o.fMinQuality),
    fHitThreshold(o.fHitThreshold),
    fCache(o.fCache),
    fDoTiming(o.fDoTiming),
    fHTiming(o.fHTiming), 
  fMaxOutliers(o.fMaxOutliers),
  fOutlierCut(o.fOutlierCut)
{
  // 
  // Copy constructor 
  // 
  // Parameters:
  //    o Object to copy from 
  //
  DGUARD(fDebug, 3, "Copy CTOR of FMD density calculator");
  TIter    next(&o.fRingHistos);
  TObject* obj = 0;
  while ((obj = next())) fRingHistos.Add(obj);
}

//____________________________________________________________________
AliFMDDensityCalculator::~AliFMDDensityCalculator()
{
  // 
  // Destructor 
  //
  DGUARD(fDebug, 3, "DTOR of FMD density calculator");
  // fRingHistos.Delete();
}

//____________________________________________________________________
AliFMDDensityCalculator&
AliFMDDensityCalculator::operator=(const AliFMDDensityCalculator& o)
{
  // 
  // Assignement operator
  // 
  // Parameters:
  //    o Object to assign from 
  // 
  // Return:
  //    Reference to this object
  //
  DGUARD(fDebug, 3, "Assignment of FMD density calculator");
  if (&o == this) return *this; 
  TNamed::operator=(o);

  fDebug              = o.fDebug;
  fMaxParticles       = o.fMaxParticles;
  fUsePoisson         = o.fUsePoisson;
  fUsePhiAcceptance   = o.fUsePhiAcceptance;
  fAccI               = o.fAccI;
  fAccO               = o.fAccO;
  fFMD1iMax           = o.fFMD1iMax;
  fFMD2iMax           = o.fFMD2iMax;
  fFMD2oMax           = o.fFMD2oMax;
  fFMD3iMax           = o.fFMD3iMax;
  fFMD3oMax           = o.fFMD3oMax;
  fMaxWeights         = o.fMaxWeights;
  fLowCuts            = o.fLowCuts;
  fEtaLumping         = o.fEtaLumping;
  fPhiLumping         = o.fPhiLumping;
  fCuts               = o.fCuts;
  fRecalculatePhi     = o.fRecalculatePhi;
  fMinQuality         = o.fMinQuality;
  fHitThreshold       = o.fHitThreshold;
  fCache              = o.fCache;
  fDoTiming           = o.fDoTiming;
  fHTiming            = o.fHTiming;
  fMaxOutliers        = o.fMaxOutliers;
  fOutlierCut         = o.fOutlierCut;

  fRingHistos.Delete();
  TIter    next(&o.fRingHistos);
  TObject* obj = 0;
  while ((obj = next())) fRingHistos.Add(obj);
  
  return *this;
}

//____________________________________________________________________
void
AliFMDDensityCalculator::SetupForData(const TAxis& axis)
{
  // Intialize this sub-algorithm 
  //
  // Parameters:
  //   etaAxis   Eta axis
  DGUARD(fDebug, 1, "Initialize FMD density calculator");
  CacheMaxWeights(axis);
 
  fCache.Init(axis);

  TIter    next(&fRingHistos);
  RingHistos* o = 0;
  while ((o = static_cast<RingHistos*>(next()))) {
    o->SetupForData(axis);
    // o->fMultCut = fCuts.GetFixedCut(o->fDet, o->fRing);
    // o->fPoisson.Init(o->fDet,o->fRing,fEtaLumping, fPhiLumping);
  }
}

//____________________________________________________________________
AliFMDDensityCalculator::RingHistos*
AliFMDDensityCalculator::GetRingHistos(UShort_t d, Char_t r) const
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
  if (idx < 0 || idx >= fRingHistos.GetEntries()) {
    AliWarning(Form("Index %d of FMD%d%c out of range", idx, d, r));
    return 0;
  }
  
  return static_cast<RingHistos*>(fRingHistos.At(idx));
}

namespace {
  Double_t Rng2Cut(UShort_t d, Char_t r, Int_t xbin, TH2* h) 
  {
    Double_t ret = 1024;
    if (xbin < 1 && xbin > h->GetXaxis()->GetNbins()) return ret;
    Int_t ybin = 0;							
    switch(d) {								
    case 1: ybin = 1; break;						
    case 2: ybin = (r=='i' || r=='I') ? 2 : 3; break;			
    case 3: ybin = (r=='i' || r=='I') ? 4 : 5; break;			
    default: return ret;
    }									
    ret = h->GetBinContent(xbin,ybin);					
    return ret;
  }
}

//____________________________________________________________________
Double_t
AliFMDDensityCalculator::GetMultCut(UShort_t d, Char_t r, Int_t ieta,
				    Bool_t /*errors*/) const
{
  // 
  // Get the multiplicity cut.  If the user has set fMultCut (via
  // SetMultCut) then that value is used.  If not, then the lower
  // value of the fit range for the enery loss fits is returned.
  // 
  // Return:
  //    Lower cut on multiplicity
  //
  return Rng2Cut(d, r, ieta, fLowCuts);
  // return fCuts.GetMultCut(d,r,ieta,errors);
}
    
//____________________________________________________________________
Double_t
AliFMDDensityCalculator::GetMultCut(UShort_t d, Char_t r, Double_t eta,
				    Bool_t /*errors*/) const
{
  // 
  // Get the multiplicity cut.  If the user has set fMultCut (via
  // SetMultCut) then that value is used.  If not, then the lower
  // value of the fit range for the enery loss fits is returned.
  // 
  // Return:
  //    Lower cut on multiplicity
  //
  Int_t ieta = fLowCuts->GetXaxis()->FindBin(eta);				
  return Rng2Cut(d, r, ieta, fLowCuts);
  // return fCuts.GetMultCut(d,r,eta,errors);
}

#ifndef NO_TIMING
# define START_TIMER(T) if (fDoTiming) T.Start(true)
# define GET_TIMER(T,V) if (fDoTiming) V = T.CpuTime()
# define ADD_TIMER(T,V) if (fDoTiming) V += T.CpuTime()
#else
# define START_TIMER(T) do {} while (false)
# define GET_TIMER(T,V) do {} while (false)
# define ADD_TIMER(T,V) do {} while (false)
#endif

//____________________________________________________________________
Bool_t
AliFMDDensityCalculator::Calculate(const AliESDFMD&        fmd,
				   AliForwardUtil::Histos& hists,
				   Bool_t                  lowFlux,
				   Double_t                /*cent*/, 
				   const TVector3&         ip)
{
  // 
  // Do the calculations 
  // 
  // Parameters:
  //    fmd      AliESDFMD object (possibly) corrected for sharing
  //    hists    Histogram cache
  //    vtxBin   Vertex bin 
  //    lowFlux  Low flux flag. 
  // 
  // Return:
  //    true on successs 
  DGUARD(fDebug, 1, "Calculate density in FMD density calculator");

  TStopwatch timer;
  TStopwatch totalT;
  
  // First measurements of timing
  //  Re-calculation      : fraction of sum 32.0%   of total 18.1%
  //  N_{particle}        : fraction of sum 15.2%   of total  8.6%
  //  Correction          : fraction of sum 26.4%   of total 14.9%
  //  #phi acceptance     : fraction of sum  0.2%   of total  0.1%
  //  Copy to cache       : fraction of sum  3.9%   of total  2.2%
  //  Poisson calculation : fraction of sum 18.7%   of total 10.6%
  //  Diagnostics         : fraction of sum  3.7%   of total  2.1%
  Double_t nPartTime  = 0;
  Double_t corrTime   = 0;
  Double_t rePhiTime  = 0;
  Double_t copyTime   = 0;
  Double_t poissonTime= 0;
  Double_t diagTime   = 0;
  // Double_t ipPhi      = TMath::ATan2(ip.Y(),ip.X());
  // Double_t ipR        = TMath::Sqrt(TMath::Power(ip.X(),2)+
  //                       TMath::Power(ip.Y(),2));
  START_TIMER(totalT);
  
  Double_t etaCache[20*512]; // Same number of strips per ring 
  Double_t phiCache[20*512]; // whether it is inner our outer. 
  // We do not use TArrayD because we do not wont a bounds check 
  // TArrayD etaCache(20*512); // Same number of strips per ring
  // TArrayD phiCache(20*512); // whether it is inner our outer. 
  
  // --- Loop over detectors -----------------------------------------
  for (UShort_t d=1; d<=3; d++) { 
    UShort_t nr = (d == 1 ? 1 : 2);
    for (UShort_t q=0; q<nr; q++) { 
      Char_t      r = (q == 0 ? 'I' : 'O');
      UShort_t    ns= (q == 0 ?  20 :  40);
      UShort_t    nt= (q == 0 ? 512 : 256);
      TH2D*       h = hists.Get(d,r);
      RingHistos* rh= GetRingHistos(d,r);
      if (!rh) { 
	AliError(Form("No ring histogram found for FMD%d%c", d, r));
	fRingHistos.ls();
	return false;
      }
      // rh->fPoisson.SetObject(d,r,vtxbin,cent);
      rh->fPoisson.Reset(0);
      rh->fTotal->Reset();
      rh->fGood->Reset();
      // rh->ResetPoissonHistos(h, fEtaLumping, fPhiLumping);

      // Reset our eta cache 
      // for (Int_t i = 0; i < 20*512; i++) 
      // etaCache[i] = phiCache[i] = AliESDFMD::kInvalidEta;
      memset(etaCache, 0, sizeof(Double_t)*20*512);
      memset(phiCache, 0, sizeof(Double_t)*20*512);
      // etaCache.Reset(AliESDFMD::kInvalidEta);
      // phiCache.Reset(AliESDFMD::kInvalidEta);

      // --- Loop over sectors and strips ----------------------------
      for (UShort_t s=0; s<ns; s++) { 
	for (UShort_t t=0; t<nt; t++) {
	  
	  Float_t  mult   = fmd.Multiplicity(d,r,s,t);
	  Double_t phi    = fmd.Phi(d,r,s,t) * TMath::DegToRad();
	  Double_t eta    = fmd.Eta(d,r,s,t);
	  Double_t oldPhi = phi;
	  Double_t oldEta = eta;
	  START_TIMER(timer);
	  if (fRecalculatePhi) {
	    // Correct for (x,y) off set of the interaction point 
	    // AliForwardUtil::GetEtaPhiFromStrip(r,t,eta,phi,ip.X(),ip.Y());
	    if (!AliForwardUtil::GetEtaPhi(d,r,s,t,ip,eta,phi) ||
		TMath::Abs(eta) < 1) {
	      AliWarningF("FMD%d%c[%2d,%3d] (%f,%f,%f) eta=%f phi=%f (%f)",
			  d, r, s, t, ip.X(), ip.Y(), ip.Z(), eta,
			  phi, oldEta);
	      eta = oldEta;
	      phi = oldPhi;
	    }
	    DMSG(fDebug, 10, "IP(x,y,z)=%f,%f,%f Eta=%f -> %f Phi=%f -> %f",
		 ip.X(), ip.Y(), ip.Z(), oldEta, eta, oldPhi, phi);
	  }
	  ADD_TIMER(timer,rePhiTime);
	  START_TIMER(timer);
	  etaCache[s*nt+t] = eta;
	  phiCache[s*nt+t] = phi;

	  // --- Check this strip ------------------------------------
	  rh->fTotal->Fill(eta);
	  if (mult == AliESDFMD::kInvalidMult) { //  || mult > 20) {
	    // Do not count invalid stuff 
	    rh->fELoss->Fill(-1);
	    // rh->fEvsN->Fill(mult,-1);
	    // rh->fEvsM->Fill(mult,-1);
	    continue;
	  }
	  if (mult > 20) 
	    AliWarningF("Raw multiplicity of FMD%d%c[%02d,%03d] = %f > 20",
			d, r, s, t, mult);
	  // --- Automatic calculation of acceptance -----------------
	  rh->fGood->Fill(eta);

	  // --- If we asked to re-calculate phi for (x,y) IP --------
	  // START_TIMER(timer);
	  // if (fRecalculatePhi) {
	  // oldPhi = phi;
	  //  phi = AliForwardUtil::GetPhiFromStrip(r, t, phi, ip.X(), ip.Y());
	  // }
	  // phiCache[s*nt+t] = phi;
	  // ADD_TIMER(timer,rePhiTime);

	  // --- Apply phi corner correction to eloss ----------------
	  if (fUsePhiAcceptance == kPhiCorrectELoss) 
	    mult *= AcceptanceCorrection(r,t);

	  // --- Get the low multiplicity cut ------------------------
	  Double_t cut  = 1024;
	  if (eta != AliESDFMD::kInvalidEta) cut = GetMultCut(d, r, eta,false);
	  else AliWarningF("Eta for FMD%d%c[%02d,%03d] is invalid: %f", 
			   d, r, s, t, eta);

	  // --- Now caluculate Nch for this strip using fits --------
	  START_TIMER(timer);
	  Double_t n   = 0;
	  if (cut > 0 && mult > cut) n = NParticles(mult,d,r,eta,lowFlux);
	  rh->fELoss->Fill(mult);
	  // rh->fEvsN->Fill(mult,n);
	  // rh->fEtaVsN->Fill(eta, n);
	  ADD_TIMER(timer,nPartTime);
	  
	  // --- Calculate correction if needed ----------------------
	  START_TIMER(timer);
	  // Temporary stuff - remove Correction call 
	  Double_t c = 1;
	  if (fUsePhiAcceptance == kPhiCorrectNch) 
	    c = AcceptanceCorrection(r,t);
	  // Double_t c = Correction(d,r,t,eta,lowFlux);
	  ADD_TIMER(timer,corrTime);
	  fCorrections->Fill(c);
	  if (c > 0) n /= c;
	  // rh->fEvsM->Fill(mult,n);
	  // rh->fEtaVsM->Fill(eta, n);
	  rh->fCorr  ->Fill(eta, c);
	  
	  // --- Accumulate Poisson statistics -----------------------
	  Bool_t hit = (n > fHitThreshold && c > 0);
	  if (hit) {
	    rh->fELossUsed->Fill(mult);
	    if (fRecalculatePhi) {
	      rh->fPhiBefore->Fill(oldPhi);
	      rh->fPhiAfter->Fill(phi);
	      rh->fEtaBefore->Fill(oldEta);
	      rh->fEtaAfter->Fill(oldEta);	      
	    }
	    rh->fSignal->Fill(eta, mult);
	  }
	  rh->fPoisson.Fill(t,s,hit,1./c);
	  h->Fill(eta,phi,n);

	  // --- If we use ELoss fits, apply now ---------------------
	  if (!fUsePoisson) rh->fDensity->Fill(eta,phi,n);
	} // for t
      } // for s 

      // --- Automatic acceptance - Calculate as an efficiency -------
      // This is very fast, so we do not bother to time it 
      rh->fGood->Divide(rh->fGood, rh->fTotal, 1, 1, "B");

      // --- Make a copy and reset as needed -------------------------
      START_TIMER(timer);
      TH2D* hclone = fCache.Get(d,r);
      // hclone->Reset();
      // TH2D* hclone = static_cast<TH2D*>(h->Clone("hclone"));
      if (!fUsePoisson) hclone->Reset();
      else { 
	for (Int_t i = 0; i <= h->GetNbinsX()+1; i++) { 
	  for (Int_t j = 0; j <= h->GetNbinsY()+1; j++) {
	    hclone->SetBinContent(i,j,h->GetBinContent(i,j));
	    hclone->SetBinError(i,j,h->GetBinError(i,j));
	  }
	}
	// hclone->Add(h); 
	h->Reset(); 
      }
      ADD_TIMER(timer,copyTime);
      
      // --- Store Poisson result ------------------------------------
      START_TIMER(timer);
      TH2D* poisson = rh->fPoisson.Result();
      for (Int_t t=0; t < poisson->GetNbinsX(); t++) { 
	for (Int_t s=0; s < poisson->GetNbinsY(); s++) { 
	  
	  Double_t poissonV = poisson->GetBinContent(t+1,s+1);
	  // Use cached eta - since the calls to GetEtaFromStrip and
	  // GetPhiFromStrip are _very_ expensive
	  Double_t  phi  = phiCache[s*nt+t];
	  Double_t  eta  = etaCache[s*nt+t]; 
	  // Double_t  phi  = fmd.Phi(d,r,s,t) * TMath::DegToRad();
	  // Double_t  eta  = fmd.Eta(d,r,s,t);
	  if (fUsePoisson) {
	    h->Fill(eta,phi,poissonV);
	    rh->fDensity->Fill(eta, phi, poissonV);
	  }
	  else
	    hclone->Fill(eta,phi,poissonV);
	}
      }
      ADD_TIMER(timer,poissonTime);
      
      // --- Make diagnostics - eloss vs poisson ---------------------
      START_TIMER(timer);
      Int_t nY = h->GetNbinsY();
      Int_t nIn  = 0; // Count non-outliers
      Int_t nOut = 0; // Count outliers
      for (Int_t ieta=1; ieta <= h->GetNbinsX(); ieta++) { 
	// Set the overflow bin to contain the phi acceptance 
	Double_t phiAcc  = rh->fGood->GetBinContent(ieta);
	Double_t phiAccE = rh->fGood->GetBinError(ieta);
	h->SetBinContent(ieta, nY+1, phiAcc);
	h->SetBinError(ieta, nY+1, phiAccE);
	Double_t eta     = h->GetXaxis()->GetBinCenter(ieta);
	rh->fPhiAcc->Fill(eta, ip.Z(), phiAcc);
	for (Int_t iphi=1; iphi<= nY; iphi++) { 
	  
	  Double_t poissonV =  0; //h->GetBinContent(,s+1);
	  Double_t eLossV =  0;
	  if(fUsePoisson) { 
	    poissonV = h->GetBinContent(ieta,iphi);
	    eLossV  = hclone->GetBinContent(ieta,iphi);
	  }
	  else { 
	    poissonV = hclone->GetBinContent(ieta,iphi);
	    eLossV  = h->GetBinContent(ieta,iphi);
	  }
	  
	  if (poissonV < 1e-12 && eLossV < 1e-12) 
	    // we do not care about trivially empty bins 
	    continue;
				      
	  Bool_t   outlier = CheckOutlier(eLossV, poissonV, fOutlierCut);
	  Double_t rel     = eLossV < 1e-12 ? 0 : (poissonV - eLossV) / eLossV;
	  if (outlier) {
	    rh->fELossVsPoissonOut->Fill(eLossV, poissonV);
	    rh->fDiffELossPoissonOut->Fill(rel);
	    nOut++;
	}
	  else {
	    rh->fELossVsPoisson->Fill(eLossV, poissonV);
	    rh->fDiffELossPoisson->Fill(rel);
	    nIn++;
	  } // if (outlier)
	} // for (iphi)
      } // for (ieta)
      Int_t    nTotal   = (nIn+nOut);
      Double_t outRatio = (nTotal > 0 ? Double_t(nOut) / nTotal : 0);
      rh->fOutliers->Fill(outRatio);
      if (outRatio < fMaxOutliers) rh->fPoisson.FillDiagnostics();
      else                         h->SetBit(AliForwardUtil::kSkipRing);
      ADD_TIMER(timer,diagTime);
      // delete hclone;
      
    } // for q
  } // for d

  if (fDoTiming) {
    // fHTiming->Fill(1,reEtaTime);
    fHTiming->Fill(2,nPartTime);
    fHTiming->Fill(3,corrTime);
    fHTiming->Fill(4,rePhiTime);
    fHTiming->Fill(5,copyTime);
    fHTiming->Fill(6,poissonTime);
    fHTiming->Fill(7,diagTime);
    fHTiming->Fill(8,totalT.CpuTime());
  }

  return kTRUE;
}

//_____________________________________________________________________
Bool_t 
AliFMDDensityCalculator::CheckOutlier(Double_t eloss, 
				      Double_t poisson,
				      Double_t cut) const
{
  if (eloss < 1e-6) return true;
  Double_t diff = TMath::Abs(poisson - eloss) / eloss;
  return diff > cut;
}
//_____________________________________________________________________
Int_t
AliFMDDensityCalculator::FindMaxWeight(const AliFMDCorrELossFit* cor,
				       UShort_t d, Char_t r, Int_t iEta) const
{
  // 
  // Find the max weight to use for FMD<i>dr</i> in eta bin @a iEta
  // 
  // Parameters:
  //    cor   Correction
  //    d     Detector 
  //    r     Ring 
  //    iEta  Eta bin 
  //
  DGUARD(fDebug, 10, "Find maximum weight in FMD density calculator");
  if(!cor) return -1; 

  AliFMDCorrELossFit::ELossFit* fit = cor->FindFit(d,r,iEta, -1);
  if (!fit) { 
    // AliWarning(Form("No energy loss fit for FMD%d%c at eta=%f", d, r, eta));
    return -1;
  }
  return fit->FindMaxWeight(2*AliFMDCorrELossFit::ELossFit::fgMaxRelError, 
			    AliFMDCorrELossFit::ELossFit::fgLeastWeight, 
			    fMaxParticles);
}
//_____________________________________________________________________
Int_t
AliFMDDensityCalculator::FindMaxWeight(const AliFMDCorrELossFit* cor,
				       UShort_t d, Char_t r, Double_t eta) const
{
  // 
  // Find the max weight to use for FMD<i>dr</i> in eta bin @a iEta
  // 
  // Parameters:
  //    cor   Correction
  //    d     Detector 
  //    r     Ring 
  //    eta   Eta
  //
  DGUARD(fDebug, 10, "Find maximum weight in FMD density calculator");
  if(!cor) return -1; 

  AliFMDCorrELossFit::ELossFit* fit = cor->FindFit(d,r,eta, -1);
  if (!fit) { 
    // AliWarning(Form("No energy loss fit for FMD%d%c at eta=%f", d, r, eta));
    return -1;
  }
  return fit->FindMaxWeight(2*AliFMDCorrELossFit::ELossFit::fgMaxRelError, 
			    AliFMDCorrELossFit::ELossFit::fgLeastWeight, 
			    fMaxParticles);
}
  
//_____________________________________________________________________
void
AliFMDDensityCalculator::CacheMaxWeights(const TAxis& axis)
{
  // 
  // Find the max weights and cache them 
  // 
  DGUARD(fDebug, 2, "Cache maximum weights in FMD density calculator");
  AliForwardCorrectionManager&  fcm = AliForwardCorrectionManager::Instance();
  const AliFMDCorrELossFit*     cor = fcm.GetELossFit();
  cor->CacheBins(fMinQuality);
  cor->Print(fDebug > 5 ? "RCS" : "C");

  TAxis eta(axis.GetNbins(),
	    axis.GetXmin(),
	    axis.GetXmax());
  
  eta.Set(cor->GetEtaAxis().GetNbins(), 
	  cor->GetEtaAxis().GetXmin(), 
	  cor->GetEtaAxis().GetXmax());

  Int_t nEta = eta.GetNbins();
  fFMD1iMax.Set(nEta);
  fFMD2iMax.Set(nEta);
  fFMD2oMax.Set(nEta);
  fFMD3iMax.Set(nEta);
  fFMD3oMax.Set(nEta);
  
  fMaxWeights->SetBins(nEta, eta.GetXmin(), eta.GetXmax(), 5, .5, 5.5);
  fMaxWeights->GetYaxis()->SetBinLabel(1, "FMD1i");
  fMaxWeights->GetYaxis()->SetBinLabel(2, "FMD2i");
  fMaxWeights->GetYaxis()->SetBinLabel(3, "FMD2o");
  fMaxWeights->GetYaxis()->SetBinLabel(4, "FMD3i");
  fMaxWeights->GetYaxis()->SetBinLabel(5, "FMD3o");

  AliInfo(Form("Get eta axis with %d bins from %f to %f",
	       nEta, eta.GetXmin(), eta.GetXmax()));
  fLowCuts->SetBins(nEta, eta.GetXmin(), eta.GetXmax(), 5, .5, 5.5);
  fLowCuts->GetYaxis()->SetBinLabel(1, "FMD1i");
  fLowCuts->GetYaxis()->SetBinLabel(2, "FMD2i");
  fLowCuts->GetYaxis()->SetBinLabel(3, "FMD2o");
  fLowCuts->GetYaxis()->SetBinLabel(4, "FMD3i");
  fLowCuts->GetYaxis()->SetBinLabel(5, "FMD3o");

  for (Int_t i = 0; i < nEta; i++) {
    Double_t leta = fMaxWeights->GetXaxis()->GetBinCenter(i+1);
    Double_t w[5];
    w[0] = fFMD1iMax[i] = FindMaxWeight(cor, 1, 'I', leta);
    w[1] = fFMD2iMax[i] = FindMaxWeight(cor, 2, 'I', leta);
    w[2] = fFMD2oMax[i] = FindMaxWeight(cor, 2, 'O', leta);
    w[3] = fFMD3iMax[i] = FindMaxWeight(cor, 3, 'I', leta);
    w[4] = fFMD3oMax[i] = FindMaxWeight(cor, 3, 'O', leta);
    for (Int_t j = 0; j < 5; j++) 
      if (w[j] > 0) fMaxWeights->SetBinContent(i+1, j+1, w[j]);
  }

  // Cache cuts in histogram
  fCuts.FillHistogram(fLowCuts);
}

//_____________________________________________________________________
Int_t
AliFMDDensityCalculator::GetMaxWeight(UShort_t d, Char_t r, Int_t iEta) const
{
  // 
  // Find the (cached) maximum weight for FMD<i>dr</i> in 
  // @f$\eta@f$ bin @a iEta
  // 
  // Parameters:
  //    d     Detector
  //    r     Ring
  //    iEta  Eta bin
  // 
  // Return:
  //    max weight or <= 0 in case of problems 
  //
  if (iEta < 0) return -1;

  const TArrayI* max  = 0;
  switch (d) { 
  case 1:  max = &fFMD1iMax;                                       break;
  case 2:  max = (r == 'I' || r == 'i' ? &fFMD2iMax : &fFMD2oMax); break;
  case 3:  max = (r == 'I' || r == 'i' ? &fFMD3iMax : &fFMD3oMax); break;
  }
  if (!max) { 
    AliWarning(Form("No array for FMD%d%c", d, r));
    return -1;
  }
  
  if (iEta >= max->fN) { 
    AliWarning(Form("Eta bin %3d out of bounds [0,%d]", 
		    iEta, max->fN-1));
    return -1;
  }

  AliDebug(30,Form("Max weight for FMD%d%c eta bin %3d: %d", d, r, iEta, 
		   max->At(iEta)));
  return max->At(iEta);
}

//_____________________________________________________________________
Int_t
AliFMDDensityCalculator::GetMaxWeight(UShort_t d, Char_t r, Float_t eta) const
{
  // 
  // Find the (cached) maximum weight for FMD<i>dr</i> iat
  // @f$\eta@f$ 
  // 
  // Parameters:
  //    d     Detector
  //    r     Ring
  //    eta   Eta bin
  // 
  // Return:
  //    max weight or <= 0 in case of problems 
  //
  AliForwardCorrectionManager&  fcm  = AliForwardCorrectionManager::Instance();
  Int_t                         iEta = fcm.GetELossFit()->FindEtaBin(eta) -1;

  return GetMaxWeight(d, r, iEta);
}

//_____________________________________________________________________
Float_t 
AliFMDDensityCalculator::NParticles(Float_t  mult, 
				    UShort_t d, 
				    Char_t   r, 
				    Float_t  eta,
				    Bool_t   lowFlux) const
{
  // 
  // Get the number of particles corresponding to the signal mult
  // 
  // Parameters:
  //    mult     Signal
  //    d        Detector
  //    r        Ring 
  //    s        Sector 
  //    t        Strip (not used)
  //    v        Vertex bin 
  //    eta      Pseudo-rapidity 
  //    lowFlux  Low-flux flag 
  // 
  // Return:
  //    The number of particles 
  //
  // if (mult <= GetMultCut()) return 0;
  DGUARD(fDebug, 3, "Calculate Nch in FMD density calculator");
  if (lowFlux) return 1;
  
  AliForwardCorrectionManager&  fcm = AliForwardCorrectionManager::Instance();
  AliFMDCorrELossFit::ELossFit* fit = fcm.GetELossFit()->FindFit(d,r,eta, -1);
  if (!fit) { 
    AliWarning(Form("No energy loss fit for FMD%d%c at eta=%f qual=%d", 
		    d, r, eta, fMinQuality));
    return 0;
  }
  
  Int_t    m   = GetMaxWeight(d,r,eta); // fit->FindMaxWeight();
  if (m < 1) { 
    AliWarning(Form("No good fits for FMD%d%c at eta=%f", d, r, eta));
    return 0;
  }
  
  UShort_t n   = TMath::Min(fMaxParticles, UShort_t(m));
  Double_t ret = fit->EvaluateWeighted(mult, n);
  
  if (fDebug > 10) {
    AliInfo(Form("FMD%d%c, eta=%7.4f, %8.5f -> %8.5f", d, r, eta, mult, ret));
  }
    
  fWeightedSum->Fill(ret);
  fSumOfWeights->Fill(ret);
  
  return ret;
}

//_____________________________________________________________________
Float_t 
AliFMDDensityCalculator::Correction(UShort_t d, 
				    Char_t   r, 
				    UShort_t t, 
				    Float_t  eta,
				    Bool_t   lowFlux) const
{
  // 
  // Get the inverse correction factor.  This consist of
  // 
  // - acceptance correction (corners of sensors) 
  // - double hit correction (for low-flux events) 
  // - dead strip correction 
  // 
  // Parameters:
  //    d        Detector
  //    r        Ring 
  //    s        Sector 
  //    t        Strip (not used)
  //    v        Vertex bin 
  //    eta      Pseudo-rapidity 
  //    lowFlux  Low-flux flag 
  // 
  // Return:
  //    
  //
  DGUARD(fDebug, 10, "Apply correction in FMD density calculator");
  Float_t correction = 1; 
  if (fUsePhiAcceptance == kPhiCorrectNch) 
    correction = AcceptanceCorrection(r,t);
  if (lowFlux) { 
    AliForwardCorrectionManager&  fcm = AliForwardCorrectionManager::Instance();

    TH1D* dblHitCor = 0;
    if (fcm.GetDoubleHit()) 
      dblHitCor = fcm.GetDoubleHit()->GetCorrection(d,r);

    if (dblHitCor) {
      Double_t dblC = dblHitCor->GetBinContent(dblHitCor->FindBin(eta));
      if (dblC > 0) correction *= dblC;
    }
    // else {
    //   AliWarning(Form("Missing double hit correction for FMD%d%c",d,r));
    // }
  }
  return correction;
}

//_____________________________________________________________________
TH1D*
AliFMDDensityCalculator::GenerateAcceptanceCorrection(Char_t r) const
{
  // 
  // Generate the acceptance corrections 
  // 
  // Parameters:
  //    r Ring to generate for 
  // 
  // Return:
  //    Newly allocated histogram of acceptance corrections
  //
  DGUARD(fDebug, 3, "Make acceptance correction in FMD density calculator");
  const Double_t ic1[] = { 4.9895, 15.3560 };
  const Double_t ic2[] = { 1.8007, 17.2000 };
  const Double_t oc1[] = { 4.2231, 26.6638 };
  const Double_t oc2[] = { 1.8357, 27.9500 };
  const Double_t* c1   = (r == 'I' || r == 'i' ? ic1      : oc1);
  const Double_t* c2   = (r == 'I' || r == 'i' ? ic2      : oc2);
  Double_t  minR       = (r == 'I' || r == 'i' ?   4.5213 :  15.4);
  Double_t  maxR       = (r == 'I' || r == 'i' ?  17.2    :  28.0);
  Int_t     nStrips    = (r == 'I' || r == 'i' ? 512      : 256);
  Int_t     nSec       = (r == 'I' || r == 'i' ?  20      :  40);
  Float_t   basearc    = 2 * TMath::Pi() / nSec;
  Double_t  rad        = maxR - minR;
  Float_t   segment    = rad / nStrips;
  Float_t   cr         = TMath::Sqrt(c1[0]*c1[0]+c1[1]*c1[1]);

  // Numbers used to find end-point of strip.
  // (See http://mathworld.wolfram.com/Circle-LineIntersection.html)
  Float_t D            = c1[0] * c2[1] - c1[1] * c2[0];
  Float_t dx           = c2[0] - c1[0];
  Float_t dy           = c2[1] - c1[1];
  Float_t dr           = TMath::Sqrt(dx*dx+dy*dy);

  TH1D* ret = new TH1D(Form("acc%c", r), 
		       Form("Acceptance correction for FMDx%c", r), 
		       nStrips, -.5, nStrips-.5);
  ret->SetXTitle("Strip");
  ret->SetYTitle("#varphi acceptance");
  ret->SetDirectory(0);
  ret->SetFillColor(r == 'I' || r == 'i' ? kRed+1 : kBlue+1);
  ret->SetFillStyle(3001);

  for (Int_t t = 0; t < nStrips; t++) { 
    Float_t   radius     = minR + t * segment;
    
    // If the radius of the strip is smaller than the radius corresponding 
    // to the first corner we have a full strip length 
    if (radius <= cr) {
      ret->SetBinContent(t+1, 1);
      continue;
    }

    // Next, we should find the end-point of the strip - that is, 
    // the coordinates where the line from c1 to c2 intersects a circle 
    // with radius given by the strip. 
    // (See http://mathworld.wolfram.com/Circle-LineIntersection.html)
    // Calculate the determinant 
    Float_t det = radius * radius * dr * dr - D*D;

    if (det <=  0) { 
      // <0 means No intersection
      // =0 means Exactly tangent
      ret->SetBinContent(t+1, 1);
      continue;
    }

    // Calculate end-point and the corresponding opening angle 
    Float_t x   = (+D * dy + dx * TMath::Sqrt(det)) / dr / dr;
    Float_t y   = (-D * dx + dy * TMath::Sqrt(det)) / dr / dr;
    Float_t th  = TMath::ATan2(x, y);

    ret->SetBinContent(t+1, th / basearc);
  }
  return ret;
}

//_____________________________________________________________________
Float_t 
AliFMDDensityCalculator::AcceptanceCorrection(Char_t r, UShort_t t) const
{
  // 
  // Get the acceptance correction for strip @a t in an ring of type @a r
  // 
  // Parameters:
  //    r  Ring type ('I' or 'O')
  //    t  Strip number 
  // 
  // Return:
  //    Inverse acceptance correction 
  //
  TH1D* acc = (r == 'I' || r == 'i' ? fAccI : fAccO);
  return acc->GetBinContent(t+1);
}

//____________________________________________________________________
void
AliFMDDensityCalculator::Terminate(const TList* dir, TList* output, 
				   Int_t nEvents)
{
  // 
  // Scale the histograms to the total number of events 
  // 
  // Parameters:
  //    dir     where to put the output
  //    nEvents Number of events 
  //
  DGUARD(fDebug, 1, "Scale histograms in FMD density calculator");
  if (nEvents <= 0) return;
  TList* d = static_cast<TList*>(dir->FindObject(GetName()));
  if (!d) return;

  TList* out = new TList;
  out->SetName(d->GetName());
  out->SetOwner();
  
  fRingHistos.Add(new RingHistos(1, 'I'));
  fRingHistos.Add(new RingHistos(2, 'I'));
  fRingHistos.Add(new RingHistos(2, 'O'));
  fRingHistos.Add(new RingHistos(3, 'I'));
  fRingHistos.Add(new RingHistos(3, 'O'));
  TIter    next(&fRingHistos);
  RingHistos* o = 0;
  THStack* sums = new THStack("sums", "sums of ring signals");
  while ((o = static_cast<RingHistos*>(next()))) {
    o->Terminate(d, nEvents);
    if (!o->fDensity) { 
      Warning("Terminate", "No density in %s", o->GetName());
      continue;
    }
    TH1D* sum = o->fDensity->ProjectionX(o->GetName(), 1, 
					 o->fDensity->GetNbinsY(),"e");
    sum->Scale(1., "width");
    sum->SetTitle(o->GetName());
    sum->SetDirectory(0);
    sum->SetYTitle("#sum N_{ch,incl}");
    sums->Add(sum);
  }
  out->Add(sums);
  output->Add(out);
}

//____________________________________________________________________
void
AliFMDDensityCalculator::CreateOutputObjects(TList* dir)
{
  // 
  // Output diagnostic histograms to directory 
  // 
  // Parameters:
  //    dir List to write in
  //  
  DGUARD(fDebug, 1, "Define output FMD density calculator");
  TList* d = new TList;
  d->SetOwner();
  d->SetName(GetName());
  dir->Add(d);
  d->Add(fWeightedSum);
  d->Add(fSumOfWeights);
  d->Add(fCorrections);
  d->Add(fAccI);
  d->Add(fAccO);
  d->Add(fMaxWeights);
  d->Add(fLowCuts);

  TParameter<int>* nFiles = new TParameter<int>("nFiles", 1);
  nFiles->SetMergeMode('+');
  
  // d->Add(sigma);
  d->Add(AliForwardUtil::MakeParameter("maxParticle",  fMaxParticles));
  d->Add(AliForwardUtil::MakeParameter("minQuality",   fMinQuality));
  d->Add(AliForwardUtil::MakeParameter("method",       fUsePoisson));
  d->Add(AliForwardUtil::MakeParameter("phiAcceptance",fUsePhiAcceptance));
  d->Add(AliForwardUtil::MakeParameter("etaLumping",   fEtaLumping));
  d->Add(AliForwardUtil::MakeParameter("phiLumping",   fPhiLumping));
  d->Add(AliForwardUtil::MakeParameter("recalcPhi",    fRecalculatePhi));
  d->Add(AliForwardUtil::MakeParameter("maxOutliers",  fMaxOutliers));
  d->Add(AliForwardUtil::MakeParameter("outlierCut",   fOutlierCut));
  d->Add(AliForwardUtil::MakeParameter("hitThreshold", fHitThreshold));
  d->Add(nFiles);
  // d->Add(nxi);
  fCuts.Output(d,"lCuts");

  fRingHistos.Add(new RingHistos(1, 'I'));
  fRingHistos.Add(new RingHistos(2, 'I'));
  fRingHistos.Add(new RingHistos(2, 'O'));
  fRingHistos.Add(new RingHistos(3, 'I'));
  fRingHistos.Add(new RingHistos(3, 'O'));
  TIter    next(&fRingHistos);
  RingHistos* o = 0;
  while ((o = static_cast<RingHistos*>(next()))) {
    o->fPoisson.SetLumping(fEtaLumping, fPhiLumping);
    o->CreateOutputObjects(d);
  }

  if (!fDoTiming) return;

  fHTiming = new TProfile("timing", "#LTt_{part}#GT", 8, .5, 8.5);
  fHTiming->SetDirectory(0);
  fHTiming->SetYTitle("#LTt_{part}#GT");
  fHTiming->SetXTitle("Part");
  fHTiming->SetFillColor(kRed+1);
  fHTiming->SetFillStyle(3001);
  fHTiming->SetMarkerStyle(20);
  fHTiming->SetMarkerColor(kBlack);
  fHTiming->SetLineColor(kBlack);
  fHTiming->SetStats(0);
  TAxis* xaxis = fHTiming->GetXaxis();
  xaxis->SetBinLabel(1, "Re-calculation of #eta");
  xaxis->SetBinLabel(2, "N_{particle}");
  xaxis->SetBinLabel(3, "Correction");
  xaxis->SetBinLabel(4, "Re-calculation of #phi");
  xaxis->SetBinLabel(5, "Copy to cache");
  xaxis->SetBinLabel(6, "Poisson calculation");
  xaxis->SetBinLabel(7, "Diagnostics");
  xaxis->SetBinLabel(8, "Total");
  d->Add(fHTiming);
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
AliFMDDensityCalculator::Print(Option_t* option) const
{
  // 
  // Print information 
  // 
  // Parameters:
  //    option Not used
  //
  AliForwardUtil::PrintTask(*this);
  gROOT->IncreaseDirLevel();

  TString phiM("none");
  switch (fUsePhiAcceptance) { 
  case kPhiNoCorrect:    phiM = "none"; break;
  case kPhiCorrectNch:   phiM = "correct Nch"; break;
  case kPhiCorrectELoss: phiM = "correct energy loss"; break;
  }

  PFV("Max(particles)",		fMaxParticles);
  PFB("Poisson method",		fUsePoisson);
  PFV("Eta lumping",		fEtaLumping);
  PFV("Phi lumping",		fPhiLumping);
  PFB("Recalculate phi",	fRecalculatePhi);
  PFB("Use phi acceptance",     phiM);
  PFV("Min(quality)",           fMinQuality);
  PFV("Threshold(hit)",         fHitThreshold);
  PFV("Max(outliers)",          fMaxOutliers);
  PFV("Cut(outlier)",           fOutlierCut);
  PFV("Lower cut", "");
  fCuts.Print();

  TString opt(option);
  opt.ToLower();
  if (!opt.Contains("nomax")) {
    PFV("Max weights", "");

    char ind[64];
    for (Int_t i = 0; i < gROOT->GetDirLevel(); i++) ind[i] = ' ';
    ind[gROOT->GetDirLevel()] = '\0';
    for (UShort_t d=1; d<=3; d++) { 
      UShort_t nr = (d == 1 ? 1 : 2);
      for (UShort_t q=0; q<nr; q++) { 
	ind[gROOT->GetDirLevel()]   = ' ';
	ind[gROOT->GetDirLevel()+1] = '\0';
	Char_t      r = (q == 0 ? 'I' : 'O');
	std::cout << ind << " FMD" << d << r << ":";
	ind[gROOT->GetDirLevel()+1] = ' ';
	ind[gROOT->GetDirLevel()+2] = '\0';
	
	const TArrayI& a = (d == 1 ? fFMD1iMax : 
			    (d == 2 ? (r == 'I' ? fFMD2iMax : fFMD2oMax) : 
			     (r == 'I' ? fFMD3iMax : fFMD3oMax)));
	Int_t j = 0;
	for (Int_t i = 0; i < a.fN; i++) { 
	  if (a.fArray[i] < 1) continue; 
	  if (j % 6 == 0)      std::cout << "\n " << ind;
	  j++;
	  std::cout << "  " << std::setw(3) << i << ": " << a.fArray[i];
	}
	std::cout << std::endl;
      }
    }
  }
  gROOT->DecreaseDirLevel();
}

//====================================================================
AliFMDDensityCalculator::RingHistos::RingHistos()
  : AliForwardUtil::RingHistos(),
    fList(0),
    // fEvsN(0), 
    // fEvsM(0), 
    // fEtaVsN(0),
    // fEtaVsM(0),
    fCorr(0),
    fSignal(0),
    fDensity(0),
    fELossVsPoisson(0),
    fDiffELossPoisson(0),
    fELossVsPoissonOut(0),
    fDiffELossPoissonOut(0),
    fOutliers(0),
    fPoisson(),
    fELoss(0),
    fELossUsed(0),
    fMultCut(0), 
    fTotal(0),
    fGood(0),
    fPhiAcc(0), 
    fPhiBefore(0),
    fPhiAfter(0),
    fEtaBefore(0),
    fEtaAfter(0)
{
  // 
  // Default CTOR
  //
}

//____________________________________________________________________
AliFMDDensityCalculator::RingHistos::RingHistos(UShort_t d, Char_t r)
  : AliForwardUtil::RingHistos(d,r),
    fList(0),
    // fEvsN(0), 
    // fEvsM(0),
    // fEtaVsN(0),
    // fEtaVsM(0),
    fCorr(0),
    fSignal(0),
    fDensity(0),
    fELossVsPoisson(0),
    fDiffELossPoisson(0),
    fELossVsPoissonOut(0),
    fDiffELossPoissonOut(0),
    fOutliers(0),
    fPoisson("ignored"),
    fELoss(0),
    fELossUsed(0),
    fMultCut(0), 
    fTotal(0),
    fGood(0),
    fPhiAcc(0), 
    fPhiBefore(0),
    fPhiAfter(0),
    fEtaBefore(0),
    fEtaAfter(0)
{
  // 
  // Constructor
  // 
  // Parameters:
  //    d detector
  //    r ring 
  //
#if 0
  fEvsN = new TH2D("elossVsNnocorr", 
		   "#Delta E/#Delta E_{mip} vs uncorrected inclusive N_{ch}",
		   250, -.5, 24.5, 251, -1.5, 24.5);
  fEvsN->SetXTitle("#Delta E/#Delta E_{mip}");
  fEvsN->SetYTitle("Inclusive N_{ch} (uncorrected)");
  fEvsN->Sumw2();
  fEvsN->SetDirectory(0);

  fEvsM = static_cast<TH2D*>(fEvsN->Clone("elossVsNcorr"));
  fEvsM->SetTitle("#Delta E/#Delta E_{mip} vs corrected inclusive N_{ch}");
  fEvsM->SetDirectory(0);

  fEtaVsN = new TProfile("etaVsNnocorr",
			 "Average inclusive N_{ch} vs #eta (uncorrected)",
			 200, -4, 6);
  fEtaVsN->SetXTitle("#eta");
  fEtaVsN->SetYTitle("#LT N_{ch,incl}#GT (uncorrected)");
  fEtaVsN->SetDirectory(0);
  fEtaVsN->SetLineColor(Color());
  fEtaVsN->SetFillColor(Color());

  fEtaVsM = static_cast<TProfile*>(fEtaVsN->Clone("etaVsNcorr"));
  fEtaVsM->SetTitle("Average inclusive N_{ch} vs #eta (corrected)");
  fEtaVsM->SetYTitle("#LT N_{ch,incl}#GT (corrected)");
  fEtaVsM->SetDirectory(0);
#endif

  fCorr = new TProfile("corr", "Average correction", 200, -4, 6);
  fCorr->SetXTitle("#eta");
  fCorr->SetYTitle("#LT correction#GT");
  fCorr->SetDirectory(0);
  fCorr->SetLineColor(Color());
  fCorr->SetFillColor(Color());

  fSignal = new TH2D("signal", "Signal distribution", 200, -4, 6, 1000, 0, 10);
  fSignal->SetXTitle("#eta");
  fSignal->SetYTitle("Signal");
  fSignal->SetDirectory(0);
  
  fDensity = new TH2D("inclDensity", "Inclusive N_{ch} density",
		      200, -4, 6, (r == 'I' || r == 'i' ? 20 : 40), 
		      0, 2*TMath::Pi());
  fDensity->SetDirectory(0);
  fDensity->Sumw2();  fDensity->SetMarkerColor(Color());
  fDensity->SetXTitle("#eta");
  fDensity->SetYTitle("#phi [radians]");
  fDensity->SetZTitle("Inclusive N_{ch} density");

  // --- Create increasing sized bins --------------------------------
  TArrayD bins;
  // bins, lowest order, higest order, return array
  const char* nchP = "N_{ch}^{Poisson}";
  const char* nchE = "N_{ch}^{#Delta}";
  AliForwardUtil::MakeLogScale(300, 0, 2, bins);
  fELossVsPoisson = new TH2D("elossVsPoisson", 
			     "N_{ch} from energy loss vs from Poisson",
			     bins.GetSize()-1, bins.GetArray(), 
			     bins.GetSize()-1, bins.GetArray());
  fELossVsPoisson->SetDirectory(0);
  fELossVsPoisson->SetXTitle(nchE);
  fELossVsPoisson->SetYTitle(nchP);
  fELossVsPoisson->SetZTitle("Correlation");
  fELossVsPoissonOut = 
    static_cast<TH2D*>(fELossVsPoisson
		       ->Clone(Form("%sOutlier", 
				    fELossVsPoisson->GetName())));
  fELossVsPoissonOut->SetDirectory(0);
  fELossVsPoissonOut->SetMarkerStyle(20);
  fELossVsPoissonOut->SetMarkerSize(0.3);
  fELossVsPoissonOut->SetMarkerColor(kBlack);
  fELossVsPoissonOut->SetTitle(Form("%s for outliers", 
				    fELossVsPoisson->GetTitle()));

  fDiffELossPoisson = new TH1D("diffElossPoisson",
			       Form("(%s-%s)/%s", nchP, nchE, nchE),
			       100, -1, 1);
  fDiffELossPoisson->SetDirectory(0);
  fDiffELossPoisson->SetXTitle(fDiffELossPoisson->GetTitle());
  fDiffELossPoisson->SetYTitle("Frequency");
  fDiffELossPoisson->SetMarkerColor(Color());
  fDiffELossPoisson->SetFillColor(Color());
  fDiffELossPoisson->SetFillStyle(3001);
  fDiffELossPoisson->Sumw2();
			       
  fDiffELossPoissonOut = 
    static_cast<TH1D*>(fDiffELossPoisson
		       ->Clone(Form("%sOutlier",fDiffELossPoisson->GetName())));
  fDiffELossPoissonOut->SetDirectory(0);
  fDiffELossPoissonOut->SetTitle(Form("%s for outliers", 
				      fDiffELossPoisson->GetTitle()));
  fDiffELossPoissonOut->SetMarkerColor(Color()-2);
  fDiffELossPoissonOut->SetFillColor(Color()-2);
  fDiffELossPoissonOut->SetFillStyle(3002);
  
  fOutliers = new TH1D("outliers", "Fraction of outliers", 100, 0, 1);
  fOutliers->SetDirectory(0);
  fOutliers->SetXTitle("N_{outlier}/(N_{outlier}+N_{inside})");
  fOutliers->SetYTitle("#sum_{events}#sum_{bins}");
  fOutliers->SetFillColor(Color());
  fOutliers->SetFillStyle(3001);
  fOutliers->SetLineColor(kBlack);

  fELoss = new TH1D("eloss", "#Delta/#Delta_{mip} in all strips", 
		    640, -1, 15);
  fELoss->SetXTitle("#Delta/#Delta_{mip} (selected)");
  fELoss->SetYTitle("P(#Delta/#Delta_{mip})");
  fELoss->SetFillColor(Color()-2);
  fELoss->SetFillStyle(3003);
  fELoss->SetLineColor(kBlack);
  fELoss->SetLineStyle(2);
  fELoss->SetLineWidth(1);
  fELoss->SetDirectory(0);

  fELossUsed = static_cast<TH1D*>(fELoss->Clone("elossUsed"));
  fELossUsed->SetTitle("#Delta/#Delta_{mip} in used strips");
  fELossUsed->SetFillStyle(3002);
  fELossUsed->SetLineStyle(1);
  fELossUsed->SetDirectory(0);

  fPhiBefore = new TH1D("phiBefore", "#phi distribution (before recalc)",
			(r == 'I' || r == 'i' ? 20 : 40), 0, 2*TMath::Pi());
  fPhiBefore->SetDirectory(0);
  fPhiBefore->SetXTitle("#phi");
  fPhiBefore->SetYTitle("Events");
  fPhiBefore->SetMarkerColor(Color());
  fPhiBefore->SetLineColor(Color());
  fPhiBefore->SetFillColor(Color());
  fPhiBefore->SetFillStyle(3001);
  fPhiBefore->SetMarkerStyle(20);

  fPhiAfter = static_cast<TH1D*>(fPhiBefore->Clone("phiAfter"));
  fPhiAfter->SetTitle("#phi distribution (after re-calc)");
  fPhiAfter->SetDirectory(0);

  fEtaBefore = new TH1D("etaBefore", "#eta distribution (before recalc)",
			200, -4, 6);
  fEtaBefore->SetDirectory(0);
  fEtaBefore->SetXTitle("#eta");
  fEtaBefore->SetYTitle("Events");
  fEtaBefore->SetMarkerColor(Color());
  fEtaBefore->SetLineColor(Color());
  fEtaBefore->SetFillColor(Color());
  fEtaBefore->SetFillStyle(3001);
  fEtaBefore->SetMarkerStyle(20);

  fEtaAfter = static_cast<TH1D*>(fEtaBefore->Clone("etaAfter"));
  fEtaAfter->SetTitle("#eta distribution (after re-calc)");
  fEtaAfter->SetDirectory(0);

}
//____________________________________________________________________
AliFMDDensityCalculator::RingHistos::RingHistos(const RingHistos& o)
  : AliForwardUtil::RingHistos(o),
    fList(o.fList), 
    // fEvsN(o.fEvsN), 
    // fEvsM(o.fEvsM),
    // fEtaVsN(o.fEtaVsN),
    // fEtaVsM(o.fEtaVsM),
    fCorr(o.fCorr),
    fSignal(o.fSignal),
    fDensity(o.fDensity),
    fELossVsPoisson(o.fELossVsPoisson),
    fDiffELossPoisson(o.fDiffELossPoisson),
    fELossVsPoissonOut(o.fELossVsPoissonOut),
    fDiffELossPoissonOut(o.fDiffELossPoissonOut),
    fOutliers(o.fOutliers),
    fPoisson(o.fPoisson),
    fELoss(o.fELoss),
    fELossUsed(o.fELossUsed),
    fMultCut(o.fMultCut), 
    fTotal(o.fTotal),
    fGood(o.fGood),
    fPhiAcc(o.fPhiAcc), 
    fPhiBefore(o.fPhiBefore),
    fPhiAfter(o.fPhiAfter),
    fEtaBefore(o.fEtaBefore),
    fEtaAfter(o.fEtaAfter)
{
  // 
  // Copy constructor 
  // 
  // Parameters:
  //    o Object to copy from 
  //
}

//____________________________________________________________________
AliFMDDensityCalculator::RingHistos&
AliFMDDensityCalculator::RingHistos::operator=(const RingHistos& o)
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
  
  // if (fEvsN)             delete fEvsN;
  // if (fEvsM)             delete fEvsM;
  // if (fEtaVsN)           delete fEtaVsN;
  // if (fEtaVsM)           delete fEtaVsM;
  if (fCorr)             delete fCorr;
  if (fSignal)           delete fSignal;
  if (fDensity)          delete fDensity;
  if (fELossVsPoisson)   delete fELossVsPoisson;
  if (fDiffELossPoisson) delete fDiffELossPoisson;
  if (fTotal)            delete fTotal;
  if (fGood)             delete fGood;
  if (fPhiAcc)           delete fPhiAcc;
  if (fPhiBefore)        delete fPhiBefore;
  if (fPhiAfter)         delete fPhiAfter;
  if (fEtaBefore)        delete fEtaBefore;
  if (fEtaAfter)         delete fEtaAfter;

  // fEvsN             = static_cast<TH2D*>(o.fEvsN->Clone());
  // fEvsM             = static_cast<TH2D*>(o.fEvsM->Clone());
  // fEtaVsN           = static_cast<TProfile*>(o.fEtaVsN->Clone());
  // fEtaVsM           = static_cast<TProfile*>(o.fEtaVsM->Clone());
  fCorr                = static_cast<TProfile*>(o.fCorr->Clone());
  fSignal              = static_cast<TH2D*>(o.fSignal->Clone());
  fDensity             = static_cast<TH2D*>(o.fDensity->Clone());
  fELossVsPoisson      = static_cast<TH2D*>(o.fELossVsPoisson->Clone());
  fDiffELossPoisson    = static_cast<TH1D*>(o.fDiffELossPoisson->Clone());
  fELossVsPoissonOut   = static_cast<TH2D*>(o.fELossVsPoisson->Clone());
  fDiffELossPoissonOut = static_cast<TH1D*>(o.fDiffELossPoisson->Clone());
  fOutliers            = static_cast<TH1D*>(o.fOutliers->Clone());
  fPoisson             = o.fPoisson;
  fELoss               = static_cast<TH1D*>(o.fELoss->Clone());
  fELossUsed           = static_cast<TH1D*>(o.fELossUsed->Clone());
  fTotal               = static_cast<TH1D*>(o.fTotal->Clone());
  fGood                = static_cast<TH1D*>(o.fGood->Clone());
  fPhiAcc              = static_cast<TH2D*>(o.fPhiAcc->Clone());
  fPhiBefore           = static_cast<TH1D*>(o.fPhiBefore->Clone());
  fPhiAfter            = static_cast<TH1D*>(o.fPhiAfter->Clone());
  fEtaBefore           = static_cast<TH1D*>(o.fEtaBefore->Clone());
  fEtaAfter            = static_cast<TH1D*>(o.fEtaAfter->Clone());
  return *this;
}
//____________________________________________________________________
AliFMDDensityCalculator::RingHistos::~RingHistos()
{
  // 
  // Destructor 
  //
}


//____________________________________________________________________
void
AliFMDDensityCalculator::RingHistos::SetupForData(const TAxis& eAxis)
{
  // Initialize 
  // This is called on first event 
  fPoisson.Init(-1,-1);
  fTotal = new TH1D("total", "Total # of strips per #eta",
		    eAxis.GetNbins(), eAxis.GetXmin(), eAxis.GetXmax());
  fTotal->SetDirectory(0);
  fTotal->SetXTitle("#eta");
  fTotal->SetYTitle("# of strips");
  fGood = static_cast<TH1D*>(fTotal->Clone("good"));
  fGood->SetTitle("# of good strips per #eta");
  fGood->SetDirectory(0);
  
  fPhiAcc = new TH2D("phiAcc", "#phi acceptance vs Ip_{z}", 
		     eAxis.GetNbins(), eAxis.GetXmin(), eAxis.GetXmax(),
		     10, -10, 10);
  fPhiAcc->SetXTitle("#eta");
  fPhiAcc->SetYTitle("v_{z} [cm]");
  fPhiAcc->SetZTitle("#phi acceptance");
  fPhiAcc->SetDirectory(0);

  if (fList) fList->Add(fPhiAcc);
}

//____________________________________________________________________
void
AliFMDDensityCalculator::RingHistos::CreateOutputObjects(TList* dir)
{
  // 
  // Make output.  This is called as part of SlaveBegin
  // 
  // Parameters:
  //    dir Where to put it 
  //
  TList* d = DefineOutputList(dir);
  // d->Add(fEvsN);
  // d->Add(fEvsM);
  // d->Add(fEtaVsN);
  // d->Add(fEtaVsM);
  d->Add(fCorr);
  d->Add(fSignal);
  d->Add(fDensity);
  d->Add(fELossVsPoisson);
  d->Add(fELossVsPoissonOut);
  d->Add(fDiffELossPoisson);
  d->Add(fDiffELossPoissonOut);
  d->Add(fOutliers);
  fPoisson.Output(d);
  fPoisson.GetOccupancy()->SetFillColor(Color());
  fPoisson.GetMean()->SetFillColor(Color());
  fPoisson.GetOccupancy()->SetFillColor(Color());
  d->Add(fELoss);
  d->Add(fELossUsed);
  d->Add(fPhiBefore);
  d->Add(fPhiAfter);
  d->Add(fEtaBefore);
  d->Add(fEtaAfter);

  TAxis x(NStrip(), -.5, NStrip()-.5);
  TAxis y(NSector(), -.5, NSector()-.5);
  x.SetTitle("strip");
  y.SetTitle("sector");
  fPoisson.Define(x, y);

  d->Add(AliForwardUtil::MakeParameter("cut", fMultCut));
  fList = d;
}

//____________________________________________________________________
void
AliFMDDensityCalculator::RingHistos::Terminate(TList* dir, Int_t nEvents)
{
  // 
  // Scale the histograms to the total number of events 
  // 
  // Parameters:
  //    dir     Where the output is 
  //    nEvents Number of events 
  //
  TList* l = GetOutputList(dir);
  if (!l) return; 

  TH2D* density = static_cast<TH2D*>(GetOutputHist(l,"inclDensity"));
  if (density) density->Scale(1./nEvents);
  fDensity = density;
}

//____________________________________________________________________
//
// EOF
//
	  


