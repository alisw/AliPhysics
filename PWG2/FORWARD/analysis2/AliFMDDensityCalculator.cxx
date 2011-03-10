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
#include <TH2D.h>
#include <TProfile.h>
#include <TROOT.h>
#include <iostream>
#include <iomanip>

ClassImp(AliFMDDensityCalculator)
#if 0
; // For Emacs
#endif 

//____________________________________________________________________
AliFMDDensityCalculator::AliFMDDensityCalculator()
  : TNamed(), 
    fRingHistos(),
    fMultCut(0),
    fSumOfWeights(0),
    fWeightedSum(0),
    fCorrections(0),
    fMaxParticles(5),
    fUsePoisson(false),
    fAccI(0),
    fAccO(0),
    fFMD1iMax(0),
    fFMD2iMax(0),
    fFMD2oMax(0),
    fFMD3iMax(0),
    fFMD3oMax(0),
    fDebug(0)
{
  // 
  // Constructor 
  //
   
}

//____________________________________________________________________
AliFMDDensityCalculator::AliFMDDensityCalculator(const char* title)
  : TNamed("fmdDensityCalculator", title), 
    fRingHistos(), 
    fMultCut(0),
    fSumOfWeights(0),
    fWeightedSum(0),
    fCorrections(0),
    fMaxParticles(5),
    fUsePoisson(false),
    fAccI(0),
    fAccO(0),
    fFMD1iMax(0),
    fFMD2iMax(0),
    fFMD2oMax(0),
    fFMD3iMax(0),
    fFMD3oMax(0),
    fDebug(0)
{
  // 
  // Constructor 
  // 
  // Parameters:
  //    name Name of object
  //
  fRingHistos.SetName(GetName());
  fRingHistos.Add(new RingHistos(1, 'I'));
  fRingHistos.Add(new RingHistos(2, 'I'));
  fRingHistos.Add(new RingHistos(2, 'O'));
  fRingHistos.Add(new RingHistos(3, 'I'));
  fRingHistos.Add(new RingHistos(3, 'O'));
  fSumOfWeights = new TH1D("sumOfWeights", "Sum of Landau weights",
			   200, 0, 20);
  fSumOfWeights->SetFillColor(kRed+1);
  fSumOfWeights->SetXTitle("#sum_{i} a_{i} f_{i}(#Delta)");
  fWeightedSum  = new TH1D("weightedSum", "Weighted sum of Landau propability",
			   200, 0, 20);
  fWeightedSum->SetFillColor(kBlue+1);
  fWeightedSum->SetXTitle("#sum_{i} i a_{i} f_{i}(#Delta)");
  fCorrections  = new TH1D("corrections", "Distribution of corrections", 
			   100, 0, 10);
  fCorrections->SetFillColor(kBlue+1);
  fCorrections->SetXTitle("correction");

  fAccI = GenerateAcceptanceCorrection('I');
  fAccO = GenerateAcceptanceCorrection('O');
}

//____________________________________________________________________
AliFMDDensityCalculator::AliFMDDensityCalculator(const 
						 AliFMDDensityCalculator& o)
  : TNamed(o), 
    fRingHistos(), 
    fMultCut(o.fMultCut),
    fSumOfWeights(o.fSumOfWeights),
    fWeightedSum(o.fWeightedSum),
    fCorrections(o.fCorrections),
    fMaxParticles(o.fMaxParticles),
    fUsePoisson(o.fUsePoisson),
    fAccI(o.fAccI),
    fAccO(o.fAccO),
    fFMD1iMax(o.fFMD1iMax),
    fFMD2iMax(o.fFMD2iMax),
    fFMD2oMax(o.fFMD2oMax),
    fFMD3iMax(o.fFMD3iMax),
    fFMD3oMax(o.fFMD3oMax),
    fDebug(o.fDebug)
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
AliFMDDensityCalculator::~AliFMDDensityCalculator()
{
  // 
  // Destructor 
  //
  fRingHistos.Delete();
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
  TNamed::operator=(o);

  fMultCut      = o.fMultCut;
  fDebug        = o.fDebug;
  fMaxParticles = o.fMaxParticles;
  fUsePoisson   = o.fUsePoisson;
  fAccI         = o.fAccI;
  fAccO         = o.fAccO;
  fFMD1iMax     = o.fFMD1iMax;
  fFMD2iMax     = o.fFMD2iMax;
  fFMD2oMax     = o.fFMD2oMax;
  fFMD3iMax     = o.fFMD3iMax;
  fFMD3oMax     = o.fFMD3oMax;

  fRingHistos.Delete();
  TIter    next(&o.fRingHistos);
  TObject* obj = 0;
  while ((obj = next())) fRingHistos.Add(obj);
  
  return *this;
}

//____________________________________________________________________
void
AliFMDDensityCalculator::Init(const TAxis& axis)
{
  // Intialize this sub-algorithm 
  //
  // Parameters:
  //   etaAxis   Not used
  CacheMaxWeights();

  TIter    next(&fRingHistos);
  RingHistos* o = 0;
  while ((o = static_cast<RingHistos*>(next())))
    o->Init(axis);
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
  if (idx < 0 || idx >= fRingHistos.GetEntries()) return 0;
  
  return static_cast<RingHistos*>(fRingHistos.At(idx));
}
    
//____________________________________________________________________
Double_t
AliFMDDensityCalculator::GetMultCut() const
{
  // 
  // Get the multiplicity cut.  If the user has set fMultCut (via
  // SetMultCut) then that value is used.  If not, then the lower
  // value of the fit range for the enery loss fits is returned.
  // 
  // Return:
  //    Lower cut on multiplicity
  //
  if (fMultCut > 0) return fMultCut;

  AliForwardCorrectionManager&  fcm = AliForwardCorrectionManager::Instance();
  AliFMDCorrELossFit* fits = fcm.GetELossFit();
  return fits->GetLowCut();
}
  
//____________________________________________________________________
Bool_t
AliFMDDensityCalculator::Calculate(const AliESDFMD&        fmd,
				   AliForwardUtil::Histos& hists,
				   UShort_t                vtxbin, 
				   Bool_t                  lowFlux)
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
  //
  for (UShort_t d=1; d<=3; d++) { 
    UShort_t nr = (d == 1 ? 1 : 2);
    for (UShort_t q=0; q<nr; q++) { 
      Char_t      r = (q == 0 ? 'I' : 'O');
      UShort_t    ns= (q == 0 ?  20 :  40);
      UShort_t    nt= (q == 0 ? 512 : 256);
      TH2D*       h = hists.Get(d,r);
      RingHistos* rh= GetRingHistos(d,r);
      rh->fEmptyStrips->Reset();
      rh->fTotalStrips->Reset();
      rh->fBasicHits->Reset();
      
      for (UShort_t s=0; s<ns; s++) { 
	for (UShort_t t=0; t<nt; t++) {
	  Float_t mult = fmd.Multiplicity(d,r,s,t);
	  Float_t phi  = fmd.Phi(d,r,s,t) / 180 * TMath::Pi();
	  Float_t eta  = fmd.Eta(d,r,s,t);
	  rh->fTotalStrips->Fill(eta, phi);
	  if (mult < GetMultCut() || mult > 20) rh->fEmptyStrips->Fill(eta,phi);
	  if (mult == AliESDFMD::kInvalidMult || mult > 20) continue;
	  
	  //if (mult == 0) { 
	  //  rh->fEmptyStrips->Fill(eta,phi);
	  //  continue;
	  // }
	  
	  Float_t n   = NParticles(mult,d,r,s,t,vtxbin,eta,lowFlux);
	  
	  rh->fEvsN->Fill(mult,n);
	  rh->fEtaVsN->Fill(eta, n);
	  
	  Float_t c = Correction(d,r,s,t,vtxbin,eta,lowFlux);
	  fCorrections->Fill(c);
	  if (c > 0) n /= c;
	  rh->fEvsM->Fill(mult,n);
	  rh->fEtaVsM->Fill(eta, n);
	  rh->fCorr->Fill(eta, c);
	  
	  //if (n < .9) rh->fEmptyStrips->Fill(eta,phi);
	  if (n > 0.9) rh->fBasicHits->Fill(eta,phi);
	    
	  h->Fill(eta,phi,n);
	  rh->fDensity->Fill(eta,phi,n);
	} // for t
      } // for s 
      //std::cout<<rh->fTotalStrips->GetEntries()<<"  "<<rh->fEmptyStrips->GetEntries()<<std::endl;
      // Loop over poisson histograms 
      for (Int_t ieta = 1; ieta <= h->GetNbinsX(); ieta++) { 
	for (Int_t iphi = 1; iphi <= h->GetNbinsY(); iphi++) { 
	  Double_t eLossV   = h->GetBinContent(ieta, iphi);
	  // Double_t eLossE   = h->GetBinError(ieta, iphi);
	  Float_t eta       = h->GetXaxis()->GetBinCenter(ieta);
	  Float_t phi       = h->GetYaxis()->GetBinCenter(iphi);
	  Double_t empty    = rh->fEmptyStrips->GetBinContent(rh->fEmptyStrips->GetXaxis()->FindBin(eta),
							      rh->fEmptyStrips->GetYaxis()->FindBin(phi)) ;
	  Double_t total    = rh->fTotalStrips->GetBinContent(rh->fTotalStrips->GetXaxis()->FindBin(eta),
							      rh->fTotalStrips->GetYaxis()->FindBin(phi));
	  
	  Double_t corr     = rh->fCorr->GetBinContent(rh->fCorr->GetXaxis()->FindBin(eta));
	  
	  
	  Double_t hits     = rh->fBasicHits->GetBinContent(ieta,iphi);
	  
	  Double_t poissonV =  (total <= 0 || empty <= 0 ? 0 : 
				-TMath::Log(empty / total));
	  if(poissonV > 0)
	    poissonV = (hits * poissonV) / (1 - TMath::Exp(-1*poissonV));
	  if(corr > 0) {
	    poissonV = poissonV / corr;
	  }
	  //else poissonV = 0;
	  // if ( eLossV > 0)
	  //  std::cout<<"event : "<<total<<"  "<<empty<<"  "<<hits<<"  "<<poissonV<<"  "<<eLossV<<std::endl;
	  Double_t poissonE = 0 ;
	  if(poissonE > 0) poissonE = TMath::Sqrt(poissonV);
	  

	  
	  rh->fELossVsPoisson->Fill(eLossV, poissonV);
	  rh->fEmptyVsTotal->Fill(total, empty);
	  if (fUsePoisson) {
	    h->SetBinContent(ieta,iphi,poissonV);
	    h->SetBinError(ieta,iphi,poissonE);
	  }
	}
      }
	
    } // for q
  } // for d
  
  return kTRUE;
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
  AliFMDCorrELossFit::ELossFit* fit = cor->GetFit(d,r,iEta);
  if (!fit) { 
    // AliWarning(Form("No energy loss fit for FMD%d%c at eta=%f", d, r, eta));
    return -1;
  }
  return fit->FindMaxWeight();
}
  
//_____________________________________________________________________
void
AliFMDDensityCalculator::CacheMaxWeights()
{
  // 
  // Find the max weights and cache them 
  // 
  AliForwardCorrectionManager&  fcm = AliForwardCorrectionManager::Instance();
  AliFMDCorrELossFit*           cor = fcm.GetELossFit();
  const TAxis&                  eta = cor->GetEtaAxis();

  Int_t nEta = eta.GetNbins();
  fFMD1iMax.Set(nEta);
  fFMD2iMax.Set(nEta);
  fFMD2oMax.Set(nEta);
  fFMD3iMax.Set(nEta);
  fFMD3oMax.Set(nEta);
  
  for (Int_t i = 0; i < nEta; i++) {
    fFMD1iMax[i] = FindMaxWeight(cor, 1, 'I', i+1);
    fFMD2iMax[i] = FindMaxWeight(cor, 2, 'I', i+1);
    fFMD2oMax[i] = FindMaxWeight(cor, 2, 'O', i+1);
    fFMD3iMax[i] = FindMaxWeight(cor, 3, 'I', i+1);
    fFMD3oMax[i] = FindMaxWeight(cor, 3, 'O', i+1);
  }
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
				    UShort_t /*s*/, 
				    UShort_t /*t*/, 
				    UShort_t /*v*/, 
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
  if (mult <= GetMultCut()) return 0;
  if (lowFlux) return 1;
  
  AliForwardCorrectionManager&  fcm = AliForwardCorrectionManager::Instance();
  AliFMDCorrELossFit::ELossFit* fit = fcm.GetELossFit()->FindFit(d,r,eta);
  if (!fit) { 
    AliWarning(Form("No energy loss fit for FMD%d%c at eta=%f", d, r, eta));
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
				    UShort_t /*s*/, 
				    UShort_t t, 
				    UShort_t /*v*/, 
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
  AliForwardCorrectionManager&  fcm = AliForwardCorrectionManager::Instance();

  Float_t correction = AcceptanceCorrection(r,t);
  if (lowFlux) { 
    TH1D* dblHitCor = 0;
    if (fcm.GetDoubleHit()) 
      dblHitCor = fcm.GetDoubleHit()->GetCorrection(d,r);

    if (dblHitCor) {
      Double_t dblC = dblHitCor->GetBinContent(dblHitCor->FindBin(eta));
      if (dblC > 0) correction *= dblC;
    }
    else {
      AliWarning(Form("Missing double hit correction for FMD%d%c",d,r));
    }
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
AliFMDDensityCalculator::ScaleHistograms(const TList* dir, Int_t nEvents)
{
  // 
  // Scale the histograms to the total number of events 
  // 
  // Parameters:
  //    dir     where to put the output
  //    nEvents Number of events 
  //
  if (nEvents <= 0) return;
  TList* d = static_cast<TList*>(dir->FindObject(GetName()));
  if (!d) return;

  TIter    next(&fRingHistos);
  RingHistos* o = 0;
  while ((o = static_cast<RingHistos*>(next())))
    o->ScaleHistograms(d, nEvents);
}

//____________________________________________________________________
void
AliFMDDensityCalculator::DefineOutput(TList* dir)
{
  // 
  // Output diagnostic histograms to directory 
  // 
  // Parameters:
  //    dir List to write in
  //  
  TList* d = new TList;
  d->SetName(GetName());
  dir->Add(d);
  d->Add(fWeightedSum);
  d->Add(fSumOfWeights);
  d->Add(fCorrections);
  d->Add(fAccI);
  d->Add(fAccO);

  TIter    next(&fRingHistos);
  RingHistos* o = 0;
  while ((o = static_cast<RingHistos*>(next()))) {
    o->Output(d);
  }
}
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
  char ind[gROOT->GetDirLevel()+3];
  for (Int_t i = 0; i < gROOT->GetDirLevel(); i++) ind[i] = ' ';
  ind[gROOT->GetDirLevel()] = '\0';
  std::cout << ind << "AliFMDDensityCalculator: " << GetName() << '\n'
	    << ind << " Multiplicity cut:       " << fMultCut << '\n'
	    << ind << " Max(particles):         " << fMaxParticles 
	    << std::endl;
  TString opt(option);
  opt.ToLower();
  if (opt.Contains("nomax")) return;
  
  std::cout << ind << " Max weights:\n";

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

//====================================================================
AliFMDDensityCalculator::RingHistos::RingHistos()
  : AliForwardUtil::RingHistos(),
    fEvsN(0), 
    fEvsM(0), 
    fEtaVsN(0),
    fEtaVsM(0),
    fCorr(0),
    fDensity(0),
    fELossVsPoisson(0),
    fTotalStrips(0),
    fEmptyStrips(0),
    fBasicHits(0),
    fEmptyVsTotal(0)
{
  // 
  // Default CTOR
  //
}

//____________________________________________________________________
AliFMDDensityCalculator::RingHistos::RingHistos(UShort_t d, Char_t r)
  : AliForwardUtil::RingHistos(d,r),
    fEvsN(0), 
    fEvsM(0),
    fEtaVsN(0),
    fEtaVsM(0),
    fCorr(0),
    fDensity(0),
    fELossVsPoisson(0),
    fTotalStrips(0),
    fEmptyStrips(0),
    fBasicHits(0),
    fEmptyVsTotal(0)
{
  // 
  // Constructor
  // 
  // Parameters:
  //    d detector
  //    r ring 
  //
  fEvsN = new TH2D(Form("%s_Eloss_N_nocorr", fName.Data()), 
		   Form("#Delta E/#Delta E_{mip} vs uncorrected inclusive "
			"N_{ch} in %s", fName.Data()), 
		   2500, -.5, 24.5, 2500, -.5, 24.5);
  fEvsM = new TH2D(Form("%s_Eloss_N_corr", fName.Data()), 
		   Form("#Delta E/#Delta E_{mip} vs corrected inclusive "
			"N_{ch} in %s", fName.Data()), 
		   2500, -.5, 24.5, 2500, -.5, 24.5);
  fEvsN->SetXTitle("#Delta E/#Delta E_{mip}");
  fEvsN->SetYTitle("Inclusive N_{ch} (uncorrected)");
  fEvsN->Sumw2();
  fEvsN->SetDirectory(0);
  fEvsM->SetXTitle("#Delta E/#Delta E_{mip}");
  fEvsM->SetYTitle("Inclusive N_{ch} (corrected)");
  fEvsM->Sumw2();
  fEvsM->SetDirectory(0);

  fEtaVsN = new TProfile(Form("%s_Eta_N_nocorr", fName.Data()),
			 Form("Average inclusive N_{ch} vs #eta (uncorrected) "
			      "in %s", fName.Data()), 200, -4, 6);
  fEtaVsM = new TProfile(Form("%s_Eta_N_corr", fName.Data()),
			 Form("Average inclusive N_{ch} vs #eta (corrected) "
			      "in %s", fName.Data()), 200, -4, 6);
  fEtaVsN->SetXTitle("#eta");
  fEtaVsN->SetYTitle("#LT N_{ch,incl}#GT (uncorrected)");
  fEtaVsN->SetDirectory(0);
  fEtaVsN->SetLineColor(Color());
  fEtaVsN->SetFillColor(Color());
  fEtaVsM->SetXTitle("#eta");
  fEtaVsM->SetYTitle("#LT N_{ch,incl}#GT (corrected)");
  fEtaVsM->SetDirectory(0);
  fEtaVsM->SetLineColor(Color());
  fEtaVsM->SetFillColor(Color());


  fCorr = new TProfile(Form("%s_corr", fName.Data()),
			 Form("Average correction in %s", fName.Data()), 
		       200, -4, 6);
  fCorr->SetXTitle("#eta");
  fCorr->SetYTitle("#LT correction#GT");
  fCorr->SetDirectory(0);
  fCorr->SetLineColor(Color());
  fCorr->SetFillColor(Color());

  fDensity = new TH2D(Form("%s_Incl_Density", fName.Data()), 
		      Form("Inclusive N_{ch} density in %s", fName.Data()), 
		      200, -4, 6, (r == 'I' || r == 'i' ? 20 : 40), 
		      0, 2*TMath::Pi());
  fDensity->SetDirectory(0);
  fDensity->SetXTitle("#eta");
  fDensity->SetYTitle("#phi [radians]");
  fDensity->SetZTitle("Inclusive N_{ch} density");

  fELossVsPoisson = new TH2D(Form("%s_eloss_vs_poisson", fName.Data()),
			     Form("N_{ch} from energy loss vs from Poission %s",
				  fName.Data()), 100, 0, 20, 100, 0, 20);
  fELossVsPoisson->SetDirectory(0);
  fELossVsPoisson->SetXTitle("N_{ch} from #DeltaE");
  fELossVsPoisson->SetYTitle("N_{ch} from Poisson");
  fELossVsPoisson->SetZTitle("Correlation");

  fEmptyVsTotal = new TH2D(Form("%s_empty_vs_total", fName.Data()), 
			   Form("# of empty strips vs. total @ # strips in %s", 
				fName.Data()), 21, -.5, 20.5, 21, -0.5, 20.5);
  fEmptyVsTotal->SetDirectory(0);
  fEmptyVsTotal->SetXTitle("Total # strips");
  fEmptyVsTotal->SetYTitle("Empty # strips");
  fEmptyVsTotal->SetZTitle("Correlation");
  
  //Inserted by HHD
  
  //Float_t nbinlimitsFMD3I[8] =  {-3.4,-3.2,-3,-2.8,-2.6,-2.4,-2.2,-2};
  Int_t nbins = 40;
  /*if(fName.Contains("FMD1I")) nbins = 7;
  if(fName.Contains("FMD2I")) nbins = 8;
  if(fName.Contains("FMD2O")) nbins = 4;
  if(fName.Contains("FMD3I")) nbins = 8;
  if(fName.Contains("FMD3O")) nbins = 4;*/
  Float_t lowlimit = -4, highlimit = 6;
  /*if(fName.Contains("FMD1I")) {lowlimit = 3.6; highlimit = 5;}
  if(fName.Contains("FMD2I")) {lowlimit = 2.2; highlimit = 3.8;}
  if(fName.Contains("FMD2O")) {lowlimit = 1.6; highlimit = 2.4;} 
  if(fName.Contains("FMD3I")) {lowlimit = -2.4; highlimit = -1.6;}
  if(fName.Contains("FMD3O")) {lowlimit = -3.5; highlimit = -2.1;} 
  
  std::cout<<nbins<<"   "<<lowlimit<<"    "<<highlimit<<std::endl;
  */
  fTotalStrips = new TH2D(Form("total%s", fName.Data()), 
			  Form("Total number of strips in %s", fName.Data()), 
			  nbins, 
			  lowlimit,
			  highlimit, 
			  (fRing == 'I' || fRing == 'i' ? 5 : 10), 
			  0, 2*TMath::Pi());
  fEmptyStrips = new TH2D(Form("empty%s", fName.Data()), 
			  Form("Empty number of strips in %s", fName.Data()), 
			  nbins, 
			  lowlimit,
			  highlimit, 
			  (fRing == 'I' || fRing == 'i' ? 5 : 10), 
			  0, 2*TMath::Pi());
  fBasicHits   = new TH2D(Form("basichits%s", fName.Data()), 
			  Form("Basic number of hits in %s", fName.Data()), 
			  200, 
			  -4, 
			  6, 
			  (fRing == 'I' || fRing == 'i' ? 20 : 40), 
			  0, 2*TMath::Pi());
  
  fTotalStrips->SetDirectory(0);
  fEmptyStrips->SetDirectory(0);
  fBasicHits->SetDirectory(0);
  fTotalStrips->SetXTitle("#eta");
  fEmptyStrips->SetXTitle("#eta");
  fBasicHits->SetXTitle("#eta");
  fTotalStrips->SetYTitle("#varphi [radians]");
  fEmptyStrips->SetYTitle("#varphi [radians]");
  fBasicHits->SetYTitle("#varphi [radians]");
  fTotalStrips->Sumw2();
  fEmptyStrips->Sumw2();
  fBasicHits->Sumw2();
  
}
//____________________________________________________________________
AliFMDDensityCalculator::RingHistos::RingHistos(const RingHistos& o)
  : AliForwardUtil::RingHistos(o), 
    fEvsN(o.fEvsN), 
    fEvsM(o.fEvsM),
    fEtaVsN(o.fEtaVsN),
    fEtaVsM(o.fEtaVsM),
    fCorr(o.fCorr),
    fDensity(o.fDensity),
    fELossVsPoisson(o.fELossVsPoisson),
    fTotalStrips(o.fTotalStrips),
    fEmptyStrips(o.fEmptyStrips),
    fBasicHits(o.fBasicHits),
    fEmptyVsTotal(o.fEmptyVsTotal)
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
  AliForwardUtil::RingHistos::operator=(o);
  
  if (fEvsN)           delete fEvsN;
  if (fEvsM)           delete fEvsM;
  if (fEtaVsN)         delete fEtaVsN;
  if (fEtaVsM)         delete fEtaVsM;
  if (fCorr)           delete fCorr;
  if (fDensity)        delete fDensity;
  if (fELossVsPoisson) delete fELossVsPoisson;
  if (fTotalStrips)    delete fTotalStrips;
  if (fEmptyStrips)    delete fEmptyStrips;
  
  fEvsN           = static_cast<TH2D*>(o.fEvsN->Clone());
  fEvsM           = static_cast<TH2D*>(o.fEvsM->Clone());
  fEtaVsN         = static_cast<TProfile*>(o.fEtaVsN->Clone());
  fEtaVsM         = static_cast<TProfile*>(o.fEtaVsM->Clone());
  fCorr           = static_cast<TProfile*>(o.fCorr->Clone());
  fDensity        = static_cast<TH2D*>(o.fDensity->Clone());
  fELossVsPoisson = static_cast<TH2D*>(o.fELossVsPoisson);
  fTotalStrips    = static_cast<TH2D*>(o.fTotalStrips);
  fEmptyStrips    = static_cast<TH2D*>(o.fEmptyStrips);
  
  return *this;
}
//____________________________________________________________________
AliFMDDensityCalculator::RingHistos::~RingHistos()
{
  // 
  // Destructor 
  //
  if (fEvsN)           delete fEvsN;
  if (fEvsM)           delete fEvsM;
  if (fEtaVsN)         delete fEtaVsN;
  if (fEtaVsM)         delete fEtaVsM;
  if (fCorr)           delete fCorr;
  if (fDensity)        delete fDensity;
  if (fELossVsPoisson) delete fELossVsPoisson;
  if (fTotalStrips)    delete fTotalStrips;
  if (fEmptyStrips)    delete fEmptyStrips;
}

  //____________________________________________________________________
void
AliFMDDensityCalculator::RingHistos::Init(const TAxis& /*eAxis*/)
{
  /*if (fTotalStrips) delete fTotalStrips;
  if (fEmptyStrips) delete fEmptyStrips;
  if (fBasicHits) delete fBasicHits;
  
  fTotalStrips = new TH2D(Form("total%s", fName.Data()), 
			  Form("Total number of strips in %s", fName.Data()), 
			  eAxis.GetNbins()/5., 
			  eAxis.GetXmin(),
			  eAxis.GetXmax(), 
			  (fRing == 'I' || fRing == 'i' ? 5 : 5), 
			  0, 2*TMath::Pi());
  fEmptyStrips = new TH2D(Form("empty%s", fName.Data()), 
			  Form("Empty number of strips in %s", fName.Data()), 
			  eAxis.GetNbins()/5., 
			  eAxis.GetXmin(), 
			  eAxis.GetXmax(), 
			  (fRing == 'I' || fRing == 'i' ? 5 : 10), 
			  0, 2*TMath::Pi());
  fBasicHits   = new TH2D(Form("basichits%s", fName.Data()), 
			  Form("Basic number of hits in %s", fName.Data()), 
			  eAxis.GetNbins(), 
			  eAxis.GetXmin(), 
			  eAxis.GetXmax(), 
			  (fRing == 'I' || fRing == 'i' ? 20 : 40), 
			  0, 2*TMath::Pi());
  
  fTotalStrips->SetDirectory(0);
  fEmptyStrips->SetDirectory(0);
  fBasicHits->SetDirectory(0);
  fTotalStrips->SetXTitle("#eta");
  fEmptyStrips->SetXTitle("#eta");
  fBasicHits->SetXTitle("#eta");
  fTotalStrips->SetYTitle("#varphi [radians]");
  fEmptyStrips->SetYTitle("#varphi [radians]");
  fBasicHits->SetYTitle("#varphi [radians]");
  fTotalStrips->Sumw2();
  fEmptyStrips->Sumw2();
  fBasicHits->Sumw2();*/
}

//____________________________________________________________________
void
AliFMDDensityCalculator::RingHistos::Output(TList* dir)
{
  // 
  // Make output 
  // 
  // Parameters:
  //    dir Where to put it 
  //
  TList* d = DefineOutputList(dir);
  d->Add(fEvsN);
  d->Add(fEvsM);
  d->Add(fEtaVsN);
  d->Add(fEtaVsM);
  d->Add(fCorr);
  d->Add(fDensity);
  d->Add(fELossVsPoisson);
  d->Add(fEmptyVsTotal);
  d->Add(fTotalStrips);
  d->Add(fEmptyStrips);
  d->Add(fBasicHits);
  
  
}

//____________________________________________________________________
void
AliFMDDensityCalculator::RingHistos::ScaleHistograms(TList* dir, Int_t nEvents)
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

  TH1* density = GetOutputHist(l,Form("%s_Incl_Density", fName.Data()));
  if (density) density->Scale(1./nEvents);
}

//____________________________________________________________________
//
// EOF
//
	  


