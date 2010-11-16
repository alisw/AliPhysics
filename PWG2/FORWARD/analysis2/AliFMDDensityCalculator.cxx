#include "AliFMDDensityCalculator.h"
#include <AliESDFMD.h>
#include <TAxis.h>
#include <TList.h>
#include <TMath.h>
#include "AliFMDAnaParameters.h"
#include "AliLog.h"
#include <TH2D.h>

ClassImp(AliFMDDensityCalculator)
#if 0
; // For Emacs
#endif 

//____________________________________________________________________
AliFMDDensityCalculator::AliFMDDensityCalculator()
  : TNamed(), 
    fRingHistos(),
    fMultCut(0.3)
{}

//____________________________________________________________________
AliFMDDensityCalculator::AliFMDDensityCalculator(const char* title)
  : TNamed("fmdDensityCalculator", title), 
    fRingHistos(), 
    fMultCut(0.3)
{
  fRingHistos.SetName(GetName());
  fRingHistos.Add(new RingHistos(1, 'I'));
  fRingHistos.Add(new RingHistos(2, 'I'));
  fRingHistos.Add(new RingHistos(2, 'O'));
  fRingHistos.Add(new RingHistos(3, 'I'));
  fRingHistos.Add(new RingHistos(3, 'O'));
}

//____________________________________________________________________
AliFMDDensityCalculator::AliFMDDensityCalculator(const 
						 AliFMDDensityCalculator& o)
  : TNamed(o), 
    fRingHistos(), 
    fMultCut(o.fMultCut)
{
  TIter    next(&o.fRingHistos);
  TObject* obj = 0;
  while ((obj = next())) fRingHistos.Add(obj);
}

//____________________________________________________________________
AliFMDDensityCalculator::~AliFMDDensityCalculator()
{
  fRingHistos.Delete();
}

//____________________________________________________________________
AliFMDDensityCalculator&
AliFMDDensityCalculator::operator=(const AliFMDDensityCalculator& o)
{
  SetName(o.GetName());
  SetTitle(o.GetTitle());

  fMultCut = o.fMultCut;

  fRingHistos.Delete();
  TIter    next(&o.fRingHistos);
  TObject* obj = 0;
  while ((obj = next())) fRingHistos.Add(obj);
  
  return *this;
}

//____________________________________________________________________
AliFMDDensityCalculator::RingHistos*
AliFMDDensityCalculator::GetRingHistos(UShort_t d, Char_t r) const
{
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
Bool_t
AliFMDDensityCalculator::Calculate(const AliESDFMD&        fmd,
				   AliForwardUtil::Histos& hists,
				   Int_t                   vtxbin, 
				   Bool_t                  lowFlux)
{
  for (UShort_t d=1; d<=3; d++) { 
    UShort_t nr = (d == 1 ? 1 : 2);
    for (UShort_t q=0; q<nr; q++) { 
      Char_t      r = (q == 0 ? 'I' : 'O');
      UShort_t    ns= (q == 0 ?  20 :  40);
      UShort_t    nt= (q == 0 ? 512 : 256);
      TH2D*       h = hists.Get(d,r);
      RingHistos* rh= GetRingHistos(d,r);

      for (UShort_t s=0; s<ns; s++) { 
	for (UShort_t t=0; t<nt; t++) {
	  Float_t mult = fmd.Multiplicity(d,r,s,t);
	  
	  if (mult == 0 || mult > 20) continue;

	  Float_t phi = fmd.Phi(d,r,s,t) / 180 * TMath::Pi();
	  Float_t eta = fmd.Eta(d,r,s,t);
	  
	  Float_t n = NParticles(mult,d,r,s,t,vtxbin,eta,lowFlux);
	  rh->fEvsN->Fill(mult,n);

	  Float_t c = Correction(d,r,s,t,vtxbin,eta,lowFlux);
	  if (c > 0) n /= c;
	  rh->fEvsM->Fill(mult,n);

	  h->Fill(eta,phi,n);
	  rh->fDensity->Fill(eta,phi,n);
	} // for t
      } // for s 
    } // for q
  } // for d
  
  return kTRUE;
}

//_____________________________________________________________________
Float_t 
AliFMDDensityCalculator::NParticles(Float_t  mult, 
				    UShort_t d, 
				    Char_t   r, 
				    UShort_t /*s*/, 
				    UShort_t /*t*/, 
				    Int_t    /*v*/, 
				    Float_t  eta,
				    Bool_t   lowFlux) const
{
  if (mult <= fMultCut) return 0;
  if (lowFlux) return 1;

  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();

  Float_t mpv  = pars->GetMPV(d,r,eta);
  Float_t w    = pars->GetSigma(d,r,eta);
  Float_t w2   = pars->Get2MIPWeight(d,r,eta);
  Float_t w3   = pars->Get3MIPWeight(d,r,eta);
  Float_t mpv2 = 2*mpv+2*w*TMath::Log(2);
  Float_t mpv3 = 3*mpv+3*w*TMath::Log(3);
  
  Float_t sum  = (TMath::Landau(mult,mpv,w,kTRUE) +
		  w2 * TMath::Landau(mult,mpv2,2*w,kTRUE) + 
		  w3  * TMath::Landau(mult,mpv3,3*w,kTRUE));
  Float_t wsum = (TMath::Landau(mult,mpv,w,kTRUE) +
		  2*w2 * TMath::Landau(mult,mpv2,2*w,kTRUE) + 
		  3*w3 * TMath::Landau(mult,mpv3,3*w,kTRUE));
  
  return (sum > 0) ? wsum / sum : 1;
}

//_____________________________________________________________________
Float_t 
AliFMDDensityCalculator::Correction(UShort_t d, Char_t r, UShort_t /*s*/, 
				    UShort_t t, Int_t v, Float_t eta,
				    Bool_t lowFlux) const
{
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();

  Float_t correction = AcceptanceCorrection(r,t);
  if (lowFlux) { 
    TH1F* dblHitCor = pars->GetDoubleHitCorrection(d,r);
    if (dblHitCor) {
      Float_t dblC = dblHitCor->GetBinContent(dblHitCor->FindBin(eta));
      if (dblC > 0) 
	correction *= dblC;
    }
    else 
      AliWarning(Form("Missing double hit correction for FMD%d%c",d,r));
  }
  
  TH1F* deadCor = pars->GetFMDDeadCorrection(v);
  if (deadCor) { 
    Float_t deadC = deadCor->GetBinContent(deadCor->FindBin(eta));
    if (deadC > 0) 
      correction *= deadC; 
  }
  
  return correction;
}


//_____________________________________________________________________
Float_t 
AliFMDDensityCalculator::AcceptanceCorrection(Char_t r, UShort_t t) const
{
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  
  //AliFMDRing fmdring(ring);
  //fmdring.Init();
  Float_t   rad       = pars->GetMaxR(r)-pars->GetMinR(r);
  Float_t   slen      = pars->GetStripLength(r,t);
  Float_t   sblen     = pars->GetBaseStripLength(r,t);
  Float_t   nstrips   = (r == 'I' ? 512 : 256);
  Float_t   segment   = rad / nstrips;
  Float_t   radius    = pars->GetMinR(r) + segment*t;
  
  Float_t   basearea1 = 0.5*sblen*TMath::Power(radius,2);
  Float_t   basearea2 = 0.5*sblen*TMath::Power((radius-segment),2);
  Float_t   basearea  = basearea1 - basearea2;
  
  Float_t   area1     = 0.5*slen*TMath::Power(radius,2);
  Float_t   area2     = 0.5*slen*TMath::Power((radius-segment),2);
  Float_t   area      = area1 - area2;
  
  return area/basearea;
}

//____________________________________________________________________
void
AliFMDDensityCalculator::ScaleHistograms(Int_t nEvents)
{
  if (nEvents <= 0) return;

  TIter    next(&fRingHistos);
  RingHistos* o = 0;
  while ((o = static_cast<RingHistos*>(next()))) {
    o->fDensity->Scale(1. / nEvents);
  }
}

//____________________________________________________________________
void
AliFMDDensityCalculator::Output(TList* dir)
{
  TList* d = new TList;
  d->SetName(GetName());
  dir->Add(d);
  TIter    next(&fRingHistos);
  RingHistos* o = 0;
  while ((o = static_cast<RingHistos*>(next()))) {
    o->Output(d);
  }
}

//====================================================================
AliFMDDensityCalculator::RingHistos::RingHistos()
  : fDet(0),
    fRing('\0'),
    fEvsN(0), 
    fEvsM(0), 
    fDensity(0)
{}

//____________________________________________________________________
AliFMDDensityCalculator::RingHistos::RingHistos(UShort_t d, Char_t r)
  : fDet(d), 
    fRing(r),
    fEvsN(0), 
    fEvsM(0),
    fDensity(0)
{
  fEvsN = new TH2D(Form("FMD%d%c_Eloss_N_nocorr", d, r), 
		   Form("Energy loss vs uncorrected inclusive "
			"N_{ch} in FMD%d%c", d, r), 
		   100, -.5, 24.5, 100, -.5, 24.5);
  fEvsM = new TH2D(Form("FMD%d%c_Eloss_N_corr", d, r), 
		   Form("Energy loss vs corrected inclusive N_{ch} in FMD%d%c",
			d, r), 100, -.5, 24.5, 100, -.5, 24.5);
  fEvsN->SetXTitle("#Delta E/#Delta E_{mip}");
  fEvsN->SetYTitle("Inclusive N_{ch} (uncorrected)");
  fEvsN->Sumw2();
  fEvsN->SetDirectory(0);
  fEvsN->SetXTitle("#Delta E/#Delta E_{mip}");
  fEvsN->SetYTitle("Inclusive N_{ch} (corrected)");
  fEvsM->Sumw2();
  fEvsM->SetDirectory(0);

  fDensity = new TH2D(Form("FMD%d%c_Incl_Density", d, r), 
		      Form("in FMD%d%c", d, r), 
		      200, -4, 6, (r == 'I' || r == 'i' ? 20 : 40), 
		      0, 2*TMath::Pi());
  fDensity->SetDirectory(0);
  fDensity->SetXTitle("#eta");
  fDensity->SetYTitle("#phi [radians]");
  fDensity->SetZTitle("Inclusive N_{ch} density");
}
//____________________________________________________________________
AliFMDDensityCalculator::RingHistos::RingHistos(const RingHistos& o)
  : TObject(o), 
    fDet(o.fDet), 
    fRing(o.fRing), 
    fEvsN(o.fEvsN), 
    fEvsM(o.fEvsM),
    fDensity(o.fDensity)
{}

//____________________________________________________________________
AliFMDDensityCalculator::RingHistos&
AliFMDDensityCalculator::RingHistos::operator=(const RingHistos& o)
{
  fDet  = o.fDet;
  fRing = o.fRing;
  
  if (fEvsN)    delete  fEvsN;
  if (fEvsM)    delete  fEvsM;
  if (fDensity) delete fDensity;
  
  fEvsN    = static_cast<TH2D*>(o.fEvsN->Clone());
  fEvsM    = static_cast<TH2D*>(o.fEvsM->Clone());
  fDensity = static_cast<TH2D*>(o.fDensity->Clone());
  
  return *this;
}
//____________________________________________________________________
AliFMDDensityCalculator::RingHistos::~RingHistos()
{
  if (fEvsN)    delete fEvsN;
  if (fEvsM)    delete fEvsM;
  if (fDensity) delete fDensity;
}

//____________________________________________________________________
void
AliFMDDensityCalculator::RingHistos::Output(TList* dir)
{
  TList* d = new TList;
  d->SetName(Form("FMD%d%c", fDet, fRing)); 
  d->Add(fEvsN);
  d->Add(fEvsM);
  d->Add(fDensity);
  dir->Add(d);
}

//____________________________________________________________________
//
// EOF
//
	  


