#include "AliFMDMCTrackInspector.h"
#include "AliMCEvent.h"
#include "AliESDEvent.h"
#include "AliStack.h"
#include "AliESDFMD.h"
#include "AliGenEventHeader.h"
#include "AliHeader.h"
#include "AliLog.h"
#include "AliFMDEncodedEdx.h"
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TClonesArray.h>

//____________________________________________________________________
AliFMDMCTrackInspector::AliFMDMCTrackInspector()
  : AliFMDEnergyFitter(),
    fTracker(),
    fIp(),
    fNTrack(0),
    fNPrimary(0)
{
}

//____________________________________________________________________
AliFMDMCTrackInspector::AliFMDMCTrackInspector(const char* title)
  : AliFMDEnergyFitter(title),
    fTracker("tracker"),
    fIp(),
    fNTrack(0),
    fNPrimary(0)
{
  // Some defaults 
  fUseIncreasingBins = true;
  fDoFits            = true;
  fDoMakeObject      = false;
  fResidualMethod    = kNoResiduals;
}
//____________________________________________________________________
AliFMDMCTrackInspector::~AliFMDMCTrackInspector()
{}

//____________________________________________________________________
void
AliFMDMCTrackInspector::CreateOutputObjects(TList* dir)
{
  DGUARD(fDebug, 2, "Create output objects w/MC hits (%p)", dir);
  fTracker.CreateOutputObjects(dir);
  AliFMDEnergyFitter::CreateOutputObjects(dir);

  TIter    next(&fRingHistos);
  RingHistos* o = 0;
  while ((o = static_cast<RingHistos*>(next()))) {
    o->fBetaGammadEdx = static_cast<TH2*>(fTracker.GetBetaGammadEdx()->Clone());
    o->fBetaGammaEta  = static_cast<TH2*>(fTracker.GetBetaGammaEta()->Clone());
    o->fDedxEta       = static_cast<TH2*>(fTracker.GetDEdxEta()->Clone());
    o->fBetaGammadEdx->SetDirectory(0);
    o->fBetaGammaEta ->SetDirectory(0);
    o->fDedxEta      ->SetDirectory(0);
    o->fList->Add(o->fBetaGammadEdx);
    o->fList->Add(o->fBetaGammaEta );
    o->fList->Add(o->fDedxEta      );
  }

  // TList* d = static_cast<TList*>(l->FindObject(GetName()));
}
//____________________________________________________________________
Bool_t
AliFMDMCTrackInspector::PreEvent(const AliMCEvent& mcInput)
{
  DGUARD(fDebug, 5, "Reset for event w/MC hits (%p)", &mcInput);
  AliMCEvent&        mc        = const_cast<AliMCEvent&>(mcInput);
  AliHeader*         header    = mc.Header();
  AliStack*          stack     = mc.Stack();
  AliGenEventHeader* genHeader = header->GenEventHeader();
  
  genHeader->PrimaryVertex(fIp);
  fNTrack   = stack->GetNtrack();
  fNPrimary = stack->GetNprimary();

  return true;
}

//____________________________________________________________________
Bool_t
AliFMDMCTrackInspector::Event(const AliESDEvent& esdInput,
			      const AliMCEvent&  mcInput,
			      Double_t           cent)
{
  DGUARD(fDebug, 5, "Process an event for MC hits (%p,%p)",
	 &esdInput, &mcInput);
  PreEvent(mcInput);

  AliESDFMD*  esdFMD    = esdInput.GetFMDData();
  if (!esdFMD) return true;

  TVector3 ip(fIp[0], fIp[1], fIp[2]);
  fTracker.Calculate(*esdFMD, mcInput, ip, cent);

  return PostEvent();
}


//____________________________________________________________________
Bool_t AliFMDMCTrackInspector::PostEvent()
{
  DGUARD(fDebug, 5, "Fill MC hit energy loss");

  // AliESDFMD*  esdFMD    = input.GetFMDData();
  for (UShort_t d = 1; d <= 3; d++) { 
    UShort_t nQ = (d == 1 ? 1 : 2);
    for (UShort_t q = 0; q < nQ; q++) { 
      UShort_t    nS = (q == 0 ?  20 :  40);
      UShort_t    nT = (q == 0 ? 512 : 256);
      Char_t      r  = (q == 0 ? 'I' : 'O');
      RingHistos* h  = 
	static_cast<AliFMDMCTrackInspector::RingHistos*>(GetRingHistos(d,r));
      if (!h) continue;

      for (UShort_t s = 0; s < nS; s++) { 
	for (UShort_t t = 0; t < nT; t++) {
	  Float_t  totPrim = fTracker.GetPrimaries()(d, r, s, t);
	  Float_t  totSec  = fTracker.GetSecondaries()(d, r, s, t);
	  Float_t  totAll  = fTracker.GetAll()(d, r, s, t);
	  Double_t esdEta  = fTracker.GetEta()(d, r, s, t);
	  if (totAll  > 0) h->FillMC(0,  esdEta, totAll);
	  if (totPrim > 0) h->FillMC(1,  esdEta, totPrim);
	  if (totSec  > 0) h->FillMC(2,  esdEta, totSec);
	}
      }
    }
  }
  TClonesArray* hits = fTracker.GetHits();
  TIter next(hits);
  AliFMDMCTrackELoss::Hit* hit = 0;
  while ((hit = static_cast<AliFMDMCTrackELoss::Hit*>(next()))) {
    RingHistos* h  = 
      static_cast<AliFMDMCTrackInspector::RingHistos*>(GetRingHistos(hit->D(),
								     hit->R()));
    if (!h) continue;
    h->fBetaGammadEdx->Fill(hit->BetaGamma(), hit->DeDx());
    h->fBetaGammaEta->Fill(hit->Eta(), hit->BetaGamma());
    h->fDedxEta->Fill(hit->Eta(), hit->DeDx());
    
  }
  return true;
}

//____________________________________________________________________
AliFMDEnergyFitter::RingHistos*
AliFMDMCTrackInspector::CreateRingHistos(UShort_t d, Char_t r) const
{ 
  // DGUARD(fDebug, 1, "Create Ring cache of MC hit energy loss (FMD%d%c)",d,r);
  return new AliFMDMCTrackInspector::RingHistos(d,r);
}

//====================================================================
AliFMDMCTrackInspector::RingHistos::RingHistos()
  : AliFMDEnergyFitter::RingHistos(), 
    fPrimary(0), 
    fSecondary(0),
    fBetaGammadEdx(0),
    fBetaGammaEta(0),
    fDedxEta(0)
{}

//____________________________________________________________________
AliFMDMCTrackInspector::RingHistos::RingHistos(UShort_t d, Char_t r)
  : AliFMDEnergyFitter::RingHistos(d,r), 
    fPrimary(0),
    fSecondary(0),
    fBetaGammadEdx(0),
    fBetaGammaEta(0),
    fDedxEta(0)
{}
//____________________________________________________________________
TArrayD
AliFMDMCTrackInspector::RingHistos::MakeIncreasingAxis(Int_t, 
							Double_t,
							Double_t) const
{
  // Make an increasing axis for ELoss distributions. 
  // 
  // We use the service function of AliFMDEncodedEdx to do this 
  TArrayD ret;
  const AliFMDEncodedEdx::Spec& s = AliFMDEncodedEdx::GetdEAxis();
  s.FillBinArray(ret);

  return ret;
}
//____________________________________________________________________
void
AliFMDMCTrackInspector::RingHistos::SetupForData(const TAxis& eAxis, 
						  const TAxis& /*cAxis*/, 
						  Double_t maxDE, 
						  Int_t    nDEbins, 
						  Bool_t   useIncrBin)
{
  // AliFMDEnergyFitter::RingHistos::SetupForData(eAxis, cAxis, maxDE, nDEbins, 
  // useIncrBin);

  fHist      = Make("eloss", "#sum#Delta_{true} of all", 
		    eAxis, maxDE, nDEbins, useIncrBin);
  fPrimary   = Make("primary", "#sum#Delta_{true} of primaries", 
		    eAxis, maxDE, nDEbins, useIncrBin);
  fPrimary->SetMarkerStyle(24);
  fPrimary->SetMarkerSize(fPrimary->GetMarkerSize()*1.2);
  fSecondary = Make("secondary","#sum#Delta_{true} of secondaries",
 		    eAxis, maxDE, nDEbins, useIncrBin);
  fSecondary->SetMarkerStyle(25);
  fSecondary->SetMarkerSize(fSecondary->GetMarkerSize()*1.2);
  fHist->SetXTitle("#Delta_{true}");
  fSecondary->SetXTitle("#Delta_{true}");
  fPrimary->SetXTitle("#Delta_{true}");

  fList->Add(fHist);
  fList->Add(fPrimary);
  fList->Add(fSecondary);
}

//____________________________________________________________________
void
AliFMDMCTrackInspector::RingHistos::FillMC(UShort_t flag, Double_t eta, 
					    Double_t mult)
{
  switch (flag) { 
    // case 0: AliFMDEnergyFitter::RingHistos::Fill(false, eta, 0, mult); break;
  case 0: fHist->Fill(eta, mult); break;
  case 1: fPrimary->Fill(eta, mult); break;
  case 2: fSecondary->Fill(eta, mult); break;
  }
}

//____________________________________________________________________
void
AliFMDMCTrackInspector::RingHistos::Scale(TH1* dist) const
{
  // First scale to bin width 
  AliFMDEnergyFitter::RingHistos::Scale(dist);
  // Then smoothen the histogram 
  dist->Smooth(2);
}
  
//____________________________________________________________________
TObjArray*
AliFMDMCTrackInspector::RingHistos::Fit(TList*           dir, 
					Double_t         lowCut, 
					UShort_t         nParticles,
					UShort_t         minEntries,
					UShort_t         minusBins, 
					Double_t         relErrorCut, 
					Double_t         chi2nuCut,
					Double_t         minWeight,
					Double_t         regCut,
					EResidualMethod  residuals) const
{
  TObjArray* all  = FitSlices(dir, "eloss", lowCut, nParticles, minEntries, 
			      minusBins, relErrorCut, chi2nuCut, minWeight, 
			      regCut, residuals, false, 0);
  TObjArray* prim = FitSlices(dir, "primary", lowCut, nParticles, minEntries, 
			      minusBins, relErrorCut, chi2nuCut, minWeight, 
			      regCut, residuals, false, 0);
  TObjArray* sec  = FitSlices(dir, "secondary", lowCut, nParticles, 
			      minEntries, minusBins, relErrorCut, chi2nuCut, 
			      minWeight, regCut, residuals, false, 0);
  if (!all || !prim || !sec) {
    AliWarningF("Failed to get results for all(%p), primary(%p), and/or "
		"secondary(%p)", all, prim, sec);
    return 0;
  }
  // Already added to the sub-folders
  // dir->Add(all);
  // dir->Add(prim);
  // dir->Add(sec);

  Int_t      nPrim = prim->GetEntriesFast();
  TObjArray* out   = new TObjArray;
  for (Int_t i = 0; i < nPrim; i++) { 
    TH1* h = static_cast<TH1*>(prim->At(i));
    if (!h) continue;

    TAxis* xAxis  = h->GetXaxis();
    TH2*   hh     = 0;
    if (xAxis->GetXbins()->GetArray()) 
      hh = new TH2D(h->GetName(), h->GetTitle(), 
		    xAxis->GetNbins(), xAxis->GetXbins()->GetArray(),
		    3, 0, 3);
    else 
      hh = new TH2D(h->GetName(), h->GetTitle(), 
		    xAxis->GetNbins(), xAxis->GetXmin(), xAxis->GetXmax(),
		    3, 0, 3);
    hh->GetYaxis()->SetBinLabel(1, "Primary");
    hh->GetYaxis()->SetBinLabel(2, "Secondary");
    hh->GetYaxis()->SetBinLabel(3, "All");
    hh->GetXaxis()->SetTitle("#eta");
    out->Add(hh);
  }
  for (Int_t i = 0; i < nPrim; i++) { 
    TH2* h = static_cast<TH2*>(out->At(i));
    if (!h) continue;
    
    TH1* hp = static_cast<TH1*>(prim->At(i));
    TH1* hs = static_cast<TH1*>(sec->At(i));
    TH1* ha = static_cast<TH1*>(all->At(i));
    TH1* hh[] = { hp, hs, ha };
    for (Int_t j = 0; j < 3; j++) { 
      TH1* ph = hh[j];
      if (!ph) continue;

      for (Int_t k = 1; k <= ph->GetNbinsX(); k++) { 
	Double_t c = ph->GetBinContent(k);
	Double_t e = ph->GetBinError(k);
	h->SetBinContent(k, j+1, c);
	h->SetBinError(k, j+1, e);
      }
    }
  }
  TList* l = GetOutputList(dir);
  if (!l) return 0; 

  out->SetName("comparative");
  out->SetOwner();
  l->Add(out);
  return out;
}
// 
// EOF
//
