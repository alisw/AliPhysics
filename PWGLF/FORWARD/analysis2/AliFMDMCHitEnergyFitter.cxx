#include "AliFMDMCHitEnergyFitter.h"
#include "AliMCEvent.h"
#include "AliESDEvent.h"
#include "AliStack.h"
#include "AliFMDHit.h"
#include "AliESDFMD.h"
#include "AliMCAuxHandler.h"
#include "AliGenEventHeader.h"
#include "AliLog.h"
#include "AliFMDEncodedEdx.h"
#include <TH1.h>
#include <TH2.h>
#include <TNtuple.h>
#include <TMath.h>


//____________________________________________________________________
AliFMDMCHitEnergyFitter::AliFMDMCHitEnergyFitter()
  : AliFMDEnergyFitter(),
    fSumPrimary(0),
    fSumSecondary(0),
    fIp(),
    fNTrack(0),
    fNPrimary(0),
    fTuple(0)
{
}

//____________________________________________________________________
AliFMDMCHitEnergyFitter::AliFMDMCHitEnergyFitter(const char* title,
						 Bool_t useTuple)
  : AliFMDEnergyFitter(title),
    fSumPrimary(0),
    fSumSecondary(0),
    fIp(),
    fNTrack(0),
    fNPrimary(0),
    fTuple(0)
{
  // Some defaults 
  fUseIncreasingBins = true;
  fDoFits            = true;
  fDoMakeObject      = false;
  fResidualMethod    = kNoResiduals;

  if (!useTuple) return;

  fTuple = new TNtuple("hits", "Hits", 
		       "primary:eta:phi:edep:length:"
		       "detector:ring:sector:strip:ipz");

}
//____________________________________________________________________
AliFMDMCHitEnergyFitter::~AliFMDMCHitEnergyFitter()
{}

//____________________________________________________________________
void
AliFMDMCHitEnergyFitter::CreateOutputObjects(TList* dir)
{
  DGUARD(fDebug, 2, "Create output objects w/MC hits (%p)", dir);
  AliFMDEnergyFitter::CreateOutputObjects(dir);
  // TList* d = static_cast<TList*>(l->FindObject(GetName()));
}
//____________________________________________________________________
Bool_t
AliFMDMCHitEnergyFitter::PreEvent(const AliMCEvent& mcInput)
{
  DGUARD(fDebug, 5, "Reset for event w/MC hits (%p)", &mcInput);
  fSumPrimary.Reset(0);
  fSumSecondary.Reset(0);

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
AliFMDMCHitEnergyFitter::Event(const AliESDEvent& esdInput,
			       const AliMCEvent&  mcInput,
			       AliMCAuxHandler&   handler)
{
  DGUARD(fDebug, 5, "Process an event for MC hits (%p,%p,%p)",
	 &esdInput, &mcInput, &handler);
  PreEvent(mcInput);

  Int_t  nEntries = handler.GetNEntry();
  DMSG(fDebug,5, "We got a total of %d particles", nEntries);

  for (Int_t i = 0; i < nEntries; i++) {
    TClonesArray* hits = handler.GetEntryArray(i);
    if (!hits) continue;

    AccumulateHits(mcInput, *hits);
  }
  return PostEvent(esdInput);
}

//____________________________________________________________________
Bool_t
AliFMDMCHitEnergyFitter::AccumulateHits(const AliMCEvent& mcInput,
					const TClonesArray& hits)
{
  DGUARD(fDebug, 5, "Accumulate MC hit energy loss (%p,%p,%d)",
	 &mcInput, &hits, hits.GetEntries());

  AliStack*  stack = const_cast<AliMCEvent&>(mcInput).Stack();
  Float_t    cache[10];
  Int_t      nHit  = hits.GetEntries();
  Int_t      oldTr = -1;
  UShort_t   oldD  = 0;
  Char_t     oldR  = '\0';
  for (Int_t j = 0; j < nHit; j++) { 
	
    // Check the hit 
    AliFMDHit* hit  = static_cast<AliFMDHit*>(hits.At(j));
    if (!hit) continue;

    // Zero fill array 
    if (fTuple) for (Int_t k = 0; k < 10; k++) cache[k] = 0;

    // Get the track number 
    Int_t iTr = hit->Track();
    if (iTr < 0 || iTr >= fNTrack) 
      AliWarningF("Track # %d seems out of bounds [0,%d]", 
		  iTr, fNTrack-1);

    // Get the track 
    AliMCParticle* particle = 
      static_cast<AliMCParticle*>(mcInput.GetTrack(iTr));
	
    // Particle type 0: unknown, 1: primary, 2: secondary
    // Int_t flag = 0; - not used at the moment 

    // Check if this charged and a primary 
    if (particle) {
      // Only charged tracks 
      if (!particle->Charge() != 0) continue;
    }
    Bool_t   isPrimary = stack->IsPhysicalPrimary(iTr) && iTr < fNPrimary;
    Double_t eloss     = hit->Edep();
    Double_t length    = hit->Length();
    // Double_t eta  = particle->Eta();
    Double_t x         = hit->X() - fIp[0];
    Double_t y         = hit->Y() - fIp[1];
    Double_t z         = hit->Z() - fIp[2];
    Double_t r         = TMath::Sqrt(x*x + y*y);
    Double_t phi       = TMath::ATan2(y, x);
    Double_t theta     = TMath::ATan2(r, z);
    Double_t eta       = -TMath::Log(TMath::Tan(theta/2));
    UShort_t det       = hit->Detector();
    Char_t   rng       = hit->Ring();
    UShort_t sec       = hit->Sector();
    UShort_t str       = hit->Strip();

    Bool_t isDistinct = (iTr != oldTr) || (det != oldD) || (rng != oldR);
    if (isDistinct) {
      RingHistos* h = static_cast<RingHistos*>(GetRingHistos(det, rng));
      h->fKind->Fill(eta, isPrimary ? .5 : 1.5);
    }
    oldTr = iTr;
    oldD  = det;
    oldR  = rng;

    AliFMDFloatMap& sum     =  (isPrimary ? fSumPrimary : fSumSecondary);
    // Float_t         old     =  sum(det, rng, sec, str);
    sum(det, rng, sec, str) += eloss;
#if 0
    Info("", "FMD%d%c[%02d,%03d]-%s: Adding %10.7f to %10.7f -> %10.7f",
	 det, rng, sec, str, (isPrimary ? "1st" : "2nd"), 
	 eloss, old, sum(det, rng, sec, str));
#endif
    if (!fTuple) continue;

    // Leaves primary:eta:phi:edep:length:detector:ring:sector:strip:ipz
    cache[0] = isPrimary ? 1 : 0;
    cache[1] = eta;
    cache[2] = phi;
    cache[3] = eloss;
    cache[4] = length;
    cache[5] = hit->Detector();
    cache[6] = Int_t(hit->Ring());
    cache[7] = hit->Sector();
    cache[8] = hit->Strip();
    cache[9] = fIp[2];
    fTuple->Fill(cache);
  }
  return true;
}

//____________________________________________________________________
Bool_t AliFMDMCHitEnergyFitter::PostEvent(const AliESDEvent& input)
{
  DGUARD(fDebug, 5, "Convert MC hit energy loss (%p)", &input);

  AliESDFMD*  esdFMD    = input.GetFMDData();
  for (UShort_t d = 1; d <= 3; d++) { 
    UShort_t nQ = (d == 1 ? 1 : 2);
    for (UShort_t q = 0; q < nQ; q++) { 
      UShort_t    nS = (q == 0 ?  20 :  40);
      UShort_t    nT = (q == 0 ? 512 : 256);
      Char_t      r  = (q == 0 ? 'I' : 'O');
      RingHistos* h  = 
	static_cast<AliFMDMCHitEnergyFitter::RingHistos*>(GetRingHistos(d,r));
      if (!h) continue;

      for (UShort_t s = 0; s < nS; s++) { 
	for (UShort_t t = 0; t < nT; t++) {
	  Float_t  totPrim = fSumPrimary(d, r, s, t);
	  Float_t  totSec  = fSumSecondary(d, r, s, t);
	  Float_t  totAll  = totPrim + totSec;
	  Double_t esdEta  = esdFMD->Eta(d, r, s, t);
	  if (totAll  > 0) h->FillMC(0,  esdEta, totAll);
	  if (totPrim > 0) h->FillMC(1,  esdEta, totPrim);
	  if (totSec  > 0) h->FillMC(2,  esdEta, totSec);
	}
      }
    }
  }
  return true;
}

//____________________________________________________________________
AliFMDEnergyFitter::RingHistos*
AliFMDMCHitEnergyFitter::CreateRingHistos(UShort_t d, Char_t r) const
{ 
  // DGUARD(fDebug, 1, "Create Ring cache of MC hit energy loss (FMD%d%c)",d,r);
  return new AliFMDMCHitEnergyFitter::RingHistos(d,r);
}

//====================================================================
AliFMDMCHitEnergyFitter::RingHistos::RingHistos()
  : AliFMDEnergyFitter::RingHistos(), 
    fPrimary(0), 
    fSecondary(0),
    fKind(0)
{}

//____________________________________________________________________
AliFMDMCHitEnergyFitter::RingHistos::RingHistos(UShort_t d, Char_t r)
  : AliFMDEnergyFitter::RingHistos(d,r), 
    fPrimary(0),
    fSecondary(0),
    fKind(0)
{}
//____________________________________________________________________
TArrayD
AliFMDMCHitEnergyFitter::RingHistos::MakeIncreasingAxis(Int_t, 
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
AliFMDMCHitEnergyFitter::RingHistos::SetupForData(const TAxis& eAxis, 
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
  fSecondary = Make("secondary","#sum#Delta_{true} of secondaries",
		    eAxis, maxDE, nDEbins, useIncrBin);
  fHist->SetXTitle("#Delta_{true}");
  fSecondary->SetXTitle("#Delta_{true}");
  fPrimary->SetXTitle("#Delta_{true}");

  fKind       = 0;
  if (eAxis.GetXbins()->GetArray()) 
    fKind = new TH2I("kind", "Particle types", 
		     eAxis.GetNbins(), eAxis.GetXbins()->GetArray(), 
		     2, 0, 2);
  else 
    fKind = new TH2I("kind", "Particle types", 
		     eAxis.GetNbins(), eAxis.GetXmin(), eAxis.GetXmax(),
		     2, 0, 2);
  fKind->SetXTitle("#eta");
  fKind->GetYaxis()->SetBinLabel(1, "Primary");
  fKind->GetYaxis()->SetBinLabel(2, "Secondary");
  fKind->SetDirectory(0);
  fKind->Sumw2();

  fList->Add(fHist);
  fList->Add(fPrimary);
  fList->Add(fSecondary);
  fList->Add(fKind);
}

//____________________________________________________________________
void
AliFMDMCHitEnergyFitter::RingHistos::FillMC(UShort_t flag, Double_t eta, 
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
TObjArray*
AliFMDMCHitEnergyFitter::RingHistos::Fit(TList*           dir, 
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

  out->SetName("compartive");
  out->SetOwner();
  l->Add(out);
  return out;
}
