#include "AliForwardMultDists.h" 
#include <AliForwardUtil.h> 
#include <AliAODForwardMult.h>
#include <AliAODCentralMult.h>
#include <AliAODEvent.h>
#include <AliLog.h>
#include <TH1.h>
#include <TH2.h>
#include <TVector2.h>
#include <THStack.h>
#include <TParameter.h>
#include <TError.h>
#include <iostream>
#include <iomanip> 

//====================================================================
namespace {
  /**
   * Marker styles 
   */
  enum { 
    kSolid        = 0x000, 
    kHollow       = 0x001, 
    kCircle       = 0x002,
    kSquare       = 0x004, 
    kUpTriangle   = 0x006, 
    kDownTriangle = 0x008, 
    kDiamond      = 0x00a,
    kCross        = 0x00c,
    kStar         = 0x00e
  };
  /** 
   * Get a ROOT marker style from bit options 
   * 
   * @param bits Bit mask of options 
   * 
   * @return ROOT marker style
   */
  Int_t GetMarkerStyle(UShort_t bits)
  {
    Int_t  base   = bits & (0xFE);
    Bool_t hollow = bits & kHollow;
    switch (base) { 
    case kCircle:       return (hollow ? 24 : 20);
    case kSquare:       return (hollow ? 25 : 21);
    case kUpTriangle:   return (hollow ? 26 : 22);
    case kDownTriangle: return (hollow ? 32 : 23);
    case kDiamond:      return (hollow ? 27 : 33); 
    case kCross:        return (hollow ? 28 : 34); 
    case kStar:         return (hollow ? 30 : 29); 
    }
    return 1;
  }
#if 0
  /** 
   * Get the marker option bits from a ROOT style 
   * 
   * @param style ROOT style 
   * 
   * @return Pattern of marker options
   */
  UShort_t GetMarkerBits(Int_t style) 
  { 
    UShort_t bits = 0;
    switch (style) { 
    case 24: case 25: case 26: case 27: case 28: case 30: case 32: 
      bits |= kHollow; break;
    }
    switch (style) { 
    case 20: case 24: bits |= kCircle;       break;
    case 21: case 25: bits |= kSquare;       break;
    case 22: case 26: bits |= kUpTriangle;   break;
    case 23: case 32: bits |= kDownTriangle; break;
    case 27: case 33: bits |= kDiamond;      break;
    case 28: case 34: bits |= kCross;        break;
    case 29: case 30: bits |= kStar;         break;
    }
    return bits;
  }
#endif 
  static Int_t GetIndexMarker(Int_t idx)
  {
    const UShort_t nMarkers = 7;
    UShort_t markers[] = { 
      kCircle, 
      kSquare, 
      kUpTriangle, 
      kDownTriangle, 
      kDiamond,
      kCross,
      kStar 
    };

    return markers[idx % nMarkers];
  }
#if 0
  /** 
   * Flip the 'hollow-ness' of a marker 
   * 
   * @param style ROOT Style 
   * 
   * @return ROOT style
   */
  Int_t FlipHollowStyle(Int_t style) 
  {
    UShort_t bits = GetMarkerBits(style);
    Int_t    ret  = GetMarkerStyle(bits ^ kHollow);
    return ret;
  }
#endif
}

//====================================================================
AliForwardMultDists::BinSpec::BinSpec(Double_t etaMin, 
				      Double_t etaMax, 
				      Double_t nchLow)
  : fEtaMin(etaMin), fEtaMax(etaMax), fLow(nchLow), fN(0), fD(0), fAxis(1,0,1)
{
}
//____________________________________________________________________
void
AliForwardMultDists::BinSpec::Push(UShort_t n, Double_t d)
{
  Int_t l = fN.GetSize();
  fN.Set(l+1); // Not terribly efficient, but that's ok because we do
  fD.Set(l+1); // this infrequently 
  fN[l] = n;
  fD[l] = d;
}
//____________________________________________________________________
const TAxis&
AliForwardMultDists::BinSpec::Axis() const
{
  if (fAxis.GetNbins() > 1) return fAxis;
  if (fN.GetSize() <= 0) return fAxis;
  if (fN.GetSize() == 1) { 
    // Exactly one spec, 
    fAxis.Set(fN[0], fLow, fN[0]*fD[0]+fLow);
  }
  else { 
    Int_t n = Int_t(fN.GetSum());
    Int_t l = fN.GetSize();
    Int_t i = 0;
    TArrayD b(n+1);
    b[0]    = fLow;
    for (Int_t j = 0; j < l; j++) { 
      for (Int_t k = 0; k < fN[j]; k++) {
	b[i+1] = b[i] + fD[j];
	i++;
      }
    }
    if (i != n) { 
      ::Warning("Axis", "Assigned bins, and number of bins not equal");
      n = TMath::Min(i, n);
    }
    fAxis.Set(n, b.GetArray());
  }
  return fAxis;
}

//====================================================================
AliForwardMultDists::AliForwardMultDists()
  : AliBaseAODTask(), 
    fBins(), 
    fSymmetric(0),
    fNegative(0), 
    fPositive(0),
    fMCVertex(0),
    fDiag(0),
    fForwardCache(0),
    fCentralCache(0),
    fMCCache(0),
    fUsePhiAcc(true),
    fIsSelected(false)
{}

//____________________________________________________________________
AliForwardMultDists::AliForwardMultDists(const char* name)
  : AliBaseAODTask(name,"AliForwardMultDists"),
    fBins(), 
    fSymmetric(0),
    fNegative(0), 
    fPositive(0),
    fMCVertex(0),
    fDiag(0),
    fForwardCache(0),
    fCentralCache(0),
    fMCCache(0),
    fUsePhiAcc(true),
    fIsSelected(false)
{
  /* 
   * User constructor 
   * 
   * @param name Name of the task 
   */
}


//____________________________________________________________________
Bool_t
AliForwardMultDists::Book()
{
  /* 
   * Create output objects - called at start of job in slave 
   * 
   */
  return true;
}
//____________________________________________________________________
Bool_t
AliForwardMultDists::PreData()
{
  /* 
   * Set-up internal structures on first seen event 
   * 
   */
  AliAODEvent*       aod     = AliForwardUtil::GetAODEvent(this);
  AliAODForwardMult* forward = GetForward(*aod);
  const TH2D&        hist    = forward->GetHistogram();
  Bool_t             useMC   = GetPrimary(*aod) != 0;

  fSums->Add(AliForwardUtil::MakeParameter("minIpZ", fIPzAxis.GetXmin()));
  fSums->Add(AliForwardUtil::MakeParameter("maxIpZ", fIPzAxis.GetXmax()));

  const TAxis* xaxis = hist.GetXaxis();
  if (!xaxis->GetXbins() || xaxis->GetXbins()->GetSize() <= 0) {
    fForwardCache = new TH1D("forwardCache", "Projection of Forward data", 
			     xaxis->GetNbins(), xaxis->GetXmin(), 
			     xaxis->GetXmax());
    fCentralCache = new TH1D("centralCache", "Projection of Central data", 
			     xaxis->GetNbins(), xaxis->GetXmin(), 
			     xaxis->GetXmax());
    if (useMC)
      fMCCache = new TH1D("mcCache", "Projection of Mc data", 
			  xaxis->GetNbins(), xaxis->GetXmin(), 
			  xaxis->GetXmax());
  }
  else { 
    fForwardCache = new TH1D("forwardCache", "Projection of Forward data", 
			     xaxis->GetNbins(),xaxis->GetXbins()->GetArray());
    fCentralCache = new TH1D("centralCache", "Projection of Central data", 
			     xaxis->GetNbins(),xaxis->GetXbins()->GetArray());
    if (useMC) 
      fMCCache = new TH1D("mcCache", "Projection of Mc data", 
			  xaxis->GetNbins(),xaxis->GetXbins()->GetArray());
  }
  fForwardCache->SetDirectory(0);
  fForwardCache->GetXaxis()->SetTitle("#eta");
  fForwardCache->GetYaxis()->SetTitle("#sum#frac{d^{2}N}{d#etadphi}");
  fForwardCache->Sumw2();
  fCentralCache->SetDirectory(0);
  fCentralCache->GetXaxis()->SetTitle("#eta");
  fCentralCache->GetYaxis()->SetTitle("#sum#frac{d^{2}N}{d#etadphi}");
  fCentralCache->Sumw2();

  if (useMC) {
    fMCCache->SetDirectory(0);
    fMCCache->GetXaxis()->SetTitle("#eta");
    fMCCache->GetYaxis()->SetTitle("#sum#frac{d^{2}N}{d#etadphi}");
    fMCCache->Sumw2();

    fMCVertex = static_cast<TH1*>(fVertex->Clone("mcVertex"));
    fMCVertex->SetTitle("Vertex distribution from MC");
    fMCVertex->SetDirectory(0);
    fMCVertex->SetFillColor(kBlue+2);
    fSums->Add(fMCVertex);
  }

  UShort_t xBase = kTrigger-1;
  UShort_t yBase = kAnalysis-1;
  fDiag = new TH2I("diagnostics", "Event selection", 
		   kTriggerVertex-xBase, kTrigger-.5,  kTriggerVertex+.5, 
		   kTriggerVertex-yBase, kAnalysis-.5, kTriggerVertex+.5);
  fDiag->SetDirectory(0);
  fDiag->GetXaxis()->SetTitle("Selection");
  fDiag->GetXaxis()->SetBinLabel(kTrigger      -xBase, "Trigger");
  fDiag->GetXaxis()->SetBinLabel(kVertex       -xBase, "Vertex");
  fDiag->GetXaxis()->SetBinLabel(kTriggerVertex-xBase, "Trigger+Vertex");
  fDiag->GetYaxis()->SetTitle("Type/MC selection");
  fDiag->GetYaxis()->SetBinLabel(kAnalysis     -yBase, "Analysis");
  fDiag->GetYaxis()->SetBinLabel(kMC           -yBase, "MC");
  fDiag->GetYaxis()->SetBinLabel(kTrigger      -yBase, "MC Trigger");
  fDiag->GetYaxis()->SetBinLabel(kVertex       -yBase, "MC Vertex");
  fDiag->GetYaxis()->SetBinLabel(kTriggerVertex-yBase, "MC Trigger+Vertex");
  fDiag->SetMarkerSize(3);
  fDiag->SetMarkerColor(kWhite);
  fSums->Add(fDiag);

  TIter   next(&fBins);
  EtaBin* bin = 0;
  while ((bin = static_cast<EtaBin*>(next()))) {
    bin->SetupForData(fSums, hist, useMC);
  }
  return true;
}
//____________________________________________________________________
Bool_t
AliForwardMultDists::Event(AliAODEvent& aod)
{
  /* 
   * Analyse a single event 
   * 
   * @param aod Input event
   */
  AliAODForwardMult* forward = GetForward(aod, false, true);
  if (!forward) return false;
    
  AliAODCentralMult* central = GetCentral(aod, false, true);
  if (!central) return false;
    
  TH2* primary = GetPrimary(aod);

  const TH2& forwardData = forward->GetHistogram();
  const TH2& centralData = central->GetHistogram();

  Double_t vz    = forward->GetIpZ();
  Bool_t   trg   = forward->IsTriggerBits(fTriggerMask);
  Bool_t   vtx   = (vz <= fIPzAxis.GetXmax() && vz >= fIPzAxis.GetXmin());
  Bool_t   ok    = true;  // Should bins process this event?
  Bool_t   mcTrg = false;
  Bool_t   mcVtx = false;
  // If we have MC inpunt
  if (primary) {
    // For MC, we need to check if we should process the event for
    // MC information or not.
    if ((fTriggerMask & (AliAODForwardMult::kV0AND|AliAODForwardMult::kNSD))
	&& !forward->IsTriggerBits(AliAODForwardMult::kMCNSD)) 
      // Bail out if this event is not MC NSD event
      ok = false;
    else 
      mcTrg = true;
    Double_t mcVz = primary->GetBinContent(0,0);
    fMCVertex->Fill(mcVz);
    if (mcVz > fIPzAxis.GetXmax() || mcVz < fIPzAxis.GetXmin()) 
      // Bail out if this event was not in the valid range 
      ok = false;
    else 
      mcVtx = true;
  }
  // Fill diagnostics 
  if (trg) { 
    fDiag->Fill(kTrigger, kAnalysis);
    if (mcTrg)          fDiag->Fill(kTrigger,       kTrigger);
    if (mcVtx)          fDiag->Fill(kTrigger,       kVertex);
    if (mcTrg && mcVtx) fDiag->Fill(kTrigger,       kTriggerVertex);
  }
  if (vtx) { 
    fDiag->Fill(kVertex, kAnalysis);
    if (mcTrg)          fDiag->Fill(kVertex,        kTrigger);
    if (mcVtx)          fDiag->Fill(kVertex,        kVertex);
    if (mcTrg && mcVtx) fDiag->Fill(kVertex,        kTriggerVertex);
  }
  if (trg && vtx) {
    fDiag->Fill(kTriggerVertex, kAnalysis);
    if (mcTrg)          fDiag->Fill(kTriggerVertex, kTrigger);
    if (mcVtx)          fDiag->Fill(kTriggerVertex, kVertex);
    if (mcTrg && mcVtx) fDiag->Fill(kTriggerVertex, kTriggerVertex);
  }
  if (primary) {
    if (mcTrg)          fDiag->Fill(kTrigger,       kMC);
    if (mcVtx)          fDiag->Fill(kVertex,        kMC);
    if (mcTrg && mcVtx) fDiag->Fill(kTriggerVertex, kMC);
  }
  if (!fIsSelected && !ok) {
    // Printf("Event is neither accepted by analysis nor by MC");
    return false; 
  }

  // forward->Print();
  // Info("", "Event vz=%f", forward->GetIpZ());
  
  ProjectX(forwardData, *fForwardCache, fUsePhiAcc);
  ProjectX(centralData, *fCentralCache, fUsePhiAcc);
  ProjectX(primary, fMCCache);
    
  TIter   next(&fBins);
  EtaBin* bin = 0;
  while ((bin = static_cast<EtaBin*>(next()))) {
    bin->Process(*fForwardCache, *fCentralCache,
		 forwardData,    centralData, 
		 fIsSelected,	 fMCCache);
  }

  return true;
}

//_____________________________________________________________________
Bool_t 
AliForwardMultDists::CheckEvent(const AliAODForwardMult& fwd)
{
  fIsSelected = AliBaseAODTask::CheckEvent(fwd);
  // We always return true, so that we can process MC truth;
  return true;
}

//____________________________________________________________________
Bool_t
AliForwardMultDists::Finalize()
{
  /*
   * Called at the end of the final processing of the job on the
   * full data set (merged data)
   */
  fResults->Add(fSums->FindObject("triggers")->Clone());
  fResults->Add(fSums->FindObject("status")->Clone());
  fResults->Add(fSums->FindObject("diagnostics")->Clone());

  THStack* sym    = new THStack("all", "All distributions");
  THStack* pos    = new THStack("all", "All distributions");
  THStack* neg    = new THStack("all", "All distributions");
  THStack* oth    = new THStack("all", "All distributions");
  THStack* sta    = 0;
  EtaBin*  bin    = 0;
  TIter    next(&fBins);
  while ((bin = static_cast<EtaBin*>(next()))) {
    bin->Terminate(fSums, fResults);

    sta = oth;
    if      (bin->IsSymmetric()) sta = sym;
    else if (bin->IsNegative())  sta = neg;
    else if (bin->IsPositive())  sta = pos;

    TH1*  hh[] = { bin->fSum, bin->fTruth, bin->fTruthAccepted, 0 };
    TH1** ph   = hh;

    Int_t idx     = sta->GetHists() ? sta->GetHists()->GetEntries() : 0;
    if (bin->fTruth) idx /= 2;
    
    Int_t mkrBits = GetIndexMarker(idx);
    while (*ph) { 
      Int_t factor = Int_t(TMath::Power(10, idx));
      TH1* h = static_cast<TH1*>((*ph)->Clone());
      h->SetDirectory(0);
      h->Scale(factor);
      h->SetTitle(Form("%s (#times%d)", h->GetTitle(), Int_t(factor)));
      h->SetMarkerStyle(GetMarkerStyle(mkrBits));
      mkrBits ^= kHollow;
      
      sta->Add(h, "p");
      ph++;
    }
  }
  TList* lsym = static_cast<TList*>(fResults->FindObject("symmetric"));
  TList* lneg = static_cast<TList*>(fResults->FindObject("negative"));
  TList* lpos = static_cast<TList*>(fResults->FindObject("positive"));
  TList* loth = static_cast<TList*>(fResults->FindObject("other"));
  
  if (lsym) lsym->Add(sym);
  if (lneg) lneg->Add(neg);
  if (lpos) lpos->Add(pos);
  if (loth) loth->Add(oth);

  return true;
}

//____________________________________________________________________
void
AliForwardMultDists::ProjectX(const TH2& input, TH1& cache, Bool_t usePhiAcc)
{
  /*
   * Project a 2D histogram into a 1D histogram taking care to use
   * either the @f$\phi@f$ acceptance stored in the overflow bins, or
   * the @f$\eta@f$ coverage stored in the underflow bins.
   * 
   * @param input      2D histogram to project 
   * @param cache      1D histogram to project into 
   * @param usePhiAcc  If true, use the @f$\phi@f$ acceptance stored in
   * the overflow bins, or if false the @f$\eta@f$ coverage stored in
   * the underflow bins.
   */
  cache.Reset();
    
  Int_t nPhi = input.GetNbinsY();
  Int_t nEta = input.GetNbinsX();

  for (Int_t iEta = 1; iEta <= nEta; iEta++) { 
    Double_t phiAcc = input.GetBinContent(iEta, nPhi+1);
    Double_t etaCov = input.GetBinContent(iEta, 0);
    Double_t factor = usePhiAcc ? phiAcc : etaCov;
      
    if (factor < 1e-3) continue; 
    Double_t sum    = 0;
    Double_t e2sum  = 0;
    for (Int_t iPhi = 1; iPhi <= nPhi; iPhi++) {
      Double_t c = input.GetBinContent(iEta, iPhi);
      Double_t e = input.GetBinError(iEta, iPhi);
	
      sum   += c;
      e2sum += e * e;
    }
    cache.SetBinContent(iEta, factor * sum);
    cache.SetBinError(iEta, factor * TMath::Sqrt(e2sum));
  }
}
//____________________________________________________________________
void
AliForwardMultDists::ProjectX(const TH2* input, TH1* cache)
{
  /*
   * Project on @f$\eta@f$ axis.  If any of the pointers passed is
   * zero, do nothing.
   * 
   * @param input 
   * @param cache 
   */
  if (!input || !cache) return;
  cache->Reset();
    
  Int_t nPhi = input->GetNbinsY();
  Int_t nEta = input->GetNbinsX();

  for (Int_t iEta = 1; iEta <= nEta; iEta++) { 
    Double_t sum    = 0;
    Double_t e2sum  = 0;
    for (Int_t iPhi = 1; iPhi <= nPhi; iPhi++) {
      Double_t c = input->GetBinContent(iEta, iPhi);
      Double_t e = input->GetBinError(iEta, iPhi);
	
      sum   += c;
      e2sum += e * e;
    }
    cache->SetBinContent(iEta, sum);
    cache->SetBinError(iEta, TMath::Sqrt(e2sum));
  }
}
//____________________________________________________________________
void
AliForwardMultDists::AddBin(const BinSpec& spec)
{
  AddBin(spec.fEtaMin, spec.fEtaMax, spec.Axis());
}

//____________________________________________________________________
void
AliForwardMultDists::AddBin(Double_t etaLow, Double_t etaMax, 
			    const TAxis& mAxis) 
{
  /*
   * Add an @f$\eta@f$ bin
   * 
   * @param etaLow Low cut on @f$\eta@f$
   * @param etaMax High cut on @f$\eta@f$
   */
  // Symmetric bin
  if (etaMax >= kInvalidEta) { 
    etaLow = -TMath::Abs(etaLow);
    etaMax = +TMath::Abs(etaLow);
  }
  EtaBin* bin = new EtaBin(etaLow, etaMax, mAxis);
  // AliInfoF("Adding bin %f, %f: %s", etaLow, etaMax, bin->GetName());
  fBins.Add(bin);
}
//____________________________________________________________________
void
AliForwardMultDists::AddBin(Double_t etaLow, Double_t etaMax, 
			    UShort_t nMax, UShort_t nDiv)
{
  /* 
   * Add an @f$\eta@f$ bin
   * 
   * @param etaLow Low cut on @f$\eta@f$
   * @param etaMax High cut on @f$\eta@f$
   */
  // Symmetric bin
  if (etaMax >= kInvalidEta) { 
    etaLow = -TMath::Abs(etaLow);
    etaMax = +TMath::Abs(etaLow);
  }
  TAxis mAxis((nMax+1)/nDiv, -.5, nMax+.5);
  EtaBin* bin = new EtaBin(etaLow, etaMax, mAxis);
  // AliInfoF("Adding bin %f, %f: %s", etaLow, etaMax, bin->GetName());
  fBins.Add(bin);
}
#define PFV(N,VALUE)					\
  do {							\
    AliForwardUtil::PrintName(N);			\
    std::cout << (VALUE) << std::endl; } while(false)
#define PFB(N,FLAG)				\
  do {									\
    AliForwardUtil::PrintName(N);					\
    std::cout << std::boolalpha << (FLAG) << std::noboolalpha << std::endl; \
  } while(false)

//____________________________________________________________________
void
AliForwardMultDists::Print(Option_t* option) const 
{
  /* 
   * Print this task 
   * 
   * @param option Not used
   */
  AliBaseAODTask::Print(option);
  gROOT->IncreaseDirLevel();
  PFB("Use phi acceptance", fUsePhiAcc);
  // AliForwardUtil::PrintName("Bins");

  Int_t first = true;
  TIter   next(&fBins);
  EtaBin* bin = 0;
  while ((bin = static_cast<EtaBin*>(next()))) {
    PFV(first ? "Bins" : "", bin->GetName());
    first = false;
  }
  gROOT->DecreaseDirLevel();
}

//====================================================================
AliForwardMultDists::EtaBin::EtaBin() 
  : fName(""),
    fMAxis(1,0,0),
    fTAxis(1,0,0),
    fMinEta(0),
    fMaxEta(0),
    fMinBin(0), 
    fMaxBin(0), 
    fSum(0),
    fCorr(0),
    fResponse(0), 
    fTruth(0),
    fTruthAccepted(0),
    fCoverage(0)
{
  /* 
   * I/O constructor
   */
}

//____________________________________________________________________
AliForwardMultDists::EtaBin::EtaBin(Double_t minEta, Double_t maxEta,
				    const TAxis& mAxis) 
  : fName(""),
    fMAxis(1,0,0),
    fTAxis(1,0,0),
    fMinEta(minEta), 
    fMaxEta(maxEta),
    fMinBin(0), 
    fMaxBin(0), 
    fSum(0),
    fCorr(0),
    fResponse(0), 
    fTruth(0),
    fTruthAccepted(0),
    fCoverage(0)
{
  /*
   * User constructor 
   * 
   * @param minEta Least @f$\eta@f$ to consider 
   * @param maxEta Largest @f$\eta@f$ to consider 
   */
  fName = TString::Format("%+05.2f_%+05.2f", fMinEta, fMaxEta);
  fName.ReplaceAll("-", "m");
  fName.ReplaceAll("+", "p");
  fName.ReplaceAll(".", "d");

  // Copy to other our object
  mAxis.Copy(fMAxis);
  
  if (mAxis.GetXbins() && mAxis.GetXbins()->GetArray()) {
    const TArrayD& mA  = *(mAxis.GetXbins());
    TArrayD        tA(mA.GetSize());
    Int_t          j   = 0;
    Double_t       min = mA[0];
    tA[0]              = min;
    for (Int_t i = 1; i < mA.GetSize(); i++) { 
      Double_t d = mA[i] - min;
      if (d < 1) 
	// Not full integer bin
	continue;  
      tA[j+1] = tA[j] + Int_t(d);
      min     = tA[j+1];
      j++;
    }
    fTAxis.Set(j, tA.GetArray());
  }
  else {
    // Rounded down maximum and minimum
    Int_t max = Int_t(mAxis.GetXmax()); 
    Int_t min = Int_t(mAxis.GetXmin());
    fTAxis.Set((max-min)+1, min-.5, max+.5);
  }
}
//____________________________________________________________________
AliForwardMultDists::EtaBin::EtaBin(const EtaBin& o) 
  : TObject(o),
    fName(o.fName),
    fMAxis(o.fMAxis), 
    fTAxis(o.fTAxis),
    fMinEta(o.fMinEta),
    fMaxEta(o.fMaxEta),
    fMinBin(o.fMinBin), 
    fMaxBin(o.fMaxBin), 
    fSum(o.fSum),
    fCorr(o.fCorr),
    fResponse(o.fResponse), 
    fTruth(o.fTruth),
    fTruthAccepted(o.fTruthAccepted),
    fCoverage(o.fCoverage)
{}
//____________________________________________________________________
AliForwardMultDists::EtaBin&
AliForwardMultDists::EtaBin::operator=(const EtaBin& o) 
{
  if (&o == this) return *this;
  
  fName		 = o.fName;
  fMAxis         = o.fMAxis;
  fTAxis         = o.fTAxis;
  fMinEta	 = o.fMinEta;
  fMaxEta	 = o.fMaxEta;
  fMinBin	 = o.fMinBin; 
  fMaxBin	 = o.fMaxBin; 
  fSum		 = o.fSum;
  fCorr		 = o.fCorr;
  fResponse	 = o.fResponse; 
  fTruth	 = o.fTruth;
  fTruthAccepted = o.fTruthAccepted;
  fCoverage      = o.fCoverage;

  return *this;
}
//____________________________________________________________________
Bool_t
AliForwardMultDists::EtaBin::IsSymmetric() const
{
  return TMath::Abs(fMaxEta + fMinEta) < 1e-6;
}
//____________________________________________________________________
Bool_t
AliForwardMultDists::EtaBin::IsNegative() const
{
  return TMath::Abs(fMaxEta) < 1e-6 && fMinEta < 0;
}
//____________________________________________________________________
Bool_t
AliForwardMultDists::EtaBin::IsPositive() const
{
  return TMath::Abs(fMinEta) < 1e-6 && fMaxEta > 0;
}
//____________________________________________________________________
const char*
AliForwardMultDists::EtaBin::ParentName() const
{
  if      (IsSymmetric()) return "symmetric";
  else if (IsNegative())  return "negative";
  else if (IsPositive())  return "positive";
  return "other";
}
//____________________________________________________________________
TList*
AliForwardMultDists::EtaBin::FindParent(TList* l, Bool_t create) const
{
  const char* parent = ParentName();
  TObject*    op     = l->FindObject(parent);

  if (op) return static_cast<TList*>(op);
  if (!create) return 0;

  // Info("FindParent", "Parent %s not found in %s, creating and adding",
  //      parent, l->GetName());
  TList* p = new TList;
  p->SetName(parent);
  p->SetOwner();
  l->Add(p);

  TList* lowEdges = new TList;
  lowEdges->SetName("lowEdges");
  lowEdges->SetOwner();
  p->Add(lowEdges);

  TList* highEdges = new TList;
  highEdges->SetName("highEdges");
  highEdges->SetOwner();
  p->Add(highEdges);
  
  return p;
}

//____________________________________________________________________
TH1*
AliForwardMultDists::EtaBin::CreateH1(const char*  name, 
				      const char*  title, 
				      const TAxis& xAxis)
{
  TH1* ret = 0;
  
  if (xAxis.GetXbins() && xAxis.GetXbins()->GetArray()) 
    ret = new TH1D(name,title,xAxis.GetNbins(), xAxis.GetXbins()->GetArray());
  else 
    ret = new TH1D(name,title,xAxis.GetNbins(),xAxis.GetXmin(),xAxis.GetXmax());
  return ret;
}
//____________________________________________________________________
TH2*
AliForwardMultDists::EtaBin::CreateH2(const char*  name, 
				      const char*  title, 
				      const TAxis& xAxis,
				      const TAxis& yAxis)
{
  TH2* ret = 0;
  
  if (xAxis.GetXbins() && xAxis.GetXbins()->GetArray()) {
    if (yAxis.GetXbins() && yAxis.GetXbins()->GetArray()) 
      // Variable variable
      ret = new TH2D(name,title, 
		     xAxis.GetNbins(), xAxis.GetXbins()->GetArray(),
		     yAxis.GetNbins(), yAxis.GetXbins()->GetArray());
    else 
      // variable fixed
      ret = new TH2D(name,title, 
		     xAxis.GetNbins(), xAxis.GetXbins()->GetArray(),
		     yAxis.GetNbins(),yAxis.GetXmin(),yAxis.GetXmax());
  }
  else {
    if (yAxis.GetXbins() && yAxis.GetXbins()->GetArray()) 
      // fixed variable 
      ret = new TH2D(name,title, 
		     xAxis.GetNbins(), xAxis.GetXmin(), xAxis.GetXmax(),
		     yAxis.GetNbins(), yAxis.GetXbins()->GetArray());
    else 
      // fixed fixed
      ret = new TH2D(name,title, 
		     xAxis.GetNbins(), xAxis.GetXmin(), xAxis.GetXmax(),
		     yAxis.GetNbins(), yAxis.GetXmin(), yAxis.GetXmax());
  }
  return ret;
}

//____________________________________________________________________
void
AliForwardMultDists::EtaBin::SetupForData(TList* list, const TH2& hist, 
					  Bool_t useMC)
{
  /*
   * Set-up internal structures on first event. 
   * 
   * @param list  List to add information to
   * @param hist  Template histogram 
   * @param max   Maximum number of particles 
   * @param useMC Whether to set-up for MC input 
   */
  TList* p = FindParent(list, true);
  TList* l = new TList;
  l->SetName(GetName());
  l->SetOwner();
  p->Add(l);
  TList* le = static_cast<TList*>(p->FindObject("lowEdges"));
  TList* he = static_cast<TList*>(p->FindObject("highEdges"));
  if (!le || !he) { 
    AliError("Failed to get bin edge lists from parent");
    return;
  }
  else {
    Int_t n = le->GetEntries();
    TParameter<double>* lp = 
      new TParameter<double>(Form("minEta%02d", n), fMinEta);
    TParameter<double>* hp = 
      new TParameter<double>(Form("maxEta%02d", n), fMaxEta);
    lp->SetMergeMode('f');
    hp->SetMergeMode('f');
    le->Add(lp);
    he->Add(hp);
  }

  const TAxis * xax = hist.GetXaxis();
  fMinBin = xax->FindFixBin(fMinEta);
  fMaxBin = xax->FindFixBin(fMaxEta-.00001);

  TString t(Form("%+5.2f<#eta<%+5.2f", fMinEta, fMaxEta));
  fSum = CreateH1("rawDist",Form("Raw P(#it{N}_{ch}) in %s",t.Data()),fMAxis);
  fSum->SetDirectory(0);
  fSum->GetXaxis()->SetTitle("#it{N}_{ch}");
  fSum->GetYaxis()->SetTitle("Raw P(#it{N}_{ch})");
  fSum->SetFillColor(kRed+1);
  fSum->SetFillStyle(0);
  // fSum->SetOption("hist e");
  fSum->SetMarkerStyle(20);
  fSum->SetMarkerColor(kRed+1);
  fSum->SetLineColor(kBlack);
  fSum->Sumw2();
  l->Add(fSum);
  
  fCorr = CreateH2("corr",Form("C_{SPD,FMD} in %s", t.Data()),fTAxis,fTAxis);
  fCorr->SetDirectory(0);
  fCorr->GetXaxis()->SetTitle("Forward");
  fCorr->GetYaxis()->SetTitle("Central");
  fCorr->SetOption("colz");
  l->Add(fCorr);

  fCoverage = new TH1D("coverage", "Fraction of bins with coverage",
		       101, -.5, 100.5);
  fCoverage->SetDirectory(0);
  fCoverage->SetXTitle("Fraction of bins [%]");
  fCoverage->SetFillStyle(3001);
  fCoverage->SetFillColor(kGreen+1);
  l->Add(fCoverage);

  if (!useMC) return;
  fResponse = CreateH2("response", Form("Reponse matrix in %s", t.Data()),
		       fMAxis, fTAxis);
  fResponse->SetDirectory(0);
  fResponse->SetYTitle("MC");
  fResponse->SetXTitle("Signal");
  fResponse->SetOption("colz");
  l->Add(fResponse);
  
  fTruth = CreateH1("truth",Form("True P(#it{N}_{ch}) in %s",t.Data()),fTAxis);
  fTruth->SetXTitle(fSum->GetXaxis()->GetTitle());
  fTruth->SetYTitle(fSum->GetYaxis()->GetTitle());
  fTruth->SetLineColor(kBlack);
  fTruth->SetFillColor(kBlue+1);
  fTruth->SetFillStyle(0);
  fTruth->SetDirectory(0);
  /// fTruth->SetOption("");
  fTruth->SetMarkerColor(kBlue+1);
  fTruth->SetMarkerStyle(24);
  fTruth->Sumw2();
  l->Add(fTruth);

  fTruthAccepted = static_cast<TH1D*>(fTruth->Clone("truthAccepted"));
  fTruthAccepted->SetTitle(Form("True (accepted) P(#it{N}_{ch}) in %s", 
				t.Data()));
  fTruthAccepted->SetLineColor(kGray+2);
  fTruthAccepted->SetFillColor(kOrange+2);
  fTruthAccepted->SetDirectory(0);
  /// fTruth->SetOption("");
  fTruthAccepted->SetMarkerColor(kOrange+2);
  l->Add(fTruthAccepted);
}
//____________________________________________________________________
void
AliForwardMultDists::EtaBin::Process(const TH1& sumForward, 
				     const TH1& sumCentral,
				     const TH2& forward,   
				     const TH2& central,
				     Bool_t     accepted, 
				     const TH1* mc)
{
  /*
   * Process a single event 
   * 
   * @param sumForward  Projection of forward data
   * @param sumCentral  Projection of the central data
   * @param forward     The original forward data 
   * @param central     The original central data
   */
  if (!mc && !accepted) 
    // If we're not looking at MC data, and the event wasn't selected,
    // we bail out - this check is already done, but we do it again
    // for safety.
    return;

  Double_t sum        = 0;
  Double_t e2sum      = 0;
  Int_t    covered    = 0;
  Double_t fsum       = -1;
  Double_t csum       = -1;
  Double_t mcsum      = 0;
  Double_t mce2sum    = 0;
  for (Int_t iEta = fMinBin; iEta <= fMaxBin; iEta++) { 
    if (mc) {
      Double_t cM = mc->GetBinContent(iEta);
      Double_t eM = mc->GetBinError(iEta);
      mcsum   += cM;
      mce2sum += eM * eM;
    }
    if (!accepted) 
      // Event wasn't selected, but we still need to get the rest of
      // the MC data.
      continue;

    Double_t cF = sumForward.GetBinContent(iEta);
    Double_t eF = sumForward.GetBinError(iEta);
    Bool_t   bF = forward.GetBinContent(iEta, 0) > 0;
    Double_t cC = sumCentral.GetBinContent(iEta);
    Double_t eC = sumCentral.GetBinError(iEta);
    Bool_t   bC = central.GetBinContent(iEta, 0) > 0;
    Double_t c  = 0;
    Double_t e  = 0;

    // If we have an overlap - as given by the eta-coverage, 
    // calculate the mean 
    if (bF & bC) { 
      c = (cF + cC) / 2;
      e = TMath::Sqrt(eF*eF + eC*eC);
      covered++;
    }
    // Else, if we have coverage from forward, use that 
    else if (bF) { 
      c = cF;
      e = eF;
      covered++;
    }
    // Else, if we have covrage from central, use that 
    else if (bC) { 
      c = cC;
      e = eC;
      covered++;
    }
    // Else, we have incomplete coverage 

    if (bF) { 
      if (fsum < 0) fsum = 0;
      fsum += cF;
    }
    if (bC) { 
      if (csum < 0) csum = 0;
      csum += cC;
    }
	
    sum   += c;
    e2sum += e*e;
  }
      
  if (accepted) {
    // Only update the histograms if the event was triggered. 
    // Fill with weight 
    Double_t w = 1; // e2sum > 0 ? 1/TMath::Sqrt(e2sum) : 1
    fSum->Fill(sum, w);
    fCorr->Fill((fsum<0?0:fsum), (csum<0?0:csum));
    
    Int_t nTotal = fMaxBin - fMinBin + 1;
    fCoverage->Fill(100*float(covered) / nTotal);
  }

  if (mc) {
    Double_t w = 1; // mce2sum > 0 ? 1/TMath::Sqrt(mce2sum) : 1
    if (fTruth) {
      fTruth->Fill(mcsum, w);
    }
    if (accepted) {
      if (fResponse) 
	// Only fill response matrix for accepted events
	fResponse->Fill(sum, mcsum);
      if (fTruthAccepted) 
	fTruthAccepted->Fill(mcsum, w);
    }
  }
}
//____________________________________________________________________
void
AliForwardMultDists::EtaBin::Terminate(TList* in, TList* out)
{
  /*
   * Called at the end of the final processing of the job on the
   * full data set (merged data)
   * 
   * @param in    Input list
   * @param out   Output list 
   * @param maxN  Maximum number of @f$N_{ch}@f$ to consider
   */
  TList* ip = FindParent(in, false);
  if (!ip) { 
    AliErrorF("Parent folder %s not found in input", ParentName());
    return;
  }

  TList* i = dynamic_cast<TList*>(ip->FindObject(GetName()));
  if (!i) { 
    AliErrorF("List %s not found in input", GetName());
    return;
  }
      
  TList* op = FindParent(out, true);
  TList* o  = static_cast<TList*>(i->Clone());
  o->SetOwner();
  op->Add(o);

  fSum           = static_cast<TH1*>(o->FindObject("rawDist"));
  fTruth         = static_cast<TH1*>(o->FindObject("truth"));
  fTruthAccepted = static_cast<TH1*>(o->FindObject("truthAccepted"));

  TH1*  hists[] = { fSum, fTruth, fTruthAccepted, 0 };
  TH1** phist   = hists;
  while (*phist) { 
    TH1* h = *phist;
    if (h) { 
      Int_t    maxN = h->GetNbinsX();
      Double_t intg = h->Integral(1, maxN);
      h->Scale(1. / intg, "width");
    }
    phist++;
  }

  if (fTruth && fTruthAccepted) {
    TH1*  trgVtx  = static_cast<TH1*>(fTruthAccepted->Clone("triggerVertex"));
    TString tit(fTruth->GetTitle());
    tit.ReplaceAll("True P(#it{N}_{ch})", "C_{trigger,vertex}");
    trgVtx->SetTitle(tit);
    trgVtx->SetYTitle("P_{MC}(#it{N}_{ch})/P_{MC,accepted}(#it{N}_{ch})");
    trgVtx->Divide(fTruth);
    trgVtx->SetDirectory(0);
    o->Add(trgVtx);
  }
  
}
//====================================================================
// 
// EOF
//
