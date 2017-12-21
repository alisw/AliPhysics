//
// Class to do the energy correction of FMD ESD data
//
#include "AliFMDEnergyFitter.h"
#include "AliForwardUtil.h"
#include "AliLandauGausFitter.h"
#include <AliESDFMD.h>
#include <TAxis.h>
#include <TList.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TMath.h>
#include <AliLog.h>
#include <TClonesArray.h>
#include <TFitResult.h>
#include <THStack.h>
#include <TROOT.h>
#include <iostream>
#include <iomanip>

ClassImp(AliFMDEnergyFitter)
#if 0
; // This is for Emacs
#endif 
namespace {
  const char* fgkEDistFormat = "%s_etabin%03d";
}


//____________________________________________________________________
AliFMDEnergyFitter::AliFMDEnergyFitter()
  : TNamed(), 
    fRingHistos(),
    fLowCut(0.4),
    fNParticles(5),
    fMinEntries(10000),
    fFitRangeBinWidth(4),
    fDoFits(false),
    fDoMakeObject(false),
    fEtaAxis(),
    fCentralityAxis(),
    fMaxE(10),
    fNEbins(300), 
    fUseIncreasingBins(true),
    fMaxRelParError(AliFMDCorrELossFit::ELossFit::fgMaxRelError),
    fMaxChi2PerNDF(AliFMDCorrELossFit::ELossFit::fgMaxChi2nu), 
    fMinWeight(AliFMDCorrELossFit::ELossFit::fgLeastWeight),
    fDebug(0),
    fResidualMethod(kNoResiduals),
    fSkips(0),
    fRegularizationCut(3e6)
{
  // 
  // Default Constructor - do not use 
  //
  // DGUARD(fDebug, 1, "Default CTOR of AliFMDEnergyFitter");
  fRingHistos.SetOwner();
}

//____________________________________________________________________
AliFMDEnergyFitter::AliFMDEnergyFitter(const char* title)
  : TNamed("fmdEnergyFitter", title), 
    fRingHistos(), 
    fLowCut(0.4),
    fNParticles(5),
    fMinEntries(10000),
    fFitRangeBinWidth(4),
    fDoFits(false),
    fDoMakeObject(false),
    fEtaAxis(0,0,0),
    fCentralityAxis(0,0,0),
    fMaxE(10),
    fNEbins(300), 
    fUseIncreasingBins(true),
    fMaxRelParError(AliFMDCorrELossFit::ELossFit::fgMaxRelError),
    fMaxChi2PerNDF(AliFMDCorrELossFit::ELossFit::fgMaxChi2nu), 
    fMinWeight(AliFMDCorrELossFit::ELossFit::fgLeastWeight),
    fDebug(0),
    fResidualMethod(kNoResiduals),
    fSkips(0),
    fRegularizationCut(3e6)
{
  // 
  // Constructor 
  // 
  // Parameters:
  //    title Title of object  - not significant 
  //
  // DGUARD(fDebug, 1, "Named CTOR of AliFMDEnergyFitter: %s", title);
  fEtaAxis.SetName("etaAxis");
  fEtaAxis.SetTitle("#eta");
  fCentralityAxis.SetName("centralityAxis");
  fCentralityAxis.SetTitle("Centrality [%]");
}

//____________________________________________________________________
AliFMDEnergyFitter::~AliFMDEnergyFitter()
{
  // 
  // Destructor
  //
  DGUARD(fDebug, 1, "DTOR of AliFMDEnergyFitter");
  // fRingHistos.Delete();
}

//____________________________________________________________________
AliFMDEnergyFitter::RingHistos*
AliFMDEnergyFitter::CreateRingHistos(UShort_t d, Char_t r) const
{
  return new RingHistos(d,r);
}

//____________________________________________________________________
AliFMDEnergyFitter::RingHistos*
AliFMDEnergyFitter::GetRingHistos(UShort_t d, Char_t r) const
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
AliFMDEnergyFitter::Init()
{
  // Create the ring histograms.  
  // 
  // Should be called from task initialization. 
  DGUARD(1,fDebug, "Creating histogram caches for each ring");
  fRingHistos.Add(CreateRingHistos(1, 'I'));
  fRingHistos.Add(CreateRingHistos(2, 'I'));
  fRingHistos.Add(CreateRingHistos(2, 'O'));
  fRingHistos.Add(CreateRingHistos(3, 'I'));
  fRingHistos.Add(CreateRingHistos(3, 'O'));
  TIter    next(&fRingHistos);
  RingHistos* o = 0;
  while ((o = static_cast<RingHistos*>(next()))) {
    o->fDebug = fDebug;
  }
}

//____________________________________________________________________
void
AliFMDEnergyFitter::CreateOutputObjects(TList* dir)
{
  // 
  // Define the output histograms.  These are put in a sub list of the
  // passed list.   The histograms are merged before the parent task calls 
  // AliAnalysisTaskSE::Terminate.  Called on task initialization on slaves.  
  // 
  // Parameters:
  //    dir Directory to add to 
  //
  DGUARD(fDebug, 1, "Define output in AliFMDEnergyFitter");
  TList* d = new TList;
  d->SetName(GetName());
  d->SetOwner(true);
  dir->Add(d);

  // Store eta axis as a histogram, since that can be merged!
  TH1* hEta = 0;
  if (fEtaAxis.GetXbins()->GetArray()) 
    hEta = new TH1I(fEtaAxis.GetName(), fEtaAxis.GetTitle(), 
		    fEtaAxis.GetNbins(), fEtaAxis.GetXbins()->GetArray());
  else 
    hEta = new TH1I(fEtaAxis.GetName(), fEtaAxis.GetTitle(), 
		    fEtaAxis.GetNbins(),fEtaAxis.GetXmin(),fEtaAxis.GetXmax());
  hEta->SetXTitle("#eta");
  hEta->SetYTitle("Nothing");
  hEta->SetDirectory(0);

  d->Add(hEta);
  d->Add(AliForwardUtil::MakeParameter("lowCut",        fLowCut));
  d->Add(AliForwardUtil::MakeParameter("nParticles",    fNParticles));
  d->Add(AliForwardUtil::MakeParameter("minEntries",    fMinEntries));
  d->Add(AliForwardUtil::MakeParameter("subtractBins",  fFitRangeBinWidth));
  d->Add(AliForwardUtil::MakeParameter("doFits",        fDoFits));
  d->Add(AliForwardUtil::MakeParameter("doObject",      fDoMakeObject));
  d->Add(AliForwardUtil::MakeParameter("maxE",          fMaxE));
  d->Add(AliForwardUtil::MakeParameter("nEbins",        fNEbins));
  d->Add(AliForwardUtil::MakeParameter("increasingBins",fUseIncreasingBins));
  d->Add(AliForwardUtil::MakeParameter("maxRelPerError",fMaxRelParError));
  d->Add(AliForwardUtil::MakeParameter("maxChi2PerNDF", fMaxChi2PerNDF));
  d->Add(AliForwardUtil::MakeParameter("minWeight",     fMinWeight));
  d->Add(AliForwardUtil::MakeParameter("regCut",        fRegularizationCut));
  d->Add(AliForwardUtil::MakeParameter("deltaShift", 
				       AliLandauGaus::EnableSigmaShift()));

  if (fRingHistos.GetEntries() <= 0) { 
    AliFatal("No ring histograms where defined - giving up!");
    return;
  }
  TIter    next(&fRingHistos);
  RingHistos* o = 0;
  while ((o = static_cast<RingHistos*>(next()))) {
    o->fDebug = fDebug;
    o->CreateOutputObjects(d);
  }
}

//____________________________________________________________________
void
AliFMDEnergyFitter::SetupForData(const TAxis& eAxis, UShort_t sys)
{
  // 
  // Initialise the task - called at first event 
  // 
  // Parameters:
  //    etaAxis The eta axis to use.  Note, that if the eta axis
  // has already been set (using SetEtaAxis), then this parameter will be 
  // ignored
  //
  DGUARD(fDebug, 1, "Initialize of AliFMDEnergyFitter");
  switch (sys) {
  case 1: fNParticles = 3; break; // pp 
  case 2: // Fall through
  case 3: // Fall through
  case 4: // Fall through
  case 5: fNParticles = 5; break; // pA, Ap, PbPb, XeXe
  default: break; // Do nothing. 
  }
  
  if (fEtaAxis.GetNbins() == 0 || 
      TMath::Abs(fEtaAxis.GetXmax() - fEtaAxis.GetXmin()) < 1e-7) 
    SetEtaAxis(eAxis);
  if (fCentralityAxis.GetNbins() == 0) {
    UShort_t n = 12;
    Double_t bins[] = {  0.,  5., 10., 15., 20., 30., 
			 40., 50., 60., 70., 80., 100. };
    SetCentralityAxis(n, bins);
  }
  TIter    next(&fRingHistos);
  RingHistos* o = 0;
  while ((o = static_cast<RingHistos*>(next())))
    o->SetupForData(fEtaAxis, fCentralityAxis, fMaxE, 
		    fNEbins, fUseIncreasingBins);
}  
//____________________________________________________________________
void
AliFMDEnergyFitter::SetEtaAxis(const TAxis& eAxis)
{
  // 
  // Set the eta axis to use.  This will force the code to use this
  // eta axis definition - irrespective of whatever axis is passed to
  // the Init member function.  Therefore, this member function can be
  // used to force another eta axis than one found in the correction
  // objects. 
  // 
  // Parameters:
  //    etaAxis Eta axis to use 
  //
  SetEtaAxis(eAxis.GetNbins(),eAxis.GetXmin(),eAxis.GetXmax());
}
//____________________________________________________________________
void
AliFMDEnergyFitter::SetEtaAxis(Int_t nBins, Double_t etaMin, Double_t etaMax)
{
  // 
  // Set the eta axis to use.  This will force the code to use this
  // eta axis definition - irrespective of whatever axis is passed to
  // the Init member function.  Therefore, this member function can be
  // used to force another eta axis than one found in the correction
  // objects. 
  // 
  // Parameters:
  //    nBins  Number of bins 
  //    etaMin Minimum of the eta axis 
  //    etaMax Maximum of the eta axis 
  //
  fEtaAxis.Set(nBins,etaMin,etaMax);
}
//____________________________________________________________________
void
AliFMDEnergyFitter::SetCentralityAxis(UShort_t n, Double_t* bins)
{
  fCentralityAxis.Set(n-1, bins);
}
//____________________________________________________________________
void
AliFMDEnergyFitter::SetEnableDeltaShift(Bool_t use) 
{
  AliLandauGaus::EnableSigmaShift(use ? 1 : 0);
}

//____________________________________________________________________
Bool_t
AliFMDEnergyFitter::Accumulate(const AliESDFMD& input,
			       Double_t         cent,
			       Bool_t           empty)
{
  // 
  // Fitter the input AliESDFMD object
  // 
  // Parameters:
  //    input     Input 
  //    cent      Centrality 
  //    empty     Whether the event is 'empty'
  // 
  // Return:
  //    True on success, false otherwise 
  DGUARD(fDebug, 5, "Accumulate statistics in AliFMDEnergyFitter - cholm");
  Int_t icent = fCentralityAxis.FindBin(cent);
  if (icent < 1 || icent > fCentralityAxis.GetNbins()) icent = 0;

  UShort_t nFills = 0;
  for(UShort_t d = 1; d <= 3; d++) {
    Int_t nRings = (d == 1 ? 1 : 2);
    for (UShort_t q = 0; q < nRings; q++) {
      Char_t      r    = (q == 0 ? 'I' : 'O');
      UShort_t    nsec = (q == 0 ?  20 :  40);
      UShort_t    nstr = (q == 0 ? 512 : 256);

      RingHistos* histos = GetRingHistos(d, r);
      if (!histos) {
	AliWarningF("No histograms for FMD%d%c", d, r);
	continue;
      }
      
      for(UShort_t s = 0; s < nsec;  s++) {
	for(UShort_t t = 0; t < nstr; t++) {
	  Float_t mult = input.Multiplicity(d,r,s,t);
	  
	  // Keep dead-channel information. 
	  if (mult == AliESDFMD::kInvalidMult || mult > 10 || mult <= 0) {
	    DMSG(fDebug,4,"FMD%d%c[%2d,%3d]=%f is invalid, too large, or <0",
		 d, r, s, t);
	    continue;
	  }

	  // Get the pseudo-rapidity 
	  Double_t eta1 = input.Eta(d,r,s,t);

	  // Int_t ieta = fEtaAxis.FindBin(eta1);
	  // if (ieta < 1 || ieta >  fEtaAxis.GetNbins()) continue; 

	  // histos->Fill(empty, ieta-1, icent, mult);
	  histos->Fill(empty, eta1, icent, mult);
	  nFills++;
	} // for strip
      } // for sector
    } // for ring 
  } // for detector

  DMSG(fDebug, 3, "Found a total of %d signals for c=%f, and %sempty event", 
       nFills, cent, (empty ? "" : "non-"));
  return kTRUE;
}

//____________________________________________________________________
void
AliFMDEnergyFitter::Fit(const TList* dir)
{
  // 
  // Scale the histograms to the total number of events 
  // 
  // Parameters:
  //    dir Where the histograms are  
  //
  DGUARD(fDebug, 1, "Fit distributions in AliFMDEnergyFitter");
  if (!fDoFits) {
    AliInfo("Not asked to do fits, returning");
    return;
  }

  TList* d = static_cast<TList*>(dir->FindObject(GetName()));
  if (!d) {
    AliWarningF("No list named %s found in %s", GetName(), dir->GetName());
    return;
  }    

  // +1          for chi^2
  // +3          for 1 landau 
  // +1          for N 
  // +fNParticles-1 for weights 
  Int_t nStack = kN+fNParticles;
  THStack* stack[nStack]; 
  stack[0]         = new THStack("chi2",   "#chi^{2}/#nu");
  stack[kC     +1] = new THStack("c",      "Constant");
  stack[kDelta +1] = new THStack("delta",  "#Delta_{p}");
  stack[kXi    +1] = new THStack("xi",     "#xi");
  stack[kSigma +1] = new THStack("sigma",  "#sigma");
  stack[kSigmaN+1] = new THStack("sigman", "#sigma_{n}");
  stack[kN     +1] = new THStack("n",      "# Particles");
  for (Int_t i = 2; i <= fNParticles; i++) 
    stack[kN+i] = new THStack(Form("a%d", i), Form("a_{%d}", i));
  for (Int_t i = 0; i < nStack; i++) 
    d->Add(stack[i]);

  // If we have no ring histograms, re-init. 
  if (fRingHistos.GetEntries() <= 0) Init();

  AliInfoF("Will do fits for %d rings", fRingHistos.GetEntries());
  TIter    next(&fRingHistos);
  RingHistos* o = 0;
  while ((o = static_cast<RingHistos*>(next()))) {
    AliInfoF("Fill fit for FMD%d%c", o->fDet, o->fRing);
    if (CheckSkip(o->fDet, o->fRing, fSkips)) {
      AliWarningF("Skipping FMD%d%c for fitting", o->fDet, o->fRing);
      continue;
    }
    
    TObjArray* l = o->Fit(d, fLowCut, fNParticles,
			  fMinEntries, fFitRangeBinWidth,
			  fMaxRelParError, fMaxChi2PerNDF,
			  fMinWeight, fRegularizationCut,
			  fResidualMethod);
    if (!l) continue;
    for (Int_t i = 0; i < l->GetEntriesFast()-1; i++) { // Last is status 
      stack[i % nStack]->Add(static_cast<TH1*>(l->At(i))); 
    }
  }

  if (!fDoMakeObject) return;

  MakeCorrectionsObject(d);
}

//____________________________________________________________________
void
AliFMDEnergyFitter::MakeCorrectionsObject(TList* d)
{
  // 
  // Generate the corrections object 
  // 
  // Parameters:
  //    dir List to analyse 
  //
  DGUARD(fDebug, 1, "Make the correction objec in AliFMDEnergyFitter");
    
  AliFMDCorrELossFit* obj = new AliFMDCorrELossFit;
  obj->SetEtaAxis(fEtaAxis);
  obj->SetLowCut(fLowCut);
  if (AliLandauGaus::EnableSigmaShift()) 
    obj->SetBit(AliFMDCorrELossFit::kHasShift);

  TIter    next(&fRingHistos);
  RingHistos* o = 0;
  while ((o = static_cast<RingHistos*>(next()))) {
    if (CheckSkip(o->fDet, o->fRing, fSkips)) {
      AliWarningF("Skipping FMD%d%c for correction object", o->fDet, o->fRing);
      continue;
    }
    
    o->FindBestFits(d, *obj, fEtaAxis);
  }
  obj->IsGood();
  d->Add(obj, "elossFits");
}

//____________________________________________________________________
void
AliFMDEnergyFitter::SetDebug(Int_t dbg)
{
  // 
  // Set the debug level.  The higher the value the more output 
  // 
  // Parameters:
  //    dbg Debug level 
  //
  fDebug = dbg;
  TIter    next(&fRingHistos);
  RingHistos* o = 0;
  while ((o = static_cast<RingHistos*>(next())))
    o->fDebug = dbg;
}  

//____________________________________________________________________
namespace {
  template <typename T>
  void GetParam(Bool_t& ret, const TCollection* col, 
		const TString& name, T& var)
  {
    TObject* o = col->FindObject(name);	      
    if (o) AliForwardUtil::GetParameter(o,var); 
    else   ret = false;				
  }
}

//____________________________________________________________________
Bool_t
AliFMDEnergyFitter::ReadParameters(const TCollection* col)
{
  // Read parameters of this object from a collection
  //
  // Parameters:
  //    col   Collection to read parameters from 
  // 
  // Return value:
  //   true on success, false otherwise 
  //
  if (!col) return false;
  Bool_t ret  = true;
  TH1*   hist = static_cast<TH1*>(col->FindObject("etaAxis"));
  if (!hist) ret = false;
  else {
    TAxis* axis = hist->GetXaxis();
    if (axis->GetXbins()->GetArray()) 
      fEtaAxis.Set(axis->GetNbins(), axis->GetXbins()->GetArray());
   else 
      fEtaAxis.Set(axis->GetNbins(), axis->GetXmin(), axis->GetXmax());
  } 
  GetParam(ret,col,"lowCut",        fLowCut);
  GetParam(ret,col,"nParticles",    fNParticles);
  GetParam(ret,col,"minEntries",    fMinEntries);
  GetParam(ret,col,"subtractBins",  fFitRangeBinWidth);
  GetParam(ret,col,"doFits",        fDoFits);
  GetParam(ret,col,"doObject",      fDoMakeObject);
  GetParam(ret,col,"maxE",          fMaxE);
  GetParam(ret,col,"nEbins",        fNEbins);
  GetParam(ret,col,"increasingBins",fUseIncreasingBins);
  GetParam(ret,col,"maxRelPerError",fMaxRelParError);
  GetParam(ret,col,"maxChi2PerNDF", fMaxChi2PerNDF);
  GetParam(ret,col,"minWeight",     fMinWeight);
  Bool_t dummy;
  GetParam(dummy,col,"regCut",      fRegularizationCut);

  return ret;
}

//____________________________________________________________________
Bool_t
AliFMDEnergyFitter::CheckSkip(UShort_t d, Char_t r, UShort_t skips) 
{
  UShort_t q  = (r == 'I' || r == 'i' ? 0 : 1);
  UShort_t c = 1 << (d-1);
  UShort_t t = 1 << (c+q-1);

  return (t & skips) == t;
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
AliFMDEnergyFitter::Print(Option_t*) const
{
  // 
  // Print information
  // 
  // Parameters:
  //    option Not used 
  //
  AliForwardUtil::PrintTask(*this);

  gROOT->IncreaseDirLevel();
  PFV("Low cut [E/E_mip]",	fLowCut);
  PFV("Max(particles)",	        fNParticles);
  PFV("Min(entries)",	        fMinEntries);
  PFV("Fit range [bins]",       fFitRangeBinWidth);
  PFB("Make fits",              fDoFits);
  PFB("Make object",            fDoMakeObject);
  PFV("Max E/E_mip",	        fMaxE);
  PFV("N bins",	                fNEbins);
  PFB("Increasing bins",        fUseIncreasingBins);
  PFV("max(delta p/p)",    	fMaxRelParError);
  PFV("max(chi^2/nu)",	        fMaxChi2PerNDF);
  PFV("min(a_i)",	        fMinWeight);
  PFV("Regularization cut",     fRegularizationCut);
  TString r = "";
  switch (fResidualMethod) { 
  case kNoResiduals:              r = "None";       break;
  case kResidualDifference:       r = "Difference"; break;
  case kResidualScaledDifference: r = "Scaled difference"; break;
  case kResidualSquareDifference: r = "Square difference"; break;
  }
  PFV("Residuals", r);
  gROOT->DecreaseDirLevel();
}
  
//====================================================================
AliFMDEnergyFitter::RingHistos::RingHistos()
  : AliForwardUtil::RingHistos(), 
    fEDist(0), 
    fEmpty(0),
    // fEtaEDists(0), 
    fHist(0),
    fList(0),
    fBest(0),
    fFits("AliFMDCorrELossFit::ELossFit", 200),
    fDebug(0)
{
  // 
  // Default CTOR
  //
  DGUARD(fDebug, 3, "Default CTOR AliFMDEnergyFitter::RingHistos");
  fBest.Expand(0);
}

//____________________________________________________________________
AliFMDEnergyFitter::RingHistos::RingHistos(UShort_t d, Char_t r)
  : AliForwardUtil::RingHistos(d,r), 
    fEDist(0), 
    fEmpty(0),
    // fEtaEDists(0), 
    fHist(0),
    fList(0),
    fBest(0),
    fFits("AliFMDCorrELossFit::ELossFit", 200),
    fDebug(0)
{
  // 
  // Constructor
  // 
  // Parameters:
  //    d detector
  //    r ring 
  //
  DGUARD(fDebug, 3, "Named CTOR AliFMDEnergyFitter::RingHistos: FMD%d%c",
	 d, r);
  fBest.Expand(0);
}
//____________________________________________________________________
AliFMDEnergyFitter::RingHistos::~RingHistos()
{
  // 
  // Destructor 
  //
  DGUARD(fDebug, 3, "DTOR of AliFMDEnergyFitter::RingHistos");
  // if (fEDist) delete fEDist;
  // fEtaEDists.Delete();
}
//__________________________________________________________________
TArrayD
AliFMDEnergyFitter::RingHistos::MakeIncreasingAxis(Int_t n, Double_t min, 
						   Double_t max) const
{
  // 
  // Make an axis with increasing bins 
  // 
  // Parameters:
  //    n    Number of bins 
  //    min  Minimum 
  //    max  Maximum
  // 
  // Return:
  //    An axis with quadratically increasing bin size 
  //

  TArrayD bins(n+1);
  Double_t dx = (max-min) / n;
  bins[0]     = min;
  Int_t    i  = 1;
  for (i = 1; i < n+1; i++) {
    Double_t dI   = float(i)/n;
    Double_t next = bins[i-1] + dx + dI * dI * dx;
    bins[i]       = next;
    if (next > max) break;
  }
  bins.Set(i+1);
  return bins;
}

//____________________________________________________________________
TH2*
AliFMDEnergyFitter::RingHistos::Make(const char*  name, 
				     const char*  title, 
				     const TAxis& eAxis, 
				     Double_t     deMax, 
				     Int_t        nDeBins, 
				     Bool_t       incr)
{
  // 
  // Make E/E_mip histogram 
  // 
  // Parameters:
  //    deMax   Maximum energy loss 
  //    nDeBins Number energy loss bins 
  //    incr    Whether to make bins of increasing size
  //
  TH2* h = 0;
  TAxis mAxis;
  if (incr) {
    TArrayD deAxis = MakeIncreasingAxis(nDeBins, 0, deMax);
    mAxis.Set(deAxis.GetSize()-1, deAxis.GetArray());
  }
  else 
    mAxis.Set(nDeBins, 0, deMax);
  
  if (mAxis.GetXbins()->GetArray()) {
    // Variable bin since in Delta 
    if (eAxis.GetXbins()->GetArray()) {
      // Variadic bin size in eta 
      h = new TH2D(name, title, 
		   eAxis.GetNbins(), eAxis.GetXbins()->GetArray(),
		   mAxis.GetNbins(), mAxis.GetXbins()->GetArray());
    }
    else { 
      // Evenly spaced bins in eta
      h = new TH2D(name, title, 
		   eAxis.GetNbins(), eAxis.GetXmin(), eAxis.GetXmax(), 
		   mAxis.GetNbins(), mAxis.GetXbins()->GetArray());
    }
  }
  else { 
    // Make increasing bins axis 
    if (eAxis.GetXbins()->GetArray()) {
      // Variable size eta bins 
      h = new TH2D(name, title, 
		   eAxis.GetNbins(), eAxis.GetXbins()->GetArray(),
		   mAxis.GetNbins(), mAxis.GetXmin(), mAxis.GetXmax());
    }
    else {
      // Fixed size eta bins 
      h = new TH2D(name, title, 
		   eAxis.GetNbins(), eAxis.GetXmin(), eAxis.GetXmax(), 
		   mAxis.GetNbins(), mAxis.GetXmin(), mAxis.GetXmax());
    }
  }
  h->SetDirectory(0);
  h->SetYTitle("#sum#DeltaE/#DeltaE_{mip}");	
  h->SetXTitle("#eta");	
  h->SetFillColor(Color());
  h->SetMarkerColor(Color());
  h->SetLineColor(Color());
  h->SetFillStyle(3001);
  h->SetMarkerStyle(20);
  h->Sumw2();

  return h;
}
//____________________________________________________________________
void
AliFMDEnergyFitter::RingHistos::CreateOutputObjects(TList* dir)
{
  // 
  // Define outputs
  // 
  // Parameters:
  //    dir 
  //
  DGUARD(fDebug, 2, "Define output in AliFMDEnergyFitter::RingHistos");
  fList = DefineOutputList(dir);
  
  // fEtaEDists = new TList;
  // fEtaEDists->SetOwner();
  // fEtaEDists->SetName("EDists");
  // 
  // fList->Add(fEtaEDists);
}
//____________________________________________________________________
void
AliFMDEnergyFitter::RingHistos::SetupForData(const TAxis& eAxis, 
					     const TAxis& /* cAxis */,
					     Double_t     maxDE, 
					     Int_t        nDEbins, 
					     Bool_t       useIncrBin)
{
  // 
  // Initialise object 
  // 
  // Parameters:
  //    eAxis      Eta axis
  //    maxDE      Max energy loss to consider 
  //    nDEbins    Number of bins 
  //    useIncrBin Whether to use an increasing bin size 
  //
  DGUARD(fDebug, 2, "Initialize in AliFMDEnergyFitter::RingHistos");
  fEDist = new TH1D(Form("%s_edist", fName.Data()), 
		    Form("#sum#DeltaE/#DeltaE_{mip} for %s", fName.Data()), 
		    200, 0, 6);
  fEDist->SetXTitle("#sum#Delta/#Delta_{mip}");
  fEDist->SetFillColor(Color());
  fEDist->SetLineColor(Color());
  fEDist->SetMarkerColor(Color());
  fEDist->SetFillStyle(3001);
  fEDist->SetMarkerStyle(20);
  fEDist->Sumw2();
  fEDist->SetDirectory(0);

  fEmpty = static_cast<TH1D*>(fEDist->Clone(Form("%s_empty", fName.Data())));
  fEmpty->SetTitle(Form("#sum#DeltaE/#DeltaE_{mip} for %s (empty events)", 
			fName.Data()));
  fEmpty->SetDirectory(0);

  fList->Add(fEDist);
  fList->Add(fEmpty);
  fHist = Make("eloss", "#sum#Delta/#Delta_{mip}", eAxis, 
	       maxDE, nDEbins, useIncrBin);
  fList->Add(fHist);
  // fList->ls();
  // fEtaEDists.ls();
}

//____________________________________________________________________
void
AliFMDEnergyFitter::RingHistos::Fill(Bool_t   empty, 
				     Double_t eta, 
				     Int_t /* icent */,  
				     Double_t mult)
{
  // 
  // Fill histogram 
  // 
  // Parameters:
  //    empty  True if event is empty
  //    ieta   Eta bin (0-based)
  //    icent  Centrality bin (1-based)
  //    mult   Signal 
  //
  DGUARD(fDebug, 10, "Filling in AliFMDEnergyFitter::RingHistos");
  if (empty) { 
    fEmpty->Fill(mult);
    return;
  }
  fEDist->Fill(mult);

  // if (!fEtaEDists) { 
  if (!fHist) {
    Warning("Fill", "No list of E dists defined");
    return;
  }
  fHist->Fill(eta, mult);
  // Here, we should first look up in a centrality array, and then in
  // that array look up the eta bin - or perhaps we should do it the
  // other way around?
  // if (ieta < 0 || ieta >= fEtaEDists->GetEntries()) return;
  
  // TH1D* h = static_cast<TH1D*>(fEtaEDists->At(ieta));
  // if (!h) return;

  // h->Fill(mult);
}

//____________________________________________________________________
TH1* 
AliFMDEnergyFitter::RingHistos::MakePar(const char* name, 
					const char* title, 
					const TAxis& eta) const
{
  // 
  // Make a parameter histogram
  // 
  // Parameters:
  //    name   Name of histogram.
  //    title  Title of histogram. 
  //    eta    Eta axis 
  // 
  // Return:
  //    
  //
  TH1D* h = new TH1D(Form("%s_%s", fName.Data(), name), 
		     Form("%s for %s", title, fName.Data()), 
		     eta.GetNbins(), eta.GetXmin(), eta.GetXmax());
  h->SetXTitle("#eta");
  h->SetYTitle(title);
  h->SetDirectory(0);
  h->SetFillColor(Color());
  h->SetMarkerColor(Color());
  h->SetLineColor(Color());
  h->SetFillStyle(3001);
  h->Sumw2();

  return h;
}
//____________________________________________________________________
TH1*
AliFMDEnergyFitter::RingHistos::MakeTotal(const char* name, 
					  const char* title, 
					  const TAxis& eta, 
					  Int_t low, 
					  Int_t high, 
					  Double_t val, 
					  Double_t err) const
{
  // 
  // Make a histogram that contains the results of the fit over the full ring 
  // 
  // Parameters:
  //    name  Name 
  //    title Title
  //    eta   Eta axis 
  //    low   Least bin
  //    high  Largest bin
  //    val   Value of parameter 
  //    err   Error on parameter 
  // 
  // Return:
  //    The newly allocated histogram 
  //
  Double_t xlow  = eta.GetBinLowEdge(low);
  Double_t xhigh = eta.GetBinUpEdge(high);
  TH1D* h = new TH1D(Form("%s_%s", fName.Data(), name), 
		     Form("%s for %s", title, fName.Data()), 
		     1, xlow, xhigh);
  h->SetBinContent(1, val);
  h->SetBinError(1, err);
  h->SetXTitle("#eta");
  h->SetYTitle(title);
  h->SetDirectory(0);
  h->SetFillColor(0);
  h->SetMarkerColor(Color());
  h->SetLineColor(Color());
  h->SetFillStyle(0);
  h->SetLineStyle(2);
  h->SetLineWidth(2);

  return h;
}
		     
//____________________________________________________________________
TObjArray*
AliFMDEnergyFitter::RingHistos::Fit(TList*           dir, 
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
  return FitSlices(dir, "eloss", lowCut, nParticles, minEntries, 
		   minusBins, relErrorCut, chi2nuCut, minWeight, regCut, 
		   residuals, true, &fBest);
}

//____________________________________________________________________
TObjArray*
AliFMDEnergyFitter::RingHistos::FitSlices(TList*           dir, 
					  const char*      name, 
					  Double_t         lowCut, 
					  UShort_t         nParticles,
					  UShort_t         minEntries,
					  UShort_t         minusBins, 
					  Double_t         relErrorCut, 
					  Double_t         chi2nuCut,
					  Double_t         minWeight,
					  Double_t         regCut,
					  EResidualMethod  residuals,
					  Bool_t           scaleToPeak,
					  TObjArray*       best) const
{
  // 
  // Fit each histogram to up to @a nParticles particle responses.
  // 
  // Parameters:
  //    dir         Output list 
  //    eta         Eta axis 
  //    lowCut      Lower cut 
  //    nParticles  Max number of convolved landaus to fit
  //    minEntries  Minimum number of entries 
  //    minusBins   Number of bins from peak to subtract to 
  //                    get the fit range 
  //    relErrorCut Cut applied to relative error of parameter. 
  //                    Note, for multi-particle weights, the cut 
  //                    is loosend by a factor of 2 
  //    chi2nuCut   Cut on @f$ \chi^2/\nu@f$ - 
  //                    the reduced @f$\chi^2@f$ 
  //
  DGUARD(fDebug, 2, "Fit in AliFMDEnergyFitter::RingHistos");

  // Get our ring list from the output container 
  TList* l = GetOutputList(dir);
  if (!l) return 0; 

  TList* dists = 0;
  // Get the 2D histogram 
  TH2* h = static_cast<TH2*>(l->FindObject(name));
  if (!h) { 
    AliWarningF("Didn't find 2D histogram '%s' in %s", name, l->GetName());
    // Get the energy distributions from the output container 
    dists = static_cast<TList*>(l->FindObject("EDists"));
    if (!dists) { 
      AliWarningF("Didn't find EtaEDists (%s) in %s", 
		  fName.Data(), l->GetName());
      l->ls();
      return 0;
    }
  }
  if (!h && !dists) return 0;

  const TAxis* pEta = (h ? h->GetXaxis() : 
			    static_cast<TAxis*>(dir->FindObject("etaAxis")));
  if (!pEta) { 
    AliWarningF("Didn't find the eta axis - either from histogram %p or "
		"list %p (%s)", h, dir, (dir ? dir->GetName() : "-"));
    return 0;
  }
  const TAxis& eta = *pEta;

  // Create an output list for the fitted distributions 
  TList* out = new TList;
  out->SetOwner();
  out->SetName(Form("%sDists", name));
  l->Add(out);

  // Optional container for residuals 
  TList* resi = 0;
  if (residuals != kNoResiduals) {
    resi = new TList();
    resi->SetName(Form("%sResiduals", name));
    resi->SetOwner();
    l->Add(resi);
  }

  // Container of the fit results histograms 
  TObjArray* pars  = new TObjArray(3+nParticles+1);
  pars->SetName(Form("%sResults", name));
  l->Add(pars);

  // Result objects stored in separate list on the output 
  TH1* hChi2   = 0;
  TH1* hC      = 0;
  TH1* hDelta  = 0;
  TH1* hXi     = 0;
  TH1* hSigma  = 0;
  TH1* hSigmaN = 0;
  TH1* hN      = 0;
  TH1* hA[nParticles-1];
  pars->Add(hChi2   = MakePar("chi2",   "#chi^{2}/#nu", eta));
  pars->Add(hC      = MakePar("c",      "Constant",     eta));
  pars->Add(hDelta  = MakePar("delta",  "#Delta_{p}",   eta));
  pars->Add(hXi     = MakePar("xi",     "#xi",          eta));
  pars->Add(hSigma  = MakePar("sigma",  "#sigma",       eta));
  pars->Add(hSigmaN = MakePar("sigman", "#sigma_{n}",   eta));
  pars->Add(hN      = MakePar("n",      "N", eta));
  for (UShort_t i = 1; i < nParticles; i++) 
    pars->Add(hA[i-1] = MakePar(Form("a%d",i+1), Form("a_{%d}",i+1), eta));

  
  Int_t nDists = h ? h->GetNbinsX() : dists->GetEntries();
  Int_t low    = nDists;
  Int_t high   = 0;
  Int_t nEmpty = 0;
  Int_t nLow   = 0;
  Int_t nFitted= 0;
  if (best) {
    best->Expand(nDists+1);
    best->Clear();
    best->SetOwner(false);
  }
  for (Int_t i = 0; i < nDists; i++) { 
    // Ignore empty histograms altoghether 
    Int_t b    = i+1;
    TH1D* dist = (h ? h->ProjectionY(Form(fgkEDistFormat,GetName(),b),b,b,"e") 
		  : static_cast<TH1D*>(dists->At(i)));
    if (!dist) { 
      // If we got the null pointer, return 0
      nEmpty++;
      continue;
    }
    // Then releasing the histogram from the it's directory
    dist->SetDirectory(0);
    // Set a meaningful title
    dist->SetTitle(Form("#Delta/#Delta_{mip} for %s in %6.2f<#eta<%6.2f",
			GetName(), eta.GetBinLowEdge(b),
			eta.GetBinUpEdge(b)));

    // Now fit 
    UShort_t    status1 = 0;
    ELossFit_t* res     = FitHist(dist,
				  lowCut, 
				  nParticles,
				  minEntries,
				  minusBins,   
				  relErrorCut,
				  chi2nuCut,
				  minWeight,
				  regCut,
				  scaleToPeak,
				  status1);
    if (!res) {
      switch (status1) { 
      case 1: nEmpty++; break;
      case 2: nLow++;   break;
      }
      // Only clean up if we have no input list 
      if (h) delete dist;
      continue;
    }
      
    // Now count as fitted, store as best fits, and add distribution
    // to the dedicated output list
    nFitted++;
    res->fBin = b; // We only have the bin information here 
    if (best) best->AddAt(res, b);
    out->Add(dist);
    // dist->GetListOfFunctions()->Add(res);
    
    // If asked to calculate residuals, do so, and store result on the
    // dedicated output list
    if (residuals != kNoResiduals && resi) 
      CalculateResiduals(residuals, lowCut, dist, res, resi);

    // Store eta limits 
    low   = TMath::Min(low,b);
    high  = TMath::Max(high,b);

    // Get the reduced chi square
    Double_t chi2 = res->fChi2; // GetChisquare();
    Int_t    ndf  = res->fNu;   // GetNDF();
    
    // Store results of best fit in output histograms 
    // res->SetLineWidth(4);
    hChi2   ->SetBinContent(b, ndf > 0 ? chi2 / ndf : 0);
    hC      ->SetBinContent(b, res->fC); 
    hDelta  ->SetBinContent(b, res->fDelta); 
    hXi     ->SetBinContent(b, res->fXi); 
    hSigma  ->SetBinContent(b, res->fSigma); 
    hSigmaN ->SetBinContent(b, res->fSigmaN); 
    hN      ->SetBinContent(b, res->fN); 

    hC     ->SetBinError(b, res->fEC); 
    hDelta ->SetBinError(b, res->fEDelta);
    hXi    ->SetBinError(b, res->fEXi); 
    hSigma ->SetBinError(b, res->fESigma); 
    hSigmaN->SetBinError(b, res->fESigmaN); 
    // hN     ->SetBinError(b, res->fGetParError(kN));

    for (UShort_t j = 0; j < nParticles-1; j++) {
      hA[j]->SetBinContent(b, res->fA[j]); 
      hA[j]->SetBinError(b, res->fEA[j]); 
    }
  }
  printf("%s: Out of %d histograms, %d where empty, %d had too little data,"
	 "leaving %d to be fitted, of which %d succeeded\n",  
	 GetName(), nDists, nEmpty, nLow, nDists-nEmpty-nLow, nFitted);

  // Fit the full-ring histogram 
  TH1*        total   = GetOutputHist(l, Form("%s_edist", fName.Data()));
  if (total) {
    UShort_t    statusT = 0;
    ELossFit_t* resT    = FitHist(total,
				  lowCut, 
				  nParticles,
				  minEntries, 
				  minusBins,
				  relErrorCut,
				  chi2nuCut,
				  minWeight,
				  regCut,
				  scaleToPeak,
				  statusT);
    if (resT) { 
      // Make histograms for the result of this fit 
      Double_t chi2 = resT->GetChi2();
      Int_t    ndf  = resT->GetNu();
      pars->Add(MakeTotal("t_chi2",   "#chi^{2}/#nu", eta, low, high,
			  ndf > 0 ? chi2/ndf : 0, 0));
      pars->Add(MakeTotal("t_c",      "Constant",     eta, low, high,
			  resT->GetC(), resT->GetEC()));
      pars->Add(MakeTotal("t_delta",  "#Delta_{p}",   eta, low, high,
			  resT->GetDelta(), resT->GetEDelta()));
      pars->Add(MakeTotal("t_xi",     "#xi",          eta, low, high,
			  resT->GetXi(), resT->GetEXi()));
      pars->Add(MakeTotal("t_sigma",  "#sigma",       eta, low, high,
			  resT->GetSigma(), resT->GetESigma()));
      pars->Add(MakeTotal("t_sigman", "#sigma_{n}",   eta, low, high,
			  resT->GetSigmaN(), resT->GetESigmaN()));
      pars->Add(MakeTotal("t_n",    "N",              eta, low, high,
			  resT->GetN(), 0));
      for (UShort_t j = 0; j < nParticles-1; j++) 
	pars->Add(MakeTotal(Form("t_a%d",j+2), 
			    Form("a_{%d}",j+2), eta, low, high,
			    resT->GetA(j), resT->GetEA(j)));
    }
  }

  TH1* status = new TH1I(Form("%sStatus",name), "Status of Fits", 5, 0, 5);
  status->GetXaxis()->SetBinLabel(1, "Total");
  status->GetXaxis()->SetBinLabel(2, "Empty");
  status->GetXaxis()->SetBinLabel(3, Form("<%d", minEntries));
  status->GetXaxis()->SetBinLabel(4, "Candidates");
  status->GetXaxis()->SetBinLabel(5, "Fitted");
  status->SetXTitle("Status");
  status->SetYTitle("# of #Delta distributions");
  status->SetBinContent(1, nDists);
  status->SetBinContent(2, nEmpty);
  status->SetBinContent(3, nLow);
  status->SetBinContent(4, nDists-nLow-nEmpty);
  status->SetBinContent(5, nFitted);
  status->SetFillColor(Color());
  status->SetFillStyle(3001);
  status->SetLineColor(Color());
  status->SetDirectory(0);
  status->SetStats(0);
  pars->Add(status);
  return pars;
}


//____________________________________________________________________
void
AliFMDEnergyFitter::RingHistos::Scale(TH1* dist) const
{
  // Scale to the bin-width
  dist->Scale(1., "width");
}  


//____________________________________________________________________
AliFMDEnergyFitter::RingHistos::ELossFit_t*
AliFMDEnergyFitter::RingHistos::FitHist(TH1*      dist,
					Double_t  lowCut, 
					UShort_t  nParticles, 
					UShort_t  minEntries,
					UShort_t  minusBins, 
					Double_t  relErrorCut, 
					Double_t  chi2nuCut,
					Double_t  minWeight,
					Double_t  regCut,
					Bool_t    scaleToPeak,
					UShort_t& status) const
{
  // 
  // Fit a signal histogram.  First, the bin @f$ b_{min}@f$ with
  // maximum bin content in the range @f$ [E_{min},\infty]@f$ is
  // found.  Then the fit range is set to the bin range 
  // @f$ [b_{min}-\Delta b,b_{min}+2\Delta b]@f$, and a 1 
  // particle signal is fitted to that.  The parameters of that fit 
  // is then used as seeds for a fit of the @f$ N@f$ particle response 
  // to the data in the range 
  // @f$ [b_{min}-\Delta b,N(\Delta_1+\xi_1\log(N))+2N\xi@f$
  // 
  // Parameters:
  //    dist        Histogram to fit 
  //    lowCut      Lower cut @f$ E_{min}@f$ on signal 
  //    nParticles  Max number @f$ N@f$ of convolved landaus to fit
  //    minusBins   Number of bins @f$ \Delta b@f$ from peak to 
  //                    subtract to get the fit range 
  //    relErrorCut Cut applied to relative error of parameter. 
  //                    Note, for multi-particle weights, the cut 
  //                    is loosend by a factor of 2 
  //    chi2nuCut   Cut on @f$ \chi^2/\nu@f$ - 
  //                    the reduced @f$\chi^2@f$ 
  // 
  // Return:
  //    The best fit function 
  //
  DGUARD(fDebug, 2, "Fit histogram in AliFMDEnergyFitter::RingHistos: %s",
	 dist->GetName());
  Double_t maxRange = 10;


  if (dist->GetEntries() <= 0) { 
    status = 1; // `empty'
    return 0;
  }
  Scale(dist);
  
  // Narrow search window for the peak 
  Int_t    cutBin  = TMath::Max(dist->GetXaxis()->FindBin(lowCut),3);
  Int_t    maxBin  = TMath::Min(dist->GetXaxis()->FindBin(10),
				dist->GetNbinsX());
  dist->GetXaxis()->SetRange(cutBin, maxBin);
    
  // Get the bin with maximum 
  Int_t    peakBin = dist->GetMaximumBin();
    
  // Normalize peak to 1 
  // Double_t max = dist->GetMaximum(); 
  Double_t max = dist->GetBinContent(peakBin); // Maximum(); 
  if (max <= 0) {
    status = 1; // `empty'
    return 0;
  }
  if (scaleToPeak) dist->Scale(1/max);
  DMSG(fDebug,5,"max(%s) -> %f", dist->GetName(), max);

  // Check that we have enough entries 
  Double_t nEntries = dist->GetEntries();
  if (nEntries <= minEntries) { 
    AliWarning(Form("Histogram at %s has too few entries (%f <= %d)",
		    dist->GetName(), nEntries, minEntries));
    status = 2;
    return 0;
  }

  // Create a fitter object 
  AliLandauGausFitter f(lowCut, maxRange, minusBins); 
  f.Clear();
  f.SetDebug(fDebug > 3); 

  // regularization cut - should be a parameter of the class 
  if (dist->GetEntries() > regCut) { 
    // We should rescale the errors 
    Double_t s = TMath::Sqrt(dist->GetEntries() / regCut);
    if (fDebug > 2) printf("Error scale: %f ", s);
    for (Int_t i = 1; i <= dist->GetNbinsX(); i++) {
      Double_t e = dist->GetBinError(i);
      dist->SetBinError(i, e * s);
    }
  }
  // If we are only asked to fit a single particle, return this fit, 
  // no matter what. 
  if (nParticles == 1) {
    TF1* r = f.Fit1Particle(dist, 0);
    if (!r) {
      status = 3; // No-fit
      return 0;
    }
    TF1* ff = new TF1(*r);
    dist->GetListOfFunctions()->Add(ff);
    ELossFit_t* ret = new ELossFit_t(0, *ff);
    ret->CalculateQuality(chi2nuCut, relErrorCut, minWeight);
    status = 0; // OK
    return ret;
  }

  // Fit from 2 upto n particles  
  for (Int_t i = 2; i <= nParticles; i++) f.FitNParticle(dist, i, 0);

  // Now, we need to select the best fit 
  Int_t nFits = f.GetFitResults().GetEntriesFast();
  for (Int_t i = nFits-1; i >= 0; i--) { 
    TF1* ff = static_cast<TF1*>(f.GetFunctions().At(i));
    // ff->SetRange(0, maxRange);
    dist->GetListOfFunctions()->Add(new TF1(*ff));
  }
  status = 0; // OK

  // Here, we use the real quality assesor instead of the old
  // `CheckResult' to ensure consitency in all output.
  ELossFit_t* ret = FindBestFit(dist, relErrorCut, chi2nuCut, minWeight);
  if (!ret) status = 3;
  return ret;
}

//__________________________________________________________________
AliFMDEnergyFitter::RingHistos::ELossFit_t* 
AliFMDEnergyFitter::RingHistos::FindBestFit(const TH1* dist,
					    Double_t   relErrorCut, 
					    Double_t   chi2nuCut,
					    Double_t   minWeightCut)  const
{
  // 
  // Find the best fit 
  // 
  // Parameters:
  //    dist           Histogram 
  //    relErrorCut    Cut applied to relative error of parameter. 
  //                       Note, for multi-particle weights, the cut 
  //                       is loosend by a factor of 2 
  //    chi2nuCut      Cut on @f$ \chi^2/\nu@f$ - 
  //                       the reduced @f$\chi^2@f$ 
  //    minWeightCut   Least valid @f$ a_i@f$ 
  // 
  // Return:
  //    Best fit or null
  //
  TList* funcs = dist->GetListOfFunctions();
  TF1*   func  = 0;
  Int_t  i     = 0;
  TIter  next(funcs);
  fFits.Clear(); // This is only ever used here

  if (fDebug) printf("Find best fit for %s ... ", dist->GetName());
  if (fDebug > 2) printf("\n");

  // Loop over all functions stored in distribution, 
  // and calculate the quality 
  while ((func = static_cast<TF1*>(next()))) { 
    ELossFit_t* fit = new(fFits[i++]) ELossFit_t(0,*func);
    fit->fDet  = fDet;
    fit->fRing = fRing;
    // fit->fBin  = b;
    fit->CalculateQuality(chi2nuCut, relErrorCut, minWeightCut);
    if (fDebug > 2) 
      Printf("%10s: %3d (chi^2/nu: %6.3f)", 
	     func->GetName(), fit->fQuality, 
	     (fit->fNu > 0 ? fit->fChi2 / fit->fNu : 999));
  }

  // Sort all the found fit objects in increasing quality 
  fFits.Sort();
  if (fDebug > 2) fFits.Print("s");

  // Get the top-most fit
  ELossFit_t* ret = static_cast<ELossFit_t*>(fFits.At(i-1));
  if (!ret) {
    AliWarningF("No fit found for %s", GetName());
    return 0;
  }
  if (ret && fDebug > 0) {
    if (fDebug > 1) printf(" %d: ", i-1);
    ret->Print("s");
  }
  // We have to make a copy here, because other wise the clones array
  // will overwrite the address
  return new ELossFit_t(*ret);
}

//____________________________________________________________________
void
AliFMDEnergyFitter::RingHistos::CalculateResiduals(EResidualMethod mode, 
						   Double_t        lowCut,
						   TH1*            dist, 
						   ELossFit_t*     fit, 
						   TCollection*    out) const
{

  // Clone the input, and reset
  TH1*    resi = static_cast<TH1*>(dist->Clone());
  TString tit(resi->GetTitle());
  tit.ReplaceAll("#DeltaE/#DeltaE_{mip}", "Residuals");
  resi->SetTitle(tit);
  resi->SetDirectory(0);

  // Set title on Y axis
  switch (mode) { 
  case kResidualDifference:       
    resi->SetYTitle("h_{i}-f(#Delta_{i}) #pm #delta_{i}");
    break;
  case kResidualScaledDifference:  
    resi->SetYTitle("[h_{i}-f(#Delta_{i})]/#delta_{i}"); break;
  case kResidualSquareDifference:  
    resi->SetYTitle("#chi_{i}^{2}=[h_{i}-f(#Delta_{i})]^{2}/#delta^{2}_{i}");
    break;
  default: 
    resi->SetYTitle("Unknown");
    break;
  }
  out->Add(resi);

  // Try to find the function 
  Double_t highCut = dist->GetXaxis()->GetXmax();
  TString funcName("landau1");
  if (fit->GetN() > 1) 
    funcName = Form("nlandau%d", fit->GetN());
  TF1* func = dist->GetFunction(funcName);
  if (func) func->GetRange(lowCut, highCut);
  resi->Reset("ICES");
  resi->GetListOfFunctions()->Clear();
  resi->SetUniqueID(mode);

    // Reset histogram
  Int_t nX = resi->GetNbinsX();
  for (Int_t i  = 1; i <= nX; i++) { 
    Double_t x  = dist->GetBinCenter(i);
    if (x < lowCut)  continue;
    if (x > highCut) break;

    Double_t h  = dist->GetBinContent(i);
    Double_t e  = dist->GetBinError(i);
    Double_t r  = 0;
    Double_t er = 0;
    if (h > 0 && e > 0) { 
      Double_t f = fit->GetC() * fit->Evaluate(x);
      if (f > 0) { 
	r  = h-f;
	switch (mode) { 
	case kResidualDifference: er = e; break;
	case kResidualScaledDifference:  r /= e; break;
	case kResidualSquareDifference:  r *= r; r /= (e*e); break;
	default: r = 0; break;
	}
      }
    }
    resi->SetBinContent(i, r);
    resi->SetBinError(i, er);
  }  
}

//__________________________________________________________________
void
AliFMDEnergyFitter::RingHistos::FindBestFits(const TList*        d, 
					     AliFMDCorrELossFit& obj,
					     const TAxis&        eta)
{
  // 
  // Find the best fits 
  // 
  // Parameters:
  //    d              Parent list
  //    obj            Object to add fits to
  //    eta            Eta axis 
  //    relErrorCut    Cut applied to relative error of parameter. 
  //                       Note, for multi-particle weights, the cut 
  //                       is loosend by a factor of 2 
  //    chi2nuCut      Cut on @f$ \chi^2/\nu@f$ - 
  //                       the reduced @f$\chi^2@f$ 
  //    minWeightCut   Least valid @f$ a_i@f$ 
  //

  // Get our ring list from the output container 
  TList* l = GetOutputList(d);
  if (!l) return; 

  // Get the energy distributions from the output container 
  TList* dists = static_cast<TList*>(l->FindObject("elossDists"));
  if (!dists) { 
    AliWarningF("Didn't find elossDists in %s", l->GetName());
    l->ls();
    return;
  }
  Int_t nBin = eta.GetNbins();
  if (fBest.GetEntriesFast() <= 0) { 
    AliWarningF("No fits found for %s", GetName());
    return;
  }
  
  for (Int_t b = 1; b <= nBin; b++) { 
    ELossFit_t* best = static_cast<ELossFit_t*>(fBest.At(b));
    if (!best) { 
      // AliErrorF("No best fit found @ %d for %s", b, GetName());
      continue;
    }
    // FindBestFit(b, dist, relErrorCut, chi2nuCut, minWeightCut);
    best->fDet  = fDet; 
    best->fRing = fRing;
    best->fBin  = b; // 
    if (fDebug > 0) {
      printf("Bin # %3d: ", b); 
      best->Print("s");
    }
    // Double_t eta = fAxis->GetBinCenter(b);
    obj.SetFit(fDet, fRing, b, best); 
    // new ELossFit_t(*best));
  }
}



//____________________________________________________________________
//
// EOF
//
