//
// Class to do the energy correction of FMD ESD data
//
#include "AliFMDEnergyFitter.h"
#include "AliForwardCorrectionManager.h"
#include <AliESDFMD.h>
#include <TAxis.h>
#include <TList.h>
#include <TH1.h>
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
    fMinEntries(1000),
    fFitRangeBinWidth(4),
    fDoFits(false),
    fDoMakeObject(false),
    fEtaAxis(),
    fCentralityAxis(),
    fMaxE(10),
    fNEbins(300), 
    fUseIncreasingBins(true),
    fMaxRelParError(.25),
    fMaxChi2PerNDF(20), 
    fMinWeight(1e-7),
    fDebug(0)
{
  // 
  // Default Constructor - do not use 
  //
}

//____________________________________________________________________
AliFMDEnergyFitter::AliFMDEnergyFitter(const char* title)
  : TNamed("fmdEnergyFitter", title), 
    fRingHistos(), 
    fLowCut(0.4),
    fNParticles(5),
    fMinEntries(1000),
    fFitRangeBinWidth(4),
    fDoFits(false),
    fDoMakeObject(false),
    fEtaAxis(0,0,0),
    fCentralityAxis(0,0,0),
    fMaxE(10),
    fNEbins(300), 
    fUseIncreasingBins(true),
    fMaxRelParError(.25),
    fMaxChi2PerNDF(20),  
    fMinWeight(1e-7),
    fDebug(3)
{
  // 
  // Constructor 
  // 
  // Parameters:
  //    title Title of object  - not significant 
  //
  fEtaAxis.SetName("etaAxis");
  fEtaAxis.SetTitle("#eta");
  fCentralityAxis.SetName("centralityAxis");
  fCentralityAxis.SetTitle("Centrality [%]");
  fRingHistos.Add(new RingHistos(1, 'I'));
  fRingHistos.Add(new RingHistos(2, 'I'));
  fRingHistos.Add(new RingHistos(2, 'O'));
  fRingHistos.Add(new RingHistos(3, 'I'));
  fRingHistos.Add(new RingHistos(3, 'O'));
}

//____________________________________________________________________
AliFMDEnergyFitter::AliFMDEnergyFitter(const AliFMDEnergyFitter& o)
  : TNamed(o), 
    fRingHistos(), 
    fLowCut(o.fLowCut),
    fNParticles(o.fNParticles),
    fMinEntries(o.fMinEntries),
    fFitRangeBinWidth(o.fFitRangeBinWidth),
    fDoFits(o.fDoFits),
    fDoMakeObject(o.fDoMakeObject),
    fEtaAxis(o.fEtaAxis),
    fCentralityAxis(o.fCentralityAxis),
    fMaxE(o.fMaxE),
    fNEbins(o.fNEbins), 
    fUseIncreasingBins(o.fUseIncreasingBins),
    fMaxRelParError(o.fMaxRelParError),
    fMaxChi2PerNDF(o.fMaxChi2PerNDF), 
    fMinWeight(o.fMinWeight),
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
AliFMDEnergyFitter::~AliFMDEnergyFitter()
{
  // 
  // Destructor
  //
  fRingHistos.Delete();
}

//____________________________________________________________________
AliFMDEnergyFitter&
AliFMDEnergyFitter::operator=(const AliFMDEnergyFitter& o)
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

  fLowCut        = o.fLowCut;
  fNParticles       = o.fNParticles;
  fMinEntries    = o.fMinEntries;
  fFitRangeBinWidth= o.fFitRangeBinWidth;
  fDoFits        = o.fDoFits;
  fDoMakeObject  = o.fDoMakeObject;
  fEtaAxis.Set(o.fEtaAxis.GetNbins(),
	       o.fEtaAxis.GetXmin(),
	       o.fEtaAxis.GetXmax());
  if (o.fCentralityAxis.GetXbins()) {
    const TArrayD* bins = o.fCentralityAxis.GetXbins();
    fCentralityAxis.Set(bins->GetSize()-1,bins->GetArray());
  }
  else 
    fCentralityAxis.Set(o.fCentralityAxis.GetNbins(),
			o.fCentralityAxis.GetXmin(),
			o.fCentralityAxis.GetXmax());
  fDebug         = o.fDebug;
  fMaxRelParError= o.fMaxRelParError;
  fMaxChi2PerNDF = o.fMaxChi2PerNDF;
  fMinWeight     = o.fMinWeight;

  fRingHistos.Delete();
  TIter    next(&o.fRingHistos);
  TObject* obj = 0;
  while ((obj = next())) fRingHistos.Add(obj);
  
  return *this;
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
AliFMDEnergyFitter::Init(const TAxis& eAxis)
{
  // 
  // Initialise the task
  // 
  // Parameters:
  //    etaAxis The eta axis to use.  Note, that if the eta axis
  // has already been set (using SetEtaAxis), then this parameter will be 
  // ignored
  //
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
    o->Init(fEtaAxis, fCentralityAxis, fMaxE, fNEbins, fUseIncreasingBins);
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
  //
  Int_t icent = fCentralityAxis.FindBin(cent);
  if (icent < 1 || icent > fCentralityAxis.GetNbins()) icent = 0;

  for(UShort_t d = 1; d <= 3; d++) {
    Int_t nRings = (d == 1 ? 1 : 2);
    for (UShort_t q = 0; q < nRings; q++) {
      Char_t      r    = (q == 0 ? 'I' : 'O');
      UShort_t    nsec = (q == 0 ?  20 :  40);
      UShort_t    nstr = (q == 0 ? 512 : 256);

      RingHistos* histos = GetRingHistos(d, r);
      
      for(UShort_t s = 0; s < nsec;  s++) {
	for(UShort_t t = 0; t < nstr; t++) {
	  Float_t mult = input.Multiplicity(d,r,s,t);
	  
	  // Keep dead-channel information. 
	  if(mult == AliESDFMD::kInvalidMult || mult > 10 || mult <= 0) 
	    continue;

	  // Get the pseudo-rapidity 
	  Double_t eta1 = input.Eta(d,r,s,t);
	  Int_t ieta = fEtaAxis.FindBin(eta1);
	  if (ieta < 1 || ieta >  fEtaAxis.GetNbins()) continue; 

	  histos->Fill(empty, ieta-1, icent, mult);

	} // for strip
      } // for sector
    } // for ring 
  } // for detector

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
  if (!fDoFits) return;

  TList* d = static_cast<TList*>(dir->FindObject(GetName()));
  if (!d) return;

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

  TIter    next(&fRingHistos);
  RingHistos* o = 0;
  while ((o = static_cast<RingHistos*>(next()))) {
    TObjArray* l = o->Fit(d, fEtaAxis, fLowCut, fNParticles,
			  fMinEntries, fFitRangeBinWidth,
			  fMaxRelParError, fMaxChi2PerNDF);
    if (!l) continue;
    for (Int_t i = 0; i < l->GetEntriesFast(); i++) { 
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
  AliForwardCorrectionManager& mgr = AliForwardCorrectionManager::Instance();
    
  AliFMDCorrELossFit* obj = new AliFMDCorrELossFit;
  obj->SetEtaAxis(fEtaAxis);
  obj->SetLowCut(fLowCut);

  TIter    next(&fRingHistos);
  RingHistos* o = 0;
  while ((o = static_cast<RingHistos*>(next()))) {
    o->FindBestFits(d, *obj, fEtaAxis, fMaxRelParError, 
		    fMaxChi2PerNDF, fMinWeight);
  }
  
  TString oName(mgr.GetObjectName(AliForwardCorrectionManager::kELossFits));
  TString fileName(mgr.GetFilePath(AliForwardCorrectionManager::kELossFits));
  AliInfo(Form("Object %s created in output - should be extracted and copied "
	       "to %s", oName.Data(), fileName.Data()));
  d->Add(new TNamed("filename", fileName.Data()));
  d->Add(obj, oName.Data());
}

//____________________________________________________________________
void
AliFMDEnergyFitter::DefineOutput(TList* dir)
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
  d->Add(&fEtaAxis);
  TIter    next(&fRingHistos);
  RingHistos* o = 0;
  while ((o = static_cast<RingHistos*>(next()))) {
    o->Output(d);
  }
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
void
AliFMDEnergyFitter::Print(Option_t*) const
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
	    << ind << " Low cut:                " << fLowCut << " E/E_mip\n"
	    << ind << " Max(particles):         " << fNParticles << '\n'
	    << ind << " Min(entries):           " << fMinEntries << '\n'
	    << ind << " Fit range:              " 
	    << fFitRangeBinWidth << " bins\n"
	    << ind << " Make fits:              " 
	    << (fDoFits ? "yes\n" : "no\n")
	    << ind << " Make object:            " 
	    << (fDoMakeObject ? "yes\n" : "no\n")
	    << ind << " Max E/E_mip:            " << fMaxE << '\n'
	    << ind << " N bins:                 " << fNEbins << '\n'
	    << ind << " Increasing bins:        " 
	    << (fUseIncreasingBins ?"yes\n" : "no\n")
	    << ind << " max(delta p/p):         " << fMaxRelParError << '\n'
	    << ind << " max(chi^2/nu):          " << fMaxChi2PerNDF << '\n'
	    << ind << " min(a_i):               " << fMinWeight << std::endl;
}
  
//====================================================================
AliFMDEnergyFitter::RingHistos::RingHistos()
  : AliForwardUtil::RingHistos(), 
    fEDist(0), 
    fEmpty(0),
    fEtaEDists(), 
    fList(0),
    fFits("AliFMDCorrELossFit::ELossFit"),
    fDebug(0)
{
  // 
  // Default CTOR
  //
}

//____________________________________________________________________
AliFMDEnergyFitter::RingHistos::RingHistos(UShort_t d, Char_t r)
  : AliForwardUtil::RingHistos(d,r), 
    fEDist(0), 
    fEmpty(0),
    fEtaEDists(), 
    fList(0),
    fFits("AliFMDCorrELossFit::ELossFit"),
    fDebug(0)
{
  // 
  // Constructor
  // 
  // Parameters:
  //    d detector
  //    r ring 
  //
  fEtaEDists.SetName("EDists");
}
//____________________________________________________________________
AliFMDEnergyFitter::RingHistos::RingHistos(const RingHistos& o)
  : AliForwardUtil::RingHistos(o), 
    fEDist(o.fEDist), 
    fEmpty(o.fEmpty),
    fEtaEDists(), 
    fList(0),
    fFits("AliFMDCorrELossFit::ELossFit"),
    fDebug(0)
{
  // 
  // Copy constructor 
  // 
  // Parameters:
  //    o Object to copy from 
  //
  fFits.Clear();
  TIter next(&o.fEtaEDists);
  TObject* obj = 0;
  while ((obj = next())) fEtaEDists.Add(obj->Clone());
  if (o.fList) {
    fList = new TList;
    fList->SetName(fName);
    TIter next2(o.fList);
    while ((obj = next2())) fList->Add(obj->Clone());
  }
}

//____________________________________________________________________
AliFMDEnergyFitter::RingHistos&
AliFMDEnergyFitter::RingHistos::operator=(const RingHistos& o)
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
  
  if (fEDist) delete fEDist;
  if (fEmpty) delete fEmpty;
  fEtaEDists.Delete();
  if (fList) fList->Delete();
  fFits.Clear();

  fEDist = static_cast<TH1D*>(o.fEDist->Clone());
  fEmpty = static_cast<TH1D*>(o.fEmpty->Clone());
  
  TIter next(&o.fEtaEDists);
  TObject* obj = 0;
  while ((obj = next())) fEtaEDists.Add(obj->Clone());

  if (o.fList) {
    fList = new TList;
    fList->SetName(fName);
    TIter next2(o.fList);
    while ((obj = next2())) fList->Add(obj->Clone());
  }

  return *this;
}
//____________________________________________________________________
AliFMDEnergyFitter::RingHistos::~RingHistos()
{
  // 
  // Destructor 
  //
  if (fEDist) delete fEDist;
  fEtaEDists.Delete();
}

//____________________________________________________________________
void
AliFMDEnergyFitter::RingHistos::Fill(Bool_t empty, Int_t ieta, 
				     Int_t /* icent */,  Double_t mult)
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
  if (empty) { 
    fEmpty->Fill(mult);
    return;
  }
  fEDist->Fill(mult);
  
  // Here, we should first look up in a centrality array, and then in
  // that array look up the eta bin - or perhaps we should do it the
  // other way around?
  if (ieta < 0 || ieta >= fEtaEDists.GetEntries()) return;
  
  TH1D* h = static_cast<TH1D*>(fEtaEDists.At(ieta));
  if (!h) return;

  h->Fill(mult);
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
void
AliFMDEnergyFitter::RingHistos::Make(Int_t    ieta, 
				     Double_t emin, 
				     Double_t emax,
				     Double_t deMax, 
				     Int_t    nDeBins, 
				     Bool_t   incr)
{
  // 
  // Make E/E_mip histogram 
  // 
  // Parameters:
  //    ieta    Eta bin
  //    eMin    Least signal
  //    eMax    Largest signal 
  //    deMax   Maximum energy loss 
  //    nDeBins Number energy loss bins 
  //    incr    Whether to make bins of increasing size
  //
  TH1D* h = 0;
  TString name  = Form(fgkEDistFormat, fName.Data(), ieta);
  TString title = Form("#DeltaE/#DeltaE_{mip} for %s %+5.3f#leq#eta<%+5.3f "
		       "(#eta bin %d)", fName.Data(), emin, emax, ieta);
  if (!incr) 
    h = new TH1D(name.Data(), title.Data(), nDeBins, 0, deMax);
  else { 
    TArrayD deAxis = MakeIncreasingAxis(nDeBins, 0, deMax);
    h = new TH1D(name.Data(), title.Data(), deAxis.fN-1, deAxis.fArray);
  }
    
  h->SetDirectory(0);
  h->SetXTitle("#DeltaE/#DeltaE_{mip}");	
  h->SetFillColor(Color());
  h->SetMarkerColor(Color());
  h->SetLineColor(Color());
  h->SetFillStyle(3001);
  h->Sumw2();
  
  fEtaEDists.AddAt(h, ieta-1);
}
//____________________________________________________________________
TH1D* AliFMDEnergyFitter::RingHistos::MakePar(const char* name, 
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
TH1D*
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
void
AliFMDEnergyFitter::RingHistos::Init(const TAxis& eAxis, 
				     const TAxis& /* cAxis */,
				     Double_t maxDE, Int_t nDEbins, 
				     Bool_t useIncrBin)
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
  fEDist = new TH1D(Form("%s_edist", fName.Data()), 
		    Form("#DeltaE/#DeltaE_{mip} for %s", fName.Data()), 
		    200, 0, 6);
  fEmpty = new TH1D(Form("%s_empty", fName.Data()), 
		    Form("#DeltaE/#DeltaE_{mip} for %s (empty events)", 
			 fName.Data()), 200, 0, 6);
  fList->Add(fEDist);
  fList->Add(fEmpty);
  // Here, we should make an array of centrality bins ...
  // In that case, the structure will be 
  // 
  //    -+- Ring1 -+- Centrality1 -+- Eta1
  //     |         |               +- Eta2
  //     ...       ...             ...
  //     |         |               +- EtaM
  //     |         +- Centrality2 -+- Eta1
  //     ...       ...             ...
  //     |         |               +- EtaM
  //     ...       ...
  //     |         +- CentralityN -+- Eta1
  //     ...       ...             ...
  //     |                         +- EtaM
  //    -+- Ring2 -+- Centrality1 -+- Eta1
  //     |         |               +- Eta2
  //     ...       ...             ...
  //     |         |               +- EtaM
  //     |         +- Centrality2 -+- Eta1
  //     ...       ...             ...
  //     |         |               +- EtaM
  //     ...       ...
  //     |         +- CentralityN -+- Eta1
  //     ...       ...             ...
  //     |                         +- EtaM
  //     ...       ...             ...
  //    -+- RingO -+- Centrality1 -+- Eta1
  //               |               +- Eta2
  //               ...             ...
  //               |               +- EtaM
  //               +- Centrality2 -+- Eta1
  //               ...             ...
  //               |               +- EtaM
  //               ...
  //               +- CentralityN -+- Eta1
  //               ...             ...
  //                               +- EtaM
  // 
  // fEtaEDists.Expand(eAxis.GetNbins());
  for (Int_t i = 1; i <= eAxis.GetNbins(); i++) { 
    Double_t min = eAxis.GetBinLowEdge(i);
    Double_t max = eAxis.GetBinUpEdge(i);
    // Or may we should do it here?  In that case we'd have the
    // following structure:
    //     Ring1 -+- Eta1 -+- Centrality1 
    //            |        +- Centrality2
    //            ...      ...
    //            |        +- CentralityN
    //            +- Eta2 -+- Centrality1 
    //            |        +- Centrality2
    //            ...      ...
    //            |        +- CentralityN
    //            ...
    //            +- EtaM -+- Centrality1 
    //                     +- Centrality2
    //                     ...
    //                     +- CentralityN
    Make(i, min, max, maxDE, nDEbins, useIncrBin);
  }
  fList->Add(&fEtaEDists);
  // fEtaEDists.ls();
}
//____________________________________________________________________
TObjArray*
AliFMDEnergyFitter::RingHistos::Fit(TList*           dir, 
				    const TAxis&     eta,
				    Double_t         lowCut, 
				    UShort_t         nParticles,
				    UShort_t         minEntries,
				    UShort_t         minusBins, 
				    Double_t         relErrorCut, 
				    Double_t         chi2nuCut) const
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

  // Get our ring list from the output container 
  TList* l = GetOutputList(dir);
  if (!l) return 0; 

  // Get the energy distributions from the output container 
  TList* dists = static_cast<TList*>(l->FindObject("EDists"));
  if (!dists) { 
    AliWarning(Form("Didn't find %s_EtaEDists in %s", 
		    fName.Data(), l->GetName()));
    l->ls();
    return 0;
  }

  // Container of the fit results histograms 
  TObjArray* pars  = new TObjArray(3+nParticles+1);
  pars->SetName("FitResults");
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

  
  Int_t nDists = dists->GetEntries();
  Int_t low    = nDists;
  Int_t high   = 0;
  Int_t nEmpty = 0;
  Int_t nLow   = 0;
  Int_t nFitted= 0;
  for (Int_t i = 0; i < nDists; i++) { 
    TH1D* dist = static_cast<TH1D*>(dists->At(i));
    // Ignore empty histograms altoghether 
    if (!dist || dist->GetEntries() <= 0) { 
      nEmpty++;
      continue;
    }

    // Scale to the bin-width
    dist->Scale(1., "width");
    
    // Normalize peak to 1 
    Double_t max = dist->GetMaximum(); 
    if (max <= 0) continue;
    dist->Scale(1/max);
    
    // Check that we have enough entries 
    Int_t nEntries = Int_t(dist->GetEntries());
    if (nEntries <= minEntries) { 
      AliWarning(Form("Histogram at %3d (%s) has too few entries (%d <= %d)",
		      i, dist->GetName(), nEntries, minEntries));
      nLow++;
      continue;
    }

    // Now fit 
    TF1* res = FitHist(dist,lowCut,nParticles,minusBins,
		       relErrorCut,chi2nuCut);
    if (!res) continue;
    nFitted++;
    // dist->GetListOfFunctions()->Add(res);

    // Store eta limits 
    low   = TMath::Min(low,i+1);
    high  = TMath::Max(high,i+1);

    // Get the reduced chi square
    Double_t chi2 = res->GetChisquare();
    Int_t    ndf  = res->GetNDF();
    
    // Store results of best fit in output histograms 
    res->SetLineWidth(4);
    hChi2   ->SetBinContent(i+1, ndf > 0 ? chi2 / ndf : 0);
    hC      ->SetBinContent(i+1, res->GetParameter(kC));   
    hDelta  ->SetBinContent(i+1, res->GetParameter(kDelta)); 
    hXi     ->SetBinContent(i+1, res->GetParameter(kXi));   
    hSigma  ->SetBinContent(i+1, res->GetParameter(kSigma));   
    hSigmaN ->SetBinContent(i+1, res->GetParameter(kSigmaN));   
    hN      ->SetBinContent(i+1, res->GetParameter(kN));   

    hC     ->SetBinError(i+1, res->GetParError(kC));
    hDelta ->SetBinError(i+1, res->GetParError(kDelta));
    hXi    ->SetBinError(i+1, res->GetParError(kXi));
    hSigma ->SetBinError(i+1, res->GetParError(kSigma));
    hSigmaN->SetBinError(i+1, res->GetParError(kSigmaN));
    hN     ->SetBinError(i+1, res->GetParError(kN));

    for (UShort_t j = 0; j < nParticles-1; j++) {
      hA[j]->SetBinContent(i+1, res->GetParameter(kA+j));
      hA[j]->SetBinError(i+1, res->GetParError(kA+j));
    }
  }
  printf("%s: Out of %d histograms, %d where empty, %d had too little data,"
	 "leaving %d to be fitted, of which %d succeeded\n",  
	 GetName(), nDists, nEmpty, nLow, nDists-nEmpty-nLow, nFitted);

  TH1* status = new TH1I("status", "Status of Fits", 5, 0, 5);
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

  // Fit the full-ring histogram 
  TH1* total = GetOutputHist(l, Form("%s_edist", fName.Data()));
  if (total && total->GetEntries() >= minEntries) { 

    // Scale to the bin-width
    total->Scale(1., "width");
    
    // Normalize peak to 1 
    Double_t max = total->GetMaximum(); 
    if (max > 0) total->Scale(1/max);

    TF1* res = FitHist(total,lowCut,nParticles,minusBins,
		       relErrorCut,chi2nuCut);
    if (res) { 
      // Make histograms for the result of this fit 
      Double_t chi2 = res->GetChisquare();
      Int_t    ndf  = res->GetNDF();
      pars->Add(MakeTotal("t_chi2",   "#chi^{2}/#nu", eta, low, high,
			  ndf > 0 ? chi2/ndf : 0, 0));
      pars->Add(MakeTotal("t_c",      "Constant",     eta, low, high,
			  res->GetParameter(kC),
			  res->GetParError(kC)));
      pars->Add(MakeTotal("t_delta",  "#Delta_{p}",   eta, low, high,
			  res->GetParameter(kDelta),
			  res->GetParError(kDelta)));
      pars->Add(MakeTotal("t_xi",     "#xi",          eta, low, high,
			  res->GetParameter(kXi),
			  res->GetParError(kXi)));
      pars->Add(MakeTotal("t_sigma",  "#sigma",       eta, low, high,
			  res->GetParameter(kSigma),
			  res->GetParError(kSigma)));
      pars->Add(MakeTotal("t_sigman", "#sigma_{n}",   eta, low, high,
			  res->GetParameter(kSigmaN),
			  res->GetParError(kSigmaN)));
      pars->Add(MakeTotal("t_n",    "N",              eta, low, high,
			  res->GetParameter(kN),0));
      for (UShort_t j = 0; j < nParticles-1; j++) 
	pars->Add(MakeTotal(Form("t_a%d",j+2), 
			    Form("a_{%d}",j+2), eta, low, high,
			    res->GetParameter(kA+j), 
			    res->GetParError(kA+j)));
    }
  }
    
  // Clean up list of histogram.  Histograms with no entries or 
  // no functions are deleted.  We have to do this using the TObjLink 
  // objects stored in the list since ROOT cannot guaranty the validity 
  // of iterators when removing from a list - tsck.  Should just implement
  // TIter::Remove(). 
  TObjLink* lnk = dists->FirstLink();
  while (lnk) {
    TH1* h = static_cast<TH1*>(lnk->GetObject());
    bool remove = false;
    if (h->GetEntries() <= 0) { 
      // AliWarning(Form("No entries in %s - removing", h->GetName()));
      remove = true;
    }
    else if (h->GetListOfFunctions()->GetEntries() <= 0) {
      // AliWarning(Form("No fuctions associated with %s - removing",
      //            h->GetName()));
      remove = true;
    }
    if (remove) {
      TObjLink* keep = lnk->Next();
      dists->Remove(lnk);
      lnk = keep;
      continue;
    }
    lnk = lnk->Next();
  }

  return pars;
}

//____________________________________________________________________
TF1*
AliFMDEnergyFitter::RingHistos::FitHist(TH1*     dist,
					Double_t lowCut, 
					UShort_t nParticles, 
					UShort_t minusBins, 
					Double_t relErrorCut, 
					Double_t chi2nuCut) const
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
  Double_t maxRange = 10;

  // Create a fitter object 
  AliForwardUtil::ELossFitter f(lowCut, maxRange, minusBins); 
  f.Clear();
  
  
  // If we are only asked to fit a single particle, return this fit, 
  // no matter what. 
  if (nParticles == 1) {
    TF1* r = f.Fit1Particle(dist, 0);
    if (!r) return 0;
    TF1* ret = new TF1(*r);
    dist->GetListOfFunctions()->Add(ret);
    return ret;
  }

  // Fit from 2 upto n particles  
  for (Int_t i = 2; i <= nParticles; i++) f.FitNParticle(dist, i, 0);


  // Now, we need to select the best fit 
  Int_t nFits = f.GetFitResults().GetEntriesFast();
  TF1*  good[nFits];
  for (Int_t i = nFits-1; i >= 0; i--) { 
    good[i] = 0;
    TF1* ff = static_cast<TF1*>(f.GetFunctions().At(i));
    // ff->SetLineColor(Color());
    ff->SetRange(0, maxRange);
    dist->GetListOfFunctions()->Add(new TF1(*ff));
    if (CheckResult(static_cast<TFitResult*>(f.GetFitResults().At(i)),
		    relErrorCut, chi2nuCut)) {
      good[i] = ff;
      ff->SetLineWidth(2);
      // f.fFitResults.At(i)->Print("V");
    }
  }
  // If all else fails, use the 1 particle fit 
  TF1* ret = static_cast<TF1*>(f.GetFunctions().At(0));

  // Find the fit with the most valid particles 
  for (Int_t i = nFits-1; i > 0; i--) {
    if (!good[i]) continue;
    if (fDebug > 1) {
      AliInfo(Form("Choosing fit with n=%d", i+1));
      f.GetFitResults().At(i)->Print();
    }
    ret = good[i];
    break;
  }
  // Give a warning if we're using fall-back 
  if (ret == f.GetFunctions().At(0)) {
    AliWarning("Choosing fall-back 1 particle fit");
  }
  // Copy our result and return (the functions are owned by the fitter)
  TF1* fret = new TF1(*ret);
  return fret;
}

//____________________________________________________________________
Bool_t
AliFMDEnergyFitter::RingHistos::CheckResult(TFitResult* r,
					    Double_t relErrorCut, 
					    Double_t chi2nuCut) const
{
  // 
  // Check the result of the fit. Returns true if 
  // - @f$ \chi^2/\nu < \max{\chi^2/\nu}@f$
  // - @f$ \Delta p_i/p_i < \delta_e@f$ for all parameters.  Note, 
  //   for multi-particle fits, this requirement is relaxed by a 
  //   factor of 2
  // - @f$ a_{n} > 10^{-7}@f$ when fitting to an @f$ n@f$ 
  //   particle response 
  // 
  // Parameters:
  //    r           Result to check
  //    relErrorCut Cut @f$ \delta_e@f$ applied to relative error 
  //                    of parameter.  
  //    chi2nuCut   Cut @f$ \max{\chi^2/\nu}@f$ 
  // 
  // Return:
  //    true if fit is good. 
  //
  if (fDebug > 10) r->Print();
  TString  n    = r->GetName();
  n.ReplaceAll("TFitResult-", "");
  Double_t chi2 = r->Chi2();
  Int_t    ndf  = r->Ndf();
  // Double_t prob = r.Prob();
  Bool_t ret = kTRUE;
  
  // Check that the reduced chi square isn't larger than 5
  if (ndf <= 0 || chi2 / ndf > chi2nuCut) { 
    if (fDebug > 2) {
      AliWarning(Form("%s: chi^2/ndf=%12.5f/%3d=%12.5f>%12.5f", 
		      n.Data(), chi2, ndf, (ndf<0 ? 0 : chi2/ndf),
		      chi2nuCut)); }
    ret = kFALSE;
  }
    
  // Check each parameter 
  UShort_t nPar = r->NPar();
  for (UShort_t i = 0; i < nPar; i++) { 
    if (i == kN) continue;  // Skip the number parameter 
    
    // Get value and error and check value 
    Double_t v  = r->Parameter(i);
    Double_t e  = r->ParError(i);
    if (v == 0) continue;

    // Calculate the relative error and check it 
    Double_t re  = e / v;
    Double_t cut = relErrorCut * (i < kN ? 1 : 2);
    if (re > cut) { 
      if (fDebug > 2) {
	AliWarning(Form("%s: Delta %s/%s=%9.5f/%9.5f=%5.1f%%>%5.1f%%",
			n.Data(), r->ParName(i).c_str(), 
			r->ParName(i).c_str(), e, v, 
			100*(v == 0 ? 0 : e/v),
			100*(cut))); }
      ret = kFALSE;
    }
  }

  // Check if we have scale parameters 
  if (nPar > kN) { 
    
    // Check that the last particle has a significant contribution 
    Double_t lastScale = r->Parameter(nPar-1);
    if (lastScale <= 1e-7) { 
      if (fDebug) {
	AliWarning(Form("%s: %s=%9.6f<1e-7", 
			n.Data(), r->ParName(nPar-1).c_str(), lastScale)); }
      ret = kFALSE;
    }
  }
  return ret;
}


//__________________________________________________________________
void
AliFMDEnergyFitter::RingHistos::FindBestFits(const TList*        d, 
					     AliFMDCorrELossFit& obj,
					     const TAxis&        eta,     
					     Double_t            relErrorCut, 
					     Double_t            chi2nuCut,
					     Double_t            minWeightCut)
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
  TList* dists = static_cast<TList*>(l->FindObject("EDists"));
  if (!dists) { 
    AliWarning(Form("Didn't find %s_EtaEDists in %s", 
		    fName.Data(), l->GetName()));
    l->ls();
    return;
  }
  Int_t nBin = eta.GetNbins();

  for (Int_t b = 1; b <= nBin; b++) { 
    TString n(Form(fgkEDistFormat, fName.Data(), b));
    TH1D*   dist = static_cast<TH1D*>(dists->FindObject(n));
    // Ignore empty histograms altoghether 
    if (!dist || dist->GetEntries() <= 0) continue; 
    
    AliFMDCorrELossFit::ELossFit* best = FindBestFit(dist,
						     relErrorCut,
						     chi2nuCut,  
						     minWeightCut);
    best->fDet  = fDet; 
    best->fRing = fRing;
    best->fBin  = b; // 
    // Double_t eta = fAxis->GetBinCenter(b);
    obj.SetFit(fDet, fRing, b, new AliFMDCorrELossFit::ELossFit(*best));
  }
}

//__________________________________________________________________
AliFMDCorrELossFit::ELossFit* 
AliFMDEnergyFitter::RingHistos::FindBestFit(const TH1* dist,
					    Double_t relErrorCut, 
					    Double_t chi2nuCut,
					    Double_t minWeightCut) 
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
  //    Best fit 
  //
  TList* funcs = dist->GetListOfFunctions();
  TIter  next(funcs);
  TF1*   func = 0;
  fFits.Clear();
  Int_t  i = 0;
  // Info("FindBestFit", "%s", dist->GetName());
  while ((func = static_cast<TF1*>(next()))) { 
    AliFMDCorrELossFit::ELossFit* fit = 
      new(fFits[i++]) AliFMDCorrELossFit::ELossFit(0,*func);
    fit->CalculateQuality(chi2nuCut, relErrorCut, minWeightCut);
  }
  fFits.Sort(false);
  // fFits.Print();
  return static_cast<AliFMDCorrELossFit::ELossFit*>(fFits.At(0));
}


//____________________________________________________________________
void
AliFMDEnergyFitter::RingHistos::Output(TList* dir)
{
  // 
  // Define outputs
  // 
  // Parameters:
  //    dir 
  //
  fList = DefineOutputList(dir);
}

//____________________________________________________________________
//
// EOF
//
