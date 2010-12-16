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
    fLowCut(0.3),
    fNParticles(3),
    fMinEntries(100),
    fFitRangeBinWidth(4),
    fDoFits(false),
    fDoMakeObject(false),
    fEtaAxis(),
    fMaxE(10),
    fNEbins(300), 
    fUseIncreasingBins(true),
    fMaxRelParError(.25),
    fMaxChi2PerNDF(20), 
    fMinWeight(1e-7),
    fDebug(0)
{}

//____________________________________________________________________
AliFMDEnergyFitter::AliFMDEnergyFitter(const char* title)
  : TNamed("fmdEnergyFitter", title), 
    fRingHistos(), 
    fLowCut(0.3),
    fNParticles(3),
    fMinEntries(100),
    fFitRangeBinWidth(4),
    fDoFits(false),
    fDoMakeObject(false),
    fEtaAxis(0,0,0),
    fMaxE(10),
    fNEbins(300), 
    fUseIncreasingBins(true),
    fMaxRelParError(.25),
    fMaxChi2PerNDF(20),  
    fMinWeight(1e-7),
    fDebug(3)
{
  fEtaAxis.SetName("etaAxis");
  fEtaAxis.SetTitle("#eta");
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
    fMaxE(o.fMaxE),
    fNEbins(o.fNEbins), 
    fUseIncreasingBins(o.fUseIncreasingBins),
    fMaxRelParError(o.fMaxRelParError),
    fMaxChi2PerNDF(o.fMaxChi2PerNDF), 
    fMinWeight(o.fMinWeight),
    fDebug(o.fDebug)
{
  TIter    next(&o.fRingHistos);
  TObject* obj = 0;
  while ((obj = next())) fRingHistos.Add(obj);
}

//____________________________________________________________________
AliFMDEnergyFitter::~AliFMDEnergyFitter()
{
  fRingHistos.Delete();
}

//____________________________________________________________________
AliFMDEnergyFitter&
AliFMDEnergyFitter::operator=(const AliFMDEnergyFitter& o)
{
  TNamed::operator=(o);

  fLowCut        = o.fLowCut;
  fNParticles       = o.fNParticles;
  fMinEntries    = o.fMinEntries;
  fFitRangeBinWidth= o.fFitRangeBinWidth;
  fDoFits        = o.fDoFits;
  fDoMakeObject  = o.fDoMakeObject;
  fEtaAxis.Set(o.fEtaAxis.GetNbins(),o.fEtaAxis.GetXmin(),o.fEtaAxis.GetXmax());
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
  if (fEtaAxis.GetNbins() == 0 || 
      fEtaAxis.GetXmin() == fEtaAxis.GetXmax()) 
    SetEtaAxis(eAxis);
  TIter    next(&fRingHistos);
  RingHistos* o = 0;
  while ((o = static_cast<RingHistos*>(next())))
    o->Init(fEtaAxis, fMaxE, fNEbins, fUseIncreasingBins);
}  
//____________________________________________________________________
void
AliFMDEnergyFitter::SetEtaAxis(const TAxis& eAxis)
{
  SetEtaAxis(eAxis.GetNbins(),eAxis.GetXmin(),eAxis.GetXmax());
}
//____________________________________________________________________
void
AliFMDEnergyFitter::SetEtaAxis(Int_t nBins, Double_t etaMin, Double_t etaMax)
{
  fEtaAxis.Set(nBins,etaMin,etaMax);
}

//____________________________________________________________________
Bool_t
AliFMDEnergyFitter::Accumulate(const AliESDFMD& input, 
			       Bool_t           empty)
{
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

	  histos->Fill(empty, ieta-1, mult);

	} // for strip
      } // for sector
    } // for ring 
  } // for detector

  return kTRUE;
}

//____________________________________________________________________
void
AliFMDEnergyFitter::Fit(TList* dir)
{
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
  d->Add(obj, oName.Data());
}

//____________________________________________________________________
void
AliFMDEnergyFitter::DefineOutput(TList* dir)
{
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
  char ind[gROOT->GetDirLevel()+1];
  for (Int_t i = 0; i < gROOT->GetDirLevel(); i++) ind[i] = ' ';
  ind[gROOT->GetDirLevel()] = '\0';
  std::cout << ind << "AliFMDEnergyFitter: " << GetName() << '\n'
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
{}

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
  if (fEDist) delete fEDist;
  fEtaEDists.Delete();
}

//____________________________________________________________________
void
AliFMDEnergyFitter::RingHistos::Fill(Bool_t empty, Int_t ieta, Double_t mult)
{
  if (empty) { 
    fEmpty->Fill(mult);
    return;
  }
  fEDist->Fill(mult);
  
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
  // Service function to define a logarithmic axis. 
  // Parameters: 
  //   n    Number of bins 
  //   min  Minimum of axis 
  //   max  Maximum of axis 
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
				     Double_t maxDE, Int_t nDEbins, 
				     Bool_t useIncrBin)
{
  fEDist = new TH1D(Form("%s_edist", fName.Data()), 
		    Form("#DeltaE/#DeltaE_{mip} for %s", fName.Data()), 
		    200, 0, 6);
  fEmpty = new TH1D(Form("%s_empty", fName.Data()), 
		    Form("#DeltaE/#DeltaE_{mip} for %s (empty events)", 
			 fName.Data()), 200, 0, 6);
  fList->Add(fEDist);
  fList->Add(fEmpty);
  // fEtaEDists.Expand(eAxis.GetNbins());
  for (Int_t i = 1; i <= eAxis.GetNbins(); i++) { 
    Double_t min = eAxis.GetBinLowEdge(i);
    Double_t max = eAxis.GetBinUpEdge(i);
    Make(i, min, max, maxDE, nDEbins, useIncrBin);
  }
  fList->Add(&fEtaEDists);
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
  for (Int_t i = 0; i < nDists; i++) { 
    TH1D* dist = static_cast<TH1D*>(dists->At(i));
    // Ignore empty histograms altoghether 
    if (!dist || dist->GetEntries() <= 0) continue; 

    // Scale to the bin-width
    dist->Scale(1., "width");
    
    // Normalize peak to 1 
    Double_t max = dist->GetMaximum(); 
    if (max <= 0) continue;
    dist->Scale(1/max);
    
    // Check that we have enough entries 
    if (dist->GetEntries() <= minEntries) { 
      AliWarning(Form("Histogram at %3d (%s) has too few entries (%d <= %d)",
		      i, dist->GetName(), Int_t(dist->GetEntries()), 
		      minEntries));
      continue;
    }

    // Now fit 
    TF1* res = FitHist(dist,lowCut,nParticles,minusBins,
		       relErrorCut,chi2nuCut);
    if (!res) continue;
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
  Double_t maxRange = 10;

  // Create a fitter object 
  AliForwardUtil::ELossFitter f(lowCut, maxRange, minusBins); 
  f.Clear();
  
  
  // If we are only asked to fit a single particle, return this fit, 
  // no matter what. 
  if (nParticles == 1) {
    TF1* r = f.Fit1Particle(dist, 0);
    if (!r) return 0;
    return new TF1(*r);
  }

  // Fit from 2 upto n particles  
  for (Int_t i = 2; i <= nParticles; i++) f.FitNParticle(dist, i, 0);


  // Now, we need to select the best fit 
  Int_t nFits = f.fFitResults.GetEntriesFast();
  TF1*  good[nFits];
  for (Int_t i = nFits-1; i >= 0; i--) { 
    good[i] = 0;
    TF1* ff = static_cast<TF1*>(f.fFunctions.At(i));
    // ff->SetLineColor(Color());
    ff->SetRange(0, maxRange);
    dist->GetListOfFunctions()->Add(new TF1(*ff));
    if (CheckResult(static_cast<TFitResult*>(f.fFitResults.At(i)),
		    relErrorCut, chi2nuCut)) {
      good[i] = ff;
      ff->SetLineWidth(2);
      // f.fFitResults.At(i)->Print("V");
    }
  }
  // If all else fails, use the 1 particle fit 
  TF1* ret = static_cast<TF1*>(f.fFunctions.At(0));

  // Find the fit with the most valid particles 
  for (Int_t i = nFits-1; i > 0; i--) {
    if (!good[i]) continue;
    if (fDebug > 1) {
      AliInfo(Form("Choosing fit with n=%d", i+1));
      f.fFitResults.At(i)->Print();
    }
    ret = good[i];
    break;
  }
  // Give a warning if we're using fall-back 
  if (ret == f.fFunctions.At(0)) {
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
AliFMDEnergyFitter::RingHistos::FindBestFits(TList*              d, 
					     AliFMDCorrELossFit& obj,
					     const TAxis&        eta,     
					     Double_t            relErrorCut, 
					     Double_t            chi2nuCut,
					     Double_t            minWeightCut)
{
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
AliFMDEnergyFitter::RingHistos::FindBestFit(TH1* dist,
					    Double_t relErrorCut, 
					    Double_t chi2nuCut,
					    Double_t minWeightCut) 
{
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
  fList = DefineOutputList(dir);
}

//____________________________________________________________________
//
// EOF
//
