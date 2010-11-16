
#include "AliFMDHistCollector.h"
#include <AliESDFMD.h>
#include <TAxis.h>
#include <TList.h>
#include <TMath.h>
#include "AliFMDAnaParameters.h"
#include "AliLog.h"
#include <TH2D.h>
#include <TH1I.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TArrayI.h>

ClassImp(AliFMDHistCollector)
#if 0
; // For Emacs
#endif 

//____________________________________________________________________
AliFMDHistCollector::AliFMDHistCollector()
  : TNamed(), 
    fList(), 
    fNEvents(0),
    fNCutBins(0),
    fCorrectionCut(0),
    fFirstBins(),
    fLastBins(), 
    fUseEtaFromData(kFALSE),
    fEtaNorm(0), 
    fOutput()
{}

//____________________________________________________________________
AliFMDHistCollector::AliFMDHistCollector(const char* title)
  : TNamed("fmdHistCollector", title), 
    fList(), 
    fNEvents(0),
    fNCutBins(1), 
    fCorrectionCut(0.5),
    fFirstBins(1),
    fLastBins(1), 
    fUseEtaFromData(kFALSE),
    fEtaNorm(0), 
    fOutput()
{
  fList.SetName(GetName());
  fOutput.SetName(GetName());
}

//____________________________________________________________________
AliFMDHistCollector::AliFMDHistCollector(const AliFMDHistCollector& o)
  : TNamed(o), 
    fList(), 
    fNEvents(o.fNEvents),
    fNCutBins(o.fNCutBins), 
    fCorrectionCut(o.fCorrectionCut),
    fFirstBins(o.fFirstBins),
    fLastBins(o.fLastBins), 
    fUseEtaFromData(o.fUseEtaFromData),
    fEtaNorm(o.fEtaNorm), 
    fOutput()
{

  TObject* obj = 0;
  TIter nextL(&o.fList);
  while ((obj = nextL())) fList.Add(obj->Clone());
  TIter nextO(&o.fOutput);
  while ((obj = nextO())) fOutput.Add(obj->Clone());
}

//____________________________________________________________________
AliFMDHistCollector::~AliFMDHistCollector()
{
  fList.Delete();
  fOutput.Delete();
  if (fEtaNorm) delete fEtaNorm;
}

//____________________________________________________________________
AliFMDHistCollector&
AliFMDHistCollector::operator=(const AliFMDHistCollector& o)
{
  SetName(o.GetName());
  SetTitle(o.GetTitle());
  if (fEtaNorm) delete fEtaNorm;
  fList.Delete();
  fOutput.Delete();

  fNEvents        = o.fNEvents;
  fNCutBins       = o.fNCutBins;
  fCorrectionCut  = o.fCorrectionCut;
  fFirstBins      = o.fFirstBins;
  fLastBins       = o.fLastBins;
  fUseEtaFromData = o.fUseEtaFromData;

  TObject* obj = 0;
  TIter nextL(&o.fList);
  while ((obj = nextL())) fList.Add(obj->Clone());
  TIter nextO(&o.fOutput);
  while ((obj = nextO())) fOutput.Add(obj->Clone());
  
  return *this;
}

//____________________________________________________________________
void
AliFMDHistCollector::Init(const TAxis& vtxAxis, const TAxis& etaAxis)
{
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();

  Int_t nVz = vtxAxis.GetNbins();
  fFirstBins.Set(5*nVz);
  fLastBins.Set(5*nVz);

  // Find the eta bin ranges 
  for (Int_t iVz = 0; iVz < nVz; iVz++) {
    AliForwardUtil::Histos* h = new AliForwardUtil::Histos;
    h->Init(etaAxis);
    fList.AddAt(h, iVz);

    // Find the first and last eta bin to use for each ring for 
    // each vertex bin.   This is instead of using the methods 
    // provided by AliFMDAnaParameters 
    for (Int_t iIdx = 0; iIdx < 5; iIdx++) {
      UShort_t d = 0;
      Char_t   r = 0;
      GetDetRing(iIdx, d, r);
      
      // Get the background object 
      TH2F* bg    = pars->GetBackgroundCorrection(d,r,iVz);
      Int_t nEta  = bg->GetNbinsX();
      Int_t first = nEta+1;
      Int_t last  = 0;

      // Loop over the eta bins 
      for (Int_t ie = 1; ie <= nEta; ie++) { 
	if (bg->GetBinContent(ie,1) < fCorrectionCut) continue;

	// Loop over the phi bins to make sure that we 
	// do not have holes in the coverage 
	bool ok = true;
	for (Int_t ip = 1; ip <= bg->GetNbinsY(); ip++) { 
	  if (bg->GetBinContent(ie,ip) < 0.001) {
	    ok = false;
	    continue;
	  }
	}
	if (!ok) continue;

	first = TMath::Min(ie, first);
	last  = TMath::Max(ie, last);
      }
      
      // Store the result for later use 
      fFirstBins[iVz*5+iIdx] = first;
      fLastBins[iVz*5+iIdx]  = last;
    } // for j 
  }

  // Next, using the ranges found above, do the eta acceptance
  // normalisation.
  fEtaNorm = new TH1F("etaAcceptance", "Acceptance in eta", 
		      etaAxis.GetNbins(), 
		      etaAxis.GetXmin(), 
		      etaAxis.GetXmax());
  fEtaNorm->SetDirectory(0);
  fEtaNorm->SetFillColor(kRed+1);
  fEtaNorm->SetFillStyle(3001);
  
  TH1F* tmp = new TH1F("tmp", "tmp", 
		      etaAxis.GetNbins(), 
		      etaAxis.GetXmin(), 
		      etaAxis.GetXmax());
  tmp->SetDirectory(0);
  TH2D* bg = new TH2D("bg", "Background map", 
		      etaAxis.GetNbins(), etaAxis.GetXmin(),etaAxis.GetXmax(),
		      20, 0, 2*TMath::Pi());
  bg->SetDirectory(0);

  Int_t colors[] = { kRed,   kPink, kMagenta, kViolet, kBlue, 
		     kAzure, kCyan, kTeal,    kGreen,  kSpring };


  // Loop over the vertex bins 
  for (Int_t iVz = 0; iVz < nVz; iVz++) { 
    // Clear cache 
    tmp->Reset();
    tmp->SetFillColor(colors[iVz % 10] + iVz/10);
    tmp->SetFillStyle(3001);

    bg->Reset();
    TList* vzList = new TList;
    vzList->SetName(Form("vtxbin%02d", iVz));
    fOutput.AddAt(vzList, iVz);

    // Loop over rings 
    for (Int_t d = 1; d <= 3; d++) { 
      Int_t nq = (d == 1 ? 1 : 2);
      for (Int_t q = 0; q < nq; q++) { 
	Char_t r = (q == 0 ? 'I' : 'O'); 
	
	// Get the first and last bin to use 
	Int_t first, last;
	GetFirstAndLast(d, r, iVz, first, last);
	
	TH2F* ibg = static_cast<TH2F*>(pars->GetBackgroundCorrection(d,r,iVz)
				       ->Clone("bg"));

	// Zero out-side of selected range 
	if (q == 1) { ibg->RebinY(2); ibg->Scale(0.5); }
	for (Int_t iEta = 1; iEta < first; iEta++) 
	  for (Int_t iPhi = 1; iPhi <= 20; iPhi++) 
	    ibg->SetBinContent(iEta,iPhi,0);
	for (Int_t iEta = last+1; iEta <= etaAxis.GetNbins(); iEta++) 
	  for (Int_t iPhi = 1; iPhi <= 20; iPhi++) 
	    ibg->SetBinContent(iEta,iPhi,0);
	// Add to output 
	bg->Add(ibg);
	delete ibg;

	// Loop over the eta bins 
	for (Int_t iEta = first; iEta <= last; iEta++) {
	  Int_t   overlap = GetOverlap(d,r,iEta,iVz);
	  Float_t oc      = tmp->GetBinContent(iEta);
	  Float_t nc      = overlap >= 0? 0.5 : 1;
	  tmp->SetBinContent(iEta, oc + nc);
	}
      }
    }
    vzList->Add(tmp->Clone("etaAcceptance"));
    vzList->Add(bg->Clone("secondaryMap"));
    fEtaNorm->Add(tmp);
  }
  fEtaNorm->Scale(1. / nVz);
  fOutput.Add(fEtaNorm);
  delete tmp;
  delete bg;

}

//____________________________________________________________________
Int_t
AliFMDHistCollector::GetIdx(UShort_t d, Char_t r) const
{
  Int_t idx = -1;
  switch (d) { 
  case 1: idx = 0; break;
  case 2: idx = 1 + (r == 'I' || r == 'i' ? 0 : 1); break;
  case 3: idx = 3 + (r == 'I' || r == 'i' ? 0 : 1); break;
  }
  return idx;
}
//____________________________________________________________________
void
AliFMDHistCollector::GetDetRing(Int_t idx, UShort_t& d, Char_t& r) const
{
  d = 0; 
  r = '\0';
  switch (idx) { 
  case 0: d = 1; r = 'I'; break;
  case 1: d = 2; r = 'I'; break;
  case 2: d = 2; r = 'O'; break;
  case 3: d = 3; r = 'I'; break;
  case 4: d = 3; r = 'O'; break;
  }
}

//____________________________________________________________________
void
AliFMDHistCollector::GetFirstAndLast(Int_t idx, Int_t vtxbin, 
				     Int_t& first, Int_t& last) const
{
  first = 0; 
  last  = 0;
  
  if (idx < 0) return;
  idx += vtxbin * 5;
      
  if (idx < 0 || idx >= fFirstBins.GetSize()) return;
  
  first = fFirstBins.At(idx)+fNCutBins;  
  last  = fLastBins.At(idx)-fNCutBins;
}

//____________________________________________________________________
Int_t
AliFMDHistCollector::GetFirst(Int_t idx, Int_t v) const 
{
  Int_t f, l;
  GetFirstAndLast(idx,v,f,l);
  return f;
}


//____________________________________________________________________
Int_t
AliFMDHistCollector::GetLast(Int_t idx, Int_t v) const 
{
  Int_t f, l;
  GetFirstAndLast(idx,v,f,l);
  return l;
}

//____________________________________________________________________
Int_t 
AliFMDHistCollector::GetOverlap(UShort_t d, Char_t r, 
				Int_t bin, Int_t vtxbin) const
{
  // Get the possibly overlapping histogram 
  Int_t other = -1;
  if (d == 1) {
    if (bin <= GetLast(2,'I',vtxbin)) other = GetIdx(2,'I');
  }
  else if (d == 2 && r == 'I') {
    if      (bin <= GetLast(2,  'O', vtxbin)) other = GetIdx(2, 'O');
    else if (bin >= GetFirst(1, 'I', vtxbin)) other = GetIdx(1, 'I');
  }
  else if (d == 2 && r == 'O') {
    if (bin >= GetFirst(2, 'I', vtxbin)) other = GetIdx(2,'I');
  }
  else if (d == 3 && r == 'O') {
    if (bin <= GetLast(3, 'I', vtxbin)) other = GetIdx(3, 'I');
  }
  else if (d == 3 && r == 'I') {
    if (bin >= GetFirst(3, 'O', vtxbin)) other = GetIdx(3, 'O');
  }
  // AliInfo(Form("FMD%d%c (%d) -> %d", d, r, GetIdx(d,r), other));
  return other;
}
//____________________________________________________________________
Int_t
AliFMDHistCollector::GetOverlap(Int_t idx, Int_t bin, Int_t vtxbin) const
{
  UShort_t d = 0; 
  Char_t   r = '\0';
  GetDetRing(idx, d, r);
  if (d == 0 || r == '\0') return 0;

  return GetOverlap(d, r, bin, vtxbin);
}
  
  
//____________________________________________________________________
Bool_t
AliFMDHistCollector::Collect(AliForwardUtil::Histos& hists,
			     Int_t                   vtxbin, 
			     TH2D&                   out)
{
  if (fList.GetEntries() <= vtxbin) {
    AliWarning(Form("Vertex bin %d out of range [0,%d]", 
		    vtxbin, fList.GetEntries()));
    return kFALSE;
  }
  
  Merge(hists, vtxbin, out);
  Store(hists, vtxbin);

  return kTRUE;
}

//____________________________________________________________________
void
AliFMDHistCollector::Store(AliForwardUtil::Histos& hists,
			   Int_t                   vtxbin)
{
  AliForwardUtil::Histos* histos = 
    static_cast<AliForwardUtil::Histos*>(fList.At(vtxbin));
  if (!histos) { 
    AliWarning(Form("No histogram for vertex bin %d", vtxbin));
    return;
  }
  if (!histos->fFMD1i || 
      !histos->fFMD2i || 
      !histos->fFMD2o || 
      !histos->fFMD3i || 
      !histos->fFMD3o) { 
    AliWarning(Form("Histograms for vertex bin %d not initialised "
		    "(%p,%p,%p,%p,%p)", vtxbin, 
		    histos->fFMD1i, histos->fFMD2i, histos->fFMD2o,
		    histos->fFMD3i, histos->fFMD3o));
    return;
  }
  if (!hists.fFMD1i || 
      !hists.fFMD2i || 
      !hists.fFMD2o || 
      !hists.fFMD3i || 
      !hists.fFMD3o) { 
    AliWarning(Form("Cache histograms not initialised (%p,%p,%p,%p,%p)",
		    hists.fFMD1i, hists.fFMD2i, hists.fFMD2o,
		    hists.fFMD3i, hists.fFMD3o));
    return;
  }

  histos->fFMD1i->Add(hists.fFMD1i);
  histos->fFMD2i->Add(hists.fFMD2i);
  histos->fFMD2o->Add(hists.fFMD2o);
  histos->fFMD3i->Add(hists.fFMD3i);
  histos->fFMD3o->Add(hists.fFMD3o);
}


//____________________________________________________________________
void
AliFMDHistCollector::Merge(AliForwardUtil::Histos& hists,
			   Int_t                   vtxbin, 
			   TH2D&                   out)
{
  // AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  Bool_t isProfile = (out.InheritsFrom(TProfile2D::Class()));

  for (UShort_t d=1; d<=3; d++) { 
    UShort_t nr = (d == 1 ? 1 : 2);
    for (UShort_t q=0; q<nr; q++) { 
      Char_t      r = (q == 0 ? 'I' : 'O');
      TH2D*       h = hists.Get(d,r);
      TH2D*       t = static_cast<TH2D*>(h->Clone(Form("FMD%d%c_tmp",d,r)));
      
      // Outer rings have better phi segmentation - rebin to same as inner. 
      if (q == 1) t->RebinY(2);

      Int_t first = 0;
      Int_t last  = 0;
      GetFirstAndLast(d, r, vtxbin, first, last);

      // Now update profile output 
      for (Int_t iEta = first; iEta <= last; iEta++) { 

	// Get the possibly overlapping histogram 
	Int_t overlap = GetOverlap(d,r,iEta,vtxbin);

	// Fill eta acceptance for this event into the phi underlow bin
	Float_t ooc      = out.GetBinContent(iEta,0);
	Float_t noc      = overlap >= 0? 0.5 : 1;
	out.SetBinContent(iEta, 0, ooc + noc);

	for (Int_t iPhi = 1; iPhi <= h->GetNbinsY(); iPhi++) { 
	  Double_t c = t->GetBinContent(iEta,iPhi);
	  Double_t e = t->GetBinError(iEta,iPhi);

	  // If there's no signal, continue 
	  // if (e <= 0) continue;
	  if (c <= 0 || e <= 0)     continue;
	  
	  if (isProfile) {
	    static_cast<TProfile2D&>(out).Fill(h->GetXaxis()->GetBinCenter(iEta), 
					       h->GetYaxis()->GetBinCenter(iPhi), 
					       c, e);
	    continue;
	  }
	  
	  // If there's no overlapping histogram (ring), then 
	  // fill in data and continue to the next phi bin 
	  if (overlap < 0) { 
	    out.SetBinContent(iEta,iPhi,c);
	    out.SetBinError(iEta,iPhi,e);
	    continue;
	  }

	  // Get the current bin content and error 
	  Double_t oc = out.GetBinContent(iEta,iPhi);
	  Double_t oe = out.GetBinError(iEta,iPhi);

#define USE_STRAIGHT_MEAN
// #define USE_STRAIGHT_MEAN_NONZERO
// #define USE_WEIGHTED_MEAN
// #define USE_MOST_CERTAIN
#if defined(USE_STRAIGHT_MEAN)
	  // calculate the average of old value (half the original), 
	  // and this value, as well as the summed squared errors 
	  // of the existing content (sqrt((e_1/2)^2=sqrt(e_1^2/4)=e_1/2) 
	  // and half the error of this.   
	  //
	  // So, on the first overlapping histogram we get 
	  // 
	  //    c = c_1 / 2
	  //    e = sqrt((e_1 / 2)^2) = e_1/2
	  // 
	  // On the second we get 
	  // 
	  //    c' = c_2 / 2 + c = c_2 / 2 + c_1 / 2 = (c_1+c_2)/2 
	  //    e' = sqrt(e^2 + (e_2/2)^2) 
	  //       = sqrt(e_1^2/4 + e_2^2/4) 
	  //       = sqrt(1/4 * (e_1^2+e_2^2)) 
	  //       = 1/2 * sqrt(e_1^2 + e_2^2)
	  out.SetBinContent(iEta,iPhi,oc + c/2);
	  out.SetBinError(iEta,iPhi,TMath::Sqrt(oe*oe+(e*e)/4));
#elif defined(USE_STRAIGHT_MEAN_NONZERO)
# define ZERO_OTHER
	  // If there's data in the overlapping histogram, 
	  // calculate the average and add the errors in 
	  // quadrature.  
	  // 
	  // If there's no data in the overlapping histogram, 
	  // then just fill in the data 
	  if (oe > 0) {
	    out.SetBinContent(iEta,iPhi,(oc + c)/2);
	    out.SetBinError(iEta,iPhi,TMath::Sqrt(oe*oe + e*e)/2);
	  }
	  else {
	    out.SetBinContent(iEta,iPhi,c);
	    out.SetBinError(iEta,iPhi,e);
	  }	    
#elif defined(USE_WEIGHTED_MEAN) 
	  // Calculate the weighted mean 
	  Double_t w  = 1/(e*e);
	  Double_t sc = w * c;
	  Double_t sw = w;
	  if (oe > 0) {
	    Double_t ow = 1/(oe*oe);
	    sc          += ow * oc;
	    sw          += ow;
	  }
	  Double_t nc = sc / sw;
	  Double_t ne = TMath::Sqrt(1 / sw);
	  out.SetBinContent(iEta,iPhi,nc,ne);
#elif defined(USE_MOST_CERTAIN)
# define ZERO_OTHER
	  if (e < oe) {
	    out.SetBinContent(iEta,iPhi,c);
	    out.SetBinError(iEta,iPhi,e);
	  }
	  else {
	    out.SetBinContent(iEta,iPhi,oc);
	    out.SetBinError(iEta,iPhi,oe);
	  }
#else 
#         error No method for defining content of overlapping bins defined
#endif
#if defined(ZERO_OTHER)
	  // Get the content of the overlapping histogram, 
	  // and zero the content so that we won't use it 
	  // again 
	  UShort_t od; Char_t oR; 
	  GetDetRing(overlap, od, oR);
	  TH2D* other = hists.Get(od,oR);
	  other->SetBinContent(iEta,iPhi,0);
	  other->SetBinError(iEta,iPhi,0);
#endif
	}
      }
      // Remove temporary histogram 
      delete t;
    } // for r
  } // for d 

  // Scale the histogram to the bin size 
  // Int_t    nEta = out.GetNbinsX();
  // Int_t    nPhi = out.GetNbinsY();
  // Double_t rEta = out.GetXaxis()->GetXmax()-out.GetXaxis()->GetXmin();
  // Double_t rPhi = out.GetYaxis()->GetXmax()-out.GetYaxis()->GetXmin();
  // out.Scale(1. / (rEta/nEta*rPhi/nPhi));
  // out.Scale(1., "width"); //
}

//____________________________________________________________________
void
AliFMDHistCollector::ScaleHistograms(const TH1I& nEvents)
{
  fNEvents = &nEvents;

#if 0
  TIter    next(&fList);
  AliForwardUtil::Histos* histos = 0;
  Int_t i = 1;
  while ((histos = static_cast<AliForwardUtil::Histos*>(next()))) {
    Int_t nev = nEvents.GetBinContent(i++);
    if (nev <= 0) continue;
    histos->fFMD1i->Scale(1. / nev);
    histos->fFMD2i->Scale(1. / nev);
    histos->fFMD2o->Scale(1. / nev);
    histos->fFMD3i->Scale(1. / nev);
    histos->fFMD3o->Scale(1. / nev);
  }
#endif 
}

//____________________________________________________________________
void
AliFMDHistCollector::Output(TList* dir)
{
  dir->Add(&fOutput);

  // The axis 
  TAxis* eAxis = static_cast<AliForwardUtil::Histos*>(fList.At(0))
    ->fFMD1i->GetXaxis();
  TAxis* vAxis = fNEvents->GetXaxis();
  Int_t  nVtx  = vAxis->GetNbins();

  // Create profiles of each vertex bin 
  TList  mults;
  for (Int_t i = 0; i < nVtx; i++) {
    TProfile* mult = new TProfile(Form("dndeta_vtx%02d", i), 
				  Form("dN_{ch}/d#eta %f<v_{z}<%f", 
				       fNEvents->GetXaxis()->GetBinLowEdge(i+1),
				       fNEvents->GetXaxis()->GetBinUpEdge(i+1)),
				  eAxis->GetNbins(),eAxis->GetXmin(),
				  eAxis->GetXmax());
    mult->Sumw2();
    mults.AddAt(mult, i);
    fOutput.Add(mult);
  }
  // Loop over vertex bins 
  for (Int_t iVtx = 0; iVtx < nVtx; iVtx++) {
    AliForwardUtil::Histos* histos = 
      static_cast<AliForwardUtil::Histos*>(fList.At(iVtx));
      // Output list for this vertex 
    TList* l = 
      static_cast<TList*>(fOutput.FindObject(Form("vtxbin%02d", iVtx)));

    // Number of events in this vertex bin 
    Int_t nev = fNEvents->GetBinContent(iVtx+1);
    if (nev <= 0) { 
      continue;
    }
    
    // One of the outputs 
    TProfile* mult = static_cast<TProfile*>(mults.At(iVtx));
    if (!mult) { 
      AliWarning(Form("No multiplicity for vertex %d", iVtx));
      continue;
    }

    // Array of histograms and iterator 
    TH2D* hh[] = { histos->fFMD1i, histos->fFMD2i, histos->fFMD2o,
		   histos->fFMD3i, histos->fFMD3o, 0 };
    TH2D** pp = hh;
    Int_t idx = 0;
    while (*pp) { 
      // (eta,phi) histogram of the ring in this vertex 
      (*pp)->SetName(Form("d2ndetadphi_%s", (*pp)->GetName()));
      // Projection (sum) over phi of each vertex histogram 
      TH1D* h = (*pp)->ProjectionX(Form("%s_proj", (*pp)->GetName()));
      // Scale to number of events and bin width 
      h->Scale(1./nev, "width");
      // add to output 
      l->Add(h);

      Int_t firstNonZero = h->GetNbinsX();
      Int_t lastNonZero  = 0;
      if (!fUseEtaFromData) {
	GetFirstAndLast(idx, iVtx, firstNonZero, lastNonZero);
	// firstNonZero += 2;
	// lastNonZero -= 2;
      }
      else {
	// Find first and las non-zerio bins 
	for (Int_t j = 1; j <= h->GetNbinsX(); j++) { 
	  if (h->GetBinContent(j) != 0) { 
	    firstNonZero = TMath::Min(j, firstNonZero);
	    lastNonZero  = TMath::Max(j, lastNonZero);
	  }
	}
	firstNonZero += fNCutBins;
	lastNonZero  -= fNCutBins;
      }

      // Fill profile histogram 
      UShort_t d; Char_t r;
      GetDetRing(idx, d, r);
      for (Int_t iEta = firstNonZero; iEta <= lastNonZero; iEta++) {
	Double_t c = h->GetBinContent(iEta);
	Double_t e = h->GetBinError(iEta);
	if (e <= 0) { 
	  AliWarning(Form("error=%f in bin %d of %s/vtx=%d", 
			  e, iEta, h->GetName(), iVtx));
	  continue;
	}
	mult->Fill(h->GetBinCenter(iEta), c, e);
      }
      idx++;
      pp++;
    } // while (*pp)

    // Loop again to output 2D hists
    pp = hh;
    while (*pp) { 
      // (*pp)->Scale(1./nev, "width");
      l->Add(*pp);
      pp++;
    }
  } // for histos 

  // Result 
  TProfile* total = new TProfile("dndeta", "#frac{1}{N}#frac{dN_{ch}}{d#eta}", 
				 eAxis->GetNbins(),eAxis->GetXmin(),
				 eAxis->GetXmax());

  
  for (Int_t iVtx = 0; iVtx < vAxis->GetNbins(); iVtx++) { 
    TProfile* mult = static_cast<TProfile*>(mults.At(iVtx));
    total->Add(mult);
  }
  fOutput.Add(total);
}

//____________________________________________________________________
//
// EOF
//
	  


