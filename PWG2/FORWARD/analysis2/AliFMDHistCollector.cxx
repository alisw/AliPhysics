// 
// This class collects the event histograms into single histograms, 
// one for each ring in each vertex bin.  
//
// Input:
//   - AliESDFMD object possibly corrected for sharing
//
// Output:
//   - 5 RingHistos objects - each with a number of vertex dependent 
//     2D histograms of the inclusive charge particle density 
// 
// HistCollector used: 
//   - AliFMDCorrSecondaryMap
//
#include "AliFMDHistCollector.h"
#include <AliESDFMD.h>
#include <TAxis.h>
#include <TList.h>
#include <TMath.h>
#include "AliForwardCorrectionManager.h"
#include "AliLog.h"
#include <TH2D.h>
#include <TH1I.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TArrayI.h>
#include <TROOT.h>
#include <iostream>
#include <iomanip>

ClassImp(AliFMDHistCollector)
#if 0
; // For Emacs
#endif 

//____________________________________________________________________
AliFMDHistCollector::AliFMDHistCollector()
  : fNCutBins(0), 
    fCorrectionCut(0), 
    fFirstBins(), 
    fLastBins(), 
    fDebug(0),
    fList(0),
    fSumRings(0),
    fCoverage(0),
    fMergeMethod(kStraightMean),
    fFiducialMethod(kByCut)
{}

//____________________________________________________________________
AliFMDHistCollector::AliFMDHistCollector(const char* title)
  : TNamed("fmdHistCollector", title), 
    fNCutBins(2), 
    fCorrectionCut(0.5), 
    fFirstBins(1), 
    fLastBins(1), 
    fDebug(0),
    fList(0),
    fSumRings(0),
    fCoverage(0),
    fMergeMethod(kStraightMean),
    fFiducialMethod(kByCut)
{
}
//____________________________________________________________________
AliFMDHistCollector::AliFMDHistCollector(const AliFMDHistCollector& o)
  : TNamed(o), 
    fNCutBins(o.fNCutBins), 
    fCorrectionCut(o.fCorrectionCut),
    fFirstBins(o.fFirstBins), 
    fLastBins(o.fLastBins), 
    fDebug(o.fDebug),
    fList(o.fList),
    fSumRings(o.fSumRings),
    fCoverage(o.fCoverage),
    fMergeMethod(o.fMergeMethod),
    fFiducialMethod(o.fFiducialMethod)
{}

//____________________________________________________________________
AliFMDHistCollector::~AliFMDHistCollector()
{ 
  if (fList) delete fList;
}
//____________________________________________________________________
AliFMDHistCollector&
AliFMDHistCollector::operator=(const AliFMDHistCollector& o)
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
  if (&o == this) return *this; 
  TNamed::operator=(o);

  fNCutBins       = o.fNCutBins;
  fCorrectionCut  = o.fCorrectionCut;
  fFirstBins      = o.fFirstBins;
  fLastBins       = o.fLastBins;
  fDebug          = o.fDebug;
  fList           = o.fList;
  fSumRings       = o.fSumRings;
  fCoverage       = o.fCoverage;
  fMergeMethod    = o.fMergeMethod;
  fFiducialMethod = o.fFiducialMethod;

  return *this;
}

//____________________________________________________________________
void
AliFMDHistCollector::Init(const TAxis& vtxAxis,
			  const TAxis& etaAxis)
{
  // 
  // Intialise 
  // 
  // Parameters:
  //    vtxAxis  Vertex axis 
  //  

  AliForwardCorrectionManager& fcm = AliForwardCorrectionManager::Instance();

  fSumRings = new TH2D("sumRings", "Sum in individual rings",
		       etaAxis.GetNbins(), etaAxis.GetXmin(), etaAxis.GetXmax(),
		       5, 1, 6);
  fSumRings->Sumw2();
  fSumRings->SetDirectory(0);
  fSumRings->SetXTitle("#eta");
  fSumRings->GetYaxis()->SetBinLabel(1,"FMD1i");
  fSumRings->GetYaxis()->SetBinLabel(2,"FMD2i");
  fSumRings->GetYaxis()->SetBinLabel(3,"FMD2o");
  fSumRings->GetYaxis()->SetBinLabel(4,"FMD3i");
  fSumRings->GetYaxis()->SetBinLabel(5,"FMD3o");
  fList->Add(fSumRings);

  fCoverage = new TH2D("coverage", "#eta coverage per v_{z}", 
		       etaAxis.GetNbins(),etaAxis.GetXmin(),etaAxis.GetXmax(),
		       vtxAxis.GetNbins(),vtxAxis.GetXmin(),vtxAxis.GetXmax());
  fCoverage->SetDirectory(0);
  fCoverage->SetXTitle("#eta");
  fCoverage->SetYTitle("v_{z} [cm]");
  fCoverage->SetZTitle("n_{bins}");
  fList->Add(fCoverage);
		       
  UShort_t nVz = vtxAxis.GetNbins();
  fFirstBins.Set(5*nVz);
  fLastBins.Set(5*nVz);

  // Find the eta bin ranges 
  for (UShort_t iVz = 1; iVz <= nVz; iVz++) {
    TList*   vtxList = new TList;
    Double_t vMin    = vtxAxis.GetBinLowEdge(iVz);
    Double_t vMax    = vtxAxis.GetBinUpEdge(iVz);
    vtxList->SetName(Form("%c%02d_%c%02d", 
			  vMin < 0 ? 'm' : 'p', int(TMath::Abs(vMin)),
			  vMax < 0 ? 'm' : 'p', int(TMath::Abs(vMax))));
    fList->Add(vtxList);

    // Find the first and last eta bin to use for each ring for 
    // each vertex bin.   This is instead of using the methods 
    // provided by AliFMDAnaParameters 
    for (Int_t iIdx = 0; iIdx < 5; iIdx++) {
      UShort_t d = 0;
      Char_t   r = 0;
      GetDetRing(iIdx, d, r);
      
      // Get the background object 
      // TH2F* bg    = pars->GetBackgroundCorrection(d,r,iVz);
      TH2D* bg    = fcm.GetSecondaryMap()->GetCorrection(d,r,iVz);
      Int_t nEta  = bg->GetNbinsX();
      Int_t first = nEta+1;
      Int_t last  = 0;

      // Loop over the eta bins 
      for (Int_t ie = 1; ie <= nEta; ie++) { 
	// Loop over the phi bins to make sure that we 
	// do not have holes in the coverage 
	bool ok = true;
	for (Int_t ip = 1; ip <= bg->GetNbinsY(); ip++) { 
	  if (!CheckCorrection(bg, ie, ip)) {
	    ok = false;
	    continue;
	  }
	}
	if (!ok) continue;

	first = TMath::Min(ie, first);
	last  = TMath::Max(ie, last);
      }
      
      // Store the result for later use 
      fFirstBins[(iVz-1)*5+iIdx] = first;
      fLastBins[(iVz-1)*5+iIdx]  = last;
      TH2D* obg = static_cast<TH2D*>(bg->Clone(Form("secMapFMD%d%c", d, r)));
      obg->SetDirectory(0);
      obg->Reset();
      vtxList->Add(obg);
      
      TH2D* hitmap = static_cast<TH2D*>(bg->Clone(Form("hitMapFMD%d%c", d, r)));
      if(r == 'O') hitmap->RebinY(2);
      hitmap->SetDirectory(0);
      hitmap->GetZaxis()->SetTitle("");
      hitmap->Reset();
      vtxList->Add(hitmap);

      // Fill diagnostics histograms 
      for (Int_t ie = first+fNCutBins; ie <= last-fNCutBins; ie++) {
	Double_t old = fCoverage->GetBinContent(ie, iVz);
	fCoverage->SetBinContent(ie, iVz, old+1);
	for (Int_t ip = 1; ip <= bg->GetNbinsY(); ip++) {
	  obg->SetBinContent(ie, ip, bg->GetBinContent(ie, ip));
	  obg->SetBinError(ie, ip, bg->GetBinError(ie, ip));
	}
      }
    } // for j 
  }
}

//____________________________________________________________________
Bool_t
AliFMDHistCollector::CheckCorrection(const TH2D* bg, Int_t ie, Int_t ip) const
{
  // 
  // Check if we should include the bin in the data range 
  // 
  // Parameters:
  //    bg Secondary map histogram
  //    ie Eta bin
  //    ip Phi bin
  // 
  // Return:
  //    True if to be used
  //
  Double_t c = bg->GetBinContent(ie,ip);
  switch (fFiducialMethod) { 
  case kByCut:
    return c >= fCorrectionCut;
  case kDistance: 
    if (2 * c < bg->GetBinContent(ie+1,ip) ||
	2 * c < bg->GetBinContent(ie-1,ip)) return false;
    return true;
  default: 
    AliError("No fiducal cut method defined");
  }
  return false;
}
    
//____________________________________________________________________
void
AliFMDHistCollector::DefineOutput(TList* dir)
{
  // 
  // Output diagnostic histograms to directory 
  // 
  // Parameters:
  //    dir List to write in
  //  
  fList = new TList;
  fList->SetOwner();
  fList->SetName(GetName());
  dir->Add(fList);

}


//____________________________________________________________________
Int_t
AliFMDHistCollector::GetIdx(UShort_t d, Char_t r) const
{
  // 
  // Get the ring index from detector number and ring identifier 
  // 
  // Parameters:
  //    d Detector
  //    r Ring identifier 
  // 
  // Return:
  //    ring index or -1 in case of problems 
  //
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
  // 
  // Get the detector and ring from the ring index 
  // 
  // Parameters:
  //    idx Ring index 
  //    d   On return, the detector or 0 in case of errors 
  //    r   On return, the ring id or '0' in case of errors 
  //
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
AliFMDHistCollector::GetFirstAndLast(Int_t idx, UShort_t vtxbin, 
				     Int_t& first, Int_t& last) const
{
  // 
  // Get the first and last eta bin to use for a given ring and vertex 
  // 
  // Parameters:
  //    idx      Ring index as given by GetIdx
  //    vtxBin   Vertex bin (1 based) 
  //    first    On return, the first eta bin to use 
  //    last     On return, the last eta bin to use 
  //
  first = 0; 
  last  = 0;
  
  if (idx    <  0) return;
  if (vtxbin <= 0) return;
  idx += (vtxbin-1) * 5;
      
  if (idx < 0 || idx >= fFirstBins.GetSize()) return;
  
  first = fFirstBins.At(idx)+fNCutBins;  
  last  = fLastBins.At(idx)-fNCutBins;
}

//____________________________________________________________________
Int_t
AliFMDHistCollector::GetFirst(Int_t idx, UShort_t v) const 
{
  // 
  // Get the first eta bin to use for a given ring and vertex 
  // 
  // Parameters:
  //    idx Ring index as given by GetIdx
  //    v vertex bin (1 based)
  // 
  // Return:
  //    First eta bin to use, or -1 in case of problems 
  //  
  Int_t f, l;
  GetFirstAndLast(idx,v,f,l);
  return f;
}


//____________________________________________________________________
Int_t
AliFMDHistCollector::GetLast(Int_t idx, UShort_t v) const 
{
  // 
  // Get the last eta bin to use for a given ring and vertex 
  // 
  // Parameters:
  //    idx Ring index as given by GetIdx
  //    v vertex bin (1 based)
  // 
  // Return:
  //    Last eta bin to use, or -1 in case of problems 
  //  
  Int_t f, l;
  GetFirstAndLast(idx,v,f,l);
  return l;
}

//____________________________________________________________________
Int_t 
AliFMDHistCollector::GetOverlap(UShort_t d, Char_t r, 
				Int_t bin,  UShort_t vtxbin) const
{
  // 
  // Get the possibly overlapping histogram of eta bin @a e in 
  // detector and ring 
  // 
  // Parameters:
  //    d Detector
  //    r Ring 
  //    e Eta bin
  //    v Vertex bin (1 based)
  //
  // Return:
  //    Overlapping histogram index or -1
  //

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
AliFMDHistCollector::GetOverlap(Int_t idx, Int_t bin, UShort_t vtxbin) const
{
  // 
  // Get the possibly overlapping histogram of eta bin @a e in 
  // detector and ring 
  // 
  // Parameters:
  //    i Ring index
  //    e Eta bin
  //    v Vertex bin (1 based)
  //
  // Return:
  //    Overlapping histogram index or -1
  //
  UShort_t d = 0; 
  Char_t   r = '\0';
  GetDetRing(idx, d, r);
  if (d == 0 || r == '\0') return 0;

  return GetOverlap(d, r, bin, vtxbin);
}
  
  
//____________________________________________________________________
void
AliFMDHistCollector::MergeBins(Double_t c,   Double_t e, 
			       Double_t oc,  Double_t oe,
			       Double_t& rc, Double_t& re) const
{
  // 
  // Merge bins accoring to set method
  // 
  // Parameters:
  //    c   Current content
  //    e   Current error
  //    oc  Old content
  //    oe  Old error
  //    rc  On return, the new content
  //    re  On return, tne new error
  //
  rc = re = 0;
  switch (fMergeMethod) { 
  case kStraightMean:
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
    rc = oc + c/2;
    re = TMath::Sqrt(oe*oe+(e*e)/4);
    break;
  case kStraightMeanNoZero:
    // If there's data in the overlapping histogram, 
    // calculate the average and add the errors in 
    // quadrature.  
    // 
    // If there's no data in the overlapping histogram, 
    // then just fill in the data 
    if (oe > 0) {
      rc = (oc + c)/2;
      re = TMath::Sqrt(oe*oe + e*e)/2;
    }
    else {
      rc = c;
      re = e;
    }	    
    break;
  case kWeightedMean: {
    // Calculate the weighted mean 
    Double_t w  = 1/(e*e);
    Double_t sc = w * c;
    Double_t sw = w;
    if (oe > 0) {
      Double_t ow =  1/(oe*oe);
      sc          += ow * oc;
      sw          += ow;
    }
    rc = sc / sw;
    re = TMath::Sqrt(1 / sw);
  }
    break;
  case kLeastError:
    if (e < oe) {
      rc = c;
      re = e;
    }
    else {
      rc = oc;
      re = oe;
    }
    break;
  default:
    AliError("No method for defining content of overlapping bins defined");
    return;
  }
}

//____________________________________________________________________
Bool_t
AliFMDHistCollector::Collect(const AliForwardUtil::Histos& hists,
			     AliForwardUtil::Histos& sums,
			     UShort_t                vtxbin, 
			     TH2D&                   out)
{
  // 
  // Do the calculations 
  // 
  // Parameters:
  //    hists    Cache of histograms 
  //    vtxBin   Vertex bin (1 based)
  //    out      Output histogram
  // 
  // Return:
  //    true on successs 
  //
  AliForwardCorrectionManager& fcm = AliForwardCorrectionManager::Instance();
  const TAxis* vtxAxis = fcm.GetVertexAxis();
  Double_t vMin    = vtxAxis->GetBinLowEdge(vtxbin);
  Double_t vMax    = vtxAxis->GetBinUpEdge(vtxbin);
  TList* vtxList 
    = static_cast<TList*>(fList->FindObject(Form("%c%02d_%c%02d", 
						vMin < 0 ? 'm' : 'p', 
						int(TMath::Abs(vMin)), 
						vMax < 0 ? 'm' : 'p', 
						int(TMath::Abs(vMax)))));
  
  
  for (UShort_t d=1; d<=3; d++) { 
    UShort_t nr = (d == 1 ? 1 : 2);
    for (UShort_t q=0; q<nr; q++) { 
      Char_t      r = (q == 0 ? 'I' : 'O');
      TH2D*       h = hists.Get(d,r);
      TH2D*       t = static_cast<TH2D*>(h->Clone(Form("FMD%d%c_tmp",d,r)));
      Int_t       i = (d == 1 ? 1 : 2*d + (q == 0 ? -2 : -1));
      TH2D*       o = sums.Get(d, r);

      // Get valid range 
      Int_t first = 0;
      Int_t last  = 0;
      GetFirstAndLast(d, r, vtxbin, first, last);
      
      // Zero outside valid range 
      for (Int_t iPhi = 0; iPhi <= t->GetNbinsY()+1; iPhi++) { 
	// Lower range 
	for (Int_t iEta = 1; iEta < first; iEta++) { 
	  t->SetBinContent(iEta,iPhi,0);
	  t->SetBinError(iEta,iPhi,0);
	}
	for (Int_t iEta = last+1; iEta <= t->GetNbinsX(); iEta++) {
	  t->SetBinContent(iEta,iPhi,0);
	  t->SetBinError(iEta,iPhi,0);
	}
      }
      for (Int_t iEta = first; iEta <= last; iEta++)
	t->SetBinContent(iEta,0,1);
      // Add to our per-ring sum 
      o->Add(t);
      
      // Outer rings have better phi segmentation - rebin to same as inner. 
      if (q == 1) t->RebinY(2);

      // Now update profile output 
      for (Int_t iEta = first; iEta <= last; iEta++) { 

	// Get the possibly overlapping histogram 
	Int_t overlap = GetOverlap(d,r,iEta,vtxbin);

	// Fill eta acceptance for this event into the phi underlow bin
	Float_t ooc      = out.GetBinContent(iEta,0);
	Float_t noc      = overlap >= 0? 0.5 : 1;
	out.SetBinContent(iEta, 0, ooc + noc);

	// Should we loop over h or t Y bins - I think it's t
	for (Int_t iPhi = 1; iPhi <= t->GetNbinsY(); iPhi++) { 
	  Double_t c  = t->GetBinContent(iEta,iPhi);
	  Double_t e  = t->GetBinError(iEta,iPhi);
	  Double_t ee = t->GetXaxis()->GetBinCenter(iEta);
	  fSumRings->Fill(ee, i, c);

	  // If there's no signal, continue 
	  // if (e <= 0) continue;
	  if (c <= 0 || e <= 0)     continue;
	  
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

	  Double_t rc, re;
	  MergeBins(c, e, oc, oe, rc, re);
	  out.SetBinContent(iEta,iPhi, rc);
	  out.SetBinError(iEta,iPhi, re);
	}
      }
      // Remove temporary histogram 
      TH2D* hRingSumVtx 
	= static_cast<TH2D*>(vtxList->FindObject(Form("hitMapFMD%d%c", d, r)));
      hRingSumVtx->Add(t);
      delete t;
    } // for r
  } // for d 
  return kTRUE;
}

//____________________________________________________________________
void
AliFMDHistCollector::Print(Option_t* /* option */) const
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
	    << ind << " # of cut bins:          " << fNCutBins << '\n'
	    << ind << " Correction cut:         " << fCorrectionCut << '\n'
	    << ind << " Merge method:           ";
  switch (fMergeMethod) {
  case kStraightMean:       std::cout << "straight mean\n"; break;
  case kStraightMeanNoZero: std::cout << "straight mean (no zeros)\n"; break;
  case kWeightedMean:       std::cout << "weighted mean\n"; break;
  case kLeastError:         std::cout << "least error\n"; break;
  }
    
  std::cout << ind << " Bin ranges:\n" << ind << "  rings  ";
  Int_t nVz = fFirstBins.fN / 5;
  for (Int_t iIdx = 0; iIdx < 5; iIdx++) {
    UShort_t d = 0;
    Char_t   r = 0;
    GetDetRing(iIdx, d, r);
    std::cout << ind << " | FMD" << d << r << " ";
  }
  std::cout << '\n' << ind << "  /vz_bin ";
  for (Int_t iIdx = 0; iIdx < 5; iIdx++) 
    std::cout << "-+--------";
  std::cout << std::endl;

  for (UShort_t iVz = 1; iVz <= nVz; iVz++) {
    std::cout << " " << std::setw(7) << iVz << "   ";
    for (Int_t iIdx = 0; iIdx < 5; iIdx++) {
      UShort_t d = 0;
      Char_t   r = 0;
      GetDetRing(iIdx, d, r);
    
      Int_t first, last;
      GetFirstAndLast(iIdx, iVz, first, last);
      std::cout << " | " << std::setw(3) << first << "-" 
		<< std::setw(3) << last;
    }
    std::cout << std::endl;
  }
}

//____________________________________________________________________
//
// EOF
//
	  


