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
#include "AliFMDCorrSecondaryMap.h"
#include "AliLog.h"
#include <TH2D.h>
#include <TH3D.h>
#include <TH1I.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TObjArray.h>
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
    fDebug(0),
    fList(0),
    fSumRings(0),
    fCoverage(0),
    fSkipped(0),
    fMergeMethod(kStraightMean),
    fFiducialMethod(kByCut),
    fSkipFMDRings(0),
    fBgAndHitMaps(false),
    fVtxList(0), 
    fByCent(0),
    fDoByCent(false)
{
  DGUARD(fDebug, 3, "Default CTOR of AliFMDHistCollector");
}

//____________________________________________________________________
AliFMDHistCollector::AliFMDHistCollector(const char* title)
  : TNamed("fmdHistCollector", title), 
    fNCutBins(2), 
    fCorrectionCut(0.5), 
    fDebug(0),
    fList(0),
    fSumRings(0),
    fCoverage(0),
    fSkipped(0),
    fMergeMethod(kStraightMean),
    fFiducialMethod(kByCut),
    fSkipFMDRings(0),
    fBgAndHitMaps(false),
    fVtxList(0), 
    fByCent(0),
    fDoByCent(false)
{
  DGUARD(fDebug, 3, "Named CTOR of AliFMDHistCollector: %s", title);
}
//____________________________________________________________________
AliFMDHistCollector::AliFMDHistCollector(const AliFMDHistCollector& o)
  : TNamed(o), 
    fNCutBins(o.fNCutBins), 
    fCorrectionCut(o.fCorrectionCut),
    fDebug(o.fDebug),
    fList(o.fList),
    fSumRings(o.fSumRings),
    fCoverage(o.fCoverage),
    fSkipped(o.fSkipped),
    fMergeMethod(o.fMergeMethod),
    fFiducialMethod(o.fFiducialMethod),
    fSkipFMDRings(o.fSkipFMDRings),
    fBgAndHitMaps(o.fBgAndHitMaps),
    fVtxList(o.fVtxList), 
    fByCent(o.fByCent),
    fDoByCent(o.fDoByCent)
{
  DGUARD(fDebug, 3, "Copy CTOR of AliFMDHistCollector");
}

//____________________________________________________________________
AliFMDHistCollector::~AliFMDHistCollector()
{ 
  DGUARD(fDebug, 3, "DTOR of AliFMDHistCollector");
  // if (fList) delete fList;
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
  DGUARD(fDebug, 3, "Assignment of AliFMDHistCollector");
  if (&o == this) return *this; 
  TNamed::operator=(o);

  fNCutBins       = o.fNCutBins;
  fCorrectionCut  = o.fCorrectionCut;
  fDebug          = o.fDebug;
  fList           = o.fList;
  fSumRings       = o.fSumRings;
  fCoverage       = o.fCoverage;
  fSkipped        = o.fSkipped;
  fMergeMethod    = o.fMergeMethod;
  fFiducialMethod = o.fFiducialMethod;
  fSkipFMDRings   = o.fSkipFMDRings;
  fBgAndHitMaps   = o.fBgAndHitMaps;
  fVtxList        = o.fVtxList; 
  fByCent         = o.fByCent;  
  fDoByCent       = o.fDoByCent;
  return *this;
}

//____________________________________________________________________
void
AliFMDHistCollector::SetupForData(const TAxis& vtxAxis,
				  const TAxis& etaAxis)
{
  // 
  // Intialise 
  // 
  // Parameters:
  //    vtxAxis  Vertex axis 
  //  
  DGUARD(fDebug, 1, "Initialization of AliFMDHistCollector");

  // AliForwardCorrectionManager& fcm = AliForwardCorrectionManager::Instance();

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
		       
  fSkipped = new TH1D("skipped", "Rings skipped", 5, 1, 6);
  fSkipped->SetDirectory(0);
  fSkipped->SetFillColor(kRed+1);
  fSkipped->SetFillStyle(3001);
  fSkipped->SetYTitle("Events");
  fSkipped->GetXaxis()->SetBinLabel(1,"FMD1i");
  fSkipped->GetXaxis()->SetBinLabel(2,"FMD2i");
  fSkipped->GetXaxis()->SetBinLabel(3,"FMD2o");
  fSkipped->GetXaxis()->SetBinLabel(4,"FMD3i");
  fSkipped->GetXaxis()->SetBinLabel(5,"FMD3o");
  fList->Add(fSkipped);

  // --- Add parameters to output ------------------------------------
  fList->Add(AliForwardUtil::MakeParameter("nCutBins",fNCutBins));
  fList->Add(AliForwardUtil::MakeParameter("skipRings",fSkipFMDRings));
  fList->Add(AliForwardUtil::MakeParameter("bgAndHits",fBgAndHitMaps));
  fList->Add(AliForwardUtil::MakeParameter("fiducial",Int_t(fFiducialMethod)));
  fList->Add(AliForwardUtil::MakeParameter("fiducialCut",fCorrectionCut));
  fList->Add(AliForwardUtil::MakeParameter("merge",Int_t(fMergeMethod)));

  UShort_t nVz = vtxAxis.GetNbins();
  fVtxList = new TObjArray(nVz, 1);
  fVtxList->SetName("histCollectorVtxBins");
  fVtxList->SetOwner();
  
  // Find the eta bin ranges 
  for (UShort_t iVz = 1; iVz <= nVz; iVz++) {
    Double_t vMin    = vtxAxis.GetBinLowEdge(iVz);
    Double_t vMax    = vtxAxis.GetBinUpEdge(iVz);
    VtxBin*  bin     = new VtxBin(iVz, vMin, vMax, fNCutBins);
    fVtxList->AddAt(bin, iVz);

    bin->SetupForData(fCoverage, fSkipFMDRings, fFiducialMethod, 
		      fCorrectionCut, fList, etaAxis, 
		      fBgAndHitMaps, fBgAndHitMaps);
  }

  if (!fDoByCent) return;

  fByCent = new TList;
  fByCent->SetName("byCentrality");
  fByCent->SetOwner();
  fList->Add(fByCent);

  Int_t    nCent   = 101;
  Double_t minCent = -.5;
  Double_t maxCent = 100.5;
  for (Int_t i = 0; i < 5; i++) { 
    UShort_t d;
    Char_t   r;
    GetDetRing(i, d, r);
    
    TH3* h = new TH3D(Form("FMD%d%c", d, r),
		      Form("dN/d#eta per centrality for FMD%d%c", d, r),
		      etaAxis.GetNbins(), etaAxis.GetXmin(), etaAxis.GetXmax(), 
		      nCent, minCent, maxCent, 1, 0, 1);
    h->SetXTitle("#eta");
    h->SetYTitle("Centrality [%]");
    h->SetZTitle("dN/d#eta");
    h->SetDirectory(0);
    h->SetMarkerColor(AliForwardUtil::RingColor(d, r));
    h->SetMarkerStyle(20);
    fByCent->Add(h);
  }
}
//____________________________________________________________________
Bool_t
AliFMDHistCollector::CheckCorrection(FiducialMethod m, 
				     Double_t       cut, 
				     const TH2D*    bg, 
				     Int_t          ie, 
				     Int_t          ip) 
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
  switch (m) { 
  case kByCut:
    return c >= cut;
  case kDistance: 
    if (2 * c < bg->GetBinContent(ie+1,ip) ||
	2 * c < bg->GetBinContent(ie-1,ip)) return false;
    return true;
  default: 
    AliErrorClass("No fiducal cut method defined");
  }
  return false;
}
    
//____________________________________________________________________
void
AliFMDHistCollector::CreateOutputObjects(TList* dir)
{
  // 
  // Output diagnostic histograms to directory 
  // 
  // Parameters:
  //    dir List to write in
  //  
  DGUARD(fDebug, 1, "Define output of AliFMDHistCollector");
  fList = new TList;
  fList->SetOwner();
  fList->SetName(GetName());
  dir->Add(fList);

}

//____________________________________________________________________
Bool_t
AliFMDHistCollector::CheckSkip(UShort_t d, Char_t r, UShort_t skips) 
{
  UShort_t q  = (r == 'I' || r == 'i' ? 0 : 1);
  UShort_t c = 1 << (d-1);
  UShort_t t = 1 << (c+q-1);

  return (t & skips) == t;
}

//____________________________________________________________________
Int_t
AliFMDHistCollector::GetIdx(UShort_t d, Char_t r)
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
AliFMDHistCollector::GetDetRing(Int_t idx, UShort_t& d, Char_t& r)
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
AliFMDHistCollector::VtxBin*
AliFMDHistCollector::GetVtxBin(Int_t ivtx)
{
  // Parameters:
  //    vtxBin   Vertex bin (1 based) 
  if (!fVtxList) return 0;
  if (ivtx < 1 || ivtx > fVtxList->GetEntriesFast()) return 0;
  VtxBin* bin = static_cast<VtxBin*>(fVtxList->At(ivtx));
  return bin;
}
//____________________________________________________________________
const AliFMDHistCollector::VtxBin*
AliFMDHistCollector::GetVtxBin(Int_t ivtx) const
{
  // Parameters:
  //    vtxBin   Vertex bin (1 based) 
  if (!fVtxList) return 0;
  if (ivtx < 1 || ivtx > fVtxList->GetEntriesFast()) return 0;
  VtxBin* bin = static_cast<VtxBin*>(fVtxList->At(ivtx));
  return bin;
}

//____________________________________________________________________
void
AliFMDHistCollector::MergeBins(MergeMethod   m,
			       Double_t c,   Double_t e, 
			       Double_t oc,  Double_t oe,
			       Double_t& rc, Double_t& re)
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
  switch (m) { 
  case kStraightMean:
  case kPreferInner:  // We only get these two when there's an overlap
  case kPreferOuter:  // between like rings, and we should do a straight mean 
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
  case kSum:
    rc = c + oc;
    re = TMath::Sqrt(oe * oe + e * e);//Add in quadarature 
    break;
  default:
    AliErrorClass("No method for defining content of overlapping bins defined");
    return;
  }
}

//____________________________________________________________________
Bool_t
AliFMDHistCollector::Collect(const AliForwardUtil::Histos& hists,
			     AliForwardUtil::Histos& sums,
			     UShort_t                vtxbin, 
			     TH2D&                   out,
			     Double_t 		     cent,
			     Bool_t                  eta2phi,
			     Bool_t                  add)
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
  DGUARD(fDebug, 1, "Collect final histogram of AliFMDHistCollector");
  // AliForwardCorrectionManager& fcm = AliForwardCorrectionManager::Instance();
  // const TAxis* vtxAxis = fcm.GetVertexAxis();
  // Double_t vMin    = vtxAxis->GetBinLowEdge(vtxbin);
  // Double_t vMax    = vtxAxis->GetBinUpEdge(vtxbin);
  VtxBin*  bin     = GetVtxBin(vtxbin);
  if (!bin) return false;
  Bool_t   ret     = bin->Collect(hists, sums, out, fSumRings, fSkipped, cent, 
				  fMergeMethod, fSkipFMDRings,
				  fByCent, eta2phi, add);

  return ret;
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
AliFMDHistCollector::Print(Option_t* /* option */) const
{
  // 
  // Print information 
  // 
  // Parameters:
  //    option Not used
  //
  TString merge("unknown");
  switch (fMergeMethod) {
  case kStraightMean:       merge = "straight mean"; break;
  case kStraightMeanNoZero: merge = "straight mean (no zeros)"; break;
  case kWeightedMean:       merge = "weighted mean"; break;
  case kLeastError:         merge = "least error"; break;
  case kSum:                merge = "straight sum"; break;
  case kPreferInner:        merge = "prefer inners"; break;
  case kPreferOuter:        merge = "prefer outers"; break;
  }

  AliForwardUtil::PrintTask(*this);
  gROOT->IncreaseDirLevel();
  PFV("# of cut bins",	fNCutBins );
  PFV("Fiducal method", (fFiducialMethod == kByCut ? "cut" : "distance"));
  PFV("Fiducial cut",	fCorrectionCut );
  PFV("Merge method",   merge);
    
  if (!fVtxList) {
    gROOT->DecreaseDirLevel();
    return;
  }
  char ind[64];
  for (Int_t i = 0; i < gROOT->GetDirLevel(); i++) ind[i] = ' ';
  ind[gROOT->GetDirLevel()] = '\0';

  std::cout << ind << "Bin ranges:\n" << ind << "  rings   |   Range  ";
  Int_t nVz = fVtxList->GetEntriesFast();
  for (Int_t iIdx = 0; iIdx < 5; iIdx++) {
    UShort_t d = 0;
    Char_t   r = 0;
    GetDetRing(iIdx, d, r);
    std::cout << ind << " | FMD" << d << r << " ";
  }
  std::cout << '\n' << ind << "  /vz_bin |-----------";
  for (Int_t iIdx = 0; iIdx < 5; iIdx++) 
    std::cout << "-+--------";
  std::cout << std::endl;

  for (UShort_t iVz = 1; iVz <= nVz; iVz++) {
    const VtxBin* bin = GetVtxBin(iVz);
    if (!bin) continue;
    std::cout << "    " << std::right << std::setw(6) << iVz << " | "
	      << std::setw(3) << bin->fLow << " - " << std::left 
	      << std::setw(3) << bin->fHigh << " ";
    for (Int_t iIdx = 0; iIdx < 5; iIdx++) {
      Int_t first, last;
      bin->GetFirstAndLast(iIdx, first, last);
      std::cout << " | " << std::setw(3) << first << "-" 
		<< std::setw(3) << last;
    }
    std::cout << std::endl;
  }
    gROOT->DecreaseDirLevel();
}

//____________________________________________________________________
AliFMDHistCollector::VtxBin::VtxBin(Int_t idx, Double_t minIpZ, Double_t maxIpZ,
				    Int_t nCutBins)
  : fIndex(idx), 
    fLow(minIpZ), 
    fHigh(maxIpZ),
    fHitMap(0), 
    fFirstBin(1), 
    fLastBin(1), 
    fNCutBins(nCutBins)
{
}
//____________________________________________________________________
AliFMDHistCollector::VtxBin::VtxBin(const VtxBin& o)
  : TObject(o), 
    fIndex(o.fIndex), 
    fLow(o.fLow), 
    fHigh(o.fHigh),
    fHitMap(o.fHitMap), 
    fFirstBin(o.fFirstBin), 
    fLastBin(o.fLastBin),
    fNCutBins(o.fNCutBins)
{
}
//____________________________________________________________________
AliFMDHistCollector::VtxBin&
AliFMDHistCollector::VtxBin::operator=(const VtxBin& o)
{
  if (&o == this) return *this;
  fIndex    = o.fIndex;
  fLow      = o.fLow;
  fHigh     = o.fHigh;
  fHitMap   = o.fHitMap;
  fFirstBin = o.fFirstBin;
  fLastBin  = o.fLastBin;
  fNCutBins = o.fNCutBins;
  return *this;
}
//____________________________________________________________________
const Char_t*
AliFMDHistCollector::VtxBin::GetName() const
{
  return Form("%c%03d_%c%03d", 
	      (fLow >= 0 ? 'p' : 'm'), Int_t(TMath::Abs(fLow)), 
	      (fHigh >= 0 ? 'p' : 'm'), Int_t(TMath::Abs(fHigh)));
}
//____________________________________________________________________
void
AliFMDHistCollector::VtxBin::SetupForData(TH2*           coverage,
					  UShort_t       skips,
					  FiducialMethod fiducial, 
					  Double_t       cut,
					  TList*         l,
					  const TAxis&   etaAxis,
					  Bool_t         doHitMaps, 
					  Bool_t         storeSecMap)
{
  TList* out = 0;
  if (doHitMaps || storeSecMap) {
    out = new TList;
    out->SetName(GetName());
    out->SetOwner();
    l->Add(out);
  }
  if (doHitMaps) { 
    fHitMap = new AliForwardUtil::Histos();
    fHitMap->Init(etaAxis);
  }
  fFirstBin.Set(5);
  fLastBin.Set(5);

  AliForwardCorrectionManager& fcm = AliForwardCorrectionManager::Instance();
  
  for (Int_t iIdx = 0; iIdx < 5; iIdx++) {
    UShort_t d = 0;
    Char_t   r = 0;
    GetDetRing(iIdx, d, r);
    
    // Skipping selected FMD rings 
    if (CheckSkip(d, r, skips)) continue;

    // Get the background object 
    TH2D* bg    = fcm.GetSecondaryMap()->GetCorrection(d,r,UShort_t(fIndex));
    Int_t nEta  = bg->GetNbinsX();
    Int_t first = nEta+1;
    Int_t last  = 0;
    
    // Loop over the eta bins 
    for (Int_t ie = 1; ie <= nEta; ie++) { 
      // Loop over the phi bins to make sure that we 
      // do not have holes in the coverage 
      bool ok = true;
      for (Int_t ip = 1; ip <= bg->GetNbinsY(); ip++) { 
	if (!CheckCorrection(fiducial, cut, bg, ie, ip)) {
	  ok = false;
	  continue;
	}
      }
      if (!ok) continue;
      
      first = TMath::Min(ie, first);
      last  = TMath::Max(ie, last);      
    }
    // Store result of first/last bin for this ring
    fFirstBin[iIdx] = first;
    fLastBin[iIdx]  = last;

    if (doHitMaps) { 
      TH2* h = fHitMap->Get(d, r);
      h->SetDirectory(0);
      h->SetName(Form("hitMapFMD%d%c", d, r));
      // if (r == 'O') h->RebinY(2);
      out->Add(h);
    }

    TH2D* obg=0;
    if(storeSecMap) {
      obg = static_cast<TH2D*>(bg->Clone(Form("secMapFMD%d%c", d, r)));
      obg->SetDirectory(0);
      obg->Reset();
      out->Add(obg);
    }
    // Fill diagnostics histograms 
    for (Int_t ie = first+fNCutBins; ie <= last-fNCutBins; ie++) {
      Double_t old = coverage->GetBinContent(ie, fIndex);
      coverage->SetBinContent(ie, fIndex, old+1);
      if(obg) {
	for (Int_t ip = 1; ip <= bg->GetNbinsY(); ip++) {
	  obg->SetBinContent(ie, ip, bg->GetBinContent(ie, ip));
	  obg->SetBinError(ie, ip, bg->GetBinError(ie, ip));
	} // for (ip)
      } // if (doSecHits)
    } // for (ie)
  } // for (iIdx)  
}
  
//____________________________________________________________________
void
AliFMDHistCollector::VtxBin::GetFirstAndLast(Int_t  idx, 
					     Int_t& first, 
					     Int_t& last) const
{
  // Get the first and last eta bin to use for a given ring and vertex 
  // 
  // Parameters:
  //    idx      Ring index as given by GetIdx
  //    first    On return, the first eta bin to use 
  //    last     On return, the last eta bin to use 
  //
  first = 0; 
  last  = 0;
  
  if (idx < 0 || idx >= fFirstBin.GetSize()) return;
  
  first = fFirstBin.At(idx)+fNCutBins;  
  last  = fLastBin.At(idx)-fNCutBins;
}
//____________________________________________________________________
Int_t
AliFMDHistCollector::VtxBin::GetFirst(Int_t  idx) const
{
  Int_t first, last;
  GetFirstAndLast(idx, first , last);
  return first;
}
//____________________________________________________________________
Int_t
AliFMDHistCollector::VtxBin::GetLast(Int_t  idx) const
{
  Int_t first, last;
  GetFirstAndLast(idx, first , last);
  return last;
}

#define PRINT_OVERFLOW(D,R,T,H) do {					\
  printf("Content of FMD%d%c overflow %s rebinning", D, R, T);		\
  Int_t i = 0;								\
  for (Int_t ix = 1; ix <= t->GetNbinsX(); ix++) {			\
  Double_t c = t->GetBinContent(ix, t->GetNbinsY()+1);			\
  if (c <= 1e-9) continue;						\
  if ((i % 10) == 0) printf("\n  ");					\
  printf("%3d: %5.2f ", ix, c);						\
  i++;									\
  }									\
  printf("\n");								\
  } while (false)

//____________________________________________________________________
Bool_t
AliFMDHistCollector::VtxBin::Collect(const AliForwardUtil::Histos& hists, 
				     AliForwardUtil::Histos&       sums, 
				     TH2D&                         out,
				     TH2D*                         sumRings,
				     TH1D*                         skipped,
				     Double_t                      cent,
				     MergeMethod                   m,
				     UShort_t                      skips,
				     TList*                        byCent,
				     Bool_t                        eta2phi,
				     Bool_t                        add)
{
  for (UShort_t d=1; d<=3; d++) { 
    UShort_t nr = (d == 1 ? 1 : 2);
    for (UShort_t q=0; q<nr; q++) { 
      Char_t      r = (q == 0 ? 'I' : 'O');
      Int_t       i = (d == 1 ? 1 : 2*d + (q == 0 ? -2 : -1));
      TH2D*       h = hists.Get(d,r);
      if (CheckSkip(d, r, skips) || !h || 
	  h->TestBit(AliForwardUtil::kSkipRing)) {
	// Skipping a ring - either because disable, not there, or
	// because of flagged (too many outliers, ...)
	skipped->Fill(i);
	continue;
      }
      TH2D*       o = sums.Get(d, r);
      TH2D*       t = static_cast<TH2D*>(h->Clone(Form("FMD%d%c_tmp",d,r)));
      
      // Get valid range 
      Int_t first = 0;
      Int_t last  = 0;
      GetFirstAndLast(d, r, first, last);
      
      // Zero outside valid range 
      Int_t nY = t->GetNbinsY();
      Int_t nX = t->GetNbinsX();
      for (Int_t iPhi = 0; iPhi <= nY+1; iPhi++) { 
	// Lower range 
	for (Int_t iEta = 1; iEta < first; iEta++) { 
	  t->SetBinContent(iEta,iPhi,0);
	  t->SetBinError(iEta,iPhi,0);
	}
	for (Int_t iEta = last+1; iEta <= nX; iEta++) {
	  t->SetBinContent(iEta,iPhi,0);
	  t->SetBinError(iEta,iPhi,0);
	}
      }
      // Fill under-flow bins with eta coverage
      for (Int_t iEta = first; iEta <= last; iEta++) 
	t->SetBinContent(iEta,0,1);
      if (eta2phi) {
	for (Int_t iEta = first; iEta <= last; iEta++) 
	  t->SetBinContent(iEta,nY+1,1);
      }
      
      if (add) {
	// Add to our per-ring sum 
	o->Add(t);

	// If we store hit maps, update here If "eta2phi" is true,
	// then we are here for MC, and so we do not update the hit
	// maps - ever!
	if (fHitMap && !eta2phi) fHitMap->Get(d, r)->Add(t);
      

	if (byCent) { 
	  TH3* dNdetaCent = static_cast<TH3*>(byCent->At(i-1));
	  if (cent >= 0 && dNdetaCent) { 
	    Int_t iCent = dNdetaCent->GetYaxis()->FindBin(cent);
	  
	    if (iCent > 0 && iCent <= dNdetaCent->GetNbinsY()) { 
	      // Make a projection of data 
	      TH1* proj = static_cast<TH1*>(t->ProjectionX("tmp", 1, nY));
	      proj->SetDirectory(0);
	      for (Int_t iEta = 1; iEta <= nX; iEta++) {
		Double_t v1 = proj->GetBinContent(iEta);
		Double_t e1 = proj->GetBinError(iEta);
		Double_t v2 = dNdetaCent->GetBinContent(iEta, iCent, 1);
		Double_t e2 = dNdetaCent->GetBinError(iEta, iCent, 1);
		dNdetaCent->SetBinContent(iEta,iCent,1, v1+v2);
		dNdetaCent->SetBinError(iEta,iCent,1, TMath::Sqrt(e1*e1+e2*e2));
	      
		// Check under-/overflow bins
		Double_t uF = t->GetBinContent(iEta, 0);
		Double_t oF = t->GetBinContent(iEta, nY+1);
		if (uF > 0) {
		  Double_t old = dNdetaCent->GetBinContent(iEta, iCent, 0);
		  dNdetaCent->SetBinContent(iEta, iCent, 0, old + uF);
		}
		if (oF > 0) {
		  Double_t old = dNdetaCent->GetBinContent(iEta, iCent, 2);
		  dNdetaCent->SetBinContent(iEta, iCent, 2, old + oF);
		}
	      } // for(iEta)
	      delete proj;
	    } // if(iCent)
	  } // if (cent)
	} // if (byCent)
      } // if (add)

      // Outer rings have better phi segmentation - rebin to same as inner. 
      if (q == 1) {
	// PRINT_OVERFLOW(d, r, "before", t);
	t->RebinY(2);	
	// PRINT_OVERFLOW(d, r, "after", t);
      }

      nY = t->GetNbinsY();

      // Now update profile output 
      for (Int_t iEta = first; iEta <= last; iEta++) { 

	// Get the possibly overlapping histogram 
	Int_t overlap = GetOverlap(d,r,iEta);

	// Get factor 
	MergeMethod mm  = m; // Possibly override method locally
	Float_t     fac = 1;
	if (m != kSum && overlap >= 0) {
	  // Default is to average 
	  fac = 0.5;
	  if (m == kPreferInner) {
	    // Current one is an outer overlapping an inner 
	    if ((r == 'o' || r == 'O') && 
		(overlap == 0 || overlap == 1 || overlap == 3))
	      // Do not use this signal 
	      fac = 0;
	    // Current one is inner overlapping an outer
	    else if ((r == 'i' || r == 'I') && (overlap == 2 || overlap == 4))
	      // Prefer this one 
	      fac = 1;
	    else 
	      // In case of two overlapping inners 
	      mm = kStraightMean;
	  }
	  else if (m == kPreferOuter) {
	    // Current one is an inner overlapping an outer 
	    if ((r == 'i' || r == 'I') && (overlap == 2 || overlap == 4))
	      // Do not use this signal 
	      fac = 0;
	    else if ((r == 'O' || r == 'o') && 
		     (overlap == 0 || overlap == 1 || overlap == 3))
	      fac = 1;
	    else 
	      // In case of two overlapping outers
	      mm = kStraightMean;
	  }
	}

	// Fill eta acceptance for this event into the phi underflow bin
	Float_t ooc      = out.GetBinContent(iEta,0);
	out.SetBinContent(iEta, 0, ooc + fac);

	// Fill phi acceptance for this event into the phi overflow bin
	Float_t oop      = out.GetBinContent(iEta,nY+1);
	Float_t nop      = t->GetBinContent(iEta,nY+1);
#if 0
	Info("", "etaBin=%3d Setting phi acceptance to %f(%f+%f)=%f", 
	     iEta, fac, oop, nop, fac*(oop+nop));
#endif
	out.SetBinContent(iEta, nY+1, fac * nop + oop);

	// Should we loop over h or t Y bins - I think it's t
	for (Int_t iPhi = 1; iPhi <= nY; iPhi++) { 
	  Double_t c  = t->GetBinContent(iEta,iPhi);
	  Double_t e  = t->GetBinError(iEta,iPhi);
	  Double_t ee = t->GetXaxis()->GetBinCenter(iEta);
	  sumRings->Fill(ee, i, c);

	  // If there's no signal or the signal was ignored because we
	  // prefer the inners/outers, continue if (e <= 0) continue;
	  if (fac <= 0 || c <= 0 || e <= 0)     continue;
	  
	  // If there's no overlapping histogram (ring) or if we
	  // prefer inner/outer, then fill in data and continue to the
	  // next phi bin
	  if (overlap < 0 || fac >= 1) { 
	    out.SetBinContent(iEta,iPhi,c);
	    out.SetBinError(iEta,iPhi,e);
	    continue;
	  }

	  // Get the current bin content and error 
	  Double_t oc = out.GetBinContent(iEta,iPhi);
	  Double_t oe = out.GetBinError(iEta,iPhi);

	  Double_t rc, re;
	  MergeBins(mm, c, e, oc, oe, rc, re);
	  out.SetBinContent(iEta,iPhi, rc);
	  out.SetBinError(iEta,iPhi, re);
	}
      }
      // Remove temporary histogram 
      delete t;
    } // for r
  } // for d 
  return true;
}
//____________________________________________________________________
Int_t 
AliFMDHistCollector::VtxBin::GetOverlap(UShort_t d, Char_t r, 
					Int_t bin) const
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
    if (bin <= GetLast(2,'I')) other = GetIdx(2,'I');
  }
  else if (d == 2 && r == 'I') {
    if      (bin <= GetLast(2,  'O')) other = GetIdx(2, 'O');
    else if (bin >= GetFirst(1, 'I')) other = GetIdx(1, 'I');
  }
  else if (d == 2 && r == 'O') {
    if (bin >= GetFirst(2, 'I'))      other = GetIdx(2,'I');
  }
  else if (d == 3 && r == 'O') {
    if (bin <= GetLast(3, 'I'))       other = GetIdx(3, 'I');
  }
  else if (d == 3 && r == 'I') {
    if (bin >= GetFirst(3, 'O'))      other = GetIdx(3, 'O');
  }
  // AliInfo(Form("FMD%d%c (%d) -> %d", d, r, GetIdx(d,r), other));
  return other;
}
  

//____________________________________________________________________
//
// EOF
//
	  


