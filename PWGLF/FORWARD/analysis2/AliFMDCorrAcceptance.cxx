//
// This class contains the acceptance correction due to dead channels 
// 
//
#include "AliFMDCorrAcceptance.h"
#include <TBrowser.h>
#include <TH2D.h>
#include <AliLog.h>
#include <iostream>

//____________________________________________________________________
AliFMDCorrAcceptance::AliFMDCorrAcceptance()
  : fRingArray(), 
    fCache(0),
    fVertexAxis(0,0,0),
    fHasOverflow(false)
{
  // 
  // Default constructor 
  //
  fRingArray.SetOwner(kTRUE);
  fRingArray.SetName("rings");
  fVertexAxis.SetName("vtxAxis");
  fVertexAxis.SetTitle("v_{z} [cm]");
  
}
//____________________________________________________________________
AliFMDCorrAcceptance::AliFMDCorrAcceptance(const 
					       AliFMDCorrAcceptance& o)
  : TObject(o), 
    fRingArray(o.fRingArray), 
    fCache(o.fCache),
    fVertexAxis(o.fVertexAxis.GetNbins(), o.fVertexAxis.GetXmin(), 
		o.fVertexAxis.GetXmax()),
    fHasOverflow(o.fHasOverflow)
{
  // 
  // Copy constructor 
  // 
  // Parameters:
  //    o Object to copy from 
  //
  fVertexAxis.SetName("vtxAxis");
  fVertexAxis.SetTitle("v_{z} [cm]");
}
//____________________________________________________________________
AliFMDCorrAcceptance::~AliFMDCorrAcceptance()
{
  //
  // Destructor 
  // 
  //
  fRingArray.Clear();
  if (fCache) fCache->Clear();
}
//____________________________________________________________________
AliFMDCorrAcceptance&
AliFMDCorrAcceptance::operator=(const AliFMDCorrAcceptance& o)
{
  // 
  // Assignment operator 
  // 
  // Parameters:
  //    o Object to assign from 
  // 
  // Return:
  //    Reference to this object 
  //
  if (&o == this) return *this;

  fRingArray        = o.fRingArray;
  fCache            = o.fCache;
  fHasOverflow      = o.fHasOverflow;
  SetVertexAxis(o.fVertexAxis);

  return *this;
}
//____________________________________________________________________
TH2D*
AliFMDCorrAcceptance::GetCorrection(UShort_t d, Char_t r, Double_t v) const
{
  // 
  // Get the acceptance correction @f$ a_{r,v}@f$ 
  // 
  // Parameters:
  //    d  Detector number (1-3)
  //    r  Ring identifier (I or O)
  //    v  Primary interaction point @f$z@f$ coordinate
  // 
  // Return:
  //    The correction @f$ a_{r,v}@f$ 
  //
  Int_t b = FindVertexBin(v);
  if (b <= 0) return 0;
  return GetCorrection(d, r, UShort_t(b));
}
//____________________________________________________________________
TH2D*
AliFMDCorrAcceptance::GetCorrection(UShort_t d, Char_t r, UShort_t b) const
{
  // 
  // Get the acceptance correction @f$ a_{r,v}@f$ 
  // 
  // Parameters:
  //    d  Detector number (1-3)
  //    r  Ring identifier (I or O)
  //    b  Bin corresponding to the primary interaction point 
  //           @f$z@f$ coordinate (1 based)
  // 
  // Return:
  //    The correction @f$ a_{r,v}@f$ 
  //
  return static_cast<TH2D*>(GetObject(fRingArray, d, r, b));
}
//____________________________________________________________________
TH1D*
AliFMDCorrAcceptance::GetPhiAcceptance(UShort_t d, Char_t r, Double_t v) const
{
  // 
  // Get the acceptance correction @f$ a_{r,v}@f$ 
  // 
  // Parameters:
  //    d  Detector number (1-3)
  //    r  Ring identifier (I or O)
  //    v  Primary interaction point @f$z@f$ coordinate
  // 
  // Return:
  //    The correction @f$ a_{r,v}@f$ 
  //
  Int_t b = FindVertexBin(v);
  if (b <= 0) return 0;
  return GetPhiAcceptance(d, r, UShort_t(b));
}
//____________________________________________________________________
TH1D*
AliFMDCorrAcceptance::GetPhiAcceptance(UShort_t d, Char_t r, UShort_t b) const
{
  // 
  // Get the acceptance correction @f$ a_{r,v}@f$ 
  // 
  // Parameters:
  //    d  Detector number (1-3)
  //    r  Ring identifier (I or O)
  //    b  Bin corresponding to the primary interaction point 
  //           @f$z@f$ coordinate (1 based)
  // 
  // Return:
  //    The correction @f$ a_{r,v}@f$ 
  //
  if (!fHasOverflow) return 0;
  if (!fCache) FillCache();
  return static_cast<TH1D*>(GetObject(*fCache, d, r, b));
}
  
//____________________________________________________________________
Int_t
AliFMDCorrAcceptance::FindVertexBin(Double_t v) const
{
  // 
  // Find the vertex bin that corresponds to the passed vertex 
  // 
  // Parameters:
  //    vertex The interaction points @f$z@f$-coordinate 
  // 
  // Return:
  //    Vertex bin in @f$[1,N_{\mbox{vertex}}]@f$ or negative if 
  // out of range 
  //
  if (fVertexAxis.GetNbins() <= 0) { 
    AliWarning("No vertex array defined");
    return 0;
  }
  Int_t bin = const_cast<TAxis&>(fVertexAxis).FindBin(v);
  if (bin <= 0 || bin > fVertexAxis.GetNbins()) { 
    AliWarning(Form("vertex %+8.4f out of range [%+8.4f,%+8.4f]",
		    v, fVertexAxis.GetXmin(), fVertexAxis.GetXmax()));
    return 0;
  }
  return bin;
}
//____________________________________________________________________
Int_t
AliFMDCorrAcceptance::GetRingIndex(UShort_t d, Char_t r) const
{
  // Get the index corresponding to the given ring 
  // 
  // Parameters:
  //    d Detector
  //    r Ring 
  // 
  // Return:
  //    Index (0 based) or negative in case of errors
  //
  switch (d) {
  case 1:  return 0;
  case 2:  return (r == 'I' || r == 'i' ? 1 : 2); break;  
  case 3:  return (r == 'I' || r == 'i' ? 3 : 4); break;  

  }
  AliWarning(Form("Index for FMD%d%c not found", d, r));
  return -1;
}
//____________________________________________________________________
TObject*
AliFMDCorrAcceptance::GetObject(const TObjArray& m, UShort_t d, 
				Char_t r, UShort_t b) const
{
  // 
  // Get the object @f$ a_{r,v}@f$ 
  // 
  // Parameters:
  //    m  Mother list
  //    d  Detector number (1-3)
  //    r  Ring identifier (I or O)
  //    b  Bin corresponding to the primary interaction point 
  //           @f$z@f$ coordinate (1 based)
  // 
  // Return:
  //    The correction @f$ a_{r,v}@f$ 
  //
  TObjArray* ringArray = GetRingArray(m, d, r);
  if (!ringArray) return 0;

  if (b <= 0 || b > ringArray->GetEntriesFast()) {
    AliWarning(Form("vertex bin %d out of range [1,%d]", 
		    b, ringArray->GetEntriesFast()));
    return 0;
  }

  TObject* o = ringArray->At(b-1);
  if (o) return o;

  AliWarning(Form("No dead channels map found for FMD%d%c in vertex bin %d",
		    d,r,b));
  return 0;
}
//____________________________________________________________________
TObjArray*
AliFMDCorrAcceptance::GetRingArray(const TObjArray& m, 
				   UShort_t d, Char_t r) const
{
  // 
  // Get the ring array corresponding to the specified ring
  // 
  // Parameters:
  //    d Detector 
  //    r Ring 
  // 
  // Return:
  //    Pointer to ring array, or null in case of problems
  //
  Int_t idx = GetRingIndex(d,r);
  if (idx < 0) return 0;
  
  TObject* o = m.At(idx);
  if (!o) { 
    AliWarning(Form("No array found for FMD%d%c", d, r));
    return 0;
  }

  return static_cast<TObjArray*>(o);
}
//____________________________________________________________________
TObjArray*
AliFMDCorrAcceptance::GetOrMakeRingArray(TObjArray& m, 
					 UShort_t d, Char_t r) const
{
  // 
  // Get the ring array corresponding to the specified ring
  // 
  // Parameters:
  //    d Detector 
  //    r Ring 
  // 
  // Return:
  //    Pointer to ring array, or newly created container 
  //
  Int_t idx = GetRingIndex(d,r);
  if (idx < 0) return 0;
  
  TObject* o = m.At(idx);
  if (!o) { 
    TObjArray* a = new TObjArray(fVertexAxis.GetNbins());
    a->SetName(Form("FMD%d%c", d, r));
    a->SetOwner(kTRUE);
    m.AddAtAndExpand(a, idx);
    return a;
  }

  return static_cast<TObjArray*>(m.At(idx));
}

//____________________________________________________________________
Bool_t
AliFMDCorrAcceptance::SetCorrection(UShort_t d, Char_t r, 
				    UShort_t b, TH2D*  h) 
{
  // 
  // Set the acceptance correction @f$ a_{r,v}(\eta)@f$ 
  // Note, that the object takes ownership of the passed pointer.
  // 
  // Parameters:
  //    d    Detector number (1-3)
  //    r    Ring identifier (I or O)
  //    b    Bin corresponding to the primary interaction point 
  //             @f$z@f$ coordinate  (1 based)
  //    h    @f$ a_{r,v}(\eta)@f$ 
  // 
  // Return:
  //    true if operation succeeded 
  //
  TObjArray* ringArray = GetOrMakeRingArray(fRingArray, d, r);
  if (!ringArray) return false;
  
  if (b <= 0 || b > fVertexAxis.GetNbins()) { 
    AliWarning(Form("Vertex bin %3d out of range [1,%3d]", 
		    b, fVertexAxis.GetNbins()));
    return false;
  }
  h->SetName(Form("FMD%d%c_vtxbin%03d", d, r, b));
  h->SetTitle(Form("Acceptance correction for FMD%d%c "
		   "in vertex bin %d [%+8.4f,%+8.4f]", 
		   d, r, b, fVertexAxis.GetBinLowEdge(b), 
		   fVertexAxis.GetBinUpEdge(b)));
  h->SetXTitle("#eta");
  h->SetYTitle("N_{strips,OK}/N_{strips}");
  h->SetFillStyle(3001);
  h->SetDirectory(0);
  h->SetStats(0);
  ringArray->AddAtAndExpand(h, b-1);
  return kTRUE;
}
//____________________________________________________________________
Bool_t
AliFMDCorrAcceptance::SetCorrection(UShort_t d, Char_t r, 
				    Double_t v, TH2D*  h) 
{
  // 
  // Set the acceptance correction @f$ a_{r,v}(\eta)@f$.
  // Note, that the object takes ownership of the passed pointer.
  // 
  // Parameters:
  //    d    Detector number (1-3)
  //    r    Ring identifier (I or O)
  //    v    Primary interaction point @f$z@f$ coordinate  
  //    h    @f$ a_{r,v}(\eta)@f$ 
  // 
  // Return:
  //    true if operation succeeded 
  //
  Int_t b = FindVertexBin(v);
  if (b <= 0 || b > fVertexAxis.GetNbins()) { 
    AliWarning(Form("Vertex %+8.4f out of range [%+8.4f,%+8.4f]", 
		    v, fVertexAxis.GetXmin(), fVertexAxis.GetXmax()));
    return false;
  }
  return SetCorrection(d, r, UShort_t(b), h);
}

//____________________________________________________________________
void
AliFMDCorrAcceptance::FillCache() const
{
  if (fCache) return;

  fCache = new TObjArray;
  fCache->SetOwner(kTRUE);
  fCache->SetName("cache");

  Int_t nV = fVertexAxis.GetNbins();
  for (UShort_t v = 1; v <= nV; v++) {
    for(UShort_t d = 1; d <= 3;d++) { 
      UShort_t nR = (d == 1 ? 1 : 2);
      for (UShort_t q = 0; q < nR; q++) { 
	Char_t   r  = (q == 0 ? 'I' : 'O');

	TObjArray* a = GetOrMakeRingArray(*fCache, d, r);

	TH2D* corr = GetCorrection(d, r, v);
	if (!corr) continue;

	Int_t nY = corr->GetNbinsY();
	TH1D* h  = corr->ProjectionX("tmp", nY+1, nY+1, "");
	h->SetName(Form("FMD%d%c_vtxbin%03d", d, r, v));
	h->SetTitle(Form("#phi acceptance correction for FMD%d%c "
			   "in vertex bin %d [%+8.4f,%+8.4f]", 
			   d, r, v, fVertexAxis.GetBinLowEdge(v), 
			   fVertexAxis.GetBinUpEdge(v)));
	h->SetXTitle("#eta");
	h->SetYTitle("N_{strips}/N_{strips,OK}");
	h->SetFillStyle(3001);
	h->SetDirectory(0);
	h->SetStats(0);
	a->AddAtAndExpand(h,v-1);

	if (fHasOverflow) continue;

	// Construct the overflow bin from 
	Int_t nX = corr->GetNbinsX();
	for (Int_t eta = 1; eta <= nX; eta++) { 
	  Double_t sum = 0;
	  for (Int_t phi = 1; phi <= nY; phi++) 
	    sum += corr->GetBinContent(eta, phi);
	  if (nY <= 0) continue;
	  h->SetBinContent(eta, nY/sum);
	} // for eta 
      } // for q 
    } // for d 
  } // for v 
}
//____________________________________________________________________
void
AliFMDCorrAcceptance::Browse(TBrowser* b)
{
  // 
  // Browse this object in the browser
  // 
  // Parameters:
  //    b 
  //
  b->Add(&fRingArray);
  if (fCache) b->Add(fCache);
  b->Add(&fVertexAxis);
}
//____________________________________________________________________
void
AliFMDCorrAcceptance::Print(Option_t* option) const
{
  // 
  // Print this object 
  // 
  // Parameters:
  //    option 
  //  
  std::cout << "Acceptance correction due to dead channels" << std::endl;
  fRingArray.Print(option);
  fVertexAxis.Print(option);
}
    
//____________________________________________________________________
//
// EOF
//
