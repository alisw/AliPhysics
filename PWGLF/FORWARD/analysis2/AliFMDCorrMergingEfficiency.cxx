//
// This class contains the secondary correction and the double hit
// correction used in low-flux events.
//
#include "AliFMDCorrMergingEfficiency.h"
#include <TBrowser.h>
#include <TH1D.h>
#include <AliLog.h>
#include <iostream>

//____________________________________________________________________
AliFMDCorrMergingEfficiency::AliFMDCorrMergingEfficiency()
  : fRingArray(), 
    fVertexAxis(0,0,0)
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
AliFMDCorrMergingEfficiency::AliFMDCorrMergingEfficiency(const 
					       AliFMDCorrMergingEfficiency& o)
  : TObject(o), 
    fRingArray(o.fRingArray), 
    fVertexAxis(o.fVertexAxis.GetNbins(), o.fVertexAxis.GetXmin(), 
		o.fVertexAxis.GetXmax())
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
AliFMDCorrMergingEfficiency::~AliFMDCorrMergingEfficiency()
{
  //
  // Destructor 
  // 
  //
  fRingArray.Clear();
}
//____________________________________________________________________
AliFMDCorrMergingEfficiency&
AliFMDCorrMergingEfficiency::operator=(const AliFMDCorrMergingEfficiency& o)
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
  fRingArray        = o.fRingArray;
  SetVertexAxis(o.fVertexAxis);

  return *this;
}
//____________________________________________________________________
TH1D*
AliFMDCorrMergingEfficiency::GetCorrection(UShort_t d, Char_t r, Double_t v) const
{
  // 
  // Get the secondary correction @f$ c_{r,v}@f$ 
  // 
  // Parameters:
  //    d  Detector number (1-3)
  //    r  Ring identifier (I or O)
  //    v  Primary interaction point @f$z@f$ coordinate
  // 
  // Return:
  //    The correction @f$ c_{r,v}@f$ 
  //
  Int_t b = FindVertexBin(v);
  if (b <= 0) return 0;
  return GetCorrection(d, r, UShort_t(b));
}
//____________________________________________________________________
TH1D*
AliFMDCorrMergingEfficiency::GetCorrection(UShort_t d, Char_t r, UShort_t b) const
{
  // 
  // Get the secondary correction @f$ c_{r,v}@f$ 
  // 
  // Parameters:
  //    d  Detector number (1-3)
  //    r  Ring identifier (I or O)
  //    b  Bin corresponding to the primary interaction point 
  //           @f$z@f$ coordinate (1 based)
  // 
  // Return:
  //    The correction @f$ c_{r,v}@f$ 
  //
  TObjArray* ringArray = GetRingArray(d, r);
  if (!ringArray) return 0;

  if (b <= 0 || b > ringArray->GetEntriesFast()) {
    AliWarning(Form("vertex bin %d out of range [1,%d]", 
		    b, ringArray->GetEntriesFast()));
    return 0;
  }

  TObject* o = ringArray->At(b-1);
  if (!o) { 
    AliWarning(Form("No secondary map found for FMD%d%c in vertex bin %d",
		    d,r,b));
    return 0;
  }
  return static_cast<TH1D*>(o);
}
  
//____________________________________________________________________
Int_t
AliFMDCorrMergingEfficiency::FindVertexBin(Double_t v) const
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
AliFMDCorrMergingEfficiency::GetRingIndex(UShort_t d, Char_t r) const
{
  // 
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
TObjArray*
AliFMDCorrMergingEfficiency::GetRingArray(UShort_t d, Char_t r) const
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
  
  TObject* o = fRingArray.At(idx);
  if (!o) { 
    AliWarning(Form("No array found for FMD%d%c", d, r));
    return 0;
  }

  return static_cast<TObjArray*>(o);
}
//____________________________________________________________________
TObjArray*
AliFMDCorrMergingEfficiency::GetOrMakeRingArray(UShort_t d, Char_t r)
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
  
  TObject* o = fRingArray.At(idx);
  if (!o) { 
    TObjArray* a = new TObjArray(fVertexAxis.GetNbins());
    a->SetName(Form("FMD%d%c", d, r));
    a->SetOwner(kTRUE);
    fRingArray.AddAtAndExpand(a, idx);
    return a;
  }

  return static_cast<TObjArray*>(fRingArray.At(idx));
}

//____________________________________________________________________
Bool_t
AliFMDCorrMergingEfficiency::SetCorrection(UShort_t d, Char_t r, 
					   UShort_t b, TH1D*  h) 
{
  // 
  // Set the secondary map correction @f$ m_{r,v}(\eta)@f$ 
  // Note, that the object takes ownership of the passed pointer.
  // 
  // Parameters:
  //    d    Detector number (1-3)
  //    r    Ring identifier (I or O)
  //    b    Bin corresponding to the primary interaction point 
  //             @f$z@f$ coordinate  (1 based)
  //    h    @f$ m_{r,v}(\eta)@f$ 
  // 
  // Return:
  //    true if operation succeeded 
  //
  TObjArray* ringArray = GetOrMakeRingArray(d, r);
  if (!ringArray) return false;
  
  if (b <= 0 || b > fVertexAxis.GetNbins()) { 
    AliWarning(Form("Vertex bin %3d out of range [1,%3d]", 
		    b, fVertexAxis.GetNbins()));
    return false;
  }
  h->SetName(Form("FMD%d%c_vtxbin%03d", d, r, b));
  h->SetTitle(Form("Secondary map correction for FMD%d%c "
		   "in vertex bin %d [%+8.4f,%+8.4f]", 
		   d, r, b, fVertexAxis.GetBinLowEdge(b), 
		   fVertexAxis.GetBinUpEdge(b)));
  h->SetXTitle("#eta");
  h->SetYTitle("dN_{ch}/d#eta / sum_i N_{ch,i}");
  h->SetFillStyle(3001);
  h->SetDirectory(0);
  h->SetStats(0);
  ringArray->AddAtAndExpand(h, b-1);
  return kTRUE;
}
//____________________________________________________________________
Bool_t
AliFMDCorrMergingEfficiency::SetCorrection(UShort_t d, Char_t r, 
					   Double_t v, TH1D*  h) 
{
  // 
  // Set the secondary map correction @f$ m_{r,v}(\eta)@f$.
  // Note, that the object takes ownership of the passed pointer.
  // 
  // Parameters:
  //    d    Detector number (1-3)
  //    r    Ring identifier (I or O)
  //    v    Primary interaction point @f$z@f$ coordinate  
  //    h    @f$ m_{r,v}(\eta)@f$ 
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
AliFMDCorrMergingEfficiency::Browse(TBrowser* b)
{
  // 
  // Browse this object in the browser
  // 
  // Parameters:
  //    b 
  //
  b->Add(&fRingArray);
  b->Add(&fVertexAxis);
}
//____________________________________________________________________
void
AliFMDCorrMergingEfficiency::Print(Option_t* option) const
{
  // 
  // Print this object 
  // 
  // Parameters:
  //    option 
  //  
  std::cout << "Merging efficiency correction" << std::endl;
  fRingArray.Print(option);
  fVertexAxis.Print(option);
}
    
//____________________________________________________________________
//
// EOF
//
