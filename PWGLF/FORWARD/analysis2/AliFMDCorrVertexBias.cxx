//
// This class contains the secondary correction and the double hit
// correction used in low-flux events.
//
#include "AliFMDCorrVertexBias.h"
#include <TBrowser.h>
#include <TH2D.h>
#include <AliLog.h>
#include <iostream>

//____________________________________________________________________
AliFMDCorrVertexBias::AliFMDCorrVertexBias()
  : fVertexArray(), 
    fVertexAxis(0,0,0)
{
  // 
  // Default constructor 
  //
  fVertexArray.SetOwner(kTRUE);
  fVertexArray.SetName("rings");
  fVertexAxis.SetName("vtxAxis");
  fVertexAxis.SetTitle("v_{z} [cm]");
  
}
//____________________________________________________________________
AliFMDCorrVertexBias::AliFMDCorrVertexBias(const AliFMDCorrVertexBias& o)
  : TObject(o), 
    fVertexArray(o.fVertexArray), 
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
AliFMDCorrVertexBias::~AliFMDCorrVertexBias()
{
  //
  // Destructor 
  // 
  //
  fVertexArray.Clear();
}
//____________________________________________________________________
AliFMDCorrVertexBias&
AliFMDCorrVertexBias::operator=(const AliFMDCorrVertexBias& o)
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
  fVertexArray        = o.fVertexArray;
  SetVertexAxis(o.fVertexAxis);

  return *this;
}
//____________________________________________________________________
TH2D*
AliFMDCorrVertexBias::GetCorrection(Char_t r, Double_t v) const
{
  // 
  // Get the vertex bias correction @f$ b_{v}@f$ 
  // 
  // Parameters:
  //    r  Ring type (I or O)
  //    v  Primary interaction point @f$z@f$ coordinate
  // 
  // Return:
  //    The correction @f$ b_{v}@f$ 
  //
  Int_t b = FindVertexBin(v);
  if (b <= 0) return 0;
  return GetCorrection(r, UShort_t(b));
}
//____________________________________________________________________
TH2D*
AliFMDCorrVertexBias::GetCorrection(Char_t r, UShort_t b) const
{
  // 
  // Get the vertex bias correction @f$ b_{v}@f$ 
  // 
  // Parameters:
  //    r  Ring type (I or O)
  //    b  Bin corresponding to the primary interaction point 
  //           @f$z@f$ coordinate (1 based)
  // 
  // Return:
  //    The correction @f$ b_{v}@f$ 
  //
  TObjArray* vertexarray = GetVertexArray(b);
  if (!vertexarray) return 0;

  Int_t idx = -1;
  switch (r) { 
  case 'i': case 'I': idx = 0; break;
  case 'o': case 'O': idx = 1; break;
  }
  if (idx < 0) {
    AliWarning(Form("Unknown ting type %c, not one of [iIoO]", r));
    return 0;
  }

  TObject* o = vertexarray->At(idx);
  if (!o) { 
    AliWarning(Form("No vertex bias found for ring type %c in vertex bin %d",
		    r,b));
    return 0;
  }
  return static_cast<TH2D*>(o);
}
  
//____________________________________________________________________
Int_t
AliFMDCorrVertexBias::FindVertexBin(Double_t v) const
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
TObjArray*
AliFMDCorrVertexBias::GetVertexArray(UShort_t v) const
{
  // 
  // Get the vertex array corresponding to the specified ring
  // 
  // Parameters:
  //    v vertex bin (1 based)
  // 
  // Return:
  //    Pointer to vertex array, or null in case of problems
  //
  if (v <= 0 || v > fVertexAxis.GetNbins()) {
    AliWarning(Form("vertex bin %d out of range [1,%d]",
		    v, fVertexAxis.GetNbins()));
    return 0;
  }
  
  TObject* o = fVertexArray.At(v-1);
  if (!o) { 
    AliWarning(Form("No array found for vertex bin %d", v));
    return 0;
  }

  return static_cast<TObjArray*>(o);
}
//____________________________________________________________________
TObjArray*
AliFMDCorrVertexBias::GetOrMakeVertexArray(UShort_t v)
{
  // 
  // Get the vertex array corresponding to the specified ring
  // 
  // Parameters:
  //    v vertex bin (1 based)
  // 
  // Return:
  //    Pointer to vertex array, or newly created container 
  //
  if (v <= 0 || v > fVertexAxis.GetNbins()) {
    AliWarning(Form("vertex bin %d out of range [1,%d]",
		    v, fVertexAxis.GetNbins()));
    return 0;
  }
    
  TObject* o = fVertexArray.At(v-1);
  if (!o) { 
    TObjArray* a = new TObjArray(fVertexAxis.GetNbins());
    a->SetName(Form("vertexbin%02d", v));
    a->SetOwner(kTRUE);
    fVertexArray.AddAtAndExpand(a, v-1);
    return a;
  }

  return static_cast<TObjArray*>(fVertexArray.At(v-1));
}

//____________________________________________________________________
Bool_t
AliFMDCorrVertexBias::SetCorrection(Char_t r, UShort_t b, TH2D*  h) 
{
  // 
  // Set the vertex bias correction @f$ b_{v}(\eta,\varphi)@f$ 
  // Note, that the object takes ownership of the passed pointer.
  // 
  // Parameters:
  //    r    Ring type (I or O)
  //    b    Bin corresponding to the primary interaction point 
  //             @f$z@f$ coordinate  (1 based)
  //    h    @f$ b_{v}(\eta,\varphi)@f$ 
  // 
  // Return:
  //    true if operation succeeded 
  //
  TObjArray* vertexarray = GetOrMakeVertexArray(b);
  if (!vertexarray) return false;
  

  Int_t idx = -1;
  switch (r) { 
  case 'i': case 'I': idx = 0; break;
  case 'o': case 'O': idx = 1; break;
  }
  if (idx < 0) {
    AliWarning(Form("Unknown ting type %c, not one of [iIoO]", r));
    return false;
  }
  h->SetName(Form("FMDX%c", r));
  h->SetTitle(Form("Vertex bias correction for %c rings "
		   "in vertex bin %d [%+8.4f,%+8.4f]", 
		   r, b, fVertexAxis.GetBinLowEdge(b), 
		   fVertexAxis.GetBinUpEdge(b)));
  h->SetXTitle("#eta");
  h->SetYTitle("#phi [radians]");
  h->SetZTitle("1/N_{t}#sum_{i}^{N_{tv}} N_{ch,i,primary} / "
	       "1/N_{v}#sum_{i}^{N_{v}} N_{ch,i,primary}");
  h->SetDirectory(0);
  h->SetStats(0);

  vertexarray->AddAtAndExpand(h, idx);
  return kTRUE;
}
//____________________________________________________________________
Bool_t
AliFMDCorrVertexBias::SetCorrection(Char_t r, Double_t v, TH2D*  h) 
{
  // 
  // Set the vertex bias correction @f$ b_{v}(\eta,\varphi)@f$.
  // Note, that the object takes ownership of the passed pointer.
  // 
  // Parameters:
  //    r    Ring type (I or O)
  //    v    Primary interaction point @f$z@f$ coordinate  
  //    h    @f$ b_{v}(\eta,\varphi)@f$ 
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
  return SetCorrection(r, UShort_t(b), h);
}
//____________________________________________________________________
void
AliFMDCorrVertexBias::Browse(TBrowser* b)
{
  // 
  // Browse this object in the browser
  // 
  // Parameters:
  //    b 
  //
  b->Add(&fVertexArray);
  b->Add(&fVertexAxis);
}
//____________________________________________________________________
void
AliFMDCorrVertexBias::Print(Option_t* option) const
{
  // 
  // Print this object 
  // 
  // Parameters:
  //    option 
  //  
  std::cout << "Vertex bias correction" << std::endl;
  fVertexArray.Print(option);
  fVertexAxis.Print(option);
}
    
//____________________________________________________________________
//
// EOF
//
