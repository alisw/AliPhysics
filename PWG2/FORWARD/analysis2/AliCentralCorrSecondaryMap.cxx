//
// This class contains the secondary correction 
// for the central region
//
#include "AliCentralCorrSecondaryMap.h"
#include <TBrowser.h>
#include <TH2D.h>
#include <AliLog.h>
#include <iostream>

//____________________________________________________________________
AliCentralCorrSecondaryMap::AliCentralCorrSecondaryMap()
  : fArray(), 
    fVertexAxis(0,0,0)
{
  // 
  // Default constructor 
  //
  fArray.SetOwner(kTRUE);
  fArray.SetName("rings");
  fVertexAxis.SetName("vtxAxis");
  fVertexAxis.SetTitle("v_{z} [cm]");
  
}
//____________________________________________________________________
AliCentralCorrSecondaryMap::AliCentralCorrSecondaryMap(const 
					       AliCentralCorrSecondaryMap& o)
  : TObject(o), 
    fArray(o.fArray), 
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
AliCentralCorrSecondaryMap::~AliCentralCorrSecondaryMap()
{
  //
  // Destructor 
  // 
  //
  fArray.Clear();
}
//____________________________________________________________________
AliCentralCorrSecondaryMap&
AliCentralCorrSecondaryMap::operator=(const AliCentralCorrSecondaryMap& o)
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
  fArray        = o.fArray;
  SetVertexAxis(o.fVertexAxis);

  return *this;
}
//____________________________________________________________________
TH2D*
AliCentralCorrSecondaryMap::GetCorrection(Double_t v) const
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
  return GetCorrection(UShort_t(b));
}
//____________________________________________________________________
TH2D*
AliCentralCorrSecondaryMap::GetCorrection(UShort_t b) const
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
  
  TObject* o = fArray.At(b-1);
  if (!o) { 
    AliWarning(Form("No dead channels map found for SPD in vertex bin %d",
		    b));
    return 0;
  }
  return static_cast<TH2D*>(o);
}
  
//____________________________________________________________________
Int_t
AliCentralCorrSecondaryMap::FindVertexBin(Double_t v) const
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
Bool_t
AliCentralCorrSecondaryMap::SetCorrection(UShort_t b, TH2D*  h) 
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
  
  if (b <= 0 || b > fVertexAxis.GetNbins()) { 
    AliWarning(Form("Vertex bin %3d out of range [1,%3d]", 
		    b, fVertexAxis.GetNbins()));
    return false;
  }
  h->SetName(Form("SPD_vtxbin%03d", b));
  h->SetTitle(Form("SecondaryMap correction for SPD "
		   "in vertex bin %d [%+8.4f,%+8.4f]", 
		   b, fVertexAxis.GetBinLowEdge(b), 
		   fVertexAxis.GetBinUpEdge(b)));
  h->SetXTitle("#eta");
  h->SetYTitle("dN_{ch}/d#eta / sum_i N_{ch,i}");
  h->SetFillStyle(3001);
  h->SetDirectory(0);
  h->SetStats(0);
  fArray.AddAtAndExpand(h, b-1);
  return kTRUE;
}
//____________________________________________________________________
Bool_t
AliCentralCorrSecondaryMap::SetCorrection(Double_t v, TH2D*  h) 
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
  return SetCorrection( UShort_t(b), h);
}
//____________________________________________________________________
void
AliCentralCorrSecondaryMap::Browse(TBrowser* b)
{
  // 
  // Browse this object in the browser
  // 
  // Parameters:
  //    b 
  //
  b->Add(&fArray);
  b->Add(&fVertexAxis);
}
//____________________________________________________________________
void
AliCentralCorrSecondaryMap::Print(Option_t* option) const
{
  // 
  // Print this object 
  // 
  // Parameters:
  //    option 
  //  
  std::cout << "  SecondaryMap correction" << std::endl;
  std::cout << "   " << std::flush;  
  fArray.Print(option);
  std::cout << "   " << std::flush;
  fVertexAxis.Print(option);
}
    
//____________________________________________________________________
//
// EOF
//
