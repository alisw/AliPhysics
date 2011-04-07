//
// This class contains the acceptance correction due to dead channels 
// 
//
#include "AliCentralCorrAcceptance.h"
#include <TBrowser.h>
#include <TH1D.h>
#include <AliLog.h>
#include <iostream>

//____________________________________________________________________
AliCentralCorrAcceptance::AliCentralCorrAcceptance()
  : fArray(), 
    fVertexAxis(0,0,0)
{
  // 
  // Default constructor 
  //
  fArray.SetOwner(kTRUE);
  fArray.SetName("acceptance");
  fVertexAxis.SetName("vtxAxis");
  fVertexAxis.SetTitle("v_{z} [cm]");
  
}
//____________________________________________________________________
AliCentralCorrAcceptance::AliCentralCorrAcceptance(const 
					       AliCentralCorrAcceptance& o)
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
AliCentralCorrAcceptance::~AliCentralCorrAcceptance()
{
  //
  // Destructor 
  // 
  //
  fArray.Clear();
}
//____________________________________________________________________
AliCentralCorrAcceptance&
AliCentralCorrAcceptance::operator=(const AliCentralCorrAcceptance& o)
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
TH1D*
AliCentralCorrAcceptance::GetCorrection(Double_t v) const
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
TH1D*
AliCentralCorrAcceptance::GetCorrection(UShort_t b) const
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
  
  // TObjArray* ringArray = GetRingArray(d, r);
  //if (!ringArray) return 0;

  if (b <= 0 || b > fArray.GetEntriesFast()) {
    AliWarning(Form("vertex bin %d out of range [1,%d]", 
		    b, fArray.GetEntriesFast()));
    return 0;
  }

  TObject* o = fArray.At(b-1);
  if (!o) { 
    AliWarning(Form("No dead channels map found for SPD in vertex bin %d",
		    b));
    return 0;
  }
  return static_cast<TH1D*>(o);
}
  
//____________________________________________________________________
Int_t
AliCentralCorrAcceptance::FindVertexBin(Double_t v) const
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
AliCentralCorrAcceptance::SetCorrection(UShort_t b, TH1D*  h) 
{
  // 
  // Set the acceptance correction @f$ a_{r,v}(\eta)@f$ 
  // Note, that the object takes ownership of the passed pointer.
  // 
  // Parameters:
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
  
  h->SetName(Form("SPD_vtxbin%03d",  b));
  h->SetTitle(Form("Acceptance correction for SPD "
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
AliCentralCorrAcceptance::SetCorrection(Double_t v, TH1D*  h) 
{
  // 
  // Set the acceptance correction @f$ a_{r,v}(\eta)@f$.
  // Note, that the object takes ownership of the passed pointer.
  // 
  // Parameters
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
  return SetCorrection(UShort_t(b), h);
}
//____________________________________________________________________
void
AliCentralCorrAcceptance::Browse(TBrowser* b)
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
AliCentralCorrAcceptance::Print(Option_t* option) const
{
  // 
  // Print this object 
  // 
  // Parameters:
  //    option 
  //  
  std::cout << "  Acceptance correction due to dead channels" << std::endl;
  std::cout << "   " << std::flush;  
  fArray.Print(option);
  std::cout << "   " << std::flush;  
  fVertexAxis.Print(option);
}
    
//____________________________________________________________________
//
// EOF
//
