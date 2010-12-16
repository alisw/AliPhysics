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
  fVertexAxis.SetName("vtxAxis");
  fVertexAxis.SetTitle("v_{z} [cm]");
}
//____________________________________________________________________
AliFMDCorrMergingEfficiency::~AliFMDCorrMergingEfficiency()
{
  fRingArray.Clear();
}
//____________________________________________________________________
AliFMDCorrMergingEfficiency&
AliFMDCorrMergingEfficiency::operator=(const AliFMDCorrMergingEfficiency& o)
{
  fRingArray        = o.fRingArray;
  SetVertexAxis(o.fVertexAxis);

  return *this;
}
//____________________________________________________________________
TH1D*
AliFMDCorrMergingEfficiency::GetCorrection(UShort_t d, Char_t r, Double_t v) const
{
  Int_t b = FindVertexBin(v);
  if (b <= 0) return 0;
  return GetCorrection(d, r, UShort_t(b));
}
//____________________________________________________________________
TH1D*
AliFMDCorrMergingEfficiency::GetCorrection(UShort_t d, Char_t r, UShort_t b) const
{
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
  b->Add(&fRingArray);
  b->Add(&fVertexAxis);
}
//____________________________________________________________________
void
AliFMDCorrMergingEfficiency::Print(Option_t* option) const
{
  std::cout << "Merging efficiency correction" << std::endl;
  fRingArray.Print(option);
  fVertexAxis.Print(option);
}
    
//____________________________________________________________________
//
// EOF
//
