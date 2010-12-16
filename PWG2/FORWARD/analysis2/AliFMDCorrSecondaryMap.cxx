#include "AliFMDCorrSecondaryMap.h"
#include <TBrowser.h>
#include <TH2D.h>
#include <AliLog.h>
#include <iostream>

//____________________________________________________________________
AliFMDCorrSecondaryMap::AliFMDCorrSecondaryMap()
  : fRingArray(), 
    fVertexAxis(0,0,0),
    fEtaAxis(0,0,0)
{
  fRingArray.SetOwner(kTRUE);
  fRingArray.SetName("rings");
  fVertexAxis.SetName("vertexAxis");
  fVertexAxis.SetTitle("v_{z} [cm]");
  fEtaAxis.SetName("etaAxis");
  fEtaAxis.SetTitle("#eta");
  
}
//____________________________________________________________________
AliFMDCorrSecondaryMap::AliFMDCorrSecondaryMap(const 
					       AliFMDCorrSecondaryMap& o)
  : TObject(o), 
    fRingArray(o.fRingArray), 
    fVertexAxis(o.fVertexAxis.GetNbins(), o.fVertexAxis.GetXmin(), 
		o.fVertexAxis.GetXmax()),
    fEtaAxis(o.fEtaAxis.GetNbins(), o.fEtaAxis.GetXmin(), 
	     o.fEtaAxis.GetXmax())
{
  fVertexAxis.SetName("vertexAxis");
  fVertexAxis.SetTitle("v_{z} [cm]");
  fEtaAxis.SetName("etaAxis");
  fEtaAxis.SetTitle("v_{z} [cm]");
}
//____________________________________________________________________
AliFMDCorrSecondaryMap::~AliFMDCorrSecondaryMap()
{
  fRingArray.Clear();
}
//____________________________________________________________________
AliFMDCorrSecondaryMap&
AliFMDCorrSecondaryMap::operator=(const AliFMDCorrSecondaryMap& o)
{
  fRingArray        = o.fRingArray;
  SetVertexAxis(o.fVertexAxis);
  SetEtaAxis(o.fEtaAxis);

  return *this;
}
//____________________________________________________________________
TH2D*
AliFMDCorrSecondaryMap::GetCorrection(UShort_t d, Char_t r, Double_t v) const
{
  Int_t b = FindVertexBin(v);
  if (b <= 0) return 0;
  return GetCorrection(d, r, UShort_t(b));
}
//____________________________________________________________________
TH2D*
AliFMDCorrSecondaryMap::GetCorrection(UShort_t d, Char_t r, UShort_t b) const
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
  return static_cast<TH2D*>(o);
}
  
//____________________________________________________________________
Int_t
AliFMDCorrSecondaryMap::FindVertexBin(Double_t v) const
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
AliFMDCorrSecondaryMap::GetRingIndex(UShort_t d, Char_t r) const
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
AliFMDCorrSecondaryMap::GetRingArray(UShort_t d, Char_t r) const
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
AliFMDCorrSecondaryMap::GetOrMakeRingArray(UShort_t d, Char_t r)
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
AliFMDCorrSecondaryMap::SetCorrection(UShort_t d, Char_t r, 
				      UShort_t b, TH2D*  h) 
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
  h->SetYTitle("#phi [radians]");
  h->SetZTitle("#sum_{i} N_{ch,i,primary} / #sum_{i} N_{ch,i,FMD}");
  h->SetDirectory(0);
  h->SetStats(0);
  ringArray->AddAtAndExpand(h, b-1);
  return kTRUE;
}
//____________________________________________________________________
Bool_t
AliFMDCorrSecondaryMap::SetCorrection(UShort_t d, Char_t r, 
				      Double_t v, TH2D*  h) 
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
AliFMDCorrSecondaryMap::Browse(TBrowser* b)
{
  b->Add(&fRingArray);
  b->Add(&fVertexAxis);
  b->Add(&fEtaAxis);
}
//____________________________________________________________________
void
AliFMDCorrSecondaryMap::Print(Option_t* option) const
{
  std::cout << "Secondary correction map" << std::endl;
  fRingArray.Print(option);
  fVertexAxis.Print(option);
  fEtaAxis.Print(option);
}
    
//____________________________________________________________________
//
// EOF
//
