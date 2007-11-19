// $Header$

#include "PointSet.h"

#include <Reve/ReveManager.h>
#include <Reve/NLTProjector.h>

#include <TTree.h>
#include <TTreePlayer.h>
#include <TF3.h>

#include <TColor.h>
#include <TCanvas.h>

using namespace Reve;

//______________________________________________________________________
// PointSet
//
// PointSet is a render-element holding a collection of 3D points with
// optional per-point TRef and an arbitrary number of integer ids (to
// be used for signal, volume-id, track-id, etc).
//
// 3D point representation is implemented in base-class TPolyMarker3D.
// Per-point TRef is implemented in base-class TPointSet3D.
//
// By using the TPointSelector the points and integer ids can be
// filled directly from a TTree holding the source data.
// Setting of per-point TRef's is not supported.
//
// PointSet is a NLTProjectable: it can be projected by using the
// NLTProjector class.

ClassImp(PointSet)

//______________________________________________________________________________
PointSet::PointSet(Int_t n_points, TreeVarType_e tv_type) :
  RenderElement(fMarkerColor),
  TPointSet3D(n_points),
  TPointSelectorConsumer(tv_type),

  fTitle          (),
  fIntIds         (0),
  fIntIdsPerPoint (0)
{
  // Constructor.

  fMarkerStyle = 20;
}

//______________________________________________________________________________
PointSet::PointSet(const Text_t* name, Int_t n_points, TreeVarType_e tv_type) :
  RenderElement(fMarkerColor),
  TPointSet3D(n_points),
  TPointSelectorConsumer(tv_type),

  fTitle          (),
  fIntIds         (0),
  fIntIdsPerPoint (0)
{
  // Constructor.

  fMarkerStyle = 20;
  SetName(name);
}

//______________________________________________________________________________
PointSet::~PointSet()
{
  // Destructor.

  delete fIntIds;
}

/**************************************************************************/

//______________________________________________________________________________
void PointSet::ComputeBBox()
{
  // Override of virtual method from TAttBBox.

  TPointSet3D::ComputeBBox();
  AssertBBoxExtents(0.1);
}

//______________________________________________________________________________
void PointSet::Reset(Int_t n_points, Int_t n_int_ids)
{
  // Drop all data and set-up the data structures to recive new data.
  // n_points   specifies the initial size of the arrays.
  // n_int_ids  specifies the number of integer ids per point.

  delete [] fP; fP = 0;
  fN = n_points;
  if(fN) fP = new Float_t [3*fN];
  memset(fP, 0, 3*fN*sizeof(Float_t));
  fLastPoint = -1;
  ClearIds();
  delete fIntIds; fIntIds = 0;
  fIntIdsPerPoint = n_int_ids;
  if (fIntIdsPerPoint > 0) fIntIds = new TArrayI(fIntIdsPerPoint*fN);
  ResetBBox();
}

//______________________________________________________________________________
Int_t PointSet::GrowFor(Int_t n_points)
{
  // Resizes internal array to allow additional n_points to be stored.
  // Returns the old size which is also the location where one can
  // start storing new data.
  // The caller is *obliged* to fill the new point slots.

  Int_t old_size = Size();
  Int_t new_size = old_size + n_points;
  SetPoint(new_size - 1, 0, 0, 0);
  if (fIntIds)
    fIntIds->Set(fIntIdsPerPoint * new_size);
  return old_size;
}

/**************************************************************************/

//______________________________________________________________________________
inline void PointSet::AssertIntIdsSize()
{
  // Assert that size of IntId array is compatible with the size of
  // the point array.

  Int_t exp_size = GetN()*fIntIdsPerPoint;
  if (fIntIds->GetSize() < exp_size)
    fIntIds->Set(exp_size);
}

//______________________________________________________________________________
Int_t* PointSet::GetPointIntIds(Int_t p) const
{
  // Return a pointer to integer ids of point with index p.
  // Existence of integer id array is checked, 0 is returned if it
  // does not exist.
  // Validity of p is *not* checked.

  if (fIntIds)
    return fIntIds->GetArray() + p*fIntIdsPerPoint;
  return 0;
}

//______________________________________________________________________________
Int_t PointSet::GetPointIntId(Int_t p, Int_t i) const
{
  // Return i-th integer id of point with index p.
  // Existence of integer id array is checked, kMinInt is returned if
  // it does not exist.
  // Validity of p and i is *not* checked.

  if (fIntIds)
    return * (fIntIds->GetArray() + p*fIntIdsPerPoint + i);
  return kMinInt;
}

//______________________________________________________________________________
void PointSet::SetPointIntIds(Int_t* ids)
{
  // Set integer ids for the last point that was registerd (most
  // probably via TPolyMarker3D::SetNextPoint(x,y,z)).

  SetPointIntIds(fLastPoint, ids);
}

//______________________________________________________________________________
void PointSet::SetPointIntIds(Int_t n, Int_t* ids)
{
  // Set integer ids for point with index n.

  if (!fIntIds) return;
  AssertIntIdsSize();
  Int_t* x = fIntIds->GetArray() + n*fIntIdsPerPoint;
  for (Int_t i=0; i<fIntIdsPerPoint; ++i)
    x[i] = ids[i];
}

/**************************************************************************/

//______________________________________________________________________________
void PointSet::SetRnrElNameTitle(const Text_t* name, const Text_t* title)
{
  // Set name and title of point-set.
  // Virtual in RenderElement.

  SetName(name);
  SetTitle(title);
}

/******************************************************************************/

//______________________________________________________________________________
void PointSet::Paint(Option_t* option)
{
  // Paint point-set.

  if(fRnrSelf == kFALSE) return;

  TPointSet3D::Paint(option);
}

/**************************************************************************/

//______________________________________________________________________________
void PointSet::InitFill(Int_t subIdNum)
{
  // Initialize point-set for new filling.
  // subIdNum gives the number of integer ids that can be assigned to
  // each point.

  if (subIdNum > 0) {
    fIntIdsPerPoint = subIdNum;
    if (!fIntIds)
      fIntIds = new TArrayI(fIntIdsPerPoint*GetN());
    else
      fIntIds->Set(fIntIdsPerPoint*GetN());
  } else {
    delete fIntIds; fIntIds = 0;
    fIntIdsPerPoint = 0;
  }
}

//______________________________________________________________________________
void PointSet::TakeAction(TPointSelector* sel)
{
  // Called from TPointSelector when internal arrays of the tree-selector
  // are filled up and need to be processed.
  // Virtual from TPointSelectorConsumer.

  static const Exc_t eH("PointSet::TakeAction ");

  if(sel == 0)
    throw(eH + "selector is <null>.");

  Int_t    n = sel->GetNfill();
  Int_t  beg = GrowFor(n);

  // printf("PointSet::TakeAction beg=%d n=%d size=%d nsubid=%d dim=%d\n",
  //        beg, n, Size(), sel->GetSubIdNum(), sel->GetDimension());

  Double_t *vx = sel->GetV1(), *vy = sel->GetV2(), *vz = sel->GetV3();
  Float_t  *p  = fP + 3*beg;

  switch(fSourceCS) {
  case TVT_XYZ:
    while(n-- > 0) {
      p[0] = *vx; p[1] = *vy; p[2] = *vz;
      p += 3;
      ++vx; ++vy; ++vz;
    }
    break;
  case TVT_RPhiZ:
    while(n-- > 0) {
      p[0] = *vx * TMath::Cos(*vy); p[1] = *vx * TMath::Sin(*vy); p[2] = *vz;
      p += 3;
      ++vx; ++vy; ++vz;
    }
    break;
  default:
    throw(eH + "unknown tree variable type.");
  }

  if (fIntIds) {
    Double_t** subarr = new Double_t* [fIntIdsPerPoint];
    for (Int_t i=0; i<fIntIdsPerPoint; ++i) {
      subarr[i] = sel->GetVal(sel->GetDimension() - fIntIdsPerPoint + i);
      if (subarr[i] == 0)
        throw(eH + "sub-id array not available.");
    }
    Int_t* ids = fIntIds->GetArray() + fIntIdsPerPoint*beg;
    n = sel->GetNfill();
    while (n-- > 0) {
      for (Int_t i=0; i<fIntIdsPerPoint; ++i) {
        ids[i] = TMath::Nint(*subarr[i]);
        ++subarr[i];
      }
      ids += fIntIdsPerPoint;
    }
    delete [] subarr;
  }
}

/**************************************************************************/

//______________________________________________________________________________
TClass* PointSet::ProjectedClass() const
{
  // Virtual from NLTProjectable, returns NLTPointSet class.

  return NLTPointSet::Class();
}


/**************************************************************************/
/**************************************************************************/

//______________________________________________________________________________
// PointSetArray
//
// An array of point-sets with each point-set playing a role of a bin
// in a histogram. When a new point is added to a PointSetArray, an
// additional separating quantity needs to be specified: it determines
// into which PointSet (bin) the point will actually be stored.
//
// By using the TPointSelector the points and the separating
// quantities can be filled directly from a TTree holding the source
// data.
// Setting of per-point TRef's is not supported.
//
// After the filling, the range of separating variable can be
// controlled with a slider to choose a sub-set of PointSets that are
// actually shown.
//

ClassImp(PointSetArray)

//______________________________________________________________________________
PointSetArray::PointSetArray(const Text_t* name,
			     const Text_t* title) :
  RenderElement(fMarkerColor),
  TNamed(name, title),

  fBins(0), fDefPointSetCapacity(128), fNBins(0), fLastBin(-1),
  fMin(0), fCurMin(0), fMax(0), fCurMax(0),
  fBinWidth(0),
  fQuantName()
{
  // Constructor.
}

//______________________________________________________________________________
PointSetArray::~PointSetArray()
{
  // Destructor: deletes the fBins array. Actual removal of
  // elements done by RenderElement.

  // printf("PointSetArray::~PointSetArray()\n");
  delete [] fBins; fBins = 0;
}

//______________________________________________________________________________
void PointSetArray::Paint(Option_t* option)
{
  // Paint the subjugated PointSet's.

  if (fRnrSelf) {
    for (List_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
      if ((*i)->GetRnrSelf())
	(*i)->GetObject()->Paint(option);
    }
  }
}

//______________________________________________________________________________
void PointSetArray::RemoveElementLocal(RenderElement* el)
{
  // Virtual from RenderElement, provide bin management.

  for (Int_t i=0; i<fNBins; ++i) {
    if (fBins[i] == el) {
      fBins[i] = 0;
      break;
    }
  }
}

//______________________________________________________________________________
void PointSetArray::RemoveElementsLocal()
{
  // Virtual from RenderElement, provide bin management.

  delete [] fBins; fBins = 0; fLastBin = -1;
}

/**************************************************************************/

//______________________________________________________________________________
void PointSetArray::SetMarkerColor(Color_t tcolor)
{
  for (List_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
    TAttMarker* m = dynamic_cast<TAttMarker*>((*i)->GetObject());
    if (m && m->GetMarkerColor() == fMarkerColor)
      m->SetMarkerColor(tcolor);
  }
  TAttMarker::SetMarkerColor(tcolor);
}

//______________________________________________________________________________
void PointSetArray::SetMarkerStyle(Style_t mstyle)
{
  for (List_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
    TAttMarker* m = dynamic_cast<TAttMarker*>((*i)->GetObject());
    if (m && m->GetMarkerStyle() == fMarkerStyle)
      m->SetMarkerStyle(mstyle);
  }
  TAttMarker::SetMarkerStyle(mstyle);
}

//______________________________________________________________________________
void PointSetArray::SetMarkerSize(Size_t msize)
{
  for (List_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
    TAttMarker* m = dynamic_cast<TAttMarker*>((*i)->GetObject());
    if (m && m->GetMarkerSize() == fMarkerSize)
      m->SetMarkerSize(msize);
  }
  TAttMarker::SetMarkerSize(msize);
}

/**************************************************************************/

//______________________________________________________________________________
void PointSetArray::TakeAction(TPointSelector* sel)
{
  // Called from TPointSelector when internal arrays of the tree-selector
  // are filled up and need to be processed.
  // Virtual from TPointSelectorConsumer.

  static const Exc_t eH("PointSetArray::TakeAction ");

  if (sel == 0)
    throw(eH + "selector is <null>.");

  Int_t n = sel->GetNfill();

  // printf("PointSetArray::TakeAction n=%d\n", n);

  Double_t *vx = sel->GetV1(), *vy = sel->GetV2(), *vz = sel->GetV3();
  Double_t *qq = sel->GetV4();

  if(qq == 0)
    throw(eH + "requires 4-d varexp.");

  switch(fSourceCS) {
  case TVT_XYZ:
    while(n-- > 0) {
      Fill(*vx, *vy, *vz, *qq);
      ++vx; ++vy; ++vz; ++qq;
    }
    break;
  case TVT_RPhiZ:
    while(n-- > 0) {
      Fill(*vx * TMath::Cos(*vy), *vx * TMath::Sin(*vy), *vz, *qq);
      ++vx; ++vy; ++vz; ++qq;
    }
    break;
  default:
    throw(eH + "unknown tree variable type.");
  }
}

/**************************************************************************/

//______________________________________________________________________________
void PointSetArray::InitBins(const Text_t* quant_name,
			     Int_t nbins, Double_t min, Double_t max,
			     Bool_t addRe)
{
  static const Exc_t eH("PointSetArray::InitBins ");

  if (nbins < 1) throw(eH + "nbins < 1.");
  if (min > max) throw(eH + "min > max.");

  RemoveElements();

  fQuantName = quant_name;
  fNBins     = nbins;
  fLastBin   = -1;
  fMin = fCurMin = min;
  fMax = fCurMax = max;
  fBinWidth  = (fMax - fMin)/fNBins;

  fBins = new Reve::PointSet*[fNBins];
  for (Int_t i=0; i<fNBins; ++i) {
    fBins[i] = new Reve::PointSet
      (Form("Slice %d [%4.3lf, %4.3lf]", i, fMin + i*fBinWidth, fMin + (i+1)*fBinWidth),
       fDefPointSetCapacity);
    fBins[i]->SetMarkerColor(fMarkerColor);
    fBins[i]->SetMarkerStyle(fMarkerStyle);
    fBins[i]->SetMarkerSize(fMarkerSize);
    if (addRe)
      gReve->AddRenderElement(fBins[i], this);
    else
      AddElement(fBins[i]);
  }
}

//______________________________________________________________________________
void PointSetArray::Fill(Double_t x, Double_t y, Double_t z, Double_t quant)
{
  fLastBin = Int_t( (quant - fMin)/fBinWidth );
  if (fLastBin >= 0 && fLastBin < fNBins && fBins[fLastBin] != 0)
    fBins[fLastBin]->SetNextPoint(x, y, z);
  else
    fLastBin = -1;
}

//______________________________________________________________________________
void PointSetArray::SetPointId(TObject* id)
{
  if (fLastBin >= 0)
    fBins[fLastBin]->SetPointId(id);
}

//______________________________________________________________________________
void PointSetArray::CloseBins()
{
  for (Int_t i=0; i<fNBins; ++i) {
    if (fBins[i] != 0) {
      // HACK! PolyMarker3D does half-management of array size.
      // In fact, the error is mine, in pointset3d(gl) i use fN instead of Size().
      // Fixed in my root, but not elsewhere.
      fBins[i]->fN = fBins[i]->fLastPoint;

      fBins[i]->ComputeBBox();
    }
  }
  fLastBin = -1;
}

/**************************************************************************/

//______________________________________________________________________________
void PointSetArray::SetOwnIds(Bool_t o)
{
  for (Int_t i=0; i<fNBins; ++i)
  {
    if (fBins[i] != 0)
      fBins[i]->SetOwnIds(o);
  }
}

/**************************************************************************/

//______________________________________________________________________________
void PointSetArray::SetRange(Double_t min, Double_t max)
{
  using namespace TMath;

  fCurMin = min; fCurMax = max;
  Int_t  low_b = (Int_t) Max(Double_t(0),       Floor((min-fMin)/fBinWidth));
  Int_t high_b = (Int_t) Min(Double_t(fNBins-1), Ceil((max-fMin)/fBinWidth));
  for (Int_t i=0; i<fNBins; ++i) {
    if (fBins[i] != 0)
      fBins[i]->SetRnrSelf(i>=low_b && i<=high_b);
  }
}


/******************************************************************************/
/******************************************************************************/

//______________________________________________________________________________
// NLTPointSet
//

ClassImp(NLTPointSet)

//______________________________________________________________________________
NLTPointSet::NLTPointSet() :
  PointSet     (),
  NLTProjected ()
{
  // Default contructor.
}

//______________________________________________________________________________
void NLTPointSet::SetProjection(NLTProjector* proj, NLTProjectable* model)
{
  NLTProjected::SetProjection(proj, model);

  * (TAttMarker*)this = * dynamic_cast<TAttMarker*>(fProjectable);
}

//______________________________________________________________________________
void NLTPointSet::UpdateProjection()
{
  NLTProjection& proj = * fProjector->GetProjection();
  PointSet     & ps   = * dynamic_cast<PointSet*>(fProjectable);

  Int_t n = ps.GetN();
  Reset(n);
  Float_t *o = ps.GetP(), *p = GetP();
  for (Int_t i = 0; i < n; ++i, o+=3, p+=3)
  {
    p[0] = o[0]; p[1] = o[1]; p[2] = o[2];
    proj.ProjectPoint(p[0], p[1], p[2]);
    p[2] = fDepth;
  }
  fLastPoint = n - 1;
}
