// $Header$

#include "PointSet.h"

#include <Reve/RGTopFrame.h>

#include <TTree.h>
#include <TF3.h>

#include <TColor.h>
#include <TCanvas.h>

using namespace Reve;

//______________________________________________________________________
// PointSet
//

ClassImp(PointSet)

PointSet::PointSet(Int_t n_points) :
  TPointSet3D(n_points),
  RenderElement(fMarkerColor)
{
  fMarkerStyle = 20;
}

PointSet::PointSet(const Text_t* name, Int_t n_points) :
  TPointSet3D(n_points),
  RenderElement(fMarkerColor)
{
  fMarkerStyle = 20;
  SetName(name);
}

PointSet::PointSet(const Text_t* name, TTree* tree,
		   TreeVarType_e tv_type) :
  TPointSet3D(tree->GetSelectedRows()),
  RenderElement(fMarkerColor)
{
  static const Exc_t eH("PointSet::PointSet ");

  fMarkerStyle = 20;
  SetName(name);
  Double_t *vx = tree->GetV1(), *vy = tree->GetV2(), *vz = tree->GetV3();
  Long64_t nr = tree->GetSelectedRows();

  switch(tv_type) {
  case TVT_XYZ:
    while(nr-- > 0) {
      SetNextPoint(*vx, *vy, *vz);
      ++vx; ++vy; ++vz;
    }
    break;
  case TVT_RPhiZ:
    while(nr-- > 0) {
      SetNextPoint(*vx * TMath::Cos(*vy), *vx * TMath::Sin(*vy), *vz);
      ++vx; ++vy; ++vz;
    }
    break;
  default:
    throw(eH + "unknown tree variable type.");
  }
}

/**************************************************************************/

void PointSet::Reset(Int_t n_points)
{
  delete [] fP; fP = 0;
  fN = n_points;
  if(fN) fP = new Float_t [3*fN];
  memset(fP, 0, 3*fN*sizeof(Float_t));
  fLastPoint = -1;
  ResetBBox();
}

/**************************************************************************/

void PointSet::Paint(Option_t* option)
{
  if(fRnrElement == false) return;

  TPointSet3D::Paint(option);
}

/**************************************************************************/
/**************************************************************************/

//______________________________________________________________________
// PointSetArray
//

ClassImp(PointSetArray)

PointSetArray::PointSetArray(const Text_t* name,
			     const Text_t* title) :
  TNamed(name, title), RenderElementListBase(fMarkerColor),
  fBins(0), fDefPointSetCapacity(128), fNBins(0)
{}

PointSetArray::~PointSetArray()
{
  DeleteBins();
}

/**************************************************************************/

void PointSetArray::SetMarkerColor(Color_t tcolor)
{
  for(lpRE_i i=fList.begin(); i!=fList.end(); ++i) {
    TAttMarker* m = dynamic_cast<TAttMarker*>((*i)->GetObject());
    if(m && m->GetMarkerColor() == fMarkerColor)
      m->SetMarkerColor(tcolor);
  }
  TAttMarker::SetMarkerColor(tcolor);
}

void PointSetArray::SetMarkerStyle(Style_t mstyle)
{
  for(lpRE_i i=fList.begin(); i!=fList.end(); ++i) {
    TAttMarker* m = dynamic_cast<TAttMarker*>((*i)->GetObject());
    if(m && m->GetMarkerStyle() == fMarkerStyle)
      m->SetMarkerStyle(mstyle);
  }
  TAttMarker::SetMarkerStyle(mstyle);
}

void PointSetArray::SetMarkerSize(Size_t msize)
{
  for(lpRE_i i=fList.begin(); i!=fList.end(); ++i) {
    TAttMarker* m = dynamic_cast<TAttMarker*>((*i)->GetObject());
    if(m && m->GetMarkerSize() == fMarkerSize)
      m->SetMarkerSize(msize);
  }
  TAttMarker::SetMarkerSize(msize);
}

/**************************************************************************/

void PointSetArray::InitBins(TGListTreeItem* tree_item, const Text_t* quant_name,
			     Int_t nbins, Double_t min, Double_t max)
{
  static const Exc_t eH("PointSetArray::InitBins ");

  if(nbins < 1) throw(eH + "nbins < 1.");
  if(min > max) throw(eH + "min > max.");

  DeleteBins();

  fQuantName = quant_name;
  fNBins     = nbins;
  fMin = fCurMin = min;
  fMax = fCurMax = max;
  fBinWidth  = (fMax - fMin)/fNBins;

  fBins = new Reve::PointSet*[fNBins];
  for(Int_t i=0; i<fNBins; ++i) {
    fBins[i] = new Reve::PointSet
      (Form("Slice %d [%4.3lf, %4.3lf]", i, fMin + i*fBinWidth, fMin + (i+1)*fBinWidth),
       fDefPointSetCapacity);
    fBins[i]->SetMarkerColor(fMarkerColor);
    fBins[i]->SetMarkerStyle(fMarkerStyle);
    fBins[i]->SetMarkerSize(fMarkerSize);
    AddElement(fBins[i]);
    if(tree_item)
      gReve->AddRenderElement(tree_item, fBins[i]);
  }
}

void PointSetArray::DeleteBins()
{
  if(fBins) {
    for(Int_t i=0; i<fNBins; ++i)
      delete fBins[i];
    delete [] fBins;
    fBins = 0; fNBins = 0;
  }
  RemoveElements();
}

void PointSetArray::Fill(Double_t quant, Double_t x, Double_t y, Double_t z)
{
  Int_t bin    = Int_t( (quant - fMin)/fBinWidth );
  if(bin >= 0 && bin < fNBins)
    fBins[bin]->SetNextPoint(x, y, z);
}

void PointSetArray::Fill(TF3* , TTree* , TreeVarType_e )
{

}

void PointSetArray::CloseBins()
{
  for(Int_t i=0; i<fNBins; ++i) {
    fBins[i]->fN = fBins[i]->fLastPoint; // HACK! PolyMarker3D does half-management of array size.
    fBins[i]->ComputeBBox();
  }
}

/**************************************************************************/

void PointSetArray::SetRange(Double_t min, Double_t max)
{
  using namespace TMath;

  fCurMin = min; fCurMax = max;
  Int_t  low_b = (Int_t) Max(Double_t(0),       Floor((min-fMin)/fBinWidth));
  Int_t high_b = (Int_t) Min(Double_t(fNBins-1), Ceil((max-fMin)/fBinWidth));
  for(Int_t i=0; i<fNBins; ++i) {
    fBins[i]->SetRnrElement(i>=low_b && i<=high_b);
  }
}

/**************************************************************************/

#include <TGFrame.h>
#include <TGDoubleSlider.h>
#include <TGXYLayout.h>

void PointSetArray::MakeScrollbar()
{
  TGMainFrame* mf = new TGMainFrame(gClient->GetRoot(), 320, 60);

  TGDoubleHSlider* hs = new TGDoubleHSlider(mf);
  hs->SetRange(fMin, fMax);
  hs->SetPosition(fMin, fMax);
  hs->Resize(300, 25);
  mf->AddFrame(hs, new TGLayoutHints(kLHintsCenterX, 10, 10, 10, 10));

  hs->Connect("PositionChanged()", "Reve::PointSetArray",
	      this, "HandleScrollEvent()");

  mf->SetWindowName(fQuantName + " Selector");
  mf->MapSubwindows();
  mf->Resize(mf->GetDefaultSize()); // this is used here to init layout algorithm
  mf->MapWindow();
}

void PointSetArray::HandleScrollEvent()
{
  TGDoubleHSlider* hs = (TGDoubleHSlider*)gTQSender;

  Float_t min = hs->GetMinPosition(), max = hs->GetMaxPosition();
  printf("hslidor min=%f max=%f\n", min, max);
  SetRange(min, max);
}
