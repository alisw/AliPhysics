// $Header$

#ifndef REVE_GuiPointSet_H
#define REVE_GuiPointSet_H

#include <Reve/PODs.h>
#include <Reve/RenderElement.h>
#include <Reve/TTreeTools.h>

#include <TPointSet3D.h>

class TTree;
class TF3;
class TGListTreeItem;

namespace Reve {

class PointSet : public RenderElement,
                 public TPointSet3D,
                 public TPointSelectorConsumer
{
  friend class PointSetArray;

protected:
  TString fTitle;

public:
  PointSet(Int_t n_points=0, TreeVarType_e tv_type=TVT_XYZ);
  PointSet(const Text_t* name, Int_t n_points=0, TreeVarType_e tv_type=TVT_XYZ);
  PointSet(const Text_t* name, TTree* tree, TreeVarType_e tv_type=TVT_XYZ);

  void  Reset(Int_t n_points=0);
  Int_t GrowFor(Int_t n_points);

  virtual const Text_t* GetTitle() const          { return fTitle; }
  virtual void          SetTitle(const Text_t* t) { fTitle = t; }

  virtual void SetMarkerColor(Color_t col)
  { SetMainColor(col); }

  virtual void Paint(Option_t* option="");

  virtual void TakeAction(TSelectorDraw*);

  ClassDef(PointSet, 1);
}; // endclass PointSet

/**************************************************************************/

class PointSetArray : public RenderElementListBase,
                      public TNamed,
                      public TAttMarker,
                      public TPointSelectorConsumer
{
  friend class PointSetArrayEditor;

  PointSetArray(const PointSetArray&);            // Not implemented
  PointSetArray& operator=(const PointSetArray&); // Not implemented

protected:
  PointSet**   fBins;
  Int_t        fDefPointSetCapacity;
  Int_t        fNBins;
  Double_t     fMin, fCurMin;
  Double_t     fMax, fCurMax;
  Double_t     fBinWidth;
  TString      fQuantName;

public:
  PointSetArray(const Text_t* name="PointSetArray", const Text_t* title="");
  virtual ~PointSetArray();

  virtual void RemoveElementLocal(RenderElement* el);
  virtual void RemoveElements();

  virtual void Paint(Option_t* option="") { PaintElements(option); }

  virtual void SetMarkerColor(Color_t tcolor=1);
  virtual void SetMarkerStyle(Style_t mstyle=1);
  virtual void SetMarkerSize(Size_t msize=1);

  virtual void TakeAction(TSelectorDraw*);


  void InitBins(const Text_t* quant_name, Int_t nbins, Double_t min, Double_t max,
		Bool_t addRe=kTRUE);
  void Fill(Double_t x, Double_t y, Double_t z, Double_t quant);

  void CloseBins();

  Int_t GetDefPointSetCapacity() const  { return fDefPointSetCapacity; }
  void  SetDefPointSetCapacity(Int_t c) { fDefPointSetCapacity = c; }

  Int_t     GetNBins()        const { return fNBins; }
  PointSet* GetBin(Int_t bin) const { return fBins[bin]; }

  Double_t GetMin()    const { return fMin; }
  Double_t GetCurMin() const { return fCurMin; }
  Double_t GetMax()    const { return fMax; }
  Double_t GetCurMax() const { return fCurMax; }

  void SetRange(Double_t min, Double_t max);

  ClassDef(PointSetArray, 1);
};

}

#endif
