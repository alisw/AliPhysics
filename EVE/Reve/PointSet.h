// $Header$

#ifndef REVE_GuiPointSet_H
#define REVE_GuiPointSet_H

#include <Reve/PODs.h>
#include <Reve/RenderElement.h>
#include <Reve/NLTBases.h>
#include <Reve/TTreeTools.h>

#include <TPointSet3D.h>
#include <TArrayI.h>

class TTree;
class TF3;
class TGListTreeItem;

namespace Reve {

class PointSet : public RenderElement,
                 public TPointSet3D,
                 public TPointSelectorConsumer,
                 public NLTProjectable
{
  friend class PointSetArray;

protected:
  TString  fTitle;           // Title/tooltip of the PointSet.
  TArrayI *fIntIds;          // Optional array of integer ideices.
  Int_t    fIntIdsPerPoint;  // Number of integer indices assigned to each point.

  void AssertIntIdsSize();

public:
  PointSet(Int_t n_points=0, TreeVarType_e tv_type=TVT_XYZ);
  PointSet(const Text_t* name, Int_t n_points=0, TreeVarType_e tv_type=TVT_XYZ);
  virtual ~PointSet();

  virtual void ComputeBBox();

  void  Reset(Int_t n_points=0, Int_t n_int_ids=0);
  Int_t GrowFor(Int_t n_points);

  virtual const Text_t* GetTitle() const          { return fTitle; }
  virtual void          SetTitle(const Text_t* t) { fTitle = t; }

  Int_t  GetIntIdsPerPoint() const { return fIntIdsPerPoint; }
  Int_t* GetPointIntIds(Int_t p) const;
  Int_t  GetPointIntId(Int_t p, Int_t i) const;

  void   SetPointIntIds(Int_t* ids);
  void   SetPointIntIds(Int_t n, Int_t* ids);

  virtual void SetRnrElNameTitle(const Text_t* name, const Text_t* title);

  virtual void SetMarkerColor(Color_t col)
  { SetMainColor(col); }

  virtual void Paint(Option_t* option="");

  virtual void InitFill(Int_t subIdNum);
  virtual void TakeAction(TPointSelector*);

  virtual const TGPicture* GetListTreeIcon() { return RenderElement::fgListTreeIcons[3]; }

  virtual TClass* ProjectedClass() const;

  ClassDef(PointSet, 1); // Render element containing an array of 3D points.
}; // endclass PointSet

/**************************************************************************/

class PointSetArray : public RenderElement,
                      public TNamed,
                      public TAttMarker,
                      public TPointSelectorConsumer
{
  friend class PointSetArrayEditor;

  PointSetArray(const PointSetArray&);            // Not implemented
  PointSetArray& operator=(const PointSetArray&); // Not implemented

protected:
  PointSet**   fBins;                 //  Pointers to subjugated PointSet's.
  Int_t        fDefPointSetCapacity;  //  Default capacity of subjugated PointSet's.
  Int_t        fNBins;                //  Number of subjugated PointSet's.
  Int_t        fLastBin;              //! Index of the last filled PointSet.
  Double_t     fMin, fCurMin;         //  Overall and current minimum value of the separating quantity.
  Double_t     fMax, fCurMax;         //  Overall and current maximum value of the separating quantity.
  Double_t     fBinWidth;             //  Separating quantity bin-width.
  TString      fQuantName;            //  Name of the separating quantity.

public:
  PointSetArray(const Text_t* name="PointSetArray", const Text_t* title="");
  virtual ~PointSetArray();

  virtual void RemoveElementLocal(RenderElement* el);
  virtual void RemoveElementsLocal();

  virtual void Paint(Option_t* option="");

  virtual void SetMarkerColor(Color_t tcolor=1);
  virtual void SetMarkerStyle(Style_t mstyle=1);
  virtual void SetMarkerSize(Size_t msize=1);

  virtual void TakeAction(TPointSelector*);


  void InitBins(const Text_t* quant_name, Int_t nbins, Double_t min, Double_t max,
		Bool_t addRe=kTRUE);
  void Fill(Double_t x, Double_t y, Double_t z, Double_t quant);
  void SetPointId(TObject* id);
  void CloseBins();

  void SetOwnIds(Bool_t o);

  Int_t GetDefPointSetCapacity() const  { return fDefPointSetCapacity; }
  void  SetDefPointSetCapacity(Int_t c) { fDefPointSetCapacity = c; }

  Int_t     GetNBins()        const { return fNBins; }
  PointSet* GetBin(Int_t bin) const { return fBins[bin]; }

  Double_t GetMin()    const { return fMin; }
  Double_t GetCurMin() const { return fCurMin; }
  Double_t GetMax()    const { return fMax; }
  Double_t GetCurMax() const { return fCurMax; }

  void SetRange(Double_t min, Double_t max);

  ClassDef(PointSetArray, 1); // Array of centrally managed PointSet's.
};

/**************************************************************************/

class NLTPointSet : public PointSet,
                    public NLTProjected
{
private:
  NLTPointSet(const NLTPointSet&);            // Not implemented
  NLTPointSet& operator=(const NLTPointSet&); // Not implemented

protected:

public:
  NLTPointSet();
  virtual ~NLTPointSet() {}

  virtual void SetProjection(NLTProjector* proj, NLTProjectable* model);

  virtual void UpdateProjection();

  ClassDef(NLTPointSet, 1); // NLT projected PointSet.
}; // endclass NLTPointSet

}

#endif
