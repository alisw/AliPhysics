// $Header$

#ifndef REVE_ReveGeom_H
#define REVE_ReveGeom_H

#include <Reve/RenderElement.h>
#include <Reve/ZTrans.h>
#include <Reve/NLTBases.h>

class TGeoVolume;
class TGeoNode;
class TGeoHMatrix;
class TGeoManager;

class TGeoShape;
class TGeoShapeExtract;

namespace Reve {

class GeoNodeRnrEl : public RenderElement,
                     public TObject
{
  friend class GeoNodeRnrElEditor;

  GeoNodeRnrEl(const GeoNodeRnrEl&);            // Not implemented
  GeoNodeRnrEl& operator=(const GeoNodeRnrEl&); // Not implemented

protected:
  TGeoNode *fNode;

public:
  GeoNodeRnrEl(TGeoNode* node);

  virtual const Text_t* GetName()  const;
  virtual const Text_t* GetTitle() const;

  TGeoNode* GetNode() const { return fNode; }

  virtual Int_t ExpandIntoListTree(TGListTree* ltree, TGListTreeItem* parent);

  virtual Bool_t CanEditRnrElement() { return false; }
  virtual void SetRnrSelf(Bool_t rnr);
  virtual void SetRnrChildren(Bool_t rnr);
  virtual void SetRnrState(Bool_t rnr);

  virtual Bool_t CanEditMainColor()  { return true; }
  virtual void SetMainColor(Color_t color);
  virtual void SetMainColor(Pixel_t pixel);

  void UpdateNode(TGeoNode* node);
  void UpdateVolume(TGeoVolume* volume);

  virtual void Draw(Option_t* option="");

  ClassDef(GeoNodeRnrEl, 1);
}; // endclass GeoNodeRnrEl

//----------------------------------------------------------------

class GeoTopNodeRnrEl : public GeoNodeRnrEl
{
  GeoTopNodeRnrEl(const GeoTopNodeRnrEl&);            // Not implemented
  GeoTopNodeRnrEl& operator=(const GeoTopNodeRnrEl&); // Not implemented

protected:
  TGeoManager* fManager;
  ZTrans       fGlobalTrans;
  Int_t        fVisOption;
  Int_t        fVisLevel;  

public:
  GeoTopNodeRnrEl(TGeoManager* manager, TGeoNode* node, Int_t visopt=1, Int_t vislvl=3);
  virtual ~GeoTopNodeRnrEl();

  virtual Bool_t  CanEditMainHMTrans() { return  kTRUE; }
  virtual ZTrans* PtrMainHMTrans()     { return &fGlobalTrans; }

  ZTrans&      RefGlobalTrans() { return fGlobalTrans; }
  void         SetGlobalTrans(const TGeoHMatrix* m);
  void         UseNodeTrans();

  Int_t GetVisOption() const { return fVisOption; }
  void  SetVisOption(Int_t visopt);
  Int_t GetVisLevel()  const { return fVisLevel; }
  void  SetVisLevel(Int_t vislvl);

  virtual Bool_t CanEditRnrElement() { return true; }
  virtual void SetRnrSelf(Bool_t rnr);

  virtual void Draw(Option_t* option="");
  virtual void Paint(Option_t* option="");

  // Signals from GeoManager.
  // These are not available any more ... colors in list-tree not refreshed
  // properly.
  void VolumeVisChanged(TGeoVolume* volume);
  void VolumeColChanged(TGeoVolume* volume);
  void NodeVisChanged(TGeoNode* node);

  ClassDef(GeoTopNodeRnrEl, 1);
}; // endclass GeoTopNodeRnrEl


//----------------------------------------------------------------
//----------------------------------------------------------------

class GeoShapeRnrEl : public RenderElement,
		      public TNamed,
                      public NLTGeoProjectable
{
  GeoShapeRnrEl(const GeoShapeRnrEl&);            // Not implemented
  GeoShapeRnrEl& operator=(const GeoShapeRnrEl&); // Not implemented

protected:
  ZTrans            fHMTrans;
  Color_t           fColor;
  UChar_t           fTransparency;
  TGeoShape*        fShape;

  static GeoShapeRnrEl* SubImportShapeExtract(TGeoShapeExtract* gse, RenderElement* parent);

public:
  GeoShapeRnrEl(const Text_t* name="GeoShapeRnrEl", const Text_t* title=0);
  virtual ~GeoShapeRnrEl();

  virtual Bool_t CanEditMainColor() { return kTRUE; }

  virtual Bool_t  CanEditMainTransparency()      { return kTRUE; }
  virtual UChar_t GetMainTransparency() const    { return fTransparency; }
  virtual void    SetMainTransparency(UChar_t t) { fTransparency = t; }  

  virtual Bool_t  CanEditMainHMTrans() { return  kTRUE; }
  virtual ZTrans* PtrMainHMTrans()     { return &fHMTrans; }

  ZTrans& RefHMTrans() { return fHMTrans; }
  void SetTransMatrix(Double_t* carr)        { fHMTrans.SetFrom(carr); }
  void SetTransMatrix(const TGeoMatrix& mat) { fHMTrans.SetFrom(mat);  }

  Color_t     GetColor()        { return fColor; }
  TGeoShape*  GetShape()        { return fShape; }

  virtual void Paint(Option_t* option="");

  static GeoShapeRnrEl* ImportShapeExtract(TGeoShapeExtract* gse, RenderElement* parent);
  
  // NLTGeoProjectable
  virtual TBuffer3D*           MakeBuffer3D();

  ClassDef(GeoShapeRnrEl, 1);
};

}

#endif
