// $Header$

#ifndef REVE_ReveGeom_H
#define REVE_ReveGeom_H

#include <Reve/RenderElement.h>

class TGeoVolume;
class TGeoNode;
class TGeoHMatrix;
class TGeoManager;

namespace Reve {

class GeoNodeRnrEl : public RenderElementListBase,
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

  virtual void UpdateItems();

  virtual Bool_t CanEditRnrElement() { return false; }
  virtual void SetRnrElement(Bool_t rnr);

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
  TGeoHMatrix* fGlobalTrans;
  Bool_t       fUseNodeTrans;
  Int_t        fVisOption;
  Int_t        fVisLevel;  

public:
  GeoTopNodeRnrEl(TGeoManager* manager, TGeoNode* node, Int_t visopt=1, Int_t vislvl=3);
  virtual ~GeoTopNodeRnrEl();

  TGeoHMatrix *GetGlobalTrans()  const { return fGlobalTrans; }
  void         SetGlobalTrans(TGeoHMatrix* m);
  Bool_t       GetUseNodeTrans() const { return fUseNodeTrans; }
  void         SetUseNodeTrans(Bool_t u=kTRUE);

  Int_t GetVisOption() const { return fVisOption; }
  void  SetVisOption(Int_t visopt);
  Int_t GetVisLevel()  const { return fVisLevel; }
  void  SetVisLevel(Int_t vislvl);

  virtual void UpdateItems();

  virtual Bool_t CanEditRnrElement() { return true; }
  virtual void SetRnrElement(Bool_t rnr);

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

}

#endif
