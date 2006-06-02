// $Header$

#ifndef REVE_RenderElement_H
#define REVE_RenderElement_H

#include <Reve/Reve.h>

#include <TNamed.h>
#include <TRef.h>

#include <Gtypes.h>
#include <GuiTypes.h>

#include <list>
#include <set>

class TGListTree;
class TGListTreeItem;

namespace Reve {

class RenderElement
{
  friend class RGTopFrame;

protected:
  // TRef     fSource;

  Bool_t   fRnrElement;
  Color_t* fMainColorPtr;

  struct ListTreeInfo {
    TGListTree*     fTree;
    TGListTreeItem* fItem;

    ListTreeInfo() {}
    ListTreeInfo(TGListTree* lt, TGListTreeItem* lti) : fTree(lt), fItem(lti) {}
    virtual ~ListTreeInfo() {}

    bool operator==(const ListTreeInfo& x) const
    { return fTree == x.fTree && fItem == x.fItem; }
    bool operator<(const ListTreeInfo& x) const
    { return fTree == x.fTree ? fItem < x.fItem : fTree < x.fTree; }

    ClassDef(ListTreeInfo, 0);
  };
  typedef std::set<ListTreeInfo>           sLTI_t;
  typedef std::set<ListTreeInfo>::iterator sLTI_i;

  sLTI_t fItems;

public:
  RenderElement();
  RenderElement(Color_t& main_color);
  virtual ~RenderElement();

  virtual TObject* GetObject(Reve::Exc_t eh="RenderElement::GetObject ");

  /*
    TRef&    GetSource() { return fSource; }
    TObject* GetSourceObject() const { return fSource.GetObject(); }
    void SetSourceObject(TObject* o) { fSource.SetObject(o); }

    void DumpSourceObject();    // *MENU*
    void InspectSourceObject(); // *MENU*
  */

  virtual TGListTreeItem* AddIntoListTree(TGListTree* ltree,
					  TGListTreeItem* parent);

  virtual void FullUpdate();

  void SpawnEditor();                  // *MENU*
  void ExportToCINT(Text_t* var_name); // *MENU*

  virtual Bool_t CanEditRnrElement()   { return kTRUE; }
  virtual Bool_t GetRnrElement() const { return fRnrElement; }
  virtual void   SetRnrElement(Bool_t rnr);

  virtual Bool_t CanEditMainColor()        { return kFALSE; }
  Color_t* GetMainColorPtr()               { return fMainColorPtr; }
  void     SetMainColorPtr(Color_t* color) { fMainColorPtr = color; }

  virtual Color_t GetMainColor() const { return fMainColorPtr ? *fMainColorPtr : 0; }
  virtual void    SetMainColor(Color_t color);
          void    SetMainColorByPixel(Pixel_t pixel);

  ClassDef(RenderElement, 1);
}; // endclass RenderElement

/**************************************************************************/

class RenderElementObjPtr : public RenderElement
{
protected:
  TObject* fObject;

public:
  RenderElementObjPtr(TObject* obj);
  RenderElementObjPtr(TObject* obj, Color_t& mainColor);
  virtual ~RenderElementObjPtr();

  virtual TObject* GetObject(Reve::Exc_t eh="RenderElementObjPtr::GetObject ");

  virtual void SetRnrElement(Bool_t rnr);

  ClassDef(RenderElementObjPtr, 1);
}; // endclass RenderElementObjPtr

/**************************************************************************/

class RenderElementListBase : public RenderElement
{
protected:
  typedef std::list<RenderElement*>            lpRE_t;
  typedef std::list<RenderElement*>::iterator  lpRE_i;

  lpRE_t fList;

  void PaintElements(Option_t* option="");

public:
  RenderElementListBase() {}
  RenderElementListBase(Color_t& col) : RenderElement(col) {}
  virtual ~RenderElementListBase() {}

  virtual void AddElement(RenderElement* el) { fList.push_back(el); }

  virtual Int_t ExpandIntoListTree(TGListTree* ltree, TGListTreeItem* parent);

  void EnableListElements();   // *MENU*
  void DisableListElements();  // *MENU*

  virtual void SetMainColor(Color_t color);
  virtual void SetMainColor(Pixel_t pixel);

  ClassDef(RenderElementListBase, 1);
};

/**************************************************************************/

class RenderElementList : public TNamed, public RenderElementListBase
{
protected:
  Color_t   fColor;

public:
  RenderElementList(const Text_t* n="RenderElementList", const Text_t* t="") :
    TNamed(n, t), RenderElementListBase(fColor), fColor(0)
  {}
  virtual ~RenderElementList() {}

  virtual Bool_t CanEditMainColor()  { return kTRUE; }

  virtual void Paint(Option_t* option = "") { PaintElements(option); }

  ClassDef(RenderElementList, 1);
};

}

#endif
