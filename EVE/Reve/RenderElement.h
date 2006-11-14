// $Header$

#ifndef REVE_RenderElement_H
#define REVE_RenderElement_H

#include <Reve/Reve.h>

#include <TNamed.h>
#include <TRef.h>

#include <list>
#include <set>

class TGListTree;
class TGListTreeItem;

namespace Reve {

class RenderElement
{
  friend class RGTopFrame;

  RenderElement(const RenderElement&);            // Not implemented
  RenderElement& operator=(const RenderElement&); // Not implemented

public:
  class ListTreeInfo
  {
  public:
    TGListTree*     fTree;
    TGListTreeItem* fItem;

    ListTreeInfo() : fTree(0), fItem(0) {}
    ListTreeInfo(TGListTree* lt, TGListTreeItem* lti) : fTree(lt), fItem(lti) {}
    ListTreeInfo(const ListTreeInfo& l) : fTree(l.fTree), fItem(l.fItem) {}
    virtual ~ListTreeInfo() {}

    ListTreeInfo& operator=(const ListTreeInfo& l)
    { fTree = l.fTree; fItem = l.fItem; return *this; }

    bool operator==(const ListTreeInfo& x) const
    { return fTree == x.fTree && fItem == x.fItem; }
    bool operator<(const ListTreeInfo& x) const
    { return fTree == x.fTree ? fItem < x.fItem : fTree < x.fTree; }

    ClassDef(ListTreeInfo, 0);
  };
  typedef std::set<ListTreeInfo>               sLTI_t;
  typedef std::set<ListTreeInfo>::iterator     sLTI_i;

  typedef std::list<RenderElement*>            lpRE_t;
  typedef std::list<RenderElement*>::iterator  lpRE_i;

protected:
  // TRef     fSource;

  Bool_t   fRnrElement;
  Color_t* fMainColorPtr;

  sLTI_t fItems;
  lpRE_t fParents;

  Bool_t fDestroyOnZeroRefCnt;
  Bool_t fDenyDestroy;

public:
  RenderElement();
  RenderElement(Color_t& main_color);
  virtual ~RenderElement();

  virtual void AddParent(RenderElement* re);
  virtual void RemoveParent(RenderElement* re);

  Bool_t GetDestroyOnZeroRefCnt() const   { return fDestroyOnZeroRefCnt; }
  void   SetDestroyOnZeroRefCnt(Bool_t d) { fDestroyOnZeroRefCnt = d; }

  Bool_t GetDenyDestroy() const   { return fDenyDestroy; }
  void   SetDenyDestroy(Bool_t d) { fDenyDestroy = d; }

  virtual TObject* GetObject(Reve::Exc_t eh="RenderElement::GetObject ");

  /*
    TRef&    GetSource() { return fSource; }
    TObject* GetSourceObject() const { return fSource.GetObject(); }
    void SetSourceObject(TObject* o) { fSource.SetObject(o); }

    void DumpSourceObject();    // *MENU*
    void InspectSourceObject(); // *MENU*
  */

  virtual TGListTreeItem* AddIntoListTree(TGListTree* ltree,
					  TGListTreeItem* parent_lti);
  virtual TGListTreeItem* AddIntoListTree(TGListTree* ltree,
					  RenderElement* parent);
  virtual Bool_t          RemoveFromListTree(TGListTree* ltree);
  virtual Bool_t          RemoveFromListTree(TGListTree* ltree,
					     TGListTreeItem* parent_lti);

  virtual sLTI_i          FindItem(TGListTree* ltree);
  virtual sLTI_i          FindItem(TGListTree* ltree,
				   TGListTreeItem* parent_lti);
  virtual TGListTreeItem* FindListTreeItem(TGListTree* ltree);
  virtual TGListTreeItem* FindListTreeItem(TGListTree* ltree,
					   TGListTreeItem* parent_lti);

  virtual void UpdateItems();

  void SpawnEditor();                          // *MENU*
  virtual void ExportToCINT(Text_t* var_name); // *MENU*

  virtual void Destroy();                      // *MENU*

  virtual Bool_t CanEditRnrElement()   { return kTRUE; }
  virtual Bool_t GetRnrElement() const { return fRnrElement; }
  virtual void   SetRnrElement(Bool_t rnr);

  virtual Bool_t CanEditMainColor()        { return kFALSE; }
  Color_t* GetMainColorPtr()               { return fMainColorPtr; }
  void     SetMainColorPtr(Color_t* color) { fMainColorPtr = color; }

  virtual Color_t GetMainColor() const { return fMainColorPtr ? *fMainColorPtr : 0; }
  virtual void    SetMainColor(Color_t color);
          void    SetMainColor(Pixel_t pixel);

  ClassDef(RenderElement, 1);
}; // endclass RenderElement

/**************************************************************************/

class RenderElementObjPtr : public RenderElement,
                            public TObject
{
  RenderElementObjPtr(const RenderElementObjPtr&);            // Not implemented
  RenderElementObjPtr& operator=(const RenderElementObjPtr&); // Not implemented

protected:
  TObject* fObject;
  Bool_t   fOwnObject;

public:
  RenderElementObjPtr(TObject* obj, Bool_t own=kTRUE);
  RenderElementObjPtr(TObject* obj, Color_t& mainColor, Bool_t own=kTRUE);
  virtual ~RenderElementObjPtr();

  virtual TObject* GetObject(Exc_t eh="RenderElementObjPtr::GetObject ");
  virtual void     ExportToCINT(Text_t* var_name);

  Bool_t GetOwnObject() const   { return fOwnObject; }
  void   SetOwnObject(Bool_t o) { fOwnObject = o; }

  ClassDef(RenderElementObjPtr, 1);
}; // endclass RenderElementObjPtr

/**************************************************************************/

class RenderElementListBase : public RenderElement
{
protected:
  lpRE_t fChildren;

  void PaintElements(Option_t* option="");

public:
  RenderElementListBase() : RenderElement(), fChildren() {}
  RenderElementListBase(Color_t& col) : RenderElement(col), fChildren() {}
  virtual ~RenderElementListBase();

  virtual void AddElement(RenderElement* el);
  virtual void RemoveElement(RenderElement* el);
  virtual void RemoveElementLocal(RenderElement* el);
  virtual void RemoveElements();

  virtual void DestroyElements();

  virtual Int_t ExpandIntoListTree(TGListTree* ltree, TGListTreeItem* parent);
  virtual Int_t DestroyListSubTree(TGListTree* ltree, TGListTreeItem* parent);

  void EnableListElements();   // *MENU*
  void DisableListElements();  // *MENU*

  virtual void SetMainColor(Color_t color);
  virtual void SetMainColor(Pixel_t pixel);

  ClassDef(RenderElementListBase, 1);
};

/**************************************************************************/

class RenderElementList : public RenderElementListBase,
                          public TNamed
{
protected:
  Color_t   fColor;
  Bool_t    fDoColor;

public:
  RenderElementList(const Text_t* n="RenderElementList", const Text_t* t="",
		    Bool_t doColor=kFALSE);
  virtual ~RenderElementList() {}

  virtual Bool_t CanEditMainColor()  { return fDoColor; }

  virtual void Paint(Option_t* option = "") { PaintElements(option); }

  ClassDef(RenderElementList, 1);
};


/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

class ReferenceBackPtr : public ReferenceCount
{
protected:
  std::list<RenderElement*> fBackRefs;

public:
  ReferenceBackPtr();
  virtual ~ReferenceBackPtr();

  ReferenceBackPtr(const ReferenceBackPtr&);
  ReferenceBackPtr& operator=(const ReferenceBackPtr&);

  void IncRefCount(RenderElement* re);
  void DecRefCount(RenderElement* re);

  ClassDef(ReferenceBackPtr, 0);
};

}

#endif
