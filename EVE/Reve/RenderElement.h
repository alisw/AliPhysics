// $Header$

#ifndef REVE_RenderElement_H
#define REVE_RenderElement_H

#include <Reve/Reve.h>

#include <TNamed.h>
#include <TRef.h>

class TGListTree;
class TGListTreeItem;
class TGPicture;

namespace Reve {

class ZTrans;

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

  static const TGPicture*                      fgRnrIcons[4];

  typedef std::set<ListTreeInfo>               sLTI_t;
  typedef std::set<ListTreeInfo>::iterator     sLTI_i;

  typedef std::list<RenderElement*>            List_t;
  typedef std::list<RenderElement*>::iterator  List_i;

protected:
  // TRef     fSource;

  Bool_t   fRnrSelf;
  Bool_t   fRnrChildren;
  Color_t* fMainColorPtr;

  sLTI_t fItems;
  List_t fParents;

  Bool_t fDestroyOnZeroRefCnt;
  Bool_t fDenyDestroy;

  List_t fChildren;

public:
  RenderElement();
  RenderElement(Color_t& main_color);
  virtual ~RenderElement();

  virtual void AddParent(RenderElement* re);
  virtual void RemoveParent(RenderElement* re);

  List_i BeginParents() { return fParents.begin(); }
  List_i EndParents()   { return fParents.end();   }
  Int_t  GetNParents() const { return fParents.size(); }

  List_i BeginChildren() { return fChildren.begin(); }
  List_i EndChildren()   { return fChildren.end();   }
  Int_t  GetNChildren() const { return fChildren.size(); }

  virtual Int_t ExpandIntoListTree(TGListTree* ltree, TGListTreeItem* parent);
  virtual Int_t DestroyListSubTree(TGListTree* ltree, TGListTreeItem* parent);

  void EnableListElements();   // *MENU*
  void DisableListElements();  // *MENU*

  Bool_t GetDestroyOnZeroRefCnt() const   { return fDestroyOnZeroRefCnt; }
  void   SetDestroyOnZeroRefCnt(Bool_t d) { fDestroyOnZeroRefCnt = d; }

  Bool_t GetDenyDestroy() const   { return fDenyDestroy; }
  void   SetDenyDestroy(Bool_t d) { fDenyDestroy = d; }

  virtual void PadPaint(Option_t* option);

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


  virtual void AddElement(RenderElement* el);
  virtual void RemoveElement(RenderElement* el);
  virtual void RemoveElementLocal(RenderElement* el);
  virtual void RemoveElements();

  virtual void Destroy();                      // *MENU*
  virtual void DestroyElements();

  virtual Bool_t CanEditRnrElement()   { return kTRUE; }
  virtual Bool_t GetRnrSelf() const { return fRnrSelf; }
  virtual Bool_t GetRnrChildren() const { return fRnrChildren; }
  virtual void   SetRnrSelf(Bool_t rnr);
  virtual void   SetRnrChildren(Bool_t rnr);
  virtual void   SetRnrState(Bool_t rnr);

  virtual Bool_t CanEditMainColor()        { return kFALSE; }
  Color_t* GetMainColorPtr()               { return fMainColorPtr; }
  void     SetMainColorPtr(Color_t* color) { fMainColorPtr = color; }

  virtual Color_t GetMainColor() const { return fMainColorPtr ? *fMainColorPtr : 0; }
  virtual void    SetMainColor(Color_t color);
  void    SetMainColor(Pixel_t pixel);

  virtual Bool_t  CanEditMainTransparency()    { return kFALSE; }
  virtual UChar_t GetMainTransparency() const  { return 0; }
  virtual void    SetMainTransparency(UChar_t) {}

  virtual Bool_t  CanEditMainHMTrans() { return kFALSE; }
  virtual ZTrans* PtrMainHMTrans()     { return 0; }

  static const TGPicture* GetCheckBoxPicture(Bool_t rnrElement, Bool_t rnrDaughter);

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

class RenderElementList : public RenderElement,
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

  ClassDef(RenderElementList, 1);
};


/**************************************************************************/

class PadPrimitive : public RenderElement,
		     public TNamed
{
public:
  PadPrimitive(const Text_t* n="PadPrimitive", const Text_t* t="");
  virtual ~PadPrimitive() {}

  virtual void Paint(Option_t* option = "");

  // virtual void PaintRenderElement(RenderElement* re, Option_t* option);

  ClassDef(PadPrimitive, 1);
};


}

#endif
