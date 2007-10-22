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
  friend class ReveManager;

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
  static const TGPicture*                      fgListTreeIcons[8];

  typedef std::set<ListTreeInfo>               sLTI_t;
  typedef sLTI_t::iterator                     sLTI_i;
  typedef sLTI_t::reverse_iterator             sLTI_ri;

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
  Int_t  fDenyDestroy;

  List_t fChildren;

public:
  RenderElement();
  RenderElement(Color_t& main_color);
  virtual ~RenderElement();

  virtual void SetRnrElNameTitle(const Text_t* name, const Text_t* title);
  virtual const Text_t* GetRnrElName()  const;
  virtual const Text_t* GetRnrElTitle() const;

  virtual void AddParent(RenderElement* re);
  virtual void RemoveParent(RenderElement* re);
  virtual void CheckReferenceCount(const Reve::Exc_t& eh="RenderElement::CheckReferenceCount ");
  virtual void CollectSceneParents(List_t& scenes);
  virtual void CollectSceneParentsFromChildren(List_t& scenes, RenderElement* parent);

  List_i BeginParents() { return fParents.begin(); }
  List_i EndParents()   { return fParents.end();   }
  Int_t  GetNParents() const { return fParents.size(); }

  List_i BeginChildren() { return fChildren.begin(); }
  List_i EndChildren()   { return fChildren.end();   }
  Int_t  GetNChildren() const { return fChildren.size(); }

  void EnableListElements (Bool_t rnr_self=kTRUE,  Bool_t rnr_children=kTRUE);  // *MENU*
  void DisableListElements(Bool_t rnr_self=kFALSE, Bool_t rnr_children=kFALSE); // *MENU*

  Bool_t GetDestroyOnZeroRefCnt() const   { return fDestroyOnZeroRefCnt; }
  void   SetDestroyOnZeroRefCnt(Bool_t d) { fDestroyOnZeroRefCnt = d; }

  Int_t  GetDenyDestroy() const { return fDenyDestroy; }
  void   IncDenyDestroy()       { ++fDenyDestroy; }
  void   DecDenyDestroy()       { if (--fDenyDestroy <= 0) CheckReferenceCount("RenderElement::DecDenyDestroy "); }

  virtual void PadPaint(Option_t* option);

  virtual TObject* GetObject(Reve::Exc_t eh="RenderElement::GetObject ") const;
  virtual TObject* GetEditorObject() const { return GetObject(); }
  /*
    TRef&    GetSource() { return fSource; }
    TObject* GetSourceObject() const { return fSource.GetObject(); }
    void SetSourceObject(TObject* o) { fSource.SetObject(o); }

    void DumpSourceObject();    // *MENU*
    void InspectSourceObject(); // *MENU*
  */

  // --------------------------------

  virtual Int_t ExpandIntoListTree(TGListTree* ltree, TGListTreeItem* parent);
  virtual Int_t DestroyListSubTree(TGListTree* ltree, TGListTreeItem* parent);

  virtual TGListTreeItem* AddIntoListTree(TGListTree* ltree,
					  TGListTreeItem* parent_lti);
  virtual TGListTreeItem* AddIntoListTree(TGListTree* ltree,
					  RenderElement* parent);
  virtual TGListTreeItem* AddIntoListTrees(RenderElement* parent);

  virtual Bool_t          RemoveFromListTree(TGListTree* ltree,
					     TGListTreeItem* parent_lti);
  virtual Int_t           RemoveFromListTrees(RenderElement* parent);

  virtual sLTI_i          FindItem(TGListTree* ltree);
  virtual sLTI_i          FindItem(TGListTree* ltree,
				   TGListTreeItem* parent_lti);
  virtual TGListTreeItem* FindListTreeItem(TGListTree* ltree);
  virtual TGListTreeItem* FindListTreeItem(TGListTree* ltree,
					   TGListTreeItem* parent_lti);

  virtual Int_t GetNItems() const { return fItems.size(); }
  virtual void  UpdateItems();

  void SpawnEditor();                          // *MENU*
  virtual void ExportToCINT(Text_t* var_name); // *MENU*

  virtual Bool_t AcceptRenderElement(RenderElement* /*el*/) { return kTRUE; }

  virtual TGListTreeItem* AddElement(RenderElement* el);
  virtual void RemoveElement(RenderElement* el);
  virtual void RemoveElementLocal(RenderElement* el);
  virtual void RemoveElements();
  virtual void RemoveElementsLocal();

  virtual void Destroy();                      // *MENU*
  virtual void DestroyElements();              // *MENU*

  virtual Bool_t HandleElementPaste(RenderElement* el);
  virtual void   ElementChanged(Bool_t update_scenes=kTRUE, Bool_t redraw=kFALSE);

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

  static  const TGPicture* GetCheckBoxPicture(Bool_t rnrElement, Bool_t rnrDaughter);
  virtual const TGPicture* GetListTreeIcon() { return fgListTreeIcons[0]; }

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

  virtual TObject* GetObject(Exc_t eh="RenderElementObjPtr::GetObject ") const;
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
  TClass   *fChildClass;

public:
  RenderElementList(const Text_t* n="RenderElementList", const Text_t* t="",
		    Bool_t doColor=kFALSE);
  virtual ~RenderElementList() {}

  virtual Bool_t CanEditMainColor()  { return fDoColor; }

  TClass* GetChildClass() const { return fChildClass; }
  void SetChildClass(TClass* c) { fChildClass = c; }

  virtual Bool_t AcceptRenderElement(RenderElement* el);

  ClassDef(RenderElementList, 1);
};

}

#endif
