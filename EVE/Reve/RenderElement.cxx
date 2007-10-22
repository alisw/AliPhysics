// $Header$

#include "RenderElement.h"
#include "ReveManager.h"
#include "RGEditor.h"

#include <TColor.h>
#include <TCanvas.h>
#include <TGListTree.h>
#include <TGPicture.h>

#include <algorithm>

using namespace Reve;

//______________________________________________________________________
// RenderElement
//
//

ClassImp(RenderElement)

const TGPicture* RenderElement::fgRnrIcons[4] = { 0 };
const TGPicture* RenderElement::fgListTreeIcons[8] = { 0 };

RenderElement::RenderElement() :
  fRnrSelf             (kTRUE),
  fRnrChildren         (kTRUE),
  fMainColorPtr        (0),
  fItems               (),
  fParents             (),
  fDestroyOnZeroRefCnt (kTRUE),
  fDenyDestroy         (0),
  fChildren            ()
{}

RenderElement::RenderElement(Color_t& main_color) :
  fRnrSelf             (kTRUE),
  fRnrChildren         (kTRUE),
  fMainColorPtr        (&main_color),
  fItems               (),
  fParents             (),
  fDestroyOnZeroRefCnt (kTRUE),
  fDenyDestroy         (0),
  fChildren            ()
{}

RenderElement::~RenderElement()
{
  static const Exc_t _eh("RenderElement::RenderElement ");

  RemoveElements();

  for(List_i p=fParents.begin(); p!=fParents.end(); ++p)
  {
    (*p)->RemoveElementLocal(this);
    (*p)->fChildren.remove(this);
  }
  fParents.clear();

  for(sLTI_i i=fItems.begin(); i!=fItems.end(); ++i) {
    i->fTree->DeleteItem(i->fItem);
  }
}

/**************************************************************************/

void RenderElement::SetRnrElNameTitle(const Text_t* name, const Text_t* title)
{
  TNamed* named = dynamic_cast<TNamed*>(GetObject());
  if (named)
    named->SetNameTitle(name, title);
}

/**************************************************************************/
const Text_t* RenderElement::GetRnrElName() const
{
  TObject* named = dynamic_cast<TObject*>(GetObject());
  return named ? named->GetName() : "<noo-nom>";
}

/**************************************************************************/
const Text_t*  RenderElement::GetRnrElTitle() const
{
  TObject* named = dynamic_cast<TObject*>(GetObject());
  return named ? named->GetTitle() : "<noo-tit>";
}
/**************************************************************************/

void RenderElement::AddParent(RenderElement* re)
{
  // Adding parent is subordinate to adding an element.
  // Consider this an internal function.

  fParents.push_back(re);
}

void RenderElement::RemoveParent(RenderElement* re)
{
  // Removing parent is subordinate to removing an element.
  // Consider this an internal function.

  static const Exc_t eH("RenderElement::RemoveParent ");

  fParents.remove(re);
  CheckReferenceCount(eH);
}

void RenderElement::CheckReferenceCount(const Reve::Exc_t& eh)
{
  if(fParents.empty()   &&  fItems.empty()         &&
     fDenyDestroy <= 0  &&  fDestroyOnZeroRefCnt)
  {
    if (gDebug > 0)
       Info(eh, Form("auto-destructing '%s' on zero reference count.", GetRnrElName()));

    gReve->PreDeleteRenderElement(this);
    delete this;
  }
}

void RenderElement::CollectSceneParents(List_t& scenes)
{
  // Collect all parents of class Reve::Scene. This is needed to
  // automatically detect which scenes need to be updated.
  //
  // Overriden in Reve::Scene to include itself and return.

  for(List_i p=fParents.begin(); p!=fParents.end(); ++p)
    (*p)->CollectSceneParents(scenes);
}

void RenderElement::CollectSceneParentsFromChildren(List_t& scenes, RenderElement* parent)
{
  // Collect scene-parents from all children. This is needed to
  // automatically detect which scenes need to be updated during/after
  // a full sub-tree update.
  // Argument parent specifies parent in traversed hierarchy for which we can
  // skip the upwards search.

  for(List_i p=fParents.begin(); p!=fParents.end(); ++p)
  {
    if (*p != parent) (*p)->CollectSceneParents(scenes);
  }

  for(List_i c=fChildren.begin(); c!=fChildren.end(); ++c)
  {
    (*c)->CollectSceneParentsFromChildren(scenes, this);
  }
}

/**************************************************************************/
// List-tree stuff
/**************************************************************************/

Int_t RenderElement::ExpandIntoListTree(TGListTree* ltree,
					TGListTreeItem* parent)
{
  // Populates parent with elements.
  // parent must be an already existing representation of *this*.
  // Returns number of inserted elements.
  // If parent already has children, it does nothing.
  //
  // RnrEl can be inserted in a list-tree several times, thus we can not
  // search through fItems to get parent here.
  // Anyhow, it is probably known as it must have been selected by the user.

  if(parent->GetFirstChild() != 0)
    return 0;
  Int_t n = 0;
  for(List_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
    (*i)->AddIntoListTree(ltree, parent);
    ++n;
  }
  return n;
}

Int_t RenderElement::DestroyListSubTree(TGListTree* ltree,
					TGListTreeItem* parent)
{
  Int_t n = 0;
  TGListTreeItem* i = parent->GetFirstChild();
  while (i != 0)
  {
    //n += DestroyListSubTree(ltree, i);
    RenderElement* re = (RenderElement*) i->GetUserData();
    i = i->GetNextSibling();
    re->RemoveFromListTree(ltree, parent);
  }
  return n;
}

/**************************************************************************/

TGListTreeItem* RenderElement::AddIntoListTree(TGListTree* ltree,
					       TGListTreeItem* parent_lti)
{
  // Add this render element into ltree to already existing item
  // parent_lti.

  static const Exc_t eH("RenderElement::AddIntoListTree ");

  TObject* tobj = GetObject(eH);
  TGListTreeItem* item = ltree->AddItem(parent_lti, tobj->GetName(), this,
					0, 0, kTRUE);
  item->SetCheckBoxPictures(GetCheckBoxPicture(1, fRnrChildren),
			    GetCheckBoxPicture(0, fRnrChildren));

  item->SetPictures(GetListTreeIcon(),GetListTreeIcon());
  item->CheckItem(fRnrSelf);

  if(fMainColorPtr != 0) item->SetColor(GetMainColor());
  item->SetTipText(tobj->GetTitle());

  fItems.insert(ListTreeInfo(ltree, item));
  ltree->ClearViewPort();

  return item;
}

TGListTreeItem* RenderElement::AddIntoListTree(TGListTree* ltree,
					       RenderElement* parent)
{
  // Add this render element into ltree to all items belonging to
  // parent. Returns list-tree-item from the first register entry (but
  // we use a set for that so it can be anything).

  TGListTreeItem* lti = 0;
  if (parent == 0) {
    lti = AddIntoListTree(ltree, (TGListTreeItem*) 0);
  } else {
    for(sLTI_ri i = parent->fItems.rbegin(); i != parent->fItems.rend(); ++i)
    {
      if(i->fTree == ltree)
	lti = AddIntoListTree(ltree, i->fItem);
    }
  }
  return lti;
}

TGListTreeItem* RenderElement::AddIntoListTrees(RenderElement* parent)
{
  // Add this render element into all list-trees and all items
  // belonging to parent. Returns list-tree-item from the first
  // register entry (but we use a set for that so it can be anything).

  TGListTreeItem* lti = 0;
  for(sLTI_ri i = parent->fItems.rbegin(); i != parent->fItems.rend(); ++i)
  {
    lti = AddIntoListTree(i->fTree, i->fItem);
  }
  return lti;
}

Bool_t RenderElement::RemoveFromListTree(TGListTree* ltree,
					 TGListTreeItem* parent_lti)
{
  static const Exc_t eH("RenderElement::RemoveFromListTree ");

  sLTI_i i = FindItem(ltree, parent_lti);
  if(i != fItems.end()) {
    DestroyListSubTree(ltree, i->fItem);
    ltree->DeleteItem(i->fItem);
    ltree->ClearViewPort();
    fItems.erase(i);
    if (parent_lti == 0) CheckReferenceCount(eH);
    return kTRUE;
  } else {
    return kFALSE;
  }
}

Int_t RenderElement::RemoveFromListTrees(RenderElement* parent)
{
  Int_t count = 0;

  sLTI_i i  = fItems.begin();
  while (i != fItems.end())
  {
    sLTI_i j = i++;
    TGListTreeItem *plti = j->fItem->GetParent();
    if (plti != 0 && (RenderElement*) plti->GetUserData() == parent)
    {
      DestroyListSubTree(j->fTree, j->fItem);
      j->fTree->DeleteItem(j->fItem);
      j->fTree->ClearViewPort();
      fItems.erase(j);
      ++count;
    }
  }

  return count;
}

/**************************************************************************/

RenderElement::sLTI_i RenderElement::FindItem(TGListTree* ltree)
{
  for(sLTI_i i = fItems.begin(); i != fItems.end(); ++i)
    if(i->fTree == ltree)
      return i;
  return fItems.end();
}

RenderElement::sLTI_i RenderElement::FindItem(TGListTree* ltree,
					      TGListTreeItem* parent_lti)
{
  for(sLTI_i i = fItems.begin(); i != fItems.end(); ++i)
    if(i->fTree == ltree && i->fItem->GetParent() == parent_lti)
      return i;
  return fItems.end();
}

TGListTreeItem* RenderElement::FindListTreeItem(TGListTree* ltree)
{
  for(sLTI_i i = fItems.begin(); i != fItems.end(); ++i)
    if(i->fTree == ltree)
      return i->fItem;
  return 0;
}

TGListTreeItem* RenderElement::FindListTreeItem(TGListTree* ltree,
						TGListTreeItem* parent_lti)
{
  for(sLTI_i i = fItems.begin(); i != fItems.end(); ++i)
    if(i->fTree == ltree && i->fItem->GetParent() == parent_lti)
      return i->fItem;
  return 0;
}

/**************************************************************************/

void RenderElement::UpdateItems()
{
  static const Exc_t eH("RenderElement::UpdateItems ");

  TObject* tobj = GetObject(eH);

  for(sLTI_i i=fItems.begin(); i!=fItems.end(); ++i) {
    i->fItem->Rename(tobj->GetName());
    i->fItem->SetTipText(tobj->GetTitle());
    i->fItem->CheckItem(fRnrSelf);
    if(fMainColorPtr != 0) i->fItem->SetColor(GetMainColor());
    i->fTree->ClearViewPort();
  }
}

/**************************************************************************/

TObject* RenderElement::GetObject(Exc_t eh) const
{
  const TObject* obj = dynamic_cast<const TObject*>(this);
  if(obj == 0)
    throw(eh + "not a TObject.");
  return const_cast<TObject*>(obj);
}

void RenderElement::SpawnEditor()
{
  gReve->EditRenderElement(this);
}

void RenderElement::ExportToCINT(Text_t* var_name)
{
  static const Exc_t eH("RenderElement::ExportToCINT ");

  const char* cname = IsA()->GetName();
  gROOT->ProcessLine(Form("%s* %s = (%s*)0x%lx;", cname, var_name, cname, this));
}

/**************************************************************************/

void RenderElement::PadPaint(Option_t* option)
{
  if (GetRnrSelf() && GetObject()) 
    GetObject()->Paint(option);
  

  if (GetRnrChildren()) {
    for(List_i i=BeginChildren(); i!=EndChildren(); ++i) {
      (*i)->PadPaint(option);
    }
  }
}

/**************************************************************************/

void RenderElement::SetRnrSelf(Bool_t rnr)
{
  if(rnr != fRnrSelf)
  {
    fRnrSelf = rnr;
    
    for(sLTI_i i=fItems.begin(); i!=fItems.end(); ++i) 
    {
      if(i->fItem->IsChecked() != rnr) {
        i->fItem->SetCheckBoxPictures(GetCheckBoxPicture(1, fRnrChildren),
				      GetCheckBoxPicture(0, fRnrChildren));
        i->fItem->CheckItem(fRnrSelf);
        i->fTree->ClearViewPort();
      }
    }
  }
}

void RenderElement::SetRnrChildren(Bool_t rnr)
{
  if(rnr != fRnrChildren)
  {
    fRnrChildren = rnr; 

    for(sLTI_i i=fItems.begin(); i!=fItems.end(); ++i) 
    {
      i->fItem->SetCheckBoxPictures(GetCheckBoxPicture(fRnrSelf, fRnrChildren),
				    GetCheckBoxPicture(fRnrSelf, fRnrChildren));
      i->fTree->ClearViewPort();
    }
  }
}

void RenderElement::SetRnrState(Bool_t rnr)
{
  if (fRnrSelf != rnr || fRnrChildren != rnr)
  {
    fRnrSelf = fRnrChildren = rnr; 

    for (sLTI_i i=fItems.begin(); i!=fItems.end(); ++i) 
    {
      i->fItem->SetCheckBoxPictures(GetCheckBoxPicture(1,1), GetCheckBoxPicture(0,0));
      i->fItem->CheckItem(fRnrSelf);
      i->fTree->ClearViewPort();
    }
  }
}

/**************************************************************************/

void RenderElement::SetMainColor(Color_t color)
{
  Color_t oldcol = GetMainColor();
  for(List_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
    if((*i)->GetMainColor() == oldcol) (*i)->SetMainColor(color);
  }

  if (fMainColorPtr) {
    *fMainColorPtr = color;
    for(sLTI_i i=fItems.begin(); i!=fItems.end(); ++i) {
      if(i->fItem->GetColor() != color) {
        i->fItem->SetColor(GetMainColor());
	i->fTree->ClearViewPort();
      }
    }
  }
}

void RenderElement::SetMainColor(Pixel_t pixel)
{
  SetMainColor(Color_t(TColor::GetColor(pixel)));
}

/**************************************************************************/

TGListTreeItem* RenderElement::AddElement(RenderElement* el)
{
  el->AddParent(this);
  fChildren.push_back(el);
  TGListTreeItem* ret = el->AddIntoListTrees(this);
  ElementChanged();
  return ret;
}

void RenderElement::RemoveElement(RenderElement* el)
{
  // Remove el from the list of children.

  el->RemoveFromListTrees(this);
  RemoveElementLocal(el);
  el->RemoveParent(this);
  fChildren.remove(el);
  ElementChanged();
}

void RenderElement::RemoveElementLocal(RenderElement* /*el*/)
{
  // Perform additional local removal of el.
  // Called from RemoveElement() which does whole untangling.
  // Put into special function as framework-related handling of
  // element removal should really be common to all classes and
  // clearing of local structures happens in between removal
  // of list-tree-items and final untangling.
  // If you override this, you should also override
  // RemoveElementsLocal().
}

void RenderElement::RemoveElements()
{
  // Remove all elements. This assumes removing of all elements can be
  // done more efficiently then looping over them and removing one by
  // one.

  for(sLTI_i i=fItems.begin(); i!=fItems.end(); ++i)
  {
    DestroyListSubTree(i->fTree, i->fItem);
  }
  RemoveElementsLocal();
  for(List_i i=fChildren.begin(); i!=fChildren.end(); ++i)
  {
    (*i)->RemoveParent(this);
  }
  fChildren.clear();
  ElementChanged();
}

void RenderElement::RemoveElementsLocal()
{
  // Perform additional local removal of all elements.
  // See comment to RemoveelementLocal(RenderElement*).
}

/**************************************************************************/

void RenderElement::EnableListElements(Bool_t rnr_self,  Bool_t rnr_children)
{
  for(List_i i=fChildren.begin(); i!=fChildren.end(); ++i)
  {
    (*i)->SetRnrSelf(rnr_self);
    (*i)->SetRnrChildren(rnr_children);
  }

  ElementChanged(kTRUE, kTRUE);
}

void RenderElement::DisableListElements(Bool_t rnr_self,  Bool_t rnr_children)
{
  for(List_i i=fChildren.begin(); i!=fChildren.end(); ++i)
  {
    (*i)->SetRnrSelf(rnr_self);
    (*i)->SetRnrChildren(rnr_children);
  }

  ElementChanged(kTRUE, kTRUE);
}

/**************************************************************************/

void RenderElement::Destroy()
{
  static const Exc_t eH("RenderElement::Destroy ");

  if (fDenyDestroy > 0)
    throw(eH + "this element '%s' is protected against destruction.", GetRnrElName());

  gReve->PreDeleteRenderElement(this);
  delete this;
  gReve->Redraw3D();
}

void RenderElement::DestroyElements()
{
  static const Exc_t eH("RenderElement::DestroyElements ");

  while( ! fChildren.empty()) {
    RenderElement* c = fChildren.front();
    if (c->fDenyDestroy <= 0)
    {
      try {
        c->Destroy();
      }
      catch(Exc_t exc) {
        Warning(eH, Form("element destruction failed: '%s'.", exc.Data()));
        RemoveElement(c);
      }
    }
    else
    {
      if (gDebug > 0)
        Info(eH, Form("element '%s' is protected agains destruction, removin locally.", c->GetRnrElName()));

      RemoveElement(c);
    }
  }
}

/**************************************************************************/

Bool_t RenderElement::HandleElementPaste(RenderElement* el)
{
  // React to element being pasted or dnd-ed.
  // Return true if redraw is needed.

  gReve->AddRenderElement(el, this);
  return kTRUE;
}

void RenderElement::ElementChanged(Bool_t update_scenes, Bool_t redraw)
{
  if (update_scenes)
    gReve->RenderElementChanged(this);
  if (redraw)
    gReve->Redraw3D();
}

/**************************************************************************/
// Statics
/**************************************************************************/

const TGPicture*
RenderElement::GetCheckBoxPicture(Bool_t rnrSelf, Bool_t rnrDaughters)
{
  Int_t idx = 0;
  if(rnrSelf)       idx = 2;
  if(rnrDaughters ) idx++;

  return fgRnrIcons[idx];
}


/**************************************************************************/
/**************************************************************************/

//______________________________________________________________________
// Reve::RenderElementObjPtr
//

ClassImp(RenderElementObjPtr)

RenderElementObjPtr::RenderElementObjPtr(TObject* obj, Bool_t own) :
  RenderElement(),
  fObject(obj),
  fOwnObject(own)
{}

RenderElementObjPtr::RenderElementObjPtr(TObject* obj, Color_t& mainColor, Bool_t own) :
  RenderElement(mainColor),
  fObject(obj),
  fOwnObject(own)
{}

TObject* RenderElementObjPtr::GetObject(Reve::Exc_t eh) const
{
  if(fObject == 0)
    throw(eh + "fObject not set.");
  return fObject;
}

void RenderElementObjPtr::ExportToCINT(Text_t* var_name)
{
  static const Exc_t eH("RenderElementObjPtr::ExportToCINT ");

  TObject* obj = GetObject(eH);
  const char* cname = obj->IsA()->GetName();
  gROOT->ProcessLine(Form("%s* %s = (%s*)0x%lx;", cname, var_name, cname, obj));
}

RenderElementObjPtr::~RenderElementObjPtr()
{
  if(fOwnObject)
    delete fObject;
}


/**************************************************************************/
/**************************************************************************/

//______________________________________________________________________
// Reve::RenderElementList
//

// !!! should have two ctors (like in RenderElement), one with Color_t&
// and set fDoColor automatically, based on which ctor is called.

ClassImp(RenderElementList)

RenderElementList::RenderElementList(const Text_t* n, const Text_t* t, Bool_t doColor) :
  RenderElement(),
  TNamed(n, t),
  fColor(0),
  fDoColor(doColor),
  fChildClass(0)
{
  if(fDoColor) {
    SetMainColorPtr(&fColor);
  }
}

Bool_t RenderElementList::AcceptRenderElement(RenderElement* el)
{
  if (fChildClass && ! el->IsA()->InheritsFrom(fChildClass))
    return kFALSE;
  return kTRUE;
}
