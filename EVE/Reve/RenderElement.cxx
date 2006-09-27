// $Header$

#include "RenderElement.h"
#include "RGTopFrame.h"

#include <TColor.h>
#include <TCanvas.h>
#include <TGListTree.h>
#include <THashList.h>

using namespace Reve;

//______________________________________________________________________
// RenderElement
//
//

ClassImp(RenderElement)

RenderElement::RenderElement() :
  fRnrElement          (kTRUE),
  fMainColorPtr        (0),
  fItems               (),
  fParents             (),
  fDestroyOnZeroRefCnt (kTRUE),
  fDenyDestroy         (kFALSE)
{}

RenderElement::RenderElement(Color_t& main_color) :
  fRnrElement          (kTRUE),
  fMainColorPtr        (&main_color),
  fItems               (),
  fParents             (),
  fDestroyOnZeroRefCnt (kTRUE),
  fDenyDestroy         (kFALSE)
{}

RenderElement::~RenderElement()
{
  static const Exc_t _eh("RenderElement::RenderElement ");

  for(lpRE_i p=fParents.begin(); p!=fParents.end(); ++p) {
    RenderElementListBase* l = dynamic_cast<RenderElementListBase*>(*p);
    if(l)
      l->RemoveElementLocal(this);
  }
  fParents.clear();

  for(sLTI_i i=fItems.begin(); i!=fItems.end(); ++i) {
    i->fTree->DeleteItem(i->fItem);
    gClient->NeedRedraw(i->fTree);
  }
}

void RenderElement::Destroy()
{
  static const Exc_t eH("RenderElement::Destroy ");

  if(fDenyDestroy)
    throw(eH + "this object is protected against destruction.");

  delete this;
  gReve->Redraw3D();
}

/**************************************************************************/

void RenderElement::AddParent(RenderElement* re)
{
  fParents.push_back(re);
}

void RenderElement::RemoveParent(RenderElement* re)
{
  static const Exc_t eH("RenderElement::RemoveParent ");

  fParents.remove(re);
  if(fParents.empty() && fDestroyOnZeroRefCnt) {
    // TObject* tobj = GetObject(eH);
    // Warning(eH.Data(), Form("auto-destructing '%s' on zero reference count.", tobj->GetName()));
    delete this;
  }
}

/**************************************************************************/

TObject* RenderElement::GetObject(Exc_t eh)
{
  TObject* obj = dynamic_cast<TObject*>(this);
  if(obj == 0)
    throw(eh + "not a TObject.");
  return obj;
}

/**************************************************************************/

TGListTreeItem* RenderElement::AddIntoListTree(TGListTree* ltree,
					       TGListTreeItem* parent_lti)
{
  static const Exc_t eH("RenderElement::AddIntoListTree ");

  TObject* tobj = GetObject(eH);
  TGListTreeItem* item = ltree->AddItem(parent_lti, tobj->GetName(), this,
					0, 0, kTRUE);
  gClient->NeedRedraw(ltree);
  item->CheckItem(GetRnrElement());
  
  if(fMainColorPtr != 0) item->SetColor(GetMainColor());
  item->SetTipText(tobj->GetTitle());

  fItems.insert(ListTreeInfo(ltree, item));

  return item;
}

TGListTreeItem* RenderElement::AddIntoListTree(TGListTree* ltree,
					       RenderElement* parent)
{
  TGListTreeItem* parent_lti = parent ? parent->FindListTreeItem(ltree) : 0;
  return AddIntoListTree(ltree, parent_lti);
}

Bool_t RenderElement::RemoveFromListTree(TGListTree* ltree)
{
  sLTI_i i = FindItem(ltree);
  if(i != fItems.end()) {
    ltree->DeleteItem(i->fItem);
    fItems.erase(i);
    gClient->NeedRedraw(ltree);
    return kTRUE;
  } else {
    return kFALSE;
  }
}

Bool_t RenderElement::RemoveFromListTree(TGListTree* ltree,
					 TGListTreeItem* parent_lti)
{
  sLTI_i i = FindItem(ltree, parent_lti);
  if(i != fItems.end()) {
    ltree->DeleteItem(i->fItem);
    fItems.erase(i);
    gClient->NeedRedraw(ltree);
    return kTRUE;
  } else {
    return kFALSE;
  }
}

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
    i->fItem->CheckItem(fRnrElement);
    if(fMainColorPtr != 0) i->fItem->SetColor(GetMainColor());
    gClient->NeedRedraw(i->fTree);
  }
  gReve->Redraw3D(); // This will go away once editor can notify ALL changes.
}

/**************************************************************************/

void RenderElement::SpawnEditor()
{
  gReve->EditRenderElement(this);
}

void RenderElement::ExportToCINT(Text_t* var_name)
{
  static const Exc_t eH("RenderElement::ExportToCINT ");

  TObject* obj = GetObject(eH);
  const char* cname = obj->IsA()->GetName();
  gROOT->ProcessLine(Form("%s* %s = (%s*) %p;", cname, var_name, cname, obj));
}

/**************************************************************************/
/*
void RenderElement::DumpSourceObject()
{
  TObject* o = GetSourceObject();
  if(o == 0) {
    printf("Source object not set.\n");
  } else {
    o->Dump();
  }
}

void RenderElement::InspectSourceObject()
{
  TObject* o = GetSourceObject();
  if(o == 0) {
    printf("Source object not set.\n");
  } else {
    o->Inspect();
  }
}
*/
/**************************************************************************/

void RenderElement::SetRnrElement(Bool_t rnr)
{
  if(rnr != fRnrElement) {
    fRnrElement = rnr;
    UpdateItems();
  }
}

/**************************************************************************/

void RenderElement::SetMainColor(Color_t color)
{
  if (fMainColorPtr) {
    *fMainColorPtr = color;
    UpdateItems();
  }
}

void RenderElement::SetMainColor(Pixel_t pixel)
{
  SetMainColor(Color_t(TColor::GetColor(pixel)));
}

/**************************************************************************/
/**************************************************************************/

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

RenderElementObjPtr::~RenderElementObjPtr()
{
  if(fOwnObject)
    delete fObject;
}

/**************************************************************************/

TObject* RenderElementObjPtr::GetObject(Reve::Exc_t eh)
{
  if(fObject == 0)
    throw(eh + "fObject not set.");
  return fObject;
}

/**************************************************************************/
/**************************************************************************/

ClassImp(RenderElementListBase)

RenderElementListBase::~RenderElementListBase()
{
  // Collapse all sub-trees
  for(sLTI_i i=fItems.begin(); i!=fItems.end(); ++i) {
    DestroyListSubTree(i->fTree, i->fItem);
    // My own items removed in ~RenderElement
  }
  RemoveElements();
}

/**************************************************************************/

void RenderElementListBase::AddElement(RenderElement* el)
{
  fChildren.push_back(el);
  el->AddParent(this);
}

void RenderElementListBase::RemoveElement(RenderElement* el)
{
  el->RemoveParent(this);
  RemoveElementLocal(el);
}

void RenderElementListBase::RemoveElementLocal(RenderElement* el)
{
  fChildren.remove(el);
}

void RenderElementListBase::RemoveElements()
{
  for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
    (*i)->RemoveParent(this);
  }
  fChildren.clear();
}

/**************************************************************************/

void RenderElementListBase::DestroyElements()
{
  static const Exc_t eH("RenderElementListBase::DestroyElements ");

  while( ! fChildren.empty()) {
    RenderElement* c = fChildren.front();
    try {
      c->Destroy();
    }
    catch(Exc_t exc) {
      Warning(eH.Data(), Form("element destruction failed: '%s'.", exc.Data()));
      RemoveElement(c);
    }
  }
}

/**************************************************************************/

Int_t RenderElementListBase::ExpandIntoListTree(TGListTree* ltree,
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
  for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
    (*i)->AddIntoListTree(ltree, parent);
    ++n;
  }
  return n;
}

Int_t RenderElementListBase::DestroyListSubTree(TGListTree* ltree,
						TGListTreeItem* parent)
{
  Int_t n = 0;
  TGListTreeItem* i = parent->GetFirstChild();
  while(i != 0) {
    n += DestroyListSubTree(ltree, i);
    RenderElement* re = (RenderElement*) i->GetUserData();
    i = i->GetNextSibling();
    re->RemoveFromListTree(ltree, parent);
  }
  return n;
}

/**************************************************************************/

void RenderElementListBase::EnableListElements()
{
  for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i)
    (*i)->SetRnrElement(kTRUE);
}

void RenderElementListBase::DisableListElements()
{
  for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i)
    (*i)->SetRnrElement(kFALSE);
}

/**************************************************************************/

void RenderElementListBase::SetMainColor(Color_t col)
{
  Color_t oldcol = GetMainColor();
  for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
    if((*i)->GetMainColor() == oldcol) (*i)->SetMainColor(col);
  }
  RenderElement::SetMainColor(col);
}

void RenderElementListBase::SetMainColor(Pixel_t pixel)
{
  // This one needed for proper calling via CINT (signals).

  SetMainColor(Color_t(TColor::GetColor(pixel)));
}

/**************************************************************************/

void RenderElementListBase::PaintElements(Option_t* option)
{
  if(fRnrElement) {
    for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
      if((*i)->GetRnrElement())
	(*i)->GetObject()->Paint(option);
    }
  }
}

/**************************************************************************/
/**************************************************************************/

ClassImp(RenderElementList)

RenderElementList::RenderElementList(const Text_t* n, const Text_t* t, Bool_t doColor) :
  RenderElementListBase(),
  TNamed(n, t),
  fColor(0),
  fDoColor(doColor)
{
  if(fDoColor) {
    SetMainColorPtr(&fColor);
  }
}
