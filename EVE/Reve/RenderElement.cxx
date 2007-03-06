// $Header$

#include "RenderElement.h"
#include "RGTopFrame.h"
#include "RGEditor.h"

#include <TColor.h>
#include <TCanvas.h>
#include <TGListTree.h>
#include <TGPicture.h>
#include <THashList.h>

#include <algorithm>

using namespace Reve;

//______________________________________________________________________
// RenderElement
//
//

ClassImp(RenderElement)

const TGPicture* RenderElement::fgRnrIcons[4] = { 0 };

RenderElement::RenderElement() :
  fRnrSelf          (kTRUE),
  fRnrChildren        (kTRUE),
  fMainColorPtr        (0),
  fItems               (),
  fParents             (),
  fDestroyOnZeroRefCnt (kTRUE),
  fDenyDestroy         (kFALSE),
  fChildren(0)
{}

RenderElement::RenderElement(Color_t& main_color) :
  fRnrSelf          (kTRUE),
  fRnrChildren        (kTRUE),
  fMainColorPtr        (&main_color),
  fItems               (),
  fParents             (),
  fDestroyOnZeroRefCnt (kTRUE),
  fDenyDestroy         (kFALSE),
  fChildren(0)
{}

RenderElement::~RenderElement()
{
  static const Exc_t _eh("RenderElement::RenderElement ");

  for(sLTI_i i=fItems.begin(); i!=fItems.end(); ++i) {
    DestroyListSubTree(i->fTree, i->fItem);
  }
  RemoveElements();

  for(List_i p=fParents.begin(); p!=fParents.end(); ++p) {
    (*p)->RemoveElementLocal(this);
  }
  fParents.clear();

  for(sLTI_i i=fItems.begin(); i!=fItems.end(); ++i) {
    i->fTree->DeleteItem(i->fItem);
    gClient->NeedRedraw(i->fTree);
  }
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
    gReve->PreDeleteRenderElement(this);
    delete this;
  }
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
  item->CheckItem(GetRnrSelf());
  
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
    gClient->NeedRedraw(i->fTree);
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
  if ( GetRnrSelf() && GetObject() ) 
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
      if (i->fItem->IsChecked() != rnr) {
        i->fItem->SetCheckBoxPictures(GetCheckBoxPicture(1,fRnrChildren), GetCheckBoxPicture(0,fRnrChildren));
#if ROOT_VERSION_CODE <= ROOT_VERSION(5,14,0)
        // cover undesired behavior TGListTree::UpdateChecked
	if(i->fItem->GetParent()) {
	  RenderElement* p = (RenderElement*) i->fItem->GetParent()->GetUserData();
	  if(p) 
	    i->fItem->GetParent()->SetCheckBoxPictures(GetCheckBoxPicture(1, p->GetRnrChildren()),
						       GetCheckBoxPicture(0, p->GetRnrChildren()));
	}
#endif
        i->fItem->CheckItem(fRnrSelf);
        gClient->NeedRedraw(i->fTree);
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
      i->fItem->SetCheckBoxPictures(GetCheckBoxPicture(fRnrSelf, fRnrChildren), GetCheckBoxPicture(fRnrSelf,fRnrChildren));
      gClient->NeedRedraw(i->fTree);
    }
  }
}

void RenderElement::ToggleRnrState()
{
  if (fRnrSelf || fRnrChildren)
    fRnrSelf = fRnrChildren = kFALSE;  
  else
    fRnrSelf = fRnrChildren = kTRUE; 

  for (sLTI_i i=fItems.begin(); i!=fItems.end(); ++i) 
  {
    i->fItem->SetCheckBoxPictures(GetCheckBoxPicture(1,1), GetCheckBoxPicture(0,0));
    i->fItem->CheckItem(fRnrSelf);

#if ROOT_VERSION_CODE <= ROOT_VERSION(5,14,0)
    // cover undesired behavior TGListTree::UpdateChecked
    if(i->fItem->GetParent()) {
      RenderElement* p = (RenderElement*) i->fItem->GetParent()->GetUserData();
      if(p) 
	i->fItem->GetParent()->SetCheckBoxPictures(GetCheckBoxPicture(1, p->GetRnrChildren()),
						   GetCheckBoxPicture(0, p->GetRnrChildren()));
    }
#endif
    gClient->NeedRedraw(i->fTree);
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
        gClient->NeedRedraw(i->fTree);
      }
    }
  }
}

void RenderElement::SetMainColor(Pixel_t pixel)
{
  SetMainColor(Color_t(TColor::GetColor(pixel)));
}

/**************************************************************************/

void RenderElement::AddElement(RenderElement* el)
{
  fChildren.push_back(el);
  el->AddParent(this);
}

void RenderElement::RemoveElement(RenderElement* el)
{
   el->RemoveParent(this);
   RemoveElementLocal(el);
}

void RenderElement::RemoveElementLocal(RenderElement* el)
{
  fChildren.remove(el);
}

void RenderElement::RemoveElements()
{
  for(List_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
    (*i)->RemoveParent(this);
  }
  fChildren.clear();
}

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

/**************************************************************************/

void RenderElement::EnableListElements()
{
  for(List_i i=fChildren.begin(); i!=fChildren.end(); ++i)
    (*i)->SetRnrSelf(kTRUE);

  gReve->Redraw3D();
}

void RenderElement::DisableListElements()
{
  for(List_i i=fChildren.begin(); i!=fChildren.end(); ++i)
    (*i)->SetRnrSelf(kFALSE);

  gReve->Redraw3D();
}

/**************************************************************************/

Int_t RenderElement::DestroyListSubTree(TGListTree* ltree,
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

void RenderElement::Destroy()
{
  static const Exc_t eH("RenderElement::Destroy ");

  if(fDenyDestroy)
    throw(eH + "this object is protected against destruction.");

  gReve->PreDeleteRenderElement(this);
  delete this;
  gReve->Redraw3D();
}

void RenderElement::DestroyElements()
{
  static const Exc_t eH("RenderElement::DestroyElements ");

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

const TGPicture*  RenderElement::GetCheckBoxPicture(Bool_t rnrSelf, Bool_t rnrDaughters)
{
  Int_t idx = 0;
  if(rnrSelf) idx = 2;
  if( rnrDaughters ) idx ++;

  return fgRnrIcons[idx];
}

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

TObject* RenderElementObjPtr::GetObject(Reve::Exc_t eh)
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

ClassImp(RenderElementList)

RenderElementList::RenderElementList(const Text_t* n, const Text_t* t, Bool_t doColor) :
  RenderElement(),
  TNamed(n, t),
  fColor(0),
  fDoColor(doColor)
{
  if(fDoColor) {
    SetMainColorPtr(&fColor);
  }
}

/**************************************************************************/

ClassImp(PadPrimitive)

PadPrimitive::PadPrimitive(const Text_t* n, const Text_t* t) :
  RenderElement(),
  TNamed(n, t)
{}

void PadPrimitive::Paint(Option_t* option)
{
  if (fRnrChildren)
  {
    for(List_i i=fChildren.begin(); i!=fChildren.end(); ++i)
      (*i)->PadPaint(option);
  }
}

