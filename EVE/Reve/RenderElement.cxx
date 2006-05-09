// $Header$

#include "RenderElement.h"
#include "RGTopFrame.h"

#include <TColor.h>
#include <TCanvas.h>
#include <TGListTree.h>
#include <THashList.h>

using namespace Reve;
using namespace Reve;

//______________________________________________________________________
// Reve
//

ClassImp(RenderElement)

RenderElement::RenderElement()
{
  fRnrElement    = true;
  fMainColorPtr  = 0;
}

RenderElement::RenderElement(Color_t& main_color) : fMainColorPtr(&main_color)
{
  fRnrElement    = true;
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
					       TGListTreeItem* parent)
{
  static const Exc_t eH("RenderElement::AddIntoListTree ");

  TObject* tobj = GetObject(eH);
  TGListTreeItem* item = ltree->AddItem(parent, tobj->GetName(), tobj,
					 0, 0, kTRUE);
  item->CheckItem(GetRnrElement());
  item->SetColor(GetMainColor());
  // printf("%s setting title for %s, '%s'\n", eH.Data(), tobj->GetName(), tobj->GetTitle());
  item->SetTipText(tobj->GetTitle());

  fItems.insert(ListTreeInfo(ltree, item));

  return item;
}

/**************************************************************************/

void RenderElement::FullUpdate()
{
  for(sLTI_i i=fItems.begin(); i!=fItems.end(); ++i) {
    // Setup name and title/tooltip? Need update calls from setname/title as well.
    i->fItem->CheckItem(fRnrElement);
    i->fItem->SetColor(GetMainColor());
  }
  gReve->Redraw3D();
  gReve->NotifyBrowser();
}

/**************************************************************************/

void RenderElement::SpawnEditor()
{
  // Here spawn a sub-class of TGedFrame with special UpdateMethod.
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
    FullUpdate();
  }
}
/**************************************************************************/

void RenderElement::SetMainColor(Color_t color)
{
  if (fMainColorPtr) {
    *fMainColorPtr = color;
    FullUpdate();
  }
}

void RenderElement::SetMainColorByPixel(Pixel_t pixel)
{
  SetMainColor(Color_t(TColor::GetColor(pixel)));
}

/**************************************************************************/
/**************************************************************************/

ClassImp(RenderElementListBase)

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
  for(lpRE_i i=fList.begin(); i!=fList.end(); ++i) {
    (*i)->AddIntoListTree(ltree, parent);
    ++n;
  }
  return n;
}

/**************************************************************************/

void RenderElementListBase::EnableListElements()
{
  for(lpRE_i i=fList.begin(); i!=fList.end(); ++i)
    (*i)->SetRnrElement(true);
}

void RenderElementListBase::DisableListElements()
{
  for(lpRE_i i=fList.begin(); i!=fList.end(); ++i)
    (*i)->SetRnrElement(false);
}

/**************************************************************************/

void RenderElementListBase::SetMainColor(Color_t col)
{
  Color_t oldcol = GetMainColor();
  for(lpRE_i i=fList.begin(); i!=fList.end(); ++i) {
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
    for(lpRE_i i=fList.begin(); i!=fList.end(); ++i) {
      if((*i)->GetRnrElement())
	(*i)->GetObject()->Paint(option);
    }
  }
}

/**************************************************************************/
/**************************************************************************/

ClassImp(RenderElementList)
