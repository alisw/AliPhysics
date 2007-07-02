// $Header$

#include "GeoNode.h"
#include <Reve/RGTopFrame.h>

#include "TGeoShapeExtract.h"

#include <TPad.h>
#include <TBuffer3D.h>
#include <TVirtualViewer3D.h>
#include <TColor.h>

#include <TGeoShape.h>
#include <TGeoVolume.h>
#include <TGeoNode.h>
#include <TGeoManager.h>
#include <TVirtualGeoPainter.h>

using namespace Reve;

//______________________________________________________________________
// GeoNodeRnrEl
//

ClassImp(GeoNodeRnrEl)

GeoNodeRnrEl::GeoNodeRnrEl(TGeoNode* node) :
  RenderElement(),
  TObject(),
  fNode(node)
{
  // Hack!! Should use cint to retrieve TAttLine::fLineColor offset.
  char* l = (char*) dynamic_cast<TAttLine*>(node->GetVolume());
  SetMainColorPtr((Color_t*)(l + sizeof(void*)));

  fRnrSelf      = fNode->TGeoAtt::IsVisible();
}

const Text_t* GeoNodeRnrEl::GetName()  const { return fNode->GetName(); }
const Text_t* GeoNodeRnrEl::GetTitle() const { return fNode->GetTitle(); }

/**************************************************************************/

Int_t GeoNodeRnrEl::ExpandIntoListTree(TGListTree* ltree,
				       TGListTreeItem* parent)
{
  // Checks if child-nodes have been imported ... imports them if not.
  // Then calls RenderElement::ExpandIntoListTree.

  if(fChildren.empty() && fNode->GetVolume()->GetNdaughters() > 0) {
    TIter next(fNode->GetVolume()->GetNodes());
    TGeoNode* dnode;
    while((dnode = (TGeoNode*) next()) != 0) {
      GeoNodeRnrEl* node_re = new GeoNodeRnrEl(dnode);
      AddElement(node_re);
    }
  }
  return RenderElement::ExpandIntoListTree(ltree, parent);
}

/**************************************************************************/

void GeoNodeRnrEl::SetRnrSelf(Bool_t rnr)
{
  RenderElement::SetRnrSelf(rnr);
  fNode->SetVisibility(rnr);
}

void GeoNodeRnrEl::SetRnrChildren(Bool_t rnr)
{
  RenderElement::SetRnrChildren(rnr);
  fNode->VisibleDaughters(rnr);
}

void GeoNodeRnrEl::SetRnrState(Bool_t rnr)
{
  RenderElement::SetRnrState(rnr);
  fNode->SetVisibility(rnr);
  fNode->VisibleDaughters(rnr);
}

/**************************************************************************/

void GeoNodeRnrEl::SetMainColor(Color_t color)
{
  fNode->GetVolume()->SetLineColor(color);
  UpdateItems();
}

void GeoNodeRnrEl::SetMainColor(Pixel_t pixel)
{
  // This one needed for proper calling via CINT (signals).

  SetMainColor(Color_t(TColor::GetColor(pixel)));
}

/**************************************************************************/

void GeoNodeRnrEl::UpdateNode(TGeoNode* node)
{
  // Updates all reve-browsers having the node in their contents.
  // All 3D-pads updated if any change found.
  //
  // Should (could?) be optimized with some assumptions about
  // volume/node structure (search for parent, know the same node can not
  // reoccur on lower level once found).

  static const Exc_t eH("GeoNodeRnrEl::UpdateNode ");

  // printf("%s node %s %p\n", eH.Data(), node->GetName(), node);

  if(fNode == node)
    UpdateItems();

  for(List_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
    ((GeoNodeRnrEl*)(*i))->UpdateNode(node);
  }

}

void GeoNodeRnrEl::UpdateVolume(TGeoVolume* volume)
{
  // Updates all reve-browsers having the volume in their contents.
  // All 3D-pads updated if any change found.
  //
  // Should (could?) be optimized with some assumptions about
  // volume/node structure (search for parent, know the same node can not
  // reoccur on lower level once found).

  static const Exc_t eH("GeoNodeRnrEl::UpdateVolume ");

  // printf("%s volume %s %p\n", eH.Data(), volume->GetName(), volume);

  if(fNode->GetVolume() == volume)
    UpdateItems();

  for(List_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
    ((GeoNodeRnrEl*)(*i))->UpdateVolume(volume);
  }
}

/**************************************************************************/

void GeoNodeRnrEl::Draw(Option_t* option)
{
  TString opt("SAME");
  opt += option;
  fNode->GetVolume()->Draw(opt);
}

/**************************************************************************/
//______________________________________________________________________
// GeoTopNodeRnrEl
//
// A wrapper over a TGeoNode, possibly displaced with a global
// trasformation fGlobalTrans (the matrix is owned by this class).
/**************************************************************************/

ClassImp(GeoTopNodeRnrEl)

GeoTopNodeRnrEl::GeoTopNodeRnrEl(TGeoManager* manager, TGeoNode* node,
				 Int_t visopt, Int_t vislvl) :
  GeoNodeRnrEl (node),
  fGlobalTrans (),
  fManager     (manager),
  fVisOption   (visopt),
  fVisLevel    (vislvl)
{
  fRnrSelf = true;
}

GeoTopNodeRnrEl::~GeoTopNodeRnrEl()
{}

/**************************************************************************/

void GeoTopNodeRnrEl::SetGlobalTrans(const TGeoHMatrix* m)
{
  fGlobalTrans.SetFrom(*m);
}

void GeoTopNodeRnrEl::UseNodeTrans()
{
  fGlobalTrans.SetFrom(*fNode->GetMatrix());
}

/**************************************************************************/

void GeoTopNodeRnrEl::SetVisOption(Int_t visopt)
{
  fVisOption = visopt;
  gReve->Redraw3D();
}

void GeoTopNodeRnrEl::SetVisLevel(Int_t vislvl)
{
  fVisLevel = vislvl;
  gReve->Redraw3D(); 
}

/**************************************************************************/

void GeoTopNodeRnrEl::SetRnrSelf(Bool_t rnr)
{
  // Revert from GeoNode to back to standard behaviour.
  RenderElement::SetRnrSelf(rnr);
}

/**************************************************************************/

void GeoTopNodeRnrEl::Draw(Option_t* option)
{
  AppendPad(option);
}

void GeoTopNodeRnrEl::Paint(Option_t* option)
{
  if(fRnrSelf) {
    gGeoManager = fManager;
    TVirtualPad* pad = gPad;
    gPad = 0;
    TGeoVolume* top_volume = fManager->GetTopVolume();
    fManager->SetVisOption(fVisOption);
    fManager->SetVisLevel(fVisLevel);
    fManager->SetTopVolume(fNode->GetVolume());
    gPad = pad;
    TVirtualGeoPainter* vgp = fManager->GetGeomPainter();
    if(vgp != 0) {
#if ROOT_VERSION_CODE > ROOT_VERSION(5,11,6)
      TGeoHMatrix geomat;
      fGlobalTrans.SetGeoHMatrix(geomat);
      vgp->PaintNode(fNode, option, &geomat);
#else
      vgp->PaintNode(fNode, option);
#endif
    }
    fManager->SetTopVolume(top_volume);
  }
}

/**************************************************************************/

void GeoTopNodeRnrEl::VolumeVisChanged(TGeoVolume* volume)
{
  static const Exc_t eH("GeoTopNodeRnrEl::VolumeVisChanged ");
  printf("%s volume %s %p\n", eH.Data(), volume->GetName(), (void*)volume);
  UpdateVolume(volume);
}

void GeoTopNodeRnrEl::VolumeColChanged(TGeoVolume* volume)
{
  static const Exc_t eH("GeoTopNodeRnrEl::VolumeColChanged ");
  printf("%s volume %s %p\n", eH.Data(), volume->GetName(), (void*)volume);
  UpdateVolume(volume);
}

void GeoTopNodeRnrEl::NodeVisChanged(TGeoNode* node)
{
  static const Exc_t eH("GeoTopNodeRnrEl::NodeVisChanged ");
  printf("%s node %s %p\n", eH.Data(), node->GetName(), (void*)node);
  UpdateNode(node);
}


/**************************************************************************/
//______________________________________________________________________
// GeoShapeRnrEl
//
// Minimal shape-wrapper allowing import of stuff from gled and retaining
// user-set visibility, colors and transparency.
/**************************************************************************/

ClassImp(GeoShapeRnrEl)

GeoShapeRnrEl::GeoShapeRnrEl(const Text_t* name, const Text_t* title) :
  RenderElement (fColor),
  TNamed        (name, title),
  fHMTrans      (),
  fColor        (0),
  fTransparency (0),
  fShape        (0)
{}

GeoShapeRnrEl::~GeoShapeRnrEl()
{
  if (fShape) {
    fShape->SetUniqueID(fShape->GetUniqueID() - 1);
    if (fShape->GetUniqueID() == 0)
      delete fShape;
  }
}

/**************************************************************************/

void GeoShapeRnrEl::Paint(Option_t* /*option*/)
{
  if (fShape == 0)
    return;

  TBuffer3D& buff = (TBuffer3D&) fShape->GetBuffer3D
    (TBuffer3D::kCore, false);

  buff.fID           = this;
  buff.fColor        = fColor;
  buff.fTransparency = fTransparency;
  fHMTrans.SetBuffer3D(buff);

  fShape->GetBuffer3D(TBuffer3D::kBoundingBox | TBuffer3D::kShapeSpecific, true);

  Int_t reqSec = gPad->GetViewer3D()->AddObject(buff);

  if (reqSec != TBuffer3D::kNone) {
    fShape->GetBuffer3D(reqSec, true);
    reqSec = gPad->GetViewer3D()->AddObject(buff);
  }

  if (reqSec != TBuffer3D::kNone)
    printf("spooky reqSec=%d for %s\n", reqSec, GetName());
}

/**************************************************************************/

Int_t GeoShapeRnrEl::ImportShapeExtract(TGeoShapeExtract * gse,
					RenderElement    * parent)
{
  gReve->DisableRedraw();
  Int_t n = SubImportShapeExtract(gse, parent);
  printf ("GeoShapeRnrEl::ImportShapeExtract imported %d elements\n", n);
  gReve->EnableRedraw();
  return n;
}

Int_t GeoShapeRnrEl::SubImportShapeExtract(TGeoShapeExtract * gse,
					   RenderElement    * parent)
{
  Int_t ncreated = 1;

  GeoShapeRnrEl* gsre = new GeoShapeRnrEl(gse->GetName(), gse->GetTitle());
  gsre->fHMTrans.SetFromArray(gse->GetTrans());
  const Float_t* rgba = gse->GetRGBA();
  gsre->fColor        = TColor::GetColor(rgba[0], rgba[1], rgba[2]);
  gsre->fTransparency = (UChar_t) (100.0f*(1.0f - rgba[3]));
  gsre->SetRnrSelf(gse->GetRnrSelf());
  gsre->SetRnrChildren(gse->GetRnrElements());
  gsre->fShape = gse->GetShape();
  if (gsre->fShape)
    gsre->fShape->SetUniqueID(gsre->fShape->GetUniqueID() + 1);

  if (parent)
    gReve->AddGlobalRenderElement(parent, gsre);
  else
    gReve->AddGlobalRenderElement(gsre);

  if (gse->HasElements())
  {
    TIter next(gse->GetElements());
    TGeoShapeExtract* chld;
    while ((chld = (TGeoShapeExtract*) next()) != 0)
      ncreated += SubImportShapeExtract(chld, gsre);
  }

  return ncreated;
}
