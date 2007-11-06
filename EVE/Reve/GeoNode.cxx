// $Header$

#include "GeoNode.h"
#include <Reve/ReveManager.h>
#include <Reve/NLTPolygonSet.h>

#include "TGeoShapeExtract.h"

#include <TPad.h>
#include <TBuffer3D.h>
#include <TVirtualViewer3D.h>
#include <TColor.h>
#include <TFile.h>

#include <TGeoShape.h>
#include <TGeoVolume.h>
#include <TGeoNode.h>
#include <TGeoShapeAssembly.h>
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

void GeoNodeRnrEl::Save(const char* file, const char* name)
{
  TGeoShapeExtract* gse = DumpShapeTree(this, 0, 0);
  
  TFile f(file, "RECREATE");
  gse->Write(name);
  f.Close();
}

/**************************************************************************/

TGeoShapeExtract* GeoNodeRnrEl::DumpShapeTree(GeoNodeRnrEl* geon, TGeoShapeExtract* parent, Int_t level)
{
  printf("dump_shape_tree %s \n", geon->GetName());
  TGeoNode*   tnode   = 0;
  TGeoVolume* tvolume = 0;
  TGeoShape*  tshape  = 0;

  tnode = geon->GetNode();
  if(tnode == 0) {
    printf("Null node for %s; assuming it's a holder and descending.\n", geon->GetName());
    goto do_dump;
  }

  tvolume = tnode->GetVolume();
  if(tvolume == 0) {
    printf("Null volume for %s; skipping.\n", geon->GetName());
    return 0;
  }

  tshape  = tvolume->GetShape();

do_dump:
  // transformation
  ZTrans trans;
  if (parent) if (parent) trans.SetFromArray(parent->GetTrans());
  TGeoMatrix* gm =  tnode->GetMatrix();
  const Double_t* rm = gm->GetRotationMatrix();
  const Double_t* tv = gm->GetTranslation();
  ZTrans t;
  t(1,1) = rm[0]; t(1,2) = rm[1]; t(1,3) = rm[2];
  t(2,1) = rm[3]; t(2,2) = rm[4]; t(2,3) = rm[5];
  t(3,1) = rm[6]; t(3,2) = rm[7]; t(3,3) = rm[8];
  t(1,4) = tv[0]; t(2,4) = tv[1]; t(3,4) = tv[2];
  trans *= t;

  TGeoShapeExtract* gse = new TGeoShapeExtract(geon->GetName(), geon->GetTitle());
  gse->SetTrans(trans.Array());
  Int_t ci = 0;
  if(tvolume) ci = tvolume->GetLineColor();
  TColor* c = gROOT->GetColor(ci);
  Float_t rgba[4] = {1, 0, 0, 1};
  if (c) {
    rgba[0] = c->GetRed();
    rgba[1] = c->GetGreen();
    rgba[2] = c->GetBlue();
  } 
  gse->SetRGBA(rgba);
  Bool_t rnr = geon->GetRnrSelf();
  if(level > gGeoManager->GetVisLevel())
    rnr = kFALSE;
  gse->SetRnrSelf(rnr);
  gse->SetRnrElements(geon->GetRnrChildren());

  if(dynamic_cast<TGeoShapeAssembly*>(tshape)){
    //    printf("<TGeoShapeAssembly \n");
    tshape = 0;
  }
  gse->SetShape(tshape);
  level ++;
  if ( geon->GetNChildren())
  {
    TList* ele = new TList();
    gse->SetElements(ele);
    gse->GetElements()->SetOwner(true);

    RenderElement::List_i i = geon->BeginChildren();
    while (i != geon->EndChildren()) {
      GeoNodeRnrEl* l = dynamic_cast<GeoNodeRnrEl*>(*i);
      DumpShapeTree(l, gse, level+1);
      i++;
    }
  }

  if(parent)
    parent->GetElements()->Add(gse);
   
  return gse;
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
  fManager     (manager),
  fGlobalTrans (),
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
      TGeoHMatrix geomat;
      fGlobalTrans.SetGeoHMatrix(geomat);
      vgp->PaintNode(fNode, option, &geomat);
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
  RenderElement(),
  TNamed        (name, title),
  fHMTrans      (),
  fColor        (0),
  fTransparency (0),
  fShape        (0)
{
  fMainColorPtr = &fColor;
}

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
  buff.fLocalFrame   = kTRUE; // Always enforce local frame (no geo manager).

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

void GeoShapeRnrEl::Save(const char* file, const char* name)
{
  TGeoShapeExtract* gse = DumpShapeTree(this, 0);
  
  TFile f(file, "RECREATE");
  gse->Write(name);
  f.Close();
}

/**************************************************************************/

TGeoShapeExtract* GeoShapeRnrEl::DumpShapeTree(GeoShapeRnrEl* gsre, TGeoShapeExtract* parent)
{
  //  printf("dump_shape_tree %s \n", gsre->GetName());
  TGeoShapeExtract* she = new TGeoShapeExtract(gsre->GetName(), gsre->GetTitle());
  she->SetTrans(gsre->RefHMTrans().Array());
  Int_t ci = gsre->GetColor();
  TColor* c = gROOT->GetColor(ci);
  Float_t rgba[4] = {1, 0, 0, 1 - gsre->GetMainTransparency()/100.};
  if (c)
 {
    rgba[0] = c->GetRed();
    rgba[1] = c->GetGreen();
    rgba[2] = c->GetBlue();
  } 
  she->SetRGBA(rgba);
  she->SetRnrSelf(gsre->GetRnrSelf());
  she->SetRnrElements(gsre->GetRnrChildren());
  she->SetShape(gsre->GetShape());
  if ( gsre->GetNChildren())
  {
    TList* ele = new TList();
    she->SetElements(ele);
    she->GetElements()->SetOwner(true);
    RenderElement::List_i i = gsre->BeginChildren();
    while (i != gsre->EndChildren()) {
      GeoShapeRnrEl* l = dynamic_cast<GeoShapeRnrEl*>(*i);
      DumpShapeTree(l, she);
      i++;
    }
  }
  if(parent)
    parent->GetElements()->Add(she);

  return she; 
}

GeoShapeRnrEl* GeoShapeRnrEl::ImportShapeExtract(TGeoShapeExtract * gse,
					RenderElement    * parent)
{
  gReve->DisableRedraw();
  GeoShapeRnrEl* gsre = SubImportShapeExtract(gse, parent);
  gsre->ElementChanged();
  gReve->EnableRedraw();
  return gsre;
}


GeoShapeRnrEl* GeoShapeRnrEl::SubImportShapeExtract(TGeoShapeExtract * gse,
					   RenderElement    * parent)
{
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

  gReve->AddGlobalRenderElement(gsre, parent);

  if (gse->HasElements())
  {
    TIter next(gse->GetElements());
    TGeoShapeExtract* chld;
    while ((chld = (TGeoShapeExtract*) next()) != 0)
     SubImportShapeExtract(chld, gsre);
  }

  return gsre;
}

/**************************************************************************/

TClass* GeoShapeRnrEl::ProjectedClass() const
{
   return NLTPolygonSet::Class();
} 

/**************************************************************************/

TBuffer3D* GeoShapeRnrEl::MakeBuffer3D()
{
  if(fShape == 0) return 0;

  if(dynamic_cast<TGeoShapeAssembly*>(fShape)){
    // !!!! TGeoShapeAssembly makes a bad TBuffer3D
    return 0;
  }

  TBuffer3D* buff  = fShape->MakeBuffer3D();
  if (fHMTrans.GetUseTrans())
  {
    Reve::ZTrans& mx = RefHMTrans();
    Int_t N = buff->NbPnts();
    Double_t* pnts = buff->fPnts;
    for(Int_t k=0; k<N; k++)
    {
      mx.MultiplyIP(&pnts[3*k]);
    }
  }
  return buff;
}



