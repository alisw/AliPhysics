// $Header$

#include "GeoNode.h"
#include <Reve/RGTopFrame.h>

#include <TPad.h>
#include <TColor.h>

#include <TGeoVolume.h>
#include <TGeoNode.h>
#include <TGeoManager.h>

using namespace Reve;

//______________________________________________________________________
// GeoNodeRnrEl
//

ClassImp(GeoNodeRnrEl)

GeoNodeRnrEl::GeoNodeRnrEl(TGeoNode* node) :
  RenderElementListBase(), 
  fNode(node)
{
  // Hack!! Should use cint to retrieve TAttLine::fLineColor offset.
  char* l = (char*) dynamic_cast<TAttLine*>(node->GetVolume());
  SetMainColorPtr((Color_t*)(l + sizeof(void*)));

  fRnrElement      = fNode->TGeoAtt::IsVisible();
}

const Text_t* GeoNodeRnrEl::GetName()  const { return fNode->GetName(); }
const Text_t* GeoNodeRnrEl::GetTitle() const { return fNode->GetTitle(); }

/**************************************************************************/

Int_t GeoNodeRnrEl::ExpandIntoListTree(TGListTree* ltree,
				       TGListTreeItem* parent)
{
  // Checks if child-nodes have been imported ... imports them if not.
  // Then calls RenderElementListBase::ExpandIntoListTree.

  if(fList.empty() && fNode->GetVolume()->GetNdaughters() > 0) {
    TIter next(fNode->GetVolume()->GetNodes());
    TGeoNode* dnode;
    while((dnode = (TGeoNode*) next()) != 0) {
      GeoNodeRnrEl* node_re = new GeoNodeRnrEl(dnode);
      AddElement(node_re);
    }
  }
  return RenderElementListBase::ExpandIntoListTree(ltree, parent);
}

/**************************************************************************/

void GeoNodeRnrEl::FullUpdate()
{
  fRnrElement      = fNode->TGeoAtt::IsVisible(); 
  RenderElementListBase::FullUpdate();
}

/**************************************************************************/

void GeoNodeRnrEl::SetRnrElement(Bool_t rnr)
{
  fNode->SetVisibility(rnr);
  FullUpdate();
}

/**************************************************************************/

void GeoNodeRnrEl::SetMainColor(Color_t color)
{
  fNode->GetVolume()->SetLineColor(color);
  FullUpdate();
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
    FullUpdate();

  for(lpRE_i i=fList.begin(); i!=fList.end(); ++i) {
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
    FullUpdate();

  for(lpRE_i i=fList.begin(); i!=fList.end(); ++i) {
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
/**************************************************************************/

ClassImp(GeoTopNodeRnrEl)

GeoTopNodeRnrEl::GeoTopNodeRnrEl(TGeoManager* manager, TGeoNode* node,
				 Int_t visopt, Int_t vislvl) :
  GeoNodeRnrEl(node),
  fManager(manager),
  fVisOption(visopt), fVisLevel(vislvl)
{
  fRnrElement = true;
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

void GeoTopNodeRnrEl::FullUpdate()
{
  RenderElementListBase::FullUpdate();
}

/**************************************************************************/

void GeoTopNodeRnrEl::SetRnrElement(Bool_t rnr)
{
  // Revert from GeoNode to back to standard behaviour.
  RenderElementListBase::SetRnrElement(rnr);
}

/**************************************************************************/

void GeoTopNodeRnrEl::Draw(Option_t* option)
{
  AppendPad(option);
}

void GeoTopNodeRnrEl::Paint(Option_t* option)
{
  if(fRnrElement) {
    gGeoManager = fManager;
    TVirtualPad* pad = gPad;
    gPad = 0;
    TGeoVolume* top_volume = fManager->GetTopVolume();
    fManager->SetVisOption(fVisOption);
    fManager->SetVisLevel(fVisLevel);
    fManager->SetTopVolume(fNode->GetVolume());
    gPad = pad;
    fNode->Paint(option);
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
