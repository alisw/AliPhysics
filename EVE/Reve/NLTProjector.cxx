#include "NLTProjector.h"
#include "ReveManager.h"
#include "NLTBases.h"

#include "TBuffer3D.h"
#include "TBuffer3DTypes.h"
#include <TVirtualPad.h>
#include <TVirtualViewer3D.h>

#include <list>

using namespace Reve;

//______________________________________________________________________________
// NLTProjector
//
// Recursively projects RenderElement and draws axis in the projected scene.
// It enables to interactivly set NLTProjection parameters and updates
// projected scene accordingly.

ClassImp(NLTProjector)

NLTProjector::NLTProjector():
  RenderElementList("NLTProjector",""),

  fProjection (0),

  fDrawCenter(kFALSE),
  fDrawOrigin(kFALSE),

  fSplitInfoMode(0),
  fSplitInfoLevel(1),
  fAxisColor(0),

  fCurrentDepth(0)
{
  // Constructor.

  fProjection  = new NLTCircularFishEye(fCenter);
  UpdateName();
}

//______________________________________________________________________________
NLTProjector::~NLTProjector()
{
  // Destructor.

  if(fProjection) delete fProjection;
}

//______________________________________________________________________________
void NLTProjector::UpdateName()
{
  // Updates name to have consitent information with prjection.

  SetName(Form ("%s (%3.1f)", fProjection->GetName(), fProjection->GetDistortion()*1000));
  UpdateItems();
}

//______________________________________________________________________________
void NLTProjector::SetProjection(NLTProjection::PType_e type, Float_t distort)
{
  // Set projection type and distortion.

  static const Exc_t eH("NLTProjector::SetProjection ");

  delete fProjection;
  fProjection = 0;

  switch (type)
  {
    case NLTProjection::PT_CFishEye:
    {
      fProjection  = new NLTCircularFishEye(fCenter);
      break;
    }
    case NLTProjection::PT_RhoZ:
    {
      fProjection  = new NLTRhoZ(fCenter);
      break;
    }
    default:
      throw(eH + "projection type not valid.");
      break;
  }
  fProjection->SetDistortion(distort);
  UpdateName();
}
//______________________________________________________________________________
void NLTProjector::SetCenter(Float_t x, Float_t y, Float_t z)
{
  // Set projection center and rebuild projected scene.

  fCenter.Set(x, y, z);
  fProjection->SetCenter(fCenter);
  ProjectChildren();
}

//______________________________________________________________________________
Bool_t NLTProjector::HandleElementPaste(RenderElement* el)
{
  // React to element being pasted or dnd-ed.
  // Return true if redraw is needed (virtual method).

  size_t n_children  = fChildren.size();
  ImportElements(el);
  return n_children != fChildren.size();
}

//______________________________________________________________________________
Bool_t NLTProjector::ShouldImport(RenderElement* rnr_el)
{
  // Returns true if rnr_el or any of its children is NTLProjectable.

  if (rnr_el->IsA()->InheritsFrom(NLTProjectable::Class()))
    return kTRUE;
  for (List_i i=rnr_el->BeginChildren(); i!=rnr_el->EndChildren(); ++i)
    if (ShouldImport(*i))
      return kTRUE;
  return kFALSE;
}

//______________________________________________________________________________
void NLTProjector::ImportElementsRecurse(RenderElement* rnr_el, RenderElement* parent)
{
  // If rnr_el is NLTProjectable add projected instance else add plain RenderElementList
  // to parent. Call same function on rnr_el children.

  if (ShouldImport(rnr_el))
  {
    RenderElement  *new_re = 0;
    NLTProjected   *new_pr = 0;
    NLTProjectable *pble   = dynamic_cast<NLTProjectable*>(rnr_el);
    if (pble)
    {
      new_re = (RenderElement*) pble->ProjectedClass()->New();
      new_pr = dynamic_cast<NLTProjected*>(new_re);
      new_pr->SetProjection(this, pble);
      new_pr->SetDepth(fCurrentDepth);
    }
    else
    {
      new_re = new RenderElementList;
    }
    TObject *tobj   = rnr_el->GetObject();
    new_re->SetRnrElNameTitle(Form("NLT %s", tobj->GetName()),
			      tobj->GetTitle());
    new_re->SetRnrSelf     (rnr_el->GetRnrSelf());
    new_re->SetRnrChildren(rnr_el->GetRnrChildren());
    gReve->AddRenderElement(new_re, parent);

    for (List_i i=rnr_el->BeginChildren(); i!=rnr_el->EndChildren(); ++i)
      ImportElementsRecurse(*i, new_re);
  }
}

//______________________________________________________________________________
void NLTProjector::ImportElements(RenderElement* rnr_el)
{
  // Recursively import elements and update projection on the projected objects.

  ImportElementsRecurse(rnr_el, this);
  ProjectChildren();
}

//______________________________________________________________________________
void NLTProjector::ProjectChildrenRecurse(RenderElement* rnr_el)
{
  // Go recursively through rnr_el tree and call UpdateProjection() on NLTProjected.

  NLTProjected* pted = dynamic_cast<NLTProjected*>(rnr_el);
  if (pted)
  {
    pted->UpdateProjection();
    TAttBBox* bb = dynamic_cast<TAttBBox*>(pted);
    if(bb)
    {
      Float_t* b = bb->AssertBBox();
      BBoxCheckPoint(b[0], b[2], b[4]);
      BBoxCheckPoint(b[1], b[3], b[5]);
    }
    rnr_el->ElementChanged(kFALSE);
  }

  for (List_i i=rnr_el->BeginChildren(); i!=rnr_el->EndChildren(); ++i)
    ProjectChildrenRecurse(*i);
}

//______________________________________________________________________________
void NLTProjector::ProjectChildren()
{
  // Project children recursevly, update BBox and notify ReveManger
  // the scenes have chenged.

  BBoxZero();
  ProjectChildrenRecurse(this);
  AssertBBoxExtents(0.1);
  {
    using namespace TMath;
    fBBox[0] = 10.0f * Floor(fBBox[0]/10.0f);
    fBBox[1] = 10.0f * Ceil (fBBox[1]/10.0f);
    fBBox[2] = 10.0f * Floor(fBBox[2]/10.0f);
    fBBox[3] = 10.0f * Ceil (fBBox[3]/10.0f);
  }

  List_t scenes;
  CollectSceneParentsFromChildren(scenes, 0);
  gReve->ScenesChanged(scenes);
}

//______________________________________________________________________________
void NLTProjector::Paint(Option_t* /*option*/)
{
  // Paint this object. Only direct rendering is supported.

  static const Exc_t eH("NLTProjector::Paint ");
  TBuffer3D buff(TBuffer3DTypes::kGeneric);

  // Section kCore
  buff.fID           = this;
  buff.fColor        = fAxisColor;
  buff.fTransparency = 0;
  buff.SetSectionsValid(TBuffer3D::kCore);

  Int_t reqSections = gPad->GetViewer3D()->AddObject(buff);
  if (reqSections != TBuffer3D::kNone)
    Error(eH, "only direct GL rendering supported.");
}


//______________________________________________________________________________
void NLTProjector::ComputeBBox()
{
  // Fill bounding-box information in base-class TAttBBox (virtual method).

  static const Exc_t eH("NLTProjector::ComputeBBox ");

  if(GetNChildren() == 0) {
    BBoxZero();
    return;
  }

  BBoxInit();
}
