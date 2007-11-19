#include "NLTProjector.h"
#include "ReveManager.h"
#include "NLTBases.h"

#include "TBuffer3D.h"
#include "TBuffer3DTypes.h"
#include <TVirtualPad.h>
#include <TVirtualViewer3D.h>

#include <list>

using namespace Reve;

Float_t NLTProjection::fgEps = 0.005f;

//______________________________________________________________________________
// NLTProjection
//
// Base-class for non-linear projection of 3D point.
// Enables to define an external center of distortion and a scale to
// fixate a bounding box of a projected point.


ClassImp(Reve::NLTProjection)

//______________________________________________________________________________
NLTProjection::NLTProjection(Vector& center) :
  fType(PT_Unknown),
  fGeoMode(GM_Unknown),
  fName(0),
  fCenter(center.x, center.y, center.z),
  fDistortion(0.0f),
  fFixedRadius(300),
  fScale(1.0f)
{
  // Constructor.
}

//______________________________________________________________________________
void NLTProjection::ProjectVector(Vector& v)
{
  // Project Reve::Vector.

  ProjectPoint(v.x, v.y, v.z);
}

//______________________________________________________________________________
void NLTProjection::UpdateLimit()
{
  // Update convergence in +inf and -inf.

  if ( fDistortion == 0.0f )
    return;

  Float_t lim =  1.0f/fDistortion + fFixedRadius;
  Float_t* c = GetProjectedCenter();
  fUpLimit.Set(lim + c[0], lim + c[1], c[2]);
  fLowLimit.Set(-lim + c[0], -lim + c[1], c[2]);
}

//______________________________________________________________________________
void NLTProjection::SetDistortion(Float_t d)
{
  // Set distortion.

  fDistortion=d;
  fScale = 1+fFixedRadius*fDistortion;
  UpdateLimit();
}

//______________________________________________________________________________
void NLTProjection::SetFixedRadius(Float_t r)
{
  // Set fixed radius.

  fFixedRadius=r;
  fScale = 1 + fFixedRadius*fDistortion;
  UpdateLimit();
}

//______________________________________________________________________________
void NLTProjection::SetDirectionalVector(Int_t screenAxis, Vector& vec)
{
  // Get vector for axis in a projected space.

  for (Int_t i=0; i<3; i++)
  {
    vec[i] = (i==screenAxis) ? 1. : 0.;
  }
}

//______________________________________________________________________________
Float_t NLTProjection::GetValForScreenPos(Int_t i, Float_t sv)
{
  // Inverse projection.

  static const Exc_t eH("NLTProjector::GetValForScreenPos ");

  Float_t xL, xM, xR;
  Vector V, DirVec;
  SetDirectionalVector(i, DirVec);
  if (fDistortion > 0.0f && ((sv > 0 && sv > fUpLimit[i]) || (sv < 0 && sv < fLowLimit[i])))
    throw(eH + Form("screen value '%f' out of limit '%f'.", sv, sv > 0 ? fUpLimit[i] : fLowLimit[i]));

  Vector zero; ProjectVector(zero);
  // search from -/+ infinity according to sign of screen value
  if (sv > zero[i])
  {
    xL = 0; xR = 1000;
    while (1)
    {
      V.Mult(DirVec, xR); ProjectVector(V);
      // printf("positive projected %f, value %f,xL, xR ( %f, %f)\n", V[i], sv, xL, xR);
      if (V[i] > sv || V[i] == sv) break;
      xL = xR; xR *= 2;
    }
  }
  else if (sv < zero[i])
  {
    xR = 0; xL = -1000;
    while (1)
    {
      V.Mult(DirVec, xL); ProjectVector(V);
      // printf("negative projected %f, value %f,xL, xR ( %f, %f)\n", V[i], sv, xL, xR);
      if (V[i] < sv || V[i] == sv) break;
      xR = xL; xL *= 2;
    }
  }
  else
  {
    return 0.0f;
  }

  do
  {
    xM = 0.5f * (xL + xR);
    V.Mult(DirVec, xM);
    ProjectVector(V);
    if (V[i] > sv)
      xR = xM;
    else
      xL = xM;
  } while(TMath::Abs(V[i] - sv) >= fgEps);

  return xM;
}

//______________________________________________________________________________
Float_t NLTProjection::GetScreenVal(Int_t i, Float_t x)
{
  // Project point on given axis and return projected value.

  Vector dv;
  SetDirectionalVector(i, dv); dv = dv*x;
  ProjectVector(dv);
  return dv[i];
}

//______________________________________________________________________________
//
// NLTRhoZ
//
// Transformation from 3D to 2D. X axis represent Z coordinate. Y axis have value of
// radius with a sign of Y coordinate.


ClassImp(Reve::NLTRhoZ)

//______________________________________________________________________________
void NLTRhoZ::SetCenter(Vector& v)
{
  // Set center of distortion (virtual method).

  fCenter = v;

  Float_t R = TMath::Sqrt(v.x*v.x+v.y*v.y);
  fProjectedCenter.x = fCenter.z;
  fProjectedCenter.y = TMath::Sign(R, fCenter.y);
  fProjectedCenter.z = 0;
  UpdateLimit();
}

//______________________________________________________________________________
void NLTRhoZ::ProjectPoint(Float_t& x, Float_t& y, Float_t& z,  PProc_e proc )
{
  // Project point.

  using namespace TMath;

  if(proc == PP_Plane || proc == PP_Full)
  {
    // project
    y = Sign((Float_t)Sqrt(x*x+y*y), y);
    x = z;
  }
  if(proc == PP_Distort || proc == PP_Full)
  {
    // move to center
    x -= fProjectedCenter.x;
    y -= fProjectedCenter.y;
    // distort
    y = (y*fScale) / (1.0f + Abs(y)*fDistortion);
    x = (x*fScale) / (1.0f + Abs(x)*fDistortion);
    // move back from center
    x += fProjectedCenter.x;
    y += fProjectedCenter.y;
  }
  z = 0.0f;
}

//______________________________________________________________________________
void NLTRhoZ::SetDirectionalVector(Int_t screenAxis, Vector& vec)
{
  // Get direction in the unprojected space for axis index in the projected space.
  // This is virtual method from base-class NLTProjection.

  if(screenAxis == 0)
    vec.Set(0., 0., 1);
  else if (screenAxis == 1)
    vec.Set(0., 1., 0);

}
//______________________________________________________________________________
Bool_t NLTRhoZ::AcceptSegment(Vector& v1, Vector& v2, Float_t tolerance)
{
  // Check if segment of two projected points is valid.

  Float_t a = fProjectedCenter.y;
  Bool_t val = kTRUE;
  if((v1.y <  a && v2.y > a) || (v1.y > a && v2.y < a))
  {
    val = kFALSE;
    if (tolerance > 0)
    {
      Float_t a1 = TMath::Abs(v1.y - a), a2 = TMath::Abs(v2.y - a);
      if (a1 < a2)
      {
	if (a1 < tolerance) { v1.y = a; val = kTRUE; }
      }
      else
      {
	if (a2 < tolerance) { v2.y = a; val = kTRUE; }
      }
    }
  }
  return val;
}

//______________________________________________________________________________
//
// NLTCircularFishEye
//
// XY projection with distortion around given center.

ClassImp(Reve:: NLTCircularFishEye)

//______________________________________________________________________________
void  NLTCircularFishEye::ProjectPoint(Float_t& x, Float_t& y, Float_t& z, PProc_e proc)
{
  // Project point.

  using namespace TMath;

  if(proc !=  PP_Plane)
  {
    x -= fCenter.x;
    y -= fCenter.y;
    Float_t phi = x == 0.0 && y == 0.0 ? 0.0 : ATan2(y,x);
    Float_t R = Sqrt(x*x+y*y);
    // distort
    Float_t NR = (R*fScale) / (1.0f + R*fDistortion);
    x = NR*Cos(phi) + fCenter.x;
    y = NR*Sin(phi) + fCenter.y;
  }
  z = 0.0f;
}

//______________________________________________________________________________
//
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

