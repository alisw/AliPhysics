#include "NLTProjector.h"
#include "ReveManager.h"
#include "NLTPolygonSet.h"
#include "PointSet.h"
#include "Track.h"

#include "TQObject.h"
#include "TBuffer3D.h"
#include "TBuffer3DTypes.h"
#include "TColor.h"
#include "TPointSet3D.h"
#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include <TVirtualPad.h>
#include <TVirtualViewer3D.h>

#include <list>

using namespace Reve;

ClassImp(Reve::NLTProjection)

Float_t NLTProjection::fgEps = 0.01f;

//______________________________________________________________________________
void NLTProjection::ProjectVector(Vector& v)
{
  ProjectPoint(v.x, v.y, v.z);
}

//______________________________________________________________________________
Vector* NLTProjection::Project(Vector* origPnts, Int_t Npnts, Bool_t copy)
{
  Vector* pnts = 0; 
  if(copy) 
  {
    pnts = new Vector[Npnts];
    memcpy(pnts, origPnts, Npnts*sizeof(Vector));
  }
  else
  { 
    pnts =  origPnts;
  }

  for(Int_t i = 0; i<Npnts; i++)
  {
    ProjectPoint(pnts[i].x, pnts[i].y, pnts[i].z);
  }
  return pnts;
}

//______________________________________________________________________________
void NLTProjection::SetDistortion(Float_t d)
{
  // Prevent scaling down whole projected scene.
  // Point at fixe distance shold be on constant screen coorinate. 

  fDistortion=d; 
  fScale = 1+fFixedRadius*fDistortion;
}

//______________________________________________________________________________
void NLTProjection::SetDirectionalVector(Int_t screenAxis, Vector& vec)
{
  for(Int_t i=0; i<3; i++)
  {
    vec[i] = (i==screenAxis) ? 1. : 0.;
  }
}

//______________________________________________________________________________
Float_t NLTProjection::GetValForScreenPos(Int_t i, Float_t sv)
{
  // move on unprojected axis, find when projected value is geiven screen values
  Float_t xL, xM, xR;

  Vector V, DirVec;
  SetDirectionalVector(i, DirVec);
  
  // search from -/+ infinity according to sign of screen value
  if (sv > 0)
  {
    xL = 0; xR = 1000;
    while (1)
    {
      V.Mult(DirVec, xR); ProjectVector(V);
      if (V[i] > sv) break;
      xL = xR; xR *= 2;
    }
  }
  else if (sv < 0)
  {
    xR = 0; xL = -1000;
    while (1)
    {
      V.Mult(DirVec, xL); ProjectVector(V);
      if (V[i] < sv) break;
      xR = xL; xL *= 2;
    }
  }
  else
  {
    return 0;
  }

  //  printf("start bisection %f, %f , SCREEN VALUE %f \n", xL, xR, sv);
  do 
  {  
    xM = 0.5f * (xL + xR);
    V.Mult(DirVec, xM);
    ProjectVector(V);
    if (V[i] > sv) 
      xR = xM; 
    else
      xL = xM;
  } while(TMath::Abs(V[i] - sv) > fgEps);

  return xM;
}

//______________________________________________________________________________
Float_t NLTProjection::GetScreenVal(Int_t i, Float_t x)
{
  Vector dv; 
  SetDirectionalVector(i, dv); dv = dv*x;
  ProjectVector(dv);
  return dv[i];
}


/**************************************************************************/
/**************************************************************************/

ClassImp(Reve::RhoZ)

void RhoZ::ProjectPoint(Float_t& x, Float_t& y, Float_t& z)
{
    Float_t R = TMath::Sqrt(x*x+y*y);
    Float_t NR = R / (1.0f + R*fDistortion);
    y = (y > 0) ? NR : -NR;
    y *= fScale;
    x = (z*fScale) / (1.0f + TMath::Abs(z)*fDistortion);
    z = 0;
}


//______________________________________________________________________________
void RhoZ::SetDirectionalVector(Int_t screenAxis, Vector& vec)
{
  if(screenAxis == 0) 
    vec.Set(0., 0., 1);
  else 
    vec.Set(0., 1., 0);

}
//______________________________________________________________________________
Bool_t RhoZ::AcceptSegment(Vector& v1, Vector& v2, Float_t tolerance) 
{
  if((v1.y < 0 && v2.y > 0) || (v1.y > 0 && v2.y < 0))
  {
    if (tolerance > 0)
    {
      Float_t a1 = TMath::Abs(v1.y), a2 = TMath::Abs(v2.y);
      if (a1 < a2)
      {
	if (a1 < tolerance) { v1.y = 0; return kTRUE; }
      }
      else
      {
	if (a2 < tolerance) { v2.y = 0; return kTRUE; }
      }
    }
    return kFALSE;
  }
  return kTRUE;
}

/**************************************************************************/
/**************************************************************************/
ClassImp(Reve:: CircularFishEye)
//______________________________________________________________________________
void  CircularFishEye::ProjectPoint(Float_t& x, Float_t& y, Float_t& z) 
{
  Float_t R = TMath::Sqrt(x*x+y*y);
  Float_t phi = x == 0.0 && y == 0.0 ? 0.0 : TMath::ATan2(y,x);
  Float_t NR = (R*fScale) / (1.0f + R*fDistortion);
  x = NR*TMath::Cos(phi);
  y = NR*TMath::Sin(phi);
  z = 0;
}

/**************************************************************************/
/**************************************************************************/
//______________________________________________________________________
// NLTProjector
//

ClassImp(NLTProjector)

NLTProjector::NLTProjector():
  RenderElementList("NLTProjector",""),

  fProjection (0),

  fSplitInfoMode(0),
  fSplitInfoLevel(1),
  fAxisColor(0),

  fCurrentDepth(0)
{
  SetProjection(NLTProjection::PT_CFishEye, 0);
}

//______________________________________________________________________________
NLTProjector::~NLTProjector()
{
  if(fProjection) delete fProjection;
}

//______________________________________________________________________________
void NLTProjector::UpdateName()
{
  SetName(Form ("%s (%3.1f)", fProjection->GetName(), fProjection->GetDistortion()*1000));
  UpdateItems();
}

//______________________________________________________________________________
void NLTProjector::SetProjection(NLTProjection::PType_e type, Float_t distort)
{
  static const Exc_t eH("NLTProjector::SetProjection ");

  delete fProjection;
  fProjection = 0;

  switch (type)
  {
    case NLTProjection::PT_CFishEye:
    {
      fProjection  = new CircularFishEye();
      break;
    }
    case NLTProjection::PT_RhoZ:
    {
      fProjection  = new RhoZ;
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
void NLTProjector::SetProjection(NLTProjection* p)
{
  delete fProjection;
  fProjection = p;
  UpdateName();
}

//______________________________________________________________________________
Bool_t NLTProjector::HandleElementPaste(RenderElement* el)
{
  size_t n_children  = fChildren.size();
  ImportElements(el);
  return n_children != fChildren.size();
}

//______________________________________________________________________________
Bool_t NLTProjector::ShouldImport(RenderElement* rnr_el)
{
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
  ImportElementsRecurse(rnr_el, this);
  ProjectChildren();
}

//______________________________________________________________________________
void NLTProjector::ProjectChildrenRecurse(RenderElement* rnr_el)
{
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

