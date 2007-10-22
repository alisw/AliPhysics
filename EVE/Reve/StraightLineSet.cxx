// $Header$

#include "StraightLineSet.h"

#include <TBuffer3D.h>
#include <TBuffer3DTypes.h>
#include <TVirtualPad.h>
#include <TVirtualViewer3D.h>

#include <TRandom.h>
#include <Reve/NLTProjector.h>

using namespace Reve;

//______________________________________________________________________
// StraightLineSet
//

ClassImp(StraightLineSet)

StraightLineSet::StraightLineSet(const Text_t* n, const Text_t* t):
  RenderElement (),
  TNamed        (n, t),

  fLinePlex      (sizeof(Line), 4),
  fMarkerPlex    (sizeof(Marker), 8),
  fOwnLinesIds   (kFALSE),
  fOwnMarkersIds (kFALSE),
  fRnrMarkers    (kTRUE),
  fRnrLines      (kTRUE),
  fLastLine      (0),
  fTrans         (kFALSE),
  fHMTrans       ()
{
  fMainColorPtr = &fLineColor;
  fLineColor    = 4;
  fMarkerColor  = 2;
  fMarkerStyle  = 20;
}

/**************************************************************************/

void StraightLineSet::AddLine(Float_t x1, Float_t y1, Float_t z1, 
			      Float_t x2, Float_t y2, Float_t z2)
{
  fLastLine = new (fLinePlex.NewAtom()) Line(x1, y1, z1, x2, y2, z2);
}

/**************************************************************************/

void StraightLineSet::AddMarker(Int_t line, Float_t pos)
{
  /*Marker* marker = */new (fMarkerPlex.NewAtom()) Marker(line, pos);
}

/**************************************************************************/

void StraightLineSet::ComputeBBox()
{
  static const Exc_t eH("StraightLineSet::ComputeBBox ");
  if(fLinePlex.Size() == 0) {
    BBoxZero();
    return;
  }

  BBoxInit();

  VoidCPlex::iterator li(fLinePlex);
  while (li.next()) {
    BBoxCheckPoint(((Line*)li())->fV1);
    BBoxCheckPoint(((Line*)li())->fV2);
  }
}

/**************************************************************************/

void StraightLineSet::Paint(Option_t* /*option*/)
{
  static const Exc_t eH("StraightLineSet::Paint ");

  TBuffer3D buff(TBuffer3DTypes::kGeneric);

  // Section kCore
  buff.fID           = this;
  buff.fColor        = fLineColor;
  buff.fTransparency = 0;
  buff.fLocalFrame   = kFALSE;
  fHMTrans.SetBuffer3D(buff);
  buff.SetSectionsValid(TBuffer3D::kCore);

  Int_t reqSections = gPad->GetViewer3D()->AddObject(buff);
  if (reqSections != TBuffer3D::kNone)
    Error(eH, "only direct GL rendering supported.");
}

/**************************************************************************/

TClass* StraightLineSet::ProjectedClass() const
{
  return NLTSLineSet::Class();
}

//______________________________________________________________________
// NLTSLineSet
//

ClassImp(NLTSLineSet)

NLTSLineSet::NLTSLineSet() : StraightLineSet(), NLTProjected ()
{}

/**************************************************************************/

void NLTSLineSet::SetProjection(NLTProjector* proj, NLTProjectable* model)
{
  NLTProjected::SetProjection(proj, model);

  // copy line and marker attributes
  * (TAttMarker*)this = * dynamic_cast<TAttMarker*>(fProjectable);
  * (TAttLine*)this   = * dynamic_cast<TAttLine*>(fProjectable);
}
/**************************************************************************/

void NLTSLineSet::UpdateProjection()
{
  NLTProjection&   proj   = * fProjector->GetProjection();
  StraightLineSet& orig  = * dynamic_cast<StraightLineSet*>(fProjectable);

  // lines
  Int_t NL = orig.GetLinePlex().Size();
  fLinePlex.Reset(sizeof(Line), NL);
  Line* l;
  Float_t p1[3];
  Float_t p2[3];
  VoidCPlex::iterator li(orig.GetLinePlex());

  Double_t s1, s2, s3;
  orig.RefHMTrans().GetScale(s1, s2, s3);
  ZTrans mx; mx.Scale(s1, s2, s3);
  while (li.next()) 
  {
    l = (Line*) li();
    p1[0] = l->fV1[0];  p1[1] = l->fV1[1]; p1[2] = l->fV1[2];
    p2[0] = l->fV2[0];  p2[1] = l->fV2[1]; p2[2] = l->fV2[2];
    mx.MultiplyIP(p1);
    mx.MultiplyIP(p2);
    proj.ProjectPointFv(p1);
    proj.ProjectPointFv(p2);
    AddLine(p1[0], p1[1], p1[2], p2[0], p2[1], p2[2]);
  }

  // markers
  Int_t NM = orig.GetMarkerPlex().Size();
  fMarkerPlex.Reset(sizeof(Marker), NM);
  Marker* m;
  VoidCPlex::iterator mi(orig.GetMarkerPlex());
  while (mi.next()) 
  {
    m = (Marker*) mi();
    AddMarker(m->fLineID, m->fPos);
  }

  // set position
  Float_t pos[3];
  orig.RefHMTrans().GetPos(pos);
  proj.ProjectPointFv(pos);
  fHMTrans.SetPos(pos);
}
