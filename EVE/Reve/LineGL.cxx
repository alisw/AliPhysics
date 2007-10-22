// $Header$

#include "LineGL.h"
#include <Reve/Line.h>

#include <TGLRnrCtx.h>
#include <Reve/GLUtilNS.h>
#include <TGLIncludes.h>

using namespace Reve;

//______________________________________________________________________
// LineGL
//

ClassImp(LineGL)

LineGL::LineGL() : TPointSet3DGL(), fM(0)
{
  // fDLCache = false; // Disable display list.
}

LineGL::~LineGL()
{}

/**************************************************************************/

Bool_t LineGL::SetModel(TObject* obj, const Option_t* /*opt*/)
{
  // TPointSet3DGL::SetModel(obj);
  if(SetModelCheckClass(obj, Line::Class())) {
    fM = dynamic_cast<Line*>(obj);
    return kTRUE;
  }
  return kFALSE;
}

/**************************************************************************/

void LineGL::DirectDraw(TGLRnrCtx & rnrCtx) const
{
  // Direct GL rendering for Line.

  // printf("LineGL::DirectDraw Style %d, LOD %d\n", rnrCtx.Style(), rnrCtx.LOD());

  if (rnrCtx.DrawPass() == TGLRnrCtx::kPassOutlineLine)
    return;

  Line& q = *fM;
  if (q.Size() <= 0) return;

  if (q.fRnrLine)
    GLUtilNS::RenderLine(q, q.GetP(), q.Size());

  if (q.fRnrPoints)
    GLUtilNS::RenderPolyMarkers(q, q.GetP(), q.Size());
}
