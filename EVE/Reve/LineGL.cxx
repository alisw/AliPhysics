// $Header$

#include "LineGL.h"
#include <Reve/Line.h>

#include <TGLDrawFlags.h>

#ifdef WIN32
#include "Windows4root.h"
#endif
#include <GL/gl.h>
#include <GL/glu.h>

using namespace Reve;

//______________________________________________________________________
// LineGL
//

ClassImp(LineGL)

LineGL::LineGL() : TPointSet3DGL(), fM(0)
{
  // fCached = false; // Disable display list.
}

LineGL::~LineGL()
{}

/**************************************************************************/

Bool_t LineGL::SetModel(TObject* obj)
{
  // TPointSet3DGL::SetModel(obj);
  if(SetModelCheckClass(obj, Line::Class())) {
    fM = dynamic_cast<Line*>(obj);
    return kTRUE;
  }
  return kFALSE;
}

/**************************************************************************/

void LineGL::DirectDraw(const TGLDrawFlags & flags) const
{
  // printf("LineGL::DirectDraw Style %d, LOD %d\n", flags.Style(), flags.LOD());

  Line& q = *fM;
  if (q.Size() <= 0) return;

  glPushAttrib(GL_POINT_BIT | GL_LINE_BIT | GL_ENABLE_BIT);
  glDisable(GL_LIGHTING);

  UChar_t color[4];

  if (q.fRnrLine)
  {
    glPushAttrib(GL_LINE_BIT);
    ColorFromIdx(q.GetLineColor(), color);
    glColor4ubv(color);
    glLineWidth(q.GetLineWidth());
    if (q.GetLineStyle() > 1) {
      Int_t    fac = 1;
      UShort_t pat = 0xffff;
      switch (q.GetLineStyle()) {
      case 2:  pat = 0x3333; break;
      case 3:  pat = 0x5555; break;
      case 4:  pat = 0xf040; break;
      case 5:  pat = 0xf4f4; break;
      case 6:  pat = 0xf111; break;
      case 7:  pat = 0xf0f0; break;
      case 8:  pat = 0xff11; break;
      case 9:  pat = 0x3fff; break;
      case 10: pat = 0x08ff; fac = 2; break;
      }

      glLineStipple(1, pat);
      glEnable(GL_LINE_STIPPLE);
    }

    const Float_t* p = q.GetP();
    const Int_t    n = q.Size();
    glBegin(GL_LINE_STRIP);
    for (Int_t i=0; i<n; ++i, p+=3)
      glVertex3fv(p);
    glEnd();
    glPopAttrib();
  }

  if (q.fRnrPoints)
  {
    ColorFromIdx(q.GetMarkerColor(), color);
    glColor4ubv(color);
    Int_t ms = q.GetMarkerStyle();
    if (ms != 2 && ms != 3 && ms != 5 && ms != 28)
      RenderPoints(flags);
    else
      RenderCrosses(flags);
  }

  glPopAttrib();
}
