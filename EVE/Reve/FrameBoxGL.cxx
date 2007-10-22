// $Header$

#include "FrameBoxGL.h"
#include <Reve/FrameBox.h>

#include <TGLIncludes.h>

#include <TMath.h>

using namespace Reve;

//______________________________________________________________________
// FrameBoxGL
//
// A class encapsulating GL rendering of Reve::FrameBox via a static
// meber function.

ClassImp(FrameBoxGL)

void FrameBoxGL::RenderFrame(const FrameBox& b, Bool_t fillp)
{
  const Float_t*  p =  b.fFramePoints;

  if (b.fFrameType == FrameBox::FT_Quad)
  {
    glBegin(fillp ? GL_POLYGON : GL_LINE_LOOP);
    glVertex3fv(p);       glVertex3fv(p + 3);
    glVertex3fv(p + 6);   glVertex3fv(p + 9);
    glEnd();
  }
  else if (b.fFrameType == FrameBox::FT_Box)
  {
    // !!! frame-fill not implemented for 3D frame.
    glBegin(GL_LINE_STRIP);
    glVertex3fv(p);       glVertex3fv(p + 3);
    glVertex3fv(p + 6);   glVertex3fv(p + 9);
    glVertex3fv(p);
    glVertex3fv(p + 12);  glVertex3fv(p + 15);
    glVertex3fv(p + 18);  glVertex3fv(p + 21);
    glVertex3fv(p + 12);
    glEnd();
    glBegin(GL_LINES);
    glVertex3fv(p + 3);   glVertex3fv(p + 15);
    glVertex3fv(p + 6);   glVertex3fv(p + 18);
    glVertex3fv(p + 9);   glVertex3fv(p + 21);
    glEnd();
  }
}

void FrameBoxGL::Render(const FrameBox* box)
{
  const FrameBox& b = *box;

  glPushAttrib(GL_POLYGON_BIT | GL_LINE_BIT | GL_ENABLE_BIT);

  glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
  glEnable(GL_COLOR_MATERIAL);
  glDisable(GL_CULL_FACE);

  if (b.fFrameType == FrameBox::FT_Quad && b.fDrawBack)
  {
    GLboolean lmts;
    glGetBooleanv(GL_LIGHT_MODEL_TWO_SIDE, &lmts);
    if (!lmts) glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(2, 2);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    const Float_t*  p =  b.fFramePoints;
    Float_t normal[3];
    TMath::Normal2Plane(p, p+3, p+6, normal);
    glNormal3fv(normal);

    glColor4ubv(b.fBackRGBA);
    RenderFrame(b, kTRUE);

    if (!lmts) glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
  }

  glDisable(GL_LIGHTING);
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_LINE_SMOOTH);

  glLineWidth(b.fFrameWidth);
  glColor4ubv(b.fFrameRGBA);
  RenderFrame(b, b.fFrameFill);

  glPopAttrib();
}
