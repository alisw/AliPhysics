// $Header$

#include "FrameBoxGL.h"
#include <Reve/FrameBox.h>

#include <TGLDrawFlags.h>

#ifdef WIN32
#include "Windows4root.h"
#endif
#include <GL/gl.h>
#include <GL/glu.h>

using namespace Reve;

//______________________________________________________________________
// FrameBoxGL
//
// A class encapsulating GL rendering of Reve::FrameBox via a static
// meber function.

ClassImp(FrameBoxGL)

void FrameBoxGL::Render(const FrameBox* box)
{
  GLboolean lightp;
  glGetBooleanv(GL_LIGHTING, &lightp);
  if (lightp) glDisable(GL_LIGHTING);

  const FrameBox& b = *box;
  const Float_t*  p =  b.fFramePoints;
  glColor4ubv(b.fFrameRGBA);
  if (b.fFrameType == FrameBox::FT_Quad)
  {
    glBegin(GL_LINE_LOOP);
    glVertex3fv(p);       glVertex3fv(p + 3);
    glVertex3fv(p + 6);   glVertex3fv(p + 9);
    glEnd();
  }
  else if (b.fFrameType == FrameBox::FT_Box)
  {
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

  if (lightp) glEnable(GL_LIGHTING);
}
