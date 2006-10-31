// $Header$

#include "QuadSetGL.h"
#include <Reve/QuadSet.h>

#include <TGLDrawFlags.h>

#include <GL/gl.h>
#include <GL/glu.h>

using namespace Reve;

//______________________________________________________________________
// OldQuadSetGL
//

ClassImp(OldQuadSetGL)

/**************************************************************************/

OldQuadSetGL::OldQuadSetGL() : TGLObject()
{
  // fCached = false; // Disable DL.
}

OldQuadSetGL::~OldQuadSetGL()
{}

/**************************************************************************/

Bool_t OldQuadSetGL::SetModel(TObject* obj)
{
  return SetModelCheckClass(obj, Reve::OldQuadSet::Class());
}

void OldQuadSetGL::SetBBox()
{
  SetAxisAlignedBBox(((OldQuadSet*)fExternalObj)->AssertBBox());
}

/**************************************************************************/

void OldQuadSetGL::DirectDraw(const TGLDrawFlags & ) const
{
  // printf("OldQuadSetGLRenderer::DirectDraw Style %d, LOD %d\n", flags.Style(), flags.LOD());

  OldQuadSet& Q = * (OldQuadSet*) fExternalObj;

  glPushAttrib(GL_ENABLE_BIT | GL_POLYGON_BIT);

  glDisable(GL_LIGHTING);
  glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
  glEnable(GL_COLOR_MATERIAL);
  glPolygonMode(GL_FRONT, GL_LINE);
  glPolygonMode(GL_BACK,  GL_LINE);
  glDisable(GL_CULL_FACE);

  Float_t c[4]; glGetFloatv(GL_CURRENT_COLOR, c);
  //  UChar_t alpha = (UChar_t)(255*c[3]);

  glBegin(GL_QUADS);
  for(std::vector<Quad>::iterator q=Q.fQuads.begin(); q!=Q.fQuads.end(); ++q) {
    UChar_t* c = (UChar_t*) &q->color;
    //glColor4ub(c[0], c[1], c[2], (c[3]*alpha) >> 8);
    glColor3ub(c[0], c[1], c[2]);
    glVertex3fv(q->vertices);
    glVertex3fv(q->vertices + 3);
    glVertex3fv(q->vertices + 6);
    glVertex3fv(q->vertices + 9);
  }
  glEnd();

  glPopAttrib();
}
