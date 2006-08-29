// $Header$

#include "QuadSetGL.h"
#include <Reve/QuadSet.h>

#include <TGLDrawFlags.h>

#include <GL/gl.h>
#include <GL/glu.h>

using namespace Reve;

//______________________________________________________________________
// QuadSetGL
//

ClassImp(QuadSetGL)

/**************************************************************************/

QuadSetGL::QuadSetGL() : TGLObject()
{
  // fCached = false; // Disable DL.
}

QuadSetGL::~QuadSetGL()
{}

/**************************************************************************/

Bool_t QuadSetGL::SetModel(TObject* obj)
{
#if ROOT_VERSION_CODE <= ROOT_VERSION(5,11,2)
  return set_model(obj, "Reve::QuadSet");
#elif ROOT_VERSION_CODE <= ROOT_VERSION(5,13,0)
  return SetModelCheckClass(obj, "Reve::QuadSet");
#else
  return SetModelCheckClass(obj, Reve::QuadSet::Class());
#endif
}

void QuadSetGL::SetBBox()
{
#if ROOT_VERSION_CODE <= ROOT_VERSION(5,11,2)
  set_axis_aligned_bbox(((QuadSet*)fExternalObj)->AssertBBox());
#else
  SetAxisAlignedBBox(((QuadSet*)fExternalObj)->AssertBBox());
#endif
}

/**************************************************************************/

void QuadSetGL::DirectDraw(const TGLDrawFlags & ) const
{
  // printf("QuadSetGLRenderer::DirectDraw Style %d, LOD %d\n", flags.Style(), flags.LOD());

  QuadSet& Q = * (QuadSet*) fExternalObj;

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
