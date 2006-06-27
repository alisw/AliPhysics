// $Header$

#include "BoxSetGL.h"
#include <Reve/BoxSet.h>

#include <TGLDrawFlags.h>
#ifdef WIN32
#include "Windows4root.h"
#endif
#include <GL/gl.h>
#include <GL/glu.h>

using namespace Reve;

//______________________________________________________________________
// BoxSetGL
//

ClassImp(BoxSetGL)

BoxSetGL::BoxSetGL()
{
  // fCached = false; // Disable display list.
}

BoxSetGL::~BoxSetGL()
{}

/**************************************************************************/

Bool_t BoxSetGL::SetModel(TObject* obj)
{
#if ROOT_VERSION_CODE <= ROOT_VERSION(5,11,2)
  return set_model(obj, "Reve::BoxSet");
#else
  return SetModelCheckClass(obj, "Reve::BoxSet");
#endif
}

void BoxSetGL::SetBBox()
{
#if ROOT_VERSION_CODE <= ROOT_VERSION(5,11,2)
  set_axis_aligned_bbox(((BoxSet*)fExternalObj)->AssertBBox());
#else
  SetAxisAlignedBBox(((BoxSet*)fExternalObj)->AssertBBox());
#endif
}

/**************************************************************************/

void BoxSetGL::DirectDraw(const TGLDrawFlags& /*flags*/) const
{
  BoxSet& mB = * (BoxSet*) fExternalObj;
  // printf("BoxSetGL::DirectDraw N boxes %d\n", mB.fBoxes.size());
  if(mB.fBoxes.size() == 0)
    return;

  glPushAttrib(GL_ENABLE_BIT | GL_POLYGON_BIT);
  glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
  glEnable(GL_COLOR_MATERIAL);
  glDisable(GL_CULL_FACE);

  Float_t c[4]; glGetFloatv(GL_CURRENT_COLOR, c);

  glBegin(GL_QUADS);
  for(std::vector<Box>::iterator q=mB.fBoxes.begin(); q!=mB.fBoxes.end(); ++q) {
    UChar_t* c = (UChar_t*) &q->color;
    glColor3ub(c[0], c[1], c[2]);

    // bottom: 3210
    glNormal3f(0, 0, -1);
    glVertex3fv(q->vertices + 9);  glVertex3fv(q->vertices + 6);
    glVertex3fv(q->vertices + 3);  glVertex3fv(q->vertices);
    // top:   4567
    glNormal3f(0, 0, 1);
    glVertex3fv(q->vertices + 12); glVertex3fv(q->vertices + 15);
    glVertex3fv(q->vertices + 18); glVertex3fv(q->vertices + 21);
    // front: 0154
    glNormal3f(0, -1, 0);
    glVertex3fv(q->vertices);      glVertex3fv(q->vertices + 3);
    glVertex3fv(q->vertices + 15); glVertex3fv(q->vertices + 12);
    // back:  7623
    glNormal3f(0, 1, 0);
    glVertex3fv(q->vertices + 21); glVertex3fv(q->vertices + 18);
    glVertex3fv(q->vertices + 6);  glVertex3fv(q->vertices + 9);
    // left:  4730
    glNormal3f(-1, 0, 0);
    glVertex3fv(q->vertices + 12); glVertex3fv(q->vertices + 21);
    glVertex3fv(q->vertices + 9);  glVertex3fv(q->vertices);
    // right: 5126
    glNormal3f(1, 0, 0);
    glVertex3fv(q->vertices + 15); glVertex3fv(q->vertices + 3);
    glVertex3fv(q->vertices + 6);  glVertex3fv(q->vertices + 18);
  }
  glEnd();

  glPopAttrib();
}
