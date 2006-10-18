// $Header$

#include "TriangleSetGL.h"
#include "TriangleSet.h"
#include <TVector3.h>

#include <TGLDrawFlags.h>

#ifdef WIN32
#include "Windows4root.h"
#endif
#include <GL/gl.h>

//______________________________________________________________________
// TriangleSetGL
//

using namespace Reve;

ClassImp(TriangleSetGL)

TriangleSetGL::TriangleSetGL() : TGLObject(), fM(0)
{
  // fCached = false; // Disable display list.
}

TriangleSetGL::~TriangleSetGL()
{}

/**************************************************************************/

Bool_t TriangleSetGL::SetModel(TObject* obj)
{
  if(SetModelCheckClass(obj, TriangleSet::Class())) {
    fM = dynamic_cast<TriangleSet*>(obj);
    return kTRUE;
  }
  return kFALSE;
}

void TriangleSetGL::SetBBox()
{
  // !! This ok if master sub-classed from TAttBBox
  SetAxisAlignedBBox(((TriangleSet*)fExternalObj)->AssertBBox());
}

/**************************************************************************/

void TriangleSetGL::DirectDraw(const TGLDrawFlags& /*flags*/) const
{
  TriangleSet& TS = *fM;
  Bool_t isScaled = TS.fHMTrans.IsScale();

  GLint ex_shade_model;
  glGetIntegerv(GL_SHADE_MODEL, &ex_shade_model);
  glShadeModel(GL_FLAT);

  glPushAttrib(GL_ENABLE_BIT | GL_POLYGON_BIT);

  glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
  glEnable(GL_COLOR_MATERIAL);
  glPolygonMode(GL_FRONT, GL_FILL);
  glPolygonMode(GL_BACK,  GL_LINE);
  glDisable(GL_CULL_FACE);
  if (isScaled) glEnable(GL_NORMALIZE);
  glPushClientAttrib(GL_CLIENT_VERTEX_ARRAY_BIT);
  glVertexPointer(3, GL_FLOAT, 0, TS.fVerts);
  glEnableClientState(GL_VERTEX_ARRAY);

  Int_t*   T = TS.fTrings;
  Float_t* N = TS.fTringNorms;
  UChar_t* C = TS.fTringCols;

  TVector3 e1, e2, n;

  glBegin(GL_TRIANGLES);
  for(Int_t t=0; t<TS.fNTrings; ++t) {
    if (N) {
      glNormal3fv(N); N += 3;
    } else {
      Float_t* v0 = TS.Vertex(T[0]);
      Float_t* v1 = TS.Vertex(T[1]);
      Float_t* v2 = TS.Vertex(T[2]);
      e1.SetXYZ(v1[0]-v0[0], v1[1]-v0[1], v1[2]-v0[2]);
      e2.SetXYZ(v2[0]-v0[0], v2[1]-v0[1], v2[2]-v0[2]);
      n = e1.Cross(e2);
      if (!isScaled) n.SetMag(1);
      glNormal3d(n.x(), n.y(), n.z());
    }
    if (C) {
      glColor3ubv(C);  C += 3;
    }
    glArrayElement(T[0]);
    glArrayElement(T[1]);
    glArrayElement(T[2]);
    T += 3;
  }
  glEnd();

  glPopClientAttrib();
  glPopAttrib();
  glShadeModel(ex_shade_model);
}
