// $Header$

#include "TPCSector3DGL.h"
#include <Alieve/TPCSector3D.h>

#include <Reve/BoxSetGL.h>

#include <TGLDrawFlags.h>

#ifdef WIN32
#include "Windows4root.h"
#endif
#include <GL/gl.h>
#include <GL/glu.h>

using namespace Reve;
using namespace Alieve;

//______________________________________________________________________
// TPCSector3DGL
//

ClassImp(TPCSector3DGL)

  TPCSector3DGL::TPCSector3DGL() : fSector(0), fBoxRnr(0)
{
  // fCached = false; // Disable display list.
}

TPCSector3DGL::~TPCSector3DGL()
{
  delete fBoxRnr;
}

/**************************************************************************/

Bool_t TPCSector3DGL::SetModel(TObject* obj)
{
#if ROOT_VERSION_CODE <= ROOT_VERSION(5,11,2)
  if(set_model(obj, "Alieve::TPCSector3D")) {
#else
  if(SetModelCheckClass(obj, "Alieve::TPCSector3D")) {
#endif
    fSector = (TPCSector3D*) fExternalObj;
    if(fBoxRnr == 0) {
      fBoxRnr = new BoxSetGL;
      fBoxRnr->SetModel(&fSector->fBoxSet);
    }
    return kTRUE;
  }
  return kFALSE;
}

void TPCSector3DGL::SetBBox()
{
#if ROOT_VERSION_CODE <= ROOT_VERSION(5,11,2)
  set_axis_aligned_bbox(((TPCSector3D*)fExternalObj)->AssertBBox());
#else
  SetAxisAlignedBBox(((TPCSector3D*)fExternalObj)->AssertBBox());
#endif
}

/**************************************************************************/

void TPCSector3DGL::DirectDraw(const TGLDrawFlags & flags) const
{
  // printf("TPCSector3DGL::DirectDraw Style %d, LOD %d\n", flags.Style(), flags.LOD());

  fBoxRnr->Render(flags);

  if(fSector->fRnrFrame) {
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    glPushAttrib(GL_CURRENT_BIT | GL_ENABLE_BIT);

    glDisable(GL_LIGHTING);
    glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);

    UChar_t col[4];
    ColorFromIdx(fSector->fFrameColor, col);
    glColor4ubv(col);

    if(fSector->fRnrInn)
      DrawSegmentFrame(TPCSectorData::GetInnSeg());
    if(fSector->fRnrOut1)
      DrawSegmentFrame(TPCSectorData::GetOut1Seg());
    if(fSector->fRnrOut2)
      DrawSegmentFrame(TPCSectorData::GetOut2Seg());

    glPopAttrib();
  }
}

void TPCSector3DGL::DrawSegmentFrame(const TPCSectorData::SegmentInfo& s) const
{
  Float_t xl, xh, yl, yh, zl, zh;
  xl = 0.5 * TPCSectorData::GetNPadsInRow(s.GetFirstRow()) * s.GetPadWidth();
  xh = 0.5 * TPCSectorData::GetNPadsInRow(s.GetLastRow())  * s.GetPadWidth();
  yl = s.GetRLow();
  yh = yl + s.GetNRows()*s.GetPadHeight();
  zl = -0.5;
  zh = 250.5;

  glBegin(GL_LINE_LOOP);
  glVertex3f( xl, yl, zl);  glVertex3f( xh, yh, zl);
  glVertex3f(-xh, yh, zl);  glVertex3f(-xl, yl, zl);
  glEnd();
  glBegin(GL_LINE_LOOP);
  glVertex3f( xl, yl, zh);  glVertex3f( xh, yh, zh);
  glVertex3f(-xh, yh, zh);  glVertex3f(-xl, yl, zh);
  glEnd();
  glBegin(GL_LINES);
  glVertex3f( xl, yl, zl);  glVertex3f( xl, yl, zh);
  glVertex3f( xh, yh, zl);  glVertex3f( xh, yh, zh);
  glVertex3f(-xh, yh, zl);  glVertex3f(-xh, yh, zh);
  glVertex3f(-xl, yl, zl);  glVertex3f(-xl, yl, zh);
  glEnd();
}
