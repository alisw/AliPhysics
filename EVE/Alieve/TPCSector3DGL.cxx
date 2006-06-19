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

  if(fSector->fPointSetOn) {
    UChar_t col[4];

    glPushAttrib(GL_POINT_BIT | GL_ENABLE_BIT);
    glDisable(GL_LIGHTING);
    glEnable(GL_BLEND);
    glEnable(GL_POINT_SMOOTH);
    glPointSize(fSector->GetMarkerSize());

    glPushClientAttrib(GL_CLIENT_VERTEX_ARRAY_BIT);
    glEnableClientState(GL_VERTEX_ARRAY);

    const Reve::PointSetArray& psa = fSector->fPointSetArray;
    for(Int_t b=0; b<psa.GetNBins(); ++b) {
      Reve::PointSet* ps = psa.GetBin(b);
      ColorFromIdx(ps->GetMarkerColor(), col);
      glColor4ubv(col);

      glVertexPointer(3, GL_FLOAT, 0, ps->GetP());
      glDrawArrays(GL_POINTS, 0, ps->GetN());
    }

    glPopClientAttrib();
    glPopAttrib();
  }

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
      DrawSegmentFrame(TPCSectorData::GetInnSeg(),  0, 2);
    if(fSector->fRnrOut1)
      DrawSegmentFrame(TPCSectorData::GetOut1Seg(), 2, 1);
    if(fSector->fRnrOut2)
      DrawSegmentFrame(TPCSectorData::GetOut2Seg(), 2, 2);

    glPopAttrib();
  }
}

void TPCSector3DGL::DrawSegmentFrame(const TPCSectorData::SegmentInfo& s,
                                     Int_t botExtraPads, Int_t topExtraPads) const
{
  Float_t xl, xh, yl, yh, zl, zh;
  xl = 0.5*s.GetPadWidth()*(TPCSectorData::GetNPadsInRow(s.GetFirstRow()) + botExtraPads);
  xh = 0.5*s.GetPadWidth()*(TPCSectorData::GetNPadsInRow(s.GetLastRow())  + topExtraPads);
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
