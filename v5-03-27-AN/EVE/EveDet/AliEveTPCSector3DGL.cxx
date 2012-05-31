// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveTPCSector3DGL.h"
#include <EveDet/AliEveTPCSector3D.h>

#include <TEveBoxSetGL.h>

#include <TGLIncludes.h>
#include <TGLRnrCtx.h>
#include <TGLSelectRecord.h>

//______________________________________________________________________________
//
// GL renderer for AliEveTPCSector3D.

ClassImp(AliEveTPCSector3DGL)

AliEveTPCSector3DGL::AliEveTPCSector3DGL() :
  TGLObject(),
  fSector(0), fBoxRnr(0),
  fRTS(0)
{
  // Constructor.

  fDLCache = false; // Disable display list.
}

AliEveTPCSector3DGL::~AliEveTPCSector3DGL()
{
  // Destructor.

  delete fBoxRnr;
}

/******************************************************************************/

//______________________________________________________________________________
Short_t AliEveTPCSector3DGL::QuantizeShapeLOD(Short_t shapeLOD, Short_t combiLOD) const
{
   // Factor in scene/viewer LOD and quantize.

   Int_t lod = ((Int_t)shapeLOD * (Int_t)combiLOD) / 100;

   if (lod >= 100)
     return 100;
   else
     return (Short_t)(10 * TMath::Nint(0.1*lod));
}

/******************************************************************************/

Bool_t AliEveTPCSector3DGL::SetModel(TObject* obj, const Option_t* /*opt*/)
{
  // Set model object.

  if(SetModelCheckClass(obj, AliEveTPCSector3D::Class())) {
    fSector = (AliEveTPCSector3D*) fExternalObj;
    if(fBoxRnr == 0) {
      fBoxRnr = new TEveBoxSetGL;
      fBoxRnr->SetModel(&fSector->fBoxSet);
    }
    return kTRUE;
  }
  return kFALSE;
}

void AliEveTPCSector3DGL::SetBBox()
{
  // Set bounding-box.

  SetAxisAlignedBBox(((AliEveTPCSector3D*)fExternalObj)->AssertBBox());
}

/******************************************************************************/

void AliEveTPCSector3DGL::DirectDraw(TGLRnrCtx & rnrCtx) const
{
  // Render object.

  // printf("AliEveTPCSector3DGL::DirectDraw Style %d, LOD %d\n", rnrCtx.Style(), rnrCtx.LOD());

  if(fRTS < fSector->fRTS) {
    fSector->UpdateBoxesAndPoints();
    fRTS = fSector->fRTS;
  }

  if (rnrCtx.SecSelection()) glPushName(0);

  Bool_t hasData = (fSector->GetSectorData() != 0);

  if(hasData)
  {
    if (rnrCtx.SecSelection()) glLoadName(9999);
    fBoxRnr->Render(rnrCtx);
  }

  glPushAttrib(GL_CURRENT_BIT | GL_POINT_BIT | GL_ENABLE_BIT);
  glDisable(GL_LIGHTING);

  if(hasData && fSector->fPointSetOn)
  {
    glEnable(GL_BLEND);
    glEnable(GL_POINT_SMOOTH);
    glPointSize(fSector->fPointSize);

    glPushClientAttrib(GL_CLIENT_VERTEX_ARRAY_BIT);
    glEnableClientState(GL_VERTEX_ARRAY);

    const TEvePointSetArray& psa = fSector->fPointSetArray;
    for(Int_t b=0; b<psa.GetNBins(); ++b)
    {
      TEvePointSet* ps = psa.GetBin(b);
      if(ps->Size() > 0)
      {
	TGLUtil::Color(ps->GetMarkerColor());

        if (rnrCtx.SecSelection()) glLoadName(b + 1);
	glVertexPointer(3, GL_FLOAT, 0, ps->GetP());
	glDrawArrays(GL_POINTS, 0, ps->Size());
      }
    }

    glPopClientAttrib();
  }

  if(fSector->fRnrFrame && ! rnrCtx.SecSelection())
  {
    TGLUtil::Color(fSector->fFrameColor);

    if(fSector->fRnrInn)
      DrawSegmentFrame(AliEveTPCSectorData::GetInnSeg(),  0, 2);
    if(fSector->fRnrOut1)
      DrawSegmentFrame(AliEveTPCSectorData::GetOut1Seg(), 2, 1);
    if(fSector->fRnrOut2)
      DrawSegmentFrame(AliEveTPCSectorData::GetOut2Seg(), 2, 2);
  }

  glPopAttrib();
}

void AliEveTPCSector3DGL::DrawSegmentFrame(const AliEveTPCSectorData::SegmentInfo& s,
                                           Int_t botExtraPads, Int_t topExtraPads) const
{
  // Draw frame of given segment.

  Float_t xl, xh, yl, yh, zl, zh;
  xl = 0.5*s.GetPadWidth()*(AliEveTPCSectorData::GetNPadsInRow(s.GetFirstRow()) + botExtraPads);
  xh = 0.5*s.GetPadWidth()*(AliEveTPCSectorData::GetNPadsInRow(s.GetLastRow())  + topExtraPads);
  yl = s.GetRLow();
  yh = yl + s.GetNRows()*s.GetPadHeight();
  zl = 0;
  zh = AliEveTPCSectorData::GetZLength();

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

/******************************************************************************/
/******************************************************************************/

//______________________________________________________________________________
void AliEveTPCSector3DGL::ProcessSelection(TGLRnrCtx & /*rnrCtx*/, TGLSelectRecord & rec)
{
  // Processes secondary selection from TGLViewer.
  // Calls TPointSet3D::PointSelected(Int_t) with index of selected
  // point as an argument.

  if (rec.GetN() < 3) return;

  if (rec.GetItem(1) == 9999)
  {
    printf("TPC3D Box selected idx=%u\n", rec.GetItem(2));
    return;
  }

  const TEvePointSetArray& psa = fSector->fPointSetArray;

  if (rec.GetItem(1) > 0 && rec.GetItem(1) <= (UInt_t) psa.GetNBins())
  {
    // TEvePointSet& ps = * psa.GetBin(rec.GetItem(1) - 1);
    printf("TPC3D Point selected, bin=%u, idx=%u\n", rec.GetItem(1) - 1, rec.GetItem(2));
    return;
  }
}
