// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveTPCSector2DGL.h"

#include <EveDet/AliEveTPCSector2D.h>

#include <TGLRnrCtx.h>
#include <TGLSelectRecord.h>
#include <TGLIncludes.h>

//______________________________________________________________________________
//
// GL renderer for AliEveTPCSector2D.

ClassImp(AliEveTPCSector2DGL)

// This can be optimized to non-pow-2 values once everybody has GL 1.4.

const Int_t AliEveTPCSector2DGL::fgkTextureWidth    = 256;
const Int_t AliEveTPCSector2DGL::fgkTextureHeight   = 128;
const Int_t AliEveTPCSector2DGL::fgkTextureByteSize = 4*256*128;

/******************************************************************************/

AliEveTPCSector2DGL::AliEveTPCSector2DGL() :
  TGLObject(),

  fSector     (0),
  fSectorData (0),

  fImage   (0),
  fTexture (0),
  fRTS     (0)
{
  // Constructor.
}

AliEveTPCSector2DGL::~AliEveTPCSector2DGL()
{
  // Destructor.
  // !!!! Should unregister texture via ContextIdentity!

  if (fImage)   delete [] fImage;
  if (fTexture) glDeleteTextures(1, &fTexture);
}

/******************************************************************************/

Bool_t AliEveTPCSector2DGL::SetModel(TObject* obj, const Option_t* /*opt*/)
{
  // Set model object.

  if (SetModelCheckClass(obj, AliEveTPCSector2D::Class())) {
    fSector = (AliEveTPCSector2D*) fExternalObj;
    return kTRUE;
  }
  return kFALSE;
}

void AliEveTPCSector2DGL::SetBBox()
{
  // Set bounding-box.

  SetAxisAlignedBBox(((AliEveTPCSector2D*)fExternalObj)->AssertBBox());
}

/******************************************************************************/

void AliEveTPCSector2DGL::ProcessSelection(TGLRnrCtx       & /*rnrCtx*/,
                                           TGLSelectRecord & rec)
{
  // Process selection record.
  // Determine row and pad, call PadSelected() in model object.

  if (rec.GetN() != 3) return;
  Int_t row = rec.GetItem(1);
  Int_t pad = rec.GetItem(2);
  if (row < 0 || row >= AliEveTPCSectorData::GetNAllRows())      return;
  if (pad < 0 || pad >= AliEveTPCSectorData::GetNPadsInRow(row)) return;
  fSector->PadSelected(row, pad);
}

/******************************************************************************/

void AliEveTPCSector2DGL::DirectDraw(TGLRnrCtx& rnrCtx) const
{
  // Actual GL drawing.

  // printf("AliEveTPCSector2DGL::DirectDraw \n");

  fSectorData = fSector->GetSectorData();

  if (fRTS < fSector->fRTS && fSectorData != 0) {
    CreateTexture();
    fRTS = fSector->fRTS;
  }

  glPushAttrib(GL_CURRENT_BIT | GL_COLOR_BUFFER_BIT | GL_POLYGON_BIT | GL_ENABLE_BIT);

  glDisable(GL_LIGHTING);
  glDisable(GL_CULL_FACE);

  // Display digits
  if (fSectorData != 0)
  {
    const AliEveTPCSectorData::SegmentInfo&  iSeg = AliEveTPCSectorData::GetInnSeg();
    const AliEveTPCSectorData::SegmentInfo& o1Seg = AliEveTPCSectorData::GetOut1Seg();
    const AliEveTPCSectorData::SegmentInfo& o2Seg = AliEveTPCSectorData::GetOut2Seg();

    if (rnrCtx.SecSelection())
    {

      if(fSector->fRnrInn)  DisplayNamedQuads(iSeg, 0, 0);
      if(fSector->fRnrOut1) DisplayNamedQuads(o1Seg, iSeg.GetNMaxPads(), 0);
      if(fSector->fRnrOut2) DisplayNamedQuads(o2Seg, 0, o1Seg.GetNRows());
    }
    else
    {
      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
      if (fSector->fUseTexture)
      {
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_BLEND);

	glAlphaFunc(GL_GREATER, 0.2);
	glEnable(GL_ALPHA_TEST);

	glBindTexture(GL_TEXTURE_2D, fTexture);
	glEnable(GL_TEXTURE_2D);

	if(fSector->fRnrInn)  DisplayTexture(iSeg, 0, 0);
	if(fSector->fRnrOut1) DisplayTexture(o1Seg, iSeg.GetNMaxPads(), 0);
	if(fSector->fRnrOut2) DisplayTexture(o2Seg, 0, o1Seg.GetNRows());

	glDisable(GL_TEXTURE_2D);
      }
      else
      {
	if(fSector->fRnrInn)  DisplayQuads(iSeg, 0, 0);
	if(fSector->fRnrOut1) DisplayQuads(o1Seg, iSeg.GetNMaxPads(), 0);
	if(fSector->fRnrOut2) DisplayQuads(o2Seg, 0, o1Seg.GetNRows());
      }
      glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
      DisplayFrame();
    }
  }

  glPopAttrib();
}

/******************************************************************************/
// Data import
/******************************************************************************/

void AliEveTPCSector2DGL::LoadPadrow(AliEveTPCSectorData::RowIterator& iter,
                                     Int_t row, Int_t colOff) const
{
  // Load data for one pad-row into the texture.

  Int_t    padVal;
  Int_t    time, val;

  Int_t    minTime = fSector->fMinTime;
  Int_t    maxTime = fSector->fMaxTime;
  Bool_t   halfBorderTime = ((maxTime - minTime) % 2 == 0);

  UChar_t* imgPos = GetRowCol(row, colOff);
  while (iter.NextPad()) {
    padVal = 0;

    while (iter.Next()) {
      time = iter.Time();
      val  = iter.Signal();

      if(time < minTime || time > maxTime)
	continue;

      if(fSector->fShowMax) {
        if(val > padVal) {
          padVal = val;
        }
      } else {
	if(halfBorderTime && (time == minTime || time == maxTime))
	  padVal += val/2;
	else
	  padVal += val;
      }
    }

    if (fSector->fShowMax == kFALSE && fSector->fAverage) {
      padVal = (Int_t)((Float_t)padVal / (maxTime - minTime));
    }
    padVal = TMath::Min(padVal, fSector->fMaxVal);
    if (padVal > fSector->fThreshold) {
      fSector->ColorFromArray(padVal, imgPos);
    }
    imgPos += 4;
  }
}

/******************************************************************************/

void AliEveTPCSector2DGL::CreateTexture() const
{
  // Create texture that holds pad data.

  if (fImage == 0)
  {
    fImage = new UChar_t[fgkTextureByteSize];
    glGenTextures(1, &fTexture);
  }
  memset(fImage, 0, fgkTextureByteSize);

  Int_t rowOff[3], colOff[3];
  Bool_t isOn[3];
  rowOff[0] = 0;
  rowOff[1] = rowOff[2] = -AliEveTPCSectorData::GetSeg(1).GetFirstRow();
  colOff[0] = colOff[2] = 0;
  colOff[1] = AliEveTPCSectorData::GetSeg(0).GetNMaxPads();
  isOn[0] = fSector->fRnrInn;
  isOn[1] = fSector->fRnrOut1;
  isOn[2] = fSector->fRnrOut2;

  fSector->SetupColorArray();

  // Loop over 3 main segments
  for (Int_t sId = 0; sId <= 2; ++sId)
  {
    if (isOn[sId] == kFALSE)
      continue;
    const AliEveTPCSectorData::SegmentInfo& sInfo = AliEveTPCSectorData::GetSeg(sId);
    for (Int_t row = sInfo.GetFirstRow(); row <= sInfo.GetLastRow(); ++row)
    {
      AliEveTPCSectorData::RowIterator i = fSectorData->MakeRowIterator(row);
      Int_t offset = (sInfo.GetNMaxPads() - AliEveTPCSectorData::GetNPadsInRow(row))/2;
      LoadPadrow(i, row + rowOff[sId], offset + colOff[sId]);
    }
  }

  glBindTexture  (GL_TEXTURE_2D, fTexture);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S,GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T,GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexEnvf      (GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,   GL_REPLACE);
  glTexImage2D   (GL_TEXTURE_2D, 0, GL_RGBA, fgkTextureWidth, fgkTextureHeight,
                  0, GL_RGBA, GL_UNSIGNED_BYTE, fImage);

}

/******************************************************************************/
// Data display
/******************************************************************************/

void AliEveTPCSector2DGL::DisplayTexture(const AliEveTPCSectorData::SegmentInfo& seg,
                                         Int_t startCol, Int_t startRow) const
{
  // Display segment data via one textured rectangle.

  Float_t w  = seg.GetNMaxPads()*seg.GetPadWidth()/2;
  Float_t y1 = seg.GetRLow();
  Float_t y2 = y1 + seg.GetNRows()*seg.GetPadHeight();

  Float_t u1 = (Float_t) startCol / fgkTextureWidth;
  Float_t v1 = (Float_t) startRow / fgkTextureHeight;
  Float_t u2 = u1 + (Float_t) seg.GetNMaxPads() / fgkTextureWidth;
  Float_t v2 = v1 + (Float_t) seg.GetNRows()    / fgkTextureHeight;

  glBegin(GL_QUADS);
  glTexCoord2f(u1, v1);  glVertex2f(-w, y1);
  glTexCoord2f(u1, v2);  glVertex2f(-w, y2);
  glTexCoord2f(u2, v2);  glVertex2f( w, y2);
  glTexCoord2f(u2, v1);  glVertex2f( w, y1);
  glEnd();
}

/******************************************************************************/

void AliEveTPCSector2DGL::DisplayQuads(const AliEveTPCSectorData::SegmentInfo& seg,
                                       Int_t startCol, Int_t startRow) const
{
  // Display segment data by rendering one quad per pad.

  Float_t yD, yU;
  Float_t xOff, x;
  Float_t padW = seg.GetPadWidth();
  Float_t padH = seg.GetPadHeight();

  glBegin(GL_QUADS);
  for (Int_t row = 0; row < seg.GetNRows(); row++)
  {
    yD = seg.GetRLow() + row*padH;
    yU = yD + padH;
    xOff = -seg.GetNMaxPads()*padW/2;
    Int_t tpcRow = row + seg.GetFirstRow();
    Int_t deltaPad = (seg.GetNMaxPads() - AliEveTPCSectorData::GetNPadsInRow(tpcRow))/2;
    Int_t   maxPad = seg.GetNMaxPads() - deltaPad;
    UChar_t   *pix = GetRowCol(row + startRow, startCol + deltaPad);
    for (Int_t pad = deltaPad; pad < maxPad; pad++, pix+=4)
    {
      x = xOff + pad*padW;
      if (pix[3] != 0)
      {
        TGLUtil::Color4ubv(pix);
        glVertex2f(x+padW, yD);
        glVertex2f(x,      yD);
        glVertex2f(x,      yU);
        glVertex2f(x+padW, yU);
      }
    }
  }
  glEnd();
}

void AliEveTPCSector2DGL::DisplayNamedQuads(const AliEveTPCSectorData::SegmentInfo& seg,
                                            Int_t startCol, Int_t startRow) const
{
  // Display segmen data as one quad per pad.
  // Tag the rows and pads for selection.

  Float_t yD, yU;
  Float_t xOff, x;
  Float_t padW = seg.GetPadWidth();
  Float_t padH = seg.GetPadHeight();

  glPushName(0);
  for (Int_t row = 0; row < seg.GetNRows(); row++)
  {
    yD = seg.GetRLow() + row*padH;
    yU = yD + padH;
    xOff = -seg.GetNMaxPads()*padW/2;
    Int_t tpcRow = row + seg.GetFirstRow();
    glLoadName(tpcRow);
    Int_t deltaPad = (seg.GetNMaxPads() - AliEveTPCSectorData::GetNPadsInRow(tpcRow))/2;
    Int_t   maxPad = seg.GetNMaxPads() - deltaPad;
    UChar_t   *pix = GetRowCol(row + startRow, startCol + deltaPad);
    glPushName(0);
    for (Int_t pad = deltaPad; pad < maxPad; pad++, pix+=4) 
    {
      x = xOff + pad*padW;
      if (pix[3] != 0 || fSector->fPickEmpty)
      {
	glLoadName(pad - deltaPad);
	glBegin(GL_QUADS);
        glVertex2f(x+padW, yD);
        glVertex2f(x,      yD);
        glVertex2f(x,      yU);
        glVertex2f(x+padW, yU);
	glEnd();
      }
    }
    glPopName();
  }
  glPopName();
}

/******************************************************************************/
// Frame drawing
/******************************************************************************/

void AliEveTPCSector2DGL::TraceStepsUp(const AliEveTPCSectorData::SegmentInfo& s)
{
  // Trace border of segment upwards.

  Float_t x = -(s.GetNMaxPads()*1.0/2 - s.GetNYSteps())*s.GetPadWidth();
  Float_t y  = s.GetRLow();
  glVertex2f(x, y);
  for (Int_t i = 0; i < s.GetNYSteps(); ++i)
  {
    y = s.GetYStep(i);
    glVertex2f(x, y);
    x -= s.GetPadWidth();
    glVertex2f(x, y);
  }
  y =  s.GetRLow() + s.GetNRows()*s.GetPadHeight();
  glVertex2f(-s.GetNMaxPads()*s.GetPadWidth()/2, y);
}

void AliEveTPCSector2DGL::TraceStepsDown(const AliEveTPCSectorData::SegmentInfo& s)
{
  // Trace border of segment downwards.

  Float_t x = s.GetNMaxPads()*s.GetPadWidth()/2;
  Float_t y = s.GetRLow() + s.GetNRows()*s.GetPadHeight();
  glVertex2f(x, y);
  for (Int_t i = s.GetNYSteps() - 1; i >= 0; --i)
  {
    y =  s.GetYStep(i);
    glVertex2f(x, y);
    x -= s.GetPadWidth();
    glVertex2f(x, y);
  }
  y = s.GetRLow();
  glVertex2f((0.5*s.GetNMaxPads() - s.GetNYSteps())*s.GetPadWidth(), y);
}

void AliEveTPCSector2DGL::DisplayFrame() const
{
  // Display frame of the sector.
  // Each segment's frame is drawn only if its data is drawn, too.

  TGLUtil::Color(fSector->fFrameColor);

  if(fSector->fRnrInn)
  {
    glBegin(GL_POLYGON);
    TraceStepsUp  (AliEveTPCSectorData::GetInnSeg());
    TraceStepsDown(AliEveTPCSectorData::GetInnSeg());
    glEnd();
  }
  if(fSector->fRnrOut1)
  {
    glBegin(GL_POLYGON);
    TraceStepsUp  (AliEveTPCSectorData::GetOut1Seg());
    TraceStepsDown(AliEveTPCSectorData::GetOut1Seg());
    glEnd();
  }
  if(fSector->fRnrOut2)
  {
    glBegin(GL_POLYGON);
    TraceStepsUp  (AliEveTPCSectorData::GetOut2Seg());
    TraceStepsDown(AliEveTPCSectorData::GetOut2Seg());
    glEnd();
  }
}
