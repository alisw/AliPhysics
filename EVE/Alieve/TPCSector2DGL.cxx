// $Header$

#include "TPCSector2DGL.h"

#include <Alieve/TPCData.h>

#include <TGLDrawFlags.h>
#include <GL/gl.h>

using namespace Reve;
using namespace Alieve;
using namespace std;

//______________________________________________________________________
// TPCSector2DGL
//

ClassImp(TPCSector2DGL)

// This can be optimized to non-pow-2 values once everybody has GL 1.4.

const Int_t TPCSector2DGL::fgkTextureWidth    = 256;
const Int_t TPCSector2DGL::fgkTextureHeight   = 128;
const Int_t TPCSector2DGL::fgkTextureByteSize = 4*256*128;

/**************************************************************************/

TPCSector2DGL::TPCSector2DGL() :
  TGLObject(),
  
  fSector     (0),
  fSectorData (0),

  fImage   (0),
  fTexture (0),
  fRTS     (0)
{}

TPCSector2DGL::~TPCSector2DGL()
{
  if(fImage)   delete fImage;
  if(fTexture) glDeleteTextures(1, &fTexture);
}

/**************************************************************************/

Bool_t TPCSector2DGL::SetModel(TObject* obj)
{
#if ROOT_VERSION_CODE <= ROOT_VERSION(5,11,2)
  if(set_model(obj, "Alieve::TPCSector2D")) {
#elif ROOT_VERSION_CODE <= ROOT_VERSION(5,13,0)
  if(SetModelCheckClass(obj, "Alieve::TPCSector2D")) {
#else
  if(SetModelCheckClass(obj, Alieve::TPCSector2D::Class())) {
#endif
    fSector = (TPCSector2D*) fExternalObj;
    return kTRUE;
  }
  return kFALSE;
}

void TPCSector2DGL::SetBBox()
{
#if ROOT_VERSION_CODE <= ROOT_VERSION(5,11,2)
  set_axis_aligned_bbox(((TPCSector2D*)fExternalObj)->AssertBBox());
#else
  SetAxisAlignedBBox(((TPCSector2D*)fExternalObj)->AssertBBox());
#endif
}

/**************************************************************************/

void TPCSector2DGL::ProcessSelection(UInt_t* ptr, TGLViewer*, TGLScene*)
{
  if (ptr[0] != 3) return;
  ptr += 3; // skip n, zmin, zmax
  Int_t row = ptr[1];
  Int_t pad = ptr[2];
  if (row < 0 || row >= TPCSectorData::GetNAllRows())      return;
  if (pad < 0 || pad >= TPCSectorData::GetNPadsInRow(row)) return;
  fSector->PadSelected(row, pad);
}

/**************************************************************************/

void TPCSector2DGL::DirectDraw(const TGLDrawFlags& flags) const
{
  // Actual GL drawing.

  // printf("TPCSector2DGL::DirectDraw \n");

  fSectorData = fSector->GetSectorData();

  if(fRTS < fSector->fRTS && fSectorData != 0) {
    CreateTexture();
    fRTS = fSector->fRTS;
  }

  glPushAttrib(GL_CURRENT_BIT      | GL_DEPTH_BUFFER_BIT |
               GL_COLOR_BUFFER_BIT | GL_ENABLE_BIT       | GL_POLYGON_BIT);

  glDisable(GL_LIGHTING);
  glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
  glEnable(GL_COLOR_MATERIAL);
  glDisable(GL_CULL_FACE);

  // Display digits
  if(fSectorData != 0) {

    const TPCSectorData::SegmentInfo&  iSeg = TPCSectorData::GetInnSeg();
    const TPCSectorData::SegmentInfo& o1Seg = TPCSectorData::GetOut1Seg();
    const TPCSectorData::SegmentInfo& o2Seg = TPCSectorData::GetOut2Seg();

    if(flags.SecSelection()) {

      if(fSector->fRnrInn)  DisplayNamedQuads(iSeg, 0, 0);
      if(fSector->fRnrOut1) DisplayNamedQuads(o1Seg, iSeg.GetNMaxPads(), 0);
      if(fSector->fRnrOut2) DisplayNamedQuads(o2Seg, 0, o1Seg.GetNRows());

    } else {

      if(fSector->fUseTexture) {
	//texture
	glEnable(GL_BLEND);
	glDepthMask(GL_FALSE);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glPolygonOffset(2,2);
	glEnable(GL_POLYGON_OFFSET_FILL);

	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glBindTexture(GL_TEXTURE_2D, fTexture);
	glEnable(GL_TEXTURE_2D);

	if(fSector->fRnrInn)  DisplayTexture(iSeg, 0, 0);
	if(fSector->fRnrOut1) DisplayTexture(o1Seg, iSeg.GetNMaxPads(), 0);
	if(fSector->fRnrOut2) DisplayTexture(o2Seg, 0, o1Seg.GetNRows());

	glDisable(GL_TEXTURE_2D);
      } else {
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

/**************************************************************************/
// Data import
/**************************************************************************/

void TPCSector2DGL::LoadPadrow(TPCSectorData::RowIterator& iter,
			       Int_t row, Int_t col_off) const
{
  Int_t    padVal;
  Int_t    time, val;   

  Int_t    minTime = fSector->fMinTime;
  Int_t    maxTime = fSector->fMaxTime;
  Bool_t   halfBorderTime = ((maxTime - minTime) % 2 == 0);

  UChar_t* img_pos = GetRowCol(row, col_off);
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

    if(fSector->fShowMax == kFALSE && fSector->fAverage) {
      padVal = (Int_t)((Float_t)padVal / (maxTime - minTime));
    }
    padVal = TMath::Min(padVal, fSector->fMaxVal);
    if(padVal > fSector->fThreshold)
      fSector->ColorFromArray(padVal, img_pos);
    img_pos += 4;
  }
}

/**************************************************************************/

void TPCSector2DGL::CreateTexture() const
{
  if (fImage == 0 ) {
    fImage = new UChar_t[fgkTextureByteSize];
    glGenTextures(1, &fTexture);
  }
  memset(fImage, 0, fgkTextureByteSize);

  Int_t rowOff[3], colOff[3];
  Bool_t isOn[3];
  rowOff[0] = 0;
  rowOff[1] = rowOff[2] = -TPCSectorData::GetSeg(1).GetFirstRow();
  colOff[0] = colOff[2] = 0;
  colOff[1] = TPCSectorData::GetSeg(0).GetNMaxPads();
  isOn[0] = fSector->fRnrInn;
  isOn[1] = fSector->fRnrOut1;
  isOn[2] = fSector->fRnrOut2;

  fSector->SetupColorArray();

  // Loop over 3 main segments
  for (Int_t sId = 0; sId <= 2; ++sId) {
    if(isOn[sId] == kFALSE)
      continue;
    const TPCSectorData::SegmentInfo& sInfo = TPCSectorData::GetSeg(sId);
    for (Int_t row=sInfo.GetFirstRow(); row<=sInfo.GetLastRow(); ++row) {
      TPCSectorData::RowIterator i = fSectorData->MakeRowIterator(row);
      Int_t offset = (sInfo.GetNMaxPads() - TPCSectorData::GetNPadsInRow(row))/2;
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

/**************************************************************************/
// Data display
/**************************************************************************/

void TPCSector2DGL::DisplayTexture(const TPCSectorData::SegmentInfo& seg,
                                   Int_t startCol, Int_t startRow) const
{
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

/**************************************************************************/

void TPCSector2DGL::DisplayQuads(const TPCSectorData::SegmentInfo& seg,
                                 Int_t startCol, Int_t startRow) const
{
  Float_t y_d, y_u;
  Float_t x_off, x;
  Float_t padW = seg.GetPadWidth();
  Float_t padH = seg.GetPadHeight();

  glBegin(GL_QUADS);
  for (Int_t row=0; row<seg.GetNRows(); row++) {
    y_d = seg.GetRLow() + row*padH;
    y_u = y_d + padH;
    x_off = -seg.GetNMaxPads()*padW/2;
    Int_t tpcRow = row + seg.GetFirstRow();
    Int_t deltaPad = (seg.GetNMaxPads() - TPCSectorData::GetNPadsInRow(tpcRow))/2;
    Int_t   maxPad = seg.GetNMaxPads() - deltaPad;
    UChar_t   *pix = GetRowCol(row + startRow, startCol + deltaPad);
    for (Int_t pad=deltaPad; pad<maxPad; pad++, pix+=4) {
      x = x_off + pad*padW;
      if (pix[3] != 0) {
        glColor4ubv(pix);
        glVertex2f(x+padW, y_d);
        glVertex2f(x,      y_d);
        glVertex2f(x,      y_u);
        glVertex2f(x+padW, y_u);
      }
    }
  }
  glEnd();
}

void TPCSector2DGL::DisplayNamedQuads(const TPCSectorData::SegmentInfo& seg,
				      Int_t startCol, Int_t startRow) const
{
  Float_t y_d, y_u;
  Float_t x_off, x;
  Float_t padW = seg.GetPadWidth();
  Float_t padH = seg.GetPadHeight();

  glPushName(0);
  for (Int_t row=0; row<seg.GetNRows(); row++) {
    y_d = seg.GetRLow() + row*padH;
    y_u = y_d + padH;
    x_off = -seg.GetNMaxPads()*padW/2;
    Int_t tpcRow = row + seg.GetFirstRow();
    glLoadName(tpcRow);
    Int_t deltaPad = (seg.GetNMaxPads() - TPCSectorData::GetNPadsInRow(tpcRow))/2;
    Int_t   maxPad = seg.GetNMaxPads() - deltaPad;
    UChar_t   *pix = GetRowCol(row + startRow, startCol + deltaPad);
    glPushName(0);
    for (Int_t pad=deltaPad; pad<maxPad; pad++, pix+=4) {
      x = x_off + pad*padW;
      if (pix[3] != 0 || fSector->fPickEmpty) {
	glLoadName(pad - deltaPad);
	glBegin(GL_QUADS);
        glVertex2f(x+padW, y_d);
        glVertex2f(x,      y_d);
        glVertex2f(x,      y_u);
        glVertex2f(x+padW, y_u);
	glEnd();
      }
    }
    glPopName();
  }
  glPopName();
}

/**************************************************************************/
// Frame drawing
/**************************************************************************/

void TPCSector2DGL::TraceStepsUp(const TPCSectorData::SegmentInfo& s)
{
  Float_t x = -(s.GetNMaxPads()*1.0/2 - s.GetNYSteps())*s.GetPadWidth();
  Float_t y  = s.GetRLow();
  glVertex2f(x, y);
  for (Int_t i=0; i<s.GetNYSteps(); ++i) {
    y = s.GetYStep(i);
    glVertex2f(x, y);
    x -= s.GetPadWidth();
    glVertex2f(x, y);
  }
  y =  s.GetRLow() + s.GetNRows()*s.GetPadHeight();
  glVertex2f(-s.GetNMaxPads()*s.GetPadWidth()/2, y);
}

void TPCSector2DGL::TraceStepsDown(const TPCSectorData::SegmentInfo& s) 
{
  Float_t x = s.GetNMaxPads()*s.GetPadWidth()/2;
  Float_t y = s.GetRLow() + s.GetNRows()*s.GetPadHeight();
  glVertex2f(x, y);
  for (Int_t i=s.GetNYSteps() - 1; i>=0; --i) {
    y =  s.GetYStep(i);
    glVertex2f(x, y);
    x -= s.GetPadWidth();
    glVertex2f(x, y);
  }
  y = s.GetRLow();
  glVertex2f((0.5*s.GetNMaxPads() - s.GetNYSteps())*s.GetPadWidth(), y);
}

void TPCSector2DGL::DisplayFrame() const
{
  UChar_t col[4];
  ColorFromIdx(fSector->fFrameColor, col);
  glColor4ubv(col);

  if(fSector->fRnrInn) {
    glBegin(GL_POLYGON);
    TraceStepsUp  (TPCSectorData::GetInnSeg());
    TraceStepsDown(TPCSectorData::GetInnSeg());
    glEnd();
  }
  if(fSector->fRnrOut1) {
    glBegin(GL_POLYGON);
    TraceStepsUp  (TPCSectorData::GetOut1Seg());
    TraceStepsDown(TPCSectorData::GetOut1Seg());
    glEnd();
  }
  if(fSector->fRnrOut2) {
    glBegin(GL_POLYGON);
    TraceStepsUp  (TPCSectorData::GetOut2Seg());
    TraceStepsDown(TPCSectorData::GetOut2Seg());
    glEnd();
  }
}
