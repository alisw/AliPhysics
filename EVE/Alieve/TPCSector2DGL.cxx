// $Header$

#include "TPCSector2DGL.h"

#include <Alieve/TPCData.h>

#include <TStyle.h>
#include <TColor.h>
#include <TStopwatch.h>

#include <GL/gl.h>

using namespace Reve;
using namespace Alieve;
using namespace std;

  // This can be optimized to non-pow-2 values once everybody has GL 1.4.

const Int_t TPCSector2DGL::fgkTextureWidth    = 256;
const Int_t TPCSector2DGL::fgkTextureHeight   = 128;
const Int_t TPCSector2DGL::fgkTextureByteSize = 4*256*128;

/**************************************************************************/

TPCSector2DGL::TPCSector2DGL() : TGLObject()
{
  fSector     = 0;
  fSectorData = 0;

  fImage   = 0;
  fTexture = 0;
  fRTS     = 0;
}

TPCSector2DGL::~TPCSector2DGL()
{
  if(fImage)   delete fImage;
  if(fTexture) glDeleteTextures(1, &fTexture);
}

/**************************************************************************/

Bool_t TPCSector2DGL::SetModel(TObject* obj)
{
  if (set_model(obj, "Alieve::TPCSector2D")) {
    fSector = (TPCSector2D*) fExternalObj;
    return true;
  }
  return false;
}

void TPCSector2DGL::SetBBox()
{
  set_axis_aligned_bbox(((TPCSector2D*)fExternalObj)->AssertBBox());
}

/**************************************************************************/

void TPCSector2DGL::SetCol(Float_t z, UChar_t* pixel) const
{
 
  Int_t n_col = gStyle->GetNumberOfColors();

  Int_t ci = gStyle->GetColorPalette
    (TMath::Min(n_col - 1,
                Int_t((n_col*(z - fSector->fthreshold))/(fSector->fMaxVal - fSector->fthreshold))));

  TColor* c = gROOT->GetColor(ci);

  if(c) {
    //    UChar_t *x = (UChar_t*) &c;
    pixel[0] = (UChar_t)(255*c->GetRed());
    pixel[1] = (UChar_t)(255*c->GetGreen());
    pixel[2] = (UChar_t)(255*c->GetBlue());
    pixel[3] = 255;
  }
}

/**************************************************************************/

void TPCSector2DGL::DirectDraw(const TGLDrawFlags& /*flags*/) const
{
  // Actual GL drawing.

  // printf("TPCSector2DGL::DirectDraw \n");

  if(fSector->fTPCData == 0)
    fSectorData = 0;
  else
    fSectorData = fSector->fTPCData->GetSectorData(fSector->fSectorID);

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

    if(fSector->fUseTexture) {
      //texture
      glEnable(GL_BLEND);
      glDepthMask(GL_FALSE);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      glShadeModel(GL_FLAT);

      glBindTexture  (GL_TEXTURE_2D, fTexture);

      glPolygonOffset(2,2);
      glEnable(GL_POLYGON_OFFSET_FILL);

      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
      glBindTexture(GL_TEXTURE_2D, fTexture);
      glEnable(GL_TEXTURE_2D);

      DisplayTexture(iSeg.GetPadWidth(),  iSeg.GetPadHeight(),  iSeg.GetRLow(),
		     iSeg.GetNMaxPads(),  iSeg.GetNRows(),
		     0,                   0);
      DisplayTexture(o1Seg.GetPadWidth(), o1Seg.GetPadHeight(), o1Seg.GetRLow(),
		     o1Seg.GetNMaxPads(), o1Seg.GetNRows(),
		     iSeg.GetNMaxPads(),  0);
      DisplayTexture(o2Seg.GetPadWidth(), o2Seg.GetPadHeight(), o2Seg.GetRLow(),
		     o2Seg.GetNMaxPads(), o2Seg.GetNRows(),
		     0,                   o1Seg.GetNRows());

      glDisable(GL_TEXTURE_2D);
    } else {
      DisplayQuads(iSeg.GetPadWidth(),  iSeg.GetPadHeight(),  iSeg.GetRLow(),
		   iSeg.GetNMaxPads(),  iSeg.GetNRows(),
		   0,                   0);
      DisplayQuads(o1Seg.GetPadWidth(), o1Seg.GetPadHeight(), o1Seg.GetRLow(),
		   o1Seg.GetNMaxPads(), o1Seg.GetNRows(),
		   iSeg.GetNMaxPads(),    0);
      DisplayQuads(o2Seg.GetPadWidth(), o2Seg.GetPadHeight(), o2Seg.GetRLow(),
		   o2Seg.GetNMaxPads(), o2Seg.GetNRows(),
		   0,                   o1Seg.GetNRows());
    }
  }

  // Display frame
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  DisplayFrame();

  glPopAttrib();
}

/**************************************************************************/
// Data import
/**************************************************************************/

void TPCSector2DGL::LoadPadrow(TPCSectorData::RowIterator& iter,
			       Int_t row, Int_t col_off) const
{
  Int_t    pad_var;
  Int_t    time, val;   

  Int_t    min_time = fSector->fMinTime;
  Int_t    max_time = fSector->fMaxTime;
  Bool_t   half_border_time = ((fSector->fMaxTime - fSector->fMinTime) % 2 == 0);

  UChar_t* img_pos = GetRowCol(row, col_off);
  while (iter.NextPad()) {
    pad_var = 0; 

    while (iter.Next()) {
      time = iter.Time();
      val  = iter.Signal();

      if(fSector->fShowMax) {
        if(val > pad_var) {
          pad_var = val;
        }
      } else {
        // Integrate int max_val.
        if(time >= min_time && time <= max_time) {
          if(half_border_time && (time == min_time || time == max_time))
            pad_var += val/2;
          else
            pad_var += val;
        }
      }
    }

    pad_var = TMath::Min(pad_var, fSector->fMaxVal);
    if(pad_var > fSector->fthreshold)
      SetCol(pad_var, img_pos);
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
  rowOff[0] = 0; rowOff[1] = rowOff[2] = -TPCSectorData::GetSeg(1).GetFirstRow();
  colOff[0] = colOff[2] = 0; colOff[1] =  TPCSectorData::GetSeg(0).GetNMaxPads();

  // Loop over 3 main segments
  for (Int_t sId = 0; sId <= 2; ++sId) {
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
  // glTexEnvf      (GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,   GL_MODULATE); // Lightning is off anyway.
  glTexImage2D   (GL_TEXTURE_2D, 0, GL_RGBA, fgkTextureWidth, fgkTextureHeight,
                  0, GL_RGBA, GL_UNSIGNED_BYTE, fImage);

}

/**************************************************************************/
// Data display
/**************************************************************************/

void TPCSector2DGL::DisplayTexture(Float_t padW,     Float_t padH, Float_t startR,
                                   Int_t numMaxPads, Int_t numRows, 
                                   Int_t startCol,   Int_t startRow) const
{
  Float_t w  = numMaxPads*padW/2;
  Float_t u1 = (Float_t) startCol / fgkTextureWidth;
  Float_t v1 = (Float_t) startRow / fgkTextureHeight;
  Float_t u2 = u1 + (Float_t) numMaxPads / fgkTextureWidth;
  Float_t v2 = v1 + (Float_t) numRows    / fgkTextureHeight;

  glBegin(GL_QUADS);  
  glTexCoord2f(u1, v1);  glVertex2f(-w, startR);
  glTexCoord2f(u1, v2);  glVertex2f(-w, startR + numRows*padH);
  glTexCoord2f(u2, v2);  glVertex2f( w, startR + numRows*padH);
  glTexCoord2f(u2, v1);  glVertex2f( w, startR);
  glEnd();
}

/**************************************************************************/

void TPCSector2DGL::DisplayQuads(Float_t padW,     Float_t padH, Float_t startR,
                                 Int_t numMaxPads, Int_t numRows, 
                                 Int_t startCol,   Int_t startRow) const
{
  UChar_t *pix;
  Float_t y_d, y_u;
  Float_t x_off, x;

  glBegin(GL_QUADS);
  for (Int_t row=0; row<numRows; row++) {
    y_d = startR + row*padH;
    y_u = y_d + padH;
    x_off = -numMaxPads*padW/2;
    pix = GetRowCol(row + startRow, startCol);
    for (Int_t pad=0; pad<numMaxPads; pad++, pix+=4) {
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
  TColor* c = gROOT->GetColor(fSector->fFrameCol);
  if(c == 0) return; 

  // x[0] = (UChar_t)(255*c->GetRed());  x[1] = (UChar_t)(255*c->GetGreen());
  // x[2] = (UChar_t)(255*c->GetBlue()); x[3] = 255;

  glColor3ub((UChar_t)(255*c->GetRed()), 
             (UChar_t)(255*c->GetGreen()),
             (UChar_t)(255*c->GetBlue()));

  glBegin(GL_LINE_LOOP);
  TraceStepsUp  (TPCSectorData::GetInnSeg());
  TraceStepsDown(TPCSectorData::GetInnSeg());
  glEnd();

  glBegin(GL_LINE_LOOP);
  TraceStepsUp  (TPCSectorData::GetOut1Seg());
  TraceStepsDown(TPCSectorData::GetOut1Seg());
  glEnd();

  glBegin(GL_LINE_STRIP);
  TraceStepsUp  (TPCSectorData::GetOut2Seg());
  TraceStepsDown(TPCSectorData::GetOut2Seg());
  glEnd();
}
