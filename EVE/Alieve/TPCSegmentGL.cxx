// $Header$

#include "TPCSegmentGL.h"

#include <TStyle.h>
#include <TColor.h>
#include <TStopwatch.h>

using namespace Reve;
using namespace Alieve;
using namespace std;

/**************************************************************************/

TPCSegmentGL::TPCSegmentGL() : TGLObject()
{
  fImage   = 0;
  fTexture = 0;
  fRTS     = 0;
}

TPCSegmentGL::~TPCSegmentGL()
{
  if(fImage) delete fImage;
  if(fTexture) glDeleteTextures(1, &fTexture);
}
 

/**************************************************************************/

Bool_t TPCSegmentGL::SetModel(TObject* obj)
{
  if( set_model(obj, "Alieve::TPCSegment")){
    fSegment = (TPCSegment*) fExternalObj;
    return true;
  }
  return false;
}

void TPCSegmentGL::SetBBox()
{
  set_axis_aligned_bbox(((TPCSegment*)fExternalObj)->AssertBBox());
}

/**************************************************************************/

void TPCSegmentGL::DirectDraw(const TGLDrawFlags & ) const
{
  // printf("TPCSegmentGL::DirectDraw \n");

  if(fSegment->fInfo == 0 ) return;
  TPCDigitsInfo* info = fSegment->fInfo;

  if(fRTS < fSegment->fRTS) {
    CreateTexture(info);
    fRTS = fSegment->fRTS;
  }

  glPushAttrib(GL_CURRENT_BIT      | GL_DEPTH_BUFFER_BIT |
	       GL_COLOR_BUFFER_BIT | GL_ENABLE_BIT       | GL_POLYGON_BIT);

  glDisable(GL_LIGHTING);
  glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
  glEnable(GL_COLOR_MATERIAL);
  glDisable(GL_CULL_FACE);


  // rnr digits
  if(fSegment->fUseTexture) {
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

    DisplayTexture(info->fInnSeg.fPadWidth, info->fInnSeg.fPadLength, info->fInnSeg.fRlow,
		    info->fInnSeg.fNMaxPads, info->fInnSeg.fNRows, 0,0);
    DisplayTexture(info->fOut1Seg.fPadWidth, info->fOut1Seg.fPadLength, info->fOut1Seg.fRlow,
    		    info->fOut1Seg.fNMaxPads,info->fOut1Seg.fNRows,info->fInnSeg.fNMaxPads,0);
    DisplayTexture(info->fOut2Seg.fPadWidth, info->fOut2Seg.fPadLength, info->fOut2Seg.fRlow,
    		    info->fOut2Seg.fNMaxPads, info->fOut2Seg.fNRows,0,info->fOut1Seg.fNRows);

    glDisable(GL_TEXTURE_2D);
  } else {
    DisplayQuads(info->fInnSeg.fPadWidth, info->fInnSeg.fPadLength, info->fInnSeg.fRlow,
		  info->fInnSeg.fNMaxPads, info->fInnSeg.fNRows, 0,0);
    DisplayQuads(info->fOut1Seg.fPadWidth, info->fOut1Seg.fPadLength, info->fOut1Seg.fRlow,
		  info->fOut1Seg.fNMaxPads,info->fOut1Seg.fNRows,info->fInnSeg.fNMaxPads,0);
    DisplayQuads(info->fOut2Seg.fPadWidth, info->fOut2Seg.fPadLength, info->fOut2Seg.fRlow,
		  info->fOut2Seg.fNMaxPads, info->fOut2Seg.fNRows,0,info->fOut1Seg.fNRows);
  }

  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  DisplayFrame(info);

  glPopAttrib();
}

/**************************************************************************/

void TPCSegmentGL::LoadPadrow(Int_t row, Int_t col_off) const
{
  AliSimDigits *digit =  &fSegment->fInfo->fSimDigits;

  Int_t    pad_var = 0;
  Int_t    time, pad, val;   
  UChar_t* img_pos;

  Int_t  min_time = fSegment->fMinTime;
  Int_t max_time = fSegment->fMaxTime;
  Bool_t half_border_time = ((fSegment->fMaxTime - fSegment->fMinTime) % 2 == 0);

  Bool_t done_p = false;
  Bool_t save_p = false;

  digit->First();
  do {

    time = digit->CurrentRow();
    pad  = digit->CurrentColumn();
    val  = digit->CurrentDigit();

    Bool_t use_digit = true;
    if(use_digit) {
      if(fSegment->fShowMax) {
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

    if(digit->Next()) {
      if(pad != digit->CurrentColumn())
	save_p = true;
    } else {
      done_p = true;
      save_p = true;
    }

    if(save_p) {
      pad_var = TMath::Min(pad_var, fSegment->fMaxVal);
      if(pad_var > fSegment->fTreshold) {
	img_pos = GetRowCol(row, pad + col_off);
	SetCol(pad_var, img_pos);
      }
      pad_var = 0; 
    }

  } while (!done_p);
}
 
/**************************************************************************/

void TPCSegmentGL::SetCol(Float_t z, UChar_t* pixel) const
{
 
  Int_t n_col = gStyle->GetNumberOfColors();

  Int_t ci = gStyle->GetColorPalette
    (TMath::Min(n_col - 1,
		Int_t((n_col*(z - fSegment->fTreshold))/(fSegment->fMaxVal - fSegment->fTreshold))));

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

void TPCSegmentGL::CreateTexture(TPCDigitsInfo* info) const
{
  AliSimDigits *digit = &info->fSimDigits;
  AliTPCParam  *par   = info->fParameter;

  TTree* t = info->fTree;
  Int_t s, row, ent, off;

  TStopwatch* sw = new  TStopwatch();
  sw->Start();

  // init fImage table
  if (fImage == 0 ){
    fImage = new UChar_t[ImageWidth*ImageHeight*4];
    glGenTextures(1, &fTexture);
  }
  memset(fImage, 0, ImageWidth*ImageHeight*4);

  ent = info->fSegEnt[fSegment->fID];
  if(ent != -1) {
    row=0;
    // info->fInnSeg.Dump();
    while(ent < t->GetEntriesFast()) {
      t->GetEntry(ent);
      par->AdjustSectorRow(digit->GetID(),s,row);
      // printf("AdjustSectorRow DigitID %d sector %d row %d \n",digit->GetID(),s,row );
      if(s != fSegment->fID) break;
      off =  (info->fInnSeg.fNMaxPads - par->GetNPadsLow(row))/2;
      LoadPadrow(row, off);
      ent++;
    }
  }
  ent = info->fSegEnt[fSegment->fID + 36];
  if(ent != -1) {
    row=0;
    // info->fOut1Seg.Dump();
    while(ent < t->GetEntriesFast()) {
      t->GetEntry(ent);
      par->AdjustSectorRow(digit->GetID(),s,row);
      // printf("AdjustSectorRow DigitID %d sector %d row %d \n",digit->GetID(),s,row );
      if(s != (fSegment->fID+36)) break;

      if(row < par->GetNRowUp1()) {
	off =  (info->fOut1Seg.fNMaxPads - par->GetNPadsUp(row))/2;
	LoadPadrow(row, off + info->fInnSeg.fNMaxPads);
      } else {
	off =  (info->fOut2Seg.fNMaxPads - par->GetNPadsUp(row))/2;
	LoadPadrow(row, off); // info->fInnSeg.fNRows - info->fOut1Seg.fNRows);
      }
      ent++;
    }
  }
  sw->Stop();
  // printf("TPCSegment::CreateTexture timer %f\n", sw->RealTime());
  // sw->Dump();

  glBindTexture  (GL_TEXTURE_2D, fTexture);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S,GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T,GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  // glTexEnvf      (GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,   GL_MODULATE); // Lightning is off anyway.
  glTexImage2D   (GL_TEXTURE_2D, 0, GL_RGBA, ImageWidth, ImageHeight,
		  0, GL_RGBA, GL_UNSIGNED_BYTE, fImage);

}

void TPCSegmentGL::DisplayTexture(Float_t pw, Float_t pl, Float_t vR,
				  Int_t fNMaxPads, Int_t fNRows, 
				  Int_t startCol, Int_t startRow) const
{
  Float_t w  = fNMaxPads*pw/2;
  Float_t v1 = 1.0*startRow/ImageHeight;
  Float_t v2 = v1 + 1.0*fNRows/ImageHeight;
  Float_t u1 = 1.0 *startCol/ImageWidth;
  Float_t u2 = u1 + 1.0 *fNMaxPads/ImageWidth;
  // printf("tex coord u1,v1: (%f, %f), v2,u2: (%f,%f) \n", v1, u1, v2, u2);
  // printf("vertex coord >>> nPads %d pw %f, w: %f, y1: %f, y2: %f \n",fNMaxPads,pw, w, vR, vR+fNRows*pl);
  // glColor4f(1.,1.,1.,fSegment->mAlpha); 
  // glColor4f(1.,1.,1.,1); 
  glBegin (GL_QUADS);
  
  glTexCoord2f(u1, v1);  glVertex3f (-w, vR,           0.0);
  glTexCoord2f(u1, v2);  glVertex3f (-w, vR+fNRows*pl, 0.0);
  glTexCoord2f(u2, v2);  glVertex3f ( w, vR+fNRows*pl, 0.0);
  glTexCoord2f(u2, v1);  glVertex3f ( w, vR,           0.0);

  glEnd();
}

/**************************************************************************/

void TPCSegmentGL::DisplayQuads(Float_t pw, Float_t pl, Float_t vR,
				Int_t fNMaxPads, Int_t fNRows, 
				Int_t startCol, Int_t startRow) const
{
  UChar_t *pix;
  Float_t y_d, y_u;
  Float_t x_off, x;

  glBegin(GL_QUADS);
  for(Int_t row = 0; row<fNRows; row++){
    y_d = vR + row*pl;
    y_u = y_d + pl;
    x_off = -fNMaxPads*pw/2;
    for(Int_t pad = 0; pad<fNMaxPads; pad++){
      pix = GetRowCol(row + startRow, pad + startCol);
      x = x_off + pad*pw;
      if(pix[3] != 0){
	glColor4ubv(pix);
	glVertex3f(x+pw, y_d, 0);
	glVertex3f(x,    y_d, 0);
	glVertex3f(x,    y_u, 0);
	glVertex3f(x+pw, y_u, 0);
      }
    }
  }
  glEnd();
}

/**************************************************************************/

void TPCSegmentGL::DisplayFrame(TPCDigitsInfo* info) const
{
  TPCSeg* seg;

  TColor* c = gROOT->GetColor(fSegment->fFrameCol);
  if(c == 0) return; 

  //  x[0] = (UChar_t)(255*c->GetRed());  x[1] = (UChar_t)(255*c->GetGreen());
  // x[2] = (UChar_t)(255*c->GetBlue()); x[3] = 255;

  glColor3ub((UChar_t)(255*c->GetRed()), 
	     (UChar_t)(255*c->GetGreen()),
	     (UChar_t)(255*c->GetBlue()));

  seg = &info->fInnSeg;
  glBegin(GL_LINE_LOOP);
  LoopStepsUp(seg);
  LoopStepsDown(seg);
  glEnd();

  // outer1 segment
  seg = &info->fOut1Seg;
  glBegin(GL_LINE_LOOP);
  LoopStepsUp(seg);
  LoopStepsDown(seg);
  glEnd();

  // outer2 segment
  seg = &info->fOut2Seg;
  glBegin(GL_LINE_STRIP);
  LoopStepsUp(seg);
  LoopStepsDown(seg);
  glEnd();
}
