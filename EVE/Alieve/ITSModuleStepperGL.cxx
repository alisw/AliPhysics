// $Header$

#include "ITSModuleStepperGL.h"

#include <Reve/GLTextNS.h>
#include <Reve/GLUtilNS.h>
#include <Reve/RGTopFrame.h>
#include <Reve/RGEditor.h>
#include <Reve/RGBAPalette.h>

#include <Alieve/ITSModuleStepper.h>
#include <Alieve/ITSScaledModule.h>

#include <TGLRnrCtx.h>
#include <TGLSelectRecord.h>

using namespace Reve;
using namespace Alieve;

//______________________________________________________________________
// ITSModuleStepperGL
//

ClassImp(ITSModuleStepperGL)

ITSModuleStepperGL::ITSModuleStepperGL() : TGLObject(), fM(0)
{
  fDLCache = false; // Disable display list.
}

ITSModuleStepperGL::~ITSModuleStepperGL()
{}

/**************************************************************************/

Bool_t ITSModuleStepperGL::SetModel(TObject* obj, const Option_t* /*opt*/)
{
  if(SetModelCheckClass(obj, ITSModuleStepper::Class())) {
    fM = dynamic_cast<ITSModuleStepper*>(obj);
    return kTRUE;
  }
  return kFALSE;
}

void ITSModuleStepperGL::SetBBox()
{
  // !! This ok if master sub-classed from TAttBBox
  SetAxisAlignedBBox(fM->AssertBBox());
}

/**************************************************************************/

void ITSModuleStepperGL::DirectDraw(TGLRnrCtx & rnrCtx) const
{
  // printf("ITSModuleStepperGL::DirectDraw Style %d, LOD %d\n", rnrCtx.Style(), rnrCtx.LOD()); 

  ITSModuleStepper& MS = *fM;
  Int_t W = Int_t(MS.fStepper->Dx*MS.fStepper->Nx);
  Int_t H = Int_t(MS.fStepper->Dy*MS.fStepper->Ny);
  Float_t dx = W*MS.fWWidth;
  Float_t dy = 6; // H*MS.fWHeight;

  GLboolean lightp;
  glGetBooleanv(GL_LIGHTING, &lightp);
  if (lightp) glDisable(GL_LIGHTING);
   
  UChar_t color[4];
  ColorFromIdx(MS.fWColor, color);
  glColor4ubv(color);

  // render frame of grid stepper
  if (MS.fRnrFrame)
  {
    glBegin(GL_LINE_LOOP);
    glVertex2f(-1,  -1);       
    glVertex2f(W+1, -1);
    glVertex2f(W+1,  H+1);
    glVertex2f(-1 ,  H+1);
    glEnd();
  }

  // triangles
  glPushAttrib(GL_ENABLE_BIT | GL_POLYGON_BIT);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glDisable(GL_CULL_FACE);
  glEnable(GL_BLEND);
  glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
    
  Float_t sx =0 ,sy = 0;
  switch(MS.fWCorner) {
    case ITSModuleStepper::PT_BottomLeft:
      sy = -dy;
      break;
    case ITSModuleStepper::PT_BottomRight:
      sy = -dy;
      break;
    case ITSModuleStepper::PT_TopLeft:
      sy = H;
      break;
    case ITSModuleStepper::PT_TopRight:
      sy = H;
      break;
    default:
      sy = dy;
      break;
  }
  
  if (rnrCtx.SecSelection()) glPushName(0);
  glPushMatrix();
  glTranslatef(sx, sy, 0.);

  // pager
  if (rnrCtx.SecSelection()) glLoadName(2);
  RenderSymbol(dx, dy, 2);
  glTranslatef(dx, 0, 0);
  if (rnrCtx.SecSelection()) glLoadName(1);
  RenderSymbol(dx*1.2, dy, 1);
  glTranslatef(dx, 0, 0);
  RenderString(Form(" %d/%d ", MS.GetCurrentPage(), MS.GetPages()), dy);
  if (rnrCtx.SecSelection()) glLoadName(3);
  RenderSymbol(dx*1.2, dy, 3);
  glTranslatef(dx, 0, 0);
  if (rnrCtx.SecSelection()) glLoadName(4);
  RenderSymbol(dx, dy, 4);
  glTranslatef(2*dx, 0, 0);
  
  // scale info
  Int_t cnx = 0, cnz = 0;
  ITSDigitsInfo* di = MS.fDigitsInfo;
  Int_t scale = fM->fScaleInfo->GetScale() - 1;
  ITSScaledModule* sm = dynamic_cast<ITSScaledModule*>(*fM->BeginChildren());
  switch(sm->GetSubDetID())
  {
    case 0: 
      cnx = di->fSPDScaleX[scale], cnz = di->fSPDScaleZ[scale];
      break;
    case 1: 
      cnx = di->fSDDScaleX[scale], cnz = di->fSDDScaleZ[scale];
      break;
    case 2:
      cnx = di->fSSDScale[scale], cnz = 1;
      break;
  }
  if (rnrCtx.SecSelection()) glLoadName(0);
  RenderString(Form("Scale"), dy);
  glTranslatef(0.07*dx, 0, 0);
  // up down arrows 
  if (rnrCtx.SecSelection()) glLoadName(6);
  RenderSymbol(dx*1.2, dy*0.9, 5);

  if (rnrCtx.SecSelection()) glLoadName(7);
  RenderSymbol(dx*1.2, dy*0.9, 6);

  glTranslatef(1*dx, 0, 0);
  if (rnrCtx.SecSelection()) glLoadName(0);
  RenderString(Form(" %dx%d ", cnx, cnz), dy, kFALSE);

  glPopMatrix();
  if (rnrCtx.SecSelection()) glLoadName(5);
  glPushMatrix();
  glTranslatef(W+2, 0, 0);
  RenderPalette(H, 4);
  glPopMatrix();

  if (rnrCtx.SecSelection()) glPopName();

  glPopAttrib();

  if (lightp) glEnable(GL_LIGHTING);
}

/**************************************************************************/

void ITSModuleStepperGL::RenderPalette(Float_t dx, Float_t dy) const
{
  ITSModule* qs = dynamic_cast<ITSModule*>(*fM->BeginChildren());
  RGBAPalette* p = qs->GetPalette();
  Float_t xs = dx/(p->GetMaxVal()- p->GetMinVal());
  Float_t ys = dy;

  Float_t x  = 0;
  glBegin(GL_QUAD_STRIP);
  for(Int_t i=p->GetMinVal(); i<=p->GetMaxVal(); i++) 
  {
    glColor4ubv(p->ColorFromValue(i + p->GetMinVal()));
    glVertex2f(0,  x);
    glVertex2f(ys, x);
    x+=xs;
  }
  glEnd();
}

/**************************************************************************/

void ITSModuleStepperGL::RenderSymbol(Float_t dx, Float_t dy, Int_t id) const
{
  Float_t xs = dx/4, ys = dy/4;

  if(id == 0) {
    glBegin(GL_QUADS);
    glVertex2f(0,ys); glVertex2f(0, ys*3); 
    glVertex2f(dx, ys*3); glVertex2f(dx, ys);
    glEnd();
    return;
  }
  

  glBegin(GL_TRIANGLES);
  switch (id) {
    case 1:
    {
      glVertex2f(xs*2.5, ys*3); glVertex2f(xs*1.5, ys*2); glVertex2f(xs*2.5, ys);
      break;
    }
    case 2:
    {
      glVertex2f(xs*2, ys*3); glVertex2f(xs, ys*2); glVertex2f(xs*2, ys);
      glVertex2f(xs*3, ys*3); glVertex2f(xs*2, ys*2); glVertex2f(xs*3, ys);
      break;
    }
    case 3:
    {
      glVertex2f(xs*1.5, ys);  glVertex2f(xs*2.5, ys*2); glVertex2f(xs*1.5, ys*3);
      break;
    }
    case 4:
    {
      glVertex2f(xs, ys);  glVertex2f(xs*2, ys*2); glVertex2f(xs, ys*3);
      glVertex2f(xs*2, ys); glVertex2f(xs*3, ys*2); glVertex2f(xs*2, ys*3);
      break;
    }
    case 5:
    {
      glVertex2f(xs, ys*2.5);  glVertex2f(xs*2, ys*3.5); glVertex2f(xs*3, ys*2.5);
      break;
    }
    case 6:
    {
      glVertex2f(xs, ys*1.5);  glVertex2f(xs*2, ys*0.5); glVertex2f(xs*3, ys*1.5);
      break;
    }
   
    default:
      break;
  }
  glEnd();
}

/**************************************************************************/
void ITSModuleStepperGL::RenderString(TString info, Float_t dy, Bool_t trans) const
{
  Float_t movex = 0; 
  
  GLUtilNS::GL_Capability_Switch texure_on(GL_TEXTURE_2D, true);
  GLTextNS::txfBindFontTexture(GLTextNS::fgDefaultFont);

  glPushMatrix();
  glTranslatef(0, dy*0.25, 0);
  Float_t s = (dy)/ (GLTextNS::fgDefaultFont->max_height());
  Float_t sx = s*0.75; Float_t sy = s*0.8;
  glScalef(sx, sy, 1);
  txfRenderString(GLTextNS::fgDefaultFont, info.Data(), info.Length());
  Int_t w, ma, md;
  txfGetStringMetrics(GLTextNS::fgDefaultFont,info.Data(), info.Length() , w, ma, md);
  movex = w*sx;
  glPopMatrix();

  if(trans)
    glTranslatef(movex, 0, 0);
}

/**************************************************************************/

void ITSModuleStepperGL::ProcessSelection(TGLRnrCtx       & /*rnrCtx*/,
					  TGLSelectRecord & rec)
{
  // Processes secondary selection from TGLViewer.
  // Calls TPointSet3D::PointSelected(Int_t) with index of selected
  // point as an argument.

  if (rec.GetN() < 2) return;

  switch (rec.GetItem(1)) {
    case 1:
      fM->Previous();
      break;
    case 2:
      fM->Start();
      break;
    case 3:
      fM->Next();
      break;
    case 4:
      fM->End();
      break;
    case 5:
      gReve->GetEditor()->DisplayRenderElement(*fM->BeginChildren());
      break;
    case 6:
    {
      DigitScaleInfo* si = fM->fScaleInfo;
      if(si->fScale < 5) 
      {
	si->ScaleChanged(si->fScale + 1);
	gReve->GetEditor()->DisplayObject(gReve->GetEditor()->GetModel());
	gReve->Redraw3D();
      }
      break;
    }
    case 7:
    {
      DigitScaleInfo* si = fM->fScaleInfo;
      if(si->fScale > 1) 
      {
	si->ScaleChanged(si->GetScale() - 1);
	gReve->GetEditor()->DisplayObject(gReve->GetEditor()->GetModel());
	gReve->Redraw3D();
      }
      break;
    }
    default:
      break;
  }
}
