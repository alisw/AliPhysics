// $Header$

#include "ITSModuleStepperGL.h"

#include <Reve/GLTextNS.h>
#include <Reve/GLUtilNS.h>
#include <Alieve/ITSModuleStepper.h>

#include <TGLDrawFlags.h>

using namespace Reve;
using namespace Alieve;

//______________________________________________________________________
// ITSModuleStepperGL
//

ClassImp(ITSModuleStepperGL)

ITSModuleStepperGL::ITSModuleStepperGL() : TGLObject(), fM(0)
{
  // fCached = false; // Disable display list.
}

ITSModuleStepperGL::~ITSModuleStepperGL()
{}

/**************************************************************************/

Bool_t ITSModuleStepperGL::SetModel(TObject* obj)
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

void ITSModuleStepperGL::DirectDraw(const TGLDrawFlags & flags) const
{
  // printf("ITSModuleStepperGL::DirectDraw Style %d, LOD %d\n", flags.Style(), flags.LOD()); 

  ITSModuleStepper& MS = *fM;
  Int_t W = Int_t(MS.fStepper->Dx*MS.fStepper->Nx);
  Int_t H = Int_t(MS.fStepper->Dy*MS.fStepper->Ny);
  Float_t dx = W*MS.fWWidth;
  Float_t dy = H*MS.fWHeight;

  GLboolean lightp;
  glGetBooleanv(GL_LIGHTING, &lightp);
  if (lightp) glDisable(GL_LIGHTING);
   
  // render frame of grid stepper
  if (MS.fRnrFrame)
  {
    UChar_t color[4];
    ColorFromIdx(MS.fFrameColor, color);
    glColor4ubv(color);

    glBegin(GL_LINE_LOOP);
    glVertex2f(0., 0.);       
    glVertex2f(W , 0.);
    glVertex2f(W , H );
    glVertex2f(0 , H );
    glEnd();
  }

  // triangles
  glPushAttrib(GL_ENABLE_BIT | GL_POLYGON_BIT);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  // glEnable(GL_POLYGON_OFFSET_FILL);

  UChar_t color[4];
  ColorFromIdx(MS.fWColor, color);
  glColor4ubv(color);
    
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
  
  if (flags.SecSelection()) glPushName(0);

  glPushMatrix();
  // traslate according to orentation 
  glTranslatef(sx, sy, 0.);

  if (flags.SecSelection()) glLoadName(2);
  RenderTriangle(dx, dy, 2);

  glTranslatef(dx, 0, 0);
  if (flags.SecSelection()) glLoadName(1);
  RenderTriangle(dx, dy, 1);
  glTranslatef(dx, 0, 0);

  Float_t movex = 0; 
  {
    GLUtilNS::GL_Capability_Switch texure_on(GL_TEXTURE_2D, true);
    GLTextNS::txfBindFontTexture(GLTextNS::fgDefaultFont);
    glPushMatrix();
    glTranslatef(0, dy*0.25, 0);
    Float_t s = (dy)/ (GLTextNS::fgDefaultFont->max_height());
    glScalef(s, s, 1);
    TString info = Form(" %d / %d ", MS.GetCurrentPage(),MS.GetPages());
    txfRenderString(GLTextNS::fgDefaultFont, info.Data(), info.Length());
    Int_t w, ma, md;
    txfGetStringMetrics(GLTextNS::fgDefaultFont,info.Data(), info.Length() , w, ma, md);
    movex = w*s;
    glPopMatrix();
  }
  glTranslatef(movex, 0, 0); 

  if (flags.SecSelection()) glLoadName(3);
  RenderTriangle(dx, dy, 3);

  glTranslatef(dx, 0, 0);
  if (flags.SecSelection()) glLoadName(4);
  RenderTriangle(dx, dy, 4);

  glPopMatrix();

  if (flags.SecSelection()) glPopName();

  glPopAttrib();

  if (lightp) glEnable(GL_LIGHTING);
      
}
/**************************************************************************/
void ITSModuleStepperGL::RenderTriangle(Float_t dx, Float_t dy, Int_t id) const
{
  Float_t xs = dx/4, ys = dy/4;

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
    default:
      break;
  }
  glEnd();

  /*
    glBegin(GL_LINE_LOOP);
    glVertex2f(0,0);
    glVertex2f(0,dy);
    glVertex2f(dx,dy);
    glVertex2f(dx,0);
    glEnd();
  */
}

/**************************************************************************/

void ITSModuleStepperGL::ProcessSelection(UInt_t* ptr, TGLViewer*, TGLScene*)
{
  // Processes secondary selection from TGLViewer.
  // Calls TPointSet3D::PointSelected(Int_t) with index of selected
  // point as an argument.

  if (ptr[0] < 2) return;
  UInt_t id = ptr[4]; 

  if(id == 1)
    fM->Previous();
  if(id == 2)
    fM->Start();
  if(id == 3)
    fM->Next();
  if(id == 4)
    fM->End();
}
