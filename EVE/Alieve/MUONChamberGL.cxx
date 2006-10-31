#include "MUONChamberGL.h"

#include <Alieve/MUONChamber.h>
#include <Alieve/MUONChamberData.h>

#include <TGLDrawFlags.h>
#include <Reve/QuadSetGL.h>
#include <GL/gl.h>
#include <GL/glu.h>

using namespace Reve;
using namespace Alieve;

//______________________________________________________________________
// MUONChamberGL
//

ClassImp(MUONChamberGL)

//______________________________________________________________________
MUONChamberGL::MUONChamberGL() :
  TGLObject(),
  fChamber(0),
  fRTS(0)
{
  //
  // constructor
  //

}

//______________________________________________________________________
MUONChamberGL::~MUONChamberGL()
{
  //
  // destructor
  //

}

//______________________________________________________________________
Bool_t MUONChamberGL::SetModel(TObject* obj)
{
  //
  // ...
  //

#if ROOT_VERSION_CODE <= ROOT_VERSION(5,11,2)
  if(set_model(obj, "Alieve::MUONChamber")) {
#elif ROOT_VERSION_CODE <= ROOT_VERSION(5,13,0)
  if(SetModelCheckClass(obj, "Alieve::MUONChamber")) {
#else
  if(SetModelCheckClass(obj, Alieve::MUONChamber::Class())) {
#endif
    
    fChamber = (MUONChamber*) fExternalObj;
    
    return kTRUE;

  }

  return kFALSE;

}

//______________________________________________________________________
void MUONChamberGL::SetBBox()
{
  //
  // ...
  //

#if ROOT_VERSION_CODE <= ROOT_VERSION(5,11,2)
  set_axis_aligned_bbox(((MUONChamber*)fExternalObj)->AssertBBox());
#else
  SetAxisAlignedBBox(((MUONChamber*)fExternalObj)->AssertBBox());
#endif

}

//______________________________________________________________________
void MUONChamberGL::DirectDraw(const TGLDrawFlags& /*flags*/) const
{
  //
  // Actual GL drawing.
  //

  glDisable(GL_LIGHTING);

  //Double_t width = 10;
  //glOrtho(-width,+width,-width,+width,-width,+width);

  if(fRTS < fChamber->fRTS) {
    fChamber->UpdateQuads();
    fRTS = fChamber->fRTS;
  }
  
  Bool_t hasData = (fChamber->GetChamberData() != 0);
  
  if(hasData) {

    DrawQuads();
  
  }

  DrawChamberFrame();

}

//______________________________________________________________________
void MUONChamberGL::DrawQuads() const
{
  //
  // draw the digits as GL_QUADS
  //

  glPushAttrib(GL_ENABLE_BIT | GL_POLYGON_BIT);

  glDisable(GL_LIGHTING);
  glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
  glEnable(GL_COLOR_MATERIAL);
  glDisable(GL_CULL_FACE);

  //Float_t c[4]; glGetFloatv(GL_CURRENT_COLOR, c);

  glPolygonMode(GL_FRONT, GL_FILL);
  glPolygonMode(GL_BACK,  GL_LINE);

  glBegin(GL_QUADS);
  for(std::vector<Quad>::iterator q=fChamber->fQuadSet1.Quads().begin(); q!=fChamber->fQuadSet1.Quads().end(); ++q) {
    UChar_t* c = (UChar_t*) &q->color;
    glColor4ubv(c);
    glVertex3fv(q->vertices);
    glVertex3fv(q->vertices + 3);
    glVertex3fv(q->vertices + 6);
    glVertex3fv(q->vertices + 9);
  }
  glEnd();

  glPolygonMode(GL_FRONT, GL_LINE);
  glPolygonMode(GL_BACK,  GL_FILL);

  glBegin(GL_QUADS);
  for(std::vector<Quad>::iterator q=fChamber->fQuadSet2.Quads().begin(); q!=fChamber->fQuadSet2.Quads().end(); ++q) {
    UChar_t* c = (UChar_t*) &q->color;
    glColor4ubv(c);
    glVertex3fv(q->vertices);
    glVertex3fv(q->vertices + 3);
    glVertex3fv(q->vertices + 6);
    glVertex3fv(q->vertices + 9);
  }
  glEnd();

  glPopAttrib();

}

//______________________________________________________________________
void MUONChamberGL::DrawChamberFrame() const
{
  //
  // draw the chamber frame as GL_LINE_LOOP
  // 

  MUONChamberData* chamberData = fChamber->GetChamberData();
  Int_t nDetElem = chamberData->GetNDetElem();
  Float_t *frameCoord;
  Float_t xOrig, yOrig, xRad, yRad, x, y, z;

  UChar_t pix[4];
  pix[0] = 255;
  pix[1] =   0;
  pix[2] =   0;
  pix[3] = 255;

  glColor4ubv(pix);
  
  for (Int_t id = 0; id < nDetElem; id++) {

    frameCoord = chamberData->GetFrameCoord(id);

    if (fChamber->GetID() < 4) {

      xOrig = frameCoord[0];
      yOrig = frameCoord[1];
      xRad  = frameCoord[2];
      yRad  = frameCoord[3];
      z     = frameCoord[4];
      
      xRad += 0.0;
      yRad += 0.0;

      glBegin(GL_LINE_LOOP);

      glVertex3f(xOrig,yOrig,z);
      
      Int_t nstep = 100;
      Float_t dstep = TMath::Pi()/2.0 / (Float_t)nstep;
      Float_t d;
      for (Int_t istep = 0; istep < nstep; istep++) {

	d = istep * dstep;
	x = xOrig + xRad * TMath::Cos(d);
	y = yOrig + yRad * TMath::Sin(d);

	glVertex3f(x,y,z);

      }
      
      glVertex3f(xOrig,yOrig,z);

      glEnd();

    } else {

      glBegin(GL_LINE_LOOP);
      glVertex3f(frameCoord[0],frameCoord[1],frameCoord[4]);
      glVertex3f(frameCoord[0],frameCoord[3],frameCoord[4]);
      glVertex3f(frameCoord[2],frameCoord[3],frameCoord[4]);
      glVertex3f(frameCoord[2],frameCoord[1],frameCoord[4]);
      glVertex3f(frameCoord[0],frameCoord[1],frameCoord[4]);
      glEnd();

    }

  }

}
