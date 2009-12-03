// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel & Bogdan Vulpescu: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveMUONChamberGL.h"

#include <EveDet/AliEveMUONChamber.h>
#include <EveDet/AliEveMUONChamberData.h>

#include <TGLRnrCtx.h>
#include <TGLIncludes.h>
#include <TMath.h>


//______________________________________________________________________________
// AliEveMUONChamberGL
//

ClassImp(AliEveMUONChamberGL)

//______________________________________________________________________________
AliEveMUONChamberGL::AliEveMUONChamberGL() :
  TGLObject(),
  fChamber(0),
  fQS1(), fQS2(),
  fRTS(0)
{
  //
  // constructor
  //

}

//______________________________________________________________________________
AliEveMUONChamberGL::~AliEveMUONChamberGL()
{
  //
  // destructor
  //

}

//______________________________________________________________________________
Bool_t AliEveMUONChamberGL::SetModel(TObject* obj, const Option_t* /*opt*/)
{
  //
  // ...
  //

  if(SetModelCheckClass(obj, AliEveMUONChamber::Class())) {

    fChamber = (AliEveMUONChamber*) fExternalObj;
    fQS1.SetModel(&fChamber->fQuadSet1);
    fQS2.SetModel(&fChamber->fQuadSet2);
    return kTRUE;

  }

  return kFALSE;

}

//______________________________________________________________________________
void AliEveMUONChamberGL::SetBBox()
{
  //
  // ...
  //

  SetAxisAlignedBBox(((AliEveMUONChamber*)fExternalObj)->AssertBBox());

}

//______________________________________________________________________________
void AliEveMUONChamberGL::DirectDraw(TGLRnrCtx& rnrCtx) const
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

    DrawQuads(rnrCtx);
    DrawPoints();

  }

  DrawChamberFrame();

}

//______________________________________________________________________________
void AliEveMUONChamberGL::DrawQuads(TGLRnrCtx& rnrCtx) const
{
  //
  // draw the digits as GL_QUADS
  //

  glPushAttrib(GL_ENABLE_BIT | GL_POLYGON_BIT);

  glDisable(GL_LIGHTING);
  glDisable(GL_CULL_FACE);

  //Float_t c[4]; glGetFloatv(GL_CURRENT_COLOR, c);

  glPolygonMode(GL_FRONT, GL_FILL);
  glPolygonMode(GL_BACK,  GL_LINE);

  fQS1.DirectDraw(rnrCtx);

  glPolygonMode(GL_FRONT, GL_LINE);
  glPolygonMode(GL_BACK,  GL_FILL);

  fQS2.DirectDraw(rnrCtx);

  glPopAttrib();

}

//______________________________________________________________________________
void AliEveMUONChamberGL::DrawPoints() const
{
  //
  // draw the clusters as GL_QUADS
  //

  Float_t x, y, z;

  glDisable(GL_LIGHTING);
  glLineWidth(1.0);

  TGLUtil::Color3f(1.0,1.0,1.0);

  glBegin(GL_LINES);

  // clusters

  Int_t clsSize = fChamber->fClusterSize;

  if (clsSize > 1) {

    for (Int_t i = 0; i < fChamber->fPointSet1.GetN(); i++) {

      fChamber->fPointSet1.GetPoint(i,x,y,z);

      glVertex3f(x-clsSize,y+clsSize,z);
      glVertex3f(x+clsSize,y-clsSize,z);

      glVertex3f(x-clsSize,y-clsSize,z);
      glVertex3f(x+clsSize,y+clsSize,z);

    }

  }

  // hits

  Int_t hitSize = fChamber->fHitSize;

  if (hitSize > 1) {

    for (Int_t i = 0; i < fChamber->fPointSet2.GetN(); i++) {

      fChamber->fPointSet2.GetPoint(i,x,y,z);

      glVertex3f(x-hitSize,y,z);
      glVertex3f(x+hitSize,y,z);

      glVertex3f(x,y-hitSize,z);
      glVertex3f(x,y+hitSize,z);

    }

  }

  glEnd();

}

//______________________________________________________________________________
void AliEveMUONChamberGL::DrawChamberFrame() const
{
  //
  // draw the chamber frame as GL_LINE_LOOP
  //

  AliEveMUONChamberData* chamberData = fChamber->GetChamberData();
  Int_t nDetElem = chamberData->GetNDetElem();
  Float_t *frameCoord;
  Float_t xOrig, yOrig, xRad, yRad, x, y, z;

  TGLUtil::Color4ub(255, 0, 0, 255);

  for (Int_t id = 0; id < nDetElem; id++) {

    frameCoord = chamberData->GetFrameCoord(id);

    if (fChamber->GetID() < 4) {

      xOrig = frameCoord[0];
      yOrig = frameCoord[1];
      xRad  = frameCoord[2];
      yRad  = frameCoord[3];
      z     = frameCoord[4];

      if (fChamber->GetID() < 2) {
	xRad += TMath::Sign(15.0,(Double_t)xRad);
	yRad += TMath::Sign(15.0,(Double_t)yRad);
      } else {
	xRad += TMath::Sign( 5.0,(Double_t)xRad);
	yRad += TMath::Sign( 5.0,(Double_t)yRad);
      }

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
