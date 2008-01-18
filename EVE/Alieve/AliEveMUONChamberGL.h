// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
#ifndef ALIEVE_MUONChamberGL_H
#define ALIEVE_MUONChamberGL_H

#include <TGLObject.h>
#include <TEveQuadSetGL.h>

class TEveQuadSetGL;


class AliEveMUONChamber;

class AliEveMUONChamberGL : public TGLObject
{

  AliEveMUONChamberGL(const AliEveMUONChamberGL&);            // Not implemented
  AliEveMUONChamberGL& operator=(const AliEveMUONChamberGL&); // Not implemented

 protected:

  virtual void DirectDraw(TGLRnrCtx & rnrCtx) const;
  void DrawChamberFrame() const;
  void DrawQuads(TGLRnrCtx& rnrCtx) const;
  void DrawPoints() const;

  AliEveMUONChamber*             fChamber; // fModel dynamic-casted to AliEveMUONChamberGL
  TEveQuadSetGL            fQS1;
  TEveQuadSetGL            fQS2;

  mutable UInt_t           fRTS;     // render time stamp

 public:

  AliEveMUONChamberGL();
  virtual ~AliEveMUONChamberGL();

  virtual Bool_t SetModel(TObject* obj, const Option_t* opt=0);
  virtual void   SetBBox();

  ClassDef(AliEveMUONChamberGL,1);   // the GL drawing class of one chamber

};

#endif
