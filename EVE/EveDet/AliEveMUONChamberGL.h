// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel & Bogdan Vulpescu: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
#ifndef AliEveMUONChamberGL_H
#define AliEveMUONChamberGL_H

// #include <TGLObject.h>
#include <TEveQuadSetGL.h>

class TEveQuadSetGL;

class AliEveMUONChamber;

class AliEveMUONChamberGL : public TGLObject
{
public:
  AliEveMUONChamberGL();
  virtual ~AliEveMUONChamberGL();

  virtual Bool_t SetModel(TObject* obj, const Option_t* opt=0);
  virtual void   SetBBox();

protected:
  virtual void DirectDraw(TGLRnrCtx & rnrCtx) const;
  void DrawChamberFrame() const;
  void DrawQuads(TGLRnrCtx& rnrCtx) const;
  void DrawPoints() const;

  AliEveMUONChamber       *fChamber; // Model object.
  TEveQuadSetGL            fQS1;
  TEveQuadSetGL            fQS2;

  mutable UInt_t           fRTS;     // render time stamp

private:
  AliEveMUONChamberGL(const AliEveMUONChamberGL&);            // Not implemented
  AliEveMUONChamberGL& operator=(const AliEveMUONChamberGL&); // Not implemented

  ClassDef(AliEveMUONChamberGL, 0);   // the GL drawing class of one chamber
};

#endif
