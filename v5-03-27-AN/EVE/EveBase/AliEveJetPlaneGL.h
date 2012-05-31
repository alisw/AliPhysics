// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveJetPlaneGL_H
#define AliEveJetPlaneGL_H

#include <TGLObject.h>

class TGLViewer;
class TGLScene;

class AliEveJetPlane;

//==============================================================================
//
// AliEveJetPlaneGL
//
// GL rendering code for AliEveJetPlane class.

class AliEveJetPlaneGL : public TGLObject
{
public:
  AliEveJetPlaneGL();
  virtual ~AliEveJetPlaneGL() {}

  virtual Bool_t SetModel(TObject* obj, const Option_t* opt=0);
  virtual void   SetBBox();

  // To support two-level selection
  // virtual Bool_t SupportsSecondarySelect() const { return kTRUE; }
  // virtual void ProcessSelection(TGLRnrCtx & rnrCtx, TGLSelectRecord & rec);

protected:
  AliEveJetPlane* fM; // Model object.

  virtual void DirectDraw(TGLRnrCtx & rnrCtx) const;

private:
  AliEveJetPlaneGL(const AliEveJetPlaneGL&);            // Not implemented
  AliEveJetPlaneGL& operator=(const AliEveJetPlaneGL&); // Not implemented

  ClassDef(AliEveJetPlaneGL, 0); // GL renderer for AliEveJetPlane.
}; // endclass AliEveJetPlaneGL

#endif
