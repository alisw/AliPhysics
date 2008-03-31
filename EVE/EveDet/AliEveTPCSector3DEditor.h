// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveTPCSector3DEditor_H
#define AliEveTPCSector3DEditor_H

#include <TGedFrame.h>

class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;

class TEveGValuator;
class TEveGDoubleValuator;

class AliEveTPCSector3D;

//------------------------------------------------------------------------------
// AliEveTPCSector3DEditor
//
// Editor for AliEveTPCSector3D.
//

class AliEveTPCSector3DEditor : public TGedFrame
{

public:
  AliEveTPCSector3DEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
		    UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~AliEveTPCSector3DEditor() {}

  virtual void SetModel(TObject* obj);

  void DoRnrFrame();
  void DoDriftVel();

  void DoPointFrac();
  void DoPointSize();

protected:
  AliEveTPCSector3D   *fM;         // Model object.

  TGCheckButton       *fRnrFrame;  // Check-box for frame rendering.
  TEveGValuator       *fDriftVel;  // Valuator for drift velocity.

  TEveGValuator       *fPointFrac; // Valuator for signal fraction displayed as points.
  TEveGValuator       *fPointSize; // Size of point in GL.

private:
  AliEveTPCSector3DEditor(const AliEveTPCSector3DEditor&);            // Not implemented
  AliEveTPCSector3DEditor& operator=(const AliEveTPCSector3DEditor&); // Not implemented

  ClassDef(AliEveTPCSector3DEditor, 0); // Editor for AliEveTPCSector3D.
};

#endif
