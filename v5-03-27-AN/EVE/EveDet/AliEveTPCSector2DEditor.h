// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveTPCSector2DEditor_H
#define AliEveTPCSector2DEditor_H

#include <TGedFrame.h>

class TGCheckButton;
class TGComboBox;

class AliEveTPCSector2D;

//------------------------------------------------------------------------------
// AliEveTPCSector2DEditor
//
// GUI editor for AliEveTPCSector2D.
//

class AliEveTPCSector2DEditor : public TGedFrame
{
public:
  AliEveTPCSector2DEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
		    UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~AliEveTPCSector2DEditor() {}

  virtual void SetModel(TObject* obj);

  void DoShowMax();
  void DoAverage();
  void SetupAverage();

  void DoUseTexture();
  void DoPickEmpty();
  void DoPickMode(Int_t mode);

protected:
  AliEveTPCSector2D *fM;            // Model object.

  TGCheckButton     *fShowMax;      // Check to show maximum signal.
  TGCheckButton     *fAverage;      // Check-box to show average of the signal.

  TGCheckButton     *fUseTexture;   // Check-box to use texture.
  TGCheckButton     *fPickEmpty;    // Check-box for picking of empty pads.
  TGComboBox        *fPickMode;     // Selector of pick-mode.

private:
  AliEveTPCSector2DEditor(const AliEveTPCSector2DEditor&);            // Not implemented
  AliEveTPCSector2DEditor& operator=(const AliEveTPCSector2DEditor&); // Not implemented

  ClassDef(AliEveTPCSector2DEditor, 0); // Editor for AliEveTPCSector2D.
};

#endif
