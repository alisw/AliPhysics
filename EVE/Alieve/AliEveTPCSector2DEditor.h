// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 * 
 **************************************************************************/

#ifndef ALIEVE_TPCSector2DEditor_H
#define ALIEVE_TPCSector2DEditor_H

#include <TGedFrame.h>

class TGCheckButton;
class TGComboBox;


class AliEveTPCSector2D;

class AliEveTPCSector2DEditor : public TGedFrame
{
  AliEveTPCSector2DEditor(const AliEveTPCSector2DEditor&);            // Not implemented
  AliEveTPCSector2DEditor& operator=(const AliEveTPCSector2DEditor&); // Not implemented

protected:
  AliEveTPCSector2D* fM; // fModel dynamic-casted to AliEveTPCSector2DEditor

  TGCheckButton*   fShowMax;
  TGCheckButton*   fAverage;

  TGCheckButton*   fUseTexture;
  TGCheckButton*   fPickEmpty;
  TGComboBox*      fPickMode;

public:
  AliEveTPCSector2DEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
		    UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  ~AliEveTPCSector2DEditor();

  virtual void SetModel(TObject* obj);

  void DoShowMax();
  void DoAverage();
  void SetupAverage();

  void DoUseTexture();
  void DoPickEmpty();
  void DoPickMode(Int_t mode);

  ClassDef(AliEveTPCSector2DEditor, 0); // Editor for AliEveTPCSector2D
}; // endclass AliEveTPCSector2DEditor

#endif
