/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveVZEROModuleEditor_H
#define AliEveVZEROModuleEditor_H

#include <TGedFrame.h>

class TEveGValuator;

class AliEveVZEROModule;

//------------------------------------------------------------------------------
// AliEveVZEROModuleEditor
//
// Editor for AliEveVZEROModule.

class AliEveVZEROModuleEditor : public TGedFrame
{
public:
  AliEveVZEROModuleEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
		     UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~AliEveVZEROModuleEditor() {}

  virtual void SetModel(TObject* obj);

  void DoSampleIndex();

protected:
  AliEveVZEROModule    *fM;            // Model dynamic-casted to AliEveVZEROModuleEditor

  TEveGValuator        *fSampleIndex;  // Widget for ADC sample index.

private:
  AliEveVZEROModuleEditor(const AliEveVZEROModuleEditor&);            // Not implemented
  AliEveVZEROModuleEditor& operator=(const AliEveVZEROModuleEditor&); // Not implemented

  ClassDef(AliEveVZEROModuleEditor, 0) // Editor for AliEveVZEROModule
};

#endif
