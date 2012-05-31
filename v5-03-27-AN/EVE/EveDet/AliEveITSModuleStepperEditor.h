// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveITSModuleStepperEditor_H
#define AliEveITSModuleStepperEditor_H

#include <TGedFrame.h>

class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;

class TEveGridStepperSubEditor;


class AliEveITSModuleStepper;

class AliEveITSModuleStepperEditor : public TGedFrame
{
public:
  AliEveITSModuleStepperEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30, UInt_t options = kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~AliEveITSModuleStepperEditor() {}

  virtual void SetModel(TObject* obj);

  void         UpdateStepper();

protected:
  AliEveITSModuleStepper   *fM;         // Model object.

  TEveGridStepperSubEditor *fStepper;   // GUI component for grid-stepper control.

private:
  AliEveITSModuleStepperEditor(const AliEveITSModuleStepperEditor&);            // Not implemented
  AliEveITSModuleStepperEditor& operator=(const AliEveITSModuleStepperEditor&); // Not implemented

  ClassDef(AliEveITSModuleStepperEditor, 0); // Editor for AliEveITSModuleStepper.
};

#endif
