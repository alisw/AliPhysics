// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveITSScaledModuleEditor_H
#define AliEveITSScaledModuleEditor_H

#include <TGedFrame.h>

class TGNumberEntry;
class TGColorSelect;
class TGComboBox;

class TEveGValuator;
class TEveGDoubleValuator;
class TEveRGBAPalette;

class AliEveDigitScaleInfo;
class AliEveITSScaledModule;
class AliITSsegmentation;

/******************************************************************************/

class AliEveITSScaledModuleEditor : public TGedFrame
{
public:
  AliEveITSScaledModuleEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30, UInt_t options = kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~AliEveITSScaledModuleEditor() {}

  virtual void   SetModel(TObject* obj);

  void DoScale();
  void DoStatType(Int_t t);

protected:
  AliEveITSScaledModule *fModule;       // Model object.

  TGNumberEntry         *fScale;        // Number-entry for digit-scale.
  TGComboBox            *fStatistic;    // Selection of scaling algorithm.

  TGVerticalFrame       *fInfoFrame;    // Frame in "Info" tab.
  TGLabel               *fInfoLabel0;   // Info text.
  TGLabel               *fInfoLabel1;   // Info text.

private:
  void CreateInfoFrame();

  AliEveITSScaledModuleEditor(const AliEveITSScaledModuleEditor&);            // Not implemented
  AliEveITSScaledModuleEditor& operator=(const AliEveITSScaledModuleEditor&); // Not implemented

  ClassDef(AliEveITSScaledModuleEditor, 0); // Editor for AliEveITSScaledModule.
};

#endif
