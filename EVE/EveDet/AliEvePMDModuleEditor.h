// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEvePMDModuleEditor_H
#define AliEvePMDModuleEditor_H

#include <TGedFrame.h>

class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;

class AliEvePMDModule;

class AliEvePMDModuleEditor : public TGedFrame
{
public:
  AliEvePMDModuleEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30, UInt_t options = kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~AliEvePMDModuleEditor() {}

  virtual void SetModel(TObject* obj);
  void DisplayHistos();

protected:
  AliEvePMDModule* fM; // Model object.

  TGVerticalFrame*  fInfoFrame; // Top frame for info labels.

  TGLabel*   fInfoLabel0; // label
  TGLabel*   fInfoLabel1; // label
  TGLabel*   fInfoLabel2; // label
  TGLabel*   fInfoLabel3; // label
  TGLabel*   fInfoLabel4; // label
  TGLabel*   fInfoLabel5; // label

private:
  void CreateInfoFrame();

  AliEvePMDModuleEditor(const AliEvePMDModuleEditor&);            // Not implemented
  AliEvePMDModuleEditor& operator=(const AliEvePMDModuleEditor&); // Not implemented

  ClassDef(AliEvePMDModuleEditor, 0); // Editor for AliEvePMDModule
};

#endif
