// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 * 
 **************************************************************************/

#ifndef ALIEVE_PMDModuleEditor_H
#define ALIEVE_PMDModuleEditor_H

#include <TGedFrame.h>

class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;


class AliEvePMDModule;

class AliEvePMDModuleEditor : public TGedFrame
{
private:
  AliEvePMDModuleEditor(const AliEvePMDModuleEditor&);            // Not implemented
  AliEvePMDModuleEditor& operator=(const AliEvePMDModuleEditor&); // Not implemented

  void CreateInfoFrame();

protected:
  AliEvePMDModule* fM; // fModel dynamic-casted to AliEvePMDModuleEditor

  TGVerticalFrame*  fInfoFrame;

  TGLabel*   fInfoLabel0;
  TGLabel*   fInfoLabel1;
  TGLabel*   fInfoLabel2;
  TGLabel*   fInfoLabel3;
  TGLabel*   fInfoLabel4;
  TGLabel*   fInfoLabel5;

public:
  AliEvePMDModuleEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30, UInt_t options = kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~AliEvePMDModuleEditor();

  virtual void SetModel(TObject* obj);
  void DisplayHistos();
  //  void PrintADC();



  // Declare callback/slot methods
  // void DoXYZZ();

  ClassDef(AliEvePMDModuleEditor, 0); // Editor for AliEvePMDModule
}; // endclass AliEvePMDModuleEditor

#endif
