// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveMacroEditor_H
#define AliEveMacroEditor_H

#include "TGedFrame.h"

class AliEveMacro;

class TGCheckButton;
class TGTextEntry;
class TGComboBox;

//______________________________________________________________________________
// Short description of AliEveMacroEditor
//

class AliEveMacroEditor : public TGedFrame
{
public:
  AliEveMacroEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
		    UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~AliEveMacroEditor() {}

  virtual void SetModel(TObject* obj);

  void DoSources(Int_t v);
  void DoTags();
  void DoMacro();
  void DoFunc();
  void DoArgs();
  void DoActive();

protected:
  AliEveMacro           *fM; // Model object.

  TGComboBox            *fSources;
  TGTextEntry           *fTags;
  TGTextEntry           *fMacro;
  TGTextEntry           *fFunc;
  TGTextEntry           *fArgs;
  TGCheckButton         *fActive;

  TGHorizontalFrame* MkHFrame(TGCompositeFrame* p=0);
  TGLabel*           MkLabel (TGCompositeFrame* p, const char* txt, Int_t width,
			      Int_t lo=0, Int_t ro=0, Int_t to=2, Int_t bo=0);

private:
  AliEveMacroEditor(const AliEveMacroEditor&);            // Not implemented
  AliEveMacroEditor& operator=(const AliEveMacroEditor&); // Not implemented

  ClassDef(AliEveMacroEditor, 0); // GUI editor for AliEveMacro.
};

#endif
