#ifndef ALIMONITORDIALOG_H
#define ALIMONITORDIALOG_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TGFrame.h>
#include <TGButton.h>
#include <TGLabel.h>
#include <RQ_OBJECT.h>


class AliMonitorDialog : public TObject {

RQ_OBJECT("AliMonitorDialog")

public:
  AliMonitorDialog(TGFrame* main, Int_t width = 300, Int_t height = 80,
		   Bool_t cancelBtn = kTRUE);
  virtual ~AliMonitorDialog();

  void               CloseWindow();
  void               DoOk();
  virtual void       OnOkClicked() {};
  void               DoCancel();
  virtual void       OnCancelClicked() {};

protected:
  TGTransientFrame*  fMain;
  TGLayoutHints*     fFrameLayout;
  TGHorizontalFrame* fFrame;
  TGLayoutHints*     fButtonFrameLayout;
  TGHorizontalFrame* fButtonFrame;
  TGLayoutHints*     fButtonLayout;
  TGTextButton*      fOkButton;
  TGTextButton*      fCancelButton;

  ClassDef(AliMonitorDialog, 0)   // base class for dialogs with ok and cancel button
};
 

#endif









