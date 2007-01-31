#ifndef ALIMONITORDIALOG_H
#define ALIMONITORDIALOG_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <RQ_OBJECT.h>
#include <TObject.h>

class TGFrame;
class TGTransientFrame;
class TGLayoutHints;
class TGHorizontalFrame;
class TGTextButton;


class AliMonitorDialog : public TObject {

RQ_OBJECT("AliMonitorDialog")

public:
  AliMonitorDialog(TGFrame* main, Int_t width = 300, Int_t height = 80,
		   Bool_t cancelBtn = kTRUE);
  virtual ~AliMonitorDialog();

  void               CloseWindow() const;
  void               DoOk();
  virtual void       OnOkClicked() {};
  void               DoCancel();
  virtual void       OnCancelClicked() {};

protected:
  TGTransientFrame*  fMain;                // the main window
  TGLayoutHints*     fFrameLayout;         // layout of the main frame
  TGHorizontalFrame* fFrame;               // the main frame
  TGLayoutHints*     fButtonFrameLayout;   // layout of the buttons frame
  TGHorizontalFrame* fButtonFrame;         // the frame for buttons
  TGLayoutHints*     fButtonLayout;        // layout of the buttons
  TGTextButton*      fOkButton;            // the Ok button
  TGTextButton*      fCancelButton;        // the cancel button

 private:
  AliMonitorDialog(const AliMonitorDialog& dlg);
  AliMonitorDialog& operator = (const AliMonitorDialog& dlg);

  ClassDef(AliMonitorDialog, 0)   // base class for dialogs with ok and cancel button
};
 

#endif









