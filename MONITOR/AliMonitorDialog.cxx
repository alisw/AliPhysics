/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  This is a base class for dialogs with ok and cancel button               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "AliMonitorDialog.h"
#include "AliLog.h"
#include <TGFrame.h>
#include <TGButton.h>


ClassImp(AliMonitorDialog) 


//_____________________________________________________________________________
AliMonitorDialog::AliMonitorDialog(TGFrame* main, Int_t width, Int_t height,
				   Bool_t cancelBtn) :
  TObject(),
  fQObject(),
  fMain(new TGTransientFrame(gClient->GetRoot(), main, width, height)),
  fFrameLayout(new TGLayoutHints(kLHintsCenterX | kLHintsCenterY, 2, 2, 2, 2)),
  fFrame(new TGHorizontalFrame(fMain, 0, 0)),
  fButtonFrameLayout(new TGLayoutHints(kLHintsCenterX | kLHintsCenterY, 2, 2, 2, 2)),
  fButtonFrame(new TGHorizontalFrame(fMain, 0, 0)),
  fButtonLayout(new TGLayoutHints(kLHintsCenterX | kLHintsCenterY, 15, 15, 2, 2)),
  fOkButton(new TGTextButton(fButtonFrame, "    &Ok    ", 1)),
  fCancelButton(NULL)
{
// create a dialog with an ok and a cancel button

  fMain->AddFrame(fFrame, fFrameLayout);

  fMain->AddFrame(fButtonFrame, fButtonFrameLayout);

  fButtonFrame->AddFrame(fOkButton, fButtonLayout);
  fOkButton->Connect("Clicked()", "AliMonitorDialog", this, "DoOk()");
  if (cancelBtn) {
    fCancelButton = new TGTextButton(fButtonFrame, " &Cancel ", 2);
    fButtonFrame->AddFrame(fCancelButton, fButtonLayout);
    fCancelButton->Connect("Clicked()", "AliMonitorDialog", this, 
			   "DoCancel()");
  }

  fMain->Connect("CloseWindow()", "AliMonitorDialog", this, "CloseWindow()");
  Int_t x;
  Int_t y;
  Window_t child;
  gVirtualX->TranslateCoordinates(main->GetId(), gClient->GetRoot()->GetId(),
				  (main->GetWidth() - width) >> 1,
				  (main->GetHeight() - height) >> 1, 
				  x, y, child);
  if (x < 0) x = 0;
  if (y < 0) y = 0;
  fMain->Move(x, y);
  fMain->SetWMPosition(x, y);
  fMain->SetWMSize(width, height);
  fMain->SetWMSizeHints(width, height, width, height, 0, 0);
  fMain->MapSubwindows();
  fMain->Layout();
  fMain->MapWindow();
}

//_____________________________________________________________________________
AliMonitorDialog::~AliMonitorDialog()
{
// clean up

  delete fOkButton;
  delete fCancelButton;
  delete fButtonLayout;
  delete fButtonFrame;
  delete fButtonFrameLayout;
  delete fFrame;
  delete fFrameLayout;

  delete fMain;
}

//_____________________________________________________________________________
void AliMonitorDialog::CloseWindow() const
{
// called when the window is closed

  delete this;
}

//_____________________________________________________________________________
void AliMonitorDialog::DoOk()
{
// called when the ok button is clicked

  OnOkClicked();
  fMain->SendCloseMessage();
}

//_____________________________________________________________________________
void AliMonitorDialog::DoCancel()
{
// called when the cancel button is clicked

  OnCancelClicked();
  fMain->SendCloseMessage();
}
