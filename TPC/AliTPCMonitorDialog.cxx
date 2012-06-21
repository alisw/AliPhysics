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

/*
$Log$
Revision 1.1  2007/09/17 10:23:31  cvetan
New TPC monitoring package from Stefan Kniege. The monitoring package can be started by running TPCMonitor.C macro located in macros folder.

*/ 

////////////////////////////////////////////////////////////////////////
////
//// AliTPCMonitorDialog class
////
//// Class to handle dialogs for settings of files and configurations 
//// for the AliTPCMonitor
//// 
//// The dialog will be called by an AliTPCMonitor object from the macro TPCMonitor.C.
//// Depending on the version number passed when creating an Object of this class
//// a certain dialog window (see constructor)  for the TPCMonitor will be opened.
//// The data inserted will be passed to the configuration handler and read out
//// by the monitor object or a the monitor member function will be directly called.
////  
//// Author: Stefan Kniege, IKF, Frankfurt
////       
////
/////////////////////////////////////////////////////////////////////////




#include <cstdlib>
#include "AliTPCMonitorDialog.h"
#include "TGTab.h"
#include "TGButton.h"
#include "TGLabel.h"
#include "TGListBox.h"
#include "TGLayout.h"
#include "TGTextBuffer.h"
#include "TGTextEntry.h" 
#include "TGWindow.h"
#include "TVirtualPadEditor.h"
#include "TTimer.h"
#include <Riostream.h>
#include "Rtypes.h"
#include "AliLog.h"

using std::string;

ClassImp(AliTPCMonitorDialog)
//_____________________________________________________________________________________________
AliTPCMonitorDialog::AliTPCMonitorDialog(TGWindow *p, TGWindow *main, UInt_t w,
                                         UInt_t h, UInt_t options, Int_t version,AliTPCMonitor* monitor):
  //fQObject(0),
fFrameMain(new TGTransientFrame(p, main, w, h, options)),
fFrameComp(0),
fFrameHor(new TGHorizontalFrame(fFrameMain, 60, 20, kFixedWidth)),
fFrameGroup(0),
fOkButton(new TGTextButton(fFrameHor, "&Ok", 1)),
fListBox(0),
fTab(new TGTab(fFrameMain, 300, 300)),
fLayout1(new TGLayoutHints(kLHintsTop    | kLHintsLeft | kLHintsExpandX, 2, 2, 2, 2)),
fLayout2(new TGLayoutHints(kLHintsBottom | kLHintsRight                , 2, 2, 5, 1)),
fLayout3(new TGLayoutHints(kLHintsTop    | kLHintsLeft, 5, 5, 5, 5)),
fMonitor(monitor),
fWidth(w),
fHeight(h),
fOptions(options),
fVersion(version),
fClient(p),
fMain(main)
{
  // Constructor for Dialog window.
  // Create a dialog window.depending on the version it is called with..
  // Verrion 0: Choose DATA Format
  // Version 1: Choose FEC components to display
  // Version 2: Choose Ranges for base and max adc calculation
  
  fFrameMain->Connect("CloseWindow()", "AliTPCMonitorDialog", this, "DoClose()");
  fFrameMain->DontCallClose();
  fFrameMain->SetCleanup(kDeepCleanup);
  fOkButton->Connect("Clicked()", "AliTPCMonitorDialog", this, "DoOK()");
  fFrameHor->AddFrame(fOkButton, fLayout1);
  fFrameHor->Resize(150, fOkButton->GetDefaultHeight());
  fFrameMain->AddFrame(    fFrameHor, fLayout2);
  fTab->Connect("Selected(Int_t)", "AliTPCMonitorDialog", this, "DoTab(Int_t)");
  CreateDialogVersion(version);
}

//_____________________________________________________________________________
AliTPCMonitorDialog::AliTPCMonitorDialog(const AliTPCMonitorDialog &dialog) :
TNamed(dialog.GetName(),dialog.GetTitle()),
//fQObject(0),
fFrameMain(new TGTransientFrame(dialog.GetClient(), dialog.GetMainFrame(), dialog.GetWidth(), dialog.GetHeight())),
fFrameComp(0),
fFrameHor(new TGHorizontalFrame(fFrameMain, 60, 20, kFixedWidth)),
fFrameGroup(0),
fOkButton(new TGTextButton(fFrameHor, "&Ok", 1)),
fListBox(0),
fTab(new TGTab(fFrameMain, 300, 300)),
fLayout1(new TGLayoutHints(kLHintsTop    | kLHintsLeft | kLHintsExpandX, 2, 2, 2, 2)),
fLayout2(new TGLayoutHints(kLHintsBottom | kLHintsRight                , 2, 2, 5, 1)),
fLayout3(new TGLayoutHints(kLHintsTop    | kLHintsLeft, 5, 5, 5, 5)),
fMonitor(dialog.GetMonitor()),
fWidth(dialog.GetWidth()),
fHeight(dialog.GetHeight()),
fOptions(dialog.GetOptions()),
fVersion(dialog.GetVersion()),
fClient(dialog.GetClient()),
fMain(dialog.GetMainFrame())
{
  // copy constructor (actually none forseen for this class)
  fFrameMain->Connect("CloseWindow()", "AliTPCMonitorDialog", this, "DoClose()");
  fFrameMain->DontCallClose();
  fFrameMain->SetCleanup(kDeepCleanup);
  fOkButton->Connect("Clicked()", "AliTPCMonitorDialog", this, "DoOK()");
  fFrameHor->AddFrame(fOkButton, fLayout1);
  fFrameHor->Resize(150, fOkButton->GetDefaultHeight());
  fFrameMain->AddFrame(    fFrameHor, fLayout2);
  fTab->Connect("Selected(Int_t)", "AliTPCMonitorDialog", this, "DoTab(Int_t)");
  CreateDialogVersion(dialog.GetVersion());
  
}

 //_____________________________________________________________________________
AliTPCMonitorDialog &AliTPCMonitorDialog::operator =(const AliTPCMonitorDialog& dialog) 
{
  // assignement operator (actually none forseen for this class)
  if(this!=&dialog)
  {
    fFrameMain  = new TGTransientFrame(dialog.GetClient(), dialog.GetMainFrame(), dialog.GetWidth(), dialog.GetHeight());
    fFrameComp  = 0;
    fFrameHor   = new TGHorizontalFrame(fFrameMain, 60, 20, kFixedWidth);
    fFrameGroup = 0;
    fOkButton   =  new TGTextButton(fFrameHor, "&Ok", 1);
    fListBox    = 0;
    fTab        = new TGTab(fFrameMain, 300, 300);
    fLayout1    = new TGLayoutHints(kLHintsTop    | kLHintsLeft | kLHintsExpandX, 2, 2, 2, 2);
    fLayout2    = new TGLayoutHints(kLHintsBottom | kLHintsRight                , 2, 2, 5, 1);
    fLayout3    = new TGLayoutHints(kLHintsTop    | kLHintsLeft, 5, 5, 5, 5);
    fMonitor    = dialog.GetMonitor();
    fWidth      = dialog.GetWidth();
    fHeight     = dialog.GetHeight();
    fOptions    = dialog.GetOptions();
    fClient     = dialog.GetClient();
    fVersion    = dialog.GetVersion();
    fMain       = dialog.GetMainFrame();
    fFrameMain->Connect("CloseWindow()", "AliTPCMonitorDialog", this, "DoClose()");
    fFrameMain->DontCallClose();
    fFrameMain->SetCleanup(kDeepCleanup);
    fOkButton->Connect("Clicked()", "AliTPCMonitorDialog", this, "DoOK()");
    fFrameHor->AddFrame(fOkButton, fLayout1);
    fFrameHor->Resize(150, fOkButton->GetDefaultHeight());
    fFrameMain->AddFrame(    fFrameHor, fLayout2);
    fTab->Connect("Selected(Int_t)", "AliTPCMonitorDialog", this, "DoTab(Int_t)");
    CreateDialogVersion(dialog.GetVersion());
  }
  return *this;
}



//_____________________________________________________________________________________________
void AliTPCMonitorDialog::CreateDialogVersion(Int_t version)
{
  // Create a dialog window depending on the version passed
  
  if(version==1)
  {
      // Choose DATA Format ////////////////////////////////////////////////////
    TGTextButton *bt;
    TGCompositeFrame* tf = fTab->AddTab("Select Format");
    fFrameComp = new TGCompositeFrame(tf, 60, 20, kVerticalFrame);
    fFrameComp->AddFrame(bt = new TGTextButton(fFrameComp, "Select Entry", 1), fLayout3);
    
    bt->Connect("Clicked()", "AliTPCMonitorDialog", this, "HandleButtons()");
    fFrameComp->AddFrame(fListBox = new TGListBox(fFrameComp, 1), fLayout3);
    tf->AddFrame(fFrameComp, fLayout3);
    fListBox->AddEntry("ROOT  file        ", 0);
    fListBox->AddEntry("DATE  file        ", 1);
    fListBox->AddEntry("DATE  file/stream ", 2);
    
    fListBox->Resize(200, 80);
    
  }
  else if(version==2)
  {
      // choose component to be displayed /////////////////////////////////////
    fFrameGroup = new TGGroupFrame(fFrameMain, "Components", kVerticalFrame);
    fFrameGroup->SetTitlePos(TGGroupFrame::kRight); // right aligned
    fFrameMain->AddFrame(fFrameGroup, fLayout3);
    fFrameGroup->SetLayoutManager(new TGMatrixLayout(fFrameGroup, 0, 2, 10));
    
    const Char_t names[6][40] ={
      "FEC (global)",
        "RCU (patch) ",
        "Branch      ",
        "FEC (local) ",
        "Connector   ",
        "Chip        "
    };
    for (Int_t j = 0; j < 6; j++)
    {
      fFrameGroup->AddFrame(new TGLabel(fFrameGroup, new TGHotString(names[j])));
      
      fBuf[j] = new TGTextBuffer(10);
      fBuf[j]->AddText(0, "all");
      
      fEnt[j] = new TGTextEntry(fFrameMain, fBuf[j]);
      fEnt[j]->Resize(50, fEnt[j]->GetDefaultHeight());
      fEnt[j]->SetFont("-adobe-courier-bold-r-*-*-14-*-*-*-*-*-iso8859-1");
      fFrameGroup->AddFrame(fEnt[j]);
      
    }
    TGTextButton* bt2 =0;
    fFrameGroup->AddFrame(bt2 = new TGTextButton(fFrameGroup, "Select Config", 2), fLayout3);
    bt2->Connect("Clicked()", "AliTPCMonitorDialog", this, "HandleButtons()");
    fFrameGroup->Resize();
  }
  else if(version==3)
  {
      // choose component to be displayed /////////////////////////////////////
    fFrameGroup = new TGGroupFrame(fFrameMain, "Ranges BASE/MAX/SUM ", kVerticalFrame);
    fFrameGroup->SetTitlePos(TGGroupFrame::kRight);
    fFrameMain->AddFrame(fFrameGroup, fLayout3);
    fFrameGroup->SetLayoutManager(new TGMatrixLayout(fFrameGroup, 0, 2, 10));
    
    const Char_t names[3][40] = {
      "baseline",
        "maximum adc",
        "sum adc ",
    };
    
    for (Int_t j = 0; j < 3; j++) {
      fFrameGroup->AddFrame(new TGLabel(fFrameGroup, new TGHotString(names[j])));
      
      fBuf[j] = new TGTextBuffer(30);
      if(j==0) fBuf[j]->AddText(0,Form("%i,%i",fMonitor->GetRangeBaseMin()  , fMonitor->GetRangeBaseMax()));
      if(j==1) fBuf[j]->AddText(0,Form("%i,%i",fMonitor->GetRangeMaxAdcMin(), fMonitor->GetRangeMaxAdcMax()));
      if(j==2) fBuf[j]->AddText(0,Form("%i,%i",fMonitor->GetRangeSumMin()   , fMonitor->GetRangeSumMax()));
      
      fEnt[j] = new TGTextEntry(fFrameMain, fBuf[j]);
      fEnt[j]->Resize(100, fEnt[j]->GetDefaultHeight());
      fEnt[j]->SetFont("-adobe-courier-bold-r-*-*-14-*-*-*-*-*-iso8859-1");
      fFrameGroup->AddFrame(fEnt[j]);
      
    }
    TGTextButton* bt2 =0;
    fFrameGroup->AddFrame(bt2 = new TGTextButton(fFrameGroup, "Select Config", 3), fLayout3);
    bt2->Connect("Clicked()", "AliTPCMonitorDialog", this, "HandleButtons()");
    fFrameGroup->Resize(100,300);
  }
  
  TGLayoutHints *fL5 = new TGLayoutHints(kLHintsBottom | kLHintsExpandX |
                                         kLHintsExpandY, 2, 2, 5, 1);
  fFrameMain->AddFrame(fTab, fL5);
  fFrameMain->MapSubwindows();
  fFrameMain->Resize();
  
  fFrameMain->CenterOnParent();
  if(     version==3)  fFrameMain->SetWindowName("Select Ranges");
  else if(version==1)  fFrameMain->SetWindowName("Select Component");
  else                 fFrameMain->SetWindowName("Select Format");
  fFrameMain->MapWindow();
}

//_____________________________________________________________________________________________
AliTPCMonitorDialog::~AliTPCMonitorDialog()
{
  // Delete test dialog widgets.
  fFrameMain->DeleteWindow();
}

//_____________________________________________________________________________________________
void AliTPCMonitorDialog::HandleButtons(Int_t id)
{
  // Handle different buttons.
  if (id == -1) {
    TGButton *btn = (TGButton*) gTQSender;
    id = btn->WidgetId();
  }
  
  switch (id)
  {
  case 1:
    {
  // Select Format
      Int_t fsel = fListBox->GetSelected();
      if(fsel==-1) {  AliError("Warning: No Format selected "); }
      TGLBEntry* libentr = (TGLBEntry*)fListBox->GetEntry(fsel);
      AliInfo(Form("Selected Format : %i : %s",fsel,libentr->GetTitle()));
      fMonitor->SetFormat(fsel);
      
      break;
    }
  case 2:
    {
  // Select FEE Component
      Int_t compsel=1;
      Int_t fComponents[10];
      for(Int_t j = 0; j < 6; j++)
      {
        TString s1 =  fEnt[j]->GetDisplayText();
        if(s1=="all"){  fComponents[j] = -1;	      }
        else         {  fComponents[j] = (int)atoi(s1.Data()); }
        if(
            (j==0 && (fComponents[j] < -1 || fComponents[j]>5000)) ||
            (j==1 && (fComponents[j] < -1 || fComponents[j]>5   )) ||
            (j==2 && (fComponents[j] < -1 || fComponents[j]>1   )) ||
            (j==3 && (fComponents[j] < -1 || fComponents[j]>25  )) ||
            (j==4 && (fComponents[j] < -1 || fComponents[j]>6   )) ||
            (j==5 && (fComponents[j] < -1 || fComponents[j]>8   ))    )
        {
          compsel =0;
          AliError("Settings out of range ( version 2) ");
        }
      }
      if(compsel==1)
      {
        if(fMonitor!=0 && fMonitor->GetEventProcessed() )fMonitor->ShowSel((Int_t*)fComponents);
        else                                             AliError("No event processed up to now");
      }
      break;
    }
  case 3:
    {
  // Set Configuration (base, max, sum)
      Int_t fConfigArr[6];
      
      string s1((TString)fEnt[0]->GetDisplayText().Data());
      string s2((TString)fEnt[1]->GetDisplayText().Data());
      string s3((TString)fEnt[2]->GetDisplayText().Data());
      
      Int_t error = 1;
      
      if(s1.find(",")!=string::npos){
        string sub11 = s1.substr(0,s1.find(","))               ;fConfigArr[0] = atoi(sub11.data());
        string sub12 = s1.substr(s1.find(",")+1,s1.length() )  ;fConfigArr[1] = atoi(sub12.data());
        if(fConfigArr[0]<fConfigArr[1] && fConfigArr[1]<1024){
          fMonitor->SetRangeBase(fConfigArr[0],fConfigArr[1]);
          error=0;
        }
      }
      if(s2.find(",")!=string::npos){
        string sub21 = s2.substr(0,s2.find(","))               ;fConfigArr[2] = atoi(sub21.data());
        string sub22 = s2.substr(s2.find(",")+1,s2.length() )  ;fConfigArr[3] = atoi(sub22.data());
        if(fConfigArr[2]<fConfigArr[3] && fConfigArr[3]<1024){
          fMonitor->SetRangeMax(fConfigArr[2],fConfigArr[3]);
          error=0;
        }
      }
      
      if(s3.find(",")!=string::npos){
        string sub31 = s3.substr(0,s3.find(","))               ;fConfigArr[4] = atoi(sub31.data());
        string sub32 = s3.substr(s3.find(",")+1,s3.length() )  ;fConfigArr[5] = atoi(sub32.data());
        if(fConfigArr[4]<fConfigArr[5] && fConfigArr[5]<1024){
          fMonitor->SetRangeSum(fConfigArr[4],fConfigArr[5]);
          error=0;
        }
      }
      if(error!=0) AliError("Settings out of range (version 3)");
      break;
    }
  default:
    break;
  }
}

//_____________________________________________________________________________________________
void AliTPCMonitorDialog::DoTab(Int_t id) 
{
  printf("Tab item %d activated\n", id);
}

//_____________________________________________________________________________________________
void AliTPCMonitorDialog::DoClose() const
{
  // Close the dialog window
  printf("\nTerminating dialog: via window manager\n");
  CloseWindow();
  // Close the Ged editor if it was activated.
  if (TVirtualPadEditor::GetPadEditor(kFALSE) != 0)
    TVirtualPadEditor::Terminate();
}

//_____________________________________________________________________________________________
void AliTPCMonitorDialog::CloseWindow() const
{
  // Called when window is closed via the window manager.
  delete this;
}

//_____________________________________________________________________________________________
void AliTPCMonitorDialog::DoOK()
{
  // Call from OK button
  printf("\nTerminating dialog: OK pressed\n");
  
  // The same effect can be obtained by using a singleshot timer:
  TTimer::SingleShot(150, "AliTPCMonitorDialog", this, "CloseWindow()");
  
  // Close the Ged editor if it was activated.
  if (TVirtualPadEditor::GetPadEditor(kFALSE) != 0)
    TVirtualPadEditor::Terminate();
}


//_____________________________________________________________________________________________
void AliTPCMonitorDialog::DoCancel()
{
  // Cancel operation and close the dialog window
  printf("\nTerminating dialog: Cancel pressed\n");
  TTimer::SingleShot(150, "AliTPCMonitorDialog", this, "CloseWindow()");
  
  // Close the Ged editor if it was activated.
  if (TVirtualPadEditor::GetPadEditor(kFALSE) != 0)
    TVirtualPadEditor::Terminate();
}

