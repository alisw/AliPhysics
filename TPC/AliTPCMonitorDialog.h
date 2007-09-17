#ifndef ALITPCMONITORDIALOG_H
#define ALITPCMONITORDIALOG_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////
//
// AliTPCMonitorDialog class
//
// Class to handle dialogs for settings of files and configurations 
// for the AliTPCMonitor
// 
// Author: Stefan Kniege, IKF, Frankfurt
//       
//
/////////////////////////////////////////////////////////////////////////



#include <iostream>
#include "TGWindow.h"
#include "TRootGuiBuilder.h"
#include "TGMenu.h"
#include "TGButtonGroup.h"
#include "TGDockableFrame.h"
#include "TGToolBar.h"
#include "TGButton.h"
#include "TGToolTip.h"
#include "TGuiBldDragManager.h"
#include "TGMdiMainFrame.h"
#include "TGLabel.h"
#include "TG3DLine.h"
#include "TGNumberEntry.h"
#include "TGuiBldHintsButton.h"
#include "TGuiBldHintsEditor.h"
#include "TGuiBldEditor.h"
#include "TGShutter.h"
#include "TGCanvas.h"
#include "TGStatusBar.h"
#include "TGTab.h"
#include "TGSplitter.h"
#include "TGMdiMenu.h"
#include "TGMdiDecorFrame.h"
#include "TGMdiFrame.h"
#include "TGListBox.h"
#include "TGListView.h"
#include "TGFSContainer.h"
#include "TGFSComboBox.h"
#include "TGFileDialog.h"
#include "TGScrollBar.h"
#include "TGMsgBox.h"
#include "TGLayout.h"  
#include "TGFrame.h"  
#include "TGTextEdit.h"
#include "RQ_OBJECT.h"
#include "TRootEmbeddedCanvas.h"
#include "TH2F.h"
#include "TRandom.h"
#include "TVirtualPadEditor.h"
#include "TSystem.h"
#include "Rtypes.h"
#include "AliTPCMonitor.h" 
#include "AliLog.h"

using namespace std;

 
class TGHorizontalFrame; 
class AliLog;
class AliTPCMonitorDialog : public TNamed{
    
    RQ_OBJECT("AliTPCMonitorDialog")
    
    
 public:
    AliTPCMonitorDialog(const TGWindow *p, const TGWindow *main, UInt_t w, UInt_t h, UInt_t options = kVerticalFrame, Int_t version =1, AliTPCMonitor* monitor =0);
    virtual ~AliTPCMonitorDialog();
    
    void DoClose();
    void CloseWindow();
    void DoOK();
    void DoCancel();
    void DoTab(Int_t id);
    void HandleButtons(Int_t id = -1);
    void CreateDialogVersion(Int_t version);
    
 private:
    
    TGTransientFrame*    fFrameMain;    // Frames for dialog windows      
    TGCompositeFrame*    fFrameComp;    // Frames for dialog windows 
    TGHorizontalFrame*   fFrameHor;     // Frames for dialog windows 
    TGGroupFrame*        fFrameGroup;   // Frames for dialog windows 
    TGButton*            fOkButton;     // Ok button for windows
    TGButton*            fCancelButton; // Cancel button for windows
    
    TGListBox*           fListBox;      // ListBox for entries to be selected
    TGTab*               fTab;          // Tabs for several pages in one window 
    TGLayoutHints*       fLayout1;      // Layout for window version 1 
    TGLayoutHints*       fLayout2;      // Layout for window versoin 2  
    TGLayoutHints*       fLayout3;      // Layout for window version 3 
    
    
    TGTextBuffer*        fBuf[7];       // Text buffer for GroupFrame
    TGTextEntry*         fEnt[7];       // Text entries for GroupFrame
    
    AliTPCMonitor*       fMonitor;      // Pointer to AliTPCMonitor to be called
    
    ClassDef(AliTPCMonitorDialog,1);
    
};



#endif
