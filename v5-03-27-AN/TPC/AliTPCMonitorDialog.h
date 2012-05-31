#ifndef ALITPCMONITORDIALOG_H
#define ALITPCMONITORDIALOG_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////
////
//// AliTPCMonitorDialog class
////
//// Class to handle dialogs for settings of files and configurations 
//// for the AliTPCMonitor
//// 
//// Author: Stefan Kniege, IKF, Frankfurt
////       
////
/////////////////////////////////////////////////////////////////////////


#include "TGFrame.h"  
#include "AliTPCMonitor.h" 
class TGWindow;
class TGTransientFrame;
class TGLayoutHints;
class TGListBox;
class AliTPCMonitor;
class TGTab;
class TGButton;
class TGHorizontalFrame; 
class TGTextBuffer;
class TGTextEntry;
class TGTransientFrame;
class TGCompositeFrame;
class TGGroupFrame;
class AliTPCMonitorDialog : public TNamed
{
    
public:
    AliTPCMonitorDialog(TGWindow *p, TGWindow *main, UInt_t w, UInt_t h, UInt_t options = kVerticalFrame, Int_t version =1, AliTPCMonitor* monitor =0);
    AliTPCMonitorDialog(const  AliTPCMonitorDialog &dialog);
    AliTPCMonitorDialog& operator= (const AliTPCMonitorDialog& dialog);
    virtual ~AliTPCMonitorDialog();
    
    void           DoClose()     const;
    void           CloseWindow() const;
    void           DoOK();
    void           DoCancel(); 
    void           DoTab(Int_t id);
    void           HandleButtons(Int_t id = -1);
    void           CreateDialogVersion(Int_t version);
    UInt_t         GetWidth()     const { return fWidth;   }
    UInt_t         GetHeight()    const { return fHeight;  }
    UInt_t         GetOptions()   const { return fOptions; }
    AliTPCMonitor* GetMonitor()   const { return fMonitor; }
    TGWindow*      GetClient()	  const { return fClient;  }
    TGWindow*      GetMainFrame() const { return fMain;    }
    Int_t          GetVersion()   const { return fVersion; }
 private:
    
    TGTransientFrame*    fFrameMain;    // Frames for dialog windows      
    TGCompositeFrame*    fFrameComp;    // Frames for dialog windows 
    TGHorizontalFrame*   fFrameHor;     // Frames for dialog windows 
    TGGroupFrame*        fFrameGroup;   // Frames for dialog windows 
    TGButton*            fOkButton;     // Ok button for windows
    TGListBox*           fListBox;      // ListBox for entries to be selected
    TGTab*               fTab;          // Tabs for several pages in one window 
    TGLayoutHints*       fLayout1;      // Layout for window version 1 
    TGLayoutHints*       fLayout2;      // Layout for window versoin 2  
    TGLayoutHints*       fLayout3;      // Layout for window version 3 
    
  

    TGTextBuffer*        fBuf[7];       // Text buffer for GroupFrame
    TGTextEntry*         fEnt[7];       // Text entries for GroupFrame
    
    AliTPCMonitor*       fMonitor;      // Pointer to AliTPCMonitor to be called
    
    UInt_t               fWidth;        // Width of the window
    UInt_t               fHeight;       // Height of the window 
    UInt_t               fOptions;      // Frame options 
    Int_t                fVersion;       // Version of the dialog
    TGWindow*            fClient;       // TGClient of Dialog
    TGWindow*            fMain;         // Main fraem for Dialog window
    ClassDef(AliTPCMonitorDialog,1);
    
};



#endif
