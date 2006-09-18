#ifndef ALIMONITORCONTROL_H
#define ALIMONITORCONTROL_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TObject.h>
#include <RQ_OBJECT.h>
#include "AliMonitorDialog.h"

class AliMonitorProcess;
class TTimer;
class TGMainFrame;
class TGLayoutHints;
class TGPopupMenu;
class TGMenuBar;
class TGVerticalFrame;
class TGHorizontalFrame;
class TGLabel;
class TGTextEntry;
class TGTextButton;
class TGNumberEntry;
class TGTextView;


class AliMonitorControl : public TObject {

RQ_OBJECT("AliMonitorControl")

public:
  AliMonitorControl(AliMonitorProcess* process);
  virtual ~AliMonitorControl();

  void               HandleMenu(Int_t id);
  void               DoReset();
  void               DoStartStop();

private:
  AliMonitorControl(const AliMonitorControl& control);
  AliMonitorControl& operator = (const AliMonitorControl& control);

  virtual Bool_t     HandleTimer(TTimer* timer);

  void               UpdateStatus();

  AliMonitorProcess* fMonitorProcess;       // the controlled monitor process

  ULong_t            fColorStatus;          // color for status info
  ULong_t            fColorStart;           // color for start button
  ULong_t            fColorStop;            // color for stop button

  TGMainFrame*       fMain;                 // the main window

  TGLayoutHints*     fMenuBarLayout;        // layout of the menu bar
  TGLayoutHints*     fMenuBarItemLayout;    // layout of the menu items
  TGLayoutHints*     fMenuBarHelpLayout;    // layout of the help menu
  TGPopupMenu*       fMenuFile;             // the file menu
  TGPopupMenu*       fMenuOptions;          // the options menu
  TGPopupMenu*       fMenuHelp;             // the help menu
  TGMenuBar*         fMenuBar;              // the menu bar

  TGLayoutHints*     fFrameLayout;          // layout of the main frame
  TGVerticalFrame*   fFrame;                // the main frame
  TGLayoutHints*     fStatusLayout;         // layout of status info
  TGLayoutHints*     fStatusFrameLayout;    // layout of status frames

  TGHorizontalFrame* fStatus1Frame;         // frame for run/event number
  TGLabel*           fRunNumberLabel;       // label for run number
  TGTextEntry*       fRunNumber;            // run number display
  TGLabel*           fEventNumberLabel;     // label for event number
  TGTextEntry*       fEventNumber;          // event number display

  TGHorizontalFrame* fStatus2Frame;         // frame for current status
  TGLabel*           fStatusLabel;          // label for status
  TGTextEntry*       fStatus;               // current status display

  TGHorizontalFrame* fStatus3Frame;         // frame for number of event/clients
  TGLabel*           fEventsLabel;          // label for number of events
  TGTextEntry*       fEvents;               // number of monitored events display
  TGLabel*           fClientsLabel;         // label for number of clients
  TGTextEntry*       fClients;              // number of clients display

  TGLayoutHints*     fButtonFrameLayout;    // layout of frame with buttons
  TGHorizontalFrame* fButtonFrame;          // frame for buttons
  TGLayoutHints*     fButtonLayout;         // layout of buttons
  TGTextButton*      fResetButton;          // the rest button
  TGTextButton*      fStartStopButton;      // the start/stop button
  Bool_t             fStartButtonStatus;    // current status of the start/stop button

  Bool_t             fTerminating;          // true if program will be terminated

  TTimer*            fTimer;                // timer for X update


  class AliMonitorBufferDlg : public AliMonitorDialog {

  public:
    AliMonitorBufferDlg(Int_t& size, TGFrame* main);
    virtual ~AliMonitorBufferDlg();

    virtual void       OnOkClicked();

  private:
    AliMonitorBufferDlg(const AliMonitorBufferDlg& dlg);
    AliMonitorBufferDlg& operator = (const AliMonitorBufferDlg& dlg);

    TGLayoutHints*     fBufferLayout;       // layout of buffer entry
    TGLabel*           fBufferLabel;        // label for buffer entry
    TGNumberEntry*     fBufferEntry;        // buffer number entry

    Int_t&             fSize;               // result
  };


  class AliMonitorClientsDlg : public AliMonitorDialog {

  public:
    AliMonitorClientsDlg(TObjArray* clients, TGFrame* main);
    virtual ~AliMonitorClientsDlg();

  private:
    AliMonitorClientsDlg(const AliMonitorClientsDlg& dlg);
    AliMonitorClientsDlg& operator = (const AliMonitorClientsDlg& dlg);

    TGLayoutHints*     fClientsLayout;      // layout of clients list
    TGTextView*        fClients;            // list of clients
  };


  ClassDef(AliMonitorControl, 0)   // class for controlling the AliMonitorProcess
};
 

#endif









