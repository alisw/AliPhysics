#ifndef ALIMONITORCONTROL_H
#define ALIMONITORCONTROL_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TGFrame.h>
#include <TGMenu.h>
#include <TGButton.h>
#include <TGLabel.h>
#include <TGTextEntry.h>
#include <RQ_OBJECT.h>
#include <TSysEvtHandler.h>
#include <TTimer.h>
#include "AliMonitorProcess.h"


class AliMonitorControl : public TObject {

RQ_OBJECT("AliMonitorControl")

public:
  AliMonitorControl(AliMonitorProcess* process);
  virtual ~AliMonitorControl();

  void               HandleMenu(Int_t id);
  void               DoReset();
  void               DoStartStop();

private:
  virtual Bool_t     HandleTimer(TTimer* timer);

  void               UpdateStatus();

  AliMonitorProcess* fMonitorProcess;

  ULong_t            fColorStatus;
  ULong_t            fColorStart;
  ULong_t            fColorStop;

  TGMainFrame*       fMain;

  TGLayoutHints*     fMenuBarLayout;
  TGLayoutHints*     fMenuBarItemLayout;
  TGLayoutHints*     fMenuBarHelpLayout;
  TGPopupMenu*       fMenuFile;
  TGPopupMenu*       fMenuOptions;
  TGPopupMenu*       fMenuHelp;
  TGMenuBar*         fMenuBar;

  TGLayoutHints*     fFrameLayout;
  TGVerticalFrame*   fFrame;
  TGLayoutHints*     fStatusLayout;
  TGLayoutHints*     fStatusFrameLayout;

  TGHorizontalFrame* fStatus1Frame;
  TGLabel*           fRunNumberLabel;
  TGTextEntry*       fRunNumber;
  TGLabel*           fEventNumberLabel;
  TGTextEntry*       fEventNumber;

  TGHorizontalFrame* fStatus2Frame;
  TGLabel*           fStatusLabel;
  TGTextEntry*       fStatus;

  TGHorizontalFrame* fStatus3Frame;
  TGLabel*           fEventsLabel;
  TGTextEntry*       fEvents;
  TGLabel*           fClientsLabel;
  TGTextEntry*       fClients;

  TGLayoutHints*     fButtonFrameLayout;
  TGHorizontalFrame* fButtonFrame;
  TGLayoutHints*     fButtonLayout;
  TGTextButton*      fResetButton;
  TGTextButton*      fStartStopButton;
  Bool_t             fStartButtonStatus;

  Bool_t             fTerminating;

  TTimer*            fTimer;

  ClassDef(AliMonitorControl, 0)   // class for controlling the AliMonitorProcess
};
 

#endif









