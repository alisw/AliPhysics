// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliOnlineReco_H
#define AliOnlineReco_H

#include "TObject.h"
#include "TString.h"
#include <TGFrame.h>

#include <map>

class AliDimIntNotifier;

class TTimer;

class TGTextButton;
class TGCheckButton;
class TGListBox;

//______________________________________________________________________________
// Short description of AliOnlineReco
//


class AliOnlineReco : public TGMainFrame
{
public:
  AliOnlineReco();
  virtual ~AliOnlineReco();

  AliDimIntNotifier* GetSOR(Int_t i) const { return fSOR[i]; }
  AliDimIntNotifier* GetEOR(Int_t i) const { return fEOR[i]; }

  Int_t  GetLastRun() const;

  Bool_t GetAutoRunMode() const;
  void   SetAutoRunMode(Bool_t ar);

  void SetTestMode() { fTestMode = kTRUE; }

  //------------------------------------------------------------------------------
  // Handlers of DIM signals.
  //------------------------------------------------------------------------------

  void StartOfRun(Int_t run);
  void EndOfRun(Int_t run);

  //------------------------------------------------------------------------------
  // Handlers of OS signals.
  //------------------------------------------------------------------------------

  void ChildProcTerm(Int_t pid, Int_t status);

  //------------------------------------------------------------------------------
  // Handlers of button signals.
  //------------------------------------------------------------------------------

  void DoAutoRun();
  void DoStart();
  void DoStop();
  void DoExit();

  virtual void CloseWindow();

  Int_t RetrieveGRP(UInt_t run, TString &gdc);

private:
  AliOnlineReco(const AliOnlineReco&);            // Not implemented
  AliOnlineReco& operator=(const AliOnlineReco&); // Not implemented

  // GUI components.
  TGListBox     *fRunList;
  TGCheckButton *fAutoRun;
  TGTextButton  *fStartButt;
  TGTextButton  *fStopButt;
  TGTextButton  *fExitButt;

  // DIM interface. Could do without members and just leak them ...
  AliDimIntNotifier *fSOR[5];
  AliDimIntNotifier *fEOR[5];

  // AutoRun state and timer
  TTimer        *fAutoRunTimer;
  Int_t          fAutoRunScheduled;
  Int_t          fAutoRunRunning;

  // Run-state, process mngmnt
  typedef std::map<Int_t, Int_t> mIntInt_t; // value should be struct { pid, state, ... };
  typedef mIntInt_t::iterator    mIntInt_i;

  mIntInt_t      fRun2PidMap;

  Bool_t         fTestMode;

  mIntInt_i FindMapEntryByPid(Int_t pid);

  void StartAliEve(mIntInt_i& mi);
  void KillPid(Int_t pid);

  void StartAutoRunTimer(Int_t run);
  void StopAutoRunTimer();

  // Things that should be private but have to be public for CINT.
public:
  void AutoRunTimerTimeout();

  ClassDef(AliOnlineReco, 0);
};

#endif
