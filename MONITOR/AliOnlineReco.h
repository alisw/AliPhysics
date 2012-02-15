// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef ALIONLINERECO_H
#define ALIONLINERECO_H

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

  void ChildProcTerm(Int_t pid, Int_t status); // *SIGNAL*
  void ExitLoopChildProcTerm();

  //------------------------------------------------------------------------------
  // Handlers of button signals.
  //------------------------------------------------------------------------------

  void DoAutoRun();
  void DoStart();
  void DoStop();
  void DoExit();

  virtual void CloseWindow();

  Int_t RetrieveGRP(UInt_t run, TString &gdc);

  // Things that should be private but have to be public for CINT.
  void AutoRunTimerTimeout();

private:
  AliOnlineReco(const AliOnlineReco&);            // Not implemented
  AliOnlineReco& operator=(const AliOnlineReco&); // Not implemented

  // GUI components.
  TGListBox     *fRunList;    // List-box for listing current runs.
  TGCheckButton *fAutoRun;    // Check-box toggling auto-run when a new run starts.
  TGTextButton  *fStartButt;  // Start for selected run.
  TGTextButton  *fStopButt;   // Stop for selected run.
  TGTextButton  *fExitButt;   // Exit button.

  // DIM interface. Could do without members and just leak them ...
  AliDimIntNotifier *fSOR[5]; // DIM listeners for SOR.
  AliDimIntNotifier *fEOR[5]; // DIM listeners for EOR.

  // AutoRun state and timer
  TTimer        *fAutoRunTimer;     // Timer for auto-run on new run.
  Int_t          fAutoRunScheduled; // Run for which auto-run is scheduled.
  Int_t          fAutoRunRunning;   // Run for which auto-run was executed.

  // Run-state, process mngmnt
  typedef std::map<Int_t, Int_t> mIntInt_t; // value should be struct { pid, state, ... };
  typedef mIntInt_t::iterator    mIntInt_i;

  mIntInt_t      fRun2PidMap;  // Map from run-number to process id.

  Bool_t         fTestMode;    // Flag for test mode (run alitestproc instead of alieve).
  Bool_t         fDoExit;     // Flag for exit mode

  mIntInt_i FindMapEntryByPid(Int_t pid);

  void StartAliEve(mIntInt_i& mi);
  void KillPid(Int_t pid);

  void StartAutoRunTimer(Int_t run);
  void StopAutoRunTimer();

  ClassDef(AliOnlineReco, 0);
};

#endif
