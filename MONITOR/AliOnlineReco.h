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
#include <TGFrame.h>

#include <map>

class AliDimIntNotifier;

class TGTextButton;
class TGListBox;

//______________________________________________________________________________
// Short description of AliOnlineReco
//


class AliOnlineReco : public TGMainFrame
{
public:
  AliOnlineReco();
  virtual ~AliOnlineReco() {}

  AliDimIntNotifier* GetSOR() const { return fSOR; }
  AliDimIntNotifier* GetEOR() const { return fEOR; }

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

  void DoStart();
  void DoStop();
  void DoXyzz();

  virtual void CloseWindow();

private:
  AliOnlineReco(const AliOnlineReco&);            // Not implemented
  AliOnlineReco& operator=(const AliOnlineReco&); // Not implemented

  // GUI components.
  TGListBox     *fRunList;
  TGTextButton  *fStartButt;
  TGTextButton  *fStopButt;
  TGTextButton  *fXyzzButt;

  // DIM interface. Could do without ...
  AliDimIntNotifier *fSOR;
  AliDimIntNotifier *fEOR;

  // Run-state, process mngmnt
  typedef std::map<Int_t, Int_t> mIntInt_t; // value should be struct { pid, state, ... };
  typedef mIntInt_t::iterator    mIntInt_i;

  mIntInt_t      fRun2PidMap;

  Bool_t         fTestMode;

  mIntInt_i FindMapEntryByPid(Int_t pid);

  ClassDef(AliOnlineReco, 0);
};

#endif
