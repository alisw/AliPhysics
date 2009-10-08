// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliDimIntNotifier_H
#define AliDimIntNotifier_H

#include <TObject.h>
#include <TQObject.h>

#include <TMutex.h>
#include <TCondition.h>
#include <TTimer.h>

#ifdef ALI_DIM
#include <dic.hxx>
#else
class DimUpdatedInfo
{
public:
  DimUpdatedInfo(const Char_t*, Int_t) {}

  Bool_t getData() { return kFALSE; }
  Int_t  getInt()  { return -1; }
};
#endif

//______________________________________________________________________________
// Short description of AliDimIntNotifier
//

class AliDimIntNotifier : public TObject,
			  public TQObject,
			  public DimUpdatedInfo
{
public:

  AliDimIntNotifier(const TString& service);

  virtual ~AliDimIntNotifier() {}

  void infoHandler();
  void infoHandlerTest(Int_t fake);

  void DimMessage(Int_t=-1); // *SIGNAL*

  static void SetMainThreadId();

private:
  AliDimIntNotifier(const AliDimIntNotifier&);            // Not implemented
  AliDimIntNotifier& operator=(const AliDimIntNotifier&); // Not implemented

  void StartTimer();
  void StopTimer();

  TTimer fReThreader;
  TMutex fNotifyLck;

  Int_t  fLastMessage;

  static Long_t fgMainThreadId;

  ClassDef(AliDimIntNotifier, 0);
};

#endif
