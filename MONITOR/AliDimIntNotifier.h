// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliDimIntNotifier_H
#define AliDimIntNotifier_H

#include <TQObject.h>

#ifdef ALI_DIM
#include <dic.hxx>
#else
class DimUpdatedInfo
{
public:
  DimUpdatedInfo(const Char_t*, Int_t) {}
  virtual ~DimUpdatedInfo() {}

  Bool_t getData() { return kFALSE; }
  Int_t  getInt()  { return -1; }
};
#endif

//______________________________________________________________________________
// Short description of AliDimIntNotifier
//

class AliDimIntNotifier :  public TQObject,
			  public DimUpdatedInfo
{
public:

  AliDimIntNotifier(const TString& service);

  virtual ~AliDimIntNotifier() {}

  void infoHandler();
  void DimMessage(Int_t mess=-1); // *SIGNAL*

private:
  AliDimIntNotifier(const AliDimIntNotifier&);            // Not implemented
  AliDimIntNotifier& operator=(const AliDimIntNotifier&); // Not implemented

  Int_t  fLastMessage;

  ClassDef(AliDimIntNotifier, 0);
};

#endif
