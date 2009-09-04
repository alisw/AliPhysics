// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliChildProcTerminator_H
#define AliChildProcTerminator_H

#include "TObject.h"
#include "TQObject.h"

//______________________________________________________________________________
// Short description of AliChildProcTerminator
//

class AliChildProcTerminator : public TObject,
		       public TQObject
{
public:
  void ChildProcTerm(Int_t pid, Int_t status); // *SIGNAL*

  static AliChildProcTerminator* Instance();

private:
  AliChildProcTerminator();
  virtual ~AliChildProcTerminator() {}

  AliChildProcTerminator(const AliChildProcTerminator&);            // Not implemented
  AliChildProcTerminator& operator=(const AliChildProcTerminator&); // Not implemented

  static void sig_handler(int sig);

  static AliChildProcTerminator* fgTheOne;

  ClassDef(AliChildProcTerminator, 0);
};

#endif
