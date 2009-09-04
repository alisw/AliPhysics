// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliChildReaper_H
#define AliChildReaper_H

#include "TObject.h"
#include "TQObject.h"

//______________________________________________________________________________
// Short description of AliChildReaper
//

class AliChildReaper : public TObject,
		       public TQObject
{
public:
  void ChildDeath(Int_t pid, Int_t status); // *SIGNAL*

  static AliChildReaper* Instance();

private:
  AliChildReaper();
  virtual ~AliChildReaper() {}

  AliChildReaper(const AliChildReaper&);            // Not implemented
  AliChildReaper& operator=(const AliChildReaper&); // Not implemented

  static void sig_handler(int sig);

  static AliChildReaper* fgTheOne;

  ClassDef(AliChildReaper, 0);
};

#endif
