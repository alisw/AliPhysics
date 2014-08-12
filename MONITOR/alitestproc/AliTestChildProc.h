// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliTestChildProc_H
#define AliTestChildProc_H

#include <TGFrame.h>

//______________________________________________________________________________
// Short description of AliTestChildProc
//

class AliTestChildProc : public TGMainFrame
{

public:
  AliTestChildProc(Int_t run);

  void DoExit();
  void DoCrash();
  void DoXyzz();
  
private:
  AliTestChildProc(const AliTestChildProc&);            // Not implemented
  AliTestChildProc& operator=(const AliTestChildProc&); // Not implemented

  ClassDef(AliTestChildProc, 0);
};

#endif
