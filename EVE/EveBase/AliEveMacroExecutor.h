// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveMacroExecutor_H
#define AliEveMacroExecutor_H

#include "TObject.h"

class AliEveMacro;

//______________________________________________________________________________
// Short description of AliEveMacroExecutor
//

class AliEveMacroExecutor : public TObject
{

public:
  AliEveMacroExecutor();
  virtual ~AliEveMacroExecutor();

  void         AddMacro(AliEveMacro* mac);
  AliEveMacro* FindMacro(const TString& func);

  void ExecMacros();
  void RemoveMacros();
//  void SaveAddedMacros();

protected:
  TList*   fMacros;

private:
  AliEveMacroExecutor(const AliEveMacroExecutor&);            // Not implemented
  AliEveMacroExecutor& operator=(const AliEveMacroExecutor&); // Not implemented

  ClassDef(AliEveMacroExecutor, 0); // Container for and executor of AliEveMacros.
};

#endif
