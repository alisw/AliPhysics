// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveMacroExecutor.h"
#include "AliEveMacro.h"
#include "AliEveEventManager.h"

#include <TEveUtil.h>
#include <TList.h>
#include <TROOT.h>

//______________________________________________________________________________
//
// Contains a list of AliEveMacros.
// The macros are added via AddMacro() and are owned by the executor.
// The macros can be executed via ExecMacros().
// They are executed in order in which they are registered.

ClassImp(AliEveMacroExecutor)

//______________________________________________________________________________
AliEveMacroExecutor::AliEveMacroExecutor() :
  TObject(),
  fMacros(new TList)
{
  // Constructor.

  fMacros->SetOwner(kTRUE);
}

//______________________________________________________________________________
AliEveMacroExecutor::~AliEveMacroExecutor()
{
  // Destructor.

  delete fMacros;
}

/******************************************************************************/

void AliEveMacroExecutor::AddMacro(AliEveMacro* mac)
{
  // Add a new macro. Ownership transfered to the executor.

  static const TEveException kEH("AliEveMacroExecutor::AddMacro ");

  const TString mname = mac->GetMacro();
  if ( ! mname.IsNull() && TEveUtil::CheckMacro(mname) == kFALSE)
  {
    TEveUtil::LoadMacro(mname);
  }
  fMacros->Add(mac);
}

AliEveMacro* AliEveMacroExecutor::FindMacro(const TString& func)
{
  // Find macro with given function name (it is supposed to be unique).
  // Returns 0 if not found.

  TIter next(fMacros);
  AliEveMacro* mac;
  while ((mac = (AliEveMacro*) next()))
  {
    if (mac->GetFunc() == func)
      return mac;
  }
  return 0;
}

/******************************************************************************/

#include "Api.h"
#include "TInterpreter.h"

void AliEveMacroExecutor::ExecMacros()
{
  // Execute registered macros.

  TIter next(fMacros);
  AliEveMacro* mac;
  while ((mac = (AliEveMacro*) next()))
  {
    // printf ("macro '%s'; func '%s'; args '%s'\n", mac->GetMacro().Data(), mac->GetFunc().Data(), mac->GetArgs().Data());

    mac->ResetExecState();
    
    if (mac->GetActive() == kFALSE || mac->GetFunc().IsNull())
    {
      continue;
    }

    if ((mac->RequiresRunLoader() && ! AliEveEventManager::HasRunLoader()) ||
        (mac->RequiresESD()       && ! AliEveEventManager::HasESD())       ||
        (mac->RequiresESDfriend() && ! AliEveEventManager::HasESDfriend()) ||
        (mac->RequiresRawReader() && ! AliEveEventManager::HasRawReader()))
    {
      mac->SetExecNoData();
      continue;
    }

    TString cmd(mac->FormForExec());
    try
    {
      Long_t                   result = 0;
      TInterpreter::EErrorCode error  = TInterpreter::kNoError;

      result = gInterpreter->ProcessLine(cmd, &error);

      // Try to fix broken cint state? Code taken form pyroot.
      if (G__get_return(0) > G__RETURN_NORMAL)
      {
	printf ("*** FIXING CINT STATE AFTER RETURN ***\n");
	G__security_recover(0);
      }

      if (error)
      {
        mac->SetExecError();
        Error("ExecMacros", "Executing %s::%s, CINT error ... hopefully recovered.",
              mac->GetMacro().Data(), cmd.Data());
      }
      else
      {
        TEveElement *el  = (TEveElement*) result;
        TObject     *obj = dynamic_cast<TObject*>(el);
        if (el != 0 && obj == 0)
        {
          Warning("ExecMacros", "Executing %s::%s, returned TEveElement seems bad, setting it to 0.",
                  mac->GetMacro().Data(), cmd.Data());
          el = 0;
        }
        mac->SetExecOK(el);
      }
    }
    catch(TEveException& exc)
    {
      mac->SetExecException(exc);

      // Try to fix broken cint state? Code taken form pyroot.
      if (G__get_return(0) > G__RETURN_NORMAL)
      {
	printf ("*** FIXING CINT STATE AFTER EXCEPTION ***\n");
	G__security_recover(0);
      }

      Error("ExecMacros", "Executing %s::%s, caught exception: '%s'.",
	    mac->GetMacro().Data(), cmd.Data(), exc.Data());
    }
  }
}
