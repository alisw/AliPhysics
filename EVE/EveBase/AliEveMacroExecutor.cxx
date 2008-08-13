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
// Full description of AliEveMacroExecutor
//

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

    if (mac->GetActive() == kFALSE || mac->GetFunc().IsNull())
    {
      continue;
    }

    switch (mac->GetSources())
    {
      case AliEveMacro::kRunLoader:
	if ( ! AliEveEventManager::HasRunLoader())
	  continue;
	break;
      case AliEveMacro::kESD:
	if ( ! AliEveEventManager::HasESD())
	  continue;
	break;
      case AliEveMacro::kESDfriend:
	if ( ! AliEveEventManager::HasESDfriend())
	  continue;
	break;
      case AliEveMacro::kRawReader:
	if ( ! AliEveEventManager::HasRawReader())
	  continue;
	break;
      default:
	break;
    }

    TString cmd(mac->FormForExec());
    try
    {
      gInterpreter->ProcessLine(cmd);
      // Try to fix broken cint state? Code taken form pyroot.
      if ( G__get_return( 0 ) > G__RETURN_NORMAL )
      {
	printf ("***INFIXING***\n");
	G__security_recover( 0 );    // 0 ensures silence
      }
    }
    catch(TEveException& exc)
    {
      Error("ExecMacros", "Executing %s::%s, caught exception: '%s'.",
	    mac->GetMacro().Data(), cmd.Data(), exc.Data());
    }

    // Try to fix broken cint state? Code taken form pyroot.
    if ( G__get_return( 0 ) > G__RETURN_NORMAL )
    {
      printf ("***POSTFIXING****\n");
      G__security_recover( 0 );    // 0 ensures silence
    }
  }
}
