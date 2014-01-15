/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
/* $Id: AliADTrigger.cxx 49869 2011-05-18 04:49:51Z hristov $ */
// ---------------------
// Class AliADTrigger
// ---------------------
// Top class to simulate the AD trigger response
// This class is only used for interface with AliTriggerDetector
// Its create and Set  Inputs of the CTP
//


#include <TClonesArray.h>
#include <TTree.h>

#include "AliRun.h"
#include "AliLoader.h"
#include "AliLog.h"
#include "AliRunLoader.h"
#include "AliTriggerInput.h"

#include "AliADdigit.h"
#include "AliADTrigger.h"

//______________________________________________________________________
ClassImp(AliADTrigger)

//______________________________________________________________________

AliADTrigger::AliADTrigger():AliTriggerDetector()
{
   SetName("AD");
   CreateInputs();
}
//______________________________________________________________________
void AliADTrigger::CreateInputs()
{
   // Do not create inputs again!!
   if( fInputs.GetEntriesFast() > 0 ) return;

   fInputs.AddLast( new AliTriggerInput( "AD_ADA", "AD", 0 ) );
   fInputs.AddLast( new AliTriggerInput( "AD_ADD","AD", 0 ) );
   fInputs.AddLast( new AliTriggerInput( "AD_ADA2", "AD", 0 ) );
   fInputs.AddLast( new AliTriggerInput( "AD_ADD2","AD", 0 ) );
}

//______________________________________________________________________
void AliADTrigger::Trigger()
{
  //  ********** Get run loader for the current event **********
   AliRunLoader* runLoader = AliRunLoader::Instance();

   AliLoader* loader = runLoader->GetLoader( "ADLoader" );

   if(!loader) {
      AliError("Can not get AD loader");
      return;
   }
   loader->LoadDigits("update");
   TTree* vzeroDigitsTree = loader->TreeD();

   if (!vzeroDigitsTree) {
      AliError("Can not get the AD digit tree");
      return;
   }
   TClonesArray* vzeroDigits = NULL;
   TBranch* digitBranch = vzeroDigitsTree->GetBranch("ADDigit");
   digitBranch->SetAddress(&vzeroDigits);
   
   // Check trigger contitions
   // .... Ex. number of digit over threshold
   //
   
   loader->UnloadDigits();     


   //   if(  )  SetInput( "AD_ADA" );
   //   if(  )  SetInput( "AD_ADD" );

  return;
}



