/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                         *
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

/* $Id$ */
/***************************************************************
 * T0 trigger class for T0 trigger signals:
 *  - T0A
 *  - T0C
 *  - T0vertex
 *  - T0 semi central event for ions
 *  - T0 central            for ions
 ****************************************************************/


#include "AliLog.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliTriggerInput.h"

#include "AliT0.h"
#include "AliT0digit.h"
#include "AliT0Trigger.h"

//----------------------------------------------------------------------
ClassImp(AliT0Trigger)

//----------------------------------------------------------------------
AliT0Trigger::AliT0Trigger()
  : AliTriggerDetector(),
    fT0(0x0),
    fDigits(0x0)
{
   SetName("T0");
   CreateInputs();
}

//----------------------------------------------------------------------
void AliT0Trigger::CreateInputs()
{
   // inputs 
   
   // Do not create inputs again!!
   if( fInputs.GetEntriesFast() > 0 ) return;
   
   fInputs.AddLast( new AliTriggerInput( "T0_A_L0", "T0",  0 ) );
   fInputs.AddLast( new AliTriggerInput( "T0_C_L0", "T0", 0 ) );
   fInputs.AddLast( new AliTriggerInput( "T0_Vertex_L0", "T0", 0 ) );
   fInputs.AddLast( new AliTriggerInput( "T0_Centr_L0", "T0",  0 ) );
   fInputs.AddLast( new AliTriggerInput( "T0_SemiCentral_L0", "T0",  0 ) );
}

//----------------------------------------------------------------------
void AliT0Trigger::Trigger()
{
  // trigger input

   AliRunLoader* runLoader = AliRunLoader::Instance();
   AliLoader * fT0Loader = runLoader->GetLoader("T0Loader");
   //   AliT0digit *fDigits; 
   fT0Loader->LoadDigits("READ");
   // Creating T0 data container

   TTree* treeD = fT0Loader->TreeD();
  if (!treeD) {
    AliError("no digits tree");
    return;
  }
  fDigits = new AliT0digit();

  TBranch *brDigits = treeD->GetBranch("T0");
  if (brDigits) {
    brDigits->SetAddress(&fDigits);
  }else{
    AliError("Branch T0 DIGIT not found");
    return;
  } 
  brDigits->GetEntry(0);
  Int_t   besttimeA = fDigits->BestTimeA();
  Int_t   besttimeC = fDigits->BestTimeC();
  Int_t   timeDiff = fDigits->TimeDiff();
  Int_t    sumMult=   fDigits->SumMult();

  if (besttimeA > 0 && besttimeA <99999  )  SetInput("T0_A_L0");
  if (besttimeC>0  && besttimeC<99999)   SetInput("T0_C_L0"); 
  //6093 corrsponds to vertex -20cm, 6202 vertex +20 with delay 150nc eqalized on the TVDC unit 
  if (timeDiff >6090 && timeDiff < 6210)       SetInput("T0_Vertex_L0");
  if (sumMult > 175)                           SetInput("T0_Centr_L0");
  if (sumMult>155 && sumMult <= 175)           SetInput("T0_SemiCentral_L0");;

   
}
