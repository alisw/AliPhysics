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

/* $Id: AliEmcalTenderTask.cxx 48689 2011-03-23 21:35:50Z agheata $ */

#include <TChain.h>
#include <TFile.h>

#include "TenderSupplies/AliEMCALTenderSupply.h"
#include "AliAnalysisManager.h"

#include "AliEmcalTenderTask.h"

ClassImp(AliEmcalTenderTask)

//______________________________________________________________________________
AliEmcalTenderTask::AliEmcalTenderTask():
           AliAnalysisTaskSE(),
           fEMCALTender(NULL)
{
// Dummy constructor
}

//______________________________________________________________________________
AliEmcalTenderTask::AliEmcalTenderTask(const char* name):
           AliAnalysisTaskSE(name),
           fEMCALTender(NULL)
{
// Default constructor
  //DefineOutput(1,  AliESDEvent::Class());
}

//______________________________________________________________________________
AliEmcalTenderTask::~AliEmcalTenderTask()
{
// Destructor
  if (fEMCALTender)
    fEMCALTender->Delete();
}

//______________________________________________________________________________
void AliEmcalTenderTask::SetEMCALTenderSupply(AliEMCALTenderSupply *supply)
{
  fEMCALTender = supply;
  supply->SetTask(this);
}
   
//______________________________________________________________________________
void AliEmcalTenderTask::ConnectInputData(Option_t *option)
{
  AliAnalysisTaskSE::ConnectInputData(option);
  fEMCALTender->Init();
}

//______________________________________________________________________________
void AliEmcalTenderTask::UserCreateOutputObjects()
{
// Nothing for the moment, but we may need ESD/AOD event replication here.
 
}

//______________________________________________________________________________
void AliEmcalTenderTask::UserExec(Option_t* option)
{
  fEMCALTender->ProcessEvent();
}
