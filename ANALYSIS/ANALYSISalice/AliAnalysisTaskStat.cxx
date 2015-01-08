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

/* $Id$ */
// Author: Andrei Gheata, 20/12/2010

//==============================================================================
// AliAnalysisTaskStat - basic task that attaches a AliAnalysisTaskstatistics 
// object to the analysis manager. Use: AliAnalysisManager::AddStatisticsTask
// to attach to a train.
//==============================================================================

#include "AliAnalysisTaskStat.h"

#include <TList.h>
#include "AliVEvent.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisStatistics.h"

ClassImp(AliAnalysisTaskStat)

//______________________________________________________________________________
AliAnalysisTaskStat::AliAnalysisTaskStat(const char *name)
                    :AliAnalysisTaskSE(name),
                     fStatistics(0),
                     fOutputList(0)
{
// Named constructor.
  DefineOutput(1, TList::Class());
  fBranchNames = "ESD:AliESDHeader. AOD:header";
  fStatistics = new AliAnalysisStatistics("MgrStat");
}

//______________________________________________________________________________
AliAnalysisTaskStat::~AliAnalysisTaskStat()
{
// Destructor.
  if (fOutputList) {
    if (!AliAnalysisManager::GetAnalysisManager()->IsProofMode()) delete fOutputList;
  } else {
    if (fStatistics) delete fStatistics;
  }
} 
    
//______________________________________________________________________________
AliAnalysisTaskStat *AliAnalysisTaskStat::AddToManager(UInt_t offlineMask)
{
// Add this task to the analysis manager. By default it selects MB events.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AliAnalysisTaskStat::AddToManager", "You need a manager first");
    return 0;
  }
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  if (!cinput) {
    ::Error("AliAnalysisTaskStat::AddToManager", "Attach first the input handler");
    return 0;
  }  
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("MgrStat", TList::Class(), AliAnalysisManager::kOutputContainer,
                mgr->GetCommonFileName());
  AliAnalysisTaskStat *taskStatistics = new AliAnalysisTaskStat("MgrStat");
  mgr->AddTask(taskStatistics);
  AliAnalysisStatistics *stat = taskStatistics->GetStatistics();
  stat->SetOfflineMask(offlineMask);
  mgr->SetStatistics(stat);
  taskStatistics->SelectCollisionCandidates(offlineMask);
  mgr->ConnectInput(taskStatistics, 0, cinput);
  mgr->ConnectOutput(taskStatistics, 1, coutput);
  return taskStatistics;
}

//______________________________________________________________________________
void  AliAnalysisTaskStat::UserCreateOutputObjects()
{
// Create the output list.
  if (!fStatistics) {
    Fatal("UserCreateOutputObjects", "You are not allowed to create this task using the dummy constructor. Use the named one.");
  }
  fOutputList = new TList();
  fOutputList->SetOwner();
  if (fStatistics) fOutputList->Add(fStatistics);
  PostData(1, fOutputList);
}   

//______________________________________________________________________________
void  AliAnalysisTaskStat::UserExec(Option_t *)
{
// Event loop.
  fStatistics->AddAccepted();
}

//______________________________________________________________________________
void  AliAnalysisTaskStat::Terminate(Option_t *)
{
// Get the statistics from its container and copy to manager.
  fOutputList = dynamic_cast<TList*>(GetOutputData(1));
  if (!fOutputList) {
    Error("Terminate", "Cannot get output list from container");
    return;
  }
  AliAnalysisStatistics *stat = dynamic_cast<AliAnalysisStatistics*>(fOutputList->At(0));
  if (!stat) {
    Error("Terminate", "Statistics object not found in list");
    return;
  }
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (stat != fStatistics) {
    // Non-local mode
    fStatistics->AddInput(stat->GetNinput());
    fStatistics->AddProcessed(stat->GetNprocessed());
    fStatistics->AddFailed(stat->GetNfailed());
    fStatistics->AddAccepted(stat->GetNaccepted());
    mgr->SetStatistics(fStatistics);
  }
  fStatistics->Print();
}  
