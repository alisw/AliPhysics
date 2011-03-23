/**
 * @file   AddTaskCopyHeader.C
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 12:13:43 2011
 * 
 * @brief  
 * 
 * 
 * @ingroup pwg2_forward_scripts_tasks
 */
#include "AliAnalysisTaskSE.h"
#include "AliInputEventHandler.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliCentrality.h"
#include "AliAnalysisManager.h"
#include <TROOT.h>

/**
 * Task to copy header from ESD to AOD 
 * 
 * @ingroup pwg2_forward_scripts_tasks
 * @ingroup pwg2_forward_aod
 */
class CopyHeaderTask : public AliAnalysisTaskSE
{
public:
  CopyHeaderTask(const char* name="header") 
    : AliAnalysisTaskSE(name)
  {}
  CopyHeaderTask(const CopyHeaderTask& other) 
    : AliAnalysisTaskSE(other)
  {}
  virtual ~CopyHeaderTask() {}
  CopyHeaderTask& operator=(const CopyHeaderTask& other) 
  {
    AliAnalysisTaskSE::operator=(other);
    return *this;
  }
  /** 
   * @{ 
   * @name Implementation of interface methods
   */
  virtual void   UserCreateOutputObjects() {}
  virtual void   Init() {}
  virtual void   LocalInit() {Init();}
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *option);
  /* @} */

  ClassDef(CopyHeaderTask,1);
};

void
CopyHeaderTask::UserExec(Option_t*)
{
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
  AliAODEvent* aod = dynamic_cast<AliAODEvent*>(AODEvent());

  if (!esd) { 
    AliWarning("Missing ESD event");
    return;
  }
  if (!aod) { 
    AliWarning("Missing AOD event");
    return;
  }
  
  AliAODHeader* aodHeader = aod->GetHeader();
  if (!aodHeader) { 
    AliWarning("Missing AOD header");
    return;
  }

  aodHeader->SetRunNumber(esd->GetRunNumber());
  aodHeader->SetOfflineTrigger(fInputHandler->IsEventSelected());
  aodHeader->SetBunchCrossNumber(esd->GetBunchCrossNumber());
  aodHeader->SetOrbitNumber(esd->GetOrbitNumber());
  aodHeader->SetPeriodNumber(esd->GetPeriodNumber());
  aodHeader->SetEventType(esd->GetEventType());
  aodHeader->SetEventNumberESDFile(esd->GetHeader()->GetEventNumberInFile());
  if(esd->GetCentrality())
    aodHeader->SetCentrality(new AliCentrality(*(esd->GetCentrality())));
  else
    aodHeader->SetCentrality(0);

  aodHeader->SetFiredTriggerClasses(esd->GetFiredTriggerClasses());
  aodHeader->SetTriggerMask(esd->GetTriggerMask()); 
  aodHeader->SetTriggerCluster(esd->GetTriggerCluster());
  aodHeader->SetL0TriggerInputs(esd->GetHeader()->GetL0TriggerInputs());    
  aodHeader->SetL1TriggerInputs(esd->GetHeader()->GetL1TriggerInputs());    
  aodHeader->SetL2TriggerInputs(esd->GetHeader()->GetL2TriggerInputs());    
  
  aodHeader->SetMagneticField(esd->GetMagneticField());
  aodHeader->SetMuonMagFieldScale(esd->GetCurrentDip()/6000.);
  aodHeader->SetZDCN1Energy(esd->GetZDCN1Energy());
  aodHeader->SetZDCP1Energy(esd->GetZDCP1Energy());
  aodHeader->SetZDCN2Energy(esd->GetZDCN2Energy());
  aodHeader->SetZDCP2Energy(esd->GetZDCP2Energy());
  aodHeader->SetZDCEMEnergy(esd->GetZDCEMEnergy(0),esd->GetZDCEMEnergy(1));
}

void
CopyHeaderTask::Terminate(Option_t*)
{}

/** 
 * 
 * 
 * @ingroup pwg2_forward_aod
 */
void
AddTaskCopyHeader()
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskCopyHeader", "No analysis manager to connect to.");
    return;
  }   
  
  gROOT->Macro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
  
  CopyHeaderTask* task = new CopyHeaderTask;
  mgr->AddTask(task);
  
  // AliAnalysisDataContainer* histOut = 
  //   mgr->CreateContainer("Forward", TList::Class(), 
  // 		 AliAnalysisManager::kOutputContainer,outputfile);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  // mgr->ConnectOutput(task, 1, histOut);
}
//
// EOF
//
