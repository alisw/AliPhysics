///////////////////////////////////////////////////////////////////////////
///\file AddTaskEMCALNeutralMesons.C
///\brief Configuration of AliAnalysisTaskNeutralMesons
///
/// Neutral mesons analysis in EMCal
///
/// \author Paraskevi Ganoti <Paraskevi.Ganoti@cern.ch>, Athens 
///////////////////////////////////////////////////////////////////////////


AliAnalysisTaskNeutralMesons *AddTaskEMCALNeutralMesons(TString trigger = "MC",
							TString geoname = "EMCAL_COMPLETEV1",
							Int_t mctype = 1,
							Int_t mcpart = 111,
							Double_t lowEnrgCut = 0.3,
							Bool_t ncellsInCl = "kFALSE",
							TString clusterBranch = "newEMCALClusters",
							TString lCustomName=""
							)
{

/// analysis manager

// Get the pointer to the existing analysis manager via the static access method.

  //trigger: MC, EMC, INT7, MB

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskEMCALNeutralMesons", "No analysis manager to connect to.");
    return NULL;
  }
 // Check the analysis type using the event handlers connected to the analysis manager.
 //=====================================================================================
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskEMCALNeutralMesons", "This task requires an input event handler");
    return NULL;
  }

	TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
 	 if(type.Contains("AOD"))
   	  {	 
    	    ::Error("AddTaskKinkpp13TeVMultMC", "This task requires to run on ESD");
   	       return NULL;
    }
    
   AliAnalysisTaskNeutralMesons *task = new AliAnalysisTaskNeutralMesons("taskMesons7TeV");
   task->SetMCtype(mctype);
   task->SetMCprt(mcpart);
   task->SetGeoName(geoname);
   task->SetLowEnergyCut(lowEnrgCut);
   task->SetNcellsPerClusterCut(ncellsInCl);
   task->SetNewClBranchName(clusterBranch);
   
   if(trigger=="EMC7") task->SelectCollisionCandidates(AliVEvent::kEMC7);
   if(trigger=="INT7") task->SelectCollisionCandidates(AliVEvent::kINT7);
   if(trigger=="MB") task->SelectCollisionCandidates(AliVEvent::kMB);

   mgr->AddTask(task);

   AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
   mgr->ConnectInput(task,0,cinput);

  TString outputFileName = Form("%s:PWGGANeutralMesons7TeV", AliAnalysisManager::GetCommonFileName());
  
  TString lContainerName="PWGGANeutralMesons7TeV";
  lContainerName.Append(lCustomName);
  AliAnalysisDataContainer *coutput1= mgr->CreateContainer(lContainerName.Data(),TList::Class(), AliAnalysisManager::kOutputContainer,outputFileName);
  mgr->ConnectOutput(task, 1, coutput1);

   return task; 
}
