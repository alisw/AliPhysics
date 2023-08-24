///////////////////////////////////////////////////////////////////
//                                                               //
// AddTaskUeSpectraRT                                                 //
// Author: Aditya Nath Mishra, Aditya.Nath.Mishra@cern.ch          //
//                                                               //
///////////////////////////////////////////////////////////////////


class AliAnalysisDataContainer;

AliAnalysisTask* AddTaskUeSpectraRT (Bool_t AnalysisMC = kFALSE,
		                        const Char_t* taskname = "UeSpectraRT")
{
  // set authomatically to true if MC
  // Bool_t AnalysisMC= (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
  // set to true if you wish to build corrections (memory consuming)
  Bool_t AnalysisCorr = kTRUE;

  // Get the pointer to the existing analysis manager via the static access method.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if(!mgr){
    ::Error("AddTaskUeSpectraRT", "No analysis manager to connect to.");
      return NULL;
    }

  // set the MC online spectra weights task

  //AliMCSpectraWeights* fMCSpectraWeights = new AliMCSpectraWeights("pp", "fMCSpectraWeights", (AliMCSpectraWeights::SysFlag) 0);
  // TODO: VZ this is still not implemented properly
  //fMCSpectraWeights->SetMCSpectraFile("alien:///alice/cern.ch/user/a/alitrain/PWGLF/LF_pp_MC/1154_20190521-2046/merge/AnalysisResults.root" );
  //fMCSpectraWeights->Init();
  //if (!fMCSpectraWeights) printf("no fMCSpectraWeights object!");//

// Check the analysis type using the event handlers connected to the analysis manager.
//==============================================================================
  if(!mgr->GetInputEventHandler()){
   ::Error("AddTaskUeSpectraRT", "This task requires an input event handler");
      return NULL;
   }


   AliAnalysisTaskUeSpectraRT * task = new AliAnalysisTaskUeSpectraRT("AnalysisUeSpectraRT");
   if(!task) return 0x0;
   TString kInputDataType = mgr->GetInputEventHandler()->GetDataType();

   Bool_t is13TeV = kTRUE;
   UInt_t trigSel = AliVEvent::kINT7;

   task->SetTrigger(trigSel);
   task->SetAnalysisMC(AnalysisMC);
   task->SetAnalysisCorr(AnalysisCorr);
   task->SetAnalysisType(kInputDataType);
   task->SetDebugLevel(0);
   task->SetEtaCut(0.8);
   task->SetVtxCut(10.0);
   task->SetPtLeadMin(5.0);
   task->SetPileUpRej(kTRUE);
   task->SetAveMultiInTrans(4.939);
   task->SetAveRecMultiInTrans(4.895); // reco MC PYTHIA 8
   task->SetAveGenMultiInTrans(7.392); // true MC PYTHIA 8
   //task->SetAveRecMultiInTrans(5.003); // reco MC EPOS LHC
   //task->SetAveGenMultiInTrans(7.72); // true MC EPOS LHC
   //task->SetMCSpectraWeightObject(fMCSpectraWeights);
   mgr->AddTask(task);

   TString fileName = AliAnalysisManager::GetCommonFileName();
   TString kName  = "outputUeSpectraRT";

   mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task,1,mgr->CreateContainer(kName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
   return task;
}
